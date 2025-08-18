
use anyhow::{anyhow, bail, ensure, Context};
use indicatif::ParallelProgressIterator;
use log::{debug, info, trace};
use noodles::core::region::Interval;
use noodles::core::Region;
use noodles::vcf;
use noodles::vcf::variant::record::samples::keys::key as vcf_key;
use noodles_util::variant::io::IndexedReader as VcfReader;
use noodles_util::variant::io::indexed_reader::Builder as VcfBuilder;
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use std::collections::{BTreeMap, HashMap, VecDeque};
use std::fs::File;
use std::path::{Path, PathBuf};
use rayon::prelude::*;

use crate::data_types::coordinates::Coordinates;
use crate::data_types::multi_region::MultiRegion;
use crate::data_types::phase_enums::PhasedZygosity;
use crate::data_types::variants::{Variant, VariantType};
use crate::parsing::noodles_helper::LoadedBed;
use crate::util::progress_bar::get_progress_style;

type VcfChromKey = (usize, usize);
pub struct RegionIterator {
    /// Tracks the next block index
    next_region_id: u64,
    /// The file paths to all input VCFs
    vcf_paths: Vec<PathBuf>,
    /// Index of the samples in the VCFs
    vcf_sample_index: Vec<usize>,
    /// The current chromosome index
    chrom_index: usize,
    /// Reader for the high confidence BED file
    hc_bed_regions: LoadedBed,
    /// Contains regions we parsed out of the previous HC region
    region_queue: std::collections::VecDeque<MultiRegion>,
    /// Amount of flanks to apply to each variant
    flank_size: usize,
    /// Allows variants to get trimmed, removing extra bases from REF/ALT
    enable_trimming: bool,
    /// Tracks the chromosome lengths
    chrom_lengths: HashMap<String, usize>,
    /// If true, changes output messages slightly
    is_compare_iterator: bool,
    /// Lookup from (vcf_index, chrom_index) to the list of variants with phases
    preloaded_variants: Option<BTreeMap<VcfChromKey, Vec<(Variant, PhasedZygosity)>>>
}

impl RegionIterator {
    /// Creates a region iterator from the provided input files.
    /// # Arguments
    /// * `truth_vcf_fn` - filepath for the truth VCF, multiple formats supported
    /// * `truth_sample` - sample name to read from in the truth VCF
    /// * `query_vcf_fn` - filepath for the query VCF, multiple formats supported
    /// * `query_sample` - sample name to read from in the query VCF
    /// * `confidence_regions` - filepath for the confidence region, BED expected (Optional)
    /// * `reference_genome` - the pre-loaded reference genome dictionary
    /// * `flank_size` - the window to merge variants
    /// * `enable_trimming` - allows variants to get trimmed, removing extra bases from REF/ALT
    #[allow(clippy::too_many_arguments)]
    pub fn new_compare_iterator(
        truth_vcf_fn: &Path,
        truth_sample: &str,
        query_vcf_fn: &Path,
        query_sample: &str,
        confidence_regions: Option<&Path>,
        reference_genome: &ReferenceGenome,
        flank_size: usize,
        enable_trimming: bool
    ) -> anyhow::Result<Self> {
        // Open the VCF files
        let mut truth_vcf_reader = VcfBuilder::default()
            .build_from_path(truth_vcf_fn)
            .with_context(|| format!("Error while opening {truth_vcf_fn:?} (or associated index):"))?;
        let mut query_vcf_reader = VcfBuilder::default()
            .build_from_path(query_vcf_fn)
            .with_context(|| format!("Error while opening {query_vcf_fn:?} (or associated index):"))?;

        // get the headers also
        let truth_vcf_header = truth_vcf_reader.read_header()
            .with_context(|| format!("Error while reading header of {truth_vcf_fn:?}:"))?;
        let query_vcf_header = query_vcf_reader.read_header()
            .with_context(|| format!("Error while reading header of {query_vcf_fn:?}:"))?;

        // make sure the provided sample names are in the VCF files
        let truth_index = truth_vcf_header.sample_names().get_index_of(truth_sample)
            .ok_or(anyhow!("Sample name {truth_sample:?} was not found in {truth_vcf_fn:?}"))?;
        let query_index = query_vcf_header.sample_names().get_index_of(query_sample)
            .ok_or(anyhow!("Sample name {query_sample:?} was not found in {query_vcf_fn:?}"))?;

        let hc_bed_regions = if let Some(cr_fn) = confidence_regions {
            LoadedBed::preload_bed_file(cr_fn)?
        } else {
            bail!("High confidence regions are currently required.");
        };

        let chrom_lengths = reference_genome.contig_keys().iter()
            .map(|k| (k.clone(), reference_genome.get_full_chromosome(k).len()))
            .collect();

        // TODO: do we want to add a check here verifying non-overlapping intervals?
        //       right now, it does an implicit prune since we send each variant to at most one interval

        // everything opened fine at least
        Ok(Self {
            next_region_id: 0,
            vcf_paths: vec![truth_vcf_fn.to_path_buf(), query_vcf_fn.to_path_buf()],
            vcf_sample_index: vec![truth_index, query_index],
            chrom_index: 0, // start with w/e the first chromosome is
            hc_bed_regions,
            region_queue: Default::default(),
            flank_size,
            enable_trimming,
            chrom_lengths,
            is_compare_iterator: true,
            preloaded_variants: None
        })
    }

    /// Creates a region iterator from the provided input files.
    /// # Arguments
    /// * `vcf_filenames` - filepath for the input VCFs, multiple formats supported
    /// * `sample_names` - sample name to read from in the provided VCFs
    /// * `confidence_regions` - filepath for the confidence region, BED expected (Optional)
    /// * `reference_genome` - the pre-loaded reference genome dictionary
    /// * `flank_size` - the window to merge variants
    /// * `enable_trimming` - allows variants to get trimmed, removing extra bases from REF/ALT
    pub fn new_merge_iterator<P: AsRef<Path> + std::fmt::Debug, S: AsRef<str> + std::fmt::Debug>(
        vcf_filenames: &[P],
        sample_names: &[S],
        confidence_regions: Option<&Path>,
        reference_genome: &ReferenceGenome,
        flank_size: usize,
        enable_trimming: bool
    ) -> anyhow::Result<Self> {
        // sanity checks
        ensure!(vcf_filenames.len() == sample_names.len(), "vcf_filenames and sample_names must be equal length");
        ensure!(!vcf_filenames.is_empty(), "Must provide at least 1 VCF to iterate on");

        // Prep our outputs
        let mut vcf_paths = vec![];
        let mut vcf_sample_index = vec![];

        // we want to open each file and do some sanity checks before going into iteration mode
        for (p, s) in vcf_filenames.iter().zip(sample_names.iter()) {
            // open the VCF reader
            let mut vr = VcfBuilder::default()
                .build_from_path(p)
                .with_context(|| format!("Error while opening {p:?} (or associated index):"))?;

            // parse the header
            let vh = vr.read_header()
                .with_context(|| format!("Error while reading header of {p:?}:"))?;

            // make sure we have the sample of interest in the header
            let si = vh.sample_names().get_index_of(s.as_ref())
                .ok_or(anyhow!("Sample name {s:?} was not found in {p:?}"))?;

            // just save the VCF and the index in the VCF
            vcf_paths.push(p.as_ref().to_path_buf());
            vcf_sample_index.push(si);
        }

        let hc_bed_regions = if let Some(cr_fn) = confidence_regions {
            LoadedBed::preload_bed_file(cr_fn)?
        } else {
            bail!("High confidence regions are currently required.");
        };

        let chrom_lengths = reference_genome.contig_keys().iter()
            .map(|k| (k.clone(), reference_genome.get_full_chromosome(k).len()))
            .collect();

        // TODO: do we want to add a check here verifying non-overlapping intervals?
        //       right now, it does an implicit prune since we send each variant to at most one interval

        // everything opened fine at least
        Ok(Self {
            next_region_id: 0,
            vcf_paths,
            vcf_sample_index,
            chrom_index: 0, // start with w/e the first chromosome is
            hc_bed_regions,
            region_queue: Default::default(),
            flank_size,
            enable_trimming,
            chrom_lengths,
            is_compare_iterator: false,
            preloaded_variants: None
        })
    }

    /// Triggers the region iterator to pre-load all of the variants into memory in parallel.
    /// This is most useful in the typical run, where each chromosome and file can be loaded independently, and then mixed later.
    /// # Errors
    /// * if there are any VCF parsing errors
    /// * if there are any variant/zygosity creation issues
    pub fn preload_all_variants(&mut self) -> anyhow::Result<()> {
        ensure!(self.chrom_index == 0, "Can only call preload if no chromosomes have been parsed");
        ensure!(self.preloaded_variants.is_none(), "Cannot preload variants twice");

        // figure out what the loads are, the key here is (vcf_index, chrom_index)
        let mut load_keys: Vec<VcfChromKey> = vec![];
        let mut full_chrom_regions: Vec<Region> = vec![];
        let mut chrom_index = 0;
        loop {
            // check if there are more chromosomes
            let opt_next_chrom = self.hc_bed_regions.get_index(chrom_index);
            if opt_next_chrom.is_none() {
                break;
            }

            // we have more, add each
            let (chrom, intervals) = opt_next_chrom.unwrap();
            ensure!(self.chrom_lengths.contains_key(chrom), "Chromosome {chrom} was not found in reference genome");

            // add a load entry for each VCF
            for i in 0..self.vcf_paths.len() {
                load_keys.push((i, chrom_index));
            }

            // add an entry for the chromosome
            let full_start = intervals[0].start().unwrap();
            let full_end = intervals.last().unwrap().end().unwrap();

            // build the region noodles will recognize
            let chrom_region = Region::new(
                chrom.clone(),
                full_start..=full_end // Position here is 1-based
            );
            full_chrom_regions.push(chrom_region);

            // increment so we get the next chromosome in the next loop
            chrom_index += 1;
        }

        // sort them so they are clustered by VCF index
        load_keys.sort();

        // now pre-load all the variants in parallel
        let style = get_progress_style();
        self.preloaded_variants = Some(load_keys.into_par_iter()
            .progress_with_style(style)
            .map(|(vi, chrom_index)| {
                // get the full chrom region
                let full_chrom_region = &full_chrom_regions[chrom_index];

                // open the reader, get the header
                let vcf_path = self.vcf_paths[vi].as_path();
                let mut vcf_reader = VcfBuilder::default()
                    .build_from_path(vcf_path)
                    .with_context(|| format!("Error while opening {vcf_path:?} (or associated index):"))?;
                let vcf_header = vcf_reader.read_header()
                    .with_context(|| format!("Error while reading header of {vcf_path:?}:"))?;

                // we pre-identified the sample index
                let vcf_index = self.vcf_sample_index[vi];

                // load the variants into memory
                match load_variants_in_region(
                    &vcf_header, &mut vcf_reader, vcf_index, full_chrom_region, self.enable_trimming
                ).with_context(|| format!("Error while pre-loading variants from input #{vi} in {full_chrom_region}:")) {
                    Ok(v) => Ok(
                        ((vi, chrom_index), v)
                    ),
                    Err(e) => Err(e)
                }
            })
            .collect::<anyhow::Result<_>>()?);

        Ok(())
    }
}

impl Iterator for RegionIterator {
    type Item = anyhow::Result<MultiRegion>;

    /// This iterator works by loading chunks of intervals at a time.
    /// Specifically, it parses all intervals on a chromosome, which is *much* faster because we only do one tabix lookup per chrom.
    fn next(&mut self) -> Option<Self::Item> {
        while self.region_queue.is_empty() {
            let opt_next_chrom = self.hc_bed_regions.get_index(self.chrom_index);
            let (chrom, intervals) = opt_next_chrom?;
            let chrom_length = match self.chrom_lengths.get(chrom) {
                Some(&cl) => cl,
                None => {
                    return Some(Err(anyhow!("Chromosome {chrom} was not found in reference genome")));
                }
            };
            self.chrom_index += 1;

            let num_input_vcfs = self.vcf_paths.len();
            let full_start = intervals[0].start().unwrap();
            let full_end = intervals.last().unwrap().end().unwrap();

            // build the region noodles will recognize
            let full_chrom_region = Region::new(
                chrom.clone(),
                full_start..=full_end // Position here is 1-based
            );

            let joint_results: Vec<(usize, Vec<(Variant, PhasedZygosity)>)> = if let Some(preloaded_variants) = self.preloaded_variants.as_mut() {
                // we have pre-loaded variants, pop them off of our preloaded Map and into this Vec for sorting
                debug!("Collecting pre-loaded variants in {full_chrom_region}...");
                let mut jr = vec![];
                let chrom_index = self.chrom_index - 1; // already incremented above, so decrement here
                for vcf_index in 0..self.vcf_paths.len() {
                    jr.push(
                        (vcf_index, preloaded_variants.remove(&(vcf_index, chrom_index)).unwrap())
                    );
                }
                jr
            } else {
                // we did not pre-load, so the variants get loaded in parallel here
                // in general, this is slower because we are only parallelizing by file
                debug!("Loading all variants in {full_chrom_region}...");

                // this is a parallelized variant loading loop
                // we found that we can pass around the VCF file paths, but the VCF readers themselves do not Send
                match self.vcf_paths.par_iter()
                    .enumerate()
                    .map(|(i, vcf_path)| {
                        // open the reader, get the header
                        let mut vcf_reader = VcfBuilder::default()
                            .build_from_path(vcf_path)
                            .with_context(|| format!("Error while opening {vcf_path:?} (or associated index):"))?;
                        let vcf_header = vcf_reader.read_header()
                            .with_context(|| format!("Error while reading header of {vcf_path:?}:"))?;

                        // we pre-identified the sample index
                        let vcf_index = self.vcf_sample_index[i];

                        // load the variants into memory
                        match load_variants_in_region(
                            &vcf_header, &mut vcf_reader, vcf_index, &full_chrom_region, self.enable_trimming
                        ).with_context(|| format!("Error while parsing variants from input #{i} in {full_chrom_region}:")) {
                            Ok(v) => Ok((i, v)),
                            Err(e) => Err(e)
                        }
                    })
                    .collect() {
                        Ok(v) => v,
                        Err(e) => {
                            // translate the error into an option for main
                            return Some(Err(e));
                        }
                    }
            };

            // variants were returned with an associated index, we want to push that onto each variant for sorting
            let mut joint_vec: Vec<(usize, Variant, PhasedZygosity)> = vec![];
            for (i, pvariants) in joint_results.into_iter() {
                // output the number of variants found while we go
                if self.is_compare_iterator {
                    if i == 0 {
                        info!("Loaded {} truth variants on {chrom}.", pvariants.len());
                    } else {
                        info!("Loaded {} query variants on {chrom}.", pvariants.len());
                    }
                } else {
                    info!("Loaded {} variants from input #{i} on {chrom}.", pvariants.len());
                }

                // extend our joint loop with the vcf index + variant + zygosity
                joint_vec.extend(
                    pvariants.into_iter()
                        .map(|(v, p)| (i, v, p))
                );
            }

            // ...which we sort by position
            joint_vec.sort_by_key(|(_b, v, _p)| v.position());

            // ... and then convert into groups of variants for problem solving
            // convert to deque and iterate
            let mut joint_deque: VecDeque<(usize, Variant, PhasedZygosity)> = joint_vec.into();
            for interval in intervals.iter() {
                // TODO: this tracking seems a little messy to me, but I'm not sure we can do it better
                //       maybe we push it into logic inside CompareRegion, something similar to a builder; not sure that's any better
                let mut variants = vec![vec![]; num_input_vcfs];
                let mut zygosities = vec![vec![]; num_input_vcfs];
                let mut window_start = None;
                let mut window_end = None;

                // iterate until we clear out the deque
                while let Some((vcf_index, variant, zygosity)) = joint_deque.pop_front() {
                    match get_variant_containment(interval, &variant) {
                        Containment::Before => {
                            // this variant is BEFORE, so just skip it
                            // no-op
                        },
                        Containment::Contained => {
                            // contained, add it to this tracking set
                            // check if the position of this variant should be within the current block, or if it should start a new one
                            let pos = variant.position() as usize;
                            if pos >= window_end.unwrap_or(usize::MAX) {
                                // we need to create a new block, this one is too far away to be included in the current
                                let coordinates = Coordinates::new(
                                    chrom.clone(), window_start.unwrap() as u64, window_end.unwrap() as u64
                                );
                                let cr = match MultiRegion::new(
                                    self.next_region_id, coordinates, variants, zygosities
                                ) {
                                    Ok(cr) => cr,
                                    Err(e) => return Some(Err(e))
                                };
                                self.region_queue.push_back(cr);
                                self.next_region_id += 1;

                                // now recent all our sentinels also
                                variants = vec![vec![]; num_input_vcfs];
                                zygosities = vec![vec![]; num_input_vcfs];
                                window_start = None;
                                // window_end = None;
                            }

                            // now do any updates to the current block (which may be empty/new)
                            if window_start.is_none() {
                                // this must be the first variant, set the window start
                                window_start = Some((variant.position() as usize).saturating_sub(self.flank_size));
                            }
                            // adjust the end overlaps
                            let ref_len = variant.ref_len();
                            let var_flank_end = (pos + ref_len + self.flank_size).min(chrom_length);
                            window_end = Some(match window_end {
                                Some(current_end) => current_end.max(var_flank_end), // it is possible that an earlier variant is longer than this one
                                None => var_flank_end
                            });

                            // save the actual variant/zygosity
                            variants[vcf_index].push(variant);
                            zygosities[vcf_index].push(zygosity);
                        },
                        Containment::After => {
                            // this means we reach the end of variants within this interval
                            // we need to put this variant back in front in case the next interval wants it
                            joint_deque.push_front((vcf_index, variant, zygosity));
                            break;
                        },
                        Containment::Overlapping => {
                            // we are assuming that the BED regions are non-overlapping
                            // so if it partially overlaps this BED region, then it cannot be used elsewhere
                            // no-op
                        },
                    };
                }

                // check if we have a block to save
                if let (Some(ws), Some(we)) = (window_start, window_end) {
                    // make sure we have a variant also
                    assert!(variants.iter().any(|v| !v.is_empty()));

                    // we need to create a new block, this one is too far away to be included in the current
                    let coordinates = Coordinates::new(
                        chrom.clone(), ws as u64, we as u64
                    );

                    let cr = match MultiRegion::new(
                        self.next_region_id, coordinates, variants, zygosities
                    ) {
                        Ok(cr) => cr,
                        Err(e) => return Some(Err(e))
                    };
                    self.region_queue.push_back(cr);
                    self.next_region_id += 1;
                } else {
                    // empty interval, nothing to do here
                }
            }

            info!("Found {} segments on {chrom}.", self.region_queue.len());
        }

        // we have something in the queue to pop
        let front_region = self.region_queue.pop_front().unwrap();
        Some(Ok(front_region))
    }
}

/// This will load all variants in a given region into a Deque, order is the same as in the VCF file.
/// Variants are pre-parsed into the `Variant` and `PhasedZygosity` types.
/// # Arguments
/// * `vcf_header` - pre-loaded VCF header from noodles
/// * `vcf_reader` - dynamic typed index VCF reader
/// * `sample_index` - index of the sample in the VCF file, usually 0 in our case
/// * `region` - the full region to load variants from
/// * `enable_trimming` - Allows variants to get trimmed, removing extra bases from REF/ALT
fn load_variants_in_region(
    vcf_header: &vcf::Header,
    vcf_reader: &mut VcfReader<noodles::bgzf::io::Reader<File>>,
    sample_index: usize,
    region: &Region,
    enable_trimming: bool
) -> anyhow::Result<Vec<(Variant, PhasedZygosity)>> {
    // first, check if the sequence name is in the index
    // this is does outside the Err(e) because of borrow checking
    let index = vcf_reader.index();
    let seqname_absent = if let Some(index_header) = index.header() {
        let seq_names = index_header.reference_sequence_names();
        !seq_names.contains(region.name())
    } else {
        false
    };

    // now try to run the query
    let mut ret: Vec<(Variant, PhasedZygosity)> = Default::default();
    let query = match vcf_reader.query(vcf_header, region) {
        Ok(q) => q,
        Err(e) => {
            // check if this is an absent sequence name in the tabix index
            if seqname_absent {
                // the index does not have the seq name, so it is likely empty
                trace!("Error while fetching region, likely empty: {e:?}");
                return Ok(Default::default());
            } else {
                // seq name is present, this is likely some other error to propagate upstream
                return Err(e.into());
            }
        }
    };

    for result in query {
        let record: Box<dyn vcf::variant::Record> = result?;
        let record_buf = vcf::variant::RecordBuf::try_from_variant_record(vcf_header, record.as_ref())?;
    
        let variants = parse_variant(&record_buf, sample_index, enable_trimming)
            .with_context(|| {
                format!("Error parsing variants in {record_buf:?}:")
            })?;

        // add them to the list
        trace!("\tFound {variants:?}");

        for (variant, zygosity) in variants.into_iter() {
            if is_variant_contained(region, &variant) {
                ret.push((variant, zygosity));
            }
        }
    }

    Ok(ret)
}

/// Given a pre-parsed variant record, this will convert it into `Variant` and `PhasedZygosity` types for a sample.
/// Note that this method currently split multi-ALT calls into multiple entries (e.g., 2|1 -> 2|0 and 0|1).
/// # Arguments
/// * `record` - the record to parse
/// * `sample_index` - index of the sample to pull genotypes from
/// * `enable_trimming` - Allows variants to get trimmed, removing extra bases from REF/ALT
fn parse_variant(
    record: &vcf::variant::RecordBuf,
    sample_index: usize,
    enable_trimming: bool
) -> anyhow::Result<Vec<(Variant, PhasedZygosity)>> {
    // variant level column
    let chrom = record.reference_sequence_name();
    let pos = record.variant_start().ok_or(anyhow!("Missing POS"))?; // 1-based
    let ref_seq = record.reference_bases();
    let alts = record.alternate_bases().as_ref();
    
    // sample specific information
    let all_samples = record.samples();
    let sample = all_samples.get_index(sample_index).unwrap();
    let gt = match sample.get(vcf_key::GENOTYPE)
        .ok_or(anyhow!("Missing GT"))? {
        Some(gt) => gt,
        None => {
            // this entry's GT = '.', which is a no-op for us
            return Ok(Default::default());
        }
    };

    trace!("{chrom}\t{pos}\t{ref_seq:?}\t{alts:?}\tGT={gt:?}");

    // parse the genotype
    let genotypes = parse_genotype(gt)?;
    
    // go through the genotypes we found and convert them into Variants
    let mut ret = vec![];
    for (alt_index, zygosity) in genotypes.into_iter() {
        if alts[alt_index-1].as_bytes() == b"*" {
            // this is effectively a reference allele
            continue;
        }

        let mut ref_sequence = ref_seq.as_bytes().to_vec();
        let mut alt_sequence = alts[alt_index-1].as_bytes().to_vec();
        if alt_sequence[0] == b'<' {
            // e.g. <INV>; we need sequence-resolved
            continue;
        }
        let position = (pos.get() - 1) as u64; // convert to 0-based

        // trim off any extra bases at the end
        while enable_trimming && ref_sequence.len() > 1 && alt_sequence.len() > 1 && ref_sequence.last().unwrap() == alt_sequence.last().unwrap() {
            assert!(ref_sequence.pop().is_some());
            assert!(alt_sequence.pop().is_some());
        }

        // TODO: do we want to parameterize this?
        let allele_size_limit = 10000;
        if ref_sequence.len() > allele_size_limit || alt_sequence.len() > allele_size_limit {
            debug!("Ignoring variant at {chrom}:{pos} {}>{}", ref_sequence.len(), alt_sequence.len());
            // TODO: do we care to track the number of variant ignored due to length?
            continue;
        }

        let variant_type = get_variant_type(record, &ref_sequence, &alt_sequence)?;
        let variant = match variant_type {
            VariantType::Snv => Variant::new_snv(0, position, ref_sequence, alt_sequence),
            VariantType::Insertion => Variant::new_insertion(0, position, ref_sequence, alt_sequence),
            VariantType::Deletion => Variant::new_deletion(0, position, ref_sequence, alt_sequence),
            VariantType::Indel => Variant::new_indel(0, position, ref_sequence, alt_sequence),
            VariantType::TrContraction => Variant::new_tr_contraction(0, position, ref_sequence, alt_sequence),
            VariantType::TrExpansion => Variant::new_tr_expansion(0, position, ref_sequence, alt_sequence),
            VariantType::SvDeletion => Variant::new_sv_deletion(0, position, ref_sequence, alt_sequence),
            VariantType::SvInsertion => Variant::new_sv_insertion(0, position, ref_sequence, alt_sequence),
            // here are the ones we explicitly are not supporting at this time
            VariantType::SvBreakend |
            VariantType::SvDuplication => {
                continue;
            },
            // Whenever you add a new VariantType, update `SUPPORTED_VARIANT_TYPES` in waffle_solver.rs
            _ => bail!("No impl for {variant_type:?}")
        }?;

        ret.push((variant, zygosity));
    }
    Ok(ret)
}

/// Parses the GT field of a record and returns a Vec of ALT allele index with genotype.
/// Multiple returns values are possible if multiple ALTs are part of the genotype.
/// # Arguments
/// * `gt` - the GT field from the record
fn parse_genotype(gt: &vcf::variant::record_buf::samples::sample::Value) -> anyhow::Result<Vec<(usize, PhasedZygosity)>> {
    use vcf::variant::record::samples::series::value::genotype::Phasing;

    let mut ret = vec![];
    if let vcf::variant::record_buf::samples::sample::Value::Genotype(genotype) = gt {
        let alleles = genotype.as_ref();

        // figure out the indices
        let (a1, a2) = match alleles.len() {
            // treats hemi-zygous as homozygous
            // TODO: do we want a special hemizygous state? how do we score that?
            1 => (alleles[0].position(), alleles[0].position()),
            2 => (alleles[0].position(), alleles[1].position()),
            _ => bail!("allele.len() != [1, 2]: {gt:?}")
        };

        // if it's None ('.'), we treat it as a reference ('0') call
        let i1 = a1.unwrap_or(0);
        let i2 = a2.unwrap_or(0);

        // both have an index set
        if i1 == i2 {
            // homozygous path
            if i1 == 0 {
                // homozygous reference
            } else {
                // homozygous alternate
                ret.push((i1, PhasedZygosity::HomozygousAlternate))
            }
        } else {
            // heterozygous, figure out if they're phased
            let p1 = alleles[0].phasing();
            let p2 = alleles[1].phasing();
            let is_phased = p1 == Phasing::Phased || p2 == Phasing::Phased;

            // now build the phase enum for the variant
            let ap1 = if is_phased { PhasedZygosity::PhasedHet10 } else { PhasedZygosity::UnphasedHeterozygous };
            let ap2 = if is_phased { PhasedZygosity::PhasedHet01 } else { PhasedZygosity::UnphasedHeterozygous };

            // save them if non-reference
            if i1 != 0 {
                ret.push((i1, ap1));
            }
            if i2 != 0 {
                ret.push((i2, ap2));
            }
        }
    } else {
        // it's not a Genotype; should we throw an error?
        bail!("Non-genotype value provided to parse_genotype: {gt:?}");
    }
    Ok(ret)
}

/// Given a record, this will extract the type of the variant contained
/// # Arguments
/// * `record` - the original loaded record, which may have additional tags to parse
/// * `ref_sequence` - the reference sequence, which may be adjusted from what is in the record
/// * `ref_sequence` - the alt sequence, which may be adjusted from what is in the record
fn get_variant_type(record: &vcf::variant::RecordBuf, ref_sequence: &[u8], alt_sequence: &[u8]) -> anyhow::Result<VariantType> {
    use vcf::variant::record::info::field::key as info_key;

    // check for SVs, some of which are supported
    let opt_sv_type = record.info().get(info_key::SV_TYPE);
    if let Some(Some(sv_type)) = opt_sv_type {
        let vt = match sv_type {
            vcf::variant::record_buf::info::field::Value::String(s) => {
                match s.as_str() {
                    "BND" => VariantType::SvBreakend,
                    "DEL" => VariantType::SvDeletion,
                    "DUP" => VariantType::SvDuplication,
                    "INS" => VariantType::SvInsertion,
                    _ => bail!("Unsupported SVTYPE detected: {sv_type:?}")
                }
            },
            _ => bail!("Unsupported SVTYPE detected: {sv_type:?}")
        };
        return Ok(vt);
    }

    // check for STRs, which are split into TrContraction and TrExpansion base on the length of the REF/ALT sequences
    let opt_trid = record.info().get("TRID");
    let vt = if let Some(Some(_trid)) = opt_trid {
        if alt_sequence.len() < ref_sequence.len() {
            VariantType::TrContraction
        } else {
            VariantType::TrExpansion
        }
    } else {
        match (ref_sequence.len(), alt_sequence.len()) {
            (0, _) | (_, 0) => bail!("cannot have alleles with 0 length"),
            (1, 1) => VariantType::Snv,
            (1, _) => VariantType::Insertion,
            (_, 1) => VariantType::Deletion,
            (_, _) => VariantType::Indel
        }
    };
    Ok(vt)
}

/// Checks if the given region _fully_ contains the described variant.
/// # Arguments
/// * `region` - a defined BED region
/// * `variant` - the variant to check
fn is_variant_contained(region: &Region, variant: &Variant) -> bool {
    // get the region interval, which is 1-based and convert to 0-based
    let interval = region.interval();
    let zb_start = interval.start().unwrap().get() - 1;
    let zb_end = interval.end().unwrap().get();
    let zb_interval = zb_start..zb_end;

    let variant_start = variant.position() as usize;
    let last_contained_pos = variant_start + variant.ref_len() - 1; // subtract one since this is exclusive range

    // if it contains both the first and last base position, then we're golden
    zb_interval.contains(&variant_start) &&
        zb_interval.contains(&last_contained_pos)
}

enum Containment {
    /// Indicates the variant start is before the region start
    Before,
    /// Indicates the variant is fully contained within the region
    Contained,
    /// Indicates the variant starts inside the region, but ends outside the region
    Overlapping,
    /// Indicates the variant is fully after the region
    After
}

/// Checks if the given region _fully_ contains the described variant.
/// # Arguments
/// * `interval` - an interval from a BED region
/// * `variant` - the variant to check
fn get_variant_containment(interval: &Interval, variant: &Variant) -> Containment {
    // get the region interval, which is 1-based and convert to 0-based
    let zb_start = interval.start().unwrap().get() - 1;
    let zb_end = interval.end().unwrap().get();

    let variant_start = variant.position() as usize;
    let variant_end = variant_start + variant.ref_len();

    if variant_start < zb_start {
        Containment::Before
    } else if variant_start >= zb_end {
        // the start is after the end, so definitely full after
        Containment::After
    } else if variant_end <= zb_end {
        Containment::Contained
    } else {
        Containment::Overlapping
    }
}

#[cfg(test)]
mod tests {
    // use super::*;

    /*
    TODO: we likely need to add tests for all the various parsing functions above
    I think this will be easier to just make a mock VCF and then check the results here, we can do this later though; I'm less concerned with that right now
     */
}