
use std::path::{Path, PathBuf};

/// Wrapper that indexes a VCF file using noodles indexer, generating a .tbi file.
/// # Arguments
/// * `vcf_fn` - the filename to index
pub fn index_vcf(vcf_fn: &Path) -> anyhow::Result<()> {
    // first, build the index
    let index = noodles::vcf::fs::index(vcf_fn)?;

    // add the .tbi extension to the filename
    let mut tbi_fn = vcf_fn.to_owned().into_os_string();
    tbi_fn.push(".tbi");
    let tbi_fn = PathBuf::from(tbi_fn);
    
    // write the index out to file
    noodles::tabix::fs::write(&tbi_fn, &index)?;

    Ok(())
}

/// Wrapper that indexes a BED file using noodles indexer, generating a .tbi file.
/// # Arguments
/// * `bed_fn` - the filename to index
pub fn index_bed(bed_fn: &Path) -> anyhow::Result<()> {
    // first, build the index
    let index = noodles::bed::fs::index(bed_fn)?;

    // add the .tbi extension to the filename
    let mut tbi_fn = bed_fn.to_owned().into_os_string();
    tbi_fn.push(".tbi");
    let tbi_fn = PathBuf::from(tbi_fn);

    // write the index out to file
    noodles::tabix::fs::write(&tbi_fn, &index)?;

    Ok(())
}
