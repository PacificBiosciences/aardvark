
use anyhow::Context;
use std::io::{BufWriter, Write};
use std::fs::File;
use std::path::Path;

/// Helper function that loads a file into some type, helpful generic
/// # Arguments
/// * `filename` - the file path to open and parse
/// # Errors
/// * if the file does not open properly
/// * if the deserialization throws errors
pub fn load_json<T: serde::de::DeserializeOwned>(filename: &Path) -> anyhow::Result<T> {
    let fp: Box<dyn std::io::Read> = if filename.extension().unwrap_or_default() == "gz" {
        Box::new(
            flate2::read::MultiGzDecoder::new(
                File::open(filename)?
            )
        )
    } else {
        Box::new(File::open(filename)?)
    };
    let result: T = serde_json::from_reader(fp)
        .with_context(|| format!("Error while deserializing {filename:?}:"))?;
    Ok(result)
}

/// This will save a generic serializable struct to JSON.
/// # Arguments
/// * `data` - the data in memory
/// * `out_filename` - user provided path to write to 
/// # Errors
/// * if opening or writing to the file throw errors
/// * if JSON serialization throws errors
pub fn save_json<T: serde::Serialize>(data: &T, out_filename: &Path) -> anyhow::Result<()> {
    let file: Box<dyn std::io::Write> = if out_filename.extension().unwrap_or_default() == "gz" {
        Box::new(
            flate2::write::GzEncoder::new(
                File::create(out_filename)?,
                flate2::Compression::best()
            )
        )
    } else {
        Box::new(File::create(out_filename)?)
    };
    let mut writer = BufWriter::new(file);
    serde_json::to_writer_pretty(&mut writer, data)
        .with_context(|| format!("Error while serializing {out_filename:?}:"))?;
    writer.flush()
        .with_context(|| format!("Error while flushing output to {out_filename:?}:"))?;
    Ok(())
}
