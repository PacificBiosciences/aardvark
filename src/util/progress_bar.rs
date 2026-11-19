
use indicatif::{ProgressState, ProgressStyle};

/// Shared function to pull our progress bar styling
pub fn get_progress_style() -> ProgressStyle {
    ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({percent}); ETA: {eta_precise}; Speed: {per_sec} {msg}")
        .unwrap()
        .with_key("percent", |state: &ProgressState, w: &mut dyn std::fmt::Write| write!(w, "{:.1}%", state.fraction()*100.0).unwrap())
        .with_key("per_sec", |state: &ProgressState, w: &mut dyn std::fmt::Write| write!(w, "{:.0}/s", state.per_sec()).unwrap())
        .progress_chars("##-")
}

/// Progress style for a bar counted in basepairs.
/// The counter and speed are shown in megabasepairs.
pub fn get_basepair_progress_style() -> ProgressStyle {
    ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {mbp_pos}/{mbp_len} Mbp ({percent}); ETA: {eta_precise}; Speed: {per_sec} {msg}")
        .unwrap()
        .with_key("percent", |state: &ProgressState, w: &mut dyn std::fmt::Write| write!(w, "{:.1}%", state.fraction()*100.0).unwrap())
        .with_key("mbp_pos", |state: &ProgressState, w: &mut dyn std::fmt::Write| write!(w, "{:.2}", state.pos() as f64 / 1_000_000.0).unwrap())
        .with_key("mbp_len", |state: &ProgressState, w: &mut dyn std::fmt::Write| write!(w, "{:.2}", state.len().unwrap_or(0) as f64 / 1_000_000.0).unwrap())
        .with_key("per_sec", |state: &ProgressState, w: &mut dyn std::fmt::Write| write!(w, "{:.2} Mbp/s", state.per_sec() / 1_000_000.0).unwrap())
        .progress_chars("##-")
}
