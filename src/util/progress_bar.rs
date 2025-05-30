
use indicatif::{ProgressState, ProgressStyle};

/// Shared function to pull our progress bar styling
pub fn get_progress_style() -> ProgressStyle {
    ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({percent}); ETA: {eta_precise}; Speed: {per_sec} {msg}")
        .unwrap()
        .with_key("percent", |state: &ProgressState, w: &mut dyn std::fmt::Write| write!(w, "{:.1}%", state.fraction()*100.0).unwrap())
        .with_key("per_sec", |state: &ProgressState, w: &mut dyn std::fmt::Write| write!(w, "{:.0}/s", state.per_sec()).unwrap())
        .progress_chars("##-")
}
