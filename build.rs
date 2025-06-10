
use std::error::Error;
use vergen_gitcl::{Emitter, GitclBuilder};

/// This function will emit the git instructions.
/// # Errors
/// * if `git` is not installed
/// * if there is no .git folder, which happens when we pull a version from GitHub on bioconda
fn emit_git() -> Result<(), Box<dyn Error>> {
    let gitcl = GitclBuilder::default()
        .all()
        .describe(false, true, Some("ThisPatternShouldNotMatchAnythingEver"))
        .build()?;

    Emitter::default()
        .fail_on_error()
        .add_instructions(&gitcl)?
        .emit()?;
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    match emit_git() {
        Ok(()) => {
            // normal path when we have a git clone or similar
        },
        Err(_e) => {
            // allow for user override if the data is missing, otherwise use "unknown"
            let git_desc = option_env!("CUSTOM_VERGEN_GIT_DESCRIBE")
                .unwrap_or("unknown");

            // we need to emit a custom one
            println!("cargo:rustc-env=VERGEN_GIT_DESCRIBE={git_desc}");
        }
    }
        
    // emit build handles the git configuration and build.rs, but we also need to track the toml and src folder 
    let rerun_if_changed = "cargo:rerun-if-changed=Cargo.toml
cargo:rerun-if-changed=src";
    println!("{rerun_if_changed}");

    Ok(())
}