use std::error::Error;
use vergen_gitcl::{Emitter, GitclBuilder};

fn main() -> Result<(), Box<dyn Error>> {
    let gitcl = GitclBuilder::default()
        .all()
        .describe(false, true, Some("ThisPatternShouldNotMatchAnythingEver"))
        .build()?;

    Emitter::default()
        .fail_on_error()
        .add_instructions(&gitcl)?
        .emit()?;
        
    // emit build handles the git configuration and build.rs, but we also need to track the toml and src folder 
    let rerun_if_changed = "cargo:rerun-if-changed=Cargo.toml
cargo:rerun-if-changed=src";
    println!("{rerun_if_changed}");

    Ok(())
}