fn main() {
    // Extract cf1-rs version from Cargo.lock so we can embed it in the version JSON.
    let lock_contents =
        std::fs::read_to_string("Cargo.lock").expect("Failed to read Cargo.lock");
    let lock: toml::Value = lock_contents.parse().expect("Failed to parse Cargo.lock");
    let version = lock["package"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["name"].as_str() == Some("cf1-rs"))
        .and_then(|p| p["version"].as_str())
        .unwrap_or("unknown");
    println!("cargo:rustc-env=CF1_RS_VERSION={version}");
    println!("cargo:rerun-if-changed=Cargo.lock");
}
