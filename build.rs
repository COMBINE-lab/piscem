use cmake::Config;
use std::env;

fn main() {
    let custom_cc = env::var("CC");
    let custom_cxx = env::var("CXX");
    let conda_build = env::var("CONDA_BUILD");
    let nopie_build = env::var("NOPIE");
    let zlib_ng_path = env::var("ZLIB_NG_PATH");

    let is_conda_build = match conda_build {
        Ok(val) => match val.to_uppercase().as_str() {
            "TRUE" | "1" | "YES" => true,
            "FALSE" | "0" | "NO" => false,
            _ => true,
        },
        Err(_e) => false,
    };

    println!("cargo:rerun-if-changed=cuttlefish/CMakeLists.txt");
    println!("cargo:rerun-if-env-changed=ZLIB_NG_PATH");

    let mut cfg_cf = Box::new(Config::new("cuttlefish"));

    (*cfg_cf).define("INSTANCE_COUNT", "32");
    if let Ok(cc_var) = custom_cc {
        (*cfg_cf).define("CMAKE_C_COMPILER", cc_var);
    }

    if let Ok(cxx_var) = custom_cxx {
        (*cfg_cf).define("CMAKE_CXX_COMPILER", cxx_var);
    }

    if is_conda_build {
        (*cfg_cf).define("CONDA_BUILD", "TRUE");
        (*cfg_cf).define("CMAKE_OSX_DEPLOYMENT_TARGET", "10.15");
        (*cfg_cf).define("MACOSX_SDK_VERSION", "10.15");
    }

    (*cfg_cf).always_configure(false);

    let dst_cf = (*cfg_cf).build();

    if let Ok(nopie) = nopie_build {
        match nopie.as_str() {
            "1" | "TRUE" | "true" | "True" => {
                println!("cargo:rustc-link-arg=-no-pie");
            }
            _ => {}
        }
    }

    println!(
        "cargo:rustc-link-search=native={}",
        dst_cf.join("lib").display()
    );

    // For some reason, if we are using
    // *some* linux distros (and on conda) and are
    // building for the linux target;
    // things get put in the lib64 directory
    // rather than lib... So, we add that here
    println!(
        "cargo:rustc-link-search=native={}",
        dst_cf.join("lib64").display()
    );

    println!("cargo:rustc-link-lib=static=kmc_core");

    // --- zlib-ng (compat mode) linking ---
    // Tier 1: user-provided path via ZLIB_NG_PATH
    if let Ok(zlib_path) = zlib_ng_path {
        let zlib_path = std::path::PathBuf::from(&zlib_path);
        println!("cargo:rustc-link-search=native={}", zlib_path.display());
        println!("cargo:rustc-link-lib=static=z");
    } else {
        // Tier 2: try pkg-config
        let pkg_config_result = pkg_config::Config::new()
            .statik(true)
            .cargo_metadata(false)
            .probe("zlib-ng");

        match pkg_config_result {
            Ok(lib) => {
                // Emit the search paths and link the static library ourselves
                // (cargo_metadata is false so we control the output)
                for path in &lib.link_paths {
                    println!("cargo:rustc-link-search=native={}", path.display());
                }
                println!("cargo:rustc-link-lib=static=z");
            }
            Err(_) => {
                // Tier 3: build vendored zlib-ng submodule with CMake
                println!("cargo:rerun-if-changed=zlib-ng/CMakeLists.txt");
                let dst_zlib = Config::new("zlib-ng")
                    .define("ZLIB_COMPAT", "ON")
                    .define("BUILD_SHARED_LIBS", "OFF")
                    .define("ZLIB_ENABLE_TESTS", "OFF")
                    .define("ZLIBNG_ENABLE_TESTS", "OFF")
                    .always_configure(false)
                    .build();

                // zlib-ng installs into lib or lib64 depending on platform
                println!(
                    "cargo:rustc-link-search=native={}",
                    dst_zlib.join("lib").display()
                );
                println!(
                    "cargo:rustc-link-search=native={}",
                    dst_zlib.join("lib64").display()
                );
                // In compat mode, the library is named libz (same as zlib)
                println!("cargo:rustc-link-lib=static=z");
            }
        }
    }

    #[cfg(target_os = "linux")]
    println!("cargo:rustc-link-lib=dylib=stdc++");

    #[cfg(target_os = "macos")]
    println!("cargo:rustc-link-lib=dylib=c++");
}
