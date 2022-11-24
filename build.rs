use cmake::Config;
use std::env;

fn main() {
    let custom_cc = env::var("CC");
    let custom_cxx = env::var("CXX");
    let conda_build = env::var("CONDA_BUILD");
    let mut is_conda_build = false;

    println!("cargo:rerun-if-changed=cuttlefish/CMakeLists.txt");
    println!("cargo:rerun-if-changed=piscem-cpp/CMakeLists.txt");

    let mut cfg_piscem_cpp = Box::new(Config::new("piscem-cpp"));
    let mut cfg_cf = Box::new(Config::new("cuttlefish"));

    (*cfg_cf).define("INSTANCE_COUNT", "32");
    if let Ok(cc_var) = custom_cc {
        (*cfg_piscem_cpp).define("CMAKE_C_COMPILER", cc_var.clone());
        (*cfg_cf).define("CMAKE_C_COMPILER", cc_var);
    }

    if let Ok(cxx_var) = custom_cxx {
        (*cfg_piscem_cpp).define("CMAKE_CXX_COMPILER", cxx_var.clone());
        (*cfg_cf).define("CMAKE_CXX_COMPILER", cxx_var);
    }

    if let Ok(_conda_build) = conda_build {
        (*cfg_cf).define("CONDA_BUILD", "TRUE");
        is_conda_build = true;
    }

    (*cfg_piscem_cpp).always_configure(false);
    (*cfg_cf).always_configure(false);

    let dst_piscem_cpp = (*cfg_piscem_cpp).build();
    let dst_cf = (*cfg_cf).build();

    println!(
        "cargo:rustc-link-search=native={}",
        dst_cf.join("lib").display()
    );
    println!(
        "cargo:rustc-link-search=native={}",
        dst_piscem_cpp.join("lib").display()
    );

    println!("cargo:rustc-link-lib=static=kmc_core");
    println!("cargo:rustc-link-lib=static=pesc_static");
    println!("cargo:rustc-link-lib=static=build_static");
    println!("cargo:rustc-link-lib=static=sshash_static");
    println!("cargo:rustc-link-lib=static=z");
    println!("cargo:rustc-link-lib=static=bz2");

    #[cfg(target_os = "linux")]
    println!("cargo:rustc-link-lib=dylib=stdc++");

    #[cfg(target_os = "macos")]
    println!("cargo:rustc-link-lib=dylib=c++");

    if is_conda_build {
        // if we are on OSX, building on conda
        // the filesystem support is borked and
        // we have to jump through some hoops.
        #[cfg(target_os = "macos")]
        println!("cargo:rustc-link-lib=dylib=c++experimental");
    }
}
