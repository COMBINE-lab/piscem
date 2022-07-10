# piscem

Rust wrapper for the next generation index + mapper (still currently in C++)

This repository currently pulls in (as a submodule) `piscem-cpp`, and then uses 
cargo + the cmake crate to build the C++ code.  It then calls out to the main C++ 
function from `main.rs`.  The idea is that in this framework, code can slowly be 
migrated from C++ -> Rust in a piecemeal fashion, while the overall top-level repo
remains a "rust build".
