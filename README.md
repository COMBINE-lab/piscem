# piscem

`piscem` is a rust wrapper for a next-generation index + mapper tool (still currently written in C++17).

This repository currently pulls in (as submodules) `piscem-cpp` and `cuttlefish`, and then uses cargo + the cmake crate to build the C++ code.  It then calls out to the main C++ function from `main.rs`.  The idea is that in this framework, code can slowly be migrated from C++ -> Rust in a piecemeal fashion, while the overall top-level repo remains a "rust build".  Importantly, `piscem` unifies and simplifies the other tools, making them runnable via a single executable and providing an improved command line interface.

Usage
=====

```
piscem
Indexing and mapping to compacted colored de Bruijn graphs

USAGE:
    piscem <SUBCOMMAND>

OPTIONS:
    -h, --help    Print help information

SUBCOMMANDS:
    build       Index a reference sequence
    help        Print this message or the help of the given subcommand(s)
    map-bulk    map bulk reads
    map-sc      map sc reads
```

`piscem` has several sub-commands; `build`, `map-sc` and `map-bulk` described below.

build
-----

The build subcommand indexes one or more reference sequences, building a piscem index over them.  The usage is as so:

```
piscem-build
Index a reference sequence

USAGE:
    piscem build [OPTIONS] --references <REFERENCES> --klen <KLEN> --mlen <MLEN> --threads <THREADS> --output <OUTPUT>

OPTIONS:
    -h, --help                       Print help information
    -k, --klen <KLEN>                length of k-mer to use
    -m, --mlen <MLEN>                length of minimizer to use
    -o, --output <OUTPUT>            output file stem
    -q                               be quiet during the indexing phase (no effect yet for cDBG building)
    -r, --references <REFERENCES>    reference FASTA location
    -t, --threads <THREADS>          number of threads to use
```

The parameters should be reasonably self-expalanatory.  The `-k` parameter is the k-mer size for the underlying colored compacted de Bruijn graph, and the `-m` parameter is the minimizer size used to build the [`sshash`](https://github.com/jermp/sshash) data structure.  The quiet `-q` flag applies to the `sshash` indexing step (not yet the CdBG construction step) and will prevent extra output being written to `stderr`.

Finally, the `-r` argument takes a list of `FASTA` format files containing the references to be indexed.  Here, if there is more than one reference, they should be provided to `-r` in the form of a `,` separated list.  For example, if you wish to index `ref1.fa`, `ref2.fa`, `ref3.fa` then your invocation should include `-r ref1.fa,ref2.fa,ref3.fa`.  The references present within all of the `FASTA` files will be indexed by the `build` command.


map-sc
------

The `map-sc` command maps single-cell sequencing reads against a piscem index, and produces a RAD format output file that can be processed by [`alevin-fry`](https://github.com/COMBINE-lab/alevin-fry).  The usage is as so:

```
piscem-map-sc
map sc reads

USAGE:
    piscem map-sc [OPTIONS] --index <INDEX> --geometry <GEOMETRY> --threads <THREADS> --output <OUTPUT>

OPTIONS:
    -1, --read1 <READ1>          path to list of read 1 files
    -2, --read2 <READ2>          path to list of read 1 files
    -g, --geometry <GEOMETRY>    geometry of barcode, umi and read
    -h, --help                   Print help information
    -i, --index <INDEX>          input index prefix
    -o, --output <OUTPUT>        path to output directory
    -q                           be quiet during mapping
    -t, --threads <THREADS>      number of threads to use
```

Here, you can provide multiple files to `-1` and `-2` as a `,` separated list just like the `-r` argument to the `build` command. Of course, it is important to ensure that you provide that information in the same order to the `-1` and `-2` flags.  The `--geometry` flag specifies the geometry of the UMIs and cell barcodes for the reads; you can find a description [here]().

