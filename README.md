# piscem

`piscem` is a rust wrapper for a next-generation index + mapper tool (still currently written in C++17).

Building
========

This repository currently pulls in (as submodules) `piscem-cpp` and `cuttlefish`, and then uses cargo + the cmake crate to build the C++ code.  It then calls out to the main C++ function from `main.rs`.  The idea is that in this framework, code can slowly be migrated from C++ -> Rust in a piecemeal fashion, while the overall top-level repo remains a "rust build".  Importantly, `piscem` unifies and simplifies the other tools, making them runnable via a single executable and providing an improved command line interface.

To build `piscem`, first check out this repository **with recursive dependencies**:

```
git clone --recursive https://github.com/COMBINE-lab/piscem.git
```

if you have accidentally checked out the repo without the `--recursive` flag, you can change into the top-level directory and run:

```
git submodule update --init --recursive
```

Once you have checked out the repository, you can build `piscem` with 

```
cargo build --release
```

It is worth noting that the build process respects the `CC` and `CXX` environment flags, so if you wish to use a specific C++ compiler, you can run:

```
CC=<path_to_c_compiler> CXX=<path_to_cxx_compiler> cargo build --release
```


Usage
=====

```
piscem 0.1.0
Indexing and mapping to compacted colored de Bruijn graphs

USAGE:
    piscem <SUBCOMMAND>

OPTIONS:
    -h, --help       Print help information
    -V, --version    Print version information

SUBCOMMANDS:
    build       Index a reference sequence
    help        Print this message or the help of the given subcommand(s)
    map-bulk    map reads for bulk processing
    map-sc      map reads for single-cell processing
```

`piscem` has several sub-commands; `build`, `map-sc` and `map-bulk` described below.

Info for different sub-commands


build
-----

The build subcommand indexes one or more reference sequences, building a piscem index over them.  The usage is as so:

```
piscem-build
Index a reference sequence

USAGE:
    piscem build [OPTIONS] --klen <KLEN> --mlen <MLEN> --threads <THREADS> --output <OUTPUT> <--ref-seqs <REF_SEQS>|--ref-lists <REF_LISTS>|--ref-dirs <REF_DIRS>>

OPTIONS:
    -d, --ref-dirs <REF_DIRS>      ',' separated list of directories (all FASTA files in each directory will be indexed, but not recursively)
    -h, --help                     Print help information
    -k, --klen <KLEN>              length of k-mer to use
    -l, --ref-lists <REF_LISTS>    ',' separated list of files (each listing input FASTA files)
    -m, --mlen <MLEN>              length of minimizer to use
    -o, --output <OUTPUT>          output file stem
    -q                             be quiet during the indexing phase (no effect yet for cDBG building)
    -s, --ref-seqs <REF_SEQS>      ',' separated list of reference FASTA files
    -t, --threads <THREADS>        number of threads to use
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

Here, you can provide multiple files to `-1` and `-2` as a `,` separated list just like the `-r` argument to the `build` command. Of course, it is important to ensure that you provide that information in the same order to the `-1` and `-2` flags.  The `--geometry` flag specifies the geometry of the UMIs and cell barcodes for the reads; you can find a description [here](https://github.com/COMBINE-lab/piscem/blob/main/README.md#geometry).

map-bulk
--------

The `map-bulk` command maps bulk sequencing reads against a piscem index. The tool performs _non-spliced_ alignment, and therefore is applicable to e.g. metagenomic reads against a set of metagenomes, DNA-seq alignment against one or more references, or RNA-seq alignment against a transcriptome (but not a genome). The program and produces a bulk RAD format output file that can be processed by [`piscem-infer`](https://github.com/COMBINE-lab/alevin-fry](https://github.com/COMBINE-lab/piscem-infer) to estimate the abundances of all references in the index given the mapped reads.  The usage is as so:

```
piscem-map-bulk
map bulk reads

USAGE:
    piscem map-bulk [OPTIONS] --index <INDEX> --read1 <READ1> --read2 <READ2> --threads <THREADS> --output <OUTPUT>

OPTIONS:
    -1, --read1 <READ1>        path to list of read 1 files
    -2, --read2 <READ2>        path to list of read 1 files
    -h, --help                 Print help information
    -i, --index <INDEX>        input index prefix
    -o, --output <OUTPUT>      path to output directory
    -q                         be quiet during mapping
    -t, --threads <THREADS>    number of threads to use
```

Here, you can provide multiple files to `-1` and `-2` as a `,` separated list just like the `-r` argument to the `build` command. Of course, it is important to ensure that you provide that information in the same order to the `-1` and `-2` flags.

geometry
--------

The geometry parameter `--geometry|-g` can take either a specific geometry name, or a generic specifier string.  The current valid names are `chromium_v2` and `chromium_v3` for 10x Genomics Chromium v2 and v3 protocols respectively. The custom format is as follows: you must specify the content of read 1 and read 2 in terms of the barcode, UMI, and mappable read sequence. A specification looks like this:

```
1{b[16]u[12]x:}2{r:}
```

In particular, this is how one would specify the 10x Chromium v3 geometry using the custom syntax.  The format string says that the read pair should be interpreted as read 1 `1{...}` followed by read 2 `2{...}`.  The syntax inside the `{}` says how the read should be interpreted.  Here `b[16]u[12]x:` means that the first 16 bases constitute the barcode, the next 12 constitute the UMI, and anything that comes after that (if it exists) until the end of read 1 should be discarded (`x`).  For read 2, we have `2{r:}`, meaning that we should interpret read 2, in it's full length, as biological sequence.

It is possible to have pieces of geometry repeated, in which case they will be extracted and concatenated together.  For example, `1{b[16]u[12]b[4]x:}` would mean that we should obtain the barcode by extracting bases 1-16 (1-based indexing) and 29-32 and concatenating them togehter to obtain the full barcode.  A specification that is followed by a specific length (i.e. a number in `[]` like `b[10]` or `x[4]` is said to be *bounded*).  The specification string can have many bounded pieces, but only one *unbounded* piece (and unbounded piece is a specifier like `r` or `x`, followed by `:`).  Likewise, since the `:` specifier means to extract this piece until the end of the string, the unbounded specifier must be the last specifier in the description of each read (_if it occurs_).
