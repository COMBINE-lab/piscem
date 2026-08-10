# piscem

`piscem` is a next-generation index + mapper tool consisting of a compacted de Bruijn graph
constructor (with tiling information), and a very efficient index and read mapper. It is a thin
command-line wrapper around [`piscem-rs`](https://github.com/COMBINE-lab/piscem-rs) (indexing and
mapping) and [`cf1-rs`](https://github.com/COMBINE-lab/cf1-rs) (cDBG construction); both are pure
Rust, so building `piscem` needs no C or C++ toolchain.

Notes
=====

 - If you are primarily interested in simply *using* `piscem`, you can obtain pre-compiled binaries
from the GitHub [releases](https://github.com/COMBINE-lab/piscem/releases) page for linux x86-64,
linux ARM, OSX x86-64, and OSX ARM. Likewise, `piscem` can be installed using
[bioconda](https://bioconda.github.io/recipes/piscem/README.html). The instructions below are
primarily for those who need to *build* `piscem` from source.

 - **Please ensure that the user file handle limit is set to 2048**. This may already be set (and
should be fine already on OSX), but you can accomplish this by executing:

```
$ ulimit -n 2048
```

before running `piscem`.

Building
========

```
git clone https://github.com/COMBINE-lab/piscem.git
cd piscem
cargo build --release
```

The checked-in `.cargo/config.toml` selects a *portable* baseline per target (`x86-64-v3` + AVX2 on
x86-64, `neoverse-n1` on linux ARM, `apple-a14` on OSX ARM). It deliberately does **not** use
`target-cpu=native`, because source packagers such as bioconda honour this file, and a binary built
for the build machine's exact CPU crashes with `SIGILL` elsewhere. If you are building only for the
machine you are on and want the last few percent, override it locally:

```
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

Usage
=====

```
Fast, space-efficient k-mer-based read mapper for bulk RNA-seq, scRNA-seq, and scATAC-seq

Usage: piscem [OPTIONS] <COMMAND>

Commands:
  build        Index a reference sequence
  map-sc       map reads for single-cell processing
  map-bulk     map reads for bulk processing
  map-sc-atac  map reads for scAtac processing
  help         Print this message or the help of the given subcommand(s)

Options:
  -q, --quiet    be quiet (no effect yet for cDBG building phase of indexing)
  -h, --help     Print help
  -V, --version  Print version
```

Each subcommand is summarised below. Run `piscem <subcommand> --help` for the authoritative and
complete option list — the options are numerous enough that reproducing them here only invites
drift.

build
-----

Indexes one or more reference sequences, building a piscem index over them. References are supplied
through exactly one of `--ref-seqs` (a `,` separated list of FASTA files), `--ref-lists` (files that
each list FASTA paths), or `--ref-dirs` (directories whose FASTA files are all indexed,
non-recursively).

The `-k` parameter is the k-mer size for the underlying colored compacted de Bruijn graph (must be
odd and `<= 31`), and `-m` is the minimizer size used to build the
[`sshash`](https://github.com/jermp/sshash) data structure (must be `< k`). Passing `--decoy-paths`
adds poison k-mer information from the given decoy sequences to the index.

Two directories are involved, and they are not the same:

  - `--work-dir` (default `./workdir.noindex`) holds the temporary files for **cDBG construction**.
  - `--tmp-dir` holds SSHash's **external minimizer-sort** scratch, a later phase. Unset, SSHash
    uses `sshash_tmp` in the current directory.

`--ram-limit-gib` caps the RAM that external minimizer sort may use (default 8 GiB); a smaller value
spills to disk sooner. Both are useful on HPC systems where `/tmp` is small or the job runs under a
cgroup memory limit. The index build also detects available memory automatically, walking the cgroup
chain so a container limit is respected rather than the host's total.

> **Note**
> If you are running on an Apple silicon machine, it is highly recommended that you set `-t` to be
> less than or equal to the number of **high performance** cores that you have (rather than the
> total number of cores including efficiency cores), as using efficiency cores in the `piscem build`
> step has been observed to severely degrade performance.

map-sc
------

Maps single-cell sequencing reads against a piscem index and produces a RAD format output file that
can be processed by [`alevin-fry`](https://github.com/COMBINE-lab/alevin-fry).

You can provide multiple files to `-1` and `-2` as a `,` separated list; it is important to provide
them in the same order to both flags. The `--geometry` flag specifies the geometry of the UMIs and
cell barcodes for the reads; see [geometry](#geometry) below.

map-bulk
--------

Maps bulk sequencing reads against a piscem index. The tool performs _non-spliced_ alignment, and is
therefore applicable to e.g. metagenomic reads against a set of metagenomes, DNA-seq alignment
against one or more references, or RNA-seq alignment against a transcriptome (but not a genome). It
produces a bulk RAD format output file that can be processed by
[`piscem-infer`](https://github.com/COMBINE-lab/piscem-infer) to estimate the abundances of all
references in the index given the mapped reads.

Reads are given either as paired `-1`/`-2` lists or as a single `--reads` list.

map-sc-atac
-----------

Maps single-cell ATAC-seq reads against a piscem index, producing a RAD file for
`alevin-fry atac`. Reads may be paired (`-1`/`-2`) or single (`--reads`), and the cell barcodes are
supplied separately through `--barcode`.

threading and decompression
---------------------------

For all three mapping subcommands, `-t` is a **single execution-slot budget shared between mapping
and gzip decoding** — not a mapping-thread count. The effective budget is capped by the parallelism
actually available, so a cpuset or cgroup limit is respected; both the requested and effective
values are recorded in `map_info.json`.

`--decoder` selects how gzip input is decompressed:

| value | behaviour |
| --- | --- |
| `auto` (default) | opens seekable gzip through the shared parallel decoder pool and lets a live thread broker adapt the mapping/decode split during the run |
| `serial` | disables the parallel pool and gives mapping the full budget |
| `parallel` | forces the parallel path wherever the input supports positional reads, still adapting the split |
| `parallel=N` | fixes `N` decode slots *per decoder-capable input* and disables adaptation |

Non-regular inputs — FIFOs, process substitutions such as `-1 <(zcat …)` — always use the serial
path, because the parallel decoder requires positional reads. The downgrade is per file, so one FIFO
among the inputs does not cost the regular files their parallel decoder.

`--thread-policy <FILE>` points at a JSON document overriding thread and decoder policy. Every field
is optional and falls back to a measured default, so a file need only name what it changes; an
unknown field is a hard error rather than a silent no-op. Currently:

```json
{ "parallel_decode": { "min_threads_per_stream": 8 } }
```

That value is how many threads must be available per gzip input before the parallel decoder is
engaged at all. Below it, the serial decoder is faster, because it decodes inline on threads that
also map and so is never idle, while a dedicated decode slot cannot map.

Note that with more than one mapping thread the byte layout of the RAD file varies from run to run,
because records are grouped into chunks by whichever worker finishes first. The mapping results are
unaffected. Use `-t 1` if you need byte-identical output to compare against.

For the full rationale and the measurements behind these defaults, see
[`docs/threading-and-decompression.md`](https://github.com/COMBINE-lab/piscem-rs/blob/main/docs/threading-and-decompression.md)
in `piscem-rs`.

geometry
--------

The geometry parameter `--geometry|-g` can take either a specific geometry name, or a generic
specifier string. The current valid names are `chromium_v2`, `chromium_v2_5p`, `chromium_v3`,
`chromium_v3_5p`, and `chromium_v4_3p`. The custom format is as follows: you must specify the
content of read 1 and read 2 in terms of the barcode, UMI, and mappable read sequence. A
specification looks like this:

```
1{b[16]u[12]x:}2{r:}
```

In particular, this is how one would specify the 10x Chromium v3 geometry using the custom syntax.
The format string says that the read pair should be interpreted as read 1 `1{...}` followed by read
2 `2{...}`. The syntax inside the `{}` says how the read should be interpreted. Here `b[16]u[12]x:`
means that the first 16 bases constitute the barcode, the next 12 constitute the UMI, and anything
that comes after that (if it exists) until the end of read 1 should be discarded (`x`). For read 2,
we have `2{r:}`, meaning that we should interpret read 2, in its full length, as biological
sequence.

It is possible to have pieces of geometry repeated, in which case they will be extracted and
concatenated together. For example, `1{b[16]u[12]b[4]x:}` would mean that we should obtain the
barcode by extracting bases 1-16 (1-based indexing) and 29-32 and concatenating them together to
obtain the full barcode. A specification that is followed by a specific length (i.e. a number in
`[]` like `b[10]` or `x[4]`) is said to be *bounded*. The specification string can have many bounded
pieces, but only one *unbounded* piece (an unbounded piece is a specifier like `r` or `x`, followed
by `:`). Likewise, since the `:` specifier means to extract this piece until the end of the string,
the unbounded specifier must be the last specifier in the description of each read (_if it occurs_).
