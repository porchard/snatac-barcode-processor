# Barcode preprocessing

This is a small preprocessor that can be used for the following:

1. Extracting barcodes (e.g., 10X cell barcodes) from sequencing reads. Depending on the experimental workflow used, barcodes are sometimes embedded in reads containing sequence beyond just the barcode, and the barcode may be reverse complemented. The `parse-barcodes` option takes a fastq file and a barcode whitelist, and tries to infer the location and orientation (reverse complemented or not) of the barcode in the read, and write a new fastq file containing the parsed and properly-oriented barcode.

2. Correcting barcodes, using an algorithm similar to that employed in CellRanger's ATAC workflow (described below).

## Installation

You must have Rust installed. You can then simply run `cargo build --release` to compile, which will create a binary at `./target/release/barcodes`.

## Usage