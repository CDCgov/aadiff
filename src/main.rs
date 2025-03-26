#![feature(let_chains)]
#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]

use clap::Parser;
use either::Either;
use std::{
    fs::OpenOptions,
    io::{BufReader, BufWriter, Write, stdin, stdout},
    ops::Range,
    path::PathBuf,
};
use zoe::{data::fasta::FastaNT, prelude::*};

// TODO: we want to change the delimiter

/// Holds arguments for the APD CLI.
#[derive(Debug, Parser)]
#[command(about = "Tool for calculating amino acid difference tables")]
pub struct APDArgs {
    #[arg(short = 'i', long)]
    /// Optional input fasta
    input_fasta: Option<PathBuf>,

    #[arg(short = 'o', long)]
    /// Optional output delimited
    output_tsv: Option<PathBuf>,

    #[arg(short = 't', long)]
    /// Turns on simultaneous multi-threading
    multi_threaded: bool,

    // This appeared to be available in the legacy program.
    #[arg(short = 'r', long)]
    /// Restrict to non-ambiguous alignable regions, pairwise.
    restrict_to_pairwise_alignable: bool,

    #[arg(short = 'e', long)]
    /// Use unix line-endings instead of Windows ones
    unix_line_endings: bool,
}

/// CLI entry point
fn main() {
    let args = APDArgs::parse();
    let line_ending = if args.unix_line_endings { "" } else { "\r" };

    let mut reader = if let Some(ref file_path) = args.input_fasta {
        FastaReader::new(BufReader::new(Either::Left(
            OpenOptions::new().read(true).open(file_path).expect("File opening error"),
        )))
    } else {
        FastaReader::new(BufReader::new(Either::Right(stdin())))
    };

    let mut writer = if let Some(ref file_path) = args.output_tsv {
        BufWriter::new(Either::Left(
            OpenOptions::new()
                .write(true)
                .create(true)
                .truncate(true)
                .open(file_path)
                .expect("File write error"),
        ))
    } else {
        BufWriter::new(Either::Right(stdout()))
    };

    let Some(Ok(dna_reference)) = reader.next() else {
        eprintln!("No first record available!");
        std::process::exit(1);
    };
    // TODO: recode to IUPAC if not in IUPAC.
    let reference = {
        let mut r = dna_reference.translate();
        r.sequence.replace_all_bytes(b'~', b'X');
        r
    };
    let ref_range = get_valid_range(&reference.sequence, args.restrict_to_pairwise_alignable);

    let other_sequences = reader
        .map(|record|
            // TODO: recode to IUPAC
            // TODO: don't translate, instead defer until later
            record.map(|r| {
                let FastaNT { name, sequence } = r.into();
                let mut residues: AminoAcids = sequence.to_aa_iter().collect();
                residues.replace_all_bytes(b'~',b'X');
                let valid_range = get_valid_range(&residues, args.restrict_to_pairwise_alignable);

                ValidSeq {
                    name, residues,
                    _codons: sequence,
                    valid_range
                }
              }))
        .collect::<Result<Vec<_>, _>>()
        .unwrap_or_die("Could not process other data.");

    let mut buffer = format!(",{}", reference.name);
    for query_header in other_sequences.iter().map(|f| f.name.as_str()) {
        buffer.push(',');
        buffer.push_str(query_header);
    }
    writeln!(&mut writer, "{buffer}{line_ending}").unwrap_or_fail();

    for (i, &ref_aa) in reference.sequence[ref_range].iter().enumerate() {
        let mut differences_found = false;
        buffer.clear();

        // TODO: grab codon in the future instead and consider not resolvable ambiguations
        // TODO: make sure you use uppercase for both
        for (query_aa, valid_range) in other_sequences.iter().map(|f| (f.residues[i], &f.valid_range)) {
            if valid_range.contains(&i) && ref_aa != query_aa {
                buffer.push_str(r#",""#);

                if query_aa == b'-' {
                    buffer.push_str("del");
                } else {
                    buffer.push(query_aa as char);
                }
                buffer.push('"');
                differences_found = true;
            } else {
                buffer.push(',');
            }
        }

        if differences_found {
            writeln!(&mut writer, "{p},{aa}{buffer}{line_ending}", p = i + 1, aa = ref_aa as char).unwrap_or_fail();
        }
    }

    writer.flush().unwrap_or_fail();
}

struct ValidSeq {
    name:        String,
    residues:    AminoAcids,
    _codons:     Nucleotides,
    valid_range: std::ops::Range<usize>,
}

fn get_valid_range(aa: &AminoAcids, restrict: bool) -> Range<usize> {
    if restrict {
        let (Some(s), Some(e)) = (
            aa.iter().position(|&aa| aa != b'X' && aa != b'-'),
            aa.iter().rposition(|&aa| aa != b'X' && aa != b'-'),
        ) else {
            eprintln!("Sequence doesn't contain valid data for comparison.");
            std::process::exit(1);
        };

        s..e + 1
    } else {
        0..aa.len()
    }
}
