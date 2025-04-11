#![feature(let_chains)]
#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]

pub(crate) mod data;

use crate::data::GC3;
use clap::Parser;
use either::Either;
use std::{
    fs::OpenOptions,
    io::{BufReader, BufWriter, Write, stdin, stdout},
    ops::Range,
    path::PathBuf,
};
use zoe::{
    data::{
        fasta::{FastaAA, FastaNT},
        nucleotides::TranslatedNucleotidesIter,
    },
    prelude::*,
};

#[derive(Debug, Parser)]
#[command(about = "Tool for calculating amino acid difference tables")]
pub struct APDArgs {
    #[arg(short = 'i', long)]
    /// Optional input fasta
    input_fasta: Option<PathBuf>,

    #[arg(short = 'o', long)]
    /// Optional output delimited file
    output_xsv: Option<PathBuf>,

    #[arg(short = 'r', long)]
    /// Restrict to non-ambiguous alignable regions, pairwise.
    restrict_to_pairwise_alignable: bool,

    #[arg(short = 'e', long)]
    /// Use unix line-endings instead of Windows ones
    unix_line_endings: bool,

    #[arg(short = 'd', long)]
    /// Use the provider delimiter for separating fields. Default is ','
    output_delimiter: Option<char>,
}

fn main() {
    let args = APDArgs::parse();
    let line_ending = if args.unix_line_endings { "" } else { "\r" };
    let delim = args.output_delimiter.unwrap_or(',');

    let mut reader = if let Some(ref file_path) = args.input_fasta {
        FastaReader::new(BufReader::new(Either::Left(
            OpenOptions::new().read(true).open(file_path).expect("File opening error"),
        )))
    } else {
        FastaReader::new(BufReader::new(Either::Right(stdin())))
    };

    let mut writer = if let Some(ref file_path) = args.output_xsv {
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

    let reference = {
        let mut r = dna_reference.recode_to_dna();
        FastaAA {
            name:     r.name,
            sequence: r.sequence.to_aa_iter_with(b'X'),
        }
    };
    let ref_range = get_valid_range(&reference.sequence, args.restrict_to_pairwise_alignable);

    let other_sequences = reader
        .map(|record|
            // TODO: don't translate, instead defer until later
            record.map(|r| {
                let FastaNT { name, sequence } = r.recode_to_dna();
                let mut residues = sequence.to_aa_iter_with(b'X').collect();
                let valid_range = get_valid_range(&residues, args.restrict_to_pairwise_alignable);

                ValidSeq {
                    name, residues,
                    codons: sequence,
                    valid_range
                }
              }))
        .collect::<Result<Vec<_>, _>>()
        .unwrap_or_die("Could not process other data.");

    let mut buffer = format!("{delim}{name}", name = reference.name);
    for query_header in other_sequences.iter().map(|f| f.name.as_str()) {
        buffer.push(delim);
        buffer.push_str(query_header);
    }
    writeln!(&mut writer, "{buffer}{line_ending}").unwrap_or_fail();

    for (i, &ref_aa) in reference.sequence[ref_range].iter().enumerate() {
        let mut differences_found = false;
        buffer.clear();

        for ValidSeq {
            name: _,
            residues,
            codons,
            valid_range,
        } in other_sequences.iter()
        {
            let aa = residues[i];
            let codon = [codons[i * 3], codons[i * 3 + 1], codons[i * 3 + 2]];

            if valid_range.contains(&i) && ref_aa != aa {
                buffer.push_str(&format!("{delim}\""));

                if aa == b'-' {
                    buffer.push_str("del");
                } else if aa == b'X'
                    && let Some(degen_aa) = GC3.get(&codon)
                {
                    // We currently support degeneracy up to 3 distinct as beyond that it is kind of useless.
                    buffer.push_str(degen_aa);
                } else {
                    buffer.push(aa as char);
                }
                buffer.push('"');
                differences_found = true;
            } else {
                buffer.push(delim);
            }
        }

        if differences_found {
            writeln!(
                &mut writer,
                "{p}{delim}{aa}{buffer}{line_ending}",
                p = i + 1,
                aa = ref_aa as char
            )
            .unwrap_or_fail();
        }
    }

    writer.flush().unwrap_or_fail();
}

struct ValidSeq {
    name:        String,
    residues:    AminoAcids,
    codons:      Nucleotides,
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
