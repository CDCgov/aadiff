# AA Diff

A port of work by James Smagala for acceleration of RainbowTree and other use cases.

## Legacy algorithm

### Requirements

- First sequence is used as the reference
- Must be DNA
- Must be aligned
- Optionally: only alignable portions can be compared, e.g., the range containing the first and last non-ambiguous amino acid residue.

### Disambiguation

- Must be IUPAC
- Reference is not disambiguated
- Change: we only report up to 3 ambiguous translations

### Open questions

- Should we scrub delimiters from strain names? Otherwise we could fail.