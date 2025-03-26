# AA Diff

A port of work by James Smagala for acceleration of RainbowTree and other use cases.

## Legacy algorithm

### Requirements

- First sequence is used as the reference
- Must be DNA
- Must be aligned
- Optionally: only alignable portions can be compared, e.g., the range containing the first and last non-ambiguous amino acid residue.
- To be confirmed: deletion variants may be denoted as "del"

### Disambiguation

- Must be IUPAC
- Reference is not disambiguated
- All possible amino acid translations will be reported.

Note: would it be better to limit the translations to 3 AA? Beyond that it is kind of silly.