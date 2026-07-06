# rbin

Last change: 2026/07/02-13:12:57.

Hilofumi Yamamoto, Ph.D. Institute of Science Tokyo

## Rank Bin: rank/bin-based frequency distribution

[![DOI](https://zenodo.org/badge/1286672754.svg)](https://doi.org/10.5281/zenodo.21124284)

`rbin` is a small command-line utility for making rank/bin-based frequency distributions.

## Binning rules

- `-T0`: Sturges rule, default
- `-T1`: Scott rule
- `-T2`: Freedman-Diaconis rule

## Example

```sh
printf "1\n2\n3\n4\n5\n" | ./rbin -c -T0
printf "1\n2\n3\n4\n5\n" | ./rbin -c -T1
printf "1\n2\n3\n4\n5\n" | ./rbin -c -T2
```

## Compile

```sh
make
```

## License

This software is released under the MIT License, see LICENSE.txt for details.
