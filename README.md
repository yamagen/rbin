# rbin

Last change: 2026/07/07-08:28:49.

Hilofumi Yamamoto, Ph.D.
Institute of Science Tokyo

## Rank Bin: rank/bin-based frequency distribution

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21229729.svg)](https://doi.org/10.5281/zenodo.21229729)

`rbin` is a small command-line utility for making rank/bin-based frequency distributions from numeric input.

It reads numeric values from standard input or files and prints a frequency distribution using rank/bin classes. With the `-s` option, it also reports basic descriptive statistics and distribution diagnostics.

## Binning rules

- `-T0`: Sturges rule, default
- `-T1`: Scott rule
- `-T2`: Freedman-Diaconis rule

For Scott and Freedman-Diaconis, `rbin` first calculates a method-specific bin width and then converts it into the actual number of bins used for output.

## Statistical output

With `-s`, `rbin` prints:

- `N`, `min`, `mean`, `max`
- `SD`, `SE`, `CV`
- reliable interval
- skewness
- kurtosis
- Jarque-Bera statistic
- Jarque-Bera p-value
- normality judgment
- rank/bin count
- actual bin width, `delta`
- input accuracy
- binning method
- method-specific bin width for Scott and Freedman-Diaconis

The normality judgment is based on the Jarque-Bera statistic. `normality=not-rejected` means that normality is not rejected by this test; it does not prove that the data are normally distributed.

## Example

```sh
printf "1\n2\n3\n4\n5\n" | ./rbin -c -T0
printf "1\n2\n3\n4\n5\n" | ./rbin -c -T1
printf "1\n2\n3\n4\n5\n" | ./rbin -c -T2
```

With statistical output:

```sh
printf "1\n2\n3\n4\n5\n" | ./rbin -s -T2
```

Example output includes lines such as:

```text
skew=0.228 kurtosis=2.170 JB=0.935 JB_p=0.6265 normality=not-rejected
rank=4 delta=1.250 accuracy=1.000000
method=Freedman-Diaconis method_bin_width=1.368
```

## Compile

```sh
make
```

## Install

```sh
make install
```

To install under a custom prefix:

```sh
make PREFIX=$HOME/.local install
```

## Clean

```sh
make clean
```

## Related tools

`rbin` is part of a small family of command-line exploratory data filters, together with `cw`, SLP (stem leaf plotter), and KSF (kanji statistics filter).

## License

This software is released under the MIT License. See `LICENSE.txt` for details.
