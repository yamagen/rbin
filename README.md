# rbin

Last change: 2026/11/04-17:26:02.

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

## Empty-bin hold

With the `-H` option, an empty bin inherits the displayed frequency of the preceding bin instead of being shown as `0`.

This changes only the displayed frequency distribution. The underlying bin counts and statistical calculations are not modified.

For example:

```sh
printf "1\n1\n1\n5\n5\n" | ./rbin -c
printf "1\n1\n1\n5\n5\n" | ./rbin -c -H
```

Without `-H`:

```text
    1.00  3
    2.33  0
    3.67  2
```

With `-H`:

```text
    1.00  3
    2.33  3
    3.67  2
```

This option is intended for exploratory plots in which an empty interval is displayed by carrying the preceding observed level forward.

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

### Vertical output

With `-V`, statistical output is printed as one item per line:

```sh
./rbin -nsV data.txt
```

For example:

```text
N	1111
min	0.000
mean	18.320
max	224.000
SD	13.233
SE	0.397
CV	72.237%
```

### JSON output

With `-j`, statistical output is printed in JSON format:

```sh
./rbin -nsj data.txt
```

For example:

```json
{
  "N": 1111,
  "min": 0.0,
  "mean": 18.32,
  "max": 224.0,
  "SD": 13.233,
  "SE": 0.397,
  "CV": 72.237
}
```

JSON output can be used directly with tools such as `jq`:

```sh
./rbin -nsj data.txt | jq '.mean'
```

```text
18.320
```

The `-n` option suppresses the histogram, allowing the JSON output to be passed directly to another command.

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
