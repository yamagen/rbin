# Sturges rule calculation

## Purpose

Sturges rule estimates the number of bins from the number of observations.

It is useful for obtaining a rough overview of the distribution.

## Formula

```text
k = 1 + log2(n)
```

where:

- `k` is the number of bins
- `n` is the number of observations

## Required values

Only `n` is required.

## Calculation steps

1. Count the number of observations.
2. Calculate `1 + log2(n)`.
3. Round the result to an integer.
4. Use this value as the number of bins.

## Behavior

Sturges rule tends to produce a small number of bins.

It is useful for readable overview figures, but it may be too coarse for large datasets or strongly skewed distributions.

## Use in rbin

```sh
rbin -c -T0
```

`-T0` is the default binning rule.

## Notes

In the cw workflow, Sturges rule is useful for showing the general shape of the cooccurrence weight distribution without flooding the reader with too much detail.
