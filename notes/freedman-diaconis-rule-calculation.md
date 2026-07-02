# Freedman-Diaconis rule calculation

## Purpose

Freedman-Diaconis rule estimates the bin width from the interquartile range.

It is useful for observing distribution shape while reducing the effect of extreme values.

## Formula

```text
h = 2 * IQR * n^(-1/3)
```

where:

- `h` is the bin width
- `IQR` is the interquartile range
- `n` is the number of observations

The number of bins is then calculated from the range and bin width:

```text
k = ceil((max - min) / h)
```

## Required values

- `n`
- `min`
- `max`
- first quartile
- third quartile
- interquartile range

## Calculation steps

1. Store all observations.
2. Sort the observations.
3. Calculate the first quartile.
4. Calculate the third quartile.
5. Calculate `IQR = Q3 - Q1`.
6. Calculate the bin width using Freedman-Diaconis rule.
7. Divide the range by the bin width.
8. Round up to obtain the number of bins.

## Behavior

Freedman-Diaconis rule is based on the middle half of the data.

It is less affected by extreme values than Scott rule. It often gives a finer view of the central part and the tail of the distribution.

## Use in rbin

```sh
rbin -c -T2
```

## Notes

In the cw workflow, Freedman-Diaconis rule is useful for checking whether the general distribution shape seen with Sturges rule is stable.

It may produce many bins, so it is often better for inspection than for printed figures.
