# Scott rule calculation

## Purpose

Scott rule estimates the bin width from the standard deviation and the number of observations.

It is based on the spread of the data.

## Formula

```text
h = 3.49 * SD * n^(-1/3)
```

where:

- `h` is the bin width
- `SD` is the sample standard deviation
- `n` is the number of observations

The number of bins is then calculated from the range and bin width:

```text
k = ceil((max - min) / h)
```

## Required values

- `n`
- `min`
- `max`
- standard deviation

## Calculation steps

1. Count the number of observations.
2. Calculate the standard deviation.
3. Calculate the bin width using Scott rule.
4. Divide the range by the bin width.
5. Round up to obtain the number of bins.

## Behavior

Scott rule uses the standard deviation, so it is sensitive to the overall spread of the data.

If the distribution has a long tail, Scott rule may produce very wide bins and the resulting histogram may become too coarse.

## Use in rbin

```sh
rbin -c -T1
```

## Notes

In the cw workflow, Scott rule is useful as a comparison, but it may be too coarse for z-value distributions with a long right tail.
