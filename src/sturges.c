#include "config.h"
#include <math.h>
#include <stddef.h>

int sturges(double count) { return (int)(1 + log(count) / log(2)); }

int sturges_two(double count) { return (int)(1 + 3.32 * log10(count)); }

double scott_bin_width(size_t n, double s) {
  /*     Scott(1979)
         h = 3.5 s n^{-1/3}
         Where, h is the width of a class, s is a standard deviation,
         n is the number of samples.
  */

  return (int)(3.5 * s * pow(n, -1.0 / 3.0));

  if (n <= 1 || s <= 0.0)
    return 0.0;

  return 3.5 * s * pow((double)n, -1.0 / 3.0);
}
