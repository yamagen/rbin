#include "config.h"
#include <math.h>

int sturges(double count) { return (int)(1 + log(count) / log(2)); }

int sturges_two(double count) { return (int)(1 + 3.32 * log10(count)); }

int scott(double n, double s) {

  /*     Scott(1979)
         h = 3.5 s n^{-1/3}
         Where, h is the width of a class, s is a standard deviation,
         n is the number of samples.

  */
  return (int)(3.5 * s * pow(n, -1 / 3));
}
