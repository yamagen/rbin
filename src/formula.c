#include "formula.h"
#include "config.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

/* Freedman-Diaconis bin width */
double fd_bin_width(double *data, size_t n);

/* IQR */
double calculate_iqr(double *data, size_t n);

/* qsort comparison */
int compare(const void *a, const void *b);

double skew(int n, double s2, double s3) {
  /* skewness */
  return s3 / (n * sqrt(s2 / n) * sqrt(s2 / n) * sqrt(s2 / n));
}

double kurt(int n, double s2, double s4) {
  /* kurtosis */
  return s4 / (n * sqrt(s2 / n) * sqrt(s2 / n) * sqrt(s2 / n) * sqrt(s2 / n));
}

double sdev(int n, double s2) {
  /* standard deviation */
  return sqrt(s2 / (n - 1));
}

double serr(int n, double s2) {
  /* standard error */
  return sdev(n, s2) / sqrt(n);
}

double cvar(int n, double s1, double s2) {
  /*
           coefficient of variation -- Mizutani (1983:44)
           正しくは、変異係数=母分散/母平均*100で、％表示が慣例である。
  */
  if (s1 != 0)
    return 100 * sdev(n, s2) / s1;
  else
    return 0;
}

double tvar(int n, double s1, double s2) {
  /* tasseido of variation -- Mizutani (1983:64) 作業中 */
  if (s1 != 0)
    return 100 * serr(n, s2) / s1;
  else
    return 0;
}

double rlo(int n, double s1, double s2, double k) {
  /* reliable interval which is lower (Mizuno, Hitoshi)*/
  if (k == 0)
    k = 2.6; /* default */
  return ((s1 - serr(n, s2) * k) * 100) / 100;
}

double rup(int n, double s1, double s2, double k) {
  /* reliable interval which is upper (Mizuno, Hitoshi: Tchevichev)*/
  if (k == 0)
    k = 2.6; /* default */
  return ((s1 + serr(n, s2) * k) * 100) / 100;
}

//
/*
 * The following functions are used to calculate the bin width and the number of bins for a
 * histogram. The bin width is the width of each bin in the histogram, and the number of bins is the
 * number of bins in the histogram. The bin width is calculated using the Wand bin width formula,
 * and the number of bins is calculated using the Freedman-Diaconis rule. The Wand bin width formula
 * is used to calculate the bin width based on the standard deviation of the data and the number of
 * data points. The Freedman-Diaconis rule is used to calculate the number of bins based on the
 * interquartile range of the data and the number of data points. The bin width and the number of
 * bins are used to create a histogram that accurately represents the data. In statistics, the
 * Freedman芯iaconis rule can be used to select the width of the bins to be used in a histogram.[1]
 * It is named after David A. Freedman and Persi Diaconis.
 */

// double fd_class_width(int n, double IQ) {
//     // Freedman-Diaconis
//     if (n <= 0) {
//         return 0.0;
//     }
//     return 2.0 * IQ / pow(n, 1.0 / 3.0);
// }

double wand_bin_width(int n, double sigma) {
  if (n <= 0 || sigma <= 0.0) {
    return 0.0;
  }
  return 3.49 * sigma * pow(n, -1.0 / 3.0);
}

int calculate_bins(double min, double max, double bin_width) {
  if (bin_width <= 0) {
    return 1;
  }
  return (int)ceil((max - min) / bin_width);
}

// 関数プロトタイプ
double calculate_bin_width_fd(double *data, size_t n, int use_iqr);

// IQR（四分位範囲）を計算する関数
double calculate_iqr(double *data, size_t n);

// 比較関数（qsort用）
int compare(const void *a, const void *b);

#ifdef __test
int main(void) {
  double data[] = {1.1, 2.2, 3.5, 5.6, 7.8, 3.3, 4.4, 6.6, 8.8, 10.1};
  size_t n = sizeof(data) / sizeof(data[0]);

  printf("FD bin width: %f\n", fd_bin_width(data, n));

  return 0;
}
#endif

// Freedman-Diaconisの計算関数
/* Freedman-Diaconis bin width */
double fd_bin_width(double *data, size_t n) {
  double *sorted_data;
  double iqr;
  double h;

  if (n <= 1) {
    fprintf(stderr, "Error: not enough data points.\n");
    return -1.0;
  }

  sorted_data = malloc(n * sizeof(double));
  if (sorted_data == NULL) {
    fprintf(stderr, "Error: memory allocation failed.\n");
    return -1.0;
  }

  for (size_t i = 0; i < n; i++)
    sorted_data[i] = data[i];

  qsort(sorted_data, n, sizeof(double), compare);

  iqr = calculate_iqr(sorted_data, n);
  free(sorted_data);

  if (iqr <= 0.0) {
    fprintf(stderr, "Error: IQR is zero, cannot calculate bin width.\n");
    return -1.0;
  }

  h = 2.0 * iqr / cbrt((double)n);

  return h;
}

// 四分位範囲（IQR）の計算
double calculate_iqr(double *data, size_t n) {
  double q1, q3;
  size_t mid = n / 2;

  // 第一四分位数
  if (mid % 2 == 0) {
    q1 = (data[mid / 2 - 1] + data[mid / 2]) / 2.0;
  } else {
    q1 = data[mid / 2];
  }

  // 第三四分位数
  if (mid % 2 == 0) {
    q3 = (data[n - mid / 2 - 1] + data[n - mid / 2]) / 2.0;
  } else {
    q3 = data[n - mid / 2 - 1];
  }

  return q3 - q1;
}

// qsort用の比較関数
int compare(const void *a, const void *b) {
  double diff = *(double *)a - *(double *)b;
  return (diff > 0) - (diff < 0);
}

// jarque_bera統計量の計算
double jarque_bera(int n, double s2, double s3, double s4) {
  double g1, g2;
  double jb;

  if (n < 8 || s2 <= 0.0)
    return -1.0;

  g1 = skew(n, s2, s3);       /* skewness */
  g2 = kurt(n, s2, s4) - 3.0; /* excess kurtosis */

  jb = (n / 6.0) * (g1 * g1 + (g2 * g2) / 4.0);

  return jb;
}

double jarque_bera_pvalue(double jb) {
  if (jb < 0.0)
    return -1.0;

  /*
   * Jarque-Bera statistic is asymptotically chi-square with df=2.
   * For chi-square df=2, upper-tail probability is exp(-x / 2).
   */
  return exp(-jb / 2.0);
}

int sturges(double count) { return (int)(1 + log(count) / log(2)); }

int sturges_two(double count) { return (int)(1 + 3.32 * log10(count)); }

double scott_bin_width(size_t n, double s) {
  /*
   * Scott (1979)
   * h = 3.5 * s * n^(-1/3)
   * h: bin width
   * s: standard deviation
   * n: number of samples
   */

  if (n <= 1 || s <= 0.0)
    return 0.0;

  return 3.5 * s * pow((double)n, -1.0 / 3.0);
}
