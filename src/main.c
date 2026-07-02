#include "config.h"
#include <stdio.h>
#include <stdlib.h> /* EXIT_FAILURE,exit */
#include <string.h>

double fd_bin_width(double *data, size_t n);
int calculate_bins(double min, double max, double bin_width);

int digitstr(char *str);
int count_freq(FILE *fp);
int sturges(double count);
int scott(double n, double s);
double check_precision(char *str);
int main(int argc, char **argv);
double cvar(int n, double s1, double s2);
double kurt(int n, double s2, double s4);
double sdev(int n, double s2);
double serr(int n, double s2);
double skew(int n, double s2, double s3);
double rlo(int n, double s1, double s2, double k);
double rup(int n, double s1, double s2, double k);
double tvar(int n, double s1, double s2);

int screen_w = 60; /* screen width default */
int append_return = 0;
int only_num = 0;
int with_stat = 0;
int out_format = 0;
int print_rank = 0;
int err_f = 0;
int int_f = 0;
double prec_level = 1;
int rank_type = 0; /* default is sturges. */

int count_freq(FILE *fp) {

  char buf[BUF];
  char p[BUF];
  int n = 0;
  int i = 0, k = 0;
  int rank = 0;
  int out_type = 0;

  double x = 0.0;
  double max = 0.0;
  double min = 0.0;
  double rmax = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  double s3 = 0.0;
  double s4 = 0.0;
  double delta = 0.0;

  double *data = NULL;
  size_t data_cap = 0;
  double bin_width = 0.0;

  double prec = 1;
  ST_RANGE *r;

  FILE *fs = tmpfile();
  if (fs == NULL)
    exit(EXIT_FAILURE);

  /* read from file */
  while (fgets(buf, BUF, fp) != NULL) {
    if (buf[0] != '#' && (out_type = digitstr(buf))) {
      out_format = MAX(out_format, out_type);

      /* 測定精度 */
      prec = check_precision(buf);
      prec_level = MIN(prec_level, prec);
      x = atof(buf);

      if ((size_t)n >= data_cap) {
        size_t new_cap = data_cap ? data_cap * 2 : 1024;
        double *tmp = realloc(data, new_cap * sizeof(double));
        if (tmp == NULL) {
          free(data);
          fclose(fs);
          exit(EXIT_FAILURE);
        }
        data = tmp;
        data_cap = new_cap;
      }

      data[n] = x;

      fprintf(fs, "%f\n", x);
      if (n == 0)
        min = max = x; /* init variable min and max */
      max = MAX(max, x);
      min = MIN(min, x);
      n++;                       /* N */
      x -= s1;                   /* gap with temp average */
      s1 += x / n;               /* average */
      s2 += (n - 1) * x * x / n; /* square  */
    } else if (!err_f) {
      for (i = 0; buf[i] != '\0' && buf[i] != '\n'; i++)
        p[i] = buf[i];
      k++;
      fprintf(stderr, "%3d: error? %s in line %d\n", k, p, n);
    }
  }

  /* number of bins/ranks */
  switch (rank_type) {
  case RANK_STURGES:
    rank = sturges(n);
    break;

  case RANK_SCOTT:
    rank = scott(n, sdev(n, s2));
    break;

  case RANK_FD:
    bin_width = fd_bin_width(data, (size_t)n);
    if (bin_width <= 0.0)
      rank = sturges(n);
    else
      rank = calculate_bins(min, max, bin_width);
    break;

  default:
    rank = sturges(n);
    break;
  }

  /* calc rank width: delta */
  if (int_f) {
    delta = (int)((max - min) / rank);
    rank = (int)((max - min) / delta) + 1;
  } else
    delta = (max - min) / rank;

  rewind(fs);

  /* get memory for rank distribution */
  if ((r = (ST_RANGE *)malloc((rank + 1) * sizeof(ST_RANGE))) == NULL)
    exit(EXIT_FAILURE);

  /* initialize */
  for (i = 0; i < rank; i++) {
    r[i].cnt = 0;
    if (int_f)
      r[i].range = (int)(min + i * delta);
    else
      r[i].range = min + i * delta;
  }

  /* read from tmpfile for rank distribution */
  while (fgets(buf, BUF, fs) != NULL) {
    x = atof(buf);
    s3 += (x - s1) * (x - s1) * (x - s1);            /* 歪度計算 */
    s4 += (x - s1) * (x - s1) * (x - s1) * (x - s1); /* 尖度計算 */
    i = (int)((x - min) / delta);                    /* 階級へふりわけ */
    if (i == rank) /* 指定階級数と一致した場合には,一段下の最大階級へ */
      i--;
    r[i].cnt++;
    rmax = MAX(rmax, r[i].cnt);
  }
  fclose(fs); /* end of tmpfile */

  if (with_stat) {
    printf("N=%d ", n);
    printf("min=%.3f ", min);
    printf("mean=%.3f ", s1);
    printf("max=%.3f ", max);
    printf("SD=%.3f ", sdev(n, s2));
    printf("SE=%.3f ", serr(n, s2));
    printf("CV=%.3f%% ", cvar(n, s1, s2));
    printf("\n");
    printf("Reliable interval low - upper = %0.3f - %0.3f (Mizutani/Chebyshev; k=2.6)\n",
           rlo(n, s1, s2, 2.6), rup(n, s1, s2, 2.6));
    printf("skew=%.3f ", skew(n, s2, s3));
    printf("kurtosis=%.3f ", kurt(n, s2, s4));
    printf("rank=%d ", rank);
    printf("delta=%.3f ", delta);
    printf("accuracy=%f\n", prec_level);
  }

  /* output */
  for (i = 0; i < rank; i++) {
    if (print_rank)
      printf("%3d ", i + 1);

    if (out_format > 1 || !int_f)
      printf("%8.2f ", r[i].range);
    else
      printf("%5d ", (int)r[i].range);

    if (!only_num) {
      printf("- :");
      for (k = 0; k < r[i].cnt / (rmax / screen_w); k++)
        putchar('*');
    }
    printf(" %d\n", r[i].cnt);
  }

  if (append_return)
    printf("\n");

  free(r);
  return 0;
}

int main(int argc, char **argv) {
  FILE *fp;
  char *op;
  int ai;

  for (ai = 1; ai < argc && argv[ai][0] == '-' && argv[ai][1]; ai++) {
    for (op = &argv[ai][1]; *op; op++) {
      switch (*op) {
      case 'a':
        append_return = 1;
        break;
      case 'c':
        only_num = 1;
        break;
      case 'e':
        err_f = 1;
        break;
      case 'i':
        int_f = 1;
        break;
      case 's':
        with_stat = 1;
        break;
      case 'r':
        print_rank = 1;
        break;
      case 'v':
        fprintf(stderr, "%s\n", REVISION);
        exit(0);
      case 'h':
        fprintf(stderr, "%s\n", USAGE);
        exit(0);

      case 'T':
        if (!digitstr(op + 1)) {
          fprintf(stderr, "%s\n", USAGE);
          exit(0);
        }
        rank_type = atoi(op + 1);
        while (*(++op + 1))
          ;
        break;

      case 'W':
        if (!digitstr(op + 1)) {
          fprintf(stderr, "%s\n", USAGE);
          exit(0);
        }
        screen_w = atoi(op + 1);
        while (*(++op + 1))
          ;
        break;

      default:
        fprintf(stderr, "%% %s: invalid option -%c\n", PROG_NAME, *op);
        fprintf(stderr, "%s\n", USAGE);
        exit(0);
      }
    }
  }

  if (ai >= argc) {
    count_freq(stdin);
  } else {
    for (; ai < argc; ai++) {
      if ((fp = fopen(argv[ai], "r")) == NULL) {
        fprintf(stderr, "can't open %s!\n", argv[ai]);
        exit(1);
      }
      count_freq(fp);
      fclose(fp);
    }
  }

  return 0;
}
