#include "config.h"
#include "formula.h"
#include "util.h"
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> /* EXIT_FAILURE,exit */
#include <string.h>

static int parse_number(const char *s, double *value);
int count_freq(FILE *fp);
int main(int argc, char **argv);

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


static int parse_number(const char *s, double *value) {
  char *end;

  while (isspace((unsigned char)*s))
    s++;

  if (*s == '\0')
    return 0;

  errno = 0;
  *value = strtod(s, &end);

  if (s == end || errno == ERANGE || !isfinite(*value))
    return 0;

  while (isspace((unsigned char)*end))
    end++;

  return *end == '\0';
}

int count_freq(FILE *fp) {

  char buf[BUF];
  int n = 0;
  int i = 0, k = 0;
  int rank = 0;
  long line_no = 0;
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
    line_no++;

    if (!parse_number(buf, &x)) {
      size_t len = strcspn(buf, "\r\n");
      buf[len] = '\0';

      if (err_f)
        continue;

      fprintf(stderr, "%% %s: input format error at line %ld: %s\n",
              PROG_NAME, line_no, buf);
      free(data);
      fclose(fs);
      return EXIT_FAILURE;
    }

    out_type = digitstr(buf);
    if (out_type == 0)
      out_type = 2; /* strtod() accepted a real number such as 1e-3 */
    out_format = MAX(out_format, out_type);

    /* 測定精度 */
    prec = check_precision(buf);
    prec_level = MIN(prec_level, prec);

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

    fprintf(fs, "%.17g\n", x);
    if (n == 0)
      min = max = x; /* init variable min and max */
    max = MAX(max, x);
    min = MIN(min, x);
    n++;                       /* N */
    x -= s1;                   /* gap with temp average */
    s1 += x / n;               /* average */
    s2 += (n - 1) * x * x / n; /* square  */
  }

  if (n == 0) {
    fprintf(stderr, "%% %s: no numeric input\n", PROG_NAME);
    free(data);
    fclose(fs);
    return EXIT_FAILURE;
  }

  /* number of bins/ranks */
  switch (rank_type) {
  case RANK_STURGES:
    rank = sturges(n);
    break;

  case RANK_SCOTT:
    bin_width = scott_bin_width((size_t)n, sdev(n, s2));
    if (bin_width <= 0.0)
      rank = sturges(n);
    else
      rank = calculate_bins(min, max, bin_width);
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

  if (rank < 1)
    rank = 1;

  /* calc rank width: delta */
  if (max == min) {
    rank = 1;
    delta = 1.0;
  } else if (int_f) {
    delta = (int)((max - min) / rank);
    if (delta < 1.0)
      delta = 1.0;
    rank = (int)((max - min) / delta) + 1;
  } else {
    delta = (max - min) / rank;
  }

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
    i = (int)((x - min) / delta); /* 階級へふりわけ */
    if (i < 0)
      i = 0;
    else if (i >= rank)
      i = rank - 1;
    r[i].cnt++;
    rmax = MAX(rmax, r[i].cnt);
  }
  fclose(fs); /* end of tmpfile */

  if (with_stat) {
    double jb;
    double jbp;

    jb = jarque_bera(n, s2, s3, s4);
    jbp = jarque_bera_pvalue(jb);

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

    if (jb >= 0.0) {
      printf("JB=%.3f ", jb);
      printf("JB_p=%.4f ", jbp);

      if (jbp < 0.05)
        printf("normality=rejected ");
      else
        printf("normality=not-rejected ");
    } else {
      printf("JB=NA JB_p=NA normality=NA ");
    }

    printf("\n");
    printf("rank=%d ", rank);
    printf("delta=%.3f ", delta);
    printf("accuracy=%f ", prec_level);
    printf("\n");

    if (rank_type == RANK_SCOTT || rank_type == RANK_FD) {
      printf("method=%s ", rank_type == RANK_SCOTT ? "Scott" : "Freedman-Diaconis");
      printf("method_bin_width=%.3f ", bin_width);

      if (bin_width <= 0.0)
        printf("fallback=Sturges ");

    } else {
      printf("method=Sturges ");
    }
    printf("\n");
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
    return count_freq(stdin);
  }

  for (; ai < argc; ai++) {
    if ((fp = fopen(argv[ai], "r")) == NULL) {
      fprintf(stderr, "can't open %s!\n", argv[ai]);
      return EXIT_FAILURE;
    }
    if (count_freq(fp) != 0) {
      fclose(fp);
      return EXIT_FAILURE;
    }
    fclose(fp);
  }

  return EXIT_SUCCESS;
}
