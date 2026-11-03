#include "config.h"
#include "formula.h"
#include "util.h"
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> /* EXIT_FAILURE,exit */
#include <string.h>

typedef struct {
  double q;
  double marker[5];
  double pos[5];
  double desired[5];
  double increment[5];
  double init[5];
  int count;
} ST_QUANTILE;

typedef struct {
  int n;
  double min;
  double max;
  double s1;
  double s2;
} ST_STAT;

void quantile_init(ST_QUANTILE *qs, double q);
void quantile_add(ST_QUANTILE *qs, double x);
double quantile_value(const ST_QUANTILE *qs);

void stat_add(ST_STAT *st, double x);
int read_data(FILE *fp, FILE *fs, ST_STAT *st, ST_QUANTILE *q1,
              ST_QUANTILE *q3);
static int parse_number(const char *s, double *value);
int count_freq(FILE *fp);
int main(int argc, char **argv);

int screen_w = 60; /* screen width default */
int append_return = 0;
int only_num = 0;
int with_stat = 0;
int no_plot = 0;
int out_format = 0;
int print_rank = 0;
int hold_empty = 0;
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

void quantile_init(ST_QUANTILE *qs, double q) {
  memset(qs, 0, sizeof(*qs));
  qs->q = q;
}

void stat_add(ST_STAT *st, double x) {
  if (st->n == 0)
    st->min = st->max = x;

  st->min = MIN(st->min, x);
  st->max = MAX(st->max, x);

  st->n++;
  x -= st->s1;
  st->s1 += x / st->n;
  st->s2 += (st->n - 1) * x * x / st->n;
}

static int compare_double(const void *a, const void *b) {
  double x = *(const double *)a;
  double y = *(const double *)b;

  if (x < y)
    return -1;
  if (x > y)
    return 1;
  return 0;
}

double quantile_value(const ST_QUANTILE *qs) {
  if (qs->count < 5)
    return NAN;

  return qs->marker[2];
}

void quantile_add(ST_QUANTILE *qs, double x) {
  if (qs->count < 5) {
    qs->init[qs->count++] = x;

    if (qs->count == 5) {
      qsort(qs->init, 5, sizeof(double), compare_double);

      for (int i = 0; i < 5; i++)
        qs->marker[i] = qs->init[i];

      qs->pos[0] = 1;
      qs->pos[1] = 2;
      qs->pos[2] = 3;
      qs->pos[3] = 4;
      qs->pos[4] = 5;

      qs->desired[0] = 1;
      qs->desired[1] = 1 + 2 * qs->q;
      qs->desired[2] = 1 + 4 * qs->q;
      qs->desired[3] = 3 + 2 * qs->q;
      qs->desired[4] = 5;

      qs->increment[0] = 0;
      qs->increment[1] = qs->q / 2;
      qs->increment[2] = qs->q;
      qs->increment[3] = (1 + qs->q) / 2;
      qs->increment[4] = 1;
    }

    return;
  }

  int k;

  if (x < qs->marker[0]) {
    qs->marker[0] = x;
    k = 0;
  } else if (x < qs->marker[1]) {
    k = 0;
  } else if (x < qs->marker[2]) {
    k = 1;
  } else if (x < qs->marker[3]) {
    k = 2;
  } else if (x <= qs->marker[4]) {
    k = 3;
  } else {
    qs->marker[4] = x;
    k = 3;
  }

  for (int i = k + 1; i < 5; i++)
    qs->pos[i]++;

  for (int i = 0; i < 5; i++)
    qs->desired[i] += qs->increment[i];

  for (int i = 1; i < 4; i++) {
    double d = qs->desired[i] - qs->pos[i];

    if ((d >= 1 && qs->pos[i + 1] - qs->pos[i] > 1) ||
        (d <= -1 && qs->pos[i - 1] - qs->pos[i] < -1)) {
      double direction = (d > 0) ? 1 : -1;

      double next =
          qs->marker[i] + direction / (qs->pos[i + 1] - qs->pos[i - 1]) *
                              ((qs->pos[i] - qs->pos[i - 1] + direction) *
                                   (qs->marker[i + 1] - qs->marker[i]) /
                                   (qs->pos[i + 1] - qs->pos[i]) +
                               (qs->pos[i + 1] - qs->pos[i] - direction) *
                                   (qs->marker[i] - qs->marker[i - 1]) /
                                   (qs->pos[i] - qs->pos[i - 1]));

      if (qs->marker[i - 1] < next && next < qs->marker[i + 1]) {
        qs->marker[i] = next;
      } else {
        int j = i + (int)direction;

        qs->marker[i] += direction * (qs->marker[j] - qs->marker[i]) /
                         (qs->pos[j] - qs->pos[i]);
      }
      qs->pos[i] += direction;
    }
  }
}

int read_data(FILE *fp, FILE *fs, ST_STAT *st, ST_QUANTILE *q1,
              ST_QUANTILE *q3) {

  char buf[BUF];
  long line_no = 0;
  int out_type = 0;
  double prec = 1;
  double x = 0.0;

  while (fgets(buf, BUF, fp) != NULL) {
    line_no++;

    if (!parse_number(buf, &x)) {
      size_t len = strcspn(buf, "\r\n");
      buf[len] = '\0';

      if (err_f)
        continue;

      fprintf(stderr, "%% %s: input format error at line %ld: %s\n", PROG_NAME,
              line_no, buf);
      return EXIT_FAILURE;
    }

    out_type = digitstr(buf);
    if (out_type == 0)
      out_type = 2;
    out_format = MAX(out_format, out_type);

    prec = check_precision(buf);
    prec_level = MIN(prec_level, prec);

    fprintf(fs, "%.17g\n", x);
    stat_add(st, x);
    quantile_add(q1, x);
    quantile_add(q3, x);
  }

  if (st->n == 0) {
    fprintf(stderr, "%% %s: no numeric input\n", PROG_NAME);
    return EXIT_FAILURE;
  }

  return 0;
}

double fd_bin_width_iqr(size_t n, double iqr) {
  if (n <= 1 || iqr <= 0.0)
    return 0.0;

  return 2.0 * iqr / cbrt((double)n);
}

int count_freq(FILE *fp) {

  ST_STAT st = {0};

  ST_QUANTILE q1;
  ST_QUANTILE q3;

  quantile_init(&q1, 0.25);
  quantile_init(&q3, 0.75);

  char buf[BUF];

  int i = 0, k = 0;
  int rank = 0;
  long line_no = 0;
  int out_type = 0;

  double x = 0.0;
  double rmax = 0.0;
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

      fprintf(stderr, "%% %s: input format error at line %ld: %s\n", PROG_NAME,
              line_no, buf);
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

    if ((size_t)st.n >= data_cap) {
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

    data[st.n] = x;

    fprintf(fs, "%.17g\n", x);

    stat_add(&st, x);
    quantile_add(&q1, x);
    quantile_add(&q3, x);
  }

  if (st.n == 0) {
    fprintf(stderr, "%% %s: no numeric input\n", PROG_NAME);
    free(data);
    fclose(fs);
    return EXIT_FAILURE;
  }

  double sketch_q1 = quantile_value(&q1);
  double sketch_q3 = quantile_value(&q3);
  double sketch_iqr = sketch_q3 - sketch_q1;

  double sketch_bin_width = fd_bin_width_iqr((size_t)st.n, sketch_iqr);
  int sketch_rank = calculate_bins(st.min, st.max, sketch_bin_width);

  /* number of bins/ranks */
  switch (rank_type) {
  case RANK_STURGES:
    rank = sturges(st.n);
    break;

  case RANK_SCOTT:
    bin_width = scott_bin_width((size_t)st.n, sdev(st.n, st.s2));
    if (bin_width <= 0.0)
      rank = sturges(st.n);
    else
      rank = calculate_bins(st.min, st.max, bin_width);
    break;

  case RANK_FD:
    bin_width = fd_bin_width(data, (size_t)st.n);
    if (bin_width <= 0.0)
      rank = sturges(st.n);
    else
      rank = calculate_bins(st.min, st.max, bin_width);
    break;

  default:
    rank = sturges(st.n);
    break;
  }

  if (rank < 1)
    rank = 1;

  /* calc rank width: delta */
  if (st.max == st.min) {
    rank = 1;
    delta = 1.0;
  } else if (int_f) {
    delta = (int)((st.max - st.min) / rank);
    if (delta < 1.0)
      delta = 1.0;
    rank = (int)((st.max - st.min) / delta) + 1;
  } else {
    delta = (st.max - st.min) / rank;
  }

  rewind(fs);

  /* get memory for rank distribution */
  if ((r = (ST_RANGE *)malloc((rank + 1) * sizeof(ST_RANGE))) == NULL)
    exit(EXIT_FAILURE);

  /* initialize */
  for (i = 0; i < rank; i++) {
    r[i].cnt = 0;
    if (int_f)
      r[i].range = (int)(st.min + i * delta);
    else
      r[i].range = st.min + i * delta;
  }

  /* read from tmpfile for rank distribution */
  while (fgets(buf, BUF, fs) != NULL) {
    x = atof(buf);
    s3 += (x - st.s1) * (x - st.s1) * (x - st.s1);               /* 歪度計算 */
    s4 += (x - st.s1) * (x - st.s1) * (x - st.s1) * (x - st.s1); /* 尖度計算 */
    i = (int)((x - st.min) / delta); /* 階級へふりわけ */
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

    jb = jarque_bera(st.n, st.s2, s3, s4);
    jbp = jarque_bera_pvalue(jb);

    printf("N=%d ", st.n);
    printf("min=%.3f ", st.min);
    printf("mean=%.3f ", st.s1);
    printf("max=%.3f ", st.max);
    printf("SD=%.3f ", sdev(st.n, st.s2));
    printf("SE=%.3f ", serr(st.n, st.s2));
    printf("CV=%.3f%% ", cvar(st.n, st.s1, st.s2));
    printf("\n");

    printf("Reliable interval low - upper = %0.3f - %0.3f (Mizutani/Chebyshev; "
           "k=2.6)\n",
           rlo(st.n, st.s1, st.s2, 2.6), rup(st.n, st.s1, st.s2, 2.6));

    printf("skew=%.3f ", skew(st.n, st.s2, s3));
    printf("kurtosis=%.3f ", kurt(st.n, st.s2, s4));

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

    printf("sketch_Q1=%.6f sketch_Q3=%.6f sketch_IQR=%.6f "
           "sketch_bin_width=%.6f sketch_rank=%d\n",
           sketch_q1, sketch_q3, sketch_iqr, sketch_bin_width, sketch_rank);

    printf("\n");
    printf("rank=%d ", rank);
    printf("delta=%.3f ", delta);
    printf("accuracy=%f ", prec_level);
    printf("\n");

    if (rank_type == RANK_SCOTT || rank_type == RANK_FD) {
      printf("method=%s ",
             rank_type == RANK_SCOTT ? "Scott" : "Freedman-Diaconis");
      printf("method_bin_width=%.3f ", bin_width);

      if (bin_width <= 0.0)
        printf("fallback=Sturges ");

    } else {
      printf("method=Sturges ");
    }
    printf("\n");
  }

  /* output */
  if (!no_plot) {
    int previous_count = 0;

    for (i = 0; i < rank; i++) {
      int display_count = r[i].cnt;

      if (hold_empty && display_count == 0 && i > 0)
        display_count = previous_count;

      if (print_rank)
        printf("%3d ", i + 1);

      if (out_format > 1 || !int_f)
        printf("%8.2f ", r[i].range);
      else
        printf("%5d ", (int)r[i].range);

      if (!only_num) {
        printf("- :");
        for (k = 0; k < display_count / (rmax / screen_w); k++)
          putchar('*');
      }
      printf(" %d\n", display_count);
      previous_count = display_count;
    }
  }

  if (append_return)
    printf("\n");

  free(r);
  return 0;
}

int process(FILE *fp) { return count_freq(fp); }

int main(int argc, char **argv) {
  FILE *fp;
  char *op;

  argc--;
  argv++;

  while (*argv != NULL && **argv == '-' && *(*argv + 1)) {
    for (op = *argv + 1; *op; op++) {
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
      case 'H':
        hold_empty = 1;
        break;
      case 'i':
        int_f = 1;
        break;
      case 'n':
        no_plot = 1;
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

    argc--;
    argv++;
  }

  if (argc == 0)
    return process(stdin);

  if (*argv == NULL)
    return count_freq(stdin);

  for (; *argv != NULL; argv++) {
    if ((fp = fopen(*argv, "r")) == NULL) {
      fprintf(stderr, "can't open %s!\n", *argv);
      return EXIT_FAILURE;
    }

    if (process(fp) != 0) {
      fclose(fp);
      return EXIT_FAILURE;
    }

    fclose(fp);
  }

  return EXIT_SUCCESS;
}
