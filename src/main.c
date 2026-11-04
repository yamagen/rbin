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
  int n;
  double min;
  double mean;
  double max;
  double sd;
  double se;
  double cv;

  double reliable_low;
  double reliable_upper;

  double skew;
  double kurtosis;
  double jb;
  double jb_p;

  int rank;
  double delta;
  double accuracy;

  int method;
  double bin_width;
} ST_OUTPUT;

typedef struct {
  int n;
  double min;
  double max;
  double s1;
  double s2;
} ST_STAT;

void stat_json(ST_OUTPUT *out);
void stat_vertical(ST_OUTPUT *out);
void stat_default(ST_OUTPUT *out);
void prep_stat(ST_OUTPUT *out, ST_STAT *st, int rank, double delta,
               double bin_width, double s3, double s4);
void stat(ST_STAT *st, int rank, double delta, double bin_width, double s3,
          double s4);
int plot(FILE *fs, ST_STAT *st, int rank, double delta);
void bin2delta_rank(double min, double max, double bin_width, double *delta,
                    int *rank);
double scott(size_t n, double s2);
double fd(FILE *fs, size_t n);
double count_bins(FILE *fs, ST_RANGE *r, int rank, double min, double delta);
void calc_skew_kurt(FILE *fs, double mean, double *s3, double *s4);
void stat_add(ST_STAT *st, double x);
int read_data(FILE *fp, FILE *fs, ST_STAT *st);
static int parse_number(const char *s, double *value);
int process(FILE *fp);
int main(int argc, char **argv);

int screen_w = 60;  /* screen width default */
int blank_line = 0; /* print blank line after each file */
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

#define STAT_DEFAULT 0
#define STAT_VERTICAL 1
#define STAT_JSON 2
int stat_format = STAT_DEFAULT;

void stat_json(ST_OUTPUT *out) {
  printf("{\n");
  printf("  \"N\": %d,\n", out->n);
  printf("  \"min\": %.3f,\n", out->min);
  printf("  \"mean\": %.3f,\n", out->mean);
  printf("  \"max\": %.3f,\n", out->max);
  printf("  \"SD\": %.3f,\n", out->sd);
  printf("  \"SE\": %.3f,\n", out->se);
  printf("  \"CV\": %.3f,\n", out->cv);
  printf("  \"reliable_low\": %.3f,\n", out->reliable_low);
  printf("  \"reliable_upper\": %.3f,\n", out->reliable_upper);
  printf("  \"skew\": %.3f,\n", out->skew);
  printf("  \"kurtosis\": %.3f,\n", out->kurtosis);

  if (out->jb >= 0.0) {
    printf("  \"JB\": %.3f,\n", out->jb);
    printf("  \"JB_p\": %.4f,\n", out->jb_p);
    printf("  \"normality\": \"%s\",\n",
           out->jb_p < 0.05 ? "rejected" : "not-rejected");
  } else {
    printf("  \"JB\": null,\n");
    printf("  \"JB_p\": null,\n");
    printf("  \"normality\": null,\n");
  }

  printf("  \"rank\": %d,\n", out->rank);
  printf("  \"delta\": %.3f,\n", out->delta);
  printf("  \"accuracy\": %f,\n", out->accuracy);

  if (out->method == RANK_SCOTT || out->method == RANK_FD) {
    printf("  \"method\": \"%s\",\n",
           out->method == RANK_SCOTT ? "Scott" : "Freedman-Diaconis");
    printf("  \"method_bin_width\": %.3f\n", out->bin_width);
  } else {
    printf("  \"method\": \"Sturges\"\n");
  }

  printf("}\n");
}

void prep_stat(ST_OUTPUT *out, ST_STAT *st, int rank, double delta,
               double bin_width, double s3, double s4) {
  out->n = st->n;
  out->min = st->min;
  out->mean = st->s1;
  out->max = st->max;

  out->sd = sdev(st->n, st->s2);
  out->se = serr(st->n, st->s2);
  out->cv = cvar(st->n, st->s1, st->s2);

  out->reliable_low = rlo(st->n, st->s1, st->s2, 2.6);
  out->reliable_upper = rup(st->n, st->s1, st->s2, 2.6);

  out->skew = skew(st->n, st->s2, s3);
  out->kurtosis = kurt(st->n, st->s2, s4);

  out->jb = jarque_bera(st->n, st->s2, s3, s4);
  out->jb_p = jarque_bera_pvalue(out->jb);

  out->rank = rank;
  out->delta = delta;
  out->accuracy = prec_level;

  out->method = rank_type;
  out->bin_width = bin_width;
}

double scott(size_t n, double s2) { return scott_bin_width(n, sdev(n, s2)); }

double fd(FILE *fs, size_t n) {
  double iqr = calc_iqr(fs, n);

  return fd_bin_width(n, iqr);
}

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

int read_data(FILE *fp, FILE *fs, ST_STAT *st) {

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
  }

  if (st->n == 0) {
    fprintf(stderr, "%% %s: no numeric input\n", PROG_NAME);
    return EXIT_FAILURE;
  }

  return 0;
}

double count_bins(FILE *fs, ST_RANGE *r, int rank, double min, double delta) {
  char buf[BUF];
  double x;
  double rmax = 0.0;
  int i;

  rewind(fs);

  while (fgets(buf, BUF, fs) != NULL) {
    x = atof(buf);

    i = (int)((x - min) / delta);
    if (i < 0)
      i = 0;
    else if (i >= rank)
      i = rank - 1;

    r[i].cnt++;
    rmax = MAX(rmax, r[i].cnt);
  }

  return rmax;
}

void calc_skew_kurt(FILE *fs, double mean, double *s3, double *s4) {
  char buf[BUF];
  double x;

  *s3 = 0.0;
  *s4 = 0.0;

  rewind(fs);

  while (fgets(buf, BUF, fs) != NULL) {
    x = atof(buf);
    *s3 += (x - mean) * (x - mean) * (x - mean);
    *s4 += (x - mean) * (x - mean) * (x - mean) * (x - mean);
  }
}

void bin2delta_rank(double min, double max, double bin_width, double *delta,
                    int *rank) {
  if (max == min) {
    *rank = 1;
    *delta = 1.0;
  } else {
    *delta = bin_width;

    if (int_f) {
      *delta = (int)*delta;
      if (*delta < 1.0)
        *delta = 1.0;
    }

    *rank = calculate_bins(min, max, *delta);
  }
}

int plot(FILE *fs, ST_STAT *st, int rank, double delta) {
  ST_RANGE *r;
  int i, k;
  double rmax;

  if ((r = malloc((rank + 1) * sizeof(ST_RANGE))) == NULL)
    return EXIT_FAILURE;

  /* initialize */
  for (i = 0; i < rank; i++) {
    r[i].cnt = 0;

    if (int_f)
      r[i].range = (int)(st->min + i * delta);
    else
      r[i].range = st->min + i * delta;
  }

  rmax = count_bins(fs, r, rank, st->min, delta);

  /* output */
  {
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

  free(r);

  return EXIT_SUCCESS;
}

void stat(ST_STAT *st, int rank, double delta, double bin_width, double s3,
          double s4) {
  ST_OUTPUT out = {0};

  prep_stat(&out, st, rank, delta, bin_width, s3, s4);

  switch (stat_format) {
  case STAT_VERTICAL:
    stat_vertical(&out);
    break;

  case STAT_JSON:
    stat_json(&out);
    break;

  default:
    stat_default(&out);
    break;
  }
}

void stat_vertical(ST_OUTPUT *out) {

  printf("N\t%d\n", out->n);
  printf("min\t%.3f\n", out->min);
  printf("mean\t%.3f\n", out->mean);
  printf("max\t%.3f\n", out->max);
  printf("SD\t%.3f\n", out->sd);
  printf("SE\t%.3f\n", out->se);
  printf("CV\t%.3f%%\n", out->cv);

  printf("reliable_low\t%.3f\n", out->reliable_low);
  printf("reliable_upper\t%.3f\n", out->reliable_upper);

  printf("skew\t%.3f\n", out->skew);
  printf("kurtosis\t%.3f\n", out->kurtosis);

  if (out->jb >= 0.0) {
    printf("JB\t%.3f\n", out->jb);
    printf("JB_p\t%.4f\n", out->jb_p);
    printf("normality\t%s\n", out->jb_p < 0.05 ? "rejected" : "not-rejected");
  } else {
    printf("JB\tNA\n");
    printf("JB_p\tNA\n");
    printf("normality\tNA\n");
  }

  printf("rank\t%d\n", out->rank);
  printf("delta\t%.3f\n", out->delta);
  printf("accuracy\t%f\n", out->accuracy);

  if (out->method == RANK_SCOTT || out->method == RANK_FD) {
    printf("method\t%s\n",
           out->method == RANK_SCOTT ? "Scott" : "Freedman-Diaconis");
    printf("method_bin_width\t%.3f\n", out->bin_width);
  } else {
    printf("method\tSturges\n");
  }
}

void stat_default(ST_OUTPUT *out) {

  printf("N=%d ", out->n);
  printf("min=%.3f ", out->min);
  printf("mean=%.3f ", out->mean);
  printf("max=%.3f ", out->max);
  printf("SD=%.3f ", out->sd);
  printf("SE=%.3f ", out->se);
  printf("CV=%.3f%% ", out->cv);
  printf("\n");

  printf("Reliable interval low - upper = %0.3f - %0.3f "
         "(Mizutani/Chebyshev; k=2.6)\n",
         out->reliable_low, out->reliable_upper);

  printf("skew=%.3f ", out->skew);
  printf("kurtosis=%.3f ", out->kurtosis);

  if (out->jb >= 0.0) {
    printf("JB=%.3f ", out->jb);
    printf("JB_p=%.4f ", out->jb_p);

    if (out->jb_p < 0.05)
      printf("normality=rejected ");
    else
      printf("normality=not-rejected ");
  } else {
    printf("JB=NA JB_p=NA normality=NA ");
  }

  printf("\n");
  printf("rank=%d ", out->rank);
  printf("delta=%.3f ", out->delta);
  printf("accuracy=%f ", out->accuracy);
  printf("\n");

  if (out->method == RANK_SCOTT || out->method == RANK_FD) {
    printf("method=%s ",
           out->method == RANK_SCOTT ? "Scott" : "Freedman-Diaconis");
    printf("method_bin_width=%.3f ", out->bin_width);
  } else {
    printf("method=Sturges ");
  }

  printf("\n");
}

int process(FILE *fp) {

  ST_STAT st = {0};

  int rank = 0;
  double s3 = 0.0;
  double s4 = 0.0;
  double delta = 0.0;

  double bin_width = 0.0;

  FILE *fs = tmpfile();
  if (fs == NULL)
    exit(EXIT_FAILURE);

  if (read_data(fp, fs, &st) != 0) {
    fclose(fs);
    return EXIT_FAILURE;
  }

  /* number of bins/ranks */
  switch (rank_type) {
  case RANK_STURGES:
    bin_width = sturges(st.n, st.min, st.max);
    break;

  case RANK_SCOTT:
    bin_width = scott((size_t)st.n, st.s2);
    break;

  case RANK_FD:
    bin_width = fd(fs, (size_t)st.n);
    break;

  default:
    break;
  }

  if (bin_width <= 0.0)
    bin_width = sturges(st.n, st.min, st.max);

  /* calc rank width: delta */
  bin2delta_rank(st.min, st.max, bin_width, &delta, &rank);

  if (with_stat)
    calc_skew_kurt(fs, st.s1, &s3, &s4);

  /* Everything concerning r exists only while plotting. */
  if (!no_plot) {
    if (plot(fs, &st, rank, delta) != EXIT_SUCCESS) {
      fclose(fs);
      return EXIT_FAILURE;
    }
  }

  fclose(fs); /* end of tmpfile */

  if (with_stat)
    stat(&st, rank, delta, bin_width, s3, s4);

  if (blank_line)
    printf("\n");

  return 0;
}

int main(int argc, char **argv) {
  FILE *fp;
  char *op;

  argc--;
  argv++;

  while (*argv != NULL && **argv == '-' && *(*argv + 1)) {
    for (op = *argv + 1; *op; op++) {
      switch (*op) {
      case 'b':
        blank_line = 1;
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
      case 'j':
        stat_format = STAT_JSON;
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
      case 'V':
        stat_format = STAT_VERTICAL;
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
    return process(stdin);

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
