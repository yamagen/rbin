#ifndef FORMULA_H
#define FORMULA_H

#include <stddef.h>

int sturges(double count);
int sturges_two(double count);

double scott_bin_width(size_t n, double s);

double fd_bin_width(double *data, size_t n);
int calculate_bins(double min, double max, double bin_width);

double jarque_bera(int n, double s2, double s3, double s4);
double jarque_bera_pvalue(double jb);

double cvar(int n, double s1, double s2);
double kurt(int n, double s2, double s4);
double sdev(int n, double s2);
double serr(int n, double s2);
double skew(int n, double s2, double s3);
double rlo(int n, double s1, double s2, double k);
double rup(int n, double s1, double s2, double k);
double tvar(int n, double s1, double s2);

#endif
