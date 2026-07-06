#ifndef UTIL_H
#define UTIL_H

#define TRUE 1
#define FALSE 0

extern char delimiter[];

void swap(char *v[], int i, int j);
int isdelimiter(unsigned char c);
char *get_nth_field(int n, char *str);
int dstrlen(char *str);

int digitstr(char *str);
double check_precision(char *str);

#endif
