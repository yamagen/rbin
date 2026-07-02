#include "util.h"
#include "config.h"
#include <ctype.h>
#include <stdlib.h>

double check_precision(char *str) {
  double precision_level = 1;
  int decimal_point = FALSE;

  for (; *str != '\0'; str++) {
    if (decimal_point && isdigit(*str)) {
      precision_level *= 0.1;
    }
    if (*str == '.')
      decimal_point = TRUE;
  }
  return precision_level;
}

int digitstr(char *str) {
  int sign_num = 0;
  int float_num = 0;

  if (*str == '\0')
    return FALSE;
  for (; *str != '\0'; str++) {
    /* if (!isdigit(*str) && *str !='\n' && *str !='.' && *str !='-')
       return FALSE;*/
    if (isalpha(*str))
      return FALSE;
    if (*str == '-' && sign_num == 0)
      sign_num = 1;
    if (*str == '.' && float_num == 0)
      float_num = 2;
  }
  return TRUE + sign_num + float_num;
}

int cfnc(const void **a, const void **b) { return strcmp((char *)*a, (char *)*b); }

int dstrlen(char *str) {

  int n;

  for (n = 0; *str != '\0'; n++, str++)
    if (isdelimiter(*str))
      return (n);
  return (n); /* maybe not reached */
}

int isdelimiter(unsigned char c) {

  char *d;

  for (d = delimiter; *d != '\0'; d++)
    if (*d == c)
      return (TRUE);
  return (FALSE);
}

char *get_nth_field(int n, char *str) {
  char *p, *l;
  char *field;
  int f;
  size_t len;

  if (n) {
    f = 0;
  } else {
    f = 1;
    n = 999;
    /*
     * There are no lines which include
     * such a lot of fields, I would say.
     */
  }

  l = p = str;

  while (isdelimiter(*p))
    p++; /* skip delimiters */

  for (;;) {
    if (--n < 1)
      break;

    l = p;
    p += dstrlen(p);

    while (isdelimiter(*p))
      p++; /* skip delimiters */

    if (*p == '\0') {
      if (f) {
        p = l;
        break;
      } else {
        return NULL;
      }
    }
  }

  len = dstrlen(p);

  if ((field = malloc(len + 1)) == NULL) {
    fprintf(stderr, "malloc fail\n");
    exit(1);
  }

  memcpy(field, p, len);
  field[len] = '\0';

  return field; /* free() needed after using the pointer */
}

/* swap: v[i]: v[j] */
void swap(char *v[], int i, int j) {

  char *temp;

  temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}
