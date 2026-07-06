#define PROG_NAME "rbin"
#define REVISION PROG_NAME " - Last change: 2026/07/02-08:21:15."
#define AUTHOR "H. Yamamoto, Ph.D. (Institute of Science Tokyo)"

#define VERSION REVISION

#define TRUE 1
#define FALSE 0

#define ALL_DIGIT 1
#define INC_SIGNED 2
#define INC_FLOAT 3
#define INC_FLOAT_SIGNED 4

#define BUF 4096

/* bin/rank type */
#define RANK_STURGES 0
#define RANK_SCOTT 1
#define RANK_FD 2
#define RANK_DEFAULT RANK_STURGES

#define USAGE                                                                  \
  "% " PROG_NAME " [-aschiertv] [-Tn] [-Wnum] file...\n"                       \
  "  -a   append a return code (data for plot)\n"                              \
  "  -c   show only range and frequency count\n"                               \
  "  -e   print no error messages\n"                                           \
  "  -h   print this help\n"                                                   \
  "  -i   round delta to integer value\n"                                      \
  "  -r   print rank\n"                                                        \
  "  -s   print statistics\n"                                                  \
  "  -Tn  set bin/rank rule:\n"                                                \
  "       n=0  Sturges (1926; default): sample-size based, rough overview\n"    \
  "       n=1  Scott (1979): SD-based, sensitive to long tails/outliers\n"      \
  "       n=2  Freedman-Diaconis (1981): IQR-based, robust to outliers\n"       \
  "  -v   print " PROG_NAME " version\n"                                       \
  "  -Wn  set screen width to n (default: n=60)\n"                             \
  "(c) " AUTHOR ": " REVISION "\n"

#define SQU(x) ((x) * (x))                /* Square */
#define CUBE(x) ((x) * (x) * (x))         /* Cube */
#define MIN(x, y) ((x) < (y) ? (x) : (y)) /* Return smallest */
#define MAX(x, y) ((x) > (y) ? (x) : (y)) /* Return largest */
#define DEG(x) ((x) / (PI) * 180.0)       /* Radians -> Degrees */
#define RAD(x) ((x) / 180.0 * (PI))       /* Degrees -> Radians */

typedef struct /* def of range */
{
  double range;
  int cnt;
} ST_RANGE;
