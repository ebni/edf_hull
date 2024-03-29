#ifndef _TS_LIB_H
#define _TS_LIB_H
#include <stdio.h>

/* TASK SET DATA STRUCTURE */

/*
 * Task set  related data. It  is preferred to store  homogeneous data
 * (periods,  deadline,...) in  unique vector,  rather than  storing a
 * vector  of  "task  structures",  where each  element  contains  the
 * parameters. This is for efficiency.
 */
typedef struct {
  int num;             /* number of tasks */
  double *per;         /* periods */
  double *dl;          /* rel. deadline (even > period) */
  double *phi;         /* offsets (in [0,period)) */
  int has_phi;         /* 1 if offset; 0 if not */
  double h_per;        /* hyperperiod = lcm of periods */
  double h_per_tol;    /* h_per tolerance. If 0, exact */
  double max_d, max_o; /* max deadline and offset */
  double eps;          /* (t_1-t_0)/t_1 < eps => t_1, t_0 are same */
} ts_t;

/* Enumerated types for the methods of period and deadline generation*/
typedef enum {
  per_null = 0, /* no method */
  per_unif      /* integer periods uniform between min and max */
} per_method_t;

typedef enum {
  dl_null = 0, /* no method*/
  dl_unif      /* double deadlines uniformly distributed around period average
                  with variance dl_var */
} dl_method_t;

/* Structure reserved for values attributed to options when the random mode is
 * activated*/
typedef struct {
  per_method_t per_m; /* method for period generation */
  dl_method_t dl_m;   /*method for deadline generation*/
  int num;            /* number of tasks to be generated */
  unsigned int seed;  /* seed for random num generator */
  double per_min;     /* stored in double, may be integers */
  double per_max;     /* stored in double, may be integers */
  double norm_dl_avg; /* normalized deadline: average */
  double norm_dl_var; /* norm_dl is in [avg-var,avg+var] */
  int phasing;        /* 1 if offset; 0 if not */
  double eps;         /* (t_1-t_0)/t_1 < eps => t_1, t_0 are same */
  int n_repeat;       /*number of tasks set randomly generated given the same
                         specifications*/
  FILE *fp_stats_csv; /*file pointer for the csv file listing at
                      each row the task set specifications and the
                      corresponding minimal number of constraints */
} ts_rand_t;

/*
 * Set the initial values to zero WITHOUT allocating
 */
void ts_set_zero(ts_t *cur_ts);

/*
 * Generate (and allocate if needed) a random set of num tasks based
 * on the settings in (ts_rand_t*)settings.
 */
void ts_rand(ts_t *cur_ts, const ts_rand_t *settings);

/*
 * Allocate arrays for task set. The size is taken from
 * cur_ts->number. If cur_ts->number == 0, then the pointers MUST be
 * NULL
 */
void ts_realloc(ts_t *cur_ts);

/*
 * Initialize the task set using data read from stdin
 */
void ts_read_alloc(ts_t *cur_ts, char *input_filename);

/*
 * Print the data of the task set
 */
void ts_print(const ts_t *cur_ts);

/*
 * Free all the dynamically allocated memory
 */
void ts_free(ts_t *cur_ts);

#endif /* _TS_LIB_H */
