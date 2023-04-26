#include "ts_lib.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Helper macros to compact a bit the code. They are relative to a
 * data structure of type ts_t, with name cur_ts
 */
#define N (cur_ts->num)
#define T(i) (cur_ts->per[(i)])
#define D(i) (cur_ts->dl[(i)])
#define Phi(i) (cur_ts->phi[(i)])
#define IS_PHI (cur_ts->has_phi)
#define H (cur_ts->h_per)
#define EPS (cur_ts->eps)

/*
 * It computes the least common multiple H of all periods (often
 * called hyperperiod). The algorithm guarantees exactness by at least
 * EPS: there is certainly a multiple of any period in the interval
 * [H, H*(1+EPS)). If (cur_ts->h_per_tol == 0) then the hyperperiod is
 * exact.
 */
static void compute_hyperperiod(ts_t *cur_ts) {
  unsigned long *jobs;
  double *multi;
  int id_max = 0, id_min = 0, i;

  /* Number of jobs per task within (current) hyperperiod  */
  jobs = malloc(sizeof(*jobs) * N);
  for (i = 0; i < N; i++)
    jobs[i] = 1;
  /* Array of jobs[i]*T(i) */
  multi = malloc(sizeof(*multi) * N);
  memcpy(multi, cur_ts->per, sizeof(*multi) * N);
  /* searching for min and max */
  for (i = 1; i < N; i++) {
    if (T(i) > T(id_max))
      id_max = i;
    if (T(i) < T(id_min))
      id_min = i;
  }
  while (multi[id_max] / multi[id_min] - 1 > EPS) {
    if (++jobs[id_min] >= ULONG_MAX) {
      fprintf(stderr, "%s:%d: Exceeded precision\n", __FILE__, __LINE__);
      H = 0;
      return;
    }
    multi[id_min] = jobs[id_min] * T(id_min);
    if (multi[id_min] > multi[id_max])
      id_max = id_min;
    for (i = 0; i < N; i++) {
      if (multi[i] < multi[id_min])
        id_min = i;
    }
  }
  H = multi[id_min];
  cur_ts->h_per_tol = multi[id_max] / multi[id_min] - 1;
  free(jobs);
  free(multi);
}

/*
 * Compute the maximum deadline and offset
 */
static void store_max_dl_phi(ts_t *cur_ts) {
  int i;

  /* Computing max deadline and max offset */
  cur_ts->max_d = D(0);
  cur_ts->max_o = Phi(0);
  for (i = 1; i < N; i++) {
    if (D(i) > cur_ts->max_d)
      cur_ts->max_d = D(i);
    if (Phi(i) > cur_ts->max_o)
      cur_ts->max_o = Phi(i);
  }
}

void ts_set_zero(ts_t *cur_ts) {
  cur_ts->num = 0;
  cur_ts->per = NULL;
  cur_ts->dl = NULL;
  cur_ts->phi = NULL;
}

void ts_realloc(ts_t *cur_ts) {
  /* allocate the data structure */
  cur_ts->per = (double *)realloc(cur_ts->per, sizeof(double) * N);
  cur_ts->dl = (double *)realloc(cur_ts->dl, sizeof(double) * N);
  cur_ts->phi = (double *)realloc(cur_ts->phi, sizeof(double) * N);
}

void ts_read_alloc(ts_t *cur_ts, char *input_filename) {
  int i, num;
  double cur_T, cur_D, cur_O;
  int cur_phasing;

  FILE *file_ptr = fopen(input_filename, "r");
  if (NULL == file_ptr) {
    fprintf(stderr, "%s:%d: Input file can't be opened.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  /* initially set no phasing */
  cur_phasing = 0;

  /* read the number of tasks */
  fscanf(file_ptr, "%d", &num);
  if (N != num) {
    N = num;
    /* Now alloc proper amount */
    ts_realloc(cur_ts);
  } else {
    /*Nothing, already allocated */
  }

  fscanf(file_ptr, "%lf", &(EPS));

  for (i = 0; i < N; i++) {
    /* read the task set parameters */
    fscanf(file_ptr, "%lf\t%lf\t%lf", &cur_T, &cur_D, &cur_O);
    T(i) = cur_T;
    D(i) = cur_D;
    cur_O = cur_O - floor(cur_O / cur_T) * cur_T;
    Phi(i) = cur_O;
    /* update the phasing */
    cur_phasing |= cur_O != 0;
  }

  /* Closing input file once the task set has been read*/
  fclose(file_ptr);

  IS_PHI = cur_phasing;
  compute_hyperperiod(cur_ts);
  store_max_dl_phi(cur_ts);
}

void ts_rand(ts_t *cur_ts, const ts_rand_t *settings) {
  int i;

  switch (settings->per_m) {
  case per_null:
    fprintf(stderr,
            "%s:%d: Null method specified. Please specify a method for period "
            "generation\n",
            __FILE__, __LINE__);
    return;
  case per_unif:
    srand(settings->seed);
    if (N != settings->num) {
      N = settings->num;
      ts_realloc(cur_ts);
    }
    IS_PHI = settings->phasing;
    for (i = 0; i < N; i++) {
      T(i) = floor(settings->per_min +
                   (double)rand() * (settings->per_max - settings->per_min) /
                       (double)RAND_MAX);
      if (IS_PHI) {
        Phi(i) = (double)rand() / (double)RAND_MAX * T(i);
      }
    }
    break;
  default:
    fprintf(stderr, "%s:%d: Unknown method specified for period generation\n",
            __FILE__, __LINE__);
    return;
  }

  /*Generating deadlines */
  switch (settings->dl_m) {
  case dl_null:
    fprintf(stderr,
            "%s:%d: Null method specified. Please specify a method for "
            "deadline generation\n",
            __FILE__, __LINE__);
    return;
  case dl_unif:
    for (i = 0; i < N; i++) {
      D(i) = T(i) *
             (settings->norm_dl_avg - settings->norm_dl_var +
              ((double)rand() / (double)RAND_MAX) * 2 * settings->norm_dl_var);
    }
    break;
  default:
    fprintf(stderr, "%s:%d: Unknown method for deadline generation specified\n",
            __FILE__, __LINE__);
    return;
  }

  compute_hyperperiod(cur_ts);
  store_max_dl_phi(cur_ts);
}

void ts_print(const ts_t *cur_ts) {
  int i;

  printf("\nThe task set is %s phasing\n", IS_PHI ? "WITH" : "WITHOUT");
  printf(" i: (      Ti,       Di,       Oi)\n");
  for (i = 0; i < N; i++) {
    printf("%2i: (%8.1f, %8f, %8f)\n", i, T(i), D(i), Phi(i));
  }
  printf("Hyperperiod: %f, computed with tolerance %g (if == 0, then exact)\n",
         H, cur_ts->h_per_tol);
  fflush(stdout);
}

void ts_free(ts_t *cur_ts) {
  free(cur_ts->per);
  free(cur_ts->dl);
  free(cur_ts->phi);
}
