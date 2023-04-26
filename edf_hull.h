#ifndef _EDF_HULL_H
#define _EDF_HULL_H

#include <strings.h>

#include "modules/qhull/src/libqhull_r/libqhull_r.h"
#include "modules/qhull/src/libqhull_r/mem_r.h"
#include "modules/qhull/src/libqhull_r/qhull_ra.h"
#include "ts_lib.h" /* task set basic routines */

/* Qhull global data structure necessary to cut points */
extern qhT qh_qh;

/*
 * Data structure  storing the  interval info.  For each  interval the
 * test of the demand is performed. The test performed is:
 *
 *                sum_i vec_p[i]*C[i] <= t1[i]-t0[i]
 *
 * where C[i] denotes the computation time of the i-th task.
 */
typedef struct {
  int num_points;   /* number of points */
  int alloc_points; /* number of allocated points in the array */
  int num_tasks;    /* number of tasks */
  double *t0;       /* starting point */
  double *t1;       /* finishing point */
  double *vec_p;    /* array of coefficients */
  double interior;  /* the point interior*[1,1,...,1] is in
           the interior */
  coordT *qh_vec_p; /* array given to qhull */
  int num_sel;      /* number of selected points (by qhull) */
  int *vec_sel;     /* indexes of selected points */
} edf_points_t;

/*
 * Set the initial values to zero WITHOUT allocating
 */
void edf_set_zero(edf_points_t *cur_points);

/*
 * Starting from the task set given in cur_task_set, it returns all
 * the EDF points, for all the pairs t0 < t1, according to the EDF
 * schedulability [Baruah 1990]
 */
void edf_create_points(const ts_t *cur_task_set, edf_points_t *cur_points);

/*
 * Print all EDF points
 */
void edf_print_points(const edf_points_t *cur_points);

/*
 * Prints the selection in the following format
 *
 * t0    t1    eta_1    eta_2    ...    eta_N
 *
 * For each row the test that is necessary (and sufficient) to run
 * is
 *
 * eta_1*C_1 +eta_2*C_2 + ... +eta_N*C_N <= t_1-t_0
 */
void edf_print_constraints_C(const edf_points_t *cur_points);

/*
 * Prints the selection in the following format
 *
 * a_1    a_2    ...    a_N   t0    t1
 *
 * For each row the test that is necessary (and sufficient) to run
 * is
 *
 * a_1*U_1 +a_2*U_2 + ... +a_N*U_N <= t_1-t_0
 */
void edf_print_constraints_U(const ts_t *cur_task_set,
                             const edf_points_t *cur_points);

/*
 * It selects the only necessary points using a geometric
 * reasoning. The algorithm for computing the convex hull is invoked through
 * a function internal to the API .
 */
double edf_qhull_points(edf_points_t *cur_points);

/*
 * Free the data structure storing the points
 */
void edf_free_points(edf_points_t *cur_points);

#endif /* _EDF_HULL_H */
