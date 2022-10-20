#include "edf_hull.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define IND (cur_points->vec_sel[i])

void edf_set_zero(edf_points_t* cur_points) {
    cur_points->alloc_points = 0;
    cur_points->t0 = NULL;
    cur_points->t1 = NULL;
    cur_points->vec_p = NULL;
    cur_points->qh_vec_p = NULL;
    cur_points->num_sel = 0;
    cur_points->vec_sel = NULL;
}

void edf_create_points(const ts_t* cur_task_set, edf_points_t* cur_points) {
#define N (cur_task_set->num)
#define T(i) (cur_task_set->per[i])
#define D(i) (cur_task_set->dl[i])
#define O(i) (cur_task_set->phi[i])

    unsigned long* max_job;
    int i, j, k, reserve_p, s_i, s_j;
    int num_deadline;
    int point_count;
    double cur_t, cur_s, coef, sum_rate;
    double big_enough;

    /* store the number of tasks */
    cur_points->num_tasks = N;

    /* Set the limit until which computing the points */
    if (cur_task_set->has_phi) {
        big_enough = 2 * cur_task_set->h_per + cur_task_set->max_d + cur_task_set->max_o;
    } else {
        big_enough = cur_task_set->h_per + cur_task_set->max_d;
    }

    /*
     * Maximum job of the tasks. Also compute the sum of max_job[i],
     * because this is also the number of considered deadlines.
     */
    max_job = malloc(sizeof(*max_job) * N);
    num_deadline = 0;
    for (i = 0; i < N; i++) {
        max_job[i] = (unsigned long)floor((big_enough - O(i) - D(i)) / T(i));
        num_deadline += max_job[i] + 1;
    }

    /*
     * Start computing the points. Reserve the initial point for a
     * special purpose. The number of reserved points is reserv_p.
     */
    /* Reserve one point for sum_i U_i <= 1, N points for the Ci >= 0 */
    reserve_p = N + 1;

    /* if we have offset the points can be really many */
    if (cur_task_set->has_phi) {
        /* the number of points (overestimating it) */
        cur_points->num_points = reserve_p + num_deadline * num_deadline;

        if (cur_points->num_points > cur_points->alloc_points) {
            /* allocate the t0, t1 data structure for the points */
            cur_points->t0 =
                (double*)realloc(cur_points->t0, sizeof(double) * cur_points->num_points);
            cur_points->t1 =
                (double*)realloc(cur_points->t1, sizeof(double) * cur_points->num_points);
            cur_points->vec_p =
                (double*)realloc(cur_points->vec_p, sizeof(double) * cur_points->num_points * N);
            /* The points to be used are in cur_points->qh_vec_p */
            cur_points->qh_vec_p =
                (coordT*)realloc(cur_points->qh_vec_p, sizeof(coordT) * cur_points->num_points * N);
            cur_points->alloc_points = cur_points->num_points;
        }

        /* initialize the point counter */
        point_count = 0;

        /* store the -Ci <= 0 constraint */
        for (i = 0; i < N; i++) {
            cur_points->t0[point_count] = 0;
            cur_points->t1[point_count] = 0;
            for (j = 0; j < N; j++) {
                cur_points->vec_p[point_count * N + j] = ((i == j) ? -1 : 0);
            }

            /* increment the point counter */
            point_count++;
        }

        /* store the sum Ui <= 1 constraint */
        cur_points->t0[point_count] = 0;
        cur_points->t1[point_count] = 1;
        for (i = 0; i < N; i++) {
            cur_points->vec_p[point_count * N + i] = 1 / T(i);
        }
        point_count++;

        /*
         * compute the points corresponding to the (start
         * time, deadline) pairs
         */
        /* loop on the absolute deadlines */
        for (i = 0; i < N; i++) {
            for (j = 0; j <= max_job[i]; j++) {
                /* t1 is equal to the j-th deadline of task i */
                cur_t = O(i) + j * T(i) + D(i);

                /* loop on the start times */
                for (s_i = 0; s_i < N; s_i++) {
                    for (s_j = 0; s_j <= max_job[s_i]; s_j++) {
                        /* the current start time */
                        cur_s = O(s_i) + s_j * T(s_i);
                        if (cur_s >= cur_t) /* start time is ahead the deadline */
                            break;
                        cur_points->t1[point_count] = cur_t;
                        cur_points->t0[point_count] = cur_s;
                        for (k = 0; k < N; k++) {
                            /* there may be numerical problems... */
                            coef = floor((cur_t - O(k) - D(k)) / T(k)) -
                                   ceil((cur_s - O(k)) / T(k)) + 1;
                            /* must be non-negative */
                            if (coef < 0) coef = 0;
                            /* store it */
                            cur_points->vec_p[point_count * N + k] = coef;
                        }
                        /* increment the point counter */
                        point_count++;
                    }
                }
            }
        }
        /* update the num_points */
        cur_points->num_points = point_count;
    } else { /* no offset here */

        /* the number of points */
        cur_points->num_points = reserve_p + num_deadline;

        if (cur_points->num_points > cur_points->alloc_points) {
            /* allocate the t0, t1 data structure for the points */
            cur_points->t0 =
                (double*)realloc(cur_points->t0, sizeof(double) * cur_points->num_points);
            cur_points->t1 =
                (double*)realloc(cur_points->t1, sizeof(double) * cur_points->num_points);
            cur_points->vec_p =
                (double*)realloc(cur_points->vec_p, sizeof(double) * cur_points->num_points * N);
            /* The points to be used are in cur_points->qh_vec_p */
            cur_points->qh_vec_p =
                (coordT*)realloc(cur_points->qh_vec_p, sizeof(coordT) * cur_points->num_points * N);
            cur_points->alloc_points = cur_points->num_points;
        }

        /* initialize the point counter */
        point_count = 0;

        /* store the -Ci <= 0 constraint */
        for (i = 0; i < N; i++) {
            cur_points->t0[point_count] = 0;
            cur_points->t1[point_count] = 0;
            for (j = 0; j < N; j++) {
                cur_points->vec_p[point_count * N + j] = ((i == j) ? -1 : 0);
            }

            /* increment the point counter */
            point_count++;
        }

        /* store the sum Ui <= 1 constraint */
        cur_points->t0[point_count] = 0;
        cur_points->t1[point_count] = 1;
        for (i = 0; i < N; i++) {
            cur_points->vec_p[point_count * N + i] = 1 / T(i);
        }
        point_count++;

        /* compute the points corresponding to deadlines */
        for (i = 0; i < N; i++) {
            for (j = 0; j <= max_job[i]; j++) {
                /* t0 is zero, when no offset */
                cur_points->t0[point_count] = 0;
                /* t1 is equal to the j-th deadline of task i */
                cur_points->t1[point_count] = cur_t = j * T(i) + D(i);
                for (k = 0; k < N; k++) {
                    /* there may be numerical problems... */
                    coef = floor((cur_t - D(k)) / T(k)) + 1;
                    /* must be non-negative */
                    if (coef < 0) coef = 0;
                    /* store it */
                    cur_points->vec_p[point_count * N + k] = coef;
                }

                /* increment the point counter */
                point_count++;
            }
        }
    }
    /* compute the interior for traslation purpose */
    sum_rate = 0;
    for (i = 0; i < N; i++) {
        sum_rate += 1 / (T(i) < D(i) ? T(i) : D(i));
    }
    cur_points->interior = 1 / (2 * sum_rate);

    /* free the max_job which is only a local data structure */
    free(max_job);
#undef N
#undef T
#undef D
#undef O
}

/*
 * Translate the points so that the origin belongs to the interior.
 * The translated vectors are written to the array
 * cur_points->qh_vec_p
 */
static void trans_points(edf_points_t* cur_points) {
    int i, j;
    double sum_coef, sum_rate, mul_factor;

    /* all the points */
    for (i = 0; i < cur_points->num_points; i++) {
        /* compute the sum of all the coefficients of the
         * current vector */
        sum_coef = 0;
        for (j = 0; j < cur_points->num_tasks; j++) {
            sum_coef += cur_points->vec_p[i * cur_points->num_tasks + j];
        }

        /*
         * the  vector  coefficients  are  multiplied  by  the
         * following factor.  Due to  the proper  selection on
         * the interior  point, this  factor should  ALWAYS be
         * GREATER than zero .
         */
        mul_factor = 1 / (cur_points->t1[i] - cur_points->t0[i] - cur_points->interior * sum_coef);
        for (j = 0; j < cur_points->num_tasks; j++) {
            cur_points->qh_vec_p[i * cur_points->num_tasks + j] =
                cur_points->vec_p[i * cur_points->num_tasks + j] * mul_factor;
        }
    }
}

/*
 * Prints the selection in the following format
 *
 * t0    t1    a_1    a_2    ...    a_N
 *
 * For each row the test that is necessary (and sufficient) to run
 * is
 *
 * a_1*C_1 +a_2*C_2 + ... +a_N*C_N <= t_1-t_0
 */
void edf_print_constraints_C(const edf_points_t* cur_points) {
    int i, j;
    printf("\nMinimal set of constraints are written in the following form\n\n");
    printf("    eta_1*C_1 +eta_2*C_2 + ... +eta_N*C_N <= t_1-t_0\n\n");
    for (j = 0; j < cur_points->num_tasks; j++) {
        printf("eta_%d\t", j + 1);
    }
    printf("t1\tt0\n");
    for (i = 0; i < cur_points->num_sel; i++) {
        for (j = 0; j < cur_points->num_tasks; j++) {
            printf("%.4f\t", cur_points->vec_p[IND * cur_points->num_tasks + j]);
        }
        printf("%.0f\t%.0f\n", cur_points->t1[IND], cur_points->t0[IND]);
    }
}

void edf_print_constraints_U(const ts_t* cur_task_set, const edf_points_t* cur_points) {
    int i, j;
    printf("\nOr, alternatively, minimal set of constraints can also be written as\n\n");
    printf("    eta_1*U_1 +eta_2*U_2 + ... +eta_N*U_N <= t_1-t_0\n\n");
    for (j = 0; j < cur_points->num_tasks; j++) {
        printf("eta_%d\t", j + 1);
    }
    printf("t1\tt0\n");

#define T(j) (cur_task_set->per[j])
    for (i = 0; i < cur_points->num_sel; i++) {
        for (j = 0; j < cur_points->num_tasks; j++) {
            if (i < cur_points->num_tasks) {
                printf("%.4f\t", cur_points->vec_p[IND * cur_points->num_tasks + j] * (1 / T(j)));
            } else {
                printf("%.4f\t", (cur_points->vec_p[IND * cur_points->num_tasks + j] * T(j)) /
                                     cur_points->t1[IND]);
            }
        }
        if (cur_points->t1[IND] != 0)
            printf("%.0f\t%.0f\n", cur_points->t1[IND] / cur_points->t1[IND], cur_points->t0[IND]);
        else
            printf("%.0f\t%.0f\n", cur_points->t1[IND], cur_points->t0[IND]);

#undef T
    }
}

void edf_print_points(const edf_points_t* cur_points) {
    int i, j;

    printf("Number of tasks:\t%d\n", cur_points->num_tasks);
    printf("Number of points:\t%d\n", cur_points->num_points);
    printf("Interior:\t\t%f\n", cur_points->interior);
    printf("Points:\n");
    for (i = 0; i < cur_points->num_points; i++) {
        printf("t0=%.3f\tt1=%.3f", cur_points->t0[i], cur_points->t1[i]);
        for (j = 0; j < cur_points->num_tasks; j++) {
            printf("\t%.3f", cur_points->vec_p[i * cur_points->num_tasks + j]);
        }
        printf("\n");
    }
}

/**
 * Fills array vec_sel with indexes of the selected points. Makes use of
 * a set data structure to sort vertices by id.
 */
void edf_array_indexes(qhT* qh, int* arr) {
    setT *vertices, *points;
    pointT* point;
    vertexT *vertex, **vertexp;
    int id, count_points = 0;
    int numpoints = 0, point_i, point_n;
    int allpoints = qh->num_points + qh_setsize(qh, qh->other_points);

    points = qh_settemp(qh, allpoints);
    qh_setzero(qh, points, 0, allpoints);
    vertices = qh_facetvertices(qh, qh->facet_list, NULL, 0);
    FOREACHvertex_(vertices) {
        id = qh_pointid(qh, vertex->point);
        if (id >= 0) {
            SETelem_(points, id) = vertex->point;
            numpoints++;
        }
    }
    qh_settempfree(qh, &vertices);
    FOREACHpoint_i_(qh, points) {
        if (point) {
            arr[count_points] = point_i;
            count_points++;
        }
    }
    qh_settempfree(qh, &points);
}

/*
 * Uses the API of qhull to calculate the convex hull.
 * Returns the CPU computing time to obtain the output.
 */
double edf_qhull_points(edf_points_t* cur_points) {
#define DIM (cur_points->num_tasks)
#define TOTpoints (cur_points->num_points)
#define SIZEobj (1 >> DIM)
    int dim = DIM; /* dimension of points */
    int numpoints; /* number of points */
    /* coordT points[(DIM + 1) * TOTpoints]; array of coordinates for each point */
    coordT** rows;
    boolT ismalloc = True;  /* True if qhull should free points in qh_freeqhull() or reallocation */
    char flags[250];        /* option flags for qhull, see qh-quick.htm */
    FILE* outfile = NULL;   /* output from qh_produce_output()
                                 use NULL to skip qh_produce_output() */
    FILE* errfile = stderr; /* error messages from qhull code */
    int exitcode;           /* 0 if no error from qhull */
    facetT* facet;          /* set by FORALLfacets */
    int curlong, totlong;   /* memory remaining after qh_memfreeshort, used if !qh_NOmem  */
    int i;
    double total_t; /* CPU computing time to build the convex hull*/

    qhT qh_qh; /* Qhull's data structure.  First argument of most calls */
    qhT* qh = &qh_qh;

    rows = malloc(TOTpoints * sizeof(coordT*));

    QHULL_LIB_CHECK

    qh_zero(qh, errfile); /*initialize Qhull memory*/

    sprintf(flags, "qhull Fx");

    /* Translate the points to have the origin within the interior */
    trans_points(cur_points);

    /* Position points in a matrix , an input format required by qh_new_qhull */
    numpoints = cur_points->num_points;
    for (i = numpoints; i--;) {
        rows[i] = cur_points->qh_vec_p + dim * i;
    }
    exitcode =
        qh_new_qhull(qh, dim, numpoints, cur_points->qh_vec_p, ismalloc, flags, outfile, errfile);

    /*
     * At this point the solution is stored as follows (more
     * details in qhull-src/src/libqhull_r/libqhull_r.h for the
     * detailed description of the struct):
     *   qh->num_vertices: number of vertices in the hull
     *   qh->vertex_list: pointer to the head of a double linked list
     *   qh->vertex_tail: pointer to a dummy node after the last one
     *
     * For each vertex in the list pointed by p (for example
     * starting from qh->vertex_list), the coordinates are listed
     * starting from:
     *   p->point
     */
    if (!exitcode) { /* if no error */
        cur_points->num_sel = qh->num_vertices;
        cur_points->vec_sel = (int*)malloc(sizeof(int) * (cur_points->num_sel));
        edf_array_indexes(qh, cur_points->vec_sel);
    }

    /*Qhull records number of clocks to complete computation in qh->hulltime*/
    total_t = (double)(qh->hulltime) / CLOCKS_PER_SEC;

    qh_freeqhull(qh, !qh_ALL);               /* free long memory  */
    qh_memfreeshort(qh, &curlong, &totlong); /* free short memory and memory allocator */
    if (curlong || totlong)
        fprintf(errfile,
                "qhull internal warning: did not free %d bytes of long memory (%d "
                "pieces)\n",
                totlong, curlong);
    free(rows);
    return total_t;
}

void edf_free_points(edf_points_t* cur_points) {
    free(cur_points->t0);
    free(cur_points->t1);
    free(cur_points->vec_p);
    free(cur_points->vec_sel);
}
