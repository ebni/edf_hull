#include "edf_hull.h"

#include <glpk.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define IND (cur_points->vec_sel[i])

/**
 * [AD] Define this macro if you want to keep in the Convex Hull all the points
 * for which the objective value is equal to the RHS of relative constraint.
 */
/* #define KEEP_EQUALS */

/**
 * [AD] Define this macro as "glp_simplex" if you want traditional simplex, or
 * as "glp_exact" if you prefer GLPK exact rational solver.
 */
#define SOLVE_METHOD glp_simplex
/* #define SOLVE_METHOD glp_exact */

void edf_set_zero(edf_points_t *cur_points) {
  cur_points->alloc_points = 0;
  cur_points->t0 = NULL;
  cur_points->t1 = NULL;
  cur_points->vec_p = NULL;
  cur_points->num_sel = 0;
  cur_points->vec_sel = NULL;
}

void edf_create_points(const ts_t *cur_task_set, edf_points_t *cur_points) {
#define N (cur_task_set->num)
#define T(i) (cur_task_set->per[i])
#define D(i) (cur_task_set->dl[i])
#define O(i) (cur_task_set->phi[i])

  unsigned long *max_job;
  unsigned long i, j, k, s_i, s_j;
  int reserve_p;
  int num_deadline;
  int point_count;
  double cur_t, cur_s, coef;
  double big_enough;

  /* store the number of tasks */
  cur_points->num_tasks = N;

  /* Reserve one point for sum_i U_i <= 1, N points for the Ci >= 0 */
  reserve_p = N + 1;

  /* Limit until which computing the points */
  if (cur_task_set->has_phi) {
    big_enough =
        2 * cur_task_set->h_per + cur_task_set->max_d + cur_task_set->max_o;
  } else {
    big_enough = cur_task_set->h_per + cur_task_set->max_d;
  }

  /* Maximum job of the tasks, and number of tot abs deadlines */
  max_job = malloc(sizeof(*max_job) * N);
  num_deadline = 0;
  for (i = 0; i < N; i++) {
    max_job[i] = (unsigned long)floor((big_enough - O(i) - D(i)) / T(i));
    num_deadline += max_job[i] + 1;
  }

  /* Number of absolute deadlines  */
  if (cur_task_set->has_phi) {
    cur_points->num_points = reserve_p + num_deadline * num_deadline;
  } else {
    cur_points->num_points = reserve_p + num_deadline;
  }

  /* allocate the t0, t1 data structure for the points */
  if (cur_points->num_points > cur_points->alloc_points) {
    cur_points->t0 =
      (double *)realloc(cur_points->t0, sizeof(double) * cur_points->num_points);
    cur_points->t1 =
      (double *)realloc(cur_points->t1, sizeof(double) * cur_points->num_points);
    cur_points->vec_p =
      (double *)realloc(cur_points->vec_p, sizeof(double) * cur_points->num_points * N);
    cur_points->alloc_points = cur_points->num_points;
  }

  /* Start computing the points. */
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
  cur_points->t1[point_count] = cur_task_set->h_per;
  for (i = 0; i < N; i++) {
    cur_points->vec_p[point_count * N + i] = (int)(cur_task_set->h_per / T(i));
  }
  point_count++;

  /* if we have offset the points can be really many */
  if (cur_task_set->has_phi) {
    /* compute the points corresponding to the (start time, deadline) pairs */
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
              if (coef < 0)
                coef = 0;
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
      cur_points->t0 = (double *)realloc(
          cur_points->t0, sizeof(double) * cur_points->num_points);
      cur_points->t1 = (double *)realloc(
          cur_points->t1, sizeof(double) * cur_points->num_points);
      cur_points->vec_p = (double *)realloc(
          cur_points->vec_p, sizeof(double) * cur_points->num_points * N);
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
    cur_points->t1[point_count] = cur_task_set->h_per;
    for (i = 0; i < N; i++) {
      cur_points->vec_p[point_count * N + i] = (int)(cur_task_set->h_per / T(i));
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
          if (coef < 0)
            coef = 0;
          /* store it */
          cur_points->vec_p[point_count * N + k] = coef;
        }

        /* increment the point counter */
        point_count++;
      }
    }
  }

  /* free the max_job which is only a local data structure */
  free(max_job);
#undef N
#undef T
#undef D
#undef O
}

void edf_print_constraints_C(const edf_points_t *cur_points) {
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

unsigned int edf_linprog_points(edf_points_t *cur_points, unsigned int disp) {
   unsigned int i, j, index, rows, cols, util;
   char buf[256];
   glp_prob* lp;
   glp_smcp params;
   int* idxs;
   double* coefs, e;
   unsigned int* idx_sel;

   /* Absolute tolerance to check if constraint has to be added or not. */
   e = 1e-6;

   /* List of flags for selected points. */
   idx_sel = malloc(sizeof(*idx_sel) * cur_points->num_points);

   /* List of column indices, used by GLPK to index coefficients in rows. */
   idxs = malloc(sizeof(*idxs) * (cur_points->num_tasks + 1));
   for (j = 1; j <= cur_points->num_tasks; ++j)
      idxs[j] = j;

   /* List of row coefficients, used when adding a new row to the problem. */
   coefs = malloc(sizeof(*coefs) * (cur_points->num_tasks + 1));

   /* Setting up the problem, turning out GLPK screen output. */
   glp_term_out(GLP_OFF);
   lp = glp_create_prob();
   glp_set_prob_name(lp, "Convex Hull");

   /* Setting GLPK parameters. */
   glp_init_smcp(&params);
   /* params.meth = GLP_DUAL; */

   /* The problem is a MAXIMIZATION one. */
   glp_set_obj_dir(lp, GLP_MAX);
   glp_set_obj_name(lp, "Check");

   /**
    * FIRST PASS - Adding the constraints one by one.
    */
   if(disp) {
      printf("\n[%s]\n", ALGO_NAME);
      printf("--- FIRST PASS ---\n");
   }

   /* Adding variable C_j for each task j with LB set as 0. */
   index = 0;
   cols = 1;
   for (j = 0; j < cur_points->num_tasks; ++j) {
      glp_add_cols(lp, 1);

      sprintf(buf, "C_%d", j+1);
      glp_set_col_name(lp, cols, buf);

      /* The LB equal to 0 enforces the constraint [-C_j <= 0]. */
      glp_set_col_bnds(lp, cols, GLP_LO, 0.0, 0.0);

      idx_sel[index] = 1;
      cur_points->num_sel++;

      cols++;
      index++;
   }

   /* Adding the mandatory [\sum U_j <= 1] constraint. */
   rows = 1;
   glp_add_rows(lp, 1);

   sprintf(buf, "Util");
   glp_set_row_name(lp, rows, buf);

   glp_set_row_bnds(lp, rows, GLP_UP, 0,
      cur_points->t1[index]-cur_points->t0[index]);

   for (j = 0; j < cur_points->num_tasks; ++j) {
      coefs[j+1] = cur_points->vec_p[index * cur_points->num_tasks + j];
   }
   glp_set_mat_row(lp, rows, cur_points->num_tasks, idxs, coefs);

   idx_sel[index] = 1;
   cur_points->num_sel++;

   rows++;
   index++;

   /* Checking all constraints if binding by setting them as objective. */
   for (; index < cur_points->num_points; ++index) {
      for (j = 0; j < cur_points->num_tasks; ++j) {
         glp_set_obj_coef(lp, j+1,
            cur_points->vec_p[index * cur_points->num_tasks + j]);
      }

      SOLVE_METHOD(lp, &params);

      /** FIXME
       * Ora il controllo butta via i vincoli che han OBJ == RHS.
       * Casomai servisse conservarli, è sufficiente cambiare il valore della
       * costante KEEP_EQUALS nella relativa #define.
       */
#ifdef KEEP_EQUALS
      if (glp_get_obj_val(lp) < cur_points->t1[index]-cur_points->t0[index]) {
         /* Skipped if the objective value is lower than RHS of constraint. */
         idx_sel[index] = 0;

         if(disp) {
            printf("[Facet_%d] OBJ: %.2f --- ", index+1, glp_get_obj_val(lp));
            printf("RHS: %.2f\n", cur_points->t1[index]-cur_points->t0[index]);
         }

         continue;
      } else {
         /* Taken if the objective value breaks the RHS of constraint. */
         idx_sel[index] = 1;
         cur_points->num_sel++;
      }
#else
      if (glp_get_obj_val(lp) > cur_points->t1[index]-cur_points->t0[index] + e) {
         /* Taken if the objective value breaks the RHS of constraint. */
         idx_sel[index] = 1;
         cur_points->num_sel++;
      } else {
         /* Skipped if the objective value is lower than RHS of constraint. */
         idx_sel[index] = 0;

         if(disp) {
            printf("[Facet_%d] OBJ: %.2f --- ", index+1, glp_get_obj_val(lp));
            printf("RHS: %.2f\n", cur_points->t1[index]-cur_points->t0[index]);
         }

         continue;
      }
#endif

      glp_add_rows(lp, 1);

      sprintf(buf, "Facet_%d", index+1);
      glp_set_row_name(lp, rows, buf);

      glp_set_row_bnds(lp, rows, GLP_UP, 0,
         cur_points->t1[index]-cur_points->t0[index]);

      for (j = 0; j < cur_points->num_tasks; ++j) {
         coefs[j+1] = cur_points->vec_p[index * cur_points->num_tasks + j];
      }
      glp_set_mat_row(lp, rows, cur_points->num_tasks, idxs, coefs);

      rows++;
   }

   /* Clearing objective function coefficients. */
   for (j = 0; j < cur_points->num_tasks; ++j) {
      glp_set_obj_coef(lp, j+1, 0);
   }

   /* Saving the problem as a CPLEX LP human-readable file. */
   if(disp) {
      glp_write_lp(lp, NULL, "before.lp");
   }

   /**
    * SECOND PASS - Checking and eventually removing all inserted constraints.
    */
   if(disp) {
      printf("--- SECOND PASS ---\n");
   }

   /* The non-negativity constraints are mandatory. */
   index = cur_points->num_tasks;

   /* Checking, one by one, the inserted constraints if they are binding. */
   i = 1;
   util = 1;
   for(; index < cur_points->num_points; ++index) {
      if (!idx_sel[index]) {
         continue;
      }

      /* Set the constraint coefficients as objective value coefficients. */
      for (j = 0; j < cur_points->num_tasks; ++j) {
         glp_set_obj_coef(lp, j+1,
            cur_points->vec_p[index * cur_points->num_tasks + j]);
      }

      /* Remove the RHS of tested constraint. */
      glp_set_row_bnds(lp, i, GLP_FR, 0, 0);

      SOLVE_METHOD(lp, &params);

      /** FIXME
       * Ora il controllo butta via i vincoli che han OBJ == RHS.
       * Casomai servisse conservarli, è sufficiente cambiare il valore della
       * costante KEEP_EQUALS nella relativa #define.
       */
#ifdef KEEP_EQUALS
      if (glp_get_obj_val(lp) < cur_points->t1[index]-cur_points->t0[index]) {
         /* Removed if the objective value is lower than RHS of constraint. */
         idx_sel[index] = 0;
         cur_points->num_sel--;

         if (i == 1) {
            util = 0;
         }

         if(disp) {
            printf("[%s] OBJ: %.2f --- ", glp_get_row_name(lp, i),
               glp_get_obj_val(lp));
            printf("RHS: %.2f\n", cur_points->t1[index]-cur_points->t0[index]);
         }
      } else {
         /* Add again the correct RHS of tested constraint. */
         glp_set_row_bnds(lp, i, GLP_UP, 0,
         cur_points->t1[index]-cur_points->t0[index]);
      }
#else
      if (glp_get_obj_val(lp) > cur_points->t1[index]-cur_points->t0[index] + e
            || glp_get_status(lp) == GLP_UNBND) {
         /* Add again the correct RHS of tested constraint. */
         glp_set_row_bnds(lp, i, GLP_UP, 0,
         cur_points->t1[index]-cur_points->t0[index]);
      } else {
         /* Removed if the objective value is lower than RHS of constraint. */
         idx_sel[index] = 0;
         cur_points->num_sel--;

         if (i == 1) {
            util = 0;
         }

         if(disp) {
            printf("[%s] OBJ: %.2f --- ", glp_get_row_name(lp, i),
               glp_get_obj_val(lp));
            printf("RHS: %.2f\n", cur_points->t1[index]-cur_points->t0[index]);
         }
      }
#endif

      i++;
   }

   /* Allocating the list of selected points. */
   cur_points->vec_sel = malloc(
      sizeof(*cur_points->vec_sel) * cur_points->num_sel);

   /* Filling the list of selected points, with correct indexes. */
   index = 0;
   for (i = 0; i < cur_points->num_points; ++i) {
      if (idx_sel[i]) {
         cur_points->vec_sel[index] = i;
         index++;
      }
   }

   /* Clearing objective function coefficients. */
   for (j = 0; j < cur_points->num_tasks; ++j) {
      glp_set_obj_coef(lp, j+1, 0);
   }

   /* Saving the problem as a CPLEX LP human-readable file. */
   if(disp) {
      glp_write_lp(lp, NULL, "after.lp");
   }

   /* Cleaning up stuff. */
   glp_delete_prob(lp);

   free(coefs);
   free(idxs);
   free(idx_sel);

   if(disp) {
      printf("--- END ---\n");
      printf("[%s]\n", ALGO_NAME);
   }

   return util;
}

void edf_free_points(edf_points_t *cur_points) {
  free(cur_points->t0);
  free(cur_points->t1);
  free(cur_points->vec_p);
  free(cur_points->vec_sel);
}

void edf_print_constraints_U(const ts_t *cur_task_set,
                             const edf_points_t *cur_points) {
  int i, j;
  printf("\nOr, alternatively, minimal set of constraints can also be written "
         "as\n\n");
  printf("    a_1*U_1 +a_2*U_2 + ... +a_N*U_N <= t_1-t_0\n\n");
  for (j = 0; j < cur_points->num_tasks; j++) {
    printf("a_%d\t", j + 1);
  }
  printf("t1\tt0\n");

#define T(j) (cur_task_set->per[j])
  for (i = 0; i < cur_points->num_sel; i++) {
    for (j = 0; j < cur_points->num_tasks; j++) {
      if (i < cur_points->num_tasks) {
	/* Special case of -Ci <= 0 constraint */
        printf("%.4f\t", cur_points->vec_p[IND * cur_points->num_tasks + j]);
      } else {
	/* Other constraints: the 1st one is tot util<=1 */
        printf("%.4f\t",
               (cur_points->vec_p[IND * cur_points->num_tasks + j] * T(j)) /
                   cur_points->t1[IND]);
      }
    }
    if (i < cur_points->num_tasks) {
      printf("0\t%.0f\n", cur_points->t0[IND]);
    } else {
      printf("1\t%.0f\n", cur_points->t0[IND]);
    }

#undef T
  }
}

void edf_print_points(const edf_points_t *cur_points) {
  int i, j;

  printf("Number of tasks:\t%d\n", cur_points->num_tasks);
  printf("Number of points:\t%d\n", cur_points->num_points);
  printf("Points:\n");
  for (i = 0; i < cur_points->num_points; i++) {
    printf("t0=%.3f\tt1=%.3f", cur_points->t0[i], cur_points->t1[i]);
    for (j = 0; j < cur_points->num_tasks; j++) {
      printf("\t%.3f", cur_points->vec_p[i * cur_points->num_tasks + j]);
    }
    printf("\n");
  }
}
