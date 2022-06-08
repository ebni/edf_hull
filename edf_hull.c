#include      <stdlib.h>
#include      <stdio.h>
#include      <math.h>
#include      "edf_hull.h"

void edf_set_zero(edf_points_t* cur_points)
{
	cur_points->alloc_points = 0;
	cur_points->t0 = NULL;
	cur_points->t1 = NULL;
	cur_points->vec_p = NULL;
	cur_points->qh_vec_p = NULL;
	cur_points->num_sel = 0;
	cur_points->vec_sel = NULL;
}


void edf_create_points(const ts_t* cur_task_set,
		       edf_points_t* cur_points)
{
#define     N          (cur_task_set->num)
#define     T(i)       (cur_task_set->per[i])
#define     D(i)       (cur_task_set->dl[i])
#define     O(i)       (cur_task_set->phi[i])
	
	unsigned long * max_job;
	int i, j, k, reserve_p, s_i, s_j;
	int num_deadline;
	int point_count;
	double cur_t, cur_s, coef, sum_rate;
	double big_enough;
	
	/* store the number of tasks */
	cur_points->num_tasks = N;

	/* Set the limit until which computing the points */
	if (cur_task_set->has_phi) {
		big_enough = 2*cur_task_set->h_per+cur_task_set->max_d
			+cur_task_set->max_o;
	} else {
		big_enough = cur_task_set->h_per+cur_task_set->max_d;
	}
	
	/* 
	 * Maximum job of the tasks. Also compute the sum of max_job[i],
	 * because this is also the number of considered deadlines.
	 */
	max_job = malloc(sizeof(*max_job)*N);
	num_deadline = 0;
	for(i=0; i<N; i++) {
		max_job[i] = (unsigned long)floor((big_enough-O(i)-D(i))/T(i));
		num_deadline += max_job[i]+1;
	}
	
	/*
	 * Start computing the points. Reserve the initial point for a
	 * special purpose. The number of reserved points is reserv_p.
	 */
	/* Reserve one point for sum_i U_i <= 1, N points for the Ci >= 0 */ 
	reserve_p = N+1;
	
	/* if we have offset the points can be really many */
	if (cur_task_set->has_phi) {
		/* the number of points (overestimating it) */
		cur_points->num_points = reserve_p + num_deadline*num_deadline;
		
		if (cur_points->num_points > cur_points->alloc_points) {
			/* allocate the t0, t1 data structure for the points */
			cur_points->t0 = (double*)realloc(cur_points->t0,sizeof(double)*cur_points->num_points);
			cur_points->t1 = (double*)realloc(cur_points->t1,sizeof(double)*cur_points->num_points);
			cur_points->vec_p = (double*)realloc(cur_points->vec_p,sizeof(double)*cur_points->num_points*N);
			/* The points to be used are in cur_points->qh_vec_p */
			cur_points->qh_vec_p = (coordT*)realloc(cur_points->qh_vec_p,sizeof(coordT)*cur_points->num_points*N);
			cur_points->alloc_points = cur_points->num_points;
		}

		/* initialize the point counter */
		point_count = 0;
		
		/* store the -Ci <= 0 constraint */
		for(i=0; i<N; i++) {
			cur_points->t0[point_count] = 0;
			cur_points->t1[point_count] = 0;
			for(j=0; j<N; j++) {
				cur_points->vec_p[point_count*N+j] = ((i==j) ? -1 : 0);
			}
			
			/* increment the point counter */
			point_count++;
		}
		
		/* store the sum Ui <= 1 constraint */
		cur_points->t0[point_count] = 0;
		cur_points->t1[point_count] = 1;
		for(i=0; i<N; i++) {
			cur_points->vec_p[point_count*N+i] = 1/T(i);
		}
		point_count++;
		
		/* 
		 * compute the points corresponding to the (start
		 * time, deadline) pairs
		 */
		/* loop on the absolute deadlines */
		for(i=0; i<N; i++) {
			for(j=0; j<=max_job[i]; j++) {
				/* t1 is equal to the j-th deadline of task i */
				cur_t = O(i)+j*T(i)+D(i);
				
				/* loop on the start times */
				for(s_i=0; s_i<N; s_i++) {
					for(s_j=0; s_j<=max_job[s_i]; s_j++) {
						/* the current start time */
						cur_s = O(s_i)+s_j*T(s_i);
						if (cur_s >= cur_t) /* start time is ahead the deadline */
							break;
						cur_points->t1[point_count] = cur_t;    
						cur_points->t0[point_count] = cur_s;
						for(k=0; k<N; k++) {
							/* there may be numerical problems... */
							coef = floor((cur_t-O(k)-D(k))/T(k))-ceil((cur_s-O(k))/T(k))+1;
							/* must be non-negative */
							if (coef < 0)
								coef = 0;
							/* store it */
							cur_points->vec_p[point_count*N+k] = coef;
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
			cur_points->t0 = (double*)realloc(cur_points->t0,sizeof(double)*cur_points->num_points);
			cur_points->t1 = (double*)realloc(cur_points->t1,sizeof(double)*cur_points->num_points);
			cur_points->vec_p = (double*)realloc(cur_points->vec_p,sizeof(double)*cur_points->num_points*N);
			/* The points to be used are in cur_points->qh_vec_p */
			cur_points->qh_vec_p = (coordT*)realloc(cur_points->qh_vec_p,sizeof(coordT)*cur_points->num_points*N);
			cur_points->alloc_points = cur_points->num_points;
		}
		
		/* initialize the point counter */
		point_count = 0;
		
		/* store the -Ci <= 0 constraint */
		for(i=0; i<N; i++) {
			cur_points->t0[point_count] = 0;
			cur_points->t1[point_count] = 0;
			for(j=0; j<N; j++) {
				cur_points->vec_p[point_count*N+j] = ((i==j) ? -1 : 0);
			}
			
			/* increment the point counter */
			point_count++;
		}
		
		/* store the sum Ui <= 1 constraint */
		cur_points->t0[point_count] = 0;
		cur_points->t1[point_count] = 1;
		for(i=0; i<N; i++) {
			cur_points->vec_p[point_count*N+i] = 1/T(i);
		}
		point_count++;
		
		/* compute the points corresponding to deadlines */
		for(i=0; i<N; i++) {
			for(j=0; j<=max_job[i]; j++) {
				/* t0 is zero, when no offset */
				cur_points->t0[point_count] = 0;
				/* t1 is equal to the j-th deadline of task i */
				cur_points->t1[point_count] = cur_t = O(i)+j*T(i)+D(i);
				for(k=0; k<N; k++) {
					/* there may be numerical problems... */
					coef = floor((cur_t-O(k)-D(k))/T(k))+1;
					/* must be non-negative */
					if (coef < 0)
						coef = 0;
					/* store it */
					cur_points->vec_p[point_count*N+k] = coef;
				}
				
				/* increment the point counter */
				point_count++;
			}
		}
	}
	/* compute the interior for traslation purpose */
	sum_rate = 0;
	for(i=0; i<N; i++) {
		sum_rate += 1/(T(i) < D(i) ? T(i) : D(i));
	}
	cur_points->interior = 1/(2*sum_rate);
	
	/* free the max_job which is only a local data structure */
	free(max_job);
#undef   N
#undef   T
#undef   D
#undef   O
}

/*
 * Translate the points so that the origin belongs to the interior.
 * The translated vectors are written to the array
 * cur_points->qh_vec_p
 */
static void trans_points(edf_points_t* cur_points)
{
	int i,j;
	double sum_coef, sum_rate, mul_factor;
	
	/* all the points */
	for(i=0; i<cur_points->num_points; i++) {
		/* compute the sum of all the coefficients of the
		 * current vector */
		sum_coef = 0;
		for(j=0; j<cur_points->num_tasks; j++) {
			sum_coef += cur_points->vec_p[i*cur_points->num_tasks+j];
		}
		
		/*
		 * the  vector  coefficients  are  multiplied  by  the
		 * following factor.  Due to  the proper  selection on
		 * the interior  point, this  factor should  ALWAYS be
		 * GREATER than zero .
		 */
		mul_factor = 1/(cur_points->t1[i]-cur_points->t0[i]-cur_points->interior*sum_coef);
		for(j=0; j<cur_points->num_tasks; j++) {
			cur_points->qh_vec_p[i*cur_points->num_tasks+j] = 
				cur_points->vec_p[i*cur_points->num_tasks+j]*mul_factor;
		}
	}
}

void edf_print_points(const edf_points_t* cur_points)
{
	int i,j;
	
	printf("Number of tasks:\n");
	printf("%d\n", cur_points->num_tasks);
	printf("Number of points:\n");
	printf("%d\n", cur_points->num_points);
	printf("Interior:\n");
	printf("%f\n", cur_points->interior);
	printf("Points:\n");
	for(i=0; i<cur_points->num_points; i++) {
		printf("t0=%f\tt1=%f",cur_points->t0[i], cur_points->t1[i]);
		for(j=0; j<cur_points->num_tasks; j++) {
			printf("\t%f",cur_points->vec_p[i*cur_points->num_tasks+j]);
		}
		printf("\n");
	}
	/* 
	 * Print the selection in the following format
	 * 
	 * t0    t1    a_1    a_2    ...    a_N
	 *
	 * For each row the test that is necessary (and sufficient) to run
	 * is
	 * 
	 * a_1*C_1 +a_2*C_2 + ... +a_N*C_N <= t_1-t_0
	 */
	printf("Minimal set of constraints are written in the following form\n\n");
	printf("    eta_1*C_1 +eta_2*C_2 + ... +eta_N*C_N <= t_1-t_0\n\n");
	printf("t0\tt1");
	for(j=0; j<cur_points->num_tasks; j++) {
		printf("\teta_%d", j+1);
	}
	printf("\n");
	for(i=0; i<cur_points->num_sel; i++) {
#define IND (cur_points->vec_sel[i])
		printf("%.0f\t%.0f",cur_points->t0[IND], cur_points->t1[IND]);
		for(j=0; j<cur_points->num_tasks; j++) {
			printf("\t%.4f",cur_points->vec_p[IND*cur_points->num_tasks+j]);
		}
		printf("\n");
#undef IND
	}
	

}

#ifdef USE_QHULL_LIB
/*
 * Necessary initialization of the qhull data structure. DONT TOUCH
 * unless you know what you are doing. (EB: I knew very little when I
 * did this and it seemed to work by miracle)
 */
static void init_qhull_struct(qhT * qh)
{
	int seed, exitcode;
  
	qh->fin = stdin;
	qh->fout = stdout;
	qh->ferr = stderr;
	qh->qhmem.ferr = stderr;
	if ((qh->MINdenom_1 = 1.0/REALmax) < REALmin) {
		qh->MINdenom_1 = REALmin;
	}
	seed= (int)time(NULL);
	qh_RANDOMseed_(qh, seed);
	qh->run_id= qh_RANDOMint;
	exitcode = setjmp(qh->errexit); /* simple statement for CRAY J916 */
	qh_initstatistics(qh);
}
#endif /* USE_QHULL_LIB */

/*
 * If the USE_QHULL_LIB is defined, then the API of qhull is
 * used. This choice is faster but still experimental. If the
 * USE_QHULL_LIB is *not* defined, then qhull is invoked via command
 * line interface and by properly writing/reading text files.
 */
void edf_qhull_points(edf_points_t * cur_points)
{
	int i,j;
#ifdef USE_QHULL_LIB
	/* We aim at simulating to "qhull Fx TI data_in.txt TO data_out.txt" */
	qhT * qh;
#else
	/* Using text file as interface with the qhull executable */
	FILE * points_file;
	FILE * selection_file;
#endif
	
	/* Translate the points to have the origin within the interior */
	trans_points(cur_points);
#ifdef USE_QHULL_LIB
	/* 
	 * All the qhull data is stored in the qh_qh data structure,
	 * global variable declared in qh_init_struct.c. The
	 * initialization of the stuct qh_qh is made:
	 *   (1) at variable initialization time in qh_init_struct.c
	 *   (2) by invoking the function init_qhull_struct(qh)
	 *   (3) in the next few lines
	 */
	/* setting the pointer to the qhull data structure (global variable) */
	init_qhull_struct(&qh_qh);

	qh = &qh_qh;
	qh->normal_size = cur_points->num_tasks * sizeof(coordT);
	/* maxline is length of some string. 500 is MAGIC number
	 * copied from qhull code */
	qh->maxline = 500 > (cur_points->num_tasks*(qh_REALdigits + 5)) ? 500 : (cur_points->num_tasks*(qh_REALdigits + 5));
	qh_init_B(qh,
		  cur_points->qh_vec_p,      /* array of coordinates */
		  cur_points->num_points,   /* number of vectors */
		  cur_points->num_tasks,    /* dimension of vectors */
		  1);                       /* true is vec_p is allocated */

	/* Finally invoking the power of qhull */
	qh_qhull(qh);
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
	cur_points->num_sel = qh->num_vertices;
	/*
	 * Una cosa buona da fare sembra invocare 
	 *
	 *   qh_printextremes(qh, fp, facetlist, facets, printall);
	 *
	 * definita in io_r.c: 1808--1833, in particolare sfruttare il loop
	 *
	 *   FOREACHvertex_(vertices)
	 *     id= qh_pointid(qh, vertex->point);
	 *  
	 * che ho verificato che ritorna gli ID dei punti passati
	 */
	/* FIXME, TODO: to get from the qhull data struct the indices
	 * of the vertices. Below NULL should be replaced by a proper
	 * array with the right indices of vertices */
	cur_points->vec_sel = NULL;
#else  /* not USE_QHULL_LIB */
	/* open/create the file in write mode */
	points_file = fopen("data_in.txt","w");
	
	/* number of tasks */
	fprintf(points_file,"%d\n", cur_points->num_tasks);
	
	/* number of points */
	fprintf(points_file,"%d\n", cur_points->num_points);
	
	/* all the points */
	for(i=0; i<cur_points->num_points; i++) {
		for(j=0; j<cur_points->num_tasks; j++) {
			fprintf(points_file,"%.24f ", cur_points->qh_vec_p[i*cur_points->num_tasks+j]);
		}
		fprintf(points_file,"\n");
	}
	
	/* close the file */
	fclose(points_file);
	
	/*
	 * At this point we run the following command
	 *
	 *   qhull Fx TI data_in.txt TO data_out.txt
	 *
	 * After the execution is finished, the file data_out.txt contains
	 * the indexes of the constraints, which describe the necessary and
	 * sufficient condition.
	 *
	 * This list of constraints should ALWAYS contain the constraint
	 * indexes from 0 to N-1, because they correspond to the condition
	 * Ci >=0.
	 */
	/*
	 * An old version of qhull was suggesting to add the options
	 * Q12 Q15 upon an error message. At the time of committing
	 * this change these options don't work anymore, so removing
	 * them 
	 */
	system("./qhull Fx TI data_in.txt TO data_out.txt");
	
	/* open the file produced by qhull in read mode */
	selection_file = fopen("data_out.txt","r");
	
	/* number of points */
	fscanf(selection_file,"%d\n", &(cur_points->num_sel));
	cur_points->vec_sel = (int*)malloc(sizeof(int)*(cur_points->num_sel));
	
	/* load the index of the points */
	for(i=0; i<cur_points->num_sel; i++) {
	fscanf(selection_file,"%d\n", &(cur_points->vec_sel[i]));
}
	
	/* 
	 * Print the selection in the following format
	 * 
	 * t0    t1    a_1    a_2    ...    a_N
	 *
	 * For each row the test that is necessary (and sufficient) to run
	 * is
	 * 
	 * a_1*C_1 +a_2*C_2 + ... +a_N*C_N <= t_1-t_0
	 */
	printf("Minimal set of constraints are written in the following form\n\n");
	printf("    eta_1*C_1 +eta_2*C_2 + ... +eta_N*C_N <= t_1-t_0\n\n");
	for(j=0; j<cur_points->num_tasks; j++) {
		printf("eta_%d\t", j+1);
	}
	printf("t1\tt0\n");
	for(i=0; i<cur_points->num_sel; i++) {
#define IND (cur_points->vec_sel[i])
		for(j=0; j<cur_points->num_tasks; j++) {
			printf("%.4f\t",cur_points->vec_p[IND*cur_points->num_tasks+j]);
		}
		printf("%.0f\t%.0f\n", cur_points->t1[IND], cur_points->t0[IND]);
#undef IND
	}
	
	/* close the file */
	fclose(selection_file);
#endif /* USE_QHULL_LIB */
}

void edf_free_points(edf_points_t * cur_points)
{
	free(cur_points->t0);
	free(cur_points->t1);
	free(cur_points->vec_p);
	free(cur_points->vec_sel);
	free(cur_points->qh_vec_p);
}
