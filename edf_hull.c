#include      <stdlib.h>
#include      <stdio.h>
#include      <math.h>

#define    MAX(a,b)      (a>b ? a : b)
#define    MIN(a,b)      (a<b ? a : b)

/*
 * Task set  related data. It  is preferred to store  homogeneous data
 * (periods,  deadline,...) in  unique vector,  rather than  storing a
 * vector  of  "task  structures",  where each  element  contains  the
 * parameters. This is for efficiency.
 */
typedef struct {
	int        number;          /* number of tasks */
	double*    period;          /* periods */
	double*    rel_deadline;    /* rel. deadline (even > period) */
	double*    offset;          /* offsets (in [0,period)) */
	int        phasing;         /* 1 if offset; 0 if not */
} task_set_t;

/*
 * Data structure  storing the  interval info.  For each  interval the
 * test of the demand is performed. The test performed is:
 *
 *                sum_i vec_p[i]*C[i] <= t1[i]-t0[i]
 *
 * where C[i] denotes the computation time of the i-th task.
 */
typedef struct {
	int        num_points;      /* number of points */
	int        num_tasks;       /* number of tasks */
	double*    t0;              /* starting point */
	double*    t1;              /* finishing point */
	double*    vec_p;           /* matrix of coefficients */
	double     interior;        /* the point interior*[1,1,...,1] is in
				       the interior */
	int        num_sel;         /* number of selected points (by qhull) */
	int*       vec_sel;         /* indexes of selected points */
} points_t;

/* The global data structure of the task set */
task_set_t  my_task_set;

/* The global data structure of the scheduling points */
points_t    my_points;

/*
 * It extends the greatest common divisor algorithm to doubles. It stops
 * when the difference is smaller than some small-enough value. To
 * simplify the algorithm a MUST be greater than of equal to b.
 */
#define   EPS      1.0e-6  /* two double closer than EPS are
			    * considered to be equal*/
double gcd_double(double a, double b)
{
	double diff;
	
	diff = a-b;
	if (diff < EPS)
		return a;
	else
		return gcd_double(MAX(diff,b), MIN(diff,b));
}

/*
 * It extends the least common multiple algorithm to doubles. It stops
 * when the difference is smaller than some small-enough value. To
 * simplify the algorithm a MUST be greater than of equal to b.
 */
double lcm_double(double a, double b)
{
	return a*b/gcd_double(MAX(a,b),MIN(a,b));
}

/*
 * It returns a good-enough time limit of a task set.
 */
double time_limit(const task_set_t* cur_task_set)
{
#define     N          (cur_task_set->number)
#define     T(i)       (cur_task_set->period[i])
#define     D(i)       (cur_task_set->rel_deadline[i])
#define     O(i)       (cur_task_set->offset[i])

	double hyperp, max_d, max_o;
	int i;
	
	/* compute the hyperperiod, max deadline, max offset */
	hyperp = T(0);
	max_d = D(0);
	max_o = O(0);
	for(i=1; i<N; i++) {
		hyperp = lcm_double(MAX(hyperp, T(i)), MIN(hyperp, T(i)));
		max_d = MAX(max_d, D(i));
		max_o = MAX(max_o, O(i));
	}
	
	/* the time limit value */
	return 2*hyperp+max_d+max_o;
}


/*
 * Initialize the task set reading from stdin
 */
void init_task_set(task_set_t* cur_task_set)
{
	
	int i, num;
	double cur_T, cur_D, cur_O;
	int cur_phasing;
	
	/* initially set no phasing */
	cur_phasing = 0;
	
	/* read the number of tasks */
	printf("Please insert the number of tasks\n");
	scanf("%d", &num);
	cur_task_set->number = num;
	
	/* allocate the data structure */
	cur_task_set->period =        (double*)malloc(sizeof(double)*num);
	cur_task_set->rel_deadline =  (double*)malloc(sizeof(double)*num);
	cur_task_set->offset =        (double*)malloc(sizeof(double)*num);
	
	/* read the task set parameters */
	printf("For each task, please insert\n");
	printf("<period> <rel. deadline> <offset>\n");
	for(i=0; i<cur_task_set->number; i++) {
		/* read... */
		printf("Task %2i\n", i);
		scanf("%lf %lf %lf", &cur_T, &cur_D, &cur_O);
    
		/* store the period */
		cur_task_set->period[i] = cur_T;
		
		/* store the deadline */
		cur_task_set->rel_deadline[i] = cur_D;
		
		/* compute the remainder of the phase and store it */
		cur_O = cur_O - floor(cur_O/cur_T)*cur_T;
		cur_task_set->offset[i] = cur_O;
		
		/* update the phasing */
		cur_phasing |= cur_O != 0;
	}
	cur_task_set->phasing = cur_phasing;
	
	/* if no task phasing the offset structure can be deallocated */
	/*   if (!cur_task_set->phasing) */
	/*     free(cur_task_set->offset); */
}

/*
 * Print the data of the task set
 */
void print_task_set(const task_set_t* cur_task_set)
{
	int i;
	
	printf("The task set is %s phasing\n", 
	       cur_task_set->phasing ? "WITH" : "WITHOUT");
	printf(" i: (      Ti,       Di,       Oi)\n");
	for(i=0; i<cur_task_set->number; i++) {
		printf("%2i: (%8.1f, %8f, %8f)\n", i,
		       cur_task_set->period[i], 
		       cur_task_set->rel_deadline[i], 
		       cur_task_set->offset[i]);
	}
}

/*
 * Starting from the task set given in cur_task_set, the time limit
 * given in big_enough, it returns all the points, for all the pairs
 * t0 < t1.
 */
void create_points(const task_set_t* cur_task_set,
		   const double big_enough,
		   points_t* cur_points)
{
#define     N          (cur_task_set->number)
#define     T(i)       (cur_task_set->period[i])
#define     D(i)       (cur_task_set->rel_deadline[i])
#define     O(i)       (cur_task_set->offset[i])
	
	int* max_job;
	int i, j, k, reserve_p, s_i, s_j;
	int num_deadline;
	int point_count;
	double cur_t, cur_s, coef, sum_rate;
	
#ifdef      DEBUG
	printf("Big enough is %f\n", big_enough);
#endif
	
	/* store the number of tasks */
	cur_points->num_tasks = N;
	
	/* 
	 * Maximum job of the tasks. Also compute the sum of max_job[i],
	 * because this is also the number of considered deadlines.
	 */
	max_job = (int*)malloc(sizeof(int)*N);
	num_deadline = 0;
	for(i=0; i<N; i++) {
		max_job[i] = (int)floor((big_enough-O(i)-D(i))/T(i));
#ifdef      DEBUG
		printf("max_job[%d] = %d\n", i, max_job[i]);
#endif
		num_deadline += max_job[i]+1;
	}
	
	/*
	 * Start computing the points. Reserve the initial point for a
	 * special purpose. The number of reserved points is reserv_p.
	 */
	/* Reserve one point for sum_i U_i <= 1, N points for the Ci >= 0 */ 
	reserve_p = N+1;
	
	/* if we have offset the points can be really many */
	if (cur_task_set->phasing) {
		/* the number of points (overestimating it) */
		cur_points->num_points = reserve_p + num_deadline*num_deadline;
		
		/* allocate the t0, t1 data structure for the points */
		cur_points->t0 = (double*)malloc(sizeof(double)*cur_points->num_points);
		cur_points->t1 = (double*)malloc(sizeof(double)*cur_points->num_points);
		cur_points->vec_p = (double*)malloc(sizeof(double)*cur_points->num_points*N);
		
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
		 * compute the points corresponding to the (start time, deadline)
		 * pairs
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
#ifdef DEBUG
						printf("i=%d, j=%d, s_i=%d, s_j=%d\n", i, j, s_i, s_j);
#endif
						cur_points->t1[point_count] = cur_t;    
						cur_points->t0[point_count] = cur_s;
						for(k=0; k<N; k++) {
							/* there may be numerical problems... */
							coef = floor((cur_t-O(k)-D(k))/T(k))-ceil((cur_s-O(k))/T(k))+1;
							/* must be non-negative */
							coef = MAX(0, coef);
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
		
		/* allocate the data structure for the points */
		cur_points->t0 = (double*)malloc(sizeof(double)*cur_points->num_points);
		cur_points->t1 = (double*)malloc(sizeof(double)*cur_points->num_points);
		cur_points->vec_p = (double*)malloc(sizeof(double)*cur_points->num_points*N);
		
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
					coef = MAX(0, coef);
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
		sum_rate += 1/MIN(T(i),D(i));
	}
	cur_points->interior = 1/(2*sum_rate);
	
	/* free the max_job which is only a local data structure */
	free(max_job);
#undef   N
#undef   T
#undef   D
#undef   O
}

void print_points(const points_t* cur_points)
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
}

/*
 * It selects the only necessary points using a geometric
 * reasoning. The algorithm for computing the convex hull is invoked
 * by the external executable "qhull" (see below for the syntax of
 * invocation. The data is passed to qhull by the exchange of text
 * files.
 */
void select_qhull_points(points_t * cur_points)
{
	FILE* points_file;
	FILE* selection_file;
	int i,j;
	double sum_coef, mul_factor;
	
	/* open/create the file in write mode */
	points_file = fopen("data_in.txt","w");
	
	/* number of tasks */
	fprintf(points_file,"%d\n", cur_points->num_tasks);
	
	/* number of points */
	fprintf(points_file,"%d\n", cur_points->num_points);
	
	/* all the points */
	for(i=0; i<cur_points->num_points; i++) {
		/* compute the sum of all the coefficients of the current vector */
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
			fprintf(points_file,"%.16f ",
				cur_points->vec_p[i*cur_points->num_tasks+j]*mul_factor);
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
	
	/* close the file */
	fclose(selection_file);
}

/*
 * Free all the dynamically allocated memory
 */
void free_all()
{
	free(my_task_set.period);
	free(my_task_set.rel_deadline);
	free(my_task_set.offset);
	free(my_points.t0);
	free(my_points.t1);
	free(my_points.vec_p);
	free(my_points.vec_sel);
}


int main(int argc, char* argv[])
{
	
	init_task_set(&my_task_set);
#ifdef DEBUG
	print_task_set(&my_task_set);
#endif
	create_points(&my_task_set, time_limit(&my_task_set), &my_points);
#ifdef DEBUG
	print_points(&my_points);
#endif
	select_qhull_points(&my_points);
	
	free_all();
}
