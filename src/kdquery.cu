/* 
 * Copyright (C) 2016-2017 Andrew Pellegrino
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//#include <pthread.h>
extern "C" {
	#include <string.h>
	#include <stddef.h>
	#include <stdio.h>
	#include <stdlib.h>
	#include <cuda.h>
	#include "kdtree.h"
}
#define HANDLE_ERROR( err ) ( HandleError( err, __FILE__, __LINE__ ) )


static void HandleError( cudaError_t err, const char *file, int line )
{
	    if (err != cudaSuccess)
			    {
					    printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
								            file, line );
						    exit( EXIT_FAILURE );
							    }
}

void *__gxx_personality_v0;
static int need_alloc = 1;
__device__ double * _data_query;
__device__ int * _field_query;
__device__ node_t * _tree_data;
__device__ double rad;
__device__ int tree_node_count;

/*
static double * _data_query;
static int *_field_query;
static node_t * _tree_data;
static double rad;
*/

__device__ void count_node(double n, double rsq, field_counter_t * counter,
                       int this_field, int query_field, int this_weight)
{
    if (n < rsq) {
        counter->array[query_field * counter->num_fields + this_field]+=this_weight;
    }
    return;
}
/*
static inline void count_node(double n, double rsq, field_counter_t * counter,
                       int this_field, int query_field, int this_weight)
{
    if (n < rsq) {
        counter->array[query_field * counter->num_fields + this_field] += this_weight;
    }
    return;
}
*/
__device__ double norm2(datum_t * a, datum_t * b)
{
    double sum = 0.0;
    int i;
    for(i=0; i<NDIM; i++) {
        sum += (a->value[i] - b->value[i])*(a->value[i] - b->value[i]);
    }
    return sum;
}
/*
static inline double norm2(datum_t * a, datum_t * b)
{
    double sum = 0.0;
    int i;
    for(i=0; i<NDIM; i++) {
        sum += (a->value[i] - b->value[i])*(a->value[i] - b->value[i]);
    }
    return sum;
}
*/

__device__ int left_child(int p)
{
    return 2*p;
}

__device__ int right_child(int p)
{
    return 2*p+1;
}
/*
static inline int left_child(int p)
{
    return 2*p;
}

static inline int right_child(int p)
{
    return 2*p+1;
}
*/
/*
 * Query how many points in the tree with head p lie within radius r of point q
 */

__device__ void radius(int pi, int qi, field_counter_t * counter, int tid)
{
    node_t p = *(_tree_data+pi);
    int this_field = p.data.field;
    int query_field = 0;
    if(_field_query != NULL) {
        query_field = _field_query[qi];
    }

    double val;

    double rsq;
    double pos_upper = 0.0, pos_lower = 0.0, point = 0.0;

    rsq = rad*rad;

    datum_t datum;
    int i;
    for(i=0;i<NDIM;i++) {
        datum.value[i] = _data_query[NDIM*qi+i];
    }
    double n = norm2(&p.data, &datum);


    //if (tid == 1) {printf("Tree Size: %i\n", tree_node_count); }
    int buff_size = int(log2f(tree_node_count)); // size of our node_id_buffer
    buff_size = int(0.4472135954999579*powf(tree_node_count,0.694242));//
    int* node_id_buffer = (int*)malloc(buff_size*sizeof(int)); // holds values of pi for later
    /* holy shit why have I not been doing this? */
    if (node_id_buffer == NULL) { printf("Failed to Allocate Memory\n"); 
	    return;}/// radius(pi, qi, counter, tid); }
    node_id_buffer[0] = pi; // init the buffer
    int stack_size = 1; //number of elements in the stack, also index of next unalloacted element
    i = 0; // stack pointer
    int stack_max = 0; // for debugging
    int n_iter = 0; // for debugging
    // could do this with pointers rather than ints too
    while (i < stack_size) {

	    n_iter++;
	    if (stack_size >= buff_size) {
		    printf("Error! stack_size > %i\n", buff_size);
	    }
	    if ( i > buff_size) {
		    printf("Error 2: i > %i\n", buff_size);
	    }
	    pi = node_id_buffer[i];
    	    node_t p = *(_tree_data+pi);
    	    n = norm2(&p.data, &datum);
	    this_field = p.data.field;

	    if( !p.has_lchild ) { 
	        count_node(n, rsq, counter, this_field, query_field, p.data.weight);
		++i;
	    } else if ( !p.has_rchild ) { 
	        count_node(n, rsq, counter, this_field, query_field, p.data.weight);
		node_id_buffer[i] = left_child(pi);
	    } else {
		val = _data_query[NDIM*qi+p.dim];
		pos_upper = val + rad;
		pos_lower = val - rad;
		point = p.data.value[p.dim];
		if (pos_upper < point) {
		    node_id_buffer[i] = left_child(pi);
		} else if (pos_lower > point) {
		    node_id_buffer[i] = right_child(pi);
		} else {
		    count_node(n, rsq, counter, this_field, query_field, p.data.weight);
		    if ( i ) {
			if (i > 1) {
			    node_id_buffer[i] = node_id_buffer[--stack_size];
			    node_id_buffer[--i] = right_child(pi);
			    node_id_buffer[--i] = left_child(pi);
			}

			else {
		            node_id_buffer[i] = right_child(pi);
			    node_id_buffer[--i] = left_child(pi);
			}
		    } else {
		        node_id_buffer[i] = left_child(pi);
		    	node_id_buffer[stack_size++] = right_child(pi);
                        if (stack_size > stack_max) { stack_max = stack_size; }
		    }
		}
	    }
    }
    free(node_id_buffer);
    //if (tid == 0) {
    //    printf("Buff Size: %i\n", buff_size);
    //	printf("i: %i, stack_size: %i\n", i, stack_size);
    // 	printf("n_iter: %i, stack_max: %i\n", n_iter, stack_max);
    //}
    return;
	
} 

/*
__device__ void radius(int pi, int qi, field_counter_t * counter, int tid)
{
    printf("%i: Entering Radius Function\n", tid);
    node_t p = *(_tree_data+pi);
    int this_field = p.data.field;

    int query_field = 0;
    if(_field_query != NULL) {
        query_field = _field_query[qi];
    }

    double val;

    double rsq;
    double pos_upper = 0.0, pos_lower = 0.0, point = 0.0;

    rsq = rad*rad;

    datum_t datum;
    int i;
    printf("%i: A\n", tid);
    for(i=0;i<NDIM;i++) {
        datum.value[i] = _data_query[NDIM*qi+i];
    }
    double n = norm2(&p.data, &datum);

    if( !p.has_lchild ) { 
        count_node(n, rsq, counter, this_field, query_field, p.data.weight);

    } else if ( !p.has_rchild ) { 
        count_node(n, rsq, counter, this_field, query_field, p.data.weight);
        radius(left_child(pi),qi,counter, tid);

    } else {
	val = _data_query[NDIM*qi+p.dim];
        pos_upper = val + rad;
        pos_lower = val - rad;
        point = p.data.value[p.dim];
        if (pos_upper < point) {
            radius(left_child(pi),qi,counter, tid);
	} else if (pos_lower > point) {
            radius(right_child(pi),qi,counter, tid);
        } else {
	    count_node(n, rsq, counter, this_field, query_field, p.data.weight);
            radius(left_child(pi),qi,counter, tid);
            radius(right_child(pi),qi,counter, tid);
        }
    }

    return;
} 
*/

/*
 * Query how many points in the tree with head p lie within radius r of point q
 */
/*
static void radius(int pi, int qi, field_counter_t * counter)
{
    node_t p = *(_tree_data+pi);
    int this_field = p.data.field;

    int query_field = 0;
    if(_field_query != NULL) {
        query_field = _field_query[qi];
    }

    double val;

    double rsq;
    double pos_upper = 0.0, pos_lower = 0.0, point = 0.0;

    rsq = rad*rad;

    datum_t datum;
    int i;
    for(i=0;i<NDIM;i++) {
        datum.value[i] = _data_query[NDIM*qi+i];
    }
    double n = norm2(&p.data, &datum);

    if( !p.has_lchild ) { / * no children * /
        count_node(n, rsq, counter, this_field, query_field, p.data.weight);
        return;

    } else if ( !p.has_rchild ) { / * one child * /
        count_node(n, rsq, counter, this_field, query_field, p.data.weight);
        radius(left_child(pi),qi,counter);
        return;

    } else { / * two children * /
        val = _data_query[NDIM*qi+p.dim];
        pos_upper = val + rad;
        pos_lower = val - rad;
        point = p.data.value[p.dim];

        if (pos_upper < point) {
            radius(left_child(pi),qi,counter);
        } else if (pos_lower > point) {
            radius(right_child(pi),qi,counter);
        } else {
            count_node(n, rsq, counter, this_field, query_field, p.data.weight);
            radius(left_child(pi),qi,counter);
            radius(right_child(pi),qi,counter);
        }
        return;
    }
    
}
*/
typedef struct thread_args {
    int rank;
    int start;
    int stop; // actually one after the last index to stop at
    int * finished;
    field_counter_t * counter;
} thread_args_t;

__device__ int _max(int a, int b)
{
    return a > b ? a : b;
}

__device__ int _min(int a, int b)
{
    return a < b ? a : b;
}
/*
static inline int max(int a, int b)
{
    return a > b ? a : b;
}

static inline int min(int a, int b)
{
    return a < b ? a : b;
}
*/
field_counter_t * init_field_counter(int qsize, int tsize)
{
    field_counter_t * c = (field_counter_t*)malloc(sizeof(field_counter_t));
    c->array = (long long *)calloc(qsize*tsize, sizeof(long long));
    c->num_fields = qsize;
    return c;
}

void free_field_counter(field_counter_t * c)
{
    free(c->array);
    free(c);
}

// for now I'll just rewrite this to run on mutiple GPU threads exactly as the pthreads run
// at least until I understand the code well enough

__global__ void static_allocator(node_t* d_tree_data, double* d_data_query, int* d_field_query) {

	_tree_data = d_tree_data;
	_data_query = d_data_query;
	_field_query = d_field_query;
}

__global__ void radius_on_device(int step, int length, int qsize, int tsize, long long* array) {

	int bid = blockIdx.x;
	int tid = threadIdx.x;
	int idx = bid*gridDim.x+tid;
	if (idx < length) {
		field_counter_t counter;// = (field_counter_t*)malloc(sizeof(field_counter_t));
		counter.num_fields = qsize;
        	counter.array = array+(qsize*tsize*tid);
                radius(1, idx, &counter, tid);
	}

	//Note: I could use barriers to collapse the counter arrays at the end
}
/*
void * radius_on_thread(void *voidargs)
{
    int i;
    thread_args_t *targs;
    targs = (thread_args_t *) voidargs;

    for(i=targs->start; i<targs->stop; i++) {
        radius(1, i, targs->counter);
    }

    *(targs->finished) = 1;
    return NULL;
}
*/
/* 
 * Compute the sum of the number of data points in the tree which lie within
 * radius r of each of the points in the array
 */
void status(cudaError_t cstatus, const char* msg) {
	if (cstatus != cudaSuccess) {printf("*Failed at: %s\n", msg); printf("*--%s\n",cudaGetErrorString(cstatus));}
}

void alloc_globals(node_t* h_tree_data, size_t ts, double* h_data_query, size_t ds, int* h_field_query, size_t fs) {
	node_t* d_tree_data;
	double* d_data_query;
	int* d_field_query;


	HANDLE_ERROR(cudaMalloc((void**)&d_tree_data, ts));
	HANDLE_ERROR(cudaMalloc((void**)&d_data_query, ds));
	HANDLE_ERROR(cudaMalloc((void**)&d_field_query, fs));

	HANDLE_ERROR(cudaMemcpy(d_tree_data, h_tree_data, ts, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_data_query, h_data_query, ds, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_field_query, h_field_query, fs, cudaMemcpyHostToDevice));

	static_allocator<<<1,1>>>(d_tree_data, d_data_query, d_field_query);
	
	cudaDeviceSynchronize();
	printf("Finished Allocating Data\n");
}

extern "C"
long long * pair_count(kdtree_t tree, double * data,
                                int * fields, int length, int num_fields,
                                double r, int num_threads)
{

    if (need_alloc) {
        need_alloc = 0; // no longer need globals allocated in the future 
        HANDLE_ERROR(cudaDeviceSetLimit(cudaLimitMallocHeapSize, size_t(2147483648)));
        alloc_globals(tree.node_data, tree.memsize*sizeof(node_t), data, 3*length*sizeof(double), fields, length*sizeof(int));
        HANDLE_ERROR(cudaMemcpyToSymbol(rad, &r, sizeof(double), size_t(0), cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpyToSymbol(tree_node_count, &tree.memsize, sizeof(int), size_t(0), cudaMemcpyHostToDevice));
    }

    int i, j;


    int step = length/num_threads;


    size_t fca_size = num_fields*tree.num_fields*sizeof(long long); //size of a field counter array

    field_counter_t* counters = (field_counter_t*)malloc(num_threads*sizeof(field_counter_t));
    long long* arrays = (long long*)calloc(num_fields*tree.num_fields*num_threads, sizeof(long long));
    // A flattened 2d array
    printf("---Creating the counters...\n");
    for(i=0; i<num_threads; i++) { // init the counters (could copying the data be faster?)
   	counters[i].array = arrays + i*num_fields*tree.num_fields;
   	counters[i].num_fields = num_fields;
    }

    printf("---Allocating Device Data... %lu B\n", fca_size*num_threads);
    long long* dev_arrays;
    HANDLE_ERROR(cudaMalloc((void**)&dev_arrays, fca_size*num_threads));


    HANDLE_ERROR(cudaMemcpy(dev_arrays, arrays, fca_size*num_threads, cudaMemcpyHostToDevice));

    printf("---Executing the Kernel...\n");
    printf("----Num Threads: %i...\n", num_threads);
    printf("----Num Fields: %i...\n", num_fields);
    printf("----Num Tree Fields: %i...\n", tree.num_fields);
    printf("----Data Length: %i...\n", length);
    printf("----Grid Size: %i...\n", step);
    printf("----Kernel Args: step: %i, length: %i, qsize: %i, tsize: %i\n", step, length, num_fields, tree.num_fields);

    dim3 grid(step+1);


    //execute radius kernel
    radius_on_device<<<grid,num_threads>>>(step, length, num_fields, tree.num_fields, dev_arrays);
    cudaDeviceSynchronize();



    printf("---Recovering the Device Data... %lu B\n", fca_size*num_threads);


    //This is cool, because pointers are awesome, our data is already in our counters!
    HANDLE_ERROR(cudaMemcpy(arrays, dev_arrays, fca_size*num_threads, cudaMemcpyDeviceToHost)); 

    printf("----counters[0].array[0] = %lli\n", counters[0].array[0]);
    printf("---Summing results...\n");
    // sum the array for each thread into one array
    long long * results = (long long*)calloc(num_fields*tree.num_fields, sizeof(long long));

    for(i=0; i<num_threads; i++) {
        for(j=0; j<num_fields*tree.num_fields; j++) {
            results[j] += counters[i].array[j];
        }
    }
	printf("--Result: %lli\n", results[0]);

	printf("---Freeing Memory\n");
    free(counters);
    free(arrays);
    cudaFree(dev_arrays);
    printf("--Done\n");
    return results;
}


/* 
 * Compute the sum of the number of data points in the tree which lie within
 * radius r of each of the points in the array
 */
/*
long long * pair_count(kdtree_t tree, double * data,
                                int * fields, int length, int num_fields,
                                double r, int num_threads)
{
    int i, j;
    int next_assign = 0;
    _tree_data = tree.node_data;
    _data_query = data;
    _field_query = fields;
    rad = r;

    if(num_threads > length / 20) {
        printf("More than %d threads!\n", length / 20);
        exit(EXIT_FAILURE);
    }

    int step = max(1000, length/num_threads);
    
    pthread_t threads[num_threads];
    thread_args_t targs[num_threads];
    int * started = calloc(num_threads, sizeof(int));
    int * finished = calloc(num_threads, sizeof(int));
    for(i=0; i<num_threads; i++) {
        targs[i].counter = init_field_counter(num_fields, tree.num_fields);

        // give each thread a unique index, like an MPI worker
        targs[i].rank = i;
        targs[i].start = next_assign;
        targs[i].stop = min(length, next_assign+step);
        next_assign = min(length, next_assign+step);

        targs[i].finished = finished+i;
    }

    // create initial threads
    for(i=0; i<num_threads; i++) {
        pthread_create(threads+i, NULL, radius_on_thread, targs+i);
        started[i] = 1;
    }

    // while there are queries to do, create new threads for `step` queries
    while(1) {

        for(i=0; i<num_threads; i++) {
            if(finished[i]) {
                pthread_join(threads[i], NULL);
                finished[i] = 0;
                started[i] = 0;

                if (next_assign >= length) break;

                targs[i].start = next_assign;
                targs[i].stop = min(length, next_assign+step);
                next_assign = min(length, next_assign+step);

                started[i] = 1;
                pthread_create(threads+i, NULL, radius_on_thread, targs+i);
            }
        }

        if (next_assign >= length) break;
        usleep(500);
    }

    // join the last of the threads
    for(i=0; i<num_threads; i++) {
        if(started[i]) {
            pthread_join(threads[i], NULL);
        }
    }

    // sum the array for each thread into one array
    long long * results = calloc(num_fields * tree.num_fields, sizeof(long long));

    for(i=0; i<num_threads; i++) {
        for(j=0; j<num_fields*tree.num_fields; j++) {
            results[j] += ((targs[i].counter)->array)[j];
        }
        free_field_counter(targs[i].counter);
    }

    free(started);
    free(finished);
    return results;
}
*/
