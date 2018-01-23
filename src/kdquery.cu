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
	#include "kdtree.h"
}

void *__gxx_personality_v0;

__device__ double * _data_query;
__device__ int *_field_query;
__device__ node_t * _tree_data;
__device__ double rad;

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
        counter->array[query_field * counter->num_fields + this_field] += this_weight;
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
__device__ void radius(int pi, int qi, field_counter_t * counter)
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

    if( !p.has_lchild ) { /* no children */
        count_node(n, rsq, counter, this_field, query_field, p.data.weight);
        return;

    } else if ( !p.has_rchild ) { /* one child */
        count_node(n, rsq, counter, this_field, query_field, p.data.weight);
        radius(left_child(pi),qi,counter);
        return;

    } else { /* two children */
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

__global__ void radius_on_device(int step, int length, field_counter_t* field_counters) {

	//I need the counter (struct), rank (or similar), start, stop, finished (ints)
	//counter: pair counter (qsize is number of fields, array is counter for each field 2d array)
	//rank: thread rank duh
	//I need to do the thread initialization here too
	//the field counters should be input to the function
	//have to partition the data to a smaller set of threads

	int rank = threadIdx.x;
	field_counter_t* counter = field_counters+rank; 
	int start = rank*step;
	int stop = _min((rank+1)*step, length);
	int i;
	for (i=start; i < stop; ++i) { //got to make sure we dont run out of threads
		radius(1, i, counter);
	}
	counter->array[0] = 50;
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

extern "C"
long long * pair_count(kdtree_t tree, double * data,
                                int * fields, int length, int num_fields,
                                double r, int num_threads)
{
    int i, j;
    //int next_assign = 0;

    //allocate the gloabal variables
    cudaError_t cudaStatus;
		cudaStatus = cudaMemcpyToSymbol("_tree_data", &tree.node_data[0], tree.memsize, size_t(0), cudaMemcpyHostToDevice);
		status(cudaStatus, "Copying Symbols");
    cudaMemcpyToSymbol("_data_query", &data[0], 3*length*sizeof(double), size_t(0), cudaMemcpyHostToDevice); 
    //I think the '3' in the previous line is NDIM, but not sure, good to know if I generalize the code
    cudaMemcpyToSymbol("_field_query", &fields[0], sizeof(int)*num_fields, size_t(0), cudaMemcpyHostToDevice);
    cudaStatus = cudaMemcpyToSymbol(rad, &r, sizeof(double), size_t(0), cudaMemcpyHostToDevice);
		status(cudaStatus, "Copying rad");

    /*
    _tree_data = tree.node_data;
    _data_query = data;
    _field_query = fields;
    rad = r;
    */

    if(num_threads > length / 20) {
        printf("More than %d threads!\n", length / 20);
        exit(EXIT_FAILURE);
    }

    printf("--Beginning my Edits\n");
    int step = length/num_threads;
    size_t fca_size = num_fields*tree.num_fields*sizeof(long long); //size of a field counter array

    field_counter_t* counters = (field_counter_t*)malloc(num_threads*sizeof(field_counter_t));

    printf("---Creating the counters...\n");
    for(i=0; i<num_threads; i++) { // init the counters (could copying the data be faster?)
   	counters[i].array = (long long *)calloc(num_fields*tree.num_fields, sizeof(long long));
   	counters[i].num_fields = num_fields;
    }

    printf("---Allocating Device Data...\n");
    field_counter_t* dev_counter_arrays = (field_counter_t*)malloc(num_threads*sizeof(field_counter_t));
    memcpy(dev_counter_arrays, counters, num_threads*sizeof(field_counter_t)); //copy the matrix data

    for(i=0; i<num_threads; i++) {
    		cudaStatus = cudaMalloc(&(dev_counter_arrays[i].array), fca_size); // allocate on the GPU the array
				status(cudaStatus, "Malloc This Thing");
      	cudaStatus = cudaMemcpy(dev_counter_arrays[i].array, counters[i].array, fca_size, cudaMemcpyHostToDevice);	
				status(cudaStatus, "Memcpy the array");
    }

    field_counter_t* dev_counters; //the full counter data to be sent to the GPU
    cudaMalloc(&dev_counters, num_threads*sizeof(field_counter_t));
    cudaMemcpy(dev_counters, dev_counter_arrays, num_threads*sizeof(field_counter_t),cudaMemcpyDeviceToDevice); 

    printf("---Executing the Kernel...\n");
    printf("----Num Threads: %i...\n", num_threads);
    printf("----Num Fields: %i...\n", length);
    //execute radius kernel
    radius_on_device<<<num_threads,1>>>(step, length, dev_counters);

    printf("---Recovering the Device Data...\n");

    cudaMemcpy(dev_counter_arrays, dev_counters, num_threads*sizeof(field_counter_t), cudaMemcpyDeviceToDevice); 
    for(i=0; i<num_threads; i++) {
      	cudaMemcpy(counters[i].array, dev_counter_arrays[i].array, fca_size, cudaMemcpyDeviceToHost);	
    }
		printf("----counters[0].array[0] = %lli\n", counters[0].array[0]);
    printf("---Summing results...\n");
    // sum the array for each thread into one array
    long long * results = (long long*)calloc(num_fields*tree.num_fields, sizeof(long long));

    for(i=0; i<num_threads; i++) {
        for(j=0; j<num_fields*tree.num_fields; j++) {
            results[j] += counters[i].array[j];
        }
				free(counters[i].array);
        //free_field_counter(counters[i]);
    }
	printf("--Result: %lli\n", results[0]);

	printf("---Freeing Memory\n");
    free(counters);
    cudaFree(dev_counter_arrays);
    cudaFree(dev_counters);
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
