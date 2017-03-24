#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#ifndef KDTREE_H
#include "kdtree.h"
#endif

static double *_x_query, *_y_query, *_z_query;

inline int min(int a, int b)
{
    return a < b ? a : b;
}

static inline FLOAT norm2(FLOAT a, FLOAT b, FLOAT c)
{
    return a*a+b*b+c*c;
}

void destroy(node_t *p)
{
    if(p->lchild != NULL) destroy(p->lchild);
    if(p->rchild != NULL) destroy(p->rchild);
    free(p);
}

/*
 * Query how many points in the tree with head p lie within radius r of point
 * (x, y, z). Recursive.
 */
static int radius(node_t *p, enum dim d, int index, FLOAT r)
{

    d=d%3;

    int i;
    double x, y, z;

    x = _x_query[index]; y = _y_query[index]; z = _z_query[index];

    FLOAT rsq, dx, dy, dz;
    FLOAT pos_upper, pos_lower, point;

    rsq = r*r;
    dx = p->x - x; dy = p->y - y; dz = p->z - z;

    FLOAT n = norm2(dx,dy,dz);

    if(p->lchild == NULL && p->rchild == NULL) { /* no children */
        return (n < rsq);
    } else if (p->rchild == 0) { /* one child */
        return (n < rsq) + radius(p->lchild,d+1,index,r);
    } else { /* two children */

        switch(d) {
        case X:
            pos_upper = x+r;
            pos_lower = x-r;
            point = p->x;
            break;
        case Y:
            pos_upper = y+r;
            pos_lower = y-r;
            point = p->y;
            break;
        case Z:
            pos_upper = z+r;
            pos_lower = z-r;
            point = p->z;
            break;
        }

        if (pos_upper < point) {
            i=radius(p->lchild,d+1,index,r);
        } else if (pos_lower > point) {
            i=radius(p->rchild,d+1,index,r);
        } else {
            i = (n < rsq) +
                radius(p->lchild,d+1,index,r) +
                radius(p->rchild,d+1,index,r);
        }
        return i;
    }
    
}

typedef struct shared_args {
    kdtree_t tree;
    FLOAT *x, *y, *z;
    int n;
    int node_start;
    int node_stop;
    int num_threads;
    long long *sum;
    FLOAT r;
} shared_args_t;

typedef struct thread_args {
    int thread_rank;
    struct shared_args * shared_args_p;
} thread_args_t;


void assign_idx(int this_thread, int num_threads, int n_array, int * start,
                                                                int * stop)
{

    int mpi_rank, mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    
    int size = mpi_size * num_threads;
    int rank = mpi_rank * num_threads + this_thread;

    int n_array_local = n_array / size + ( n_array % size > rank );

    *start = rank * ( n_array / size ) + min( rank, n_array % size );
    *stop = *start + n_array_local;

}

void * twopoint_wrap(void *voidargs)
{

    int start, stop, i;
    thread_args_t *targs;
    targs = (thread_args_t *) voidargs;
    shared_args_t * args;
    args = targs->shared_args_p;
    
    int rank = targs->thread_rank;
    assign_idx(rank, args->num_threads, args->n, &start, &stop);

    args->sum[rank] = 0;

    for(i=start; i<stop; i++) {
        args->sum[rank] += radius(args->tree.root, 0, i, args->r);
    }

    return NULL;
}

/* 
 * Compute the sum of the number of data points in the tree which lie within
 * radius r of each of the (x,y,z) points in the array. Result may easily
 * exceed the size of a 32-bit int, so we return a long long.
 */
long long two_point_correlation(kdtree_t tree, FLOAT x[], FLOAT y[],
                FLOAT z[], int n, FLOAT r, int num_threads, MPI_Comm comm)
{
    int i, rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    _x_query = x; _y_query = y; _z_query = z;

    shared_args_t ss;
    ss.tree = tree;
    ss.n = n;
    ss.r = r;
    ss.num_threads = num_threads;
    ss.sum = (long long *)malloc(num_threads * sizeof(long long));
    
    thread_args_t targs[num_threads];
    for(i=0; i<num_threads; i++) {
        targs[i].shared_args_p = &ss;
        targs[i].thread_rank = i;
    }

    pthread_t threads[num_threads];
    threads[0] = pthread_self();

    long long result = 0;

    double t1, t2;
    t1 = MPI_Wtime();

    for(i=1; i<num_threads; i++) 
        pthread_create(threads+i, NULL, twopoint_wrap, targs+i);

    twopoint_wrap(targs);

    for(i=1; i<num_threads; i++)
        pthread_join(threads[i], NULL);

    for (i=0; i<num_threads; i++)
        result += ss.sum[i];

    MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);

    t2 = MPI_Wtime();

    if(!rank) {
        printf("Time on rank 0: %f sec\n", t2 - t1);
        printf("Sum: %lld\n", result);
    }

    return result;
}
