/*  Andrew Pellegrino, 2017
 *
 *  K-d tree building algorithm adapted from Russell A. Brown, Building a
 *  Balanced k-d Tree in O(kn log n) Time, Journal of Computer Graphics
 *  Techniques (JCGT), vol. 4, no. 1, 50-68, 2015
 */

#ifndef KDTREE_H
#define KDTREE_H

#define NDIM 3

typedef struct datum {
    double value[NDIM];
    unsigned int field:8;
    unsigned int weight:8;
} datum_t;

typedef struct node {
    //double data[NDIM];
    datum_t data;
    unsigned int has_lchild:1;
    unsigned int has_rchild:1;
    unsigned int dim:8;
} node_t;

typedef struct field_counter {
    long long total;
    long long * array;
    int num_fields;
} field_counter_t;

typedef struct kdtree {
    node_t * node_data;
    int size;
    int memsize;
    int num_fields;
    datum_t * data;
    int * fields;
    int * args[NDIM];
} kdtree_t;

kdtree_t tree_construct(double *, int *, int, int);

long long * pair_count(kdtree_t, double *, int *, int, int, double, int);

#endif /* KDTREE_H */
