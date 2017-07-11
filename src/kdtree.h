/*  Andrew Pellegrino, 2017
 *
 *  K-d tree building algorithm adapted from Russell A. Brown, Building a
 *  Balanced k-d Tree in O(kn log n) Time, Journal of Computer Graphics
 *  Techniques (JCGT), vol. 4, no. 1, 50-68, 2015
 */

#define KDTREE_H

#define HAS_LCHILD (1<<0)
#define HAS_RCHILD (1<<1)

#define ID_MASK_MAXINT ((1<<9) - 1)
// bits 2 through 10 (indexing from 0)
#define ID_MASK (ID_MASK_MAXINT << 2)

#define NDIM 3

typedef struct datum {
    double value[NDIM];
} datum_t;

typedef struct node {
    //double data[NDIM];
    datum_t data;
    unsigned int has_lchild:1;
    unsigned int has_rchild:1;
    unsigned int id:8;
    unsigned int dim:8;
} node_t;

typedef struct field_counter {
    long long * array;
    int size;
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

int left_child(int);
int right_child(int);

kdtree_t tree_construct(datum_t *, int *, int, int);

long long * pair_count_jackknife(kdtree_t, datum_t *, int *, int, int, double, int);
long long * pair_count_ftf(kdtree_t, datum_t *, int *, int, int, double, int);
long long * pair_count_noerr(kdtree_t, datum_t *, int, double, int);
