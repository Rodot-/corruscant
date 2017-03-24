#include "kdtree.h"
#include <stdio.h>

int nodesize(void) {
    return sizeof(node_t);
}

void verify(node_t *root, enum dim d) {
	d=d%3;
    char *errmsg = "node %p should not be %s child of %p (%lf, %lf) \n";

	if (root->lchild != NULL) {
		switch(d) {
        case X:
            if(root->x < root->lchild->x) {
                fprintf(stderr,errmsg,
                (void*)root->lchild,"left",(void*)root);
            }
            break;
        case Y:
            if(root->y < root->lchild->y) {
                fprintf(stderr,errmsg,
                (void*)root->lchild,"left",(void*)root);
            }
            break;
        case Z:
            if(root->z < root->lchild->z) {
                fprintf(stderr,errmsg,
                (void*)root->lchild,"left",(void*)root);
            }
            break;
		}
		verify(root->lchild,d+1);
	}
	if (root->rchild != NULL) {
		switch(d) {
        case X:
            if(root->x > root->rchild->x) {
                fprintf(stderr,errmsg,
                (void*)root->lchild,"left",(void*)root);
            }
            break;
        case Y:
            if(root->y > root->rchild->y) {
                fprintf(stderr,errmsg,
                (void*)root->lchild,"left",(void*)root);
            }
            break;
        case Z:
            if(root->z > root->rchild->z) {
                fprintf(stderr,errmsg,
                (void*)root->lchild,"left",(void*)root);
            }
            break;
		}
		verify(root->rchild,d+1);
    }
}

int count(node_t *p) {
    if(p->rchild == NULL) return (p->lchild != NULL) + 1;
    return (count(p->lchild) + count(p->rchild) + 1);
}
