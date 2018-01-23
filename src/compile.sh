$(gcc -I/usr/lcoal/cuda/include -c kdbench.c);

$(nvcc -c kdquery.cu);

$(gcc -L/usr/local/cuda/lib64 -o test kdbuild.c kdtest.c kdbench.c kdquery.o -lcuda -lcudart -lm);
