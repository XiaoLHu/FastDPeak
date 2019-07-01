#include "stack.h"
#include<stdio.h>
#include<stdlib.h>

typedef float* point;

float complete_distance(point v1, point v2);
float distance(point v1, point v2, float upper_bound);
v_array<point > parse_points(FILE *input);
void print(point &p);
int posix_memalign(void **mptr, size_t, size_t bytes);
v_array<point > parse_points(float *input, int data_size, int dim);
