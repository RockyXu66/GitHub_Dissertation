#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <thrust/host_vector.h> 
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include "matrix.h"
#include "Vmatrix.h"
#include <iostream>

#define dev_dim 400

struct Arr_dim{
	double data[dev_dim];
};

using namespace std;
using namespace thrust;
using namespace thrust::placeholders;

void cudaAdd(double m[dev_dim], double n[dev_dim], double size, double projectionValue);

void cudaAdd2(double m[dev_dim], double n[dev_dim], double size);

host_vector<Vmatrix> cudaReconstruction(vector<Arr_dim> cuda_means, Vmatrix cuda_oneImage_scores, vector<vector<Arr_dim>> cuda_oneVecs, int cellImage_num, int dim, int pca_dim);
//void cudaReconstruction(matrix cuda_means, matrix cuda_oneImage_scores, host_vector<matrix> cuda_oneVecs, int cellImage_num, int dim, int pca_dim);

#pragma once
