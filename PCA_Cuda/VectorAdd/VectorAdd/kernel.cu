#include "kernel.h"
#define dev_dim 400
__global__ void arradd(double* d_m, double* d_n, double size, double projectionValue) {
	int myid = threadIdx.x;

	d_n[myid] = d_n[myid] + projectionValue * d_m[myid];
}
__global__ void arradd2(double* d_m, double* d_n, double size) {
	int myid = threadIdx.x;

	d_n[myid] += d_m[myid];
}

void cudaAdd(double m[dev_dim], double n[dev_dim], double size, double projectionValue) {

	double *d_m, *d_n;

	cudaMalloc(&d_m, size);
	cudaMalloc(&d_n, size);

	cudaMemcpy(d_m, m, size, cudaMemcpyHostToDevice);

	dim3 DimGrid(1, 1);
	dim3 DimBlock(dev_dim, 1);

	arradd <<< DimGrid, DimBlock >>> (d_m, d_n, size, projectionValue);
	// arradd <<< 1, 200 >>> (d_m, d_n, d_p, size);

	cudaMemcpy(n, d_n, size, cudaMemcpyDeviceToHost);

	cudaFree(d_m);
	cudaFree(d_n);
}

void cudaAdd2(double m[dev_dim], double n[dev_dim], double size) {

	double *d_m, *d_n;

	cudaMalloc(&d_m, size);
	cudaMalloc(&d_n, size);

	cudaMemcpy(d_m, m, size, cudaMemcpyHostToDevice);

	dim3 DimGrid(1, 1);
	dim3 DimBlock(dev_dim, 1);

	arradd2 << < DimGrid, DimBlock >> > (d_m, d_n, size);
	// arradd <<< 1, 200 >>> (d_m, d_n, d_p, size);

	cudaMemcpy(n, d_n, size, cudaMemcpyDeviceToHost);

	cudaFree(d_m);
	cudaFree(d_n);
}



host_vector<Vmatrix> cudaReconstruction(vector<Arr_dim> cuda_means, Vmatrix cuda_oneImage_scores, vector<vector<Arr_dim>> cuda_oneVecs, int cellImage_num, int dim, int pca_dim) {
	cout << "Start cuda now" << endl;
	host_vector<Vmatrix> whole_OriginalData;
	double size = dev_dim * sizeof(double);
	double *d_m, *d_n;

	cudaMalloc(&d_m, size);
	cudaMalloc(&d_n, size);
	dim3 DimGrid(1, 1);
	dim3 DimBlock(dev_dim, 1);

	for (int cellIndex = 0; cellIndex < cellImage_num; cellIndex++) {
		Vmatrix OriginalData(1, dim);

		double d_OriginalData[dev_dim];
		for (int i = 0; i < dim; i++) {
			d_OriginalData[i] = 0;
		}

		cudaMemcpy(d_n, d_OriginalData, size, cudaMemcpyHostToDevice);
		for (int i = 0; i < pca_dim; i++) {
			double projectionValue = cuda_oneImage_scores.data[cellIndex][i];
			
			cudaMemcpy(d_m, cuda_oneVecs[cellIndex][i].data, size, cudaMemcpyHostToDevice);
			
			arradd <<< DimGrid, DimBlock >>> (d_m, d_n, size, projectionValue);
		}

		cudaMemcpy(d_m, cuda_means[cellIndex].data, size, cudaMemcpyHostToDevice);
		arradd2 << < DimGrid, DimBlock >> > (d_m, d_n, size);
		cudaMemcpy(d_OriginalData, d_n, size, cudaMemcpyDeviceToHost);

		for (int i = 0; i < dim; i++) {
			OriginalData.data[0][i] = d_OriginalData[i];
		}
		whole_OriginalData.push_back(OriginalData);

	}

	cudaFree(d_m);
	cudaFree(d_n);

	return whole_OriginalData;
}





















