//#include <iostream>
//#include <thrust/host_vector.h> 
//#include <thrust/device_vector.h>
//
//using namespace std;
//using namespace thrust;
//
//#define dev_dim 20
//__global__ void arradd(double* d_m, double* d_n, double size, double projectionValue) {
//	int myid = threadIdx.x;
//
//	d_n[myid] += projectionValue * d_m[myid];
//}
//
//void cudaAdd(double m[dev_dim], double n[dev_dim], double size, double projectionValue) {
//
//	double *d_m, *d_n;
//
//	cudaMalloc(&d_m, size);
//	cudaMalloc(&d_n, size);
//
//	cudaMemcpy(d_m, m, size, cudaMemcpyHostToDevice);
//
//	dim3 DimGrid(1, 1);
//	dim3 DimBlock(dev_dim, 1);
//
//	arradd << < DimGrid, DimBlock >> > (d_m, d_n, size, projectionValue);
//	// arradd <<< 1, 200 >>> (d_m, d_n, d_p, size);
//
//	cudaMemcpy(n, d_n, size, cudaMemcpyDeviceToHost);
//
//	cudaFree(d_m);
//	cudaFree(d_n);
//}
//
//int main() {
//
//	vector<vector<double>> data1;
//	vector<double> data11;
//	for (int i = 0; i < 20; i++) {
//		data11.push_back(i);
//	}
//	data1.push_back(data11);
//	data1.push_back(data11);
//
//	host_vector<double> H(data11.begin(), data11.end());
//	device_vector<double> D = H;
//	
//	// extract raw pointer from device_ptr//	//double * raw_ptr = thrust::raw_pointer_cast(D);
//	size_t N = 10;//	// create a device_ptr//	thrust::device_ptr<int> dev_ptr = thrust::device_malloc<int>(N);//	// extract raw pointer from device_ptr//	int * raw_ptr = thrust::raw_pointer_cast(dev_ptr);
//	//cout << raw_ptr[10] << endl;
//	raw_ptr[10] = 100;
//	cout << D[10] << endl;
//
//
//	double projectionValue = 0.3;
//	double a[dev_dim], b[dev_dim];
//	for (int j = 0; j < 20; j++) {
//		a[j] = j;// oneVec.data[j][i];
//		b[j] = 0; // OriginalData.data[0][j];
//	}
//	double size = dev_dim * sizeof(double);
//	cudaAdd(a, b, size, projectionValue);
//
//	for (int i = 0; i < 20; i++) {
//		cout << b[i] << " ";
//	}
//	cout << endl;
//
//	int x;
//	cin >> x;
//	return 0;
//}


//#include <thrust/device_vector.h>
//#include <thrust/transform.h>
//#include <thrust/functional.h>
//#include <thrust/copy.h>
//#include <iostream>
//
//using namespace std;
//using namespace thrust;
//using namespace thrust::placeholders;
//
//// SAXPY
//struct saxpy {
//	double a;
//	saxpy (double a) : a(a){}
//	__host__ __device__ 
//		double operator()(double x) {
//		return a*x;
//	}
//};
//
//int main(void)
//{
//	// Vector addition Z = X + Y
//	thrust::device_vector<float> X(3);
//	thrust::device_vector<float> Y(3);
//	thrust::device_vector<float> Z(3);
//	X[0] = 10; X[1] = 20; X[2] = 30;
//	Y[0] = 15; Y[1] = 35; Y[2] = 10;
//	thrust::transform(X.begin(), X.end(),
//		Y.begin(),
//		Z.begin(),
//		thrust::plus<float>());
//	for (size_t i = 0; i < Z.size(); i++)
//		std::cout << "Z[" << i << "] = " << Z[i] << "\n";
//
//	// Sum of a vector   result = sum(M)
//	vector<double> test;
//	test.push_back(2);
//	test.push_back(1);
//	test.push_back(8);
//	device_vector<double> M = test;
//	//M[0] = 2; M[1] = 1; M[2] = 10;
//	double result = reduce(M.begin(), M.end());
//	cout << "sum is " << result << endl;
//
//	// SAXPY   N = a * M
//	device_vector<double> N(3);
//	double a = 3;
//	transform(M.begin(), M.end(),
//		N.begin(),
//		saxpy(a));
//	for (size_t i = 0; i < N.size(); i++) {
//		cout << "N[" << i << "] = " << N[i] <<" ";
//	}
//	cout << endl;
//	
//	for (int i = 0; i < 3; i++) {
//		transform(M.begin(), M.end(),
//			N.begin(),
//			N.begin(),
//			_1 + _2);
//	}
//	
//
//	
//	for (size_t i = 0; i < N.size(); i++) {
//		cout << "N[" << i << "] = " << N[i] << " ";
//	}
//	cout << endl;
//	host_vector<double> tmp = N;
//	thrust:: copy(N.begin(), N.end(), test.begin());
//	for (size_t i = 0; i < N.size(); i++) {
//		cout << "N[" << i << "] = " << test[i] << " ";
//	}
//	
//	int x;
//	cin >> x;
//	return 0;
//}





//#include <thrust/device_vector.h>
//
//thrust::device_vector<int> iVec;
//
//int* iArray = thrust::raw_pointer_cast(&iVec[0]);
//
//fooKernel << <x, y >> >(iArray);
//
//// Template structure to pass to kernel
//template <typename T>
//struct KernelArray
//{
//	T*  _array;
//	int _size;
//};
//
//// Function to convert device_vector to structure
//template <typename T>
//KernelArray<T> convertToKernel(thrust::device_vector<T>& dVec)
//{
//	KernelArray<T> kArray;
//	kArray._array = thrust::raw_pointer_cast(&dVec[0]);
//	kArray._size = (int)dVec.size();
//
//	return kArray;
//}
//
//thrust::device_vector<int> iVec;
//
//fooKernel << <x, y >> >(convertToKernel(iVec)); // Explicit conversion from iVec to KernelArray<int>
//
//__global__ fooKernel(KernelArray<int> inArray)
//{
//	for (int i = 0; i < inArray._size; ++i)
//		something = inArray._array[i];
//	// ...
//	return;
//}

#include <iostream>
#include <vector>
using namespace std;
#define dev_dim 10

typedef double arr_dim[dev_dim];

__global__ void arradd(double* d_m, double* d_n, double size) {
	int myid = threadIdx.x;

	d_n[myid] = d_n[myid] + d_m[myid];
}

int main() {
	//arr_dim Ta;
	//arr_dim Tb;
	//double Ta[dev_dim];
	//double Tb[dev_dim];

	//for (int i = 0; i < dev_dim; i++) {
	//	Ta[i] = 0;
	//	Tb[i] = i + 1;
	//}

	// Test
	vector<vector<double *>> test1, test2;
	vector<double *> test11, test22;
	double test111[dev_dim], test222[dev_dim];
	for (int i = 0; i < dev_dim; i++) {
		test111[i] = i;
		test222[i] = 0;
	}
	test11.push_back(test111);
	test22.push_back(test222);
	test1.push_back(test11);
	test2.push_back(test22);


	double *d_m, *d_n;
	double size = dev_dim * sizeof(double);

	cudaMalloc(&d_m, size);
	cudaMalloc(&d_n, size);

	cudaMemcpy(d_n, test2[0][0], size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_m, test1[0][0], size, cudaMemcpyHostToDevice);

	dim3 DimGrid(1, 1);
	dim3 DimBlock(dev_dim, 1);
	arradd << <DimGrid, DimBlock >> > (d_m, d_n, size);
	cudaMemcpy(test2[0][0], d_n, size, cudaMemcpyDeviceToHost);

	//for (int i = 0; i < dev_dim; i++) {
		//cout << Ta[i] << " ";
	//}
	//cout << endl;

	cudaFree(d_m);
	cudaFree(d_n);

	
	

	for (int i = 0; i < dev_dim; i++) {
		test1[0][0][i] += 5;
		cout << test1[0][0][i] << " ";
	}
	cout << endl;

	// Copy
	test2 = test1;
	for (int i = 0; i < dev_dim; i++) {
		cout << test2[0][0][i] << " ";
	}
	cout << endl;
	

	int x;
	cin >> x;
	return 0;
}

