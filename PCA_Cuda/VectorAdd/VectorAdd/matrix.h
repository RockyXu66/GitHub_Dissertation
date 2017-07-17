#pragma once

class matrix {
public:
	double **data;
	int i;
	int j;
	matrix() {}
	matrix(int i, int j) {
		this->i = i;
		this->j = j;
		data = new double*[i];
		for (int k = 0; k < i; k++) {
			data[k] = new double[j];
		}
		// Initialize matrix
		for (int x = 0; x < i; x++) {
			for (int y = 0; y < j; y++) {
				data[x][y] = 0;
			}
		}
	}
	void free() {
		for (int x = 0; x < i; x++) {
			delete[] data[x];
		}
		delete[] data;
	}
};