
#include <iostream>
#include <vector>
using namespace std;

class Vmatrix {
public:
	vector<vector<double>> data;
	Vmatrix(){}
	Vmatrix(int m, int n) {
		for (int i = 0; i < m; i++) {
			vector<double> row;
			for (int j = 0; j < n; j++) {
				row.push_back(0);
			}
			data.push_back(row);
		}
	}
	
};


#pragma once
