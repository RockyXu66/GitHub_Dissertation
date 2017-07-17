//
//  PCA.hpp
//  PCA_opengl_0522
//
//  Created by Yinghan Xu on 22/05/2017.
//  Copyright © 2017 Yinghan Xu. All rights reserved.
//

#ifndef PCA_hpp
#define PCA_hpp

#include <iostream>
#include <fstream>
#include <ctime>
#include <U:/Eigen/Eigen>
using namespace Eigen;
using namespace std;

class PCA {

public:
	int num;        // Number of images
	int cellImage_num; //Number of cells
	int dim;        // Dimension of one image
	int pca_dim;    // PCA dimension of the image
	string fileName;
	MatrixXd EigenVectors_Red, EigenVectors_Green, EigenVectors_Blue;
	MatrixXd oneImage_scores_Red, oneImage_scores_Green, oneImage_scores_Blue;
	vector<MatrixXd> wholeImage_scores_Red, wholeImage_scores_Green, wholeImage_scores_Blue;
	vector<MatrixXd> oneVecs_Red, oneVecs_Green, oneVecs_Blue;
	MatrixXd means_Red, means_Green, means_Blue;


	unsigned char* CalculateEigen(vector<unsigned char*> images, int width, int height);

	void Mean(MatrixXd data_r, MatrixXd data_g, MatrixXd data_b, unsigned char* meanImage);

	MatrixXd Eigen(MatrixXd data, unsigned char* newImage, string name, string fileName);

	MatrixXd Reconstruction(MatrixXd data, string fileName);

	vector<MatrixXd> WholeReconstruction(string fileName, string meanFileName, string scoresFileName, int imageIndex, string channel);

	void SaveHighEigen(string newFileName, string oldFileName, int highDimension);

	void Mean(vector<MatrixXd> data, string fileName);

	void CalScores(vector<MatrixXd> whole_data, string fileName, string meanFileName, string newFileName);

	unsigned char* Reconstruction_RealTime(unsigned char* image, int imageIndex, int width, int height);

	vector<MatrixXd> WholeReconstruction_RealTime(string meanFileName, string scoresFileName, int imageIndex, string channel);

private:
	// Used for sortting eigen value
	MatrixXd BubbleSort(MatrixXd eigenValues, int num);

};

#endif /* PCA_hpp */