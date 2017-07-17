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
#include <GLFW/glfw3.h>
#include <Eigen/Eigen>
//#include <opencv2/opencv.hpp>
using namespace Eigen;
using namespace std;
//using namespace cv;

class PCA_{
    
public:
    int num;        // Number of images
    int cellImage_num; //Number of cells
    int dim;        // Dimension of one image
    int pca_dim;    // PCA dimension of the image
    string fileName;
    
    unsigned char* CalculateEigen(vector<unsigned char*> images, int image_num, int width, int height);
    
    unsigned char* Mean(MatrixXd data_r, MatrixXd data_g, MatrixXd data_b, unsigned char* meanImage);
    
    void Eigen(MatrixXd data, string name, string fileName);
    
    MatrixXd Reconstruction(MatrixXd data, string fileName);
    
    vector<MatrixXd> WholeReconstruction(string fileName, string meanFileName, string scoresFileName, int imageIndex);
    
    void SaveHighEigen(string newFileName, string oldFileName, int highDimension);
    
    void Mean(vector<MatrixXd> data, string fileName);
    
    void CalScores(vector<MatrixXd> whole_data, string fileName, string meanFileName, string newFileName);
    
    void Cal(string fileName, string meanFileName, string newFileName);
    
private:
    // Used for sortting eigen value
    MatrixXd BubbleSort(MatrixXd eigenValues, int num);
    
};

#endif /* PCA_hpp */





