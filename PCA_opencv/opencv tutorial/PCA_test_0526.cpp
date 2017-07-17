////
////  PCA_test_0526.cpp
////  PCA_opencv
////
////  Created by Yinghan Xu on 26/05/2017.
////  Copyright © 2017 Yinghan Xu. All rights reserved.
////
//
//#include <iostream>
//#include <opencv2/opencv.hpp>
//
//// GLM Mathemtics
//#include <glm/glm.hpp>
//#include <glm/gtc/matrix_transform.hpp>
//#include <glm/gtc/type_ptr.hpp>
//#include <Eigen/Dense>
//#include <Eigen/Core>
//#include "ImageVector.h"
//#include "PCA_Demension.hpp"
//
//using namespace std;
//using namespace cv;
//using namespace glm;
//using namespace Eigen;
//
//float Cov(VectorXf X, VectorXf Y);
//
//int main(){
//    
//    // Load image
//    int image_num = 3;
//    Mat img1 = imread("image_face/face1.png");
//    Mat img2 = imread("image_face/face2.png");
//    Mat img3 = imread("image_face/face3.png");
//
//    
//    // Convert to grayscale
//    cvtColor(img1, img1, CV_BGR2GRAY);
//    cvtColor(img2, img2, CV_BGR2GRAY);
//    cvtColor(img3, img3, CV_BGR2GRAY);
//    Mat img_mean = img1.clone();
//    
//    // Check if image is loaded successfully
//    if(!img1.data || img1.empty())
//    {
//        cout << "Problem loading image!!!" << endl;
//        return EXIT_FAILURE;
//    }
//    
//    int width =  img1.size().width;
//    int height = img1.size().height;
//    
//    int length = width*height;
//    
//    vector<ImageVector> imageVectors;
//    VectorXf mean_vector(length);
//    
//    // Create image vectors
//    for(int i=0; i<image_num; i++){
//        ImageVector imageVector = ImageVector(width, height);
//        imageVectors.push_back(imageVector);
//    }
//    
//    // Put gray pixel value into the image vector and calculate the mean image
//    for(int i=0; i<width; i++){
//        for(int j=0; j<height; j++){
//            uchar color_temp1 = img1.at<uchar>(Point(j, i));
//            imageVectors[0].value[j*width+i] = color_temp1;
//            uchar color_temp2 = img2.at<uchar>(Point(j, i));
//            imageVectors[1].value[j*width+i] = color_temp2;
//            uchar color_temp3 = img3.at<uchar>(Point(j, i));
//            imageVectors[2].value[j*width+i] = color_temp3;
//            
//            float mean = 0;
//            for(int k=0; k<image_num; k++){
//                mean += imageVectors[k].value[j*width+i];
//            }
//            mean /= 3.0f;
//            img_mean.at<uchar>(Point(j, i)) = mean;
//            mean_vector(j*width+i) = mean;
//        }
//    }
//    
//    vector<VectorXf> vectors;
//    
//    for(int j=0; j<3; j++){
//        VectorXf i_vector(length);
//        for(int i=0; i<length; i++){
//            i_vector(i) = imageVectors[j].value[i];
//        }
//        vectors.push_back(i_vector);
//    }
//    
//    MatrixXf sub_M(image_num, length);
//    for(int i=0; i<image_num; i++){
//        for(int j=0; j<length; j++){
//            sub_M(i, j) = vectors[i][j] - mean_vector(j);
//        }
//    }
//    
////    VectorXi I(length);
////    for(int i=0; i<length; i++){
////        
////    }
////    
////    
//    MatrixXf cov_M(length, length);
//    MatrixXf cov_T(length, length);
//    for(int i=0; i<length; i++){
//        VectorXf temp(image_num);
//        for(int k=0; k<image_num; k++){
//            temp(k) = sub_M(k, i);
//        }
//        cov_M(i, i) = Cov(temp, temp);
//    }
//    
////    for(int i=1; i<length-1; i++){
////        for(int j=i+1; j<length; j++){
////            float total = 0;
////            for(int k=0; k<image_num; k++){
////                total += sub_M(k, i)*sub_M(k, j);
////            }
////            cov_T(i, j) = total/float(image_num-1);
////        }
////    }
//    
////    cout<<cov_M(100, 100)<<endl;
////    cout<<cov_T(100, 1000)<<endl;
//    
//    
//    
//    // Test
//    vector<double> X;
//    double x;
//    x=1; X.push_back(x);
//    x=1; X.push_back(x);
//    x=4; X.push_back(x);
////    double x_mean = 1.3333;
////    double result;
//    PCA_Demension pca_test;
////    pca_test.covariance(X, x_mean, X, x_mean, result);
////
////    cout<<result<<endl;
//    
//    x=2; X.push_back(x);
//    x=2; X.push_back(x);
//    x=3; X.push_back(x);
//    x=4; X.push_back(x);
//    x=9; X.push_back(x);
//    x=10; X.push_back(x);
//    x=7; X.push_back(x);
//    
//    vector<double> Y;
//    x=3; Y.push_back(x);
//    x=4; Y.push_back(x);
//    x=3; Y.push_back(x);
//    x=2; Y.push_back(x);
//    x=1; Y.push_back(x);
//    x=5; Y.push_back(x);
//    x=12; Y.push_back(x);
//    x=6; Y.push_back(x);
//    x=5; Y.push_back(x);
//    x=2; Y.push_back(x);
//    
//    vector<double> Z;
//    x=6; Z.push_back(x);
//    x=8; Z.push_back(x);
//    x=12; Z.push_back(x);
//    x=1; Z.push_back(x);
//    x=11; Z.push_back(x);
//    x=14; Z.push_back(x);
//    x=9; Z.push_back(x);
//    x=20; Z.push_back(x);
//    x=10; Z.push_back(x);
//    x=10; Z.push_back(x);
//    
//    vector<vector<double>> Feature_Data;
//    Feature_Data.push_back(X);
//    Feature_Data.push_back(Y);
//    Feature_Data.push_back(Z);
//    
//    vector<vector<double>> PCA_Features;
//    vector<vector<double>> Final_Data;
//    
//    pca_test.PCA_demension(Feature_Data, 5, PCA_Features, Final_Data);
//    
//    cout<<"test"<<endl;
//    
//    
////    m(1, 10) = 1;
////    cout<<m(1, 10)<<endl;
////    cout<<m(length-1, length-1)<<endl;
//    
//    
////    MatrixXd test(length, length);
////    cout<<test(10, 10)<<endl;
//
//    
//////    ImageVector test = ImageVector(2,2);
//////    test.value[0] = 8; test.value[1] = 9; test.value[2] = 11; test.value[3] = 12;
//////    int testMean = 0;
//////    for(int i=0; i<test.value.size(); i++){
//////        testMean += test.value[i];
//////    }
//////    testMean /= test.value.size();
//////    cout<<Cov(test, testMean, test, testMean )<<endl;
////    
////
/////* three channels */
//////    for(int i=0; i<width; i++){
//////        for(int j=0; j<height; j++){
//////            Vec3b color_temp1 = img1.at<Vec3b>(Point(i, j)); //+img2.at<Vec3b>(Point(i, j))+img3.at<Vec3b>(Point(i, j)));
//////            Vec3b color_temp2 = img2.at<Vec3b>(Point(i, j));
//////            Vec3b color_temp3 = img3.at<Vec3b>(Point(i, j));
//////            
//////            cout<<color_temp1<<endl;
////////            vec3 color_mean = vec3(color_temp1[0]+color_temp2[0]+color_temp3[0],
////////                                   color_temp1[1]+color_temp2[1]+color_temp3[1],
////////                                   color_temp1[2]+color_temp2[2]+color_temp3[2]);
////////            color_mean /= 3.0f;
////////            
////////            Vec3b color_new = Vec3b(color_mean[0], color_mean[1], color_mean[2]);
//////////            cout<<color_new<<endl;
//////////            cout<<color_mean[0]<<" "<<color_mean[1]<<" "<<color_mean[2]<<endl;
////////            img_mean.at<Vec3b>(j, i) = color_new;
//////        }
//////    }
//////    Scalar intensity;
//////    intensity = img1.at<Vec3b>(Point(20, 20));
//////    cout<<intensity<<endl;
//////    Vec3b color_temp;
//////    color_temp = img1.at<Vec3b>(Point(20, 20));
//////    color_temp[0] += 80;
//////    color_temp[1] += 80;
//////    color_temp[2] += 80;
//////    
//////    Vec3f red = color_temp;
//////    cout<<"red: "<<red[0]<<endl;
//////    
//////    img1.at<Vec3b>(20, 20) = color_temp;
//////    intensity = img1.at<Vec3b>(Point(20, 20));
//////    cout<<intensity<<endl;
//    
//    
//    imshow("scrouce image", img1);
//    imshow("mean image", img_mean);
//    
//    waitKey(0);
//    return 0;
//}
//
//// Calculate convariance
//float Cov(VectorXf X, VectorXf Y){
//    float total = 0.0f;
//    for(int i=0; i<X.size(); i++){
//        total += X[i]* Y[i];
//    }
//    return total/float(X.size()-1);
//}
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
