////
////  PCA_Get_ResultImage.cpp
////  opencv tutorial
////
////  Created by Yinghan Xu on 03/07/2017.
////  Copyright © 2017 Yinghan Xu. All rights reserved.
////
//
//#include <stdio.h>
//#include <iostream>
//#include <fstream>
//#include <string>
//
////#include "opencv2\opencv.hpp"
////#include "opencv2\stitching.hpp"
//#include "opencv2/core/core.hpp"
//#include "opencv2/stitching.hpp"
//#include "opencv2/opencv.hpp"
//#include "opencv2/imgcodecs.hpp"
//#include "opencv2/highgui.hpp"
//#include"Eigen/Eigen"
//#include <fstream>
//
//using namespace cv;
//using namespace std;
//using namespace Eigen;
//
//
//int main()
//{
//    Mat img;
//    
//    img = imread("PCAImages/1.png");
//    
////    imshow("test", img);
//    
////    Vec3b intensity = img.at<Vec3b>(1599, 2559); //at<Vec3b>(y, x)  => y is the height / x is the width
////    int blue = intensity.val[0];
////    int green = intensity.val[1];
////    int red = intensity.val[2];
////    
////    cout<<blue<<endl;
////    cout<<green<<endl;
////    cout<<red<<endl;
//    
//    MatrixXd newImage_r(1, 300*300);
//    MatrixXd newImage_g(1, 300*300);
//    MatrixXd newImage_b(1, 300*300);
//    
//    ifstream PCAImages_Red("PCAImages/PCAImages_Red.txt");
//    ifstream PCAImages_Green("PCAImages/PCAImages_Green.txt");
//    ifstream PCAImages_Blue("PCAImages/PCAImages_Blue.txt");
//    for(int i=0; i<300*300; i++){
////        newImage_r(0, i) = 255;
//        PCAImages_Red >> newImage_r(0, i);
//        PCAImages_Green >> newImage_g(0, i);
//        PCAImages_Blue >> newImage_b(0, i);
//    }
//    PCAImages_Red.close();
//    
//    Vec3b newPixel;
//    for(int i=0; i<300; i++){
//        for(int j=0; j<300; j++){
//            newPixel = Vec3b(newImage_b(0, j+i*300), newImage_g(0, j+i*300), newImage_r(0, j+i*300));
//            img.at<Vec3b>(i, j) = newPixel;
//        }
//    }
//    
//    imshow("newImage", img);
//    
//    imwrite("PCAImages/new2.png", img);
//    
//    waitKey(0);
//    cvDestroyAllWindows();
//    
//}
//
//
//
//
//
//
//
//
