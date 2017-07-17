////
////  Stitch.cpp
////  opencv tutorial
////
////  Created by Yinghan Xu on 27/03/2017.
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
//
//#ifdef _DEBUG
//#pragma comment(lib, "opencv_core246d.lib")
//#pragma comment(lib, "opencv_imgproc246d.lib")   //MAT processing
//#pragma comment(lib, "opencv_highgui246d.lib")
//#pragma comment(lib, "opencv_stitching246d.lib");
//
//#else
//#pragma comment(lib, "opencv_core246.lib")
//#pragma comment(lib, "opencv_imgproc246.lib")
//#pragma comment(lib, "opencv_highgui246.lib")
//#pragma comment(lib, "opencv_stitching246.lib");
//#endif
//
//using namespace cv;
//using namespace std;
//
//void loadImages();
//void loadImagesTest();
//
//vector<Mat> vImg;
//
//int main()
//{
//    Mat panoramaImage;
//    
////    loadImages();
////    loadImagesTest();
//    
////    vImg.push_back(imread("images/livingRoom3/result_part2.jpg"));
////    vImg.push_back(imread("images/livingRoom3/IMG_4235.jpg"));
////    vImg.push_back(imread("images/livingRoom3/IMG_4242.jpg"));
//    vImg.push_back(imread("images/livingRoom4/result_part1.jpg"));
//    vImg.push_back(imread("images/livingRoom4/result_part2.jpg"));
//    
//    Stitcher stitcher = Stitcher::createDefault();
//    
//    
//    unsigned long AAtime=0, BBtime=0; //check processing time
//    AAtime = getTickCount(); //check processing time
//    
//    Stitcher::Status status = stitcher.stitch(vImg, panoramaImage);
//    
//    BBtime = getTickCount(); //check processing time
//    printf("%.2lf sec \n",  (BBtime - AAtime)/getTickFrequency() ); //check processing time
//    
////    resize(rImg, rImg, Size(), 0.075, 0.075);
//    
//    if (Stitcher::OK == status){
////        imshow("Stitching Result",panoramaImage);
//        String path = "images/resultImage/result"+to_string((BBtime - AAtime)/getTickFrequency())+".jpg";
//        imwrite(path, panoramaImage);
//    }
//    else{
//        printf("Stitching fail.");
//    }
//    
//    waitKey(0);
//    cvDestroyAllWindows();
//    
//}
//
//void loadImages(){
//    for(int i=8; i<15; i++){
//        String path = "images/livingRoom4/IMG_42"+to_string(i+57)+".jpg";
//        vImg.push_back(imread(path));
//    }
//}
//
//void loadImagesTest(){
//    for(int i=0; i<5; i++){
//        String path = "images/school2/IMG_410"+to_string(i+5)+".JPG";
//        vImg.push_back(imread(path));
//    }
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
