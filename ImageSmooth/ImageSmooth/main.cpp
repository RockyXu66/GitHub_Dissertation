//
//  main.cpp
//  ImageSmooth
//
//  Created by Yinghan Xu on 05/08/2017.
//  Copyright © 2017 YinghanXu. All rights reserved.
//

#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "opencv2/imgproc/imgproc.hpp"

using namespace std;
using namespace cv;

int main(int argc, const char * argv[]) {
    int cell_dimension = 30;
    int num_components = 90;
    Mat img = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/ResultImages/head(720)_cell"+to_string(cell_dimension)+"_"+to_string(num_components)+"/head11.png");
    
    Mat smoothedImage = img.clone();
    Mat resultImage = img.clone();
//    for(int i=0; i<img.rows; i++){
//        for(int j=0; j<img.cols; j++){
//            smoothedImage.at<Vec3b>(i, j)[0] = 0;
//            smoothedImage.at<Vec3b>(i, j)[1] = 0;
//            smoothedImage.at<Vec3b>(i, j)[2] = 0;
//        }
//    }
    
//    for(int i=0; i<img.rows/30; i++){
//        for(int j=0; j<img.cols/30; j++){
//            resultImage.col(30*j+1) = smoothedImage.col(30*j+1).clone();
//            resultImage.col(30*j+2) = smoothedImage.col(30*j+2).clone();
//            resultImage.col(30*j+3) = smoothedImage.col(30*j+3).clone();
//            resultImage.col(30*j+4) = smoothedImage.col(30*j+4).clone();
//            resultImage.col(30*j+5) = smoothedImage.col(30*j+5).clone();
//        }
//    }
    int kernelSize = 5;
    medianBlur (smoothedImage, smoothedImage, kernelSize);
//    blur(smoothedImage, smoothedImage, Size(kernelSize, kernelSize));
    
    for(int i=1; i<img.rows-1; i++){
        for(int j=1; j<img.cols-1; j++){
            if(j%cell_dimension==0){
                resultImage.at<Vec3b>(i, j-1) = smoothedImage.at<Vec3b>(i, j-1);
                resultImage.at<Vec3b>(i, j) = smoothedImage.at<Vec3b>(i, j);
                resultImage.at<Vec3b>(i, j+1) = smoothedImage.at<Vec3b>(i, j+1);
                resultImage.at<Vec3b>(i, j+2) = smoothedImage.at<Vec3b>(i, j+2);
            }
            if(i%cell_dimension==0){
                resultImage.at<Vec3b>(i-1, j) = smoothedImage.at<Vec3b>(i-1, j);
                resultImage.at<Vec3b>(i, j) = smoothedImage.at<Vec3b>(i, j);
                resultImage.at<Vec3b>(i+1, j) = smoothedImage.at<Vec3b>(i+1, j);
                resultImage.at<Vec3b>(i+2, j) = smoothedImage.at<Vec3b>(i+2, j);
            }
        }
    }
    
    imwrite("medianBlur_cell_"+to_string(cell_dimension)+"_"+to_string(num_components)+"_"+to_string(kernelSize)+".png", resultImage);
    
//    imshow("OriginalImage", img);
//    imshow("SmoothedImage", smoothedImage);
//    imshow("ResultImage", resultImage);
//    waitKey(0);
    return 0;
}
