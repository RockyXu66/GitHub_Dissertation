////
////  OpenCV_PCA_0529.cpp
////  PCA_opencv
////
////  Created by Yinghan Xu on 29/05/2017.
////  Copyright © 2017 Yinghan Xu. All rights reserved.
////
//
//#include <opencv2/opencv.hpp>
//#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/ml/ml.hpp>
//#include"Eigen/Eigen"
//
//using namespace Eigen;
//using namespace cv;
//using namespace std;
//
//void DoPca(const Mat &_data, int dim, Mat &eigenvalues, Mat &eigenvectors);
//
//void printMat( Mat _data )
//{
//    Mat data = cv::Mat_<double>(_data);
//    for ( int i=0; i<data.rows; i++ )
//    {
//        for ( int j=0; j< data.cols; j++ )
//        {
//            cout << data.at<double>(i,j) << "  ";
//        }
//        cout << endl;
//    }
//}
//
//int main(int argc, char* argv[])
//{
//    float A[ 20 ]={
//        2.5, 2.4, 0.5, 0.7, 2.2, 2.9, 1.9, 2.2, 3.1, 3.0, 2.3, 2.7, 2.0, 1.6, 1.0, 1.1, 1.5, 1.6, 1.1, 0.9 };
//    Mat DataMat = Mat::zeros( 10, 2, CV_32F );
//    
//    Mat img1 = imread("image_face/eye_1.png");
//    Mat img2 = imread("image_face/eye_2.png");
//    Mat img3 = imread("image_face/eye_3.png");
//    // Convert to grayscale
//        cvtColor(img1, img1, CV_BGR2GRAY);
//        cvtColor(img2, img2, CV_BGR2GRAY);
//        cvtColor(img3, img3, CV_BGR2GRAY);
//    
//    int num = 3;
//    int dim = 400;
//    int dim_pca = 100;
//    int width = 20;
//    int height = 20;
//    MatrixXf data(num, dim);
//    
//    //将数组A里的数据放入DataMat矩阵中
//    for ( int i=0; i<width; i++ )
//    {
//        for ( int j=0; j<height; j++ )
//        {
////            DataMat.at<float>(i, j) = A[i * 6 + j];
//            data(0, j*width+i) = float(img1.at<uchar>(Point(j, i)));
//                        data(1, j*width+i) = float(img2.at<uchar>(Point(j, i)));
//                        data(2, j*width+i) = float(img3.at<uchar>(Point(j, i)));
//        }
//    }
//    
////    for(int i=0; i<num; i++){
////        for(int j=0; j<dim; j++){
////            DataMat.at<float>(i, j) = data(i, j);
////        }
////    }
//    
//    for(int i=0; i<10; i++){
//        for(int j=0; j<2; j++){
//            DataMat.at<float>(i, j) = A[i*10+j];
//        }
//    }
//    
//    // OPENCV PCA
//    PCA pca(DataMat, noArray(), CV_PCA_DATA_AS_ROW);
//    
//    Mat eigenvalues;//特征值
//    Mat eigenvectors;//特征向量
//    
//    DoPca(DataMat, 2, eigenvalues, eigenvectors);
//    
//    cout << "eigenvalues:" << endl;
//    printMat( eigenvalues );
//    cout << "\n" << endl;
//    cout << "eigenvectors:" << endl;
//    printMat( eigenvectors );
//    
//    
//    
////    system("pause");
//    return 0;
//}
//void DoPca(const Mat &_data, int dim, Mat &eigenvalues, Mat &eigenvectors)
//{
//    assert( dim>0 );
//    Mat data =  cv::Mat_<double>(_data);
//    
//    int R = data.rows;
//    int C = data.cols;
//    
//    if ( dim>C )
//        dim = C;
//    
//    //计算均值
//    Mat m = Mat::zeros( 1, C, data.type() );
//    
//    for ( int j=0; j<C; j++ )
//    {
//        for ( int i=0; i<R; i++ )
//        {
//            m.at<double>(0,j) += data.at<double>(i,j);
//        }
//    }
//    
//    m = m/R;
//    //求取6列数据对应的均值存放在m矩阵中，均值： [1.67、2.01、1.67、2.01、1.67、2.01]
//    
//    
//    //计算协方差矩阵
//    Mat S =  Mat::zeros( R, C, data.type() );
//    for ( int i=0; i<R; i++ )
//    {
//        for ( int j=0; j<C; j++ )
//        {
//            S.at<double>(i,j) = data.at<double>(i,j) - m.at<double>(0,j); // 数据矩阵的值减去对应列的均值
//        }
//    }
//    Mat Average = S.t() * S /(R);
//    //计算协方差矩阵的方式----(S矩阵的转置 * S矩阵)/行数
//    
//    
//    //使用opencv提供的eigen函数求特征值以及特征向量
//    eigen(Average, eigenvalues, eigenvectors);
//}
