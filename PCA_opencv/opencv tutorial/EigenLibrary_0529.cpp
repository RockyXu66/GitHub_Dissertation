//
//  EigenLibrary_0529.cpp
//  PCA_opencv
//
//  Created by Yinghan Xu on 29/05/2017.
//  Copyright © 2017 Yinghan Xu. All rights reserved.
//
#include <iostream>
#include"Eigen/Eigen"
#include<opencv2/opencv.hpp>
using namespace Eigen;
using namespace std;
using namespace cv;

const int num = 3;    // 样本数量
const int dim = 20*20;    // 样本的维度（一张图片的维度）
const int pca_dim = 2;//50;// pca dimension of each sample

//const int num = 10;    // 样本数量
//const int dim = 2;    // 样本的维度（一张图片的维度）
//const int pca_dim = 2;//50;// pca dimension of each sample

void BubbleSort(double *pData, int Count);

int main(){
    
    MatrixXd data(num, dim);
//    data<<2.5, 2.4, 0.5, 0.7, 2.2, 2.9, 1.9, 2.2, 3.1, 3.0, 2.3, 2.7, 2.0, 1.6, 1.0, 1.1, 1.5, 1.6, 1.1, 0.9;
    
    Mat img1 = imread("image_face/eye_1.png");
    Mat img2 = imread("image_face/eye_2.png");
    Mat img3 = imread("image_face/eye_3.png");

//    Mat img1 = imread("image_face/face_1.png");
//    Mat img2 = imread("image_face/face_2.png");
//    Mat img3 = imread("image_face/face_3.png");
    
    // Convert to grayscale
    cvtColor(img1, img1, CV_BGR2GRAY);
    cvtColor(img2, img2, CV_BGR2GRAY);
    cvtColor(img3, img3, CV_BGR2GRAY);
    
    int width =  img1.size().width;
    int height = img1.size().height;
    
    cout<<width*height<<endl;
    
    
    // Put gray pixel value into the data matrix
    for(int i=0; i<width; i++){
        for(int j=0; j<height; j++){
            data(0, j*width+i) = float(img1.at<uchar>(Point(j, i)));
            data(1, j*width+i) = float(img2.at<uchar>(Point(j, i)));
            data(2, j*width+i) = float(img3.at<uchar>(Point(j, i)));
        }
    }
    
    cout << "Here is the input data matrix:" << endl << data << endl << endl;
    
    MatrixXd mean =data.colwise().mean();

    Mat mean_image;
    mean_image = img1.clone();
    for(int i=0; i<width; i++){
        for(int j=0; j<height; j++){
            mean_image.at<uchar>(Point(j, i)) = mean(j*width+i);
        }
    }
    imshow("mean_image", mean_image);

    cout<<"----------------Mean------------------\n"<<mean<<endl;

    MatrixXd C(num, dim);  // data-mean
    MatrixXd C_T(dim, num);// C的转置
    
    for(int i=0; i<num; i++)
    {
        C.row(i) =data.row(i)-mean;  //C= data-mean
    }
    C_T =C.transpose();
    
    MatrixXd A =C_T*C;
    EigenSolver<MatrixXd> es(A);
    //求协方差矩阵的特征值和特征向量
    MatrixXd eigenvals_A =es.eigenvalues().real();
    MatrixXd eigenvecs_A =es.eigenvectors().real();
    cout << "The eigenvalues of A are:" << endl << eigenvals_A << endl;
    cout << "The matrix of eigenvectors, V, is:" << endl << eigenvecs_A << endl << endl;
    
    //排序
    double a[dim]={0.0};
    for(int i=0; i<dim; i++)
    {
        a[i] =eigenvals_A(i, 0);
    }
    BubbleSort(a, dim);
    MatrixXd sorted_eigenvals_A(dim, 1);
    MatrixXd sorted_eigenvecs_A(dim, dim);
    //对特征值从大到小排序
    for(int i=0; i<dim; i++)
    {
        sorted_eigenvals_A(i, 0) =a[i];
    }
    //对应的改变特征向量的排序
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            if(eigenvals_A(i, 0) == a[j])
            {
                sorted_eigenvecs_A.col(j) =eigenvecs_A.col(i);
            }
        }
    }
    
    cout << "The eigenvalues of A are:" << endl << sorted_eigenvals_A << endl;
    cout << "The matrix of eigenvectors, V, is:" << endl << sorted_eigenvecs_A << endl << endl;
    
    //求PCA的投影矩阵及降维结果
    MatrixXd extract_eigenvecs(dim, pca_dim); //投影矩阵
    MatrixXd projection(num, pca_dim);         //最终输出，维度为num*pca_dim
    for(int i=0; i<pca_dim; i++)        //抽取主成分
    {
        extract_eigenvecs.col(i) =sorted_eigenvecs_A.col(i);
    }
    
    /*
     res就是原始数据降维后的矩阵
     特别注意，这里是C*project_mat,即(data-mean)*project_mat
     在网上看到很多人这里都是直接用的data*project_mat,那样是不对的
     */
    projection = C * extract_eigenvecs;
    
    MatrixXd RowDataAdjust = C.transpose();
    MatrixXd RowFeatureVector = extract_eigenvecs.transpose();
//    cout<<"test\n"<<RowDataAdjust<<endl;
//    cout<<"test\n"<<RowFeatureVector<<endl;
    
    MatrixXd FinalData = RowFeatureVector * RowDataAdjust;
//    cout<<"FinalData: \n"<<FinalData<<endl;
    
    MatrixXd RowOriginalData(num, dim);
    RowOriginalData = RowFeatureVector.transpose() * FinalData;
    cout<<endl<<endl<<endl;
    
    cout<<"----------------Final projection data-----------------\n"<<FinalData<<endl;
    
    MatrixXd OriginalData = RowOriginalData.transpose();
    for(int i=0; i<num; i++){
        OriginalData.row(i) += mean;
    }
    
    cout<<endl<<endl;
    cout<<"Original Data: \n"<<OriginalData<<endl;
//    cout<<"OriginalData:  \n"<<RowOriginalData.transpose()<<endl;
    
//    109 110 107  100 95 90 88 87  85  84  82  80  82
//    110 107.667 101.333 94.6667
//    -1.88482    0.554253  3.51313  3.61696
    
    
//    complex<double> lambda = es.eigenvalues()[0];
//    cout << "Consider the first eigenvalue, lambda = " << lambda << endl;
//    VectorXcd v = es.eigenvectors().col(0);
//    cout << "If v is the corresponding eigenvector, then lambda * v = " << endl << lambda * v << endl;
//    cout << "... and A * v = " << endl << A.cast<complex<double> >() * v << endl << endl;
//    MatrixXcd D = es.eigenvalues().asDiagonal();
//    MatrixXcd V = es.eigenvectors();
//    cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;
}
//冒泡排序算法
void BubbleSort(double *pData, int Count)
{
    double iTemp;
    for(int i=1; i<Count; i++)
    {
        //一共进行（count-1）轮，每次得到一个最大值
        for(int j=Count-1; j>=i; j--) //每次从最后往前交换，得到最大值
        {
            if(pData[j]>pData[j-1])
            {
                iTemp      = pData[j-1];
                pData[j-1] = pData[j];
                pData[j]   = iTemp;
            }
        }
    }
}










