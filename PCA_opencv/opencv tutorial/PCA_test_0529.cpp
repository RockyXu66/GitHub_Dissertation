////
////  PCA_test_0529.cpp
////  PCA_opencv
////
////  Created by Yinghan Xu on 29/05/2017.
////  Copyright © 2017 Yinghan Xu. All rights reserved.
////
//
//#include<iostream>
//#include<fstream>
//#include"Eigen/Eigen"
//#include<opencv2/opencv.hpp>
//using namespace Eigen;
//using namespace std;
//using namespace cv;
//void BubbleSort(float *pData, int Count);//冒泡排序
//const int num = 10;    // 样本数量
//const int dim = 400;    // 样本的维度（一张图片的维度）
//const int pca_dim = 100;// pca dimension of each sample
//
//int main()
//{
//    // 原始数据：每一行代表一个样本（这里为了测试时随机生成的）
//    MatrixXf data =MatrixXf::Random(num, dim);
////    MatrixXf data(10, 2);
////    data<<2.5, 2.4, 0.5, 0.7, 2.2, 2.9, 1.9, 2.2, 3.1, 3.0, 2.3, 2.7, 2.0, 1.6, 1.0, 1.1, 1.5, 1.6, 1.1, 0.9;
//    
//    
//    // 求每一维度上的平均值，即每一列的平均值，最后求出来是一个行向量
//    MatrixXf mean =data.colwise().mean();
//    cout<<"------------------Orignal data-------------------------\n"<<data<<endl;
//    cout<<"----------------Mean------------------\n"<<mean<<endl;
//    
//    MatrixXf C(num, dim);  // data-mean
//    MatrixXf C_T(dim, num);// C的转置
//    for(int i=0; i<num; i++)
//    {
//        C.row(i) =data.row(i)-mean;  //C= data-mean
//    }
//    C_T =C.transpose();
////    cout<<"----------------Subtract matrix------------------\n"<<C<<endl;
////    cout<<"----------------Transpose of C------------------\n"<<C_T<<endl;
//    
//    /*
//     协方差矩阵A的维度为dim*dim，这里除不除dim最后结果都一样，
//     在书中的定义中是要除的，估计是有数学含义的。
//     */
//    MatrixXf A =C_T*C /float(num-1);
//    EigenSolver<MatrixXf> sol_A;
//    sol_A.compute(A);
//    //求协方差矩阵的特征值和特征向量
//    MatrixXf eigenvals_A =sol_A.eigenvalues().real();
//    MatrixXf eigenvecs_A =sol_A.eigenvectors().real();
//    
//    cout<<"----------------EigenValues------------------\n"<<eigenvals_A<<endl;
//    cout<<"----------------EigenVectors------------------\n"<<eigenvecs_A<<endl;
//    
//    //排序
//    float a[dim]={0.0};
//    for(int i=0; i<dim; i++)
//    {
//        a[i] =eigenvals_A(i, 0);
//    }
//    BubbleSort(a, dim);
//    MatrixXf sorted_eigenvals_A(dim, 1);
//    MatrixXf sorted_eigenvecs_A(dim, dim);
//    //对特征值从大到小排序
//    for(int i=0; i<dim; i++)
//    {
//        sorted_eigenvals_A(i, 0) =a[i];
//    }
//    //对应的改变特征向量的排序
//    for(int i=0; i<dim; i++)
//    {
//        for(int j=0; j<dim; j++)
//        {
//            if(eigenvals_A(i, 0) == a[j])
//            {
//                sorted_eigenvecs_A.col(j) =eigenvecs_A.col(i);
//            }
//        }
//    }
//
//    //求PCA的投影矩阵及降维结果
//    MatrixXf project_mat(dim, pca_dim); //投影矩阵
//    MatrixXf res(num, pca_dim);         //最终输出，维度为num*pca_dim
//    for(int i=0; i<pca_dim; i++)        //抽取主成分
//    {
//        project_mat.col(i) =sorted_eigenvecs_A.col(i);
//    }
//    
//    /*
//     res就是原始数据降维后的矩阵
//     特别注意，这里是C*project_mat,即(data-mean)*project_mat
//     在网上看到很多人这里都是直接用的data*project_mat,那样是不对的
//     */
//    res =C* project_mat;
//    
//    MatrixXf RowDataAdjust = C.transpose();
//    MatrixXf RowFeatureVector = project_mat.transpose();
//    cout<<"test\n"<<RowDataAdjust<<endl;
//    cout<<"test\n"<<RowFeatureVector<<endl;
//    
//    MatrixXf FinalData = RowFeatureVector * RowDataAdjust;
//    cout<<"FinalData: \n"<<FinalData<<endl;
//    
//    MatrixXf RowOriginalData(num, dim);
//    RowOriginalData = RowFeatureVector.transpose() * FinalData;
//    cout<<"RowOriginalData:  \n"<<RowOriginalData<<endl;
//    
////    //-----------------------------OpenCV中的PCA--------------------------------//
////    CvMat *pca_input         =cvCreateMat( num, dim, CV_32FC1);               // 输入
////    CvMat *pca_avg           =cvCreateMat( 1, dim, CV_32FC1);		          // 平均值
////    CvMat *pca_eigenvalue    =cvCreateMat( 1, std::min(num, dim), CV_32FC1);  // 特征值
////    CvMat *pca_eigenvector   =cvCreateMat( std::min(num, dim), dim, CV_32FC1);// 特征向量
////    CvMat *pca_eigenvector_T =cvCreateMat(dim, std::min(num, dim), CV_32FC1); // 特征向量
////    CvMat *pca_output        =cvCreateMat(num, pca_dim, CV_32FC1);            // 输出
////    // 数据导入
////    for(int i=0; i<num; i++)
////    {
////        for(int j=0; j<dim; j++)
////        {
////            cvSetReal2D(pca_input, i, j, data(i, j));
////        }
////    }
////    cvCalcPCA(pca_input, pca_avg, pca_eigenvalue, pca_eigenvector, CV_PCA_DATA_AS_ROW);
////    
////    
//    cout<<"\n---------------Eeigen Tool降维结果-----------------"<<endl;
//    cout<<res<<endl;
////
////    cout<<"---------------OpenCV Tool降维结果---------------"<<endl;
////    cvProjectPCA(pca_input, pca_avg, pca_eigenvector, pca_output);
////    // 数据输出
////    for(int i=0; i<num; i++)
////    {
////        for(int j=0; j<pca_dim; j++)
////        {
////            float temp =0.0;
////            temp =cvGetReal2D(pca_output, i, j);
////            cout<<temp<<"   ";
////        }
////        cout<<endl;
////    }
////    
////    system("pause");
//    return 0;
//}
////冒泡排序算法
//void BubbleSort(float *pData, int Count)  
//{   
//    float iTemp;  
//    for(int i=1; i<Count; i++)   
//    { 
//        //一共进行（count-1）轮，每次得到一个最大值   
//        for(int j=Count-1; j>=i; j--) //每次从最后往前交换，得到最大值  
//        {     
//            if(pData[j]>pData[j-1])    
//            {     
//                iTemp      = pData[j-1];    
//                pData[j-1] = pData[j];    
//                pData[j]   = iTemp;    
//            }   
//        }  
//    }  
//}
