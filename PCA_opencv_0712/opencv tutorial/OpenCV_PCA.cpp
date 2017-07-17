#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>

using namespace cv;
using namespace std;

// Normalizes a given image into a value range between 0 and 255.
Mat norm_0_255(const Mat& src) {
    // Create and return normalized image:
    Mat dst;
    switch(src.channels()) {
        case 1:
            cv::normalize(src, dst, 0, 255, NORM_MINMAX, CV_8UC1);
            break;
        case 3:
            cv::normalize(src, dst, 0, 255, NORM_MINMAX, CV_8UC3);
            break;
        default:
            src.copyTo(dst);
            break;
    }
    return dst;
}

// Converts the images given in src into a row matrix.
Mat asRowMatrix(const vector<Mat>& data, int rtype, double alpha = 1, double beta = 0) {
    Mat dst(static_cast<int>(data.size()), data[0].rows*data[0].cols, CV_32F);
    for(unsigned int i = 0; i < data.size(); i++)
    {
        Mat image_row = data[i].clone().reshape(1,1);
        image_row.convertTo(dst.row(i),CV_32F);
    }
    return dst;
}

void save(const string &file_name,PCA pca_)
{
    FileStorage fs(file_name,FileStorage::WRITE);
    fs << "mean" << pca_.mean;
    fs << "e_vectors" << pca_.eigenvectors;
    fs << "e_values" << pca_.eigenvalues;
    fs.release();
}

PCA load(const string &file_name, int num_components)
{
    PCA pca_;
    FileStorage fs(file_name,FileStorage::READ);
    fs["mean"] >> pca_.mean;
    cout<<"mean ";
    fs["e_vectors"] >> pca_.eigenvectors;
    cout<<"e_vectors ";
    fs["e_values"] >> pca_.eigenvalues;
    cout<<"e_values"<<endl;
    pca_.eigenvectors = pca_.eigenvectors(Rect(0,0,960*960,num_components));
    pca_.eigenvalues = pca_.eigenvalues(Rect(0,0,1,num_components));
    
    fs.release();
    return pca_;
}

// Number of components to keep for the PCA:
const int num_components = 80;

const int imageIndex = 10;

const int image_num = 120;

int main(int argc, const char *argv[]) {
    
    /*
    for(int i=0; i<120; i++){
        Mat img;
        if(i<9){
            img = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/artifix_120x_png/artifix000"+to_string(i+1)+".png");
            imwrite("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/new/artifix"+to_string(i+1)+".png", img);
        }else if(i>=9&&i<99){
            img = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/artifix_120x_png/artifix00"+to_string(i+1)+".png");
            imwrite("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/new/artifix"+to_string(i+1)+".png", img);
        }else if(i>=99){
            img = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/artifix_120x_png/artifix0"+to_string(i+1)+".png");
            imwrite("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/new/artifix"+to_string(i+1)+".png", img);
        }
    }
     */
    
    
    // Holds some images:
    vector<Mat> db_b;
    vector<Mat> db_g;
    vector<Mat> db_r;
    
    String oneImagePath = "/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/artifix_120/artifix1.png";
//    Mat grayImg_b = imread("textures/new_Model_screenShot/1.png", IMREAD_GRAYSCALE);
//    Mat grayImg_b = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/training_set_head/fig_2.bmp", IMREAD_GRAYSCALE);
    Mat grayImg_b = imread(oneImagePath, IMREAD_GRAYSCALE);
    
//    Mat c = imread("textures/new_Model_screenShot/1.png", IMREAD_GRAYSCALE);
//    Mat c = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/training_set_head/fig_2.bmp", IMREAD_GRAYSCALE);
    Mat c = imread(oneImagePath, IMREAD_GRAYSCALE);
    c = c.reshape(1, 10);
    Mat d = c.t(); // Transpose
    cout<<c.rows<<" "<<c.cols<<endl;
    cout<<d.rows<<" "<<d.cols<<endl;
//    for(int i=0; i<c.rows; i++){
//        cout<<double(c.at<uchar>(i, 0))<<endl;
//    }

    Mat test;
    cout<<"Start loading images:"<<endl;
    for(int i=0; i<image_num; i++){
//        Mat img = imread("textures/new_Model_screenShot/" + to_string(i+1)+".png", IMREAD_GRAYSCALE);
//        test = imread("textures/new_Model_screenShot/" + to_string(i+1)+".png");
//        for(int x=0; x<grayImg_b.rows; x++){
//            for(int y=0; y<grayImg_b.cols; y++){
//                grayImg_b.at<uchar>(Point(x, y)) = img.at<uchar>(Point(x, y)); // channel b
//            }
//        }
//        Mat src = imread("textures/new_Model_screenShot/" + to_string(i+1)+".png");
//        Mat src = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/training_set_head/fig_" + to_string(i+1)+".bmp");
        Mat src = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/artifix_120/artifix" + to_string(i+1)+".png");
        Mat bgr[3];   //destination array
        split(src,bgr);//split source
        db_b.push_back(bgr[0]);
        db_g.push_back(bgr[1]);
        db_r.push_back(bgr[2]);
    }
    cout<<"Finish loading images."<<endl;
    
    // Build a matrix with the observations in row:
    Mat data_b = asRowMatrix(db_b, CV_32FC1);
    Mat data_g = asRowMatrix(db_g, CV_32FC1);
    Mat data_r = asRowMatrix(db_r, CV_32FC1);
    
//    // Perform a PCA:
//    cout<<"===Start calculating PCA==="<<endl;
//    auto start1 = chrono::high_resolution_clock::now();
//    PCA pca_b(data_b, Mat(), CV_PCA_DATA_AS_ROW, num_components);
//    auto finish1 = chrono::high_resolution_clock::now();
//    chrono::duration<double> elapsed1 = finish1 - start1;
//    cout << "Blue PCA time: "<< float(elapsed1.count()) / 60.0f << " min" << endl;
//    auto start2 = chrono::high_resolution_clock::now();
//    PCA pca_g(data_g, Mat(), CV_PCA_DATA_AS_ROW, num_components);
//    auto finish2 = chrono::high_resolution_clock::now();
//    chrono::duration<double> elapsed2 = finish2 - start2;
//    cout << "Green PCA time: "<< float(elapsed2.count()) / 60.0f << " min" << endl;
//    auto start3 = chrono::high_resolution_clock::now();
//    PCA pca_r(data_r, Mat(), CV_PCA_DATA_AS_ROW, num_components);
//    auto finish3 = chrono::high_resolution_clock::now();
//    chrono::duration<double> elapsed3 = finish3 - start3;
//    cout << "Red PCA time: "<< float(elapsed3.count()) / 60.0f << " min" << endl;

//    // Save the PCA result
//    cout<<"===Save to files==="<<endl;
////    save("ResultPCA/PCA_b_Result(soilder).txt", pca_b);
////    save("ResultPCA/PCA_g_Result(soilder).txt", pca_g);
////    save("ResultPCA/PCA_r_Result(soilder).txt", pca_r);
////    save("ResultPCA/PCA_b_Result(head).txt", pca_b);
////    save("ResultPCA/PCA_g_Result(head).txt", pca_g);
////    save("ResultPCA/PCA_r_Result(head).txt", pca_r);
//    save("ResultPCA/PCA_b_Result(chest).txt", pca_b);
//    save("ResultPCA/PCA_g_Result(chest).txt", pca_g);
//    save("ResultPCA/PCA_r_Result(chest).txt", pca_r);
////    FileStorage fs("ResultPCA/PCA_b_Result(chest).txt",FileStorage::WRITE);
////    pca_b.write(fs);
//    cout<<"Finish saving to files."<<endl;
    
//    PCA pca_b_new;
//    FileStorage fs("ResultPCA/PCA_b_Result(chest).txt",FileStorage::READ);
//    FileNode fn(fs["PCA"]);
//    pca_b_new.read(fn);
    
//    // Load the PCA result
    cout<<"===Load PCA results==="<<endl;
    auto load1 = chrono::high_resolution_clock::now();
////    PCA pca_b = load("ResultPCA/PCA_b_Result(soilder).txt", num_components);
////    PCA pca_g = load("ResultPCA/PCA_g_Result(soilder).txt", num_components);
////    PCA pca_r = load("ResultPCA/PCA_r_Result(soilder).txt", num_components);
    cout<<"Blue ";
    PCA pca_b = load("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/PCA_b_Result(chest).txt", num_components);
    cout<<"Green ";
    PCA pca_g = load("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/PCA_g_Result(chest).txt", num_components);
    cout<<"Red"<<endl;
    PCA pca_r = load("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/PCA_r_Result(chest).txt", num_components);
    auto load2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_load = load2 - load1;
    cout<<"Finish loading PCA results. Duration time: "<<float(elapsed_load.count()) / 60.0f << " min" <<endl;
    
//    PCA pca_b;
//    FileStorage fs("ResultPCA/PCA_b_Result(chest).txt",FileStorage::READ);
//    fs["mean"] >> pca_b.mean;
//    cout<<"mean ";
//    fs["e_vectors"] >> pca_b.eigenvectors;
//    cout<<"e_vectors ";
//    fs["e_values"] >> pca_b.eigenvalues;
//    cout<<"e_values"<<endl;
//    pca_.eigenvectors = pca_b.eigenvectors(Rect(0,0,960*960,num_components));
//    pca_.eigenvalues = pca_b.eigenvalues(Rect(0,0,1,num_components));
    
    
////    Mat mean = pca.mean.clone();
//    Mat eigenvalues_b = pca_b.eigenvalues.clone();
//    Mat eigenvectors_b = pca_b.eigenvectors.clone();
//    cout<<"Eigen Vectors: ";
//    cout<<eigenvectors_b.rows<<" "<<eigenvectors_b.cols<<endl;
////    cout<<eigenvectors_b<<endl;
//    cout<<"Eigen Values: ";
//    cout<<eigenvalues_b.rows<<" "<<eigenvalues_b.cols<<endl;
////    cout<<eigenvalues_b<<endl;
    
    // The mean image:
    Mat m[3];
    m[0] = norm_0_255(pca_b.mean.row(0)).reshape(1, db_b[0].rows).clone();
    m[1] = norm_0_255(pca_g.mean.row(0)).reshape(1, db_b[0].rows).clone();
    m[2] = norm_0_255(pca_r.mean.row(0)).reshape(1, db_b[0].rows).clone();
    Mat mean = imread(oneImagePath);
//    Mat mean = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/training_set_head/fig_2.bmp");
    merge(m, 3, mean);
    imshow("mean", mean);

    // The first three eigen texture:
    Mat c1[3];
    c1[0] = norm_0_255(pca_b.eigenvectors.row(0)).reshape(1, db_b[0].rows).clone();
    c1[1] = norm_0_255(pca_g.eigenvectors.row(0)).reshape(1, db_b[0].rows).clone();
    c1[2] = norm_0_255(pca_r.eigenvectors.row(0)).reshape(1, db_b[0].rows).clone();
    Mat pca1 = imread(oneImagePath);
//    Mat pca1 = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/training_set_head/fig_2.bmp");
    merge(c1, 3, pca1);
    imshow("pc1", pca1);
//    imshow("pc2", norm_0_255(pca.eigenvectors.row(1)).reshape(1, db[0].rows));
//    imshow("pc3", norm_0_255(pca.eigenvectors.row(2)).reshape(1, db[0].rows));
    
    Mat dst_b = pca_b.project(data_b.row(imageIndex));
    Mat reconstruction_b = norm_0_255(pca_b.backProject(dst_b)).reshape(1, db_b[0].rows);
    Mat dst_g = pca_g.project(data_g.row(imageIndex));
    Mat reconstruction_g = norm_0_255(pca_g.backProject(dst_g)).reshape(1, db_b[0].rows);
    Mat dst_r = pca_r.project(data_r.row(imageIndex));
    Mat reconstruction_r = norm_0_255(pca_r.backProject(dst_r)).reshape(1, db_b[0].rows);
    
    // Put three channel into one bgr image
    Mat channel[3];
    channel[0] = reconstruction_b.clone();
    channel[1] = reconstruction_g.clone();
    channel[2] = reconstruction_r.clone();
//    Mat resultImage= imread("textures/new_Model_screenShot/1.png");
    Mat resultImage = imread(oneImagePath);
    merge(channel, 3, resultImage);
    imwrite("ResultImages/artifix"+to_string(imageIndex+1)+"_"+to_string(num_components)+".png", resultImage);
    cout<<"Finish reconstruction"<<endl;
    
    imshow("Reconstruction", resultImage);
    
    // Crop to cell
    Mat cell = resultImage(Rect(100,0,150,150));
    imshow("cell", cell);
    
    
    // Show the images:
    waitKey(0);
    
    // Success!
    return 0;
}





































