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

void save(const string &file_name,PCA pca_){
    FileStorage fs(file_name,FileStorage::WRITE);
    fs << "mean" << pca_.mean;
    cout<<"mean ";
    fs << "e_vectors" << pca_.eigenvectors;
    cout<<"e_vectors ";
    fs << "e_values" << pca_.eigenvalues;
    cout<<"e_values"<<endl;
    fs.release();
}

// Put three channel into one bgr image
Mat mergeChannels(Mat b, Mat g, Mat r, string oneImagePath){
    Mat channels[3];
    channels[0] = b;
    channels[1] = g;
    channels[2] = r;
    Mat result = imread(oneImagePath);
    merge(channels, 3, result);
    return result;
}

void saveScores(const string &file_name, Mat scores){
    FileStorage fs(file_name,FileStorage::WRITE);
    fs << "scores" << scores;
    cout<<"one scores"<<endl;
    fs.release();
}

PCA load(const string &file_name, int num_components){
    PCA pca_;
    FileStorage fs(file_name,FileStorage::READ);
    fs["mean"] >> pca_.mean;
    cout<<"mean ";
    fs["e_vectors"] >> pca_.eigenvectors;
    cout<<"e_vectors ";
    fs["e_values"] >> pca_.eigenvalues;
    cout<<"e_values"<<endl;
    pca_.eigenvectors = pca_.eigenvectors(Rect(0,0,720*720,num_components));
    pca_.eigenvalues = pca_.eigenvalues(Rect(0,0,1,num_components));
    
    fs.release();
    return pca_;
}

PCA loadPCA(const string &file_name){
    PCA pca_;
    FileStorage fs(file_name,FileStorage::READ);
    fs["mean"] >> pca_.mean;
    cout<<"mean ";
    fs["e_vectors"] >> pca_.eigenvectors;
    cout<<"e_vectors ";
    fs["e_values"] >> pca_.eigenvalues;
    cout<<"e_values"<<endl;
    if(pca_.mean.empty()){
        cout<<"Cannot load file "<<file_name<<endl;
    }
    fs.release();
    return pca_;
}

Mat loadScores(const string &file_name){
    Mat scores;
    FileStorage fs(file_name,FileStorage::READ);
    fs["scores"] >> scores;
    if(scores.empty()){
        cout<<"Cannot load file "<<file_name<<endl;
    }
    fs.release();
    cout<<"Finish loading one scores file"<<endl;
    return scores;
}

// Number of components to keep for the PCA:
const int num_components = 10;

const int imageIndex = 10;

const int image_num = 899;

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
    
    String oneImagePath = "/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution720)/head2.png";
    Mat oneImage = imread(oneImagePath);

//    cout<<"Start loading images:"<<endl;
//    for(int i=1; i<image_num+1; i++){
//        Mat src = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution720)/head" + to_string(i+1)+".png");
//        Mat bgr[3];   //destination array
//        split(src,bgr);//split source
//        db_b.push_back(bgr[0]);
//        db_g.push_back(bgr[1]);
//        db_r.push_back(bgr[2]);
//    }
//    cout<<"Finish loading images."<<endl;
//    
//    // Build a matrix with the observations in row:
//    Mat data_b = asRowMatrix(db_b, CV_32FC1);
//    Mat data_g = asRowMatrix(db_g, CV_32FC1);
//    Mat data_r = asRowMatrix(db_r, CV_32FC1);
//    
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
//
//    // Save the PCA result
//    cout<<"===Save PCA results to files==="<<endl;
//    save("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_whole/PCA_b(head720).txt", pca_b);
//    save("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_whole/PCA_g(head720).txt", pca_g);
//    save("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_whole/PCA_r(head720).txt", pca_r);
////    FileStorage fs("ResultPCA/PCA_b_Result(chest).txt",FileStorage::WRITE);
////    pca_b.write(fs);
//    cout<<"Finish saving to files."<<endl;
//    
//    // Save scores
//    int pca_dim = pca_b.eigenvectors.rows;
//    Mat scores_b(image_num, pca_dim, CV_32F);
//    Mat scores_g(image_num, pca_dim, CV_32F);
//    Mat scores_r(image_num, pca_dim, CV_32F);
//    for(int i=0; i<image_num; i++){
//        Mat score_b = pca_b.project(data_b.row(i));
//        Mat score_g = pca_g.project(data_g.row(i));
//        Mat score_r = pca_r.project(data_r.row(i));
//        score_b.copyTo(scores_b(Rect(0,i,pca_dim,1)));
//        score_g.copyTo(scores_g(Rect(0,i,pca_dim,1)));
//        score_r.copyTo(scores_r(Rect(0,i,pca_dim,1)));
//        cout<<i+1<<" ";
//    }
//    cout<<endl;
//    cout<<"===Save scores to files==="<<endl;
////    Mat dst_b = pca_b.project(data_b.row(imageIndex));
//    saveScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_whole/Scores_b(head720).txt", scores_b);
////    Mat dst_g = pca_g.project(data_g.row(imageIndex));
//    saveScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_whole/Scores_g(head720).txt", scores_g);
////    Mat dst_r = pca_r.project(data_r.row(imageIndex));
//    saveScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_whole/Scores_r(head720).txt", scores_r);
//    cout<<"Finish saving scores to files."<<endl;
    
    
    // Load PCA results and scores
    cout<<"===Start loading PCA results and scores==="<<endl;
    PCA pca_b = load("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_whole/PCA_b(head720).txt", num_components);
    PCA pca_g = load("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_whole/PCA_g(head720).txt", num_components);
    PCA pca_r = load("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_whole/PCA_r(head720).txt", num_components);
    Mat scores_b = loadScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_whole/Scores_b(head720).txt");
    Mat scores_g = loadScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_whole/Scores_g(head720).txt");
    Mat scores_r = loadScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_whole/Scores_r(head720).txt");
    cout<<"Load successfully."<<endl;
    
    // Reconstruction
    cout<<"===Start reconstructing==="<<endl;
    auto re1 = chrono::high_resolution_clock::now();
    Mat scoreB = scores_b(Rect(0,imageIndex,num_components,1));
    Mat scoreG = scores_g(Rect(0,imageIndex,num_components,1));
    Mat scoreR = scores_r(Rect(0,imageIndex,num_components,1));
    Mat reconstruction_b = norm_0_255(pca_b.backProject(scoreB)).reshape(1, oneImage.rows);
    Mat reconstruction_g = norm_0_255(pca_g.backProject(scoreG)).reshape(1, oneImage.rows);
    Mat reconstruction_r = norm_0_255(pca_r.backProject(scoreR)).reshape(1, oneImage.rows);
    Mat resultImage = mergeChannels(reconstruction_b, reconstruction_g, reconstruction_r, oneImagePath);
    auto re2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_re = re2 - re1;
    cout<<"Reconstruction time: "<<elapsed_re.count()<<" s"<<endl;
    
    imwrite("ResultImages/head_"+to_string(num_components)+".png", resultImage);
    imshow("Reconstruction", resultImage);
    cout<<"Finish reconstruction."<<endl;
    
    
    
//    PCA pca_b_new;
//    FileStorage fs("ResultPCA/PCA_b_Result(chest).txt",FileStorage::READ);
//    FileNode fn(fs["PCA"]);
//    pca_b_new.read(fn);
    
////    // Load the PCA result
//    cout<<"===Load PCA results==="<<endl;
//    auto load1 = chrono::high_resolution_clock::now();
//////    PCA pca_b = load("ResultPCA/PCA_b_Result(soilder).txt", num_components);
//////    PCA pca_g = load("ResultPCA/PCA_g_Result(soilder).txt", num_components);
//////    PCA pca_r = load("ResultPCA/PCA_r_Result(soilder).txt", num_components);
//    cout<<"Blue ";
//    PCA pca_b = load("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/PCA_b_Result(chest).txt", num_components);
//    cout<<"Green ";
//    PCA pca_g = load("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/PCA_g_Result(chest).txt", num_components);
//    cout<<"Red"<<endl;
//    PCA pca_r = load("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/PCA_r_Result(chest).txt", num_components);
//    auto load2 = chrono::high_resolution_clock::now();
//    chrono::duration<double> elapsed_load = load2 - load1;
//    cout<<"Finish loading PCA results. Duration time: "<<float(elapsed_load.count()) / 60.0f << " min" <<endl;
    
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
    
//    // The mean image:
//    Mat m[3];
//    m[0] = norm_0_255(pca_b.mean.row(0)).reshape(1, db_b[0].rows).clone();
//    m[1] = norm_0_255(pca_g.mean.row(0)).reshape(1, db_b[0].rows).clone();
//    m[2] = norm_0_255(pca_r.mean.row(0)).reshape(1, db_b[0].rows).clone();
//    Mat mean = imread(oneImagePath);
////    Mat mean = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/training_set_head/fig_2.bmp");
//    merge(m, 3, mean);
//    imshow("mean", mean);
//
//    // The first three eigen texture:
//    Mat c1[3];
//    c1[0] = norm_0_255(pca_b.eigenvectors.row(0)).reshape(1, db_b[0].rows).clone();
//    c1[1] = norm_0_255(pca_g.eigenvectors.row(0)).reshape(1, db_b[0].rows).clone();
//    c1[2] = norm_0_255(pca_r.eigenvectors.row(0)).reshape(1, db_b[0].rows).clone();
//    Mat pca1 = imread(oneImagePath);
////    Mat pca1 = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/training_set_head/fig_2.bmp");
//    merge(c1, 3, pca1);
//    imshow("pc1", pca1);
////    imshow("pc2", norm_0_255(pca.eigenvectors.row(1)).reshape(1, db[0].rows));
////    imshow("pc3", norm_0_255(pca.eigenvectors.row(2)).reshape(1, db[0].rows));
//    
    
//    Mat reconstruction_b = norm_0_255(pca_b.backProject(dst_b)).reshape(1, db_b[0].rows);
    
//    Mat reconstruction_g = norm_0_255(pca_g.backProject(dst_g)).reshape(1, db_b[0].rows);
//    Mat reconstruction_r = norm_0_255(pca_r.backProject(dst_r)).reshape(1, db_b[0].rows);
//    
//    // Put three channel into one bgr image
//    Mat channel[3];
//    channel[0] = reconstruction_b.clone();
//    channel[1] = reconstruction_g.clone();
//    channel[2] = reconstruction_r.clone();
////    Mat resultImage= imread("textures/new_Model_screenShot/1.png");
//    Mat resultImage = imread(oneImagePath);
//    merge(channel, 3, resultImage);
//    imwrite("ResultImages/artifix"+to_string(imageIndex+1)+"_"+to_string(num_components)+".png", resultImage);
//    cout<<"Finish reconstruction"<<endl;
//    
//    imshow("Reconstruction", resultImage);
//    
//    // Crop to cell
//    Mat cell = resultImage(Rect(100,0,150,150));
//    imshow("cell", cell);
//    
//    
//    // Show the images:
//    waitKey(0);
    
    // Success!
    return 0;
}





































