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

void savePCA(const string &file_name,PCA pca_)
{
    FileStorage fs(file_name,FileStorage::WRITE);
    fs << "mean" << pca_.mean;
    fs << "e_vectors" << pca_.eigenvectors;
    fs << "e_values" << pca_.eigenvalues;
    fs.release();
}

void savePCA(const string &file_name,vector<PCA> pcas, int cell_num){
    FileStorage fs(file_name,FileStorage::WRITE);
    Mat means(cell_num, pcas[0].mean.cols, CV_32F);
    Mat eigenValues(pcas[0].eigenvalues.rows, cell_num, CV_32F);
    Mat eigenVectors(pcas[0].eigenvectors.rows*cell_num, pcas[0].eigenvectors.cols, CV_32F);
    for(int i=0; i<cell_num; i++){
        pcas[i].mean.copyTo(means(Rect(0, i, pcas[0].mean.cols, 1)));
        pcas[i].eigenvalues.copyTo(eigenValues(Rect(i, 0, 1, pcas[0].eigenvalues.rows)));
        pcas[i].eigenvectors.copyTo(eigenVectors(Rect(0, i*pcas[0].eigenvectors.rows, pcas[0].eigenvectors.cols, pcas[0].eigenvectors.rows)));
    }
    fs << "mean" << means;
    fs << "e_values" << eigenValues;
    fs << "e_vectors" << eigenVectors;
    fs.release();
}

PCA loadPCA(const string &file_name, int num_components)
{
    PCA pca_;
    FileStorage fs(file_name,FileStorage::READ);
    fs["mean"] >> pca_.mean;
    cout<<"mean ";
    fs["e_vectors"] >> pca_.eigenvectors;
    cout<<"e_vectors ";
    fs["e_values"] >> pca_.eigenvalues;
    cout<<"e_values"<<endl;
//    pca_.eigenvectors = pca_.eigenvectors(Rect(0,0,960*960,num_components));
//    pca_.eigenvalues = pca_.eigenvalues(Rect(0,0,1,num_components));
    
    fs.release();
    return pca_;
}



void saveScores(const string &file_name, Mat scores){
    FileStorage fs(file_name, FileStorage::WRITE);
    fs << "scores" << scores;
    fs.release();
}

Mat loadScores(const string &file_name, int num_components){
    Mat scores;
    FileStorage fs(file_name,FileStorage::READ);
    fs["scores"] >> scores;
    scores = scores(Rect(0,0,num_components,1));
    fs.release();
    return scores;
}

Mat mergeChannels(Mat b, Mat g, Mat r, int dimension, string oneImagePath){
    Mat channels[3];
    channels[0] = norm_0_255(b).reshape(1, dimension).clone();
    channels[1] = norm_0_255(g).reshape(1, dimension).clone();
    channels[2] = norm_0_255(r).reshape(1, dimension).clone();
    Mat result = imread(oneImagePath);
    merge(channels, 3, result);
    return result;
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

// Number of components to keep for the PCA:
const int num_components = 90;
const int imageIndex = 10;
const int image_num = 120;
const int image_width = 960;
const int image_height = 960;
const int cell_dimension = 20;

String oneImagePath = "/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/artifix_120/artifix1.png";

int main(int argc, const char *argv[]) {
    
    // Holds some images:
    vector<Mat> db_b;
    vector<Mat> db_g;
    vector<Mat> db_r;
    
    int cell_num = (image_width*image_height)/(cell_dimension*cell_dimension);
    
    Mat grayImg_b = imread(oneImagePath, IMREAD_GRAYSCALE);
    
    Mat c = imread(oneImagePath, IMREAD_GRAYSCALE);
    c = c.reshape(1, 10);
    Mat d = c.t(); // Transpose
    cout<<c.rows<<" "<<c.cols<<endl;
    cout<<d.rows<<" "<<d.cols<<endl;
//    for(int i=0; i<c.rows; i++){
//        cout<<double(c.at<uchar>(i, 0))<<endl;
//    }

//    cout<<"Start loading images:"<<endl;
//    for(int i=0; i<image_num; i++){
//        Mat src = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/artifix_120/artifix" + to_string(i+1)+".png");
//        Mat bgr[3];   //destination array
//        split(src,bgr);//split source
//        db_b.push_back(bgr[0]);
//        db_g.push_back(bgr[1]);
//        db_r.push_back(bgr[2]);
//    }
//    cout<<"Finish loading images."<<endl;
    
    /*
    // Crop to cells
    cout<<"Start loading images:"<<endl;
    auto t_load1 = chrono::high_resolution_clock::now();
    vector<Mat> cells_b, cells_g, cells_r;
    for(int i=0; i<image_height/cell_dimension; i++){
        for(int j=0; j<image_width/cell_dimension; j++){
            vector<Mat> cell_b, cell_g, cell_r;
            for(int n=0; n<image_num; n++){
                Mat img = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/artifix_120/artifix" + to_string(n+1)+".png");
                img = img(Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension));
                Mat bgr[3];   //destination array
                split(img,bgr);//split source

                cell_b.push_back(bgr[0]);
                cell_g.push_back(bgr[1]);
                cell_r.push_back(bgr[2]);
            }
            // Build a matrix with the observations in row:
            Mat data_b = asRowMatrix(cell_b, CV_32FC1);
            Mat data_g = asRowMatrix(cell_g, CV_32FC1);
            Mat data_r = asRowMatrix(cell_r, CV_32FC1);
            cells_b.push_back(data_b);
            cells_g.push_back(data_g);
            cells_r.push_back(data_r);
        }
        cout<<i+1<<"/"<<(image_height/cell_dimension)+1<<" ";
    }
    cout<<endl;
    auto t_load2 = chrono::high_resolution_clock::now();
    chrono::duration<double> t_load_elapsed = t_load2 - t_load1;
    cout<<"Finish loading images. Duration time:"<<float(t_load_elapsed.count()) / 60.0f << " min"<<endl;
    
    // Calculate PCA:
    cout<<"===Start calculating PCA==="<<endl;
    auto t_pca1 = chrono::high_resolution_clock::now();
    vector<PCA> pcas_b, pcas_g, pcas_r;
    for(int i=0; i<cell_num; i++){
        Mat data_b = cells_b[i];
        Mat data_g = cells_g[i];
        Mat data_r = cells_r[i];
    
        PCA pca_b(data_b, Mat(), CV_PCA_DATA_AS_ROW);
        PCA pca_g(data_g, Mat(), CV_PCA_DATA_AS_ROW);
        PCA pca_r(data_r, Mat(), CV_PCA_DATA_AS_ROW);
        pcas_b.push_back(pca_b);
        pcas_g.push_back(pca_g);
        pcas_r.push_back(pca_r);
        cout<<i+1<<"/"<<cell_num+1<<" ";
    }
    cout<<endl;
    auto t_pca2 = chrono::high_resolution_clock::now();
    chrono::duration<double> t_pca_elapsed = t_pca2 - t_pca1;
    cout << "Finish calculating PCA. Duration time: "<<float(t_pca_elapsed.count()) / 60.0f << " min"<< endl;

    // Save the PCA result
    cout<<"===Save to files==="<<endl;
    savePCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/PCA_cells_b(chest).txt", pcas_b, cell_num);
    savePCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/PCA_cells_g(chest).txt", pcas_g, cell_num);
    savePCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/PCA_cells_r(chest).txt", pcas_r, cell_num);
    cout<<"Finish saving to files."<<endl;
    
    cout<<"===Save scores to files==="<<endl;
    int pca_dim = pcas_b[0].eigenvectors.rows;
    Mat scores_cells_b(cell_num*image_num, pca_dim, CV_32F);
    Mat scores_cells_g(cell_num*image_num, pca_dim, CV_32F);
    Mat scores_cells_r(cell_num*image_num, pca_dim, CV_32F);
    for(int i=0; i<image_num; i++){
        for(int j=0; j<cell_num; j++){
            Mat score_b = pcas_b[j].project(cells_b[j].row(i));
            Mat score_g = pcas_g[j].project(cells_g[j].row(i));
            Mat score_r = pcas_r[j].project(cells_r[j].row(i));
            score_b.copyTo(scores_cells_b(Rect(0,j+i*cell_num,pca_dim,1)));
            score_g.copyTo(scores_cells_g(Rect(0,j+i*cell_num,pca_dim,1)));
            score_r.copyTo(scores_cells_r(Rect(0,j+i*cell_num,pca_dim,1)));
        }
    }
    saveScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/Scores_cells_b(chest).txt", scores_cells_b);
    saveScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/Scores_cells_g(chest).txt", scores_cells_g);
    saveScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/Scores_cells_r(chest).txt", scores_cells_r);
    cout<<"Finish saving scores to files."<<endl;
     */
    
//    PCA pca_b_new;
//    FileStorage fs("ResultPCA/PCA_b_Result(chest).txt",FileStorage::READ);
//    FileNode fn(fs["PCA"]);
//    pca_b_new.read(fn);
    
    // Load the PCA result
    cout<<"===Load PCA results==="<<endl;
    auto load1 = chrono::high_resolution_clock::now();
    PCA pca_b = loadPCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/PCA_cells_b(chest).txt", num_components);
    PCA pca_g = loadPCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/PCA_cells_g(chest).txt", num_components);
    PCA pca_r = loadPCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/PCA_cells_r(chest).txt", num_components);
//    Mat scores_b = loadScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/Scores_b_Result(chest).txt", num_components);
//    Mat scores_g = loadScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/Scores_g_Result(chest).txt", num_components);
//    Mat scores_r = loadScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/Scores_r_Result(chest).txt", num_components);
    auto load2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_load = load2 - load1;
    cout<<"Finish loading PCA results. Duration time: "<<float(elapsed_load.count()) / 60.0f << " min" <<endl;
    
//    // The mean image:
//    Mat mean = mergeChannels(pca_b.mean.row(0), pca_g.mean.row(0), pca_r.mean.row(0), dimension, oneImagePath);
//    imshow("mean", mean);
//
//    // The first three eigen texture:
//    Mat pca1 = mergeChannels(pca_b.eigenvectors.row(0), pca_g.eigenvectors.row(0), pca_r.eigenvectors.row(0), dimension, oneImagePath);
//    imshow("pc1", pca1);
    
//    Mat scores_b = pca_b.project(data_b.row(imageIndex));
//    saveScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/Scores_b_Result(chest).txt", dst_b);
//    cout<<dst_b.rows<<" "<<dst_b.cols<<endl;
//    Mat scores_g = pca_g.project(data_g.row(imageIndex));
//    saveScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/Scores_g_Result(chest).txt", dst_g);
//    Mat scores_r = pca_r.project(data_r.row(imageIndex));
//    saveScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA/Scores_r_Result(chest).txt", dst_r);
    
    
//    Mat reconstruction_b = norm_0_255(pca_b.backProject(scores_b)).reshape(1, dimension);
//    Mat reconstruction_g = norm_0_255(pca_g.backProject(scores_g)).reshape(1, dimension);
//    Mat reconstruction_r = norm_0_255(pca_r.backProject(scores_r)).reshape(1, dimension);
//    
//    // Put three channel into one bgr image
//    Mat resultImage = mergeChannels(reconstruction_b, reconstruction_g, reconstruction_r, oneImagePath);
//    imwrite("ResultImages/artifix"+to_string(imageIndex+1)+"_"+to_string(num_components)+".png", resultImage);
//    cout<<"Finish reconstruction"<<endl;
//    
//    imshow("Reconstruction", resultImage);
//    
//    // Crop to cell (Rect(x, y, width, height))
//    Mat cell = resultImage(Rect(100,0,150,150));
//    imshow("cell", cell);
    
    
    // Show the images:
//    waitKey(0);
    
    // Success!
    return 0;
}





































