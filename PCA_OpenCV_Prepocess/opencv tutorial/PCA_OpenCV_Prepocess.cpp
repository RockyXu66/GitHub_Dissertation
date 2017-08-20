#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>
#include <Eigen/Eigen>
//#include <string>

using namespace cv;
using namespace std;
using namespace Eigen;

// Number of components to keep for the PCA:
const int num_components = 100;
const int smallerNum = 100;      // Total num of components in saved files
const int cell_dimension = 30;  // Dimension of cells
const int imageIndex = 9;
const int image_num = 899;      //120
const int image_width = 540;   //960
const int image_height = image_width;  //960

//String oneImagePath = "/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/artifix_120/artifix1.png";
String oneImagePath = "/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution"+to_string(image_width)+")/head2.png";

String m_name = to_string(smallerNum)+"_cells"+to_string(cell_dimension);
String d_name = "head";
String file_PCA_b = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells"+to_string(image_width)+"/PCA"+m_name+"_b("+d_name+to_string(image_width)+").txt";
String file_PCA_g = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells"+to_string(image_width)+"/PCA"+m_name+"_g("+d_name+to_string(image_width)+").txt";
String file_PCA_r = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells"+to_string(image_width)+"/PCA"+m_name+"_r("+d_name+to_string(image_width)+").txt";
String file_Scores_b = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells"+to_string(image_width)+"/Scores"+m_name+"_b("+d_name+to_string(image_width)+").txt";
String file_Scores_g = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells"+to_string(image_width)+"/Scores"+m_name+"_g("+d_name+to_string(image_width)+").txt";
String file_Scores_r = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells"+to_string(image_width)+"/Scores"+m_name+"_r("+d_name+to_string(image_width)+").txt";

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

void savePCA(const string &file_name,PCA pca_){
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

// Only save num_components PCA eienvectors
void savePCA(const string &file_name,vector<PCA> pcas, int cell_num, int num_components){
    FileStorage fs(file_name,FileStorage::WRITE);
    Mat means(cell_num, pcas[0].mean.cols, CV_32F);
    Mat eigenValues(pcas[0].eigenvalues.rows, num_components, CV_32F);
    Mat eigenVectors(pcas[0].eigenvectors.rows*num_components, pcas[0].eigenvectors.cols, CV_32F);
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

PCA loadPCA(const string &file_name){
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
    if(pca_.mean.empty()){
        cout<<"Cannot load file "<<file_name<<endl;
    }
    fs.release();
    return pca_;
}

void saveScores(const string &file_name, Mat scores){
    FileStorage fs(file_name, FileStorage::WRITE);
    fs << "scores" << scores;
    fs.release();
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

void saveSmallerScores(const string &file_name, const string &new_fileName, int num){
    Mat scores;
    FileStorage fs(file_name,FileStorage::READ);
    fs["scores"] >> scores;
    if(scores.empty()){
        cout<<"Cannot load file "<<file_name<<endl;
        return;
    }
    fs.release();

    Mat new_scores = scores(Rect(0,0,num,scores.rows));
    FileStorage nfs(new_fileName, FileStorage::WRITE);
    nfs << "scores" << new_scores;
    nfs.release();
    cout<<"Save new scores successfully"<<endl;
}

void saveSmallerPCA(const string &file_name, const string &new_fileName, int num, int cell_num, int originalNum){
    PCA pca_;
    FileStorage fs(file_name,FileStorage::READ);
    fs["mean"] >> pca_.mean;
    if(pca_.mean.empty()){
        cout<<"Cannot load file "<<file_name<<endl;
        return;
    }
    fs["e_values"] >> pca_.eigenvalues;
    fs["e_vectors"] >> pca_.eigenvectors;
    fs.release();
    
    pca_.eigenvalues = pca_.eigenvalues(Rect(0,0,pca_.eigenvalues.cols,num));
    Mat tmpVectors = pca_.eigenvectors;
    pca_.eigenvectors = pca_.eigenvectors(Rect(0,0,pca_.eigenvectors.cols,num*cell_num));
    for(int i=0; i<cell_num; i++){
        Mat tmp = tmpVectors(Rect(0,i*originalNum,pca_.eigenvectors.cols,num));
        tmp.copyTo(pca_.eigenvectors(Rect(0,i*num,pca_.eigenvectors.cols,num)));
    }
    FileStorage nfs(new_fileName, FileStorage::WRITE);
    nfs << "mean" << pca_.mean;
    nfs << "e_values" << pca_.eigenvalues;
    nfs << "e_vectors" << pca_.eigenvectors;
    nfs.release();
    cout<<"Save new pca successfully"<<endl;
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

void calculatePCA(int cell_num){
    // Crop to cells
    cout<<"Start loading images:"<<endl;
    auto t_load1 = chrono::high_resolution_clock::now();
    vector<Mat> cells_b, cells_g, cells_r;
    for(int i=0; i<image_height/cell_dimension; i++){
        auto approximate1 = chrono::high_resolution_clock::now();
        for(int j=0; j<image_width/cell_dimension; j++){
            vector<Mat> cell_b, cell_g, cell_r;
            for(int n=1; n<image_num+1; n++){
                Mat img = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900(resolution"+to_string(image_width)+")/head" + to_string(n+1)+".png");
//                Mat img = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/artifix_120/artifix" + to_string(n+1)+".png");
                //                Mat img = imread("/Users/yinghanxu/Study/GitHub_Dissertation/PCA_opencv_0712/opencv tutorial/textures/new_Model_screenShot/"+to_string(n+1)+".png");
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
        auto approximate2 = chrono::high_resolution_clock::now();
        if(i==0){
            chrono::duration<double> t_approximate = approximate2 - approximate1;
            cout<<"Single loop time: "<<float(t_approximate.count())/60.0f<<" min"<<endl;
            cout<<"Approximate total time: "<<(float(t_approximate.count())/60.0f)*(image_height/cell_dimension)<<" min"<<endl;
        }
        cout<<i+1<<"/"<<(image_height/cell_dimension)<<" ";
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
        auto approximate1 = chrono::high_resolution_clock::now();
        Mat data_b = cells_b[i];
        Mat data_g = cells_g[i];
        Mat data_r = cells_r[i];
        
        PCA pca_b(data_b, Mat(), CV_PCA_DATA_AS_ROW, num_components);
        PCA pca_g(data_g, Mat(), CV_PCA_DATA_AS_ROW, num_components);
        PCA pca_r(data_r, Mat(), CV_PCA_DATA_AS_ROW, num_components);
        pcas_b.push_back(pca_b);
        pcas_g.push_back(pca_g);
        pcas_r.push_back(pca_r);
        auto approximate2 = chrono::high_resolution_clock::now();
        if(i==0){
            chrono::duration<double> t_approximate = approximate2 - approximate1;
            cout<<"Single loop time: "<<float(t_approximate.count())/60.0f<<" min"<<endl;
            cout<<"Approximate total time: "<<(float(t_approximate.count())/60.0f)*cell_num<<" min"<<endl;
        }
        cout<<i+1<<"/"<<cell_num+1<<" ";
    }
    cout<<endl;
    auto t_pca2 = chrono::high_resolution_clock::now();
    chrono::duration<double> t_pca_elapsed = t_pca2 - t_pca1;
    cout << "Finish calculating PCA. Duration time: "<<float(t_pca_elapsed.count()) / 60.0f << " min"<< endl;
    
    // Save the PCA result
    cout<<"===Save to files==="<<endl;
    savePCA(file_PCA_b, pcas_b, cell_num);
    savePCA(file_PCA_g, pcas_g, cell_num);
    savePCA(file_PCA_r, pcas_r, cell_num);
    cout<<"Finish saving to files."<<endl;
    
    cout<<"===Calculate scores and save to files==="<<endl;
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
        cout<<i+1<<" ";
    }
    cout<<endl;
//    saveScores(file_Scores_b1, scores_cells_b(Rect(0,0,pca_dim,(cell_num*image_num)/2)));
//    saveScores(file_Scores_b2, scores_cells_b(Rect(0,(cell_num*image_num)/2,pca_dim,(cell_num*image_num)/2)));
//    saveScores(file_Scores_g1, scores_cells_g(Rect(0,0,pca_dim,(cell_num*image_num)/2)));
//    saveScores(file_Scores_g2, scores_cells_g(Rect(0,(cell_num*image_num)/2,pca_dim,(cell_num*image_num)/2)));
//    saveScores(file_Scores_r1, scores_cells_r(Rect(0,0,pca_dim,(cell_num*image_num)/2)));
//    saveScores(file_Scores_r2, scores_cells_r(Rect(0,(cell_num*image_num)/2,pca_dim,(cell_num*image_num)/2)));
    saveScores(file_Scores_b, scores_cells_b);
    saveScores(file_Scores_g, scores_cells_g);
    saveScores(file_Scores_r, scores_cells_r);
    cout<<"Calculate scores and save to files."<<endl;
}

void reconstruction(int cell_num, Mat bgr[3], Mat mean_bgr[3]){
    // Load the PCA result
    cout<<"===Load PCA results==="<<endl;
    auto load1 = chrono::high_resolution_clock::now();
    PCA pca_b = loadPCA(file_PCA_b);
    PCA pca_g = loadPCA(file_PCA_g);
    PCA pca_r = loadPCA(file_PCA_r);
    cout<<"===Load scores==="<<endl;
    Mat scores_b = loadScores(file_Scores_b);
    Mat scores_g = loadScores(file_Scores_g);
    Mat scores_r = loadScores(file_Scores_r);
    auto load2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_load = load2 - load1;
    cout<<"Finish loading PCA results and scores. (Duration time: "<<elapsed_load.count()<<" s)" <<endl;
    
    
    // Reconstruction
    cout<<"===Start reconstruction==="<<endl;
    auto re1 = chrono::high_resolution_clock::now();
    
    int x = image_height/cell_dimension;
    int y = image_width/cell_dimension;
    
//    x = y = 1;
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            // Get cell image pca and score
//            auto tt1 = chrono::high_resolution_clock::now();
            PCA pcaB, pcaG, pcaR;
            pcaB.mean = pca_b.mean(Rect(0,j+i*y,cell_dimension*cell_dimension,1));
            pcaB.eigenvalues = pca_b.eigenvalues(Rect(j+i*y,0,1,num_components));
            pcaB.eigenvectors = pca_b.eigenvectors(Rect(0,(j+i*y)*smallerNum,cell_dimension*cell_dimension,num_components));
            pcaG.mean = pca_g.mean(Rect(0,j+i*y,cell_dimension*cell_dimension,1));
            pcaG.eigenvalues = pca_g.eigenvalues(Rect(j+i*y,0,1,num_components));
            pcaG.eigenvectors = pca_g.eigenvectors(Rect(0,(j+i*y)*smallerNum,cell_dimension*cell_dimension,num_components));
            pcaR.mean = pca_r.mean(Rect(0,j+i*y,cell_dimension*cell_dimension,1));
            pcaR.eigenvalues = pca_r.eigenvalues(Rect(j+i*y,0,1,num_components));
            pcaR.eigenvectors = pca_r.eigenvectors(Rect(0,(j+i*y)*smallerNum,cell_dimension*cell_dimension,num_components));
            Mat scoreB = scores_b(Rect(0,imageIndex*cell_num+j+i*y,num_components,1));
            Mat scoreG = scores_g(Rect(0,imageIndex*cell_num+j+i*y,num_components,1));
            Mat scoreR = scores_r(Rect(0,imageIndex*cell_num+j+i*y,num_components,1));
//            auto tt2 = chrono::high_resolution_clock::now();
//            chrono::duration<double> e_tt = tt2 - tt1;
//            cout<<"cout time: "<<e_tt.count()<<endl;
            
            
//            Mat Result(1,400,CV_32F);
//            for(int i=0; i<400; i++){
//                Result.at<float>(0,i) = 0;
//            }
            // Do backproject without using pca function
//            auto rr1 = chrono::high_resolution_clock::now();
//            for(int i=0; i<20; i++){
//                Result = Result + scoreB.at<float>(0,i) * pcaB.eigenvectors.row(i);
//            }
//            Result = Result + pcaB.mean;
//            Result = Result.reshape(1,20);
//            cout<<Result.rows<<" "<<Result.cols<<endl;
//            auto rr2 = chrono::high_resolution_clock::now();
//            chrono::duration<double> e_rr = rr2 - rr1;
//            cout<<"without opencl speedup time: "<<e_rr.count()<<endl;
            
//            auto t1 = chrono::high_resolution_clock::now();
            Mat b = pcaB.backProject(scoreB).reshape(1,cell_dimension);
//            auto t2 = chrono::high_resolution_clock::now();
//            chrono::duration<double> e_t = t2 - t1;
//            cout<<"pca with opencl speedup time: "<<e_t.count()<<endl;
//            cout<<"OpenCL speed up times: "<<e_rr.count()/e_t.count()<<endl;
            Mat g = pcaG.backProject(scoreG).reshape(1,cell_dimension);
            Mat r = pcaR.backProject(scoreR).reshape(1,cell_dimension);
            
            
            
//            cout<<norm_0_255(Result.reshape(1,cell_dimension))<<endl;;
//            cout<<norm_0_255(b)<<endl;
            
//            auto t2 = chrono::high_resolution_clock::now();
//            chrono::duration<double> e_t = t2 - t1;
//            cout<<"pca time: "<<e_t.count()<<endl;
//            auto t2 = chrono::high_resolution_clock::now();
//            chrono::duration<double> e_t = t2 - t1;
//            cout<<"pca time: "<<e_t.count()<<endl;
            
            // Mean image
//            pcaB.mean.reshape(1,cell_dimension).copyTo(mean_bgr[0](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
//            pcaG.mean.reshape(1,cell_dimension).copyTo(mean_bgr[1](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
//            pcaR.mean.reshape(1,cell_dimension).copyTo(mean_bgr[2](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            
            // Reconstruction image
            b.copyTo(bgr[0](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            g.copyTo(bgr[1](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            r.copyTo(bgr[2](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            
        }
    }
//        Mat reconstruction_b = norm_0_255(mean_bgr[0]);
//        Mat reconstruction_g = norm_0_255(mean_bgr[1]);
//        Mat reconstruction_r = norm_0_255(mean_bgr[2]);
//    auto n1 = chrono::high_resolution_clock::now();
    Mat reconstruction_b = norm_0_255(bgr[0]);
    Mat reconstruction_g = norm_0_255(bgr[1]);
    Mat reconstruction_r = norm_0_255(bgr[2]);
//    auto n2 = chrono::high_resolution_clock::now();
//    chrono::duration<double> e_n = n2 - n1;
//    cout<<"normalize time: "<<e_n.count()<<endl;
    
//     Put three channel into one bgr image
    Mat resultImage = mergeChannels(reconstruction_b, reconstruction_g, reconstruction_r, oneImagePath);
    auto re2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_re = re2 - re1;
    cout<<"Reconstruction time: "<<elapsed_re.count() <<" s"<<endl;
//    imwrite("ResultImages/head"+to_string(imageIndex+1+1)+"("+to_string(image_width)+")_cell"+to_string(cell_dimension)+"_"+to_string(num_components)+".png", resultImage);
//    imshow("Reconstruction", resultImage);
//    
//    cout<<"Finish reconstruction. "<<endl;
//    waitKey(0);
}

int main(int argc, const char *argv[]) {
    
    int cell_num = (image_width*image_height)/(cell_dimension*cell_dimension);
    
    Mat test(200,10,CV_32F);
    Mat newT = test(Rect(0,0,10,100));
    Mat newTT = test(Rect(0,100,10,100));
    
    Mat c = imread(oneImagePath, IMREAD_GRAYSCALE);
    c = c.reshape(1, 10);
    Mat d = c.t(); // Transpose
    cout<<c.rows<<" "<<c.cols<<endl;
    cout<<d.rows<<" "<<d.cols<<endl;
    
    calculatePCA(cell_num);
    
//    Mat bgr[3];
//    split(imread(oneImagePath), bgr);
//    Mat mean_bgr[3];
//    split(imread(oneImagePath), mean_bgr);
//    reconstruction(cell_num, bgr, mean_bgr);
    
//  //   Save smaller PCA and scores
//    int small = 10;
//    saveSmallerPCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA"+to_string(smallerNum)+"_cells"+to_string(cell_dimension)+"_b(head"+to_string(image_width)+").txt", "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA"+to_string(small)+"_cells"+to_string(cell_dimension)+"_b(head"+to_string(image_width)+").txt", small, cell_num, smallerNum);
//    saveSmallerPCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA"+to_string(smallerNum)+"_cells"+to_string(cell_dimension)+"_g(head"+to_string(image_width)+").txt", "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA"+to_string(small)+"_cells"+to_string(cell_dimension)+"_g(head"+to_string(image_width)+").txt", small, cell_num, smallerNum);
//    saveSmallerPCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA"+to_string(smallerNum)+"_cells"+to_string(cell_dimension)+"_r(head"+to_string(image_width)+").txt", "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA"+to_string(small)+"_cells"+to_string(cell_dimension)+"_r(head"+to_string(image_width)+").txt", small, cell_num, smallerNum);
//    saveSmallerScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+to_string(smallerNum)+"_cells"+to_string(cell_dimension)+"_b(head"+to_string(image_width)+").txt", "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+to_string(small)+"_cells"+to_string(cell_dimension)+"_b(head"+to_string(image_width)+").txt", small);
//    saveSmallerScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+to_string(smallerNum)+"_cells"+to_string(cell_dimension)+"_g(head"+to_string(image_width)+").txt", "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+to_string(small)+"_cells"+to_string(cell_dimension)+"_g(head"+to_string(image_width)+").txt", small);
//    saveSmallerScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+to_string(smallerNum)+"_cells"+to_string(cell_dimension)+"_r(head"+to_string(image_width)+").txt", "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+to_string(small)+"_cells"+to_string(cell_dimension)+"_r(head"+to_string(image_width)+").txt", small);
    
    
    /*
    string line;
    ifstream myfile ("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA_cells_b(head)_test.txt");
//    ofstream newfile ("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA_cells_b(head)_test.txt");
    if (myfile.is_open())
    {
//        for ( int i=0; i<600000; i++ )//0-289843-mean 289844-597162-evalues 597162
//        {
//            getline (myfile,line);
////            cout << line << '\n';
//            if(i>=289844&&i<=597162){
//                newfile << line<<'\n';
//            }
//        }
        int i=0;
        while(getline(myfile, line)){
            if(i>5){
                //newfile << line<<'\n';
                break;
            }
            std::size_t pos = line.find("       ");
            if(i==4){
                string str3 = line.substr (pos+7);
                cout<<str3<<endl;
            }
            
            i++;
        }
        myfile.close();
//        newfile.close();
    }
    
    else cout << "Unable to open file";
     */
    
    return 0;
}





































