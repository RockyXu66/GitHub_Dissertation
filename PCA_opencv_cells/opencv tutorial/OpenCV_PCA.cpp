#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>
//#include <string>

using namespace cv;
using namespace std;

// Number of components to keep for the PCA:
const int num_components = 30;
const int smallerNum = 80;
const int cell_dimension = 20;
const int imageIndex = 10;
const int image_num = 900;      //120
const int image_width = 1080;   //960
const int image_height = 1080;  //960

//String oneImagePath = "/Users/yinghanxu/Study/GitHub_Dissertation/PCA_opencv_0712/opencv tutorial/textures/new_Model_screenShot/1.png";
//String oneImagePath = "/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/artifix_120/artifix1.png";
String oneImagePath = "/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900/head2.png";

String m_name = to_string(smallerNum)+"_cells"+to_string(cell_dimension);
String file_PCA_b = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA"+m_name+"_b(head).txt";
String file_PCA_g = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA"+m_name+"_g(head).txt";
String file_PCA_r = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA"+m_name+"_r(head).txt";
//String file_Scores_b = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+m_name+"_b(head).txt";
//String file_Scores_g = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+m_name+"_g(head).txt";
//String file_Scores_r = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+m_name+"_r(head).txt";

String file_Scores_b1 = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+m_name+"_b1(head).txt";
String file_Scores_b2 = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+m_name+"_b2(head).txt";
String file_Scores_g1 = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+m_name+"_g1(head).txt";
String file_Scores_g2 = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+m_name+"_g2(head).txt";
String file_Scores_r1 = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+m_name+"_r1(head).txt";
String file_Scores_r2 = "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores"+m_name+"_r2(head).txt";

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
    cout<<"Finish load one scores file"<<endl;
    return scores;
}

void saveSmallerScores(const string &file_name, const string &new_fileName, int samllerNum){
    Mat scores;
    FileStorage fs(file_name,FileStorage::READ);
    fs["scores"] >> scores;
    if(scores.empty()){
        cout<<"Cannot load file "<<file_name<<endl;
    }
    fs.release();

    scores = scores(Rect(0,0,samllerNum,scores.rows));
    FileStorage nfs(new_fileName, FileStorage::WRITE);
    nfs << "scores" << scores;
    nfs.release();
}

void saveSmallerPCA(const string &file_name, const string &new_fileName, int samllerNum, int cell_num, int image_num){
    PCA pca_;
    FileStorage fs(file_name,FileStorage::READ);
    fs["mean"] >> pca_.mean;
    if(pca_.mean.empty()){
        cout<<"Cannot load file "<<file_name<<endl;
    }
    fs["e_values"] >> pca_.eigenvalues;
    fs["e_vectors"] >> pca_.eigenvectors;
    fs.release();
    
    pca_.eigenvalues = pca_.eigenvalues(Rect(0,0,pca_.eigenvalues.cols,samllerNum));
    Mat tmpVectors = pca_.eigenvectors;
    pca_.eigenvectors = pca_.eigenvectors(Rect(0,0,pca_.eigenvectors.cols,samllerNum*cell_num));
    for(int i=0; i<cell_num; i++){
        Mat tmp = tmpVectors(Rect(0,i*image_num,pca_.eigenvectors.cols,samllerNum));
        tmp.copyTo(pca_.eigenvectors(Rect(0,i*samllerNum,pca_.eigenvectors.cols,samllerNum)));
    }
    FileStorage nfs(new_fileName, FileStorage::WRITE);
    nfs << "mean" << pca_.mean;
    nfs << "e_values" << pca_.eigenvalues;
    nfs << "e_vectors" << pca_.eigenvectors;
    nfs.release();
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
        for(int j=0; j<image_width/cell_dimension; j++){
            vector<Mat> cell_b, cell_g, cell_r;
            for(int n=0; n<image_num; n++){
                Mat img = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900/head" + to_string(n+1)+".png");
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
        cout<<endl;
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
        
        PCA pca_b(data_b, Mat(), CV_PCA_DATA_AS_ROW, num_components);
        PCA pca_g(data_g, Mat(), CV_PCA_DATA_AS_ROW, num_components);
        PCA pca_r(data_r, Mat(), CV_PCA_DATA_AS_ROW, num_components);
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
    savePCA(file_PCA_b, pcas_b, cell_num);
    savePCA(file_PCA_g, pcas_g, cell_num);
    savePCA(file_PCA_r, pcas_r, cell_num);
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
        cout<<i+1<<" ";
    }
    cout<<endl;
    saveScores(file_Scores_b1, scores_cells_b(Rect(0,0,pca_dim,(cell_num*image_num)/2)));
    saveScores(file_Scores_b2, scores_cells_b(Rect(0,(cell_num*image_num)/2,pca_dim,(cell_num*image_num)/2)));
    saveScores(file_Scores_g1, scores_cells_g(Rect(0,0,pca_dim,(cell_num*image_num)/2)));
    saveScores(file_Scores_g2, scores_cells_g(Rect(0,(cell_num*image_num)/2,pca_dim,(cell_num*image_num)/2)));
    saveScores(file_Scores_r1, scores_cells_r(Rect(0,0,pca_dim,(cell_num*image_num)/2)));
    saveScores(file_Scores_r2, scores_cells_r(Rect(0,(cell_num*image_num)/2,pca_dim,(cell_num*image_num)/2)));
//    saveScores(file_Scores_b, scores_cells_b);
//    saveScores(file_Scores_g, scores_cells_g);
//    saveScores(file_Scores_r, scores_cells_r);
    cout<<"Finish saving scores to files."<<endl;
}

void reconstruction(int cell_num){
    // Load the PCA result
    cout<<"===Load PCA results==="<<endl;
    auto load1 = chrono::high_resolution_clock::now();
    PCA pca_b = loadPCA(file_PCA_b);
    PCA pca_g = loadPCA(file_PCA_g);
    PCA pca_r = loadPCA(file_PCA_r);
    cout<<"===Load scores==="<<endl;
    Mat scores_b1 = loadScores(file_Scores_b1);
    Mat scores_b2 = loadScores(file_Scores_b2);
    Mat scores_g1 = loadScores(file_Scores_g1);
    Mat scores_g2 = loadScores(file_Scores_g2);
    Mat scores_r1 = loadScores(file_Scores_r1);
    Mat scores_r2 = loadScores(file_Scores_r2);
    Mat scores_b(cell_num*image_num, smallerNum, CV_32F);
    Mat scores_g(cell_num*image_num, smallerNum, CV_32F);
    Mat scores_r(cell_num*image_num, smallerNum, CV_32F);
    scores_b1.copyTo(scores_b(Rect(0,0,smallerNum,(cell_num*image_num)/2)));
    scores_b2.copyTo(scores_b(Rect(0,(cell_num*image_num)/2,smallerNum,(cell_num*image_num)/2)));
    scores_g1.copyTo(scores_g(Rect(0,0,smallerNum,(cell_num*image_num)/2)));
    scores_g2.copyTo(scores_g(Rect(0,(cell_num*image_num)/2,smallerNum,(cell_num*image_num)/2)));
    scores_r1.copyTo(scores_r(Rect(0,0,smallerNum,(cell_num*image_num)/2)));
    scores_r2.copyTo(scores_r(Rect(0,(cell_num*image_num)/2,smallerNum,(cell_num*image_num)/2)));
    auto load2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_load = load2 - load1;
    cout<<"Finish loading PCA results and scores. (Duration time: "<<float(elapsed_load.count()) / 60.0f << " min)" <<endl;
    
    
    // Reconstruction
    cout<<"===Start reconstruction==="<<endl;
    auto re1 = chrono::high_resolution_clock::now();
    // Get all pcas
    vector<PCA> pcas_b, pcas_g, pcas_r;
    
    Mat mean_bgr[3], bgr[3];
    split(imread(oneImagePath), mean_bgr);
    split(imread(oneImagePath), bgr);
    int x = image_height/cell_dimension;
    int y = image_width/cell_dimension;
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            // Get cell image pca and score
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
            Mat b = pcaB.backProject(scoreB).reshape(1,cell_dimension);
            Mat g = pcaG.backProject(scoreG).reshape(1,cell_dimension);
            Mat r = pcaR.backProject(scoreR).reshape(1,cell_dimension);
            
            // Mean image
            pcaB.mean.reshape(1,cell_dimension).copyTo(mean_bgr[0](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            pcaG.mean.reshape(1,cell_dimension).copyTo(mean_bgr[1](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            pcaR.mean.reshape(1,cell_dimension).copyTo(mean_bgr[2](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            
            // Reconstruction image
            b.copyTo(bgr[0](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            g.copyTo(bgr[1](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            r.copyTo(bgr[2](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
        }
    }
    //    Mat reconstruction_b = norm_0_255(mean_bgr[0]);
    //    Mat reconstruction_g = norm_0_255(mean_bgr[1]);
    //    Mat reconstruction_r = norm_0_255(mean_bgr[2]);
    Mat reconstruction_b = norm_0_255(bgr[0]);
    Mat reconstruction_g = norm_0_255(bgr[1]);
    Mat reconstruction_r = norm_0_255(bgr[2]);
    
    // Put three channel into one bgr image
    Mat resultImage = mergeChannels(reconstruction_b, reconstruction_g, reconstruction_r, oneImagePath);
    imwrite("ResultImages/head"+to_string(imageIndex+1)+"_cell"+to_string(cell_dimension)+"_"+to_string(num_components)+".png", resultImage);
    imshow("Reconstruction", resultImage);
    auto re2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_re = re2 - re1;
    cout<<"Finish reconstruction. (Duration time: "<<float(elapsed_re.count()) / 60.0f << " min)" <<endl;
    waitKey(0);
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
    
//    calculatePCA(cell_num);
    /*
    // Crop to cells
    cout<<"Start loading images:"<<endl;
    auto t_load1 = chrono::high_resolution_clock::now();
    vector<Mat> cells_b, cells_g, cells_r;
    for(int i=0; i<image_height/cell_dimension; i++){
        for(int j=0; j<image_width/cell_dimension; j++){
            vector<Mat> cell_b, cell_g, cell_r;
            for(int n=0; n<image_num; n++){
//                Mat img = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/head_900/head" + to_string(n+1)+".png");
                Mat img = imread("/Users/yinghanxu/Study/Dissertation_ResultData/Data_Set/artifix_120/artifix" + to_string(n+1)+".png");
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
    
        PCA pca_b(data_b, Mat(), CV_PCA_DATA_AS_ROW, num_components);
        PCA pca_g(data_g, Mat(), CV_PCA_DATA_AS_ROW, num_components);
        PCA pca_r(data_r, Mat(), CV_PCA_DATA_AS_ROW, num_components);
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
    savePCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA50_cells30_b(chest).txt", pcas_b, cell_num);
    savePCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA50_cells30_g(chest).txt", pcas_g, cell_num);
    savePCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA50_cells30_r(chest).txt", pcas_r, cell_num);
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
        cout<<i+1<<" ";
    }
    cout<<endl;
    saveScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores50_cells30_b(chest).txt", scores_cells_b);
    saveScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores50_cells30_g(chest).txt", scores_cells_g);
    saveScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores50_cells30_r(chest).txt", scores_cells_r);
    cout<<"Finish saving scores to files."<<endl;
     */
    
    
    reconstruction(cell_num);
    /*
    // Load the PCA result
    cout<<"===Load PCA results==="<<endl;
    auto load1 = chrono::high_resolution_clock::now();
    PCA pca_b = loadPCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA50_cells30_b(chest).txt");
    PCA pca_g = loadPCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA50_cells30_g(chest).txt");
    PCA pca_r = loadPCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA50_cells30_r(chest).txt");
    cout<<"===Load scores==="<<endl;
    Mat scores_b = loadScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores50_cells30_b(chest).txt");
    Mat scores_g = loadScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores50_cells30_g(chest).txt");
    Mat scores_r = loadScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores50_cells30_r(chest).txt");
    auto load2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_load = load2 - load1;
    cout<<"Finish loading PCA results and scores. (Duration time: "<<float(elapsed_load.count()) / 60.0f << " min)" <<endl;
    
    
    // Reconstruction
    cout<<"===Start reconstruction==="<<endl;
    auto re1 = chrono::high_resolution_clock::now();
    // Get all pcas
    vector<PCA> pcas_b, pcas_g, pcas_r;
    
    Mat mean_bgr[3], bgr[3];
    split(imread(oneImagePath), mean_bgr);
    split(imread(oneImagePath), bgr);
    int x = image_height/cell_dimension;
    int y = image_width/cell_dimension;
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            // Get cell image pca and score
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
            Mat b = pcaB.backProject(scoreB).reshape(1,cell_dimension);
            Mat g = pcaG.backProject(scoreG).reshape(1,cell_dimension);
            Mat r = pcaR.backProject(scoreR).reshape(1,cell_dimension);
            
            // Mean image
            pcaB.mean.reshape(1,cell_dimension).copyTo(mean_bgr[0](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            pcaG.mean.reshape(1,cell_dimension).copyTo(mean_bgr[1](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            pcaR.mean.reshape(1,cell_dimension).copyTo(mean_bgr[2](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            
            // Reconstruction image
            b.copyTo(bgr[0](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            g.copyTo(bgr[1](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            r.copyTo(bgr[2](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
        }
    }
//    Mat reconstruction_b = norm_0_255(mean_bgr[0]);
//    Mat reconstruction_g = norm_0_255(mean_bgr[1]);
//    Mat reconstruction_r = norm_0_255(mean_bgr[2]);
    Mat reconstruction_b = norm_0_255(bgr[0]);
    Mat reconstruction_g = norm_0_255(bgr[1]);
    Mat reconstruction_r = norm_0_255(bgr[2]);
    auto re2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_re = re2 - re1;
    cout<<"Finish reconstruction. (Duration time: "<<float(elapsed_re.count()) / 60.0f << " min)" <<endl;
    
    // Put three channel into one bgr image
    Mat resultImage = mergeChannels(reconstruction_b, reconstruction_g, reconstruction_r, oneImagePath);
    imwrite("ResultImages/artifix"+to_string(imageIndex+1)+"_cell"+to_string(cell_dimension)+"_"+to_string(num_components)+".png", resultImage);
    imshow("Reconstruction", resultImage);
    waitKey(0);
    */
    
//    saveSmallerScores("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores_cells_g(chest).txt", "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/Scores40_cells_g(chest).txt", 40);
//    saveSmallerPCA("/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA_cells_g(chest).txt", "/Users/yinghanxu/Study/Dissertation_ResultData/ResultPCA_cells/PCA40_cells_g(chest).txt", 40, cell_num, image_num);
    
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





































