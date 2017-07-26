//
//  PCA_.cpp
//  PCA_opengl_0522
//
//  Created by Yinghan Xu on 22/05/2017.
//  Copyright © 2017 Yinghan Xu. All rights reserved.
//

#include "PCA.hpp"

const int imageIndex = 10;

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

void PCA_::load(){
    int imageIndex = 10;
    
    int cell_num = (image_width*image_height)/(cell_dimension*cell_dimension);
    
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
    cout<<"Finish loading PCA results and scores. (Duration time: "<<float(elapsed_load.count()) / 60.0f << " min)" <<endl;
    
    int x = image_height/cell_dimension;
    int y = image_width/cell_dimension;
    vector<PCA> pcas_b, pcas_g, pcas_r;
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
            
            pcas_b.push_back(pcaB);
            pcas_g.push_back(pcaG);
            pcas_r.push_back(pcaR);
        }
    }
//    // Put three channel into one bgr image
//    Mat resultImage = mergeChannels(reconstruction_b, reconstruction_g, reconstruction_r, oneImagePath);
//    imwrite("ResultImages/head"+to_string(imageIndex+1)+"_"+to_string(num_components)+".png", resultImage);
//    
//    string p = "ResultImages/head"+to_string(imageIndex+1)+"_"+to_string(num_components)+".png";
////    char* path ="ResultImages/head"+to_string(imageIndex+1)+"_"+to_string(num_components)+".png";
//    char* path = new char[p.size() + 1];
//    copy(p.begin(), p.end(), path);
//    path[p.size()] = '\0';
    
    // Save PCA results and scores for real time
    this->pcas_b = pcas_b;
    this->pcas_g = pcas_g;
    this->pcas_r = pcas_r;
    this->scores_b = scores_b;
    this->scores_g = scores_g;
    this->scores_r = scores_r;
    
//    return path;
}

unsigned char* PCA_::reconstruct(int imageIndex, unsigned char* newImage, Mat bgr[3], int cell_num, int x, int y){
    
    // Reconstruction
//    cout<<"===Start reconstruction==="<<endl;
//    auto re1 = chrono::high_resolution_clock::now();
    
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            // Reconstruction image
            int index = j+i*y;
            Mat scoreB = scores_b(Rect(0,imageIndex*cell_num+j+i*y,num_components,1));
            Mat scoreG = scores_g(Rect(0,imageIndex*cell_num+j+i*y,num_components,1));
            Mat scoreR = scores_r(Rect(0,imageIndex*cell_num+j+i*y,num_components,1));
            Mat b = pcas_b[index].backProject(scoreB).reshape(1,cell_dimension);
            Mat g = pcas_g[index].backProject(scoreG).reshape(1,cell_dimension);
            Mat r = pcas_r[index].backProject(scoreR).reshape(1,cell_dimension);
            
            b.copyTo(bgr[0](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            g.copyTo(bgr[1](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            r.copyTo(bgr[2](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
        }
    }
    Mat reconstruction_b = norm_0_255(bgr[0]);
    Mat reconstruction_g = norm_0_255(bgr[1]);
    Mat reconstruction_r = norm_0_255(bgr[2]);
    // Put three channel into one bgr image
    Mat resultImage = mergeChannels(reconstruction_b, reconstruction_g, reconstruction_r, oneImagePath);
    
//    auto re2 = chrono::high_resolution_clock::now();
//    chrono::duration<double> elapsed_re = re2 - re1;
//    cout<<"Finish reconstruction. (Duration time: "<<elapsed_re.count()<<")"<<endl;
    
//    imwrite("ResultImages/head"+to_string(imageIndex+1)+"_"+to_string(num_components)+".png", resultImage);
//    string p = "ResultImages/head"+to_string(imageIndex+1)+"_"+to_string(num_components)+".png";
//    char* path = new char[p.size() + 1];
//    copy(p.begin(), p.end(), path);
//    path[p.size()] = '\0';
//    return path;
    
//    auto p1 = chrono::high_resolution_clock::now();
    // Put data matrix into each channel of the image
    for (int x = 0; x<image_width; x++) {
        for(int y=0; y<image_height; y++){
            newImage[(y+x*image_height) * 3] = resultImage.at<cv::Vec3b>(x,y)[2];
            newImage[(y+x*image_height) * 3 + 1] = resultImage.at<cv::Vec3b>(x,y)[1];
            newImage[(y+x*image_height) * 3 + 2] = resultImage.at<cv::Vec3b>(x,y)[0];
        }
    }
//    auto p2 = chrono::high_resolution_clock::now();
//    chrono::duration<double> elapsed_p = p2 - p1;
//    cout<<"Put into image. (Duration time: "<<elapsed_p.count()<<")"<<endl;
    
    return newImage;
    
}

unsigned char* PCA_::Reconstruction_RealTime(unsigned char* image, int imageIndex, int width, int height){
    unsigned char* newImage = image;
    
    string meanFileName, scoresFileName;
    meanFileName = "eigenVectors_cell/mean_Red.txt";
    scoresFileName = "eigenVectors_cell/scores40_Red.txt";
    vector<MatrixXd> ReconstructedData_Red = WholeReconstruction_RealTime(meanFileName, scoresFileName, imageIndex, "Red");
    meanFileName = "eigenVectors_cell/mean_Green.txt";
    scoresFileName = "eigenVectors_cell/scores40_Green.txt";
    vector<MatrixXd> ReconstructedData_Green = WholeReconstruction_RealTime(meanFileName, scoresFileName, imageIndex, "Green");
    meanFileName = "eigenVectors_cell/mean_Blue.txt";
    scoresFileName = "eigenVectors_cell/scores40_Blue.txt";
    vector<MatrixXd> ReconstructedData_Blue = WholeReconstruction_RealTime(meanFileName, scoresFileName, imageIndex, "Blue");
    
    // Put cells together into one matrix
    MatrixXd data_r_R(1, width*height);
    MatrixXd data_g_R(1, width*height);
    MatrixXd data_b_R(1, width*height);
    for (int i = 0; i < height/20; i++) {
        for (int j = 0; j < width / 20; j++) {
            MatrixXd data_r_cell = ReconstructedData_Red[j+(height/20)*i];
            MatrixXd data_g_cell = ReconstructedData_Green[j+(height/20)*i];
            MatrixXd data_b_cell = ReconstructedData_Blue[j+(height/20)*i];
            
            for (int x = 0; x < 20; x++) {
                for (int y = 0; y < 20; y++) {
                    int piIndex = j * 20 + i * 300 * 20 + y + x*300;
                    data_r_R(0, piIndex) = data_r_cell(0, y + 20 * x);
                    data_g_R(0, piIndex) = data_g_cell(0, y + 20 * x);
                    data_b_R(0, piIndex) = data_b_cell(0, y + 20 * x);
                }
            }
        }
    }
    
    // Put data matrix into each channel of the image
    for (int i = 0; i<width*height; i++) {
        newImage[i * 3] = data_r_R(0, i);
        newImage[i * 3 + 1] = data_g_R(0, i);
        newImage[i * 3 + 2] = data_b_R(0, i);
    }
    
    return newImage;
}

unsigned char* PCA_::CalculateEigen(vector<unsigned char*> images, int width, int height) {
    
    this->num = images.size();
    this->dim = 20 * 20;// width*height;
    this->cellImage_num = (width*height)/(20*20);
    
    this->pca_dim = 40;
    
    // Show which image
    int imageIndex = 6;
    
    // Put each channel values into the data matrix
    /*
//    MatrixXd data_r(num, width*height);
//    MatrixXd data_g(num, width*height);
//    MatrixXd data_b(num, width*height);
//    for (int i = 0; i<images.size(); i++) {
//        unsigned char* image = images[i];
//        for (int j = 0; j<width*height; j++) {
//            int r = image[j * 3];
//            int g = image[j * 3 + 1];
//            int b = image[j * 3 + 2];
//            data_r(i, j) = r;
//            data_g(i, j) = g;
//            data_b(i, j) = b;
//        }
//    }
//    
//    // Crop to cell => 20*20 pixels
//    vector<MatrixXd> data_r_cells, data_g_cells, data_b_cells;
//    for (int i = 0; i < height/20; i++) {
//        for (int j = 0; j < width / 20; j++) {
//            MatrixXd data_r_cell(num, dim);
//            MatrixXd data_g_cell(num, dim);
//            MatrixXd data_b_cell(num, dim);
//            for (int k = 0; k < num; k++) {
//                
//                for (int x = 0; x < 20; x++) {
//                    for (int y = 0; y < 20; y++) {
//                        int piIndex = j * 20 + i * 300 * 20 + y + x*300;
//                        data_r_cell(k, y + 20 * x) = data_r(k, piIndex);
//                        data_g_cell(k, y + 20 * x) = data_g(k, piIndex);
//                        data_b_cell(k, y + 20 * x) = data_b(k, piIndex);
//                    }
//                }
//            }
//            
//            data_r_cells.push_back(data_r_cell);
//            data_g_cells.push_back(data_g_cell);
//            data_b_cells.push_back(data_b_cell);
//        }
//    }
    
    // Calculate mean image
    //unsigned char* meanImage = images[0];
    //Mean(data_r, data_g, data_b, meanImage);
    
    
    unsigned char* newImage = images[0];
    time_t currentTime;
    struct tm *localTime;
     */
    string fileName;
    string meanFileName;
    string scoresFileName;
    
    // Calculate Eigen values and Eigen vectors
    /*
//    time(&currentTime);                   // Get the current time
//    localTime = localtime(&currentTime);
//    fileName = "EigenVector/EigenVector_Red_" + to_string(localTime->tm_hour)
//    + to_string(localTime->tm_min) + to_string(localTime->tm_sec) + ".txt";
//    ofstream outFile1(fileName);
//    
//    for (int i = 0; i < (width*height)/(20*20); i++) {
//        cout << "Cell image " << i << endl;
//        MatrixXd eigenVectors_r = Eigen(data_r_cells[i], newImage, "Red", fileName);
//    }
//    
//    
//    time(&currentTime);                   // Get the current time
//    localTime = localtime(&currentTime);
//    fileName = "EigenVector/EigenVector_Green_" + to_string(localTime->tm_hour)
//    + to_string(localTime->tm_min) + to_string(localTime->tm_sec) + ".txt";
//    ofstream outFile2(fileName);
//    for (int i = 0; i < (width*height) / (20 * 20); i++) {
//        cout << "Cell image " << i << endl;
//        MatrixXd eigenVectors_g = Eigen(data_g_cells[i], newImage, "Green", fileName);
//    }
    
//    time(&currentTime);                   // Get the current time
//    localTime = localtime(&currentTime);
//    fileName = "EigenVector/EigenVector_Blue_" + to_string(localTime->tm_hour)
//    + to_string(localTime->tm_min) + to_string(localTime->tm_sec) + ".txt";
//    ofstream outFile3(fileName);
//    for (int i = 0; i < (width*height) / (20 * 20); i++) {
//        cout << "Cell image " << i << endl;
//        MatrixXd eigenVectors_b = Eigen(data_b_cells[i], newImage, "Blue", fileName);
//    }
     */
    
    unsigned char* eigenTexture = images[0];
    
    // Reconstructoin one cell
    /*
//    fileName = "eigenVectors_cell/EigenVector_Red_0-150.txt";
//    MatrixXd OriginalData_Red = Reconstruction(data_r_cells[149], fileName);
//    fileName = "eigenVectors_cell/EigenVector_Green_0-150.txt";
//    MatrixXd OriginalData_Green = Reconstruction(data_g_cells[149], fileName);
//    fileName = "eigenVectors_cell/EigenVector_Blue_0-150.txt";
//    MatrixXd OriginalData_Blue = Reconstruction(data_b_cells[149], fileName);
//
//    for (int i = 0; i<dim; i++) {
//        eigenTexture[i * 3] = OriginalData_Red(0, i);
//        eigenTexture[i * 3 + 1] = OriginalData_Green(0, i);
//        eigenTexture[i * 3 + 2] = OriginalData_Blue(0, i);
//    }
     */
    
    // Save 40 high-dimensional eigenspace
    /*
    fileName = "eigenVectors_cell/EigenVector_Red_0-224.txt";
    SaveHighEigen("eigenVectors_cell/highEigen"+to_string(pca_dim)+"_Red.txt", fileName, pca_dim);
    fileName = "eigenVectors_cell/EigenVector_Green_0-224.txt";
    SaveHighEigen("eigenVectors_cell/highEigen"+to_string(pca_dim)+"_Green.txt", fileName, pca_dim);
    fileName = "eigenVectors_cell/EigenVector_Blue_0-224.txt";
    SaveHighEigen("eigenVectors_cell/highEigen"+to_string(pca_dim)+"_Blue.txt", fileName, pca_dim);
     */
    
    // Calculate mean image and save to files
    /*
    fileName = "eigenVectors_cell/mean_Red.txt";
    Mean(data_r_cells, fileName);
    fileName = "eigenVectors_cell/mean_Green.txt";
    Mean(data_g_cells, fileName);
    fileName = "eigenVectors_cell/mean_Blue.txt";
    Mean(data_b_cells, fileName);
     */
    
    // Calculate scores for eigentextures and save to files
    /*
    fileName = "eigenVectors_cell/highEigen40_Red.txt";
    meanFileName = "eigenVectors_cell/mean_Red.txt";
    scoresFileName = "eigenVectors_cell/scores40_Red.txt";
    CalScores(data_r_cells, fileName, meanFileName, scoresFileName);
    fileName = "eigenVectors_cell/highEigen40_Green.txt";
    meanFileName = "eigenVectors_cell/mean_Green.txt";
    scoresFileName = "eigenVectors_cell/scores40_Green.txt";
    CalScores(data_g_cells, fileName, meanFileName, scoresFileName);
    fileName = "eigenVectors_cell/highEigen40_Blue.txt";
    meanFileName = "eigenVectors_cell/mean_Blue.txt";
    scoresFileName = "eigenVectors_cell/scores40_Blue.txt";
    CalScores(data_b_cells, fileName, meanFileName, scoresFileName);
     */
    
    // Reconstruction the whole image
    
    fileName = "eigenVectors_cell/highEigen40_Red.txt";
    meanFileName = "eigenVectors_cell/mean_Red.txt";
    scoresFileName = "eigenVectors_cell/scores40_Red.txt";
    vector<MatrixXd> ReconstructedData_Red = WholeReconstruction(fileName, meanFileName, scoresFileName, imageIndex, "Red");
    fileName = "eigenVectors_cell/highEigen40_Green.txt";
    meanFileName = "eigenVectors_cell/mean_Green.txt";
    scoresFileName = "eigenVectors_cell/scores40_Green.txt";
    vector<MatrixXd> ReconstructedData_Green = WholeReconstruction(fileName, meanFileName, scoresFileName, imageIndex, "Green");
    fileName = "eigenVectors_cell/highEigen40_Blue.txt";
    meanFileName = "eigenVectors_cell/mean_Blue.txt";
    scoresFileName = "eigenVectors_cell/scores40_Blue.txt";
    vector<MatrixXd> ReconstructedData_Blue = WholeReconstruction(fileName, meanFileName, scoresFileName, imageIndex, "Blue");
    
    // Put cells together into one matrix
    MatrixXd data_r_R(1, width*height);
    MatrixXd data_g_R(1, width*height);
    MatrixXd data_b_R(1, width*height);
    for (int i = 0; i < height/20; i++) {
        for (int j = 0; j < width / 20; j++) {
            MatrixXd data_r_cell = ReconstructedData_Red[j+(height/20)*i];
            MatrixXd data_g_cell = ReconstructedData_Green[j+(height/20)*i];
            MatrixXd data_b_cell = ReconstructedData_Blue[j+(height/20)*i];
            
            for (int x = 0; x < 20; x++) {
                for (int y = 0; y < 20; y++) {
                    int piIndex = j * 20 + i * 300 * 20 + y + x*300;
                    data_r_R(0, piIndex) = data_r_cell(0, y + 20 * x);
                    data_g_R(0, piIndex) = data_g_cell(0, y + 20 * x);
                    data_b_R(0, piIndex) = data_b_cell(0, y + 20 * x);
                }
            }
        }
    }
    
    // Put data matrix into each channel of the image
    for (int i = 0; i<width*height; i++) {
        eigenTexture[i * 3] = data_r_R(0, i);
        eigenTexture[i * 3 + 1] = data_g_R(0, i);
        eigenTexture[i * 3 + 2] = data_b_R(0, i);
    }
    
//    return images[0];
    //return meanImage;
    return eigenTexture;
}

void PCA_::Mean(MatrixXd data_r, MatrixXd data_g, MatrixXd data_b, unsigned char* meanImage) {
    MatrixXd mean_r = data_r.colwise().mean();
    MatrixXd mean_g = data_g.colwise().mean();
    MatrixXd mean_b = data_b.colwise().mean();
    
    for (int i = 0; i<dim; i++) {
        meanImage[i * 3] = mean_r(0, i);
        meanImage[i * 3 + 1] = mean_g(0, i);
        meanImage[i * 3 + 2] = mean_b(0, i);
    }
}

void PCA_::SaveHighEigen(string newFileName, string oldFileName, int highDimension){
    
    ifstream oldfile(oldFileName);
    MatrixXd EigenVectors(cellImage_num*dim, dim);
    for (int i = 0; i < cellImage_num*dim; i++) {
        for (int j = 0; j < dim; j++) {
            oldfile >> EigenVectors(i, j);
        }
    }
    oldfile.close();
    
    ofstream newfile(newFileName);
    for (int i = 0; i < cellImage_num*dim; i++) {
        for (int j = 0; j < highDimension; j++) {
            newfile << EigenVectors(i, j) << " ";
        }
        newfile << endl;
    }
    newfile.close();
    
}

void PCA_::Mean(vector<MatrixXd> data, string fileName){
    vector<MatrixXd> means;
    
    for(int i=0; i<cellImage_num; i++){
        MatrixXd mean = data[i].colwise().mean();
        means.push_back(mean);
    }
    
    ofstream newFile(fileName);
    for(int i=0; i<cellImage_num; i++){
        for(int j=0; j<dim; j++){
            newFile << means[i](0, j)<<" ";
        }
        newFile << endl;
    }
    newFile.close();
}

void PCA_::CalScores(vector<MatrixXd> whole_data, string fileName, string meanFileName, string newFileName){
    
    // Get Eigen textures
    MatrixXd EigenVectors(cellImage_num*dim, 40);
    ifstream myfile(fileName);
    for (int i = 0; i < cellImage_num*dim; i++) {
        for (int j = 0; j < 40; j++) {
            myfile >> EigenVectors(i, j);
        }
    }
    myfile.close();
    
    vector<MatrixXd> whole_OriginalData;
    
    // Get mean image
    MatrixXd means(cellImage_num, dim);
    ifstream meanfile(meanFileName);
    for (int i = 0; i < cellImage_num; i++) {
        for (int j = 0; j < dim; j++) {
            meanfile >> means(i, j);
        }
    }
    meanfile.close();
    
    ofstream scoresFile(newFileName);
    
    for(int cellIndex=0; cellIndex<cellImage_num; cellIndex++){
        MatrixXd data = whole_data[cellIndex];
        
        MatrixXd mean = means.row(cellIndex);   //data.colwise().mean();  // mean  mean =
        MatrixXd dataAdjust(num, dim);  // data-mean
        for (int i = 0; i<num; i++) {
            dataAdjust.row(i) = data.row(i) - mean;  // dataAdjust = data-mean
        }
        
        MatrixXd oneVec(dim, 40);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < 40; j++) {
                oneVec(i, j) = EigenVectors(cellIndex*dim + i, j);
            }
        }
        
        // Calculate scores and save to the file
        for(int k=0; k<num; k++){
            for(int i=0; i<40; i++){
                double projectionValue = dataAdjust.row(k) * oneVec.col(i);
                scoresFile << projectionValue << " ";
            }
            scoresFile <<endl;
        }
    }
    
    scoresFile.close();
}


vector<MatrixXd> PCA_::WholeReconstruction(string fileName, string meanFileName, string scoresFileName, int imageIndex, string channel){
    cout<<"Reconstruction data: "<<endl<<endl;
    
    cout<<"Get Eigen Textures:"<<endl<<endl;
    // Get Eigen textures
    MatrixXd EigenVectors(cellImage_num*dim, 40);
    ifstream myfile(fileName);
    for (int i = 0; i < cellImage_num*dim; i++) {
        for (int j = 0; j < 40; j++) {
            myfile >> EigenVectors(i, j);
        }
    }
    myfile.close();
    
    // Save eigen vectors for real-time rendering
    if(channel == "Red"){
        this->EigenVectors_Red = EigenVectors;
    }else if(channel == "Green"){
        this->EigenVectors_Green = EigenVectors;
    }else if(channel == "Blue"){
        this->EigenVectors_Blue = EigenVectors;
    }
    
    
    vector<MatrixXd> whole_OriginalData;
    
    cout<<"Get mean image:"<<endl<<endl;
    // Get mean image
    MatrixXd means(cellImage_num, dim);
    ifstream meanfile(meanFileName);
    for (int i = 0; i < cellImage_num; i++) {
        for (int j = 0; j < dim; j++) {
            meanfile >> means(i, j);
        }
    }
    meanfile.close();
    
    // Save mean image for real-time rendering
    if(channel == "Red"){
        this->means_Red = means;
    }else if(channel == "Green"){
        this->means_Green = means;
    }else if(channel == "Blue"){
        this->means_Blue = means;
    }
    
    cout<<"Get scores:"<<endl<<endl;
    // Get scores
    MatrixXd whole_scores(cellImage_num*num, 40);
    ifstream scoresFile(scoresFileName);
    for (int i = 0; i < cellImage_num*num; i++) {
        for (int j = 0; j < 40; j++) {
            scoresFile >> whole_scores(i, j);
        }
    }
    scoresFile.close();

    cout<<"Get one scores:"<<endl<<endl;
    MatrixXd oneImage_scores(cellImage_num, 40);
    vector<MatrixXd> wholeImage_scores;
    for(int i=0; i<cellImage_num; i++){
        for(int j=0; j<40; j++){
            oneImage_scores(i, j) = whole_scores(imageIndex+i*num, j);
        }
    }
    
    // Save all score for real-time rendering
    for(int k=0; k<36; k++){
        for(int i=0; i<cellImage_num; i++){
            for(int j=0; j<40; j++){
                oneImage_scores(i, j) = whole_scores(k+i*num, j);
            }
        }
        wholeImage_scores.push_back(oneImage_scores);
    }
    
    
    if(channel == "Red"){
        this->wholeImage_scores_Red = wholeImage_scores;
        this->oneImage_scores_Red = oneImage_scores;
    }else if(channel == "Green"){
        this->wholeImage_scores_Green = wholeImage_scores;
        this->oneImage_scores_Green = oneImage_scores;
    }else if(channel == "Blue"){
        this->wholeImage_scores_Blue = wholeImage_scores;
        this->oneImage_scores_Blue = oneImage_scores;
    }
    
    // Save oneVectors
    vector<MatrixXd> oneVecs;
    
    for(int cellIndex=0; cellIndex<cellImage_num; cellIndex++){
//        MatrixXd data = whole_data[cellIndex];
//        
        MatrixXd mean = means.row(cellIndex);   //data.colwise().mean();  // mean  mean =
//        MatrixXd dataAdjust(num, dim);  // data-mean
//        for (int i = 0; i<num; i++) {
//            dataAdjust.row(i) = data.row(i) - mean;  // dataAdjust = data-mean
//        }
        
        MatrixXd oneVec(dim, pca_dim);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < pca_dim; j++) {
                oneVec(i, j) = EigenVectors(cellIndex*dim + i, j);
            }
        }
        oneVecs.push_back(oneVec);
        
        MatrixXd OriginalData = MatrixXd::Constant(1, dim, 0.0);
        for(int i=0; i<pca_dim; i++){
//            double projectionValue = dataAdjust.row(imageIndex) * oneVec.col(i);
            double projectionValue = oneImage_scores(cellIndex, i);
            MatrixXd eigenTexture(1, dim);
            eigenTexture = projectionValue * oneVec.col(i).transpose();
            OriginalData += eigenTexture;
        }
        
        OriginalData += mean;
        
        whole_OriginalData.push_back(OriginalData);
    }
    
    if(channel == "Red"){
        this->oneVecs_Red = oneVecs;
    }else if(channel == "Green"){
        this->oneVecs_Green = oneVecs;
    }else if(channel == "Blue"){
        this->oneVecs_Blue = oneVecs;
    }
    
    return whole_OriginalData;
}

vector<MatrixXd> PCA_::WholeReconstruction_RealTime(string meanFileName, string scoresFileName, int imageIndex, string channel){
    cout<<"Reconstruction data: "<<endl;
    
    
    vector<MatrixXd> whole_OriginalData;
    
//    cout<<"Get mean image:"<<endl;
//    // Get mean image
//    MatrixXd means(cellImage_num, dim);
//    ifstream meanfile(meanFileName);
//    for (int i = 0; i < cellImage_num; i++) {
//        for (int j = 0; j < dim; j++) {
//            meanfile >> means(i, j);
//        }
//    }
//    meanfile.close();
//    cout<<"Finish get mean image"<<endl;
    
//    cout<<"Get scores:"<<endl;
//    // Get scores
//    MatrixXd whole_scores(cellImage_num*num, 40);
//    ifstream scoresFile(scoresFileName);
//    for (int i = 0; i < cellImage_num*num; i++) {
//        for (int j = 0; j < 40; j++) {
//            scoresFile >> whole_scores(i, j);
//        }
//    }
//    scoresFile.close();
//    
//    cout<<"Get one scores:"<<endl;
//    MatrixXd oneImage_scores(cellImage_num, 40);
//    for(int i=0; i<cellImage_num; i++){
//        for(int j=0; j<40; j++){
//            oneImage_scores(i, j) = whole_scores(imageIndex+i*num, j);
//        }
//    }
//    cout<<"Finish get scores"<<endl;
    
    cout<<"Get eigenvectors and scores"<<endl;
    MatrixXd EigenVectors, oneImage_scores, means;
    vector<MatrixXd> oneVecs;
    if(channel == "Red"){
        EigenVectors = EigenVectors_Red;
        oneImage_scores = wholeImage_scores_Red[imageIndex];
        oneVecs = oneVecs_Red;
        means = means_Red;
    }else if(channel == "Green"){
        EigenVectors = EigenVectors_Green;
        oneImage_scores = wholeImage_scores_Green[imageIndex];
        oneVecs = oneVecs_Green;
        means = means_Green;
    }else if(channel == "Blue"){
        EigenVectors = EigenVectors_Blue;
        oneImage_scores = wholeImage_scores_Blue[imageIndex];
        oneVecs = oneVecs_Blue;
        means = means_Blue;
    }
    cout<<"Finish get eigenvectors and scores"<<endl;
    
    
    
    cout<<"Start reconstruction"<<endl;
    for(int cellIndex=0; cellIndex<cellImage_num; cellIndex++){
        MatrixXd mean = means.row(cellIndex);   //data.colwise().mean();  // mean  mean =
        
        MatrixXd oneVec = oneVecs[cellIndex];
        
        MatrixXd OriginalData = MatrixXd::Constant(1, dim, 0.0);
        for(int i=0; i<pca_dim; i++){
            double projectionValue = oneImage_scores(cellIndex, i);
            MatrixXd eigenTexture(1, dim);
            eigenTexture = projectionValue * oneVec.col(i).transpose();
            OriginalData += eigenTexture;
        }
        
        OriginalData += mean;
        
        whole_OriginalData.push_back(OriginalData);
    }
    cout<<"Finish reconstruction"<<endl;
    
    return whole_OriginalData;
}

// Reconstruct one cell
MatrixXd PCA_::Reconstruction(MatrixXd data, string fileName) {
//    cout << data << endl;
    MatrixXd mean = data.colwise().mean();

    MatrixXd dataAdjust(num, dim);  // data-mean
    for (int i = 0; i<num; i++) {
        dataAdjust.row(i) = data.row(i) - mean;  // dataAdjust = data-mean
    }
    
//    int cellImage_num = (300*300)/(20*20);
    MatrixXd EigenVectors(cellImage_num*dim, dim);
    ifstream myfile(fileName);
    for (int i = 0; i < cellImage_num*dim; i++) {
        for (int j = 0; j < dim; j++) {
            myfile >> EigenVectors(i, j);
        }
    }
    myfile.close();

    MatrixXd oneVec(dim, dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            oneVec(i, j) = EigenVectors(149*dim + i, j);
        }
    }
    
    cout<<oneVec<<endl;
    
    pca_dim = 40;
    
    int imageIndex = 0;
    
    MatrixXd OriginalData = MatrixXd::Constant(1, dim, 0.0);
    for(int i=0; i<pca_dim; i++){
        double projectionValue = dataAdjust.row(imageIndex) * oneVec.col(i);
        MatrixXd eigenTexture(1, dim);
        eigenTexture = projectionValue * oneVec.col(i).transpose();
        OriginalData += eigenTexture;
    }
    
    OriginalData += mean;
    
    cout<<"Reconstruction Data: "<<endl;
    cout << OriginalData << endl;
    return OriginalData;// EigenVectors.col(0).transpose() * 400;
}

MatrixXd PCA_::Eigen(MatrixXd data, unsigned char* newImage, string name, string fileName) {
    //cout << endl << endl;
    cout << name << endl;
    //cout<<data<<endl;
    MatrixXd mean = data.colwise().mean();
    
    MatrixXd dataAdjust(num, dim);  // data-mean
    MatrixXd dataAdjust_Transpose(dim, num);// transpose of dataAdjust
    
    for (int i = 0; i<num; i++) {
        dataAdjust.row(i) = data.row(i) - mean;  // dataAdjust = data-mean
    }
    dataAdjust_Transpose = dataAdjust.transpose();
    
    // Covariance matrix
    MatrixXd Covariance = (1 / (double)num) * dataAdjust_Transpose * dataAdjust;
    
    EigenSolver<MatrixXd> eigenSolver(Covariance);
    
    // Calculate eigen value and eigen vector from covariance matrix
    MatrixXd eigenValues = eigenSolver.eigenvalues().real();
    MatrixXd eigenVectors = eigenSolver.eigenvectors().real();
    
    // Sort eigen values
    MatrixXd sorted_eigenValues = BubbleSort(eigenValues, dim);
    
    MatrixXd sorted_eigenVectors(dim, dim);
    for (int i = 0; i<dim; i++)
    {
        for (int j = 0; j<dim; j++)
        {
            if (eigenValues(i, 0) == sorted_eigenValues(j, 0))
            {
                sorted_eigenVectors.col(j) = eigenVectors.col(i);
            }
        }
    }
    
    /*
     pca_dim = 2;
     
     // Extrace pca_dim eigenvectors which represent the original eigenspace adequately
     MatrixXd extract_eigenVectors(dim, pca_dim);
     for (int i = 0; i<pca_dim; i++)
     {
     extract_eigenVectors.col(i) = sorted_eigenVectors.col(i);
     }
     
     // The projection of each cell image can be described as a num * pca_dim matrix
     // Calculated by dataAdjust(data-mean) * extract_eigenVectors
     MatrixXd projection(num, pca_dim);
     projection = dataAdjust * extract_eigenVectors;
     
     MatrixXd OriginalData(num, dim);
     
     OriginalData = (extract_eigenVectors * projection.transpose()).transpose();
     for (int i = 0; i<num; i++) {
     OriginalData.row(i) += mean;
     }
     //cout<<endl<<endl<<endl<<endl;
     //cout<<"Original Data: \n"<<OriginalData<<endl;
     */
    
    ofstream outFile;
    outFile.open(fileName, ios::app);
    
    for (int i = 0; i < sorted_eigenVectors.rows(); i++) {
        for (int j = 0; j < sorted_eigenVectors.cols(); j++) {
            outFile << sorted_eigenVectors(i, j) << " ";
        }
        outFile << endl;
    }
    outFile.close();
    
    
    return sorted_eigenVectors;
    
}

MatrixXd PCA_::BubbleSort(MatrixXd eigenValues, int num) {
    double temp;
    for (int i = 1; i<num; i++) {
        for (int j = num - 1; j >= i; j--) {
            if (eigenValues(j, 0)>eigenValues(j - 1, 0)) {
                temp = eigenValues(j - 1, 0);
                eigenValues(j - 1, 0) = eigenValues(j, 0);
                eigenValues(j, 0) = temp;
            }
        }
    }
    return eigenValues;
}





















