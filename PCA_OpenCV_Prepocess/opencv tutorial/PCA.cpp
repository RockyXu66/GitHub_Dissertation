

//
//  PCA.cpp
//  PCA_opengl_0522
//
//  Created by Yinghan Xu on 22/05/2017.
//  Copyright © 2017 Yinghan Xu. All rights reserved.
//

#include "PCA.hpp"

Mat PCA_::CalculateEigen(vector<Mat> images, int image_num, int width, int height) {
    
    this->num = images.size();
    this->dim = 20 * 20;// width*height;
    this->cellImage_num = (width*height)/(20*20);
    
    this->pca_dim = 100;
    
    int img_dimension = 300;
    
    // Show which image
    int imageIndex = 7;
    
//    Mat img = imread("../../Data_Set/artifix_120x_png/artifix0001.png");
//    imshow("img", img);
    
    time_t currentTime;
    struct tm *localTime;
    string fileName;
    string fileNameRed, fileNameGreen, fileNameBlue;
    string meanFileName;
    string scoresFileName;
    
    int testTime1, testTime2, duration;
    
    time(&currentTime);                   // Get the current time
    localTime = localtime(&currentTime);
    cout<<"Start get each channel values into the data matrix"<<endl;
    // Put each channel values into the data matrix
    MatrixXd data_r(num, width*height);
    MatrixXd data_g(num, width*height);
    MatrixXd data_b(num, width*height);
    for (int i = 0; i<images.size(); i++) {
//        for (int j = 0; j<width*height; j++) {
//            int r = images[i].at<cv::Vec3b>(y,x)[0][j * 3];
//            int g = images[i][j * 3 + 1];
//            int b = images[i][j * 3 + 2];
//            data_r(i, j) = r;
//            data_g(i, j) = g;
//            data_b(i, j) = b;
//        }
        for(int x=0; x<width; x++){
            for(int y=0; y<height; y++){
                int r = images[i].at<Vec3b>(y,x)[2];
                int g = images[i].at<Vec3b>(y,x)[1];
                int b = images[i].at<Vec3b>(y,x)[0];
                data_r(i, y+x*height) = r;
                data_g(i, y+x*height) = g;
                data_b(i, y+x*height) = b;
            }
        }
        cout<<"Matrix "<<i+1<<endl;
    }
    cout<<"Finish get each channel values into the data matrix"<<endl;
    cout<<"Put into data matrix duration time: "<<endl;


////    MatrixXd data_r = MatrixXd::Constant(num, width*height, 1.0);
////    MatrixXd data_g = MatrixXd::Constant(num, width*height, 2.0);
////    MatrixXd data_b = MatrixXd::Constant(num, width*height, 3.0);
    
//    // Calculate mean image
//    unsigned char* meanImage = images[0];
//    meanImage = Mean(data_r, data_g, data_b, meanImage);
//    return meanImage;
    
    cout<<"Start crop to cell"<<endl;
    // Crop to cell => 20*20 pixels
    vector<MatrixXd> data_g_cells, data_b_cells, data_r_cells;
    for (int i = 0; i < height/20; i++) {
        for (int j = 0; j < width/20; j++) {
            MatrixXd data_r_cell(num, dim);
            MatrixXd data_g_cell(num, dim);
            MatrixXd data_b_cell(num, dim);
            for (int k = 0; k < num; k++) {
                
                for (int x = 0; x < 20; x++) {
                    for (int y = 0; y < 20; y++) {
                        int piIndex = j * 20 + i * img_dimension * 20 + y + x*img_dimension;
                        data_r_cell(k, y + 20 * x) = data_r(k, piIndex);
                        data_g_cell(k, y + 20 * x) = data_g(k, piIndex);
                        data_b_cell(k, y + 20 * x) = data_b(k, piIndex);
                    }
                }
            }
            
            data_r_cells.push_back(data_r_cell);
            data_g_cells.push_back(data_g_cell);
            data_b_cells.push_back(data_b_cell);
        }
    }
    cout<<"Finish crop to cell"<<endl;
    cout<<"Crop to cell duration time: "<<endl;
    
    Mat eigenTexture = images[3];
    
    
    // Calculate Eigen values and Eigen vectors
    fileName = "EigenVector_artifix/EigenVector_Red.txt";
    ofstream outFile1(fileName);
    for (int i = 0; i < (width*height)/(20*20); i++) {
        cout << "Cell image " << i << endl;
        Eigen(data_r_cells[i], data_g_cells[i], data_b_cells[i], "Red", fileName);
    }
//
//    fileName = "EigenVector_artifix/EigenVector_Green.txt";
//    ofstream outFile2(fileName);
//    for (int i = 0; i < (width*height) / (20 * 20); i++) {
//        cout << "Cell image " << i << endl;
//        Eigen(data_g_cells[i], "Green", fileName);
//    }
//    
//    fileName = "EigenVector_artifix/EigenVector_Blue.txt";
//    ofstream outFile3(fileName);
//    for (int i = 0; i < (width*height) / (20 * 20); i++) {
//        cout << "Cell image " << i << endl;
//        Eigen(data_b_cells[i], "Blue", fileName);
//    }
//
    
//
//    // Save high-dimensional eigenspace
//    cout<<"Start saving 150 high-dimensional eigenspace"<<endl;
//    fileName = "EigenVector_artifix/EigenVector_Red.txt";
//    SaveHighEigen("EigenVector_artifix/highEigen"+to_string(pca_dim)+"_Red.txt", fileName, pca_dim);
//    fileName = "EigenVeEigenVector_artifixctor/EigenVector_Green.txt";
//    SaveHighEigen("EigenVector_artifix/highEigen"+to_string(pca_dim)+"_Green.txt", fileName, pca_dim);
//    fileName = "EigenVector_artifix/EigenVector_Blue.txt";
//    SaveHighEigen("EigenVector_artifix/highEigen"+to_string(pca_dim)+"_Blue.txt", fileName, pca_dim);
//    cout<<"Finish saving 150 high-dimensional eigenspace"<<endl;
//    
//    // Calculate mean image and save to files
//    cout<<"Start calculating mean image and saving to files"<<endl;
//    fileName = "EigenVector_artifix/mean_Red.txt";
//    Mean(data_r_cells, fileName);
//    fileName = "EigenVector_artifix/mean_Green.txt";
//    Mean(data_g_cells, fileName);
//    fileName = "EigenVector_artifix/mean_Blue.txt";
//    Mean(data_b_cells, fileName);
//    cout<<"Finish calculating mean image and saving to files"<<endl;
//    
////    MatrixXd EigenVectors(cellImage_num*dim, 150);
////    MatrixXd means(cellImage_num, dim);
//    
//    
//    // Calculate scores for eigentextures and save to files
//    cout<<"Start calculating scores for eigentextures and saving to files"<<endl;
//    fileName = "EigenVector_artifix/highEigen150_Red.txt";
//    meanFileName = "EigenVector_artifix/mean_Red.txt";
//    scoresFileName = "EigenVector_artifix/scores150_Red.txt";
//    CalScores(data_r_cells, fileName, meanFileName, scoresFileName);
//    fileName = "EigenVector_artifix/highEigen150_Green.txt";
//    meanFileName = "EigenVector_artifix/mean_Green.txt";
//    scoresFileName = "EigenVector_artifix/scores150_Green.txt";
//    CalScores(data_g_cells, fileName, meanFileName, scoresFileName);
//    fileName = "EigenVector_artifix/highEigen150_Blue.txt";
//    meanFileName = "EigenVector_artifix/mean_Blue.txt";
//    scoresFileName = "EigenVector_artifix/scores150_Blue.txt";
//    CalScores(data_b_cells, fileName, meanFileName, scoresFileName);
    /*
//    cout<<"Start calculating scores for eigentextures and saving to files"<<endl;
//    fileName = "EigenVector/highEigen150_Red.txt";
//    meanFileName = "EigenVector/mean_Red.txt";
//    scoresFileName = "EigenVector/scores150_Red.txt";
////    Cal(fileName, meanFileName, scoresFileName);
////    CalScores(data_r_cells, fileName, meanFileName, scoresFileName);
//    
//    // Get Eigen textures
//    cout<<"Start getting eigen textures"<<endl;
//    ifstream myfile_red(fileName);
//    for (int i = 0; i < cellImage_num*dim; i++) {
//        for (int j = 0; j < 150; j++) {
//            myfile_red >> EigenVectors(i, j);
//        }
//    }
//    myfile_red.close();
//    cout<<"Finish getting eigen textures"<<endl;
//    // Get mean image
//    cout<<"Start getting mean image"<<endl;
//    ifstream meanfile_red(meanFileName);
//    for (int i = 0; i < cellImage_num; i++) {
//        for (int j = 0; j < dim; j++) {
//            meanfile_red >> means(i, j);
//        }
//    }
//    meanfile_red.close();
//    cout<<"Finish getting mean image"<<endl;
//    ofstream scoresFile_red(scoresFileName);
//    cout<<"Start calculating"<<endl;
//    for(int cellIndex=0; cellIndex<cellImage_num; cellIndex++){
//        MatrixXd data = data_r_cells[cellIndex];
//        
//        MatrixXd mean = means.row(cellIndex);   //data.colwise().mean();  // mean  mean =
//        MatrixXd dataAdjust(num, dim);  // data-mean
//        for (int i = 0; i<num; i++) {
//            dataAdjust.row(i) = data.row(i) - mean;  // dataAdjust = data-mean
//        }
//        
//        MatrixXd oneVec(dim, 150);
//        for (int i = 0; i < dim; i++) {
//            for (int j = 0; j < 150; j++) {
//                oneVec(i, j) = EigenVectors(cellIndex*dim + i, j);
//            }
//        }
//        
//        // Calculate scores and save to the file
//        for(int k=0; k<num; k++){
//            for(int i=0; i<150; i++){
//                double projectionValue = dataAdjust.row(k) * oneVec.col(i);
//                scoresFile_red << projectionValue << " ";
//            }
//            scoresFile_red <<endl;
//        }
//        cout<<cellIndex+1<<" ";
//    }
//    cout<<endl;
//    scoresFile_red.close();
//    
//    fileName = "EigenVector/highEigen150_Green.txt";
//    meanFileName = "EigenVector/mean_Green.txt";
//    scoresFileName = "EigenVector/scores150_Green.txt";
//    // Get Eigen textures
//    cout<<"Start getting eigen textures"<<endl;
//    ifstream myfile_green(fileName);
//    for (int i = 0; i < cellImage_num*dim; i++) {
//        for (int j = 0; j < 150; j++) {
//            myfile_green >> EigenVectors(i, j);
//        }
//    }
//    myfile_green.close();
//    cout<<"Finish getting eigen textures"<<endl;
//    // Get mean image
//    cout<<"Start getting mean image"<<endl;
//    ifstream meanfile_green(meanFileName);
//    for (int i = 0; i < cellImage_num; i++) {
//        for (int j = 0; j < dim; j++) {
//            meanfile_green >> means(i, j);
//        }
//    }
//    meanfile_green.close();
//    cout<<"Finish getting mean image"<<endl;
//    ofstream scoresFile_green(scoresFileName);
//    cout<<"Start calculating"<<endl;
//    for(int cellIndex=0; cellIndex<cellImage_num; cellIndex++){
//        MatrixXd data = data_g_cells[cellIndex];
//        
//        MatrixXd mean = means.row(cellIndex);   //data.colwise().mean();  // mean  mean =
//        MatrixXd dataAdjust(num, dim);  // data-mean
//        for (int i = 0; i<num; i++) {
//            dataAdjust.row(i) = data.row(i) - mean;  // dataAdjust = data-mean
//        }
//        
//        MatrixXd oneVec(dim, 150);
//        for (int i = 0; i < dim; i++) {
//            for (int j = 0; j < 150; j++) {
//                oneVec(i, j) = EigenVectors(cellIndex*dim + i, j);
//            }
//        }
//        
//        // Calculate scores and save to the file
//        for(int k=0; k<num; k++){
//            for(int i=0; i<150; i++){
//                double projectionValue = dataAdjust.row(k) * oneVec.col(i);
//                scoresFile_green << projectionValue << " ";
//            }
//            scoresFile_green <<endl;
//        }
//        cout<<cellIndex+1<<" ";
//    }
//    cout<<endl;
//    scoresFile_green.close();
////    CalScores(data_g_cells, fileName, meanFileName, scoresFileName);

    
//    fileName = "EigenVector/highEigen150_Blue.txt";
//    meanFileName = "EigenVector/mean_Blue.txt";
//    scoresFileName = "EigenVector/scores150_Blue.txt";
//    // Get Eigen textures
//    cout<<"Start getting eigen textures"<<endl;
//    ifstream myfile_blue(fileName);
//    for (int i = 0; i < cellImage_num*dim; i++) {
//        for (int j = 0; j < 150; j++) {
//            myfile_blue >> EigenVectors(i, j);
//        }
//    }
//    myfile_blue.close();
//    cout<<"Finish getting eigen textures"<<endl;
//    // Get mean image
//    cout<<"Start getting mean image"<<endl;
//    ifstream meanfile_blue(meanFileName);
//    for (int i = 0; i < cellImage_num; i++) {
//        for (int j = 0; j < dim; j++) {
//            meanfile_blue >> means(i, j);
//        }
//    }
//    meanfile_blue.close();
//    cout<<"Finish getting mean image"<<endl;
//    ofstream scoresFile_blue(scoresFileName);
//    cout<<"Start calculating"<<endl;
//    for(int cellIndex=0; cellIndex<cellImage_num; cellIndex++){
//        MatrixXd data = data_b_cells[cellIndex];
//        
//        MatrixXd mean = means.row(cellIndex);   //data.colwise().mean();  // mean  mean =
//        MatrixXd dataAdjust(num, dim);  // data-mean
//        for (int i = 0; i<num; i++) {
//            dataAdjust.row(i) = data.row(i) - mean;  // dataAdjust = data-mean
//        }
//        
//        MatrixXd oneVec(dim, 150);
//        for (int i = 0; i < dim; i++) {
//            for (int j = 0; j < 150; j++) {
//                oneVec(i, j) = EigenVectors(cellIndex*dim + i, j);
//            }
//        }
//        
//        // Calculate scores and save to the file
//        for(int k=0; k<num; k++){
//            for(int i=0; i<150; i++){
//                double projectionValue = dataAdjust.row(k) * oneVec.col(i);
//                scoresFile_blue << projectionValue << " ";
//            }
//            scoresFile_blue <<endl;
//        }
//        cout<<cellIndex+1<<" ";
//    }
//    cout<<endl;
//    scoresFile_blue.close();
////    CalScores(data_b_cells, fileName, meanFileName, scoresFileName);
//
     */
    cout<<"Finish calculating scores for eigentextures and saving to files"<<endl;
    
    // Reconstruction the whole image
    /*
     fileName = "EigenVector/highEigen100_Red.txt";
    meanFileName = "EigenVector/mean_Red.txt";
    scoresFileName = "EigenVector/scores150_Red.txt";
    cout<<"Red"<<endl;
    vector<MatrixXd> ReconstructedData_Red = WholeReconstruction(fileName, meanFileName, scoresFileName, imageIndex);
    fileName = "EigenVector/highEigen100_Green.txt";
    meanFileName = "EigenVector/mean_Green.txt";
    scoresFileName = "EigenVector/scores150_Green.txt";
    cout<<"Green"<<endl;
    vector<MatrixXd> ReconstructedData_Green = WholeReconstruction(fileName, meanFileName, scoresFileName, imageIndex);
    fileName = "EigenVector/highEigen100_Blue.txt";
    meanFileName = "EigenVector/mean_Blue.txt";
    scoresFileName = "EigenVector/scores150_Blue.txt";
    cout<<"Blue"<<endl;
    vector<MatrixXd> ReconstructedData_Blue = WholeReconstruction(fileName, meanFileName, scoresFileName, imageIndex);
    
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
                    int piIndex = j * 20 + i * 1080 * 20 + y + x*1080;
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
     */
    
//    return images[0];
    //return meanImage;
    return eigenTexture;
}

static Mat formatImagesForPCA(const vector<Mat> &data)
{
    Mat dst(static_cast<int>(data.size()), data[0].rows*data[0].cols, CV_32F);
    for(unsigned int i = 0; i < data.size(); i++)
    {
        Mat image_row = data[i].clone().reshape(1,1);
        Mat row_i = dst.row(i);
        image_row.convertTo(row_i,CV_32F);
    }
    return dst;
}

unsigned char* PCA_::Mean(MatrixXd data_r, MatrixXd data_g, MatrixXd data_b, unsigned char* meanImage) {
    MatrixXd mean_r = data_r.colwise().mean();
    MatrixXd mean_g = data_g.colwise().mean();
    MatrixXd mean_b = data_b.colwise().mean();
    
    for (int i = 0; i<dim*cellImage_num; i++) {
        meanImage[i * 3] = mean_r(0, i);
        meanImage[i * 3 + 1] = mean_g(0, i);
        meanImage[i * 3 + 2] = mean_b(0, i);
    }
    return meanImage;
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
        cout<<i+1<<" ";
    }
    
    ofstream newFile(fileName);
    for(int i=0; i<cellImage_num; i++){
        for(int j=0; j<dim; j++){
            newFile << means[i](0, j)<<" ";
        }
        newFile << endl;
        cout<<i+1<<" ";
    }
    cout<<endl;
    newFile.close();
}

void PCA_::CalScores(vector<MatrixXd> whole_data, string fileName, string meanFileName, string newFileName){
    int score_dim = 150;
    // Get Eigen textures
    cout<<"Start getting eigen textures"<<endl;
    MatrixXd EigenVectors(cellImage_num*dim, score_dim);
    ifstream myfile(fileName);
    for (int i = 0; i < cellImage_num*dim; i++) {
        for (int j = 0; j < score_dim; j++) {
            myfile >> EigenVectors(i, j);
        }
    }
    myfile.close();
    cout<<"Finish getting eigen textures"<<endl;
    
    // Get mean image
    cout<<"Start getting mean image"<<endl;
    MatrixXd means(cellImage_num, dim);
    ifstream meanfile(meanFileName);
    for (int i = 0; i < cellImage_num; i++) {
        for (int j = 0; j < dim; j++) {
            meanfile >> means(i, j);
        }
    }
    meanfile.close();
    cout<<"Finish getting mean image"<<endl;
    
    ofstream scoresFile(newFileName);
    
    cout<<"Start calculating"<<endl;
    for(int cellIndex=0; cellIndex<cellImage_num; cellIndex++){
        MatrixXd data = whole_data[cellIndex];
        
        MatrixXd mean = means.row(cellIndex);   //data.colwise().mean();  // mean  mean =
        MatrixXd dataAdjust(num, dim);  // data-mean
        cout<<"1"<<endl;
        for (int i = 0; i<num; i++) {
            dataAdjust.row(i) = data.row(i) - mean;  // dataAdjust = data-mean
        }
        cout<<"Finish getting dataAdjust"<<endl;
        
        MatrixXd oneVec(dim, score_dim);
        for (int i = 0; i < dim; i++) {
            cout<<"2";
            for (int j = 0; j < score_dim; j++) {
                oneVec(i, j) = EigenVectors(cellIndex*dim + i, j);
            }
        }
        
        // Calculate scores and save to the file
        cout<<endl;
        cout<<"Start saving"<<endl;
        for(int k=0; k<num; k++){
            for(int i=0; i<score_dim; i++){
                double projectionValue = dataAdjust.row(k) * oneVec.col(i);
                scoresFile << projectionValue << " ";
            }
            scoresFile <<endl;
        }
        cout<<cellIndex+1<<" ";
    }
    cout<<endl;
    
    scoresFile.close();
}

void PCA_::Cal(string fileName, string meanFileName, string newFileName){
    
    // Get Eigen textures
    cout<<"1"<<endl;
    MatrixXd EigenVectors(cellImage_num*dim, 100);
    ifstream myfile(fileName);
    for (int i = 0; i < cellImage_num*dim; i++) {
        for (int j = 0; j < 100; j++) {
            myfile >> EigenVectors(i, j);
        }
    }
    myfile.close();
    cout<<"2"<<endl;
    
    // Get mean image
    cout<<"3"<<endl;
    MatrixXd means(cellImage_num, dim);
    ifstream meanfile(meanFileName);
    for (int i = 0; i < cellImage_num; i++) {
        for (int j = 0; j < dim; j++) {
            meanfile >> means(i, j);
        }
    }
    meanfile.close();
    cout<<"4"<<endl;
    
    ofstream scoresFile(newFileName);
    
    for(int cellIndex=0; cellIndex<cellImage_num; cellIndex++){
//        MatrixXd data = whole_data[cellIndex];
        
        MatrixXd mean = means.row(cellIndex);   //data.colwise().mean();  // mean  mean =
//        MatrixXd dataAdjust(num, dim);  // data-mean
//        for (int i = 0; i<num; i++) {
//            dataAdjust.row(i) = data.row(i) - mean;  // dataAdjust = data-mean
//        }
        
        MatrixXd oneVec(dim, 100);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < 100; j++) {
                oneVec(i, j) = EigenVectors(cellIndex*dim + i, j);
            }
        }
        
        // Calculate scores and save to the file
//        for(int k=0; k<num; k++){
//            for(int i=0; i<100; i++){
//                double projectionValue = dataAdjust.row(k) * oneVec.col(i);
//                scoresFile << projectionValue << " ";
//            }
//            scoresFile <<endl;
//        }
        cout<<cellIndex+1<<" ";
    }
    cout<<endl;
    
    scoresFile.close();
}


vector<MatrixXd> PCA_::WholeReconstruction(string fileName, string meanFileName, string scoresFileName, int imageIndex){
    cout<<"Reconstruction data: "<<endl;
    
    // Get Eigen textures
    cout<<"Start getting eigen textures"<<endl;
    MatrixXd EigenVectors(cellImage_num*dim, pca_dim);
    ifstream myfile(fileName);
    for (int i = 0; i < cellImage_num*dim; i++) {
        for (int j = 0; j < pca_dim; j++) {
            myfile >> EigenVectors(i, j);
        }
    }
    myfile.close();
    cout<<"Finish getting eigen textures"<<endl;
    
    vector<MatrixXd> whole_OriginalData;
    
    // Get mean image
    cout<<"Start getting mean image"<<endl;
    MatrixXd means(cellImage_num, dim);
    ifstream meanfile(meanFileName);
    for (int i = 0; i < cellImage_num; i++) {
        for (int j = 0; j < dim; j++) {
            meanfile >> means(i, j);
        }
    }
    meanfile.close();
    cout<<"Finish getting mean image"<<endl;
    
    // Get scores
    cout<<"Start getting scores"<<endl;
    MatrixXd whole_scores(cellImage_num*num, pca_dim);
    ifstream scoresFile(scoresFileName);
    for (int i = 0; i < cellImage_num*num; i++) {
        for (int j = 0; j < pca_dim; j++) {
            scoresFile >> whole_scores(i, j);
        }
    }
    scoresFile.close();
    cout<<"Finish getting scores"<<endl;

    MatrixXd oneImage_scores(cellImage_num, pca_dim);
    for(int i=0; i<cellImage_num; i++){
        for(int j=0; j<pca_dim; j++){
            oneImage_scores(i, j) = whole_scores(imageIndex+i*num, j);
        }
    }
    
    cout<<"Cell: "<<endl;
    for(int cellIndex=0; cellIndex<cellImage_num; cellIndex++){
        MatrixXd mean = means.row(cellIndex);   //data.colwise().mean();  // mean  mean =
        
        MatrixXd oneVec(dim, pca_dim);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < pca_dim; j++) {
                oneVec(i, j) = EigenVectors(cellIndex*dim + i, j);
            }
        }

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
    cout<<endl;
    
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

void PCA_::Eigen(MatrixXd r_data, MatrixXd g_data, MatrixXd b_data, string name, string fileName) {
    cout << name << endl;
    //cout<<data<<endl;
//    MatrixXd mean = data.colwise().mean();
//    
//    MatrixXd dataAdjust(num, dim);  // data-mean
//    MatrixXd dataAdjust_Transpose(dim, num);// transpose of dataAdjust
//    
//    for (int i = 0; i<num; i++) {
//        dataAdjust.row(i) = data.row(i) - mean;  // dataAdjust = data-mean
//    }
//    for(int i=0; i<num; i++){
//        for(int j=0; j<dim; j++){
//            if(dataAdjust(i, j)<1){
//                dataAdjust(i, j) = 0;
//            }
//        }
//    }
//    dataAdjust_Transpose = dataAdjust.transpose();
//    
//    // Covariance matrix
//    MatrixXd Covariance = (1 / (double)num) * dataAdjust_Transpose * dataAdjust;
//    
//    cout<<"Start eigen calculation"<<endl;
//    EigenSolver<MatrixXd> eigenSolver(Covariance);
//    
//    // Calculate eigen value and eigen vector from covariance matrix
//    MatrixXd eigenValues = eigenSolver.eigenvalues().real();
//    MatrixXd eigenVectors = eigenSolver.eigenvectors().real();
//    
//    // Sort eigen values
//    MatrixXd sorted_eigenValues = BubbleSort(eigenValues, dim);
//    
//    MatrixXd sorted_eigenVectors(dim, dim);
//    for (int i = 0; i<dim; i++){
//        for (int j = 0; j<dim; j++){
//            if (eigenValues(i, 0) == sorted_eigenValues(j, 0)){
//                sorted_eigenVectors.col(j) = eigenVectors.col(i);
//            }
//        }
//    }
    
    
//    vector<Mat> images;
//    for(int x=0; x<120; x++){
//        Mat image(20, 20, CV_32F);
//        for(int i=0; i<20; i++){
//            for(int j=0; j<20; j++){
//                Vec3b color = Vec3b(b_data(x, j+i*20), g_data(x, j+i*20), r_data(x, j+i*20));
//                image.at<Vec3b>(Point(j,i)) = color;
//            }
//        }
//        images.push_back(image);
//    }
//    
//    Mat data = formatImagesForPCA(images);
    
    Mat mat(r_data.rows(), r_data.cols(), CV_64FC1);
    for(int i=0; i<r_data.rows(); i++){
        for(int j=0; j<r_data.cols(); j++){
            mat.at<double>(i, j) = r_data(i, j);
        }
    }
    // Perform a PCA:
    PCA pca(mat, Mat(), CV_PCA_DATA_AS_ROW);
    Mat eigenvalues = pca.eigenvalues.clone();
    Mat eigenvectors = pca.eigenvectors.clone();
    
//    for (int i = 0; i < 6; ++i)
//        std::cout << eigenvalues.at<float>(i, 0) << std::endl;
    
    ofstream outFile;
    outFile.open(fileName, ios::app);
    
    for (int i = 0; i < eigenvectors.rows; i++) {
        for (int j = 0; j < eigenvectors.cols; j++) {
            outFile << eigenvectors.at<float>(i, j) << " ";
        }
        outFile << endl;
    }
    outFile.close();
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





















