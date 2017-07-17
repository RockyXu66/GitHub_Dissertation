/*
 //
//  PCA.cpp
//  PCA_opengl_0522
//
//  Created by Yinghan Xu on 22/05/2017.
//  Copyright © 2017 Yinghan Xu. All rights reserved.
//

#include "PCA.hpp"

unsigned char* PCA::CalculateEigen(vector<unsigned char*> images, int width, int height){
    
    this->num = images.size();
    this->dim = 20*20;
    
    // Put each channel values into the data matrix
    MatrixXd data_r(num, dim);
    MatrixXd data_g(num, dim);
    MatrixXd data_b(num, dim);
    
    for(int i=0; i<images.size(); i++){
        unsigned char* image = images[i];
        for(int j=0; j<width*height; j++){
            int r = image[j*3];
            int g = image[j*3+1];
            int b = image[j*3+2];
            
            data_r(i, j) = r;
            data_g(i, j) = g;
            data_b(i, j) = b;
        }
    }
    
    // Crop to cell => 20*20 pixels
    vector<MatrixXd> data_r_cells, data_g_cells, data_b_cells;
    for (int i = 0; i < height/20; i++) {
        for (int j = 0; j < width / 20; j++) {
            MatrixXd data_r_cell(num, dim);
            MatrixXd data_g_cell(num, dim);
            MatrixXd data_b_cell(num, dim);
            for (int k = 0; k < num; k++) {
                
                for (int x = 0; x < 20; x++) {
                    for (int y = 0; y < 20; y++) {
                        int piIndex = j * 20 + i * 300 * 20 + y + x*300;
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
    
    // Calculate mean image
//    unsigned char* meanImage = images[1];
//    Mean(data_r, data_g, data_b, meanImage);
    
    // Calculate Eigen values and Eigen vectors
    unsigned char* newImage = images[0];
    time_t currentTime;
    struct tm *localTime;
    string fileName;

//    time(&currentTime);                   // Get the current time
//    localTime = localtime(&currentTime);
//    fileName = "EigenVector_Red_" + to_string(localTime->tm_hour)
//    + to_string(localTime->tm_min) + to_string(localTime->tm_sec) + ".txt";
//    ofstream outFile1(fileName);
//    MatrixXd eigenVectors_r = Eigen(data_r, newImage, "Red", fileName);

//    time(&currentTime);                   // Get the current time
//    localTime = localtime(&currentTime);
//    fileName = "EigenVector_Green_" + to_string(localTime->tm_hour)
//    + to_string(localTime->tm_min) + to_string(localTime->tm_sec) + ".txt";
//    ofstream outFile2(fileName);
//    MatrixXd eigenVectors_g = Eigen(data_g, newImage, "Green", fileName);
    time(&currentTime);                   // Get the current time
    localTime = localtime(&currentTime);
    fileName = "EigenVector_Blue_" + to_string(localTime->tm_hour)
    + to_string(localTime->tm_min) + to_string(localTime->tm_sec) + ".txt";
    ofstream outFile3(fileName);
    for (int i = 0; i < (width*height) / (20 * 20); i++) {
        cout << "Cell image " << i << endl;
        MatrixXd eigenVectors_b = Eigen(data_b_cells[i], newImage, "Blue", fileName);
    }
    
    // Reconstructoin
//    fileName = "EigenVector_Red_233738.txt";
//    MatrixXd eigenTexture_Red = Reconstruction(data_r, newImage, fileName);
//    
//    fileName = "EigenVector_Green_233754.txt";
//    MatrixXd eigenTexture_Green = Reconstruction(data_g, newImage, fileName);
//    
//    fileName = "EigenVector_Blue_23388.txt";
//    MatrixXd eigenTexture_Blue = Reconstruction(data_b, newImage, fileName);
//    
//    unsigned char* eigenTexture = images[0];
//    for(int i=0; i<dim; i++){
//        eigenTexture[i*3] = eigenTexture_Red(0, i);
//        eigenTexture[i*3+1] = eigenTexture_Green(0, i);
//        eigenTexture[i*3+2] = eigenTexture_Blue(0, i);
//    }
    
//    fileName = "EigenVector_Red_cellTest.txt";
//    MatrixXd original_Red = Reconstruction(data_r, newImage, fileName);
//    fileName = "EigenVector_Green_cellTest.txt";
//    MatrixXd original_Green = Reconstruction(data_g, newImage, fileName);
//    fileName = "EigenVector_Blue_cellTest.txt";
//    MatrixXd original_Blue = Reconstruction(data_b, newImage, fileName);
//    unsigned char* originalTexture = images[0];
//    for(int i=0; i<dim; i++){
//        originalTexture[i*3] = original_Red(0, i);
//        originalTexture[i*3+1] = original_Green(0, i);
//        originalTexture[i*3+2] = original_Blue(0, i);
//    }
    
//    fileName = "EigenVector_Red_outTest.txt";
//    MatrixXd eigenTexture_Red = Reconstruction(data_r, newImage, fileName);
    
    
    return images[2];
//    return meanImage;
//    return eigenTexture;
//    return originalTexture;
}

void PCA::Mean(MatrixXd data_r, MatrixXd data_g, MatrixXd data_b, unsigned char* meanImage){
    MatrixXd mean_r =data_r.colwise().mean();
    MatrixXd mean_g =data_g.colwise().mean();
    MatrixXd mean_b =data_b.colwise().mean();
    
    for(int i=0; i<dim; i++){
        meanImage[i*3] = mean_r(0, i);
        meanImage[i*3+1] = mean_g(0, i);
        meanImage[i*3+2] = mean_b(0, i);
    }
}

MatrixXd PCA::Reconstruction(MatrixXd data, unsigned char* newImage, string fileName) {
    cout << data.row(0) << endl;
    MatrixXd mean = data.colwise().mean();
    
    MatrixXd dataAdjust(num, dim);  // data-mean
    for (int i = 0; i<num; i++) {
        dataAdjust.row(i) = data.row(i) - mean;  // dataAdjust = data-mean
    }
    
    MatrixXd EigenVectors(5*dim, dim);
    ifstream myfile(fileName);
    for (int i = 0; i < 5*dim; i++) {
        for (int j = 0; j < dim; j++) {
            myfile >> EigenVectors(i, j);
        }
    }
    myfile.close();
    
    MatrixXd oneVec(dim, dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            oneVec(i, j) = EigenVectors(i, j);
        }
    }
    
//    cout<<oneVec<<endl;
    
    pca_dim = 40;
    
    int imageIndex = 1;
    
    MatrixXd OriginalData = MatrixXd::Constant(1, dim, 0.0);
    for(int i=0; i<pca_dim; i++){
        double projectionValue = dataAdjust.row(imageIndex) * oneVec.col(i);
        MatrixXd eigenTexture(1, dim);
        eigenTexture = projectionValue * oneVec.col(i).transpose();
        OriginalData += eigenTexture;
    }
//    double projectionValue1 = dataAdjust.row(0) * oneVec.col(0);
//    MatrixXd eigenTexture1(1, dim);
//    eigenTexture1 = projectionValue1 * oneVec.col(0).transpose();
//    
//    double projectionValue2 = dataAdjust.row(0) * oneVec.col(1);
//    MatrixXd eigenTexture2(1, dim);
//    eigenTexture2 = projectionValue2 * oneVec.col(1).transpose();
    
    OriginalData += mean;
    cout<<"Reconstruction Data: "<<endl;
    cout<<OriginalData<<endl;
    
    
//    // Extrace pca_dim eigenvectors which represent the original eigenspace adequately
//    MatrixXd extract_eigenVectors(dim, pca_dim);
//    for (int i = 0; i<pca_dim; i++){
//        extract_eigenVectors.col(i) = oneVec.col(i);
//    }
//    
//    MatrixXd projection(num, pca_dim);
//    projection = dataAdjust * extract_eigenVectors;
//    
////    double projectionTest = dataAdjust.row(0) * extract_eigenVectors.col(0);
//    
//    MatrixXd All_OriginalData(num, dim);
//    
//    All_OriginalData = (extract_eigenVectors * projection.transpose()).transpose();
//    for (int i = 0; i<num; i++) {
//        All_OriginalData.row(i) += mean;
//    }
//    
//    cout << endl;
//    cout<<"Reconstruction Data:"<<endl;
//    cout<<All_OriginalData.row(imageIndex)<<endl;
//    return EigenVectors.col(0).transpose()*400;
    return OriginalData;
}

MatrixXd PCA::Eigen(MatrixXd data, unsigned char* newImage, string name, string fileName){
    cout<<endl<<endl;
    cout<<name<<endl;
    //cout<<data<<endl;
    MatrixXd mean =data.colwise().mean();
    
    MatrixXd dataAdjust(num, dim);  // data-mean
    MatrixXd dataAdjust_Transpose(dim, num);// transpose of dataAdjust
    
    for(int i=0; i<num; i++){
        dataAdjust.row(i) =data.row(i)-mean;  // dataAdjust = data-mean
    }
    dataAdjust_Transpose = dataAdjust.transpose();
    
    // Covariance matrix
    MatrixXd Covariance = (1/(double)num) * dataAdjust_Transpose * dataAdjust;
    
    EigenSolver<MatrixXd> eigenSolver(Covariance);
    
    // Calculate eigen value and eigen vector from covariance matrix
    MatrixXd eigenValues =eigenSolver.eigenvalues().real();
    MatrixXd eigenVectors =eigenSolver.eigenvectors().real();
    
    // Sort eigen values
    MatrixXd sorted_eigenValues = BubbleSort(eigenValues, dim);
    
    MatrixXd sorted_eigenVectors(dim, dim);
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            if(eigenValues(i, 0) == sorted_eigenValues(j, 0))
            {
                sorted_eigenVectors.col(j) =eigenVectors.col(i);
            }
        }
    }
 
//    pca_dim = 2;
//    
//    // Extrace pca_dim eigenvectors which represent the original eigenspace adequately
//    MatrixXd extract_eigenVectors(dim, pca_dim);
//    for(int i=0; i<pca_dim; i++)
//    {
//        extract_eigenVectors.col(i) =sorted_eigenVectors.col(i);
//    }
//    
//    // The projection of each cell image can be described as a num * pca_dim matrix
//    // Calculated by dataAdjust(data-mean) * extract_eigenVectors
//    MatrixXd projection(num, pca_dim);
//    projection = dataAdjust * extract_eigenVectors;
//    
//    MatrixXd OriginalData(num, dim);
//    
//    OriginalData = (extract_eigenVectors * projection.transpose()).transpose();
//    for(int i=0; i<num; i++){
//        OriginalData.row(i) += mean;
//    }
    //cout<<endl<<endl<<endl<<endl;
    //cout<<"Original Data: \n"<<OriginalData<<endl;
    
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



MatrixXd PCA::BubbleSort(MatrixXd eigenValues, int num){
    double temp;
    for(int i=1; i<num; i++){
        for(int j=num-1; j>=i; j--) {
            if(eigenValues(j, 0)>eigenValues(j-1, 0)){
                temp      = eigenValues(j-1, 0);
                eigenValues(j-1, 0) = eigenValues(j, 0);
                eigenValues(j, 0)   = temp;
            }
        }
    }
    return eigenValues;
}
*/



//
//  PCA.cpp
//  PCA_opengl_0522
//
//  Created by Yinghan Xu on 22/05/2017.
//  Copyright © 2017 Yinghan Xu. All rights reserved.
//

#include "PCA.hpp"

unsigned char* PCA::CalculateEigen(vector<unsigned char*> images, int width, int height) {
    
    this->num = images.size();
    this->dim = 20 * 20;// width*height;
    
    // Put each channel values into the data matrix
    MatrixXd data_r(num, width*height);
    MatrixXd data_g(num, width*height);
    MatrixXd data_b(num, width*height);
    
    for (int i = 0; i<images.size(); i++) {
        unsigned char* image = images[i];
        for (int j = 0; j<width*height; j++) {
            int r = image[j * 3];
            int g = image[j * 3 + 1];
            int b = image[j * 3 + 2];
            data_r(i, j) = r;
            data_g(i, j) = g;
            data_b(i, j) = b;
        }
    }
    
    // Crop to cell => 20*20 pixels
    vector<MatrixXd> data_r_cells, data_g_cells, data_b_cells;
    for (int i = 0; i < height/20; i++) {
        for (int j = 0; j < width / 20; j++) {
            MatrixXd data_r_cell(num, dim);
            MatrixXd data_g_cell(num, dim);
            MatrixXd data_b_cell(num, dim);
            for (int k = 0; k < num; k++) {
                
                for (int x = 0; x < 20; x++) {
                    for (int y = 0; y < 20; y++) {
                        int piIndex = j * 20 + i * 300 * 20 + y + x*300;
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
    
    
    
    
    // Calculate mean image
    //unsigned char* meanImage = images[0];
    //Mean(data_r, data_g, data_b, meanImage);
    
    // Calculate Eigen values and Eigen vectors
    unsigned char* newImage = images[0];
    time_t currentTime;
    struct tm *localTime;
    string fileName;
    
    // Copy one file to another file
    ifstream infile("eigenVectors_cell copy/EigenVector_Blue_180-224.txt");
    ofstream outfile("eigenVectors_cell copy/EigenVector_Blue_0-150.txt", ios::app);
    string content = "";
    int i;
    
    for(i=0 ; infile.eof()!=true ; i++) // get content of infile
        content += infile.get();
    
    i--;
    content.erase(content.end()-1);     // erase last character
    
    cout << i << " characters read...\n";
    infile.close();
    
    outfile << content;                 // output
    outfile.close();
    
    
//    time(&currentTime);                   // Get the current time
//    localTime = localtime(&currentTime);
//    fileName = "EigenVector/EigenVector_Red_" + to_string(localTime->tm_hour)
//    + to_string(localTime->tm_min) + to_string(localTime->tm_sec) + ".txt";
//    ofstream outFile1(fileName);
//    MatrixXd eigenVectors_r = Eigen(data_r_cells[0], newImage, "Red", fileName);
    
//    for (int i = 180; i < (width*height)/(20*20); i++) {
//        cout << "Cell image " << i << endl;
//        MatrixXd eigenVectors_r = Eigen(data_r_cells[i], newImage, "Red", fileName);
//    }

    
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
//    for (int i = 180; i < (width*height) / (20 * 20); i++) {
//        cout << "Cell image " << i << endl;
//        MatrixXd eigenVectors_b = Eigen(data_b_cells[i], newImage, "Blue", fileName);
//    }
    
    // Reconstructoin
    /*fileName = "EigenVector_Red_233738.txt";
     MatrixXd OriginalData_Red = Reconstruction(data_r, fileName);
     fileName = "EigenVector_Green_233754.txt";
     MatrixXd OriginalData_Green = Reconstruction(data_g, fileName);
     fileName = "EigenVector_Blue_23388.txt";
     MatrixXd OriginalData_Blue = Reconstruction(data_b, fileName);
     unsigned char* eigenTexture = images[0];
     for (int i = 0; i<dim; i++) {
     eigenTexture[i * 3] = OriginalData_Red(0, i);
     eigenTexture[i * 3 + 1] = OriginalData_Green(0, i);
     eigenTexture[i * 3 + 2] = OriginalData_Blue(0, i);
     }*/
    
    return images[0];
    //return meanImage;
    //return eigenTexture;
}

void PCA::Mean(MatrixXd data_r, MatrixXd data_g, MatrixXd data_b, unsigned char* meanImage) {
    MatrixXd mean_r = data_r.colwise().mean();
    MatrixXd mean_g = data_g.colwise().mean();
    MatrixXd mean_b = data_b.colwise().mean();
    
    for (int i = 0; i<dim; i++) {
        meanImage[i * 3] = mean_r(0, i);
        meanImage[i * 3 + 1] = mean_g(0, i);
        meanImage[i * 3 + 2] = mean_b(0, i);
    }
}

MatrixXd PCA::Reconstruction(MatrixXd data, string fileName) {
    //    cout << data << endl;
    MatrixXd mean = data.colwise().mean();
    
    MatrixXd dataAdjust(num, dim);  // data-mean
    for (int i = 0; i<num; i++) {
        dataAdjust.row(i) = data.row(i) - mean;  // dataAdjust = data-mean
    }
    
    MatrixXd EigenVectors(dim, dim);
    ifstream myfile(fileName);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            myfile >> EigenVectors(i, j);
        }
    }
    myfile.close();
    
    double projectionValue1 = dataAdjust.row(1) * EigenVectors.col(0);
    MatrixXd eigenTexture1(1, dim);
    eigenTexture1 = projectionValue1 * EigenVectors.col(0).transpose();
    
    double projectionValue2 = dataAdjust.row(1) * EigenVectors.col(1);
    MatrixXd eigenTexture2(1, dim);
    eigenTexture2 = projectionValue2 * EigenVectors.col(1).transpose();
    MatrixXd OriginalData(1, dim);
    OriginalData = eigenTexture1 + eigenTexture2 + mean;
    
    
    //    pca_dim = 2;
    //
    //    // Extrace pca_dim eigenvectors which represent the original eigenspace adequately
    //    MatrixXd extract_eigenVectors(dim, pca_dim);
    //    for (int i = 0; i<pca_dim; i++){
    //        extract_eigenVectors.col(i) = EigenVectors.col(i);
    //    }
    //
    //    MatrixXd projection(num, pca_dim);
    //    projection = dataAdjust * extract_eigenVectors;
    //
    //    double projectionTest = dataAdjust.row(0) * extract_eigenVectors.col(0);
    //
    //    MatrixXd OriginalData(num, dim);
    //
    //    OriginalData = (extract_eigenVectors * projection.transpose()).transpose();
    //    for (int i = 0; i<num; i++) {
    //        OriginalData.row(i) += mean;
    //    }
    
    cout << endl;
    //    cout << projectionTest/400.0f << endl;
    //    cout <<projection/400.0f<< endl;
    //    cout<<projectionValue/400.0f<<endl;
    cout << "Original Data: " << endl;
    cout << projectionValue1 / 400.0f << endl;
    cout << projectionValue2 / 400.0f << endl;
    cout << OriginalData << endl;
    return OriginalData;// EigenVectors.col(0).transpose() * 400;
}

MatrixXd PCA::Eigen(MatrixXd data, unsigned char* newImage, string name, string fileName) {
    //cout << endl << endl;
    cout << name << endl;
    cout<<data<<endl;
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
    
    
//    return sorted_eigenVectors;
    return mean;
    
}



MatrixXd PCA::BubbleSort(MatrixXd eigenValues, int num) {
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





















