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
    cout<<"Finish loading PCA results and scores. (Duration time: "<<elapsed_load.count()<< " s)" <<endl;
    
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
    
    // Save PCA results and scores for real time
    this->pcas_b = pcas_b;
    this->pcas_g = pcas_g;
    this->pcas_r = pcas_r;
    this->scores_b = scores_b;
    this->scores_g = scores_g;
    this->scores_r = scores_r;
    
//    return path;
}

void PCA_::windowLoopLoad(int loop_num_components){
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
    cout<<"Finish loading PCA results and scores. (Duration time: "<<elapsed_load.count()<< " s)" <<endl;
    
    int x = image_height/cell_dimension;
    int y = image_width/cell_dimension;
    vector<PCA> pcas_b, pcas_g, pcas_r;
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            // Get cell image pca and score
            PCA pcaB, pcaG, pcaR;
            pcaB.mean = pca_b.mean(Rect(0,j+i*y,cell_dimension*cell_dimension,1));
            pcaB.eigenvalues = pca_b.eigenvalues(Rect(j+i*y,0,1,loop_num_components));
            pcaB.eigenvectors = pca_b.eigenvectors(Rect(0,(j+i*y)*smallerNum,cell_dimension*cell_dimension,loop_num_components));
            pcaG.mean = pca_g.mean(Rect(0,j+i*y,cell_dimension*cell_dimension,1));
            pcaG.eigenvalues = pca_g.eigenvalues(Rect(j+i*y,0,1,loop_num_components));
            pcaG.eigenvectors = pca_g.eigenvectors(Rect(0,(j+i*y)*smallerNum,cell_dimension*cell_dimension,loop_num_components));
            pcaR.mean = pca_r.mean(Rect(0,j+i*y,cell_dimension*cell_dimension,1));
            pcaR.eigenvalues = pca_r.eigenvalues(Rect(j+i*y,0,1,loop_num_components));
            pcaR.eigenvectors = pca_r.eigenvectors(Rect(0,(j+i*y)*smallerNum,cell_dimension*cell_dimension,loop_num_components));
            Mat scoreB = scores_b(Rect(0,imageIndex*cell_num+j+i*y,loop_num_components,1));
            Mat scoreG = scores_g(Rect(0,imageIndex*cell_num+j+i*y,loop_num_components,1));
            Mat scoreR = scores_r(Rect(0,imageIndex*cell_num+j+i*y,loop_num_components,1));
//            Mat b = pcaB.backProject(scoreB).reshape(1,cell_dimension);
//            Mat g = pcaG.backProject(scoreG).reshape(1,cell_dimension);
//            Mat r = pcaR.backProject(scoreR).reshape(1,cell_dimension);
            
            pcas_b.push_back(pcaB);
            pcas_g.push_back(pcaG);
            pcas_r.push_back(pcaR);
        }
    }
    
    // Release the original pcas data
//    for(int i=0; i<x*y; i++){
//        this->pcas_b[i].mean.deallocate();
//        this->pcas_b[i].eigenvalues.release();
//        this->pcas_b[i].eigenvectors.release();
//    }
//    this->pcas_b.clear();
//    this->pcas_g.clear();
//    this->pcas_r.clear();
//    this->scores_b.release();
//    this->scores_g.release();
//    this->scores_r.release();
    
    // Save PCA results and scores for real time
    this->pcas_b = pcas_b;
    this->pcas_g = pcas_g;
    this->pcas_r = pcas_r;
    this->scores_b = scores_b;
    this->scores_g = scores_g;
    this->scores_r = scores_r;
    
    //    return path;
}

vector<Mat> bs, gs,rs;


unsigned char* PCA_::reconstruct(int imageIndex, unsigned char* newImage, Mat bgr[3], int cell_num, int x, int y, bool isBlur){
    
    if(CPU_only){
        bs.clear();
        gs.clear();
        rs.clear();
        for(int i=0; i<x*y; i++){
            Mat b(1,cell_dimension*cell_dimension,CV_32F);
            Mat g(1,cell_dimension*cell_dimension,CV_32F);
            Mat r(1,cell_dimension*cell_dimension,CV_32F);
            for(int k=0; k<cell_dimension*cell_dimension; k++){
                b.at<float>(0,k) = 0;
                g.at<float>(0,k) = 0;
                r.at<float>(0,k) = 0;
            }
            bs.push_back(b);
            gs.push_back(g);
            rs.push_back(r);
        }
    }
    
    
    
//    cout<<"===Start reconstruction==="<<endl;
//    auto re1 = chrono::high_resolution_clock::now();
    
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            // Reconstruction image
            int index = j+i*y;
            Mat scoreB = scores_b(Rect(0,imageIndex*cell_num+j+i*y,num_components,1));
            Mat scoreG = scores_g(Rect(0,imageIndex*cell_num+j+i*y,num_components,1));
            Mat scoreR = scores_r(Rect(0,imageIndex*cell_num+j+i*y,num_components,1));
            
            
            if(CPU_only){
                // Reconstruction without using opencv pca function
                Mat b = bs[index];
                Mat g = gs[index];
                Mat r = rs[index];
//                auto re1 = chrono::high_resolution_clock::now();
                for(int k=0; k<num_components; k++){
                    bs[index] = bs[index] + scoreB.at<float>(0,k) * pcas_b[index].eigenvectors.row(k);
                    gs[index] = gs[index] + scoreG.at<float>(0,k) * pcas_g[index].eigenvectors.row(k);
                    rs[index] = rs[index] + scoreR.at<float>(0,k) * pcas_r[index].eigenvectors.row(k);
                }
                bs[index] = bs[index] + pcas_b[index].mean;
                gs[index] = gs[index] + pcas_g[index].mean;
                rs[index] = rs[index] + pcas_r[index].mean;
//                auto re2 = chrono::high_resolution_clock::now();
//                chrono::duration<double> elapsed_re = re2 - re1;
//                this->cccTime += elapsed_re.count();
//                this->ccc++;
                bs[index] = bs[index].reshape(1,cell_dimension);
                gs[index] = gs[index].reshape(1,cell_dimension);
                rs[index] = rs[index].reshape(1,cell_dimension);
                bs[index].copyTo(bgr[0](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
                gs[index].copyTo(bgr[1](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
                rs[index].copyTo(bgr[2](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            }else{
//                auto re1 = chrono::high_resolution_clock::now();
                Mat b = pcas_b[index].backProject(scoreB).reshape(1,cell_dimension);
                Mat g = pcas_g[index].backProject(scoreG).reshape(1,cell_dimension);
                Mat r = pcas_r[index].backProject(scoreR).reshape(1,cell_dimension);
//                auto re2 = chrono::high_resolution_clock::now();
//                chrono::duration<double> elapsed_re = re2 - re1;
//                this->cccTime += elapsed_re.count();
//                this->ccc++;
                b.copyTo(bgr[0](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
                g.copyTo(bgr[1](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
                r.copyTo(bgr[2](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            }
            
            
//            Mat b = pcas_b[index].backProject(scoreB).reshape(1,cell_dimension);
//            Mat g = pcas_g[index].backProject(scoreG).reshape(1,cell_dimension);
//            Mat r = pcas_r[index].backProject(scoreR).reshape(1,cell_dimension);
            
            
            
        }
    }
//    Mat reconstruction_b = norm_0_255(bgr[0]);
//    Mat reconstruction_g = norm_0_255(bgr[1]);
//    Mat reconstruction_r = norm_0_255(bgr[2]);
    // Put three channel into one bgr image
    Mat resultImage = mergeChannels(norm_0_255(bgr[0]), norm_0_255(bgr[1]), norm_0_255(bgr[2]), oneImagePath);
    
    if(isBlur){
        Mat smoothedImage = resultImage.clone();
        int kernelSize = 3;
        medianBlur (smoothedImage, smoothedImage, kernelSize);
        for(int i=1; i<resultImage.rows-1; i++){
            for(int j=1; j<resultImage.cols-1; j++){
                if(j%cell_dimension==0){
                    resultImage.at<Vec3b>(i, j-1) = smoothedImage.at<Vec3b>(i, j-1);
                    resultImage.at<Vec3b>(i, j) = smoothedImage.at<Vec3b>(i, j);
                    resultImage.at<Vec3b>(i, j+1) = smoothedImage.at<Vec3b>(i, j+1);
                    resultImage.at<Vec3b>(i, j+2) = smoothedImage.at<Vec3b>(i, j+2);
                }
                if(i%cell_dimension==0){
                    resultImage.at<Vec3b>(i-1, j) = smoothedImage.at<Vec3b>(i-1, j);
                    resultImage.at<Vec3b>(i, j) = smoothedImage.at<Vec3b>(i, j);
                    resultImage.at<Vec3b>(i+1, j) = smoothedImage.at<Vec3b>(i+1, j);
                    resultImage.at<Vec3b>(i+2, j) = smoothedImage.at<Vec3b>(i+2, j);
                }
            }
        }
    }
    
//    auto re2 = chrono::high_resolution_clock::now();
//    chrono::duration<double> elapsed_re = re2 - re1;
//    cout<<"Reconstruction time: "<<elapsed_re.count()<<" s"<<endl;
    
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

unsigned char* PCA_::windowLoopReconstruct(int imageIndex, unsigned char* newImage, Mat bgr[3], int cell_num, int x, int y, bool isBlur, int loop_num_components){
    
    if(CPU_only){
        bs.clear();
        gs.clear();
        rs.clear();
        for(int i=0; i<x*y; i++){
            Mat b(1,cell_dimension*cell_dimension,CV_32F);
            Mat g(1,cell_dimension*cell_dimension,CV_32F);
            Mat r(1,cell_dimension*cell_dimension,CV_32F);
            for(int k=0; k<cell_dimension*cell_dimension; k++){
                b.at<float>(0,k) = 0;
                g.at<float>(0,k) = 0;
                r.at<float>(0,k) = 0;
            }
            bs.push_back(b);
            gs.push_back(g);
            rs.push_back(r);
        }
    }
    
    
    
    //    cout<<"===Start reconstruction==="<<endl;
    //    auto re1 = chrono::high_resolution_clock::now();
    
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            // Reconstruction image
            int index = j+i*y;
            Mat scoreB = scores_b(Rect(0,imageIndex*cell_num+j+i*y,loop_num_components,1));
            Mat scoreG = scores_g(Rect(0,imageIndex*cell_num+j+i*y,loop_num_components,1));
            Mat scoreR = scores_r(Rect(0,imageIndex*cell_num+j+i*y,loop_num_components,1));
            
            
            if(CPU_only){
                // Reconstruction without using opencv pca function
                Mat b = bs[index];
                Mat g = gs[index];
                Mat r = rs[index];
                
                for(int k=0; k<loop_num_components; k++){
                    bs[index] = bs[index] + scoreB.at<float>(0,k) * pcas_b[index].eigenvectors.row(k);
                    gs[index] = gs[index] + scoreG.at<float>(0,k) * pcas_g[index].eigenvectors.row(k);
                    rs[index] = rs[index] + scoreR.at<float>(0,k) * pcas_r[index].eigenvectors.row(k);
                }
                bs[index] = bs[index] + pcas_b[index].mean;
                gs[index] = gs[index] + pcas_g[index].mean;
                rs[index] = rs[index] + pcas_r[index].mean;
                bs[index] = bs[index].reshape(1,cell_dimension);
                gs[index] = gs[index].reshape(1,cell_dimension);
                rs[index] = rs[index].reshape(1,cell_dimension);
                bs[index].copyTo(bgr[0](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
                gs[index].copyTo(bgr[1](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
                rs[index].copyTo(bgr[2](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            }else{
                Mat b = pcas_b[index].backProject(scoreB).reshape(1,cell_dimension);
                Mat g = pcas_g[index].backProject(scoreG).reshape(1,cell_dimension);
                Mat r = pcas_r[index].backProject(scoreR).reshape(1,cell_dimension);
                
                b.copyTo(bgr[0](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
                g.copyTo(bgr[1](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
                r.copyTo(bgr[2](Rect(i*cell_dimension,j*cell_dimension,cell_dimension,cell_dimension)));
            }
        }
    }
    // Put three channel into one bgr image
    Mat resultImage = mergeChannels(norm_0_255(bgr[0]), norm_0_255(bgr[1]), norm_0_255(bgr[2]), oneImagePath);
    
    if(isBlur){
        Mat smoothedImage = resultImage.clone();
        int kernelSize = 3;
        medianBlur (smoothedImage, smoothedImage, kernelSize);
        for(int i=1; i<resultImage.rows-1; i++){
            for(int j=1; j<resultImage.cols-1; j++){
                if(j%cell_dimension==0){
                    resultImage.at<Vec3b>(i, j-1) = smoothedImage.at<Vec3b>(i, j-1);
                    resultImage.at<Vec3b>(i, j) = smoothedImage.at<Vec3b>(i, j);
                    resultImage.at<Vec3b>(i, j+1) = smoothedImage.at<Vec3b>(i, j+1);
                    resultImage.at<Vec3b>(i, j+2) = smoothedImage.at<Vec3b>(i, j+2);
                }
                if(i%cell_dimension==0){
                    resultImage.at<Vec3b>(i-1, j) = smoothedImage.at<Vec3b>(i-1, j);
                    resultImage.at<Vec3b>(i, j) = smoothedImage.at<Vec3b>(i, j);
                    resultImage.at<Vec3b>(i+1, j) = smoothedImage.at<Vec3b>(i+1, j);
                    resultImage.at<Vec3b>(i+2, j) = smoothedImage.at<Vec3b>(i+2, j);
                }
            }
        }
    }
    
    //    auto re2 = chrono::high_resolution_clock::now();
    //    chrono::duration<double> elapsed_re = re2 - re1;
    //    cout<<"Reconstruction time: "<<elapsed_re.count()<<" s"<<endl;
    
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





















