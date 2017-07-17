//
//  PCA.cpp
//  PCA_opengl_0522
//
//  Created by Yinghan Xu on 22/05/2017.
//  Copyright © 2017 Yinghan Xu. All rights reserved.
//

#include "PCA.hpp"
#define dev_cellImage_num 225


unsigned char* PCA::Reconstruction_RealTime(unsigned char* image, int imageIndex, int width, int height) {
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
	for (int i = 0; i < height / 20; i++) {
		for (int j = 0; j < width / 20; j++) {
			MatrixXd data_r_cell = ReconstructedData_Red[j + (height / 20)*i];
			MatrixXd data_g_cell = ReconstructedData_Green[j + (height / 20)*i];
			MatrixXd data_b_cell = ReconstructedData_Blue[j + (height / 20)*i];

			for (int x = 0; x < 20; x++) {
				for (int y = 0; y < 20; y++) {
					int piIndex = j * 20 + i * 300 * 20 + y + x * 300;
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

unsigned char* PCA::CalculateEigen(vector<unsigned char*> images, int width, int height) {

	this->num = images.size();
	this->dim = 20 * 20;// width*height;
	this->cellImage_num = (width*height) / (20 * 20);

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
	for (int i = 0; i < height / 20; i++) {
		for (int j = 0; j < width / 20; j++) {
			MatrixXd data_r_cell = ReconstructedData_Red[j + (height / 20)*i];
			MatrixXd data_g_cell = ReconstructedData_Green[j + (height / 20)*i];
			MatrixXd data_b_cell = ReconstructedData_Blue[j + (height / 20)*i];

			for (int x = 0; x < 20; x++) {
				for (int y = 0; y < 20; y++) {
					int piIndex = j * 20 + i * 300 * 20 + y + x * 300;
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

void PCA::SaveHighEigen(string newFileName, string oldFileName, int highDimension) {

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

void PCA::Mean(vector<MatrixXd> data, string fileName) {
	vector<MatrixXd> means;

	for (int i = 0; i<cellImage_num; i++) {
		MatrixXd mean = data[i].colwise().mean();
		means.push_back(mean);
	}

	ofstream newFile(fileName);
	for (int i = 0; i<cellImage_num; i++) {
		for (int j = 0; j<dim; j++) {
			newFile << means[i](0, j) << " ";
		}
		newFile << endl;
	}
	newFile.close();
}

void PCA::CalScores(vector<MatrixXd> whole_data, string fileName, string meanFileName, string newFileName) {

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

	for (int cellIndex = 0; cellIndex<cellImage_num; cellIndex++) {
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
		for (int k = 0; k<num; k++) {
			for (int i = 0; i<40; i++) {
				double projectionValue = dataAdjust.row(k) * oneVec.col(i);
				scoresFile << projectionValue << " ";
			}
			scoresFile << endl;
		}
	}

	scoresFile.close();
}


vector<MatrixXd> PCA::WholeReconstruction(string fileName, string meanFileName, string scoresFileName, int imageIndex, string channel) {
	cout << "Reconstruction data: " << endl << endl;

	cout << "Get Eigen Textures:" << endl << endl;
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
	if (channel == "Red") {
		this->EigenVectors_Red = EigenVectors;
	}
	else if (channel == "Green") {
		this->EigenVectors_Green = EigenVectors;
	}
	else if (channel == "Blue") {
		this->EigenVectors_Blue = EigenVectors;
	}


	vector<MatrixXd> whole_OriginalData;

	cout << "Get mean image:" << endl << endl;
	// Get mean image
	MatrixXd means(cellImage_num, dim);
	//Vmatrix cuda_means(dev_cellImage_num, dev_dim);
	vector<Arr_dim> cuda_means;
	ifstream meanfile(meanFileName);
	for (int i = 0; i < cellImage_num; i++) {
		Arr_dim tmp;
		for (int j = 0; j < dim; j++) {
			meanfile >> means(i, j);
			tmp.data[j] = means(i, j);
			//cuda_means.data[i][j] = means(i, j);
		}
		cuda_means.push_back(tmp);
	}
	meanfile.close();


	// Save mean image for real-time rendering
	if (channel == "Red") {
		this->means_Red = means;
		this -> cuda_means_Red = cuda_means;
	}
	else if (channel == "Green") {
		this->means_Green = means;
		this->cuda_means_Green = cuda_means;
	}
	else if (channel == "Blue") {
		this->means_Blue = means;
		this->cuda_means_Blue = cuda_means;
	}

	cout << "Get scores:" << endl << endl;
	// Get scores
	MatrixXd whole_scores(cellImage_num*num, 40);
	ifstream scoresFile(scoresFileName);
	for (int i = 0; i < cellImage_num*num; i++) {
		for (int j = 0; j < 40; j++) {
			scoresFile >> whole_scores(i, j);
		}
	}
	scoresFile.close();

	cout << "Get one scores:" << endl << endl;
	MatrixXd oneImage_scores(cellImage_num, 40);
	vector<MatrixXd> wholeImage_scores;
	Vmatrix cuda_oneImage_scores(cellImage_num, 40);
	vector<Vmatrix> cuda_wholeImage_scores;
	/*for (int i = 0; i<cellImage_num; i++) {
		for (int j = 0; j<40; j++) {
			oneImage_scores(i, j) = whole_scores(imageIndex + i*num, j);
		}
	}*/

	// Save all score for real-time rendering
	for (int k = 0; k<36; k++) {
		for (int i = 0; i<cellImage_num; i++) {
			for (int j = 0; j<40; j++) {
				oneImage_scores(i, j) = whole_scores(k + i*num, j);
				cuda_oneImage_scores.data[i][j] = oneImage_scores(i, j);
			}
		}
		wholeImage_scores.push_back(oneImage_scores);
		cuda_wholeImage_scores.push_back(cuda_oneImage_scores);
	}


	if (channel == "Red") {
		this->wholeImage_scores_Red = wholeImage_scores;
		this->oneImage_scores_Red = oneImage_scores;
		this->cuda_wholeImage_scores_Red = cuda_wholeImage_scores;
		this->cuda_oneImage_scores_Red = cuda_oneImage_scores;
	}
	else if (channel == "Green") {
		this->wholeImage_scores_Green = wholeImage_scores;
		this->oneImage_scores_Green = oneImage_scores;
		this->cuda_wholeImage_scores_Green = cuda_wholeImage_scores;
		this->cuda_oneImage_scores_Green = cuda_oneImage_scores;
	}
	else if (channel == "Blue") {
		this->wholeImage_scores_Blue = wholeImage_scores;
		this->oneImage_scores_Blue = oneImage_scores;
		this->cuda_wholeImage_scores_Blue = cuda_wholeImage_scores;
		this->cuda_oneImage_scores_Blue = cuda_oneImage_scores;
	}

	// Save oneVectors
	vector<MatrixXd> oneVecs;
	vector<vector<Arr_dim>> cuda_oneVecs;
	for (int cellIndex = 0; cellIndex<cellImage_num; cellIndex++) {
		//        MatrixXd data = whole_data[cellIndex];
		//        
		MatrixXd mean = means.row(cellIndex);   //data.colwise().mean();  // mean  mean =
												//        MatrixXd dataAdjust(num, dim);  // data-mean
												//        for (int i = 0; i<num; i++) {
												//            dataAdjust.row(i) = data.row(i) - mean;  // dataAdjust = data-mean
												//        }

		MatrixXd oneVec(pca_dim, dim);
		vector<Arr_dim> cuda_oneVec;
		/*typedef double arr[dev_dim];
		vector<arr> cuda_oneVec;*/

		for (int i = 0; i < pca_dim; i++) {
			Arr_dim tmp;
			for (int j = 0; j < dim; j++) {
				oneVec(i, j) = EigenVectors(cellIndex*dim + j, i);
				tmp.data[j] = oneVec(i, j);
			}
			cuda_oneVec.push_back(tmp);
		}

		oneVecs.push_back(oneVec);
		cuda_oneVecs.push_back(cuda_oneVec);

		MatrixXd OriginalData = MatrixXd::Constant(1, dim, 0.0);
		for (int i = 0; i<pca_dim; i++) {
			//            double projectionValue = dataAdjust.row(imageIndex) * oneVec.col(i);
			double projectionValue = oneImage_scores(cellIndex, i);
			MatrixXd eigenTexture(1, dim);
			eigenTexture = projectionValue * oneVec.row(i);// .transpose();
			OriginalData += eigenTexture;
		}

		OriginalData += mean;

		whole_OriginalData.push_back(OriginalData);
	}

	if (channel == "Red") {
		this->oneVecs_Red = oneVecs;
		this->cuda_oneVecs_Red = cuda_oneVecs;
	}
	else if (channel == "Green") {
		this->oneVecs_Green = oneVecs;
		this->cuda_oneVecs_Green = cuda_oneVecs;
	}
	else if (channel == "Blue") {
		this->oneVecs_Blue = oneVecs;
		this->cuda_oneVecs_Blue = cuda_oneVecs;
	}

	return whole_OriginalData;
}

vector<MatrixXd> PCA::WholeReconstruction_RealTime(string meanFileName, string scoresFileName, int imageIndex, string channel) {
	cout << "Reconstruction data: " << endl;


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

	cout << "Get eigenvectors and scores" << endl;
	MatrixXd EigenVectors, oneImage_scores, means;
	vector<MatrixXd> oneVecs;
	Vmatrix cuda_oneImage_scores(cellImage_num, pca_dim);//cuda_means(cellImage_num, dim), 
	vector<Arr_dim> cuda_means;
	vector<vector<Arr_dim>> cuda_oneVecs;
	if (channel == "Red") {
		EigenVectors = EigenVectors_Red;
		oneImage_scores = wholeImage_scores_Red[imageIndex];
		oneVecs = oneVecs_Red;
		means = means_Red;
		cuda_means = cuda_means_Red;
		cuda_oneImage_scores = cuda_wholeImage_scores_Red[imageIndex];
		cuda_oneVecs = cuda_oneVecs_Red;
	}
	else if (channel == "Green") {
		EigenVectors = EigenVectors_Green;
		oneImage_scores = wholeImage_scores_Green[imageIndex];
		oneVecs = oneVecs_Green;
		means = means_Green;
		cuda_means = cuda_means_Green;
		cuda_oneImage_scores = cuda_wholeImage_scores_Green[imageIndex];
		cuda_oneVecs = cuda_oneVecs_Green;
	}
	else if (channel == "Blue") {
		EigenVectors = EigenVectors_Blue;
		oneImage_scores = wholeImage_scores_Blue[imageIndex];
		oneVecs = oneVecs_Blue;
		means = means_Blue;
		cuda_means = cuda_means_Blue;
		cuda_oneImage_scores = cuda_wholeImage_scores_Blue[imageIndex];
		cuda_oneVecs = cuda_oneVecs_Blue;
	}
	cout << "Finish get eigenvectors and scores" << endl;



	cout << "Start reconstruction" << endl;
	//// Use Cuda to speed up
	

	// Copy Eigen format data to cuda format data
	/*for (int i = 0; i < cellImage_num; i++) {
		for (int j = 0; j < dim; j++) {
			cuda_means.data[i][j] = means(i, j);
		}
	}
	for (int i = 0; i < cellImage_num; i++) {
		for (int j = 0; j < pca_dim; j++) {
			cuda_oneImage_scores.data[i][j] = oneImage_scores(i, j);
		}
	}
	for (int i = 0; i < cellImage_num; i++) {
		matrix cuda_oneVec(dim, pca_dim);
		for (int j = 0; j < dim; j++) {
			for (int k = 0; k < pca_dim; k++) {
				cuda_oneVec.data[j][k] = oneVecs[i](j, k);
			}
		}
		cuda_oneVecs.push_back(cuda_oneVec);
	}*/
	// Cuda acceleration
	
	cout << "Pass data to cuda" << endl;
	host_vector<Vmatrix> cuda_OriginalData = 
	cudaReconstruction(cuda_means, cuda_oneImage_scores, cuda_oneVecs,
		cellImage_num, dim, pca_dim);
	cout << "Get data back from cuda" << endl;
	// Copy cuda format data to Eigen format data
	for (int cellIndex = 0; cellIndex < cellImage_num; cellIndex++) {
		MatrixXd OriginalData = MatrixXd::Constant(1, dim, 0.0);
		for (int i = 0; i < dim; i++) {
			OriginalData(0, i) = cuda_OriginalData[cellIndex].data[0][i];
		}
		whole_OriginalData.push_back(OriginalData);
	}
	//// Free memory
	//cuda_means.free();
	//cuda_oneImage_scores.free();
	//for (int cellIndex = 0; cellIndex < cellImage_num; cellIndex++) {
	//	cuda_oneVecs[cellIndex].free();
	//	cuda_OriginalData[cellIndex].free();
	//}

	/*for (int cellIndex = 0; cellIndex<cellImage_num; cellIndex++) {
		MatrixXd mean = means.row(cellIndex); 

		MatrixXd oneVec = oneVecs[cellIndex];

		MatrixXd OriginalData = MatrixXd::Constant(1, dim, 0.0);
		for (int i = 0; i<pca_dim; i++) {
			double projectionValue = oneImage_scores(cellIndex, i);
			MatrixXd eigenTexture(1, dim);
			eigenTexture = projectionValue * oneVec.col(i).transpose();
			OriginalData += eigenTexture;
		}

		OriginalData += mean;

		whole_OriginalData.push_back(OriginalData);
	}*/
	cout << "Finish reconstruction" << endl;

	return whole_OriginalData;
}

// Reconstruct one cell
MatrixXd PCA::Reconstruction(MatrixXd data, string fileName) {
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
			oneVec(i, j) = EigenVectors(149 * dim + i, j);
		}
	}

	cout << oneVec << endl;

	pca_dim = 40;

	int imageIndex = 0;

	MatrixXd OriginalData = MatrixXd::Constant(1, dim, 0.0);
	for (int i = 0; i<pca_dim; i++) {
		double projectionValue = dataAdjust.row(imageIndex) * oneVec.col(i);
		MatrixXd eigenTexture(1, dim);
		eigenTexture = projectionValue * oneVec.col(i).transpose();
		OriginalData += eigenTexture;
	}

	OriginalData += mean;

	cout << "Reconstruction Data: " << endl;
	cout << OriginalData << endl;
	return OriginalData;// EigenVectors.col(0).transpose() * 400;
}

MatrixXd PCA::Eigen(MatrixXd data, unsigned char* newImage, string name, string fileName) {
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












