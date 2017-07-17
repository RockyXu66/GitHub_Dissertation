//
//  ImageVectors.h
//  PCA_opencv
//
//  Created by Yinghan Xu on 27/05/2017.
//  Copyright © 2017 Yinghan Xu. All rights reserved.
//

#ifndef ImageVectors_h
#define ImageVectors_h

#include <iostream>
using namespace std;

class ImageVector{
public:
//    int imageVector
    vector<int> value;
    ImageVector(int width, int height){
        for(int i=0; i<width*height; i++){
            value.push_back(0);
        }
    }
};


#endif /* ImageVectors_h */
