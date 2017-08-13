////
////  simpleOpenCL.cpp
////  OpenCL_Test
////
////  Created by Yinghan Xu on 10/08/2017.
////  Copyright © 2017 Yinghan Xu. All rights reserved.
////
//
//#include <iostream>
//#include <OpenCL/opencl.h>
//
//using namespace std;
//
//int main(){
//    //get all platforms (drivers)
//    std::vector<cl::Platform> all_platforms;
//    cl::Platform::get(&all_platforms);
//    if(all_platforms.size()==0){
//        std::cout<<" No platforms found. Check OpenCL installation!\n";
//        exit(1);
//    }
//    cl::Platform default_platform=all_platforms[0];
//    std::cout << "Using platform: "<<default_platform.getInfo<CL_PLATFORM_NAME>()<<"\n";
//}
//
