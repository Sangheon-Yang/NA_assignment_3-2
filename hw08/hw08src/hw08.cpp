//  hw08.cpp
//  NA_hw08_OPENCV
//
//  Created by 양상헌 on 11/11/2019.
//  Copyright © 2019 양상헌. All rights reserved.
//

//file I/O and resize with interpolation
#include <opencv2/opencv.hpp>
#include <iostream>
using namespace std;
using namespace cv;

//the place where edition is needed.
void BilinearInterpolation(cv::Mat &src, cv::Mat &dst){
    // bilinear interpolation
    // using built-in-function in OPENCV;
    cv::resize(src, dst, cv::Size( dst.cols, dst.rows ), 0, 0, cv::INTER_LINEAR);
}

// implementing Bilinear interpolation NOT using built-in OPENCV function
void BilinearInterpolationMyown(cv::Mat &src, cv::Mat &dst){
    // bilinear Interpolation
    // implementing MyOwn Resize function
    
    //-------------IndexError-------------------//
    //double rate_of_Height = dst.rows/src.rows;
    //double rate_of_Width = dst.cols/src.cols;
    
    //-----------Solved Index Error Problem---//
    double rate_of_Height = (double)dst.rows/(double)src.rows;
    double rate_of_Width = (double)dst.cols/(double)src.cols;
    
    for(int y = 0;y < dst.rows;y++){
        for(int x = 0; x< dst.cols;x++){
            
            //src image's index of Mat src[px][py] (in order to get pixel value)
            int px = (int)(x / rate_of_Width);
            int py = (int)(y / rate_of_Height);
            
            // calculating a, (1-a), b, (1-b)
            double fx0 = (double)x/(double)rate_of_Width - (double)px;
            double fx1 = 1 - fx0;
            double fy0 = (double)y/(double)rate_of_Height - (double)py;
            double fy1 = 1 - fy0;
            
            double area00 = fx1*fy1; //(1-a)*(1-b)
            double area01 = fx0*fy1; //a*(1-b)
            double area10 = fx1*fy0; //(1-a)*b
            double area11 = fx0*fy0; //a*b
            
            if(src.channels() == 1){//when src image is gray image
                unsigned char P00 = src.at<unsigned char>(py,px);
                unsigned char P01 = src.at<unsigned char>(py,px+1);
                unsigned char P10 = src.at<unsigned char>(py+1,px);
                unsigned char P11 = src.at<unsigned char>(py+1,px+1);
                dst.at<unsigned char>(y,x)= area00*P00+area01*P01+area10*P10+area11*P11;
            }
            
            else{ //case of RGB image
                Vec3b P00 = src.at<cv::Vec3b>(py,px);
                Vec3b P01 = src.at<Vec3b>(py,px+1);
                Vec3b P10 = src.at<Vec3b>(py+1,px);
                Vec3b P11 = src.at<Vec3b>(py+1,px+1);
                dst.at<Vec3b>(y,x)= area00*P00+area01*P01+area10*P10+area11*P11;
            }
        }
    }
}
cv::Mat bi_dst; // for output of built-in function
cv::Mat bi_dst2; // for output of Self-Implemented MyOwn function

int main()
{
    cv::Mat src;
    src = cv::imread("example.jpeg", 1);
    if (src.empty())
    {
        std::cout << "Cannot find an image" << std::endl;
        return -1;
    }
    cv::imshow("SourceImage", src);
    cv::waitKey(0);

    //¿ÃπÃ¡ˆ∏¶ 2πË ≈∞ø¸¿ª ∂ß¿« Bilinear interpolation
    int height = src.rows;
    int width = src.cols;
    
    //------------- edited --------------//
    printf("\nOriginal Image size : %d x %d \n", height , width );
    
    int newHeight, newWidth;
    std::cout << "\nDesired scale of height : "<< endl;
    std::cin >> newHeight;
    std::cout << "\nDesired scale of width : " << endl;
    std::cin >> newWidth;
    printf("\ndesired size of picture: %d x %d\n",newHeight, newWidth);
    
    bi_dst = cv::Mat(newHeight, newWidth, src.type(), cv::Scalar(0));
    bi_dst2 = cv::Mat(newHeight, newWidth, src.type(), cv::Scalar(0));
    
    //bilinearInterpolation
    BilinearInterpolation(src, bi_dst);
    BilinearInterpolationMyown(src, bi_dst2);
    cv::imshow("BIImage_of_built-in-function", bi_dst);
    cv::imshow("BIImage_of_Myown_function", bi_dst2);
    cv::waitKey(0);
    return 0;
}


