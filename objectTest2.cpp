//Massey University
//ID:15145234 Name:Kyuwon Shim
//Mail:kyu.shim88@gmail.com

#include "opencv2/opencv.hpp"
#include <stdio.h>
#include <chrono>
#include <ctime>
#include <iostream>
#include <bits/stdc++.h>
#include <opencv2/core/core.hpp>
#include <vector>
#include <string>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;
using namespace std;
using namespace chrono;

//g++ -std=c++11 objectTest2.cpp -o objectTest2 `pkg-config --cflags --libs opencv`
/*--------------- SubtractionTest ---------------*/


#define Mpixel(image,x,y) ((uchar *)((image).data +(y)*((image).step)))[(x)*((image).channels())]//gray color space

#define pixelB(image,x,y) image.data[image.step[0]*y+image.step[1]*x]	//Blue color space
#define pixelG(image,x,y) image.data[image.step[0]*y+image.step[1]*x+1]	//Green color space
#define pixelR(image,x,y) image.data[image.step[0]*y+image.step[1]*x+2]	//Red color space


void Integral_test(Mat image1, double **integral_image){
//This is the function for printing out 1/20 size of original metrices values and summed-area table values.
//Checking whether properly summed or not. 	

	/*
		1/20 size of original metrices values
	*/
	cout<<"The First image"<<endl;
	for(int i=0;i<image1.cols/20;i++){
		cout<<"[ ";
		for(int j=0;j<image1.rows/20;j++){
			cout<<(int)Mpixel(image1,i,j)<<" ";
		}
		cout<<" ]"<<endl;
	}
	///////////////////////////////////////////

	/*
		summed-area table values
	*/
	cout<<"The Second image"<<endl;
	for(int i=0;i<image1.cols;i++){
		cout<<"[ ";
		for(int j=0;j<image1.rows;j++){
			cout<<integral_image[i][j]<<" ";
		}
		cout<<" ]"<<endl;
	}
	//////////////////////////////////////////
}

void Integral_Gray_Initialize(Mat image1, double **integral_image){
	//This is the function for creating summed-area table in gray image.
	for(int i=0;i<image1.cols;i++){
		for(int j=0;j<image1.rows;j++){
			if((i==0)||(j==0)){//The value of the first row and column is zero. 
				integral_image[i][j]=0;
			}else{
				integral_image[i][j]=integral_image[i-1][j]+integral_image[i][j-1]-integral_image[i-1][j-1]+(int)Mpixel(image1,i-1,j-1);
			}
			
		}
	}
}

void Integral_Color_Image(Mat3b image1, double **integral_image_B, double **integral_image_G, double **integral_image_R){
	//This is the function for creating summed-area table in color image.
	for(int i=0;i<image1.cols;i++){
		for(int j=0;j<image1.rows;j++){
			if((i==0)||(j==0)){//The value of the first row and column is zero.
				integral_image_B[i][j]=0;
				integral_image_G[i][j]=0;
				integral_image_R[i][j]=0;
			}else{
				integral_image_B[i][j]=integral_image_B[i-1][j]+integral_image_B[i][j-1]-integral_image_B[i-1][j-1]+(int)pixelB(image1,i-1,j-1);
				integral_image_G[i][j]=integral_image_G[i-1][j]+integral_image_G[i][j-1]-integral_image_G[i-1][j-1]+(int)pixelG(image1,i-1,j-1);
				integral_image_R[i][j]=integral_image_R[i-1][j]+integral_image_R[i][j-1]-integral_image_R[i-1][j-1]+(int)pixelR(image1,i-1,j-1);
			}
			
		}
	}
}



void Filter_Gray_Integral(Mat image1, Mat image2, double** integral_image, int window_size){
	//This is function for applying the Kuwahara filter to image2, it must require summed-area table.

	//Mat image1 is the source image with gray scale. 
	//Mat image2 is the final output which is applied to Kuwahara filter in gray scale image.
	//double** integral_image is 2-dimensional array which stores the summed-area table.
	//int window_size is the size of window it is big more blurly when the value is high.

	int picture_x_size=image1.cols;//the x-axis length of source image.
	int picture_y_size=image1.rows;//the y-axis length of source image.
	int window_x_size=window_size;//the x-axis length of window.
	int window_y_size=window_size;//the y-axis length of window.


	/* The function of Kuhawara filter*/
			
	// mk=1(n+1)×(n+1)×∑(x,y)∈θkφ(f(x,y)

	// k∈{0,1,2,3},fis the source image function,
	// f(x,y)is the value of the pixel at coordinates(x,y),
	// φis a function calculating the value of a particular pixel,
	// 1/(n+1)×(n+1)is the number of pixels in the current area,
	// n is the value obtained directly from the filter windowsize.

	// reference by 
	// https://link-springer-com.ezproxy.massey.ac.nz/content/pdf/10.1007%2Fs11760-015-0791-3.pdf
	// page 665
	/**********************************/

	for(int i=(window_x_size/2);i<picture_x_size-(window_x_size/2+1);i++){
		for(int j=(window_y_size/2);j<picture_y_size-(window_y_size/2+1);j++){
			
			int f[4];
			int small_window_size=(window_x_size/2)*(window_y_size/2);

			int i1=i+1;
			int j1=j+1;

			int i2=i+1+(window_x_size/2);
			int j2=j+1;

			int i3=i+1;
			int j3=j+1+(window_y_size/2);

			int i4=i+1+(window_x_size/2);
			int j4=j+1+(window_y_size/2);

			



			f[0]=integral_image[i1][j1]-integral_image[i1-(window_x_size/2)][j1]
			-integral_image[i1][j1-(window_y_size/2)]+integral_image[i1-(window_x_size/2)][j1-(window_y_size/2)];
			f[0]=f[0]/small_window_size;


			f[1]=integral_image[i2][j2]-integral_image[i2-(window_x_size/2)][j2]
			-integral_image[i2][j2-(window_y_size/2)]+integral_image[i2-(window_x_size/2)][j2-(window_y_size/2)];
			f[1]=f[1]/small_window_size;

			f[2]=integral_image[i3][j3]-integral_image[i3-(window_x_size/2)][j3]
			-integral_image[i3][j3-(window_y_size/2)]+integral_image[i3-(window_x_size/2)][j3-(window_y_size/2)];
			f[2]=f[2]/small_window_size;

			f[3]=integral_image[i4][j4]-integral_image[i4-(window_x_size/2)][j4]
			-integral_image[i4][j4-(window_y_size/2)]+integral_image[i4-(window_x_size/2)][j4-(window_y_size/2)];
			f[3]=f[3]/small_window_size;

			
			/*
				The process finding the smallest value between f[0] and f[3].
			*/
			int final=9999;
			for(int l=0;l<4;l++){
				if(final>f[l]){
					final=f[l];
				}
			}

			Mpixel(image2,i,j)=final;
			//////////////////////////// 

		}
	}

}


void Filter_Color_Integral(Mat3b image1, Mat3b image2, double** integral_image_B, double** integral_image_G, double** integral_image_R, int window_size){
	// This is function for applying the Kuwahara filter to image2, it must require summed-area table.

	// Mat image1 is the source image with RGB scale. 
	// Mat image2 is the final output which is applied to Kuwahara filter in RGB image.
	// double** integral_image_B is 2-dimensional array which stores the summed-area table of Blue color space.
	// double** integral_image_G is 2-dimensional array which stores the summed-area table of Green color space.
	// double** integral_image_R is 2-dimensional array which stores the summed-area table of Red color space.
	// int window_size is the size of window it is big more blurly when the value is high.

	// This is color version of kuwahara filter, the process is mostly similar with gray version.

	int picture_x_size=image1.cols;
	int picture_y_size=image1.rows;
	int window_x_size=window_size;
	int window_y_size=window_size;

	for(int i=(window_x_size/2);i<picture_x_size-(window_x_size/2+1);i++){
		for(int j=(window_y_size/2);j<picture_y_size-(window_y_size/2+1);j++){
			int f_B[4];
			int f_G[4];
			int f_R[4];

			int small_window_size=(window_x_size/2)*(window_y_size/2);


			int i1=i+1;
			int j1=j+1;

			int i2=i+1+(window_x_size/2);
			int j2=j+1;

			int i3=i+1;
			int j3=j+1+(window_y_size/2);

			int i4=i+1+(window_x_size/2);
			int j4=j+1+(window_y_size/2);


			f_B[0]=integral_image_B[i1][j1]-integral_image_B[i1-(window_x_size/2)][j1]
			-integral_image_B[i1][j1-(window_y_size/2)]+integral_image_B[i1-(window_x_size/2)][j1-(window_y_size/2)];
			f_B[0]=f_B[0]/small_window_size;

			f_B[1]=integral_image_B[i2][j2]-integral_image_B[i2-(window_x_size/2)][j2]
			-integral_image_B[i2][j2-(window_y_size/2)]+integral_image_B[i2-(window_x_size/2)][j2-(window_y_size/2)];
			f_B[1]=f_B[1]/small_window_size;

			f_B[2]=integral_image_B[i3][j3]-integral_image_B[i3-(window_x_size/2)][j3]
			-integral_image_B[i3][j3-(window_y_size/2)]+integral_image_B[i3-(window_x_size/2)][j3-(window_y_size/2)];
			f_B[2]=f_B[2]/small_window_size;

			f_B[3]=integral_image_B[i4][j4]-integral_image_B[i4-(window_x_size/2)][j4]
			-integral_image_B[i4][j4-(window_y_size/2)]+integral_image_B[i4-(window_x_size/2)][j4-(window_y_size/2)];
			f_B[3]=f_B[3]/small_window_size;


			f_G[0]=integral_image_G[i1][j1]-integral_image_G[i1-(window_x_size/2)][j1]
			-integral_image_G[i1][j1-(window_y_size/2)]+integral_image_G[i1-(window_x_size/2)][j1-(window_y_size/2)];
			f_G[0]=f_G[0]/small_window_size;

			f_G[1]=integral_image_G[i2][j2]-integral_image_G[i2-(window_x_size/2)][j2]
			-integral_image_G[i2][j2-(window_y_size/2)]+integral_image_G[i2-(window_x_size/2)][j2-(window_y_size/2)];
			f_G[1]=f_G[1]/small_window_size;

			f_G[2]=integral_image_G[i3][j3]-integral_image_G[i3-(window_x_size/2)][j3]
			-integral_image_G[i3][j3-(window_y_size/2)]+integral_image_G[i3-(window_x_size/2)][j3-(window_y_size/2)];
			f_G[2]=f_G[2]/small_window_size;

			f_G[3]=integral_image_G[i4][j4]-integral_image_G[i4-(window_x_size/2)][j4]
			-integral_image_G[i4][j4-(window_y_size/2)]+integral_image_G[i4-(window_x_size/2)][j4-(window_y_size/2)];
			f_G[3]=f_G[3]/small_window_size;


			f_R[0]=integral_image_R[i1][j1]-integral_image_R[i1-(window_x_size/2)][j1]
			-integral_image_R[i1][j1-(window_y_size/2)]+integral_image_R[i1-(window_x_size/2)][j1-(window_y_size/2)];
			f_R[0]=f_R[0]/small_window_size;

			f_R[1]=integral_image_R[i2][j2]-integral_image_R[i2-(window_x_size/2)][j2]
			-integral_image_R[i2][j2-(window_y_size/2)]+integral_image_R[i2-(window_x_size/2)][j2-(window_y_size/2)];
			f_R[1]=f_R[1]/small_window_size;

			f_R[2]=integral_image_R[i3][j3]-integral_image_R[i3-(window_x_size/2)][j3]
			-integral_image_R[i3][j3-(window_y_size/2)]+integral_image_R[i3-(window_x_size/2)][j3-(window_y_size/2)];
			f_R[2]=f_R[2]/small_window_size;

			f_R[3]=integral_image_R[i4][j4]-integral_image_R[i4-(window_x_size/2)][j4]
			-integral_image_R[i4][j4-(window_y_size/2)]+integral_image_R[i4-(window_x_size/2)][j4-(window_y_size/2)];
			f_R[3]=f_R[3]/small_window_size;


			
			/*
				The process finding the smallest value of Blue color space between f_B[0] and f_B[3].
			*/			
			int final=9999;
			for(int l=0;l<4;l++){
				if(final>f_B[l]){
					final=f_B[l];
				}
			}

			pixelB(image2,i,j)=final;
			/////////////////////////

			/*
				The process finding the smallest value of Green color space between f_G[0] and f_G[3].
			*/

			final=9999;
			for(int l=0;l<4;l++){
				if(final>f_G[l]){
					final=f_G[l];
				}
			}

			pixelG(image2,i,j)=final;
			/////////////////////////

			/*
				The process finding the smallest value of Red color space between f_R[0] and f_R[3].
			*/

			final=9999;
			for(int l=0;l<4;l++){
				if(final>f_R[l]){
					final=f_R[l];
				}
			}

			pixelR(image2,i,j)=final;
			/////////////////////////

		}
	}

}
void Filter_Gray(Mat image1, Mat image2, int window_size){
	// This is function for applying the Kuwahara filter to gray-image2 (without summed-table).

	// Mat image1 is the source image with Gray scale. 
	// Mat image2 is the final output which is applied to Kuwahara filter in Gray image.

	// int window_size is the size of window it is big more blurly when the value is high.


	int picture_x_size=image1.cols;//the x-axis length of source image.
	int picture_y_size=image1.rows;//the y-axis length of source image.
	int window_x_size=window_size;//the x-axis length of window. 
	int window_y_size=window_size;//the y-axis length of window.

	for(int i=(window_x_size/2);i<picture_x_size-(window_x_size/2);i++){
		for(int j=(window_y_size/2);j<picture_y_size-(window_y_size/2);j++){
			int value[4];
			int f[4];
			int small_window_size=((window_x_size/2)+1)*((window_y_size/2)+1);

			//	This is 3*3 window
			//  | | | |        
			//	| | | |
			//	| | | |



			//the section of f[0]
			//  |*|*| |        
			//	|*|*| |
			//	| | | |

			value[0]=0;
			for(int l=i-(window_x_size/2);l<=i;l++){
				for(int k=j-(window_y_size/2);k<=j;k++){
					value[0]=value[0]+(int)Mpixel(image1,l,k);
				}
			}

			f[0]=value[0]/small_window_size;


			//the section of f[1]
			//  | |*|*|        
			//	| |*|*|
			//	| | | |

			value[1]=0;
			for(int l=i;l<=i+(window_x_size/2);l++){
				for(int k=j-(window_y_size/2);k<=j;k++){
					value[1]=value[1]+(int)Mpixel(image1,l,k);
				}
			}

			f[1]=value[1]/small_window_size;


			//the section of f[2]
			//  | | | |        
			//	|*|*| |
			//	|*|*| |

			value[2]=0;
			for(int l=i-(window_x_size/2);l<=i;l++){
				for(int k=j;k<=j+(window_y_size/2);k++){
					value[2]=value[2]+(int)Mpixel(image1,l,k);
				}
			}

			f[2]=value[2]/small_window_size;


			//the section of f[3]
			//  | | | |        
			//	| |*|*|
			//	| |*|*|

			value[3]=0;
			for(int l=i;l<=i+(window_x_size/2);l++){
				for(int k=j;k<=j+(window_x_size/2);k++){
					value[3]=value[3]+(int)Mpixel(image1,l,k);
				}
			}

			f[3]=value[3]/small_window_size;

			
			/*
				The process finding the smallest value between f[0] and f[3].
			*/
			int final=9999;
			for(int l=0;l<4;l++){
				if(final>f[l]){
					final=f[l];
				}
			}

			Mpixel(image2,i,j)=final;
			/////////////////////////

		}
	}

}


void Filter(Mat3b image1, Mat3b image2, int window_size){
	// This is function for applying the Kuwahara filter to color-image2 (without summed-table).

	// Mat3b image1 is the source image with color. 
	// Mat3b image2 is the final output which is applied to Kuwahara filter in color image.

	// int window_size is the size of window it is big more blurly when the value is high.
	
	int picture_x_size=image1.cols;//the x-axis length of source image.
	int picture_y_size=image1.rows;//the y-axis length of source image.
	int window_x_size=window_size;//the x-axis length of window.
	int window_y_size=window_size;//the y-axis length of window.

	for(int i=(window_x_size/2);i<picture_x_size-(window_x_size/2);i++){
		for(int j=(window_y_size/2);j<picture_y_size-(window_y_size/2);j++){
			
			int valueB[4];
			int valueR[4];
			int valueG[4];
			int f_B[4];
			int f_R[4];
			int f_G[4];

			int small_window_size=((window_x_size/2)+1)*((window_y_size/2)+1);


			//	This is 3*3 window
			//  | | | |        
			//	| | | |
			//	| | | |



			//the section of f[0]
			//  |*|*| |        
			//	|*|*| |
			//	| | | |

			valueB[0]=0;
			valueR[0]=0;
			valueG[0]=0;
			for(int l=i-(window_x_size/2);l<=i;l++){
				for(int k=j-(window_y_size/2);k<=j;k++){
					valueB[0]=valueB[0]+(int)pixelB(image1,l,k);
					valueR[0]=valueR[0]+(int)pixelR(image1,l,k);
					valueG[0]=valueG[0]+(int)pixelG(image1,l,k);
				}
			}

			f_B[0]=valueB[0]/small_window_size;
			f_R[0]=valueR[0]/small_window_size;
			f_G[0]=valueG[0]/small_window_size;

			
			//the section of f[1]
			//  | |*|*|        
			//	| |*|*|
			//	| | | |

			valueB[1]=0;
			valueR[1]=0;
			valueG[1]=0;
			for(int l=i;l<=i+(window_x_size/2);l++){
				for(int k=j-(window_y_size/2);k<=j;k++){
					valueB[1]=valueB[1]+(int)pixelB(image1,l,k);
					valueR[1]=valueR[1]+(int)pixelR(image1,l,k);
					valueG[1]=valueG[1]+(int)pixelG(image1,l,k);
				}
			}

			f_B[1]=valueB[1]/small_window_size;
			f_R[1]=valueR[1]/small_window_size;
			f_G[1]=valueG[1]/small_window_size;

			//the section of f[2]
			//  | | | |        
			//	|*|*| |
			//	|*|*| |

			valueB[2]=0;
			valueR[2]=0;
			valueG[2]=0;
			for(int l=i-(window_x_size/2);l<=i;l++){
				for(int k=j;k<=j+(window_y_size/2);k++){
					valueB[2]=valueB[2]+(int)pixelB(image1,l,k);
					valueR[2]=valueR[2]+(int)pixelR(image1,l,k);
					valueG[2]=valueG[2]+(int)pixelG(image1,l,k);
				}
			}

			f_B[2]=valueB[2]/small_window_size;
			f_R[2]=valueR[2]/small_window_size;
			f_G[2]=valueG[2]/small_window_size;


			//the section of f[3]
			//  | | | |        
			//	| |*|*|
			//	| |*|*|

			valueB[3]=0;
			valueR[3]=0;
			valueG[3]=0;
			for(int l=i;l<=i+(window_x_size/2);l++){
				for(int k=j;k<=j+(window_x_size/2);k++){
					valueB[3]=valueB[3]+(int)pixelB(image1,l,k);
					valueR[3]=valueR[3]+(int)pixelR(image1,l,k);
					valueG[3]=valueG[3]+(int)pixelG(image1,l,k);
				}
			}

			f_B[3]=valueB[3]/small_window_size;
			f_R[3]=valueR[3]/small_window_size;
			f_G[3]=valueG[3]/small_window_size;


			/*
				The process finding the smallest value between f[0] and f[3].
			*/

			int final_r=9999;
			int final_g=9999;
			int final_b=9999;
			for(int l=0;l<4;l++){
				if(final_b>f_B[l]){
					final_b=f_B[l];
				}

				if(final_g>f_G[l]){
					final_g=f_G[l];
				}

				if(final_r>f_R[l]){
					final_r=f_R[l];
				}
			}

			pixelR(image2,i,j)=final_r;
			pixelG(image2,i,j)=final_g;
			pixelB(image2,i,j)=final_b;
			///////////////////////////
			

		}
	}

}


int main(int argc,char **argv){
	double fps=0.0;	
	if(argc==3){//static version require 2 images.
		/*The First image*/
        Mat3b image1;
        Mat gray_image1;
	   	
	   	image1=imread(argv[1],1);
	   	if(!image1.data){printf("Could not open the file\n"); exit(0);}
		cvtColor(image1,gray_image1, CV_BGR2GRAY);//color image1 to gray scale
	   	/*****************/

	   	/*The Second image process*/
   		Mat3b image2;
   		Mat gray_image2;
   		Mat output;
   		Mat3b final_output;

   		image2=imread(argv[2],1);
   		if(!image2.data){printf("Could not open the file\n"); exit(0);}
   		cvtColor(image2,gray_image2, CV_BGR2GRAY);//color image2 to gray scale
		
		
		output=Mat::zeros(gray_image2.size(),CV_LOAD_IMAGE_GRAYSCALE);//initialize the value of output metrices to zero
		final_output=Mat::zeros(gray_image2.size(),CV_8UC3);//initialize the value of final_output metrices to zero

		//////////////////////////////////////
		//Kuwahara filter using Summed-table//
		//////////////////////////////////////

		/*Memory Allocation*/
		double** integral_image=new double*[gray_image2.cols+1];
		for(int i = 0; i < gray_image2.cols+1; ++i)
			integral_image[i] = new double[gray_image2.rows+1];
		/*********************/

		Integral_Gray_Initialize(gray_image2,integral_image);//create summed-table to integral_image array.
		Filter_Gray_Integral(gray_image2,output,integral_image,35);//Applying kuwahara filter to output using integral_image.
		
		/*Memory deallocation*/
		for(int i = 0; i < gray_image1.cols+1; ++i) {
			delete [] integral_image[i];
		}
		delete [] integral_image;
		/***************/
		
		///////////////////////////////////////



		////////////////////////////////////////
		//Kuwahara filter without Summed-table//
		////////////////////////////////////////
		
		// Filter_Gray(gray_image2,output,15);
		
		///////////////////////////////////////
		

		/*subtraction process between The first image and the second image*/
		for(int x=0; x<gray_image2.cols;x++){
			for(int y=0; y<gray_image2.rows;y++){
				// Mpixel(output,x,y)=(float)Mpixel(gray_image1,x,y)-(float)Mpixel(gray_image2,x,y);
				Mpixel(output,x,y)=(float)Mpixel(gray_image2,x,y)-(float)Mpixel(gray_image1,x,y);
				// cout<<"x:"<<x<<" y: "<<y<<endl;

				if(Mpixel(output,x,y)<100){//thresholding
					Mpixel(output,x,y)=0;
				}
			}
		}

		/*****************************************************************/

		/*Load the subtracted area to final_output from source*/
		for(int x=0; x<gray_image2.cols;x++){
			for(int y=0; y<gray_image2.rows;y++){
				if((int)Mpixel(output,x,y)!=0){
					pixelB(final_output,x,y)=pixelB(image2,x,y);
					pixelG(final_output,x,y)=pixelG(image2,x,y);
					pixelR(final_output,x,y)=pixelR(image2,x,y);
				}

			}
		}
		/*******************************************************/

    	/*final_output show*/
    	imshow("final_output" ,final_output);
    	waitKey(0);
    	/*************/

	}else{//dynamic mode with camera
		cout<<"dynamic mode"<<endl;

		/*Camera setup*/
		VideoCapture cap;
	  	cap.open(0);
	 	if (!cap.isOpened()){
	        cout << "Failed to open camera" << endl;
	        return 0;
	    }
	    cout << "Opened camera" << endl;
	   	cap.set(CV_CAP_PROP_FRAME_WIDTH, 640);
        cap.set(CV_CAP_PROP_FRAME_HEIGHT, 480);
        /**************/

        /*The First image*/
        Mat3b image1;
        Mat gray_image1;
	   	int key=0;
	   	cap >> image1;
	   	if(!image1.data){printf("Could not open the file\n"); exit(0);}
	   	cvtColor(image1,gray_image1, CV_BGR2GRAY);//copy camera color image to gray scale
	   	/*****************/

	   	/*The Second image process*/
	   	while (1){
	   		system_clock::time_point start = system_clock::now();//start clock
	   		
	   		
	   		Mat3b image2;
	   		Mat gray_image2;
	   		Mat output;
	   		Mat3b final_output;
	   		cap >> image2;
			
			cvtColor(image2,gray_image2, CV_BGR2GRAY);//copy camera color image to gray scale
			
			if(!image2.data){printf("Could not open the file\n"); exit(0);}
			output=Mat::zeros(gray_image2.size(),CV_LOAD_IMAGE_GRAYSCALE);//initialize the value of output metrices to zero
    		final_output=Mat::zeros(gray_image2.size(),CV_8UC3);//initialize the value of final_output metrices to zero

    		//////////////////////////////////////
    		//Kuwahara filter using Summed-table//
    		//////////////////////////////////////

			/*Memory Allocation*/
			double** integral_image=new double*[gray_image2.cols+1];
			for(int i = 0; i < gray_image2.cols+1; ++i)
    			integral_image[i] = new double[gray_image2.rows+1];
    		/*********************/

    		Integral_Gray_Initialize(gray_image2,integral_image);//create summed-table to integral_image array.
    		Filter_Gray_Integral(gray_image2,output,integral_image,35);//Applying kuwahara filter to output using integral_image.
    		
    		/*Memory deallocation*/
			for(int i = 0; i < gray_image1.cols+1; ++i) {
    			delete [] integral_image[i];
			}
			delete [] integral_image;
			/***************/
			
			///////////////////////////////////////



			////////////////////////////////////////
    		//Kuwahara filter without Summed-table//
    		////////////////////////////////////////
    		
    		// Filter_Gray(gray_image2,output,15);
    		
    		///////////////////////////////////////
    		

			/*subtraction process between The first image and the second image*/
			for(int x=0; x<gray_image2.cols;x++){
				for(int y=0; y<gray_image2.rows;y++){
					Mpixel(output,x,y)=(float)Mpixel(gray_image1,x,y)-(float)Mpixel(gray_image2,x,y);
					// Mpixel(output,x,y)=(float)Mpixel(gray_image2,x,y)-(float)Mpixel(gray_image1,x,y);
					// cout<<"x:"<<x<<" y: "<<y<<endl;

					if(Mpixel(output,x,y)<100){//thresholding
						Mpixel(output,x,y)=0;
					}
				}
			}

			/*****************************************************************/

			/*Load the subtracted area to final_output from source*/
			for(int x=0; x<gray_image2.cols;x++){
				for(int y=0; y<gray_image2.rows;y++){
					if((int)Mpixel(output,x,y)!=0){
						pixelB(final_output,x,y)=pixelB(image2,x,y);
						pixelG(final_output,x,y)=pixelG(image2,x,y);
						pixelR(final_output,x,y)=pixelR(image2,x,y);
					}

				}
			}
			/*******************************************************/

	    	/*final_output show*/
	    	imshow("final_output" ,final_output);
	    	/*************/
	    	
	    	/*program termination*/
	    	key=waitKey(1);
	       	if(key==113 || key==27) return 0;//either esc or 'q'
	       	/****************/

	       	/*Caculate performance of program*/
	       	system_clock::time_point end = system_clock::now();
	       	double seconds = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	     	fps = 1000000/seconds;
	     	cout << "frames " << fps << " seconds " << seconds << endl;
	     	/*********************************/
	    }
	    /************/
	}
}