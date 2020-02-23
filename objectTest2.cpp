//g++ -std=c++11 objectTest2.cpp -o objectTest2 `pkg-config --cflags --libs opencv`
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


/*--------------- SubtractionTest ---------------*/


#define Mpixel(image,x,y) ((uchar *)((image).data +(y)*((image).step)))[(x)*((image).channels())]//gray color space

#define pixelB(image,x,y) image.data[image.step[0]*y+image.step[1]*x]	//Blue color space
#define pixelG(image,x,y) image.data[image.step[0]*y+image.step[1]*x+1]	//Green color space
#define pixelR(image,x,y) image.data[image.step[0]*y+image.step[1]*x+2]	//Red color space

void Grey_to_Color(Mat3b source_image, Mat filtered_image,Mat3b output_image){
	for(int x=0; x<source_image.cols;x++){
		for(int y=0; y<source_image.rows;y++){
			if(Mpixel(filtered_image,x,y)!=0){//when pixel is not zero
				pixelB(output_image,x,y)=pixelB(source_image,x,y);
				pixelG(output_image,x,y)=pixelG(source_image,x,y);
				pixelR(output_image,x,y)=pixelR(source_image,x,y);
			}
		}	
	}
}
void Determining_ROI_Size(Mat source_image, int *small_x, int *small_y,int *large_x,int *large_y){
	//This is function for finding coordinates (p1 and p2) of rectangle before cropping Region of Interest area

	// * p1(small_x,small_y)
	// @ p2(large_x,large_y)
	//  |*| | |        
	//	| | | |
	//	| | |@|


	*small_x=9999;
	*small_y=9999;
	*large_x=-9999;
	*large_y=-9999;
	for(int x=0; x<source_image.cols;x++){
		for(int y=0; y<source_image.rows;y++){
			if(Mpixel(source_image,x,y)!=0){//when pixel is not zero
				if(x>*large_x){
					*large_x=x;
				}
				if(y>*large_y){
					*large_y=y;
				}
				if(x<*small_x){
					*small_x=x;
				}
				if(y<*small_y){
					*small_y=y;
				}
			}
		}	
	}
}


void Image_stitching(Mat source_image,Mat filtered_image,Mat output_image){
	
	//vertical way
	int sm[source_image.cols];
	int lar[source_image.cols];
	for(int x=0; x<source_image.cols;x++){
		sm[x]=9999;
		lar[x]=-9999;
		for(int y=0; y<source_image.rows;y++){
			if(Mpixel(filtered_image,x,y)<30){//thresholding
											  //this is better to do seperate but I put on here for performancing
				Mpixel(filtered_image,x,y)=0;
			}

			if((int)Mpixel(filtered_image,x,y)!=0){//The process for detecting the beginning point of pixel and the last point of pixel of each column. 
				if(y<sm[x]){
					sm[x]=y;
				}
				if(y>lar[x]){
					lar[x]=y;
				}
			}
			if(x!=0){
				if( (y>=sm[x-1])&&(y<=lar[x-1])){//The vertical way copy all previous line pixels from the point of the first detected pixcel to the last.  
					Mpixel(output_image,x-1,y)=Mpixel(source_image,x-1,y);
				}
			}
			if( (x==source_image.cols-1)&&(y==source_image.rows-1) ){//The same process of upper process but the last one is exception.
				for(int a=sm[x]; a<lar[x];a++){
					Mpixel(output_image,x,a)=Mpixel(source_image,x,a);	
				}						
			}

		}
	}
	//horizontal way
	int sm2[source_image.rows];
	int lar2[source_image.rows];
	for(int x=0; x<source_image.rows;x++){
		sm2[x]=9999;
		lar2[x]=-9999;
		for(int y=0; y<source_image.cols;y++){
			if((int)Mpixel(output_image,y,x)!=0){//The process for detecting the beginning point of pixel and the last point of pixel of each row.
				if(y<sm2[x]){
					sm2[x]=y;
				}
				if(y>lar2[x]){
					lar2[x]=y;
				}
			}
			if(x!=0){
				if( (y>sm2[x-1])&&(y<lar2[x-1])  ){//The horizontal way copy all previous line pixels from the point of the first detected pixcel to the last.  
					Mpixel(output_image,y,x-1)=Mpixel(source_image,y,x-1);	
				}							
			}
			if( (x==source_image.rows-1)&&(y==source_image.cols-1) ){//the case of last
				for(int a=sm2[x]; a<lar2[x];a++){
					Mpixel(output_image,a,x)=Mpixel(source_image,a,x);	
				}						
			}
		}
	}
}
void median_filter(Mat image1,Mat image2,int window_size){
    //picture elements(pels)
    int function_size_input=window_size;//this is for window size

    int picture_x_size=image1.cols;
    int picture_y_size=image1.rows;
    int mdn;//median value in a window
    int ltmdn=0;// number of pels having gray levels less than mdn in a window
    int window_x_size=function_size_input;
    int window_y_size=function_size_input;
    //int hist[window_x_size*window_y_size];
    int index_for_array=0;
    int count_of_less_than_median=0;
    int median_value;
    int left_column[window_y_size];
    int right_column[window_y_size];
    // int left_column[window_y_size-1];
    // int right_column[window_y_size-1];
    int hist_index=0;
    int g1;

    int th=(window_x_size*window_y_size/2);
    
    // for(int i=0;i<picture_y_size;i++){
    //      Mpixel(image2,637,i)=255;
    //  }
    
    for(int i=(window_x_size/2);i<picture_x_size-(window_x_size/2);i++){
        
        int hist[256];
        for(int q=0;q<256;q++){
            hist[q]=0;
        }

        int index_for_hist=0;
        int y_sizeof=i+window_y_size;
        int x_sizeof=(window_x_size/2)+window_x_size;
        for(int a=i;a<y_sizeof;a++){
            for(int b=(window_x_size/2);b<x_sizeof;b++){
                index_for_hist=(int)Mpixel(image1,a,b);
                hist[index_for_hist]=hist[index_for_hist]+1;
            }
        }


        int counter_for_find_median=(window_x_size*window_y_size/2)+1;
        int counter_for_find_less_than_median=0;

        for(int z=0;z<256;z++){
            if(hist[z]!=0){ 
                counter_for_find_median=counter_for_find_median-hist[z];
                if(counter_for_find_median<=0){
                    median_value=z;
                    mdn=median_value;
                    break;
                }else{
                    counter_for_find_less_than_median
                    =counter_for_find_less_than_median+hist[z];
                }
            }       
        }

        ltmdn=counter_for_find_less_than_median;

        //Mpixel(image2,i,(window_y_size/2))=mdn;

        for(int j=(window_y_size/2)+1;j<picture_y_size-(window_y_size/2);j++){//j indicates picture column number

    
            int index_for_left_column=0;
            int index_for_right_column=0;

            for(int l=i;l<i+window_x_size;l++){
                left_column[index_for_left_column]=(int)Mpixel(image1,l,j);
                index_for_left_column++;

                right_column[index_for_right_column]=(int)Mpixel(image1,l,j+window_x_size-1);               
                index_for_right_column++;   
            }
            
                    
            

            for(int k=0;k<window_y_size;k++){
                g1=left_column[k];
                hist[g1]=hist[g1]-1;
                if(g1<mdn){
                    ltmdn=ltmdn-1;
                }
                g1=right_column[k];
                hist[g1]=hist[g1]+1;
                if(g1<mdn){
                    ltmdn=ltmdn+1;
                }

                if(ltmdn>th){
                    while(true){
                        mdn=mdn-1;
                        ltmdn=ltmdn-hist[mdn];
                        if(ltmdn<=th){
                            break;
                        }
                    }
                }else{
                    while(ltmdn+hist[mdn]<=th){
                        ltmdn=ltmdn+hist[mdn];
                        mdn=mdn+1;
                    }
                }
            }
            Mpixel(image2,i,j)=mdn; //determine pixel
        }
    }

    for(int j=0;j<window_size/2+1;j++){ 
        for(int i=0;i<picture_y_size;i++){
            Mpixel(image2,j,i)=0;
            Mpixel(image2,picture_x_size-1-j,i)=0;
        }
    }
    for(int j=0;j<window_size/2+1;j++){ 
        for(int i=window_size/2;i<picture_x_size-(window_size/2);i++){
            Mpixel(image2,i,j)=0;
            Mpixel(image2,i,picture_y_size-1-j)=0;
        }
    }

}

int FindTheLargestContour(std::vector<vector<Point>>contours){
    int largestcontour=0;
    long int largestsize=0;
    for(int i=0;i<contours.size();i++){
        if(largestsize<contours[i].size()){
            largestsize=contours[i].size();
            largestcontour=i;
        }
    }
    return largestcontour;
}


void Integral_test(Mat image1, double **integral_image, double **squared_integral_image){
//This is the function for printing out 1/20 size of original metrices values and summed-area table values.
//Checking whether properly summed or not. 	

	/*
		1/20 size of original metrices values
	*/
	cout<<"The First image"<<endl;
	for(int i=0;i<image1.cols;i++){
		cout<<"[ ";
		for(int j=0;j<image1.rows;j++){
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


	cout<<"The squared image"<<endl;
	for(int i=0;i<image1.cols/20;i++){
		cout<<"[ ";
		for(int j=0;j<image1.rows/20;j++){
			cout<<squared_integral_image[i][j]<<" ";
		}
		cout<<" ]"<<endl;
	}
	//////////////////////////////////////////
}

void Integral_Gray_Initialize(Mat image1, double **integral_image,double **squared_integral_image){
	//This is the function for creating summed-area table in gray image.
	for(int i=0;i<image1.cols+1;i++){
		for(int j=0;j<image1.rows+1;j++){
			if((i==0)||(j==0)){//The value of the first row and column is zero. 
				integral_image[i][j]=0;
				squared_integral_image[i][j]=0;
			}else{
				integral_image[i][j]=integral_image[i-1][j]+integral_image[i][j-1]-integral_image[i-1][j-1]+(int)Mpixel(image1,i-1,j-1);
				squared_integral_image[i][j]=squared_integral_image[i-1][j]+squared_integral_image[i][j-1]-squared_integral_image[i-1][j-1]
				+(int)Mpixel(image1,i-1,j-1)*(int)Mpixel(image1,i-1,j-1);
			}
			
		}
	}
}



void Integral_Color_Initialize(Mat3b image1, double **integral_image_B, double **integral_image_G, double **integral_image_R, double **squared_integral_image_B, double **squared_integral_image_G, double **squared_integral_image_R){
	//This is the function for creating summed-area table in color image.
	for(int i=0;i<image1.cols+1;i++){
		for(int j=0;j<image1.rows+1;j++){
			if((i==0)||(j==0)){//The value of the first row and column is zero.
				integral_image_B[i][j]=0;
				integral_image_G[i][j]=0;
				integral_image_R[i][j]=0;

				squared_integral_image_B[i][j]=0;
				squared_integral_image_G[i][j]=0;
				squared_integral_image_R[i][j]=0;
			}else{
				integral_image_B[i][j]=integral_image_B[i-1][j]+integral_image_B[i][j-1]-integral_image_B[i-1][j-1]+(int)pixelB(image1,i-1,j-1);
				integral_image_G[i][j]=integral_image_G[i-1][j]+integral_image_G[i][j-1]-integral_image_G[i-1][j-1]+(int)pixelG(image1,i-1,j-1);
				integral_image_R[i][j]=integral_image_R[i-1][j]+integral_image_R[i][j-1]-integral_image_R[i-1][j-1]+(int)pixelR(image1,i-1,j-1);

				squared_integral_image_B[i][j]=squared_integral_image_B[i-1][j]+squared_integral_image_B[i][j-1]-squared_integral_image_B[i-1][j-1]+(int)pixelB(image1,i-1,j-1)*(int)pixelB(image1,i-1,j-1);
				squared_integral_image_G[i][j]=squared_integral_image_G[i-1][j]+squared_integral_image_G[i][j-1]-squared_integral_image_G[i-1][j-1]+(int)pixelG(image1,i-1,j-1)*(int)pixelG(image1,i-1,j-1);
				squared_integral_image_R[i][j]=squared_integral_image_R[i-1][j]+squared_integral_image_R[i][j-1]-squared_integral_image_R[i-1][j-1]+(int)pixelR(image1,i-1,j-1)*(int)pixelR(image1,i-1,j-1);
			}
			
		}
	}
}


void Filter_Gray_Integral(Mat image1, Mat image2, double** integral_image,double** squared_integral_image, int window_size){
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
			int s_f[4];
			int result[4];
			int small_window_size=(window_x_size/2)*(window_y_size/2);

			int i1=i;
			int j1=j;

			int i2=i+(window_x_size/2);
			int j2=j;

			int i3=i;
			int j3=j+(window_y_size/2);

			int i4=i+(window_x_size/2);
			int j4=j+(window_y_size/2);

			int mean_f[4];

			



			// f[0]=integral_image[i1][j1]-integral_image[i1-(window_x_size/2)][j1]
			// -integral_image[i1][j1-(window_y_size/2)]+integral_image[i1-(window_x_size/2)][j1-(window_y_size/2)];


			f[0]=integral_image[i1-1][j1]-integral_image[i1][j1-1]
			-integral_image[i1-1][j1-1]+(int)Mpixel(image1,i1-1,j1-1);


			// s_f[0]=squared_integral_image[i1][j1]-squared_integral_image[i1-(window_x_size/2)][j1]
			// -squared_integral_image[i1][j1-(window_y_size/2)]+squared_integral_image[i1-(window_x_size/2)][j1-(window_y_size/2)];


			cout<<"s_f: "<<s_f[0]<<endl;
			cout<<"f[0]: "<<f[0]<<endl;
			// f[0]=f[0]/small_window_size;
			// cout<<"s_f: "<<s_f[0]<<endl;
			// cout<<"(f[0]*f[0])/small_window_size): "<<((f[0]*f[0])/small_window_size)<<endl;
			
			// result[0]=(f[0]*f[0]-((f[0]*f[0])/small_window_size))/small_window_size;

			result[0]=(s_f[0]-((f[0]*f[0])/small_window_size));
			// result[0]=sqrt(result[0]);


			f[1]=integral_image[i2][j2]-integral_image[i2-(window_x_size/2)][j2]
			-integral_image[i2][j2-(window_y_size/2)]+integral_image[i2-(window_x_size/2)][j2-(window_y_size/2)];

			s_f[1]=squared_integral_image[i2][j2]-squared_integral_image[i2-(window_x_size/2)][j2]
			-squared_integral_image[i2][j2-(window_y_size/2)]+squared_integral_image[i2-(window_x_size/2)][j2-(window_y_size/2)];
			// f[1]=f[1]/small_window_size;
			


			cout<<"s_f: "<<s_f[1]<<endl;
			cout<<"f[1]: "<<f[1]<<endl;
			// result[1]=(f[1]*f[1]-((f[1]*f[1])/small_window_size))/small_window_size;
			result[1]=(s_f[1]-((f[1]*f[1])/small_window_size));
			// result[1]=sqrt(result[1]);


			f[2]=integral_image[i3][j3]-integral_image[i3-(window_x_size/2)][j3]
			-integral_image[i3][j3-(window_y_size/2)]+integral_image[i3-(window_x_size/2)][j3-(window_y_size/2)];

			s_f[2]=squared_integral_image[i3][j3]-squared_integral_image[i3-(window_x_size/2)][j3]
			-squared_integral_image[i3][j3-(window_y_size/2)]+squared_integral_image[i3-(window_x_size/2)][j3-(window_y_size/2)];
			// f[2]=f[2]/small_window_size;


			cout<<"s_f: "<<s_f[2]<<endl;
			cout<<"f[2]: "<<f[2]<<endl;
			
			// result[2]=(f[2]*f[2]-((f[2]*f[2])/small_window_size))/small_window_size;
			result[2]=(s_f[2]-((f[2]*f[2])/small_window_size));
			// result[2]=sqrt(result[2]);


			f[3]=integral_image[i4][j4]-integral_image[i4-(window_x_size/2)][j4]
			-integral_image[i4][j4-(window_y_size/2)]+integral_image[i4-(window_x_size/2)][j4-(window_y_size/2)];

			s_f[3]=squared_integral_image[i4][j4]-squared_integral_image[i4-(window_x_size/2)][j4]
			-squared_integral_image[i4][j4-(window_y_size/2)]+squared_integral_image[i4-(window_x_size/2)][j4-(window_y_size/2)];
			// f[3]=f[3]/small_window_size;


			cout<<"s_f: "<<s_f[3]<<endl;
			cout<<"f[3]: "<<f[3]<<endl;

			
			// result[3]=(f[3]*f[3]-((f[3]*f[3])/small_window_size))/small_window_size;
			result[3]=(s_f[3]-((f[3]*f[3])/small_window_size));
			// result[3]=sqrt(result[3]);
			
			/*
				The process finding the smallest value between f[0] and f[3].
			*/
			int final=9999;
			for(int l=0;l<4;l++){
				cout<<"result: "<<result[l]<<endl;
				// if(final>f[l]){
				// 	final=f[l];
				// }
				if(final>result[l]){
					final=result[l];
				}
			}
			cout<<"final: "<<final<<endl;
			getchar();

			Mpixel(image2,i,j)=final;
			//////////////////////////// 

		}
	}

}

void Filter_Gray_Integral2(Mat image1, Mat image2, double** integral_image,double** squared_integral_image, int window_size){
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
			int s_f[4];
			int result[4];
			int small_window_size=(window_x_size/2+1)*(window_y_size/2+1);

			int i0=i;
			int j0=j;

			int i1=i+(window_x_size/2);
			int j1=j;

			int i2=i;
			int j2=j+(window_y_size/2);

			int i3=i+(window_x_size/2);
			int j3=j+(window_y_size/2);

			int mean_f[4];

			// cout<<"small_window_size: "<<small_window_size<<endl;




			// f[0]=integral_image[i1][j1]-integral_image[i1-(window_x_size/2)][j1]
			// -integral_image[i1][j1-(window_y_size/2)]+integral_image[i1-(window_x_size/2)][j1-(window_y_size/2)];


			f[0]=integral_image[i0+(window_x_size/2)][j0+(window_y_size/2)]-integral_image[i0-1][j0+(window_y_size/2)]
			-integral_image[i0+(window_y_size/2)][j0-1]+integral_image[i0-1][j0-1];

			s_f[0]=squared_integral_image[i0+(window_x_size/2)][j0+(window_y_size/2)]-squared_integral_image[i0-1][j0+(window_y_size/2)]
			-squared_integral_image[i0+(window_y_size/2)][j0-1]+squared_integral_image[i0-1][j0-1];

			// cout<<"integral_image[i1+(window_x_size/2)][j1+(window_y_size/2)]:"<<integral_image[i1+(window_x_size/2)][j1+(window_y_size/2)]<<endl;
			// cout<<"-integral_image[i1][j1+(window_y_size/2)]:"<<integral_image[i1][j1+(window_y_size/2)]<<endl;
			// cout<<"-integral_image[i1+(window_y_size/2)][j1]:"<<integral_image[i1+(window_y_size/2)][j1]<<endl;
			// cout<<"integral_image[i1][j1]:"<<integral_image[i1][j1]<<endl;

			result[0]=(s_f[0]-((f[0]*f[0])/small_window_size))/small_window_size;

			// cout<<"s_f[0]: "<<s_f[0]<<endl;
			// cout<<"f[0]: "<<	f[0]<<endl;


			f[1]=integral_image[i1+(window_x_size/2)][j1+(window_y_size/2)]-integral_image[i1-1][j1+(window_y_size/2)]
			-integral_image[i1+(window_y_size/2)][j1-1]+integral_image[i1-1][j1-1];

			s_f[1]=squared_integral_image[i1+(window_x_size/2)][j1+(window_y_size/2)]-squared_integral_image[i1-1][j1+(window_y_size/2)]
			-squared_integral_image[i1+(window_y_size/2)][j1-1]+squared_integral_image[i1-1][j1-1];

			result[1]=(s_f[1]-((f[1]*f[1])/small_window_size))/small_window_size;


			f[2]=integral_image[i2+(window_x_size/2)][j2+(window_y_size/2)]-integral_image[i2-1][j2+(window_y_size/2)]
			-integral_image[i2+(window_y_size/2)][j2-1]+integral_image[i2-1][j2-1];

			s_f[2]=squared_integral_image[i2+(window_x_size/2)][j2+(window_y_size/2)]-squared_integral_image[i2-1][j2+(window_y_size/2)]
			-squared_integral_image[i2+(window_y_size/2)][j2-1]+squared_integral_image[i2-1][j2-1];

			result[2]=(s_f[2]-((f[2]*f[2])/small_window_size))/small_window_size;


			f[3]=integral_image[i3+(window_x_size/2)][j3+(window_y_size/2)]-integral_image[i3-1][j3+(window_y_size/2)]
			-integral_image[i3+(window_y_size/2)][j3-1]+integral_image[i3-1][j3-1];

			s_f[3]=squared_integral_image[i3+(window_x_size/2)][j3+(window_y_size/2)]-squared_integral_image[i3-1][j3+(window_y_size/2)]
			-squared_integral_image[i3+(window_y_size/2)][j3-1]+squared_integral_image[i3-1][j3-1];

			result[3]=(s_f[3]-((f[3]*f[3])/small_window_size))/small_window_size;

			// cout<<"s_f[3]: "<<s_f[3]<<endl;
			// cout<<"f[3]: "<<	f[3]<<endl;
			// cout<<"small_window_size: "<<small_window_size<<endl; 
			// cout<<"result[3]: "<<result[3]<<endl;			
			// cout<<"f[3]*f[3]: "<<f[3]*f[3]<<endl;
			// cout<<"f[3]*f[3])/small_window_size: "<<(f[3]*f[3])/small_window_size<<endl;
			/*
				The process finding the smallest value between f[0] and f[3].
			*/
			int final=9999;
			for(int l=0;l<4;l++){
				// cout<<"result: "<<result[l]<<endl;
				// if(final>f[l]){
				// 	final=f[l];
				// }
				if(final>result[l]){
					final=result[l];
				}
			}
			cout<<"final: "<<final<<endl;
			
			// getchar();

			Mpixel(image2,i,j)=final;
			//////////////////////////// 

		}
	}

}

void Kuwahara_Filter_Gray_With_Sum_Table(Mat source_image, Mat output_image, double** integral_image, double** squared_integral_image, int window_size){
	//This is function for applying the Kuwahara filter to output_image, it must require summed-area table.

	//Mat source_image is the source image with gray scale. 
	//Mat output_image is the final output which is applied to Kuwahara filter in gray scale image.
	//double** integral_image is 2-dimensional array which stores the summed-area table.
	//int window_size is the size of window it is big more blurly when the value is high.

	int picture_x_size=source_image.cols;//the x-axis length of source image.
	int picture_y_size=source_image.rows;//the y-axis length of source image.
	

	for(int i=0;i<=source_image.cols-(window_size);i++){
		for(int j=0;j<=source_image.rows-(window_size);j++){
			double f[4];
			double s_f[4];
			double result[4];
			int small_window_size=(window_size/2+1)*(window_size/2+1);

			int i_col[4];
			int i_row[4];

			i_col[0]=i;
			i_row[0]=j;

			i_col[1]=i+(window_size/2);
			i_row[1]=j;

			i_col[2]=i;
			i_row[2]=j+(window_size/2);

			i_col[3]=i+(window_size/2);
			i_row[3]=j+(window_size/2);

			double mean_fa[4];

			
			for(int a=0;a<4;a++){
				f[a]=integral_image[i_col[a]+(window_size/2)+1][i_row[a]+(window_size/2)+1]-integral_image[i_col[a]][i_row[a]+(window_size/2)+1]
				-integral_image[i_col[a]+(window_size/2)+1][i_row[a]]+integral_image[i_col[a]][i_row[a]];

				s_f[a]=squared_integral_image[i_col[a]+(window_size/2)+1][i_row[a]+(window_size/2)+1]-squared_integral_image[i_col[a]][i_row[a]+(window_size/2)+1]
				-squared_integral_image[i_col[a]+(window_size/2)+1][i_row[a]]+squared_integral_image[i_col[a]][i_row[a]];

				mean_fa[a]=f[a]/small_window_size;
				// cout<<"mean_fa["<<a<<"]: "<<mean_fa[a]<<endl;
			}
			for(int a=0;a<4;a++){
				result[a]=(s_f[a]-((f[a]*f[a])/small_window_size))/small_window_size;
			}


			double final=9999;
			int min_index=0;
			for(int l=0;l<4;l++){
				if(final>result[l]){
					final=result[l];
					min_index=l;
				}
			}
			Mpixel(output_image,i+window_size/2,j+window_size/2)=(int)mean_fa[min_index];
		}
	}

}

void Kuwahara_Filter_Color_With_Sum_Table(Mat3b source_image, Mat3b output_image, double **integral_image_B, double **integral_image_G, double **integral_image_R, double **squared_integral_image_B, double **squared_integral_image_G, double **squared_integral_image_R, int window_size){
	//This is function for applying the Kuwahara filter to output_image, it must require summed-area table.

	//Mat source_image is the source image with gray scale. 
	//Mat output_image is the final output which is applied to Kuwahara filter in gray scale image.
	//double** integral_image is 2-dimensional array which stores the summed-area table.
	//int window_size is the size of window it is big more blurly when the value is high.

	int picture_x_size=source_image.cols;//the x-axis length of source image.
	int picture_y_size=source_image.rows;//the y-axis length of source image.
	

	for(int i=0;i<=source_image.cols-(window_size);i++){
		for(int j=0;j<=source_image.rows-(window_size);j++){
			double f_B[4];
			double f_G[4];
			double f_R[4];
			double s_f_B[4];
			double s_f_G[4];
			double s_f_R[4];
			double result_B[4];
			double result_G[4];
			double result_R[4];
			int small_window_size=(window_size/2+1)*(window_size/2+1);

			int i_col[4];
			int i_row[4];

			i_col[0]=i;
			i_row[0]=j;

			i_col[1]=i+(window_size/2);
			i_row[1]=j;

			i_col[2]=i;
			i_row[2]=j+(window_size/2);

			i_col[3]=i+(window_size/2);
			i_row[3]=j+(window_size/2);

			double mean_fa_B[4];
			double mean_fa_G[4];
			double mean_fa_R[4];

			for(int a=0;a<4;a++){
				f_B[a]=integral_image_B[i_col[a]+(window_size/2)+1][i_row[a]+(window_size/2)+1]-integral_image_B[i_col[a]][i_row[a]+(window_size/2)+1]
				-integral_image_B[i_col[a]+(window_size/2)+1][i_row[a]]+integral_image_B[i_col[a]][i_row[a]];

				s_f_B[a]=squared_integral_image_B[i_col[a]+(window_size/2)+1][i_row[a]+(window_size/2)+1]-squared_integral_image_B[i_col[a]][i_row[a]+(window_size/2)+1]
				-squared_integral_image_B[i_col[a]+(window_size/2)+1][i_row[a]]+squared_integral_image_B[i_col[a]][i_row[a]];

				mean_fa_B[a]=f_B[a]/small_window_size;
				// cout<<"mean_fa_B["<<a<<"]: "<<mean_fa_B[a]<<endl;
			}

			for(int a=0;a<4;a++){
				f_G[a]=integral_image_G[i_col[a]+(window_size/2)+1][i_row[a]+(window_size/2)+1]-integral_image_G[i_col[a]][i_row[a]+(window_size/2)+1]
				-integral_image_G[i_col[a]+(window_size/2)+1][i_row[a]]+integral_image_G[i_col[a]][i_row[a]];

				s_f_G[a]=squared_integral_image_G[i_col[a]+(window_size/2)+1][i_row[a]+(window_size/2)+1]-squared_integral_image_G[i_col[a]][i_row[a]+(window_size/2)+1]
				-squared_integral_image_G[i_col[a]+(window_size/2)+1][i_row[a]]+squared_integral_image_G[i_col[a]][i_row[a]];

				mean_fa_G[a]=f_G[a]/small_window_size;
				// cout<<"mean_fa_G["<<a<<"]: "<<mean_fa_G[a]<<endl;
			}
			
			for(int a=0;a<4;a++){
				f_R[a]=integral_image_R[i_col[a]+(window_size/2)+1][i_row[a]+(window_size/2)+1]-integral_image_R[i_col[a]][i_row[a]+(window_size/2)+1]
				-integral_image_R[i_col[a]+(window_size/2)+1][i_row[a]]+integral_image_R[i_col[a]][i_row[a]];

				s_f_R[a]=squared_integral_image_R[i_col[a]+(window_size/2)+1][i_row[a]+(window_size/2)+1]-squared_integral_image_R[i_col[a]][i_row[a]+(window_size/2)+1]
				-squared_integral_image_R[i_col[a]+(window_size/2)+1][i_row[a]]+squared_integral_image_R[i_col[a]][i_row[a]];

				mean_fa_R[a]=f_R[a]/small_window_size;
				// cout<<"mean_fa_R["<<a<<"]: "<<mean_fa_R[a]<<endl;
			}

			for(int a=0;a<4;a++){
				result_B[a]=(s_f_B[a]-((f_B[a]*f_B[a])/small_window_size))/small_window_size;
			}

			for(int a=0;a<4;a++){
				result_G[a]=(s_f_G[a]-((f_G[a]*f_G[a])/small_window_size))/small_window_size;
			}

			for(int a=0;a<4;a++){
				result_R[a]=(s_f_R[a]-((f_R[a]*f_R[a])/small_window_size))/small_window_size;
			}


			double final=9999;
			int min_index=0;
			for(int l=0;l<4;l++){
				if(final>result_B[l]){
					final=result_B[l];
					min_index=l;
				}
			}

			pixelB(output_image,i+window_size/2,j+window_size/2)=(int)mean_fa_B[min_index];

			final=9999;
			min_index=0;
			for(int l=0;l<4;l++){
				if(final>result_G[l]){
					final=result_G[l];
					min_index=l;
				}
			}

			pixelG(output_image,i+window_size/2,j+window_size/2)=(int)mean_fa_G[min_index];

			final=9999;
			min_index=0;

			for(int l=0;l<4;l++){
				if(final>result_R[l]){
					final=result_R[l];
					min_index=l;
				}
			}

			pixelR(output_image,i+window_size/2,j+window_size/2)=(int)mean_fa_R[min_index];
		}
	}

}




void Filter_Gray_Integral4(Mat image1, Mat image2, double** integral_image,double** squared_integral_image, int window_size){
	//This is function for applying the Kuwahara filter to image2, it must require summed-area table.

	//Mat image1 is the source image with gray scale. 
	//Mat image2 is the final output which is applied to Kuwahara filter in gray scale image.
	//double** integral_image is 2-dimensional array which stores the summed-area table.
	//int window_size is the size of window it is big more blurly when the value is high.

	int picture_x_size=image1.cols;//the x-axis length of source image.
	int picture_y_size=image1.rows;//the y-axis length of source image.
	int window_x_size=window_size;//the x-axis length of window.
	int window_y_size=window_size;//the y-axis length of window.


	for(int i=(window_x_size/2);i<picture_x_size-(window_x_size/2+1);i++){
		for(int j=(window_y_size/2);j<picture_y_size-(window_y_size/2+1);j++){
			float f[4];
			float s_f[4];
			float result[4];
			int small_window_size=(window_x_size/2+1)*(window_y_size/2+1);

			int i0=i;
			int j0=j;

			int i1=i+(window_x_size/2);
			int j1=j;

			int i2=i;
			int j2=j+(window_y_size/2);

			int i3=i+(window_x_size/2);
			int j3=j+(window_y_size/2);

			float mean_f[4];

			cout<<"small_window_size: "<<small_window_size<<endl;

			mean_f[0]=0;
			for(int l=i-(window_x_size/2);l<=i;l++){
				for(int k=j-(window_y_size/2);k<=j;k++){
					mean_f[0]=mean_f[0]+(float)Mpixel(image1,l,k);
				}
			}
			mean_f[0]=mean_f[0]/small_window_size;

			cout<<"mean_f[0]: "<<mean_f[0]<<endl;
			
			//the section of f[1]
			//  | |*|*|        
			//	| |*|*|
			//	| | | |
			mean_f[1]=0;
			for(int l=i;l<=i+(window_x_size/2);l++){
				for(int k=j-(window_y_size/2);k<=j;k++){
					mean_f[1]=mean_f[1]+(float)Mpixel(image1,l,k);
				}
			}
			mean_f[1]=mean_f[1]/small_window_size;

			cout<<"mean_f[1]: "<<mean_f[1]<<endl;

			//the section of f[2]
			//  | | | |        
			//	|*|*| |
			//	|*|*| |
			mean_f[2]=0;
			for(int l=i-(window_x_size/2);l<=i;l++){
				for(int k=j;k<=j+(window_y_size/2);k++){
					mean_f[2]=mean_f[2]+(float)Mpixel(image1,l,k);
				}
			}
			mean_f[2]=mean_f[2]/small_window_size;


			//the section of f[3]
			//  | | | |        
			//	| |*|*|
			//	| |*|*|
			mean_f[3]=0;
			for(int l=i;l<=i+(window_x_size/2);l++){
				for(int k=j;k<=j+(window_x_size/2);k++){
					mean_f[3]=mean_f[3]+(float)Mpixel(image1,l,k);
				}
			}

			mean_f[3]=mean_f[3]/small_window_size;





			//the section of f[0]
			//  |*|*| |        
			//	|*|*| |
			//	| | | |

			float temp=0;
			result[0]=0;


			for(int l=i-(window_x_size/2);l<=i;l++){
				for(int k=j-(window_y_size/2);k<=j;k++){
					temp=mean_f[0]-(float)Mpixel(image1,l,k);
					temp=temp*temp;
					result[0]=result[0]+temp;
				}
			}

			result[0]=result[0]/small_window_size;


			
			//the section of f[1]
			//  | |*|*|        
			//	| |*|*|
			//	| | | |
			result[1]=0;
			for(int l=i;l<=i+(window_x_size/2);l++){
				for(int k=j-(window_y_size/2);k<=j;k++){
					temp=mean_f[1]-(float)Mpixel(image1,l,k);
					temp=temp*temp;
					result[1]=result[1]+temp;
				}
			}
			result[1]=result[1]/small_window_size;


			//the section of f[2]
			//  | | | |        
			//	|*|*| |
			//	|*|*| |
			result[2]=0;
			for(int l=i-(window_x_size/2);l<=i;l++){
				for(int k=j;k<=j+(window_y_size/2);k++){
					temp=mean_f[2]-(float)Mpixel(image1,l,k);
					temp=temp*temp;
					result[2]=result[2]+temp;
				}
			}
			result[2]=result[2]/small_window_size;


			//the section of f[3]
			//  | | | |        
			//	| |*|*|
			//	| |*|*|
			result[3]=0;
			for(int l=i;l<=i+(window_x_size/2);l++){
				for(int k=j;k<=j+(window_x_size/2);k++){
					temp=mean_f[3]-(float)Mpixel(image1,l,k);
					temp=temp*temp;
					// cout<<"temp: "<<temp<<endl;
					result[3]=result[3]+temp;
				}
			}
			// getchar();

			result[3]=result[3]/small_window_size;












			// cout<<"s_f[3]: "<<s_f[3]<<endl;
			// cout<<"f[3]: "<<	f[3]<<endl;
			// cout<<"small_window_size: "<<small_window_size<<endl; 
			cout<<"result[0]: "<<result[0]<<endl;			
			cout<<"result[1]: "<<result[1]<<endl;			
			cout<<"result[2]: "<<result[2]<<endl;			
			cout<<"result[3]: "<<result[3]<<endl;
			cout<<"i: "<<i<<endl;
			cout<<"j: "<<j<<endl;			
			// cout<<"f[3]*f[3]: "<<f[3]*f[3]<<endl;
			// cout<<"f[3]*f[3])/small_window_size: "<<(f[3]*f[3])/small_window_size<<endl;
			/*
				The process finding the smallest value between f[0] and f[3].
			*/
			float final=9999;
			for(int l=0;l<4;l++){
				// cout<<"result: "<<result[l]<<endl;
				// if(final>f[l]){
				// 	final=f[l];
				// }
				if(final>result[l]){
					final=result[l];
				}
			}

			for(int l=0;l<4;l++){
				if(final==result[l]){
					Mpixel(image2,i,j)=mean_f[l];
					break;
				}
			}
			// cout<<"final: "<<final<<endl;
			// if(final>255){
			// 	// getchar();
			// 	Mpixel(image2,i,j)=(int)Mpixel(image1,i,j);
			// }else{
			// Mpixel(image2,i,j)=(int)final;	
			// }

			
			//////////////////////////// 

		}
	}

}

void Kuwahara_Filter_Gray_Without_Sum_Table(Mat source_image, Mat output_image, int window_size){
	//This is function for applying the Kuwahara filter to output_image, it must require summed-area table.

	//Mat source_image is the source image with gray scale. 
	//Mat output_image is the final output which is applied to Kuwahara filter in gray scale image.
	//double** integral_image is 2-dimensional array which stores the summed-area table.
	//int window_size is the size of window it is big more blurly when the value is high.

	int picture_x_size=source_image.cols;//the x-axis length of source image.
	int picture_y_size=source_image.rows;//the y-axis length of source image.


	
	for(int i=0;i<=source_image.cols-(window_size);i++){
		for(int j=0;j<=source_image.rows-(window_size);j++){
			double f[4];
			double s_f[4];
			double result[4];
			int small_window_size=(window_size/2+1)*(window_size/2+1);

			int i_col[4];
			int i_row[4];

			i_col[0]=i;
			i_row[0]=j;

			i_col[1]=i+(window_size/2);
			i_row[1]=j;

			i_col[2]=i;
			i_row[2]=j+(window_size/2);

			i_col[3]=i+(window_size/2);
			i_row[3]=j+(window_size/2);

			double mean_f[4];

			// cout<<"small_window_size: "<<small_window_size<<endl;

			mean_f[0]=0;
			for(int l=i;l<=i+(window_size/2);l++){
				for(int k=j;k<=j+(window_size/2);k++){
					mean_f[0]=mean_f[0]+(float)Mpixel(source_image,l,k);
				}
			}
			mean_f[0]=mean_f[0]/small_window_size;

			// cout<<"mean_f[0]: "<<mean_f[0]<<endl;
			
			//the section of f[1]
			//  | |*|*|        
			//	| |*|*|
			//	| | | |
			mean_f[1]=0;
			for(int l=i+(window_size/2);l<i+window_size;l++){
				for(int k=j;k<=j+(window_size/2);k++){
					mean_f[1]=mean_f[1]+(float)Mpixel(source_image,l,k);
				}
			}
			mean_f[1]=mean_f[1]/small_window_size;

			// cout<<"mean_f[1]: "<<mean_f[1]<<endl;

			//the section of f[2]
			//  | | | |        
			//	|*|*| |
			//	|*|*| |

			mean_f[2]=0;
			for(int l=i;l<=i+(window_size/2);l++){
				for(int k=j+(window_size/2);k<j+window_size;k++){
					mean_f[2]=mean_f[2]+(float)Mpixel(source_image,l,k);
				}
			}
			mean_f[2]=mean_f[2]/small_window_size;
			// cout<<"mean_f[2]: "<<mean_f[2]<<endl;


			//the section of f[3]
			//  | | | |        
			//	| |*|*|
			//	| |*|*|

			mean_f[3]=0;
			for(int l=i+(window_size/2);l<i+window_size;l++){
				for(int k=j+(window_size/2);k<j+window_size;k++){
					mean_f[3]=mean_f[3]+(float)Mpixel(source_image,l,k);
				}
			}

			mean_f[3]=mean_f[3]/small_window_size;
			// cout<<"mean_f[3]: "<<mean_f[3]<<endl;
			
			

			double temp=0;
			result[0]=0;
			for(int l=i;l<=i+(window_size/2);l++){
				for(int k=j;k<=j+(window_size/2);k++){
					temp=mean_f[0]-(double)Mpixel(source_image,l,k);
					temp=temp*temp;
					result[0]=result[0]+temp;
				}
			}

			result[0]=result[0]/small_window_size;


			
			//the section of f[1]
			//  | |*|*|        
			//	| |*|*|
			//	| | | |
			result[1]=0;
			for(int l=i+(window_size/2);l<i+window_size;l++){
				for(int k=j;k<=j+(window_size/2);k++){
					temp=mean_f[1]-(double)Mpixel(source_image,l,k);
					temp=temp*temp;
					result[1]=result[1]+temp;
				}
			}
			result[1]=result[1]/small_window_size;


			//the section of f[2]
			//  | | | |        
			//	|*|*| |
			//	|*|*| |
			result[2]=0;
			for(int l=i;l<=i+(window_size/2);l++){
				for(int k=j+(window_size/2);k<j+window_size;k++){
					temp=mean_f[2]-(double)Mpixel(source_image,l,k);
					temp=temp*temp;
					result[2]=result[2]+temp;
				}
			}
			result[2]=result[2]/small_window_size;


			//the section of f[3]
			//  | | | |        
			//	| |*|*|
			//	| |*|*|
			result[3]=0;
			for(int l=i+(window_size/2);l<i+window_size;l++){
				for(int k=j+(window_size/2);k<j+window_size;k++){
					temp=mean_f[3]-(double)Mpixel(source_image,l,k);
					temp=temp*temp;
					// cout<<"temp: "<<temp<<endl;
					result[3]=result[3]+temp;
				}
			}
			// getchar();

			result[3]=result[3]/small_window_size;

			double final=9999;
			int min_index=0;
			for(int l=0;l<4;l++){
				if(final>result[l]){
					final=result[l];
					min_index=l;
				}
			}
			Mpixel(output_image,i+window_size/2,j+window_size/2)=(int)mean_f[min_index];
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

		gray_image2=Mat::ones(gray_image2.size(),CV_LOAD_IMAGE_GRAYSCALE);//initialize the value of output metrices to zero
		Mat gray_image3=Mat::ones(gray_image2.size(),CV_LOAD_IMAGE_GRAYSCALE);//initialize the value of output metrices to zero
		gray_image2=gray_image2+gray_image3;



		// //////////////////////////////////////
		// //Kuwahara filter using Summed-table//
		// //////////////////////////////////////

		// /*Memory Allocation*/
		// double** integral_image=new double*[gray_image2.cols+1];
		// double** squared_integral_image=new double*[gray_image2.cols+1];
		// for(int i = 0; i < gray_image2.cols+1; ++i){
		// 	integral_image[i] = new double[gray_image2.rows+1];
		// 	squared_integral_image[i] = new double[gray_image2.rows+1];
		// }
		// /*********************/

		// Integral_Gray_Initialize(gray_image2,integral_image,squared_integral_image);//create summed-table to integral_image array.
		// Integral_test(gray_image2,integral_image,squared_integral_image);//create summed-table to integral_image array.
		
		// Kuwahara_Filter_Gray_Without_Sum_Table(gray_image2,output,7);//Applying kuwahara filter to output using integral_image.
		
		// /*Memory deallocation*/
		// for(int i = 0; i < gray_image1.cols+1; ++i) {
		// 	delete [] integral_image[i];
		// 	delete [] squared_integral_image[i];
		// }
		// delete [] integral_image;
		// delete [] squared_integral_image;
		// /***************/
		
		// ///////////////////////////////////////





		//////////////////////////////////////
		//Kuwahara filter using Summed-table with colour//
		//////////////////////////////////////

		/*Memory Allocation*/
		double** integral_image_R=new double*[gray_image2.cols+1];
		double** squared_integral_image_R=new double*[gray_image2.cols+1];
		double** integral_image_G=new double*[gray_image2.cols+1];
		double** squared_integral_image_G=new double*[gray_image2.cols+1];
		double** integral_image_B=new double*[gray_image2.cols+1];
		double** squared_integral_image_B=new double*[gray_image2.cols+1];

		for(int i = 0; i < gray_image2.cols+1; ++i){
			integral_image_R[i] = new double[gray_image2.rows+1];
			squared_integral_image_R[i] = new double[gray_image2.rows+1];
			integral_image_G[i] = new double[gray_image2.rows+1];
			squared_integral_image_G[i] = new double[gray_image2.rows+1];
			integral_image_B[i] = new double[gray_image2.rows+1];
			squared_integral_image_B[i] = new double[gray_image2.rows+1];
		}
		/*********************/
		Mat3b color_output;
		color_output=Mat::zeros(gray_image2.size(),CV_8UC3);//initialize the value of final_output metrices to zero
		Integral_Color_Initialize(image2,integral_image_B,integral_image_G,integral_image_R,squared_integral_image_B,squared_integral_image_G,squared_integral_image_R);//create summed-table to integral_image array.
		Kuwahara_Filter_Color_With_Sum_Table(image2,color_output,integral_image_B,integral_image_G,integral_image_R,squared_integral_image_B,squared_integral_image_G,squared_integral_image_R,9);//Applying kuwahara filter to output using integral_image.
		
		/*Memory deallocation*/
		for(int i = 0; i < gray_image1.cols+1; ++i) {
			delete [] integral_image_R[i];
			delete [] squared_integral_image_R[i];
			delete [] integral_image_G[i];
			delete [] squared_integral_image_G[i];
			delete [] integral_image_B[i];
			delete [] squared_integral_image_B[i];
		}
		delete [] integral_image_R;
		delete [] squared_integral_image_R;
		delete [] integral_image_G;
		delete [] squared_integral_image_G;
		delete [] integral_image_B;
		delete [] squared_integral_image_B;
		/***************/
		
		///////////////////////////////////////



		////////////////////////////////////////
		//Kuwahara filter without Summed-table//
		////////////////////////////////////////
		
		// Filter_Gray(gray_image2,output,15);
		
		///////////////////////////////////////
		namedWindow("image2",WINDOW_NORMAL);
		namedWindow("color_output",WINDOW_NORMAL);
		resizeWindow("color_output", 1200,1200);
		resizeWindow("image2", 1200,1200);
		imshow("color_output" ,color_output);
		imshow("image2" ,image2);
		// imwrite( "lee_eun.jpg", color_output );
		waitKey(0);
		exit(0);

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


		// Scalar color=CV_RGB(255,0,0);
		// std::vector<vector<Point>>contours;
		// findContours(output,contours,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE);

        
  //       int largestcontour=FindTheLargestContour(contours);
  //       Mat drawing=Mat::zeros(output.size(),CV_8UC3);
		// drawContours(drawing,contours,largestcontour,color,2,8);
		/*******************************************************/

    	/*final_output show*/
    	imshow("final_output" ,final_output);
    	imshow("filtered image" ,output);
    	imshow("gray_image2" ,gray_image2);
    	// imshow("drawing" ,drawing);



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
	   	int j=0;
	   	while (1){
	   		system_clock::time_point start = system_clock::now();//start clock
	   		
	   		Mat3b image2;
	   		Mat gray_image2;
	   		Mat output;
	   		Mat output1;
	   		Mat output2;
	   		Mat3b final_output;
	   		cap >> image2;
			
			cvtColor(image2,gray_image2, CV_BGR2GRAY);//copy camera color image to gray scale
			
			if(!image2.data){printf("Could not open the file\n"); exit(0);}
			output=Mat::zeros(gray_image2.size(),CV_LOAD_IMAGE_GRAYSCALE);//initialize the value of output metrices to zero
			output1=Mat::zeros(gray_image1.size(),CV_LOAD_IMAGE_GRAYSCALE);//initialize the value of output metrices to zero
			output2=Mat::zeros(gray_image2.size(),CV_LOAD_IMAGE_GRAYSCALE);//initialize the value of output metrices to zero
    		final_output=Mat::zeros(gray_image2.size(),CV_8UC3);//initialize the value of final_output metrices to zero

    		//////////////////////////////////////
    		//Kuwahara filter using Summed-table//
    		//////////////////////////////////////

			/*Memory Allocation*/
			double** integral_image1=new double*[gray_image1.cols+1];
			double** squared_integral_image1=new double*[gray_image1.cols+1];

			double** integral_image2=new double*[gray_image2.cols+1];
			double** squared_integral_image2=new double*[gray_image2.cols+1];
			for(int i = 0; i < gray_image2.cols+1; ++i){
				integral_image1[i] = new double[gray_image1.rows+1];
				squared_integral_image1[i] = new double[gray_image1.rows+1];

				integral_image2[i] = new double[gray_image2.rows+1];
				squared_integral_image2[i] = new double[gray_image2.rows+1];
			}
    		/*********************/

    		Integral_Gray_Initialize(gray_image2,integral_image2,squared_integral_image2);//create summed-table to integral_image array.
    		Kuwahara_Filter_Gray_With_Sum_Table(gray_image2,output2,integral_image2,squared_integral_image2,11);//Applying kuwahara filter to output using integral_image.
    		Integral_Gray_Initialize(gray_image1,integral_image1,squared_integral_image1);//create summed-table to integral_image array.
    		Kuwahara_Filter_Gray_With_Sum_Table(gray_image1,output1,integral_image1,squared_integral_image1,11);//Applying kuwahara filter to output using integral_image.
		
			/*Memory deallocation*/
			for(int i = 0; i < gray_image1.cols+1; ++i) {
				delete [] integral_image1[i];
				delete [] squared_integral_image1[i];

				delete [] integral_image2[i];
				delete [] squared_integral_image2[i];
			}
			delete [] integral_image1;
			delete [] squared_integral_image1;
			delete [] integral_image2;
			delete [] squared_integral_image2;
			/***************/
			
			///////////////////////////////////////



			////////////////////////////////////////
    		//Kuwahara filter without Summed-table//
    		////////////////////////////////////////
    		
    		// Filter_Gray(gray_image2,output,5);

    		
    		///////////////////////////////////////
    		

			/*subtraction process between The first image and the second image*/
			output=gray_image1-gray_image2;
			/********************************/
			
			/***Cropping Object by outline of object***/
			Mat temp_window=Mat::zeros(image2.size(),CV_LOAD_IMAGE_GRAYSCALE);//gray case
			Image_stitching(gray_image2,output,temp_window);
			/******************************************/

			/***ROI Section***/
			int small_x,small_y,large_x,large_y;//size of
			Determining_ROI_Size(temp_window,&small_x,&small_y,&large_x,&large_y);
			Mat ROI;
			if((large_x==-9999)||(large_y==-9999)||(small_x==9999)||(small_y==9999)){
				ROI=Mat::zeros(gray_image2.size(),CV_LOAD_IMAGE_GRAYSCALE);//Default size
			}else{

				Rect faces;
				faces.width=large_x;
				faces.height=large_y;
				ROI = gray_image2(faces);
				for(int x=small_x; x<large_x; x++){
					for(int y=small_y; y<large_y; y++){
						/*in case of color*/
						// pixelB(ROI,x,y)=pixelB(image2,x,y);
						// pixelG(ROI,x,y)=pixelG(image2,x,y);
						// pixelR(ROI,x,y)=pixelR(image2,x,y);

						Mpixel(ROI,x-small_x,y-small_y)=Mpixel(gray_image2,x,y);					
					}
				}
			}
			
			/***final output color***/
			Grey_to_Color(image2,temp_window,final_output);

			/*******************************************************/

	    	/*output Showing Section*/
	    	imshow("color_output" ,final_output);
	    	imshow("ROI" ,ROI);
	    	imshow("temp_window" ,temp_window);
	    	imshow("Gray_output" ,output);

	    	// imshow("image1" ,output1);
	    	// imshow("image2" ,output2);
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
	     	
	     	/***Delay for image1***/
	     	// j++;
	     	// if(j==3){
	     	// 	image1=image2;
	     	// 	cvtColor(image1,gray_image1, CV_BGR2GRAY);//copy camera color image to gray scale
	     	// 	j=0;
	     	// }
	     	/**********************/
	     	image1=image2;
	     	cvtColor(image1,gray_image1, CV_BGR2GRAY);//copy camera color image to gray scale
	    }
	    /************/
	}
}