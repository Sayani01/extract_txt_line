#include <stdio.h>
#include<conio.h>
#include<iostream>
#include "opencv2/imgproc/imgproc.hpp"
//#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "ConnectedComponent.h"

using namespace std;
using namespace cv;

int lblNum = 0 , rows = 0 , col = 0;

//find the minimum in an array...
int minimum(int a, int b, int c, int d) {
	int elements[4] = { a,b,c,d },min = 999999;
	for (int i = 0; i < 4; i++) {
		if (elements[i] == 0)
			continue;
		if (elements[i] < min)
			min = elements[i];
	}
	return(min);
}

//Relabel the connected components...
void fill(int **label, int I, int J, int replace, int current) {
	int flag = 0;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < col; j++) {
			if (i == I && j == J) {
				flag = 1;
				break;
			}
			if (label[i][j] == current)
				label[i][j] = replace;
		}
		if (flag)
			break;
	}
}

//Bound the components with recatngle...
void boundingRectangles(int **box, int **label, int img_r,int img_c, int lbl){
	lblNum = lbl, rows = img_r, col = img_c;
	for (int i = 0; i < lblNum; i++) {
		box[i][0] = i;
		box[i][1] = 0; //open/close component...
		box[i][2] = -1; //x-start
		box[i][3] = -1; //x-end
		box[i][4] = -1; //y-start
		box[i][5] = -1; //y-end
	}

	for (int j = 0; j < col ; j++) {
		int *flag = new int[lblNum];
		for (int k = 0; k < lblNum; k++)
			flag[k] = 0;
		for (int i = 0; i < rows; i++) {
			int pos = label[i][j];
			flag[pos] = 1;
			if (pos == 0 || j == col - 1) //last column or white pixels...
				continue;
			if (box[pos][1] == 0) { //new component...
				box[pos][2] = j;
				box[pos][1] = 1; //open component...
			}
		}
		if (j == col - 1) { //last column...
			for (int k = 0; k < lblNum; k++) {
				if (flag[k] == 1) {
					box[k][3] = j;
					box[k][1] = 0; //close component...
				}
			}
		}
		for (int k = 0; k < lblNum; k++) {
			if (flag[k] == 0 && box[k][1] == 1) { //component open but not encountered...
				box[k][3] = j - 1;
				box[k][1] = 0; //close component...
			}
		}
		delete[] flag;
	}

	for (int i = 0; i < rows; i++) {
		int *flag = new int[lblNum];
		for (int k = 0; k < lblNum; k++)
			flag[k] = 0;
		for (int j = 0; j < col; j++) {
			int pos = label[i][j];
			flag[pos] = 1;
			if (pos == 0 || i == 0 || i == rows - 1) //first row...
				continue;
			if (box[pos][1] == 0) { //new component...
				box[pos][4] = i;
				box[pos][1] = 1; //close component;
			}
		}
		if (i == rows - 1) { //last row...
			for (int k = 0; k < lblNum; k++) {
				if (flag[k] == 1) {
					box[k][5] = i;
					box[k][1] = 0; //close component...
				}
			}
		}
		for (int k = 0; k < lblNum; k++) {
			if (flag[k] == 0 && box[k][1] == 1) { //component open but not encountered
				box[k][5] = i - 1;
				box[k][1] = 0;  //close component...
			}
		}
		delete[] flag;
	}
	/*std::cout << "label on/off  x0   x1     y0      y1 \n";
	for (int i = 0; i < lblNum; i++) {
		std::cout << "  ";
		for (int j = 0; j < 6; j++) {
			std::cout << box[i][j] << "     ";
		}
		std::cout << "\n";
	}*/
}

//label the 8-connected components...
int labelling(int **img, int **label, int img_r,int img_c ) {
	lblNum = 0;
	rows = img_r, col = img_c;
	int W, NW, N, NE;
	int *labelcount; 
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < col; j++) {
			if (img[i][j] == 255)    //background...
				continue;
			else{
				if (i == 0 && j==0) { //first cell...
					W = 0;
					NW = 0;
					N = 0;
					NE = 0;
				}
				else if (i == 0 && j > 0) {  //first row...
					W = label[i][j - 1];
					NW = 0;
					N = 0;
					NE = 0;
				}
				else if (i > 0 && j == 0){ //first column...
					W = 0;
					NW = 0;
					N = label[i - 1][j];
					NE = label[i - 1][j + 1];

				}
				else if (i > 0 && j == col-1) { //last column...
					W = label[i][j - 1];
					NW = label[i-1][j - 1];
					N = label[i - 1][j];
					NE = 0;
				}
				else {
					W = label[i][j - 1];
					NW = label[i - 1][j - 1];
					N = label[i - 1][j];
					NE = label[i - 1][j + 1];
				}
				
				if(W==0 && NW==0 && N==0 && NE==0) //new component...
					label[i][j] = ++lblNum;
				else //label with the minimum value...
					label[i][j] = minimum(N, NW, NE, W);
				//re-label...
				if (W!=0 && W != label[i][j])
					fill(label, i, j, label[i][j], W);
				if (NW!=0 && NW != label[i][j])
					fill(label, i, j, label[i][j], NW);
				if (N!=0 && N != label[i][j])
					fill(label, i, j, label[i][j], N);
				if (NE!=0 && NE != label[i][j])
					fill(label, i, j, label[i][j], NE);				
			}			
		}
	}

	lblNum = lblNum + 1;
	labelcount = new int[lblNum];
	for (int i = 0; i < lblNum; i++)
		labelcount[i] = 0;
	int pos;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < col; j++) {
			if (label[i][j] < 0)
				std::cout << "here->" << i << "," << j << "\n";
			pos = label[i][j];
			labelcount[pos]++;   //count the no. of pixels in each label....
		}
	}

	/*for (int i = 0; i < lblNum; i++)
		std::cout << "label " << i << " = " << labelcount[i] << " pixels\n";*/

	//for (int i = 0; i <= labelNum; i++) {
	//	if(labelcount[i]<6){  //negligible component....
	//		for (int k = i; k < labelNum; k++) {
	//			labelcount[k] = labelcount[k + 1];
	//		}
	//		labelNum--;
	//		i--;
	//	}
	//}

	/*std::cout << "\nNo.of Labels = " << lblNum << "\n";
	for (int i = 0; i < lblNum; i++) {
		if (labelcount[i] != 0)
			std::cout << "label " << i << " = " << labelcount[i] << " pixels\n";
	}*/
	return(lblNum);
}

