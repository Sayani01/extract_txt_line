#include <stdio.h>
#include<conio.h>
#include<iostream>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "NICK.h"
#include "GaussBlurRotated.h"

using namespace std;
using namespace cv;

void blurROI(int **img, int h, int w, int theta) {
	//int **binimg;
	double **org_img, weight[3][3];
	org_img = new double*[h + 2];
	for (int i = 0; i < h + 2; i++) {
		org_img[i] = new double[w + 2];
		for (int j = 0; j < w + 2; j++) {
			if (i > 0 && i <= h && j > 0 && j <= w) 
				org_img[i][j] = img[i - 1][j - 1]; 
			if (i == 0 ) {
				if (j > 0 && j <= w) 
					org_img[0][j] = img[0][j - 1];
				org_img[0][0] = img[0][0];
				org_img[0][w + 1] = img[0][w - 1];
			}
			if (i == h + 1) {
				if (j > 0 && j <= w)
					org_img[h + 1][j] = img[h - 1][j - 1];
				org_img[h + 1][0] = img[h - 1][0];
				org_img[h + 1][w + 1] = img[h - 1][w - 1];
			}
			if (j == 0) {
				if (i > 0 && i <= h) 
					org_img[i][0] = img[i - 1][0];
			}
			if (j == w + 1) {
				if (i > 0 && i <= h) 
					org_img[i][0] = img[i - 1][w - 1];
			}
		}
	}

	double sigmax = 1.0, sigmay = 1.0;
	theta *= (CV_PI / 180);
	/*double a = cos(theta)*cos(theta) / (2 * (sigmax)*(sigmax)) + sin(theta)*sin(theta) / (2 * (sigmay)*(sigmay));
	double c = sin(theta)*sin(theta) / (2 * (sigmax)*(sigmax)) + cos(theta)*cos(theta) / (2 * (sigmay)*(sigmay));*/

	double x, y, u;
	double wt[3][3];
	double wtSum = 0;
	//Set weight matrix...
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			x = i*cos(theta) - j*sin(theta);
			y = i*sin(theta) + j*cos(theta);
			u = (x / sigmax) * (x / sigmax) + (y / sigmay) * (y / sigmay);
			weight[i + 1][j + 1] = exp(-u / 2);
			wtSum += weight[i + 1][j + 1];
		}
	}
	
	for (int i = 1; i <= h; i++) {
		for (int j = 1; j <= w; j++) {
			double roi[3][3];
			double sum = 0;
			//set 3x3 roi...
			for (int k = -1; k <= 1; k++) {
				for (int l = -1; l <= 1; l++) {
					roi[k + 1][l + 1] = org_img[i + k][j + l] ;
				}
			}
			//multiply by weights....
			sum = 0;
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					roi[k][l] *= weight[k][l] / wtSum;
					sum += roi[k][l];
				}
			}
			//Set gaussian value to center....
			roi[1][1] = sum;
			for (int k = -1; k <= 1; k++) {
				for (int l = -1; l <= 1; l++) {
					org_img[i + k][j + l] = roi[k + 1][l + 1];
				}
			}
		}
	}

	for (int i = 1; i <= h; i++) {
		for (int j = 1; j <= w; j++) {
			img[i - 1][j - 1] = int (255*org_img[i][j]);
			if (img[i - 1][j - 1] < 0) img[i - 1][j - 1] = 0;
			if (img[i - 1][j - 1] > 255) img[i - 1][j - 1] = 255;
		}
	}
	/*xc = center(1);
	yc = center(2);
	theta = (theta / 180)*pi;
	xm = (x - xc)*cos(theta) - (y - yc)*sin(theta);
	ym = (x - xc)*sin(theta) + (y - yc)*cos(theta);
	u = (xm / sigmax) ^ 2 + (ym / sigmay) ^ 2;
	val = offset + factor*exp(-u / 2);*/

	//Delete memory...
	for (int i = 0; i < h + 2; i++)
		delete[] org_img[i];
	delete[] org_img;

}