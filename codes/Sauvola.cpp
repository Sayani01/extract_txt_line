#include "Sauvola.h"
#include<iostream>
#include<math.h>

using namespace std;

double rowSum(double **mat, int R, int c_min,int c_max) {
	double sum = 0;
	for (int j = c_min; j <= c_max; j++)
		sum += mat[R][j];
	return(sum);
}

double colSum(double **mat, int C, int r_min, int r_max) {
	double sum = 0;
	for (int i = r_min; i <= r_max; i++)
		sum += mat[i][C];
	return(sum);
}

void normalizeImage(int **im, int r, int c, double **norm_im)
{
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			norm_im[i][j] = im[i][j] / 255.0;
}

void meanImage(double **mat, int r, int c, double **mean, int fdims) {
	int lim = fdims / 2 , NP = fdims*fdims;
	for (int i = lim; i < r - lim; i++) {
		for (int j = lim; j < c - lim; j++) {
			if (i == lim && j == lim) {  //first pixel...
				double sum = 0.0;
				for (int k = j - lim; k <= j + lim; k++)
					sum += colSum(mat, k, i - lim, i + lim);
				mean[i][j] = sum / NP;
			}

			if (j > lim)  //same row next column....
				mean[i][j] = (mean[i][j - 1] * NP - colSum(mat,  j-lim-1, i - lim, i + lim) + colSum(mat, j + lim, i - lim, i + lim)) / NP;

			if (j == lim && i>lim)  //new row...
				mean[i][j] = (mean[i - 1][j] * NP - rowSum(mat, i - lim - 1, j - lim, j + lim) + rowSum(mat, i + lim, j - lim, j + lim)) / NP;
		}
	}
}

void std_devImage(double **mat, double **M, int r, int c, double **std_dev, int fdims) {
	int lim = fdims / 2, NP = fdims*fdims;
	double **D;
	D = new double*[r];
	for (int i = 0; i < r; i++) {
		D[i] = new double[c];
		for (int j = 0; j < c; j++)
			D[i][j] = pow((mat[i][j] - M[i][j]),2);			
	}
	meanImage(D, r, c, std_dev, fdims);
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			std_dev[i][j] = sqrt(std_dev[i][j]);
}

void Sauvola(int **im, int r, int c, int **bin_im) {
	int fdims = 19;
	double k = 0.5, R;
	double **norm_im,**mean, **std_dev, **T;

	norm_im = new double *[r];
	mean = new double*[r];
    std_dev = new double*[r];
	T = new double *[r];
	for (int i = 0; i < r; i++) {
		norm_im[i] = new double[c];
		mean[i] = new double[c];
		std_dev[i] = new double[c];
		T[i] = new double[c];
		for (int j = 0; j < c; j++) {
			norm_im[i][j] = 0.0;
			mean[i][j] = 0.0;
			std_dev[i][j] = 0.0;
			T[i][j] = 0.0;
		}
	}

	normalizeImage(im, r, c, norm_im);
	meanImage(norm_im, r, c, mean,fdims);
	std_devImage(norm_im, mean, r, c, std_dev, fdims);

	R = std_dev[0][0];
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			if (std_dev[i][j] > R)
				R = std_dev[i][j];
	
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			//Sauvola...
			T[i][j] = mean[i][j] * (1.0 + k*((std_dev[i][j] / R) - 1.0 ));
			//Niblack...
			//T[i][j] = mean[i][j] + k*std_dev[i][j];
			if (norm_im[i][j] > T[i][j])
				bin_im[i][j] = 255;
			else
				bin_im[i][j] = 0;
		}
	}
	
	for (int i = 0; i < r; i++) {
		delete[] norm_im[i];
		delete[] std_dev[i];
		delete[] mean[i];
		delete[] T[i];
	}
	delete[] norm_im;
	delete[] std_dev;
	delete[] mean;
	delete[] T;
}