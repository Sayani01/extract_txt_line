#include "NICK.h"
#include<iostream>
#include<math.h>

using namespace std;

double rowSum(double **mat, int R, int c_min, int c_max);

double colSum(double **mat, int C, int r_min, int r_max);

void normalizeImage(int **im, int r, int c, double **norm_im);

void meanImage(double **mat, int r, int c, double **mean, int fdims);

void NICK_dev(double **mat, double **M, int r, int c, double **std_dev, int fdims) {
	int lim = fdims / 2, NP = fdims*fdims;
	double **D;
	D = new double*[r];
	for (int i = 0; i < r; i++) {
		D[i] = new double[c];
		for (int j = 0; j < c; j++)
			D[i][j] = pow((mat[i][j] - M[i][j]), 2);
	}
	meanImage(D, r, c, std_dev, fdims);
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			std_dev[i][j] = sqrt(std_dev[i][j] + (M[i][j] * M[i][j]));
}


void NICK(int **im, int r, int c, int **bin_im) {
	int fdims = 17;
	double k = -0.1, R;
	double **norm_im, **mean, **D, **T;

	norm_im = new double *[r];
	mean = new double*[r];
	D = new double*[r];
	T = new double *[r];
	for (int i = 0; i < r; i++) {
		norm_im[i] = new double[c];
		mean[i] = new double[c];
		D[i] = new double[c];
		T[i] = new double[c];
		for (int j = 0; j < c; j++) {
			norm_im[i][j] = 0.0;
			mean[i][j] = 0.0;
			D[i][j] = 0.0;
			T[i][j] = 0.0;
		}
	}

	normalizeImage(im, r, c, norm_im);
	meanImage(norm_im, r, c, mean, fdims);
	NICK_dev(norm_im, mean, r, c, D, fdims);

	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			T[i][j] = mean[i][j] + k*D[i][j];
			if (norm_im[i][j] > T[i][j])
				bin_im[i][j] = 255;
			else
				bin_im[i][j] = 0;
		}
	}

	for (int i = 0; i < r; i++) {
		delete[] norm_im[i];
		delete[] D[i];
		delete[] mean[i];
		delete[] T[i];
	}
	delete[] norm_im;
	delete[] D;
	delete[] mean;
	delete[] T;
}