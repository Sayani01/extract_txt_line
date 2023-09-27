#include "Otsu.h"
#include<iostream>
#include<math.h>

using namespace std;

int L = 100;

void normalizeImage(int **im, int r, int c, double **norm_im);
/*{
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			norm_im[i][j] = im[i][j] / 255.0; 
}*/

void histogram(int *level, double **norm_im, int **label , int r , int c) {
	int q;
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++) {
			if (norm_im[i][j] >= 1.0) {
				level[L - 1]++;
				label[i][j] = L;
				continue;
			}
			q = int(norm_im[i][j] * L);
			level[q]++;
			label[i][j] = q + 1;
		}
}

double betweenClassVariance(double *p, int k) {
	double w0 = 0.0, w1 = 0.0, mu0 = 0.0, mu1 = 0.0, var_b;
	for (int i = 0; i <= k; i++) {
		w0 += p[i];
		mu0 += (i + 1)*p[i];
	}
	mu0 /= w0;
	for (int i = k + 1; i < L; i++) {
		w1 += p[i];
		mu1 += (i + 1)*p[i];
	}
	mu1 /= w1;
	var_b = w0*w1*(mu1 - mu0)*(mu1 - mu0);
	return var_b;
}


void binary(double **norm_im, double T, int **bin_im, int r, int c) {
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			if (norm_im[i][j] > T)
				bin_im[i][j] = 255;
			else
				bin_im[i][j] = 0;
}

double Otsu(int **im , int r , int c , int **bin_im) {
	double var_b , T , argmax_var_b = 0.0 , thresh = 0.0;
	double **norm_im, *p;
	int **label, *level;
	
	level = new int[L];
	p = new double[L];
	for (int i = 0; i < L; i++)
		level[i] = 0;

	norm_im = new double*[r];
	label = new int*[r];
	for (int i = 0; i < r; i++) {
		norm_im[i] = new double[c];
		label[i] = new int[c];
		for (int j = 0; j < c; j++) {
			label[i][j] = -1;
		}
	}

	normalizeImage(im, r, c, norm_im);
	histogram(level, norm_im, label, r, c);

	//Display histogram with 100 levels...
	/*std::cout << "\nHistogram:\n";
	for (int i = 0; i < L; i++) {
		std::cout << i / double(L) << " - " << (i + 1) / double(L) << " : " << level[i] << endl;
	}*/

	int N = r*c;
	for (int i = 0; i < L; i++)
		p[i] = level[i] / double(N);

	for (int k = 0; k < L; k++) {
		/*if (level[k] == 0) { //no occurance...
		std::cout << "\nNone\n";
		continue;
		}*/
		var_b = betweenClassVariance(p, k);
		if (var_b > argmax_var_b) {
			argmax_var_b = var_b;
			thresh = (k + 1) / double(L);
		}
	}

	T = (2 * thresh - 1 / double(L)) / 2.0;
	//std::cout << "\nThe threshold level : " << T << endl;

	binary(norm_im, T, bin_im, r, c);

	for (int i = 0; i < r; i++) {
		delete[] norm_im[i];
		delete[] label[i];
	}
	delete[] level;
	delete[] p;

	return(T);
}

