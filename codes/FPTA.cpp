#include "FPTA.h"
#include<iostream>
#include<math.h>

using namespace std;

//No.of 01 transitions in the nbd...
void transitions(int **src, int **ZO, int r, int c) {
	int *p = new int[10];
	for (int i = 1; i <= r; i++) {
		for (int j = 1; j <= c; j++) {
			p[0] = src[i][j];
			p[1] = src[i - 1][j];
			p[2] = src[i - 1][j + 1];
			p[3] = src[i][j + 1];
			p[4] = src[i + 1][j + 1];
			p[5] = src[i + 1][j];
			p[6] = src[i + 1][j - 1];
			p[7] = src[i][j - 1];
			p[8] = src[i - 1][j - 1];
			p[9] = src[i - 1][j];
			for (int k = 2; k < 10; k++) {
				if (p[k-1]==0 && p[k]==1)  //0 to 1 transition...
					ZO[i-1][j-1]++;
			}
		}
	}
}
/*____________________________________________________________________________________________*/

//No. of non-zero nbrs...
void nonZeroNbrs(int **src, int **NZ, int r, int c) {
	for (int i = 1; i <= r ; i++) {
		for (int j = 1; j <= c; j++) {
			for (int l = i - 1; l <= i + 1; l++) {
				for (int m = j - 1; m <= j + 1; m++) {
					if (l == i && m == j)
						continue;
					if (src[l][m])
						NZ[i - 1][j - 1]++;
				}
			}
		}
	}
}
/*____________________________________________________________________________________________*/

//Fast Parallel Thinning algorithm...As proposed by T.Y.Zhang and C.Y.Suen (1984)
int thinning(int **src, int r, int c, int **dest, int iter) {
	int **ZO, **NZ, **src_copy, **f1, **f2, *p;
	int filter = 3, noChange = 1, flag = 0;

	p = new int[9];
	ZO = new int*[r];
	NZ = new int*[r];
	f1 = new int*[r];
	f2 = new int*[r];
	src_copy = new int*[r + 2];
	for (int i = 0; i < r + 2; i++) {
		src_copy[i] = new int[c + 2];
		for (int j = 0; j < c + 2; j++)
			src_copy[i][j] = 0;
	}
	for (int i = 0; i < r; i++) {
		ZO[i] = new int[c];
		NZ[i] = new int[c];
		f1[i] = new int[c];
		f2[i] = new int[c];
		for (int j = 0; j < c; j++) {
			ZO[i][j] = 0;
			NZ[i][j] = 0;
			f1[i][j] = 0;
			f2[i][j] = 0;
			src_copy[i + 1][j + 1] = (255 - src[i][j]) / 255;
			dest[i][j] = src[i][j];
		}
	}
	transitions(src_copy, ZO, r, c);
	nonZeroNbrs(src_copy, NZ, r, c);
	
	for (int i = 1; i <= r; i++) {
		for (int j = 1; j <= c; j++) {
			if (src_copy[i][j] == 0)
				continue;
			p[0] = src_copy[i][j];
			p[1] = src_copy[i - 1][j];
			p[2] = src_copy[i - 1][j + 1];
			p[3] = src_copy[i][j + 1];
			p[4] = src_copy[i + 1][j + 1];
			p[5] = src_copy[i + 1][j];
			p[6] = src_copy[i + 1][j - 1];
			p[7] = src_copy[i][j - 1];
			p[8] = src_copy[i - 1][j - 1];
			if (2 <= NZ[i - 1][j - 1] && NZ[i - 1][j - 1] <= 6) {
				if (ZO[i - 1][j - 1] == 1) {
					if(iter%2 == 1){ //odd iteration....
						if (p[1] * p[3] * p[5] == 0) {
							if (p[3] * p[5] * p[7] == 0) {
								f1[i-1][j-1] = 1;
							}
						}
					}
					if(iter%2 == 0){   //even iteration....
						if (p[1] * p[3] * p[7] == 0) {
							if (p[1] * p[5] * p[7] == 0) {
								f1[i-1][j-1] = 1;
							}
						}
					}
				}
			}
		}
	}
	
	//save in dest[][]...
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			if (f1[i][j] == 1) {
				dest[i][j] = 255;
			}
			else
				dest[i][j] = 255 - 255 * src_copy[i + 1][j + 1];
		}
	}

	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			if (src[i][j] != dest[i][j]) {
				noChange = 0;
				break;
			}
		}
		if (!noChange)
			break;
	}

	for (int i = 0; i < r; i++) {
		delete[] NZ[i];
		delete[] ZO[i];
		delete[] f1[i];
		delete[] f2[i];
	}
	delete[] NZ;
	delete[] ZO;
	delete[] f1;
	delete[] f2;
	delete[] p;

	for (int i = 0; i < r+2; i++)
		delete[] src_copy[i];
	delete[] src_copy;

	return(noChange);
}
/*____________________________________________________________________________________________*/

//Function call...
void thinCC(int **roi, int row, int col, int **thin_roi) {
	int **src_im, **dest_im;
	int noChange, k = 1;
	src_im = new int *[row];
	dest_im = new int *[row];
	for (int i = 0; i < row; i++) {
		src_im[i] = new int[col];
		dest_im[i] = new int[col];
		for (int j = 0; j < col; j++) {
			dest_im[i][j] = roi[i][j];
		}
	}
	do {
		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
				src_im[i][j] = dest_im[i][j];
		noChange = thinning(src_im, row, col, dest_im, k);
		k++;
		//std::cout << "k = " << k << "\n";
	} while (!noChange);
	std::cout << "\nNo. of steps:" << k << "\n";
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			thin_roi[i][j] = dest_im[i][j];
		}
	}
}