#include "KWK.h"
#include<iostream>
#include<math.h>

using namespace std;
//int row, col;

//No.of 01 transitions in the nbd...
void KWKtransitions(int **src, int **ZO, int row, int col) {
	int *p = new int[10];
	for (int i = 1; i <= row; i++) {
		for (int j = 1; j <= col; j++) {
			p[0] = src[i][j];
			p[1] = src[i - 1][j - 1];
			p[2] = src[i - 1][j];
			p[3] = src[i - 1][j + 1];
			p[4] = src[i][j + 1];
			p[5] = src[i + 1][j + 1];
			p[6] = src[i + 1][j];
			p[7] = src[i - 1][j - 1];
			p[8] = src[i][j - 1];
			p[9] = p[1];
			for (int k = 2; k < 10; k++) {
				if (p[k - 1] == 0 && p[k] == 1)  //0 to 1 transition...
					ZO[i - 1][j - 1]++;
			}
		}
	}
}
/*____________________________________________________________________________________________*/

//No. of non-zero nbrs...
void KWKnonZeroNbrs(int **src, int **NZ, int row, int col) {
	for (int i = 1; i <= row; i++) {
		for (int j = 1; j <= col; j++) {
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

//Fast Parallel Thinning algorithm...As proposed by Kwan, Wong, Kang (2001)
int KWKthinning(int **src, int row, int col, int **dest, int iter) {
	int **ZO, **NZ, **src_copy, **f1a, **f1b, **f2, *p;
	int filter = 3, noChange = 1, flag = 0;

	p = new int[9];
	ZO = new int*[row];
	NZ = new int*[row];
	f1a = new int*[row];
	f1b = new int*[row];
	f2 = new int*[row];
	src_copy = new int*[row + 2];
	for (int i = 0; i < row + 2; i++) {
		src_copy[i] = new int[col + 2];
		for (int j = 0; j < col + 2; j++)
			src_copy[i][j] = 0;
	}
	for (int i = 0; i < row; i++) {
		ZO[i] = new int[col];
		NZ[i] = new int[col];
		f1a[i] = new int[col];
		f1b[i] = new int[col];
		f2[i] = new int[col];
		for (int j = 0; j < col; j++) {
			ZO[i][j] = 0;
			NZ[i][j] = 0;
			f1a[i][j] = 0;
			f1b[i][j] = 0;
			f2[i][j] = 0;
			src_copy[i + 1][j + 1] = (255 - src[i][j]) / 255;
			dest[i][j] = (255 - src[i][j]) / 255;
		}
	}
	KWKtransitions(src_copy, ZO, row, col);
	KWKnonZeroNbrs(src_copy, NZ, row, col);

	//1st Pass : 1st sub-iteration....
	for (int i = 1; i <= row; i++) {
		for (int j = 1; j <= col; j++) {
			if (src_copy[i][j] == 0)
				continue;
			p[0] = src_copy[i][j];
			p[1] = src_copy[i - 1][j - 1];
			p[2] = src_copy[i - 1][j];
			p[3] = src_copy[i - 1][j + 1];
			p[4] = src_copy[i][j + 1];
			p[5] = src_copy[i + 1][j + 1];
			p[6] = src_copy[i + 1][j];
			p[7] = src_copy[i - 1][j - 1];
			p[8] = src_copy[i][j - 1];
			if (2 <= NZ[i - 1][j - 1] && NZ[i - 1][j - 1] <= 6) {
				if (ZO[i - 1][j - 1] == 1) {
					if (p[2] * p[4] * p[6] == 0) {
						if (p[4] * p[6] * p[8] == 0) {
							f1a[i - 1][j - 1] = 1;
						}
					}
				}
			}
		}
	}

	//save in dest[][]...
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			if (f1a[i][j] == 1) {
				src_copy[i + 1][j + 1] = 0;
			}
		}
	}

	//1st Pass : 2nd sub-iteration....
	for (int i = 1; i <= row; i++) {
		for (int j = 1; j <= col; j++) {
			if (src_copy[i][j] == 0)
				continue;
			p[0] = src_copy[i][j];
			p[1] = src_copy[i - 1][j - 1];
			p[2] = src_copy[i - 1][j];
			p[3] = src_copy[i - 1][j + 1];
			p[4] = src_copy[i][j + 1];
			p[5] = src_copy[i + 1][j + 1];
			p[6] = src_copy[i + 1][j];
			p[7] = src_copy[i - 1][j - 1];
			p[8] = src_copy[i][j - 1];
			if (3 <= NZ[i - 1][j - 1] && NZ[i - 1][j - 1] <= 6) {
				if (ZO[i - 1][j - 1] == 1) {
					if (p[2] * p[4] * p[8] == 0) {
						if (p[2] * p[6] * p[8] == 0) {
							f1b[i - 1][j - 1] = 1;
						}
					}
				}
			}
		}
	}

	//save in dest[][]...
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			if (f1b[i][j] == 1) {
				src_copy[i + 1][j + 1] = 0;
			}
		}
	}

	//2nd Pass....
	for (int i = 1; i <= row; i++) {
		for (int j = 1; j <= col; j++) {
			if (src_copy[i][j] == 0)
				continue;
			p[0] = src_copy[i][j];
			p[1] = src_copy[i - 1][j - 1];
			p[2] = src_copy[i - 1][j];
			p[3] = src_copy[i - 1][j + 1];
			p[4] = src_copy[i][j + 1];
			p[5] = src_copy[i + 1][j + 1];
			p[6] = src_copy[i + 1][j];
			p[7] = src_copy[i - 1][j - 1];
			p[8] = src_copy[i][j - 1];
			if (p[1] * p[8] * p[6] == 1 && p[3] == 0) {
				f2[i - 1][j - 1] = 1;
				continue;
			}
			if (p[3] * p[4] * p[6] == 1 && p[1] == 0) {
				f2[i - 1][j - 1] = 1;
				continue;
			}
			if (p[5] * p[6] * p[8] == 1 && p[3] == 0) {
				f2[i - 1][j - 1] = 1;
				continue;
			}
			if (p[4] * p[6] * p[7] == 1 && p[1] == 0) {
				f2[i - 1][j - 1] = 1;
				continue;
			}
		}
	}

	//save in dest[][]...
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			if (f2[i][j] == 1) {
				src_copy[i + 1][j + 1] = 0;
			}
		}
	}

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			dest[i][j] = 255 - src_copy[i + 1][j + 1] * 255;
		}
	}

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			if (src[i][j] != dest[i][j]) {
				noChange = 0;
				break;
			}
		}
		if (!noChange)
			break;
	}

	for (int i = 0; i < row; i++) {
		delete[] NZ[i];
		delete[] ZO[i];
		delete[] f1a[i];
		delete[] f1b[i];
		delete[] f2[i];
	}
	delete[] NZ;
	delete[] ZO;
	delete[] f1a;
	delete[] f1b;
	delete[] f2;
	delete[] p;

	for (int i = 0; i < row + 2; i++)
		delete[] src_copy[i];
	delete[] src_copy;

	return(noChange);
}
/*____________________________________________________________________________________________*/

//Function call...
void KWKthinCC(int **roi, int row, int col, int **thin_roi) {
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
		noChange = KWKthinning(src_im, row, col, dest_im, k);
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