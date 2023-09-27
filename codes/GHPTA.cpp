#include "GHPTA.h"
#include<iostream>
#include<math.h>

using namespace std;

//No. of distinct components in the 8-component neighborhood....
void C(int **src, int **CP, int r, int c){
	int *p = new int[9];
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

			if (p[1]==0 && (p[2]==255 || p[3]==255))
				CP[i-1][j-1]++;
			if (p[3]==0 && (p[4]==255 || p[5]==255))
				CP[i-1][j-1]++;
			if (p[5]==0 && (p[6]==255 || p[7]==255))
				CP[i-1][j-1]++;
			if (p[7]==0 && (p[8]==255 || p[1]==255))
				CP[i-1][j-1]++;
		}
	}
}
/*____________________________________________________________________________________________*/

//N(P) = min(N1(P),N2(P))....
void N_P(int **src, int **NP, int r, int c) {
	int N1, N2;

	int *p = new int[9];
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

			N1 = 0;
			if (p[8] || p[1])
				N1++;
			if (p[2] || p[3])
				N1++;
			if (p[4] || p[5])
				N1++;
			if (p[6] || p[7])
				N1++;

			N2 = 0;
			if (p[1] || p[2])
				N2++;
			if (p[3] || p[4])
				N2++;
			if (p[5] || p[6])
				N2++;
			if (p[7] || p[8])
				N2++;

			if (N1 < N2)
				NP[i-1][j-1] = N1;
			else
				NP[i-1][j-1] = N2;
		}
	}
}
/*____________________________________________________________________________________________*/

//Guo-Hall Parallel Thinning algorithm...(1989)
int GHthinning(int **src, int r, int c, int **dest, int iter) {
	int **NP, **CP, **src_copy, **f1, **f2, *p;
	int filter = 3, noChange = 1, flag = 0;

	p = new int[9];
	NP = new int*[r];
	CP = new int*[r];
	f1 = new int*[r];
	f2 = new int*[r];
	src_copy = new int*[r + 2];
	for (int i = 0; i < r + 2; i++) {
		src_copy[i] = new int[c + 2];
		for (int j = 0; j < c + 2; j++)
			src_copy[i][j] = 0;
	}
	for (int i = 0; i < r; i++) {
		NP[i] = new int[c];
		CP[i] = new int[c];
		f1[i] = new int[c];
		f2[i] = new int[c];
		for (int j = 0; j < c; j++) {
			NP[i][j] = 0;
			CP[i][j] = 0;
			f1[i][j] = 0;
			f2[i][j] = 0;
			src_copy[i + 1][j + 1] = (255 - src[i][j]) / 255;
			dest[i][j] = src[i][j];
		}
	}
	C(src_copy, CP, r, c);
	N_P(src_copy, NP, r, c);
	
	//First sub-iteration....
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

			if (CP[i - 1][j - 1] == 1) {
				if (NP[i - 1][j - 1] == 2 || NP[i - 1][j - 1] == 3) {
					if (iter % 2 == 1) {    //odd iteration....
						std::cout << "odd\n";
						if ((p[1] == 0 || p[2] == 0 || !p[4] == 0) || p[3] == 0) {
							f1[i - 1][j - 1] = 1;  //delete pixel....							
						}
					}	
					if (iter % 2 == 0) {     //even iteration....
						std::cout << "even\n";
						if ((p[5] == 0 || p[6] == 0 || !p[8] == 0) && p[7] == 0) {
							f1[i - 1][j - 1] = 1;    //delete pixel....							
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

	//std::cout << "\nnoChange = " << noChange;

	for (int i = 0; i < r; i++) {
		delete[] CP[i];
		delete[] NP[i];
		delete[] f1[i];
		delete[] f2[i];
	}
	delete[] CP;
	delete[] NP;
	delete[] f1;
	delete[] f2;
	delete[] p;

	for (int i = 0; i < r + 2; i++)
		delete[] src_copy[i];
	delete[] src_copy;

	return(noChange);
}
/*____________________________________________________________________________________________*/

//Function call...
void GHthinCC(int **roi, int row, int col, int **thin_roi) {
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
		noChange = GHthinning(src_im, row, col, dest_im, k);
		k++;
		std::cout << "k = " << k << "\n";
	} while (!noChange);
	std::cout << "\nNo. of steps:" << k << "\n";
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			thin_roi[i][j] = dest_im[i][j];
		}
	}
}