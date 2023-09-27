#include <stdio.h>
#include<conio.h>
#include<iostream>
#include "opencv2/imgproc/imgproc.hpp"
//#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "ConnectedComponent.h"

using namespace std;
using namespace cv;

void diagonalBox(int **roi, int h, int w, int& c1, int& c2, int& c3, int& c4) {
	int flag = 0, flag2 = 0, i_min, i_max, j_min, j_max, prev_flag = 0, prev_prev_flag, c_min, c_trans;
	int **inv_roi = new int*[h];
	for (int i = 0; i < h; i++) {
		inv_roi[i] = new int[w];
		for (int j = 0; j < w; j++) {
			inv_roi[i][j] = (255 - roi[i][j]) / 255;
		}
	}

	c_min = (h < w) ? h - 1: w - 1;
	for (int c = 0; c < h + w - 1; c++) {
		prev_prev_flag = prev_flag;
		prev_flag = flag;
		flag = 0;
		i_min = 0;
		i_max = c;
		if (c > c_min && w == h) {
			i_min = c - w + 1;
			i_max = h - 1;
		}
		if (c > c_min && h > w) {
			i_min = c - w + 1;
			i_max = c;
			if (i_max >= h)
				i_max = h - 1;
		}
		if (c > c_min && w > h) {
			i_min = 0;
			i_max = h - 1;
			if (c >= w) 
				i_min = c - w + 1;
		}
		for (int i = i_max; i >= i_min; i--) {
			if (c - i >= w)
				break;
			if (inv_roi[i][c - i] == 1) {
				flag = 1;
				break;
			}
		}
		if (prev_prev_flag == 0 && prev_flag == 0 && flag == 1) {
			c1 = c - 1;
			//std::cout << "\nTop left diagonal at c = " << c - 1;
		}
		if (prev_prev_flag == 1 && prev_flag == 0 && flag == 0) {
			c2 = c - 1;
			//std::cout << "\nBottom right diagonal at c = " << c - 1;
			break;
		}
	}
	if (prev_flag == 1 && flag == 0) {
		c2 = h + w - 2;
	}
	if (flag == 1) {
		c2 = h + w - 1;
		//std::cout << "\n######################Bottom right diagonal at c = " << c2;
	}


	flag = 0; prev_flag = 0; prev_prev_flag = 0;
	c_trans = w - h;
	if (c_trans >= 0) {
		//std::cout << "\n******** w>h ";
		for (int c = w - 1; c >= c_trans; c--) {
			prev_prev_flag = prev_flag;
			prev_flag = flag;
			flag = 0;
			j_max = w - 1;
			j_min = c;
			for (int j = j_min; j <= j_max; j++) {
				if (inv_roi[j - c][j] == 1) {
					flag = 1;
					break;
				}
			}
			if (prev_prev_flag == 0 && prev_flag == 0 && flag == 1) {
				c3 = c + 1;
				//std::cout << "\nTop right diagonal at c = " << c3;
			}
			if (prev_prev_flag == 1 && prev_flag == 0 && flag == 0) {
				c4 = c - 1;
				//std::cout << "\nBottom left diagonal at c = " << c4;
				flag2 = 1;
				//std::cout << "\nFirst part";
			}
		}
		if (flag2 == 0) {
			for (int c = c_trans - 1; c >= 0; c--) {
				prev_prev_flag = prev_flag;
				prev_flag = flag;
				flag = 0;
				i_max = h - 1;
				i_min = 0;
				for (int i = i_min; i <= i_max; i++) {
					if (inv_roi[i][c + i] == 1) {
						flag = 1;
						break;
					}
				}
				if (prev_prev_flag == 0 && prev_flag == 0 && flag == 1) {
					c3 = c + 1;
					//std::cout << "\nTop right diagonal at c = " << c3;
				}
				if (prev_prev_flag == 1 && prev_flag == 0 && flag == 0) {
					c4 = c + 1;
					//std::cout << "\nBottom left diagonal at c = " << c4;
					flag2 = 1;
					//std::cout << "\n2nd part";
				}
			}
		}		
		if (flag2 == 0) {
			for (int c = -1; c >= 1 - h; c--) {
				prev_prev_flag = prev_flag;
				prev_flag = flag;
				flag = 0;
				i_max = h - 1;
				i_min = -c;
				for (int i = i_min; i <= i_max; i++) {
					if (inv_roi[i][c + i] == 1) {
						flag = 1;
						break;
					}
				}
				if (prev_prev_flag == 0 && prev_flag == 0 && flag == 1) {
					c3 = c + 1;
					//std::cout << "\nTop right diagonal at c = " << c3;
				}
				if (prev_prev_flag == 1 && prev_flag == 0 && flag == 0) {
					c4 = c + 1;
					//std::cout << "\nBottom left diagonal at c = " << c4;
					flag2 = 1;
					//std::cout << "\n3rd part";
				}
			}
		}
		if (flag2 == 0) {
			c4 = -h;
			//std::cout << "\nBottom left diagonal at c = " << c4;
		}
	}

	if (c_trans < 0) {
		//std::cout << "\n******** h>w ";
		for (int c = w - 1; c >= 0; c--) {
			prev_prev_flag = prev_flag;
			prev_flag = flag;
			flag = 0;
			j_max = w - 1;
			j_min = c;
			for (int j = j_min; j <= j_max; j++) {
				if (inv_roi[j - c][j] == 1) {
					flag = 1;
					break;
				}
			}
			if (prev_prev_flag == 0 && prev_flag == 0 && flag == 1) {
				c3 = c + 1;
				//std::cout << "\nTop right diagonal at c = " << c3;
			}
			if (prev_prev_flag == 1 && prev_flag == 0 && flag == 0) {
				c4 = c + 1;
				//std::cout << "\nBottom left diagonal at c = " << c4;
				flag2 = 1;
			}
		}
		if (flag2 == 0) {
			for (int c = -1; c >= c_trans; c--) {
				prev_prev_flag = prev_flag;
				prev_flag = flag;
				flag = 0;
				j_max = w - 1;
				j_min = 0;
				for (int j = j_min; j <= j_max; j++) {
					if (inv_roi[j - c][j] == 1) {
						flag = 1;
						break;
					}
				}
				if (prev_prev_flag == 0 && prev_flag == 0 && flag == 1) {
					c3 = c + 1;
					//std::cout << "\nTop right diagonal at c = " << c3;
				}
				if (prev_prev_flag == 1 && prev_flag == 0 && flag == 0) {
					c4 = c + 1;
					//std::cout << "\nBottom left diagonal at c = " << c4;
					flag2 = 1;
				}
			}			
		}
		if (flag2 == 0) {
			for (int c = c_trans - 1; c >= 1 - h; c--) {
				prev_prev_flag = prev_flag;
				prev_flag = flag;
				flag = 0;
				i_max = h - 1;
				i_min = -c;
				for (int i = i_min; i <= i_max; i++) {
					if (inv_roi[i][c + i] == 1) {
						flag = 1;
						break;
					}
				}
				if (prev_prev_flag == 0 && prev_flag == 0 && flag == 1) {
					c3 = c + 1;
					//std::cout << "\nTop right diagonal at c = " << c3;
				}
				if (prev_prev_flag == 1 && prev_flag == 0 && flag == 0) {
					c4 = c + 1;
					//std::cout << "\nBottom left diagonal at c = " << c4;
					flag2 = 1;
				}
			}
		}
		if (flag2 == 0) {
			c4 = -h;
			//std::cout << "\nBottom left diagonal at c = " << c4;
		}
	}
}