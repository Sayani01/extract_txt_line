//#include <opencv2\core\core.hpp>
//#include <opencv2\imgproc\imgproc.hpp>
//#include <opencv2\highgui\highgui.hpp>
//
//
//#include <stdio.h>
//#include <conio.h>
//#include <iostream>
//#include <vector>
//#include <string>
//
//#include "NICK.h"
//#include "Otsu.h"
//#include "RegLine.h"
//#include "ConnectedComponent.h"
//#include "BoundingBox.h"
//#include "GaussBlurRotated.h"
//#include "NeighbourCC.h"
//#include "Projections.h"
//#include "CentroidLine.h"
//
//using namespace std;
//using namespace cv;
//namespace cv{
//	using std::vector;
//};
//
//extern int labelNum = 0, r = 0, c = 0;
//extern double AH = 0;
//extern int line_no = 0;
//
//int minimum(int arr[3][3]) {
//	int m = arr[0][0];
//	for (int i = 0; i < 3; i++) {
//		for (int j = 0; j < 3; j++) {
//			if (arr[i][j] < m)
//				m = arr[i][j];
//		}
//	}
//	return(m);
//}
//
//void blur(int **img) {
//	int max, p[3][3];
//	int **img_copy;
//	img_copy = new int*[r];
//	for (int i = 0; i < r; i++) {
//		img_copy[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			img_copy[i][j] = 0;
//		}
//	}
//	for (int i = 1; i < r - 1; i++) {
//		for (int j = 1; j < c - 1; j++) {
//			if (img[i][j] == 0) continue;
//			for (int k = -1; k <= 1; k++) {
//				for (int l = -1; l <= 1; l++)
//					p[k + 1][l + 1] = img[i + k][j + l];
//			}
//			img_copy[i][j] = minimum(p);
//		}
//	}
//	for (int i = 1; i < r - 1; i++) {
//		for (int j = 1; j < c - 1; j++) {
//			if (img[i][j] == 0) continue;
//			img[i][j] = img_copy[i][j];
//		}
//	}
//}
//
////Matrix <--> array conversions...
//void mat2arr(Mat A, int **arr) {
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			arr[i][j] = A.at<uchar>(i, j);
//		}
//	}
//}
//void arr2mat(int **arr, Mat A) {
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			A.at<uchar>(i, j) = arr[i][j];
//		}
//	}
//}
////Overloads...
//void mat2arr(Mat A, int **arr, int row, int col) {
//	for (int i = 0; i < row; i++) {
//		for (int j = 0; j < col; j++) {
//			arr[i][j] = A.at<uchar>(i, j);
//		}
//	}
//}
//extern void arr2mat(int **arr, Mat A, int row, int col) {
//	for (int i = 0; i < row; i++) {
//		for (int j = 0; j < col; j++) {
//			A.at<uchar>(i, j) = arr[i][j];
//		}
//	}
//}
//void invImage(int **arr, int **inv_arr, int h = r, int w = c) {
//	for (int i = 0; i < h; i++) {
//		for (int j = 0; j < w; j++) {
//			inv_arr[i][j] = 255 - arr[i][j];
//		}
//	}
//}
////find histogram to calculate average height...
//double averageHeight(int *heights, int sz, int range = 5) {
//	int max = heights[0], size, pos;
//	int *hist;
//	double AH;
//
//	//get the maximum height...
//	for (int i = 1; i < sz; i++) {
//		if (heights[i] > max)
//			max = heights[i];
//	}
//	size = int(max / range) + 1;
//	hist = new int[size]();
//	for (int i = 0; i < sz; i++) {
//		int val = int(heights[i] / range);
//		hist[val]++;
//	}
//
//	std::cout << "The histogram:\n";
//	for (int i = 0; i < size; i++)
//		std::cout << range*i << " - " << range*(i + 1) << " -> " << hist[i] << "components\n";
//
//	max = hist[3];
//	pos = 3;
//	for (int i = 4; i < size; i++) {
//		if (hist[i] > max) {
//			max = hist[i];
//			pos = i;
//		}
//	}
//
//	AH = (range*(pos + 1));
//	std::cout << "average height = " << AH << "\n";
//
//	delete[] hist;
//	return(AH);
//}
//double averageDiagHeight(double *diagheights, int sz, int range = 5) {
//	int size, pos;
//	int *hist;
//	double max = diagheights[0], ah;
//
//	//get the maximum height...
//	for (int i = 1; i < sz; i++) {
//		if (diagheights[i] > max)
//			max = diagheights[i];
//	}
//	std::cout << "\nMax height:" << max;
//	size = int(max / range) + 1;
//	hist = new int[size];
//	for (int i = 0; i < size; i++) {
//		hist[i] = 0;
//	}
//	for (int i = 0; i < sz; i++) {
//		int val = int(diagheights[i] / range);
//		hist[val]++;
//	}
//
//	std::cout << "The histogram:\n";
//	for (int i = 0; i < size; i++)
//		std::cout << range*i << " - " << range*(i + 1) << " -> " << hist[i] << "components\n";
//
//	max = hist[3];
//	pos = 3;
//	for (int i = 4; i < size; i++) {
//		if (hist[i] > max) {
//			max = hist[i];
//			pos = i;
//		}
//	}
//
//	ah = double(range*(pos + 1));
//	std::cout << "average height = " << ah << "\n";
//
//	delete[] hist;
//	return(ah);
//}
//
//
///*--------------------------------------------------------------------------------------------*/
///*                  Trim CC's to exclude extensions or accent marks in words                  */
///*--------------------------------------------------------------------------------------------*/
//void trimHV(int **img, int **box, int **label, double **boundingBox, int *trimlabel) {
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1 || boundingBox[i][0] == 2 || trimlabel[i] != 1) continue;
//
//		int i_start = box[i][4], j_start = box[i][2];
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//
//		//continue trimming for larger CCs...
//		if (h < w && h <= AH) continue;
//		else if (w < h && w <= AH) continue;
//
//		int **roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				if (label[i_start + k][j_start + l] == i) {
//					roi[k][l] = img[i_start + k][j_start + l];
//				}
//				else
//					roi[k][l] = 255;
//			}
//		}
//
//		int kstart = 0, kend;
//		while (kstart < h) {
//			int sum = 0;
//			kend = kstart + 2;
//			if (kend >= h) kend = h - 1;
//			for (int stripi = kstart; stripi <= kend; stripi++)
//				for (int stripj = 0; stripj < w; stripj++)
//					if (roi[stripi][stripj] == 0)
//						sum++;
//
//			if (double(sum) / (3*w) < 0.18) {
//				for (int stripi = kstart; stripi <= kend; stripi++)
//					for (int stripj = 0; stripj < w; stripj++)
//						if (label[i_start + stripi][j_start + stripj] == i)
//							img[i_start + stripi][j_start + stripj] = 255;
//			}
//			kstart = kend + 1;			
//		}
//
//
//		//delete...
//		for (int k = 0; k < h; k++)
//			delete roi[k];
//		delete roi;
//
//	}
//}
//
//void trimV(int **img, int **box, int **label, double **boundingBox, int *trimlabel) {
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1 || boundingBox[i][0] == 2 || trimlabel[i] != 2) continue;
//
//		int i_start = box[i][4], j_start = box[i][2];
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//
//		//continue trimming for larger CCs...
//		if (h < w && h <= AH) continue;
//		else if (w < h && w <= AH) continue;
//
//		int **roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				if (label[i_start + k][j_start + l] == i) {
//					roi[k][l] = img[i_start + k][j_start + l];
//				}
//				else
//					roi[k][l] = 255;
//			}
//		}
//
//		int kstart = 0, kend;
//		while (kstart < w) {
//			int sum = 0;
//			kend = kstart + 2;
//			if (kend >= w) kend = w - 1;
//			for (int stripj = kstart; stripj <= kend; stripj++)
//				for (int stripi = 0; stripi < h; stripi++)
//					if (roi[stripi][stripj] == 0)
//						sum++;
//
//			if (double(sum) / (3 * h) < 0.18) {
//				for (int stripj = kstart; stripj <= kend; stripj++)
//					for (int stripi = 0; stripi < h; stripi++)
//						if (label[i_start + stripi][j_start + stripj] == i)
//							img[i_start + stripi][j_start + stripj] = 255;
//			}
//			kstart = kend + 1;
//		}
//
//
//		//delete...
//		for (int k = 0; k < h; k++)
//			delete roi[k];
//		delete roi;
//
//	}
//}
//
//int deletePixels(int *hist, int sz, int range) {
//	int istart = 0, iend, count = 0, countwhite = 0;
//	while (istart < sz) {
//		double sum = 0, sumwhite = 0;
//		iend = istart + 3;
//		if (iend >= sz) iend = sz - 1;
//		for (int i = istart; i <= iend; i++) {
//			sum += hist[i];
//			sumwhite += (range - hist[i]);
//		}
//		if (sum / (3 * range) < 0.18) {
//			count += sum;
//			countwhite += sumwhite;
//		}
//		istart = iend + 1;
//	}
//	return(countwhite);
//}
//
//
//
///*--------------------------------------------------------------------------------------------*/
///*       Find the longest line through centroid in CC by rotating the regression line         */
///*--------------------------------------------------------------------------------------------*/
//int crossing(int **roi, int h, int w, Point p1, Point p2, double slope, int reg) {
//	int x, y, prev = 0, curr = 0, count = 0;
//	if (reg == 1) {   //Y on X line
//		prev = 255;
//		for (int x = p1.x + 1; x <= p2.x; x++) {
//			y = (int)(p1.y + slope*(x - p1.x));
//			if (y < 0) y = 0;
//			if (y >= h) y = h - 1;
//			if (prev == 255 && roi[y][x] == 0)
//				count++;
//			prev = roi[y][x];
//		}
//	}
//	if (reg == -1) {   //X on Y line
//		prev = 255;
//		for (int y = p1.y + 1; y <= p2.y; y++) {
//			x = (int)(p1.x + (y - p1.y)/slope);
//			if (x < 0) x = 0;
//			if (x >= w) x = w - 1;
//			if (prev == 255 && roi[y][x] == 0)
//				count++;
//			prev = roi[y][x];
//		}
//	}
//	std::cout << "\ncrossing number = " << count;
//	return(count);
//}
//void rotateRegLine(int **roi, int h, int w, int center[2], double slope, int flag) {
//	double m_clk = 0, m_anticlk = 0, rot_slope, min_dist = 0, interceptyx;
//	for (int k = 0; k < 45; k += 5) {
//		Point p1, p2, p3, p4;
//		rot_slope = tan(k * CV_PI / 180);
//		m_clk = (slope - rot_slope) / (1 + slope*rot_slope);
//		interceptyx = (center[0]) - m_clk*(center[1]);
//		p1.x = 0;
//		p1.y = interceptyx;
//		if (p1.y >= h) {
//			p1.y = h - 1;
//			p1.x = (p1.y + interceptyx) / m_clk;
//		}
//		p2.x = w - 1;
//		p2.y = m_clk*(w - 1) + interceptyx;
//		if (p2.y < 0) {
//			p2.y = 0;
//			p2.x = interceptyx / m_clk;
//		}
//		m_anticlk = (slope + rot_slope) / (1 - slope*rot_slope);
//	}
//}
//
///*----------------------------------------------------------------*/
///*       Find the Regression lines through centroid in CC         */
///*----------------------------------------------------------------*/
//void centroidCC(int **roi, int h, int w, int center[2]) {
//	if (roi[center[0]][center[1]] == 0)
//		return;
//
//	int flag = 0, k = 1, x_min = 0, y_min = 0;
//	while (flag == 0) {
//		double dist = 100;
//		int i_min = center[0] - k;
//		int i_max = center[0] + k;
//		int j_min = center[1] - k;
//		int j_max = center[1] + k;
//		if (i_min < 0)
//			i_min = 0;
//		if (i_max >= h)
//			i_max = h - 1;
//		if (j_min < 0)
//			j_min = 0;
//		if (j_max >= w)
//			j_max = w - 1;
//		for (int i = i_min + 1; i < i_max; i++) {
//			if (roi[i][j_min] == 0) {
//				flag = 1;
//				double d = sqrt((i - center[0])*(i - center[0]) + (j_min - center[1])*(j_min - center[1]));
//				if (dist > d) {
//					dist = d;
//					x_min = j_min;
//					y_min = i;
//				}
//			}
//			if (roi[i][j_max] == 0) {
//				flag = 1;
//				double d = sqrt((i - center[0])*(i - center[0]) + (j_max - center[1])*(j_max - center[1]));
//				if (dist > d) {
//					dist = d;
//					x_min = j_max;
//					y_min = i;
//				}
//			}
//		}
//		for (int j = j_min; j <= j_max; j++) {
//			if (roi[i_min][j] == 0) {
//				flag = 1;
//				double d = sqrt((i_min - center[0])*(i_min - center[0]) + (j - center[1])*(j - center[1]));
//				if (dist > d) {
//					dist = d;
//					x_min = j;
//					y_min = i_min;
//				}
//			}
//			if (roi[i_max][j] == 0) {
//				flag = 1;
//				double d = sqrt((i_max - center[0])*(i_max - center[0]) + (j - center[1])*(j - center[1]));
//				if (dist > d) {
//					dist = d;
//					x_min = j;
//					y_min = i_max;
//				}
//			}
//		}
//		k++;
//	}
//	center[0] = y_min;
//	center[1] = x_min;
//}
//
//void skewLine(int **img, Mat cc_img, int **label, double **boundingBox, int **box) {
//	double **regSlope, skew;
//	int cross;
//	regSlope = new double*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		regSlope[i] = new double[2];
//		for (int j = 0; j < 2; j++)
//			regSlope[i][j] = 0;
//	}
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) //Not a significant component...
//			continue;
//		int center[2];
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//
//		//--------------------------------Select ROI and determine centroid for black pixels...
//		
//		int **roi;
//		int *X, *Y;
//		int count = 0, N = 0, i_start = box[i][4], j_start = box[i][2];
//		double  slopeyx = 0, slopexy = 0, interceptyx = 0, interceptxy = 0, spreadYX, spreadXY;
//		roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				if (label[i_start + k][j_start + l] == i) {
//					roi[k][l] = img[i_start + k][j_start + l];
//					N++;
//				}
//				else
//					roi[k][l] = 255;
//			}
//		}
//
//		X = new int[N];
//		Y = new int[N];
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				if (roi[k][l] == 0) {
//					X[count] = l;
//					Y[count] = k;
//					count++;
//				}
//			}
//		}
//		center[0] = (int)Mean(Y, N);
//		center[1] = (int)Mean(X, N);
//
//		//Locate centroid in CC...
//		//centroidCC(roi, h, w, center);
//		//std::cout << "\nCentroid for label " << i << ":" << i_start + center[0] << " , " << j_start + center[1];
//
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				if (label[i_start + k][j_start + l] == i) {
//					Vec3b color;
//					color.val[0] = 150;
//					color.val[1] = 150;
//					color.val[2] = 150;
//					cc_img.at<Vec3b>(i_start + k, j_start + l) = color;
//				}
//			}
//		}
//
//		//------------------------------------Determine regression line slope...			
//
//		regressionLine(X, Y, N, slopeyx, slopexy, interceptyx, interceptxy);
//		interceptyx = (i_start + center[0]) - slopeyx*(j_start + center[1]);
//		interceptxy = (j_start + center[1]) - slopexy*(i_start + center[0]);
//		//std::cout << "\nLabel " << i << " - Slope = " << slope << " , intercept = " << intercept << " , N = " << N;
//
//		//-------------------------------Display the regression line with longer spread through centroid...
//
//		Point p1, p2, p3, p4;
//		double m_clk = 0, m_anticlk = 0, rot_slope, min_dist = 0;
//
//		p1.x = box[i][2];
//		p1.y = slopeyx*box[i][2] + interceptyx;
//		p2.x = box[i][3];
//		p2.y = slopeyx*box[i][3] + interceptyx;
//		spreadYX = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
//		//line(cc_img, p1, p2, Scalar(255, 0, 0), 1.5, 8);
//
//		p3.y = box[i][4];
//		p3.x = slopexy*box[i][4] + interceptxy;
//		p4.y = box[i][5];
//		p4.x = slopexy*box[i][5] + interceptxy;
//		spreadXY = sqrt((p3.x - p4.x)*(p3.x - p4.x) + (p3.y - p4.y)*(p3.y - p4.y));
//		//line(cc_img, p3, p4, Scalar(255, 150, 0), 1.5, 8);
//
//		//------------------------Rotate the longer regression line about the centroid...
//		//                            ...display the line which cuts through the text maximum times...
//
//		int temp;
//		if (spreadYX > spreadXY) {
//			Point temp1, temp2;
//			double temp_spread;
//			line(cc_img, p1, p2, Scalar(255, 255, 0), 1.5, 8);
//			regSlope[i][0] = 1;
//			regSlope[i][1] = slopeyx;
//			skew = slopeyx;
//			temp1.x = p1.x - box[i][2];
//			temp1.y = p1.y - box[i][4];
//			temp2.x = p2.x - box[i][2];
//			temp2.y = p2.y - box[i][4];
//			cross = crossing(roi, h, w, temp1, temp2, slopeyx, 1);
//			for (int k = 0; k < 45; k += 5) {
//				rot_slope = tan(k * CV_PI / 180);
//
//				//Determine the clockwise rotations...
//				m_clk = (slopeyx - rot_slope) / (1 + slopeyx*rot_slope);
//				interceptyx = (i_start + center[0]) - m_clk*(j_start + center[1]);
//				p1.x = box[i][2];
//				p1.y = m_clk*box[i][2] + interceptyx;
//				if (p1.y >= box[i][5]) {
//					p1.y = box[i][5] - 1;
//					p1.x = (p1.y - interceptyx) / m_clk;
//				}
//				/*else if (p1.y < box[i][4]) {
//					p1.y = box[i][4];
//					p1.y = (p1.x - interceptyx) / m_clk;
//				}*/
//
//				p2.x = box[i][3];
//				p2.y = m_clk*box[i][3] + interceptyx;
//				if (p2.y < box[i][4]) {
//					p2.y = box[i][4];
//					p2.x = (p2.y - interceptyx) / m_clk;
//				}
//				/*else if (p2.y >= box[i][5]) {
//					p2.y = box[i][5] - 1;
//					p2.x = (p2.y - interceptyx) / m_clk;
//				}*/
//				temp1.x = p1.x - box[i][2];
//				temp1.y = p1.y - box[i][4];
//				temp2.x = p2.x - box[i][2];
//				temp2.y = p2.y - box[i][4];
//				temp_spread = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
//				temp = crossing(roi, h, w, temp1, temp2, m_clk, 1);
//				if (temp > cross) {
//					cross = temp;
//					skew = m_clk;
//					spreadYX = temp_spread;
//				}
//				if (temp == cross && temp_spread > spreadYX) {
//					skew = m_clk;
//					spreadYX = temp_spread;
//				}
//				//line(cc_img, p1, p2, Scalar(255, 150, 0), 1.5, 8);
//
//				//Determine the anti-clockwise rotation....
//				m_anticlk = (slopeyx + rot_slope) / (1 - slopeyx*rot_slope);
//				interceptyx = (i_start + center[0]) - m_anticlk*(j_start + center[1]);
//				p1.x = box[i][2];
//				p1.y = m_anticlk*box[i][2] + interceptyx; 
//				if (p1.y < box[i][4]) {
//					p1.y = box[i][4];
//					p1.x = (p1.y - interceptyx) / m_anticlk;
//				}
//				else if (p1.y >= box[i][5]) {
//					p1.y = box[i][5] - 1;
//					p1.x = (p1.y - interceptyx) / m_anticlk;
//				}
//
//				p2.x = box[i][3];
//				p2.y = m_anticlk*box[i][3] + interceptyx;
//				if (p2.y >= box[i][5]) {
//					p2.y = box[i][5] - 1;
//					p2.x = (p2.y - interceptyx) / m_anticlk;
//				}
//				else if (p2.y < box[i][4]) {
//					p2.y = box[i][4];
//					p2.x = (p2.y - interceptyx) / m_anticlk;
//				}
//				temp1.x = p1.x - box[i][2];
//				temp1.y = p1.y - box[i][4];
//				temp2.x = p2.x - box[i][2];
//				temp2.y = p2.y - box[i][4];
//				temp = crossing(roi, h, w, temp1, temp2, m_anticlk, 1);
//				temp_spread = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
//				if (temp > cross) {
//					cross = temp;
//					skew = m_anticlk;
//					spreadYX = temp_spread;
//				}
//				if (temp == cross && temp_spread > spreadYX) {
//					skew = m_anticlk;
//					spreadYX = temp_spread;
//				}
//				//line(cc_img, p1, p2, Scalar(150, 255, 0), 1.5, 8);
//			}
//
//			//------------Display the maximum spread line for Y on X box...
//
//			interceptyx = (i_start + center[0]) - skew*(j_start + center[1]);
//			p1.x = box[i][2];
//			p1.y = skew*box[i][2] + interceptyx;
//			p2.x = box[i][3];
//			p2.y = skew*box[i][3] + interceptyx;
//			line(cc_img, p1, p2, Scalar(255, 0, 255), 1.5, 8);
//		}
//		else {
//			Point temp1, temp2;
//			double temp_spread;
//			line(cc_img, p3, p4, Scalar(0, 150, 255), 1.5, 8);
//			regSlope[i][0] = -1;
//			regSlope[i][1] = slopexy;
//			skew = slopexy;
//			temp1.x = p3.x - box[i][2];
//			temp1.y = p3.y - box[i][4];
//			temp2.x = p4.x - box[i][2];
//			temp2.y = p4.y - box[i][4];
//			cross = crossing(roi, h, w, temp1, temp2, slopexy, -1);
//			for (int k = 0; k < 45; k += 5) {
//				rot_slope = tan(k * CV_PI / 180);
//
//				//Determine the clockwise rotations...
//				m_anticlk = (slopexy - rot_slope) / (1 + slopexy*rot_slope);
//				interceptxy = (j_start + center[1]) - m_anticlk*(i_start + center[0]);
//				p1.y = box[i][4];
//				p1.x = m_anticlk*box[i][4] + interceptxy;
//				if (p1.x < box[i][2]) {
//					p1.x = box[i][2];
//					p1.y = (p1.x - interceptxy) / m_anticlk;
//				}
//				else if (p1.x >= box[i][3]) {
//					p1.x = box[i][3] - 1;
//					p1.y = (p1.x - interceptxy) / m_anticlk;
//				}
//
//				p2.y = box[i][5];
//				p2.x = m_anticlk*box[i][5] + interceptxy;
//				if (p2.x >= box[i][3]) {
//					p2.x = box[i][3] - 1;
//					p2.y = (p2.x - interceptxy) / m_anticlk;
//				}
//				else if (p2.x < box[i][2]) {
//					p2.x = box[i][2];
//					p2.y = (p2.x - interceptxy) / m_anticlk;
//				}
//				temp1.x = p1.x - box[i][2];
//				temp1.y = p1.y - box[i][4];
//				temp2.x = p2.x - box[i][2];
//				temp2.y = p2.y - box[i][4];
//				temp = crossing(roi, h, w, temp1, temp2, m_anticlk, -1);
//				temp_spread = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
//				if (temp > cross) {
//					cross = temp;
//					skew = m_anticlk;
//					spreadXY = temp_spread;
//				}
//				if (temp == cross && temp_spread > spreadXY) {
//					skew = m_anticlk;
//					spreadXY = temp_spread;
//				}
//				//line(cc_img, p1, p2, Scalar(0, 150, 255), 1.5, 8);
//
//				//Determine the anti-clockwise rotation....
//				m_clk = (slopexy + rot_slope) / (1 - slopexy*rot_slope);
//				interceptxy = (j_start + center[1]) - m_clk*(i_start + center[0]);
//				p1.y = box[i][4];
//				p1.x = m_clk*box[i][4] + interceptxy;
//				if (p1.x < box[i][2]) {
//					p1.x = box[i][2];
//					p1.y = (p1.x - interceptxy) / m_clk;
//				}
//				else if (p1.x >= box[i][3]) {
//					p1.x = box[i][3] - 1;
//					p1.y = (p1.x - interceptxy) / m_clk;
//				}
//
//				p2.y = box[i][5];
//				p2.x = m_clk*box[i][5] + interceptxy;
//				if (p2.x >= box[i][3]) {
//					p2.x = box[i][3] - 1;
//					p2.y = (p2.x - interceptxy) / m_clk;
//				}
//				else if (p2.x < box[i][2]) {
//					p2.x = box[i][2];
//					p2.y = (p2.x - interceptxy) / m_clk;
//				}
//				temp1.x = p1.x - box[i][2];
//				temp1.y = p1.y - box[i][4];
//				temp2.x = p2.x - box[i][2];
//				temp2.y = p2.y - box[i][4];
//				temp = crossing(roi, h, w, temp1, temp2, m_clk, -1);
//				temp_spread = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
//				if (temp > cross) {
//					cross = temp;
//					skew = m_clk;
//					spreadXY = temp_spread;
//				}
//				if (temp == cross && temp_spread > spreadXY) {
//					skew = m_clk;
//					spreadXY = temp_spread;
//				}
//				//line(cc_img, p1, p2, Scalar(0, 255, 150), 1.5, 8);
//			}
//
//			//------------Display the maximum spread line for X on Y box...
//
//			interceptxy = (j_start + center[1]) - skew*(i_start + center[0]);
//			p3.y = box[i][4];
//			p3.x = skew*box[i][4] + interceptxy;
//			p4.y = box[i][5];
//			p4.x = skew*box[i][5] + interceptxy;
//			line(cc_img, p3, p4, Scalar(255, 0, 255), 1.5, 8);
//		}
//		
//		
//		Vec3b color;
//		color.val[0] = 0;
//		color.val[1] = 255;
//		color.val[2] = 0;
//		cc_img.at<Vec3b>(i_start + center[0], j_start + center[1]) = color;
//
//
//		//Delete ROI...
//		for (int k = 0; k < h; k++)
//			delete[] roi[k];
//		delete[] roi;
//	}
//
//	imshow("CentroidLocations", cc_img);
//	imwrite("Data/ICDAR/Data218/centroids218(new).tif", cc_img);
//}
//
///*----------------------------------------------------------------*/
///*     Re-cluster the CC wrt the nearest CC with longest line     */
///*----------------------------------------------------------------*/
//void cluster(int **img,  double **Data, int **label, int **change) {
//	double avg_length = 0;
//	int N = 0;
//	for (int i = 0; i < labelNum; i++) {
//		if (Data[i][4] == 0) continue; //length is 0...
//		avg_length += Data[i][4];
//		N++;
//	}
//	avg_length /= N;
//	std::cout << "\nAverage length of line = " << avg_length;
//	for (int i = 0; i < labelNum; i++) {
//		if (Data[i][4] == 0)   //length of line...
//			continue;
//		int **roi;
//		int h, w, centroid[2], i_min, i_max, j_min, j_max;
//		double slope, length;
//		centroid[0] = Data[i][1]; //Y co-ordinate
//		centroid[1] = Data[i][0]; //X co-ordinate
//		slope = Data[i][2];
//		length = Data[i][4];
//		i_min = centroid[0] - int(1.25*length);
//		i_max = centroid[0] + int(1.25*length);
//		j_min = centroid[1] - int(1.25*length);
//		j_max = centroid[1] + int(1.25*length);
//		if (i_min < 0) i_min = 0;
//		if (i_max >= r) i_max = r - 1;
//		if (j_min < 0) j_min = 0;
//		if (j_max >= c) j_max = c - 1;
//		h = i_max - i_min + 1;
//		w = j_max - j_min + 1;
//		roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				if (img[i_min + k][j_min + l] == 255)
//					roi[k][l] = 255;
//				else
//					roi[k][l] = 0;
//			}
//		}
//
//		//-----------------------Locate nearest neighbor-------------------------
//		int near[2], color_code, c[9] = { 0 }, label_no, max_color, max_count;
//		near[0] = centroid[0];
//		near[1] = centroid[1];
//		double line = 0, distance = 1000;
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				if (roi[k][l] == 255) continue; //background...
//				label_no = label[i_min + k][j_min + l];
//				//if (label_no == i) continue;   //reference CC...
//				c[int(Data[label_no][3])]++;
//				/*double dist = sqrt((centroid[0] - k)*(centroid[0] - k) + (centroid[1] - l)*(centroid[1] - l));
//				if (dist < distance && Data[label_no][4] > line) {
//					dist = distance;
//					line = Data[label_no][4];
//					color_code = Data[label_no][3];
//				}*/
//			}
//		}
//
//		//---------Find largest cluster---------
//		max_color = 1; max_count = c[1];
//		for (int n = 2; n < 9; n++) {
//			if (c[n] > max_count) {
//				max_count = c[n];
//				max_color = n;
//			}
//		}
//		color_code = max_color;
//		if (color_code != Data[i][3]) {
//			change[i][0] = 1;
//			change[i][1] = color_code;
//		}
//
//		//------------------Delete ROI--------------------
//		for (int k = 0; k < h; k++)
//			delete[] roi[k];
//		delete[] roi;
//	}
//
//	//--------------------Re-color CC's----------------------
//	for (int i = 0; i < labelNum; i++) {
//		if (Data[i][4] == 0) continue;
//		if (change[i][0] == 1)
//			Data[i][3] = change[i][1];
//	}
//
//}
//
///*----------------------------------------------------------------*/
///*            Blur the CC using rotated Gaussian filter           */
///*----------------------------------------------------------------*/
//void GaussBlur(int **img, int **label, int **box, double **Data) {
//	Mat blurimg;
//	blurimg.create(r, c, CV_8UC1);
//	//cvtColor(blurimg, blurimg, CV_GRAY2BGR);
//
//	int **roi = new int*[r];
//	for (int i = 0; i < r; i++) {
//		roi[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			if (Data[label[i][j]][3] == 1 || Data[label[i][j]][3] == 5) {
//				roi[i][j] = 0;
//			}
//			else 
//				roi[i][j] = 255;
//		}
//	}
//	//blurROI(roi, r, c, 45);
//	arr2mat(roi, blurimg);
//	GaussianBlur(blurimg, blurimg, Size(41, 7), 0, 0);
//	mat2arr(blurimg, roi);
//	Otsu(roi, r, c, roi);
//	arr2mat(roi, blurimg);
//	imshow("Only1&5", blurimg);
//	imwrite("Data/ICDAR/Data218/mask218.tif", blurimg);
//
//}
//
///*----------------------------------------------------------------*/
///*       Find the Partition lines through centroid in CC         */
///*----------------------------------------------------------------*/
//int verticalShift(int **roi, int h, int w, int points[4]) {
//	int count_right = 0, count_left = 0, pixels = 0, prev, count;
//	for (int j = 0; j < points[0]; j++) {
//		for (int i = 0; i < h; i++) {
//			count_left += (255 - roi[i][j]) / 255;
//		}
//	}
//	for (int j = points[0] + 1; j < w; j++) {
//		for (int i = 0; i < h; i++) {
//			count_right += (255 - roi[i][j]) / 255;
//		}
//	}
//	for (int i = 0; i < h; i++)
//		pixels += (255 - roi[i][points[0]]) / 255;
//	//std::cout << "left = " << count_left << "   right = " << count_right;
//	if (count_left == count_right) {
//		return(pixels);
//	}
//
//	if (count_left < count_right) {
//		while (count_left < count_right) {
//			points[0]++;
//			if (points[0] >= w)
//				break;
//			count_left += pixels;
//			pixels = 0;
//			for (int i = 0; i < h; i++) {
//				pixels += (255 - roi[i][points[0]]) / 255;
//			}
//			count_right -= pixels;
//		} 
//		//if (count_left > count_right) points[0]--;
//		//std::cout << "left = " << count_left << "   right = " << count_right;
//		pixels = 0;
//		for (int i = 0; i < h; i++)
//			pixels += (255 - roi[i][points[0]]) / 255;
//		return(pixels);
//	}
//	
//	if (count_left > count_right) {
//		while (count_left > count_right) {
//			points[0]--;
//			if (points[0] < 0)
//				break;
//			count_right += pixels;
//			pixels = 0;
//			for (int i = 0; i < h; i++) {
//				pixels += (255 - roi[i][points[0]]) / 255;
//			}
//			count_left -= pixels;
//		}
//		//if (count_left < count_right) points[0]++;
//		//std::cout << "left = " << count_left << "   right = " << count_right;
//		pixels = 0;
//		for (int i = 0; i < h; i++)
//			pixels += (255 - roi[i][points[0]]) / 255;
//		return(pixels);
//	}
//
//}
//
//int horizontalShift(int **roi, int h, int w, int points[4]) {
//	int count_up = 0, count_down = 0, pixels = 0, prev, count;
//	for (int i = 0; i < points[1]; i++) {
//		for (int j = 0; j < w; j++) {
//			count_up += (255 - roi[i][j]) / 255;
//		}
//	}
//	for (int i = points[1] + 1; i < h; i++) {
//		for (int j = 0; j < w; j++) {
//			count_down += (255 - roi[i][j]) / 255;
//		}
//	}
//	for (int j = 0; j < w; j++)
//		pixels += (255 - roi[points[1]][j])/255;
//	//std::cout << "up = " << count_up << "   down = " << count_down;
//	if (count_up == count_down) {
//		return(pixels);
//	}
//
//	if (count_up < count_down) {
//		while (count_up < count_down) {
//			points[1]++;
//			if (points[1] >= h)
//				break;
//			count_up += pixels;
//			pixels = 0;
//			for (int j = 0; j < w; j++) {
//				pixels += (255 - roi[points[1]][j]) / 255;
//			}
//			count_down -= pixels;
//		}
//		//if (count_up > count_down) points[1]--;
//		//std::cout << "up = " << count_up << "   down = " << count_down;
//		pixels = 0;
//		for (int j = 0; j < w; j++)
//			pixels += (255 - roi[points[1]][j]) / 255;
//		return(pixels);
//	}
//
//	if (count_up > count_down) {
//		while (count_up > count_down) {
//			points[1]--;
//			if (points[1] < 0)
//				break;
//			count_down += pixels;
//			pixels = 0;
//			for (int j = 0; j < w; j++) {
//				pixels += (255 - roi[points[1]][j]) / 255;
//			}
//			count_up -= pixels;
//		}
//		//if (count_up > count_down) points[1]++;
//		//std::cout << "up = " << count_up << "   down = " << count_down;
//		pixels = 0;
//		for (int j = 0; j < w; j++)
//			pixels += (255 - roi[points[1]][j]) / 255;
//		return(pixels);
//		/*for (int j = 0; j < w; j++) {
//			if (prev == 255 && roi[points[1]][j] == 0) count++;
//			prev = roi[points[1]][j];
//		}
//		return(count);*/
//	}
//
//}
//
//int DiagShift(int **roi, int h, int w, int points[4]) {
//	int count_full = 0, count_left = 0, count_right = 0, k = 0, pixels = 0, prev, count;
//	for (int i = 0; i < h; i++) {
//		for (int j = 0; j < w; j++)
//			count_full += (255 - roi[i][j]) / 255;
//	}
//	for (int i = points[1]; i < h; i++) {
//		for (int j = 0; j < points[0] + k; j++) {
//			count_left += (255 - roi[i][j]) / 255;
//		}
//		if (points[0] + k < w - 1) k++;
//	}
//	for (int i = points[1]; i <= points[3]; i++) {
//		pixels += (255 - roi[i][points[0] + k]) / 255;
//		k++;
//	}
//	count_right = count_full - count_left - pixels; 
//	//std::cout << "up = " << count_left << "   down = " << count_full - count_left;
//	if (count_left == count_right) {
//		return(pixels);
//	}
//	if (count_left < count_right) {
//		do {
//			count_left += pixels;
//			if (points[1] == 0) {
//				points[0]++;
//				if (points[0] >= w)	break;
//			}
//			if (points[0] == 0 ) {
//				points[1]--;
//				if (points[1] < 0) points[1] = 0;
//			}
//			if (points[2] == w - 1) {
//				points[3]--;
//				if (points[3] < 0)	break;
//			}
//			if (points[3] == h - 1) {
//				points[2]++;
//				if (points[2] >= w) points[2] = w - 1;
//			}
//			pixels = 0; k = 0;
//			for (int i = points[1]; i <= points[3]; i++) {
//				pixels += (255 - roi[i][points[0] + k]) / 255;
//				k++;
//			}
//			count_right -= pixels;
//		} while (count_left < count_right);
//		//std::cout << "up = " << count_left << "   down = " << count_full - count_left;
//		/*pixels = 0; k = 0;
//		for (int i = points[1]; i <= points[3]; i++) {
//			pixels += (255 - roi[i][points[0] + k]) / 255;
//			k++;
//		}*/
//		return(count_full - count_left - count_right);
//	}
//	if (count_left > count_right) {
//		do {
//			count_right += pixels;
//			if (points[0] == 0) {
//				points[1]++;
//				if (points[1] >= h)	break;
//			}
//			if (points[1] == 0) {
//				points[0]--;
//				if (points[0] < 0) points[0] = 0;
//			}
//			if (points[3] == h - 1) {
//				points[2]--;
//				if (points[2] < 0)	break;
//			}
//			if (points[2] == w - 1) {
//				points[3]++;
//				if (points[3] >= h) points[3] = h - 1;
//			}
//			pixels = 0; k = 0;
//			for (int i = points[1]; i <= points[3]; i++) {
//				pixels += (255 - roi[i][points[0] + k]) / 255;
//				k++;
//			}
//			count_left -= pixels;
//		} while (count_left > count_right);
//		//std::cout << "up = " << count_left << "   down = " << count_full - count_left;
//		/*pixels = 0; k = 0;
//		for (int i = points[1]; i <= points[3]; i++) {
//			pixels += (255 - roi[i][points[0] + k]) / 255;
//			k++;
//		}*/
//		return(count_full - count_left - count_right);
//	}
//	
//}
//
//int OffDiagShift(int **roi, int h, int w, int points[4]) {
//	int count_full = 0, count_left = 0, count_right = 0, k = 0, pixels = 0, prev ,count;
//	//std::cout << points[0] << " , " << points[1] << " , " << points[2] << " , " << points[3];
//	for (int i = 0; i < h; i++) {
//		for (int j = 0; j < w; j++)
//			count_full += (255 - roi[i][j]) / 255;
//	}
//	for (int i = points[3]; i >= 0; i--) {
//		for (int j = 0; j < points[2] + k; j++) {
//			count_left += (255 - roi[i][j]) / 255;
//		}
//		if (points[2] + k < w - 1) k++;
//	}
//	for (int i = points[3]; i >= points[1]; i--) {
//		pixels += (255 - roi[i][points[2] + k]) / 255;
//		k++;
//	}
//	count_right = count_full - count_left - pixels;
//	//std::cout << "up = " << count_left << "   down = " << count_full - count_left;
//	if (count_left == count_right) {
//		return(pixels);
//	}
//	if (count_left < count_right) {
//		do {
//			count_left += pixels;
//			if (points[0] == w - 1) {
//				points[1]++;
//				if (points[1] >= h)	break;
//			}
//			if (points[1] == 0) {
//				points[0]++;
//				if (points[0] >= w) points[0] = w - 1;
//			}
//			if (points[3] == h - 1) {
//				points[2]++;
//				if (points[2] >= w)	break;
//			}
//			if (points[2] == 0) {
//				points[3]++;
//				if (points[3] >= h) points[3] = h - 1;
//			}			
//			pixels = 0; k = 0;
//			for (int i = points[3]; i >= points[1]; i--) {
//				pixels += (255 - roi[i][points[2] + k]) / 255;
//				k++;
//			}
//			count_right -= pixels;
//		} while (count_left < count_right);
//		//std::cout << "up = " << count_left << "   down = " << count_full - count_left; 
//		/*pixels = 0; k = 0;
//		for (int i = points[3]; i >= points[1]; i--) {
//			pixels += (255 - roi[i][points[2] + k]) / 255;
//			k++;
//		}*/
//		return(count_full - count_left - count_right);
//	}
//	if (count_left > count_right) {
//		do {
//			count_right += pixels;
//			if (points[1] == 0) {
//				points[0]--;
//				if (points[0] < 0)	break;
//			}
//			if (points[0] == w - 1) {
//				points[1]--;
//				if (points[1] < 0 ) points[1] = 0;
//			}
//			if (points[2] == 0) {
//				points[3]--;
//				if (points[3] < 0)	break;
//			}
//			if (points[3] == h - 1) {
//				points[2]--;
//				if (points[2] < 0) points[2] = 0;
//			}
//			
//			pixels = 0; k = 0;
//			for (int i = points[3]; i >= points[1]; i--) {
//				pixels += (255 - roi[i][points[2] + k]) / 255;
//				k++;
//			}
//			count_left -= pixels;
//		} while (count_left > count_right);
//		//std::cout << "up = " << count_left << "   down = " << count_full - count_left;
//		/*pixels = 0; k = 0;
//		for (int i = points[3]; i >= points[1]; i--) {
//			pixels += (255 - roi[i][points[2] + k]) / 255;
//			k++;
//		}*/
//		return(count_full - count_left - count_right);
//	}
//
//}
//
//int crossing(int **roi, int h, int w, Point p1, Point p2, double slope, double intercept) {
//	int x, y, prev = 255, curr = 0, count = 0 ;
//	count = (255 - roi[p1.y][p1.x]) / 255;
//	for (int x = p1.x; x < p2.x; x++) {
//		y = (int)(slope*x + intercept);
//		if (y < 0) y = 0;
//		if (y >= h) y = h - 1;
//		count += (255 - roi[y][x]) / 255;
//	}
//	count+= (255 - roi[p2.y][p2.x]) / 255;
//	//std::cout << " * " << count;
//	return(count);
//}
//
//void getLineNo(int **img, int **label, int **box, double **Data, double **connect, int **lineNo) {
//	int *pre_list = new int[labelNum];
//	int *succ_list = new int[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1 || Data[i][5] == 0) //Not a significant component...
//			continue;
//		succ_list[i] = connect[i][1];
//		pre_list[i] = -1;
//	}
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1 || Data[i][5] == 0) //Not a significant component...
//			continue;
//		if (connect[i][1]>0)
//			pre_list[int(connect[i][1])] = i;
//	}
//	
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1 || Data[i][5] == 0) //Not a significant component...
//			continue;
//		int curr = i, prev = pre_list[i], next = succ_list[i];
//		if (prev == -1 && next <= 0) {   //single word... 
//			line_no++;
//			lineNo[curr][0] = line_no;
//			lineNo[curr][1] = 2;
//			continue;
//		}
//		if (!(prev == -1 && next > 0)) continue; // not beginning of line...
//		else {
//			line_no++;
//			lineNo[curr][1] = 1; //mark the first label no. of the line...
//			while (next != -1) {
//				lineNo[curr][0] = line_no;
//				prev = curr, curr = next, next = succ_list[curr];
//			}
//			lineNo[curr][0] = line_no;
//			lineNo[curr][1] = 3; //mark the last label no. of the line...
//		}
//	}
//	std::cout << "\nNo. of lines = " << line_no;
//	for (int i = 0; i < labelNum; i++)
//		std::cout << "\nLabel " << i << " : Line " << lineNo[i][0];
//}
//
//void partition(Mat scanimg, int **img, Mat cc_img, Mat bin_img_org, int **label, double **boundingBox, double **diag_box, int **box, double **Data, int **lineNo, double **connect) {
//	Mat bin_img, img_org, cluster_img, new_cluster, sub_im;
//	int **cluster_im;
//	Point p1, p2, p3, p4, p5, p6, p7, p8, p, q, P, Q;
//	Point **ends;
//	int *angles, *clusters;
//	img_org.create(r, c, CV_8UC1);
//	cvtColor(bin_img_org, bin_img, CV_GRAY2BGR);
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			Vec3b color, prev;
//			prev = bin_img.at<Vec3b>(i, j);
//			if (prev.val[0] == 255 && prev.val[1] == 255 && prev.val[2] == 255) continue;
//			color.val[0] = 150;
//			color.val[1] = 150;
//			color.val[2] = 150;
//			bin_img.at<Vec3b>(i, j) = color;
//		}
//	}
//
//	ends = new Point*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		ends[i] = new Point[2];
//		for (int j = 0; j < 2; j++) {
//			ends[i][0] = Point(0, 0);
//			ends[i][1] = Point(0, 0);
//		}
//	}
//
//	clusters = new int[labelNum];
//	angles = new int[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		angles[i] = -1;
//		clusters[i] = 0;
//	}
//	Data = new double*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		Data[i] = new double[6];
//		for (int j = 0; j < 6; j++) {
//			Data[i][j] = 0;
//		}
//	}
//	cluster_im = new int*[r];
//	for (int i = 0; i < r; i++) {
//		cluster_im[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			cluster_im[i][j] = 255;
//		}
//	}
//	cluster_img.create(r, c, CV_8UC1);
//	arr2mat(cluster_im, cluster_img);
//	cvtColor(cluster_img, cluster_img, CV_GRAY2BGR);
//
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) //Not a significant component...
//			continue;
//		int center[2];
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//
//		//--------------------------------Select ROI and determine centroid for black pixels---------------
//
//		int **roi, points[4];
//		int *X, *Y;
//		int count = 0, N = 0, i_start = box[i][4], j_start = box[i][2], max_line[2], pixels = 0;
//		double  slopeyx = 0, slopexy = 0, interceptyx = 0, interceptxy = 0, spreadYX, spreadXY, slope;
//		roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				if (label[i_start + k][j_start + l] == i) {
//					roi[k][l] = img[i_start + k][j_start + l];
//					N++;
//				}
//				else
//					roi[k][l] = 255;
//			}
//		}
//
//		X = new int[N];
//		Y = new int[N];
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				if (roi[k][l] == 0) {
//					X[count] = l;
//					Y[count] = k;
//					count++;
//				}
//			}
//		}
//		center[0] = (int)Mean(Y, N);
//		center[1] = (int)Mean(X, N);
//
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				if (label[i_start + k][j_start + l] == i) {
//					Vec3b color;
//					color.val[0] = 150;
//					color.val[1] = 150;
//					color.val[2] = 150;
//					cc_img.at<Vec3b>(i_start + k, j_start + l) = color;
//				}
//			}
//		}
//
//		//-----------------------------Determine the line with maximum text pixels--------------------
//		if (boundingBox[i][0] == 1) { //only compare 0 and 90 degrees if upright box...
//			//Set the horizontal line...
//			p1.x = box[i][2];
//			p1.y = i_start + center[0];
//			p2.x = box[i][3];
//			p2.y = i_start + center[0];
//			//line(cc_img, p1, p2, Scalar(255, 0, 255), 1.5, 8);
//			points[0] = 0;
//			points[1] = center[0];
//			points[2] = w - 1;
//			points[3] = center[0];
//			//std::cout << "\nlabel " << i << ": ";
//			max_line[0] = 1;
//			max_line[1] = horizontalShift(roi, h, w, points);
//			p1.y = i_start + points[1];
//			p2.y = i_start + points[1];
//			//line(cc_img, p1, p2, Scalar(255, 0, 255), 1.5, 8);
//			p = p1;
//			q = p2;
//			angles[i] = 0;
//			//std::cout << " horizontal = " << max_line[1];
//
//			//Set the vertical line...
//			p3.x = j_start + center[1];
//			p3.y = box[i][4];
//			p4.x = j_start + center[1];
//			p4.y = box[i][5];
//			//line(cc_img, p3, p4, Scalar(255, 0, 0), 1.5, 8);
//			points[0] = center[1];
//			points[1] = 0;
//			points[2] = center[1];
//			points[3] = h - 1;
//			//std::cout << "\nlabel " << i << ": ";
//			pixels = verticalShift(roi, h, w, points);
//			p3.x = j_start + points[0];
//			p4.x = j_start + points[0];
//			//line(cc_img, p3, p4, Scalar(255, 0, 0), 1.5, 8);
//			if (pixels > max_line[1]) {
//				max_line[1] = pixels;
//				max_line[0] = 2;
//				p = p3;
//				q = p4;
//				angles[i] = 90;
//			}
//			//std::cout << " vertical = " << pixels;
//		}
//		
//		else if (boundingBox[i][0] == 2) {
//			//Set the diagonal line...
//			interceptyx = (i_start + center[0]) - (j_start + center[1]);
//			p5.x = box[i][2];
//			p5.y = p5.x + interceptyx;
//			if (p5.y < box[i][4]) {
//				p5.y = box[i][4];
//				p5.x = p5.y - interceptyx;
//			}
//			if (p5.y > box[i][5]) {
//				p5.y = box[i][5];
//				p5.x = p5.y - interceptyx;
//			}
//			p6.x = box[i][3];
//			p6.y = p6.x + interceptyx;
//			if (p6.y < box[i][4]) {
//				p6.y = box[i][4];
//				p6.x = p6.y - interceptyx;
//			}
//			if (p6.y > box[i][5]) {
//				p6.y = box[i][5];
//				p6.x = p6.y - interceptyx;
//			}
//			//line(cc_img, p5, p6, Scalar(0, 150, 255), 1.5, 8);
//			points[0] = p5.x - box[i][2];
//			points[1] = p5.y - box[i][4];
//			points[2] = p6.x - box[i][2];
//			points[3] = p6.y - box[i][4];
//			//std::cout << "\nlabel " << i << ": ";
//			max_line[0] = 3;
//			max_line[1] = DiagShift(roi, h, w, points);
//			p5.x = j_start + points[0];
//			p5.y = i_start + points[1];
//			p6.x = j_start + points[2];
//			p6.y = i_start + points[3];
//			p = p5;
//			q = p6;
//			angles[i] = 45;
//			//line(cc_img, p5, p6, Scalar(0, 255, 150), 1.5, 8);
//			/*if (pixels > max_line[1]) {
//				max_line[1] = pixels;
//				max_line[0] = 3;
//				p = p5;
//				q = p6;
//				angles[i] = 45;
//			}*/
//			//std::cout << " diagonal = " << pixels;		
//
//			//Set the off-diagonal line...
//			interceptxy = (i_start + center[0]) + (j_start + center[1]);
//			p7.x = box[i][3];
//			p7.y = interceptxy - p7.x;
//			if (p7.y < box[i][4]) {
//				p7.y = box[i][4];
//				p7.x = interceptxy - p7.y;
//			}
//			if (p7.y > box[i][5]) {
//				p7.y = box[i][5];
//				p7.x = interceptxy - p7.y;
//			}
//			p8.x = box[i][2];
//			p8.y = interceptxy - p8.x;
//			if (p8.y < box[i][4]) {
//				p8.y = box[i][4];
//				p8.x = interceptxy - p8.y;
//			}
//			if (p8.y > box[i][5]) {
//				p8.y = box[i][5];
//				p8.x = interceptxy - p8.y;
//			}
//			//line(cc_img, p7, p8, Scalar(255, 150, 0), 1.5, 8);
//			points[0] = p7.x - box[i][2];
//			points[1] = p7.y - box[i][4];
//			points[2] = p8.x - box[i][2];
//			points[3] = p8.y - box[i][4];
//			//std::cout << "\nlabel " << i << ": ";
//			pixels = OffDiagShift(roi, h, w, points);
//			p7.x = j_start + points[0];
//			p7.y = i_start + points[1];
//			p8.x = j_start + points[2];
//			p8.y = i_start + points[3];
//			//line(cc_img, p7, p8, Scalar(255, 150, 0), 1.5, 8);
//			if (pixels > max_line[1]) {
//				max_line[1] = pixels;
//				max_line[0] = 4;
//				p = p7;
//				q = p8;
//				angles[i] = 135;
//			}
//			//std::cout << " off-diagonal = " << pixels;
//		}
//		
//
//		line(cc_img, p, q, Scalar(255, 0, 150), 1.5, 8);
//		center[0] = (int)((p.y + q.y) / 2.0) - i_start;
//		center[1] = (int)((p.x + q.x) / 2.0) - j_start;
//		if (max_line[0] == 1) {  //horizontal...
//			center[0] = p.y - i_start;
//			center[1] = (int)((p.x + q.x) / 2.0) - j_start;
//		}
//		if (max_line[0] == 2) {  //vertical...
//			center[0] = (int)((p.y + q.y) / 2.0) - i_start;
//			center[1] = p.x - j_start;
//		}
//
//
//		//----------------------------Rotate the best approximation by 45degrees on both sides for the exact one------------------
//
//		double theta_m, clk_m, anticlk_m, intercept, max_slope;
//		int max_angle = angles[i];
//		pixels = 0;
//		P = p;
//		Q = q;
//		for (int k = 5; k < 90; k += 5) {
//			theta_m = tan(k*CV_PI / 180);
//			if (max_line[0] == 1) {  //horizontal...
//				clk_m = -theta_m;
//				anticlk_m = theta_m;
//			}
//			if (max_line[0] == 2) {  //vertical...
//				clk_m = 1 / theta_m;
//				anticlk_m = -1 / theta_m;
//			}
//			if (max_line[0] == 3) {  //diagonal...
//				clk_m = (1 - theta_m) / (1 + theta_m);
//				anticlk_m = (1 + theta_m) / (1 - theta_m);
//			}
//			if (max_line[0] == 4) {  //off-diagonal...
//				clk_m = (1 + theta_m) / (theta_m - 1);
//				anticlk_m = (theta_m - 1) / (1 + theta_m);
//			}
//			//std::cout << "\nclk_m = " << clk_m << ", anticlk_m = " << anticlk_m;
//			intercept = (center[0]) - clk_m*(center[1]);
//			p.x = 0;
//			p.y = intercept;
//			//std::cout << "\n(x1, y1) = " << p.x << " , " << p.y ;
//			if (p.y < 0) {
//				p.y = 0;
//				p.x = (-intercept) / clk_m;
//			}
//			else if (p.y >= h) {
//				p.y = h - 1;
//				p.x = ((h - 1) - intercept) / clk_m;
//			}
//			q.x = w - 1;
//			q.y = (w - 1) * clk_m + intercept;
//			//std::cout << "   (x2, y2) = " << q.x << ", " << q.y;
//			if (q.y < 0) {
//				q.y = 0;
//				q.x = (-intercept) / clk_m;
//			}
//			else if (q.y >= h) {
//				q.y = h - 1;
//				q.x = ((h - 1) - intercept) / clk_m;
//			}
//			pixels = crossing(roi, h, w, p, q, clk_m, intercept);
//			//std::cout << "\nLabel " << i << " : " << pixels;
//			//std::cout << "\n(x1, y1) = " << p.x << " , " << p.y << "   (x2, y2) = " << q.x << ", " << q.y;
//			p.x += j_start;
//			p.y += i_start;
//			q.x += j_start;
//			q.y += i_start;
//			/*if (i == 166)
//				line(cc_img, p, q, Scalar(150, 0, 255), 1.5, 8);*/
//			if (pixels >= max_line[1]) {
//				max_line[1] = pixels;
//				P = p;
//				Q = q;
//				max_slope = clk_m;
//				max_angle = angles[i] - k;
//				if (max_angle < 0) max_angle = 180 + max_angle;
//			}
//			//line(cc_img, p, q, Scalar(150, 0, 255), 1.5, 8);
//
//			intercept = (center[0]) - anticlk_m*(center[1]);
//			p.x = 0;
//			p.y = intercept;
//			//std::cout << "\n(x1, y1) = " << p.x << " , " << p.y;
//			if (p.y < 0) {
//				p.y = 0;
//				p.x = (-intercept) / anticlk_m;
//			}
//			else if (p.y >= h) {
//				p.y = h - 1;
//				p.x = ((h - 1) - intercept) / anticlk_m;
//			}
//			q.x = w - 1;
//			q.y = (w - 1) * anticlk_m + intercept;
//			//std::cout << "   (x2, y2) = " << q.x << ", " << q.y;
//			if (q.y < 0) {
//				q.y = 0;
//				q.x = (-intercept) / anticlk_m;
//			}
//			else if (q.y >= h) {
//				q.y = h - 1;
//				q.x = ((h - 1) - intercept) / anticlk_m;
//			}
//			pixels = crossing(roi, h, w, p, q, anticlk_m, intercept);
//			//std::cout << "\nLabel " << i << " : " << pixels;
//			//std::cout << "\n(x1, y1) = " << p.x << " , " << p.y << "   (x2, y2) = " << q.x << ", " << q.y;
//			p.x += j_start;
//			p.y += i_start;
//			q.x += j_start;
//			q.y += i_start;
//			/*if (i == 166)
//				line(cc_img, p, q, Scalar(255, 0, 150), 1.5, 8);*/
//			if (pixels >= max_line[1]) {
//				max_line[1] = pixels;
//				P = p;
//				Q = q;
//				max_slope = anticlk_m;
//				max_angle = angles[i] + k;
//				if (max_angle > 180) max_angle = max_angle - 180;
//				//std::cout << "*this*";
//			}
//			//line(cc_img, p, q, Scalar(255, 0, 150), 1.5, 8);
//		}
//		if ((abs(max_angle - angles[i]) >= 35 && abs(max_angle - angles[i]) <= 70) || (abs(max_angle - angles[i]) >= 125 && abs(max_angle - angles[i]) <= 160)) {
//			//line(cc_img, P, Q, Scalar(255, 0, 0), 2, 8);
//			if (boundingBox[i][0] == 1) {
//				boundingBox[i][0] = 2;
//				boundingBox[i][1] = diag_box[i][0];
//				boundingBox[i][2] = diag_box[i][1];
//				boundingBox[i][3] = diag_box[i][2];
//				boundingBox[i][4] = diag_box[i][3];
//				boundingBox[i][5] = diag_box[i][4];
//				boundingBox[i][6] = diag_box[i][5];
//				boundingBox[i][7] = diag_box[i][6];
//				boundingBox[i][8] = diag_box[i][7];
//				line(cc_img, Point(boundingBox[i][2], boundingBox[i][1]), Point(boundingBox[i][6], boundingBox[i][5]), Scalar(255, 0, 0), 1.5, 8);
//				line(cc_img, Point(boundingBox[i][6], boundingBox[i][5]), Point(boundingBox[i][4], boundingBox[i][3]), Scalar(255, 0, 0), 1.5, 8);
//				line(cc_img, Point(boundingBox[i][4], boundingBox[i][3]), Point(boundingBox[i][8], boundingBox[i][7]), Scalar(255, 0, 0), 1.5, 8);
//				line(cc_img, Point(boundingBox[i][8], boundingBox[i][7]), Point(boundingBox[i][2], boundingBox[i][1]), Scalar(255, 0, 0), 1.5, 8);
//			}
//			else if (boundingBox[i][0] == 2) {
//				boundingBox[i][0] = 1;
//				boundingBox[i][1] = box[i][2]; // Left x
//				boundingBox[i][2] = box[i][3]; // Right x
//				boundingBox[i][3] = box[i][4]; // Top y
//				boundingBox[i][4] = box[i][5]; // Bottom y
//				line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(255, 0, 0), 1.5, 8);
//				line(cc_img, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(255, 0, 0), 1.5, 8);
//				line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(255, 0, 0), 1.5, 8);
//				line(cc_img, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(255, 0, 0), 1.5, 8);
//			}
//		}
//		angles[i] = max_angle;
//		/*std::cout << "\nLabel " << i << " : " << max_line[1];
//		std::cout << "\nLabel " << i << "(x1, y1) = " << P.x << " , " << P.y << "   (x2, y2) = " << Q.x << ", " << Q.y << " , slope = " << max_slope;*/
//		line(cc_img, P, Q, Scalar(0, 0, 0), 1.5, 8);
//
//		Vec3b color;
//		Point O;
//		O.x = j_start + center[1];
//		O.y = i_start + center[0];
//		if (angles[i] >= 0 && angles[i] <= 45) {
//			int lambda = angles[i] / 9;
//			clusters[i] = 1;
//			color.val[0] = int(255 * lambda / 5);
//			color.val[1] = 0;
//			color.val[2] = 255;
//		}
//		if (angles[i] >= 50 && angles[i] <= 90) {
//			int lambda = (angles[i] - 45) / 9;
//			clusters[i] = 2;
//			color.val[0] = 255;
//			color.val[1] = 0;
//			color.val[2] = int(255 * (5 - lambda) / 5);
//		}
//		if (angles[i] >= 95 && angles[i] <= 135) {
//			int lambda = (angles[i] - 90) / 9;
//			clusters[i] = 3;
//			color.val[0] = int(255 * (5 - lambda) / 5);
//			color.val[1] = int(255 * lambda / 5);
//			color.val[2] = 0;
//		}
//		if (angles[i] >= 140 && angles[i] <= 180) {
//			int lambda = (angles[i] - 135) / 9;
//			clusters[i] = 4;
//			color.val[0] = 0;
//			color.val[1] = int(127 * (5 - lambda) / 5);
//			color.val[2] = int(255 * lambda / 5);
//		}
//
//		//if (angles[i]  > 170 || angles[i] <= 5) { //around 0deg...
//		//	//std::cout << "\n slope=" << max_slope;
//		//	color.val[0] = 0;
//		//	color.val[1] = 0;
//		//	color.val[2] = 255;
//		//	clusters[i] = 1; //red --> 0 deg
//		//	cluster_img.at<Vec3b>(i_start + center[0], j_start + center[1]) = color;
//		//	circle(cluster_img, O, 1.5, color, -1, 8 , 0);
//		//}
//		//else if (angles[i] > 5 && angles[i] <= 30) { //around 22.5deg...
//		//	std::cout << "\n Label : " << i << " slope=" << angles[i];
//		//	color.val[0] = 255;
//		//	color.val[1] = 150;
//		//	color.val[2] = 0;
//		//	clusters[i] = 5; //cyan --> 22.5 deg
//		//	cluster_img.at<Vec3b>(i_start + center[0], j_start + center[1]) = color;
//		//	circle(cluster_img, O, 1.5, color, -1, 8, 0);
//		//}
//		//else if (angles[i] > 30 && angles[i]  <= 55) { //around 45deg...
//		//	//std::cout << "\n slope=" << max_slope;
//		//	color.val[0] = 255;
//		//	color.val[1] = 0;
//		//	color.val[2] = 0;
//		//	clusters[i] = 2; //blue  --> 45 deg
//		//	cluster_img.at<Vec3b>(i_start + center[0], j_start + center[1]) = color;
//		//	circle(cluster_img, O, 1.5, color, -1, 8, 0);
//		//}
//		//else if (angles[i] > 55 && angles[i] <= 80) { //around 67.5deg...
//		//	//std::cout << "\n slope=" << max_slope;
//		//	color.val[0] = 50;
//		//	color.val[1] = 0;
//		//	color.val[2] = 100;
//		//	clusters[i] = 6; //maroon --> 67.5deg
//		//	cluster_img.at<Vec3b>(i_start + center[0], j_start + center[1]) = color;
//		//	circle(cluster_img, O, 1.5, color, -1, 8, 0);
//		//}
//		//else if (angles[i] > 80 && angles[i]  <= 95) { //around 90deg...
//		//	//std::cout << "\n Label : " << i << " slope=" << angles[i];
//		//	color.val[0] = 0;
//		//	color.val[1] = 255;
//		//	color.val[2] = 0;
//		//	clusters[i] = 3; //green --> 90 deg
//		//	cluster_img.at<Vec3b>(i_start + center[0], j_start + center[1]) = color;
//		//	circle(cluster_img, O, 1.5, color, -1, 8, 0);
//		//}
//		//else if (angles[i] > 95 && angles[i] <= 120) { //around 112.5deg...
//		//	std::cout << "\n Label : " << i << " slope=" << angles[i];
//		//	color.val[0] = 255;
//		//	color.val[1] = 0;
//		//	color.val[2] = 150;
//		//	clusters[i] = 7; //violet --> 112.5 deg
//		//	cluster_img.at<Vec3b>(i_start + center[0], j_start + center[1]) = color;
//		//	circle(cluster_img, O, 1.5, color, -1, 8, 0);
//		//}
//		//else if (angles[i] > 120 && angles[i]  <= 145) { //around 135deg...
//		//	//std::cout << "\n slope=" << max_slope;
//		//	color.val[0] = 255;
//		//	color.val[1] = 0;
//		//	color.val[2] = 255;
//		//	clusters[i] = 4; //magenta --> 135 deg
//		//	cluster_img.at<Vec3b>(i_start + center[0], j_start + center[1]) = color;
//		//	circle(cluster_img, O, 1.5, color, -1, 8, 0);
//		//}
//		//else if (angles[i] > 145 && angles[i] <= 170) { //around 158.5deg...
//		//	//std::cout << "\n slope=" << max_slope;
//		//	color.val[0] = 0;
//		//	color.val[1] = 150;
//		//	color.val[2] = 255;
//		//	clusters[i] = 8; //orange --> 158.5 deg
//		//	cluster_img.at<Vec3b>(i_start + center[0], j_start + center[1]) = color;
//		//	circle(cluster_img, O, 1.5, color, -1, 8, 0);
//		//}
//		std::cout << "\nCluster " << clusters[i];
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				if (label[i_start + k][j_start + l] == i) {
//					bin_img.at<Vec3b>(i_start + k, j_start + l) = color;
//				}
//			}
//		}
//
//		//line(bin_img, P, Q, Scalar(0, 0, 0), 2.5, 8);
//
//		cluster_im[i_start + center[0]][j_start + center[1]] = 0;
//		Data[i][0] = j_start + center[1]; //centroid_x...
//		Data[i][1] = i_start + center[0]; //centroid_y...
//		Data[i][2] = angles[i]; //slope...
//		Data[i][3] = clusters[i]; //box orientation code...1:hor;2:ver;3:45deg;4:135deg...
//		Data[i][4] = sqrt((P.x - Q.x)*(P.x - Q.x) + (P.y - Q.y)*(P.y - Q.y)); //length of line...
//		Data[i][5] = 1;  //Label is included for re-clustering...
//
//		if (P.x < Q.x) {
//			ends[i][0] = P;
//			ends[i][1] = Q;
//		}
//		else if (P.x > Q.x) {
//			ends[i][0] = Q;
//			ends[i][1] = P;
//		}
//		else {
//			if (P.y < Q.y) {
//				ends[i][0] = P;
//				ends[i][1] = Q;
//			}
//			else {
//				ends[i][0] = Q;
//				ends[i][1] = P;
//			}
//		}
//
//		color.val[0] = 0;
//		color.val[1] = 255;
//		color.val[2] = 0;
//		//cc_img.at<Vec3b>(i_start + center[0], j_start + center[1]) = color;
//
//		//Delete ROI...
//		for (int k = 0; k < h; k++)
//			delete[] roi[k];
//		delete[] roi;
//	}
//
//	//	imshow("Clusters", cluster_img);
//	imshow("CentroidLocations", cc_img);
//	//	imwrite("Data/ICDAR/Data218/Clusters218.tif", cluster_img);
//	imwrite("Data/ICDAR/Data218/partitions218.tif", cc_img);
//	imwrite("Data/ICDAR/Data218/OldClusters218.tif", bin_img);
//
//	//----------------------Re-cluster based on nearness-------------------------
//
//	Mat lines;
//	lines.create(r, c, CV_8UC1);
//	new_cluster.create(r, c, CV_8UC1);
//	cvtColor(bin_img_org, new_cluster, CV_GRAY2BGR);
//	cvtColor(bin_img_org, lines, CV_GRAY2BGR);
//	int **change, change_count, prev_change_count = 0, **subimg;
//	change = new int *[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		change[i] = new int[3];
//		for (int j = 0; j < 3; j++)
//			change[i][j] = 0;
//	}
//
//	for (int i = 0; i < labelNum; i++) {
//		lineNo[i][2] = Data[i][0];
//		lineNo[i][3] = Data[i][1];
//	}
//
//	subimg = new int*[r];
//	for (int i = 0; i < r; i++) {
//		subimg[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			if (img[i][j] == 0)
//				subimg[i][j] = 187;
//			else
//				subimg[i][j] = 255;
//		}
//	}
//
//	//Remove very small and very large CCs...
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) //Not a significant component...
//			continue;
//
//		double ht, wd;
//		if (boundingBox[i][0] == 1) { //upright
//			ht = box[i][5] - box[i][4] + 1;
//			wd = box[i][3] - box[i][2] + 1;
//			if (ht <= wd) Data[i][3] = 1; //horizontal
//			else Data[i][3] = 2; //vertical
//		}
//		if (boundingBox[i][0] == 2) { //diagonal
//			ht = sqrt((boundingBox[i][5] - boundingBox[i][1])*(boundingBox[i][5] - boundingBox[i][1]) + (boundingBox[i][6] - boundingBox[i][2])*(boundingBox[i][6] - boundingBox[i][2]));
//			wd = sqrt((boundingBox[i][3] - boundingBox[i][5])*(boundingBox[i][3] - boundingBox[i][5]) + (boundingBox[i][4] - boundingBox[i][6])*(boundingBox[i][4] - boundingBox[i][6]));
//			if (ht <= wd) Data[i][3] = 4; //135deg
//			else Data[i][3] = 3; //45deg
//		}
//
//		//if (ht < 0.5*AH && wd < 0.5*AH) {
//		if (!(ht > 0.2 * AH && wd > 1.25*AH && ht <= 2.2 * AH) && !(wd > 0.2 * AH && ht > 1.25*AH && wd <= 2.2 * AH)) {
//
//		//if ((!(ht > 0.2 * AH && wd > 1.5*AH && ht <= 2 * AH) && !(wd > 0.2 * AH && ht > 1.5*AH && wd <= 2 * AH)) || abs(1 - ht / wd) < 0.2 || abs(1 - wd / ht) < 0.2) {
//			Data[i][5] = 0;
//			for (int k = box[i][4]; k <= box[i][5]; k++) {
//				for (int l = box[i][2]; l <= box[i][3]; l++) {
//					if (label[k][l] == i) {
//						subimg[k][l] = 255;
//					}
//				}
//			}
//			std::cout << "\nLabel " << i << " excluded!";
//			continue; //Exclude small pixels...
//		}
//		else {
//			for (int k = box[i][4]; k <= box[i][5]; k++) {
//				for (int l = box[i][2]; l <= box[i][3]; l++) {
//					if (label[k][l] == i) {
//						subimg[k][l] = 187;
//					}
//				}
//			}
//		}
//	}
//
//	sub_im.create(r, c, CV_8UC1);
//	arr2mat(subimg, sub_im);
//	cvtColor(sub_im, sub_im, CV_GRAY2BGR);
//	//Draw the bounding boxes...
//	//for (int i = 0; i < labelNum; i++) {
//	//	if (box[i][2] == -1 || Data[i][5] == 0) continue;
//	//	else {
//	//		if (boundingBox[i][0] == 2) { //Diagonal box....
//	//			line(sub_im, Point(boundingBox[i][2], boundingBox[i][1]), Point(boundingBox[i][6], boundingBox[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	//			line(sub_im, Point(boundingBox[i][6], boundingBox[i][5]), Point(boundingBox[i][4], boundingBox[i][3]), Scalar(0, 0, 255), 1.5, 8);
//	//			line(sub_im, Point(boundingBox[i][4], boundingBox[i][3]), Point(boundingBox[i][8], boundingBox[i][7]), Scalar(0, 0, 255), 1.5, 8);
//	//			line(sub_im, Point(boundingBox[i][8], boundingBox[i][7]), Point(boundingBox[i][2], boundingBox[i][1]), Scalar(0, 0, 255), 1.5, 8);
//	//		}
//	//		if (boundingBox[i][0] == 1) { //Upright box....
//	//			line(sub_im, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	//			line(sub_im, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	//			line(sub_im, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//	//			line(sub_im, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	//		}
//	//		//line(sub_im, ends[i][0], ends[i][1], Scalar(0, 255, 0), 2.5, 8);
//	//	}
//	//}
//
//	//Print Data...
//	std::cout << "\ncentroidx  centroidy   angle   color   span   big?";
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//		std::cout << "\nLabel " << i <<":" << Data[i][0] << "   " << Data[i][1] << "    " << Data[i][2] << "    " << Data[i][3] << "    " << Data[i][4] << "    " << Data[i][5];
//	}
//
//
//	//connectNbrs(img, label, box, boundingBox, Data, sub_im);
//	for (int it = 1; it <= 4; it++) {
//		HausdorffNbrs(subimg, label, box, boundingBox, Data, ends, sub_im);
//	}
//	for (int it = 1; it <= 4; it++) {
//		near2axis(subimg, label, box, boundingBox, Data, ends, sub_im);
//	}
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1 || Data[i][5] == 0) continue;
//		//arrowedLine(sub_im, ends[i][0], ends[i][1], Scalar(255, 0, 255), 2, 8, 0, 0.05);
//	}
//	connectHNbrs(img, label, box, boundingBox, Data, ends, connect, sub_im);
//	connectVNbrs(img, label, box, boundingBox, Data, ends, connect, sub_im);
//	connecDNbrs(img, label, box, boundingBox, Data, ends, connect, sub_im);
//	getLineNo(img, label, box, Data, connect, lineNo);
//	//Draw the connections...
//	for (int i = 1; i < labelNum; i++) {
//		if (connect[i][1] > 0) { 
//			Point c1 = Point(Data[i][0], Data[i][1]);
//			Point c2 = Point(Data[int(connect[i][1])][0], Data[int(connect[i][1])][1]);
//			line(scanimg, c1, c2, Scalar(0, 0, 255), 2, 8, 0);
//			if (lineNo[i][1] == 1) { //first word...
//				Point c1 = ends[i][0];
//				Point c2 = Point(Data[i][0], Data[i][1]);
//				line(scanimg, c1, c2, Scalar(0, 0, 255), 2, 8, 0);
//			}
//		}
//		else if (Data[i][5] == 1 && lineNo[i][1] == 3) { //last word...
//			Point c1 = Point(Data[i][0], Data[i][1]);
//			Point c2 = ends[i][1];
//			line(scanimg, c1, c2, Scalar(0, 0, 255), 2, 8, 0);
//		}
//		else if (Data[i][5] == 1 && lineNo[i][1] == 2) { //single word...
//			lineNo[i][0] = ++line_no;
//			line(scanimg, ends[i][0], ends[i][1], Scalar(0, 0, 255), 2, 8, 0);
//		}
//	}
//	imshow("MediumLabelswithPointers", scanimg);
//	imwrite("Data/ICDAR/Data218/MediumCC218.tif", scanimg);
//
//	//Re-cluster by proximity...
//	for (int it = 1; it <= 4; it++) {
//		for (int i = 0; i < labelNum; i++) {
//			for (int j = 0; j < 3; j++)
//				change[i][j] = 0;
//		}
//		
//		change_count = 0;
//		//cluster(img, Data, label, change);
//		//nbr8Label(subimg, label, box, Data, change, ends, it); 
//
//		if(it <= 4)
//			nbr8Label(subimg, label, box, boundingBox, Data, change, ends, it);
//		else
//			allNbrLabels(subimg, label, box, boundingBox, Data, change, it);
//		for (int i = 0; i < labelNum; i++) {
//			if (box[i][2] == -1) //Not a significant component...
//				continue;
//			Vec3b color;
//			int **roi, i_start = box[i][4], j_start = box[i][2];
//			int h = box[i][5] - box[i][4] + 1;
//			int w = box[i][3] - box[i][2] + 1;
//			roi = new int*[h];
//			for (int k = 0; k < h; k++) {
//				roi[k] = new int[w];
//				for (int l = 0; l < w; l++) {
//					if (label[i_start + k][j_start + l] == i) {
//						roi[k][l] = 0;
//					}
//					else
//						roi[k][l] = 255;
//				}
//			}
//
//			if (Data[i][5] == 0) {
//				color.val[0] = 150;
//				color.val[1] = 150;
//				color.val[2] = 150;
//			}
//			else {
//				int clr = int(Data[i][3]);
//				change_count += change[i][0];
//				//if (change[i][0] == 1) { //Change in clustering...
//				//	Data[i][3] = change[i][1];
//				if (clr == 1) { //0deg
//					color.val[0] = 0;
//					color.val[1] = 0;
//					color.val[2] = 255;
//				}
//				if (clr == 2) {  //45deg
//					color.val[0] = 255;
//					color.val[1] = 0;
//					color.val[2] = 0;
//				}
//				if (clr == 3) {  //90deg
//					color.val[0] = 0;
//					color.val[1] = 255;
//					color.val[2] = 0;
//				}
//				if (clr == 4) {  //135deg
//					color.val[0] = 255;
//					color.val[1] = 0;
//					color.val[2] = 255;
//				}
//				if (clr == 5) {  //22.5deg
//					color.val[0] = 255;
//					color.val[1] = 150;
//					color.val[2] = 0;
//				}
//				if (clr == 6) {  //67.5deg
//					color.val[0] = 50;
//					color.val[1] = 0;
//					color.val[2] = 100;
//				}
//				if (clr == 7) {  //112.5deg
//					color.val[0] = 255;
//					color.val[1] = 0;
//					color.val[2] = 150;
//				}
//				if (clr == 8) {  //158.5deg
//					color.val[0] = 0;
//					color.val[1] = 150;
//					color.val[2] = 255; //std::cout << i << "done$$$$";
//				}
//			}
//			
//			for (int k = 0; k < h; k++) {
//				for (int l = 0; l < w; l++) {
//					if (roi[k][l] == 0) {
//						new_cluster.at<Vec3b>(i_start + k, j_start + l) = color;
//					}
//				}
//			}
//
//			//line(cc_img, p7, p8, Scalar(255, 150, 0), 1.5, 8);
//
//			//Delete ROI...
//			for (int k = 0; k < h; k++)
//				delete[] roi[k];
//			delete[] roi;
//		}
//		std::cout << "\n iteration : " << it << " - No. of changes = " << change_count;
//		if (it == 4) {
//			imshow("NewClusters1", new_cluster);
//			imwrite("Data/ICDAR/Data218/NewClusters218(new1).tif", new_cluster);
//		}
//		if (it == 8) {
//			imshow("NewClusters2", new_cluster);
//			imwrite("Data/ICDAR/Data218/NewClusters218(new2).tif", new_cluster);
//		}
//		if (it > 15 && abs(prev_change_count - change_count) <= 1) break;
//		prev_change_count = change_count;
//	}
//	
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			Vec3b color, prev;
//			prev = lines.at<Vec3b>(i, j);
//			if (prev.val[0] == 255 && prev.val[1] == 255 && prev.val[2] == 255) continue;
//			if (Data[label[i][j]][5] == 0) {
//				color.val[0] = 225;
//				color.val[1] = 225;
//				color.val[2] = 225;
//			}
//			else {
//				color.val[0] = 187;
//				color.val[1] = 187;
//				color.val[2] = 187;
//			}
//			lines.at<Vec3b>(i, j) = color;
//		}
//	}
//
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1 || Data[i][5] == 0) //Not a significant component...
//			continue;
//		Point P, Q;
//		Scalar color;
//		double theta_m;
//		if (Data[i][3] != 3) {
//			if (Data[i][3] == 1) {
//				color = Scalar(0, 0, 255);
//				theta_m = 0;
//			}
//			if (Data[i][3] == 2) {
//				color = Scalar(255, 0, 0);
//				theta_m = tan(45 * CV_PI / 180);
//			}
//			if (Data[i][3] == 4) {
//				color = Scalar(255, 0, 255);
//				theta_m = tan(135 * CV_PI / 180);
//			}
//			if (Data[i][3] == 5) {
//				color = Scalar(255, 150, 0);
//				theta_m = tan(20 * CV_PI / 180);
//			}
//			if (Data[i][3] == 6) {
//				color = Scalar(50, 0, 100);
//				theta_m = tan(70 * CV_PI / 180);
//			}
//			if (Data[i][3] == 7) {
//				color = Scalar(255, 0, 150);
//				theta_m = tan(110 * CV_PI / 180);
//			}
//			if (Data[i][3] == 8) {
//				color = Scalar(0, 150, 255);
//				theta_m = tan(223 * CV_PI / 180);
//			}
//			
//			P.x = box[i][2];
//			P.y = Data[i][1] + tan(Data[i][2] * CV_PI / 180) * (P.x - Data[i][0]);
//			Q.x = box[i][3];
//			Q.y = Data[i][1] + tan(Data[i][2] * CV_PI / 180) * (Q.x - Data[i][0]);
//			line(lines, P, Q, color, 3, 8);
//		}
//		else {
//			color = Scalar(0, 255, 0);
//			P.x = Data[i][0];
//			P.y = box[i][4];
//			Q.x = Data[i][0];
//			Q.y = box[i][5];
//			line(lines, P, Q, color, 3, 8);
//		}
//	}
//	imshow("NewClusters", new_cluster);
//	imwrite("Data/ICDAR/Data218/NewClusters218(new).tif", new_cluster);
//	imshow("NewLines", lines);
//	imwrite("Data/ICDAR/Data218/Line218(new).tif", lines);
//
//	//GaussBlur(img, label, box, Data);
//
//	//----------------------------Delete memory---------------------------
//
//	for (int i = 0; i < r; i++) {
//		delete[] cluster_im[i];
//	}
//	delete[] cluster_im;
//
//	for (int i = 0; i < labelNum; i++) {
//		delete[] change[i];
//	}
//	delete[] change;
//	delete[] angles;
//	delete[] clusters;
//}
//
//void axisLine(Mat scanimg, int **img, Mat cc_img, Mat bin_img_org, int **label, double **boundingBox, double **diag_box, int **box, double **Data, int **lineNo, double **connect) {
//	Mat bin_img, img_org, cluster_img, new_cluster, sub_im;
//	int **cluster_im;
//	Point p1, p2, p3, p4, p5, p6, p7, p8, p, q, P, Q;
//	Point **ends;
//	int *angles, *clusters;
//	img_org.create(r, c, CV_8UC1);
//	cvtColor(bin_img_org, bin_img, CV_GRAY2BGR);
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			Vec3b color, prev;
//			prev = bin_img.at<Vec3b>(i, j);
//			if (prev.val[0] == 255 && prev.val[1] == 255 && prev.val[2] == 255) continue;
//			color.val[0] = 150;
//			color.val[1] = 150;
//			color.val[2] = 150;
//			bin_img.at<Vec3b>(i, j) = color;
//		}
//	}
//
//	ends = new Point*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		ends[i] = new Point[2];
//		for (int j = 0; j < 2; j++) {
//			ends[i][0] = Point(0, 0);
//			ends[i][1] = Point(0, 0);
//		}
//	}
//
//	clusters = new int[labelNum];
//	angles = new int[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		angles[i] = -1;
//		clusters[i] = 0;
//	}
//	
//	cluster_im = new int*[r];
//	for (int i = 0; i < r; i++) {
//		cluster_im[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			cluster_im[i][j] = 255;
//		}
//	}
//	cluster_img.create(r, c, CV_8UC1);
//	arr2mat(cluster_im, cluster_img);
//	cvtColor(cluster_img, cluster_img, CV_GRAY2BGR);
//
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) //Not a significant component...
//			continue;
//		int center[2];
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//
//		//--------------------------------Select ROI and determine centroid for black pixels---------------
//
//		int **roi, points[4];
//		int *X, *Y;
//		int count = 0, N = 0, i_start = box[i][4], j_start = box[i][2], max_line[2], pixels = 0;
//		double  slopeyx = 0, slopexy = 0, interceptyx = 0, interceptxy = 0, spreadYX, spreadXY, slope;
//		roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				if (label[i_start + k][j_start + l] == i) {
//					roi[k][l] = img[i_start + k][j_start + l];
//					N++;
//				}
//				else
//					roi[k][l] = 255;
//			}
//		}
//
//		X = new int[N];
//		Y = new int[N];
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				if (roi[k][l] == 0) {
//					X[count] = l;
//					Y[count] = k;
//					count++;
//				}
//			}
//		}
//		center[0] = (int)Mean(Y, N);
//		center[1] = (int)Mean(X, N);
//
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				if (label[i_start + k][j_start + l] == i) {
//					Vec3b color;
//					color.val[0] = 150;
//					color.val[1] = 150;
//					color.val[2] = 150;
//					cc_img.at<Vec3b>(i_start + k, j_start + l) = color;
//				}
//			}
//		}
//
//		//-----------------------------Draw the line with guess slope--------------------
//		//slope = (y-y0)/(x-x0)...
//		p.x = 0;
//		p.y = (center[0]) - Data[i][2] * center[1]; 
//		if (p.y < 0) {
//			p.y = 0;
//			p.x = (center[1]) - center[0] / Data[i][2];
//		}
//		else if (p.y >= h) {
//			p.y = h - 1;
//			p.x = (center[1]) + (p.y - center[0]) / Data[i][2];
//		}
//
//		q.x = w - 1;
//		q.y = (center[0]) + Data[i][2] * (q.x - center[1]);
//		if (q.y < 0) {
//			q.y = 0;
//			q.x = (center[1]) - center[0] / Data[i][2];
//		}
//		else if (q.y >= h) {
//			q.y = h - 1;
//			q.x = (center[1]) + (q.y - center[0]) / Data[i][2];
//		}
//
//		p.x += j_start;
//		p.y += i_start;
//		q.x += j_start;
//		q.y += i_start;
//
//		angles[i] = atan(Data[i][2]) * 180 / CV_PI;
//		std::cout << "\nLabel " << i << " has slope at " << angles[i] << " degrees";
//		line(cc_img, p, q, Scalar(255, 0, 150), 2, 8);
//
//		//----------------------------Rotate the best approximation by 45degrees on both sides for the exact one------------------
//
//		Vec3b color;
//		Point O;
//		O.x = j_start + center[1];
//		O.y = i_start + center[0];
//		P = p;
//		Q = q;
//		if (angles[i] >= 0 && angles[i] <= 45) {
//			int lambda = angles[i] / 9;
//			clusters[i] = 1;
//			color.val[0] = int(255 * lambda / 5);
//			color.val[1] = 0;
//			color.val[2] = 255;
//		}
//		if (angles[i] >= 50 && angles[i] <= 90) {
//			int lambda = (angles[i] - 45) / 9;
//			clusters[i] = 2;
//			color.val[0] = 255;
//			color.val[1] = 0;
//			color.val[2] = int(255 * (5 - lambda) / 5);
//		}
//		if (angles[i] >= 95 && angles[i] <= 135) {
//			int lambda = (angles[i] - 90) / 9;
//			clusters[i] = 3;
//			color.val[0] = int(255 * (5 - lambda) / 5);
//			color.val[1] = int(255 * lambda / 5);
//			color.val[2] = 0;
//		}
//		if (angles[i] >= 140 && angles[i] <= 180) {
//			int lambda = (angles[i] - 135) / 9;
//			clusters[i] = 4;
//			color.val[0] = 0;
//			color.val[1] = int(127 * (5 - lambda) / 5);
//			color.val[2] = int(255 * lambda / 5);
//		}
//
//		
//		std::cout << "\nCluster " << clusters[i];
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				if (label[i_start + k][j_start + l] == i) {
//					bin_img.at<Vec3b>(i_start + k, j_start + l) = color;
//				}
//			}
//		}
//
//		//line(bin_img, P, Q, Scalar(0, 0, 0), 2.5, 8);
//
//		cluster_im[i_start + center[0]][j_start + center[1]] = 0;
//		Data[i][0] = j_start + center[1]; //centroid_x...
//		Data[i][1] = i_start + center[0]; //centroid_y...
//		Data[i][2] = angles[i]; //slope...
//		Data[i][3] = clusters[i]; //box orientation code...1:hor;2:ver;3:45deg;4:135deg...
//		Data[i][4] = sqrt((P.x - Q.x)*(P.x - Q.x) + (P.y - Q.y)*(P.y - Q.y)); //length of line...
//		Data[i][5] = 1;  //Label is included for re-clustering...
//
//		if (P.x < Q.x) {
//			ends[i][0] = P;
//			ends[i][1] = Q;
//		}
//		else if (P.x > Q.x) {
//			ends[i][0] = Q;
//			ends[i][1] = P;
//		}
//		else {
//			if (P.y < Q.y) {
//				ends[i][0] = P;
//				ends[i][1] = Q;
//			}
//			else {
//				ends[i][0] = Q;
//				ends[i][1] = P;
//			}
//		}
//
//		color.val[0] = 0;
//		color.val[1] = 255;
//		color.val[2] = 0;
//		//cc_img.at<Vec3b>(i_start + center[0], j_start + center[1]) = color;
//
//		//Delete ROI...
//		for (int k = 0; k < h; k++)
//			delete[] roi[k];
//		delete[] roi;
//	}
//
//	//	imshow("Clusters", cluster_img);
//	imshow("CentroidLocations", cc_img);
//	//	imwrite("Data/ICDAR/Data218/Clusters218.tif", cluster_img);
//	imwrite("Data/ICDAR/Data218/partitions218.tif", cc_img);
//	imwrite("Data/ICDAR/Data218/OldClusters218.tif", bin_img);
//
//	//----------------------Re-cluster based on nearness-------------------------
//
//	Mat lines;
//	lines.create(r, c, CV_8UC1);
//	new_cluster.create(r, c, CV_8UC1);
//	cvtColor(bin_img_org, new_cluster, CV_GRAY2BGR);
//	cvtColor(bin_img_org, lines, CV_GRAY2BGR);
//	int **change, change_count, prev_change_count = 0, **subimg;
//	change = new int *[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		change[i] = new int[3];
//		for (int j = 0; j < 3; j++)
//			change[i][j] = 0;
//	}
//
//	for (int i = 0; i < labelNum; i++) {
//		lineNo[i][2] = Data[i][0];
//		lineNo[i][3] = Data[i][1];
//	}
//
//	subimg = new int*[r];
//	for (int i = 0; i < r; i++) {
//		subimg[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			if (img[i][j] == 0)
//				subimg[i][j] = 187;
//			else
//				subimg[i][j] = 255;
//		}
//	}
//
//	//Remove very small and very large CCs...
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) //Not a significant component...
//			continue;
//
//		double ht, wd;
//		if (boundingBox[i][0] == 1) { //upright
//			ht = box[i][5] - box[i][4] + 1;
//			wd = box[i][3] - box[i][2] + 1;
//			if (ht <= wd) Data[i][3] = 1; //horizontal
//			else Data[i][3] = 2; //vertical
//		}
//		if (boundingBox[i][0] == 2) { //diagonal
//			ht = sqrt((boundingBox[i][5] - boundingBox[i][1])*(boundingBox[i][5] - boundingBox[i][1]) + (boundingBox[i][6] - boundingBox[i][2])*(boundingBox[i][6] - boundingBox[i][2]));
//			wd = sqrt((boundingBox[i][3] - boundingBox[i][5])*(boundingBox[i][3] - boundingBox[i][5]) + (boundingBox[i][4] - boundingBox[i][6])*(boundingBox[i][4] - boundingBox[i][6]));
//			if (ht <= wd) Data[i][3] = 4; //135deg
//			else Data[i][3] = 3; //45deg
//		}
//
//		//if (ht < 0.5*AH && wd < 0.5*AH) {
//		if (!(ht > 0.2 * AH && wd > 1.25*AH && ht <= 2.2 * AH) && !(wd > 0.2 * AH && ht > 1.25*AH && wd <= 2.2 * AH)) {
//
//			//if ((!(ht > 0.2 * AH && wd > 1.5*AH && ht <= 2 * AH) && !(wd > 0.2 * AH && ht > 1.5*AH && wd <= 2 * AH)) || abs(1 - ht / wd) < 0.2 || abs(1 - wd / ht) < 0.2) {
//			Data[i][5] = 0;
//			for (int k = box[i][4]; k <= box[i][5]; k++) {
//				for (int l = box[i][2]; l <= box[i][3]; l++) {
//					if (label[k][l] == i) {
//						subimg[k][l] = 255;
//					}
//				}
//			}
//			std::cout << "\nLabel " << i << " excluded!";
//			continue; //Exclude small pixels...
//		}
//		else {
//			for (int k = box[i][4]; k <= box[i][5]; k++) {
//				for (int l = box[i][2]; l <= box[i][3]; l++) {
//					if (label[k][l] == i) {
//						subimg[k][l] = 187;
//					}
//				}
//			}
//		}
//	}
//
//	sub_im.create(r, c, CV_8UC1);
//	arr2mat(subimg, sub_im);
//	cvtColor(sub_im, sub_im, CV_GRAY2BGR);
//	//Draw the bounding boxes...
//	//for (int i = 0; i < labelNum; i++) {
//	//	if (box[i][2] == -1 || Data[i][5] == 0) continue;
//	//	else {
//	//		if (boundingBox[i][0] == 2) { //Diagonal box....
//	//			line(sub_im, Point(boundingBox[i][2], boundingBox[i][1]), Point(boundingBox[i][6], boundingBox[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	//			line(sub_im, Point(boundingBox[i][6], boundingBox[i][5]), Point(boundingBox[i][4], boundingBox[i][3]), Scalar(0, 0, 255), 1.5, 8);
//	//			line(sub_im, Point(boundingBox[i][4], boundingBox[i][3]), Point(boundingBox[i][8], boundingBox[i][7]), Scalar(0, 0, 255), 1.5, 8);
//	//			line(sub_im, Point(boundingBox[i][8], boundingBox[i][7]), Point(boundingBox[i][2], boundingBox[i][1]), Scalar(0, 0, 255), 1.5, 8);
//	//		}
//	//		if (boundingBox[i][0] == 1) { //Upright box....
//	//			line(sub_im, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	//			line(sub_im, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	//			line(sub_im, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//	//			line(sub_im, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	//		}
//	//		//line(sub_im, ends[i][0], ends[i][1], Scalar(0, 255, 0), 2.5, 8);
//	//	}
//	//}
//
//	//Print Data...
//	std::cout << "\ncentroidx  centroidy   angle   color   span   big?";
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//		std::cout << "\nLabel " << i << ":" << Data[i][0] << "   " << Data[i][1] << "    " << Data[i][2] << "    " << Data[i][3] << "    " << Data[i][4] << "    " << Data[i][5];
//	}
//
//
//	//connectNbrs(img, label, box, boundingBox, Data, sub_im);
//	for (int it = 1; it <= 4; it++) {
//		HausdorffNbrs(subimg, label, box, boundingBox, Data, ends, sub_im);
//	}
//	for (int it = 1; it <= 4; it++) {
//		near2axis(subimg, label, box, boundingBox, Data, ends, sub_im);
//	}
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1 || Data[i][5] == 0) continue;
//		//arrowedLine(sub_im, ends[i][0], ends[i][1], Scalar(255, 0, 255), 2, 8, 0, 0.05);
//	}
//	connectHNbrs(img, label, box, boundingBox, Data, ends, connect, sub_im);
//	connectVNbrs(img, label, box, boundingBox, Data, ends, connect, sub_im);
//	connecDNbrs(img, label, box, boundingBox, Data, ends, connect, sub_im);
//	getLineNo(img, label, box, Data, connect, lineNo);
//	//Draw the connections...
//	for (int i = 1; i < labelNum; i++) {
//		if (connect[i][1] > 0) {
//			Point c1 = Point(Data[i][0], Data[i][1]);
//			Point c2 = Point(Data[int(connect[i][1])][0], Data[int(connect[i][1])][1]);
//			line(scanimg, c1, c2, Scalar(0, 0, 255), 2, 8, 0);
//			if (lineNo[i][1] == 1) { //first word...
//				Point c1 = ends[i][0];
//				Point c2 = Point(Data[i][0], Data[i][1]);
//				line(scanimg, c1, c2, Scalar(0, 0, 255), 2, 8, 0);
//			}
//		}
//		else if (Data[i][5] == 1 && lineNo[i][1] == 3) { //last word...
//			Point c1 = Point(Data[i][0], Data[i][1]);
//			Point c2 = ends[i][1];
//			line(scanimg, c1, c2, Scalar(0, 0, 255), 2, 8, 0);
//		}
//		else if (Data[i][5] == 1 && lineNo[i][1] == 2) { //single word...
//			lineNo[i][0] = ++line_no;
//			line(scanimg, ends[i][0], ends[i][1], Scalar(0, 0, 255), 2, 8, 0);
//		}
//	}
//	imshow("MediumLabelswithPointers", scanimg);
//	imwrite("Data/ICDAR/Data218/MediumCC218.tif", scanimg);
//
//	//Re-cluster by proximity...
//	//for (int it = 1; it <= 4; it++) {
//	//	for (int i = 0; i < labelNum; i++) {
//	//		for (int j = 0; j < 3; j++)
//	//			change[i][j] = 0;
//	//	}
//
//	//	change_count = 0;
//	//	//cluster(img, Data, label, change);
//	//	//nbr8Label(subimg, label, box, Data, change, ends, it); 
//
//	//	if (it <= 4)
//	//		nbr8Label(subimg, label, box, boundingBox, Data, change, ends, it);
//	//	else
//	//		allNbrLabels(subimg, label, box, boundingBox, Data, change, it);
//	//	for (int i = 0; i < labelNum; i++) {
//	//		if (box[i][2] == -1) //Not a significant component...
//	//			continue;
//	//		Vec3b color;
//	//		int **roi, i_start = box[i][4], j_start = box[i][2];
//	//		int h = box[i][5] - box[i][4] + 1;
//	//		int w = box[i][3] - box[i][2] + 1;
//	//		roi = new int*[h];
//	//		for (int k = 0; k < h; k++) {
//	//			roi[k] = new int[w];
//	//			for (int l = 0; l < w; l++) {
//	//				if (label[i_start + k][j_start + l] == i) {
//	//					roi[k][l] = 0;
//	//				}
//	//				else
//	//					roi[k][l] = 255;
//	//			}
//	//		}
//
//	//		if (Data[i][5] == 0) {
//	//			color.val[0] = 150;
//	//			color.val[1] = 150;
//	//			color.val[2] = 150;
//	//		}
//	//		else {
//	//			int clr = int(Data[i][3]);
//	//			change_count += change[i][0];
//	//			//if (change[i][0] == 1) { //Change in clustering...
//	//			//	Data[i][3] = change[i][1];
//	//			if (clr == 1) { //0deg
//	//				color.val[0] = 0;
//	//				color.val[1] = 0;
//	//				color.val[2] = 255;
//	//			}
//	//			if (clr == 2) {  //45deg
//	//				color.val[0] = 255;
//	//				color.val[1] = 0;
//	//				color.val[2] = 0;
//	//			}
//	//			if (clr == 3) {  //90deg
//	//				color.val[0] = 0;
//	//				color.val[1] = 255;
//	//				color.val[2] = 0;
//	//			}
//	//			if (clr == 4) {  //135deg
//	//				color.val[0] = 255;
//	//				color.val[1] = 0;
//	//				color.val[2] = 255;
//	//			}
//	//			if (clr == 5) {  //22.5deg
//	//				color.val[0] = 255;
//	//				color.val[1] = 150;
//	//				color.val[2] = 0;
//	//			}
//	//			if (clr == 6) {  //67.5deg
//	//				color.val[0] = 50;
//	//				color.val[1] = 0;
//	//				color.val[2] = 100;
//	//			}
//	//			if (clr == 7) {  //112.5deg
//	//				color.val[0] = 255;
//	//				color.val[1] = 0;
//	//				color.val[2] = 150;
//	//			}
//	//			if (clr == 8) {  //158.5deg
//	//				color.val[0] = 0;
//	//				color.val[1] = 150;
//	//				color.val[2] = 255; //std::cout << i << "done$$$$";
//	//			}
//	//		}
//
//	//		for (int k = 0; k < h; k++) {
//	//			for (int l = 0; l < w; l++) {
//	//				if (roi[k][l] == 0) {
//	//					new_cluster.at<Vec3b>(i_start + k, j_start + l) = color;
//	//				}
//	//			}
//	//		}
//
//	//		//line(cc_img, p7, p8, Scalar(255, 150, 0), 1.5, 8);
//
//	//		//Delete ROI...
//	//		for (int k = 0; k < h; k++)
//	//			delete[] roi[k];
//	//		delete[] roi;
//	//	}
//	//	std::cout << "\n iteration : " << it << " - No. of changes = " << change_count;
//	//	if (it == 4) {
//	//		imshow("NewClusters1", new_cluster);
//	//		imwrite("Data/ICDAR/Data218/NewClusters218(new1).tif", new_cluster);
//	//	}
//	//	if (it == 8) {
//	//		imshow("NewClusters2", new_cluster);
//	//		imwrite("Data/ICDAR/Data218/NewClusters218(new2).tif", new_cluster);
//	//	}
//	//	if (it > 15 && abs(prev_change_count - change_count) <= 1) break;
//	//	prev_change_count = change_count;
//	//}
//
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			Vec3b color, prev;
//			prev = lines.at<Vec3b>(i, j);
//			if (prev.val[0] == 255 && prev.val[1] == 255 && prev.val[2] == 255) continue;
//			if (Data[label[i][j]][5] == 0) {
//				color.val[0] = 225;
//				color.val[1] = 225;
//				color.val[2] = 225;
//			}
//			else {
//				color.val[0] = 187;
//				color.val[1] = 187;
//				color.val[2] = 187;
//			}
//			lines.at<Vec3b>(i, j) = color;
//		}
//	}
//
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1 || Data[i][5] == 0) //Not a significant component...
//			continue;
//		Point P, Q;
//		Scalar color;
//		double theta_m;
//		if (Data[i][3] != 3) {
//			if (Data[i][3] == 1) {
//				color = Scalar(0, 0, 255);
//				theta_m = 0;
//			}
//			if (Data[i][3] == 2) {
//				color = Scalar(255, 0, 0);
//				theta_m = tan(45 * CV_PI / 180);
//			}
//			if (Data[i][3] == 4) {
//				color = Scalar(255, 0, 255);
//				theta_m = tan(135 * CV_PI / 180);
//			}
//			if (Data[i][3] == 5) {
//				color = Scalar(255, 150, 0);
//				theta_m = tan(20 * CV_PI / 180);
//			}
//			if (Data[i][3] == 6) {
//				color = Scalar(50, 0, 100);
//				theta_m = tan(70 * CV_PI / 180);
//			}
//			if (Data[i][3] == 7) {
//				color = Scalar(255, 0, 150);
//				theta_m = tan(110 * CV_PI / 180);
//			}
//			if (Data[i][3] == 8) {
//				color = Scalar(0, 150, 255);
//				theta_m = tan(223 * CV_PI / 180);
//			}
//
//			P.x = box[i][2];
//			P.y = Data[i][1] + tan(Data[i][2] * CV_PI / 180) * (P.x - Data[i][0]);
//			if (P.y < box[i][4]) {
//				P.y = box[i][4];
//				P.x = Data[i][0] + (P.y - Data[i][1]) / tan(Data[i][2] * CV_PI / 180);
//			}
//			else if (P.y > box[i][5]) {
//				P.y = box[i][5];
//				P.x = Data[i][0] + (P.y - Data[i][1]) / tan(Data[i][2] * CV_PI / 180);
//			}
//
//			Q.x = box[i][3];
//			Q.y = Data[i][1] + tan(Data[i][2] * CV_PI / 180) * (Q.x - Data[i][0]);
//			if (Q.y < box[i][4]) {
//				Q.y = box[i][4];
//				Q.x = Data[i][0] + (Q.y - Data[i][1]) / tan(Data[i][2] * CV_PI / 180);
//			}
//			else if (Q.y > box[i][5]) {
//				Q.y = box[i][5];
//				Q.x = Data[i][0] + (Q.y - Data[i][1]) / tan(Data[i][2] * CV_PI / 180);
//			}
//			line(lines, P, Q, color, 3, 8);
//		}
//		else {
//			color = Scalar(0, 255, 0);
//			P.x = Data[i][0];
//			P.y = box[i][4];
//			Q.x = Data[i][0];
//			Q.y = box[i][5];
//			line(lines, P, Q, color, 3, 8);
//		}
//	}
//	imshow("NewClusters", new_cluster);
//	imwrite("Data/ICDAR/Data218/NewClusters218(new).tif", new_cluster);
//	imshow("NewLines", lines);
//	imwrite("Data/ICDAR/Data218/Line218(new).tif", lines);
//
//	//GaussBlur(img, label, box, Data);
//
//	//----------------------------Delete memory---------------------------
//
//	for (int i = 0; i < r; i++) {
//		delete[] cluster_im[i];
//	}
//	delete[] cluster_im;
//
//	for (int i = 0; i < labelNum; i++) {
//		delete[] change[i];
//	}
//	delete[] change;
//	delete[] angles;
//	delete[] clusters;
//}
//
//
///*----------------------------------------------------------------*/
///*                 Invert, dilate and fill holes                  */
///*----------------------------------------------------------------*/
//
//void dilate(int **img, int h = r, int w = c) {
//	int n = 3;
//	int struc[3][3] = { {1,1,1},{1,1,1},{1,1,1} };
//	for (int i = 1; i < h - 1; i++) {
//		for (int j = 1; j < w - 1; j++) {
//			if (img[i][j] == 255) continue;
//			int sum = 0;
//			for (int k = 0; k < 3; k++) {
//				for (int l = 0; l < 3; l++) {
//					sum += img[i - 1 + k][j - 1 + l]*struc[k][l];
//				}
//			}
//			sum /= 255;
//			if (sum != 0) img[i][j] = 255;
//		}
//	}
//}
//
//void erode(int **img, int h = r, int w = c) {
//	int n = 3;
//	int **copy = new int*[r];
//	for (int i = 0; i < r; i++) {
//		copy[i] = new int[c];
//		for (int j = 0; j < c; j++)
//			copy[i][j] = img[i][j];
//	}
//	//int struc[7][7] = { { 0,0,0,0,0,0,0 },{ 1,1,1,1,1,1,1 } ,{ 1,1,1,1,1,1,1 } ,{ 1,1,1,1,1,1,1 },{ 1,1,1,1,1,1,1 },{ 1,1,1,1,1,1,1 },{ 0,0,0,0,0,0,0 } };
//	//int struc[5][5] = { {0,0,0,0,0}, { 1,1,1,1,1 },{ 1,1,1,1,1 },{ 1,1,1,1,1 },{ 0,0,0,0,0 } };
//	int struc[3][3] = { { 0,0,0 },{ 1,1,1 },{ 0,0,0 } };
//	for (int i = (n - 1) / 2; i < h - (n - 1) / 2; i++) {
//		for (int j = (n - 1) / 2; j < w - (n - 1) / 2; j++) {
//			if (copy[i][j] == 0) continue; //text pixel...
//			int sum = 0;
//			for (int k = -(n-1)/2; k <= (n - 1) / 2; k++) {
//				for (int l = -(n - 1) / 2; l <= (n - 1) / 2; l++) {
//					if (copy[i + k][j + l] == 0 && struc[(n - 1) / 2 + k][(n - 1) / 2 + l] == 1) sum++; //black nbr...
//				}
//			}
//			if (sum!=0) img[i][j] = 0;
//		}
//	}
//}
//
//void fill(int **img, int h = r, int w = c) {
//	Mat inv_img, fill_img;
//	int **inv_im;
//	inv_im = new int*[h];
//	for (int i = 0; i < h; i++) {
//		inv_im[i] = new int[w];
//		for (int j = 0; j < w; j++) {
//			inv_im[i][j] = 0;
//		}
//	}
//	invImage(img, inv_im, h, w);
//	//dilate(inv_im, h, w);
//	inv_img.create(h, w, CV_8UC1);
//	arr2mat(inv_im, inv_img, h, w);
//	//imshow("dilated", inv_img);
//	
//	//--------------------------Fill holes-------------------------
//	fill_img = inv_img.clone();
//
//	vector<vector<Point>> contours;
//
//	findContours(fill_img, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
//	drawContours(fill_img, contours, -1, Scalar(255), CV_FILLED);
//
//	mat2arr(fill_img, inv_im, h, w);
//	//imshow("filled", fill_img);
//	invImage(inv_im, img, h, w);
//}
//
///*----------------------------------------------------------------*/
///*        Find the diagonal boundary boxes for each CC            */
///*----------------------------------------------------------------*/
//void diagonalCC(int **label, int **mainCC, int **box, double **boundingBox, double **diag_box, Mat cc_img, int *trimlabel) {
//	double h1, w1, area, area_diag, sd0, sd45, sd135, sd90;
//	int count = 0, k = 0, diag_count = 0, rectangle_count = 0;
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1)
//			continue;
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//
//		int NE = 0, SE = 0, SW = 0, NW = 0;
//		int **roi, i_start = box[i][4], j_start = box[i][2];
//		roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				if (label[i_start + k][j_start + l] == i)
//					roi[k][l] = mainCC[i_start + k][j_start + l];
//				else
//					roi[k][l] = 255;
//				if (roi[k][l] == 0) count++;
//
//			}
//		}
//		std::cout << "\nLabel " << i << " : ";
//
//		int whitetotal = 0;
//		int *histH = new int[h];
//		for (int k = 0; k < h; k++)
//			histH[k] = 0;
//		sd0 = project0(roi, h, w, histH);
//		/*for (int k = 0; k < h; k++)
//			whitetotal += (w - histH[k]);
//		int count0 = deletePixels(histH, h, w);
//		double white0 = count0 / (double(whitetotal));*/
//
//		whitetotal = 0;
//		int *histV = new int[w];
//		for (int k = 0; k < w; k++)
//			histV[k] = 0;
//		sd90 = project90(roi, h, w, histV);
//		/*for (int k = 0; k < w; k++)
//			whitetotal += (h - histV[k]);
//		int count90 = deletePixels(histV, w, h);
//		double white90 = count90 / (double(whitetotal));*/
//
//		diagonalBox(roi, h, w, NW, SE, NE, SW);
//		diag_box[i][0] = box[i][4] + (NW - NE) / 2;  // N.y
//		diag_box[i][1] = box[i][2] + (NW + NE) / 2;  // N.x
//		diag_box[i][2] = box[i][4] + (SE - SW) / 2;  // S.y
//		diag_box[i][3] = box[i][2] + (SE + SW) / 2;  // S.x
//		diag_box[i][4] = box[i][4] + (SE - NE) / 2;  // E.y
//		diag_box[i][5] = box[i][2] + (SE + NE) / 2;  // E.x
//		diag_box[i][6] = box[i][4] + (NW - SW) / 2;  // W.y
//		diag_box[i][7] = box[i][2] + (NW + SW) / 2;  // W.x
//
//		
//
//		whitetotal = 0;
//		int szOD = 2*(diag_box[i][6] - diag_box[i][0] + 1);
//		int *histOD = new int[szOD]; //W.y - N.y + 1
//		for (int k = 0; k < szOD; k++)
//			histOD[k] = 0;
//		sd45 = project45(mainCC, label, diag_box[i], i, histOD);
//		/*for (int k = 0; k < szOD; k++)
//			whitetotal += ((diag_box[i][2] - diag_box[i][6] + 1) - histOD[k]);
//		int count45 = deletePixels(histOD, szOD, diag_box[i][2] - diag_box[i][6] + 1);
//		double white45 = count45 / (double(whitetotal));
//*/
//		whitetotal = 0;
//		int szD = 2*(diag_box[i][2] - diag_box[i][6] + 1);
//		int *histD = new int[szD]; //S.y - W.y + 1
//		for (int k = 0; k < szD; k++)
//			histD[k] = 0;
//		sd135 = project135(mainCC, label, diag_box[i], i, histD);
//		/*for (int k = 0; k < szD; k++)
//			whitetotal += ((diag_box[i][6] - diag_box[i][0] + 1) - histD[k]);
//		int count135 = deletePixels(histD, szD, diag_box[i][6] - diag_box[i][0] + 1);
//		double white135 = count135 / (double(whitetotal));
//
//		std::cout << "\nNo.of pixels deleted:" << count0 << "(0) , " << count90 << "(90) , " << count45 << "(45) , " << count135 << "(135)";
//		std::cout << "\nNo.of white space deleted:" << white0 << "(0) , " << white90 << "(90) , " << white45 << "(45) , " << white135 << "(135)";*/
//
//		area = ((double)h) * w;
//		h1 = sqrt((diag_box[i][0] - diag_box[i][4])*(diag_box[i][0] - diag_box[i][4]) + (diag_box[i][1] - diag_box[i][5])*(diag_box[i][1] - diag_box[i][5]));
//		w1 = sqrt((diag_box[i][2] - diag_box[i][4])*(diag_box[i][2] - diag_box[i][4]) + (diag_box[i][3] - diag_box[i][5])*(diag_box[i][3] - diag_box[i][5]));
//		area_diag = ( h1 * w1);
//		//std::cout << "\n label = " << i << "area = " << area << "  diagonal area = " << area_diag;
//
//		int length = h > w ? h : w;
//		int diag_length = h1 > w1 ? h1 : w1;
//
//		//minimum count of text pixels deleted in 0, 45, 90, 135 degrees....
//		/*int min = count0;
//		int mindir = 0;
//
//		if (count45 < min) {
//			min = count45;
//			mindir = 45;
//		}
//		if (count135 < min) {
//			min = count135;
//			mindir = 135;
//		}
//		if (count90 < min) {
//			min = count90;
//			mindir = 90;
//		}
//
//		int min = white0;
//		int mindir = 0;
//
//		if (white45 > min) {
//			min = white45;
//			mindir = 45;
//		}
//		if (white135 > min) {
//			min = white135;
//			mindir = 135;
//		}
//		if (white90 > min) {
//			min = white90;
//			mindir = 90;
//		}
//*/
//		//maximum std deviation of projection profiles in 0, 45, 90, 135 degrees....
//		double max = sd0;
//		int maxdir = 0;
//
//		if (sd45 > max) {
//			max = sd45;
//			maxdir = 45;
//		}
//		if (sd135 > max) {
//			max = sd135;
//			maxdir = 135;
//		}
//		if (sd90 > max) {
//			max = sd90;
//			maxdir = 90;
//		}
//		
//		double *sd = new double[4];
//		sd[0] = sd0;
//		sd[1] = sd90;
//		sd[2] = sd45;
//		sd[3] = sd135;
//		double std = stddev(sd, 4);
//		std::cout << "\nstd deviation = " << std;
//
//		//Decide direction of trimming the cc...
//		//if (min == 0) {
//			if (maxdir == 0)
//				trimlabel[i] = 1;
//			if (maxdir == 90)
//				trimlabel[i] = 2;
//			if (maxdir == 45)
//				trimlabel[i] = 3;
//			if (maxdir == 135)
//				trimlabel[i] = 4;
//		//}
//		/*else {
//			if (mindir == 0)
//				trimlabel[i] = 1;
//			if (mindir == 90)
//				trimlabel[i] = 2;
//			if (mindir == 45)
//				trimlabel[i] = 3;
//			if (mindir == 135)
//				trimlabel[i] = 4;
//		}
//*/
//		//Decide the boundary to display...
//		if ((std < 1.1 && area > area_diag) || (std >= 1.1 && (maxdir == 45 || maxdir == 135))) { //Diagonal box is more filling than upright box...
//			boundingBox[i][0] = 2;
//			boundingBox[i][1] = diag_box[i][0];
//			boundingBox[i][2] = diag_box[i][1];
//			boundingBox[i][3] = diag_box[i][2];
//			boundingBox[i][4] = diag_box[i][3];
//			boundingBox[i][5] = diag_box[i][4];
//			boundingBox[i][6] = diag_box[i][5];
//			boundingBox[i][7] = diag_box[i][6];
//			boundingBox[i][8] = diag_box[i][7];
//			diag_count++;
//			line(cc_img, Point(diag_box[i][1], diag_box[i][0]), Point(diag_box[i][5], diag_box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(diag_box[i][5], diag_box[i][4]), Point(diag_box[i][3], diag_box[i][2]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(diag_box[i][3], diag_box[i][2]), Point(diag_box[i][7], diag_box[i][6]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(diag_box[i][7], diag_box[i][6]), Point(diag_box[i][1], diag_box[i][0]), Scalar(0, 0, 255), 1.5, 8);
//		}
//		else if ((std < 1.1 && area < area_diag) || (std >= 1.0 && (maxdir == 0 || maxdir == 90))) {
//			boundingBox[i][0] = 1;
//			boundingBox[i][1] = box[i][2]; // Left x
//			boundingBox[i][2] = box[i][3]; // Right x
//			boundingBox[i][3] = box[i][4]; // Top y
//			boundingBox[i][4] = box[i][5]; // Bottom y
//			rectangle_count++;
//			line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//		}
//
//		for (int k = 0; k < h; k++) {
//			delete[] roi[k];
//		}
//		delete[] roi;
//		delete histD;
//		delete histH;
//		delete histV;
//		delete histOD;
//	}
//	std::cout << "\nNo. of diagonal boxes = " << diag_count;
//	std::cout << "\nNo. of upright rectangles = " << rectangle_count;
//	imshow("Components", cc_img);
//	//imwrite("Data/ICDAR/Data55/DiagBox218.tif", cc_img);
//}
//
///*----------------------------------------------------------------*/
///*           Read masked image and find the main CC's             */
///*----------------------------------------------------------------*/
//void MaskCC(Mat img, int **rawimage) {
//	Mat orgimg;
//	orgimg.create(r, c, CV_8UC1);
//	cvtColor(img, orgimg, CV_GRAY2BGR);
//
//	Mat img1, bin_img, cc_img, main_img, thin_img, lines_img;
//	int **rawimg, **image, **bin_im1, **label, **box, **mainCC, *smallCC, **thin_im, *heights, **lineNo, *trimlabel;
//	double h1, w1, area, area_diag, **diag_box, **boundingBox, *diagheights;
//	int count = 0, k = 0, diag_count = 0, rectangle_count = 0;
//
//	//----------------------Blur & Binarize image-----------------------
//	image = new int*[r];
//	rawimg = new int*[r];
//	bin_im1 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		image[i] = new int[c];
//		rawimg[i] = new int[c];
//		bin_im1[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			image[i][j] = rawimage[i][j];
//			rawimg[i][j] = rawimage[i][j];
//			bin_im1[i][j] = rawimage[i][j];
//		}
//	}
//	/*mat2arr(img, rawimg);
//	mat2arr(img, image);*/
//	/*blur(image);
//	img1.create(r, c, CV_8UC1);
//	arr2mat(image, img1);
//	imshow("Blurred", img1);*/
//	//double T = Otsu(image, img.rows, img.cols, bin_im1);
//	//NICK(image, r, c, bin_im1);
//	bin_img.create(r, c, CV_8UC1);
//	bin_img = img;
//	//arr2mat(bin_im1, bin_img);
//	cvtColor(bin_img, lines_img, CV_GRAY2BGR);
//	erode(bin_im1);
//	arr2mat(bin_im1, bin_img);
//	imwrite("Data/ICDAR/Data218/OriginalBinary1.tif", bin_img);
//	/*GaussianBlur(bin_img, bin_img, Size(21, 21), 0, 0);
//	mat2arr(bin_img, image);
//	NICK(image, r, c, bin_im1);
//	arr2mat(bin_im1, bin_img);
//	imwrite("Data/ICDAR/Data218/OriginalBinary.tif", bin_img);*/
//
//	//--------------------------Merge pixels with 1-pixel gap--------------------------
//	/*for (int i = 1; i < r - 1; i++) {
//		for (int j = 1; j < c - 1; j++) {
//			if (bin_im1[i][j] == 0)
//				continue;
//			int *p = new int[9];
//			p[0] = (255 - bin_im1[i][j]) / 255;
//			p[1] = (255 - bin_im1[i - 1][j - 1]) / 255;
//			p[2] = (255 - bin_im1[i - 1][j]) / 255;
//			p[3] = (255 - bin_im1[i - 1][j + 1]) / 255;
//			p[4] = (255 - bin_im1[i][j + 1]) / 255;
//			p[5] = (255 - bin_im1[i + 1][j + 1]) / 255;
//			p[6] = (255 - bin_im1[i + 1][j]) / 255;
//			p[7] = (255 - bin_im1[i + 1][j - 1]) / 255;
//			p[8] = (255 - bin_im1[i][j - 1]) / 255;
//
//			if (p[0] * p[2] * p[6] == 0 && p[1] + p[8] + p[7] > 1 && p[3] + p[4] + p[5] > 1) {
//				bin_im1[i][j] = 0;
//				std::cout << "\nchanged";
//			}
//			delete[] p;
//		}
//	}
//	arr2mat(bin_im1, bin_img);
//	imshow("MergedBinary.jpg", bin_img);
//*/
//	//--------------------------Connected components--------------------------
//	label = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label[i] = new int[c];
//		for (int j = 0; j < c; j++)
//			label[i][j] = 0;
//	}
//	labelNum = labelling(bin_im1, label, r, c);
//	//set_labelNum(labelling(bin_im1, label, r, c));
//	//labelNum = labelNum + 1;
//
//	//------------------Bound the components by rectangular boxes------------------
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	cvtColor(bin_img, cc_img, CV_GRAY2BGR);
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1) { //component not available...
//			count++;
//			continue;
//		}
//		std::cout << i << ":" << box[i][2] << "," << box[i][3] << "," << box[i][4] << "," << box[i][5] << "\n";
//		line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//		line(cc_img, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//		line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//		line(cc_img, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	}
//
//	imshow("Components", cc_img);
//	//imwrite("Data/ICDAR/Data25/MaskCC25.jpg", cc_img);
//
//	//-------------------Average height of components--------------------
//	int length = labelNum - count;
//	heights = new int[length];
//	for (int i = 0; i < length; i++)
//		heights[i] = 0;
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1)
//			continue;
//		heights[k++] = box[i][5] - box[i][4] + 1;
//		std::cout << "label " << i << " height " << heights[k - 1] << "\n";
//	}
//	AH = averageHeight(heights, length);
//	//set_AH(averageHeight(heights, length));
//	//AH = 35;
//	std::cout << "Average height = " << AH;
//
//	//---------------Delete CC's with small height and length------------
//	smallCC = new int[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) {
//			smallCC[i] = -1;
//			continue;
//		}
//		/*if ((box[i][5] - box[i][4]) > 0.65*AH || (box[i][3] - box[i][2]) > 0.65*AH) {
//			smallCC[i] = 1;
//		}*/
//		if ((box[i][5] - box[i][4]) < 0.5*AH || (box[i][3] - box[i][2]) < 0.5*AH) {
//			smallCC[i] = 0;
//			box[i][2] = -1;
//			box[i][3] = -1;
//			box[i][4] = -1;
//			box[i][5] = -1;
//			//std::cout << "\nLabel " << i << " deleted!";
//		}
//		else smallCC[i] = 1;
//		/*else {
//			smallCC[i] = 0;
//			box[i][2] = -1;
//			box[i][3] = -1;
//			box[i][4] = -1;
//			box[i][5] = -1;
//		}*/
//			
//	}
//
//	mainCC = new int *[r];
//	for (int i = 0; i < r; i++) {
//		mainCC[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			if (smallCC[label[i][j]] <= 0) {
//				label[i][j] = 0;
//				mainCC[i][j] = 255;
//			}
//			else
//				mainCC[i][j] = bin_im1[i][j];
//		}
//	}
//	main_img.create(r, c, CV_8UC1);
//	arr2mat(mainCC, main_img);
//	imshow("Deleting the small CC's", main_img);
//	imwrite("Data/ICDAR/Data218/Nick_bin218.tif", main_img);
//
//	//---------------Diagonal CC box---------------
//	trimlabel = new int[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		trimlabel[0] = 0;
//
//	diag_box = new double*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		diag_box[i] = new double[8];
//		for (int j = 0; j < 8; j++)
//			diag_box[i][j] = 0;
//	}
//	boundingBox = new double*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		boundingBox[i] = new double[9];
//		for (int j = 0; j < 9; j++) {
//			boundingBox[i][j] = -1;
//		}
//	}
//	cvtColor(main_img, cc_img, CV_GRAY2BGR);
//	diagonalCC(label, mainCC, box, boundingBox, diag_box, cc_img, trimlabel);
//
//	//---------------Delete diagonal CC's with small height and length------------
//	smallCC = new int[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1 || boundingBox[i][0] == 1) {
//			continue;
//		}
//		double h1 = sqrt((boundingBox[i][5] - boundingBox[i][1])*(boundingBox[i][5] - boundingBox[i][1]) + (boundingBox[i][6] - boundingBox[i][2])*(boundingBox[i][6] - boundingBox[i][2]));
//		double w1 = sqrt((boundingBox[i][5] - boundingBox[i][3])*(boundingBox[i][5] - boundingBox[i][3]) + (boundingBox[i][6] - boundingBox[i][4])*(boundingBox[i][6] - boundingBox[i][4]));
//		
//		//if (h1 > 0.65*AH || w1 > 0.65*AH) {
//		//	smallCC[i] = 1;
//		//	//std::cout << "\nLabel " << i << " deleted!";
//		//}
//		if (h1 < 0.5*AH || w1 < 0.5*AH) {
//			smallCC[i] = 0;
//			boundingBox[i][0] = 0;
//			boundingBox[i][1] = 0;
//			boundingBox[i][2] = 0;
//			boundingBox[i][3] = 0;
//			boundingBox[i][4] = 0;
//			boundingBox[i][5] = 0;
//			boundingBox[i][6] = 0;
//			boundingBox[i][7] = 0;
//			boundingBox[i][8] = 0;
//			box[i][2] = -1;
//			box[i][3] = -1;
//			box[i][4] = -1;
//			box[i][5] = -1;
//			//std::cout << "\nLabel " << i << " deleted!";
//		}
//		else smallCC[i] = 1;
//		//else {
//		//	smallCC[i] = 0;
//		//	boundingBox[i][0] = 0;
//		//	boundingBox[i][1] = 0;
//		//	boundingBox[i][2] = 0;
//		//	boundingBox[i][3] = 0;
//		//	boundingBox[i][4] = 0;
//		//	boundingBox[i][5] = 0;
//		//	boundingBox[i][6] = 0;
//		//	boundingBox[i][7] = 0;
//		//	boundingBox[i][8] = 0;
//		//	box[i][2] = -1;
//		//	box[i][3] = -1;
//		//	box[i][4] = -1;
//		//	box[i][5] = -1;
//		//	//std::cout << "\nLabel " << i << " deleted!";
//		//}
//	}
//
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (smallCC[label[i][j]] == 0) {
//				label[i][j] = 0;
//				mainCC[i][j] = 255;
//			}
//		}
//	}
//
//	//----------------Dilate and fill up holes in each CC-----------------
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//		int i_start = box[i][4], j_start = box[i][2];
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//		int **roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				if (label[i_start + k][j_start + l] == i) {
//					roi[k][l] = mainCC[i_start + k][j_start + l];
//				}
//				else
//					roi[k][l] = 255;
//			}
//
//		}
//
//		fill(roi, h, w);
//
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				if (roi[k][l] == 0 ) {
//					label[i_start + k][j_start + l] = i;
//					mainCC[i_start + k][j_start + l] = roi[k][l] ;
//				}
//			}
//		}
//
//		//Delete roi...
//		for (int k = 0; k < h; k++)
//			delete[] roi[k];
//		delete[] roi;
//	}
//
//	arr2mat(mainCC, main_img);
//	imshow("Deleting the small CC's(2)", main_img);
//	imwrite("Data/ICDAR/Data218/DeletedsmallCC.tif", main_img);
//
//	//-------------------Average height of components--------------------
//	count = 0;
//	double **Data = new double*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		Data[i] = new double[6];
//		for (int j = 0; j < 6; j++) {
//			Data[i][j] = 0;
//		}
//	}
//	cvtColor(main_img, cc_img, CV_GRAY2BGR);
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1 || boundingBox[i][0] == -1) { //component not available...
//			count++;
//			continue;
//		}
//		Point ends[4][2];
//		Centroidline(mainCC, label, i, diag_box[i], box[i], ends);
//
//		double lengths[4], maxlength = 0;
//		int dir = -1;
//		for (int k = 0; k < 4; k++) {
//			lengths[k] = norm(ends[k][0] - ends[k][1]);
//			if (lengths[k] > maxlength) {
//				maxlength = lengths[k];
//				dir = k;
//			}
//		}
//		if (dir == -1) {
//			box[i][2] = -1;
//			continue;
//		}
//		if (ends[dir][0].y == ends[dir][1].y)
//			Data[i][2] = 0; //horizontal
//		else if (ends[dir][0].x == ends[dir][1].x)
//			Data[i][2] = 9999; //vertical
//		else
//			Data[i][2] = (ends[dir][1].y - ends[dir][0].y) / double(ends[dir][1].x - ends[dir][0].x);
//		std::cout << "\nLabel " << i << " has angle " << Data[i][2];
//
//		//if (dir == 2 || dir == 3) {
//		if (boundingBox[i][0] == 2) { //Diagonal box....
//			line(cc_img, Point(diag_box[i][1], diag_box[i][0]), Point(diag_box[i][5], diag_box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(diag_box[i][5], diag_box[i][4]), Point(diag_box[i][3], diag_box[i][2]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(diag_box[i][3], diag_box[i][2]), Point(diag_box[i][7], diag_box[i][6]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(diag_box[i][7], diag_box[i][6]), Point(diag_box[i][1], diag_box[i][0]), Scalar(0, 0, 255), 1.5, 8);
//
//			//line(cc_img, ends[dir][0], ends[dir][1], Scalar(0, 120, 255), 2, 8);
//		}
//		//if (dir == 0 || dir == 1) {
//		if (boundingBox[i][0] == 1) { //Upright box....
//			line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//
//			//line(cc_img, ends[dir][0], ends[dir][1], Scalar(255, 0, 255), 2, 8);
//		}
//		if (dir == 0 || dir == 1)
//			line(cc_img, ends[dir][0], ends[dir][1], Scalar(255, 0, 255), 1.5, 8);
//		if (dir == 2 || dir == 3)
//			line(cc_img, ends[dir][0], ends[dir][1], Scalar(0, 120, 255), 1.5, 8);
//		
//		/*line(cc_img, ends[1][0], ends[1][1], Scalar(120, 255, 120), 1.5, 8);
//		line(cc_img, ends[2][0], ends[2][1], Scalar(255, 120, 120), 1.5, 8);
//		line(cc_img, ends[3][0], ends[3][1], Scalar(120, 120, 255), 1.5, 8);*/
//	}
//	imshow("Components", cc_img);
//	imwrite("Data/ICDAR/Data218/DiagBox218.tif", cc_img);
//
//	length = labelNum - count;
//	diagheights = new double[length];
//	for (int i = 0; i < length; i++)
//		diagheights[i] = 0;
//	k = 0;
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1 || boundingBox[i][0] == -1)
//			continue;
//		if(boundingBox[i][0] == 1) //upright
//			diagheights[k++] = box[i][5] - box[i][4] + 1;
//		if (boundingBox[i][0] == 2) //diagonal
//			diagheights[k++] = sqrt((boundingBox[i][5] - boundingBox[i][1])*(boundingBox[i][5] - boundingBox[i][1]) + (boundingBox[i][6] - boundingBox[i][2])*(boundingBox[i][6] - boundingBox[i][2]));
//		std::cout << "\nlabel " << i << " height " << diagheights[k - 1] << "\n";
//	}
//	AH = averageDiagHeight(diagheights, length);
//	//set_AH(averageHeight(heights, length));
//	std::cout << "Average height = " << AH;
//
//	//-----------------Find major axis for each CC and connect---------------------------
//	//skewLine(mainCC, cc_img, label, boundingBox, box);
//	lineNo = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		lineNo[i] = new int[4];
//		lineNo[i][0] = 0; //line no. for i-th label...
//		lineNo[i][1] = 0; //1 if i-th label is the starting of the line, 2 if i-th label is a single word...
//		lineNo[i][2] = 0; //centroid.x...
//		lineNo[i][3] = 0;//centroidy...
//	}
//
//	double **connect = new double*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		connect[i] = new double[3];
//		connect[i][0] = i; //label no....
//		connect[i][1] = -1;  //next label...
//		connect[i][2] = 0;  //distance...
//	}
//
//
//
//	//partition(orgimg, mainCC, cc_img, main_img, label, boundingBox, diag_box, box, Data, lineNo, connect);
//
//	axisLine(orgimg, mainCC, cc_img, main_img, label, boundingBox, diag_box, box, Data, lineNo, connect);
//
//
//	int colors[10][3] = { { 0,0,255 },{ 0,255,0 },{ 255,0,0 },{ 128,0,255 },{ 128,128,0 },{ 255,128,0 },{ 255,0,128 },{ 128,255,0 },{ 97,42,190 },{ 144,41,78 } };
//
//	std::cout << "\nTotal no. of text lines = " << line_no;
//
//
//	//---------------------assign line no. to each text pixel----------------
//	//collect the starting CC and the inclination of each line in the document...
//	double **lineData = new double*[line_no];
//	for (int i = 0; i < line_no; i++) {
//		lineData[i] = new double[3];
//		lineData[i][0] = 0; //starting label...
//		lineData[i][1] = 0; //slope of the line...
//		lineData[i][2] = 0; //no. of cc...
//	}
//
//	for (int i = 0; i < labelNum; i++) {
//		if (lineNo[i][1] == 2) { //single word...
//			int lno = lineNo[i][0] - 1;
//			lineData[lno][0] = i;
//			lineData[lno][2] = 1;
//			continue;
//		}
//
//		if (lineNo[i][1] != 1) continue; //not the starting of line...
//
//		int lno = lineNo[i][0] - 1;
//		lineData[lno][0] = i;
//		lineData[lno][2]++;
//
//		double X = lineNo[i][2], Y = lineNo[i][3], XY = X*Y, X2 = X*X;
//		int curr = i, next = connect[i][1], count = 1;
//		while (next > 0) {
//			curr = next;
//			X += lineNo[curr][2];
//			Y += lineNo[curr][3];
//			XY += (lineNo[curr][2] * lineNo[curr][3]);
//			X2 += (lineNo[curr][2] * lineNo[curr][2]);
//			lineData[lno][2]++;
//			count++;
//			next = connect[curr][1];
//		}
//		X += lineNo[curr][2];
//		Y += lineNo[curr][3];
//		XY += (lineNo[curr][2] * lineNo[curr][3]);
//		X2 += (lineNo[curr][2] * lineNo[curr][2]);
//		lineData[lno][2]++;
//		count++;
//
//		lineData[lno][1] = (count*XY - X*Y) / (count*X2 - X*X); //regression slope...
//		/*Point p1 = Point(lineNo[i][2], lineNo[i][3]);
//		Point p2 = Point(lineNo[curr][2], lineNo[i][3]+ lineData[lno][1]*(lineNo[curr][2]- lineNo[i][2]));
//		line(lines_img, p1, p2, Scalar(255, 0, 255), 2, 8, 0);*/
//
//	}
//
//	std::cout << "\nLine data::";
//	for (int i = 0; i < line_no; i++) {
//		std::cout << "\nLine " << i + 1;
//		std::cout << "\nStarts at : label " << lineData[i][0];
//		std::cout << "\nSlope =  " << lineData[i][1];
//		std::cout << "\nNo. of words = " << lineData[i][2];
//	}
//
//	//re-assign single words...
//	for (int i = 0; i < line_no; i++) {
//		if (lineData[i][2] != 1) continue;
//		double *linedist = new double[line_no];
//		double mindist = 99999;
//		int min = 0;
//		Point pt = Point(lineNo[int(lineData[i][0])][2], lineNo[int(lineData[i][0])][3]);
//		for (int k = 0; k < line_no; k++) {
//			if (lineData[k][2] < 2) {
//				linedist[k] = -1;
//				continue;
//			}
//			int startlbl = lineData[k][0];
//			Point start = Point(lineNo[startlbl][2], lineNo[startlbl][3]);
//			double slope = lineData[k][1];
//			linedist[k] = abs(slope*(pt.x - start.x) - (pt.y - start.y)) / sqrt(1 + slope*slope);
//			if (linedist[k] < mindist && linedist[k] >= 0) {
//				mindist = linedist[k];
//				min = k;
//			}
//		}
//		if (mindist < AH) {
//			lineNo[int(lineData[i][0])][0] = min + 1;
//			lineData[min][2]++;
//			lineData[i][2]--;
//			std::cout << "\nLabel " << lineData[i][0] << " added to line #" << min + 1 << " at distance " << mindist;
//		}
//		delete[] linedist;
//	}
//
//	//merge lines with same slope and intercept...
//	for (int i = 0; i < line_no; i++) {
//		if (lineData[i][2] != 2) continue;
//		Point istart = Point(lineNo[int(lineData[i][0])][2], lineNo[int(lineData[i][0])][3]); //centroid of starting label of current line...
//		double *linedist = new double[line_no];
//		double mindist = 99999;
//		int min = 0;
//		for (int k = 0; k < line_no; k++) {
//			if (lineData[k][2] < 2 || k == i) {
//				linedist[k] = -1;
//				continue;
//			}
//			int startlbl = lineData[k][0];
//			Point kstart = Point(lineNo[startlbl][2], lineNo[startlbl][3]); //centroid of other starting labels...
//			double slope = lineData[k][1];
//			linedist[k] = abs(slope*(istart.x - kstart.x) - (istart.y - kstart.y)) / sqrt(1 + slope*slope);
//			if (linedist[k] < mindist && linedist[k] >= 0) {
//				mindist = linedist[k];
//				min = k;
//			}
//		}
//		if (mindist < AH) {
//			if (lineData[i][2] < lineData[min][2]) {
//				for (int k = 0; k < labelNum; k++) {
//					if (lineNo[k][0] == i + 1) lineNo[k][0] = min + 1;
//				}
//				lineData[i][0] = 0; //starting label...
//				lineData[i][1] = 0; //slope of the line...
//				lineData[i][2] = 0; //no. of cc...
//			}
//			else {
//				for (int k = 0; k < labelNum; k++) {
//					if (lineNo[k][0] == min + 1) lineNo[k][0] = i + 1;
//				}
//				lineData[min][0] = 0; //starting label...
//				lineData[min][1] = 0; //slope of the line...
//				lineData[min][2] = 0; //no. of cc...
//			}
//		}
//
//	}
//
//	//display lines....
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			Vec3b color;
//			if (lines_img.at<Vec3b>(i, j) == Vec3b(255, 255, 255) || lineNo[label[i][j]][0] == 0) continue; //background or unmarked...
//			int rem = lineNo[label[i][j]][0] % 10;
//			color[0] = colors[rem][0];
//			color[1] = colors[rem][1];
//			color[2] = colors[rem][2];
//
//			lines_img.at<Vec3b>(i, j) = color;
//		}
//	}
//
//	imshow("TextLines", lines_img);
//	imwrite("Data/ICDAR/Data218/TxtLine218.tif", lines_img);
//
//
//
//	//---------------------------Delete memory---------------------------
//	for (int i = 0; i < r; i++) {
//		delete[] image[i];
//		delete[] bin_im1[i];
//		delete[] rawimg[i];
//		delete[] label[i];
//		delete[] mainCC[i];
//	}
//	delete[] image;
//	delete[] bin_im1;
//	delete[] rawimg;
//	delete[] label;
//	delete[] mainCC;
//
//	for (int i = 0; i < labelNum; i++) {
//		delete[] box[i];
//		delete[] boundingBox[i];
//		delete[] Data[i];
//	}
//	delete[] box;
//	delete[] boundingBox;
//	delete[] Data;
//	delete[] smallCC;
//	delete[] heights;
//	delete[] diagheights;
//	delete[] trimlabel;
//}
//
//
///*----------------------------------------------------------------*/
///*       Read image and find the boundaries of the main CC's      */
///*----------------------------------------------------------------*/
//void BoundCC(Mat img) {
//	Mat orgimg;
//	orgimg.create(r, c, CV_8UC1);
//	cvtColor(img, orgimg, CV_GRAY2BGR);
//
//	Mat img1, bin_img, cc_img, main_img, thin_img, lines_img;
//	int **rawimg, **image, **bin_im1, **label, **box, **mainCC, *smallCC, **thin_im, *heights, **lineNo, *trimlabel;
//	double h1, w1, area, area_diag, **diag_box, **boundingBox, *diagheights;
//	int count = 0, k = 0, diag_count = 0, rectangle_count = 0;
//
//	//----------------------Blur & Binarize image-----------------------
//	image = new int*[r];
//	rawimg = new int*[r];
//	bin_im1 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		image[i] = new int[c];
//		rawimg[i] = new int[c];
//		bin_im1[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			image[i][j] = 0;
//			rawimg[i][j] = 0;
//			bin_im1[i][j] = 0;
//		}
//	}
//	mat2arr(img, rawimg);
//	mat2arr(img, image);
//	/*blur(image);
//	img1.create(r, c, CV_8UC1);
//	arr2mat(image, img1);
//	imshow("Blurred", img1);*/
//	//double T = Otsu(image, img.rows, img.cols, bin_im1);
//	NICK(image, r, c, bin_im1);
//	bin_img.create(r, c, CV_8UC1);
//	arr2mat(bin_im1, bin_img);
//	cvtColor(bin_img, lines_img, CV_GRAY2BGR);
//	/*erode(bin_im1);
//	arr2mat(bin_im1, bin_img);*/
//	imwrite("Data/ICDAR/Data218/OriginalBinary1.tif", bin_img);
//	/*GaussianBlur(bin_img, bin_img, Size(21, 21), 0, 0);
//	mat2arr(bin_img, image);
//	NICK(image, r, c, bin_im1);
//	arr2mat(bin_im1, bin_img);
//	imwrite("Data/ICDAR/Data218/OriginalBinary.tif", bin_img);*/
//
//	
//	//--------------------------Connected components--------------------------
//	label = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label[i] = new int[c];
//		for (int j = 0; j < c; j++)
//			label[i][j] = 0;
//	}
//	labelNum = labelling(bin_im1, label, r, c);
//	//set_labelNum(labelling(bin_im1, label, r, c));
//	//labelNum = labelNum + 1;
//
//	//------------------Bound the components by rectangular boxes------------------
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	cvtColor(bin_img, cc_img, CV_GRAY2BGR);
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1) { //component not available...
//			count++;
//			continue;
//		}
//		std::cout << i << ":" << box[i][2] << "," << box[i][3] << "," << box[i][4] << "," << box[i][5] << "\n";
//		line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//		line(cc_img, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//		line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//		line(cc_img, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	}
//
//	imshow("Components", cc_img);
//	imwrite("Data/ICDAR/Data218/MaskCC218.tif", cc_img);
//
//	//-------------------Average height of components--------------------
//	int length = labelNum - count;
//	heights = new int[length];
//	for (int i = 0; i < length; i++)
//		heights[i] = 0;
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1)
//			continue;
//		heights[k++] = box[i][5] - box[i][4] + 1;
//		std::cout << "label " << i << " height " << heights[k - 1] << "\n";
//	}
//	AH = averageHeight(heights, length);
//	//set_AH(averageHeight(heights, length));
//	//AH = 35;
//	std::cout << "Average height = " << AH;
//
//	//---------------Delete CC's with small height and length------------
//	smallCC = new int[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) {
//			smallCC[i] = -1;
//			continue;
//		}
//		/*if ((box[i][5] - box[i][4]) > 0.65*AH || (box[i][3] - box[i][2]) > 0.65*AH) {
//			smallCC[i] = 1;
//		}*/
//		if ((box[i][5] - box[i][4]) < 0.5*AH || (box[i][3] - box[i][2]) < 0.5*AH) {
//			smallCC[i] = 0;
//			box[i][2] = -1;
//			box[i][3] = -1;
//			box[i][4] = -1;
//			box[i][5] = -1;
//			//std::cout << "\nLabel " << i << " deleted!";
//		}
//		else smallCC[i] = 1;
//		/*else {
//			smallCC[i] = 0;
//			box[i][2] = -1;
//			box[i][3] = -1;
//			box[i][4] = -1;
//			box[i][5] = -1;
//		}*/
//			
//	}
//
//	mainCC = new int *[r];
//	for (int i = 0; i < r; i++) {
//		mainCC[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			if (smallCC[label[i][j]] <= 0) {
//				label[i][j] = 0;
//				mainCC[i][j] = 255;
//			}
//			else
//				mainCC[i][j] = bin_im1[i][j];
//		}
//	}
//	main_img.create(r, c, CV_8UC1);
//	arr2mat(mainCC, main_img);
//	imshow("Deleting the small CC's", main_img);
//	imwrite("Data/ICDAR/Data218/Nick_bin218.tif", main_img);
//
//	//---------------Diagonal CC box---------------
//	trimlabel = new int[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		trimlabel[i] = -1;
//
//	diag_box = new double*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		diag_box[i] = new double[8];
//		for (int j = 0; j < 8; j++)
//			diag_box[i][j] = 0;
//	}
//	boundingBox = new double*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		boundingBox[i] = new double[9];
//		for (int j = 0; j < 9; j++) {
//			boundingBox[i][j] = -1;
//		}
//	}
//	cvtColor(main_img, cc_img, CV_GRAY2BGR);
//	diagonalCC(label, mainCC, box, boundingBox, diag_box, cc_img, trimlabel);
//
//	//---------------Delete diagonal CC's with small height and length------------
//	smallCC = new int[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1 || boundingBox[i][0] == 1) {
//			continue;
//		}
//		double h1 = sqrt((boundingBox[i][5] - boundingBox[i][1])*(boundingBox[i][5] - boundingBox[i][1]) + (boundingBox[i][6] - boundingBox[i][2])*(boundingBox[i][6] - boundingBox[i][2]));
//		double w1 = sqrt((boundingBox[i][5] - boundingBox[i][3])*(boundingBox[i][5] - boundingBox[i][3]) + (boundingBox[i][6] - boundingBox[i][4])*(boundingBox[i][6] - boundingBox[i][4]));
//		
//		//if (h1 > 0.65*AH || w1 > 0.65*AH) {
//		//	smallCC[i] = 1;
//		//	//std::cout << "\nLabel " << i << " deleted!";
//		//}
//		if (h1 < 0.5*AH || w1 < 0.5*AH) {
//			smallCC[i] = 0;
//			boundingBox[i][0] = 0;
//			boundingBox[i][1] = 0;
//			boundingBox[i][2] = 0;
//			boundingBox[i][3] = 0;
//			boundingBox[i][4] = 0;
//			boundingBox[i][5] = 0;
//			boundingBox[i][6] = 0;
//			boundingBox[i][7] = 0;
//			boundingBox[i][8] = 0;
//			box[i][2] = -1;
//			box[i][3] = -1;
//			box[i][4] = -1;
//			box[i][5] = -1;
//			//std::cout << "\nLabel " << i << " deleted!";
//		}
//		else smallCC[i] = 1;
//		//else {
//		//	smallCC[i] = 0;
//		//	boundingBox[i][0] = 0;
//		//	boundingBox[i][1] = 0;
//		//	boundingBox[i][2] = 0;
//		//	boundingBox[i][3] = 0;
//		//	boundingBox[i][4] = 0;
//		//	boundingBox[i][5] = 0;
//		//	boundingBox[i][6] = 0;
//		//	boundingBox[i][7] = 0;
//		//	boundingBox[i][8] = 0;
//		//	box[i][2] = -1;
//		//	box[i][3] = -1;
//		//	box[i][4] = -1;
//		//	box[i][5] = -1;
//		//	//std::cout << "\nLabel " << i << " deleted!";
//		//}
//	}
//
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (smallCC[label[i][j]] == 0) {
//				label[i][j] = 0;
//				mainCC[i][j] = 255;
//			}
//		}
//	}
//
//	//----------------Dilate and fill up holes in each CC-----------------
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//		int i_start = box[i][4], j_start = box[i][2];
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//		int **roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				if (label[i_start + k][j_start + l] == i) {
//					roi[k][l] = mainCC[i_start + k][j_start + l];
//				}
//				else
//					roi[k][l] = 255;
//			}
//
//		}
//
//		fill(roi, h, w);
//
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				if (roi[k][l] == 0 ) {
//					label[i_start + k][j_start + l] = i;
//					mainCC[i_start + k][j_start + l] = roi[k][l] ;
//				}
//			}
//		}
//
//		//Delete roi...
//		for (int k = 0; k < h; k++)
//			delete[] roi[k];
//		delete[] roi;
//	}
//
//	arr2mat(mainCC, main_img);
//	imshow("Deleting the small CC's(2)", main_img);
//	imwrite("Data/ICDAR/Data218/DeletedsmallCC.tif", main_img);
//
//	//-------------------Average height of components--------------------
//	count = 0;
//
//	cvtColor(main_img, cc_img, CV_GRAY2BGR);
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1 || boundingBox[i][0] == -1) { //component not available...
//			count++;
//			continue;
//		}
//		if (boundingBox[i][0] == 2) { //Diagonal box....
//			line(cc_img, Point(diag_box[i][1], diag_box[i][0]), Point(diag_box[i][5], diag_box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(diag_box[i][5], diag_box[i][4]), Point(diag_box[i][3], diag_box[i][2]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(diag_box[i][3], diag_box[i][2]), Point(diag_box[i][7], diag_box[i][6]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(diag_box[i][7], diag_box[i][6]), Point(diag_box[i][1], diag_box[i][0]), Scalar(0, 0, 255), 1.5, 8);
//		}
//		if (boundingBox[i][0] == 1) { //Upright box....
//			line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//		}
//	}
//	imshow("Components", cc_img);
//	imwrite("Data/ICDAR/Data218/DiagBox218.tif", cc_img);
//
//	length = labelNum - count;
//	diagheights = new double[length];
//	for (int i = 0; i < length; i++)
//		diagheights[i] = 0;
//	k = 0;
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1 || boundingBox[i][0] == -1)
//			continue;
//		if(boundingBox[i][0] == 1) //upright
//			diagheights[k++] = box[i][5] - box[i][4] + 1;
//		if (boundingBox[i][0] == 2) //diagonal
//			diagheights[k++] = sqrt((boundingBox[i][5] - boundingBox[i][1])*(boundingBox[i][5] - boundingBox[i][1]) + (boundingBox[i][6] - boundingBox[i][2])*(boundingBox[i][6] - boundingBox[i][2]));
//		std::cout << "\nlabel " << i << " height " << diagheights[k - 1] << "\n";
//	}
//	AH = averageDiagHeight(diagheights, length);
//	std::cout << "Average height = " << AH;
//
//	//----------------Trim components to remove accent marks and extensions-------------
//	std::cout << "\nTrimming directions :";
//	for (int i = 0; i < labelNum; i++) {
//		if (trimlabel[i] == -1)
//			continue;
//		if (trimlabel[i] == 0)
//			std::cout << "\nLabel" << i << ": undecided";
//		if (trimlabel[i] == 1)
//			std::cout << "\nLabel" << i << ": 0 degree";
//		if (trimlabel[i] == 2)
//			std::cout << "\nLabel" << i << ": 90 degree";
//		if (trimlabel[i] == 3)
//			std::cout << "\nLabel" << i << ": 45 degree";
//		if (trimlabel[i] == 4)
//			std::cout << "\nLabel" << i << ": 135 degree";
//	}
//	
//	trimHV(mainCC, box, label, boundingBox, trimlabel);
//	trimV(mainCC, box, label, boundingBox, trimlabel);
//	arr2mat(mainCC, main_img);
//	imshow("Deleting the small CC's(3)", main_img);
//	imwrite("Data/ICDAR/Data218/DeletedsmallCC(trimmed).tif", main_img);
//
//	//-----------------Re - Label--------------
//	for (int i = 0; i < labelNum; i++) {
//		delete[] box[i];
//		delete[] boundingBox[i];
//	}
//	delete[] box;
//	delete[] boundingBox;
//	delete[] smallCC;
//	delete[] heights;
//	delete[] diagheights;
//	delete[] trimlabel;
//
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++)
//			label[i][j] = 0;
//	}
//	labelNum = labelling(mainCC, label, r, c);
//
//	//------------------Bound the components by rectangular boxes------------------
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	cvtColor(main_img, cc_img, CV_GRAY2BGR);
//	count = 0;
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1) { //component not available...
//			count++;
//			continue;
//		}
//		std::cout << i << ":" << box[i][2] << "," << box[i][3] << "," << box[i][4] << "," << box[i][5] << "\n";
//		line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//		line(cc_img, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//		line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//		line(cc_img, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	}
//
//	imshow("NewComponents", cc_img);
//	imwrite("Data/ICDAR/Data218/TrimmedCC218.tif", cc_img);
//
//	//-------------------Average height of components--------------------
//	length = labelNum - count;
//	heights = new int[length];
//	for (int i = 0; i < length; i++)
//		heights[i] = 0;
//	k = 0;
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1)
//			continue;
//		heights[k++] = box[i][5] - box[i][4] + 1;
//		std::cout << "label " << i << " height " << heights[k - 1] << "\n";
//	}
//	AH = averageHeight(heights, length);
//	std::cout << "Average height = " << AH;
//
//	//---------------Delete CC's with small height and length------------
//	smallCC = new int[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) {
//			smallCC[i] = -1;
//			continue;
//		}
//		if ((box[i][5] - box[i][4]) < 0.3*AH || (box[i][3] - box[i][2]) < 0.3*AH) {
//			smallCC[i] = 0;
//			box[i][2] = -1;
//			box[i][3] = -1;
//			box[i][4] = -1;
//			box[i][5] = -1;
//			//std::cout << "\nLabel " << i << " deleted!";
//		}
//		else smallCC[i] = 1;
//	}
//
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (smallCC[label[i][j]] <= 0) {
//				label[i][j] = 0;
//				mainCC[i][j] = 255;
//			}
//		}
//	}
//	main_img.create(r, c, CV_8UC1);
//	arr2mat(mainCC, main_img);
//	imshow("Deleting the small CC's", main_img);
//	imwrite("Data/ICDAR/Data218/Nick_bin218(trimmed).tif", main_img);
//
//	//---------------Diagonal CC box---------------
//	trimlabel = new int[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		trimlabel[0] = 0;
//
//	diag_box = new double*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		diag_box[i] = new double[8];
//		for (int j = 0; j < 8; j++)
//			diag_box[i][j] = 0;
//	}
//	boundingBox = new double*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		boundingBox[i] = new double[9];
//		for (int j = 0; j < 9; j++) {
//			boundingBox[i][j] = -1;
//		}
//	}
//	cvtColor(main_img, cc_img, CV_GRAY2BGR);
//	diagonalCC(label, mainCC, box, boundingBox, diag_box, cc_img, trimlabel);
//
//	//---------------Delete diagonal CC's with small height and length------------
//	smallCC = new int[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1 || boundingBox[i][0] == 1) {
//			continue;
//		}
//		double h1 = sqrt((boundingBox[i][5] - boundingBox[i][1])*(boundingBox[i][5] - boundingBox[i][1]) + (boundingBox[i][6] - boundingBox[i][2])*(boundingBox[i][6] - boundingBox[i][2]));
//		double w1 = sqrt((boundingBox[i][5] - boundingBox[i][3])*(boundingBox[i][5] - boundingBox[i][3]) + (boundingBox[i][6] - boundingBox[i][4])*(boundingBox[i][6] - boundingBox[i][4]));
//
//		if (h1 < 0.5*AH || w1 < 0.5*AH) {
//			smallCC[i] = 0;
//			boundingBox[i][0] = 0;
//			boundingBox[i][1] = 0;
//			boundingBox[i][2] = 0;
//			boundingBox[i][3] = 0;
//			boundingBox[i][4] = 0;
//			boundingBox[i][5] = 0;
//			boundingBox[i][6] = 0;
//			boundingBox[i][7] = 0;
//			boundingBox[i][8] = 0;
//			box[i][2] = -1;
//			box[i][3] = -1;
//			box[i][4] = -1;
//			box[i][5] = -1;
//			//std::cout << "\nLabel " << i << " deleted!";
//		}
//		else smallCC[i] = 1;
//		
//	}
//
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (smallCC[label[i][j]] == 0) {
//				label[i][j] = 0;
//				mainCC[i][j] = 255;
//			}
//		}
//	}
//
//	//-------------------Average height of components--------------------
//	count = 0;
//	arr2mat(mainCC, main_img);
//	cvtColor(main_img, cc_img, CV_GRAY2BGR);
//	double **Data = new double*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		Data[i] = new double[6];
//		for (int j = 0; j < 6; j++) {
//			Data[i][j] = 0;
//		}
//	}
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1 || boundingBox[i][0] == -1) { //component not available...
//			count++;
//			continue;
//		}
//		if (boundingBox[i][0] == 2) { //Diagonal box....
//			line(cc_img, Point(diag_box[i][1], diag_box[i][0]), Point(diag_box[i][5], diag_box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(diag_box[i][5], diag_box[i][4]), Point(diag_box[i][3], diag_box[i][2]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(diag_box[i][3], diag_box[i][2]), Point(diag_box[i][7], diag_box[i][6]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(diag_box[i][7], diag_box[i][6]), Point(diag_box[i][1], diag_box[i][0]), Scalar(0, 0, 255), 1.5, 8);
//		}
//		if (boundingBox[i][0] == 1) { //Upright box....
//			line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//			line(cc_img, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//		}
//	}
//	imshow("Components", cc_img);
//	imwrite("Data/ICDAR/Data218/DiagBox218(trimmed).tif", cc_img);
//
//	length = labelNum - count;
//	diagheights = new double[length];
//	for (int i = 0; i < length; i++)
//		diagheights[i] = 0;
//	k = 0;
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1 || boundingBox[i][0] == -1)
//			continue;
//		if (boundingBox[i][0] == 1) //upright
//			diagheights[k++] = box[i][5] - box[i][4] + 1;
//		if (boundingBox[i][0] == 2) //diagonal
//			diagheights[k++] = sqrt((boundingBox[i][5] - boundingBox[i][1])*(boundingBox[i][5] - boundingBox[i][1]) + (boundingBox[i][6] - boundingBox[i][2])*(boundingBox[i][6] - boundingBox[i][2]));
//		std::cout << "\nlabel " << i << " height " << diagheights[k - 1] << "\n";
//	}
//	AH = averageDiagHeight(diagheights, length);
//	std::cout << "Average height = " << AH;
//
//	//-----------------Find major axis for each CC and connect---------------------------
//	//skewLine(mainCC, cc_img, label, boundingBox, box);
//	lineNo = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		lineNo[i] = new int[4];
//		lineNo[i][0] = 0; //line no. for i-th label...
//		lineNo[i][1] = 0; //1 if i-th label is the starting of the line, 2 if i-th label is a single word...
//		lineNo[i][2] = 0; //centroid.x...
//		lineNo[i][3] = 0;//centroidy...
//	}
//
//	double **connect = new double*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		connect[i] = new double[3];
//		connect[i][0] = i; //label no....
//		connect[i][1] = -1;  //next label...
//		connect[i][2] = 0;  //distance...
//	}
//
//	partition(orgimg, mainCC, cc_img, main_img, label, boundingBox, diag_box, box, Data, lineNo, connect);
//
//	int colors[10][3] = { {0,0,255},{0,255,0},{255,0,0},{128,0,255}, {128,128,0},{255,128,0},{255,0,128}, {128,255,0},{0,255,128}, {255,50,255} };
//
//	std::cout << "\nTotal no. of text lines = " << line_no;
//
//	//---------------------assign line no. to each text pixel----------------
//	//collect the starting CC and the inclination of each line in the document...
//	double **lineData = new double*[line_no];
//	for (int i = 0; i < line_no; i++) {
//		lineData[i] = new double[3];
//		lineData[i][0] = 0; //starting label...
//		lineData[i][1] = 0; //slope of the line...
//		lineData[i][2] = 0; //no. of cc...
//	}
//
//	for (int i = 0; i < labelNum; i++) {
//		if (lineNo[i][1] == 2) { //single word...
//			int lno = lineNo[i][0] - 1;
//			lineData[lno][0] = i;
//			lineData[lno][2] = 1;
//			continue;
//		}
//
//		if (lineNo[i][1] != 1) continue; //not the starting of line...
//		
//		int lno = lineNo[i][0] - 1;
//		lineData[lno][0] = i;
//		lineData[lno][2]++;
//
//		double X = lineNo[i][2], Y = lineNo[i][3], XY = X*Y, X2 = X*X;
//		int curr = i, next = connect[i][1], count = 1;
//		while (next > 0) {
//			curr = next;
//			X += lineNo[curr][2];
//			Y += lineNo[curr][3];
//			XY += (lineNo[curr][2] * lineNo[curr][3]);
//			X2 += (lineNo[curr][2] * lineNo[curr][2]);
//			lineData[lno][2]++;
//			count++;
//			next = connect[curr][1];
//		}
//		X += lineNo[curr][2];
//		Y += lineNo[curr][3];
//		XY += (lineNo[curr][2] * lineNo[curr][3]);
//		X2 += (lineNo[curr][2] * lineNo[curr][2]);
//		lineData[lno][2]++;
//		count++;
//
//		lineData[lno][1] = (count*XY - X*Y) / (count*X2 - X*X); //regression slope...
//		/*Point p1 = Point(lineNo[i][2], lineNo[i][3]);
//		Point p2 = Point(lineNo[curr][2], lineNo[i][3]+ lineData[lno][1]*(lineNo[curr][2]- lineNo[i][2]));
//		line(lines_img, p1, p2, Scalar(255, 0, 255), 2, 8, 0);*/
//
//	}
//
//	std::cout << "\nLine data::";
//	for (int i = 0; i < line_no; i++) {
//		std::cout << "\nLine " << i + 1;
//		std::cout << "\nStarts at : label " << lineData[i][0];
//		std::cout << "\nSlope =  " << lineData[i][1];
//		std::cout << "\nNo. of words = " << lineData[i][2];
//	}
//
//	//re-assign single words...
//	for (int i = 0; i < line_no; i++) {
//		if (lineData[i][2] != 1 ) continue;
//		double *linedist = new double[line_no];
//		double mindist = 99999;
//		int min = 0;
//		Point pt = Point(lineNo[int(lineData[i][0])][2], lineNo[int(lineData[i][0])][3]);
//		for (int k = 0; k < line_no; k++) {
//			if (lineData[k][2] < 2) {
//				linedist[k] = -1;
//				continue;
//			}
//			int startlbl = lineData[k][0];
//			Point start = Point(lineNo[startlbl][2], lineNo[startlbl][3]);
//			double slope = lineData[k][1];
//			linedist[k] = abs(slope*(pt.x - start.x) - (pt.y - start.y)) / sqrt(1 + slope*slope);
//			if (linedist[k] < mindist && linedist[k] >= 0) {
//				mindist = linedist[k];
//				min = k;
//			}
//		}
//		if (mindist < AH) {
//			lineNo[int(lineData[i][0])][0] = min + 1;
//			lineData[min][2]++;
//			lineData[i][2]--;
//			std::cout << "\nLabel " << lineData[i][0] << " added to line #" << min + 1 << " at distance " << mindist;
//		}
//		delete[] linedist;
//	}
//
//	//merge lines with same slope and intercept...
//	for (int i = 0; i < line_no; i++) {
//		if (lineData[i][2] != 2) continue;
//		Point istart = Point(lineNo[int(lineData[i][0])][2], lineNo[int(lineData[i][0])][3]); //centroid of starting label of current line...
//		double *linedist = new double[line_no];
//		double mindist = 99999;
//		int min = 0;
//		for (int k = 0; k < line_no; k++) {
//			if (lineData[k][2] < 2 || k == i) {
//				linedist[k] = -1;
//				continue;
//			}
//			int startlbl = lineData[k][0];
//			Point kstart = Point(lineNo[startlbl][2], lineNo[startlbl][3]); //centroid of other starting labels...
//			double slope = lineData[k][1];
//			linedist[k] = abs(slope*(istart.x - kstart.x) - (istart.y - kstart.y)) / sqrt(1 + slope*slope);
//			if (linedist[k] < mindist && linedist[k] >= 0) {
//				mindist = linedist[k];
//				min = k;
//			}
//		}
//		if (mindist < AH) {
//			if (lineData[i][2] < lineData[min][2]) {
//				for (int k = 0; k < labelNum; k++) {
//					if (lineNo[k][0] == i + 1) lineNo[k][0] = min + 1;
//				}
//				lineData[i][0] = 0; //starting label...
//				lineData[i][1] = 0; //slope of the line...
//				lineData[i][2] = 0; //no. of cc...
//			}
//			else {
//				for (int k = 0; k < labelNum; k++) {
//					if (lineNo[k][0] == min + 1) lineNo[k][0] = i + 1;
//				}
//				lineData[min][0] = 0; //starting label...
//				lineData[min][1] = 0; //slope of the line...
//				lineData[min][2] = 0; //no. of cc...
//			}
//		}
//
//	}
//
//	//display lines....
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			Vec3b color;
//			if (lines_img.at<Vec3b>(i, j) == Vec3b(255, 255, 255) || lineNo[label[i][j]][0] == 0) continue; //background or unmarked...
//			int rem = lineNo[label[i][j]][0] % 10;
//			color[0] = colors[rem][0];
//			color[1] = colors[rem][1];
//			color[2] = colors[rem][2];
//
//			lines_img.at<Vec3b>(i, j) = color;
//		}
//	}
//
//	imshow("TextLines", lines_img);
//	imwrite("Data/ICDAR/Data218/TxtLine218.tif", lines_img);
//
//
//	//---------------------------Delete memory---------------------------
//	for (int i = 0; i < r; i++) {
//		delete[] image[i];
//		delete[] bin_im1[i];
//		delete[] label[i];
//		delete[] mainCC[i];
//	}
//	delete[] image;
//	delete[] bin_im1;
//	delete[] label;
//	delete[] mainCC;
//
//	for (int i = 0; i < labelNum; i++) {
//		delete[] box[i];
//		delete[] boundingBox[i];
//	}
//	delete[] box;
//	delete[] boundingBox;
//	delete[] smallCC;
//	delete[] heights;
//	delete[] diagheights;
//}
//
//
////Trim CCs...
//void TrimCC(Mat img) {
//	Mat orgimg;
//	orgimg.create(r, c, CV_8UC1);
//	cvtColor(img, orgimg, CV_GRAY2BGR);
//
//	Mat img1, bin_img, cc_img, main_img, thin_img, lines_img;
//	int **rawimg, **image, **bin_im1, **label, **box, **mainCC, *smallCC, **thin_im, *heights, **lineNo, *trimlabel;
//	double h1, w1, area, area_diag, **diag_box, **boundingBox, *diagheights;
//	int count = 0, k = 0, diag_count = 0, rectangle_count = 0;
//
//	//----------------------Blur & Binarize image-----------------------
//	image = new int*[r];
//	rawimg = new int*[r];
//	bin_im1 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		image[i] = new int[c];
//		rawimg[i] = new int[c];
//		bin_im1[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			image[i][j] = 0;
//			rawimg[i][j] = 0;
//			bin_im1[i][j] = 0;
//		}
//	}
//	mat2arr(img, rawimg);
//	mat2arr(img, image);
//	NICK(image, r, c, bin_im1);
//	bin_img.create(r, c, CV_8UC1);
//	arr2mat(bin_im1, bin_img);
//	cvtColor(bin_img, lines_img, CV_GRAY2BGR);
//	imwrite("Data/ICDAR/Data218/OriginalBinary1.tif", bin_img);
//
//
//	//--------------------------Connected components--------------------------
//	label = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label[i] = new int[c];
//		for (int j = 0; j < c; j++)
//			label[i][j] = 0;
//	}
//	labelNum = labelling(bin_im1, label, r, c);
//
//	//------------------Bound the components by rectangular & diagonal boxes------------------
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	//-------------------Average height of components--------------------
//	int length = labelNum - count;
//	heights = new int[length];
//	for (int i = 0; i < length; i++)
//		heights[i] = 0;
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1)
//			continue;
//		heights[k++] = box[i][5] - box[i][4] + 1;
//		std::cout << "label " << i << " height " << heights[k - 1] << "\n";
//	}
//	AH = averageHeight(heights, length);
//	std::cout << "Average height = " << AH;
//
//	//---------------Delete CC's with small height and length------------
//	smallCC = new int[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) {
//			smallCC[i] = -1;
//			continue;
//		}
//		if ((box[i][5] - box[i][4]) < 0.25*AH || (box[i][3] - box[i][2]) < 0.25*AH) {
//			smallCC[i] = 0;
//			box[i][2] = -1;
//			box[i][3] = -1;
//			box[i][4] = -1;
//			box[i][5] = -1;
//		}
//		else smallCC[i] = 1;
//	}
//
//	mainCC = new int *[r];
//	for (int i = 0; i < r; i++) {
//		mainCC[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			if (smallCC[label[i][j]] <= 0) {
//				label[i][j] = 0;
//				mainCC[i][j] = 255;
//			}
//			else
//				mainCC[i][j] = bin_im1[i][j];
//		}
//	}
//	main_img.create(r, c, CV_8UC1);
//	arr2mat(mainCC, main_img);
//	imshow("Deleting the small CC's", main_img);
//	imwrite("Data/ICDAR/Data218/Nick_bin218.tif", main_img);
//
//	diag_box = new double*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		diag_box[i] = new double[8];
//		for (int j = 0; j < 8; j++)
//			diag_box[i][j] = 0;
//	}
//	
//	boundingBox = new double*[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		boundingBox[i] = new double[9];
//		for (int j = 0; j < 9; j++) {
//			boundingBox[i][j] = -1;
//		}
//	}
//
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1)
//			continue;
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//
//		int NE = 0, SE = 0, SW = 0, NW = 0;
//		int **roi, i_start = box[i][4], j_start = box[i][2];
//		roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				if (label[i_start + k][j_start + l] == i)
//					roi[k][l] = mainCC[i_start + k][j_start + l];
//				else
//					roi[k][l] = 255;
//				if (roi[k][l] == 0) count++;
//
//			}
//		}
//		std::cout << "\nLabel " << i << " : " << box[i][2] << "," << box[i][3] << "," << box[i][4] << "," << box[i][5] << "\n";
//
//		diagonalBox(roi, h, w, NW, SE, NE, SW);
//		diag_box[i][0] = box[i][4] + (NW - NE) / 2;  // N.y
//		diag_box[i][1] = box[i][2] + (NW + NE) / 2;  // N.x
//		diag_box[i][2] = box[i][4] + (SE - SW) / 2;  // S.y
//		diag_box[i][3] = box[i][2] + (SE + SW) / 2;  // S.x
//		diag_box[i][4] = box[i][4] + (SE - NE) / 2;  // E.y
//		diag_box[i][5] = box[i][2] + (SE + NE) / 2;  // E.x
//		diag_box[i][6] = box[i][4] + (NW - SW) / 2;  // W.y
//		diag_box[i][7] = box[i][2] + (NW + SW) / 2;  // W.x
//
//
//		fill(roi, h, w);
//
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				if (roi[k][l] == 0) {
//					label[i_start + k][j_start + l] = i;
//					mainCC[i_start + k][j_start + l] = roi[k][l];
//				}
//			}
//		}
//	}
//
//	arr2mat(mainCC, main_img);
//	imshow("Deleting the small CC's(1)", main_img);
//
//	trimCC(mainCC, label, box, diag_box);
//
//	arr2mat(mainCC, main_img);
//	imshow("Deleting the small CC's(2)", main_img);
//	imwrite("Data/ICDAR/Data218/Nick_bin218.tif", main_img);
//
//	MaskCC(main_img, mainCC);
//
//	//---------------------------Delete memory---------------------------
//	for (int i = 0; i < r; i++) {
//		delete[] image[i];
//		delete[] bin_im1[i];
//		delete[] label[i];
//		delete[] mainCC[i];
//	}
//	delete[] image;
//	delete[] bin_im1;
//	delete[] label;
//	delete[] mainCC;
//	
//	for (int i = 0; i < labelNum; i++) {
//		delete[] box[i];
//		delete[] boundingBox[i];
//	}
//	delete[] box;
//	delete[] boundingBox;
//	delete[] smallCC;
//}
//
//
///*----------------------------------------------------------------*/
///*                        Main function                           */
///*----------------------------------------------------------------*/
////int main(int argc, char** argv) {
////	//----------------------Retrieve .dat file------------------------
////	//int Ix = 1000; //Image Width (x=0...Ix-1)
////	//int Iy = 1000; //Image Height (y=0...Iy-1)
////	//unsigned int *IM_SegmResult; //Pointer to store raw data
////	//std::string SegmResultName; //File that contains the raw data
////	//FILE *f1;
////
////	//IM_SegmResult = (unsigned int *)calloc(Ix*Iy, sizeof(int));
////	//f1 = fopen("Data/data13.dat","r");
////	//fread(IM_SegmResult, Ix*Iy, sizeof(int), f1);
////	//fclose(f1);
////
////	//
////	//Mat im;
////	//im.create(Ix, Iy, CV_8UC1);
////	//int **imdat = new int*[Ix];
////	//for (int i = 0; i < Ix; i++) {
////	//	imdat[i] = new int[Iy];
////	//	for (int j = 0; j < Iy; j++)
////	//		imdat[i][j] = 255-int(IM_SegmResult[i*Ix + j]) * 100;
////	//}
////	//arr2mat(imdat, im, Ix, Iy);
////	//imshow("dat",im);
////
////	//std::cout << "\nDone!\n";
////
////
////	//-----------------------------Read image-----------------------
////	Mat img = imread("Data/ICDAR/Data218/scan218.jpg", IMREAD_GRAYSCALE);
////	imshow("Image", img);
////	std::cout << "Image dimensions:" << img.rows << "x" << img.cols << endl;
////
////	GaussianBlur(img, img, Size(3,3), 0, 0);
////	imshow("Image1", img);
////	r = img.rows;
////	c = img.cols;
////	/*set_r(img.rows);
////	set_c(img.cols);*/
////
////	//MaskCC(img);
////	//BoundCC(img);
////	TrimCC(img);;
////
////	std::cout << "\nDone!\n";
////	waitKey(0);
////	system("pause");
////	return 0;
////
////}
//
