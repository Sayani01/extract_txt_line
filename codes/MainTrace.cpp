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
//#include "FPTA.h"
//#include "GHPTA.h"
//#include "OIC.h"
//
//
//using namespace std;
//using namespace cv;
//namespace cv {
//	using std::vector;
//};
//
//extern int labelNum = 0, r = 0, c = 0;
//extern double AH = 0;
//
//int grid = 2;
//
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
//void arr2mat(int **arr, Mat A, int row, int col) {
//	for (int i = 0; i < row; i++) {
//		for (int j = 0; j < col; j++) {
//			A.at<uchar>(i, j) = arr[i][j];
//		}
//	}
//}
//
////bckground=white ; foreground=black
//int checkNbd(int **img, int h, int w, Point ref) {
//	int sum = 0;
//	for (int i = -1; i <= 1; i++) {
//		for (int j = -1; j <= 1; j++) {
//			if (ref.x + j < 0 || ref.x + j >= w || ref.y + i < 0 || ref.y + i >= h || (i == 0 && j == 0))
//				continue;
//			else {
//				if (img[ref.y + i][ref.x + j] == 0)
//					sum++;
//			}
//		}
//	}
//	return(sum);
//}
//
//int findTerminals(int **img, int h, int w, std::vector<Point> &terminal , std::vector<Point> &next) {
//	int **copy, sz = 0;
//	copy = new int*[h + 2];
//	for (int i = 0; i < h + 2; i++) {
//		copy[i] = new int[w + 2];
//		for (int j = 0; j < w + 2; j++) {
//			if (i == 0 || j == 0 || i == h + 1 || j == w + 1)
//				copy[i][j] = 0;
//			else
//				copy[i][j] = (255 - img[i - 1][j - 1]) / 255;
//		}
//	}
//	for (int i = 1; i <= h; i++) {
//		for (int j = 1; j <= w; j++) {
//			if (copy[i][j] == 0) continue;
//			int sum = checkNbd(copy, h+2, w+2, Point(j, i));
//			if (sum == 7) {
//				terminal.push_back(Point(j - 1, i - 1));
//				sz++;
//			}
//		}
//	}
//	return(sz);
//}
//
//int isTerminal(std::vector<Point> terminal, Point p) {
//	for (int i = 0; i<int(terminal.size()); i++) {
//		if (terminal[i] == p)
//			return(1);
//	}
//	return(0);
//}
//
////Forward tracing...
//int TracePWLine0(int **img, int h, int w, std::vector<Point> &sequence, std::vector<Point> terminals) {
//	int sz = 1;
//	int **copy;
//	copy = new int*[h];
//	for (int i = 0; i < h ; i++) {
//		copy[i] = new int[w];
//		for (int j = 0; j < w ; j++) {
//			copy[i][j] = (255 - img[i][j]) / 255; //1  - FG, 0 - BG
//		}
//	}
//
//	do {
//		Point nxt;
//		Point curr = sequence[sz - 1];
//		double dist = 999;
//		int flag = 0;
//		for (int k = 0; k < h; k++) {
//			if (copy[k][curr.x + 6] == 1) {
//				double d = sqrt((k - curr.y)*(k - curr.y) + 36);
//				if (d < 12) {
//					flag = 1;
//					if (d < dist) {
//						dist = d;
//						nxt = Point(curr.x + 6, k);
//					}
//				}				
//			}
//		}
//		if (flag == 1) {     //next point available...
//			/*sequence.push_back(nxt);
//			sz++;*/
//			double slope = atan2(nxt.y - curr.y, nxt.x - curr.x) * 180 / CV_PI;
//			//if (slope <= -90) slope = -(180 + slope);
//			std::cout << "\nslope = " << slope;
//			if (abs(slope) <= 90 ) {    //check for non horizontal...
//				sequence.push_back(nxt);
//				sz++;
//			}
//			else
//				flag = 0;
//		}
//		if(flag == 0){
//			int f = 0;
//			double dist = 9999;
//			for (int l = 0; l < 6; l++) {
//				for (int k = 0; k < h; k++) {
//					if (isTerminal(terminals,Point(curr.x + l,k)) ) {
//						double d = sqrt((k - curr.y)*(k - curr.y) + l*l);
//						if (d < 7) {
//							f = 1;
//							if (d < dist && d < 7) {
//								dist = d;
//								nxt = Point(curr.x + l, k);
//							}
//						}						
//					}
//				}
//			}
//			if (f == 1) {
//				sequence.push_back(nxt);
//				sz++;
//			}
//			break;
//		}
//	} while (1);
//	return(sz);
//}
//
//int TracePWLine135(int **img, int h, int w, std::vector<Point> &sequence, std::vector<Point> terminals) {
//	int sz = 1;
//	int **copy;
//	copy = new int*[h];
//	for (int i = 0; i < h; i++) {
//		copy[i] = new int[w];
//		for (int j = 0; j < w; j++) {
//			copy[i][j] = (255 - img[i][j]) / 255; //1  - FG, 0 - BG
//		}
//	}
//
//	do {
//		Point nxt;
//		Point curr = sequence[sz - 1];
//		double dist = 9999.9;
//		int flag = 0, intercept = -curr.x + curr.y - 6;
//		for (int k = 0; k < 2; k++) {
//			intercept -= k;
//			for (int y = 0; y < h; y++) {
//				int x = y - intercept;
//				if (x < 0 || x >= w)
//					continue;
//				if (copy[y][x] == 1) {
//					double d = sqrt((double(y - curr.y)*(y - curr.y)) + double((x - curr.x)*(x - curr.x)));
//					//std::cout << ".....distance = " << d;
//					if (d <= 12) {
//						flag = 1;
//						if (d < dist) {
//							dist = d;
//							nxt = Point(x, y);
//						}
//					}
//				}
//			}
//		}
//		if (flag == 1) {     //next point available...
//			/*sequence.push_back(nxt);
//			sz++;*/
//			double slope = atan2(nxt.y - curr.y, nxt.x - curr.x) * 180 / CV_PI;
//
//			if (abs(slope + 45) <= 60) {    //check for non horizontal...
//				sequence.push_back(nxt);
//				sz++;
//				//std::cout << "\n..." << nxt.y << " , " << nxt.x;
//		    }	
//			else
//				flag = 0;
//		}
//		if (flag == 0) {  //find the nearest terminal point...
//			double dist = 9999.9;
//			int f = 0;
//			for (int intercept = -curr.x + curr.y; intercept >= -curr.x + curr.y - 6; intercept--) {
//				for (int y = 0; y < h; y++) {
//					int x = y - intercept;
//					if (x < 0 || x >= w) continue;
//					if (isTerminal(terminals, Point(x, y))) {
//						double d = sqrt((double(y - curr.y)*(y - curr.y)) + double((x - curr.x)*(x - curr.x)));
//						if (d < 5) {
//							f = 1;
//							if (d < dist) {
//								dist = d;
//								nxt = Point(x, y);
//							}
//						}						
//					}
//				}					
//			}
//			if (f == 1) {
//				sequence.push_back(nxt);
//				sz++;
//			}
//			break;
//		}
//	} while (1);
//	return(sz);
//}
//
//int TracePWLine45(int **img, int h, int w, std::vector<Point> &sequence, std::vector<Point> terminals) {
//	int sz = 1;
//	int **copy;
//	copy = new int*[h];
//	for (int i = 0; i < h; i++) {
//		copy[i] = new int[w];
//		for (int j = 0; j < w; j++) {
//			copy[i][j] = (255 - img[i][j]) / 255; //1  - FG, 0 - BG
//		}
//	}
//
//	do {
//		Point nxt;
//		Point curr = sequence[sz - 1];
//		double dist = 9999.9;
//		int flag = 0, intercept = curr.x + curr.y + 6;
//		for (int k = 0; k < 2; k++) {
//			intercept += k;
//			for (int y = 0; y < h; y++) {
//				int x = -y + intercept;
//				if (x < 0 || x >= w) continue;
//				if (copy[y][x] == 1) {
//					double d = sqrt((double(y - curr.y)*(y - curr.y)) + double((x - curr.x)*(x - curr.x)));
//					//std::cout << ".....distance = " << d;
//					flag = 1;
//					if (d < dist) {
//						dist = d;
//						nxt = Point(x, y);
//					}
//				}
//			}
//
//		}
//		if (flag == 1) {     //next point available...
//							 /*sequence.push_back(nxt);
//							 sz++;*/
//			double slope = atan2(nxt.y - curr.y, nxt.x - curr.x) * 180 / CV_PI;
//
//			if (abs(slope - 45) <= 90) {    //check for non horizontal...
//				sequence.push_back(nxt);
//				sz++;
//				//std::cout << "\n..." << nxt.y << " , " << nxt.x;
//			}
//			else
//				flag = 0;
//		}
//		if (flag == 0) {  //find the nearest terminal point...
//			double dist = 9999.9;
//			int f = 0;
//			for (int intercept = curr.x + curr.y; intercept <= curr.x + curr.y + 6; intercept++) {
//				for (int y = 0; y < h; y++) {
//					int x = -y + intercept;
//					if (x < 0 || x >= w) continue;
//					if (isTerminal(terminals, Point(x, y))) {
//						double d = sqrt((double(y - curr.y)*(y - curr.y)) + double((x - curr.x)*(x - curr.x)));
//						if (d < 6) {
//							f = 1;
//							if (d < dist) {
//								dist = d;
//								nxt = Point(x, y);
//							}
//						}						
//					}
//				}
//			}
//			if (f == 1) {
//				sequence.push_back(nxt);
//				sz++;
//			}
//			break;
//		}
//	} while (1);
//	return(sz);
//}
//
////Backward tracing...
//int BackTracePWLine0(int **img, int h, int w, std::vector<Point> &sequence, std::vector<Point> terminals) {
//	int sz = 1;
//	int **copy;
//	copy = new int*[h];
//	for (int i = 0; i < h; i++) {
//		copy[i] = new int[w];
//		for (int j = 0; j < w; j++) {
//			copy[i][j] = (255 - img[i][j]) / 255; //1  - FG, 0 - BG
//		}
//	}
//
//	do {
//		Point nxt;
//		Point curr = sequence[sz - 1];
//		double dist = 999;
//		int flag = 0;
//		for (int k = 0; k < h; k++) {
//			if (copy[k][curr.x - 6] == 1) {
//				double d = sqrt((k - curr.y)*(k - curr.y) + 36);
//				if (d < 10) {
//					flag = 1;
//					if (d < dist) {
//						dist = d;
//						nxt = Point(curr.x - 6, k);
//					}
//				}
//			}
//		}
//		if (flag == 1) {     //next point available...
//			/*sequence.push_back(nxt);
//			sz++;*/
//			double slope = atan2(nxt.y - curr.y, nxt.x - curr.x) * 180 / CV_PI;
//			////if (slope <= -90) slope = -(180 + slope);
//			std::cout << "\nslope = " << slope;
//			if (180 - abs(slope) <= 90 ) {    //check for non horizontal...
//				sequence.push_back(nxt);
//				sz++;
//			}
//			else
//				flag = 0;
//		}
//		if (flag == 0) {
//			int f = 0;
//			double dist = 9999;
//			for (int l = 0; l < 6; l++) {
//				for (int k = 0; k < h; k++) {
//					if (isTerminal(terminals, Point(curr.x - l, k))) {
//						double d = sqrt((k - curr.y)*(k - curr.y) + l*l);
//						if (d < 10) {
//							f = 1;
//							if (d < dist && d < 10) {
//								dist = d;
//								nxt = Point(curr.x - l, k);
//							}
//						}
//					}
//				}
//			}
//			if (f == 1) {
//				sequence.push_back(nxt);
//				sz++;
//			}
//			break;
//		}
//	} while (1);
//	return(sz);
//}
//
//int BackTracePWLine135(int **img, int h, int w, std::vector<Point> &sequence, std::vector<Point> terminals) {
//	int sz = 1;
//	int **copy;
//	copy = new int*[h];
//	for (int i = 0; i < h; i++) {
//		copy[i] = new int[w];
//		for (int j = 0; j < w; j++) {
//			copy[i][j] = (255 - img[i][j]) / 255; //1  - FG, 0 - BG
//		}
//	}
//
//	do {
//		Point nxt;
//		Point curr = sequence[sz - 1];
//		double dist = 9999.9;
//		int flag = 0, intercept = -curr.x + curr.y + 6;
//		for (int k = 0; k < 2; k++) {
//			intercept += k;
//			for (int y = 0; y < h; y++) {
//				int x = y - intercept;
//				if (x < 0 || x >= w)
//					continue;
//				if (copy[y][x] == 1) {
//					double d = sqrt((double(y - curr.y)*(y - curr.y)) + double((x - curr.x)*(x - curr.x)));
//					//std::cout << ".....distance = " << d;
//					if (d <= 12) {
//						flag = 1;
//						if (d < dist) {
//							dist = d;
//							nxt = Point(x, y);
//						}
//					}
//				}
//			}
//		}
//		if (flag == 1) {     //next point available...
//							 sequence.push_back(nxt);
//							 sz++;
//			double slope = atan2(nxt.y - curr.y, nxt.x - curr.x) * 180 / CV_PI;
//
//			//if (abs(slope) >= 80) {    //check for non horizontal...
//			//	sequence.push_back(nxt);
//			//	sz++;
//			//	//std::cout << "\n..." << nxt.y << " , " << nxt.x;
//			//}
//			//else
//			//	flag = 0;
//		}
//		if (flag == 0) {  //find the nearest terminal point...
//			double dist = 9999.9;
//			int f = 0;
//			for (int intercept = -curr.x + curr.y; intercept <= -curr.x + curr.y + 6; intercept++) {
//				for (int y = 0; y < h; y++) {
//					int x = y - intercept;
//					if (x < 0 || x >= w) continue;
//					if (isTerminal(terminals, Point(x, y))) {
//						double d = sqrt((double(y - curr.y)*(y - curr.y)) + double((x - curr.x)*(x - curr.x)));
//						if (d < 5) {
//							f = 1;
//							if (d < dist) {
//								dist = d;
//								nxt = Point(x, y);
//							}
//						}
//					}
//				}
//			}
//			if (f == 1) {
//				sequence.push_back(nxt);
//				sz++;
//			}
//			break;
//		}
//	} while (1);
//	return(sz);
//}
//
//
//double stdDev(double *x, int sz) {
//	double X = 0, sd = 0;
//	for (int i = 0; i < sz ; i++)
//		X += x[i];
//	X /= sz;
//
//	for (int i = 0; i < sz; i++) {
//		sd += ((x[i] - X)*(x[i] - X));
//	}
//	sd = sqrt(sd / sz);
//	return(sd);
//}
//
//int isUniform0(std::vector<Point> sequence, int sz) {
//	double *Y = new double[sz];
//	for (int i = 0; i < sz; i++)
//		Y[i] = sequence[i].y;
//	std::cout << "\nStd deviation in y-direction = " << stdDev(Y, sz);
//	return(1);
//	if (stdDev(Y, sz) < 8)
//		return(1);
//	else
//		return(0);
//}
//
//double avgSlope(std::vector<Point> sequence, int sz) {
//	double *slope = new double[sz - 1];
//	double meanSlope = 0.0;
//	
//	for (int i = 0; i < sz - 1; i++) {
//		slope[i] = atan2(sequence[i + 1].y - sequence[i].y, sequence[i + 1].x - sequence[i].x) * 180 / CV_PI;
//		meanSlope += slope[i];
//	}
//	meanSlope /= (sz - 1);
//	std::cout << "\naverage = " << meanSlope;
//	return(meanSlope);
//}
//
//double regSlope(double corr, std::vector<Point> sequence, int sz) {
//	double *X = new double[sz];
//	double *Y = new double[sz];
//	for (int i = 0; i < sz; i++) {
//		X[i] = sequence[i].x;
//		Y[i] = sequence[i].y;
//	}
//	double sigma_X = stdDev(X, sz);
//	double sigma_Y = stdDev(Y, sz);
//	if (sigma_X != 0) {
//		return(corr*sigma_Y / sigma_X);
//	}
//	else
//		return(0);
//}
//
//double cor_coeff(std::vector<Point> sequence, int sz) {
//	double *X = new double[sz];
//	double *Y = new double[sz];
//	double M_xy = 0, M_x = 0, M_y = 0, M_X2 = 0, M_Y2 = 0;
//	for (int i = 0; i < sz; i++) {
//		X[i] = sequence[i].x;
//		Y[i] = sequence[i].y;
//	}
//	for (int i = 0; i < sz; i++) {
//		M_xy += X[i] * Y[i];
//		M_x += X[i];
//		M_y += Y[i];
//		M_X2 += (X[i] * X[i]);
//		M_Y2 += (Y[i] * Y[i]);
//	}
//	double coeff = (sz*M_xy - M_x*M_y) / (sqrt(sz*M_X2 - M_x*M_x)*sqrt(sz*M_Y2 - M_y*M_y));
//	std::cout << "\nCorrelation coefficient = " << coeff;
//	return(coeff);
//}
//
//void traceLine(int **img, int h, int w, std::vector<Point> &sequence) {
//	
//	int sz = 2;
//	Point next, refrnce = sequence[sz - 1];
//	
//	while (refrnce.y + 5 < w && sz<10) {
//		std::vector<Point> locus;
//		for (int j = refrnce.y; j <= refrnce.y + 5; j++) {
//			int i1 = refrnce.x - abs(sqrt(25 - (refrnce.y - j)*(refrnce.y - j)));
//			int i2 = refrnce.x + abs(sqrt(25 - (refrnce.y - j)*(refrnce.y - j)));
//			if (i1 >= 0 && i1 < h) {
//				int nbd1 = checkNbd(img, h, w, Point(j, i1));
//				if (img[i1][j] == 0 && nbd1 == 8)
//					locus.push_back(Point(j, i1));
//			}
//			if (i2 >= 0 && i2 < h) {
//				int nbd2 = checkNbd(img, h, w, Point(j, i2));
//				if (img[i2][j] == 0 && nbd2 == 8)
//					locus.push_back(Point(j, i2));
//			}
//			
//		}
//		if (int(locus.size()) == 0) 
//			break;
//		std::cout << "\nhas points";
//		std::cout << "\nPossible nbrs = " << int(locus.size());
//		if (sz >= 2) {
//			Point prev_vec = refrnce - sequence[sz - 2];
//			Point curr_vec;
//			double dot_prod, max_dot_prod = 999;
//			for (int i = 0; i < int(locus.size()); i++) {
//				std::cout << "..." << locus[i].x << "," << locus[i].y;
//				curr_vec = locus[i] - refrnce;
//				dot_prod = acos((prev_vec.dot(curr_vec)) / (norm(curr_vec)*norm(prev_vec))) * 180 / CV_PI;
//				if (dot_prod > 90) dot_prod = 180 - dot_prod;
//				std::cout << "** dot prod = " << dot_prod;
//				if (abs(dot_prod) < max_dot_prod) {
//					max_dot_prod = dot_prod;
//					next = locus[i];
//				}
//			}
//			sequence.push_back(next);
//			sz++;
//		}	
//		for (int i = 0; i < sz; i++)
//			std::cout << "->" << sequence[i].y << " , " << sequence[i].x;
//		refrnce = sequence[sz - 1];
//	}
//	std::cout << "\ntrace length = " << sz;
//}
//
//int FlowDirection(int **img, int h, int w, Point terminal) {
//	/* 0 : no flow
//	   1 : only forward
//	   2 : only backward
//	   3 : both direction      */
//	int flow, sumf = 0, sumb = 0;
//	for (int k = 0; k < h; k++) {
//		if(terminal.x + 6 < w)
//			sumf += (255 - img[k][terminal.x + 6]) / 255;
//		if(terminal.x - 6 >= 0)
//			sumb += (255 - img[k][terminal.x - 6]) / 255;
//	}
//	std::cout << "\nfwd sum = " << sumf << " , bwd sum = " << sumb;
//	if (sumf > 0 && sumb == 0)
//		flow = 1;
//	else if (sumf == 0 && sumb > 0)
//		flow = 2;
//	else if (sumf > 0 && sumb > 0)
//		flow = 3;
//	else if (sumf == 0 && sumb == 0)
//		flow = 0;
//	return(flow);
//}
//
//int FlowDirection135(int **img, int h, int w, Point terminal) {
//	/* 0 : no flow
//	1 : only forward
//	2 : only backward
//	3 : both direction      */
//	int flow, sumf = 0, sumb = 0, interceptf = -terminal.x + terminal.y - 6, interceptb = -terminal.x + terminal.y + 6;
//	for (int k = 0; k < 2; k++) {
//		interceptb += k;
//		interceptf += k;
//		for (int y = 0; y < h; y++) {
//			int xf = y - interceptf;
//			int xb = y - interceptb;
//			if (xf >= 0 && xf < w)
//				sumf += (255 - img[y][xf]) / 255;
//			if (xb >= 0 && xb < w)
//				sumb += (255 - img[y][xb]) / 255;			
//		}
//	}
//	std::cout << "\nfwd sum = " << sumf << " , bwd sum = " << sumb;
//
//	return(1);
//	if (sumf > 0 && sumb == 0)
//		flow = 1;
//	else if (sumf == 0 && sumb > 0)
//		flow = 2;
//	else if (sumf > 0 && sumb > 0)
//		flow = 3;
//	else if (sumf == 0 && sumb == 0)
//		flow = 0;
//	return(flow);
//}
//
//
//void tracing0(Mat img, Mat scancol) {
//	int **orgimg, **thin_img, **label, **box;
//	Mat colorimg, colororg;
//
//	orgimg = new int*[r];
//	thin_img = new int*[r];
//	for (int i = 0; i < r; i++) {
//		orgimg[i] = new int[c];
//		thin_img[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			orgimg[i][j] = 0;
//			thin_img[i][j] = 255;
//		}
//	}
//	mat2arr(img, orgimg);
//	colororg.create(r, c, CV_8UC1);
//	cvtColor(img, colororg, CV_GRAY2BGR);
//
//	//Label components...
//	label = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label[i] = new int[c];
//		for (int j = 0; j < c; j++)
//			label[i][j] = 0;
//	}
//	labelNum = labelling(orgimg, label, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	//Trace between the terminal points...
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//		int **roi, **thin_roi;
//		int istart = box[i][4];
//		int jstart = box[i][2];
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//
//		//thin the ROI...
//		roi = new int*[h];
//		thin_roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			thin_roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				thin_roi[k][l] = 255;
//				if (label[k + istart][l + jstart] == i)
//					roi[k][l] = 0;
//				else
//					roi[k][l] = 255;
//			}
//		}
//		thinCC(roi, h, w, thin_roi);
//
//		for (int k = istart; k < h + istart; k++) {
//			for (int l = jstart; l < w + jstart; l++) {
//				if (thin_roi[k - istart][l - jstart] == 0)
//					thin_img[k][l] = 150;
//			}
//		}
//
//		//Determine the terminal points...
//		std::vector<Point> terminals, next;
//		int sz;
//		sz = findTerminals(thin_roi, h, w, terminals, next);
//		std::cout << "\n Label " << i << " has " << sz << " terminal points";
//
//		for (int k = 0; k < sz; k++) {
//			thin_img[terminals[k].y + istart][terminals[k].x + jstart] = 0;
//			//thin_img[next[k].y + istart][next[k].x + jstart] = 1;
//			//line(colororg, terminals[k] + Point(jstart,istart), next[k] + Point(jstart, istart), Scalar(255, 0, 255), 1);
//		}
//
//		//Trace the lines from each terminal point...
//		for (int k = 0; k < int(terminals.size()); k++) {
//			int flow = FlowDirection(thin_roi, h, w, terminals[k]);
//			std::cout << "\nFlow direction = " << flow;
//			if (flow == 1) {
//				std::cout << "\nTrace in forward direction...";
//				std::vector<Point> sequence;
//				sequence.push_back(terminals[k]);
//				int sz = TracePWLine0(thin_roi, h, w, sequence, terminals);
//				std::cout << "\nlength = " << sz;
//				if (sz < 5)
//					continue;
//				double avg = avgSlope(sequence, sz);
//				std::cout << "\naverage slope = " << abs(avg);
//				if (isUniform0(sequence, sz) == 0)
//					continue;
//				if (abs(avg) >= 15)
//					continue;
//				for (int l = 0; l < sz - 1; l++) {
//					thin_img[sequence[l].y + istart][sequence[l].x + jstart] = 0;
//					line(colororg, sequence[l] + Point(jstart, istart), sequence[l + 1] + Point(jstart, istart), Scalar(0, 0, 255), 3);
//					line(scancol, sequence[l] + Point(jstart, istart), sequence[l + 1] + Point(jstart, istart), Scalar(0, 0, 255), 3);
//				}
//				thin_img[sequence[sz - 1].y + istart][sequence[sz - 1].x + jstart] = 1;
//			}
//			else if (flow == 2) {
//				std::cout << "\nTrace in backward direction...";
//				std::vector<Point> sequence;
//				sequence.push_back(terminals[k]);
//				int sz = BackTracePWLine0(thin_roi, h, w, sequence, terminals);
//				std::reverse(sequence.begin(), sequence.end());
//				if (sz < 5)
//					continue;
//				double avg = avgSlope(sequence, sz);
//				std::cout << "\naverage slope = " << abs(avg);
//				if (isUniform0(sequence, sz) == 0)
//					continue;
//				if (abs(avg) >= 15)
//					continue;
//				for (int l = 0; l < sz - 1; l++) {
//					thin_img[sequence[l].y + istart][sequence[l].x + jstart] = 0;
//					line(colororg, sequence[l] + Point(jstart, istart), sequence[l + 1] + Point(jstart, istart), Scalar(0, 0, 255), 3);
//					line(scancol, sequence[l] + Point(jstart, istart), sequence[l + 1] + Point(jstart, istart), Scalar(0, 0, 255), 3);
//				}
//				thin_img[sequence[sz - 1].y + istart][sequence[sz - 1].x + jstart] = 1;
//			}
//			
//		}
//
//
//		for (int k = 0; k < h; k++) {
//			delete[] roi[k];
//			delete[] thin_roi[k];
//		}
//		delete[] roi;
//		delete[] thin_roi;
//	}
//
//	colorimg.create(r, c, CV_8UC1);
//	arr2mat(thin_img, colorimg);
//	cvtColor(colorimg, colorimg, CV_GRAY2BGR);
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (thin_img[i][j] > 1) continue;
//			Vec3b color, color1;
//			if (thin_img[i][j] == 0) {
//				color[0] = 255;
//				color[1] = 0;
//				color[2] = 100;
//				colorimg.at<Vec3b>(i, j) = color;
//				colororg.at<Vec3b>(i, j) = color;
//			}
//			if (thin_img[i][j] == 1) {
//				color[0] = 255;
//				color[1] = 0;
//				color[2] = 0;
//				colorimg.at<Vec3b>(i, j) = color;
//				colororg.at<Vec3b>(i, j) = color;
//			}
//
//			/*color1[0] = 0;
//			color1[1] = 0;
//			color1[2] = 255;
//			colororg.at<Vec3b>(i, j) = color1;*/
//		}
//	}
//	imshow("Thin+terminal.jpg", colorimg);
//	//imwrite("Data/ICDAR/Data218/rot0Thin.tif", colorimg);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan002/rot0Thin.tif", colorimg);
//	imshow("terminal.jpg", colororg);
//	//imwrite("Data/ICDAR/Data218/rot0Terminals.tif", colororg);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan002/rot0Terminals.tif", colororg);
//}
//
//void tracing135(Mat img, Mat scancol) {
//	int **orgimg, **thin_img, **label, **box;
//	Mat colorimg , colororg;
//
//	orgimg = new int*[r];
//	thin_img = new int*[r];
//	for (int i = 0; i < r; i++) {
//		orgimg[i] = new int[c];
//		thin_img[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			orgimg[i][j] = 0;
//			thin_img[i][j] = 255;
//		}
//	}
//	mat2arr(img, orgimg);
//	colororg.create(r, c, CV_8UC1);
//	cvtColor(img, colororg, CV_GRAY2BGR);
//
//	//Label components...
//	label = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label[i] = new int[c];
//		for (int j = 0; j < c; j++)
//			label[i][j] = 0;
//	}
//	labelNum = labelling(orgimg, label, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	//Trace between the terminal points...
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//
//		line(colororg, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//		line(colororg, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//		line(colororg, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//		line(colororg, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//		
//
//		int **roi, **thin_roi;
//		int istart = box[i][4];
//		int jstart = box[i][2];
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//
//		//thin the ROI...
//		roi = new int*[h];
//		thin_roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			thin_roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				thin_roi[k][l] = 255;
//				if (label[k + istart][l + jstart] == i)
//					roi[k][l] = 0;
//				else
//					roi[k][l] = 255;
//			}
//		}
//		thinCC(roi, h, w, thin_roi);
//
//		for (int k = istart; k < h + istart; k++) {
//			for (int l = jstart; l < w + jstart; l++) {
//				if (thin_roi[k - istart][l - jstart] == 0)
//					thin_img[k][l] = 150;
//			}
//		}
//
//		//Determine the terminal points...
//		std::vector<Point> terminals , next;
//		int sz;
//		sz = findTerminals(thin_roi, h, w, terminals, next);
//		std::cout << "\n Label " << i << " has " << sz << " terminal points";
//
//		for (int k = 0; k < sz; k++) {
//			thin_img[terminals[k].y + istart][terminals[k].x + jstart] = 0;
//			//thin_img[next[k].y + istart][next[k].x + jstart] = 1;
//			//line(colororg, terminals[k] + Point(jstart,istart), next[k] + Point(jstart, istart), Scalar(255, 0, 255), 1);
//		}
//		
//		//Trace the lines from each terminal point...
//		for (int k = 0; k < int(terminals.size()); k++) {
//			int flow = FlowDirection135(thin_roi, h, w, terminals[k]);
//			std::cout << "\nFlow direction = " << flow;
//			if (flow == 1) {
//				std::cout << "\nTrace in forward direction...";
//				std::vector<Point> sequence;
//				sequence.push_back(terminals[k]);
//				int sz = TracePWLine135(thin_roi, h, w, sequence, terminals);
//				if (sz <= 6)
//					continue;
//				double corr = cor_coeff(sequence, sz);
//				double avg = atan(regSlope(corr, sequence, sz)) * 180 / CV_PI;
//				std::cout << "\naverage slope = " << avg;
//				if (abs(corr) < 0.94) {
//					std::cout << "...No" << corr;
//					continue;
//				}
//				if (!(avg <= -10 && avg >= -80))
//					continue;
//				for (int l = 0; l < sz - 1; l++) {
//					thin_img[sequence[l].y + istart][sequence[l].x + jstart] = 0;
//					line(colororg, sequence[l] + Point(jstart, istart), sequence[l + 1] + Point(jstart, istart), Scalar(0, 0, 255), 3);
//					line(scancol, sequence[l] + Point(jstart, istart), sequence[l + 1] + Point(jstart, istart), Scalar(255, 0, 255), 3);
//				}
//				thin_img[sequence[sz - 1].y + istart][sequence[sz - 1].x + jstart] = 1;
//			}
//			else if (flow == 2) {
//				std::cout << "\nTrace in backward direction...";
//				std::vector<Point> sequence;
//				sequence.push_back(terminals[k]);
//				int sz = BackTracePWLine135(thin_roi, h, w, sequence, terminals);
//				std::reverse(sequence.begin(), sequence.end());
//				if (sz <= 6)
//					continue;
//				double corr = cor_coeff(sequence, sz);
//				double avg = atan(regSlope(corr, sequence, sz)) * 180 / CV_PI;
//				std::cout << "\naverage slope = " << avg;
//				/*if (abs(corr) < 0.94) {
//					std::cout << "...No" << corr;
//					continue;
//				}
//				if (!(avg <= -10 && avg >= -80))
//					continue;*/
//				for (int l = 0; l < sz - 1; l++) {
//					thin_img[sequence[l].y + istart][sequence[l].x + jstart] = 0;
//					line(colororg, sequence[l] + Point(jstart, istart), sequence[l + 1] + Point(jstart, istart), Scalar(0, 0, 255), 3);
//					line(scancol, sequence[l] + Point(jstart, istart), sequence[l + 1] + Point(jstart, istart), Scalar(255, 0, 255), 3);
//				}
//				thin_img[sequence[sz - 1].y + istart][sequence[sz - 1].x + jstart] = 1;
//			}
//			
//		}
//		
//
//		for (int k = 0; k < h; k++) {
//			delete[] roi[k];
//			delete[] thin_roi[k];
//		}
//		delete[] roi;
//		delete[] thin_roi;
//	}
//
//	colorimg.create(r, c, CV_8UC1);	
//	arr2mat(thin_img, colorimg);
//	cvtColor(colorimg, colorimg, CV_GRAY2BGR);
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (thin_img[i][j] > 1) continue;
//			Vec3b color , color1;
//			if (thin_img[i][j] == 0) {
//				color[0] = 255;
//				color[1] = 0;
//				color[2] = 100;
//				colorimg.at<Vec3b>(i, j) = color;
//				colororg.at<Vec3b>(i, j) = color;
//			}
//			if (thin_img[i][j] == 1) {
//				color[0] = 255;
//				color[1] = 0;
//				color[2] = 0;
//				colorimg.at<Vec3b>(i, j) = color;
//				colororg.at<Vec3b>(i, j) = color;
//			}
//			
//			/*color1[0] = 0;
//			color1[1] = 0;
//			color1[2] = 255;
//			colororg.at<Vec3b>(i, j) = color1;*/
//		}
//	}
//	imshow("Thin+terminal.jpg", colorimg);
//	//imwrite("Data/ICDAR/Data218/rot135Thin.tif", colorimg);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan002/rot135Thin.tif", colorimg);
//	imshow("terminal.jpg", colororg);
//	//imwrite("Data/ICDAR/Data218/rot135Terminals.tif", colororg);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan002/rot135Terminals.tif", colororg);
//	imshow("Primary direction.tif", scancol);
//}
//
//void tracing45(Mat img, Mat scancol) {
//	int **orgimg, **thin_img, **label, **box;
//	Mat colorimg, colororg;
//
//	orgimg = new int*[r];
//	thin_img = new int*[r];
//	for (int i = 0; i < r; i++) {
//		orgimg[i] = new int[c];
//		thin_img[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			orgimg[i][j] = 0;
//			thin_img[i][j] = 255;
//		}
//	}
//	mat2arr(img, orgimg);
//	colororg.create(r, c, CV_8UC1);
//	cvtColor(img, colororg, CV_GRAY2BGR);
//
//	//Label components...
//	label = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label[i] = new int[c];
//		for (int j = 0; j < c; j++)
//			label[i][j] = 0;
//	}
//	labelNum = labelling(orgimg, label, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	//Trace between the terminal points...
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//		int **roi, **thin_roi;
//		int istart = box[i][4];
//		int jstart = box[i][2];
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//
//		//thin the ROI...
//		roi = new int*[h];
//		thin_roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			thin_roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				thin_roi[k][l] = 255;
//				if (label[k + istart][l + jstart] == i)
//					roi[k][l] = 0;
//				else
//					roi[k][l] = 255;
//			}
//		}
//		thinCC(roi, h, w, thin_roi);
//
//		for (int k = istart; k < h + istart; k++) {
//			for (int l = jstart; l < w + jstart; l++) {
//				if (thin_roi[k - istart][l - jstart] == 0)
//					thin_img[k][l] = 150;
//			}
//		}
//
//		//Determine the terminal points...
//		std::vector<Point> terminals, next;
//		int sz;
//		sz = findTerminals(thin_roi, h, w, terminals, next);
//		std::cout << "\n Label " << i << " has " << sz << " terminal points";
//
//		for (int k = 0; k < sz; k++) {
//			thin_img[terminals[k].y + istart][terminals[k].x + jstart] = 0;
//			//thin_img[next[k].y + istart][next[k].x + jstart] = 1;
//			//line(colororg, terminals[k] + Point(jstart,istart), next[k] + Point(jstart, istart), Scalar(255, 0, 255), 1);
//		}
//
//		//Trace the lines from each terminal point...
//		for (int k = 0; k < int(terminals.size()); k++) {
//			std::vector<Point> sequence;
//			sequence.push_back(terminals[k]);
//			int sz = TracePWLine45(thin_roi, h, w, sequence, terminals);
//			if (sz <= 5)
//				continue;
//			double corr = cor_coeff(sequence, sz);
//			double avg = atan(regSlope(corr, sequence, sz)) * 180 / CV_PI;
//			std::cout << "\naverage slope = " << avg;
//			if (corr < 0.92) {
//				std::cout << "...No" << corr;
//				continue;
//			}
//			if (!(avg >= 10 && avg <= 75) && !(avg <= -100 && avg >= -175))
//				continue;
//			for (int l = 0; l < sz - 1; l++) {
//				thin_img[sequence[l].y + istart][sequence[l].x + jstart] = 0;
//				line(colororg, sequence[l] + Point(jstart, istart), sequence[l + 1] + Point(jstart, istart), Scalar(0, 0, 255), 3);
//				line(scancol, sequence[l] + Point(jstart, istart), sequence[l + 1] + Point(jstart, istart), Scalar(0, 140, 255), 3);
//			}
//			thin_img[sequence[sz - 1].y + istart][sequence[sz - 1].x + jstart] = 1;
//		}
//
//
//		for (int k = 0; k < h; k++) {
//			delete[] roi[k];
//			delete[] thin_roi[k];
//		}
//		delete[] roi;
//		delete[] thin_roi;
//	}
//
//	colorimg.create(r, c, CV_8UC1);
//	arr2mat(thin_img, colorimg);
//	cvtColor(colorimg, colorimg, CV_GRAY2BGR);
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (thin_img[i][j] > 1) continue;
//			Vec3b color, color1;
//			if (thin_img[i][j] == 0) {
//				color[0] = 255;
//				color[1] = 0;
//				color[2] = 100;
//				colorimg.at<Vec3b>(i, j) = color;
//				colororg.at<Vec3b>(i, j) = color;
//			}
//			if (thin_img[i][j] == 1) {
//				color[0] = 255;
//				color[1] = 0;
//				color[2] = 0;
//				colorimg.at<Vec3b>(i, j) = color;
//				colororg.at<Vec3b>(i, j) = color;
//			}
//
//			/*color1[0] = 0;
//			color1[1] = 0;
//			color1[2] = 255;
//			colororg.at<Vec3b>(i, j) = color1;*/
//		}
//	}
//	imshow("Thin+terminal.jpg", colorimg);
//	//imwrite("Data/ICDAR/Data218/rot45Thin.tif", colorimg);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan002/rot45Thin.tif", colorimg);
//	imshow("terminal.jpg", colororg);
//	//imwrite("Data/ICDAR/Data218/rot45Terminals.tif", colororg);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan002/rot45Terminals.tif", colororg);
//	imshow("Primary direction.tif", scancol);
//}
//
////Trace the OIC component-wise with grid size = 2...
//void sort(int *arr, int n) {
//	for (int i = 0; i < n - 1; i++) {
//		for (int j = 0; j < n - i - 1; j++) {
//			if (arr[j] > arr[j + 1]) {
//				int temp = arr[j];
//				arr[j] = arr[j + 1];
//				arr[j + 1] = temp;
//			}
//		}					
//	}
//}
//void sort(std::vector<Point> ends[2], int n) {
//	for (int i = 0; i < n - 1; i++) {
//		for (int j = 0; j < n - i - 1; j++) {
//			if (ends[0][j].x > ends[0][j + 1].x) {
//				Point temp = ends[0][j];
//				ends[0][j] = ends[0][j + 1];
//				ends[0][j + 1] = temp;
//
//				temp = ends[1][j];
//				ends[1][j] = ends[1][j + 1];
//				ends[1][j + 1] = temp;
//			}
//		}
//	}
//}
//
//
//double mean(int *arr, int sz) {
//	double M = 0.0;
//	for (int i = 0; i < sz; i++) {
//		//std::cout << "\n" << arr[i]<<" M = "<<M;
//		M += (double)arr[i];
//	}
//	//std::cout << "\n mean = " << M ;
//	return(M / (double)sz);
//}
//double median(int *arr, int sz) {
//	sort(arr, sz);
//	for (int i = 0; i < sz - 2; i++)
//		std::cout << "\n" << arr[i];
//	if (sz % 2 == 1)
//		return(arr[int((sz - 1) / 2)]);
//	else
//		return((arr[int(sz / 2)] + arr[int(sz / 2) - 1]) / 2.0);
//}
//void medianFitLine(std::vector<int> X, std::vector<int> Y, Point& L, Point& M, Point& R, Mat colororg) {	
//	int n = int(X.size());
//	std::cout << "\nNo. of points = " << n;
//	int x_r, x_m, x_l, y_r, y_m, y_l;
//	int sz2 = floor(n / 2.0);
//	int sz3 = floor(n / 3.0);
//
//	//Left median
//	int *arrX = new int[sz3];
//	int *arrY = new int[sz2];
//	for (int k = 0; k < sz3; k++) {
//		arrX[k] = X[k];
//	}
//	for (int k = 0; k < sz2; k++) {
//		arrY[k] = Y[k];
//	}
//	x_l = median(arrX, sz3);
//	y_l = median(arrY, sz2);
//	delete[] arrX;
//	delete[] arrY;
//	L = Point(x_l, y_l);
//	std::cout << "\nLeft = " << L;
//
//	//Right median
//	int sz2_new = sz2 + (n % 2);
//	int sz3_new = sz3 + (n % 3);
//	arrX = new int[sz3_new];
//	arrY = new int[sz2_new];
//	for (int k = 0; k < sz3_new; k++) {
//		arrX[k] = X[k + 2 * sz3];
//	}
//	for (int k = 0; k < sz2_new; k++) {
//		arrY[k] = Y[k + sz2];
//	}
//	x_r = median(arrX, sz3_new);
//	y_r = median(arrY, sz2_new);
//	delete[] arrX;
//	delete[] arrY;
//	R = Point(x_r, y_r);
//	std::cout << "\nRight = " << R;
//
//	//Middle median
//	arrX = new int[sz3];
//	//arrY = new int[sz3];
//	for (int k = 0; k < sz3; k++) {
//		arrX[k] = X[k + sz3];
//		//arrY[k] = Y[k + sz3];
//	}
//	x_m = median(arrX, sz3);
//	//y_m = mean(arrY, sz3);
//	y_m = int(y_l + y_r) / 2;
//	delete[] arrX;
//	//delete[] arrY;
//	M = Point(x_m, y_m);
//
//
//}
//
//int checkBlack(int roi[3][5], int h, int w) {
//	for (int i = 0; i < h; i++) {
//		for (int j = 0; j< w; j++) {
//			if (roi[i][j] == 0)
//				return(1);
//		}
//	}
//	return(0);
//}
//
//void matraLine(Mat img, Mat colororg, Mat colorimg, int **orgimg, int **box, int **label, std::vector<int> OIClabel, std::vector<std::vector<Point>> OIC, double **marginData, Point **baselineData) {
//	int g = grid;
//	for (int i = 0; i < int(OIC.size()); i++) {
//		int lbl = OIClabel[i];
//
//		int **roi;
//		int istart = box[lbl][4];
//		int jstart = box[lbl][2];
//		int h = box[lbl][5] - box[lbl][4] + 1;
//		int w = box[lbl][3] - box[lbl][2] + 1;
//
//		roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				if (label[k + istart][l + jstart] == lbl)
//					roi[k][l] = 0;
//				else
//					roi[k][l] = 255;
//			}
//		}
//
//		//list the points on horizontal segments of OIC with length >= 7 (at least 3 grid points)...
//		int perimeter = 0;
//		std::vector<Point> OIC0corner[2];
//		Point origin = Point(jstart, istart);;
//		for (int k = 0; k<int(OIC[i].size()); k++) {
//			Point prev, curr, p;
//			curr = OIC[i][k] - Point(2*g, 2*g);
//			p = curr + origin;
//			if (k == 0)
//				prev = OIC[i][int(OIC[i].size()) - 1] - Point(2 * g, 2 * g);
//			else
//				prev = OIC[i][k - 1] - Point(2 * g, 2 * g);
//
//			if (prev.y != curr.y || abs(prev.x - curr.x) <= 2 * g) //not horizontal segment...
//				continue;
//			//line(colororg, curr + origin, prev + origin, Scalar(0, 255, 255), 1, 8);
//
//			//if (curr.y + 4 > ceil(h * 0.5)) //segment not in the upper strip...
//			//	continue;
//			//line(colororg, curr + origin, prev + origin, Scalar(255, 120, 0), 1, 8);
//			
//			int flag = 0;
//			for (int y = prev.y + 1; y <= prev.y + g; y++) {
//				for (int x = curr.x; x <= prev.x; x++) {
//					if (x < 0 || x >= w || y < 0 || y >= h)
//						continue;
//					if (roi[y][x] == 0) {
//						flag = 1;
//						break;
//					}
//				}
//				if (flag) break;
//			}
//			if (flag) {
//				OIC0corner[0].push_back(prev + origin);
//				OIC0corner[1].push_back(curr + origin); 
//				//line(colororg, curr + origin, prev + origin, Scalar(255, 0, 120), 1, 8);
//				perimeter += abs(curr.x - prev.x);
//			}
//			
//		}		
//
//		if (int(OIC0corner[0].size()) == 0 || perimeter <= 30) continue;
//		//List of line segments...
//		std::cout << "\nLabel " << lbl << "=>";
//		for (int k = 0; k<int(OIC0corner[0].size()); k++) {
//			std::cout << "\nfrom " << OIC0corner[0][k] << " to " << OIC0corner[1][k];
//		}
//
//		//Select only the base points...
//		int *Xbin = new int[w + 2*g];
//		for (int k = 0; k < w; k++)
//			Xbin[k] = -1;
//		std::vector<int> X, Y;
//		sort(OIC0corner, int(OIC0corner[0].size()));
//		for (int k = 0; k<int(OIC0corner[0].size()); k++) {
//			for (int l = OIC0corner[1][k].x; l <= OIC0corner[0][k].x; l++) {
//				if (Xbin[l - origin.x + g] == -1)
//					Xbin[l - origin.x + g] = OIC0corner[1][k].y;
//				else if (Xbin[l - origin.x + g] > OIC0corner[1][k].y)
//					Xbin[l - origin.x + g] = OIC0corner[1][k].y;
//			}
//		}
//		for (int k = 0; k < w; k++) {
//			if (Xbin[k] == -1) continue;
//			X.push_back(k + origin.x - g);
//			Y.push_back(Xbin[k]);
//		}
//
//		std::cout << "\nList of " << int(X.size()) << " points:";
//		for (int k = 0; k<int(X.size()); k++) {
//			Vec3b color = { 0, 255, 0 };
//			colororg.at<Vec3b>(Y[k], X[k]) = color;
//			std::cout << "\nx = " << X[k] << " , y = " << Y[k];
//		}
//
//		delete Xbin;
//		//Find the median-median line through these segments...
//		double intercept = 0;
//		double slope = 0;
//		Point L, M, R;
//		medianFitLine(X, Y, L, M, R, colororg);
//
//		Vec3b color1 = { 0, 255, 0 };
//		Vec3b color2 = { 0, 255, 120 };
//		Vec3b color3 = { 255, 120, 0 };
//		colororg.at<Vec3b>(L.y , L.x) = color1;
//		colororg.at<Vec3b>(M.y , M.x) = color2;
//		colororg.at<Vec3b>(R.y , R.x) = color3;
//		std::cout << "\n" << L << " , " << M << " , " << R;
//		//line(colororg, L, R, Scalar(0, 120, 255), 1, 8);
//		slope = (R.y - L.y) / double(R.x - L.x);
//		//intercept = ((L.y - slope*L.x) + (M.y - slope*M.x) + (R.y - slope*R.x)) / 3.0;
//		intercept = (L.y - slope*L.x);
//		std::cout << "\nSlope = " << slope << " intercept = " << intercept;
//		
//		if (abs(slope) > tan(35 * CV_PI / 180)) continue;
//
//		marginData[i][1] = intercept;
//		marginData[i][2] = slope;
//		Point p1(jstart, slope*jstart + intercept);
//		Point p2(jstart + w - 1, slope*(jstart + w - 1) + intercept);
//		baselineData[i][0] = L;
//		baselineData[i][1] = R;
//		std::cout << "\nfrom " << p1 << " to " << p2;
//		line(colororg, p1, p2, Scalar(255, 0, 255), 2.5, 8);
//		/*line(colorimg, p1, p2, Scalar(255, 0, 255), 2, 8);*/
//
//		//delete...
//		for (int k = 0; k < h; k++)
//			delete[] roi[k];
//		delete[] roi;
//		OIC0corner[0].empty();
//		OIC0corner[1].empty();
//	}
//
//	imshow("OIC.tif", colororg);
//	imwrite("Data/ICDAR/Data218/OIC.tif", colororg);
//}
//
////void baseLine1(Mat img, Mat colororg, Mat colorimg, int **orgimg, int **box, int **label, std::vector<int> OIClabel, std::vector<std::vector<Point>> OIC, double **marginData) {
////
////	for (int i = 0; i < int(OIC.size()); i++) {
////		int lbl = OIClabel[i];
////
////		int **roi;
////		int istart = box[lbl][4];
////		int jstart = box[lbl][2];
////		int h = box[lbl][5] - box[lbl][4] + 1;
////		int w = box[lbl][3] - box[lbl][2] + 1;
////
////		roi = new int*[h];
////		for (int k = 0; k < h; k++) {
////			roi[k] = new int[w];
////			for (int l = 0; l < w; l++) {
////				if (label[k + istart][l + jstart] == lbl)
////					roi[k][l] = 0;
////				else
////					roi[k][l] = 255;
////			}
////		}
////
////		//list the points on horizontal segments of OIC with length >= 7 (at least 3 grid points)...
////		std::vector<Point> OIC0corner[2];
////		Point origin = Point(jstart, istart);;
////		for (int k = 0; k<int(OIC[i].size()); k++) {
////			Point prev, curr, p;
////			curr = OIC[i][k] - Point(4, 4);
////			p = curr + origin;
////			if (k == 0)
////				prev = OIC[i][int(OIC[i].size()) - 1] - Point(4, 4);
////			else
////				prev = OIC[i][k - 1] - Point(4, 4);
////
////			if (prev.y != curr.y || abs(prev.x-curr.x) <= 4) //not horizontal segment...
////				continue;
////			//line(colororg, curr + origin, prev + origin, Scalar(0, 255, 255), 1, 8);
////
////			if (curr.y - 4 < ceil(h * 0.5)) //segment not in the bottom strip...
////				continue;
////			//line(colororg, curr + origin, prev + origin, Scalar(255, 120, 0), 1, 8);
////
////			int flag = 0;
////			for (int y = prev.y - 2; y < prev.y; y++) {
////				for (int x = prev.x; x <= curr.x; x++) {
////					if (x < 0 || x >= w || y < 0 || y >= h)
////						continue;
////					if (roi[y][x] == 0) {
////						flag = 1;
////						break;
////					}
////				}
////				if (flag) break;
////			}
////			if (flag) {
////				OIC0corner[0].push_back(curr + origin);
////				OIC0corner[1].push_back(prev + origin);
////				//line(colororg, curr + origin, prev + origin, Scalar(255, 0, 120), 1, 8);
////			}
////
////		}
////		if (int(OIC0corner[0].size()) == 0) continue;
////		//List of line segments...
////		std::cout << "\nLabel " << lbl << "=>";
////		for (int k = 0; k<int(OIC0corner[0].size()); k++) {
////			std::cout << "\nfrom " << OIC0corner[0][k] << " to " << OIC0corner[1][k];
////		}
////
////		//Find the median-median line through these segments...
////		double intercept = 0;
////		double slope = 0;
////		Point L, M, R;
////		medianFitLine(OIC0corner, L, M, R, w, origin, colororg);
////
////		Vec3b color1 = { 0, 255, 0 };
////		Vec3b color2 = { 0, 255, 120 };
////		Vec3b color3 = { 255, 120, 0 };
////		colororg.at<Vec3b>(L.y, L.x) = color1;
////		colororg.at<Vec3b>(M.y, M.x) = color2;
////		colororg.at<Vec3b>(R.y, R.x) = color3;
////		std::cout << "\n" << L << " , " << M << " , " << R;
////		//line(colororg, L, R, Scalar(0, 120, 255), 1, 8);
////		slope = (R.y - L.y) / double(R.x - L.x);
////		intercept = ((L.y - slope*L.x) + (M.y - slope*M.x) + (R.y - slope*R.x)) / 3.0;
////		std::cout << "\nSlope = " << slope << " intercept = " << intercept;
////
////		if (abs(slope) > tan(20 * CV_PI / 180)) continue;
////
////		Point p1(jstart, slope*jstart + intercept);
////		Point p2(jstart + w - 1, slope*(jstart + w - 1) + intercept);
////		std::cout << "\nfrom " << p1 << " to " << p2;
////		line(colororg, p1, p2, Scalar(255, 0, 120), 1, 8);
////		line(colorimg, p1, p2, Scalar(255, 90, 120), 2, 8);
////
////		//delete...
////		for (int k = 0; k < h; k++)
////			delete[] roi[k];
////		delete[] roi;
////		OIC0corner[0].empty();
////		OIC0corner[1].empty();
////	}
////
////	imshow("OIC.tif", colororg);
////	imwrite("Data/ICDAR/Data218/OIC.tif", colororg);
////}
////
//void baseLine(Mat img, Mat colororg, Mat colorimg, int **orgimg, int **box, int **label, std::vector<int> OIClabel, std::vector<std::vector<Point>> OIC, double **marginData, Point **baselineData) {
//	int g = grid;
//	for (int i = 0; i < int(OIC.size()); i++) {
//		int lbl = OIClabel[i];
//
//		int **roi;
//		int istart = box[lbl][4];
//		int jstart = box[lbl][2];
//		int h = box[lbl][5] - box[lbl][4] + 1;
//		int w = box[lbl][3] - box[lbl][2] + 1;
//
//		roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				if (label[k + istart][l + jstart] == lbl)
//					roi[k][l] = 0;
//				else
//					roi[k][l] = 255;
//			}
//		}
//
//		//------------------------------------------------------------------------------------------------------
//		//list the points on horizontal segments of OIC with length >= 7 (at least 3 grid points)...
//		int perimeter = 0;
//		std::vector<Point> OIC0corner[2];
//		Point origin = Point(jstart, istart);;
//		for (int k = 0; k<int(OIC[i].size()); k++) {
//			Point prev, curr, p;
//			curr = OIC[i][k] - Point(2 * g, 2 * g);
//			p = curr + origin;
//			if (k == 0)
//				prev = OIC[i][int(OIC[i].size()) - 1] - Point(2 * g, 2 * g);
//			else
//				prev = OIC[i][k - 1] - Point(2*g, 2*g);
//
//			if (prev.y != curr.y || abs(prev.x - curr.x) <= 2*g) //not horizontal segment...
//				continue;
//			//line(colororg, curr + origin, prev + origin, Scalar(0, 255, 255), 1, 8);
//			
//			//if (curr.y - 4 < ceil(h * 0.45)) //segment not in the bottom strip...
//			//	continue;
//			////line(colororg, curr + origin, prev + origin, Scalar(255, 120, 0), 1, 8);
//
//			int flag = 0;
//			for (int y = prev.y - g; y < prev.y; y++) {
//				for (int x = prev.x; x <= curr.x; x++) {
//					if (x < 0 || x >= w || y < 0 || y >= h)
//						continue;
//					if (roi[y][x] == 0) {
//						flag = 1;
//						break;
//					}
//				}
//				if (flag) break;
//			}
//			if (flag) {
//				OIC0corner[0].push_back(curr + origin);
//				OIC0corner[1].push_back(prev + origin);
//				line(colororg, curr + origin, prev + origin, Scalar(255, 0, 120), 1, 8);
//				perimeter += abs(curr.x - prev.x);
//			}
//
//		}
//		if (int(OIC0corner[0].size()) == 0 || perimeter <= 30) continue;
//		//List of line segments...
//		std::cout << "\nLabel " << lbl << "=>";
//		for (int k = 0; k<int(OIC0corner[0].size()); k++) {
//			std::cout << "\nfrom " << OIC0corner[0][k] << " to " << OIC0corner[1][k];
//		}
//
//		//------------------------------------------------------------------------------------------------------
//		//Select only the base points...
//		int *Xbin = new int[w + 2*g];
//		for (int k = 0; k < w; k++)
//			Xbin[k] = -1;
//		std::vector<int> X, Y;
//		sort(OIC0corner, int(OIC0corner[0].size()));
//		for (int k = 0; k<int(OIC0corner[0].size()); k++) {
//			for (int l = OIC0corner[1][k].x; l <= OIC0corner[0][k].x; l++) {
//				if (Xbin[l - origin.x + g] == -1)
//					Xbin[l - origin.x + g] = OIC0corner[1][k].y;
//				else if (Xbin[l - origin.x + g] < OIC0corner[1][k].y)
//					Xbin[l - origin.x + g] = OIC0corner[1][k].y;
//			}
//		}
//		for (int k = 0; k < w; k++) {
//			if (Xbin[k] == -1) continue;
//			X.push_back(k + origin.x - g);
//			Y.push_back(Xbin[k]);
//		}
//
//		std::cout << "\nList of " << int(X.size()) << " points:";
//		for (int k = 0; k<int(X.size()); k++) {
//			Vec3b color = { 0, 255, 0 };
//			colororg.at<Vec3b>(Y[k], X[k]) = color;
//			std::cout << "\nx = " << X[k] << " , y = " << Y[k];
//		}
//
//		delete Xbin;
//		
//		//------------------------------------------------------------------------------------------------------
//		//Find the median-median line through these segments...
//		double intercept = 0;
//		double slope = 0;
//		Point L, M, R;
//		medianFitLine(X, Y, L, M, R, colororg);
//
//		Vec3b color1 = { 0, 255, 0 };
//		Vec3b color2 = { 0, 120, 255 };
//		Vec3b color3 = { 255, 120, 0 };
//		colororg.at<Vec3b>(L.y, L.x) = color1;
//		colororg.at<Vec3b>(M.y, M.x) = color2;
//		colororg.at<Vec3b>(R.y, R.x) = color3;
//		std::cout << "\n" << L << " , " << M << " , " << R;
//		//line(colororg, L, R, Scalar(0, 120, 255), 1, 8);
//		slope = (R.y - L.y) / double(R.x - L.x);
//		//intercept = ((L.y - slope*L.x) + (M.y - slope*M.x) + (R.y - slope*R.x)) / 3.0;
//		intercept = (L.y - slope*L.x) ;
//		std::cout << "\nSlope = " << slope << " intercept = " << intercept;
//
//		if (abs(slope) > tan(45 * CV_PI / 180)) continue;
//
//		marginData[i][3] = intercept;
//		marginData[i][4] = slope;
//
//		Point p1(jstart, slope*jstart + intercept);
//		Point p2(jstart + w - 1, slope*(jstart + w - 1) + intercept);
//		baselineData[i][2] = L;
//		baselineData[i][3] = R;
//		std::cout << "\nfrom " << p1 << " to " << p2;
//		line(colororg, p1, p2, Scalar(255, 0, 120), 2.5, 8);
//		//line(colorimg, p1, p2, Scalar(255, 90, 120), 2, 8);
//
//		//------------------------------------------------------------------------------------------------------
//		//delete memory...
//		for (int k = 0; k < h; k++)
//			delete[] roi[k];
//		delete[] roi;
//		OIC0corner[0].empty();
//		OIC0corner[1].empty();
//		X.empty();
//		Y.empty();
//	}
//
//	imshow("OIC.tif", colororg);
//	imwrite("Data/ICDAR/Data218/OIC.tif", colororg);
//}
//
//std::vector<int> BWratio(int **img, int *box, Point p1, Point p2, Point p3, Point p4, double intercept1, double intercept2, double slope) {
//	int countBlack = 0, countWhite = 0;
//	for (int j = box[2]; j <= box[3]; j++) {
//		int imin = slope*j + intercept1;
//		int imax = slope*j + intercept2;
//		for (int i = imin; i <= imax; i++) {
//			if (i<0 || i>=r || j<0 || j>=c) continue;
//			if (img[i][j] == 0)
//				countBlack++;
//			else
//				countWhite++;
//		}		
//	}
//	std::vector<int> BWcounts;
//	BWcounts.push_back(countBlack);
//	BWcounts.push_back(countWhite);
//	return(BWcounts);
//	//return(double(countWhite) / double(countBlack + countWhite));
//}
//
//
//void displayMargins(int **img, Mat colorimg, double **marginData, Point **baselineData, int row, int col, int **box, int **label) {
//	for (int i = 0; i < row; i++) {
//		int lbl = int(marginData[i][0]);
//
//		int istart = box[lbl][4];
//		int jstart = box[lbl][2];
//		int h = box[lbl][5] - box[lbl][4] + 1;
//		int w = box[lbl][3] - box[lbl][2] + 1;
//
//		std::cout << "\nLabel " << marginData[i][0] << " : " << marginData[i][1] << " , " << marginData[i][2] << " ; " << marginData[i][3] << " , " << marginData[i][4];
//		double m1 = marginData[i][2];
//		double m2 = marginData[i][4];
//		double intercept1 = marginData[i][1];
//		double intercept2 = marginData[i][3];
//		
//		double theta1 = atan(m1) * 180 / CV_PI;
//		double theta2 = atan(m2) * 180 / CV_PI;
//		std::cout << " ; angles = " << theta1 << " , " << theta2;
//
//		if (intercept1 == 0 && intercept2 == 0) continue;
//		//if (abs(theta1 - theta2) > 15) continue;
//		
//		Point p1, p2, p3, p4, p5, p6, p7, p8, p9, p10;
//		double ratioU = -1.0, ratioL = -1.0, ratioU1 = -1.0, ratioL1 = -1.0, interceptL = 0.0, interceptR = 0.0;
//		std::vector<int> BWcounts;
//
//		//Draw upper baseline and the parallel lower baseline...
//		if (intercept1 != 0.0) {
//			p1 = Point(jstart, m1*jstart + intercept1);
//			p2 = Point(jstart + w - 1, m1*(jstart + w - 1) + intercept1);
//			interceptL = (baselineData[i][2].y - m1*baselineData[i][2].x);
//			p3 = Point(jstart, m1*jstart + interceptL);
//			p4 = Point(jstart + w - 1, m1*(jstart + w - 1) + interceptL);
//
//			BWcounts = BWratio(img, box[lbl], p1, p2, p3, p4, intercept1, interceptL, m1);
//			ratioU = double(BWcounts[1]) / BWcounts[0];  //#W/#B
//			//ratioU = double(BWcounts[1]) / ( BWcounts[0] + BWcounts[1]);          //#W/#B+W
//
//			interceptR = (baselineData[i][3].y - m1*baselineData[i][3].x);
//			p9 = Point(jstart, m1*jstart + interceptR);
//			p10 = Point(jstart + w - 1, m1*(jstart + w - 1) + interceptR);
//
//			BWcounts = BWratio(img, box[lbl], p1, p2, p9, p10, intercept1, interceptR, m1);
//			ratioU1 = double(BWcounts[1]) / BWcounts[0];  //#W/#B
//			//ratioU1 = double(BWcounts[1]) / (BWcounts[0] + BWcounts[1]);          //#W/#B+W
//
//			if (ratioU1 < ratioU) {
//				p3 = p9;
//				p4 = p10;
//				ratioU = ratioU1;
//			}
//			std::cout << "\n%age white(upper) = " << ratioU;
//		}
//		
//		//Draw lower baseline and the parallel upper baseline...
//		if (intercept2 != 0.0) {
//			p5 = Point(jstart, m2*jstart + intercept2);
//			p6 = Point(jstart + w - 1, m2*(jstart + w - 1) + intercept2);
//			interceptL = (baselineData[i][0].y - m2*baselineData[i][0].x);
//			p7 = Point(jstart, m2*jstart + interceptL);
//			p8 = Point(jstart + w - 1, m2*(jstart + w - 1) + interceptL);
//
//			BWcounts = BWratio(img, box[lbl], p7, p8, p5, p6, interceptL, intercept2, m2);
//			ratioL = double(BWcounts[1]) / BWcounts[0];  //#W/#B
//			//ratioL = double(BWcounts[1]) / (BWcounts[0] + BWcounts[1]);  //#W/#B
//
//			interceptR = (baselineData[i][1].y - m2*baselineData[i][1].x);
//			p9 = Point(jstart, m2*jstart + interceptR);
//			p10 = Point(jstart + w - 1, m2*(jstart + w - 1) + interceptR);
//
//			BWcounts = BWratio(img, box[lbl], p9, p10, p5, p6, interceptR, intercept2, m2);
//			ratioL1 = double(BWcounts[1]) / BWcounts[0];  //#W/#B
//			//ratioL1 = double(BWcounts[1]) / (BWcounts[0] + BWcounts[1]);  //#W/#B
//
//			if (ratioL1 < ratioL) {
//				p7 = p9;
//				p8 = p10;
//				ratioL = ratioL1;
//			}
//			std::cout << "\n%age white(lower) = " << ratioL;
//		}
//		
//		//Reset the baseline defining points data...
//		baselineData[i][0] = Point(0, 0);
//		baselineData[i][1] = Point(0, 0);
//		baselineData[i][2] = Point(0, 0);
//		baselineData[i][3] = Point(0, 0);
//
//		if (ratioU <= ratioL && intercept1 != 0 ) {
//			baselineData[i][0] = p1;
//			baselineData[i][1] = p2;
//			line(colorimg, p1, p2, Scalar(255, 0, 255), 2, 8);
//			if (intercept2 != 0.0) {
//				baselineData[i][2] = p3;
//				baselineData[i][3] = p4;
//				line(colorimg, p3, p4, Scalar(255, 0, 255), 2, 8);
//			}
//		}
//		else if (ratioL < 0) {
//			baselineData[i][0] = p1;
//			baselineData[i][1] = p2;
//			line(colorimg, p1, p2, Scalar(255, 0, 255), 2, 8);
//			if (m1 < 0) {
//				p3 = Point(jstart, istart + h - 1);
//				p4.x = jstart + w - 1;
//				p4.y = p3.y + m1*w;
//			}
//			else {
//				p4 = Point(jstart + w - 1, istart + h - 1);
//				p3.x = jstart;
//				p3.y = p4.y - m1*w;
//			}
//			baselineData[i][2] = p3;
//			baselineData[i][3] = p4;
//			line(colorimg, p3, p4, Scalar(255, 0, 255), 2, 8);
//		}
//
//		if (ratioU >= ratioL && intercept2 != 0) {
//			baselineData[i][2] = p5;
//			baselineData[i][3] = p6;
//			line(colorimg, p5, p6, Scalar(0, 120, 255), 2, 8);
//			if (intercept1 != 0.0) {
//				baselineData[i][0] = p7;
//				baselineData[i][1] = p8;
//				line(colorimg, p7, p8, Scalar(0, 120, 255), 2, 8);
//			}
//		}
//		else if (ratioU < 0) {
//			baselineData[i][2] = p5;
//			baselineData[i][3] = p6;
//			line(colorimg, p5, p6, Scalar(0, 120, 255), 2, 8);
//			if (m2 >= 0) {
//				p7 = Point(jstart, istart);
//				p8.x = jstart + w - 1;
//				p8.y = p7.y + m2*w;
//			}
//			else {
//				p8 = Point(jstart + w - 1, istart);
//				p7.x = jstart;
//				p7.y = p8.y - m2*w;
//			}
//			baselineData[i][0] = p7;
//			baselineData[i][1] = p8;
//			line(colorimg, p7, p8, Scalar(0, 120, 255), 2, 8);
//		}
//	}
//
//	imshow("bounds.tif", colorimg);
//	imwrite("Data/ICDAR/Data218/bounds.tif", colorimg);
//}
//
//void traceOIC(Mat img , int **orgimg, int h, int w, std::vector<Point> &sequence) {
//	Mat colimg;
//	//int grid = 3;
//	int **S, **oiclabel;
//	
//	int g = grid;
//	h = h - ((h - 1) % g);
//	w = w - ((w - 1) % g);
//
//	int R = h + 4 * g, C = w + 4 * g;
//	std::cout << "\nPadded dimension = " << R << " x " << C;
//	S = new int*[R];
//	oiclabel = new int*[R];
//	for (int i = 0; i < R; i++) {
//		S[i] = new int[C];
//		oiclabel[i] = new int[C];
//		for (int j = 0; j < C; j++) {
//			oiclabel[i][j] = 0;
//			if (i < 2*g || j < 2*g || i >= h + 2*g || j >= w + 2*g)
//				S[i][j] = 255;
//			else
//				S[i][j] = orgimg[i - 2*g][j - 2*g];
//		}
//	}
//	colimg.create(R, C, CV_8UC1);
//	arr2mat(S, colimg, R, C);
//	cvtColor(colimg, colimg, CV_GRAY2BGR);
//	imwrite("Data/ICDAR/Data218/OIC.tif", colimg);
//
//	Point p0;
//	int flag = 0;
//	for (int i = g; i < R - g - 1; i++) {
//		for (int j = g; j < C - g - 1; j++) {
//			if (S[i][j] == 0) {
//				std::cout << "\n First black = " << i << " , " << j;
//				p0 = Point(j, i);
//				flag = 1;
//				break;
//			}
//		}
//		if (flag)
//			break;
//	}
//	std::cout << "\n First black = " << p0.y << " , " << p0.x;
//
//	sequence.push_back(p0);
//	MakeOIC(S, sequence, g);
//
//	//for (int i = 0; i<int(sequence.size()) ; i++) {
//	//	oiclabel[sequence[i].y][sequence[i].x] = 1;
//	//	Vec3b color = { 0,0,255 };
//	//	colimg.at<Vec3b>(sequence[i].y, sequence[i].x) = color;
//	//	if (i<int(sequence.size()) - 1)
//	//		line(colimg, sequence[i], sequence[i + 1], Scalar(0, 0, 255), 1, 8);
//	//	else
//	//		line(colimg, sequence[i], sequence[0], Scalar(0, 0, 255), 1, 8);
//	//}
//	//imshow("OIC.tif", colimg);
//	//imwrite("Data/ICDAR/Data218/OIC.tif", colimg);
//}
//
//void ConnCompOIC(Mat img, Mat scanimg) {
//	int **orgimg, **label, **box, **oiclabel, **binscanimg;
//	double **marginData;
//	Point **baselineData;
//	Mat colorimg, colororg, colorbounds;
//	std::vector<std::vector<Point>> OIC;
//	std::vector<int> oiclabels;
//
//	orgimg = new int*[r];
//	binscanimg = new int*[r];
//	for (int i = 0; i < r; i++) {
//		orgimg[i] = new int[c];
//		binscanimg[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			orgimg[i][j] = 0;
//			binscanimg[i][j] = 0;
//		}
//	}
//	mat2arr(img, orgimg);
//	mat2arr(img, binscanimg);
//	//NICK(binscanimg, r, c, binscanimg);
//	colororg.create(r, c, CV_8UC1);
//	cvtColor(img, colororg, CV_GRAY2BGR);
//	colorbounds.create(r, c, CV_8UC1);
//	cvtColor(scanimg, colorbounds, CV_GRAY2BGR);
//
//	//Label components...
//	label = new int*[r];
//	oiclabel = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label[i] = new int[c];
//		oiclabel[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			label[i][j] = 0;
//			oiclabel[i][j] = 0;
//		}
//	}
//	labelNum = labelling(orgimg, label, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	int g = grid;
//	//Find OICs component-wise...
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//
//		Mat ROI;
//		int **roi, **roiOIC;
//		int istart = box[i][4];
//		int jstart = box[i][2];
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//
//		if (h < 10 && w < 10) continue;
//
//		roi = new int*[h];
//		roiOIC = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			roiOIC[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				roiOIC[k][l] = 0;
//				if (label[k + istart][l + jstart] == i)
//					roi[k][l] = 0;
//				else
//					roi[k][l] = 255;
//			}
//		}
//		ROI.create(h, w, CV_8UC1);
//		arr2mat(roi, ROI, h, w);
//
//		std::cout << "\nLabel = " << i;
//		std::vector<Point> sequenceOIC;
//		traceOIC(ROI, roi, h, w, sequenceOIC);
//
//		OIC.push_back(sequenceOIC);
//		oiclabels.push_back(i);
//
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				oiclabel[k + istart][l + jstart] = i*roiOIC[k][l];
//			}
//		}
//		//delete...
//		for (int k = 0; k < h; k++) {
//			delete[] roi[k];
//		}
//		delete[] roi;
//	}
//
//	//Display OICs...
//	for (int i = 0; i<int(OIC.size()); i++) {
//		for (int j = 0; j<int(OIC[i].size()); j++) {
//			Vec3b color = { 0,120,255 };
//			int istart = box[oiclabels[i]][4];
//			int jstart = box[oiclabels[i]][2];
//			Point origin = Point(jstart - 2*g, istart - 2*g);
//			Point p = OIC[i][j] + origin;
//			if (p.x >= 0 && p.x < c && p.y >= 0 && p.y < r) {
//				colororg.at<Vec3b>(p.y, p.x) = color;
//				if (j<int(OIC[i].size()) - 1)
//					line(colororg, OIC[i][j] + origin, OIC[i][j + 1] + origin, Scalar(0, 0, 255), 1, 8);
//				else
//					line(colororg, OIC[i][j] + origin, OIC[i][0] + origin, Scalar(0, 0, 255), 1, 8);
//			}
//		}
//	}
//	std::cout << "\nhere!";
//	imshow("OIC.tif", colororg);
//	imwrite("Data/ICDAR/Data218/OIC.tif", colororg);
//
//	marginData = new double*[int(oiclabels.size())];
//	for (int i = 0; i<int(oiclabels.size()); i++) {
//		marginData[i] = new double[5];
//		marginData[i][0] = oiclabels[i];
//		for (int j = 1; j < 5; j++) {
//			marginData[i][j] = 0.0;
//		}
//	}
//	baselineData = new Point*[int(oiclabels.size())];
//	for (int i = 0; i<int(oiclabels.size()); i++) {
//		baselineData[i] = new Point[4];
//		for (int j = 1; j < 4; j++) {
//			baselineData[i][j] = Point(0, 0);
//		}
//	}
//	colorimg.create(r, c, CV_8UC1);
//	cvtColor(scanimg, colorimg, CV_GRAY2BGR);
//	matraLine(img, colororg, colorimg, orgimg, box, label, oiclabels, OIC, marginData, baselineData);
//	baseLine(img, colororg, colorimg, orgimg, box, label, oiclabels, OIC, marginData, baselineData);
//
//	//Display upper and base lines...
//	displayMargins(binscanimg, colorbounds, marginData, baselineData, int(OIC.size()), 5, box, label);
//
//	//Delete...
//	for (int i = 0; i < r; i++) {
//		delete[] orgimg[i];
//	}
//	delete[] orgimg;
//}
//
//
//int main(int argc, char** argv) {
//
//	//-----------------------------Read image-----------------------
//	//Mat scanimg = imread("Data/ICDAR/Data218/218.tif", IMREAD_GRAYSCALE);
//	/*Mat scanimg = imread("Data/ICDAR/Data218/scan218.jpg", IMREAD_GRAYSCALE);
//	Mat img0 = imread("Data/ICDAR/Data218/Gauss0.tif", IMREAD_GRAYSCALE);
//	Mat img45 = imread("Data/ICDAR/Data218/Gauss45.tif", IMREAD_GRAYSCALE);
//	Mat img135 = imread("Data/ICDAR/Data218/Gauss135.tif", IMREAD_GRAYSCALE);*/
//
//	Mat scanimg = imread("Data/MyScans/Scanned Data/Scanned Data/Scan002/scan002.tif", IMREAD_GRAYSCALE);
//	Mat img0 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan002/Gauss0.tif", IMREAD_GRAYSCALE);
//	Mat img45 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan002/Gauss45.tif", IMREAD_GRAYSCALE);
//	Mat img135 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan002/Gauss135.tif", IMREAD_GRAYSCALE);
//	Mat scancol;
//
//	imshow("Image", scanimg);
//	std::cout << "Image dimensions:" << scanimg.rows << "x" << scanimg.cols << endl;
//
//	r = scanimg.rows;
//	c = scanimg.cols;
//	
//	scancol.create(r, c, CV_8UC1);
//	cvtColor(scanimg, scancol, CV_GRAY2BGR);
//	tracing0(img0, scancol);
//	tracing135(img135, scancol);
//	tracing45(img45, scancol);
//
//	imshow("Primarydirection.tif", scancol);
//	//imwrite("Data/ICDAR/Data218/GlobalDirections.tif", scancol);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan002/GlobalDirections.tif", scancol);
//
//	//ConnCompOIC(img0, scanimg);
//	//invCC(scanimg, 1);
//
//	std::cout << "\nDone!\n";
//	waitKey(0);
//	system("pause");
//	return 0;
//
//}