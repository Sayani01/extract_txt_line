#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\highgui\highgui.hpp>


#include <stdio.h>
#include <conio.h>
#include <iostream>
#include <vector>
#include <string>

#include "NICK.h"
#include "Otsu.h"
#include "RegLine.h"
#include "ConnectedComponent.h"
#include "BoundingBox.h"
#include "FPTA.h"
#include "GHPTA.h"

using namespace std;
using namespace cv;
namespace cv {
	using std::vector;
};

extern int labelNum = 0, r = 0, c = 0;
extern double AH = 0;

int grid = 2;

void mat2arr(Mat A, int **arr) {
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			arr[i][j] = A.at<uchar>(i, j);
		}
	}
}
void arr2mat(int **arr, Mat A) {
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			A.at<uchar>(i, j) = arr[i][j];
		}
	}
}
//Overloads...
void mat2arr(Mat A, int **arr, int row, int col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			arr[i][j] = A.at<uchar>(i, j);
		}
	}
}
void arr2mat(int **arr, Mat A, int row, int col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			A.at<uchar>(i, j) = arr[i][j];
		}
	}
}
double median(std::vector<int>  arr, int sz) {
	//sort(arr, sz);
	/*for (int i = 0; i < sz - 2; i++)
	std::cout << "\n" << arr[i];*/
	if (sz % 2 == 1)
		return(arr[int((sz - 1) / 2)]);
	else
		return((arr[int(sz / 2)] + arr[int(sz / 2) - 1]) / 2.0);
}

double mean(int *arr, int sz) {
	double M = 0.0;
	for (int i = 0; i < sz; i++) {
		//std::cout << "\n" << arr[i]<<" M = "<<M;
		M += (double)arr[i];
	}
	//std::cout << "\n mean = " << M ;
	return(M / (double)sz);
}

double mean(std::vector<int> arr, int sz) {
	double M = 0.0;
	for (int i = 0; i < sz; i++) {
		//std::cout << "\n" << arr[i]<<" M = "<<M;
		M += (double)arr[i];
	}
	//std::cout << "\n mean = " << M ;
	return(M / (double)sz);
}

void erode(int **img, int h = r, int w = c) {
	int n = 3;
	int **copy = new int*[r];
	for (int i = 0; i < r; i++) {
		copy[i] = new int[c];
		for (int j = 0; j < c; j++)
			copy[i][j] = img[i][j];
	}
	//int struc[7][7] = { { 0,0,0,0,0,0,0 },{ 1,1,1,1,1,1,1 } ,{ 1,1,1,1,1,1,1 } ,{ 1,1,1,1,1,1,1 },{ 1,1,1,1,1,1,1 },{ 1,1,1,1,1,1,1 },{ 0,0,0,0,0,0,0 } };
	//int struc[5][5] = { {0,0,0,0,0}, { 1,1,1,1,1 },{ 1,1,1,1,1 },{ 1,1,1,1,1 },{ 0,0,0,0,0 } };
	int struc[3][3] = { { 1,1,1 },{ 1,1,1 },{ 1,1,1 } };
	for (int i = (n - 1) / 2; i < h - (n - 1) / 2; i++) {
		for (int j = (n - 1) / 2; j < w - (n - 1) / 2; j++) {
			if (copy[i][j] == 0) continue; //text pixel...
			int sum = 0;
			for (int k = -(n - 1) / 2; k <= (n - 1) / 2; k++) {
				for (int l = -(n - 1) / 2; l <= (n - 1) / 2; l++) {
					if (copy[i + k][j + l] == 0 && struc[(n - 1) / 2 + k][(n - 1) / 2 + l] == 1) sum++; //black nbr...
				}
			}
			if (sum != 0) img[i][j] = 0;
		}
	}
}

void dilate(int **img, int h = r, int w = c) {
	int n = 3;
	int **copy = new int*[r];
	for (int i = 0; i < r; i++) {
		copy[i] = new int[c];
		for (int j = 0; j < c; j++)
			copy[i][j] = img[i][j] / 255;
	}
	int struc[3][3] = { { 1,1,1 },{ 1,1,1 },{ 1,1,1 } };
	for (int i = 1; i < h - 1; i++) {
		for (int j = 1; j < w - 1; j++) {
			if (copy[i][j] == 1) continue;
			int sum = 0;
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					sum += copy[i - 1 + k][j - 1 + l] * struc[k][l];
				}
			}
			if (sum != 0) img[i][j] = 255;
		}
	}
}


/***********************************************************************************************************
											COMBINE BOUNDING BOX
************************************************************************************************************/
//bckground=white ; foreground=black
int checkNbd(int **img, int h, int w, Point ref) {
	int sum = 0;
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			if (ref.x + j < 0 || ref.x + j >= w || ref.y + i < 0 || ref.y + i >= h || (i == 0 && j == 0))
				continue;
			else {
				if (img[ref.y + i][ref.x + j] == 0)
					sum++;
			}
		}
	}
	return(sum);
}

Point findCentroid(int **img, int h, int w) {
	std::vector<Point> textpixels;
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			if (img[i][j] == 0) {
				textpixels.push_back(Point(j, i));
			}
		}
	}
	Point centroid = Point(0, 0);
	double sx = 0, sy = 0;
	for (int i = 0; i<int(textpixels.size()); i++) {
		sx += textpixels[i].x;
		sy += textpixels[i].y;
	}
	centroid.x = int(sx / int(textpixels.size()));
	centroid.y = int(sy / int(textpixels.size()));

	return(centroid);

}

Point findMedianCentroid(int **img, int h, int w) {
	std::vector<int> X, Y;
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			if (img[i][j] == 0) {
				Y.push_back(i);
			}
		}
	}
	for (int j = 0; j < w; j++) {
		for (int i = 0; i < h; i++) {
			if (img[i][j] == 0) {
				X.push_back(j);
			}
		}
	}
	Point centroid = Point(median(X, int(X.size())), median(Y, int(Y.size())));
	return(centroid);
}


int findTerminals(int **img, int h, int w, std::vector<Point> &terminal) {
	int **copy, sz = 1;
	copy = new int*[h + 2];
	for (int i = 0; i < h + 2; i++) {
		copy[i] = new int[w + 2];
		for (int j = 0; j < w + 2; j++) {
			if (i == 0 || j == 0 || i == h + 1 || j == w + 1)
				copy[i][j] = 0;
			else
				copy[i][j] = (255 - img[i - 1][j - 1]) / 255;
		}
	}
	for (int i = 1; i <= h; i++) {
		for (int j = 1; j <= w; j++) {
			if (copy[i][j] == 0) continue;
			int sum = checkNbd(copy, h+2, w+2, Point(j, i));
			if (sum == 7) {
				terminal.push_back(Point(j - 1, i - 1));
				sz++;
			}
		}
	}
	return(sz);
}

int isTerminal(std::vector<Point> terminal, Point p) {
	for (int i = 0; i<int(terminal.size()); i++) {
		if (terminal[i] == p)
			return(1);
	}
	return(0);
}

int countCrossings(int **roi, int h, int w, Point p1, Point p2, std::vector<int> &whiteRun) {
	int count = 0;
	double slope = atan2(p2.y - p1.y, p2.x - p1.x) * 180 / CV_PI;
	if (slope < 0)
		slope += 180;

	std::cout << "\nSlope = " << slope;

	if (p1.y == p2.y) {  //0 degrees...
		int istart = p1.x < p2.x ? p1.x : p2.x;
		int istop = p1.x > p2.x ? p1.x : p2.x;
		int white = 0, prev = 255;
		int crossing;
		for (int x = istart; x <= istop; x++) {
			crossing = checkNbd(roi, h, w, Point(x, p1.y));
			if (crossing && prev == 255) { //white to black...
				count++;
				prev = 0;
				if (white != 0) {
					whiteRun.push_back(white);
					white = 0;
				}
			}
			else if (!crossing) {  //black to white...
				white++;
				prev = 255;
			}
		}
		if (white != 0)
			whiteRun.push_back(white);
	}
	else  if (p1.x == p2.x) { //90 degree...
		int jstart = p1.y < p2.y ? p1.y : p2.y;
		int jstop = p1.y > p2.y ? p1.y : p2.y;
		int white = 0, prev = 255;
		int crossing;
		for (int y = jstart; y <= jstop; y++) {
			crossing = checkNbd(roi, h, w, Point(p1.x, y));
			if (crossing && prev == 255) { //white to black...
				count++;
				prev = 0;
				if (white != 0) {
					whiteRun.push_back(white);
					white = 0;
				}
			}
			else if (!crossing) {  //black to white...
				white++;
				prev = 255;
			}
		}
		if (white != 0)
			whiteRun.push_back(white);
	}
	else if (slope >= 45 && slope <= 135) {
		int jstart = p1.y < p2.y ? p1.y : p2.y;
		int jstop = p1.y > p2.y ? p1.y : p2.y;
		int white = 0, prev = 255;
		int crossing;
		for (int y = jstart; y <= jstop; y++) {
			int x = int(p1.x + double((y - p1.y)*(p2.x - p1.x)) / (p2.y - p1.y));
			crossing = checkNbd(roi, h, w, Point(x, y));
			if (crossing && prev == 255) { //white to black...
				count++;
				prev = 0;
				if (white != 0) {
					whiteRun.push_back(white);
					white = 0;
				}
			}
			else if (!crossing) {  //black to white...
				white++;
				prev = 255;
			}
		}
		if (white != 0)
			whiteRun.push_back(white);
	}
	else if ((slope > 0 && slope < 45) || (slope > 135 && slope < 180) ) {
		int istart = p1.x < p2.x ? p1.x : p2.x;
		int istop = p1.x > p2.x ? p1.x : p2.x;
		int white = 0, prev = 255;
		int crossing;
		for (int x = istart; x <= istop; x++) {
			int y = int(p1.y + double((x - p1.x)*(p2.y - p1.y)) / (p2.x - p1.x));
			crossing = checkNbd(roi, h, w, Point(x, y));
			if (crossing && prev == 255) { //white to black...
				count++;
				prev = 0;
				if (white != 0) {
					whiteRun.push_back(white);
					white = 0;
				}
			}
			else if (!crossing) {  //black to white...
				white++;
				prev = 255;
			}
		}
		if (white != 0)
			whiteRun.push_back(white);
	}
	
	return(count);
}

std::vector<std::pair<Point, Point>> combineBoxes(int **thin_img, int **label, int **box, std::vector<std::vector<Point>> all_terminals, std::vector<std::pair<int, Point>> all_centroids, std::vector<std::pair<int, int>> &connectCC) {
	std::vector<std::pair<int, int>> matchedCC;  //store the labels that have been compared...
	  //store the labels to join...
	std::vector<std::pair<Point, Point>> connectCentroid; //store the actual points to connnect...

	int  dim = 12;
	int thresh = 30;
	std::vector<int> nbr_labels;
	for (int i = 1; i < labelNum; i++) {
		if (box[i][2] == -1) continue;
		int **roi;
		int istart = box[i][4];
		int jstart = box[i][2];
		int h = box[i][5] - box[i][4] + 1;
		int w = box[i][3] - box[i][2] + 1;

		//collect the neighbouring labels from overlapping boxes...
		std::cout << "\nLabels:";
		std::cout << "\nBox size : X = " << box[i][2] - dim << " to " << box[i][3] + dim << " , Y = " << box[i][4] - dim << " to " << box[i][5] + dim;
		for (int y = box[i][4] - dim; y <= box[i][5] + dim; y++) {
			for (int x = box[i][2] - dim; x <= box[i][3] + dim; x++) {
				if (y < 0 || y >= r || x < 0 || x >= c)
					continue;
				if (label[y][x] != i && label[y][x] != 0) {
					int flag = 0;
					//std::cout << label[y][x] << "---";
					for (int k = 0; k<int(nbr_labels.size()); k++) {
						if (nbr_labels[k] == label[y][x]) {
							flag = 1;
							break;
						}
					}
					if (!flag)
						nbr_labels.push_back(label[y][x]);
				}
			}
		}

		std::cout << "\n*****************************\nNo. of neighbors for label " << i << " = " << int(nbr_labels.size());

		double length = 0;
		double max_length, max_max_length = -1;
		Point max_max_p1, max_max_p2;
		int  max_lbl = 0, max_run = 0, min_run = 0;
		for (int k = 0; k<int(nbr_labels.size()); k++) {
			int nbr = nbr_labels[k];
			int max, min;
			Point max_p1, max_p2;
			max_length = -1;
			std::cout << "\nComparing with label " << nbr << "...";

			int flag = 0;
			for (int k = 0; k<int(connectCC.size()); k++) {
				if ((connectCC[k].first == i && connectCC[k].second == nbr) || (connectCC[k].first == nbr && connectCC[k].second == i)) {
					flag = 1;
					break;
				}
			}
			if (flag) {
				std::cout << "\nAlready compared...";
				continue;
			}

			//Create the combined CCs...
			int xstart = box[i][2] < box[nbr][2] ? box[i][2] : box[nbr][2];
			int xstop = box[i][3] > box[nbr][3] ? box[i][3] : box[nbr][3];
			int ystart = box[i][4] < box[nbr][4] ? box[i][4] : box[nbr][4];
			int ystop = box[i][5] > box[nbr][5] ? box[i][5] : box[nbr][5];
			int W = xstop - xstart + 1;
			int H = ystop - ystart + 1;
			roi = new int*[H];
			for (int y = 0; y < H; y++) {
				roi[y] = new int[W];
				for (int x = 0; x < W; x++) {
					if (thin_img[ystart + y][xstart + x] != 255 && (label[ystart + y][xstart + x] == i || label[ystart + y][xstart + x] == nbr))
						roi[y][x] = 0;
					else
						roi[y][x] = 255;
				}
			}

			int it_i, it_nbr;
			for (int it = 0; it<int(all_terminals.size()); it++) {
				if (all_terminals[it][0].x == i)
					it_i = it;
				if (all_terminals[it][0].x == nbr)
					it_nbr = it;
			}

			//lines from current's centroid to neighbor's centroid...
			Point p1 = all_centroids[it_i].second - Point(xstart, ystart);
			std::vector<int> whiteRuns;
			Point p2 = all_centroids[it_nbr].second - Point(xstart, ystart);
			int countBW = countCrossings(roi, H, W, p1, p2, whiteRuns);
			std::cout << "\nLine " << all_centroids[it_i].second << " to " << all_centroids[it_nbr].second;
			std::cout << "\nThe label " << nbr << " has " << countBW << " B/W crossings and the sequence of white run-lengths are: ";
			for (int seq = 0; seq<int(whiteRuns.size()); seq++)
				std::cout << whiteRuns[seq] << "...";

			//find min & max run length...
			max = whiteRuns[0];
			min = whiteRuns[0];
			for (int seq = 1; seq<int(whiteRuns.size()); seq++) {
				if (whiteRuns[seq] > max)
					max = whiteRuns[seq];
				else if (whiteRuns[seq] < min)
					min = whiteRuns[seq];
			}
			if (max <= thresh || min <= thresh) {
				length = double(countBW)*min / max;
				std::cout << "\nLength = " << length;
				max_length = length;
				max_run = max;
				min_run = min;
				max_p1 = p1 + Point(xstart, ystart);
				max_p2 = p2 + Point(xstart, ystart);
			}
			
			whiteRuns.clear();

			//lines from current label's centroid to neighbors...
			p1 = all_centroids[it_i].second - Point(xstart, ystart);
			std::vector<int> whiteRuns_i;
			for (int it = 1; it<int(all_terminals[it_nbr].size()); it++) {
				Point p2 = all_terminals[it_nbr][it] - Point(xstart, ystart);
				int countBW = countCrossings(roi, H, W, p1, p2, whiteRuns_i);
				std::cout << "\nLine " << all_centroids[it_i].second << " to " << all_terminals[it_nbr][it];
				std::cout << "\nThe label " << nbr << " has " << countBW << " B/W crossings and the sequence of white run-lengths are: ";
				for (int seq = 0; seq<int(whiteRuns_i.size()); seq++)
					std::cout << whiteRuns_i[seq] << "...";
				
				//find min & max run length...
				max = whiteRuns_i[0];
				min = whiteRuns_i[0];
				for (int seq = 1; seq<int(whiteRuns_i.size()); seq++) {
					if (whiteRuns_i[seq] > max)
						max = whiteRuns_i[seq];
					else if (whiteRuns_i[seq] < min)
						min = whiteRuns_i[seq];
				}
				if (max > thresh || min > thresh)
					continue;
				length = double(countBW)*min / max;
				std::cout << "\nLength = " << length;
				if (max_length < length) {
					max_length = length;
					max_run = max;
					min_run = min;
					max_p1 = p1 + Point(xstart, ystart);
					max_p2 = p2 + Point(xstart, ystart);
				}
				whiteRuns_i.clear();
			}

			//lines from neighbors' centroid to current label...
			p1 = all_centroids[it_nbr].second - Point(xstart, ystart);
			std::vector<int> whiteRuns_nbr;
			for (int it = 1; it<int(all_terminals[it_i].size()); it++) {
				Point p2 = all_terminals[it_i][it] - Point(xstart, ystart);
				int countBW = countCrossings(roi, H, W, p1, p2, whiteRuns_nbr);
				std::cout << "\nLine " << all_centroids[it_nbr].second << " to " << all_terminals[it_i][it];
				std::cout << "\nThe label " << i << " has " << countBW << " B/W crossings and the sequence of white run-lengths are: ";
				for (int seq = 0; seq<int(whiteRuns_nbr.size()); seq++)
					std::cout << whiteRuns_nbr[seq] << "...";
				
				//find min & max run length...
				max = whiteRuns_nbr[0];
				min = whiteRuns_nbr[0];
				for (int seq = 1; seq<int(whiteRuns_nbr.size()); seq++) {
					if (whiteRuns_nbr[seq] > max)
						max = whiteRuns_nbr[seq];
					else if (whiteRuns_nbr[seq] < min)
						min = whiteRuns_nbr[seq];
				}
				if (max > thresh || min > thresh)
					continue;
				length = double(countBW)*min / max;
				std::cout << "\nLength = " << length;
				if (max_length < length) {
					max_length = length;
					max_run = max;
					min_run = min;
					max_p1 = p1 + Point(xstart, ystart);
					max_p2 = p2 + Point(xstart, ystart);
				}
				whiteRuns_nbr.clear();
			}

			//lines from currentlabel's terminal points to neighbours' terminal points...
			for (int  it1= 1; it1<int(all_terminals[it_i].size()); it1++) {
				p1 = all_terminals[it_i][it1] - Point(xstart, ystart);
				std::vector<int> whiteRuns_i;

				for (int it2 = 1; it2<int(all_terminals[it_nbr].size()); it2++) {
					Point p2 = all_terminals[it_nbr][it2] - Point(xstart, ystart);
					int countBW = countCrossings(roi, H, W, p1, p2, whiteRuns_i);
					std::cout << "\nLine " << all_terminals[it_i][it1] << " to " << all_terminals[it_nbr][it2];
					std::cout << "\nThe label " << nbr << " has " << countBW << " B/W crossings and the sequence of white run-lengths are: ";
					for (int seq = 0; seq<int(whiteRuns_i.size()); seq++)
						std::cout << whiteRuns_i[seq] << "...";

					//find min & max run length...
					max = whiteRuns_i[0];
					min = whiteRuns_i[0];
					for (int seq = 1; seq<int(whiteRuns_i.size()); seq++) {
						if (whiteRuns_i[seq] > max)
							max = whiteRuns_i[seq];
						else if (whiteRuns_i[seq] < min)
							min = whiteRuns_i[seq];
					}
					if (max > thresh || min > thresh)
						continue;
					length = double(countBW)*min / max;
					std::cout << "\nLength = " << length;
					if (max_length < length) {
						max_length = length;
						max_run = max;
						min_run = min;
						max_p1 = p1 + Point(xstart, ystart);
						max_p2 = p2 + Point(xstart, ystart);
					}
					whiteRuns_i.clear();
				}
			}

			

			if (max_length > max_max_length) {
				max_max_length = max_length;
				max_lbl = nbr;
				max_max_p1 = max_p1;
				max_max_p2 = max_p2;
			}

		}


		std::cout << "\nLabel " << i << " is near to label " << max_lbl << " with length = " << max_max_length;
		if (max_lbl != 0) {
			connectCC.push_back(pair<int, int>(i, max_lbl));
			connectCentroid.push_back(pair<Point, Point>(max_max_p1, max_max_p2));
		}			

		//empty the vecctor...
		nbr_labels.clear();
	}
	//return(connectCC);
	return(connectCentroid);
}

std::vector<std::pair<Point, Point>> tracing0(Mat img, Mat colorimg, int **thin_img, std::vector<std::pair<Point, Point>> &connectCC) {
	int **orgimg, **label, **box;
	Mat thin_im;

	orgimg = new int*[r];
	for (int i = 0; i < r; i++) {
		orgimg[i] = new int[c];
		for (int j = 0; j < c; j++) {
			orgimg[i][j] = 0;
		}
	}
	mat2arr(img, orgimg);
	
	//Label components...
	label = new int*[r];
	for (int i = 0; i < r; i++) {
		label[i] = new int[c];
		for (int j = 0; j < c; j++)
			label[i][j] = 0;
	}
	labelNum = labelling(orgimg, label, r, c);

	//Bound the components by rectangular & diagonal boxes
	box = new int*[labelNum];
	for (int i = 0; i < labelNum; i++)
		box[i] = new int[6];
	boundingRectangles(box, label, r, c, labelNum);

	//Trace between the terminal points...
	std::vector<std::vector<Point>> all_terminals;
	std::vector<std::pair<int, Point>> all_centroids;
	for (int i = 1; i < labelNum; i++) {
		if (box[i][2] == -1) continue;
		int **roi, **thin_roi;
		int istart = box[i][4];
		int jstart = box[i][2];
		int h = box[i][5] - box[i][4] + 1;
		int w = box[i][3] - box[i][2] + 1;

		//thin the ROI...
		roi = new int*[h];
		thin_roi = new int*[h];
		for (int k = 0; k < h; k++) {
			roi[k] = new int[w];
			thin_roi[k] = new int[w];
			for (int l = 0; l < w; l++) {
				thin_roi[k][l] = 255;
				if (label[k + istart][l + jstart] == i)
					roi[k][l] = 0;
				else
					roi[k][l] = 255;
			}
		}
		Point centroid = findCentroid(roi, h, w) + Point(jstart, istart);
		all_centroids.push_back(std::pair<int, Point>(i, centroid));
		
		thinCC(roi, h, w, thin_roi);

		for (int k = istart; k < h + istart; k++) {
			for (int l = jstart; l < w + jstart; l++) {
				if (thin_roi[k - istart][l - jstart] == 0)
					thin_img[k][l] = 150;
			}
		}

		//Determine the terminal points...
		std::vector<Point> terminals;
		terminals.push_back(Point(i, 0)); //store the label no....
		int sz;
		sz = findTerminals(thin_roi, h, w, terminals);
		std::cout << "\n Label " << i << " has " << sz << " terminal points";

		for (int k = 1; k < sz; k++) {
			terminals[k] += Point(jstart, istart);
			thin_img[terminals[k].y][terminals[k].x] = 1;
		}

		all_terminals.push_back(terminals);

		for (int k = 0; k < h; k++) {
			delete[] roi[k];
			delete[] thin_roi[k];
		}
		delete[] roi;
		delete[] thin_roi;
	}
	
	std::vector<std::pair<int, int>> connectCC1;
	std::vector<std::pair<Point, Point>> connectCentroid = combineBoxes(thin_img, label, box, all_terminals, all_centroids, connectCC1);
	for (int k = 0; k<int(connectCC1.size()); k++) {
		Point p1, p2;
		for (int i = 0; i<int(all_centroids.size()); i++) {
			if (connectCC1[k].first == all_centroids[i].first)
				p1 = all_centroids[i].second;
			if (connectCC1[k].second == all_centroids[i].first)
				p2 = all_centroids[i].second;
		}
		connectCC.push_back(pair<Point, Point>(p1, p2));
	}

	//arr2mat(thin_img, colorimg);
	thin_im.create(r, c, CV_8UC1);

	Vec3b color, color1, color0;
	color[0] = 0;
	color[1] = 0;
	color[2] = 255;

	color1[0] = 0;
	color1[1] = 0;
	color1[2] = 0;

	color0[0] = 255;
	color0[1] = 255;
	color0[2] = 255;

	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			if (thin_img[i][j] == 150) {
				colorimg.at<Vec3b>(i, j) = color1;
				thin_im.at<uchar>(i, j) = 0;
			}			
			else if (thin_img[i][j] == 1) {
				colorimg.at<Vec3b>(i, j) = color;
				thin_im.at<uchar>(i, j) = 0;
				circle(colorimg, Point(j, i), 2, Scalar(0, 0, 255), 2);
			}
			else {
				colorimg.at<Vec3b>(i, j) = color0;
				thin_im.at<uchar>(i, j) = 255;
			}
		}
	}
	imshow("Thinbw.jpg", thin_im);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan002/ThinTeminal.tif", colorimg);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan002/ThinBW.tif", thin_im);

	return(connectCentroid);
}


void fillHoles(int **img, int **label, int **box) {
	std::vector<int> countCC[2];
	int **countPixels = new int *[labelNum];
	for (int i = 0; i < labelNum; i++) {
		countPixels[i] = new int[2];
		countPixels[i][0] = i;
		countPixels[i][1] = 0;
	}

	for (int i = 2; i < labelNum; i++) {
		if (box[i][2] == -1) continue;
		int count = 0;
		for (int k = box[i][4]; k <= box[i][5]; k++) {
			for (int l = box[i][2]; l <= box[i][3]; l++) {
				if (label[k][l] == i)
					count++;
			}
		}
		countPixels[i][1] = count;
		std::cout << "\n Label " << i << " has " << count << " pixels";
		countCC[0].push_back(count);
		countCC[1].push_back((box[i][3] - box[i][2])*(box[i][5] - box[i][4]));
	}
	int *arr = new int[int(countCC[0].size())];
	int *area = new int[int(countCC[0].size())];
	for (int i = 0; i<int(countCC[0].size()); i++) {
		arr[i] = countCC[0][i];
		area[i] = countCC[1][i];
	}
	double medCount = mean(arr, int(countCC[0].size()));
	double medArea = mean(area, int(countCC[0].size()));
	std::cout << "\nThe mean of CC size = " << medCount;
	std::cout << "\nThe mean of CC Area = " << medArea;

	for (int i = 2; i < labelNum; i++) {
		if (box[i][2] == -1) continue;
		if (countPixels[i][1] > 3 * medCount && (box[i][3] - box[i][2])*(box[i][5] - box[i][4]) > 3.5*medArea) continue;
		for (int k = box[i][4]; k <= box[i][5]; k++) {
			for (int l = box[i][2]; l <= box[i][3]; l++) {
				if (label[k][l] == i)
					img[k][l] = 255;
			}
		}
	}
}

void invCCProj(Mat img, int flag = 0) {
	int **orgimg, **label, **box, **binimg, **thin_img;
	Mat colorimg, bin_img, thin_im, near4img;

	orgimg = new int*[r];
	binimg = new int*[r];
	label = new int*[r];
	thin_img = new int*[r];
	for (int i = 0; i < r; i++) {
		orgimg[i] = new int[c];
		binimg[i] = new int[c];
		label[i] = new int[c];
		thin_img[i] = new int[c];
		for (int j = 0; j < c; j++) {
			orgimg[i][j] = 0;
			binimg[i][j] = 0;
			label[i][j] = 0;
			thin_img[i][j] = 255;
		}
	}
	mat2arr(img, orgimg);

	//----------------------------------------------------------------------------------------------------------
	//Binarize if needed...
	if (flag) {
		NICK(orgimg, r, c, binimg);
		labelNum = labelling(binimg, label, r, c);

		//Bound the components by rectangular & diagonal boxes
		box = new int*[labelNum];
		for (int i = 0; i < labelNum; i++)
			box[i] = new int[6];
		boundingRectangles(box, label, r, c, labelNum);

		for (int i = 1; i < labelNum; i++) {
			if (box[i][2] == -1) continue;
			int h = box[i][5] - box[i][4] + 1;
			int w = box[i][3] - box[i][2] + 1;
			int count = 0;
			for (int k = box[i][4]; k <= box[i][5]; k++) {
				for (int l = box[i][2]; l <= box[i][3]; l++) {
					if (binimg[k][l] == 0) count++;
				}
			}
			if (h*w < 30 || count < 100) {
				for (int k = box[i][4]; k <= box[i][5]; k++) {
					for (int l = box[i][2]; l <= box[i][3]; l++) {
						if (label[k][l] == i) {
							label[k][l] = 0;
							binimg[k][l] = 255;
						}
					}
				}
				box[i][2] = -1;
				box[i][3] = -1;
				box[i][4] = -1;
				box[i][5] = -1;
			}
		}

		bin_img.create(r, c, CV_8UC1);
		arr2mat(binimg, bin_img);
		GaussianBlur(bin_img, bin_img, Size(3, 3), 0, 0);
		mat2arr(bin_img, orgimg);
		NICK(orgimg, r, c, binimg);
		orgimg = binimg;
		arr2mat(binimg, bin_img);
		imshow("binaryImg.tif", bin_img);
		imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan002/binImg15.tif", bin_img);

		//----------------------------------------------------------------------------------------------------------
		for (int i = 0; i < labelNum; i++)
			delete[] box[i];
		delete[] box;
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++)
				label[i][j] = 0;
		}
		labelNum = 0;
	}

	//Fill holes...
	erode(orgimg);
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			orgimg[i][j] = 255 - orgimg[i][j];
		}
	}
	//dilate(orgimg);

	labelNum = labelling(orgimg, label, r, c);

	//Bound the components by rectangular & diagonal boxes
	box = new int*[labelNum];
	for (int i = 0; i < labelNum; i++)
		box[i] = new int[6];
	boundingRectangles(box, label, r, c, labelNum);

	//fillHoles(orgimg, label, box);

	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			orgimg[i][j] = 255 - orgimg[i][j];
		}
	}

	//---------------------------------------------------------------------------------------------------------
	//re-label components...
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			label[i][j] = 0;
		}
	}
	for (int i = 0; i < labelNum; i++)
		delete[] box[i];
	delete[] box;
	labelNum = 0;
	std::cout << "\nRe-labelling...";
	labelNum = labelling(orgimg, label, r, c);

	//Bound the components by rectangular & diagonal boxes
	box = new int*[labelNum];
	for (int i = 0; i < labelNum; i++)
		box[i] = new int[6];
	boundingRectangles(box, label, r, c, labelNum);

	//dilate(orgimg);
	arr2mat(orgimg, img);

	imshow("invcc.tif", img);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan002/invcc.tif", img);


	//---------------------------------------------------------------------------------------------------------
	//Exclude small components...
	for (int i = 1; i < labelNum; i++) {
		if (box[i][2] == -1) continue;
		int h = box[i][5] - box[i][4] + 1;
		int w = box[i][3] - box[i][2] + 1;
		int count = 0;
		for (int k = box[i][4]; k <= box[i][5]; k++) {
			for (int l = box[i][2]; l <= box[i][3]; l++) {
				if (orgimg[k][l] == 0) count++;
			}
		}
		if (h*w < 30 || count < 200) {
			for (int k = 0; k < 6; k++)
				box[i][k] = -1;
		}
	}

	//thin and show terminal points...
	colorimg.create(r, c, CV_8UC1);
	cvtColor(img, colorimg, CV_GRAY2BGR);
	std::vector<std::pair<Point, Point>> connectCC;
	std::vector<std::pair<Point, Point>> connectCentroid = tracing0(img, colorimg, thin_img, connectCC);

	for (int i = 0; i < labelNum; i++) {
		if (box[i][2] == -1) continue;
		int **roi;
		int h = box[i][5] - box[i][4] + 1;
		int w = box[i][3] - box[i][2] + 1;

		roi = new int*[h];
		for (int y = 0; y < h; y++) {
			roi[y] = new int[w];
			for (int x = 0; x < w; x++) {
				//if (label[box[i][4] + y][box[i][2] + x] == i && thin_img[box[i][4] + y][box[i][2] + x] != 255)
				if (label[box[i][4] + y][box[i][2] + x] == i )
					roi[y][x] = 0;
				else
					roi[y][x] = 255;
			}
		}
		Point centroid = findCentroid(roi, h, w) + Point(box[i][2], box[i][4]);
		line(colorimg, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(255, 0, 0), 2, 8);
		line(colorimg, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(255, 0, 0), 2, 8);
		line(colorimg, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(255, 0, 0), 2, 8);
		line(colorimg, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(255, 0, 0), 2, 8);
		circle(colorimg, centroid, 2, Scalar(255, 0, 255), 2);

		for (int i = 0; i < h; i++) {
			delete[] roi[i];
		}
		delete[] roi;

	}

	for (int i = 0; i<int(connectCC.size()); i++) {
		//line(colorimg, connectCC[i].first, connectCC[i].second, Scalar(255, 0, 255), 2, 8);
		arrowedLine(colorimg, connectCentroid[i].first, connectCentroid[i].second, Scalar(255, 0, 255), 1.5, 8);
		arrowedLine(colorimg, connectCC[i].first, connectCC[i].second, Scalar(0, 120, 255), 2, 8);
	}

	/*for (int i = 0; i < labelNum; i++) {
		if (box[i][2] == -1) continue;
		line(colorimg, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(255, 0, 0), 2, 8);
		line(colorimg, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(255, 0, 0), 2, 8);
		line(colorimg, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(255, 0, 0), 2, 8);
		line(colorimg, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(255, 0, 0), 2, 8);		
	}*/
	imshow("boundingBox.tif", colorimg);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan002/centroidBoundingBox.tif", colorimg);
	
	//Delete...
	for (int i = 0; i < r; i++) {
		delete[] orgimg[i];
	}
	delete[] orgimg;
}

/***********************************************************************************************************
												DETERMINE PCA
************************************************************************************************************/
void drawAxis(Mat &img, Point p, Point q, Scalar colour, const float scale = 0.2)
{
	double angle;
	double hypotenuse;
	angle = atan2((double)p.y - q.y, (double)p.x - q.x); // angle in radians
	hypotenuse = sqrt((double)(p.y - q.y) * (p.y - q.y) + (p.x - q.x) * (p.x - q.x));
	//    double degrees = angle * 180 / CV_PI; // convert radians to degrees (0-180 range)
	//    cout << "Degrees: " << abs(degrees - 180) << endl; // angle in 0-360 degrees range
	// Here we lengthen the arrow by a factor of scale
	q.x = (int)(p.x - scale * hypotenuse * cos(angle));
	q.y = (int)(p.y - scale * hypotenuse * sin(angle));
	line(img, p, q, colour, 1, CV_AA);
	// create the arrow hooks
	p.x = (int)(q.x + 9 * cos(angle + CV_PI / 4));
	p.y = (int)(q.y + 9 * sin(angle + CV_PI / 4));
	line(img, p, q, colour, 1, CV_AA);
	p.x = (int)(q.x + 9 * cos(angle - CV_PI / 4));
	p.y = (int)(q.y + 9 * sin(angle - CV_PI / 4));
	line(img, p, q, colour, 1, CV_AA);
}

double getOrientation(vector<Point> &pts, Mat img) {
	//Construct a buffer used by the pca analysis
	int sz = static_cast<int>(pts.size());
	Mat data_pts = Mat(sz, 2, CV_64FC1);
	for (int i = 0; i < data_pts.rows; ++i)
	{
		data_pts.at<double>(i, 0) = pts[i].x;
		data_pts.at<double>(i, 1) = pts[i].y;
	}
	//Perform PCA analysis
	PCA pca_analysis(data_pts, Mat(), CV_PCA_DATA_AS_ROW,2);

	//Store the center of the object
	Point cntr = Point(static_cast<int>(pca_analysis.mean.at<double>(0, 0)), static_cast<int>(pca_analysis.mean.at<double>(0, 1)));

	//Store the eigenvalues and eigenvectors
	double eig_vec[2][2], eig_val[2];

	Mat mean = pca_analysis.mean.clone();
	Mat eigenvalues = pca_analysis.eigenvalues.clone();
	Mat eigenvectors = pca_analysis.eigenvectors.clone();
	std::cout << "\n" << eigenvalues;
	std::cout << "\n" << eigenvectors; 
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			eig_vec[i][j] = eigenvectors.at<double>(i, j);
		}
		eig_val[i] = eigenvalues.at<double>(i);
	}

	// Draw the principal components
	Mat col_img = img;
	cvtColor(col_img, col_img, CV_GRAY2BGR);
	circle(col_img, cntr, 3, Scalar(255, 0, 255), 2);
	Point p1 = cntr + 0.02 * Point((eig_vec[0][0] * eig_val[0]),(eig_vec[0][1] * eig_val[0]));
	Point p2 = cntr - 0.02 * Point((eig_vec[1][0] * eig_val[1]), (eig_vec[1][1] * eig_val[1]));
	drawAxis(col_img, cntr, p1, Scalar(0, 255, 0), 1);
	drawAxis(col_img, cntr, p2, Scalar(255, 255, 0), 5);
	double angle = atan2(eig_vec[0][1], eig_vec[0][0]); // orientation in radians
	std::cout << "\nOrientation = " << angle * 180 / CV_PI;
	imshow("orientations.jpg",col_img);
	return angle;
}

void imgPCA(Mat img, int flag = 0) {
	int **orgimg, **label, **box, **binimg, **thin_img;
	Mat colorimg, bin_img, thin_im, near4img;

	orgimg = new int*[r];
	binimg = new int*[r];
	label = new int*[r];
	thin_img = new int*[r];
	for (int i = 0; i < r; i++) {
		orgimg[i] = new int[c];
		binimg[i] = new int[c];
		label[i] = new int[c];
		thin_img[i] = new int[c];
		for (int j = 0; j < c; j++) {
			orgimg[i][j] = 0;
			binimg[i][j] = 0;
			label[i][j] = 0;
			thin_img[i][j] = 255;
		}
	}
	mat2arr(img, orgimg);

	//----------------------------------------------------------------------------------------------------------
	//Binarize if needed...
	if (flag) {
		NICK(orgimg, r, c, binimg);
		labelNum = labelling(binimg, label, r, c);

		//Bound the components by rectangular & diagonal boxes
		box = new int*[labelNum];
		for (int i = 0; i < labelNum; i++)
			box[i] = new int[6];
		boundingRectangles(box, label, r, c, labelNum);

		for (int i = 1; i < labelNum; i++) {
			if (box[i][2] == -1) continue;
			int h = box[i][5] - box[i][4] + 1;
			int w = box[i][3] - box[i][2] + 1;
			int count = 0;
			for (int k = box[i][4]; k <= box[i][5]; k++) {
				for (int l = box[i][2]; l <= box[i][3]; l++) {
					if (binimg[k][l] == 0) count++;
				}
			}
			if (h*w < 30 || count < 100) {
				for (int k = box[i][4]; k <= box[i][5]; k++) {
					for (int l = box[i][2]; l <= box[i][3]; l++) {
						if (label[k][l] == i) {
							label[k][l] = 0;
							binimg[k][l] = 255;
						}
					}
				}
				box[i][2] = -1;
				box[i][3] = -1;
				box[i][4] = -1;
				box[i][5] = -1;
			}
		}

		bin_img.create(r, c, CV_8UC1);
		arr2mat(binimg, bin_img);
		GaussianBlur(bin_img, bin_img, Size(3, 3), 0, 0);
		mat2arr(bin_img, orgimg);
		NICK(orgimg, r, c, binimg);
		orgimg = binimg;
		arr2mat(binimg, bin_img);
		imshow("binaryImg.tif", bin_img);
		imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan002/binImg15.tif", bin_img);
	}

	vector<Point> textpixels;
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			if (binimg[i][j] == 0) {
				textpixels.push_back(Point(j, i));
			}
		}
	}

	getOrientation(textpixels, img);
	imshow("directions", img);
}

int main(int argc, char** argv) {
	
	//-----------------------------Read image-----------------------
	//  Mat scanimg = imread("Data/ICDAR/Data218/218.tif", IMREAD_GRAYSCALE);
	
	Mat scanimg = imread("Data/MyScans/Scanned Data/Scanned Data/Scan002/line.tif", IMREAD_GRAYSCALE);
	
	imshow("Image", scanimg);
	std::cout << "Image dimensions:" << scanimg.rows << "x" << scanimg.cols << endl;
	
	r = scanimg.rows;
	c = scanimg.cols;

	//invCCProj(scanimg, 1);
	imgPCA(scanimg, 1);
	
	std::cout << "\nDone!\n";
	waitKey(0);
	system("pause");
	return 0;
	
}