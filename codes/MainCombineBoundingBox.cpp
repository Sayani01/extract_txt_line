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
	std::sort(arr.begin(), arr.begin() + (sz - 1));
	if (sz % 2 == 1)
		return(arr[int((sz - 1) / 2)]);
	else
		return((arr[int(sz / 2)] + arr[int(sz / 2) - 1]) / 2.0);
}

double median(std::vector<double>  arr, int sz) {
	//sort(arr, sz);
	/*for (int i = 0; i < sz - 2; i++)
	std::cout << "\n" << arr[i];*/

	std::sort(arr.begin(), arr.begin() + (sz - 1));
	if (sz % 2 == 1)
		return(arr[int((sz - 1) / 2)]);
	else
		return((arr[int(sz / 2)] + arr[int(sz / 2) - 1]) / 2.0);
}

Point median(std::vector<Point>  p, int sz) {
	std::vector<int> X, Y;
	for (int i = 0; i < sz; i++) {
		X.push_back(p[i].x);
		Y.push_back(p[i].y);
	}

	sort(X.begin(), X.end());
	sort(Y.begin(), Y.end());
	Point med;
	if (sz % 2 == 1) {
		med.x = X[int((sz - 1) / 2)];
		med.y = Y[int((sz - 1) / 2)];
	}
	else {
		med.x = (X[int(sz / 2)] + X[int(sz / 2) - 1]) / 2.0;
		med.y = (Y[int(sz / 2)] + Y[int(sz / 2) - 1]) / 2.0;
	}
	return(med);
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

Point mean(std::vector<Point> P, int sz) {
	Point M = Point(0,0);
	for (int i = 0; i < sz; i++) {
		//std::cout << "\n" << P[i];
		M.x += P[i].x;
		M.y += P[i].y;
	}
	M.x = int((double)M.x / double(sz));
	M.y = int((double)M.y / double(sz));
	//std::cout << "\n mean = " << M ;
	return(M);
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

int countCrossings(int **roi, int h, int w, Point p1, Point p2, std::vector<double> &whiteRun) {
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
		double white_dist = 0;
		int crossing;
		Point prev_cross = p1.y < p2.y ? p1 : p2;
		Point end_point = p1.y > p2.y ? p1 : p2;
		for (int y = jstart; y <= jstop; y++) {
			int x = int(p1.x + double((y - p1.y)*(p2.x - p1.x)) / (p2.y - p1.y));
			crossing = checkNbd(roi, h, w, Point(x, y));
			if (crossing && prev == 255) { //white to black...
				count++;
				prev = 0;
				white_dist = norm(prev_cross - Point(x, y));
				if (white_dist > 1.5)
					whiteRun.push_back(white_dist);
				prev_cross = Point(x, y);
			}
			else if (!crossing && prev == 0) {  //black to white...
				prev = 255;
			}
		}
		white_dist = norm(prev_cross - end_point);
		if (white_dist > 1.5)
			whiteRun.push_back(white_dist);
	}
	else if ((slope > 0 && slope < 45) || (slope > 135 && slope < 180)) {
		int istart = p1.x < p2.x ? p1.x : p2.x;
		int istop = p1.x > p2.x ? p1.x : p2.x;
		int white = 0, prev = 255;
		double white_dist = 0;
		int crossing;
		Point prev_cross = p1.x < p2.x ? p1 : p2;
		Point end_point = p1.x > p2.x ? p1 : p2;
		for (int x = istart; x <= istop; x++) {
			int y = int(p1.y + double((x - p1.x)*(p2.y - p1.y)) / (p2.x - p1.x));
			crossing = checkNbd(roi, h, w, Point(x, y));
			if (crossing && prev == 255) { //white to black...
				count++;
				prev = 0;
				white_dist = norm(prev_cross - Point(x, y));
				if (white_dist > 1.5)
					whiteRun.push_back(white_dist);
				prev_cross = Point(x, y);
			}
			else if (!crossing && prev == 0) {  //black to white...
				prev = 255;
			}
		}
		white_dist = norm(prev_cross - end_point);
		if (white_dist > 1.5)
			whiteRun.push_back(white_dist);
	}

	return(count);
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

//************************************************************************************************************************
//find peaks:base ratio...

void smoothingHistogram(int *hist, int sz) {
	int win_sz = 21;
	//int win_sz = 15;
	std::cout << "\nsize of hist=" << sz;
	
	int max = hist[0];
	for (int i = 1; i < sz; i++) {
		max = hist[i] > max ? hist[i] : max;
	}

	std::cout << "\nHistogram: ";
	for (int i = 1; i < sz; i++) {
		hist[i] = hist[i] < 0.2 * max ? 0 : hist[i];
		//std::cout << hist[i] << " , ";
	}

	for (int i = 0; i < sz; i += win_sz / 2) {
		int avg_hist = 0;
		for (int j = 0; j < win_sz; j++) {
			if (i + j >= sz)
				break;
			avg_hist += hist[i + j];
		}
		avg_hist /= win_sz;
		for (int j = 0; j < win_sz; j++) {
			if (i + j >= sz)
				break;
			hist[i + j] = avg_hist;
		}
	}
	
}

std::vector<std::vector<int>> peaksAndValleys(int *hist, int sz) {
	int start, stop;
	for (int i = 0; i < sz; i++) {
		if (hist[i] > 0) {
			start = i - 1;
			break;
		}
	}
	for (int i = sz - 1; i >= 0; i--) {
		if (hist[i] > 0) {
			stop = i + 1;
			break;
		}
	}
	if (start < 0) start = 0;
	if (stop > sz - 1) stop = sz - 1;

	std::cout << "\nstart = " << start << " , stop = " << stop;

	int max = hist[0];
	int *seq = new int[sz]; //1 => possible peak , -1 => possible valley
	for (int i = 1; i < sz; i++) {
		max = hist[i] > max ? hist[i] : max;
	}
	/*double avg = 0;
	for (int i = start; i <= stop; i++) {
		avg += hist[i];
	}
	avg /= (stop - start + 1);*/

	double thresh = 0.3*max;

	for (int i = 0; i < sz; i++) {
		seq[i] = hist[i]>=thresh ? 1 : -1;
	}
	
	
	std::vector<std::vector<int>> pv;
	for (int i = start; i <= stop;) {
		if (seq[i] == 1) {
			int peak_val = 0;
			int avg = 0, count = 0;
			while (seq[i] != -1 && i <= stop) {
				peak_val = hist[i] > peak_val ? hist[i] : peak_val;
				avg += i;
				count++;
				i++;
			}
			avg /= count;
			/*if ((avg - thresh) / (max - thresh) < 0.1)
				continue;*/
			std::vector<int> peak;
			peak.push_back(1);
			peak.push_back(avg);
			peak.push_back(peak_val);
			pv.push_back(peak);
			peak.clear();
		}
		if (seq[i] == -1) {
			int valley_val = 99999;
			int avg = 0, count = 0;
			while (seq[i] != 1 && i <= stop) {
				valley_val = hist[i] < valley_val ? hist[i] : valley_val;
				avg += i;
				count++;
				i++;
			}
			avg /= count;
			/*if ((thresh - avg) / double(avg) < 0.2)
				continue;*/
			std::vector<int> valley;
			valley.push_back(-1);
			valley.push_back(avg);
			valley.push_back(valley_val);
			pv.push_back(valley);
			valley.clear();
		}
	}

	if (pv[0][0] == 1) {//starts with peak...
		std::vector<std::vector<int>> pv_new;
		std::vector<int> valley;
		valley.push_back(-1);
		valley.push_back(start);
		valley.push_back(0);
		pv_new.push_back(valley);
		pv_new.insert(pv_new.end(), pv.begin(), pv.end());
		pv.clear();
		pv = pv_new;
	}
	if (pv[int(size(pv)) - 1][0] == 1) {   //ends with peak...
		std::vector<int> valley;
		valley.push_back(-1);
		valley.push_back(stop);
		valley.push_back(0);
		pv.push_back(valley);
	}

	std::cout << "\nSequence of peaks and valleys:";
	for (int i = 0; i<int(size(pv)); i++) {
		if (pv[i][0] == 1) {
			std::cout << "\n Peak at i = " << pv[i][1] << " , height = " << pv[i][2];
		}
		if (pv[i][0] == -1) {
			std::cout << "\n Valley at i = " << pv[i][1] << " , height = " << pv[i][2];
		}
	}

	return(pv);
}

std::pair<int, double> peaks2base1(int **roi, int h, int w) {
	int *hist = new int[h];
	for (int i = 0; i < h; i++) {
		hist[i] = 0;
		int count = 0, prev = 255;
		for (int j = 0; j < w; j++) {
			if (roi[i][j] < 255) hist[i]++;
			if (prev == 255 && roi[i][j] < 255) count++;
			prev = roi[i][j];
		}
		hist[i] *= count;
	}

	smoothingHistogram(hist, h);
	std::vector<std::vector<int>> pv = peaksAndValleys(hist, h);

	int **col_roi = new int*[h];
	for (int i = 0; i < h; i++) {
		col_roi[i] = new int[w];
		for (int j = 0; j < w; j++) {
			if (j <= hist[i] && roi[i][j] == 255)
				col_roi[i][j] = 1;
			else if (j <= hist[i] && roi[i][j] < 255)
				col_roi[i][j] = 2;
			else if (j > hist[i] && roi[i][j] < 255)
				col_roi[i][j] = 0;
			else if (j > hist[i] && roi[i][j] == 255)
				col_roi[i][j] = 255;
		}
	}
	Mat colroi;
	colroi.create(h, w, CV_8UC1);
	cvtColor(colroi, colroi, CV_GRAY2BGR);
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			if(col_roi[i][j] == 0) 
				colroi.at<Vec3b>(i, j) = Vec3b(0, 0, 0);
			if (col_roi[i][j] == 255)
				colroi.at<Vec3b>(i, j) = Vec3b(255, 255, 255);
			if (col_roi[i][j] == 1)
				colroi.at<Vec3b>(i, j) = Vec3b(0, 0, 255);
			if (col_roi[i][j] == 2)
				colroi.at<Vec3b>(i, j) = Vec3b(0, 0, 150);
		}
	}

	for (int t = 0; t < int(pv.size()); t++) {
		if (pv[t][0] == -1) 
			line(colroi, Point(0, pv[t][1]), Point(w - 1, pv[t][1]), Scalar(255, 0, 0), 4, 8);
		if (pv[t][0] == 1)
			line(colroi, Point(0, pv[t][1]), Point(w - 1, pv[t][1]), Scalar(0, 255, 0), 4, 8);
	}

	std::pair<int, double> data; //base, ratio...
	double max_ratio = -1;
	int opt_base;
	for (int t = 0; t<int(pv.size()); t++) {
		if (pv[t][0] == -1 && t + 2 < int(pv.size())) { //valley...
			std::cout << "\npeak height = " << pv[t + 1][2] << " , base = " << pv[t + 2][1] - pv[t][1] << " , ratio = " << double(pv[t + 1][2]) / (pv[t + 2][1] - pv[t][1]);
			if (double(pv[t + 1][2]) / (pv[t + 2][1] - pv[t][1]) > max_ratio) {
				max_ratio = double(pv[t + 1][2]) / (pv[t + 2][1] - pv[t][1]);
				opt_base = pv[t + 2][1] - pv[t][1];
			}
		}
	}

	data.first = opt_base;
	data.second = max_ratio;

	//imwrite("Data/ICDAR/Data155/roi_HPP.tif", colroi);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/roi_HPP.tif", colroi);
	return(data);

}

//new one...
std::vector<std::vector<int>> peaks2base(int **roi, int h, int w, std::vector<int>& base_range) {
	Mat ROI;
	ROI.create(h, w, CV_8UC1);
	arr2mat(roi, ROI, h, w);
	GaussianBlur(ROI, ROI, Size(45, 9), 0, 0);
	//GaussianBlur(ROI, ROI, Size(21, 5), 0, 0); //for ICDAR
	mat2arr(ROI, roi, h, w);
	int th = 170;

	int *hist = new int[h];
	for (int i = 0; i < h; i++) {
		hist[i] = 0;
		int count = 0, prev = 255;
		for (int j = 0; j < w; j++) {
			if (roi[i][j] < th) hist[i]++;
			if (prev > th && roi[i][j] < th) count++;
			prev = roi[i][j];
		}
		//hist[i] *= count;
	}

	smoothingHistogram(hist, h);
	std::vector<std::vector<int>> pv = peaksAndValleys(hist, h);

	int range_start, range_stop;
	for (int t = 0; t < int(pv.size()) - 2; t+=2) {
		if (pv[t][0] == -1) { //valley...
			range_start = pv[t][1], range_stop = pv[t+2][1];
			while (hist[range_start] < 5) {
				range_start++;
			}
			while (hist[range_stop] < 5) {
				range_stop--;
			}
		}
		base_range.push_back(range_stop - range_start);
	}
	std::cout << "\nBase ranges: ";
	for (int i = 0; i<int(base_range.size()); i++) {
		std::cout << base_range[i] << " , ";
	}
	
	int count = 0;
	for (int t = 0; t < int(pv.size()); t++) {
		if (pv[t][0] == 1)
			count++;
	}
	std::cout << "\nNo. of peaks = " << count;

	int **col_roi = new int*[h];
	for (int i = 0; i < h; i++) {
		col_roi[i] = new int[w];
		for (int j = 0; j < w; j++) {
			if (j <= hist[i] && roi[i][j] >= th)
				col_roi[i][j] = 1;
			else if (j <= hist[i] && roi[i][j] < th)
				col_roi[i][j] = 2;
			else if (j > hist[i] && roi[i][j] < th)
				col_roi[i][j] = 0;
			else if (j > hist[i] && roi[i][j] > th)
				col_roi[i][j] = 255;
		}
	}
	Mat colroi;
	colroi.create(h, w, CV_8UC1);
	cvtColor(colroi, colroi, CV_GRAY2BGR);
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			if (col_roi[i][j] == 0)
				colroi.at<Vec3b>(i, j) = Vec3b(0, 0, 0);
			if (col_roi[i][j] == 255)
				colroi.at<Vec3b>(i, j) = Vec3b(255, 255, 255);
			if (col_roi[i][j] == 1)
				colroi.at<Vec3b>(i, j) = Vec3b(0, 0, 255);
			if (col_roi[i][j] == 2)
				colroi.at<Vec3b>(i, j) = Vec3b(0, 0, 150);
		}
	}

	for (int t = 0; t < int(pv.size()); t++) {
		if (pv[t][0] == -1)
			line(colroi, Point(0, pv[t][1]), Point(w - 1, pv[t][1]), Scalar(255, 0, 0), 4, 8);
		if (pv[t][0] == 1)
			line(colroi, Point(0, pv[t][1]), Point(w - 1, pv[t][1]), Scalar(0, 255, 0), 4, 8);
	}

	//imwrite("Data/ICDAR/Data155/roi_HPP.tif", colroi);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/roi_HPP.tif", colroi);

	return(pv);

}

double pvd(std::vector<std::vector<int>> pv) {
	double avg_pvd = 0;
	for (int i = 0; i < int(pv.size()); i+=2) {
		int p_i, p_j, v_i;
		double pvd_i = 0;
		if (i - 1 < 0)
			p_i = pv[i + 1][2];
		else
			p_i = pv[i - 1][2];
		if (i + 1 >= int(pv.size()))
			p_j = pv[i - 1][2];
		else
			p_j = pv[i + 1][2];
		v_i = pv[i][2];
		pvd_i = (p_i + p_j) / 2.0 - v_i;
		avg_pvd += pvd_i;
	}
	avg_pvd /= (int(pv.size()) + 1)*0.5;
	return(avg_pvd);
}

//***********************************************************************************************************************************************************************
//rotate image...
cv::Mat imageRotate(Mat orgimg, Mat rotimg, double angle) {
	//invert bw image...
	/*int **org = new int*[orgimg.rows];
	for (int i = 0; i < orgimg.rows; i++) {
		org[i] = new int[orgimg.cols];
		for (int j = 0; j < orgimg.cols; j++) {
			org[i][j] = 255 - orgimg.at<uchar>(i, j);
		}
	}

	arr2mat(org, orgimg, orgimg.rows, orgimg.cols);*/

	// get rotation matrix for rotating the image around its center in pixel coordinates
	cv::Point2f center((orgimg.cols - 1) / 2.0, (orgimg.rows - 1) / 2.0);
	cv::Mat rot = cv::getRotationMatrix2D(center, angle, 1.0);
	// determine bounding rectangle, center not relevant
	cv::Rect2f bbox = cv::RotatedRect(cv::Point2f(), orgimg.size(), angle).boundingRect();
	// adjust transformation matrix
	rot.at<double>(0, 2) += bbox.width / 2.0 - orgimg.cols / 2.0;
	rot.at<double>(1, 2) += bbox.height / 2.0 - orgimg.rows / 2.0;

	cv::warpAffine(orgimg, rotimg, rot, bbox.size());

	//re-invert final image...
	int **rot_im = new int*[rotimg.rows];
	for (int i = 0; i < rotimg.rows; i++) {
		rot_im[i] = new int[rotimg.cols];
		for (int j = 0; j < rotimg.cols; j++) {
			rot_im[i][j] = 255 - rotimg.at<uchar>(i, j);
		}
	}

	arr2mat(rot_im, rotimg, rotimg.rows, rotimg.cols);

	//imwrite("Data/ICDAR/Data155/rotated_im.tif", rotimg);
	cv::imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/rotated_im.tif", rotimg);
	return(rotimg);
}

cv::Mat imageRotate1(Mat orgimg, Mat rotimg, double angle) {
	
	// get rotation matrix for rotating the image around its center in pixel coordinates
	cv::Point2f center((orgimg.cols - 1) / 2.0, (orgimg.rows - 1) / 2.0);
	cv::Mat rot = cv::getRotationMatrix2D(center, angle, 1.0);
	// determine bounding rectangle, center not relevant
	cv::Rect2f bbox = cv::RotatedRect(cv::Point2f(), orgimg.size(), angle).boundingRect();
	// adjust transformation matrix
	rot.at<double>(0, 2) += bbox.width / 2.0 - orgimg.cols / 2.0;
	rot.at<double>(1, 2) += bbox.height / 2.0 - orgimg.rows / 2.0;

	cv::warpAffine(orgimg, rotimg, rot, bbox.size());

	//cv::imwrite("Data/ICDAR/Data155/rotated_im.tif", rotimg);
	cv::imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/rotated_im.tif", rotimg);
	return(rotimg);
}

//orientations of combined labels...
std::pair<int, double> findOrientation(Mat roi_img, int **roi, int h, int w, std::vector<pair<int, int>> connections, std::vector<std::pair<int, Point>> all_centroids, int& base, std::vector<Point>& lines, std::vector<std::vector<int>> connectLabelsAngles) {
	double theta, max_theta, min_max_theta;
	int min_max_range = 9999;

	//invert roi...
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			roi[i][j] = 255 - roi[i][j];
		}
	}

	arr2mat(roi, roi_img, h, w);

	double max_ratio = -1;
	std::pair<int, double> p2bratio;
	std::vector<std::vector<int>> pv_max;

	std::vector<double> thetaList;
	//list the angles between actual connections...
	for (int i = 0; i<int(connections.size()); i++) {
		for (int k = 0; k < int(connectLabelsAngles.size()); k++) {
			if ((connectLabelsAngles[k][0] == connections[i].first && connectLabelsAngles[k][1] == connections[i].second) || 
				(connectLabelsAngles[k][1] == connections[i].first && connectLabelsAngles[k][0] == connections[i].second)) {
				theta = double(connectLabelsAngles[k][2]);
				break;
			}
		}

		if (theta > 90 && theta < 180) theta = theta - 180;
		thetaList.push_back(theta);
	}

	//list the angles between centroids...
	for (int i = 0; i<int(connections.size()); i++) {
		Point p1, p2;
		int f1 = 0, f2 = 0;
		for (int k = 0; k < int(all_centroids.size()); k++) {
			if (all_centroids[k].first == connections[i].first) {
				p1 = all_centroids[k].second;
				f1 = 1;
			}
			if (all_centroids[k].first == connections[i].second) {
				p2 = all_centroids[k].second;
				f2 = 1;
			}
			if (f1 && f2) break;
		}

		theta = atan2(p2.y - p1.y, p2.x - p1.x) * 180 / CV_PI;
		if (theta > 90 && theta < 180) theta = theta - 180;
		thetaList.push_back(theta);
	}

	for (int i = 0; i<int(thetaList.size()); i++) {
		theta = thetaList[i];
		std::cout << "\nRotated by " << theta << " degrees";

		Mat rot_img;
		rot_img = imageRotate(roi_img, rot_img, int(theta));
		//imwrite("Data/ICDAR/Data155/roi_rot.tif", rot_img);
		imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/roi_rot.tif", rot_img);

		int **rotimg = new int*[rot_img.rows];
		for (int k = 0; k < rot_img.rows; k++) {
			rotimg[k] = new int[rot_img.cols];
			for (int l = 0; l < rot_img.cols; l++)
				rotimg[k][l] = 255;
		}

		mat2arr(rot_img, rotimg, rot_img.rows, rot_img.cols);
		std::vector<int> base_range;
		std::vector<std::vector<int>> pv = peaks2base(rotimg, rot_img.rows, rot_img.cols, base_range);

		//find best p:b ratio...
		std::pair<int, double> data; //base, ratio...
		data.first = 0;
		data.second = -1;

		for (int t = 0; t<int(pv.size()); t++) {
			if (pv[t][0] == -1 && t + 2 < int(pv.size())) { //valley...
				std::cout << "\npeak height = " << pv[t + 1][2] << " , base = " << base_range[t/2] << " , ratio = " << double(pv[t + 1][2]) / (base_range[t / 2]);
				/*if (double(pv[t + 1][2]) / (base_range[t / 2]) > data.second) {
					data.second = double(pv[t + 1][2]) / (base_range[t / 2]);
					data.first = base_range[t / 2];
				}*/
			}
		}

		//Peak valley distance...
		std::cout<<"\nPVD = "<<pvd(pv);

		//maximum base range and corresponding ratio...
		int max_range = base_range[0];
		int pk_ht = pv[1][2];
		if (int(base_range.size())>1) {
			for (int t = 1; t<int(base_range.size()); t++) {
				if (base_range[t] > max_range) {
					max_range = base_range[t];
					pk_ht = pv[2 * t + 1][2];
				}
			}
		}
		if (double(pk_ht) / (max_range) > data.second) {
			data.second = double(pk_ht) / (max_range);
			data.first = max_range;
		}
		std::cout << "\nMax Base Range = " << max_range;
		if (max_range < min_max_range) {
			min_max_range = max_range;
			min_max_theta = theta;
		}


		//for (int t = 0; t<int(pv.size()); t++) {
		//	if (pv[t][0] == -1 && t + 2 < int(pv.size())) { //valley...
		//		std::cout << "\npeak height = " << pv[t + 1][2] << " , base = " << pv[t + 2][1] - pv[t][1] << " , ratio = " << double(pv[t + 1][2]) / (pv[t + 2][1] - pv[t][1]);
		//		if (double(pv[t + 1][2]) / (pv[t + 2][1] - pv[t][1]) > data.second) {
		//			data.second = double(pv[t + 1][2]) / (pv[t + 2][1] - pv[t][1]);
		//			data.first = pv[t + 2][1] - pv[t][1];
		//		}
		//	}
		//}

		//compare with other rotations...
		if (data.second > max_ratio) {
			max_ratio = data.second;
			max_theta = theta;
			pv_max.clear();
			pv_max = pv;
		}

		rot_img.release();
	}

	/*int count = 0;
	for (int t = 0; t<int(pv_max.size()); t++) {
		if (pv_max[t][0] == 1)
			count++;
	}*/

	std::cout << "\n**************Best orientation is " << min_max_theta << " with range = " << min_max_range;
	
	Mat rot_img;
	rot_img = imageRotate(roi_img, rot_img, int(max_theta));
	for (int k = 0; k < rot_img.rows; k++) {
		for (int l = 0; l < rot_img.cols; l++) {
			if (rot_img.at<uchar>(k, l) == 0)
				rot_img.at<uchar>(k, l) = 255;
			else if (rot_img.at<uchar>(k, l) == 255)
				rot_img.at<uchar>(k, l) = 0;
		}
	}
	
	Mat colroi, colrot;
	colroi.create(rot_img.rows, rot_img.cols, CV_8UC1);
	cvtColor(rot_img, colroi, CV_GRAY2BGR);
	for (int t = 0; t < int(pv_max.size()); t++) {
		if (pv_max[t][0] == -1)
			line(colroi, Point(0, pv_max[t][1]), Point(rot_img.cols - 1, pv_max[t][1]), Scalar(255, 0, 0), 4, 8);
		if (pv_max[t][0] == 1)
			line(colroi, Point(0, pv_max[t][1]), Point(rot_img.cols - 1, pv_max[t][1]), Scalar(0, 255, 0), 4, 8);
	}
		
	colrot = imageRotate1(colroi, colrot, int(-max_theta));

	//crop roi...
	int x1, x2, y1, y2;
	for (int y = 0; y < colrot.rows; y++) {
		int flag = 0;
		for (int x = 0; x < colrot.cols; x++) {
			if (colrot.at<Vec3b>(y, x) == Vec3b(255, 255, 255)) {
				flag = 1;
				break;
			}
		}
		if (flag) {
			y1 = y;
			break;
		}
		
	}
	for (int y = colrot.rows - 1; y >= 0 ; y--) {
		int flag = 0;
		for (int x = 0; x < colrot.cols; x++) {
			if (colrot.at<Vec3b>(y, x) == Vec3b(255, 255, 255)) {
				flag = 1;
				break;
			}
		}
		if (flag) {
			y2 = y;
			break;
		}

	}
	for (int x = 0; x < colrot.cols; x++) {
		int flag = 0;
		for (int y = 0; y < colrot.rows; y++) {
			if (colrot.at<Vec3b>(y, x) == Vec3b(255, 255, 255)) {
				flag = 1;
				break;
			}
		}
		if (flag) {
			x1 = x;
			break;
		}

	}
	for (int x = colrot.cols - 1; x >= 0; x--) {
		int flag = 0;
		for (int y = 0; y < colrot.rows; y++) {
			if (colrot.at<Vec3b>(y, x) == Vec3b(255, 255, 255)) {
				flag = 1;
				break;
			}
		}
		if (flag) {
			x2 = x;
			break;
		}

	}
	
	//cv::Rect myROI(y1, x1, y2 - y1 + 1,  x2 - x1 + 1);
	cv::Rect myROI(x1, y1, x2 - x1 + 1, y2 - y1 + 1);
	Mat croppedImage = colrot(myROI);
	std::cout << "\nroi: " << h << "x" << w << " , final: " << croppedImage.rows << "x" << croppedImage.cols;
	
	//imwrite("Data/ICDAR/Data155/roi_final.tif", croppedImage);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/roi_final.tif", croppedImage);

	p2bratio.first = int(max_theta);
	p2bratio.second = max_ratio;

	//std::vector<Point> lines;
	int x = int((croppedImage.cols - 1) / 4.0);
	for (int y = 0; y < croppedImage.rows; y++) {
		if (croppedImage.at<Vec3b>(y, x) == Vec3b(0, 255, 0))
			lines.push_back(Point(x, y));
	}
	x = int(3*(croppedImage.cols - 1) / 4.0);
	for (int y = 0; y < croppedImage.rows; y++) {
		if (croppedImage.at<Vec3b>(y, x) == Vec3b(0, 255, 0))
			lines.push_back(Point(x, y));
	}
	int y = int((croppedImage.rows - 1) / 2.0);
	for (int x = 0; x < croppedImage.cols; x++) {
		if (croppedImage.at<Vec3b>(y, x) == Vec3b(0, 255, 0))
			lines.push_back(Point(x, y));
	}
	return(p2bratio);
}

//orientations of single labels...
std::pair<int, double> findOrientation(Mat roi_img, int **roi, int h, int w, std::vector<std::pair<int, Point>> all_centroids, std::vector<Point> roi_all_terminals, int& base, std::vector<Point>& lines) {
	double theta, max_theta, min_max_theta;
	int min_max_range = 9999;

	//invert roi...
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			roi[i][j] = 255 - roi[i][j];
		}
	}

	arr2mat(roi, roi_img, h, w);

	double max_ratio = -1;
	std::pair<int, double> p2bratio;
	std::vector<std::vector<int>> pv_max;

	//list the angles by joining the terminals pairwise...
	std::vector<double> thetaList;
	for (int i = 1; i<int(roi_all_terminals.size()); i++) {
		Point p1 = roi_all_terminals[i];
		int f1 = 0, f2 = 0;
		for (int k = i + 1; k < int(roi_all_terminals.size()); k++) {
			Point p2 = roi_all_terminals[k];
			theta = atan2(p2.y - p1.y, p2.x - p1.x) * 180 / CV_PI;
			if (theta > 90 && theta < 180) theta = theta - 180;
			thetaList.push_back(theta);
		}		
	}

	for (int i = 0; i<int(thetaList.size()); i++) {
		theta = thetaList[i];
		std::cout << "\nRotated by " << theta << " degrees";

		Mat rot_img;
		rot_img = imageRotate(roi_img, rot_img, int(theta));
		//imwrite("Data/ICDAR/Data155/roi_rot.tif", rot_img);
		imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/roi_rot.tif", rot_img);

		int **rotimg = new int*[rot_img.rows];
		for (int k = 0; k < rot_img.rows; k++) {
			rotimg[k] = new int[rot_img.cols];
			for (int l = 0; l < rot_img.cols; l++)
				rotimg[k][l] = 255;
		}

		mat2arr(rot_img, rotimg, rot_img.rows, rot_img.cols);
		std::vector<int> base_range;
		std::vector<std::vector<int>> pv = peaks2base(rotimg, rot_img.rows, rot_img.cols, base_range);

		//find best p:b ratio...
		std::pair<int, double> data; //base, ratio...
		data.first = 0;
		data.second = -1;

		for (int t = 0; t<int(pv.size()); t++) {
			if (pv[t][0] == -1 && t + 2 < int(pv.size())) { //valley...
				std::cout << "\npeak height = " << pv[t + 1][2] << " , base = " << base_range[t / 2] << " , ratio = " << double(pv[t + 1][2]) / (base_range[t / 2]);
				/*if (double(pv[t + 1][2]) / (base_range[t / 2]) > data.second) {
				data.second = double(pv[t + 1][2]) / (base_range[t / 2]);
				data.first = base_range[t / 2];
				}*/
			}
		}

		//Peak valley distance...
		std::cout << "\nPVD = " << pvd(pv);

		//maximum base range and corresponding ratio...
		int max_range = base_range[0];
		int pk_ht = pv[1][2];
		if (int(base_range.size())>1) {
			for (int t = 1; t<int(base_range.size()); t++) {
				if (base_range[t] > max_range) {
					max_range = base_range[t];
					pk_ht = pv[2 * t + 1][2];
				}
			}
		}
		if (double(pk_ht) / (max_range) > data.second) {
			data.second = double(pk_ht) / (max_range);
			data.first = max_range;
		}
		std::cout << "\nMax Base Range = " << max_range;
		if (max_range < min_max_range) {
			min_max_range = max_range;
			min_max_theta = theta;
		}


		//for (int t = 0; t<int(pv.size()); t++) {
		//	if (pv[t][0] == -1 && t + 2 < int(pv.size())) { //valley...
		//		std::cout << "\npeak height = " << pv[t + 1][2] << " , base = " << pv[t + 2][1] - pv[t][1] << " , ratio = " << double(pv[t + 1][2]) / (pv[t + 2][1] - pv[t][1]);
		//		if (double(pv[t + 1][2]) / (pv[t + 2][1] - pv[t][1]) > data.second) {
		//			data.second = double(pv[t + 1][2]) / (pv[t + 2][1] - pv[t][1]);
		//			data.first = pv[t + 2][1] - pv[t][1];
		//		}
		//	}
		//}

		//compare with other rotations...
		if (data.second > max_ratio) {
			max_ratio = data.second;
			max_theta = theta;
			pv_max.clear();
			pv_max = pv;
		}

		rot_img.release();
	}

	std::cout << "\n**************Best orientation is " << min_max_theta << " with range = " << min_max_range;

	Mat rot_img;
	rot_img = imageRotate(roi_img, rot_img, int(max_theta));
	for (int k = 0; k < rot_img.rows; k++) {
		for (int l = 0; l < rot_img.cols; l++) {
			if (rot_img.at<uchar>(k, l) == 0)
				rot_img.at<uchar>(k, l) = 255;
			else if (rot_img.at<uchar>(k, l) == 255)
				rot_img.at<uchar>(k, l) = 0;
		}
	}

	Mat colroi, colrot;
	colroi.create(rot_img.rows, rot_img.cols, CV_8UC1);
	cvtColor(rot_img, colroi, CV_GRAY2BGR);
	for (int t = 0; t < int(pv_max.size()); t++) {
		if (pv_max[t][0] == -1)
			line(colroi, Point(0, pv_max[t][1]), Point(rot_img.cols - 1, pv_max[t][1]), Scalar(255, 0, 0), 4, 8);
		if (pv_max[t][0] == 1)
			line(colroi, Point(0, pv_max[t][1]), Point(rot_img.cols - 1, pv_max[t][1]), Scalar(0, 255, 0), 4, 8);
	}

	colrot = imageRotate1(colroi, colrot, int(-max_theta));

	//crop roi...
	int x1, x2, y1, y2;
	for (int y = 0; y < colrot.rows; y++) {
		int flag = 0;
		for (int x = 0; x < colrot.cols; x++) {
			if (colrot.at<Vec3b>(y, x) == Vec3b(255, 255, 255)) {
				flag = 1;
				break;
			}
		}
		if (flag) {
			y1 = y;
			break;
		}

	}
	for (int y = colrot.rows - 1; y >= 0; y--) {
		int flag = 0;
		for (int x = 0; x < colrot.cols; x++) {
			if (colrot.at<Vec3b>(y, x) == Vec3b(255, 255, 255)) {
				flag = 1;
				break;
			}
		}
		if (flag) {
			y2 = y;
			break;
		}

	}
	for (int x = 0; x < colrot.cols; x++) {
		int flag = 0;
		for (int y = 0; y < colrot.rows; y++) {
			if (colrot.at<Vec3b>(y, x) == Vec3b(255, 255, 255)) {
				flag = 1;
				break;
			}
		}
		if (flag) {
			x1 = x;
			break;
		}

	}
	for (int x = colrot.cols - 1; x >= 0; x--) {
		int flag = 0;
		for (int y = 0; y < colrot.rows; y++) {
			if (colrot.at<Vec3b>(y, x) == Vec3b(255, 255, 255)) {
				flag = 1;
				break;
			}
		}
		if (flag) {
			x2 = x;
			break;
		}

	}

	//cv::Rect myROI(y1, x1, y2 - y1 + 1,  x2 - x1 + 1);
	cv::Rect myROI(x1, y1, x2 - x1 + 1, y2 - y1 + 1);
	Mat croppedImage = colrot(myROI);
	std::cout << "\nroi: " << h << "x" << w << " , final: " << croppedImage.rows << "x" << croppedImage.cols;

	//imwrite("Data/ICDAR/Data155/roi_final.tif", croppedImage);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/roi_final.tif", croppedImage);

	p2bratio.first = int(max_theta);
	p2bratio.second = max_ratio;

	//std::vector<Point> lines;
	int x = int((croppedImage.cols - 1) / 4.0);
	for (int y = 0; y < croppedImage.rows; y++) {
		if (croppedImage.at<Vec3b>(y, x) == Vec3b(0, 255, 0))
			lines.push_back(Point(x, y));
	}
	x = int(3 * (croppedImage.cols - 1) / 4.0);
	for (int y = 0; y < croppedImage.rows; y++) {
		if (croppedImage.at<Vec3b>(y, x) == Vec3b(0, 255, 0))
			lines.push_back(Point(x, y));
	}
	int y = int((croppedImage.rows - 1) / 2.0);
	for (int x = 0; x < croppedImage.cols; x++) {
		if (croppedImage.at<Vec3b>(y, x) == Vec3b(0, 255, 0))
			lines.push_back(Point(x, y));
	}
	return(p2bratio);


}

//orientation of connecting line pairs...
std::pair<int, double> findOrientation(Mat roi_img, int **roi, int h, int w, std::vector<int> thetaList, std::vector<Point>& lines) {
	double theta, max_theta, min_max_theta;
	int min_max_range = 9999;

	//invert roi...
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			roi[i][j] = 255 - roi[i][j];
		}
	}

	arr2mat(roi, roi_img, h, w);

	double max_ratio = -1;
	std::pair<int, double> p2bratio;
	std::vector<std::vector<int>> pv_max;

	

	for (int i = 0; i<int(thetaList.size()); i++) {
		theta = thetaList[i];
		std::cout << "\nRotated by " << theta << " degrees";

		Mat rot_img;
		rot_img = imageRotate(roi_img, rot_img, int(theta));
		//imwrite("Data/ICDAR/Data155/roi_rot.tif", rot_img);
		imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/roi_rot.tif", rot_img);

		int **rotimg = new int*[rot_img.rows];
		for (int k = 0; k < rot_img.rows; k++) {
			rotimg[k] = new int[rot_img.cols];
			for (int l = 0; l < rot_img.cols; l++)
				rotimg[k][l] = 255;
		}

		mat2arr(rot_img, rotimg, rot_img.rows, rot_img.cols);
		std::vector<int> base_range;
		std::vector<std::vector<int>> pv = peaks2base(rotimg, rot_img.rows, rot_img.cols, base_range);

		//find best p:b ratio...
		std::pair<int, double> data; //base, ratio...
		data.first = 0;
		data.second = -1;

		for (int t = 0; t<int(pv.size()); t++) {
			if (pv[t][0] == -1 && t + 2 < int(pv.size())) { //valley...
				std::cout << "\npeak height = " << pv[t + 1][2] << " , base = " << base_range[t / 2] << " , ratio = " << double(pv[t + 1][2]) / (base_range[t / 2]);
				/*if (double(pv[t + 1][2]) / (base_range[t / 2]) > data.second) {
				data.second = double(pv[t + 1][2]) / (base_range[t / 2]);
				data.first = base_range[t / 2];
				}*/
			}
		}

		//Peak valley distance...
		std::cout << "\nPVD = " << pvd(pv);

		//maximum base range and corresponding ratio...
		int max_range = base_range[0];
		int pk_ht = pv[1][2];
		if (int(base_range.size())>1) {
			for (int t = 1; t<int(base_range.size()); t++) {
				if (base_range[t] > max_range) {
					max_range = base_range[t];
					pk_ht = pv[2 * t + 1][2];
				}
			}
		}
		if (double(pk_ht) / (max_range) > data.second) {
			data.second = double(pk_ht) / (max_range);
			data.first = max_range;
		}
		std::cout << "\nMax Base Range = " << max_range;
		if (max_range < min_max_range) {
			min_max_range = max_range;
			min_max_theta = theta;
		}


		//for (int t = 0; t<int(pv.size()); t++) {
		//	if (pv[t][0] == -1 && t + 2 < int(pv.size())) { //valley...
		//		std::cout << "\npeak height = " << pv[t + 1][2] << " , base = " << pv[t + 2][1] - pv[t][1] << " , ratio = " << double(pv[t + 1][2]) / (pv[t + 2][1] - pv[t][1]);
		//		if (double(pv[t + 1][2]) / (pv[t + 2][1] - pv[t][1]) > data.second) {
		//			data.second = double(pv[t + 1][2]) / (pv[t + 2][1] - pv[t][1]);
		//			data.first = pv[t + 2][1] - pv[t][1];
		//		}
		//	}
		//}

		//compare with other rotations...
		if (data.second > max_ratio) {
			max_ratio = data.second;
			max_theta = theta;
			pv_max.clear();
			pv_max = pv;
		}

		rot_img.release();
	}

	std::cout << "\n**************Best orientation is " << min_max_theta << " with range = " << min_max_range;

	Mat rot_img;
	rot_img = imageRotate(roi_img, rot_img, int(max_theta));
	for (int k = 0; k < rot_img.rows; k++) {
		for (int l = 0; l < rot_img.cols; l++) {
			if (rot_img.at<uchar>(k, l) == 0)
				rot_img.at<uchar>(k, l) = 255;
			else if (rot_img.at<uchar>(k, l) == 255)
				rot_img.at<uchar>(k, l) = 0;
		}
	}

	Mat colroi, colrot;
	colroi.create(rot_img.rows, rot_img.cols, CV_8UC1);
	cvtColor(rot_img, colroi, CV_GRAY2BGR);
	for (int t = 0; t < int(pv_max.size()); t++) {
		if (pv_max[t][0] == -1)
			line(colroi, Point(0, pv_max[t][1]), Point(rot_img.cols - 1, pv_max[t][1]), Scalar(255, 0, 0), 4, 8);
		if (pv_max[t][0] == 1)
			line(colroi, Point(0, pv_max[t][1]), Point(rot_img.cols - 1, pv_max[t][1]), Scalar(0, 255, 0), 4, 8);
	}

	colrot = imageRotate1(colroi, colrot, int(-max_theta));

	//crop roi...
	int x1, x2, y1, y2;
	for (int y = 0; y < colrot.rows; y++) {
		int flag = 0;
		for (int x = 0; x < colrot.cols; x++) {
			if (colrot.at<Vec3b>(y, x) == Vec3b(255, 255, 255)) {
				flag = 1;
				break;
			}
		}
		if (flag) {
			y1 = y;
			break;
		}

	}
	for (int y = colrot.rows - 1; y >= 0; y--) {
		int flag = 0;
		for (int x = 0; x < colrot.cols; x++) {
			if (colrot.at<Vec3b>(y, x) == Vec3b(255, 255, 255)) {
				flag = 1;
				break;
			}
		}
		if (flag) {
			y2 = y;
			break;
		}

	}
	for (int x = 0; x < colrot.cols; x++) {
		int flag = 0;
		for (int y = 0; y < colrot.rows; y++) {
			if (colrot.at<Vec3b>(y, x) == Vec3b(255, 255, 255)) {
				flag = 1;
				break;
			}
		}
		if (flag) {
			x1 = x;
			break;
		}

	}
	for (int x = colrot.cols - 1; x >= 0; x--) {
		int flag = 0;
		for (int y = 0; y < colrot.rows; y++) {
			if (colrot.at<Vec3b>(y, x) == Vec3b(255, 255, 255)) {
				flag = 1;
				break;
			}
		}
		if (flag) {
			x2 = x;
			break;
		}

	}

	//cv::Rect myROI(y1, x1, y2 - y1 + 1,  x2 - x1 + 1);
	cv::Rect myROI(x1, y1, x2 - x1 + 1, y2 - y1 + 1);
	Mat croppedImage = colrot(myROI);
	std::cout << "\nroi: " << h << "x" << w << " , final: " << croppedImage.rows << "x" << croppedImage.cols;

	//imwrite("Data/ICDAR/Data155/roi_final.tif", croppedImage);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/roi_final.tif", croppedImage);

	p2bratio.first = int(max_theta);
	p2bratio.second = max_ratio;

	//std::vector<Point> lines;
	int x = int((croppedImage.cols - 1) / 4.0);
	for (int y = 0; y < croppedImage.rows; y++) {
		if (croppedImage.at<Vec3b>(y, x) == Vec3b(0, 255, 0))
			lines.push_back(Point(x, y));
	}
	x = int(3 * (croppedImage.cols - 1) / 4.0);
	for (int y = 0; y < croppedImage.rows; y++) {
		if (croppedImage.at<Vec3b>(y, x) == Vec3b(0, 255, 0))
			lines.push_back(Point(x, y));
	}
	int y = int((croppedImage.rows - 1) / 2.0);
	for (int x = 0; x < croppedImage.cols; x++) {
		if (croppedImage.at<Vec3b>(y, x) == Vec3b(0, 255, 0))
			lines.push_back(Point(x, y));
	}
	return(p2bratio);


}

//orientation in single direction...
std::pair<int, double> findOrientation(Mat roi_img, int **roi, int h, int w, int theta) {
	
	//invert roi...
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			roi[i][j] = 255 - roi[i][j];
		}
	}

	arr2mat(roi, roi_img, h, w);

	double max_ratio = -1;
	std::pair<int, double> p2bratio;

	std::cout << "\nRotated by " << theta << " degrees";

	Mat rot_img;
	rot_img = imageRotate(roi_img, rot_img, int(theta));
	//imwrite("Data/ICDAR/Data155/roi_rot.tif", rot_img);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/roi_rot.tif", rot_img);

	int **rotimg = new int*[rot_img.rows];
	for (int k = 0; k < rot_img.rows; k++) {
		rotimg[k] = new int[rot_img.cols];
		for (int l = 0; l < rot_img.cols; l++)
			rotimg[k][l] = 255;
	}

	mat2arr(rot_img, rotimg, rot_img.rows, rot_img.cols);
	std::vector<int> base_range;
	std::vector<std::vector<int>> pv = peaks2base(rotimg, rot_img.rows, rot_img.cols, base_range);

	//find p:b ratio...
	std::pair<int, double> data; //base, ratio...
	for (int t = 0; t<int(pv.size()); t++) {
		if (pv[t][0] == -1 && t + 2 < int(pv.size())) { //valley...
			std::cout << "\npeak height = " << pv[t + 1][2] << " , base = " << base_range[t / 2] << " , ratio = " << double(pv[t + 1][2]) / (base_range[t / 2]);
		}
	}
	//Peak valley distance...
	std::cout << "\nPVD = " << pvd(pv);

	//maximum base range and corresponding ratio...
	int max_range = base_range[0];
	int pk_ht = pv[1][2];
	if (int(base_range.size())>1) {
		for (int t = 1; t<int(base_range.size()); t++) {
			if (base_range[t] > max_range) {
				max_range = base_range[t];
				pk_ht = pv[2 * t + 1][2];
			}
		}
	}
	if (double(pk_ht) / (max_range) > data.second) {
		p2bratio.second = double(pk_ht) / (max_range);
		p2bratio.first = max_range;
	}
	std::cout << "\nMax Base Range = " << max_range;

	return(p2bratio);
}



//*********************************************************************************************************************************************************************
void groupLabels(Mat colorimg, std::vector<std::pair<int, int>> connectCC1, std::vector<std::pair<int, int>> connectCC2, int **label, int **box) {
	std::vector<int> labelList;
	std::vector<std::pair<int, int>> connectCC = connectCC1;
	connectCC.insert(connectCC.end(), connectCC2.begin(), connectCC2.end());
	std::vector<int> flag;

	for (int i = 0; i<int(connectCC.size()); i++) {
		int f = 0;

		//check if already connected...
		for (int k = 0; k<int(flag.size()); k++) {
			if (connectCC[i].first == flag[k] || connectCC[i].second == flag[k]) {
				f = 1;
				break;
			}
		}
		if (f) continue;

		//list the connected labels...
		labelList.push_back(connectCC[i].first);
		labelList.push_back(connectCC[i].second);
		int start_t = 0;
		while (start_t <= int(labelList.size()) - 1) {
			for (int t = start_t; t<int(labelList.size()); t++) {
				int c = 0;
				for (int k = 0; k<int(connectCC.size()); k++) {
					if (k == i)
						continue;
					if (connectCC[k].first == labelList[t] ) {
						labelList.push_back(connectCC[k].second);
						c++;
					}
					if (connectCC[k].second == labelList[t]) {
						labelList.push_back(connectCC[k].first);
						c++;
					}
				}
			}
			start_t = int(labelList.size()) - c;
		}	
		

		//remove repeatitions...
		for (int k = 0; k<int(labelList.size()); k++) {
			for (int t = k + 1; t<int(labelList.size()); t++) {
				if (labelList[t] == labelList[k]) {
					labelList.erase(labelList.begin() + t);
				}
			}
		}
		flag.insert(flag.end(), labelList.begin(), labelList.end());

		//draw common boundary box...
		int x0 = 99999, x1 = 0, y0 = 99999, y1 = 0;
		for (int k = 0; k<int(labelList.size()); k++) {
			int lbl = labelList[k];
			if (box[lbl][2] < x0)
				x0 = box[lbl][2];
			if (box[lbl][3] > x1)
				x1 = box[lbl][3];
			if (box[lbl][4] < y0)
				y0 = box[lbl][4];
			if (box[lbl][5] > y1)
				y1 = box[lbl][5];
		}
		/*line(colorimg, Point(x0, y0), Point(x0, y1), Scalar(0, 0, 255), 2, 8);
		line(colorimg, Point(x0, y1), Point(x1, y1), Scalar(0, 0, 255), 2, 8);
		line(colorimg, Point(x0, y0), Point(x1, y0), Scalar(0, 0, 255), 2, 8);
		line(colorimg, Point(x1, y0), Point(x1, y1), Scalar(0, 0, 255), 2, 8);*/


		labelList.clear();
	}


	//imwrite("Data/ICDAR/Data155/CommonBB.tif", colorimg);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/CommonBB.tif", colorimg);

}

//----------------------------------LEVEL 1: Combine cc's to form "words" and find orientations....
std::vector<pair<pair<int, int>, int>> combineLabels(Mat colorimg, Mat text_lines, std::vector<std::pair<int, int>> connectCC1, std::vector<std::pair<int, int>> connectCC2, 
	                                                std::vector<std::pair<int, Point>> all_centroids, int **label, int **box, std::vector<std::vector<int>> connectLabelsAngles, 
	                                                     std::vector<std::vector<int>> &TextLines, std::vector<std::pair<int, double>> &lineRatio) {

	std::vector<std::pair<int, int>> labelList;
	std::vector<std::pair<int, int>> connectCC = connectCC1;
	connectCC.insert(connectCC.end(), connectCC2.begin(), connectCC2.end());
	std::vector<int> flag;

	//list all the labels which are connected...
	for (int i = 0; i<int(connectCC.size()); i++) {
		int f1 = 0, f2 = 0;

		//check if already included...
		for (int k = 0; k<int(labelList.size()); k++) {
			if (connectCC[i].first == labelList[k].first) {
				f1 = 1;
				break;
			}
		}
		if (!f1)
			labelList.push_back(std::pair<int, int>(connectCC[i].first, 0));

		for (int k = 0; k<int(labelList.size()); k++) {
			if (connectCC[i].second == labelList[k].first) {
				f2 = 1;
				break;
			}
		}
		if (!f2)
			labelList.push_back(std::pair<int, int>(connectCC[i].second, 0));
	}

	int count = 0;
	for (int i = 0; i<int(connectCC.size()); i++) {
		int lbl1, lbl2;
		int f1 = 0, f2 = 0;
		for (int k = 0; k<int(labelList.size()); k++) {
			if (labelList[k].first == connectCC[i].first) {
				lbl1 = k;
				f1 = 1;
			}
			if (labelList[k].first == connectCC[i].second) {
				lbl2 = k;
				f2 = 1;
			}
			if (f1 && f2) break;
		}

		//if both are not labelled yet...
		if (labelList[lbl1].second == 0 && labelList[lbl2].second == 0) {
			count++;
			labelList[lbl1].second = count;
			labelList[lbl2].second = count;
			continue;
		}

		//if one is labelled...
		if (labelList[lbl1].second == 0 && labelList[lbl2].second != 0) {
			labelList[lbl1].second = labelList[lbl2].second;
			continue;
		}
		if (labelList[lbl1].second != 0 && labelList[lbl2].second == 0) {
			labelList[lbl2].second = labelList[lbl1].second;
			continue;
		}

		//if both are already labelled but differently...
		if (labelList[lbl1].second != 0 && labelList[lbl2].second != 0 && labelList[lbl1].second != labelList[lbl2].second) {
			int min = labelList[lbl1].second < labelList[lbl2].second ? labelList[lbl1].second : labelList[lbl2].second;
			int max= labelList[lbl1].second > labelList[lbl2].second ? labelList[lbl1].second : labelList[lbl2].second;
			for (int k = 0; k<int(labelList.size()); k++) {
				if (labelList[k].second == max)
					labelList[k].second = min;
			}
			continue;
		}
	}

	//draw new common boundary box...
	int **new_box = new int*[count];
	for (int i = 0; i < count; i++) {
		new_box[i] = new int[4];
		for (int j = 0; j < 4; j++)
			new_box[i][j] = 0;
	}
	std::vector<std::vector<int>> group_labels;
	for (int i = 1; i <= count; i++) {
		int x0 = 99999, x1 = 0, y0 = 99999, y1 = 0;
		std::cout << "\nJoining labels with count = " << i;
		std::vector<int> labelList_i;
		for (int k = 0; k<int(labelList.size()); k++) {
			if (labelList[k].second != i) continue;
			std::cout << labelList[k].first << " , ";
			labelList_i.push_back(labelList[k].first);
			int lbl = labelList[k].first;
			if (box[lbl][2] < x0)
				x0 = box[lbl][2];
			if (box[lbl][3] > x1)
				x1 = box[lbl][3];
			if (box[lbl][4] < y0)
				y0 = box[lbl][4];
			if (box[lbl][5] > y1)
				y1 = box[lbl][5];
		}
		group_labels.push_back(labelList_i);
		labelList_i.clear();

		new_box[i - 1][0] = x0;
		new_box[i - 1][1] = x1;
		new_box[i - 1][2] = y0;
		new_box[i - 1][3] = y1;

		line(colorimg, Point(x0, y0), Point(x0, y1), Scalar(0, 0, 255), 2, 8);
		line(colorimg, Point(x0, y1), Point(x1, y1), Scalar(0, 0, 255), 2, 8);
		line(colorimg, Point(x0, y0), Point(x1, y0), Scalar(0, 0, 255), 2, 8);
		line(colorimg, Point(x1, y0), Point(x1, y1), Scalar(0, 0, 255), 2, 8);
	}
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/CommonBB.tif", colorimg);

	//labellist...
	for (int k = 0; k<int(labelList.size()); k++) {
		std::cout << labelList[k].first << " = " << labelList[k].second << "\n";
	}

	//HPP for each new boxes...
	std::vector<std::pair<int, double>> p2bratios; //angle, ratio...
	std::vector<std::vector<int>> list_data; // label, angle, width of base, centroid_x, centroid_y...
	std::vector<std::vector<Point>> lines_data; //list of points to draw lines for each block...
	std::vector<int> ref;

	std::cout << "\n\nList of common labels:\n";
	for (int i = 0; i <= count - 1; i++) {

		std::cout << "\nLabel " << i + 1 << " : ";
		if (!group_labels[i].size()) {
			std::cout << "skip!";
			continue;
		}

		ref.push_back(i);
		for (int t = 0; t<int(group_labels[i].size()); t++) {
			std::cout << group_labels[i][t] << " , ";
		}

		std::vector<pair<int, int>> connections;
		std::cout << "\nConnections :";
		for (int t = 0; t<int(group_labels[i].size()); t++) {
			for (int k = 0; k<int(connectCC.size()); k++) {
				if (connectCC[k].first == group_labels[i][t] || connectCC[k].second == group_labels[i][t]) {
					//std::cout << connectCC[k].first << "<-->" << connectCC[k].second << " , ";
					int flag = 1;
					for (int l = 0; l<int(connections.size()); l++) {
						if ((connections[l].first == connectCC[k].first && connections[l].second == connectCC[k].second)
							|| (connections[l].first == connectCC[k].second && connections[l].second == connectCC[k].first)) {
							flag = 0;
							break;
						}
					}
					if (flag)
						connections.push_back(connectCC[k]);
				}
			}
		}
		for (int k = 0; k<int(connections.size()); k++)
			std::cout << connections[k].first << "<-->" << connections[k].second << " , ";
		std::cout << "\n-------------";

		//new ROIs...
		int h = new_box[i][3] - new_box[i][2] + 1;
		int w = new_box[i][1] - new_box[i][0] + 1;
		int **roi = new int*[h];
		for (int k = 0; k < h; k++) {
			roi[k] = new int[w];
			for (int l = 0; l < w; l++) {
				int lbl = label[new_box[i][2] + k][new_box[i][0] + l];
				int flag = 0;
				for (int t = 0; t<int(group_labels[i].size()); t++) {
					if (group_labels[i][t] == lbl) {
						flag = 1;
						break;
					}
				}
				if (flag) 
					roi[k][l] = 0;
				else
					roi[k][l] = 255;
			}
		}
		Mat roi_img;
		roi_img.create(h, w, CV_8UC1);
		//imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/roi.tif", roi_img);
		
		int base = 0;
		std::vector<Point> lines;
		std::pair<int, double> ratio = findOrientation(roi_img, roi, h, w, connections, all_centroids, base, lines, connectLabelsAngles);
		std::cout << "\n*************   Best orientation = " << ratio.first << " degrees  ***************";

		p2bratios.push_back(ratio);

		for (int k = 0; k<int(lines.size()); k++) {
			lines[k] += Point(new_box[i][0], new_box[i][2]);
			//circle(colorimg, lines[k] , 4, Scalar(255, 0, 0), 2);
		}
		lines_data.push_back(lines);

		std::vector<int> data; //label, angle, width of base, centroid_x, centroid_y...
		for (int t = 0; t<int(group_labels[i].size()); t++) {
			data.push_back(group_labels[i][t]); //label
			data.push_back(ratio.first); //angle
			data.push_back(base);   //width of base
			for (int k = 0; k<int(all_centroids.size()); k++) {
				if (all_centroids[k].first == group_labels[i][t]) { //centroid
					data.push_back(all_centroids[k].second.x);
					data.push_back(all_centroids[k].second.y);
					break;
				}
			}
		}
		list_data.push_back(data);
	}		

	std::cout << "\n________________________________\nList of data:";
	for (int i = 0; i<int(list_data.size()); i++) {
		int lbl = list_data[i][0];
		int angle = list_data[i][1];
		int base = list_data[i][2];
		int x = list_data[i][3];
		int y = list_data[i][4];

		std::cout << "\nlabel " << lbl << " is at angle " << angle << " degrees with base width " << base << " , centroid at " << Point(x, y);
		Point p1, p2;
		//p1...
		double slope = tan(angle*CV_PI / 180);
		p1.x = box[lbl][2];
		p1.y = int(y + slope*(p1.x - x));
		if (p1.y < box[lbl][4]) {
			p1.y = 0;
			p1.x = int(x - y / slope);
		}
		if (p1.y > box[lbl][5] - 1) {
			p1.y = box[lbl][5] - 1;
			p1.x = int(x - y / slope);
		}
		//p2...
		p2.x = box[lbl][3];
		p2.y = int(y + slope*(p2.x - x));
		if (p2.y < box[lbl][4]) {
			p2.y = 0;
			p2.x = int(x - y / slope);
		}
		if (p2.y > box[lbl][5] - 1) {
			p2.y = box[lbl][5] - 1;
			p2.x = int(x - y / slope);
		}
		//line(colorimg, p1, p2, Scalar(255, 0, 0), 3, 8);
	}

	
	std::vector<double> ratio_list;
	for (int i = 0; i < int(p2bratios.size()); i++) {
		ratio_list.push_back(p2bratios[i].second);
	}
	double med = median(ratio_list, int(p2bratios.size()));
	std::cout << "\n\n\nMedian of ratios = " << med;

	std::vector<int> labels_list;
	std::vector<pair<pair<int, int>, int>> label_deg; //<label, skew, line#>...
	int count_lineNo = 0, lineNo = 0;
	for (int i = 0; i < int(p2bratios.size()); i++) {
		std::cout << "\n" << i + 1 << ". orientation = " << p2bratios[i].first << " degrees , ratio = " << p2bratios[i].second;
		if (p2bratios[i].second > 0.5* med) {
			std::cout << "-----okay";
			count_lineNo++;
			
			//collect the labels that are getting marked...
			std::cout << "\n Labels with lines:" ;

			for (int k = 0; k<int(group_labels[ref[i]].size()); k++) {
				std::cout << group_labels[ref[i]][k] << " , ";
				label_deg.push_back(make_pair(make_pair(group_labels[ref[i]][k], p2bratios[i].first), count_lineNo));
			}
		
			if (p2bratios[i].first == 90 || p2bratios[i].first == -90) {
				int it = int(TextLines.size()) ;
				for (int k = 0; k<int(lines_data[i].size()); k++) {
					std::vector<int> newLine;
					int y1 = new_box[ref[i]][2] - 15;
					int y2 = new_box[ref[i]][3] + 15;
					int x = lines_data[i][k].x;
				//	line(text_lines, Point(x, y1), Point(x, y2), Scalar(255, 0, 0), 2, 8);
					line(colorimg, Point(x, y1), Point(x, y2), Scalar(255, 0, 0), 3, 8);
					circle(colorimg, lines_data[i][k], 4, Scalar(255, 0, 255), 2);
										
					newLine.push_back(++lineNo);
					if (y1 <= y2) {
						newLine.push_back(x);
						newLine.push_back(y1);
						newLine.push_back(x);
						newLine.push_back(y2);
					}
					else {
						newLine.push_back(x);
						newLine.push_back(y2);
						newLine.push_back(x);
						newLine.push_back(y1);
					}
					newLine.push_back(p2bratios[i].first);
					newLine.push_back(1);
					newLine.push_back(p2bratios[i].second);
					TextLines.push_back(newLine); // <line#, L_x, L_y, R_x, R_y, skew, type=1, ratio>...
					lineRatio.push_back(make_pair(lineNo, p2bratios[i].second)); //<line#, ratio>...
				}
			}
			//if nearly vertical then take top to bottom...
			else if ((p2bratios[i].first > 80 && p2bratios[i].first < 110) || (p2bratios[i].first > -110 && p2bratios[i].first < -80)) {
				double slope = tan(p2bratios[i].first*CV_PI / 180);

				for (int k = 0; k<int(lines_data[i].size()); k++) {
					std::vector<int> newLine;
					int x1 = new_box[ref[i]][0] - 15;
					int x2 = new_box[ref[i]][1] + 15;
					int y1 = lines_data[i][k].y + int((x1 - lines_data[i][k].x)*slope);
					if (y1 < new_box[ref[i]][2] - 15) {
						y1 = new_box[ref[i]][2] - 15;
						x1 = lines_data[i][k].x + int((y1 - lines_data[i][k].y) / slope);
					}
					if (y1 > new_box[ref[i]][3] + 15) {
						y1 = new_box[ref[i]][3] + 15;
						x1 = lines_data[i][k].x + int((y1 - lines_data[i][k].y) / slope);
					}
					int y2 = lines_data[i][k].y + int((x2 - lines_data[i][k].x)*slope);
					if (y2 < new_box[ref[i]][2] - 15) {
						y2 = new_box[ref[i]][2] - 15;
						x2 = lines_data[i][k].x + int((y2 - lines_data[i][k].y) / slope);
					}
					if (y2 > new_box[ref[i]][3] + 15) {
						y2 = new_box[ref[i]][3] + 15;
						x2 = lines_data[i][k].x + int((y2 - lines_data[i][k].y) / slope);
					}
					//line(text_lines, Point(x1, y1), Point(x2, y2), Scalar(255, 0, 0), 2, 8);
					line(colorimg, Point(x1, y1), Point(x2, y2), Scalar(255, 0, 0), 2, 8);
					circle(colorimg, lines_data[i][k], 4, Scalar(255, 0, 255), 2);

					newLine.push_back(++lineNo);
					//lines will be left to right & top to down...
					if (y1 <= y2) {
						newLine.push_back(x1);
						newLine.push_back(y1);
						newLine.push_back(x2);
						newLine.push_back(y2);
					}
					else {
						newLine.push_back(x2);
						newLine.push_back(y2);
						newLine.push_back(x1);
						newLine.push_back(y1);
					}
					newLine.push_back(p2bratios[i].first);
					newLine.push_back(1);
					newLine.push_back(p2bratios[i].second);
					TextLines.push_back(newLine); // <line#, L_x, L_y, R_x, R_y, skew, type=1, ratio>...
					lineRatio.push_back(make_pair(lineNo, p2bratios[i].second)); //<line#, ratio>...	
				}
			}
			else {
				double slope = tan(p2bratios[i].first*CV_PI / 180);

				for (int k = 0; k<int(lines_data[i].size()); k++) {
					std::vector<int> newLine;
					int x1 = new_box[ref[i]][0] - 15;
					int x2 = new_box[ref[i]][1] + 15;
					int y1 = lines_data[i][k].y + int((x1 - lines_data[i][k].x)*slope);
					if (y1 < new_box[ref[i]][2] - 15) {
						y1 = new_box[ref[i]][2] - 15;
						x1 = lines_data[i][k].x + int((y1 - lines_data[i][k].y) / slope);
					}
					if (y1 > new_box[ref[i]][3] + 15) {
						y1 = new_box[ref[i]][3] + 15;
						x1 = lines_data[i][k].x + int((y1 - lines_data[i][k].y) / slope);
					}
					int y2 = lines_data[i][k].y + int((x2 - lines_data[i][k].x)*slope);
					if (y2 < new_box[ref[i]][2] - 15) {
						y2 = new_box[ref[i]][2] - 15;
						x2 = lines_data[i][k].x + int((y2 - lines_data[i][k].y) / slope);
					}
					if (y2 > new_box[ref[i]][3] + 15) {
						y2 = new_box[ref[i]][3] + 15;
						x2 = lines_data[i][k].x + int((y2 - lines_data[i][k].y) / slope);
					}
					//line(text_lines, Point(x1, y1), Point(x2, y2), Scalar(255, 0, 0), 2, 8);
					line(colorimg, Point(x1, y1), Point(x2, y2), Scalar(255, 0, 0), 2, 8);
					circle(colorimg, lines_data[i][k], 4, Scalar(255, 0, 255), 2);

					newLine.push_back(++lineNo);
					//lines will be left to right & top to down...
					if (x1 < x2) {
						newLine.push_back(x1);
						newLine.push_back(y1);
						newLine.push_back(x2);
						newLine.push_back(y2);
					}
					else if (x1 == x2) {
						if (y1 < y2) {
							newLine.push_back(x1);
							newLine.push_back(y1);
							newLine.push_back(x2);
							newLine.push_back(y2);
						}
						else {
							newLine.push_back(x2);
							newLine.push_back(y2);
							newLine.push_back(x1);
							newLine.push_back(y1);
						}
					}
					else {
						newLine.push_back(x2);
						newLine.push_back(y2);
						newLine.push_back(x1);
						newLine.push_back(y1);
					}
					newLine.push_back(p2bratios[i].first);
					newLine.push_back(1);
					newLine.push_back(p2bratios[i].second);
					TextLines.push_back(newLine); // <line#, L_x, L_y, R_x, R_y, skew, type=1, ratio>...
					lineRatio.push_back(make_pair(lineNo, p2bratios[i].second)); //<line#, ratio>...	
				}

			}
		}
			
		else
			std::cout << "-----no";
	}
	
	//remove overlapping lines with common end points...
	for (int i = 1; i<int(TextLines.size()); i++) {
		int flg = 1;
		Point currL = Point(TextLines[i][1], TextLines[i][2]);
		Point currR = Point(TextLines[i][3], TextLines[i][4]);
		for (int k = 0; k < i; k++) {
			Point prevL = Point(TextLines[k][1], TextLines[k][2]);
			Point prevR = Point(TextLines[k][3], TextLines[k][4]);

			//both ends coincide ...
			if (abs(currL.x - prevL.x) < 10 && abs(currL.y - prevL.y) < 10 && abs(currR.x - prevR.x) < 10 && abs(currR.y - prevR.y) < 10) {
				flg = 0;
			}
			//left end coincides ...
			else if (abs(currL.x - prevL.x) < 10 && abs(currL.y - prevL.y) < 10) {
				if (norm(currL - currR) < norm(prevL - prevR) && abs(norm(currL - currR) + norm(currR - prevR) - norm(prevR - prevL)) < 5) {
					flg = 0;
				}
				else if (norm(currL - currR) >= norm(prevL - prevR) && abs(norm(currL - currR) - norm(currR - prevR) - norm(prevR - prevL)) < 5) {
					flg = 0;
				}
			}

			//right end coincides ...
			else if (abs(currR.x - prevR.x) < 10 && abs(currR.y - prevR.y) < 10) {
				if (norm(currL - currR) < norm(prevL - prevR) && abs(norm(currL - currR) + norm(currL - prevL) - norm(prevR - prevL)) < 5) {
					flg = 0;
				}
				else if (norm(currL - currR) >= norm(prevL - prevR) && abs(norm(currL - currR) - norm(currL - prevL) - norm(prevR - prevL)) < 5) {
					flg = 0;
				}
			}
			if (!flg) {
				std::cout << "\nLine " << i + 1 << " coincides with line " << k + 1;
				break;
			}
		}
		if (!flg) {		
			TextLines.erase(TextLines.begin() + i); // <line#, L_x, L_y, R_x, R_y, type=1, ratio>...
			lineRatio.erase(lineRatio.begin() + i); //<line#, ratio>...
			i--;
		}
		/*else
			line(text_lines, currL, currR, Scalar(255, 255, 0), 4, 8);	*/

	}
	std::cout << "\nList of textlines of type 1:\n";
	for (int i = 0; i<int(TextLines.size()); i++) {
		TextLines[i][0] = i + 1;
		lineRatio[i].first = i + 1;
		std::cout << "\n" << TextLines[i][0] << ". " << Point(TextLines[i][1], TextLines[i][2]) << " to " << Point(TextLines[i][3], TextLines[i][4])
			<< " skew = " << TextLines[i][5] << " -- Type " << TextLines[i][6] << "...ratio = " << lineRatio[i].second;
	}

	//display labels with orientations...
	sort(label_deg.begin(), label_deg.end());
	for (int i = 0; i < int(label_deg.size()); i++) {
		std::cout << "\nLabel " << label_deg[i].first.first << " orientation : " << label_deg[i].first.second << " Line no. : " << label_deg[i].second;
	}

	//imwrite("Data/ICDAR/Data155/CommonBB.tif", colorimg);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/CommonBB.tif", colorimg);
	
	return(label_deg);
}

//white pixels...
std::vector<std::pair<Point, Point>> combineBoxes(int **thin_img, int **label, int **box, std::vector<std::vector<Point>> all_terminals, 
	                           std::vector<std::pair<int, Point>> all_centroids, std::vector<std::pair<int, int>> &connectCC1, std::vector<std::pair<int, int>> &connectCC2) {

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
		double max_length, max_max_length = -1, max_max_length2 = -1;
		Point max_max_p1, max_max_p2;
		Point max_max2_p1, max_max2_p2;
		int  max_lbl = 0, max_lbl2 = 0, max_run = 0, min_run = 0;
		for (int k = 0; k<int(nbr_labels.size()); k++) {
			int nbr = nbr_labels[k];
			int max, min;
			Point max_p1, max_p2;
			max_length = -1;
			std::cout << "\nComparing with label " << nbr << "...";

			/*int flag = 0;
			for (int k = 0; k<int(connectCC1.size()); k++) {
				if ((connectCC1[k].first == i && connectCC1[k].second == nbr) || (connectCC1[k].first == nbr && connectCC1[k].second == i)) {
					flag = 1;
					break;
				}
			}
			if (flag) {
				std::cout << "\nAlready compared...";
				continue;
			}*/
			/*for (int k = 0; k<int(connectCC2.size()); k++) {
				if ((connectCC2[k].first == i && connectCC2[k].second == nbr) || (connectCC2[k].first == nbr && connectCC2[k].second == i)) {
					flag = 1;
					break;
				}
			}
			if (flag) {
				std::cout << "\nAlready compared...";
				continue;
			}*/

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
			else if (max_length > max_max_length2) {
				max_max_length2 = max_length;
				max_lbl2 = nbr;
				max_max2_p1 = max_p1;
				max_max2_p2 = max_p2;
			}

		}


		std::cout << "\nLabel " << i << " is near to label " << max_lbl << " with length = " << max_max_length;
		if (max_lbl != 0) {
			connectCC1.push_back(pair<int, int>(i, max_lbl));
			connectCentroid.push_back(pair<Point, Point>(max_max_p1, max_max_p2));
		}			
		if (max_lbl2 != 0) {
			connectCC2.push_back(pair<int, int>(i, max_lbl2));
			connectCentroid.push_back(pair<Point, Point>(max_max2_p1, max_max2_p2));
		}

		//empty the vecctor...
		nbr_labels.clear();
	}
	//return(connectCC);
	return(connectCentroid);
}

//white distance...
std::vector<std::pair<Point, Point>> combineBoxes_new(int **thin_img, int **label, int **box, std::vector<std::vector<Point>> all_terminals, 
	                                                  std::vector<std::pair<int, Point>> all_centroids, std::vector<std::pair<int, int>> &connectCC1, 
	                                                  std::vector<std::pair<int, int>> &connectCC2, std::vector<std::vector<int>> &connectLabelsAngles) {

	std::vector<std::pair<int, int>> matchedCC;  //store the labels that have been compared...
												 //store the labels to join...
	std::vector<std::pair<Point, Point>> connectCentroid; //store the actual points to connnect...
	//std::vector<std::vector<int>> connectLabelsAngles;  //store the orientation of the actual connections...


	int  dim = 12;
	double thresh = 0.3;
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
		double max_length, max_max_length = -1, max_max_length2 = -1;
		Point max_max_p1, max_max_p2;
		Point max_max2_p1, max_max2_p2;
		int  max_lbl = 0, max_lbl2 = 0;
		double max_run = 0, min_run = 0;
		for (int k = 0; k<int(nbr_labels.size()); k++) {
			int nbr = nbr_labels[k];
			double max, min;
			Point max_p1, max_p2;
			max_length = -1;
			std::cout << "\nComparing with label " << nbr << "...";

			/*int flag = 0;
			for (int k = 0; k<int(connectCC1.size()); k++) {
			if ((connectCC1[k].first == i && connectCC1[k].second == nbr) || (connectCC1[k].first == nbr && connectCC1[k].second == i)) {
			flag = 1;
			break;
			}
			}
			if (flag) {
			std::cout << "\nAlready compared...";
			continue;
			}*/
			/*for (int k = 0; k<int(connectCC2.size()); k++) {
			if ((connectCC2[k].first == i && connectCC2[k].second == nbr) || (connectCC2[k].first == nbr && connectCC2[k].second == i)) {
			flag = 1;
			break;
			}
			}
			if (flag) {
			std::cout << "\nAlready compared...";
			continue;
			}*/

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
			std::vector<double> whiteRuns;
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
			whiteRuns.clear();
			if (min/max >= thresh && min < 30) {
				length = double(countBW)*min / max;
				std::cout << "\nLength = " << length;
				max_length = length;
				max_run = max;
				min_run = min;
				max_p1 = p1 + Point(xstart, ystart);
				max_p2 = p2 + Point(xstart, ystart);
			}


			//lines from current label's centroid to neighbors...
			p1 = all_centroids[it_i].second - Point(xstart, ystart);
			std::vector<double> whiteRuns_i;
			for (int it = 1; it<int(all_terminals[it_nbr].size()); it++) {
				Point p2 = all_terminals[it_nbr][it] - Point(xstart, ystart);
				int countBW = countCrossings(roi, H, W, p1, p2, whiteRuns_i);
				std::cout << "\nLine " << all_centroids[it_i].second << " to " << all_terminals[it_nbr][it];
				std::cout << "\nThe label " << nbr << " has " << countBW << " B/W crossings and the sequence of white run-lengths are: ";
				for (int seq = 0; seq<int(whiteRuns_i.size()); seq++)
					std::cout << whiteRuns_i[seq] << "...";

				if (int(whiteRuns_i.size()) == 0) continue;
				//find min & max run length...
				max = whiteRuns_i[0];
				min = whiteRuns_i[0];
				for (int seq = 1; seq<int(whiteRuns_i.size()); seq++) {
					if (whiteRuns_i[seq] > max)
						max = whiteRuns_i[seq];
					else if (whiteRuns_i[seq] < min)
						min = whiteRuns_i[seq];
				}
				whiteRuns_i.clear();
				if (min/max < thresh || min > 30)
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
			}

			//lines from neighbors' centroid to current label...
			p1 = all_centroids[it_nbr].second - Point(xstart, ystart);
			std::vector<double> whiteRuns_nbr;
			for (int it = 1; it<int(all_terminals[it_i].size()); it++) {
				Point p2 = all_terminals[it_i][it] - Point(xstart, ystart);
				int countBW = countCrossings(roi, H, W, p1, p2, whiteRuns_nbr);
				std::cout << "\nLine " << all_centroids[it_nbr].second << " to " << all_terminals[it_i][it];
				std::cout << "\nThe label " << i << " has " << countBW << " B/W crossings and the sequence of white run-lengths are: ";
				for (int seq = 0; seq<int(whiteRuns_nbr.size()); seq++)
					std::cout << whiteRuns_nbr[seq] << "...";

				if (int(whiteRuns_nbr.size()) == 0) continue;
				//find min & max run length...
				max = whiteRuns_nbr[0];
				min = whiteRuns_nbr[0];
				for (int seq = 1; seq<int(whiteRuns_nbr.size()); seq++) {
					if (whiteRuns_nbr[seq] > max)
						max = whiteRuns_nbr[seq];
					else if (whiteRuns_nbr[seq] < min)
						min = whiteRuns_nbr[seq];
				}
				whiteRuns_nbr.clear();
				if (min/max < thresh || min > 30)
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
			}

			//lines from currentlabel's terminal points to neighbours' terminal points...
			for (int it1 = 1; it1<int(all_terminals[it_i].size()); it1++) {
				p1 = all_terminals[it_i][it1] - Point(xstart, ystart);
				std::vector<double> whiteRuns_i;

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
					whiteRuns_i.clear();
					if (min/max < thresh || min > 30)
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
				}
			}



			if (max_length > max_max_length) {
				max_max_length = max_length;
				max_lbl = nbr;
				max_max_p1 = max_p1;
				max_max_p2 = max_p2;
			}
			else if (max_length > max_max_length2) {
				max_max_length2 = max_length;
				max_lbl2 = nbr;
				max_max2_p1 = max_p1;
				max_max2_p2 = max_p2;
			}

		}

		std::vector<int> connectOrient;

		std::cout << "\nLabel " << i << " is near to label " << max_lbl << " with length = " << max_max_length;
		if (max_lbl != 0) {
			connectCC1.push_back(pair<int, int>(i, max_lbl));
			connectCentroid.push_back(pair<Point, Point>(max_max_p1, max_max_p2));
			connectOrient.push_back(i);
			connectOrient.push_back(max_lbl);
			connectOrient.push_back(int(atan2(max_max_p1.y - max_max_p2.y, max_max_p1.x - max_max_p2.x) * 180 / CV_PI));
			connectLabelsAngles.push_back(connectOrient);
			connectOrient.clear();
		}
		if (max_lbl2 != 0) {
			connectCC2.push_back(pair<int, int>(i, max_lbl2));
			connectCentroid.push_back(pair<Point, Point>(max_max2_p1, max_max2_p2));
			connectOrient.push_back(i);
			connectOrient.push_back(max_lbl2);
			connectOrient.push_back(int(atan2(max_max2_p1.y - max_max2_p2.y, max_max2_p1.x - max_max2_p2.x) * 180 / CV_PI));
			connectLabelsAngles.push_back(connectOrient);
			connectOrient.clear();
		}

		//empty the vecctor...
		nbr_labels.clear();
	}
	//return(connectCC);
	return(connectCentroid);
}

//--------------------------------LEVEL 1: Find orientations of "single words"....
std::vector<pair<pair<int, int>, int>> singleLabels(Mat colorimg, Mat text_lines, std::vector<std::pair<int, Point>> all_centroids, 
	                                                std::vector<std::vector<Point>> all_terminals, int **label, int **box, std::vector<int> unmarked_labels,
											        std::vector<std::vector<int>> &TextLines, std::vector<std::pair<int, double>> &lineRatio,
	                                                std::vector<pair<pair<int, int>, pair<int, double>>> &label2line) {
	std::vector<std::pair<int, double>> p2bratios; //angle, ratio...
	std::vector<std::vector<int>> list_data; // label, angle, width of base, centroid_x, centroid_y...
	std::vector<Point> lines_data; //list of points to draw lines for each block...
	
	for (int i = 1; i < labelNum; i++) {
		if (box[i][2] == -1) continue;

		int flag = 0;
		for (int k = 0; k<int(unmarked_labels.size()); k++) {
			if (i == unmarked_labels[k]) {
				flag = 1;
				break;
			}
		}
		if (!flag) continue;

		//if block has less than 3 terminal points then skip...
		int ref = 0;
		while (all_centroids[ref].first != i && ref < int(all_centroids.size())) {
			ref++;
		}
		
		if (int(all_terminals[ref].size()) <= 4 )
			continue;

		//roi...
		int **roi, **thin_roi;
		int istart = box[i][4];
		int jstart = box[i][2];
		int h = box[i][5] - box[i][4] + 1;
		int w = box[i][3] - box[i][2] + 1;

		roi = new int*[h];
		for (int k = 0; k < h; k++) {
			roi[k] = new int[w];
			for (int l = 0; l < w; l++) {
				if (label[k + istart][l + jstart] == i)
					roi[k][l] = 0;
				else
					roi[k][l] = 255;
			}
		}

		Mat roi_img;
		roi_img.create(h, w, CV_8UC1);
		//imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/roi.tif", roi_img);

		int base = 0;
		std::vector<Point> lines;
		std::pair<int, double> ratio = findOrientation(roi_img, roi, h, w, all_centroids, all_terminals[ref], base, lines);
		std::cout << "\n*************   Best orientation = " << ratio.first << " degrees  ***************";

		p2bratios.push_back(ratio);

		//for (int k = 0; k<int(lines.size()); k++) {
		//	lines[k] += Point(box[i][2], box[i][4]);
		//	//circle(colorimg, lines[k] , 4, Scalar(255, 0, 0), 2);
		//}
		lines_data.push_back(all_centroids[ref].second);

		std::vector<int> data; //label, angle, width of base, centroid_x, centroid_y...
		data.push_back(all_centroids[ref].first);//label
		data.push_back(ratio.first); //angle
		data.push_back(base);   //width of base
		data.push_back(all_centroids[ref].second.x);
		data.push_back(all_centroids[ref].second.y);
		list_data.push_back(data);


	}

	std::cout << "\n________________________________\nList of data:";
	for (int i = 0; i<int(list_data.size()); i++) {
		int lbl = list_data[i][0];
		int angle = list_data[i][1];
		int base = list_data[i][2];
		int x = list_data[i][3];
		int y = list_data[i][4];

		std::cout << "\nlabel " << lbl << " is at angle " << angle << " degrees with base width " << base << " , centroid at " << Point(x, y);
		Point p1, p2;
		//p1...
		double slope = tan(angle*CV_PI / 180);
		p1.x = box[lbl][2];
		p1.y = int(y + slope*(p1.x - x));
		if (p1.y < box[lbl][4]) {
			p1.y = 0;
			p1.x = int(x - y / slope);
		}
		if (p1.y > box[lbl][5] - 1) {
			p1.y = box[lbl][5] - 1;
			p1.x = int(x - y / slope);
		}
		//p2...
		p2.x = box[lbl][3];
		p2.y = int(y + slope*(p2.x - x));
		if (p2.y < box[lbl][4]) {
			p2.y = 0;
			p2.x = int(x - y / slope);
		}
		if (p2.y > box[lbl][5] - 1) {
			p2.y = box[lbl][5] - 1;
			p2.x = int(x - y / slope);
		}
		//line(colorimg, p1, p2, Scalar(255, 255, 0), 3, 8);
	}

	std::vector<double> ratio_list;
	for (int i = 0; i < int(p2bratios.size()); i++) {
		ratio_list.push_back(p2bratios[i].second);
	}
	double med = median(ratio_list, int(p2bratios.size()));
	std::cout << "\n\n\nMedian of ratios = " << med;

	std::vector<int> labels_list;
	std::vector<pair<pair<int, int>,int>> label_deg;
	int count_LineNo = TextLines.back()[0], lineNo = TextLines.back()[0];
	for (int i = 0; i < int(p2bratios.size()); i++) {
		std::cout << "\n" << i + 1 << ". orientation = " << p2bratios[i].first << " degrees , ratio = " << p2bratios[i].second;
		if (p2bratios[i].second > 0.5* med) {
			std::cout << "-----okay";
			std::vector<int> newLine;

			//collect the labels that are getting marked...
			std::cout << "\n Labels with lines:"; // <label, skew, line#> 
			label_deg.push_back(make_pair(make_pair(list_data[i][0], p2bratios[i].first), ++count_LineNo));
					
			int lbl = list_data[i][0];
			int angle = list_data[i][1];
			int base = list_data[i][2];
			int x = list_data[i][3];
			int y = list_data[i][4];

			if (angle == 90 || angle == -90) {
				int y1 = box[lbl][2] - 15;
				int y2 = box[lbl][3] + 15;
				int x = lines_data[i].x;
				line(colorimg, Point(x, y1), Point(x, y2), Scalar(255, 0, 0), 3, 8);
				circle(colorimg, lines_data[i], 4, Scalar(255, 0, 255), 2);
				newLine.push_back(++lineNo);
				if (y1 <= y2) {
					newLine.push_back(x);
					newLine.push_back(y1);
					newLine.push_back(x);
					newLine.push_back(y2);
				}
				else {
					newLine.push_back(x);
					newLine.push_back(y2);
					newLine.push_back(x);
					newLine.push_back(y1);
				}
				newLine.push_back(p2bratios[i].first);
				newLine.push_back(2);
				newLine.push_back(p2bratios[i].second);
				TextLines.push_back(newLine); // <line#, L_x, L_y, R_x, R_y, skew, type=2, ratio>...
				lineRatio.push_back(make_pair(lineNo, p2bratios[i].second)); //<line no. ratio>...
			}
			//if nearly vertical then take top to bottom...
			else if ((angle > 80 && angle < 110) || (angle > -110 && angle < -80)) {
				double slope = tan(angle*CV_PI / 180);

				int x1 = box[lbl][2] - 15;

				int y1 = lines_data[i].y + int((x1 - lines_data[i].x)*slope);
				if (y1 < box[lbl][4] - 15) {
					y1 = box[lbl][4] - 15;
					x1 = lines_data[i].x + int((y1 - lines_data[i].y) / slope);
				}
				if (y1 > box[lbl][5] + 15) {
					y1 = box[lbl][5] + 15;
					x1 = lines_data[i].x + int((y1 - lines_data[i].y) / slope);
				}

				int x2 = box[lbl][3] + 15;

				int y2 = lines_data[i].y + int((x2 - lines_data[i].x)*slope);
				if (y2 <  box[lbl][4] - 15) {
					y2 = box[lbl][4] - 15;
					x2 = lines_data[i].x + int((y2 - lines_data[i].y) / slope);
				}
				if (y2 >  box[lbl][5] + 15) {
					y2 = box[lbl][5] + 15;
					x2 = lines_data[i].x + int((y2 - lines_data[i].y) / slope);
				}
				line(colorimg, Point(x1, y1), Point(x2, y2), Scalar(120, 200, 0), 3, 8);
				line(text_lines, Point(x1, y1), Point(x2, y2), Scalar(120, 200, 0), 3, 8);
				circle(colorimg, lines_data[i], 4, Scalar(255, 0, 255), 2);
				newLine.push_back(++lineNo);
				//lines will be left to right & top to down...
				if (y1 <= y2) {
					newLine.push_back(x1);
					newLine.push_back(y1);
					newLine.push_back(x2);
					newLine.push_back(y2);
				}
				else {
					newLine.push_back(x2);
					newLine.push_back(y2);
					newLine.push_back(x1);
					newLine.push_back(y1);
				}
				newLine.push_back(p2bratios[i].first);
				newLine.push_back(2);
				newLine.push_back(p2bratios[i].second);
				TextLines.push_back(newLine); // <line#, L_x, L_y, R_x, R_y, skew, type=2, ratio>...
				lineRatio.push_back(make_pair(lineNo, p2bratios[i].second)); //<line no. ratio>...
			}
			else {
				double slope = tan(angle*CV_PI / 180);

				int x1 = box[lbl][2] - 15;
				
				int y1 = lines_data[i].y + int((x1 - lines_data[i].x)*slope);
				if (y1 < box[lbl][4] - 15) {
					y1 = box[lbl][4] - 15;
					x1 = lines_data[i].x + int((y1 - lines_data[i].y) / slope);
				}
				if (y1 > box[lbl][5] + 15) {
					y1 = box[lbl][5] + 15;
					x1 = lines_data[i].x + int((y1 - lines_data[i].y) / slope);
				}

				int x2 = box[lbl][3] + 15;

				int y2 = lines_data[i].y + int((x2 - lines_data[i].x)*slope);
				if (y2 <  box[lbl][4] - 15) {
					y2 = box[lbl][4] - 15;
					x2 = lines_data[i].x + int((y2 - lines_data[i].y) / slope);
				}
				if (y2 >  box[lbl][5] + 15) {
					y2 = box[lbl][5] + 15;
					x2 = lines_data[i].x + int((y2 - lines_data[i].y) / slope);
				}
				line(colorimg, Point(x1, y1), Point(x2, y2), Scalar(120, 200, 0), 3, 8);
				line(text_lines, Point(x1, y1), Point(x2, y2), Scalar(120, 200, 0), 3, 8);
				circle(colorimg, lines_data[i], 4, Scalar(255, 0, 255), 2);
				newLine.push_back(++lineNo);
				//lines will be left to right & top to down...
				if (x1 < x2) {
					newLine.push_back(x1);
					newLine.push_back(y1);
					newLine.push_back(x2);
					newLine.push_back(y2);
				}
				else if (x1 == x2) {
					if (y1 < y2) {
						newLine.push_back(x1);
						newLine.push_back(y1);
						newLine.push_back(x2);
						newLine.push_back(y2);
					}
					else {
						newLine.push_back(x2);
						newLine.push_back(y2);
						newLine.push_back(x1);
						newLine.push_back(y1);
					}
				}
				else {
					newLine.push_back(x2);
					newLine.push_back(y2);
					newLine.push_back(x1);
					newLine.push_back(y1);
				}
				newLine.push_back(p2bratios[i].first);
				newLine.push_back(2);
				newLine.push_back(p2bratios[i].second);
				TextLines.push_back(newLine); // <line#, L_x, L_y, R_x, R_y, skew, type=2, ratio>...
				lineRatio.push_back(make_pair(lineNo, p2bratios[i].second)); //<line no. ratio>...
			}
			
		}

		else
			std::cout << "-----no";
	}

	for (int i = 0; i<int(all_centroids.size()); i++) {
		circle(text_lines, all_centroids[i].second, 4, Scalar(0, 0, 255), 2);
	}

	//display labels with orientations...
	sort(label_deg.begin(), label_deg.end());
	for (int i = 0; i < int(label_deg.size()); i++) {
		std::cout << "\nLabel " << label_deg[i].first.first << " orientation : " << label_deg[i].first.second << " Line No. : " << label_deg[i].second;
		label2line.push_back(make_pair(make_pair(label_deg[i].first.first, label_deg[i].second), make_pair(label_deg[i].first.second, 0)));
	}

	//imwrite("Data/ICDAR/Data155/Orientations.tif", colorimg);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/Orientations.tif", colorimg);

	//imwrite("Data/ICDAR/Data155/text_lines.tif", text_lines);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/text_lines.tif", text_lines);

	return(label_deg);
}

double linearFit(std::vector<Point> data) {
	double X_mean = 0.0, Y_mean = 0.0, SXY = 0.0, SXX = 0.0, SYY = 0.0;
	int N = int(data.size());
	for (int i = 0; i < N; i++) {
		X_mean += data[i].x;
		Y_mean += data[i].y;
	}
	X_mean /= N;
	Y_mean /= N;
	for (int i = 0; i < N; i++) {
		SXX += data[i].x*data[i].x;
		SYY += data[i].y*data[i].y;
		SXY += data[i].x*data[i].y;
	}
	double slope = ((SYY - SXX) + sqrt((SYY - SXX)*(SYY - SXX) + 4 * SXY*SXY)) / (2 * SXY);
	double skew = atan(slope) * 180 / CV_PI;
	if (skew > 90 && skew < 180) skew = skew - 180;
	else if (skew < -90 && skew > -180) skew = skew + 180;
	std::cout << "\nSlope of linear fitting line with perpendicular offset :" << skew << " degrees";
	
	return(skew);
}

int LSM_perp_off(std::vector<Point> data, double& X_mean, double& Y_mean) {
	double SXY = 0.0, SXX = 0.0, SYY = 0.0;
	int N = int(data.size());
	for (int i = 0; i < N; i++) {
		X_mean += data[i].x;
		Y_mean += data[i].y;
	}
	X_mean /= N;
	Y_mean /= N;
	for (int i = 0; i < N; i++) {
		SXX += data[i].x*data[i].x;
		SYY += data[i].y*data[i].y;
		SXY += data[i].x*data[i].y;
	}
	double B = 0.5*((SYY - N*Y_mean*Y_mean) - (SXX - N*X_mean*X_mean)) / (N*X_mean*Y_mean - SXY);
	double slope = -B + sqrt(B*B + 1);
	int angle = atan(slope) * 180 / CV_PI;
	std::cout << "\nOrientation of LSM with perpendicular offset = " << angle;
	return(angle);
}

int medianSkew(std::vector<Point> data) {
	int med;
	std::vector<int> skew;
	for (int i = 0; i<int(data.size()); i++) {
		for (int j = i + 1; j<int(data.size()); j++) {
			skew.push_back(atan2(data[i].y - data[j].y, data[i].x - data[j].x) * 180 / CV_PI);
			if (skew.back() < 0)
				skew.back() += 180;
			std::cout << skew.back() << ",";
		}
	}
	
	return(int(median(skew, int(skew.size()))));
}

//-------------------------------LEVEL 2: Connect pair of "words" to form "line segments"....
Point linesIntersect(std::pair<int, Point> line1, std::pair<int, Point> line2) {
	double slope1 = tan(line1.first*CV_PI / 180);
	double slope2 = tan(line2.first*CV_PI / 180);
	int X = double((slope1*line1.second.x - slope2*line2.second.x) - (line1.second.y - line2.second.y)) / (slope1 - slope2);
	int Y = line1.second.y + slope1*(X - line1.second.x);
	return(Point(X, Y));
}

void connectLines(Mat text_lines, std::vector<std::pair<int, Point>> all_centroids, std::vector<std::vector<Point>> all_terminals, int **label, int **box,
	                                                   std::vector<std::vector<int>> &TextLines, std::vector<std::pair<int, double>> &lineRatio,
	                                                   std::vector<std::vector<int>> &TextLinesNew, 
	                                                   std::vector<pair<pair<int, int>, pair<int, double>>> &label2line, 
	                                                   std::vector<std::pair<std::pair<int, int>, std::pair<double, int>>> &ConnectedLines) {

	std::vector<std::pair<std::pair<int, int>, std::pair<double , double>>> line_pairs; //line #, line #, distance, ratio...
	std::vector<std::pair<std::pair<int, int>, double>> line_pairs2; //line #, line #, distance...
	std::vector<std::pair<int, double>> p2bratios; //angle, ratio...
	std::vector<std::vector<Point>> lines_data; // points to draw lines...

	int count = 0;
	int text_n = int(TextLines.size());

	int **newBox = new int*[text_n];
	for (int i = 0; i<text_n; i++) {
		newBox[i] = new int[4];
		for (int j = 0; j < 4; j++)
			newBox[i][j] = -1;
	}

	for (int i = 0; i<text_n; i++) {
		std::vector<int> collinear_labels;

		Point R = Point(TextLines[i][3], TextLines[i][4]);
		double min_dist = 9999, dist = 0;
		int min_line = 0;
		Point min_point;
		double min_dist2 = 9999, dist2 = 0;
		int min_line2 = 0;
		Point min_point2;
		for (int k = 0; k<text_n; k++) {
			if (k == i) continue;
			//print a connection from right end to nearest line's left end
			Point L_k = Point(TextLines[k][1], TextLines[k][2]);
			dist = norm(R - L_k);
			if (dist < min_dist) {
				min_dist = dist;
				min_line = k;
				min_point = L_k;
			}
			else if (dist < min_dist2) {
				min_dist2 = dist;
				min_line2 = k;
				min_point2 = L_k;
			}
			//std::cout << "\ndist between line " << i + 1 << " and line " << k + 1 << " = " << dist;
		}
		std::cout << "\nLine " << TextLines[i][0] << " is nearest to Line " << TextLines[min_line][0] << " at distance " << min_dist;
		line_pairs.push_back(make_pair(make_pair(TextLines[i][0], TextLines[min_line][0]), make_pair(min_dist, 0)));
		std::cout << "\nLine " << TextLines[i][0] << " is second nearest to Line " << TextLines[min_line2][0] << " at distance " << min_dist2;
		line_pairs2.push_back(make_pair(make_pair(TextLines[i][0], TextLines[min_line2][0]), min_dist2));
		std::cout << "\nCollinear points : ";
		std::vector<Point> data;
		data.push_back(Point(TextLines[i][1], TextLines[i][2])); //Left end of i-th line...
		data.push_back(Point(TextLines[i][3], TextLines[i][4])); //Right end of i-th line...
		data.push_back(Point(TextLines[min_line][1], TextLines[min_line][2])); //Left end of nearest line...
		data.push_back(Point(TextLines[min_line][3], TextLines[min_line][4])); //Right end of nearest line...


		//collect the combined labels...
		std::cout << "\nList of labels: ";
		for (int k = 0; k<int(label2line.size()); k++) {
			if (label2line[k].first.second == line_pairs.back().first.first || label2line[k].first.second == line_pairs.back().first.second) {
				std::cout << label2line[k].first.first << " , ";
				collinear_labels.push_back(label2line[k].first.first);
			}
		}

		//collect the centroids of combined labels...
		std::cout << "\nList of centroids: ";
		for (int k = 0; k<int(collinear_labels.size()); k++) {
			for (int t = 0; t<int(all_centroids.size()); t++) {
				if (all_centroids[t].first == collinear_labels[k]) {
					data.push_back(all_centroids[t].second);
					circle(text_lines, all_centroids[t].second, 4, Scalar(255, 0, 255), 2);
					std::cout << all_centroids[t].second << " , ";
				}
			}
		}

		//collect the skews...
		std::vector<int> skew;

		//include the original orientations of line pair in the list...
		skew.push_back(TextLines[i][0]);
		if (skew.back() > 90 && skew.back() < 180) skew.back() = skew.back() - 180;
		else if (skew.back() < -90 && skew.back() > -180) skew.back() = skew.back() + 180;

		skew.push_back(TextLines[min_line][0]);
		if (skew.back() > 90 && skew.back() < 180) skew.back() = skew.back() - 180;
		else if (skew.back() < -90 && skew.back() > -180) skew.back() = skew.back() + 180;
		if ((skew[0] + 3 > skew.back() && skew[0] - 3 < skew.back()) || (abs(skew[0] - skew.back()) > 180 - 3 && abs(skew[0] - skew.back()) < 180 + 3)) {
			skew.pop_back();
		}

		//skew between centroids of labels....
		for (int k = 0; k<int(data.size()) - 1; k++) {
			for (int kk = k + 1; kk < int(data.size()); kk++) {
				skew.push_back(atan2(data[k].y - data[kk].y, data[k].x - data[kk].x) * 180 / CV_PI);
				if (skew.back() > 90 && skew.back() < 180) skew.back() = skew.back() - 180;
				else if (skew.back() < -90 && skew.back() > -180) skew.back() = skew.back() + 180;
				for (int t = 0; t<int(skew.size()) - 1; t++) {
					if ((skew[t] + 3 > skew.back() && skew[t] - 3 < skew.back()) || (abs(skew[t] - skew.back()) > 180 - 3 && abs(skew[t] - skew.back()) < 180 + 3) ) {
						skew.pop_back();
						break;
					}
				}
			}
		}

		//skew between terminal points of labels...
		for (int k = 0; k<int(collinear_labels.size()) - 1; k++) {
			for (int t = k + 1; t<int(collinear_labels.size()); t++) {
				int lbl1, lbl2;
				for (int it = 0; it<int(all_terminals.size()); it++) {
					if (all_terminals[it][0].x == collinear_labels[k]) lbl1 = it;
					if (all_terminals[it][0].x == collinear_labels[t]) lbl2 = it;
				}
				for (int kk = 1; kk<int(all_terminals[lbl1].size()); kk++) {
					for (int tt = 1; tt < int(all_terminals[lbl2].size()); tt++) {
						skew.push_back(atan2(all_terminals[lbl1][kk].y - all_terminals[lbl2][tt].y, all_terminals[lbl1][kk].x - all_terminals[lbl2][tt].x) * 180 / CV_PI);
						if (skew.back() > 90 && skew.back() < 180) skew.back() = skew.back() - 180;
						else if (skew.back() < -90 && skew.back() > -180) skew.back() = skew.back() + 180;
						for (int t = 0; t<int(skew.size()) - 1; t++) {
							if ((skew[t] + 3 > skew.back() && skew[t] - 3 < skew.back()) || (abs(skew[t] - skew.back()) > 180 - 3 && abs(skew[t] - skew.back()) < 180 + 3)) {
								skew.pop_back();
								break;
							}
						}
					}
				}
			}
		}
		std::cout << "\nNo. of orientations to check = " << skew.size();

		std::cout << "\n**************************************************************";
		std::cout << "\nSkew of line " << TextLines[i][0] << " = " << TextLines[i][5];
		std::cout << "\nSkew of line " << TextLines[min_line][0] << " = " << TextLines[min_line][5];

		//create ROI...

		newBox[i][0] = 99999;
		newBox[i][1] = 0;
		newBox[i][2] = 99999;
		newBox[i][3] = 0;
		for (int k = 0; k<int(collinear_labels.size()); k++) {
			int lbl = collinear_labels[k];
			if (newBox[i][0] > box[lbl][2])
				newBox[i][0] = box[lbl][2];
			if (newBox[i][1] < box[lbl][3])
				newBox[i][1] = box[lbl][3];
			if (newBox[i][2] > box[lbl][4])
				newBox[i][2] = box[lbl][4]; 
			if (newBox[i][3] < box[lbl][5])
				newBox[i][3] = box[lbl][5];
		}
		int h = newBox[i][3] - newBox[i][2] + 1;
		int w = newBox[i][1] - newBox[i][0] + 1;
		int **roi = new int*[h];
		for (int k = 0; k < h; k++) {
			roi[k] = new int[w];
			for (int l = 0; l < w; l++) {
				roi[k][l] = 255;
				for (int t = 0; t<int(collinear_labels.size()); t++) {
					if (collinear_labels[t] == label[newBox[i][2] + k][newBox[i][0] + l]) {
						roi[k][l] = 0;
						break;
					}
				}
			}
		}

		Mat roi_img;
		roi_img.create(h, w, CV_8UC1);
		arr2mat(roi, roi_img, h, w);
		//imwrite("Data/ICDAR/Data155/roi_rot.tif", roi_img);
		imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/roi_img.tif", roi_img);

		//find p2bratio in possible directions...
		std::vector<Point> lines;
		p2bratios.push_back(findOrientation(roi_img, roi, h, w, skew, lines));
		line_pairs[i].second.second = p2bratios[i].second;

		for (int k = 0; k<int(lines.size()); k++) {
			lines[k] += Point(newBox[i][0], newBox[i][2]);
		}
		lines_data.push_back(lines);
		//std::cout << "\nMedian skew of all lines formed by joining terminal points: " << medianSkew(data) << " degrees";

		if (min_line != 0) {
			//line(text_lines, R, min_point, Scalar(0, 120, 120), 3, 8);
			arrowedLine(text_lines, R, min_point, Scalar(0, 255, 128), 3, 8);
		}


		//imwrite("Data/ICDAR/Data039/text_lines.tif", text_lines);
		imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/text_lines.tif", text_lines);

		//	double theta = linearFit(data); //find the orientation of the linear fitting curve through the end points and the label centroids...
		std::cout << "\n**************************************************************";

		roi_img.release();
		lines.clear();
		collinear_labels.clear();
	}

	//remove the multiple connections by keeping the best ratio connection...
	std::cout << "\nList of line pairs with distance:";
	int *counter = new int[int(line_pairs.size())]; //0: not taken , 1: best ratio...
	for (int i = 0; i<int(line_pairs.size()); i++) {
		counter[i] = -1;
	}
	for (int i = 0; i<int(line_pairs.size()); i++) {
		if (counter[i] != -1) continue;
		counter[i] = 0;
		double max_ratio = p2bratios[i].second; //current ratio
		int max_loc = i; //location
		for (int k = 0; k<int(line_pairs.size()); k++) {
			if (k == i || counter[k] != -1) continue;
			if (line_pairs[i].first.second == line_pairs[k].first.second) {
				//std::cout << "\n...another found...";
				counter[k] = 0;
				if (max_ratio < p2bratios[k].second) {
					max_ratio = p2bratios[k].second;
					max_loc = k;
				}
			}
		}
		for (int k = 0; k<int(line_pairs.size()); k++) {
			if (line_pairs[i].first.second == line_pairs[k].first.second && k != max_loc)
				counter[k] = 1;
		}
	}

	for (int i = 0; i<int(line_pairs.size()); i++) {
		std::cout << "\nline " << line_pairs[i].first.first << " -- line " << line_pairs[i].first.second << " at " << line_pairs[i].second.first << "...ratio = " << line_pairs[i].second.second;
		if (counter[i]) {
			std::cout << "...removing line";
		}
	}


	std::vector<double> ratio_list;
	for (int i = 0; i < int(p2bratios.size()); i++) {
		ratio_list.push_back(p2bratios[i].second);
	}
	double med = median(ratio_list, int(p2bratios.size()));
	std::cout << "\n\n\nMedian of ratios = " << med;

  	int lineNo = 0;
	for (int i = 0; i < int(p2bratios.size()); i++) {
		std::cout << "\n" << i + 1 << ". orientation = " << p2bratios[i].first << " degrees , ratio = " << p2bratios[i].second;
		if (counter[i]) {
			std::cout << "...removed";
			continue;
		}
		if (p2bratios[i].second > 0.5* med) {
			std::cout << "-----okay";

			if (p2bratios[i].first == 90 || p2bratios[i].first == -90) {
				int it = int(TextLines.size());
				for (int k = 0; k<int(lines_data[i].size()); k++) {
					std::vector<int> newLine;
					int y1 = newBox[i][2] - 15;
					int y2 = newBox[i][3] + 15;
					int x = lines_data[i][k].x;
					//	line(text_lines, Point(x, y1), Point(x, y2), Scalar(255, 0, 0), 2, 8);
					line(text_lines, Point(x, y1), Point(x, y2), Scalar(255, 0, 0), 3, 8);
					circle(text_lines, lines_data[i][k], 4, Scalar(255, 0, 255), 2);

					newLine.push_back(++lineNo);
					if (y1 <= y2) {
						newLine.push_back(x);
						newLine.push_back(y1);
						newLine.push_back(x);
						newLine.push_back(y2);
					}
					else {
						newLine.push_back(x);
						newLine.push_back(y2);
						newLine.push_back(x);
						newLine.push_back(y1);
					}
					newLine.push_back(p2bratios[i].first);
					newLine.push_back(3);
					newLine.push_back(p2bratios[i].second);
					TextLinesNew.push_back(newLine); // <line#, L_x, L_y, R_x, R_y, skew, type=3, ratio>...
					ConnectedLines.push_back(make_pair(line_pairs[i].first, make_pair(line_pairs[i].second.second, lineNo)));
				}
			}
			//if nearly vertical then take top to bottom...
			else if ((p2bratios[i].first > 80 && p2bratios[i].first < 110) || (p2bratios[i].first > -110 && p2bratios[i].first < -80)) {
				double slope = tan(p2bratios[i].first*CV_PI / 180);

				for (int k = 0; k<int(lines_data[i].size()); k++) {
					std::vector<int> newLine;
					int x1 = newBox[i][0] - 15;
					int x2 = newBox[i][1] + 15;
					int y1 = lines_data[i][k].y + int((x1 - lines_data[i][k].x)*slope);
					if (y1 < newBox[i][2] - 15) {
						y1 = newBox[i][2] - 15;
						x1 = lines_data[i][k].x + int((y1 - lines_data[i][k].y) / slope);
					}
					if (y1 > newBox[i][3] + 15) {
						y1 = newBox[i][3] + 15;
						x1 = lines_data[i][k].x + int((y1 - lines_data[i][k].y) / slope);
					}
					int y2 = lines_data[i][k].y + int((x2 - lines_data[i][k].x)*slope);
					if (y2 < newBox[i][2] - 15) {
						y2 = newBox[i][2] - 15;
						x2 = lines_data[i][k].x + int((y2 - lines_data[i][k].y) / slope);
					}
					if (y2 > newBox[i][3] + 15) {
						y2 = newBox[i][3] + 15;
						x2 = lines_data[i][k].x + int((y2 - lines_data[i][k].y) / slope);
					}
					//line(text_lines, Point(x1, y1), Point(x2, y2), Scalar(255, 0, 0), 2, 8);
					//line(text_lines, Point(x1, y1), Point(x2, y2), Scalar(0, 128, 255), 2, 8);
					circle(text_lines, lines_data[i][k], 4, Scalar(255, 0, 255), 2);

					newLine.push_back(++lineNo);
					//lines will be top to down...
					if (y1 <= y2) {
						newLine.push_back(x1);
						newLine.push_back(y1);
						newLine.push_back(x2);
						newLine.push_back(y2);
					}
					else {
						newLine.push_back(x2);
						newLine.push_back(y2);
						newLine.push_back(x1);
						newLine.push_back(y1);
					}
					newLine.push_back(p2bratios[i].first);
					newLine.push_back(3);
					newLine.push_back(p2bratios[i].second);
					TextLinesNew.push_back(newLine); // <line#, L_x, L_y, R_x, R_y, skew, type=3, ratio>...
					ConnectedLines.push_back(make_pair(line_pairs[i].first, make_pair(line_pairs[i].second.second, lineNo)));
				}

			}
			else {
				double slope = tan(p2bratios[i].first*CV_PI / 180);

				for (int k = 0; k<int(lines_data[i].size()); k++) {
					std::vector<int> newLine;
					int x1 = newBox[i][0] - 15;
					int x2 = newBox[i][1] + 15;
					int y1 = lines_data[i][k].y + int((x1 - lines_data[i][k].x)*slope);
					if (y1 < newBox[i][2] - 15) {
						y1 = newBox[i][2] - 15;
						x1 = lines_data[i][k].x + int((y1 - lines_data[i][k].y) / slope);
					}
					if (y1 > newBox[i][3] + 15) {
						y1 = newBox[i][3] + 15;
						x1 = lines_data[i][k].x + int((y1 - lines_data[i][k].y) / slope);
					}
					int y2 = lines_data[i][k].y + int((x2 - lines_data[i][k].x)*slope);
					if (y2 < newBox[i][2] - 15) {
						y2 = newBox[i][2] - 15;
						x2 = lines_data[i][k].x + int((y2 - lines_data[i][k].y) / slope);
					}
					if (y2 > newBox[i][3] + 15) {
						y2 = newBox[i][3] + 15;
						x2 = lines_data[i][k].x + int((y2 - lines_data[i][k].y) / slope);
					}
					//line(text_lines, Point(x1, y1), Point(x2, y2), Scalar(255, 0, 0), 2, 8);
					//line(text_lines, Point(x1, y1), Point(x2, y2), Scalar(0, 128, 255), 2, 8);
					circle(text_lines, lines_data[i][k], 4, Scalar(255, 0, 255), 2);

					newLine.push_back(++lineNo);
					//lines will be left to right...
					if (x1 < x2) {
						newLine.push_back(x1);
						newLine.push_back(y1);
						newLine.push_back(x2);
						newLine.push_back(y2);
					}
					else {
						newLine.push_back(x2);
						newLine.push_back(y2);
						newLine.push_back(x1);
						newLine.push_back(y1);
					}
					newLine.push_back(p2bratios[i].first);
					newLine.push_back(3);
					newLine.push_back(p2bratios[i].second);
					TextLinesNew.push_back(newLine); // <line#, L_x, L_y, R_x, R_y, skew, type=3, ratio>...
					ConnectedLines.push_back(make_pair(line_pairs[i].first, make_pair(line_pairs[i].second.second, lineNo)));
				}

			}

		}
		else
			std::cout << "----no";
	}

	//remove overlapping lines with common end points...
	for (int i = 1; i<int(TextLinesNew.size()); i++) {
		int flg = 1;
		Point currL = Point(TextLinesNew[i][1], TextLinesNew[i][2]);
		Point currR = Point(TextLinesNew[i][3], TextLinesNew[i][4]);
		for (int k = 0; k < i; k++) {
			Point prevL = Point(TextLinesNew[k][1], TextLinesNew[k][2]);
			Point prevR = Point(TextLinesNew[k][3], TextLinesNew[k][4]);

			//both ends coincide ...
			if (abs(currL.x - prevL.x) < 10 && abs(currL.y - prevL.y) < 10 && abs(currR.x - prevR.x) < 10 && abs(currR.y - prevR.y) < 10) {
				flg = 0;
			}
			//left end coincides ...
			else if (abs(currL.x - prevL.x) < 10 && abs(currL.y - prevL.y) < 10) {
				if (norm(currL - currR) < norm(prevL - prevR) && abs(norm(currL - currR) + norm(currR - prevR) - norm(prevR - prevL)) < 5) {
					flg = 0;
				}
				else if (norm(currL - currR) >= norm(prevL - prevR) && abs(norm(currL - currR) - norm(currR - prevR) - norm(prevR - prevL)) < 5) {
					flg = 0;
				}
			}

			//right end coincides ...
			else if (abs(currR.x - prevR.x) < 10 && abs(currR.y - prevR.y) < 10) {
				if (norm(currL - currR) < norm(prevL - prevR) && abs(norm(currL - currR) + norm(currL - prevL) - norm(prevR - prevL)) < 5) {
					flg = 0;
				}
				else if (norm(currL - currR) >= norm(prevL - prevR) && abs(norm(currL - currR) - norm(currL - prevL) - norm(prevR - prevL)) < 5) {
					flg = 0;
				}
			}
			if (!flg) {
				std::cout << "\nLine " << i + 1 << " coincides with line " << k + 1;
				break;
			}
		}
		if (!flg) {
			TextLinesNew.erase(TextLinesNew.begin() + i); // <line#, L_x, L_y, R_x, R_y, type=1>...
			ConnectedLines.erase(ConnectedLines.begin() + i);
			i--;
		}
		/*else
		line(text_lines, currL, currR, Scalar(255, 255, 0), 4, 8);	*/

	}
	//keep track of ratio of original lines...
	//include the old lines which have been removed if they have no connections...
	for (int i = 0; i<int(line_pairs.size()); i++) {
		if (counter[i] == 1) { //connection has been removed...
			int flag = 0;
			for (int k = 0; k<int(line_pairs.size()); k++) {
				if (line_pairs[k].first.second == line_pairs[i].first.first && counter[k] == 0) { //connected with other line...
					flag = 1;
					break;
				}
			}
			if (!flag) { //no other connection...
				for (int k = 0; k<int(TextLines.size()); k++) {
					if (TextLines[k][0] == line_pairs[i].first.first) {
						TextLinesNew.push_back(TextLines[k]);
						TextLinesNew[int(TextLinesNew.size()) - 1][0] = ++lineNo;
						ConnectedLines.push_back(make_pair(make_pair(line_pairs[i].first.first, 0), make_pair(lineRatio[k].second, lineNo)));
					}
				}

			}
		}
	}

	/*for (int i = 0; i<int(TextLinesNew.size()); i++) {
		TextLinesNew[i][0] = i + 1;
	}*/

	//Eliminate bad crossing lines...

	//find perp distance of left end and right end....if thr max>50...
	int *del_crossing = new int[int(ConnectedLines.size())];
	for (int i = 0; i<int(TextLinesNew.size()); i++)
		del_crossing[i] = 0;
	for (int i = 0; i<int(TextLinesNew.size()) ; i++) {
		Point L_i = Point(TextLinesNew[i][1], TextLinesNew[i][2]);
		Point R_i = Point(TextLinesNew[i][3], TextLinesNew[i][4]);
		int angle = TextLinesNew[i][5];
		std::cout << "\n********************************************************************************************";

		for (int j = 0; j < int(TextLinesNew.size()); j++) {
			if (i == j) continue;
			Point L_j = Point(TextLinesNew[j][1], TextLinesNew[j][2]);
			Point R_j = Point(TextLinesNew[j][3], TextLinesNew[j][4]);

			Point intersection = linesIntersect(make_pair(TextLinesNew[i][5], L_i), make_pair(TextLinesNew[j][5], L_j)); //<slope, (x_0,y_0)>...
			if (intersection.x > min(L_i.x, R_i.x) - 10 && intersection.x < max(L_i.x, R_i.x) + 10 &&
				intersection.y > min(L_i.y, R_i.y) - 10 && intersection.y < max(L_i.y, R_i.y) + 10 &&
				intersection.x > min(L_j.x, R_j.x) - 10 && intersection.x < max(L_j.x, R_j.x) + 10 &&
				intersection.y > min(L_j.y, R_j.y) - 10 && intersection.y < max(L_j.y, R_j.y) + 10 ) { //(i==35 && j== 37 || i==37 && j==35)

				std::cout << "\n" << TextLinesNew[i][0] <<"(ratio = "<< ConnectedLines[i].second.first << ") - "
					        << TextLinesNew[j][0] << "(ratio = " << ConnectedLines[j].second.first << ") " << " intersects at " << intersection << " , ";

				double L_dist, R_dist;
 				Point L_perp_foot, R_perp_foot;

				if (angle == 0) { //horizontal line...
								  //if (L_j.x > min(L_i.x, R_i.x) && L_j.x < max(L_i.x, R_i.x))
					L_dist = abs(L_i.y - L_j.y);
					R_dist = abs(R_i.y - R_j.y);
				}
				else if (angle == 90 || angle == -90) { //vertical line...
														// (L_j.y > min(L_i.y, R_i.y) && L_j.y < max(L_i.y, R_i.y))
					L_dist = abs(L_i.x - L_j.x);
					R_dist = abs(R_i.x - R_j.x);
				}
				else {
					double slope = tan(angle*CV_PI / 180);
					L_perp_foot.x = int(((slope*(L_j.y - L_i.y) + (slope*slope*L_i.x + L_j.x)) / (slope*slope + 1)));
					L_perp_foot.y = int(L_i.y + slope*(L_perp_foot.x - L_i.x));
					L_dist = norm(L_j - L_perp_foot);

					R_perp_foot.x = int(((slope*(R_j.y - R_i.y) + (slope*slope*R_i.x + R_j.x)) / (slope*slope + 1)));
					R_perp_foot.y = int(R_i.y + slope*(R_perp_foot.x - R_i.x));
					R_dist = norm(R_j - R_perp_foot);
				}
				std::cout << "...Left = " << L_dist << " , Right = " << R_dist;
				if (L_dist > 50 || R_dist > 50) {
					if (ConnectedLines[i].second.first < ConnectedLines[j].second.first) { //compare the ratios...
						std::cout << "\nRemove line " << TextLinesNew[i][0];
						del_crossing[i] = 1;
						for (int k = 0; k<int(line_pairs.size()); k++) {
							if (line_pairs[k].first == ConnectedLines[i].first) {
								counter[k] = 2;
								break;
							}
							else if (ConnectedLines[i].first.second == 0 && 
								(line_pairs[k].first.first == ConnectedLines[i].first.first || line_pairs[k].first.second == ConnectedLines[i].first.first)) {
								counter[k] = 2;
								break;
							}
						}
					}
					else if (ConnectedLines[i].second.first > ConnectedLines[j].second.first) {
						std::cout << "\nRemove line " << TextLinesNew[j][0];
						del_crossing[j] = 1;
						for (int k = 0; k<int(line_pairs.size()); k++) {
							if (line_pairs[k].first == ConnectedLines[j].first) {
								counter[k] = 2;
								break;
							}
							else if (ConnectedLines[j].first.second == 0 &&
								(line_pairs[k].first.first == ConnectedLines[j].first.first || line_pairs[k].first.second == ConnectedLines[j].first.first)) {
								counter[k] = 2;
								break;
							}
						}
					}					
				}
			}
		}

	}
	for (int i = 0; i<int(line_pairs.size()); i++) {
		std::cout << "\nline " << line_pairs[i].first.first << " -- line " << line_pairs[i].first.second << " at " << line_pairs[i].second.first << "...ratio = " << line_pairs[i].second.second;
		if (counter[i] == 1) {
			std::cout << "...removed earlier";
		}
		if (counter[i] == 2) {
			std::cout << "...removing line";
		}
	}

	for (int i = 0; i<int(ConnectedLines.size()); i++) {
		for (int k = 0; k<int(line_pairs.size()); k++) {
			if (line_pairs[k].first == ConnectedLines[i].first && counter[k] == 2) { //located the line no. of the deleted pair...
				for (int j = 0; j<int(TextLinesNew.size()); j++) {
					if (TextLinesNew[j][0] == ConnectedLines[i].second.second) { //delete the textline from the list...
						TextLinesNew.erase(TextLinesNew.begin() + j);
						break;
					}
				}
			}
		}
		if(del_crossing[i])
			for (int j = 0; j<int(TextLinesNew.size()); j++) {
				if (TextLinesNew[j][0] == ConnectedLines[i].second.second) { //delete the textline for bad crossing from the list...
					TextLinesNew.erase(TextLinesNew.begin() + j);
					break;
				}
			}
	}
	
	//remove the connected line corresponding to the erased textline...
	for (int i = 0; i<int(ConnectedLines.size()); i++) {
		int flg = 0;
		for (int k = 0; k<int(TextLinesNew.size()); k++) {
			if (TextLinesNew[k][0] == ConnectedLines[i].second.second) {
				flg = 1;
				break;
			}
		}
		//remove if line no. not found ...
		if (!flg) {
			ConnectedLines.erase(ConnectedLines.begin() + i);
			i--;
		}
	}
	std::cout << "\nCHECK!!!!!";



	//std::cout << "\nNew List of textlines:\n";
	for (int i = 0; i<int(TextLinesNew.size()); i++) {
		TextLinesNew[i][0] = i + 1;
		ConnectedLines[i].second.second = i + 1;
		line(text_lines, Point(TextLinesNew[i][1], TextLinesNew[i][2]), Point(TextLinesNew[i][3], TextLinesNew[i][4]), Scalar(0, 128, 255), 3.5, 8);
		std::cout << "\n" << TextLinesNew[i][0] << ". " << Point(TextLinesNew[i][1], TextLinesNew[i][2]) << " to " << Point(TextLinesNew[i][3], TextLinesNew[i][4])
			<< " skew = " << TextLinesNew[i][5] << " -- Type " << TextLinesNew[i][6] << " connecting old lines "
			<< ConnectedLines[i].first.first << " and " << ConnectedLines[i].first.second << " -- ratio = " << ConnectedLines[i].second.first;
		//imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/text_lines.tif", text_lines);
	}
	

	//imwrite("Data/ICDAR/Data155/text_lines.tif", text_lines);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/text_lines.tif", text_lines);

	delete[] counter;


}

std::vector<std::pair<int,int>> collinearSegments(std::vector<std::vector<int>> &TextLinesNew, 
	                                               std::vector<std::pair<int, Point>> all_centroids, std::vector<std::vector<Point>> all_terminals) {
	std::vector<std::vector<int>> segmentList;
	std::vector<std::pair<int, int>> line_code; //<line no, code no>...
	
	//lines that (almost)overlapp....
	for (int i = 0; i<int(TextLinesNew.size()); i++) {
		std::vector<int> collinear_seg;
		Point L_i = Point(TextLinesNew[i][1], TextLinesNew[i][2]);
		Point R_i = Point(TextLinesNew[i][3], TextLinesNew[i][4]);
		int angle = TextLinesNew[i][5];
		double slope = tan(angle*CV_PI / 180);
		collinear_seg.push_back(TextLinesNew[i][0]);

		
		std::cout << "\nLine " << TextLinesNew[i][0] << " is collinear with lines ";
		for (int k = 0; k<int(TextLinesNew.size()); k++) {
			if (k == i) continue;
			//std::cout << "here";

			Point L_k = Point(TextLinesNew[k][1], TextLinesNew[k][2]);
			Point R_k = Point(TextLinesNew[k][3], TextLinesNew[k][4]);

			//if lines intersect...
			int flg = 0;
			if (TextLinesNew[i][5] != TextLinesNew[k][5]) {
				Point intersection = linesIntersect(make_pair(TextLinesNew[i][5], L_i), make_pair(TextLinesNew[k][5], L_k)); //<slope, (x_0,y_0)>...
				if (intersection.x > min(L_i.x, R_i.x) && intersection.x<max(L_i.x, R_i.x) && 
					intersection.y>min(L_i.y, R_i.y) && intersection.y < max(L_i.y, R_i.y) && 
					intersection.x > min(L_k.x, R_k.x) && intersection.x<max(L_k.x, R_k.x) &&
					intersection.y>min(L_k.y, R_k.y) && intersection.y < max(L_k.y, R_k.y)) {
					std::cout << TextLinesNew[k][0] << "(intersects at " << intersection << ") , ";
					flg = 1;
					collinear_seg.push_back(TextLinesNew[k][0]);
				}
			}

			if(!flg) {
				double dist = -1;

				if (angle == 0) { //horizontal line...
					if (L_k.x > min(L_i.x, R_i.x) && L_k.x < max(L_i.x, R_i.x))
						dist = abs(L_i.y - L_k.y);
				}
				else if (angle == 90 || angle == -90) { //vertical line...
					if (L_k.y > min(L_i.y, R_i.y) && L_k.y < max(L_i.y, R_i.y))
						dist = abs(L_i.x - L_k.x);
				}
				else {
					Point perp_foot;
					perp_foot.x = int(((slope*(L_k.y - L_i.y) + (slope*slope*L_i.x + L_k.x)) / (slope*slope + 1)));
					perp_foot.y = int(L_i.y + slope*(perp_foot.x - L_i.x));
					dist = norm(L_k - perp_foot);
					//std::cout << "\nL = " << L_i << " R = " << R_i << " perp_foot = " << perp_foot << " dist = " << dist;
					if (abs(norm(L_i - perp_foot) + norm(R_i - perp_foot) - norm(L_i - R_i)) > 1) //only if foot is on the segment...
						continue;
					//std::cout << "\nL = " << L << " R = " << R << " perp_foot = " << perp_foot << " dist = " << (1.0 - dist / max_dist);
				}
				if (dist == -1) continue;

				if (dist < 35) { //perpendicular distance...
					std::cout << TextLinesNew[k][0] << " (" << dist << ") ,  ";
					collinear_seg.push_back(TextLinesNew[k][0]);
				}
			}
			
			////if the line is nearly vertical, try to connect in both directions....
			//if ((TextLinesNew[k][5] > -110 && TextLinesNew[k][5] < -80) || (TextLinesNew[k][5] > 80 && TextLinesNew[k][5] < 110)) {
			//	std::cout << "...checking line " << TextLinesNew[k][0] << " in reverse...";
			//	Point R_k = Point(TextLinesNew[k][1], TextLinesNew[k][2]);
			//	Point L_k = Point(TextLinesNew[k][3], TextLinesNew[k][4]);

			//	//if lines intersect...
			//	int flg = 0;
			//	if (TextLinesNew[i][5] != TextLinesNew[k][5]) {
			//		Point intersection = linesIntersect(make_pair(TextLinesNew[i][5], L_i), make_pair(TextLinesNew[k][5], L_k)); //<slope, (x_0,y_0)>...
			//		if (intersection.x > min(L_i.x, R_i.x) && intersection.x<max(L_i.x, R_i.x) &&
			//			intersection.y>min(L_i.y, R_i.y) && intersection.y < max(L_i.y, R_i.y) &&
			//			intersection.x > min(L_k.x, R_k.x) && intersection.x<max(L_k.x, R_k.x) &&
			//			intersection.y>min(L_k.y, R_k.y) && intersection.y < max(L_k.y, R_k.y)) {
			//			std::cout << TextLinesNew[k][0] << "(intersects at " << intersection << ") *in reverse* , ";
			//			flg = 1;
			//			collinear_seg.push_back(TextLinesNew[k][0]);
			//		}
			//	}

			//	if (!flg) {
			//		double dist = -1;

			//		if (angle == 0) { //horizontal line...
			//			if (L_k.x > min(L_i.x, R_i.x) && L_k.x < max(L_i.x, R_i.x))
			//				dist = abs(L_i.y - L_k.y);
			//		}
			//		else if (angle == 90 || angle == -90) { //vertical line...
			//			if (L_k.y > min(L_i.y, R_i.y) && L_k.y < max(L_i.y, R_i.y))
			//				dist = abs(L_i.x - L_k.x);
			//		}
			//		else {
			//			Point perp_foot;
			//			perp_foot.x = int(((slope*(L_k.y - L_i.y) + (slope*slope*L_i.x + L_k.x)) / (slope*slope + 1)));
			//			perp_foot.y = int(L_i.y + slope*(perp_foot.x - L_i.x));
			//			dist = norm(L_k - perp_foot);
			//			//std::cout << "\nL = " << L_i << " R = " << R_i << " perp_foot = " << perp_foot << " dist = " << dist;
			//			if (abs(norm(L_i - perp_foot) + norm(R_i - perp_foot) - norm(L_i - R_i)) > 1) //only if foot is on the segment...
			//				continue;
			//			//std::cout << "\nL = " << L << " R = " << R << " perp_foot = " << perp_foot << " dist = " << (1.0 - dist / max_dist);
			//		}
			//		if (dist == -1) continue;

			//		if (dist < 35) { //perpendicular distance...
			//			std::cout << TextLinesNew[k][0] << " (" << dist << ") ,  ";
			//			collinear_seg.push_back(TextLinesNew[k][0]);
			//		}
			//	}
			//}

		}

 		segmentList.push_back(collinear_seg);
	}


	std::cout << "\nCollinear line segments:";
	for (int i = 0; i<int(segmentList.size()); i++) {
		std::cout << "\n";
		for (int j = 0; j<int(segmentList[i].size()); j++) {
			std::cout << segmentList[i][j] << " , ";
		}
	}


	//mark the collinear line segments with same code...
	int *counter = new int[int(TextLinesNew.size())];
	int newLine = 0;
	for (int i = 0; i<int(TextLinesNew.size()); i++) {
		counter[i] = 0;
	}
	for (int i = 0; i<int(TextLinesNew.size()); i++) {
		if (counter[i] == 0) {
			counter[i] = ++newLine;
		}
		for (int k = 1; k<int(segmentList[i].size()); k++) {
			if (counter[segmentList[i][k] - 1] == counter[i])
				continue;
			else if (counter[segmentList[i][k] - 1] == 0)
				counter[segmentList[i][k] - 1] = counter[i];			
			else if (counter[segmentList[i][k] - 1] != counter[i] && counter[segmentList[i][k] - 1] != 0) {
				int temp = counter[segmentList[i][k] - 1];
				for (int j = 0; j<int(TextLinesNew.size()); j++) {
					if (counter[j] == temp)
						counter[j] = counter[i];
				}
			}
		}
		/*for (int i = 0; i<int(TextLinesNew.size()); i++) {
			std::cout << "\nLine " << TextLinesNew[i][0] << " -- Code " << counter[i];
		}*/
	}

	for (int i = 0; i<int(TextLinesNew.size()); i++) {
		std::cout << "\nLine " << TextLinesNew[i][0] << " -- Code " << counter[i];
		line_code.push_back(make_pair(TextLinesNew[i][0], counter[i]));
	}

	delete[] counter;
	return(line_code);
}

//------------------------------LEVEL 3: Connect "line segments" to get the "text lines"....
//perpendicular distance of point p1 from the line passing through p0 with slope theta...
double perpDist(Point p0, Point p1, int theta) {
	double slope = tan(theta*CV_PI / 180);
	double dist = -1;

	if (theta == 0) { //horizontal line...
		dist = abs(p0.y - p1.y);
	}
	else if (theta == 90 || theta == -90) { //vertical line...
		dist = abs(p0.x - p1.x);
	}
	else {
		Point perp_foot;
		perp_foot.x = int(((slope*(p1.y - p0.y) + (slope*slope*p0.x + p1.x)) / (slope*slope + 1)));
		perp_foot.y = int(p0.y + slope*(perp_foot.x - p0.x));
		dist = norm(p1 - perp_foot);
	}

	return(dist);
}


//void connectSegmentsOld1(int **label, int **box, std::vector<std::pair<int, int>> &label_code, std::vector<std::pair<int, int>> &line_code,
//						std::vector<std::vector<int>> &TextLinesNew, std::vector<std::pair<int, Point>> all_centroids, std::vector<std::vector<Point>> all_terminals,
//						std::vector<std::pair<std::pair<int, int>, std::pair<double, int>>> &ConnectedLines) {
//
//	std::vector<std::vector<int>> newBox; //box dimensions...
//	std::pair<int, double> p2bratios; //angle, ratio...
//	std::vector<std::pair<int, int>> code_pairs; //codes that are same line segment...
//
//	//enlist the unique line codes....
//	std::vector<std::pair<int, int>> code_list; //<code no., median skew of all lines>...
//	//std::vector<int> code_list;
//	int N = int(line_code.size());
//	
//	int *marker = new int[int(line_code.size())];
//	for (int i = 0; i<int(line_code.size()); i++) {
//		marker[i] = 0;
//	}
//
//	for (int i = 0; i<int(line_code.size()); i++) {
//		if (marker[i] != 0) continue;
//		int code = line_code[i].second;
//		std::vector<int> skew;
//		for (int j = 0; j<N; j++) {
//			if (line_code[j].second == code) {
//				skew.push_back(TextLinesNew[line_code[j].first - 1][5]); //store the skew of the textline having code...
//				marker[j] = 1;
//			}
//		}
//		code_list.push_back(make_pair(code, median(skew, int(skew.size()))));
//	}
//	sort(code_list.begin(), code_list.end());
//	std::cout << "\nList of line codes... ";
//	for (int i = 0; i<int(code_list.size()); i++)
//		std::cout << "\n" << code_list[i].first << "\t" << "median skew = " << code_list[i].second;
//
//
//	//connect segments based on line codes pairwise...
//	for (int i = 0; i<int(code_list.size()) - 1; i++) {
//		int code_i = code_list[i].first;
//		for (int j = i + 1; j<int(code_list.size()); j++) {
//			int code_j = code_list[j].first;
//
//			if (code_i == code_j) continue;
//
//			std::cout << "\n*****************************************************************************";
//			std::cout << "\nComparing codes " << code_i << "(" << code_list[i].second << ") and " << code_j << "(" << code_list[j].second << ").......";
//			if (abs(code_list[i].second - code_list[j].second) > 30) {
//				std::cout << "\nDifference = " << abs(code_list[i].second - code_list[j].second) << "....not comparable ";
//				continue;
//			}
//			else
//				std::cout << "\nDifference = " << abs(code_list[i].second - code_list[j].second) << "...comparable.....proceeding";
//			std::vector<int> collinear_lbl;
//			for (int k = 0; k<int(label_code.size()); k++) {
//				if (label_code[k].second == code_i || label_code[k].second == code_j)
//					collinear_lbl.push_back(label_code[k].first);
//			}
//			std::cout << "\nCollinear labels : ";
//			for (int k = 0; k<int(collinear_lbl.size()); k++)
//				std::cout << collinear_lbl[k] << " , ";
//
//
//			//DO THIS LATER....
//			//create ROI...
//
//			//int newBox_i[4] = { 99999,0,99999,0 };
//			//for (int k = 0; k<int(collinear_lbl.size()); k++) {
//			//	int lbl = collinear_lbl[k];
//			//	if (newBox_i[0] > box[lbl][2])
//			//		newBox_i[0] = box[lbl][2];
//			//	if (newBox_i[1] < box[lbl][3])
//			//		newBox_i[1] = box[lbl][3];
//			//	if (newBox_i[2] > box[lbl][4])
//			//		newBox_i[2] = box[lbl][4];
//			//	if (newBox_i[3] < box[lbl][5])
//			//		newBox_i[3] = box[lbl][5];
//			//}
//
//			//int h = newBox_i[3] - newBox_i[2] + 1;
//			//int w = newBox_i[1] - newBox_i[0] + 1;
//			//int **roi = new int*[h];
//			//for (int k = 0; k < h; k++) {
//			//	roi[k] = new int[w];
//			//	for (int l = 0; l < w; l++) {
//			//		roi[k][l] = 255;
//			//		for (int t = 0; t<int(collinear_lbl.size()); t++) {
//			//			if (collinear_lbl[t] == label[newBox_i[2] + k][newBox_i[0] + l]) {
//			//				roi[k][l] = 0;
//			//				break;
//			//			}
//			//		}
//			//	}
//			//}
//
//			//Mat roi_img;
//			//roi_img.create(h, w, CV_8UC1);
//			//arr2mat(roi, roi_img, h, w);
//			////imwrite("Data/ICDAR/Data039/roi_rot.tif", roi_img);
//			//imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/roi_img.tif", roi_img);
//
//
//			//collect terminal points and centroids...
//			std::vector<Point> data;
//			for (int k = 0; k<int(collinear_lbl.size()); k++) {
//				for (int t = 0; t<int(all_centroids.size()); t++) {
//					if (all_centroids[t].first == collinear_lbl[k]) {
//						data.push_back(all_centroids[t].second);
//						break;
//					}
//				}
//				for (int t = 0; t<int(all_terminals.size()); t++) {
//					if (all_terminals[t][0].x == collinear_lbl[k]) {
//						data.insert(data.end(), all_terminals[t].begin() + 1, all_terminals[t].end());
//						break;
//					}
//				}
//			}
//
//			//regression slope with perpendicular offset...
//			double X_mean = 0.0, Y_mean = 0.0;
//			int angle = LSM_perp_off(data, X_mean, Y_mean);
//			Point mean_centroid = Point(X_mean, Y_mean);
//			std::cout << "\nMean of all terminal points = " << mean_centroid;
//
//			//find the min/max distance of each point from the LSM line and the line perpendicular to LSM line...
//			double max_dist = -1, min_dist = 99999, perp_max_dist = -1, perp_min_dist = 99999;
//			std::cout << "\nDistance of " << mean_centroid << " from ";
//			for (int t = 0; t<int(data.size()); t++) {
//				double dist = perpDist(mean_centroid, data[t], angle);
//				double perp_dist = perpDist(mean_centroid, data[t], angle - 90);
//				max_dist = dist > max_dist ? dist : max_dist;
//				min_dist = dist < min_dist ? dist : min_dist;
//				perp_max_dist = perp_dist > perp_max_dist ? perp_dist : perp_max_dist;
//				perp_min_dist = perp_dist < perp_min_dist ? perp_dist : perp_min_dist;
//			}
//			std::cout << "\nLSM direction = " << angle;
//			std::cout << "\nMaximum distance = " << max_dist << "\t Minimum distance = " << min_dist;
//
//			std::cout << "\nperpendicular to LSM direction = " << angle - 90;
//			std::cout << "\nMaximum distance = " << perp_max_dist << "\t Minimum distance = " << perp_min_dist;
//			
//			if (max_dist < perp_max_dist && max_dist < 100) {
//				std::cout << "\n--------------Connecting lines at " << mean_centroid << " at " << angle << " degrees";
//				code_pairs.push_back(make_pair(code_i, code_j));
//			}
//			else if (max_dist > perp_max_dist && perp_max_dist < 100) {
//				std::cout << "\n--------------Connecting lines at " << mean_centroid << " at " << angle - 90 << " degrees";
//				code_pairs.push_back(make_pair(code_i, code_j));
//			}
//			else
//				std::cout << "\n--------Not connected";
//
//			
//			/*for (int k = 0; k<int(ConnectedLines.size()); k++) {
//				if (ConnectedLines[k].second.second == line_code[i].first) {
//					std::cout << "\nLine " << line_code[i].first << "has ratio = " << ConnectedLines[k].second.first;
//				}
//				if (ConnectedLines[k].second.second == line_code[j].first) {
//					std::cout << "\nLine " << line_code[j].first << "has ratio = " << ConnectedLines[k].second.first;
//				}
//			}*/
//			//find p2bratio in LSM directions...
//			//p2bratios = findOrientation(roi_img, roi, h, w, angle);
//			//std::cout << "\nPVD ratio in new direction = " << p2bratios.second;
//
//			/*roi_img.release();
//			for (int k = 0; k < h; k++)
//				delete[] roi[k];
//			delete[] roi;*/
//			collinear_lbl.clear();
//			data.clear();
//		}
//	}
//
//	//display list of connected codes if possible...
//	if (int(code_pairs.size())) {
//
//		//counter of codes...
//		int* counter = new int[int(code_pairs.size())];
//		for (int i = 0; i<int(code_pairs.size()); i++) {
//			counter[i] = 0;
//			std::cout << "\n" << code_pairs[i].first << " --- " << code_pairs[i].second;
//		}
//		std::vector<int> common_code;
//		int code_N = code_list[int(code_list.size()) - 1].first;
//
//		std::cout << "\nConnect the following codes...";
//		for (int i = 0; i<int(code_pairs.size()); i++) {
//			if (counter[i]) continue;
//
//			counter[i] = 1;
//			common_code.push_back(code_pairs[i].first);
//			common_code.push_back(code_pairs[i].second);
//
//			//if not the last, search for more...
//			if (i < int(code_pairs.size()) - 1) {
//				for (int k = i + 1; k<int(code_pairs.size()); k++) {
//					if (code_pairs[i].first == code_pairs[k].first) {
//						counter[k] = 1;
//						int flag = 0;
//						for (int it = 0; it < int(common_code.size()); it++) {
//							if (common_code[it] == code_pairs[k].second) {
//								flag = 1;
//								break;
//							}
//						}
//						if (!flag) 
//							common_code.push_back(code_pairs[k].second);
//					}
//					else if (code_pairs[i].first == code_pairs[k].second) {
//						counter[k] = 1;
//						int flag = 0;
//						for (int it = 0; it < int(common_code.size()); it++) {
//							if (common_code[it] == code_pairs[k].first) {
//								flag = 1;
//								break;
//							}
//						}
//						if (!flag) 
//							common_code.push_back(code_pairs[k].first);
//					}
//				}
//
//				for (int k = i + 1; k<int(code_pairs.size()); k++) {
//					if (code_pairs[i].second == code_pairs[k].first) {
//						counter[k] = 1;
//						int flag = 0;
//						for (int it = 0; it < int(common_code.size()); it++) {
//							if (common_code[it] == code_pairs[k].second) {
//								flag = 1;
//								break;
//							}
//						}
//						if (!flag)
//							common_code.push_back(code_pairs[k].second);
//					}
//					else if (code_pairs[i].second == code_pairs[k].second) {
//						counter[k] = 1;
//						int flag = 0;
//						for (int it = 0; it < int(common_code.size()); it++) {
//							if (common_code[it] == code_pairs[k].first) {
//								flag = 1;
//								break;
//							}
//						}
//						if (!flag)
//							common_code.push_back(code_pairs[k].first);
//					}
//				}
//			}
//			
//
//			int k = 2;
//			while (k<int(common_code.size())) {
//				int check = common_code[k++];
//				for (int t = i + 1; t<int(code_pairs.size()); t++) { //include code# only if not included before...
//					if (check == code_pairs[t].first) {
//						counter[t] = 1;
//						int flag = 0;
//						for (int it = 0; it < int(common_code.size()); it++) {
//							if (common_code[it] == code_pairs[t].second) {
//								flag = 1;
//								break;
//							}
//						}
//						if (!flag) 
//							common_code.push_back(code_pairs[t].second);
//					}
//					else if (check == code_pairs[t].second) {
//						counter[t] = 1;
//						int flag = 0;
//						for (int it = 0; it < int(common_code.size()); it++) {
//							if (common_code[it] == code_pairs[t].first) {
//								flag = 1;
//								break;
//							}
//						}
//						if (!flag) 
//							common_code.push_back(code_pairs[t].first);
//					}
//				}
//			} 
//
//			std::cout << "\nCommon code list: ";
//			for (int k = 0; k<int(common_code.size()); k++)
//				std::cout << "\n" << common_code[k] << " , ";
//			
//			//re-number codes...
//			code_N++;
//			for (int k = 0; k<int(line_code.size()); k++) {
//				int flag = 0;
//				for (int t = 0; t<int(common_code.size()); t++) {
//					if (common_code[t] == line_code[k].second) {
//						flag = 1;
//						break;
//					}
//				}
//				if (flag) {
//					line_code[k].second = code_N; //assign new code no. to the lines of connected pairs...
//				}
//			}
//			for (int k = 0; k<int(label_code.size()); k++) {
//				int flag = 0;
//				for (int t = 0; t<int(common_code.size()); t++) {
//					if (common_code[t] == label_code[k].second) {
//						flag = 1;
//						break;
//					}
//				}
//				if (flag) {
//					label_code[k].second = code_N; //assign new code no. to the labels of connected pairs...
//				}
//			}
//			common_code.clear();
//		}
//
//	}
//	else
//		std::cout << "\nNo codes can be combined!";
//	
//
//}
//
//
//void connectSegmentsOld2(int **label, int **box, std::vector<std::pair<int, int>> &label_code, std::vector<std::pair<int, int>> &line_code,
//	std::vector<std::vector<int>> &TextLinesNew, std::vector<std::pair<int, Point>> all_centroids, std::vector<std::vector<Point>> all_terminals,
//	std::vector<std::pair<std::pair<int, int>, std::pair<double, int>>> &ConnectedLines) {
//
//	std::vector<std::vector<int>> newBox; //box dimensions...
//	std::pair<int, double> p2bratios; //angle, ratio...
//	std::vector<std::pair<int, int>> code_pairs; //codes that are same line segment...
//
//												 //enlist the unique line codes....
//	std::vector<std::pair<int, int>> code_skew_list; //<code no., median skew of all lines>...
//	std::vector<std::pair<int, Point>> code_centroid_list; //<code no., median centroid of all terminal points of lines>...
//	std::vector<std::vector<Point>> code_ends_list;  //<code no., left end, right end>....
//
//	int *marker = new int[int(line_code.size())];
//	for (int i = 0; i<int(line_code.size()); i++) {
//		marker[i] = 0;
//	}
//
//	//median skew of code lines & left/right end points of code...
//	for (int i = 0; i<int(line_code.size()); i++) {
//		if (marker[i] != 0) continue;
//		int code = line_code[i].second;
//		std::vector<int> skew;
//		Point pL, pR;
//		pL.x = 99999, pL.y = 0;
//		pR.x = -1, pR.y = 0;
//		for (int j = 0; j<int(line_code.size()); j++) {
//			if (line_code[j].second == code) {
//				skew.push_back(TextLinesNew[line_code[j].first - 1][5]); //store the skew of the textline having code...
//				if (pL.x > TextLinesNew[line_code[j].first - 1][1])      //store the left end point with minimum x...
//					pL = Point(TextLinesNew[line_code[j].first - 1][1], TextLinesNew[line_code[j].first - 1][2]);
//				if (pR.x < TextLinesNew[line_code[j].first - 1][3])      //store the right end point with maximum x...
//					pR = Point(TextLinesNew[line_code[j].first - 1][3], TextLinesNew[line_code[j].first - 1][4]);
//				marker[j] = 1;
//			}
//		}
//		
//		std::vector<Point> temp;
//		temp.push_back(Point(code, 0));
//		temp.push_back(pL);
//		temp.push_back(pR);
//		code_skew_list.push_back(make_pair(code, median(skew, int(skew.size()))));
//		code_ends_list.push_back(temp);
//	}
//	sort(code_skew_list.begin(), code_skew_list.end());
//	std::cout << "\nList of line codes... ";
//	for (int i = 0; i<int(code_skew_list.size()); i++) {
//		std::cout << "\n" << code_skew_list[i].first << "\t" << "median skew = " << code_skew_list[i].second << " \t " << code_ends_list[i][1] << " \t " << code_ends_list[i][2];
//		int angle = atan2(code_ends_list[i][2].y - code_ends_list[i][1].y, code_ends_list[i][2].x - code_ends_list[i][1].x) * 180 / CV_PI;
//		if (angle > 90 && angle < 180) angle = angle - 180;
//		else if (angle < -90 && angle > -180) angle = angle + 180;
//		std::cout << " \t " << angle;
//	}
//	delete[] marker;
//
//	marker = new int[int(label_code.size())];
//	for (int i = 0; i<int(label_code.size()); i++) {
//		marker[i] = 0;
//	}
//
//	//median centroid of code labels...
//	std::vector<Point> points;
//	for (int i = 0; i<int(label_code.size()); i++) {
//		if (marker[i] != 0 || label_code[i].second == 0) continue;
//		int code = label_code[i].second;
//		for (int j = i; j<int(label_code.size()); j++) {
//			if (label_code[j].second == code) {
//				//collect the centroid...
//				for (int k = 0; k<int(all_centroids.size()); k++) {
//					if (all_centroids[k].first == label_code[j].first) {
//						points.push_back(all_centroids[k].second);
//						break;
//					}
//				}
//				//collect all terminal points...
//				for (int k = 0; k<int(all_terminals.size()); k++) {
//					if (all_terminals[k][0].x == label_code[j].first) {
//						points.insert(points.end(), all_terminals[k].begin() + 1, all_terminals[k].end());
//						break;
//					}
//				}
//				marker[j] = 1;
//			}
//		}
//		
//		code_centroid_list.push_back(make_pair(code, median(points, int(points.size()))));
//		points.clear();
//	}
//	sort(code_centroid_list.begin(), code_centroid_list.end());
//	std::cout << "\nList of line codes... ";
//	for (int i = 0; i<int(code_skew_list.size()); i++)
//		std::cout << "\n" << code_skew_list[i].first << "\t" << "median skew = " << code_skew_list[i].second << "\t" << "median centroid = " << code_centroid_list[i].second;
//	delete[] marker;
//
//	
//}

//void connectSegmentsOld3(int **label, int **box, std::vector<std::pair<int, int>> &label_code, std::vector<std::pair<int, int>> &line_code,
//	std::vector<std::vector<int>> &TextLinesNew, std::vector<std::pair<int, Point>> all_centroids, std::vector<std::vector<Point>> all_terminals,
//	std::vector<std::pair<std::pair<int, int>, std::pair<double, int>>> &ConnectedLines) {
//
//    //enlist the unique line codes....
//	std::vector<std::pair<int, int>> code_skew_list; //<code no., median skew of all lines>...
//	std::vector<std::pair<int, Point>> code_centroid_list; //<code no., median centroid of all terminal points of lines>...
//	std::vector<std::vector<Point>> code_ends_list;  //<code no., left end, right end>....
//
//	int *marker = new int[int(line_code.size())];
//	for (int i = 0; i<int(line_code.size()); i++) {
//		marker[i] = 0;
//	}
//
//	//median skew of code lines & left/right end points of code...
//	for (int i = 0; i<int(line_code.size()); i++) {
//		if (marker[i] != 0) continue;
//		int code = line_code[i].second;
//		std::vector<int> skew;
//		Point pL, pR;
//		pL.x = 99999, pL.y = 0;
//		pR.x = -1, pR.y = 0;
//		for (int j = 0; j<int(line_code.size()); j++) {
//			if (line_code[j].second == code) {
//				skew.push_back(TextLinesNew[line_code[j].first - 1][5]); //store the skew of the textline having code...
//				if (pL.x > TextLinesNew[line_code[j].first - 1][1])      //store the left end point with minimum x...
//					pL = Point(TextLinesNew[line_code[j].first - 1][1], TextLinesNew[line_code[j].first - 1][2]);
//				if (pR.x < TextLinesNew[line_code[j].first - 1][3])      //store the right end point with maximum x...
//					pR = Point(TextLinesNew[line_code[j].first - 1][3], TextLinesNew[line_code[j].first - 1][4]);
//				marker[j] = 1;
//			}
//		}
//
//		std::vector<Point> temp;
//		temp.push_back(Point(code, 0));
//		temp.push_back(pL);
//		temp.push_back(pR);
//		code_skew_list.push_back(make_pair(code, median(skew, int(skew.size()))));
//		code_ends_list.push_back(temp);
//	}
//	sort(code_skew_list.begin(), code_skew_list.end());
//	//sort code_ends_list...
//	std::vector<Point> temp;
//	for (int i = 0; i<int(code_ends_list.size()) - 1; i++) {
//		for (int j = i + 1; j<int(code_ends_list.size()); j++) {
//			if (code_ends_list[i][0].x > code_ends_list[j][0].x) {
//				temp = code_ends_list[i];
//				code_ends_list[i] = code_ends_list[j];
//				code_ends_list[j] = temp;
//			}
//		}
//	}
//
//	std::cout << "\nList of line codes... ";
//	for (int i = 0; i<int(code_skew_list.size()); i++) {
//		std::cout << "\n" << code_skew_list[i].first << "\t" << "median skew = " << code_skew_list[i].second << " \t " << code_ends_list[i][1] << " \t " << code_ends_list[i][2];
//		int angle = atan2(code_ends_list[i][2].y - code_ends_list[i][1].y, code_ends_list[i][2].x - code_ends_list[i][1].x) * 180 / CV_PI;
//		if (angle > 90 && angle < 180) angle = angle - 180;
//		else if (angle < -90 && angle > -180) angle = angle + 180;
//		std::cout << " \t " << angle;
//		code_skew_list[i].second = angle;
//	}
//	delete[] marker;
//
//	std::vector<std::pair<int, int>> code_pairs;
//	for (int i = 0; i<int(code_ends_list.size()); i++) {
//		double min_dist = 9999, dist, min_dist2 = 9999, dist2;
//		int min_code, min_code2;
//		for (int j = 0; j<int(code_ends_list.size()); j++) {
//			if (j == i) continue;
//			if (abs(code_skew_list[i].second - code_skew_list[j].second) > 30) continue; //skew difference should be less than 30deg to compare...
//
//			//find minimum dist = h * sin(theta)...
//
//			dist = norm(code_ends_list[i][2] - code_ends_list[j][1]);
//			if (dist < min_dist) {
//				min_dist = dist;
//				min_code = j;
//			}
//			else if (dist < min_dist2) {
//				min_dist2 = dist;
//				min_code2 = j;
//			}
//		}
//		std::cout << "\n\nLine code " << code_skew_list[i].first << " connects to " << code_skew_list[min_code].first << " distance = " << min_dist
//			<< " OR " << code_skew_list[min_code2].first << " distance = " << min_dist2;
//		int avg_angle1 = int((code_skew_list[i].second + code_skew_list[min_code].second) / 2.0);
//		int avg_angle2 = int((code_skew_list[i].second + code_skew_list[min_code2].second) / 2.0);
//
//		std::vector<Point> seed_points;
//		for (int k = 0; k<int(label_code.size()); k++) {
//			if (label_code[k].second == code_skew_list[i].first || label_code[k].second == code_skew_list[min_code].first) {
//				int lbl = label_code[k].first;
//				for (int t = 0; t<int(all_centroids.size()); t++) {
//					if (all_centroids[t].first == lbl) {
//						seed_points.push_back(all_centroids[t].second);
//						break;
//					}
//				}
//			}
//		}
//		Point meanP = mean(seed_points, int(seed_points.size()));
//		Point medP = median(seed_points, int(seed_points.size()));
//		std::cout << "\n--------draw median line through " << medP << " at angle " << avg_angle1;
//
//		//find perpendicular distance of left end, centre, right end of each line from the median line drawn....
//		double dist_i = -1, dist_min = -1;
//		Point p_i = (code_ends_list[i][1] + code_ends_list[i][2]) ;
//		p_i.x /= 2.0, p_i.y /= 2.0;
//		Point p_min = (code_ends_list[min_code][1] + code_ends_list[min_code][2]);
//		p_min.x /= 2.0, p_min.y /= 2.0;
//
//		double dist_iL = perpDist(code_ends_list[i][1], medP, code_skew_list[i].second);
//		double dist_minL = perpDist(code_ends_list[min_code][1], medP, code_skew_list[i].second);
//
//		double dist_iC = perpDist(p_i, medP, code_skew_list[i].second);
//		double dist_minC = perpDist(p_min, medP, code_skew_list[min_code].second);
//
//		double dist_iR = perpDist(code_ends_list[i][2], medP, code_skew_list[i].second);
//		double dist_minR = perpDist(code_ends_list[min_code][2], medP, code_skew_list[i].second);
//
//		if (dist_iL <= dist_iC && dist_iL <= dist_iR)
//			dist_i = dist_iL;
//
//		else if (dist_iC <= dist_iL && dist_iC <= dist_iR)
//			dist_i = dist_iC;
//
//		else
//			dist_i = dist_iR;
//
//		if (dist_minL <= dist_minC && dist_minL <= dist_minR)
//			dist_min = dist_minL;
//
//		else if (dist_minC <= dist_minL && dist_minC <= dist_minR)
//			dist_min = dist_minC;
//
//		else
//			dist_min = dist_minR;
//
//		std::cout << "\nPerpendicular distance from line code " << code_skew_list[i].first << " = " << dist_i;
//		std::cout << "\nPerpendicular distance from line code " << code_skew_list[min_code].first << " = " << dist_min;
//
//
//		//d=h*sin(theta)...
//		Point u, v1, w1, v2, w2;
//		u = code_ends_list[i][2] - code_ends_list[i][1];
//		v1 = code_ends_list[min_code][2] - code_ends_list[min_code][1];
//		w1 = code_ends_list[min_code][1] - code_ends_list[i][2];
//		v2 = code_ends_list[min_code2][2] - code_ends_list[min_code2][1];
//		w2 = code_ends_list[min_code2][1] - code_ends_list[i][2];
//
//
//		double d1_i = norm(w1)*sin(acos((u.x*w1.x + u.y*w1.y) / (norm(u)*norm(w1))));
//		double d_min1_code = norm(w1)*sin(acos((v1.x*w1.x + v1.y*w1.y) / (norm(v1)*norm(w1))));
//
//		double d2_i = norm(w2)*sin(acos((u.x*w2.x + u.y*w2.y) / (norm(u)*norm(w2))));
//		double d_min2_code = norm(w2)*sin(acos((v2.x*w2.x + v2.y*w2.y) / (norm(v2)*norm(w2))));
//
//
//		std::cout << "\n\nPerpendicular distance from line code " << code_skew_list[i].first << " = " << d1_i;
//		std::cout << "\nPerpendicular distance from line code " << code_skew_list[min_code].first << " = " << d_min1_code;
//
//		std::cout << "\n\nPerpendicular distance from line code " << code_skew_list[i].first << " = " << d2_i;
//		std::cout << "\nPerpendicular distance from line code " << code_skew_list[min_code2].first << " = " << d_min2_code;
//
//
//		if (min_dist < 75 || (d1_i < 35 && d_min1_code < 35))
//		{
//			std::cout << "\nCan connect!";
//			code_pairs.push_back(make_pair(code_skew_list[i].first, code_skew_list[min_code].first));
//		}
//		else if (min_dist2 < 75 || (d2_i < 35 && d_min2_code < 35))
//		{
//			std::cout << "\nCan connect!";
//			code_pairs.push_back(make_pair(code_skew_list[i].first, code_skew_list[min_code2].first));
//		}
//		else
//			std::cout << "\nNot possible...";
//
//		std::cout << "\n--------------------\n";
//		seed_points.clear();
//	}
//
//
//	//display list of connected codes if possible...
//	if (int(code_pairs.size())) {
//	
//		//counter of codes...
//		int* counter = new int[int(code_pairs.size())];
//		for (int i = 0; i<int(code_pairs.size()); i++) {
//			counter[i] = 0;
//			std::cout << "\n" << code_pairs[i].first << " --- " << code_pairs[i].second;
//		}
//		std::vector<int> common_code;
//		int code_N = code_skew_list[int(code_skew_list.size()) - 1].first;
//	
//		std::cout << "\nConnect the following codes...";
//		for (int i = 0; i<int(code_pairs.size()); i++) {
//			if (counter[i]) continue;
//	
//			counter[i] = 1;
//			common_code.push_back(code_pairs[i].first);
//			common_code.push_back(code_pairs[i].second);
//	
//			//if not the last, search for more...
//			if (i < int(code_pairs.size()) - 1) {
//				for (int k = i + 1; k<int(code_pairs.size()); k++) {
//					if (code_pairs[i].first == code_pairs[k].first) {
//						counter[k] = 1;
//						int flag = 0;
//						for (int it = 0; it < int(common_code.size()); it++) {
//							if (common_code[it] == code_pairs[k].second) {
//								flag = 1;
//								break;
//							}
//						}
//						if (!flag) 
//							common_code.push_back(code_pairs[k].second);
//					}
//					else if (code_pairs[i].first == code_pairs[k].second) {
//						counter[k] = 1;
//						int flag = 0;
//						for (int it = 0; it < int(common_code.size()); it++) {
//							if (common_code[it] == code_pairs[k].first) {
//								flag = 1;
//								break;
//							}
//						}
//						if (!flag) 
//							common_code.push_back(code_pairs[k].first);
//					}
//				}
//	
//				for (int k = i + 1; k<int(code_pairs.size()); k++) {
//					if (code_pairs[i].second == code_pairs[k].first) {
//						counter[k] = 1;
//						int flag = 0;
//						for (int it = 0; it < int(common_code.size()); it++) {
//							if (common_code[it] == code_pairs[k].second) {
//								flag = 1;
//								break;
//							}
//						}
//						if (!flag)
//							common_code.push_back(code_pairs[k].second);
//					}
//					else if (code_pairs[i].second == code_pairs[k].second) {
//						counter[k] = 1;
//						int flag = 0;
//						for (int it = 0; it < int(common_code.size()); it++) {
//							if (common_code[it] == code_pairs[k].first) {
//								flag = 1;
//								break;
//							}
//						}
//						if (!flag)
//							common_code.push_back(code_pairs[k].first);
//					}
//				}
//			}
//				
//	
//			int k = 2;
//			while (k<int(common_code.size())) {
//				int check = common_code[k++];
//				for (int t = i + 1; t<int(code_pairs.size()); t++) { //include code# only if not included before...
//					if (check == code_pairs[t].first) {
//						counter[t] = 1;
//						int flag = 0;
//						for (int it = 0; it < int(common_code.size()); it++) {
//							if (common_code[it] == code_pairs[t].second) {
//								flag = 1;
//								break;
//							}
//						}
//						if (!flag) 
//							common_code.push_back(code_pairs[t].second);
//					}
//					else if (check == code_pairs[t].second) {
//						counter[t] = 1;
//						int flag = 0;
//						for (int it = 0; it < int(common_code.size()); it++) {
//							if (common_code[it] == code_pairs[t].first) {
//								flag = 1;
//								break;
//							}
//						}
//						if (!flag) 
//							common_code.push_back(code_pairs[t].first);
//					}
//				}
//			} 
//	
//			std::cout << "\nCommon code list: ";
//			for (int k = 0; k<int(common_code.size()); k++)
//				std::cout << "\n" << common_code[k] << " , ";
//				
//			//re-number codes...
//			code_N++;
//			for (int k = 0; k<int(line_code.size()); k++) {
//				int flag = 0;
//				for (int t = 0; t<int(common_code.size()); t++) {
//					if (common_code[t] == line_code[k].second) {
//						flag = 1;
//						break;
//					}
//				}
//				if (flag) {
//					line_code[k].second = code_N; //assign new code no. to the lines of connected pairs...
//				}
//			}
//			for (int k = 0; k<int(label_code.size()); k++) {
//				int flag = 0;
//				for (int t = 0; t<int(common_code.size()); t++) {
//					if (common_code[t] == label_code[k].second) {
//						flag = 1;
//						break;
//					}
//				}
//				if (flag) {
//					label_code[k].second = code_N; //assign new code no. to the labels of connected pairs...
//
//				}
//			}
//			common_code.clear();
//		}
//	
//	}
//	else
//		std::cout << "\nNo codes can be combined!";
//}


void connectSegments(int **label, int **box, std::vector<std::pair<int, int>> &label_code, std::vector<std::pair<int, int>> &line_code,
	std::vector<std::vector<int>> &TextLinesNew, std::vector<std::pair<int, Point>> all_centroids, std::vector<std::vector<Point>> all_terminals,
	std::vector<std::pair<std::pair<int, int>, std::pair<double, int>>> &ConnectedLines, std::vector<std::pair<int,std::pair<Point, Point>>> &new_lines) {

	//enlist the unique line codes....
	std::vector<std::pair<int, int>> code_skew_list; //<code no., median skew of all lines>...
	std::vector<std::pair<int, Point>> code_centroid_list; //<code no., median centroid of all terminal points of lines>...
	std::vector<std::vector<Point>> code_ends_list;  //<code no., left end, right end>....

	int *marker = new int[int(line_code.size())];
	for (int i = 0; i<int(line_code.size()); i++) {
		marker[i] = 0;
	}

	//median skew of code lines & left/right end points of code...
	for (int i = 0; i<int(line_code.size()); i++) {
		if (marker[i] != 0) continue;
		int code = line_code[i].second;
		std::vector<int> skew;
		Point pL, pR;
		pL.x = 99999, pL.y = 99999;
		pR.x = -1, pR.y = -1;
		for (int j = 0; j<int(line_code.size()); j++) {
			if (line_code[j].second == code) {
				skew.push_back(TextLinesNew[line_code[j].first - 1][5]); //store the skew of the textline having code...
				if ((TextLinesNew[line_code[j].first - 1][5] > -110 && TextLinesNew[line_code[j].first - 1][5] < -80) ||
					(TextLinesNew[line_code[j].first - 1][5] > 80 && TextLinesNew[line_code[j].first - 1][5] < 110)) {
					if (pL.y > TextLinesNew[line_code[j].first - 1][2])      //store the left end point with minimum y...
						pL = Point(TextLinesNew[line_code[j].first - 1][1], TextLinesNew[line_code[j].first - 1][2]);
					if (pR.y < TextLinesNew[line_code[j].first - 1][4])      //store the right end point with maximum y...
						pR = Point(TextLinesNew[line_code[j].first - 1][3], TextLinesNew[line_code[j].first - 1][4]);
				}
				else {
					if (pL.x > TextLinesNew[line_code[j].first - 1][1])      //store the left end point with minimum x...
						pL = Point(TextLinesNew[line_code[j].first - 1][1], TextLinesNew[line_code[j].first - 1][2]);
					if (pR.x < TextLinesNew[line_code[j].first - 1][3])      //store the right end point with maximum x...
						pR = Point(TextLinesNew[line_code[j].first - 1][3], TextLinesNew[line_code[j].first - 1][4]);
				}
				marker[j] = 1;
			}
		}

		std::vector<Point> temp;
		temp.push_back(Point(code, 0));
		temp.push_back(pL);
		temp.push_back(pR);
		code_skew_list.push_back(make_pair(code, median(skew, int(skew.size()))));
		code_ends_list.push_back(temp);
	}
	sort(code_skew_list.begin(), code_skew_list.end());
	//sort code_ends_list...
	std::vector<Point> temp;
	for (int i = 0; i<int(code_ends_list.size()) - 1; i++) {
		for (int j = i + 1; j<int(code_ends_list.size()); j++) {
			if (code_ends_list[i][0].x > code_ends_list[j][0].x) {
				temp = code_ends_list[i];
				code_ends_list[i] = code_ends_list[j];
				code_ends_list[j] = temp;
			}
		}
	}

	for (int i = 0; i<int(code_ends_list.size()); i++) {
		new_lines.push_back(make_pair(code_ends_list[i][0].x,make_pair(code_ends_list[i][1], code_ends_list[i][2])));
	}

	std::cout << "\nList of line codes... ";
	for (int i = 0; i<int(code_skew_list.size()); i++) {
		std::cout << "\n" << code_skew_list[i].first << "\t" << "median skew = " << code_skew_list[i].second << " \t " << code_ends_list[i][1] << " \t " << code_ends_list[i][2];
		int angle = atan2(code_ends_list[i][2].y - code_ends_list[i][1].y, code_ends_list[i][2].x - code_ends_list[i][1].x) * 180 / CV_PI;
		/*if (angle > 90 && angle < 180) angle = angle - 180;
		else if (angle < -90 && angle > -180) angle = angle + 180;*/

		if (angle >= -180 && angle < 0) angle = angle + 180;
		else if (angle == 180) angle = 0;
		std::cout << " \t " << angle;
		code_skew_list[i].second = angle;
	}
	delete[] marker;

	std::vector<std::pair<int, int>> code_pairs;
	
	for (int i = 0; i<int(code_ends_list.size()); i++) {
		double min_dist = 9999, dist, min_dist2 = 9999, dist2;
		int min_code = -1, min_code2 = -1;
		for (int j = 0; j<int(code_ends_list.size()); j++) {
			if (j == i) continue;

			if (code_skew_list[i].second > 150 && code_skew_list[i].second < 180 && code_skew_list[j].second  > 0 && code_skew_list[j].second < 30)  {
				if (abs(180 - code_skew_list[i].second + code_skew_list[j].second) > 30) {
					std::cout << "\n***" << code_skew_list[i].first << " and " << code_skew_list[j].first << " not connected because of angle difference!";
					continue;
				}
			}
			else if (code_skew_list[i].second > 0 && code_skew_list[i].second < 30 && code_skew_list[j].second  > 150 && code_skew_list[j].second < 180) {
				if (abs(180 - code_skew_list[j].second + code_skew_list[i].second) > 30) {
					std::cout << "\n***" << code_skew_list[i].first << " and " << code_skew_list[j].first << " not connected because of angle difference!";
					continue;
				}
			}
			else if (abs(code_skew_list[i].second - code_skew_list[j].second) > 30) {	//skew difference should be less than 30deg to compare...
				std::cout << "\n***" << code_skew_list[i].first << " and " << code_skew_list[j].first << " not connected because of angle difference!";
				continue;
			}

			//find minimum dist = h * sin(theta)...
			Point u, v, w;
			u = code_ends_list[i][2] - code_ends_list[i][1];
			v = code_ends_list[j][2] - code_ends_list[j][1];
			w = code_ends_list[j][1] - code_ends_list[i][2];

			double d_i = norm(w)*sin(acos((u.x*w.x + u.y*w.y) / (norm(u)*norm(w))));
			double d_min_code = norm(w)*sin(acos((v.x*w.x + v.y*w.y) / (norm(v)*norm(w))));
			dist = norm(w);

			if (d_i > 50 || d_min_code > 50 ) continue;
			//if (d_i > 29 || d_min_code > 29) continue; //for ICDAR...
			
			if (dist < min_dist) {
				min_dist = dist;
				min_code = j;
			}
			else if (dist < min_dist2) {
				min_dist2 = dist;
				min_code2 = j;
			}
		}
		if (min_code >= 0) {
			std::cout << "\n\nLine code " << code_skew_list[i].first << " connects to " << code_skew_list[min_code].first << " distance = " << min_dist;
				//<< " OR " << code_skew_list[min_code2].first << " distance = " << min_dist2;
			int avg_angle1 = int((code_skew_list[i].second + code_skew_list[min_code].second) / 2.0);
			//int avg_angle2 = int((code_skew_list[i].second + code_skew_list[min_code2].second) / 2.0);

			if (min_code2 >= 0) {
				std::cout<< "\n\nLine code " << code_skew_list[i].first << " connects to "<< code_skew_list[min_code2].first << " distance = " << min_dist2;
				int avg_angle2 = int((code_skew_list[i].second + code_skew_list[min_code2].second) / 2.0);
			}
			code_pairs.push_back(make_pair(code_skew_list[i].first, code_skew_list[min_code].first));

			//if (min_dist < 75)
			//{
			//	std::cout << "\nCan connect with first!";
			//	code_pairs.push_back(make_pair(code_skew_list[i].first, code_skew_list[min_code].first));
			//}
			//else if (min_code2 >= 0 && min_dist2 < 75)
			//{
			//	std::cout << "\nCan connect with second!";
			//	code_pairs.push_back(make_pair(code_skew_list[i].first, code_skew_list[min_code2].first));
			//}
			//else
			//	std::cout << "\nNot possible...";
		}
		else {
			std::cout << "\nNo connections from line code " << code_skew_list[i].first;
		}

		std::cout << "\n--------------------\n";
	}


	//display list of connected codes if possible...
	if (int(code_pairs.size())) {

		//counter of codes...
		int* counter = new int[int(code_pairs.size())];
		for (int i = 0; i<int(code_pairs.size()); i++) {
			counter[i] = 0;
			std::cout << "\n" << code_pairs[i].first << " --- " << code_pairs[i].second;
		}
		std::vector<int> common_code;
		int code_N = code_skew_list[int(code_skew_list.size()) - 1].first;

		std::cout << "\nConnect the following codes...";
		for (int i = 0; i<int(code_pairs.size()); i++) {
			if (counter[i]) continue;

			counter[i] = 1;
			common_code.push_back(code_pairs[i].first);
			common_code.push_back(code_pairs[i].second);

			//if not the last, search for more...
			if (i < int(code_pairs.size()) - 1) {
				for (int k = i + 1; k<int(code_pairs.size()); k++) {
					if (code_pairs[i].first == code_pairs[k].first) {
						counter[k] = 1;
						int flag = 0;
						for (int it = 0; it < int(common_code.size()); it++) {
							if (common_code[it] == code_pairs[k].second) {
								flag = 1;
								break;
							}
						}
						if (!flag)
							common_code.push_back(code_pairs[k].second);
					}
					else if (code_pairs[i].first == code_pairs[k].second) {
						counter[k] = 1;
						int flag = 0;
						for (int it = 0; it < int(common_code.size()); it++) {
							if (common_code[it] == code_pairs[k].first) {
								flag = 1;
								break;
							}
						}
						if (!flag)
							common_code.push_back(code_pairs[k].first);
					}
				}

				for (int k = i + 1; k<int(code_pairs.size()); k++) {
					if (code_pairs[i].second == code_pairs[k].first) {
						counter[k] = 1;
						int flag = 0;
						for (int it = 0; it < int(common_code.size()); it++) {
							if (common_code[it] == code_pairs[k].second) {
								flag = 1;
								break;
							}
						}
						if (!flag)
							common_code.push_back(code_pairs[k].second);
					}
					else if (code_pairs[i].second == code_pairs[k].second) {
						counter[k] = 1;
						int flag = 0;
						for (int it = 0; it < int(common_code.size()); it++) {
							if (common_code[it] == code_pairs[k].first) {
								flag = 1;
								break;
							}
						}
						if (!flag)
							common_code.push_back(code_pairs[k].first);
					}
				}
			}


			int k = 2;
			while (k<int(common_code.size())) {
				int check = common_code[k++];
				for (int t = i + 1; t<int(code_pairs.size()); t++) { //include code# only if not included before...
					if (check == code_pairs[t].first) {
						counter[t] = 1;
						int flag = 0;
						for (int it = 0; it < int(common_code.size()); it++) {
							if (common_code[it] == code_pairs[t].second) {
								flag = 1;
								break;
							}
						}
						if (!flag)
							common_code.push_back(code_pairs[t].second);
					}
					else if (check == code_pairs[t].second) {
						counter[t] = 1;
						int flag = 0;
						for (int it = 0; it < int(common_code.size()); it++) {
							if (common_code[it] == code_pairs[t].first) {
								flag = 1;
								break;
							}
						}
						if (!flag)
							common_code.push_back(code_pairs[t].first);
					}
				}
			}

			//include the new connections...
			std::vector<std::pair<Point, Point>> commonCodes;
			int newCode = common_code[0];
			std::cout << "\nCommon code list: ";
			for (int k = 0; k<int(common_code.size()); k++) {
				std::cout << "\n" << common_code[k] << " , ";
				for (int t = 0; t<int(code_ends_list.size()); t++) {
					if (code_ends_list[t][0].x == common_code[k])
						commonCodes.push_back(make_pair(code_ends_list[t][1], code_ends_list[t][2]));
				}				
			}
			//sort the vector by ascending order of x-cordinate...
			for (int ti = 0; ti<int(commonCodes.size()) - 1; ti++) {
				for (int tj = ti + 1; tj<int(commonCodes.size()); tj++) {
					std::pair<Point, Point> temp;
					if (commonCodes[ti].first.x > commonCodes[tj].first.x || 
						    (commonCodes[ti].first.x == commonCodes[tj].first.x && commonCodes[ti].first.y > commonCodes[tj].first.y)) {
						temp = commonCodes[ti];
						commonCodes[ti] = commonCodes[tj];
						commonCodes[tj] = temp;
					}
				}
			}
			//join the new connections...
			std::cout << "\nNew connections:";
			for (int t = 0; t<int(commonCodes.size()) - 1; t++) {
				new_lines.push_back(make_pair(newCode,make_pair(commonCodes[t].second, commonCodes[t + 1].first)));
				std::cout << "\njoining " << commonCodes[t].second << " to " << commonCodes[t + 1].first;
			}

			//re-number codes...
			code_N++;
			for (int k = 0; k<int(line_code.size()); k++) {
				int flag = 0;
				for (int t = 0; t<int(common_code.size()); t++) {
					if (common_code[t] == line_code[k].second) {
						flag = 1;
						break;
					}
				}
				if (flag) {
					line_code[k].second = code_N; //assign new code no. to the lines of connected pairs...
				}
			}
			for (int k = 0; k<int(label_code.size()); k++) {
				int flag = 0;
				for (int t = 0; t<int(common_code.size()); t++) {
					if (common_code[t] == label_code[k].second) {
						flag = 1;
						break;
					}
				}
				if (flag) {
					label_code[k].second = code_N; //assign new code no. to the labels of connected pairs...

				}
			}
			for (int k = 0; k<int(new_lines.size()); k++) {
				int flag = 0;
				for (int t = 0; t<int(common_code.size()); t++) {
					if (common_code[t] == new_lines[k].first) {
						flag = 1;
						break;
					}
				}
				if (flag) {
					new_lines[k].first = code_N; //assign new code no. to the new_lines...

				}
			}
			common_code.clear();
		}

	}
	else
		std::cout << "\nNo codes can be combined!";
}


double MahalanobisDist(std::vector<Point> data, Point p) {
	double meanX = 0, meanY = 0, covXX = 0, covYY = 0, covXY = 0;
	int N = int(data.size());
	for (int i = 1; i < N; i++) {
		meanX += data[i].x;
		meanY += data[i].y;
	}
	meanX /= N;
	meanY /= N;

	Point centroid = median(data, N);

	for (int i = 1; i < N; i++) {
		covXX += (data[i].x - centroid.x)*(data[i].x - centroid.x);
		covYY += (data[i].y - centroid.y)*(data[i].y - centroid.y);
		covXY += (data[i].x - centroid.x)*(data[i].y - centroid.y);
	}
	covXX /= (N - 2);
	covYY /= (N - 2);
	covXY /= (N - 2);

	double Sigma[2][2];
	Sigma[0][0] = covXX;
	Sigma[0][1] = covXY;
	Sigma[1][0] = covXY;
	Sigma[1][1] = covYY;

	double det = covXX*covYY - covXY*covXY;
	double invSigma[2][2];
	invSigma[0][0] = covYY / det;
	invSigma[0][1] = -covXY / det;
	invSigma[1][0] = -covXY / det;
	invSigma[1][1] = covXX / det;

	double MD = sqrt(invSigma[0][0] * (p.x - centroid.x)*(p.x - centroid.x) +
		                 invSigma[1][1] * (p.y - centroid.y)*(p.y - centroid.y) + 
		                       2 * invSigma[0][1] * (p.x - centroid.x)*(p.y - centroid.y));
	//std::cout << "\n***MD^2 = " << MD*MD << " ; Sig_X = " << covXX << " , Sig_Y = " << covYY<<" , Sig_XY = " << covXY;
	return(MD);

}


double avg_MD(std::vector<Point> Cluster) {
	int N = int(Cluster.size());
	int n = N - 1;
	double sum = 0.0, avg = 0.0;
	std::vector<double> MD;
	for (int i = 1; i < N; i++) {
		Point p = Cluster[i];
		MD.push_back(MahalanobisDist(Cluster, p));
		sum += MD[i - 1];
	}
	avg = sum / n;
	double chi_sq = 0.0;
	for (int i = 0; i < n; i++) {
		chi_sq = (MD[i] - avg)*(MD[i] - avg);
	}
	chi_sq /= avg;
	return(chi_sq);
}

std::vector<std::pair<Point, Point>> tracing0(Mat img, Mat colorimg, Mat text_lines, int **thin_img, std::vector<std::pair<Point, Point>> &connectCCmax1, std::vector<std::pair<Point, Point>> &connectCCmax2) {
	int **orgimg, **label, **box;
	Mat thin_im, connections;

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

		//mark the centroids...
		for (int k = istart; k < h + istart; k++) {
			for (int l = jstart; l < w + jstart; l++) {
				if (thin_roi[k - istart][l - jstart] == 0)
					thin_img[k][l] = 150;
				for (int t = 0; t < int(all_centroids.size()); t++) {
					if (all_centroids[t].second == Point(l, k))
						thin_img[k][l] = 2;
				}
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
	
	std::vector<std::pair<int, int>> connectCC1, connectCC2;
	std::vector<std::vector<int>> connectLabelsAngles;
	std::vector<std::pair<Point, Point>> connectCentroid = combineBoxes_new(thin_img, label, box, all_terminals, all_centroids, connectCC1, connectCC2, connectLabelsAngles);
	


	for (int k = 0; k<int(connectCC1.size()); k++) {
		Point p1, p2;
		for (int i = 0; i<int(all_centroids.size()); i++) {
			if (connectCC1[k].first == all_centroids[i].first)
				p1 = all_centroids[i].second;
			if (connectCC1[k].second == all_centroids[i].first)
				p2 = all_centroids[i].second;
		}
		connectCCmax1.push_back(pair<Point, Point>(p1, p2));
	}

	for (int k = 0; k<int(connectCC2.size()); k++) {
		Point p1, p2;
		for (int i = 0; i<int(all_centroids.size()); i++) {
			if (connectCC2[k].first == all_centroids[i].first)
				p1 = all_centroids[i].second;
			if (connectCC2[k].second == all_centroids[i].first)
				p2 = all_centroids[i].second;
		}
		connectCCmax2.push_back(pair<Point, Point>(p1, p2));
	}

	//arr2mat(thin_img, colorimg);
	thin_im.create(r, c, CV_8UC1);
	cvtColor(thin_im, thin_im, CV_GRAY2BGR);

	Vec3b colorT, colorC, color1, color0;
	colorT[0] = 0;
	colorT[1] = 0;
	colorT[2] = 255;

	colorC[0] = 255;
	colorC[1] = 0;
	colorC[2] = 255;

	color1[0] = 0;
	color1[1] = 0;
	color1[2] = 0;

	color0[0] = 255;
	color0[1] = 255;
	color0[2] = 255;

	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			if (thin_img[i][j] == 150) {
				//colorimg.at<Vec3b>(i, j) = color1;
				thin_im.at<Vec3b>(i, j) = color1;
			}			
			else if (thin_img[i][j] == 1) {
				colorimg.at<Vec3b>(i, j) = colorT;
				thin_im.at<Vec3b>(i, j) = colorT;
				circle(colorimg, Point(j, i), 2, Scalar(0, 0, 255), 2);
				circle(thin_im, Point(j, i), 2, Scalar(0, 0, 255), 2);
			}
			else if (thin_img[i][j] == 2) {
				colorimg.at<Vec3b>(i, j) = colorC;
				thin_im.at<Vec3b>(i, j) = colorC;
				circle(colorimg, Point(j, i), 2, Scalar(255, 0, 255), 2);
				circle(thin_im, Point(j, i), 2, Scalar(255, 0, 255), 2);
			}
			else {
				//colorimg.at<Vec3b>(i, j) = color0;
				thin_im.at<Vec3b>(i, j) = color0;
			}
		}
	}
	imshow("Thinbw.jpg", thin_im);
	/*imwrite("Data/ICDAR/Data155/ThinTeminal.tif", colorimg);
	imwrite("Data/ICDAR/Data155/ThinBW.tif", thin_im);*/

	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/ThinTeminal.tif", colorimg);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/ThinBW.tif", thin_im);

	//************** STEP 1 : COMBINE LABELS*********************

	std::vector<std::vector<int>> TextLines; // <line no., left_x, left_y, right_x, right_y, skew, type, ratio>...
	std::vector<std::pair<int, double>> lineRatio; //<line no., ratio>...
		
	//Extract lines of type 1 connecting multiple CCs...	
	std::vector<pair<pair<int, int>, int>> label_deg_1 = combineLabels(colorimg, text_lines, connectCC1, connectCC2, all_centroids, label, box, connectLabelsAngles, TextLines, lineRatio);
	std::vector<int> unmarked_labels;
	for (int i = 1; i < labelNum; i++) {
		if (box[i][2] == -1) continue;
		int flag = 1;
		for (int k = 0; k<int(label_deg_1.size()); k++) {
			if (label_deg_1[k].first.first == i) {
				flag = 0;
				break;
			}
		}
		if (flag)
			unmarked_labels.push_back(i);
	}

	//List of labels and textlines:
	std::vector<pair<pair<int, int>, pair<int, double>>> label2line; //<label, line no., skew, distance>
	for (int i = 0; i<int(label_deg_1.size()); i++)  {
		Point centroid;
		for (int k = 0; k<int(all_centroids.size()); k++) {
			if (all_centroids[k].first == label_deg_1[i].first.first) {
				centroid = all_centroids[k].second;
				break;
			}
		}
		std::cout << "\n----------------\nLooking for nearest line to label " << label_deg_1[i].first.first << "...";
		double dist, min_dist = 9999;
		int min_dist_lineNo = 0, min_angle = 0;
		for (int k = 0; k<int(TextLines.size()); k++) {
			Point L = Point(TextLines[k][1], TextLines[k][2]);
			Point R = Point(TextLines[k][3], TextLines[k][4]);
			double skew = TextLines[k][5];

			dist = -1; 

			if (skew == 0) { //horizontal line...
				if (centroid.x > min(L.x, R.x) && centroid.x < max(L.x, R.x))
					dist = abs(centroid.y - R.y);
			}
			else if (skew == 90 || skew == -90) { //vertical line...
				if (centroid.y > min(L.y, R.y) && centroid.y < max(L.y, R.y))
					dist = abs(centroid.x - R.x);
			}
			else {
				Point perp_foot;
				double slope = tan(skew*CV_PI / 180);
				perp_foot.x = int(((slope*(centroid.y - L.y) + (slope*slope*L.x + centroid.x)) / (slope*slope + 1)));
				perp_foot.y = int(L.y + slope*(perp_foot.x - L.x));
				//std::cout << "\nL = " << L << " R = " << R << " perp_foot = " << perp_foot << " dist = " << norm(centroid - perp_foot);
				if (abs(norm(L - perp_foot) + norm(R - perp_foot) - norm(L - R)) < 1) {
					std::cout << "  ok...";
					dist = norm(centroid - perp_foot);
				}
			}

			if (dist == -1) continue;

			if (dist < min_dist) {
				min_dist = dist;
				min_dist_lineNo = TextLines[k][0];
				min_angle = TextLines[k][5];
			}
			
		}
		std::cout << "\nNearest line no. : " << min_dist_lineNo << " at distance " << min_dist << " at angle " << min_angle;
	//	label2line.push_back(make_pair(make_pair(label_deg_1[i].first.first, min_dist_lineNo), make_pair(min_angle, min_dist)));

		if (min_dist > 35) {
			label2line.push_back(make_pair(make_pair(label_deg_1[i].first.first, 0), make_pair(0, -1)));
			std::cout << "...not included";
		}
		else
			label2line.push_back(make_pair(make_pair(label_deg_1[i].first.first, min_dist_lineNo), make_pair(min_angle, min_dist)));
 	}

	
	int text_n = int(TextLines.size());
	for (int i = 0; i<int(TextLines.size()); i++) {
		int flg = 0;
		for (int k = 0; k<int(label2line.size()); k++) {
			if (label2line[k].first.second == TextLines[i][0]) {
				flg = 1;
				break;
			}
		}
		if (!flg) {
			std::cout << "\nLine " << TextLines[i][0] << " has no label"; 
			//line(text_lines, Point(TextLines[i][1], TextLines[i][2]), Point(TextLines[i][3], TextLines[i][4]), Scalar(255, 0, 255), 3, 8);
			TextLines.erase(TextLines.begin() + i);
			lineRatio.erase(lineRatio.begin() + i);
			i--;
		}
	}

	if (text_n > int(TextLines.size())) {  //some lines  have been erased...
		for (int i = 0; i<int(label2line.size()); i++) {
			std::cout << "\nLabel " << label2line[i].first.first  << " -- Line No. " << label2line[i].first.second
				<< " -- " << label2line[i].second.first << " degrees" << " -- " << label2line[i].second.second << " units";
			if (label2line[i].first.second == 0)
				std::cout << " ***";
		}

		std::vector<pair<pair<int, int>, pair<int, double>>> copy_label2line = label2line;
		for (int i = 0; i<int(TextLines.size()); i++) {
			for (int k = 0; k<int(label2line.size()); k++) {
				if (copy_label2line[k].first.second == TextLines[i][0])
					label2line[k].first.second = i + 1;
			}
			TextLines[i][0] = i + 1;
			lineRatio[i].first = i + 1;
		}
		copy_label2line.clear();
	}
	

 	for (int i = 0; i<int(label2line.size()); i++) {
		std::cout << "\nLabel " << label2line[i].first.first <<  " -- Line No. " << label2line[i].first.second
			<< " -- " << label2line[i].second.first << " degrees" << " -- " << label2line[i].second.second << " units";
		if (label2line[i].first.second == 0)
			std::cout << " ***";
	}

	//************** STEP 1 : MARK SINGLE LABELS*********************

	//Extract lines of type 2 having single CC...
	std::vector<pair<pair<int, int>, int>> label_deg_2 = singleLabels(colorimg,text_lines, all_centroids, all_terminals, label, box, unmarked_labels, TextLines, lineRatio, label2line);

	for (int i = 0; i<int(TextLines.size()); i++) {
		if(TextLines[i][6] == 1)
			line(text_lines, Point(TextLines[i][1], TextLines[i][2]), Point(TextLines[i][3], TextLines[i][4]), Scalar(255, 0, 0), 3, 8);
		if (TextLines[i][6] == 2)
			line(text_lines, Point(TextLines[i][1], TextLines[i][2]), Point(TextLines[i][3], TextLines[i][4]), Scalar(120, 200, 0), 3, 8);
	}
	for (int i = 0; i<int(all_centroids.size()); i++) {
		circle(text_lines, all_centroids[i].second, 4, Scalar(0, 0, 255), 2);
	}
	//imwrite("Data/ICDAR/Data155/text_lines.tif", text_lines);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/text_lines.tif", text_lines);
	

	//include the labels that are not of type 1 or 2...
	for (int i = 1; i < labelNum; i++) {
		if (box[i][2] == -1) continue;

		int flag = 0;
		for (int k = 0; k<int(label2line.size()); k++) {
			if (label2line[k].first.first == i) { //nearest line to label already checked...
				flag = 1;
				break;
			}
		}
		//NO
		if (!flag) { //check if not checked already...
			Point centroid;
			for (int k = 0; k<int(all_centroids.size()); k++) {
				if (all_centroids[k].first == i) {
					centroid = all_centroids[k].second;
					break;
				}
			}
			std::cout << "\n----------------\nLooking for nearest line to label " << i << "...";
			double dist, min_dist = 9999;
			int min_dist_lineNo = 0, min_angle = 0;
			for (int k = 0; k<int(TextLines.size()); k++) {
				Point L = Point(TextLines[k][1], TextLines[k][2]);
				Point R = Point(TextLines[k][3], TextLines[k][4]);
				double skew = TextLines[k][5];

				dist = -1;


				if (skew == 0) { //horizontal line...
					if (centroid.x > min(L.x, R.x) && centroid.x < max(L.x, R.x))
						dist = abs(centroid.y - R.y);
				}
				else if (skew == 90 || skew == -90) { //vertical line...
					if (centroid.y > min(L.y, R.y) && centroid.y < max(L.y, R.y))
						dist = abs(centroid.x - R.x);
				}
				else {
					Point perp_foot;
					double slope = tan(skew*CV_PI / 180);
					perp_foot.x = int(((slope*(centroid.y - L.y) + (slope*slope*L.x + centroid.x)) / (slope*slope + 1)));
					perp_foot.y = int(L.y + slope*(perp_foot.x - L.x));
					//std::cout << "\nL = " << L << " R = " << R << " perp_foot = " << perp_foot << " dist = " << norm(centroid - perp_foot);
					if (abs(norm(L - perp_foot) + norm(R - perp_foot) - norm(L - R)) < 1) {
						std::cout << "  ok...";
						dist = norm(centroid - perp_foot);
					}
				}

				if (dist == -1) continue;

				if (dist < min_dist) {
					min_dist = dist;
					min_dist_lineNo = TextLines[k][0];
					min_angle = TextLines[k][5];
				}

			}
			std::cout << "\nNearest line no. : " << min_dist_lineNo << " at distance " << min_dist << " at angle " << min_angle;
			//label2line.push_back(make_pair(make_pair(i, min_dist_lineNo), make_pair(min_angle, min_dist)));

			if (min_dist > 35) {
				label2line.push_back(make_pair(make_pair(i, 0), make_pair(0, -1)));
				std::cout << "...not included";
			}
			else
				label2line.push_back(make_pair(make_pair(i, min_dist_lineNo), make_pair(min_angle, min_dist)));
		}
	}


	sort(label2line.begin(), label2line.end());
	std::cout << "\nLabel to lines map:(many-1)\n";
	int k = 0;
	for (int i = 0; i<int(label2line.size()); i++) {
		std::cout << "\nLabel " << label2line[i].first.first << all_centroids[k++].second << " -- Line No. " << label2line[i].first.second
			<< " -- " << label2line[i].second.first << " degrees" << " -- " << label2line[i].second.second << " units";
		if (label2line[i].first.second == 0)
			std::cout << " ***";
	}

	std::cout << "\n******************************\n";
	std::cout << "\nLines to label map:(1-many)\n";
	for (int i = 0; i<int(TextLines.size()); i++) {
		int flg = 0;
		for (int k = 0; k<int(label_deg_2.size()); k++) {
			if (label_deg_2[k].second == TextLines[i][0]) {
				flg = 1;
				TextLines[i].push_back(label_deg_2[k].first.first);
			}
		}
		if (!flg)
			TextLines[i].push_back(0);
		std::cout << "\n" << TextLines[i][0] << ". " << Point(TextLines[i][1], TextLines[i][2]) << " to "
			<< Point(TextLines[i][3], TextLines[i][4]) << " skew = " << TextLines[i][5] << " -- Type " << TextLines[i][6] << " -- label " << TextLines[i][7];
	}

	

	//imwrite("Data/ICDAR/Data155/text_lines.tif", text_lines);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/text_lines.tif", text_lines);

	//************** STEP 2 : PAIR LINES*********************

	std::vector<std::vector<int>> TextLinesNew; //<line#, L_x, L_y, R_x, R_y, skew, type=3, ratio>...
	std::vector<pair<pair<int, int>, pair<int, double>>> label2lineNew; //<label, line no., skew, distance>
	std::vector<pair<int, double>> lineRatio2; //<orange line no. , ratio>...
	std::vector<std::pair<std::pair<int, int>, std::pair<double, int>>> ConnectedLines; //< <line pair 1, line pair 2>, <ratio, newline> >....

	connectLines(text_lines, all_centroids, all_terminals, label, box, TextLines, lineRatio, TextLinesNew, label2line, ConnectedLines);
	
	//List of labels and textlines:
	for (int i = 0; i<int(all_centroids.size()); i++) {
		Point centroid = all_centroids[i].second;
		std::cout << "\n----------------\nLooking for nearest line to label " << all_centroids[i].first << "...";
		double dist, min_dist = 9999;
		int min_dist_lineNo = 0, min_angle = 0;
		for (int k = 0; k<int(TextLinesNew.size()); k++) {
			Point L = Point(TextLinesNew[k][1], TextLinesNew[k][2]);
			Point R = Point(TextLinesNew[k][3], TextLinesNew[k][4]);
			double skew = TextLinesNew[k][5];

			dist = -1;

			if (skew == 0) { //horizontal line...
				if (centroid.x > min(L.x, R.x) && centroid.x < max(L.x, R.x))
					dist = abs(centroid.y - R.y);
			}
			else if (skew == 90 || skew == -90) { //vertical line...
				if (centroid.y > min(L.y, R.y) && centroid.y < max(L.y, R.y))
					dist = abs(centroid.x - R.x);
			}
			else {
				Point perp_foot;
				double slope = tan(skew*CV_PI / 180);
				perp_foot.x = int(((slope*(centroid.y - L.y) + (slope*slope*L.x + centroid.x)) / (slope*slope + 1)));
				perp_foot.y = int(L.y + slope*(perp_foot.x - L.x));
				//std::cout << "\nL = " << L << " R = " << R << " perp_foot = " << perp_foot << " dist = " << norm(centroid - perp_foot);
				if (abs(norm(L - perp_foot) + norm(R - perp_foot) - norm(L - R)) < 1) { //only if point lies in range of line...
					//std::cout << "  ok...";
					dist = norm(centroid - perp_foot);
				}
			}

			if (dist == -1) continue;

			if (dist < min_dist) {
				min_dist = dist;
				min_dist_lineNo = TextLinesNew[k][0];
				min_angle = TextLinesNew[k][5];
			}

		}
		std::cout << "\nNearest line no. : " << min_dist_lineNo << " at distance " << min_dist << " at angle " << min_angle;
		if (min_dist > 35) {
			label2lineNew.push_back(make_pair(make_pair(all_centroids[i].first, 0), make_pair(0, -1)));
			std::cout << "...not included";
		}
		else
			label2lineNew.push_back(make_pair(make_pair(all_centroids[i].first, min_dist_lineNo), make_pair(min_angle, min_dist)));		
	}

	//check for lines with no labels...
	for (int i = 0; i<int(TextLinesNew.size()); i++) {
		int flg = 0;
		for (int k = 0; k<int(label2line.size()); k++) {
			if (label2lineNew[k].first.second == TextLinesNew[i][0]) {
				flg = 1;
				break;
			}
		}
		if (!flg) {
			std::cout << "\nLine " << TextLinesNew[i][0] << " has no label";
			line(text_lines, Point(TextLines[i][1], TextLines[i][2]), Point(TextLines[i][3], TextLines[i][4]), Scalar(0, 0, 255), 3, 8);
			//TextLinesNew[i][0] = -1;
			TextLinesNew.erase(TextLinesNew.begin() + i);
			i--;
		}
	}
	//check for lines with single labels and 2 terminal points...
	for (int i = 0; i<int(TextLinesNew.size()); i++) {
		int flg = 0;
		int count = 0;
		int lbl;
		int pos = 0;
		for (int k = 0; k<int(label2lineNew.size()); k++) {
			if (label2lineNew[k].first.second == TextLinesNew[i][0]) {
				count++;
				lbl = label2lineNew[k].first.first;
				pos = k;
			}
		}
		std::cout << "\nLine " << TextLinesNew[i][0] << " has " << count << " labels";
		if (count > 1) continue;
		if (count == 1) {
			std::cout << "\nLine " << TextLinesNew[i][0] << " has " << count << " label...";
			for (int k = 0; k<int(all_terminals.size()); k++) {
				if (all_terminals[k][0].x == lbl && int(all_terminals[k].size()) <= 3) {
					flg = 1;
					break;
				}
			}
		}
		if (flg) {
			//remove line...
			std::cout << "\nLine " << TextLinesNew[i][0] << " has single label(" << lbl << ") with 2 terminal points";
			//line(text_lines, Point(TextLines[i][1], TextLines[i][2]), Point(TextLines[i][3], TextLines[i][4]), Scalar(0, 0, 255), 3, 8);
			//TextLinesNew[i][0] = -1;
			TextLinesNew.erase(TextLinesNew.begin() + i);
			i--;
			//remove label2line connection...
			label2lineNew[pos].first.second = 0;
			label2lineNew[pos].second.first = 0;
			label2lineNew[pos].second.second = -1;
		}
	}

	//re-number TextLinesNew...
 	std::vector<pair<int, int>> map_line_no;
	for (int i = 0; i<int(TextLinesNew.size()); i++) {
		map_line_no.push_back(make_pair(TextLinesNew[i][0], i + 1));
		TextLinesNew[i][0] = i + 1;
	}
	for (int i = 0; i<int(label2lineNew.size()); i++) {
		for (int k = 0; k<int(map_line_no.size()); k++) {
			if (label2lineNew[i].first.second == map_line_no[k].first) {
				label2lineNew[i].first.second = map_line_no[k].second;
				break;
			}
		}
	}
	
	//************** STEP 3 : CONNECT COLLINEAR SEGMENTS*********************
	

	//label-collinear code matching...
	std::vector<std::pair<int, int>> line_code = collinearSegments(TextLinesNew, all_centroids, all_terminals); //<line no, code>... 

	for (int i = 0; i<int(label2lineNew.size()); i++) {
		if (label2lineNew[i].first.second == 0) {
			std::cout << "\nLabel " << label2lineNew[i].first.first << all_centroids[i].second << " -- Line No. " << label2lineNew[i].first.second << "(" << 0 << ")"
				<< " -- " << label2lineNew[i].second.first << " degrees" << " -- " << label2lineNew[i].second.second << " units ***"; 
		}
		else
			std::cout << "\nLabel " << label2lineNew[i].first.first << all_centroids[i].second << " -- Line No. " << label2lineNew[i].first.second << "(" << line_code[label2lineNew[i].first.second - 1].second << ")"
			<< " -- " << label2lineNew[i].second.first << " degrees" << " -- " << label2lineNew[i].second.second << " units";
	}

	std::vector<std::pair<int, int>> label_code; //<label, code>...
	for (int i = 0; i<int(label2lineNew.size()); i++) {
		if (label2lineNew[i].first.second == 0) {
			label_code.push_back(make_pair(label2lineNew[i].first.first, 0));
		}
		else {			
			label_code.push_back(make_pair(label2lineNew[i].first.first, line_code[label2lineNew[i].first.second - 1].second));
		}
		std::cout << "\nLabel " << label_code.back().first << " -- Code " << label_code.back().second;
	}

	std::vector<std::pair<int,std::pair<Point, Point>>> new_lines;
	connectSegments(label, box, label_code, line_code, TextLinesNew, all_centroids, all_terminals, ConnectedLines, new_lines);

	std::cout << "\nNew Label-Code list:";
	for (int i = 0; i<int(label_code.size()); i++) {
		std::cout << "\nLabel " << label_code[i].first << " -- Code " << label_code[i].second;
	}

	std::cout << "\nNew Line-Code list:";
	for (int i = 0; i<int(new_lines.size()); i++) {
		std::cout << "\nLine Code " << new_lines[i].first << " :: " << new_lines[i].second.first << " -- " << new_lines[i].second.second;
	}

	//************** STEP 4 : ASSIGN REMAINING LABELS *********************


	for (int i = 0; i<int(label_code.size()); i++) {
		if (label_code[i].second == 0) {
			int lbl = label_code[i].first;
			Point p;
			for (int k = 0; k<int(all_centroids.size()); k++) {
				if (all_centroids[k].first == lbl) {
					p = all_centroids[k].second;
				}
			}
			double min_dist1 = 99999, min_dist2 = 99999;
			int min_dist_code1 = 0, min_dist_code2 = 0;
			std::cout << "\nLabel " << lbl << p << "---";
			for (int t = 0; t<int(new_lines.size()); t++) {
				Point L = new_lines[t].second.first;
				Point R = new_lines[t].second.second;

				int min_x = L.x < R.x ? L.x : R.x;
				int min_y = L.y < R.y ? L.y : R.y;

				int max_x = L.x >= R.x ? L.x : R.x;
				int max_y = L.y >= R.y ? L.y : R.y;

				int skew = int(atan2(L.y - R.y, L.x - R.x) * 180 / CV_PI);
				if (skew > 90 && skew < 180) skew = skew - 180;
				else if (skew < -90 && skew > -180) skew = skew + 180;
				double distL = norm(L - p);
				double distR = norm(R - p);
				double nearest_end = distL < distR ? distL : distR;

				Point vec1 = R - L;
				Point vec2 = p - R;
				
				
				double dist;
				Point perp_foot;

				if (skew == 0) { //horizontal line...
					dist = abs(p.y - L.y);
					perp_foot = Point(p.x, L.y);
				}
				else if (skew == 90 || skew == -90) { //vertical line...
					dist = abs(p.x - L.x);
					perp_foot = Point(p.y, L.x);
				}
				else {
					double slope = tan(skew*CV_PI / 180);
					perp_foot.x = int(((slope*(p.y - L.y) + (slope*slope*L.x + p.x)) / (slope*slope + 1)));
					perp_foot.y = int(L.y + slope*(perp_foot.x - L.x));
					dist = norm(p - perp_foot);
				}

				//std::cout << "\ncode " << new_lines[t].first << " distance = " << dist << " : base point = " << perp_foot;
				if ((perp_foot.x > min_x - 10 && perp_foot.x < max_x + 10 && perp_foot.y > min_y - 10 && perp_foot.y < max_y + 10 ) ) {
					if (dist < min_dist1) {
						min_dist1 = dist;
						min_dist_code1 = new_lines[t].first;
					}
				}
				//else {
				//	double angle_change = acos((vec1.x*vec2.x + vec1.y*vec2.y) / (norm(vec1)*norm(vec2))) * 180 / CV_PI;
				//	if (angle_change > 90 && angle_change <= 180)
				//		angle_change -= 180;
				//	std::cout << "\nAngle change = " << angle_change << " with line code " << new_lines[t].first;
				//	
				//	std::cout << "..." << sin(angle_change*CV_PI / 180)*norm(vec2);
				//	//if (abs(angle_change) <= 30  && abs(angle_change) < abs(min_dist2))
				//	if (abs(cos(angle_change*CV_PI / 180)*norm(vec2)) < abs(min_dist2)) {
				//		min_dist2 = cos(angle_change*CV_PI / 180)*norm(vec2);
				//		min_dist_code2 = new_lines[t].first;
				//	}
				//}
			}
			label_code[i].second = min_dist1 <= 60 ? min_dist_code1 : min_dist_code2;
			std::cout << "\nLabel " << label_code[i].first << " is near code " << min_dist_code1 << " at perpendicular distance = " << min_dist1;
			std::cout << "\nLabel " << label_code[i].first << " is near code " << min_dist_code2 << " at angle = " << min_dist2;
		}
	}

	//create clusters with centroids and terminal points...
	int *markCodes = new int[int(label_code.size())];
	std::vector<std::vector<Point>> lineClusters; // code-wise list of all centroids & terminal points...
	std::vector<std::vector<int>> CodeList; // code-wise list of all labels....
	for (int i = 0; i<int(label_code.size()); i++)
		markCodes[i] = 0;
	for (int i = 0; i<int(label_code.size()); i++) {
		if (label_code[i].second == 0) continue;
		if (markCodes[i]) continue;
		std::vector<Point> cluster;
		std::vector<int> code;
		cluster.push_back(Point(label_code[i].second, 0));
		code.push_back(label_code[i].second);
		for (int j = i; j<int(label_code.size()); j++) {
			if (label_code[j].second != label_code[i].second) continue;
			markCodes[j] = 1;
			code.push_back(label_code[j].first); //collect the labels in the cluster...

			//collect all the centroids and terminal points with same code...
			for (int t = 0; t<int(all_centroids.size()); t++) {
				if (all_centroids[t].first == label_code[j].first) {
					cluster.push_back(all_centroids[t].second);
					break;
				}
			}
			for (int t = 0; t<int(all_terminals.size()); t++) {
				if (all_terminals[t][0].x == label_code[j].first) {
					for (int tt = 1; tt<int(all_terminals[t].size()); tt++) {
						cluster.push_back(all_terminals[t][tt]);
					}
					break;
				}
			}
		}
		lineClusters.push_back(cluster);
		CodeList.push_back(code);

		cluster.clear();
		code.clear();
	}

	std::cout << "\nLabel list code-wise...";
	for (int i = 0; i<int(CodeList.size()); i++) {
		std::cout << "\nCode " << CodeList[i][0] << " :: ";
		for (int j = 1; j<int(CodeList[i].size()); j++) {
			std::cout << CodeList[i][j] << " , ";
		}
	}


	//Display each cluster box...
	//for (int i = 0; i<int(CodeList.size()); i++) {
	//	//find the bounding box for the cluster...
	//	int newBox[4] = { 99999,0,99999,0 };
	//	for (int j = 1; j<int(CodeList[i].size()); j++) {
	//		int lbl = CodeList[i][j];
	//		newBox[0] = newBox[0] > box[lbl][2] ? box[lbl][2] : newBox[0];
	//		newBox[1] = newBox[1] < box[lbl][3] ? box[lbl][3] : newBox[1];
	//		newBox[2] = newBox[2] > box[lbl][4] ? box[lbl][4] : newBox[2];
	//		newBox[3] = newBox[3] < box[lbl][5] ? box[lbl][5] : newBox[3];
	//	}

	//	int **roi;
	//	int istart = newBox[2];
	//	int jstart = newBox[0];
	//	int h = newBox[3] - newBox[2] + 1;
	//	int w = newBox[1] - newBox[0] + 1;

	//	//ROI...
	//	roi = new int*[h];
	//	for (int k = 0; k < h; k++) {
	//		roi[k] = new int[w];
	//		for (int l = 0; l < w; l++) {
	//			roi[k][l] = 255;
	//			int flag = 0;
	//			for (int j = 1; j<int(CodeList[i].size()); j++) {
	//				if (label[k + istart][l + jstart] == CodeList[i][j]) {
	//					flag = 1;
	//					break;
	//				}
	//			}
	//			if (flag)
	//				roi[k][l] = 0;
	//		}
	//	}
	//	Mat cluster;
	//	cluster.create(h, w, CV_8UC1);
	//	arr2mat(roi, cluster, h, w);

	//	//imwrite("Data/ICDAR/Data155/cluster.tif", cluster);
	//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/cluster.tif", cluster);
	//}

	//Remove outliers (points for which MD^2 > 9.210 ie. 99% of confidence interval) from clusters...
	std::vector<std::vector<Point>> NewlineClusters; //after removing outliers...
	for (int i = 0; i<int(lineClusters.size()); i++) {
		std::vector<Point> newCluster; 
		newCluster.push_back(lineClusters[i][0]);
		for (int k = 1; k<int(lineClusters[i].size()); k++) {
			std::vector<Point> cluster;
			for (int t = 1; t<int(lineClusters[i].size()); t++) {
				if (t == k) continue;
				cluster.push_back(lineClusters[i][t]);
			}
			double md = MahalanobisDist(cluster, lineClusters[i][k]);
			if (md*md > 9.21) {
				std::cout << "\n" << lineClusters[i][k] << " is outlier";
			}
			else {
				std::cout << "\n" << lineClusters[i][k] << " is okay";
				newCluster.push_back(lineClusters[i][k]);
			}				
		}
		NewlineClusters.push_back(newCluster);
		newCluster.clear();
	}

	//Chi-square value for each cluster...
	/*std::cout << "\n\n****Chi-square for each cluster:";
	std::vector<double> chi_sq;
	for (int i = 0; i<int(CodeList.size()); i++) {
		chi_sq.push_back(avg_MD(lineClusters[i]));
		std::cout << "\nCluster " << i + 1 << " : " << chi_sq[i];
	}*/

	//*******First level classification....
	int countlbl = 0;
	for (int i = 0; i<int(label_code.size()); i++) {
		if (label_code[i].second == 0)
			countlbl++;
	}
	std::cout << "\n\n" << "Iteration 1: " << countlbl << " labels yet to be classified...";
	//For each un-classified CC make a list of its min perp dist & MD from each cluster code....
	std::vector<std::pair<int, int>> assignedLabels;
	for (int i = 0; i<int(label_code.size()); i++) {
		if (label_code[i].second != 0) continue;

		int lbl = label_code[i].first;
		Point p;
		for (int k = 0; k<int(all_centroids.size()); k++) {
			if (all_centroids[k].first == lbl) {
				p = all_centroids[k].second;
				break;
			}
		}

		//Find MD & min perp-distance from each cluster...
		std::vector<std::pair<int, std::pair<double, double>>> all_dist; //<cluster#, <MD,perpDist>>...
		for (int k = 0; k<int(NewlineClusters.size()); k++) {

			//Find the minimum perpendicular distance from all segments of this cluster....
			double perp_dst = 9999, dst;
			for (int t = 0; t < int(new_lines.size()); t++) {
				if (new_lines[t].first != NewlineClusters[k][0].x) continue;

				int theta = int(atan2(new_lines[t].second.first.y - new_lines[t].second.second.y, new_lines[t].second.first.x - new_lines[t].second.second.x) * 180 / CV_PI);
				dst = perpDist(p, new_lines[t].second.first, theta);
				perp_dst = perp_dst > dst ? dst : perp_dst;
			}
			all_dist.push_back(make_pair(NewlineClusters[k][0].x, make_pair(MahalanobisDist(NewlineClusters[k], p), perp_dst)));
		}
		std::cout << "\n\nLabel " << lbl << p << ":";
		for (int k = 0; k<int(NewlineClusters.size()); k++) {
			std::cout << "\nCode " << all_dist[k].first << " : MD = " << all_dist[k].second.first << " , Perp-Dist = " << all_dist[k].second.second;
			if (all_dist[k].second.first*all_dist[k].second.first > 5.991 && all_dist[k].second.first*all_dist[k].second.first < 13) //chi-sq between 95% - 99.9%...
				std::cout << "---> MD critical region";
			else if (all_dist[k].second.first*all_dist[k].second.first > 13) //chi-sq above 99.9%...
				std::cout << "---> Outlier (Eliminate)";
			if (all_dist[k].second.second > 41)  //  2/3rd of 25.... 
				std::cout << "---> Perp-Distance critical region";
		}

		//find the nearest cluster wrt perp-dist, excluding outliers...
		double min_perp_dst = 9999;
		int min_perp_code = 0, min_pos = 0;
		for (int k = 0; k<int(NewlineClusters.size()); k++) {
			if (all_dist[k].second.first*all_dist[k].second.first > 13) continue; //ignore outliers...

			if (all_dist[k].second.second < min_perp_dst) {
				min_perp_dst = all_dist[k].second.second;
				min_perp_code = all_dist[k].first;
				min_pos = k;
			}
		}
		//if MD^2 is below 99.9% then assign immediately...
		if (all_dist[min_pos].second.first*all_dist[min_pos].second.first < 13) {
			label_code[i].second = min_perp_code;
			std::cout << "\n\nAssigning code " << min_perp_code;
			assignedLabels.push_back(label_code[i]);
		}
		else
		{
			std::cout << "\nNot-assigned";
		}
		
	}

	//for (int i = 0; i<int(label_code.size()); i++) {
	//	if (label_code[i].second == 0) {
	//		int lbl = label_code[i].first;
	//		std::vector<Point> p;
	//		for (int k = 0; k<int(all_centroids.size()); k++) {
	//			if (all_centroids[k].first == lbl) {
	//				p.push_back(all_centroids[k].second);
	//			}
	//		}
	//		/*for (int k = 0; k<int(all_terminals.size()); k++) {
	//			if (all_terminals[k][0].x == lbl) {
	//				for (int t = 0; t < all_terminals[k].size(); t++) {
	//					p.push_back(all_terminals[k][t]);
	//				}
	//			}
	//		}*/
	//		double min_dist1 = 99999;
	//		int min_dist_code1 = 0;
	//		std::cout << "\nLabel " << lbl << p[0] << "---";
	//		for (int t = 0; t<int(lineClusters.size()); t++) {
	//			std::vector<double> all_MD;
	//			double max_dist = -1;
	//			for (int k = 0; k<int(p.size()); k++) {
	//				all_MD.push_back(MahalanobisDist(lineClusters[t], p[k]));
	//				max_dist = all_MD[int(all_MD.size()) - 1] > max_dist ? all_MD[int(all_MD.size()) - 1] : max_dist;
	//			}

	//			std::cout << "\nat max Mahalanobis Distance = " << max_dist << " from line code " << lineClusters[t][0].x << " (chi sq = " << chi_sq[t] << ")";
	//			if (max_dist <= chi_sq[t])
	//				std::cout << "-------> okay";
	//			else
	//				std::cout << "-------> outlier";
	//			if (max_dist < min_dist1) {
	//				min_dist1 = max_dist;
	//				min_dist_code1 = lineClusters[t][0].x;
	//			}
	//		}
	//		label_code[i].second = min_dist_code1;
	//		//std::cout << "\n\n\nLabel " << label_code[i].first << " is near code " << min_dist_code1 << " at Mahalanobis distance = " << min_dist1;
	//		if (min_dist1*min_dist1 < 5.991) {  //chi-sq below 95%....
	//			label_code[i].second = min_dist_code1;
	//			std::cout << "\n\n\nLabel " << label_code[i].first << " is near code " << min_dist_code1 << " at Mahalanobis distance = " << min_dist1 << " (correct range)";
	//		}
	//		else if (min_dist1*min_dist1 < 13.816) { //chi-sq below 99.9%
	//			std::cout << "\n\n\nLabel " << label_code[i].first << " is near code " << min_dist_code1 << " at Mahalanobis distance = " << min_dist1 << " (critical range)";
	//			label_code[i].second = 0;

	//			//Find the minimum perpendicular distance of the centroid from each of the segments...
	//			double perp_dst = 9999, dst, perp_code;
	//			for (int t = 0; t < int(new_lines.size()); t++) {
	//				int theta = int(atan2(new_lines[t].second.first.y - new_lines[t].second.second.y, new_lines[t].second.first.x - new_lines[t].second.second.x) * 180 / CV_PI);
	//				dst = perpDist(p[0], new_lines[t].second.first, theta);
	//				if (perp_dst > dst) {
	//					perp_dst = dst;
	//					perp_code = new_lines[t].first;
	//				}
	//			}
	//			std::cout << "\nLabel " << label_code[i].first << " is near code " << perp_code << " at perpendicular distance = " << perp_dst;
	//			if (perp_dst < 41 && perp_code != min_dist_code1) {
	//				std::cout << "\nclosest line segment is different------re-classify!!";
	//			}
	//			if (perp_code == min_dist_code1) {
	//				std::cout << "\nClassification okay!";
	//				label_code[i].second = min_dist_code1;
	//			}
	//		}
	//		else {
	//			std::cout << "\n\n\nOutlier to the closest cluster code  " << min_dist_code1;
	//			label_code[i].second = 0;
	//		}
	//	}
	//	//check the min distance from every line segment and then check the MD from that cluster; MD^2 has to be < 9.21...
	//}

	//Include the new terminal points and centroids in the clusters...
	for (int i = 0; i<int(assignedLabels.size()); i++) {
		int lbl = assignedLabels[i].first;
		int code = assignedLabels[i].second;
		int ptr;
		for (int k = 0; k<int(all_centroids.size()); k++) {
			if (all_centroids[k].first == lbl) {
				ptr = k;
				break;
			}
		}
		for (int k = 0; k<int(NewlineClusters.size()); k++) {
			if (NewlineClusters[k][0].x == code) {
				std::cout << "\nCluster " << code << ": Previous size = " << int(NewlineClusters[k].size());
				NewlineClusters[k].push_back(all_centroids[ptr].second);
				for (int t = 1; t<int(all_terminals[ptr].size()); t++) {
					NewlineClusters[k].push_back(all_terminals[ptr][t]);
				}
				std::cout << "   New size = " << int(NewlineClusters[k].size());
				break;
			}
		}
	}

	
	//*******Second level classification....

	countlbl = 0;
	for (int i = 0; i<int(label_code.size()); i++) {
		if (label_code[i].second == 0)
			countlbl++;
	}
	std::cout << "\n\n" << "Iteration 2: " << countlbl << " labels yet to be classified...";

	//assign the remaining to the nearest cluster wrt to MD...
	assignedLabels.clear();
	for (int i = 0; i<int(label_code.size()); i++) {
		if (label_code[i].second != 0) continue;

		int lbl = label_code[i].first;
		Point p;
		for (int k = 0; k<int(all_centroids.size()); k++) {
			if (all_centroids[k].first == lbl) {
				p = all_centroids[k].second;
				break;
			}
		}
		double min_MD = 99999, min_code = -1;
		for (int k = 0; k<int(NewlineClusters.size()); k++) {
			double MD = MahalanobisDist(NewlineClusters[k], p);
			std::cout << "\n" << MD << "  Cluster " << NewlineClusters[k][0].x;
			//Exclude if more than 99%...
			if (MD * MD > 9.21)
				std::cout  << "....Outlier";
			else if (MD < min_MD) {
				min_MD = MD;
				min_code = NewlineClusters[k][0].x;
			}
		}
		std::cout << "\n\nLabel " << lbl << p << ": nearest to cluster code " << min_code << " with MD = " << min_MD;
		if (min_code != -1) {
			label_code[i].second = min_code;
			std::cout << "\n\nAssigning code " << min_code;
			assignedLabels.push_back(label_code[i]);
		}
	}

	//Include the new terminal points and centroids in the clusters...
	for (int i = 0; i<int(assignedLabels.size()); i++) {
		int lbl = assignedLabels[i].first;
		int code = assignedLabels[i].second;
		int ptr;
		for (int k = 0; k<int(all_centroids.size()); k++) {
			if (all_centroids[k].first == lbl) {
				ptr = k;
				break;
			}
		}
		for (int k = 0; k<int(NewlineClusters.size()); k++) {
			if (NewlineClusters[k][0].x == code) {
				std::cout << "\nCluster " << code << ": Previous size = " << int(NewlineClusters[k].size());
				NewlineClusters[k].push_back(all_centroids[ptr].second);
				for (int t = 1; t<int(all_terminals[ptr].size()); t++) {
					NewlineClusters[k].push_back(all_terminals[ptr][t]);
				}
				std::cout << "   New size = " << int(NewlineClusters[k].size());
				break;
			}
		}
	}


	//*******Third level classification....(small clusters)

	//re-cluster single labelled clusters...

	std::cout << "\nFinal Label-Code list:";
	for (int i = 0; i<int(label_code.size()); i++) {
		std::cout << "\nLabel " << label_code[i].first << " -- Code " << label_code[i].second;
	}

	std::cout << "\n\nTotal no of clusters(textlines) = " << int(NewlineClusters.size());
	std::cout << "\nNo of Labels in each cluster:";
	std::vector<std::pair<int, int>> clusterSize; //<code, no of labels>...
	for (int i = 0; i<int(NewlineClusters.size()); i++) {
		clusterSize.push_back(make_pair(NewlineClusters[i][0].x, 0));
		for (int t = 0; t<int(label_code.size()); t++) {
			if (label_code[t].second == clusterSize[i].first)
				clusterSize[i].second++;
		}
		std::cout << "\nCluster code " << clusterSize[i].first << " has " << clusterSize[i].second << " labels";
		//reset labels if cluster has atmost 2 labels
		if (clusterSize[i].second <= 2) {
			std::cout << "...removing";
			for (int t = 0; t<int(label_code.size()); t++) {
				if (label_code[t].second == clusterSize[i].first) {
					label_code[t].second = 0;
				}
			}
			clusterSize[i] = make_pair(-1, 0);
			for (int j = i + 1; j<int(NewlineClusters.size());j++)
				NewlineClusters[j - 1] = NewlineClusters[j];
			NewlineClusters.erase(NewlineClusters.end());
			i--;
		}
	}
	std::cout << "\n\nCurrent no of clusters(textlines) = " << int(NewlineClusters.size());

	countlbl = 0;
	for (int i = 0; i<int(label_code.size()); i++) {
		if (label_code[i].second == 0)
			countlbl++;
	}
	std::cout << "\n\n" <<"Iteration 3: " << countlbl << " labels yet to be classified...";
	
	//assign the remaining to the nearest cluster wrt to MD...
	assignedLabels.clear();
	for (int i = 0; i<int(label_code.size()); i++) {
		if (label_code[i].second != 0) continue;

		int lbl = label_code[i].first;
		Point p;
		for (int k = 0; k<int(all_centroids.size()); k++) {
			if (all_centroids[k].first == lbl) {
				p = all_centroids[k].second;
				break;
			}
		}
		double min_MD = 99999, min_code = -1;
		for (int k = 0; k<int(NewlineClusters.size()); k++) {
			if (NewlineClusters[k][0].x == -1) continue; //exclude clusters that are dropped...
			double MD = MahalanobisDist(NewlineClusters[k], p);
			std::cout << "\n" << MD << "  Cluster " << NewlineClusters[k][0].x;
			//Exclude if more than 99%...
			if (MD * MD > 9.21)
				std::cout << "....Outlier";
			else if (MD < min_MD) {
				min_MD = MD;
				min_code = NewlineClusters[k][0].x;
			}
		}
		std::cout << "\n\nLabel " << lbl << p << ": nearest to cluster code " << min_code << " with MD = " << min_MD;
		if (min_code != -1) {
			label_code[i].second = min_code;
			std::cout << "\n\nAssigning code " << min_code;
			assignedLabels.push_back(label_code[i]);
		}
	}

	//Include the new terminal points and centroids in the clusters...
	for (int i = 0; i<int(assignedLabels.size()); i++) {
		int lbl = assignedLabels[i].first;
		int code = assignedLabels[i].second;
		int ptr;
		for (int k = 0; k<int(all_centroids.size()); k++) {
			if (all_centroids[k].first == lbl) {
				ptr = k;
				break;
			}
		}
		for (int k = 0; k<int(NewlineClusters.size()); k++) {
			if (NewlineClusters[k][0].x == code) {
				std::cout << "\nCluster " << code << ": Previous size = " << int(NewlineClusters[k].size());
				NewlineClusters[k].push_back(all_centroids[ptr].second);
				for (int t = 1; t<int(all_terminals[ptr].size()); t++) {
					NewlineClusters[k].push_back(all_terminals[ptr][t]);
				}
				std::cout << "   New size = " << int(NewlineClusters[k].size());
				break;
			}
		}
	}

	//*******Fourth level classification....(final)

	//assign the remaining CCs as individual text lines...

	countlbl = 0;
	for (int i = 0; i<int(label_code.size()); i++) {
		if (label_code[i].second == 0)
			countlbl++;
	}
	std::cout << "\n\n" << "Iteration 4: " << countlbl << " labels yet to be classified...";

	int codeN = 0;
	for (int k = 0; k<int(NewlineClusters.size()); k++) {
		codeN = codeN < NewlineClusters[k][0].x ? NewlineClusters[k][0].x : codeN;
	}

	assignedLabels.clear();
	for (int i = 0; i<int(label_code.size()); i++) {
		if (label_code[i].second != 0) continue;

		int lbl = label_code[i].first;
		Point p;
		int ptr;
		for (int k = 0; k<int(all_centroids.size()); k++) {
			if (all_centroids[k].first == lbl) {
				p = all_centroids[k].second;
				ptr = k;
				break;
			}
		}
		double min_MD = 99999, min_code = -1;
		for (int k = 0; k<int(NewlineClusters.size()); k++) {
			if (NewlineClusters[k][0].x == -1) continue; //exclude clusters that are dropped...
			double MD = MahalanobisDist(NewlineClusters[k], p);
			std::cout << "\n" << MD << "  Cluster " << NewlineClusters[k][0].x;
			//Exclude if more than 99%...
			if (MD * MD > 9.21)
				std::cout << "....Outlier";
			else if (MD < min_MD) {
				min_MD = MD;
				min_code = NewlineClusters[k][0].x;
			}
		}
		std::cout << "\n\nLabel " << lbl << p << ": nearest to cluster code " << min_code << " with MD = " << min_MD;
		if (min_code != -1) {
			label_code[i].second = min_code;
			std::cout << "\n\nAssigning code " << min_code;
			assignedLabels.push_back(label_code[i]);

			//include in the cluster...
			for (int k = 0; k<int(NewlineClusters.size()); k++) {
				if (NewlineClusters[k][0].x == min_code) {
					std::cout << "\nCluster " << min_code << ": Previous size = " << int(NewlineClusters[k].size());
					NewlineClusters[k].push_back(p); //centroid...
					for (int t = 1; t<int(all_terminals[ptr].size()); t++) {
						NewlineClusters[k].push_back(all_terminals[ptr][t]);
					}
					std::cout << "   New size = " << int(NewlineClusters[k].size());
					break;
				}
			}
		}
		else {
			std::cout << "\n\nAssigning NEW code " << ++codeN;
			label_code[i].second = codeN;
			assignedLabels.push_back(label_code[i]);

			//include the NEW cluster...
			std::vector<Point> cluster;
			cluster.push_back(Point(label_code[i].second, 0));
			cluster.push_back(p); //centroid...
			for (int tt = 1; tt<int(all_terminals[ptr].size()); tt++) {
				cluster.push_back(all_terminals[ptr][tt]);
			}
			NewlineClusters.push_back(cluster);
			cluster.clear();
		}
	}


	//display the final list of clusters and labels...
	std::cout << "\n\n\n****************************************************************************";
	std::cout << "\nFinal list of clusters:";
	for (int i = 0; i<int(NewlineClusters.size()); i++) {
		std::cout << "\n " << i << ". Code " << NewlineClusters[i][0].x << " : ";
		for (int t = 0; t<int(label_code.size()); t++) {
			if (label_code[t].second == NewlineClusters[i][0].x)
				std::cout << label_code[t].first << " , ";
		}
	}
	std::cout << "\n\n" << int(NewlineClusters.size()) << " text lines have been identified!";
	std::cout << "\n****************************************************************************";


	

	connections.create(r, c, CV_8UC1);
	cvtColor(img, connections, CV_GRAY2BGR);
	
	//color the labels that are collinear...
	//define colors...
	std::vector<Vec3b> colors; //<col0,col1,col2,...col10>...
	colors.push_back(Vec3b(50, 255, 100)); //col0 green
	colors.push_back(Vec3b(255, 255, 0)); //col1 orange
	colors.push_back(Vec3b(200, 150, 255)); //col2 pink
	colors.push_back(Vec3b(50, 100, 200)); //col3 orange
	colors.push_back(Vec3b(150, 50, 200)); //col4 pink
	colors.push_back(Vec3b(255, 100, 0)); // blue
	colors.push_back(Vec3b(50, 120, 50)); //col6 green
	colors.push_back(Vec3b(200, 50, 100)); //col7 blue
	colors.push_back(Vec3b(0, 255, 255)); //col8 yellow
	colors.push_back(Vec3b(0, 100, 255)); //col9 red
	colors.push_back(Vec3b(150, 150, 150)); //col10
	std::cout << "\nLine-Color matching--";

	//color the unmarked labels black...
	std::vector<int> unmarked_lbl;
	for (int k = 0; k<int(label2lineNew.size()); k++) {
		if (label2lineNew[k].first.second == 0) {
			unmarked_lbl.push_back(label2lineNew[k].first.first);
			int lbl = label2lineNew[k].first.first;
			for (int h = box[lbl][4]; h <= box[lbl][5]; h++) {
				for (int w = box[lbl][2]; w <= box[lbl][3]; w++) {
					if (label[h][w] == lbl) {
						text_lines.at<Vec3b>(h, w) = colors[10];
						connections.at<Vec3b>(h, w) = colors[10];
					}
				}
			}
		}
	}

	int *color_counter = new int[int(label_code.size())];
	for (int i = 0; i<int(label_code.size()); i++)
		color_counter[i] = 0;
	int counter = -1;
	for (int i = 0; i<int(label_code.size()); i++) {
		if (color_counter[i] == 1 || label_code[i].second == 0) continue; //already colored...

		std::vector<int> collinear_lbl;
		std::cout << "\nLabel " << label_code[i].first << " Code "<< label_code[i].second << " = color " << ++counter;
		int col_no = counter % 10;
		//collect the labels with same code no...
		for (int k = 0; k<int(label_code.size()); k++) {
			if (label_code[k].second == label_code[i].second) {
				collinear_lbl.push_back(label_code[k].first);
				color_counter[k] = 1;
			}
		}
		int newBox[4] = { 99999,0,99999,0 };
		for (int k = 0; k<int(collinear_lbl.size()); k++) {
			int lbl = collinear_lbl[k];
			newBox[0] = newBox[0] > box[lbl][2] ? box[lbl][2] : newBox[0];
			newBox[1] = newBox[1] < box[lbl][3] ? box[lbl][3] : newBox[1];
			newBox[2] = newBox[2] > box[lbl][4] ? box[lbl][4] : newBox[2];
			newBox[3] = newBox[3] < box[lbl][5] ? box[lbl][5] : newBox[3];
		}
		for (int h = newBox[2]; h <= newBox[3]; h++) {
			for (int w = newBox[0]; w <= newBox[1]; w++) {
				for (int k = 0; k<int(collinear_lbl.size()); k++) {
					if (label[h][w] == collinear_lbl[k]) {
						text_lines.at<Vec3b>(h, w) = colors[col_no];
						connections.at<Vec3b>(h, w) = colors[col_no];
					}
				}
			}
		}
	}
	
	//mark the connections in red...
	/*for (int i = 0; i<int(new_lines.size()); i++) {
		line(connections, new_lines[i].second.first, new_lines[i].second.second, Scalar(0, 0, 255), 3, 8);
	}*/

	//imwrite("Data/ICDAR/Data155/text_lines.tif", text_lines);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/text_lines.tif", text_lines);

	//imwrite("Data/ICDAR/Data155/connections.tif", connections);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/connections.tif", connections);
	

	
	return(connectCentroid);
}

void invCCProj(Mat img, int flag = 0) {
	int **orgimg, **label, **box, **binimg, **thin_img;
	Mat colorimg, bin_img, thin_im, near4img, text_lines;

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
	text_lines.create(r, c, CV_8UC1);
	arr2mat(thin_img, text_lines);
	cvtColor(text_lines, text_lines, CV_GRAY2BGR);

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
		//imwrite("Data/ICDAR/Data155/binImg15.tif", bin_img);
		imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/binImg15.tif", bin_img);

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
	//imwrite("Data/ICDAR/Data155/invcc.tif", img);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/invcc.tif", img);

	/*Mat rotimg;
	imageRotate(img, rotimg, -30);*/

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

	std::vector<std::pair<Point, Point>> connectCC1, connectCC2;
	std::vector<std::pair<Point, Point>> connectCentroid = tracing0(img, colorimg, text_lines, thin_img, connectCC1, connectCC2);


	for(int i = 0; i<int(connectCentroid.size()); i++)
		arrowedLine(colorimg, connectCentroid[i].first, connectCentroid[i].second, Scalar(255, 0, 255), 1.5, 8);

	for (int i = 0; i<int(connectCC1.size()); i++) {
		//line(colorimg, connectCC[i].first, connectCC[i].second, Scalar(255, 0, 255), 2, 8);
		arrowedLine(colorimg, connectCC1[i].first, connectCC1[i].second, Scalar(0, 120, 255), 2.5, 8);
	}

	for (int i = 0; i<int(connectCC2.size()); i++) {
		//line(colorimg, connectCC[i].first, connectCC[i].second, Scalar(255, 0, 255), 2, 8);
		//arrowedLine(colorimg, connectCentroid[i].first, connectCentroid[i].second, Scalar(255, 0, 255), 1.5, 8);
		arrowedLine(colorimg, connectCC2[i].first, connectCC2[i].second, Scalar(120, 255, 0), 2.5, 8);
	}

	/*for (int i = 0; i < labelNum; i++) {
		if (box[i][2] == -1) continue;
		line(colorimg, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(0, 0, 255), 2, 8);
		line(colorimg, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 2, 8);
		line(colorimg, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(0, 0, 255), 2, 8);
		line(colorimg, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 2, 8);
	}*/
	imshow("boundingBox.tif", colorimg);
	//imwrite("Data/ICDAR/Data155/centroidBoundingBoxred.tif", colorimg);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/centroidBoundingBoxred.tif", colorimg);
	
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
	drawAxis(col_img, cntr, p1, Scalar(0, 0, 255), 2);
	drawAxis(col_img, cntr, p2, Scalar(255, 0, 0), 2);
	double angle = atan2(eig_vec[0][1], eig_vec[0][0]); // orientation in radians
	std::cout << "\nOrientation = " << angle * 180 / CV_PI;
	imshow("orientations.jpg",col_img);
	//imwrite("Data/ICDAR/Data039/directions.tif", col_img);
	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/directions.tif", col_img);
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
		//imwrite("Data/ICDAR/Data155/binImg.tif", bin_img);
		imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan039/binImg.tif", bin_img);
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
	//Mat scanimg = imread("Data/ICDAR/Data155/155.tif", IMREAD_GRAYSCALE);
	
	Mat scanimg = imread("Data/MyScans/Scanned Data/Scanned Data/Scan039/scan039.tif", IMREAD_GRAYSCALE);
	
	imshow("Image", scanimg);
	std::cout << "Image dimensions:" << scanimg.rows << "x" << scanimg.cols << endl;
	
	r = scanimg.rows;
	c = scanimg.cols;

	invCCProj(scanimg, 1);
	//imgPCA(scanimg, 1);
	
	std::cout << "\nDone!\n";
	waitKey(0);
	system("pause");
	return 0;
	
}