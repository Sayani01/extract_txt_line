#include <stdio.h>
#include<conio.h>
#include<iostream>
#include <vector>
#include <string>
#include "opencv2/imgproc/imgproc.hpp"
//#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "VarDeclare.h"
#include "Projections.h"
#include "ConnectedComponent.h"

using namespace std;
using namespace cv;
namespace cv {
	using std::vector;
};


//Standard deviation...
double stddev(int *x, int sz){
	double X = 0, sd = 0;
	for (int i = 0; i < sz; i++)
		X += x[i];
	X /= sz;

	for (int i = 0; i < sz; i++) {
		sd += ((x[i] - X)*(x[i] - X));
	}
	sd = sqrt(sd / sz);
	return(sd);
}

double stddev(double *x, int sz) {
	double X = 0, sd = 0;
	for (int i = 0; i < sz; i++)
		X += x[i];
	X /= sz;

	for (int i = 0; i < sz; i++) {
		sd += ((x[i] - X)*(x[i] - X));
	}
	sd = sqrt(sd / sz);
	return(sd);
}

double stddev(std::vector<int> x, int sz) {
	double X = 0, sd = 0;
	for (int i = 0; i < sz; i++)
		X += x[i];
	X /= sz;

	for (int i = 0; i < sz; i++) {
		sd += ((x[i] - X)*(x[i] - X));
	}
	sd = sqrt(sd / sz);
	return(sd);
}

double AvgOfDerivative(int *vec, int sz) {
	int s = sz - 1;
	double *diff = new double[s];
	for (int i = 1; i < sz; i++)
		diff[i - 1] = vec[i] - vec[i - 1];

	return(stddev(diff,s));
}

double AvgOfDerivative(std::vector<int> vec, int sz) {
	int s = sz - 1;
	double *diff = new double[s];
	for (int i = 1; i < sz; i++)
		diff[i - 1] = vec[i] - vec[i - 1];

	return(stddev(diff, s));
}

//find projection in horizontal direction...
double project0(int **roi, int h, int w, int *hist){
	for (int i = 0; i < h; i++) {
		int prev = 255;
		int count = 0;
		for (int j = 0; j < w; j++) {
			if (roi[i][j] == 0)
				hist[i]++;
			if (prev == 255 && roi[i][j] == 0)
				count++;
			prev = roi[i][j];
		}
		hist[i] *= count;
	}
	std::cout << "\nStd dev(0) = " << stddev(hist, h);
	return(stddev(hist, h));
}

//find projection in vertical direction...
double project90(int **roi, int h, int w, int *hist){
	for (int j = 0; j < w; j++) {
		int prev = 255;
		int count = 0;
		for (int i = 0; i < h; i++) {
			if (roi[i][j] == 0)
				hist[j]++;
			if (prev == 255 && roi[i][j] == 0)
				count++;
			prev = roi[i][j];
		}
		hist[j] *= count;
	}
	std::cout << "\nStd dev(180) = " << stddev(hist, w);
	return(stddev(hist, w));
}

//find projection in off-diagonal direction...
double project45(int **img, int **label, double *dim, int lbl, int *hist){
	Point N = Point(int(dim[1]), int(dim[0]));
	Point S = Point(int(dim[3]), int(dim[2]));
	Point E = Point(int(dim[5]), int(dim[4]));
	Point W = Point(int(dim[7]), int(dim[6]));

	int sz = 2*(W.y - N.y + 1);
	int C = W.y - W.x; //c = W.y - W.x...
	int ilimit = (S.y - W.y + 1);
	int istart = W.y, jstart = W.x, adder = 1;

	
	for (int k = 0; k < sz; k++) {
		//std::cout << "\nhere c=" << C - k << "i:" << istart << " - " << istart + ilimit;
		int prev = 255;
		int count = 0;
		for (int i = istart; i <= istart + ilimit; i++) {
			int j = i - (C - k); //y - x = c...
			if (i < 0 || i >= r || j < 0 || j >= c) {
				continue;
			}
			if (label[i][j] == lbl)
				hist[k]++;
			if (prev == 255 && label[i][j] == lbl)
				count++;
			prev = img[i][j];
		}
		adder++;
		if (adder % 2 == 0) //even run...
			istart--;
		else      //odd run...
			jstart++;
		hist[k] *= count;
	}
	/*std::cout << "Projections in 135deg:";
	for (int k = 0; k < sz; k++)
	std::cout << hist[k] << " , ";
	*/

	std::cout << "\nStd dev(45) = " << stddev(hist, sz);
	return(stddev(hist, sz));
}

//find projection in diagonal direction...
double project135(int **img, int **label, double *dim, int lbl, int *hist){
	Point N = Point(int(dim[1]), int(dim[0]));
	Point S = Point(int(dim[3]), int(dim[2]));
	Point E = Point(int(dim[5]), int(dim[4]));
	Point W = Point(int(dim[7]), int(dim[6]));
	
	int sz = 2*(S.y - W.y + 1);
	int C = N.x + N.y; //c = N.x + N.y...
	int ilimit = (W.y - N.y + 1);
	int istart = W.y, jstart = W.x, adder = 1;
	int perp = C - (S.x + S.y);
	
	for (int k = 0; k < sz; k++) {
		int prev = 255;
		int count = 0;
		for (int i = istart; i >= istart - ilimit; i--) {
			int j = (C + k) - i; //x + y = c...
			if (i < 0 || i >= r || j < 0 || j >= c) {
				continue;
			}
			if (label[i][j] == lbl) 
				hist[k]++;
			if (prev == 255 && label[i][j] == lbl)
				count++;
			prev = img[i][j];
		}
		adder++;
		if (adder % 2 == 0) //even run...
			istart++;
		else      //odd run...
			jstart++; 	
		hist[k] *= count;
		//hist[k] %= perp;
	}
	/*std::cout << "Projections in 135deg:";
	for (int k = 0; k < sz; k++)
		std::cout << hist[k] << " , ";
*/

	std::cout << "\nStd dev(135) = " << stddev(hist, sz);
	return(stddev(hist, sz));

}

int isInteriorPoint(int img[7][7]) {
	int count = 1;
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < 7; j++) {
			count *= (1 - img[i][j]);
		}
	}
	return(count);
}

void projectAllCross(int **img, int h, int w, int &inclineCrossing) {
	const int thres = 7;
	int incline = 0, maxGrad = 0;
	double slope1, slope2, intercept1, intercept2;
	double **features = new double*[36];
	for (int i = 0; i < 36; i++) {
		features[i] = new double[2];
		if (i % 2 == 0)
			features[i][0] = 5 * (i / 2);
		else
			features[i][0] = 90 + 5 * ((i - 1) / 2);
		for (int j = 1; j < 2; j++) {
			features[i][j] = -1;
		}
	}

	//project horizontally....
	int crossCount, maxCross = -1;
	for (int i = 0; i < h; i++) {
		crossCount = 0;
		int flag = 1;
		int gap = 0;
		int prev = 1, curr, whitecount[3] = {0,0,0}, x = 0;
		for (int j = 0; j < w; j++) {
			if (crossCount > 0) gap++;
			if (crossCount > 0 && img[i][j] == 1) { //in between crossings
				whitecount[x]++;
			}
			if (img[i][j] == 0 && prev == 1) { //white->black 
				if (crossCount == 0) //first crossing
					crossCount++;
				else if (whitecount[x] >= thres || gap >= 10) { //ideal crossing...reset
					crossCount++;
					gap = 0;
					whitecount[0] = 0;
					whitecount[1] = 0;
					whitecount[2] = 0;
					x = 0;
				}
				else if (x == 0 && whitecount[x] < thres) //first close crossing
					x = 1;
				else if (x == 1 && whitecount[x] < thres) //second close crossing
					x = 2;
				else if ((x == 2 && whitecount[x] < thres) || gap < 10) //third close crossing
					flag = 0;
			}
			if (!flag)
				break; //do not include line
			prev = img[i][j];
		}
		if (crossCount > maxCross && flag)
			maxCross = crossCount;
	}

	features[0][1] = maxCross;

	//project vertically...
	maxCross = -1;
	for (int j = 0; j < w; j++) {
		crossCount = 0;
		int flag = 1;
		int gap = 0;
		int prev = 1, curr, whitecount[3] = { 0,0,0 }, x = 0;
		for (int i = 1; i < h; i++) {
			if (crossCount > 0) gap++;
			if (crossCount > 0 && img[i][j] == 1) { //in between crossings
				whitecount[x]++;
			}
			if (img[i][j] == 0 && prev == 1) { //white->black 
				if (crossCount == 0) //first crossing
					crossCount++;
				else if (whitecount[x] >= thres || gap >= 10) { //ideal crossing...reset
					crossCount++;
					gap = 0;
					whitecount[0] = 0;
					whitecount[1] = 0;
					whitecount[2] = 0;
					x = 0;
				}
				else if (x == 0 && whitecount[x] < thres) //first close crossing
					x = 1;
				else if (x == 1 && whitecount[x] < thres) //second close crossing
					x = 2;
				else if ((x == 2 && whitecount[x] < thres) || gap < 10) //third close crossing
					flag = 0;
			}
			if (!flag)
				break; //do not include line
			prev = img[i][j];
		}
		if (crossCount > maxCross && flag)
			maxCross = crossCount;
	}

	features[1][1] = maxCross;

	for (int K = 1; K <= 17; K++) {
		int count = 0;
		intercept1 = -1;

		int theta = 5 * K;
		slope1 = tan((theta + 90)*CV_PI / 180);
		slope2 = tan((theta)*CV_PI / 180);
		//std::cout << "\ntheta = " << theta << " , theta - 90 = " << theta - 90 << " , slope1 = " << slope1 << " , slope2 = " << slope2;

		//projection from 95 - 175 degrees...
		maxCross = -1;
		int flag = 0;
		do {
			intercept1++;
			for (int y = 0; y < h; y++) {
				int x = int((y - intercept1) / slope1);
				if (x >= w || x < 0) continue;
				if (img[y][x] == 0) {
					flag = 1;
					break;
				}
			}
			//std::cout << "\nintercept = " << intercept1;
		} while (flag == 0);
		do {
			//std::cout << "\nintercept = " << intercept1;
			count = 0, crossCount = 0;
			int fl = 1;
			int gap = 0;
			int prev = 1, curr, whitecount[3] = { 0,0,0 }, X = 0;
			for (int y = 0; y < h; y++) {
				int x = int((y - intercept1) / slope1);
				if (x >= w || x < 0) continue;
				count += 1 - img[y][x];
				if (crossCount > 0) gap++;
				if (crossCount > 0 && img[y][x] == 1) { //in between crossings
					whitecount[X]++;
				}
				if (img[y][x] == 0 && prev == 1) { //white->black 
					if (crossCount == 0) //first crossing
						crossCount++;
					else if (whitecount[X] >= thres || gap >= 10) { //ideal crossing...reset
						crossCount++;
						gap = 0;
						whitecount[0] = 0;
						whitecount[1] = 0;
						whitecount[2] = 0;
						X = 0;
					}
					else if (X == 0 && whitecount[X] < thres) //first close crossing
						X = 1;
					else if (X == 1 && whitecount[X] < thres) //second close crossing
						X = 2;
					else if ((X == 2 && whitecount[X] < thres) || gap < 10) //second close crossing
						fl = 0;
				}
				if (!fl)
					break; //do not include line
				prev = img[y][x];
			}
			if (count == 0)
				break;
			if (crossCount > maxCross && fl)
				maxCross = crossCount;
			intercept1++;
		} while (count > 0);

		features[K * 2 + 1][1] = maxCross;

		//projection from 5 - 85 degrees...
		intercept2 = h - 1, maxCross = -1;
		flag = 0;
		do {
			intercept2--;
			for (int y = h - 1; y >= 0; y--) {
				int x = (y - intercept2) / slope2;
				if (x >= w || x < 0) continue;
				if (img[y][x] == 0) {
					flag = 1;
					break;
				}
			}
		} while (flag == 0);
		do {
			count = 0, crossCount = 0;
			int fl = 1;
			int gap = 0;
			int prev = 1, curr, whitecount[3] = { 0,0,0 }, X = 0;
			for (int y = h - 1; y >= 0; y--) {
				int x = (y - intercept2) / slope2;
				if (x >= w || x < 0) continue;
				count += 1 - img[y][x];
				if (crossCount > 0) gap++;
				if (crossCount > 0 && img[y][x] == 1) { //in between crossings
					whitecount[X]++;
				}
				if (img[y][x] == 0 && prev == 1) { //white->black 
					if (crossCount == 0) //first crossing
						crossCount++;
					else if (whitecount[X] >= thres || gap >= 10) { //ideal crossing...reset
						crossCount++;
						gap = 0;
						whitecount[0] = 0;
						whitecount[1] = 0;
						whitecount[2] = 0;
						X = 0;
					}
					else if (X == 0 && whitecount[X] < thres) //first close crossing
						X = 1;
					else if (X == 1 && whitecount[X] < thres) //second close crossing
						X = 2;
					else if ((X == 2 && whitecount[X] < thres) || gap < 10) //second close crossing
						fl = 0;
				}
				if (!fl)
					break; //do not include line
				prev = img[y][x];
			}
			if (count == 0)
				break;
			if (crossCount > maxCross && fl)
				maxCross = crossCount;
			intercept2--;
		} while (count > 0);

		features[K * 2][1] = maxCross;
	}

	double maxFeature1 = -1, maxFeature2 = -1, maxFeature3 = -1;
	std::cout << "\nFeatures List (angle/std dev/max/bw-cross):";
	for (int i = 0; i < 36; i++) {
		std::cout << "\n" << features[i][0] << "\t" << features[i][1] ;
		
		if (maxFeature3 <= features[i][1]) {
			maxFeature3 = features[i][1];
			inclineCrossing = features[i][0];
		}
	}
	std::cout << "\nInclined at " << inclineCrossing;
}

void projectAllStdDev(int **img, int h, int w, int &inclineStdDev) {
	int incline = 0, maxGrad = 0;
	double slope1, slope2, intercept1, intercept2;
	double **features = new double*[36];
	for (int i = 0; i < 36; i++) {
		features[i] = new double[2];
		if (i % 2 == 0)
			features[i][0] = 5 * (i / 2);
		else
			features[i][0] = 90 + 5 * ((i - 1) / 2);
		for (int j = 1; j < 2; j++) {
			features[i][j] = -1;
		}
	}

	//project horizontally....
	double dev;
	int *hist = new int[h];
	for (int i = 0; i < h; i++)
		hist[i] = 0;
	dev = project0(img, h, w, hist);
	features[0][1] = dev;
	delete[] hist;

	//project vertically...
	hist = new int[w];
	for (int i = 0; i < w; i++)
		hist[i] = 0;
	dev = project90(img, h, w, hist);
	features[1][1] = dev;

	for (int K = 1; K <= 17; K++) {
		int count = 0;
		intercept1 = -1;

		int theta = 5 * K;
		slope1 = tan((theta + 90)*CV_PI / 180);
		slope2 = tan((theta)*CV_PI / 180);
		//std::cout << "\ntheta = " << theta << " , theta - 90 = " << theta - 90 << " , slope1 = " << slope1 << " , slope2 = " << slope2;

		//projection from 95 - 175 degrees...
		std::vector<int> hist1;
		int flag = 0;
		do {
			intercept1++;
			for (int y = 0; y < h; y++) {
				int x = int((y - intercept1) / slope1);
				if (x >= w || x < 0) continue;
				if (img[y][x] == 0) {
					flag = 1;
					break;
				}
			}
			//std::cout << "\nintercept = " << intercept1;
		} while (flag == 0);
		do {
			//std::cout << "\nintercept = " << intercept1;
			count = 0;
			int prev = 0, curr;
			Point first, last;
			for (int y = 0; y < h; y++) {
				int x = int((y - intercept1) / slope1);
				if (x >= w || x < 0) continue;
				count += 1 - img[y][x];
			}
			if (count == 0)
				break;
			hist1.push_back(count);
			intercept1++;
		} while (count > 0);
		dev = stddev(hist1, int(hist1.size()));

		features[K * 2 + 1][1] = dev;

		//projection from 5 - 85 degrees...
		intercept2 = h - 1;
		count = 0;
		std::vector<int> hist2;
		flag = 0;
		do {
			intercept2--;
			for (int y = h - 1; y >= 0; y--) {
				int x = (y - intercept2) / slope2;
				if (x >= w || x < 0) continue;
				if (img[y][x] == 0) {
					flag = 1;
					break;
				}
			}
		} while (flag == 0);
		do {
			count = 0;
			int prev = 0, curr;
			Point first, last;
			for (int y = h - 1; y >= 0; y--) {
				int x = (y - intercept2) / slope2;
				if (x >= w || x < 0) continue;
				count += 1 - img[y][x];
			}
			if (count == 0)
				break;
			hist2.push_back(count);
			intercept2--;
		} while (count > 0);
		dev = stddev(hist2, int(hist2.size()));

		features[K * 2][1] = dev;
	}

	double maxFeature1 = -1, maxFeature2 = -1, maxFeature3 = -1;
	std::cout << "\nFeatures List (angle/std dev/max/bw-cross):";
	for (int i = 0; i < 36; i++) {
		std::cout << "\n" << features[i][0] << "\t" << features[i][1];

		if (maxFeature3 < features[i][1]) {
			maxFeature3 = features[i][1];
			inclineStdDev = features[i][0];
		}
	}
	std::cout << "\nInclined at " << inclineStdDev;
}

void projectAllRange(int **img, int h, int w, int &inclineCrossing) {
	const int thres = 8;
	int incline = 0, maxGrad = 0;
	double slope1, slope2, intercept1, intercept2;
	double **features = new double*[36];
	for (int i = 0; i < 36; i++) {
		features[i] = new double[2];
		if (i % 2 == 0)
			features[i][0] = 5 * (i / 2);
		else
			features[i][0] = 90 + 5 * ((i - 1) / 2);
		for (int j = 1; j < 2; j++) {
			features[i][j] = -1;
		}
	}

	//project horizontally....
	int crossCount, maxCross = -1;
	for (int i = 0; i < h; i++) {
		crossCount = 0;
		int flag = 1;
		int prev = 1, curr, whitecount[3] = { 0,0,0 }, x = 0;
		for (int j = 0; j < w; j++) {
			if (crossCount > 0 && img[i][j] == 1) { //in between crossings
				whitecount[x]++;
			}
			if (img[i][j] == 0 && prev == 1) { //white->black 
				if (crossCount == 0) //first crossing
					crossCount++;
				else if (whitecount[x] >= thres) { //ideal crossing...reset
					crossCount++;
					whitecount[0] = 0;
					whitecount[1] = 0;
					whitecount[2] = 0;
					x = 0;
				}
				else if (x == 0 && whitecount[x] < thres) //first close crossing
					x = 1;
				else if (x == 1 && whitecount[x] < thres) //second close crossing
					x = 2;
				else if (x == 2 && whitecount[x] < thres) //third close crossing
					flag = 0;
			}
			if (!flag)
				break; //do not include line
			prev = img[i][j];
		}
		if (crossCount > maxCross && flag)
			maxCross = crossCount;
	}

	features[0][1] = maxCross;

	//project vertically...
	maxCross = -1;
	for (int j = 0; j < w; j++) {
		crossCount = 0;
		int flag = 1;
		int prev = 1, curr, whitecount[3] = { 0,0,0 }, x = 0;
		for (int i = 1; i < h; i++) {
			if (crossCount > 0 && img[i][j] == 1) { //in between crossings
				whitecount[x]++;
			}
			if (img[i][j] == 0 && prev == 1) { //white->black 
				if (crossCount == 0) //first crossing
					crossCount++;
				else if (whitecount[x] >= thres) { //ideal crossing...reset
					crossCount++;
					whitecount[0] = 0;
					whitecount[1] = 0;
					whitecount[2] = 0;
					x = 0;
				}
				else if (x == 0 && whitecount[x] < thres) //first close crossing
					x = 1;
				else if (x == 1 && whitecount[x] < thres) //second close crossing
					x = 2;
				else if (x == 2 && whitecount[x] < thres) //third close crossing
					flag = 0;
			}
			if (!flag)
				break; //do not include line
			prev = img[i][j];
		}
		if (crossCount > maxCross && flag)
			maxCross = crossCount;
	}

	features[1][1] = maxCross;

	for (int K = 1; K <= 17; K++) {
		int count = 0;
		intercept1 = -1;

		int theta = 5 * K;
		slope1 = tan((theta + 90)*CV_PI / 180);
		slope2 = tan((theta)*CV_PI / 180);
		//std::cout << "\ntheta = " << theta << " , theta - 90 = " << theta - 90 << " , slope1 = " << slope1 << " , slope2 = " << slope2;

		//projection from 95 - 175 degrees...
		maxCross = -1;
		int flag = 0;
		do {
			intercept1++;
			for (int y = 0; y < h; y++) {
				int x = int((y - intercept1) / slope1);
				if (x >= w || x < 0) continue;
				if (img[y][x] == 0) {
					flag = 1;
					break;
				}
			}
			//std::cout << "\nintercept = " << intercept1;
		} while (flag == 0);
		do {
			//std::cout << "\nintercept = " << intercept1;
			count = 0, crossCount = 0;
			int fl = 1;
			int prev = 1, curr, whitecount[3] = { 0,0,0 }, X = 0;
			for (int y = 0; y < h; y++) {
				int x = int((y - intercept1) / slope1);
				if (x >= w || x < 0) continue;
				count += 1 - img[y][x];
				if (crossCount > 0 && img[y][x] == 1) { //in between crossings
					whitecount[X]++;
				}
				if (img[y][x] == 0 && prev == 1) { //white->black 
					if (crossCount == 0) //first crossing
						crossCount++;
					else if (whitecount[X] >= thres) { //ideal crossing...reset
						crossCount++;
						whitecount[0] = 0;
						whitecount[1] = 0;
						whitecount[2] = 0;
						X = 0;
					}
					else if (X == 0 && whitecount[X] < thres) //first close crossing
						X = 1;
					else if (X == 1 && whitecount[X] < thres) //second close crossing
						X = 2;
					else if (X == 2 && whitecount[X] < thres) //second close crossing
						fl = 0;
				}
				if (!fl)
					break; //do not include line
				prev = img[y][x];
			}
			if (count == 0)
				break;
			if (crossCount > maxCross && fl)
				maxCross = crossCount;
			intercept1++;
		} while (count > 0);

		features[K * 2 + 1][1] = maxCross;

		//projection from 5 - 85 degrees...
		intercept2 = h - 1, maxCross = -1;
		flag = 0;
		do {
			intercept2--;
			for (int y = h - 1; y >= 0; y--) {
				int x = (y - intercept2) / slope2;
				if (x >= w || x < 0) continue;
				if (img[y][x] == 0) {
					flag = 1;
					break;
				}
			}
		} while (flag == 0);
		do {
			count = 0, crossCount = 0;
			int fl = 1;
			int prev = 1, curr, whitecount[3] = { 0,0,0 }, X = 0;
			for (int y = h - 1; y >= 0; y--) {
				int x = (y - intercept2) / slope2;
				if (x >= w || x < 0) continue;
				count += 1 - img[y][x];
				if (crossCount > 0 && img[y][x] == 1) { //in between crossings
					whitecount[X]++;
				}
				if (img[y][x] == 0 && prev == 1) { //white->black 
					if (crossCount == 0) //first crossing
						crossCount++;
					else if (whitecount[X] >= thres) { //ideal crossing...reset
						crossCount++;
						whitecount[0] = 0;
						whitecount[1] = 0;
						whitecount[2] = 0;
						X = 0;
					}
					else if (X == 0 && whitecount[X] < thres) //first close crossing
						X = 1;
					else if (X == 1 && whitecount[X] < thres) //second close crossing
						X = 2;
					else if (X == 2 && whitecount[X] < thres) //second close crossing
						fl = 0;
				}
				if (!fl)
					break; //do not include line
				prev = img[y][x];
			}
			if (count == 0)
				break;
			if (crossCount > maxCross && fl)
				maxCross = crossCount;
			intercept2--;
		} while (count > 0);

		features[K * 2][1] = maxCross;
	}

	double maxFeature1 = -1, maxFeature2 = -1, maxFeature3 = -1;
	std::cout << "\nFeatures List (angle/std dev/max/bw-cross):";
	for (int i = 0; i < 36; i++) {
		std::cout << "\n" << features[i][0] << "\t" << features[i][1];

		if (maxFeature3 < features[i][1]) {
			maxFeature3 = features[i][1];
			inclineCrossing = features[i][0];
		}
	}
	std::cout << "\nInclined at " << inclineCrossing;
}


void projectAll(int **img, int h, int w, int &inclineStdDev, int &inclineExtent, int &inclineCrossing) {	
	const int Nb = 7;
	int incline = 0, maxGrad = 0;
	double extent, maxExtent;
	double slope1, slope2, intercept1, intercept2, dev, maxDev, gradient;
	double **features = new double*[36];
	for (int i = 0; i < 36; i++) {
		features[i] = new double[4];
		if (i % 2 == 0)
			features[i][0] = 5 * (i / 2);
		else
			features[i][0] = 90 + 5 * ((i - 1) / 2);
		for (int j = 1; j < 4; j++) {
			features[i][j] = -1;
		}
	}

	//project horizontally....
	int min = 9999, max = -1, crossCount, maxCross = -1;
	maxExtent = -1;
	int *hist = new int[h];
	for (int i = 0; i < h; i++)
		hist[i] = 0;
	dev = project0(img, h, w, hist);
	for (int i = 0; i < h; i++) {
		if (hist[i] < min) min = hist[i];
		if (hist[i] > max) max = hist[i];
	}
	for (int i = 0; i < h; i++) {
		crossCount = 0;
		int prev = 0, curr;
		Point first, last;
		for (int j = 0; j < w; j++) {
			int nbrs[Nb][Nb];
			for (int k = i - 3; k <= i + 3; k++) {
				for (int l = j - 3; l <= j + 3; l++) {
					if (k<0 || k>=h || l<0 || l>=w)
						nbrs[k - i + 3][l - j + 3] = 0;
					else
						nbrs[k - i + 3][l - j + 3] = img[k][l];
				}
			}
			curr = isInteriorPoint(nbrs);
			if (img[i][j] == 0 && crossCount == 0)
				first = Point(j, i);

			if(curr == 1 && prev == 0) { //white->black 
				crossCount++;
			}
			if (img[i][j] == 0)
				last = Point(j, i);
			prev = curr;
		}
		extent = norm(last - first);
		if (extent > maxExtent)
			maxExtent = extent;
		if (crossCount > maxCross)
			maxCross = crossCount;
	}

	features[0][1] = dev;
	features[0][2] = maxExtent;
	features[0][3] = maxCross;
	delete[] hist;

	//project vertically...
	min = 9999, max = -1, maxCross = -1, maxExtent = -1;
	hist = new int[w];
	for (int i = 0; i < w; i++)
		hist[i] = 0;
	dev = project90(img, h, w, hist);
	for (int i = 0; i < w; i++) {
		if (hist[i] < min) min = hist[i];
		if (hist[i] > max) max = hist[i];
	}
	for (int j = 0; j < w; j++ ) {
		crossCount = 0;
		int prev = 0, curr;
		Point first, last;
		for (int i = 1; i < h; i++) {
			int nbrs[Nb][Nb];
			for (int k = i - 3; k <= i + 3; k++) {
				for (int l = j - 3; l <= j + 3; l++) {
					if (k<0 || k >= h || l<0 || l >= w)
						nbrs[k - i + 3][l - j + 3] = 0;
					else
						nbrs[k - i + 3][l - j + 3] = img[k][l];
				}
			}
			curr = isInteriorPoint(nbrs);
			if (img[i][j] == 0 && crossCount == 0)
				first = Point(j, i);

			if (curr == 1 && prev == 0) { //white->black 
				crossCount++;
			}
			
			if (img[i][j] == 0)
				last = Point(j, i);
			prev = curr;
		}
		extent = norm(last - first);
		if (extent > maxExtent)
			maxExtent = extent;
		if (crossCount > maxCross)
			maxCross = crossCount;
	}

	features[1][1] = dev;
	features[1][2] = maxExtent;
	features[1][3] = maxCross;
	delete[] hist;

	for (int K = 1; K <= 17; K++) {
		int count = 0;
		min = 9999, max = -1;
		intercept1 = -1;

		int theta = 5 * K;
		slope1 = tan((theta + 90)*CV_PI / 180);
		slope2 = tan((theta)*CV_PI / 180);
		//std::cout << "\ntheta = " << theta << " , theta - 90 = " << theta - 90 << " , slope1 = " << slope1 << " , slope2 = " << slope2;

		//projection from 95 - 175 degrees...
		maxCross = -1, maxExtent = -1;
		std::vector<int> hist1;
		int flag = 0;
		do {
			intercept1++;
			for (int y = 0; y < h; y++) {
				int x = int((y - intercept1) / slope1);
				if (x >= w || x < 0) continue;
				if (img[y][x] == 0) {
					flag = 1;
					break;
				}
			}
			//std::cout << "\nintercept = " << intercept1;
		}while (flag == 0);
		do {
			//std::cout << "\nintercept = " << intercept1;
			count = 0, crossCount = 0;
			int prev = 0, curr;
			Point first, last;
			for (int y = 0; y < h; y++) {
				int x = int((y - intercept1) / slope1);
				if (x >= w || x < 0) continue;
				count += 1 - img[y][x];
				int nbrs[Nb][Nb];
				for (int k = y - 3; k <= y + 3; k++) {
					for (int l = x - 3; l <= x + 3; l++) {
						if (k<0 || k >= h || l<0 || l >= w)
							nbrs[k - y + 3][l - x + 3] = 0;
						else
							nbrs[k - y + 3][l - x + 3] = img[k][l];
					}
				}
				curr = isInteriorPoint(nbrs);
				if (img[y][x] == 0 && crossCount == 0)
					first = Point(x, y);

				if (curr == 1 && prev == 0) { //white->black 
					crossCount++;
				}

				if (img[y][x] == 0)
					last = Point(x, y);
				prev = curr;
			}
			if (count == 0) 
				break;
			hist1.push_back(count);
			if (count < min) {
				min = count;
			}
			if (count > max) {
				max = count;
			}
			extent = norm(last - first);
			if (extent > maxExtent)
				maxExtent = extent;
			if (crossCount > maxCross)
				maxCross = crossCount;
			intercept1++;
		}while (count > 0);
		dev = stddev(hist1, int(hist1.size()));

		features[K * 2 + 1][1] = dev;
		features[K * 2 + 1][2] = maxExtent;
		features[K * 2 + 1][3] = maxCross;

		//projection from 5 - 85 degrees...
		intercept2 = h - 1, maxCross = -1, maxExtent = -1;
		count = 0, min = 9999, max = -1;
		std::vector<int> hist2;
		flag = 0;
		do {
			intercept2--;
			for (int y = h - 1; y >= 0; y--) {
				int x = (y - intercept2) / slope2;
				if (x >= w || x < 0) continue;
				if (img[y][x] == 0) {
					flag = 1;
					break;
				}
			}
		} while (flag == 0);
		do {
			count = 0, crossCount = 0;
			int prev = 0, curr;
			Point first, last;
			for (int y = h - 1; y >= 0; y--) {
				int x = (y - intercept2) / slope2;
				if (x >= w || x < 0) continue;
				count += 1 - img[y][x];
				int nbrs[Nb][Nb];
				for (int k = y - 3; k <= y + 3; k++) {
					for (int l = x - 3; l <= x + 3; l++) {
						if (k<0 || k >= h || l<0 || l >= w)
							nbrs[k - y + 3][l - x + 3] = 0;
						else
							nbrs[k - y + 3][l - x + 3] = img[k][l];
					}
				}
				curr = isInteriorPoint(nbrs);
				if (img[y][x] == 0 && crossCount == 0)
					first = Point(x, y);

				if (curr == 1 && prev == 0) { //white->black 
					crossCount++;
				}
				if (img[y][x] == 0)
					last = Point(x, y);
				prev = curr;
			}
			if (count == 0)
				break;
			hist2.push_back(count);
			if (count < min) {
				min = count;
			}
			if (count > max) {
				max = count;
			}
			extent = norm(last - first);
			if (extent > maxExtent)
				maxExtent = extent;
			if (crossCount > maxCross)
				maxCross = crossCount;
			intercept2--;
		} while (count > 0);
		dev = stddev(hist2, int(hist2.size()));

		features[K * 2][1] = dev;
		features[K * 2][2] = maxExtent;
		features[K * 2][3] = maxCross;
	}

	double maxFeature1 = -1, maxFeature2 = -1, maxFeature3 = -1;
	std::cout << "\nFeatures List (angle/std dev/max/bw-cross):";
	for (int i = 0; i < 36; i++) {
		std::cout << "\n" << features[i][0] << "\t" << features[i][1] << "      " << features[i][2] << "\t" << features[i][3];
		if (maxFeature1 < features[i][1]) {
			maxFeature1 = features[i][1];
			inclineStdDev = features[i][0];
		}
		if (maxFeature2 < features[i][2]) {
			maxFeature2 = features[i][2];
			inclineExtent = features[i][0];
		}
		if (maxFeature3 < features[i][3]) {
			maxFeature3 = features[i][3];
			inclineCrossing = features[i][0];
		}
	}
	std::cout << "\nInclined at " << inclineStdDev << " , " << inclineExtent << " , " << inclineCrossing;
}

//collect pixels to be deleted in horizontal direction...
void trim0(int **roi, int h, int w, std::vector<Point>& deletepixels) {
	int kstart = 0, kend;
	while (kstart < h) {
		int sum = 0;
		kend = kstart + 2;
		if (kend >= h) kend = h - 1;
		for (int stripi = kstart; stripi <= kend; stripi++)
			for (int stripj = 0; stripj < w; stripj++)
				if (roi[stripi][stripj] == 0)
					sum++;

		if (double(sum) / (3 * w) < 0.2) {
			for (int stripi = kstart; stripi <= kend; stripi++)
				for (int stripj = 0; stripj < w; stripj++)
					if (roi[stripi][stripj] == 0) {     //collect pixel locations to be deleted....
						deletepixels.push_back(Point(stripj, stripi));
					}
		}
		kstart = kend + 1;
	}
	//std::cout << "\nNo. of pixels to delete = " << int(deletepixels.size());
}

//collect pixels to be deleted in vertical direction...
void trim90(int **roi, int h, int w, std::vector<Point>& deletepixels) {
	int kstart = 0, kend;
	while (kstart < w) {
		int sum = 0;
		kend = kstart + 2;
		if (kend >= w) kend = w - 1;
		for (int stripj = kstart; stripj <= kend; stripj++)
			for (int stripi = 0; stripi < h; stripi++)
				if (roi[stripi][stripj] == 0)
					sum++;

		if (double(sum) / (3 * h) < 0.2) {
			for (int stripj = kstart; stripj <= kend; stripj++)
				for (int stripi = 0; stripi < h; stripi++)
					if (roi[stripi][stripj] == 0) {     //collect pixel locations to be deleted....
						deletepixels.push_back(Point(stripj,stripi));
					}
		}
		kstart = kend + 1;
	}
	//std::cout << "\nNo. of pixels to delete = " << int(deletepixels.size());
}

//collect pixels to be deleted in off-diagonal direction...
void trim45(int **img, int **label, double *corners, int lbl, std::vector<Point>& deletepixels) {
	Point N = Point(int(corners[1]), int(corners[0]));
	Point S = Point(int(corners[3]), int(corners[2]));
	Point E = Point(int(corners[5]), int(corners[4]));
	Point W = Point(int(corners[7]), int(corners[6]));

	int sz = 2 * (W.y - N.y + 1);
	int C = W.y - W.x; //c = W.y - W.x...
	int ilimit = (S.y - W.y + 1);
	int istart = W.y, jstart = W.x, adder = 1;

	int kstart = 0, kend;
	while (kstart < sz) {
		int sum = 0;
		kend = kstart + 2;
		if (kend >= sz) kend = sz - 1;
		int copy_istart = istart, copy_jstart = jstart;
		for (int k = kstart; k <= kend; k++) {
			for (int i = istart; i <= istart + ilimit; i++) {
				int j = i - (C - k); //y - x = c...
				if (i < 0 || i >= r || j < 0 || j >= c) {
					continue;
				}
				if (label[i][j] == lbl)
					sum++;
			}
			adder++;
			if (adder % 2 == 0) //even run...
				istart--;
			else      //odd run...
				jstart++;
		}
		if (double(sum) / (3 * ilimit) < 0.2) {
			for (int k = kstart; k <= kend; k++) {
				for (int i = copy_istart; i <= copy_istart + ilimit; i++) {
					int j = i - (C - k); //y - x = c...
					if (i < 0 || i >= r || j < 0 || j >= c) {
						continue;
					}
					if (label[i][j] == lbl) {     //collect pixel locations to be deleted....
						deletepixels.push_back(Point(j,i));
					}
				}
				adder++;
				if (adder % 2 == 0) //even run...
					copy_istart--;
				else      //odd run...
					copy_istart++;
			}
		}
		kstart = kend + 1;
	}
	//std::cout << "\nNo. of pixels to delete = " << int(deletepixels.size());
}

//collect pixels to be deleted in diagonal direction...
void trim135(int **img, int **label, double *corners, int lbl, std::vector<Point>& deletepixels) {
	Point N = Point(int(corners[1]), int(corners[0]));
	Point S = Point(int(corners[3]), int(corners[2]));
	Point E = Point(int(corners[5]), int(corners[4]));
	Point W = Point(int(corners[7]), int(corners[6]));

	int sz = 2 * (S.y - W.y + 1);
	int C = N.x + N.y; //c = N.x + N.y...
	int ilimit = (W.y - N.y + 1);
	int istart = W.y, jstart = W.x, adder = 1;

	int kstart = 0, kend;
	while (kstart < sz) {
		int sum = 0;
		kend = kstart + 2;
		if (kend >= sz) kend = sz - 1;
		int copy_istart = istart, copy_jstart = jstart;
		for (int k = kstart; k <= kend; k++) {
			for (int i = istart; i >= istart - ilimit; i--) {
				int j = (C + k) - i; //x + y = c...
				if (i < 0 || i >= r || j < 0 || j >= c) {
					continue;
				}
				if (label[i][j] == lbl)
					sum++;
			}
			adder++;
			if (adder % 2 == 0) //even run...
				istart++;
			else      //odd run...
				jstart++;
		}
		if (double(sum) / (3 * ilimit) < 0.2) {
			for (int k = kstart; k <= kend; k++) {
				for (int i = istart; i >= istart - ilimit; i--) {
					int j = (C + k) - i; //x + y = c...
					if (i < 0 || i >= r || j < 0 || j >= c) {
						continue;
					}
					if (label[i][j] == lbl) {     //collect pixel locations to be deleted....
						deletepixels.push_back(Point(j,i));
					}
				}
				adder++;
				if (adder % 2 == 0) //even run...
					copy_istart++;
				else      //odd run...
					copy_istart--;
			}
		}
		kstart = kend + 1;
	}
	//std::cout << "\nNo. of pixels to delete = " << int(deletepixels.size());
}

//Trim each CC in suitable direction...
void trimCC(int **img, int **label, int **box, double **boundingBox) {
	for (int i = 0; i < labelNum; i++) {
		if (box[i][2] == -1)
			continue;
		int h = box[i][5] - box[i][4] + 1;
		int w = box[i][3] - box[i][2] + 1;

		int **roi, **labelCC, i_start = box[i][4], j_start = box[i][2], total = 0;
		roi = new int*[h];
		labelCC = new int*[h];
		for (int k = 0; k < h; k++) {
			roi[k] = new int[w];
			labelCC[k] = new int[w];
			for (int l = 0; l < w; l++) {
				if (label[i_start + k][j_start + l] == i){
					roi[k][l] = img[i_start + k][j_start + l];
					total++;
				}
				else
					roi[k][l] = 255;
			}
		}
		std::cout << "\nLabel " << i << " : " << "\nTotal no. of pixels = " << total;
		std::vector<Point> trimhor, trimver, trimdiag, trimoffdiag;
		trim0(roi, h, w, trimhor);
		trim90(roi, h, w, trimver);
		trim45(img, label, boundingBox[i], i, trimoffdiag);
		trim135(img, label, boundingBox[i], i, trimdiag);

		std::cout << "\nNo. of pixels to delete = " << int(trimhor.size());
		std::cout << "\nNo. of pixels to delete = " << int(trimver.size());
		std::cout << "\nNo. of pixels to delete = " << int(trimdiag.size());
		std::cout << "\nNo. of pixels to delete = " << int(trimoffdiag.size());
		int totalH = total - int(trimhor.size());
		int totalV = total - int(trimver.size());
		int totalD = total - int(trimdiag.size());
		int totalOD = total - int(trimoffdiag.size());

		int **subccH, **subccV, **subccOD, **subccD, count = 0;
		subccH = roi;
		subccV = roi;
		subccOD = roi;
		subccD = roi;

		subccH = new int*[h];
		subccV = new int*[h];
		subccD = new int*[h];
		subccOD = new int*[h];
		for (int k = 0; k < h; k++) {
			subccH[k] = new int[w];
			subccV[k] = new int[w];
			subccD[k] = new int[w];
			subccOD[k] = new int[w];
			for (int l = 0; l < w; l++) {
				subccH[k][l] = roi[k][l];
				subccV[k][l] = roi[k][l];
				subccD[k][l] = roi[k][l];
				subccOD[k][l] = roi[k][l];
			}
		}

		//horizontal trimming...
		for (int k = 0; k < int(trimhor.size()); k++) {
			subccH[trimhor[k].y][trimhor[k].x] = 255;
		}
		for (int k = 0; k < h; k++) {
			for (int l = 0; l < w; l++)
				labelCC[k][l] = 0;
		}
		int lblnoH = labelling(subccH, labelCC, h, w);
		int *lblcount = new int[lblnoH];
		for (int k = 0; k < lblnoH; k++) {
			lblcount[k] = 0;
		}
		for (int k = 0; k < h; k++) {
			for (int l = 0; l < w; l++)
				lblcount[labelCC[k][l]]++;
		}
		for (int k = 1; k < lblnoH; k++) {
			if (lblcount[k] != 0) count++;
		}
		lblnoH = count;
		delete[] lblcount;

		//--------------------------------------------------------------------------------------------
		//Vertical trimmming...
		count = 0;
		for (int k = 0; k < int(trimver.size()); k++) {
			subccV[trimver[k].y][trimver[k].x] = 255;
		}
		for (int k = 0; k < h; k++) {
			for (int l = 0; l < w; l++)
				labelCC[k][l] = 0;
		}
		int lblnoV = labelling(subccV, labelCC, h, w);
		lblcount = new int[lblnoV];
		for (int k = 0; k < lblnoV; k++) {
			lblcount[k] = 0;
		}
		for (int k = 0; k < h; k++) {
			for (int l = 0; l < w; l++)
				lblcount[labelCC[k][l]]++;
		}
		for (int k = 1; k < lblnoV; k++) {
			if (lblcount[k] != 0) count++;
		}
		lblnoV = count;
		delete[] lblcount;

		//--------------------------------------------------------------------------------------------
		//Diagonal trimming...
		count = 0;
		for (int k = 0; k < int(trimdiag.size()); k++) {
			subccD[trimdiag[k].y - i_start][trimdiag[k].x - j_start] = 255;
		}
		for (int k = 0; k < h; k++) {
			for (int l = 0; l < w; l++)
				labelCC[k][l] = 0;
		}
		int lblnoD = labelling(subccD, labelCC, h, w);
		lblcount = new int[lblnoD];
		for (int k = 0; k < lblnoD; k++) {
			lblcount[k] = 0;
		}
		for (int k = 0; k < h; k++) {
			for (int l = 0; l < w; l++)
				lblcount[labelCC[k][l]]++;
		}
		for (int k = 1; k < lblnoD; k++) {
			if (lblcount[k] != 0) count++;
		}
		lblnoD = count;
		delete[] lblcount;

		//--------------------------------------------------------------------------------------------
		//Off-Diagonal trimming...
		count = 0;
		for (int k = 0; k < int(trimoffdiag.size()); k++) {
			subccOD[trimoffdiag[k].y - i_start][trimoffdiag[k].x - j_start] = 255;
		}
		for (int k = 0; k < h; k++) {
			for (int l = 0; l < w; l++)
				labelCC[k][l] = 0;
		}
		int lblnoOD = labelling(subccOD, labelCC, h, w);
		lblcount = new int[lblnoOD];
		for (int k = 0; k < lblnoOD; k++) {
			lblcount[k] = 0;
		}
		for (int k = 0; k < h; k++) {
			for (int l = 0; l < w; l++)
				lblcount[labelCC[k][l]]++;
		}
		for (int k = 1; k < lblnoOD; k++) {
			if (lblcount[k] != 0) count++;
		}
		lblnoD = count;
		delete[] lblcount;


		std::cout << "\nNo.of sub-CCs : " << lblnoH << "(0) , " << lblnoV << "(90) , " << lblnoOD << "(45) , " << lblnoD << "(135)";
		std::cout << "\nNo.of pixels in sub-CCs : " << totalH << "(0) , " << totalV << "(90) , " << totalOD << "(45) , " << totalD << "(135)";

		int trimdata[4][3];
		trimdata[0][0] = 0;
		trimdata[1][0] = 90;
		trimdata[2][0] = 135;
		trimdata[3][0] = 45;
		trimdata[0][1] = lblnoH ;
		trimdata[1][1] = lblnoV;
		trimdata[2][1] = lblnoD;
		trimdata[3][1] = lblnoOD;
		trimdata[0][2] = totalH;
		trimdata[1][2] = totalV;
		trimdata[2][2] = totalD;
		trimdata[3][2] = totalOD;

		int dir = trimdata[0][0];
		int mindir = trimdata[0][1];
		int minpixel = trimdata[0][2];
		for (int k = 1; k < 4; k++) {
			if (trimdata[k][1] < mindir) {
				dir = trimdata[k][0];
				mindir = trimdata[k][1];
				minpixel = trimdata[k][2];
			}
			else if (trimdata[k][1] == mindir) {
				if (trimdata[k][2] > minpixel) {
					dir = trimdata[k][0];
					mindir = trimdata[k][1];
					minpixel = trimdata[k][2];
				}
				else if (trimdata[k][2] == minpixel) {
					dir = -1;
				}
			}
		}
		std::cout << "\nSelected direction = " << dir;
		if (dir == 0 && double(minpixel)/total > 0.65) {
			for (int k = 0; k < h; k++) {
				for (int l = 0; l < w; l++) {
					if (subccH[k][l] == 255 && label[i_start + k][j_start + l] == i) {
						img[i_start + k][j_start + l] = 255;
						label[i_start + k][j_start + l] = 0;
					}
				}
			}
		}

		if (dir == 90 && double(minpixel) / total > 0.65) {
			for (int k = 0; k < h; k++) {
				for (int l = 0; l < w; l++) {
					if (subccV[k][l] == 255 && label[i_start + k][j_start + l] == i) {
						img[i_start + k][j_start + l] = 255;
						label[i_start + k][j_start + l] = 0;
					}
				}
			}
		}

		if (dir == 45 && double(minpixel) / total > 0.65) {
			for (int k = 0; k < h; k++) {
				for (int l = 0; l < w; l++) {
					if (subccOD[k][l] == 255 && label[i_start + k][j_start + l] == i) {
						img[i_start + k][j_start + l] = 255;
						label[i_start + k][j_start + l] = 0;
					}
				}
			}
		}

		if (dir == 135 && double(minpixel) / total > 0.65) {
			for (int k = 0; k < h; k++) {
				for (int l = 0; l < w; l++) {
					if (subccD[k][l] == 255 && label[i_start + k][j_start + l] == i) {
						img[i_start + k][j_start + l] = 255;
						label[i_start + k][j_start + l] = 0;
					}
				}
			}
		}

		//delete...
		for (int k = 0; k < h; k++) {
			delete[] roi[k];
			delete[] labelCC[k];
			delete[] subccD[k];
			delete[] subccH[k];
			delete[] subccOD[k];
			delete[] subccV[k];
		}
		delete[] roi;
		delete[] labelCC;
		delete[] subccD;
		delete[] subccH;
		delete[] subccOD;
		delete[] subccV;
	}
}

