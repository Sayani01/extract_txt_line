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
//#include "Projections.h"
//#include "OIC.h"
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
//void sort(std::vector<int> arr, int n) {
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
//
//double median(std::vector<int>  arr, int sz) {
//	//sort(arr, sz);
//	/*for (int i = 0; i < sz - 2; i++)
//		std::cout << "\n" << arr[i];*/
//	if (sz % 2 == 1)
//		return(arr[int((sz - 1) / 2)]);
//	else
//		return((arr[int(sz / 2)] + arr[int(sz / 2) - 1]) / 2.0);
//}
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
//
//double mean(std::vector<int> arr, int sz) {
//	double M = 0.0;
//	for (int i = 0; i < sz; i++) {
//		//std::cout << "\n" << arr[i]<<" M = "<<M;
//		M += (double)arr[i];
//	}
//	//std::cout << "\n mean = " << M ;
//	return(M / (double)sz);
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
//	int struc[3][3] = { { 1,1,1 },{ 1,1,1 },{ 1,1,1 } };
//	for (int i = (n - 1) / 2; i < h - (n - 1) / 2; i++) {
//		for (int j = (n - 1) / 2; j < w - (n - 1) / 2; j++) {
//			if (copy[i][j] == 0) continue; //text pixel...
//			int sum = 0;
//			for (int k = -(n - 1) / 2; k <= (n - 1) / 2; k++) {
//				for (int l = -(n - 1) / 2; l <= (n - 1) / 2; l++) {
//					if (copy[i + k][j + l] == 0 && struc[(n - 1) / 2 + k][(n - 1) / 2 + l] == 1) sum++; //black nbr...
//				}
//			}
//			if (sum != 0) img[i][j] = 0;
//		}
//	}
//}
//
//void dilate(int **img, int h = r, int w = c) {
//	int n = 3;
//	int **copy = new int*[r];
//	for (int i = 0; i < r; i++) {
//		copy[i] = new int[c];
//		for (int j = 0; j < c; j++)
//			copy[i][j] = img[i][j] / 255;
//	}
//	int struc[3][3] = { { 1,1,1 },{ 1,1,1 },{ 1,1,1 } };
//	for (int i = 1; i < h - 1; i++) {
//		for (int j = 1; j < w - 1; j++) {
//			if (copy[i][j] == 1) continue;
//			int sum = 0;
//			for (int k = 0; k < 3; k++) {
//				for (int l = 0; l < 3; l++) {
//					sum += copy[i - 1 + k][j - 1 + l] * struc[k][l];
//				}
//			}
//			if (sum != 0) img[i][j] = 255;
//		}
//	}
//}
//
//
//void whiteLines(int *hist, int size, int *Y) {
//	int start = 0, stop = size - 1;
//	for (int i = 0; i < size; i++) {
//		if (hist[i] > 0) {
//			start = i;
//			break;
//		}
//	}
//	for (int i = size - 1; i >= 0; i--) {
//		if (hist[i] > 0) {
//			stop = i;
//			break;
//		}
//	}
//	Y[start] = 1;
//	Y[stop] = 1;
//	for (int i = start + 1; i < stop; i++) {
//		int count = 0;
//		if (hist[i] >= 5) continue;
//		for (int k = i; k <= stop; k++) {
//			if (hist[k] < 5)
//				count++;
//			else
//				break;
//		}
//		if (count >= 10) {
//			int k = i + (count) / 2;
//			Y[k] = 1;
//		}
//		i = i + count - 1;
//	}
//}
//
///***********************************************************************************************************
//										PROJECTION PROFILES
//************************************************************************************************************/
//
////Histogram clusters...
//std::vector<std::vector<int>> ClustersOffDiag(int *hist, int sz) {
//	std::vector<std::vector<int>> peaks;
//	for (int i = 2; i < sz; i++) {
//		if (hist[i - 1] < 5 && hist[i] > 0) { //start a cluster...
//			std::vector<int> cluster;
//			cluster.push_back(i);
//			int max = hist[i];
//			i++;
//			do {
//				if (max < hist[i])
//					max = hist[i];
//				i++;
//			} while (hist[i] > 5 && i < sz);
//			cluster.push_back(i - 1);
//			cluster.push_back(max);
//			peaks.push_back(cluster);
//		}
//	}
//	/*for (int i = 1; i<int(peaks.size()); i++) {
//		if (abs(peaks[i - 1][1] - peaks[i][0]) < 30) {
//			peaks[i - 1][1] = peaks[i][1];
//			peaks[i - 1][2] = peaks[i - 1][2] > peaks[i][2] ? peaks[i - 1][2] : peaks[i][2];
//			peaks.erase(peaks.begin() + i);
//			i--;
//		}
//	}*/
//	return(peaks);
//}
//
////void smoothingHistogram(int *hist, int sz, int p1 = 10, double p2 = 0.75, double p3 = 0.6) {
////	int win_sz = p1;
////	int avg_loc_max = 0, no_loc_max = 0;
////	for (int i = 0; i < sz ; i += win_sz / 2) {
////		int loc_max = 0;
////		for (int j = 0; j < win_sz; j++) {
////			if (i + j >= sz)
////				break;
////			if (hist[i + j] > loc_max)
////				loc_max = hist[i + j];
////		}
////		no_loc_max++;
////		avg_loc_max += loc_max;
////		for (int j = 0; j < win_sz; j++) {
////			if (i + j >= sz)
////				break;
////			hist[i + j] = int(p2*loc_max);
////		}
////		
////	}
////	avg_loc_max /= no_loc_max;
////	for (int i = 0; i < sz; i++) {
////		if (hist[i] < p3*avg_loc_max) {
////			hist[i] = 0;
////		}
////	}
////}
//
//void smoothingHistogram(int *hist, int sz, int p1, double p2) {
//	int win_sz = 15;
//	std::cout << "\nsize of hist=" << sz;
//	//int sz1 = 2 * int(sz / double(win_sz)) - 1;
//	//int **hist1 = new int*[sz1];
//	/*for (int i = 0; i < sz1; i++) {
//		hist1[i] = new int[2];
//		for (int j = 0; j < 2; j++)
//			hist1[i][j] = 0;
//	}
//	int k = 0;*/
//	for (int i = 0; i < sz; i += win_sz / 2) {
//		int avg_hist = 0;
//		for (int j = 0; j < win_sz; j++) {
//			if (i + j >= sz)
//				break;
//			avg_hist += hist[i + j];
//		}
//		avg_hist /= win_sz;
//		/*hist1[k][0] = i + int(win_sz / 2) - 1;
//		hist1[k][1] = avg_hist;*/
//		for (int j = 0; j < win_sz; j++) {
//			if (i + j >= sz)
//				break;
//			hist[i + j] = avg_hist;
//		}
//		//std::cout << "\n" << hist[i];
//		/*std::cout << "\n" << hist1[k][0] << "    " << hist1[k][1];
//		k++;*/
//	}
//}
//
//std::vector<std::vector<int>> peaksAndValleys(int *hist, int sz) {
//	int start, stop;
//	for (int i = 0; i < sz; i++) {
//		if (hist[i] > 0) {
//			start = i - 1;
//			break;
//		}
//	}
//	for (int i = sz - 1; i >= 0; i--) {
//		if (hist[i] > 0) {
//			stop = i + 1;
//			break;
//		}
//	}
//
//	int max = hist[0];
//	int *seq = new int[sz]; //1 => possible peak , -1 => possible valley
//	for (int i = 1; i < sz; i++) {
//		max = hist[i] > max ? hist[i] : max;
//	}
//	/*double avg = 0;
//	for (int i = start; i <= stop; i++) {
//		avg += hist[i];
//	}
//	avg /= (stop - start + 1);*/
//
//	double thresh = 0.3*max;
//
//	for (int i = 0; i < sz; i++) {
//		seq[i] = hist[i]>=thresh ? 1 : -1;
//	}
//	
//	
//	std::vector<std::vector<int>> pv;
//	for (int i = start; i <= stop;) {
//		if (seq[i] == 1) {
//			int peak_val = 0;
//			int avg = 0, count = 0;
//			while (seq[i] != -1 && i <= stop) {
//				peak_val = hist[i] > peak_val ? hist[i] : peak_val;
//				avg += i;
//				count++;
//				i++;
//			}
//			avg /= count;
//			/*if ((avg - thresh) / (max - thresh) < 0.1)
//				continue;*/
//			std::vector<int> peak;
//			peak.push_back(1);
//			peak.push_back(avg);
//			peak.push_back(peak_val);
//			pv.push_back(peak);
//		}
//		else if (seq[i] == -1) {
//			int valley_val = 99999;
//			int avg = 0, count = 0;
//			while (seq[i] != 1 && i <= stop + 1) {
//				valley_val = hist[i] < valley_val ? hist[i] : valley_val;
//				avg += i;
//				count++;
//				i++;
//			}
//			avg /= count;
//			/*if ((thresh - avg) / double(avg) < 0.2)
//				continue;*/
//			std::vector<int> valley;
//			valley.push_back(-1);
//			valley.push_back(avg);
//			valley.push_back(valley_val);
//			pv.push_back(valley);
//		}
//	}
//
//	std::cout << "\nSequence of peaks and valleys:";
//	for (int i = 0; i<int(size(pv)); i++) {
//		if (pv[i][0] == 1) {
//			std::cout << "\n Peak at i = " << pv[i][1] << " , height = " << pv[i][2];
//		}
//		else if (pv[i][0] == -1) {
//			std::cout << "\n Valley at i = " << pv[i][1] << " , height = " << pv[i][2];
//		}
//	}
//
//	return(pv);
//}
//
////Project horizontally for vertical strips...
//void StripVer(int **img, Mat binimg, int N) {
//	int **reg;
//	Mat colorBinary, regions, reg_mask;
//
//	reg = new int*[r];
//	for (int i = 0; i < r; i++) {
//		reg[i] = new int[c];
//		for (int j = 0; j < c; j++)
//			reg[i][j] = 255;
//	}
//
//	int firstcol = -1, lastcol = -1;
//	for (int j = 0; j < c; j++) {
//		for (int i = 0; i < r; i++) {
//			if (img[i][j] == 0) {
//				firstcol = j;
//				break;
//			}
//		}
//		if (firstcol >= 0) break;
//	}
//	for (int j = c-1; j >= 0; j--) {
//		for (int i = 0; i < r; i++) {
//			if (img[i][j] == 0) {
//				lastcol = j;
//				break;
//			}
//		}
//		if (lastcol >= 0) break;
//	}
//
//	regions.create(r, c, CV_8UC1);
//	cvtColor(binimg, regions, CV_GRAY2BGR);
//
//	colorBinary.create(r, c, CV_8UC1);
//	cvtColor(binimg, colorBinary, CV_GRAY2BGR);
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (img[i][j] == 0)
//				colorBinary.at<Vec3b>(i, j) = Vec3b(150, 150, 150);
//		}
//	}
//	//line(colorBinary, Point(firstcol, 0), Point(firstcol, r-1), Scalar(0, 0, 255), 1, 8);
//
//	//int N = (lastcol - firstcol + 1) / 10;
//	std::vector<int> bar_seq;
//	int count = 0;
//	for (int k = firstcol ; k < c; k += int(N/2)) {
//		if (k > lastcol)
//			break;
//		line(colorBinary, Point(k, 0), Point(k, r - 1), Scalar(0, 0, 255), 1, 8);
//		int **roi, *hist;
//		int w = N;
//		if (k + N >= c)
//			w = c - k;
//		roi = new int*[r];
//		hist = new int[r];
//		for (int i = 0; i < r; i++) {
//			roi[i] = new int[w];
//			hist[i] = 0;
//			for (int j = 0; j < w; j++) {
//				roi[i][j] = img[i][j + k];
//			}
//		}
//		project0(roi, r, w, hist);
//		/*int *Y = new int[r];;
//		for (int i = 0; i < r; i++)
//			Y[i] = 0;
//		whiteLines(hist, r, Y);
//*/
//		std::cout << "\nStrip " << count << " histogram->";
//		smoothingHistogram(hist, r, 10, 0.2);
//		std::vector<std::vector<int>> pv = peaksAndValleys(hist, r);
//
//		std::vector<std::vector<int>> peaks = ClustersOffDiag(hist, r);
//
//		//mark histogram for peak clusters...
//		int *mark = new int[r];
//		for (int t = 0; t < r; t++)
//			mark[t] = 0;
//		for (int t = 0; t < int(peaks.size()); t++) {
//			if (peaks[t][2] > N && peaks[t][2] > 2.5*abs(peaks[t][1]-peaks[t][0])) {
//				//mark[peaks[t][0]] = 1;
//				//mark[peaks[t][1]] = 1;
//				//if (count % 2 == 0) {
//				//	//line(regions, Point(k, peaks[t][0]), Point(k + w, peaks[t][0]), Scalar(0, 0, 255), 4, 8);
//				//	line(regions, Point(k, peaks[t][1]), Point(k + w, peaks[t][1]), Scalar(0, 0, 255), 8, 8);
//				//	//line(regions, Point(k, peaks[t][0]), Point(k , peaks[t][1]), Scalar(0, 0, 255), 4, 8);
//				//	//line(regions, Point(k + w, peaks[t][0]), Point(k + w, peaks[t][1]), Scalar(0, 0, 255), 4, 8);
//				//}
//				//else {
//				//	//line(regions, Point(k, peaks[t][0]), Point(k + w, peaks[t][0]), Scalar(255, 0, 0), 4, 8);
//				//	line(regions, Point(k, peaks[t][1]), Point(k + w, peaks[t][1]), Scalar(255, 0, 0), 4, 8);
//				//	//line(regions, Point(k, peaks[t][0]), Point(k, peaks[t][1]), Scalar(255, 0, 0), 4, 8);
//				//	//line(regions, Point(k + w, peaks[t][0]), Point(k + w, peaks[t][1]), Scalar(255, 0, 0), 4, 8);
//				//}
//				/*for (int y = peaks[t][0]; y <= peaks[t][1]; y++) {
//					for (int x = k; x <= k + w; x++) {
//						reg[y][x] = 0;
//					}
//				}*/
//			}
//		}
//
//		
//
//		for (int t = 0; t < int(pv.size()); t++) {
//			if (pv[t][0] == 1) {
//				int vleft[2], p[2], vright[2];
//				p[0] = pv[t][1];
//				p[1] = pv[t][2];
//				vleft[0] = pv[t - 1][1];
//				vleft[1] = pv[t - 1][2];
//				vright[0] = pv[t + 1][1];
//				vright[1] = pv[t + 1][2];
//				int thickness = (p[0] - vleft[0]) < (vright[0] - p[0]) ? int((p[0] - vleft[0]) / 2.0) : int((vright[0] - p[0]) / 2.0);
//				mark[pv[t][1]] = 1;
//				if ((vleft[1]) / double(p[1]) < 0.1 && (vright[1]) / double(p[1]) < 0.1 && p[1] > 10.0*double(2*thickness)) {
//					line(regions, Point(k, pv[t][1]), Point(k + w, pv[t][1]), Scalar(0, 0, 255), 8, 8);
//					bar_seq.push_back(p[0] - thickness);
//					bar_seq.push_back(p[0] + thickness);
//					for (int y = p[0] - thickness; y <= p[0] + thickness; y++) {
//						for (int x = k; x <= k + w; x++) {
//							reg[y][x] = 0;
//						}
//					}
//				}
//
//			}
//			/*else if (pv[t][0] == -1) {
//				mark[pv[t][1]] = 1;
//				line(regions, Point(k, pv[t][1]), Point(k + w, pv[t][1]), Scalar(255, 0, 0), 8, 8);
//
//			}*/
//			else
//				continue;
//		}
//
//		Vec3b color1, color2;
//		if (count % 2 == 0) {
//			color1 = { 2, 34, 228 };
//			color2 = { 4, 28, 146 };
//		}
//		else {
//			color1 = { 228, 34, 2 };
//			color2 = { 146, 28, 4 };
//		}
//		count++;
//
//		for (int i = 0; i < r; i++) {
//			for (int j = 0; j <= hist[i]; j++) {
//				if (j + k >= c) break;
//				if (img[i][j + k] == 0)
//					colorBinary.at<Vec3b>(i, k + j) = color2;
//				else
//					colorBinary.at<Vec3b>(i, k + j) = color1;
//			}
//			
//			/*if(Y[i] == 1)
//				line(colorBinary, Point(k, i), Point(k + N, i), Scalar(255, 0, 0), 2.5, 8);*/
//			//line(colorBinary, Point(k, i), Point(k + hist[i], i), Scalar(0, 120, 255), 1, 8);
//		}
//
//		std::cout << "\nSequence of masks in strip " << count;
//		for (int i = 0; i < int(bar_seq.size()) - 1; i += 2) {
//			std::cout << "\nFrom i = " << bar_seq[i] << " to i = " << bar_seq[i + 1];
//			if (i + 2 < int(bar_seq.size()) && bar_seq[i + 1] > bar_seq[i + 2])
//				std::cout << "...overlapping!";
//		}
//
//		for (int i = 0; i < r; i++)
//			delete[] roi[i];
//		delete[] roi;
//		delete[] hist;
//		bar_seq.clear();
//	}
//	line(colorBinary, Point(lastcol, 0), Point(lastcol, r - 1), Scalar(0, 0, 255), 2.5, 8);
//	imshow("stripLines.tif", colorBinary);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/VerticalStripsProj-50.tif", colorBinary);
//
//	reg_mask.create(r, c, CV_8UC1);
//	erode(reg);
//	arr2mat(reg, reg_mask);
//	cvtColor(reg_mask, reg_mask, CV_GRAY2BGR);
//
//	imshow("regions.tif", regions);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/VerticalStripsRegions-50.tif", regions);
//
//	imshow("Mask_regions.tif", reg_mask);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/VerticalStripsRegionsMask-50.tif", reg_mask);
//	//imwrite("Data/ICDAR/Data218/VerticalStrips.tif", colorBinary);
//}
//
////Project vertically for horizontal strips...
//void StripHor(int **img, Mat binimg) {
//	Mat colorBinary;
//	int firstrow = -1, lastrow = -1;
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (img[i][j] == 0) {
//				firstrow = i;
//				break;
//			}
//		}
//		if (firstrow >= 0) break;
//	}
//	for (int i = r - 1; i >= 0; i--) {
//		for (int j = 0; j < c; j++) {
//			if (img[i][j] == 0) {
//				lastrow = i;
//				break;
//			}
//		}
//		if (lastrow >= 0) break;
//	}
//	colorBinary.create(r, c, CV_8UC1);
//	cvtColor(binimg, colorBinary, CV_GRAY2BGR);
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (img[i][j] == 0)
//				colorBinary.at<Vec3b>(i, j) = Vec3b(150, 150, 150);
//		}
//	}
//	//line(colorBinary, Point(firstcol, 0), Point(firstcol, r-1), Scalar(0, 0, 255), 1, 8);
//
//	int N = (lastrow - firstrow + 1) / 10;
//	for (int k = firstrow; k < lastrow - N; k += N) {
//		line(colorBinary, Point(0, k), Point(c - 1, k), Scalar(0, 0, 255), 1, 8);
//		int **roi, *hist;
//		int h = N;
//		roi = new int*[h];
//		hist = new int[c];
//		for (int i = 0; i < h; i++) {
//			roi[i] = new int[c];
//			for (int j = 0; j < c; j++) {
//				hist[j] = 0;
//				roi[i][j] = img[i + k][j];
//			}
//		}
//		project90(roi, h, c, hist);
//		for (int j = 0; j < c; j++) {
//			Vec3b color1 = { 2, 34, 228 }, color2 = { 4, 28, 146 };
//			for (int i = 0; i <= hist[j]; i++) {
//				if (i + k >= r) break;
//				if (img[i][j] == 0)
//					colorBinary.at<Vec3b>(k + i, j) = color2;
//				else
//					colorBinary.at<Vec3b>(k + i, j) = color1;
//			}
//			//line(colorBinary, Point(k, i), Point(k + hist[i], i), Scalar(0, 120, 255), 1, 8);
//		}
//
//		for (int i = 0; i < h; i++)
//			delete[] roi[i];
//		delete[] roi;
//		delete[] hist;
//
//	}
//	line(colorBinary, Point(0, lastrow), Point(c - 1, lastrow), Scalar(0, 0, 255), 1, 8);
//	imshow("stripLines.tif", colorBinary);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/HorizontalStripsProjorg.tif", colorBinary);
//	//imwrite("Data/ICDAR/Data218/HorizontalStrips.tif", colorBinary);
//}
//
////Project diagonally for off-diagonal strips...
//void StripDiag(int **img, Mat binimg) {
//	Mat colorBinary;
//	int firstdiag = -1, lastdiag = -1;
//	Point N, S, W, E;
//
//	colorBinary.create(r, c, CV_8UC1);
//	cvtColor(binimg, colorBinary, CV_GRAY2BGR);
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (img[i][j] == 0)
//				colorBinary.at<Vec3b>(i, j) = Vec3b(150, 150, 150);
//		}
//	}
//	for (int C = 0; C <= r + c - 2; C++) {
//		for (int x = C; x >= 0; x--) {
//			int y = C - x;
//			if (y >= r || x >= c) break;
//			if (img[y][x] == 0) {
//				firstdiag = C;
//				break;
//			}
//		}
//		if (firstdiag >= 0)
//			break;
//
//	}
//
//	for (int C = r + c - 2; C >= 0; C--) {
//		for (int x = C - r + 1; x <= c - 1; x++) {
//			int y = C - x;
//			if (y >= r || x < 0) break;
//			if (img[y][x] == 0) {
//				lastdiag = C;
//				break;
//			}
//		}
//		if (lastdiag >= 0)
//			break;
//
//	}
//
//	int no = (lastdiag - firstdiag + 1) / 10;
//	for (int k = firstdiag; k <= lastdiag - no + 1; k += no) {
//		if (k < c) {
//			N = Point(k, 0);
//		}
//		else {
//			N = Point(c - 1, k - c + 1);
//		}
//		if (k < r) {
//			W = Point(0, k);
//		}
//		else {
//			W = Point(k - r + 1, r - 1);
//		}
//		line(colorBinary, N, W, Scalar(0, 0, 255), 1, 8);
//		int k1 = k + no - 1;
//		if (k1 < c) {
//			E = Point(k1, 0);
//		}
//		else {
//			E = Point(c - 1, k1 - c + 1);
//		}
//		if (k1 < r) {
//			S = Point(0, k1);
//		}
//		else {
//			S = Point(k1 - r + 1, r - 1);
//		}
//
//		int *intercepts = new int[4];
//		intercepts[0] = N.y - N.x;
//		intercepts[1] = W.y - W.x;
//		intercepts[2] = S.y - S.x;
//		intercepts[3] = E.y - E.x;
//		int min = intercepts[0], max = intercepts[0];
//		for (int i = 1; i < 4; i++) {
//			if (min > intercepts[i])
//				min = intercepts[i];
//			if (max < intercepts[i])
//				max = intercepts[i];
//		}
//		int sz = max - min + 1;
//		int *hist = new int[sz];
//		for (int i = 0; i < sz; i++)
//			hist[i] = 0;
//		double *dim = new double[8];
//		dim[0] = (k + min) / 2;
//		dim[1] = (k - min) / 2;
//		dim[2] = (k1 + max) / 2;
//		dim[3] = (k1 - max) / 2;
//		dim[4] = (k1 + min) / 2;
//		dim[5] = (k1 - min) / 2;
//		dim[6] = (k + max) / 2;
//		dim[7] = (k - max) / 2;
//		project45(img, img, dim, 0, hist);
//		int *Y = new int[sz];;
//		for (int i = 0; i < sz; i++)
//			Y[i] = 0;
//		whiteLines(hist, sz, Y);
//
//		Vec3b color1 = { 2, 34, 228 }, color2 = { 4, 28, 146 };
//		int adder = 1, ystart = (k + max) / 2;
//
//		for (int t = 0; t < sz; t++) {
//			//std::cout << "\nhere c=" << C - k << "i:" << istart << " - " << istart + ilimit;
//			for (int i = ystart; i <= ystart + hist[t]; i++) {
//				int j = i - (max - t); //y - x = c...
//				if (i < 0 || i >= r || j < 0 || j >= c) {
//					continue;
//				}
//				if (img[i][j] == 0)
//					colorBinary.at<Vec3b>(i, j) = color2;
//				else
//					colorBinary.at<Vec3b>(i, j) = color1;
//			}
//			/*if (Y[t] == 1) {
//				line(colorBinary, Point((k - max + t) / 2, (max - t + k) / 2), Point((k1 - max + t) / 2, (max - t + k1) / 2), Scalar(255, 0, 0), 2.5, 8);
//			}*/
//			adder++;
//			if (adder % 2 == 0) //even run...
//				ystart--;
//		}
//	}
//
//	if (lastdiag < c) {
//		E = Point(lastdiag, 0);
//	}
//	else {
//		E = Point(c - 1, lastdiag - c + 1);
//	}
//	if (lastdiag < r) {
//		S = Point(0, lastdiag);
//	}
//	else {
//		S = Point(lastdiag - r + 1, r - 1);
//	}
//
//	line(colorBinary, E, S, Scalar(0, 0, 255), 1, 8);
//	//int N = (lastdiag - firstdiag + 1) / 10;
//	imshow("stripLines.tif", colorBinary);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/OffDiagStripsProj.tif", colorBinary);
//	//imwrite("Data/ICDAR/Data218/DiagonalStrips.tif", colorBinary);
//}
//
////Project off-diagonally for diagonal strips...
//void StripOffDiag(int **img, Mat binimg) {
//	Mat colorBinary, regions;
//	int firstdiag = -9999, lastdiag = -9999;
//	Point N, S, W, E;
//
//	regions.create(r, c, CV_8UC1);
//	cvtColor(binimg, regions, CV_GRAY2BGR);
//	colorBinary.create(r, c, CV_8UC1);
//	cvtColor(binimg, colorBinary, CV_GRAY2BGR);
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (img[i][j] == 0)
//				colorBinary.at<Vec3b>(i, j) = Vec3b(150, 150, 150);
//		}
//	}
//	for (int C = 1 - c; C <= r - 1; C++) {
//		for (int x = -C; x <= c - 1; x++) {
//			int y = C + x;
//			if (y < 0 || y >= r || x < 0 || x >= c) continue;
//			if (img[y][x] == 0) {
//				firstdiag = C;
//				break;
//			}
//		}
//		if (firstdiag > -9999)
//			break;
//
//	}
//
//	for (int C = r - 1; C >= 1 - c; C--) {
//		for (int x = 0; x <= C + r - 1; x++) {
//			int y = C + x;
//			if (y < 0 || y >= r || x < 0 || x >= c) continue;
//			if (img[y][x] == 0) {
//				lastdiag = C;
//				break;
//			}
//		}
//		if (lastdiag > -9999)
//			break;
//
//	}
//	
//	int no = (lastdiag - firstdiag + 1) / 10;
//	int counter = 0;
//	for (int k = firstdiag; k < lastdiag - no; k += no) {
//		if (k < 0) {
//			N = Point(-k, 0);
//		}
//		else {
//			N = Point(0, k);
//		}
//		if (k < 0) {
//			E = Point(c - 1, c - 1 + k);
//		}
//		else {
//			E = Point(r - 1 - k, r - 1);
//		}
//		line(colorBinary, N, E, Scalar(0, 0, 255), 1, 8);
//		line(regions, N, E, Scalar(0, 0, 255), 1, 8);
//		int k1 = k + no - 1;
//		if (k1 < 0) {
//			W = Point(-k1, 0);
//		}
//		else {
//			W = Point(0, k1);
//		}
//		if (k1 < 0) {
//			S = Point(c - 1, c - 1 + k1);
//		}
//		else {
//			S = Point(r - 1 - k1, r - 1);
//		}
//		int *intercepts = new int[4];
//		intercepts[0] = N.y + N.x;
//		intercepts[1] = W.y + W.x;
//		intercepts[2] = S.y + S.x;
//		intercepts[3] = E.y + E.x;
//		std::cout << "\n";
//		for (int i = 0; i < 4; i++) {
//			std::cout << "intercept[" << i << "] = " << intercepts[i] << "\t";
//		}
//		int min, max;
//		min = intercepts[0] <= intercepts[1] ? intercepts[0] : intercepts[1];
//		max = intercepts[2] >= intercepts[3] ? intercepts[2] : intercepts[3];
//
//		for (int i = 1; i < 4; i++) {
//			if (min > intercepts[i])
//				min = intercepts[i];
//			if (max < intercepts[i])
//				max = intercepts[i];
//		}
//		std::cout << "\nmax = " << max << " , min = " << min;
//		int sz = max - min + 1;
//		int *hist = new int[sz];
//		for (int i = 0; i < sz; i++)
//			hist[i] = 0;
//		double *dim = new double[8];
//		dim[0] = (k + min) / 2;
//		dim[1] = -(k - min) / 2;
//		dim[2] = (k1 + max) / 2;
//		dim[3] = -(k1 - max) / 2;
//		dim[4] = (k + max) / 2;
//		dim[5] = -(k - max) / 2;
//		dim[6] = (k1 + min) / 2;
//		dim[7] = -(k1 - min) / 2;
//		project135(img, img, dim, 0, hist);
//
//		int *Y = new int[sz];;
//		for (int i = 0; i < sz; i++)
//			Y[i] = 0;
//		whiteLines(hist, sz, Y);
//
//		smoothingHistogram(hist, sz, 10, 0.2);
//
//		std::vector<std::vector<int>> peaks = ClustersOffDiag(hist, sz);
//
//		Vec3b color1, color2;
//		if (counter % 2 == 0) {
//			color1 = { 2, 34, 228 };
//			color2 = { 4, 28, 146 };
//		}
//		else {
//			color1 = { 228, 34, 2 };
//			color2 = { 146, 28, 4 };
//		}
//		int adder = 1, ystart = (k1 + min) / 2;
//
//		//mark histogram for peak clusters...
//		int *mark = new int[sz];
//		for (int t = 0; t < sz; t++)
//			mark[t] = 0;
//		for (int t = 0; t < int(peaks.size()); t++) {
//			if ( peaks[t][2] > no) {
//				mark[peaks[t][0]] = 1;
//				mark[peaks[t][1]] = 1;
//				line(regions, Point((min + peaks[t][0] - k1) / 2, (k1 + min + peaks[t][0]) / 2), Point((min + peaks[t][0] - k) / 2, (k + min + peaks[t][0]) / 2), Scalar(255, 0, 0), 4, 8);
//				line(regions, Point((min + peaks[t][1] - k1) / 2, (k1 + min + peaks[t][1]) / 2), Point((min + peaks[t][1] - k) / 2, (k + min + peaks[t][1]) / 2), Scalar(255, 0, 0), 4, 8);
//				line(regions, Point((min + peaks[t][0] - k1) / 2, (k1 + min + peaks[t][0]) / 2), Point((min + peaks[t][1] - k1) / 2, (k1 + min + peaks[t][1]) / 2), Scalar(255, 0, 0), 4, 8);
//				line(regions, Point((min + peaks[t][1] - k) / 2, (k + min + peaks[t][1]) / 2), Point((min + peaks[t][0] - k) / 2, (k + min + peaks[t][0]) / 2), Scalar(255, 0, 0), 4, 8);
//			}
//		}
//		for (int t = 0; t < sz; t++) {
//			//std::cout << "\nhere c=" << C - k << "i:" << istart << " - " << istart + ilimit;
//			for (int i = ystart; i >= ystart - hist[t]; i--) {
//				int j = (min + t) - i; //y + x = c...
//				if (i < 0 || i >= r || j < 0 || j >= c) {
//					continue;
//				}
//				if (img[i][j] == 0)
//					colorBinary.at<Vec3b>(i, j) = color2;
//				else
//					colorBinary.at<Vec3b>(i, j) = color1;
//			}
//			/*if (Y[t] == 1) {
//				line(colorBinary, Point((min + t - k1) / 2, (k1 + min + t)/2), Point((min + t - k) / 2, (k + min + t) / 2), Scalar(255, 0, 0), 2.5, 8);
//			}*/
//			adder++;
//			if (adder % 2 == 0) //even run...
//				ystart++;
//
//			/*if(mark[t] == 1)
//				line(colorBinary, Point((min + t - k1) / 2, (k1 + min + t) / 2), Point((min + t - k) / 2, (k + min + t) / 2), Scalar(255, 0, 0), 4, 8);*/
//		}
//		counter++;
//
//	}
//
//	if (lastdiag < c) {
//		E = Point(lastdiag, 0);
//	}
//	else {
//		E = Point(c - 1, lastdiag - c + 1);
//	}
//	if (lastdiag < r) {
//		S = Point(0, lastdiag);
//	}
//	else {
//		S = Point(lastdiag - r + 1, r - 1);
//	}
//
//	//line(colorBinary, E, S, Scalar(0, 0, 255), 1, 8);
//	//int N = (lastdiag - firstdiag + 1) / 10;
//	imshow("stripLines.tif", colorBinary);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/DiagStripsProj.tif", colorBinary);
//	imshow("regions.tif", regions);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/DiagStripsRegions.tif", regions);
//	//imwrite("Data/ICDAR/Data218/OffDiagonalStrips.tif", colorBinary);
//}
//
//Point CCmedian(int **roi, int h, int w) {
//	std::vector<int> X, Y;
//	std::cout << "\npixel values = \n";
//	/*for (int i = 0; i < h; i++) {
//		for (int j = 0; j < w; j++) {
//			if(roi[i][j]==1)
//				std::cout << "\t" << roi[i][j];
//		}
//	}*/
//	for (int i = 0; i < h; i++) {
//		for (int j = 0; j < w; j++) {
//			if(roi[i][j] == 0) { 
//				X.push_back(j);
//				Y.push_back(i);
//			}
//		}
//	}
//	int sz = int(X.size());
//	std::cout << "\nsize = " << sz;
//	Point med = Point(int(median(X, sz)), int(median(Y, sz)));
//	return(med);
//}
//
//int findLongestWord(Mat bin_img, int **label, int labelNum, int **box, int r, int c) {
//	int max_hist = 0, max_lbl;
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1)
//			continue;
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//		
//		if (h < 15 || w < 15)
//			continue;
//
//		if (h > w )
//			continue;
//
//		int **roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				if (label[box[i][4] + k][box[i][2] + l] == i)
//					roi[k][l] = 0;
//				else
//					roi[k][l] = 255;
//			}
//		}
//		//find horizontal projection...
//		int *hist = new int[h];
//		for (int k = 0; k < h; k++) {
//			hist[k] = 0;
//		}
//		project0(roi, h, w, hist);
//		smoothingHistogram(hist, h, 6, 0.2);
//		int max = 0;
//		for (int k = 0; k < h; k++)
//			if (hist[k]>max)
//				max = hist[k];
//
//		if (max > w) {
//			Vec3b color1 = { 2, 34, 228 };
//			Vec3b color2 = { 4, 28, 146 };
//			for (int k = 0; k < h; k++) {
//				for (int l = 0; l < hist[k]; l++) {
//					if (l + box[i][2] >= c)
//						break;
//					if (label[box[i][4] + k][box[i][2] + l] == 0)
//						bin_img.at<Vec3b>(box[i][4] + k, box[i][2] + l) = color1;
//					else
//						bin_img.at<Vec3b>(box[i][4] + k, box[i][2] + l) = color2;
//				}
//			}
//			if (max > max_hist) {
//				max_hist = max;
//				max_lbl = i;
//			}
//		}		
//
//		for (int k = 0; k < h; k++)
//			delete[] roi[k];
//		delete[] roi;
//	}
//	std::cout << "\nThe length of the strip should be " << int(1.25 * (box[max_lbl][3] - box[max_lbl][2] + 1));
//	/*imshow("binaryImgProj.tif", bin_img);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/binImgProj.tif", bin_img);*/
//
//	return(int(1.25 * (box[max_lbl][3] - box[max_lbl][2] + 1)));
//}
//
//void fillHoles(int **img, int **label, int **box) {
//	std::vector<int> countCC[2];
//	int **countPixels = new int *[labelNum];
//	for (int i = 0; i < labelNum; i++) {
//		countPixels[i] = new int[2];
//		countPixels[i][0] = i;
//		countPixels[i][1] = 0;
//	}
//
//	for (int i = 2; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//		int count = 0;
//		for (int k = box[i][4]; k <= box[i][5]; k++) {
//			for (int l = box[i][2]; l <= box[i][3]; l++) {
//				if (label[k][l] == i)
//					count++;
//			}
//		}
//		countPixels[i][1] = count;
//		std::cout << "\n Label " << i << " has " << count << " pixels";
//		countCC[0].push_back(count);
//		countCC[1].push_back((box[i][3] - box[i][2])*(box[i][5] - box[i][4]));
//	}
//	int *arr = new int[int(countCC[0].size())];
//	int *area = new int[int(countCC[0].size())];
//	for (int i = 0; i<int(countCC[0].size()); i++) {
//		arr[i] = countCC[0][i];
//		area[i] = countCC[1][i];
//	}
//	double medCount = mean(arr, int(countCC[0].size()));
//	double medArea = mean(area, int(countCC[0].size()));
//	std::cout << "\nThe mean of CC size = " << medCount;
//	std::cout << "\nThe mean of CC Area = " << medArea;
//
//	for (int i = 2; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//		if (countPixels[i][1] > 3 * medCount && (box[i][3] - box[i][2])*(box[i][5] - box[i][4]) > 3.5*medArea) continue;
//		for (int k = box[i][4]; k <= box[i][5]; k++) {
//			for (int l = box[i][2]; l <= box[i][3]; l++) {
//				if (label[k][l] == i)
//					img[k][l] = 255;
//			}
//		}
//	}
//}
//
//
//void MultiOIC(int **img, Mat binimg) {
//	int grid = 4;
//	int **S, **oiclabel;
//		
//	int g = grid;
//	int h = r - ((r - 1) % g);
//	int w = c - ((c - 1) % g);
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
//				S[i][j] = img[i - 2*g][j - 2*g];
//		}
//	}
//	Mat colimg;
//	colimg.create(R, C, CV_8UC1);
//	arr2mat(S, colimg, R, C);
//	cvtColor(colimg, colimg, CV_GRAY2BGR);
//
//	//for (int i = 0; i < R/g - 1; i++) {
//	//	line(colimg, Point(0, i*g), Point(C-1, i*g), Scalar(0, 0, 255), 1, 8);
//	//}
//	//for (int j = 0; j < C / g - 1; j++) {
//	//	line(colimg, Point(j*g, 0), Point(j*g, R - 1), Scalar(0, 0, 255), 1, 8);
//	//}
//	//imwrite("Data/ICDAR/Data218/MultiOIC.tif", colimg);
//
//	vector<Point> sequence;
//	MakeMultiOIC(colimg, S, R, C, sequence, g);
//	int i = 0;
//	do{
//		Point q = sequence[i], curr = q, next = sequence[i + 1];
//		do {
//			line(colimg, curr, next, Scalar(0, 0, 255), 1, 8);
//			curr = next;
//			next = sequence[++i];
//		} while (q != curr);	
//		i++;
//	}while (i<int(sequence.size()) - 1);
//	imshow("MultiOIC.tif", colimg);
//	imwrite("Data/ICDAR/Data218/MultiOIC.tif", colimg);
//}
//
//int isBoundary(int roi[3][3]) {
//	int sum = 0;
//	for (int i = 0; i < 3; i++) {
//		for (int j = 0; j < 3; j++) {
//			if (i == 1 && j == 1) continue;
//			if (roi[i][j] == 1) sum++;
//		}
//	}
//	return(sum);
//}
//
//void findBoundary(int **img, int h, int w) {
//	int **copy = img;
//	for (int i = 0; i < h; i++) {
//		for (int j = 0; j < w; j++) {
//			if (copy[i][j] == 1) continue;
//			int roi[3][3];
//			for (int k = -1; k <= 1; k++) {
//				for (int l = -1; l <= 1; l++) {
//					if (i + k < 0 || i + k >= h || j + l < 0 || j + l >= w)
//						roi[k + 1][l + 1] = 1;
//					else
//						roi[k + 1][l + 1] = copy[i + k][j + l];
//				}
//			}
//			if (isBoundary)
//				img[i][j] = 2;
//			else
//				img[i][j] = 1;
//		}
//	}
//}
//
//
//void invCCProj(Mat img, int flag = 0) {
//	int **orgimg, **label, **box, **binimg;
//	Mat colorimg, bin_img, thin_im, near4img;
//
//	orgimg = new int*[r];
//	binimg = new int*[r];
//	label = new int*[r];
//	for (int i = 0; i < r; i++) {
//		orgimg[i] = new int[c];
//		binimg[i] = new int[c];
//		label[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			orgimg[i][j] = 0;
//			binimg[i][j] = 0;
//			label[i][j] = 0;
//		}
//	}
//	mat2arr(img, orgimg);
//
//	//----------------------------------------------------------------------------------------------------------
//	//Binarize if needed...
//	if (flag) {
//		NICK(orgimg, r, c, binimg);
//		labelNum = labelling(binimg, label, r, c);
//
//		//Bound the components by rectangular & diagonal boxes
//		box = new int*[labelNum];
//		for (int i = 0; i < labelNum; i++)
//			box[i] = new int[6];
//		boundingRectangles(box, label, r, c, labelNum);
//
//		for (int i = 1; i < labelNum; i++) {
//			if (box[i][2] == -1) continue;
//			int h = box[i][5] - box[i][4] + 1;
//			int w = box[i][3] - box[i][2] + 1;
//			int count = 0;
//			for (int k = box[i][4]; k <= box[i][5]; k++) {
//				for (int l = box[i][2]; l <= box[i][3]; l++) {
//					if (binimg[k][l] == 0) count++;
//				}
//			}
//			if (h*w < 30 || count < 100) {
//				for (int k = box[i][4]; k <= box[i][5]; k++) {
//					for (int l = box[i][2]; l <= box[i][3]; l++) {
//						if (label[k][l] == i) {
//							label[k][l] = 0;
//							binimg[k][l] = 255;
//						}
//					}
//				}
//				box[i][2] = -1;
//				box[i][3] = -1;
//				box[i][4] = -1;
//				box[i][5] = -1;
//			}
//		}
//
//		bin_img.create(r, c, CV_8UC1);
//		arr2mat(binimg, bin_img);
//		GaussianBlur(bin_img, bin_img, Size(3, 3), 0, 0);
//		mat2arr(bin_img, orgimg);
//		NICK(orgimg, r, c, binimg);
//		orgimg = binimg;
//		arr2mat(binimg, bin_img);
//		imshow("binaryImg.tif", bin_img);
//		imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/binImg-50.tif", bin_img);
//		
//		//----------------------------------------------------------------------------------------------------------
//		for (int i = 0; i < labelNum; i++)
//			delete[] box[i];
//		delete[] box;
//		for (int i = 0; i < r; i++) {
//			for (int j = 0; j < c; j++)
//				label[i][j] = 0;
//		}
//		labelNum = 0;
//	}
//		
//	//Fill holes...
//	erode(orgimg);
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			orgimg[i][j] = 255 - orgimg[i][j];
//		}
//	}
//	//dilate(orgimg);
//
//	labelNum = labelling(orgimg, label, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	//fillHoles(orgimg, label, box);
//
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			orgimg[i][j] = 255 - orgimg[i][j];
//		}
//	}
//
//	//---------------------------------------------------------------------------------------------------------
//	//re-label components...
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			label[i][j] = 0;
//		}
//	}
//	for (int i = 0; i < labelNum; i++)
//		delete[] box[i];
//	delete[] box;
//	labelNum = 0;
//	std::cout << "\nRe-labelling...";
//	labelNum = labelling(orgimg, label, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	//dilate(orgimg);
//	arr2mat(orgimg, img);
//	colorimg.create(r, c, CV_8UC1);
//	colorimg = img;
//
//	imshow("invcc.tif", colorimg);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/invcc-50.tif", colorimg);
//
//	cvtColor(img, colorimg, CV_GRAY2BGR);
//	near4img.create(r, c, CV_8UC1);
//	cvtColor(img, near4img, CV_GRAY2BGR);
//	
//	//imwrite("Data/ICDAR/Data218/invcc.tif", img);
//
//	//---------------------------------------------------------------------------------------------------------
//	//Exclude small components...
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//		int count = 0;
//		for (int k = box[i][4]; k <= box[i][5]; k++) {
//			for (int l = box[i][2]; l <= box[i][3]; l++) {
//				if (orgimg[k][l] == 0) count++;
//			}
//		}
//		if (h*w < 30 || count < 200) {
//			for (int k = 0; k < 6; k++)
//				box[i][k] = -1;
//		}
//	}
//
//	int N = findLongestWord(colorimg, label, labelNum, box, r, c);
//
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//		line(colorimg, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(255, 0, 0), 2, 8);
//		line(colorimg, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(255, 0, 0), 2, 8);
//		line(colorimg, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(255, 0, 0), 2, 8);
//		line(colorimg, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(255, 0, 0), 2, 8);
//	}
//	imshow("binaryImgProj.tif", colorimg);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/binImgProj-50.tif", colorimg);
//
//
//	//imshow("invcc.tif", near4img);
//	//imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/Near4Network.tif", near4img);
//
//	/*imshow("boundary.tif", thin_im);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/boundary.tif", thin_im);*/
//
//	//MultiOIC(orgimg, img);
//
//	StripVer(orgimg, img, N);
//	//StripHor(orgimg, img);
//	//StripDiag(orgimg, img);
//	//StripOffDiag(orgimg, img);
//
//	//Delete...
//	for (int i = 0; i < r; i++) {
//		delete[] orgimg[i];
//	}
//	delete[] orgimg;
//}
//
//int countTextPixel(Mat orgimg, int **lbl, int l, int h, int w, Point start) {
//	int count = 0;
//	for (int i = start.y; i <= start.y + h; i++) {
//		for (int j = start.x; j <= start.x + w; j++) {
//			if (orgimg.at<uchar>(i, j) == 0 && lbl[i][j] == l)
//				count++;
//		}
//	}
//	return(count);
//}
//
//void regionMask1(Mat orgimg) {
//	Mat S1 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/0Mask.tif", IMREAD_GRAYSCALE);
//
//	int **s1;
//	int **org;
//
//	Mat allMask, allMaskCol;
//	allMask.create(r, c, CV_8UC1);
//	allMaskCol.create(r, c, CV_8UC1);
//	int **all_mask;
//
//	org = new int*[r];
//	for (int i = 0; i < r; i++) {
//		org[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			org[i][j] = 0;
//		}
//	}
//	mat2arr(orgimg, org);
//
//	//label 1...
//	s1 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		s1[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			s1[i][j] = 0;
//		}
//	}
//	mat2arr(S1, s1);
//
//	all_mask = new int*[r];
//	for (int i = 0; i < r; i++) {
//		all_mask[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			all_mask[i][j] = 255 * (s1[i][j] / 255);
//			if (orgimg.at<uchar>(i, j) == 0) {
//				all_mask[i][j] = 150;
//			}
//		}
//	}
//
//	arr2mat(all_mask, allMask);
//	cvtColor(allMask, allMaskCol, CV_GRAY2BGR);
//	Vec3b color1 = { 255, 0, 255 };
//	Vec3b color2 = { 0, 255, 255 };
//	Vec3b color3 = { 255, 255, 0 };
//	Vec3b color4 = { 0, 120, 255 };
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (all_mask[i][j] == 150 && s1[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color1;
//			}
//		}
//	}
//
//	imshow("allMasks", allMask);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/SuperimposeMasks.tif", allMask);
//	imshow("allMasks", allMaskCol);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/SuperimposeMasksColored.tif", allMaskCol);
//
//}
//
//void regionMask2(Mat orgimg) {
//	Mat S1 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/30Mask.tif", IMREAD_GRAYSCALE);
//	Mat S2 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/-50Mask.tif", IMREAD_GRAYSCALE);
//
//	int **s1;
//	int **s2;
//	int **org;
//
//	int **label1;
//	int **label2;
//
//	int labelNum1 = 0;
//	int labelNum2 = 0;
//
//	Mat allMask, allMaskCol;
//	allMask.create(r, c, CV_8UC1);
//	allMaskCol.create(r, c, CV_8UC1);
//	int **all_mask;
//
//	org = new int*[r];
//	for (int i = 0; i < r; i++) {
//		org[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			org[i][j] = 0;
//		}
//	}
//	mat2arr(orgimg, org);
//
//	//label 1...
//	std::cout << "\nLabelling s1...";
//	s1 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		s1[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			s1[i][j] = 0;
//		}
//	}
//	mat2arr(S1, s1);
//
//	label1 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label1[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			label1[i][j] = 0;
//		}
//	}
//	labelNum1 = labelling(s1, label1, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	int **box1 = new int*[labelNum1];
//	for (int i = 0; i < labelNum1; i++)
//		box1[i] = new int[6];
//	boundingRectangles(box1, label1, r, c, labelNum1);
//
//	//label S2...
//	std::cout << "\nLabelling s2...";
//	s2 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		s2[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			s2[i][j] = 0;
//		}
//	}
//	mat2arr(S2, s2);
//
//	label2 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label2[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			label2[i][j] = 0;
//		}
//	}
//	labelNum2 = labelling(s2, label2, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	int **box2 = new int*[labelNum2];
//	for (int i = 0; i < labelNum2; i++)
//		box2[i] = new int[6];
//	boundingRectangles(box2, label2, r, c, labelNum2);
//
//	//compare overlapped regions...
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			int p1 = (255 - s1[i][j]) / 255;
//			int p2 = (255 - s2[i][j]) / 255;
//
//			//no component...
//			if (p1 + p2 == 0)
//				continue;
//			//no overlapping component...
//			if ((p1 && !p2 ) || (!p1 && p2 ) ) {
//				continue;
//			}
//
//			//s1 and s2 overlapped...
//			if (p1 && p2 ) {
//				int l1 = label1[i][j];
//				int l2 = label2[i][j];
//				double diag1 = (box1[l1][5] - box1[l1][4])*(box1[l1][5] - box1[l1][4]) + (box1[l1][3] - box1[l1][2])*(box1[l1][3] - box1[l1][2]);
//				double diag2 = (box2[l2][5] - box2[l2][4])*(box2[l2][5] - box2[l2][4]) + (box2[l2][3] - box2[l2][2])*(box2[l2][3] - box2[l2][2]);
//				std::cout << "\nDiag1 = " << sqrt(diag1) << " , Diag2 = " << sqrt(diag2);
//				int c1 = countTextPixel(orgimg, label1, l1, (box1[l1][5] - box1[l1][4]), (box1[l1][3] - box1[l1][2]), Point(box1[l1][2], box1[l1][4]));
//				int c2 = countTextPixel(orgimg, label2, l2, (box2[l2][5] - box2[l2][4]), (box2[l2][3] - box2[l2][2]), Point(box2[l2][2], box2[l2][4]));
//				std::cout << "\ncount1 = " << c1 << " , count2 = " << c2;
//				if (c1 < c2) {
//					std::cout << "\nRemove cc from s1";
//					for (int k = box1[l1][4]; k <= box1[l1][5]; k++) {
//						for (int l = box1[l1][2]; l <= box1[l1][3]; l++)
//							if (k < r && l < c && label1[k][l] == l1)
//								s1[k][l] = 255;
//					}
//				}
//				else {
//					std::cout << "\nRemove cc from s2";
//					for (int k = box2[l2][4]; k <= box2[l2][5]; k++) {
//						for (int l = box2[l2][2]; l <= box2[l2][3]; l++)
//							if (k < r && l < c && label2[k][l] == l2)
//								s2[k][l] = 255;
//					}
//				}
//			}
//		}
//	}
//
//
//	all_mask = new int*[r];
//	for (int i = 0; i < r; i++) {
//		all_mask[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			all_mask[i][j] = 255 * (s1[i][j] / 255)*(s2[i][j] / 255);
//			if (orgimg.at<uchar>(i, j) == 0) {
//				all_mask[i][j] = 150;
//			}
//		}
//	}
//
//	arr2mat(all_mask, allMask);
//	cvtColor(allMask, allMaskCol, CV_GRAY2BGR);
//	Vec3b color1 = { 255, 0, 255 };
//	Vec3b color2 = { 0, 255, 255 };
//	Vec3b color3 = { 255, 255, 0 };
//	Vec3b color4 = { 0, 120, 255 };
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (all_mask[i][j] == 150 && s1[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color1;
//			}
//			else if (all_mask[i][j] == 150 && s2[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color2;
//			}
//		}
//	}
//
//	imshow("allMasks", allMask);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/SuperimposeMasks.tif", allMask);
//	imshow("allMasks", allMaskCol);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/SuperimposeMasksColored.tif", allMaskCol);
//
//}
//
//void regionMask3(Mat orgimg) {
//	Mat S1 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/0Mask.tif", IMREAD_GRAYSCALE);
//	Mat S2 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/22.5Mask.tif", IMREAD_GRAYSCALE);
//	Mat S3 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/-30Mask.tif", IMREAD_GRAYSCALE);
//
//	int **s1;
//	int **s2;
//	int **s3;
//	int **org;
//
//	int **label1;
//	int **label2;
//	int **label3;
//
//	int labelNum1 = 0;
//	int labelNum2 = 0;
//	int labelNum3 = 0;
//
//	Mat allMask, allMaskCol;
//	allMask.create(r, c, CV_8UC1);
//	allMaskCol.create(r, c, CV_8UC1);
//	int **all_mask;
//
//	org = new int*[r];
//	for (int i = 0; i < r; i++) {
//		org[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			org[i][j] = 0;
//		}
//	}
//	mat2arr(orgimg, org);
//
//	//label 1...
//	std::cout << "\nLabelling s1...";
//	s1 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		s1[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			s1[i][j] = 0;
//		}
//	}
//	mat2arr(S1, s1);
//
//	label1 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label1[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			label1[i][j] = 0;
//		}
//	}
//	labelNum1 = labelling(s1, label1, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	int **box1 = new int*[labelNum1];
//	for (int i = 0; i < labelNum1; i++)
//		box1[i] = new int[6];
//	boundingRectangles(box1, label1, r, c, labelNum1);
//
//	//label S2...
//	std::cout << "\nLabelling s2...";
//	s2 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		s2[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			s2[i][j] = 0;
//		}
//	}
//	mat2arr(S2, s2);
//
//	label2 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label2[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			label2[i][j] = 0;
//		}
//	}
//	labelNum2 = labelling(s2, label2, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	int **box2 = new int*[labelNum2];
//	for (int i = 0; i < labelNum2; i++)
//		box2[i] = new int[6];
//	boundingRectangles(box2, label2, r, c, labelNum2);
//
//	//label S3...
//	std::cout << "\nLabelling s3...";
//	s3 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		s3[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			s3[i][j] = 0;
//		}
//	}
//	mat2arr(S3, s3);
//
//	label3 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label3[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			label3[i][j] = 0;
//		}
//	}
//	labelNum3 = labelling(s3, label3, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	int **box3 = new int*[labelNum3];
//	for (int i = 0; i < labelNum3; i++)
//		box3[i] = new int[6];
//	boundingRectangles(box3, label3, r, c, labelNum3);
//
//	//compare overlapped regions...
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			int p1 = (255 - s1[i][j]) / 255;
//			int p2 = (255 - s2[i][j]) / 255;
//			int p3 = (255 - s3[i][j]) / 255;
//
//			//no component...
//			if (p1 + p2 + p3 == 0)
//				continue;
//			//no overlapping component...
//			if ((p1 && !p2 && !p3 ) || (!p1 && p2 && !p3 ) || (!p1 && !p2 && p3 ) ) {
//				continue;
//			}
//
//			//s1 and s2 overlapped...
//			if (p1 && p2 && !p3 ) {
//				int l1 = label1[i][j];
//				int l2 = label2[i][j];
//				double diag1 = (box1[l1][5] - box1[l1][4])*(box1[l1][5] - box1[l1][4]) + (box1[l1][3] - box1[l1][2])*(box1[l1][3] - box1[l1][2]);
//				double diag2 = (box2[l2][5] - box2[l2][4])*(box2[l2][5] - box2[l2][4]) + (box2[l2][3] - box2[l2][2])*(box2[l2][3] - box2[l2][2]);
//				std::cout << "\nDiag1 = " << sqrt(diag1) << " , Diag2 = " << sqrt(diag2);
//				int c1 = countTextPixel(orgimg, label1, l1, (box1[l1][5] - box1[l1][4]), (box1[l1][3] - box1[l1][2]), Point(box1[l1][2], box1[l1][4]));
//				int c2 = countTextPixel(orgimg, label2, l2, (box2[l2][5] - box2[l2][4]), (box2[l2][3] - box2[l2][2]), Point(box2[l2][2], box2[l2][4]));
//				std::cout << "\ncount1 = " << c1 << " , count2 = " << c2;
//				if (c1 < c2) {
//					std::cout << "\nRemove cc from s1";
//					for (int k = box1[l1][4]; k <= box1[l1][5]; k++) {
//						for (int l = box1[l1][2]; l <= box1[l1][3]; l++)
//							if (k < r && l < c && label1[k][l] == l1)
//								s1[k][l] = 255;
//					}
//				}
//				else {
//					std::cout << "\nRemove cc from s2";
//					for (int k = box2[l2][4]; k <= box2[l2][5]; k++) {
//						for (int l = box2[l2][2]; l <= box2[l2][3]; l++)
//							if (k < r && l < c && label2[k][l] == l2)
//								s2[k][l] = 255;
//					}
//				}
//			}
//			//s1 and s3 overlapped...
//			else if (p1 && !p2 && p3 ) {
//				int l1 = label1[i][j];
//				int l3 = label3[i][j];
//				double diag1 = (box1[l1][5] - box1[l1][4])*(box1[l1][5] - box1[l1][4]) + (box1[l1][3] - box1[l1][2])*(box1[l1][3] - box1[l1][2]);
//				double diag3 = (box3[l3][5] - box3[l3][4])*(box3[l3][5] - box3[l3][4]) + (box3[l3][3] - box3[l3][2])*(box3[l3][3] - box3[l3][2]);
//				std::cout << "\nDiag1 = " << sqrt(diag1) << " , Diag3 = " << sqrt(diag3);
//				int c1 = countTextPixel(orgimg, label1, l1, (box1[l1][5] - box1[l1][4]), (box1[l1][3] - box1[l1][2]), Point(box1[l1][2], box1[l1][4]));
//				int c3 = countTextPixel(orgimg, label3, l3, (box3[l3][5] - box3[l3][4]), (box3[l3][3] - box3[l3][2]), Point(box3[l3][2], box3[l3][4]));
//				std::cout << "\ncount1 = " << c1 << " , count3 = " << c3;
//				if (c1 < c3) {
//					std::cout << "\nRemove cc from s1";
//					for (int k = box1[l1][4]; k <= box1[l1][5]; k++) {
//						for (int l = box1[l1][2]; l <= box1[l1][3]; l++)
//							if (k < r && l < c && label1[k][l] == l1)
//								s1[k][l] = 255;
//					}
//				}
//				else {
//					std::cout << "\nRemove cc from s3";
//					for (int k = box3[l3][4]; k <= box3[l3][5]; k++) {
//						for (int l = box3[l3][2]; l <= box3[l3][3]; l++)
//							if (k < r && l < c && label3[k][l] == l3)
//								s3[k][l] = 255;
//					}
//				}
//			}
//			//s2 and s3 overlapped...
//			else if (!p1 && p2 && p3 ) {
//				int l2 = label2[i][j];
//				int l3 = label3[i][j];
//				double diag2 = (box2[l2][5] - box2[l2][4])*(box2[l2][5] - box2[l2][4]) + (box2[l2][3] - box2[l2][2])*(box2[l2][3] - box2[l2][2]);
//				double diag3 = (box3[l3][5] - box3[l3][4])*(box3[l3][5] - box3[l3][4]) + (box3[l3][3] - box3[l3][2])*(box3[l3][3] - box3[l3][2]);
//				std::cout << "\nDiag2 = " << sqrt(diag2) << " , Diag3 = " << sqrt(diag3);
//				int c2 = countTextPixel(orgimg, label2, l2, (box2[l2][5] - box2[l2][4]), (box2[l2][3] - box2[l2][2]), Point(box2[l2][2], box2[l2][4]));
//				int c3 = countTextPixel(orgimg, label3, l3, (box3[l3][5] - box3[l3][4]), (box3[l3][3] - box3[l3][2]), Point(box3[l3][2], box3[l3][4]));
//				std::cout << "\ncount2 = " << c2 << " , count3 = " << c3;
//				if (c2 < c3) {
//					std::cout << "\nRemove cc from s2";
//					for (int k = box2[l2][4]; k <= box2[l2][5]; k++) {
//						for (int l = box2[l2][2]; l <= box2[l2][3]; l++)
//							if (k < r && l < c && label2[k][l] == l2)
//								s2[k][l] = 255;
//					}
//				}
//				else {
//					std::cout << "\nRemove cc from s3";
//					for (int k = box3[l3][4]; k <= box3[l3][5]; k++) {
//						for (int l = box3[l3][2]; l <= box3[l3][3]; l++)
//							if (k < r && l < c && label3[k][l] == l3)
//								s3[k][l] = 255;
//					}
//				}
//			}
//		}
//	}
//
//
//	all_mask = new int*[r];
//	for (int i = 0; i < r; i++) {
//		all_mask[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			all_mask[i][j] = 255 * (s1[i][j] / 255)*(s2[i][j] / 255)*(s3[i][j] / 255);
//			if (orgimg.at<uchar>(i, j) == 0) {
//				all_mask[i][j] = 150;
//			}
//		}
//	}
//
//	arr2mat(all_mask, allMask);
//	cvtColor(allMask, allMaskCol, CV_GRAY2BGR);
//	Vec3b color1 = { 255, 0, 255 };
//	Vec3b color2 = { 0, 255, 255 };
//	Vec3b color3 = { 255, 255, 0 };
//	Vec3b color4 = { 0, 120, 255 };
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (all_mask[i][j] == 150 && s1[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color1;
//			}
//			if (all_mask[i][j] == 150 && s2[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color2;
//			}
//			else if (all_mask[i][j] == 150 && s3[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color3;
//			}
//		}
//	}
//
//	imshow("allMasks", allMask);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/SuperimposeMasks.tif", allMask);
//	imshow("allMasks", allMaskCol);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/SuperimposeMasksColored.tif", allMaskCol);
//
//}
//
//void regionMask4(Mat orgimg) {
//	Mat S1 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/0Mask.tif", IMREAD_GRAYSCALE);
//	Mat S2 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/15Mask.tif", IMREAD_GRAYSCALE);
//	Mat S3 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/-30Mask.tif", IMREAD_GRAYSCALE);
//	Mat S4 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/-45Mask.tif", IMREAD_GRAYSCALE);
//
//	int **s1;
//	int **s2;
//	int **s3;
//	int **s4;
//	int **org;
//
//	int **label1;
//	int **label2;
//	int **label3;
//	int **label4;
//
//	int labelNum1 = 0;
//	int labelNum2 = 0;
//	int labelNum3 = 0;
//	int labelNum4 = 0;
//
//	Mat allMask, allMaskCol;
//	allMask.create(r, c, CV_8UC1);
//	allMaskCol.create(r, c, CV_8UC1);
//	int **all_mask;
//		
//	org = new int*[r];
//	for (int i = 0; i < r; i++) {
//		org[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			org[i][j] = 0;
//		}
//	}
//	mat2arr(orgimg, org);
//
//	//label 1...
//	std::cout<<"\nLabelling s1...";
//	s1 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		s1[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			s1[i][j] = 0;
//		}
//	}
//	mat2arr(S1, s1);
//	
//	label1 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label1[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			label1[i][j] = 0;
//		}
//	}
//	labelNum1 = labelling(s1, label1, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	int **box1 = new int*[labelNum1];
//	for (int i = 0; i < labelNum1; i++)
//		box1[i] = new int[6];
//	boundingRectangles(box1, label1, r, c, labelNum1);
//
//	//label S2...
//	std::cout << "\nLabelling s2...";
//	s2 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		s2[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			s2[i][j] =  0;
//		}
//	}
//	mat2arr(S2, s2);
//
//	label2 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label2[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			label2[i][j] = 0;
//		}
//	}
//	labelNum2 = labelling(s2, label2, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	int **box2 = new int*[labelNum2];
//	for (int i = 0; i < labelNum2; i++)
//		box2[i] = new int[6];
//	boundingRectangles(box2, label2, r, c, labelNum2);
//
//	//label S3...
//	std::cout << "\nLabelling s3...";
//	s3 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		s3[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			s3[i][j] = 0;
//		}
//	}
//	mat2arr(S3, s3);
//
//	label3 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label3[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			label3[i][j] = 0;
//		}
//	}
//	labelNum3 = labelling(s3, label3, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	int **box3 = new int*[labelNum3];
//	for (int i = 0; i < labelNum3; i++)
//		box3[i] = new int[6];
//	boundingRectangles(box3, label3, r, c, labelNum3);
//
//	//label S4...
//	std::cout << "\nLabelling s4...";
//	s4 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		s4[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			s4[i][j] = 0;
//		}
//	}
//	mat2arr(S4, s4);
//
//	label4 = new int*[r];
//	for (int i = 0; i < r; i++) {
//		label4[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			label4[i][j] = 0;
//		}
//	}
//	labelNum4 = labelling(s4, label4, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	int **box4 = new int*[labelNum4];
//	for (int i = 0; i < labelNum4; i++)
//		box4[i] = new int[6];
//	boundingRectangles(box4, label4, r, c, labelNum4);
//
//	//compare overlapped regions...
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			int p1 = (255 - s1[i][j]) / 255;
//			int p2 = (255 - s2[i][j]) / 255;
//			int p3 = (255 - s3[i][j]) / 255;
//			int p4 = (255 - s4[i][j]) / 255;
//
//			//no component...
//			if (p1 + p2 + p3 + p4 == 0) 
//				continue;
//			//no overlapping component...
//			if ((p1 && !p2 && !p3 && !p4) || (!p1 && p2 && !p3 && !p4) || (!p1 && !p2 && p3 && !p4) || (!p1 && !p2 && !p3 && p4)) {
//				continue;
//			}
//
//			//s1 and s2 overlapped...
//			if (p1 && p2 && !p3 && !p4) {
//				int l1 = label1[i][j];
//				int l2 = label2[i][j];
//				double diag1 = (box1[l1][5] - box1[l1][4])*(box1[l1][5] - box1[l1][4]) + (box1[l1][3] - box1[l1][2])*(box1[l1][3] - box1[l1][2]);
//				double diag2 = (box2[l2][5] - box2[l2][4])*(box2[l2][5] - box2[l2][4]) + (box2[l2][3] - box2[l2][2])*(box2[l2][3] - box2[l2][2]);
//				int c1 = countTextPixel(orgimg, label1, l1, (box1[l1][5] - box1[l1][4]), (box1[l1][3] - box1[l1][2]), Point(box1[l1][2], box1[l1][4]));
//				int c2 = countTextPixel(orgimg, label2, l2, (box2[l2][5] - box2[l2][4]), (box2[l2][3] - box2[l2][2]), Point(box2[l2][2], box2[l2][4]));
//				std::cout << "\ncount1 = " << c1 << " , count2 = " << c2;
//				if (c1 < c2) {
//					for (int k = box1[l1][4]; k <= box1[l1][5]; k++) {
//						for (int l = box1[l1][2]; l <= box1[l1][3]; l++)
//							if (k < r && l < c && label1[k][l] == l1)
//								s1[k][l] = 255;
//					}
//				}
//				else {
//					for (int k = box2[l2][4]; k <= box2[l2][5]; k++) {
//						for (int l = box2[l2][2]; l <= box2[l2][3]; l++)
//							if (k < r && l < c && label2[k][l] == l2)
//								s2[k][l] = 255;
//					}
//				}
//			}
//			//s1 and s3 overlapped...
//			else if (p1 && !p2 && p3 && !p4) {
//				int l1 = label1[i][j];
//				int l3 = label3[i][j];
//				double diag1 = (box1[l1][5] - box1[l1][4])*(box1[l1][5] - box1[l1][4]) + (box1[l1][3] - box1[l1][2])*(box1[l1][3] - box1[l1][2]);
//				double diag3 = (box3[l3][5] - box3[l3][4])*(box3[l3][5] - box3[l3][4]) + (box3[l3][3] - box3[l3][2])*(box3[l3][3] - box3[l3][2]);
//				int c1 = countTextPixel(orgimg, label1, l1, (box1[l1][5] - box1[l1][4]), (box1[l1][3] - box1[l1][2]), Point(box1[l1][2], box1[l1][4]));
//				int c3 = countTextPixel(orgimg, label3, l3, (box3[l3][5] - box3[l3][4]), (box3[l3][3] - box3[l3][2]), Point(box3[l3][2], box3[l3][4]));
//				std::cout << "\ncount1 = " << c1 << " , count3 = " << 3;
//				if (c1 < c3) {
//					for (int k = box1[l1][4]; k <= box1[l1][5]; k++) {
//						for (int l = box1[l1][2]; l <= box1[l1][3]; l++)
//							if (k < r && l < c && label1[k][l] == l1)
//								s1[k][l] = 255;
//					}
//				}
//				else {
//					for (int k = box3[l3][4]; k <= box3[l3][5]; k++) {
//						for (int l = box3[l3][2]; l <= box3[l3][3]; l++)
//							if (k < r && l < c && label3[k][l] == l3)
//								s3[k][l] = 255;
//					}
//				}
//			}
//			//s1 and s4 overlapped...
//			else if (p1 && !p2 && !p3 && p4) {
//				int l1 = label1[i][j];
//				int l4 = label4[i][j];
//				double diag1 = (box1[l1][5] - box1[l1][4])*(box1[l1][5] - box1[l1][4]) + (box1[l1][3] - box1[l1][2])*(box1[l1][3] - box1[l1][2]);
//				double diag4 = (box4[l4][5] - box4[l4][4])*(box4[l4][5] - box4[l4][4]) + (box4[l4][3] - box4[l4][2])*(box4[l4][3] - box4[l4][2]);
//				int c1 = countTextPixel(orgimg, label1, l1, (box1[l1][5] - box1[l1][4]), (box1[l1][3] - box1[l1][2]), Point(box1[l1][2], box1[l1][4]));
//				int c4 = countTextPixel(orgimg, label4, l4, (box4[l4][5] - box4[l4][4]), (box4[l4][3] - box4[l4][2]), Point(box4[l4][2], box4[l4][4]));
//				std::cout << "\ncount1 = " << c1 << " , count4 = " << 4;
//				if (c1 < c4) {
//					for (int k = box1[l1][4]; k <= box1[l1][5]; k++) {
//						for (int l = box1[l1][2]; l <= box1[l1][3]; l++)
//							if (k < r && l < c && label1[k][l] == l1)
//								s1[k][l] = 255;
//					}
//				}
//				else {
//					for (int k = box4[l4][4]; k <= box4[l4][5]; k++) {
//						for (int l = box4[l4][2]; l <= box4[l4][3]; l++)
//							if (k < r && l < c && label4[k][l] == l4)
//								s4[k][l] = 255;
//					}
//				}
//			}
//			//s2 and s3 overlapped...
//			else if (!p1 && p2 && p3 && !p4) {
//				int l2 = label2[i][j];
//				int l3 = label3[i][j];
//				double diag2 = (box2[l2][5] - box2[l2][4])*(box2[l2][5] - box2[l2][4]) + (box2[l2][3] - box2[l2][2])*(box2[l2][3] - box2[l2][2]);
//				double diag3 = (box3[l3][5] - box3[l3][4])*(box3[l3][5] - box3[l3][4]) + (box3[l3][3] - box3[l3][2])*(box3[l3][3] - box3[l3][2]);
//				int c2 = countTextPixel(orgimg, label2, l2, (box2[l2][5] - box2[l2][4]), (box2[l2][3] - box2[l2][2]), Point(box2[l2][2], box2[l2][4]));
//				int c3 = countTextPixel(orgimg, label3, l3, (box3[l3][5] - box3[l3][4]), (box3[l3][3] - box3[l3][2]), Point(box3[l3][2], box3[l3][4]));
//				std::cout << "\ncount2 = " << c2 << " , count3 = " << 3;
//				if (c2 < c3) {
//					for (int k = box2[l2][4]; k <= box2[l2][5]; k++) {
//						for (int l = box2[l2][2]; l <= box2[l2][3]; l++)
//							if (k < r && l < c && label2[k][l] == l2)
//								s2[k][l] = 255;
//					}
//				}
//				else {
//					for (int k = box3[l3][4]; k <= box3[l3][5]; k++) {
//						for (int l = box3[l3][2]; l <= box3[l3][3]; l++)
//							if (k < r && l < c && label3[k][l] == l3)
//								s3[k][l] = 255;
//					}
//				}
//			}
//			//s2 and s4 overlapped...
//			else if (!p1 && p2 && !p3 && p4) {
//				int l2 = label2[i][j];
//				int l4 = label4[i][j];
//				double diag2 = (box2[l2][5] - box2[l2][4])*(box2[l2][5] - box2[l2][4]) + (box2[l2][3] - box2[l2][2])*(box2[l2][3] - box2[l2][2]);
//				double diag4 = (box4[l4][5] - box4[l4][4])*(box4[l4][5] - box4[l4][4]) + (box4[l4][3] - box4[l4][2])*(box4[l4][3] - box4[l4][2]);
//				int c2 = countTextPixel(orgimg, label2, l2, (box2[l2][5] - box2[l2][4]), (box2[l2][3] - box2[l2][2]), Point(box2[l2][2], box2[l2][4]));
//				int c4 = countTextPixel(orgimg, label4, l4, (box4[l4][5] - box4[l4][4]), (box4[l4][3] - box4[l4][2]), Point(box4[l4][2], box4[l4][4]));
//				std::cout << "\ncount2 = " << c2 << " , count4 = " << 4;
//				if (c2 < c4) {
//					for (int k = box2[l2][4]; k <= box2[l2][5]; k++) {
//						for (int l = box2[l2][2]; l <= box2[l2][3]; l++)
//							if (k < r && l < c && label2[k][l] == l2)
//								s2[k][l] = 255;
//					}
//				}
//				else {
//					for (int k = box4[l4][4]; k <= box4[l4][5]; k++) {
//						for (int l = box4[l4][2]; l <= box4[l4][3]; l++)
//							if (k < r && l < c && label4[k][l] == l4)
//								s4[k][l] = 255;
//					}
//				}
//			}
//			//s3 and s4 overlapped...
//			else if (!p1 && !p2 && p3 && p4) {
//				int l3 = label3[i][j];
//				int l4 = label4[i][j];
//				double diag3 = (box3[l3][5] - box3[l3][4])*(box3[l3][5] - box3[l3][4]) + (box3[l3][3] - box3[l3][2])*(box3[l3][3] - box3[l3][2]);
//				double diag4 = (box4[l4][5] - box4[l4][4])*(box4[l4][5] - box4[l4][4]) + (box4[l4][3] - box4[l4][2])*(box4[l4][3] - box4[l4][2]);
//				int c3 = countTextPixel(orgimg, label3, l3, (box3[l3][5] - box3[l3][4]), (box3[l3][3] - box3[l3][2]), Point(box3[l3][2], box3[l3][4]));
//				int c4 = countTextPixel(orgimg, label4, l4, (box4[l4][5] - box4[l4][4]), (box4[l4][3] - box4[l4][2]), Point(box4[l4][2], box4[l4][4]));
//				std::cout << "\ncount3 = " << c3 << " , count4 = " << 4;
//				if (c3 < c4) {
//					for (int k = box3[l3][4]; k <= box3[l3][5]; k++) {
//						for (int l = box3[l3][2]; l <= box3[l3][3]; l++)
//							if (k < r && l < c && label3[k][l] == l3)
//								s3[k][l] = 255;
//					}
//				}
//				else {
//					for (int k = box4[l4][4]; k <= box4[l4][5]; k++) {
//						for (int l = box4[l4][2]; l <= box4[l4][3]; l++)
//							if (k < r && l < c && label4[k][l] == l4)
//								s4[k][l] = 255;
//					}
//				}
//			}
//		}
//	}
//
//
//	all_mask = new int*[r];
//	for (int i = 0; i < r; i++) {
//		all_mask[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			all_mask[i][j] = 255 * (s1[i][j] / 255)*(s2[i][j] / 255)*(s3[i][j] / 255)*(s4[i][j] / 255);
//			if (orgimg.at<uchar>(i, j) == 0) {
//				all_mask[i][j] = 150;
//			}
//		}
//	}
//
//	arr2mat(all_mask, allMask);
//	cvtColor(allMask, allMaskCol, CV_GRAY2BGR);
//	Vec3b color1 = { 255, 0, 255 };
//	Vec3b color2 = { 0, 255, 255 };
//	Vec3b color3 = { 255, 255, 0 };
//	Vec3b color4 = { 0, 120, 255 };
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (all_mask[i][j] == 150 && s1[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color1;
//			}
//			if (all_mask[i][j] == 150 && s2[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color2;
//			}
//			else if (all_mask[i][j] == 150 && s3[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color3;
//			}
//			else if (all_mask[i][j] == 150 && s4[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color4;
//			}
//		}
//	}
//
//	imshow("allMasks", allMask);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/SuperimposeMasks.tif", allMask);
//	imshow("allMasks", allMaskCol);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/SuperimposeMasksColored.tif", allMaskCol);
//
//}
//
///*****************************************************************************************************
//											GRAPH BASED METHOD
//******************************************************************************************************/
//double HDistance(int **roi, int h, int w, int hastart, int hbstart, int hastop, int hbstop, int wastart, int wbstart, int wastop, int wbstop, int La, int Lb) {
//	double min_dAB, minmin_dAB = 99999;
//	Point arga;
//	for (int i = hastart; i < hastop; i++) {
//		for (int j = wastart; j < wastop; j++) {
//			if (roi[i][j] != La) continue;
//			arga = Point(j, i);
//			min_dAB = 99999;
//			for (int k = hbstart; k < hbstop; k++) {
//				for (int l = wbstart; l < wbstop; l++) {
//					if (roi[k][l] != Lb) continue;
//					if (min_dAB > norm(arga - Point(l, k))) {
//						min_dAB = norm(arga - Point(l, k));
//					}
//				}
//			}
//			if (min_dAB < minmin_dAB) {
//				minmin_dAB = min_dAB;
//			}
//		}
//	}
//	return(minmin_dAB);
//}
//
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
//int findTerminals(int **img, int h, int w, std::vector<Point> &terminal) {
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
//			int sum = checkNbd(copy, h + 2, w + 2, Point(j, i));
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
//int checkCorner(int **img, int h, int w, Point ref) {
//	int nbrs[10] = { 0 };
//	nbrs[0] = img[ref.y - 1][ref.x - 1];
//	nbrs[1] = img[ref.y - 1][ref.x];
//	nbrs[2] = img[ref.y - 1][ref.x + 1];
//	nbrs[3] = img[ref.y][ref.x + 1];
//	nbrs[4] = img[ref.y + 1][ref.x + 1];
//	nbrs[5] = img[ref.y + 1][ref.x];
//	nbrs[6] = img[ref.y + 1][ref.x - 1];
//	nbrs[7] = img[ref.y][ref.x - 1];
//	nbrs[8] = img[ref.y - 1][ref.x - 1];
//	nbrs[9] = img[ref.y - 1][ref.x];
//
//	int E, W, N, S, NE, NW, SE, SW;
//	NW = img[ref.y - 1][ref.x - 1];
//	N = img[ref.y - 1][ref.x];
//	NE = img[ref.y - 1][ref.x + 1];
//	E = img[ref.y][ref.x + 1];
//	SE = img[ref.y + 1][ref.x + 1];
//	S = img[ref.y + 1][ref.x];
//	SW = img[ref.y + 1][ref.x - 1];
//	W = img[ref.y][ref.x - 1];
//
//	int flag;
//	if ((S && W && (N || NE)) || (S && E && (W || NW)) || (N && E && (S || SW)) || (N && W && (E || SE)) ||
//		(S && NE && NW && !N) || (W && SE && NE && !E) || (N && SE && SW && !S) || (E && NW && SW && !W) ||
//		(SE && SW && NE && !S && !E) || (SE && SW && NW && !S && !W) || (NE && NW && SE && !N && !E) || (NE && NW && SW && !N && !W)) {
//		flag = 1;
//	}
//	else
//		flag = 0;
//
//	//for (int k = 0; k < 8; k++) {
//	//	if (nbrs[k] + nbrs[k + 1] + nbrs[k + 2] == 0) { //3 consecutive whites
//	//		std::cout << "\nNo\n" << nbrs[0] << nbrs[1] << nbrs[2] << nbrs[3] << nbrs[4] << nbrs[5] << nbrs[6] << nbrs[7];
//	//		flag = 0;
//	//		break;
//	//	}
//	//}
//	//if (flag) {
//	//	std::cout << "\nYES\n" << nbrs[0] << nbrs[1] << nbrs[2] << nbrs[3] << nbrs[4] << nbrs[5] << nbrs[6] << nbrs[7];
//	//}
//	return(flag);
//}
//
//int findCorners(int **img, int h, int w, std::vector<Point> &corners) {
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
//			int flag = checkCorner(copy, h + 2, w + 2, Point(j, i));
//			if (flag) {
//				corners.push_back(Point(j - 1, i - 1));
//				sz++;
//			}
//		}
//	}
//	return(sz);
//}
//
//int traceLine(int **roi, int h, int w, Point terminal, int theta, std::vector<Point> &sequence) {
//	Point curr = terminal, next;
//	int rho = 10, del = 15;
//	double slope = tan(theta*CV_PI / 180);
//	next.x = curr.x + 10;
//	next.y = curr.y + slope * 10;
//	int flag = 1, sz = 1;
//
//	if (next.x < 0 || next.x >= w || next.y < 0 || next.y >= h) {
//		flag = 0;
//		return(flag);
//	}
//
//	sequence.push_back(next);
//	sz++;
//	while (flag) {
//		Point p;
//		p = next - Point(rho * cos(theta*CV_PI / 180), rho * sin(theta*CV_PI / 180));
//		if (p.x < 0 || p.x >= w || p.y < 0 || p.y >= h) { //going out of block...
//			flag = 0;
//			break;
//		}
//		/*if (roi[p.y][p.x] == 1) {
//		flag = 0;
//		break;
//		}*/
//		int negflag = 1, posflag = 1;
//		int new_theta_neg, new_theta_pos;
//		int bw[31];
//		int countblack = 0;
//		//std::cout << "\n";
//		for (int k = theta - del; k <= theta + del; k++) {
//			Point p;
//			p = next - Point(rho * cos(k*CV_PI / 180), rho * sin(k*CV_PI / 180));
//			if (p.x < 0 || p.x >= w || p.y < 0 || p.y >= h) {
//				flag = 0;
//				break;
//			}
//			bw[k - theta + del] = roi[p.y][p.x];
//			if (roi[p.y][p.x] == 0) countblack++;
//			//std::cout << "   " << bw[k - theta + del];
//		}
//		if (countblack == 0) { // no black...
//			flag = 0;
//			break;
//		}
//		int next_theta;
//		int count = 0, maxcount = -1, first = 0;
//		for (int k = 0; k < 2 * del + 1; k++) {
//			if (bw[k] == 0) {
//				if (first == 0) first = 1; //start consecutive blacks...
//				count++;
//			}
//			if ((bw[k] == 1 && first == 1) || k == 2 * del) { //end of blacks...
//				if (count > maxcount) {
//					next_theta = (theta - del) + (k - int(count / 2));
//					maxcount = count;
//					count = 0;
//					first = 0;
//				}
//				else if (count == maxcount) {
//					if (abs(next_theta - theta) > abs(del - k + int(count / 2))) {
//						next_theta = (theta - del) + (k - int(count / 2));
//						maxcount = count;
//						count = 0;
//						first = 0;
//					}
//
//				}
//			}
//		}
//		if (maxcount < 0) {
//			flag = 0;
//			break;
//		}
//		//std::cout << " => next = " << next_theta << " degrees";
//
//		curr = next;
//		next = next - Point(rho * cos(next_theta*CV_PI / 180), rho * sin(next_theta*CV_PI / 180));
//		sequence.push_back(next);
//		theta = next_theta;
//		sz++;
//	}
//	return(sz);
//}
//
//int longestRay(int **roi, int h, int w, Point terminal, int theta, std::vector<Point> &sequence) {
//	int del = 7;
//	int black_count, total_count;
//	int next_theta = theta, prev_theta = theta;
//	Point curr = terminal;
//	Point next;
//	double slope = tan(theta*CV_PI) / 180;
//
//
//	/*while (next_theta <= theta + del) {
//	black_count = 0;
//	total_count = 0;
//	double slope = tan(next_theta*CV_PI / 180);
//	Point last_black;
//	for (int x = terminal.x + 1; x < w; x++) {
//	int y = terminal.y + slope*(x - terminal.x);
//	if (y < 0 || y >= h) break;
//	if (roi[y][x] != 0) break;
//	if (roi[y][x] == 0) {
//	black_count++;
//	last_black = Point(x, y);
//	}
//	}
//	int dst = black_count;
//	if (dst > distance) {
//	distance = dst;
//	last_point = last_black;
//	}
//	next_theta += del;
//	}*/
//	int x = terminal.x;
//	int y;
//	while (roi[curr.y][curr.x] == 0) {
//		sequence.push_back(curr);
//		x = x + del;
//		y = curr.y + slope * del;
//
//		if (x < 0 || x >= w || y < 0 || y >= h)
//			break;
//
//		int first = y, last = y;
//		while (roi[first][x] == 0 && first > 0)
//			first--;
//		first++;
//		while (roi[last][x] == 0 && last < h - 1)
//			last++;
//		last--;
//		next = Point(x, int(first + last) / 2);
//		int next_theta = atan2(next.y - curr.y, next.x - curr.x) * 180 / CV_PI;
//		if (next_theta < 0)
//			next_theta += 180;
//		/*if (abs(next_theta - prev_theta) > 20)
//		break;*/
//		slope = tan(next_theta * CV_PI / 180);
//		curr = next;
//		prev_theta = next_theta;
//		if (curr.x < 0 || curr.x >= w || curr.y < 0 || curr.y >= h)
//			break;
//	}
//	return(int(sequence.size()));
//}
//
//void connectLines(std::vector<int> Lines[6], int S, int theta) {
//	//Sort the list...
//	int max0 = -1, max1 = -1;
//	for (int i = 0; i < S - 1; i++) {
//		for (int j = 0; j < S - i - 1; j++) {
//			if (Lines[1][j] < Lines[1][j + 1] || (Lines[1][j] == Lines[1][j + 1] && Lines[0][j] > Lines[0][j + 1])) {
//				int *temp = new int[6];
//				for (int k = 0; k < 6; k++) {
//					temp[k] = Lines[k][j];
//					Lines[k][j] = Lines[k][j + 1];
//					Lines[k][j + 1] = temp[k];
//				}
//			}
//		}
//	}
//	int **flag = new int*[S];
//	for (int i = 0; i < S; i++) {
//		flag[i] = new int[3];
//		for (int j = 0; j < 3; j++)
//			flag[i][j] = -1;
//	}
//
//	std::vector<std::vector<int>> lineSet;
//	int count = 1;
//	for (int i = 0; i < S; i++) {
//		if (flag[i][0] > 0) continue;
//		std::vector<int> seq;
//		seq.push_back(count); //store the line no.
//		flag[i][0] = seq[0];
//		count++;
//		seq.push_back(Lines[0][i]); //store the label no.
//		Point left = Point(Lines[3][i], Lines[2][i]);
//		Point right = Point(Lines[5][i], Lines[4][i]);
//		int ref_slope = atan2(right.y - left.y, right.x - left.x) * 180 / CV_PI;
//		if (ref_slope < 0) ref_slope = 180 + ref_slope;
//		int f = 1, ref = i;
//		do {
//			int next_right = -1, next_left = -1;
//			int slope = -360;
//			for (int k = 0; k < S; k++) {
//				if (ref == k || Lines[3][k] < right.x) continue;
//				int th = atan2(Lines[2][k] - left.y, Lines[3][k] - left.x) * 180 / CV_PI;
//				if (th < 0) th = 180 + th;
//				if (abs(th - ref_slope) < abs(slope - ref_slope)) {
//					slope = th;
//					next_right = k;
//				}
//			}
//			if (next_right >= 0 && abs(slope - ref_slope) > 10)
//				f = 0;
//			else if (next_right >= 0 && flag[next_right][0] == -1) {
//				flag[next_right][0] = seq[0]; //got the next right...
//				flag[next_right][1] = slope; //left angle of next right...
//				seq.push_back(Lines[0][next_right]); //store the next label to the right...
//				left = Point(Lines[3][next_right], Lines[2][next_right]);
//				right = Point(Lines[5][next_right], Lines[4][next_right]);
//				ref = next_right;
//			}
//			else
//				f = 0;
//
//		} while (f);
//
//		lineSet.push_back(seq);
//	}
//	for (int i = 0; i<int(lineSet.size()); i++) {
//		std::cout << "\nLine no. ";
//		for (int j = 0; j<int(lineSet[i].size()); j++) {
//			std::cout << lineSet[i][j] << " -> ";
//		}
//		std::cout << "end";
//	}
//}
//
//int corrCoeff(Point p[4]) {
//	int minCoeff;
//	int cc = p[0].x*(p[1].y - p[2].y) + p[1].x*(p[2].y - p[0].y) + p[2].x*(p[0].y - p[1].y);
//	minCoeff = cc;
//	cc = p[0].x*(p[1].y - p[3].y) + p[1].x*(p[3].y - p[0].y) + p[3].x*(p[0].y - p[1].y);
//	if (cc < minCoeff)
//		minCoeff = cc;
//	cc = p[0].x*(p[2].y - p[3].y) + p[2].x*(p[3].y - p[0].y) + p[3].x*(p[0].y - p[2].y);
//	if (cc < minCoeff)
//		minCoeff = cc;
//	cc = p[1].x*(p[2].y - p[3].y) + p[2].x*(p[3].y - p[1].y) + p[3].x*(p[1].y - p[2].y);
//	if (cc < minCoeff)
//		minCoeff = cc;
//
//	return(minCoeff);
//}
//
//void drawLines(std::vector<int> Lines[6], int S, int theta) {
//	//Sort the list...
//	int max0 = -1, max1 = -1;
//	for (int i = 0; i < S - 1; i++) {
//		for (int j = 0; j < S - i - 1; j++) {
//			if (Lines[1][j] < Lines[1][j + 1] || (Lines[1][j] == Lines[1][j + 1] && Lines[0][j] > Lines[0][j + 1])) {
//				int *temp = new int[6];
//				for (int k = 0; k < 6; k++) {
//					temp[k] = Lines[k][j];
//					Lines[k][j] = Lines[k][j + 1];
//					Lines[k][j + 1] = temp[k];
//				}
//			}
//		}
//	}
//	int **flag = new int*[S];  //line# : left label : collinearity coeff : right label : collinearity coeff
//	for (int i = 0; i < S; i++) {
//		flag[i] = new int[5];
//		for (int j = 0; j < 5; j++)
//			flag[i][j] = -1;
//	}
//
//	std::vector<std::vector<int>> lineSet;
//	int count = 1;
//	for (int i = 0; i < S; i++) {
//		if (flag[i][3] > 0) continue; //has right label already...
//		std::vector<int> seq;
//		seq.push_back(count); //store the line no.
//		flag[i][0] = seq[0];
//		count++;
//		seq.push_back(Lines[0][i]); //store the label no.
//		Point ileft = Point(Lines[3][i], Lines[2][i]);
//		Point iright = Point(Lines[5][i], Lines[4][i]);
//		int ref_slope = atan2(iright.y - ileft.y, iright.x - ileft.x) * 180 / CV_PI;
//		if (ref_slope < 0) ref_slope = 180 + ref_slope;
//		int f = 1, ref = i;
//		std::cout << "\n label " << Lines[0][ref];
//
//		do {
//			int next_right = -1;
//			int roi_h1 = 0, roi_h2 = ileft.y + 1, roi_w1 = iright.x, roi_w2 = c; //select roi box...
//			std::vector<int> labels_in_roi;
//			for (int k = 0; k < S; k++) {
//				if (k == ref) continue;
//				if (roi_h1 <= Lines[2][k] && Lines[2][k] < roi_h2 && roi_w1 <= Lines[3][k] && Lines[3][k] < roi_w2)
//					labels_in_roi.push_back(k);
//			}
//			if (int(labels_in_roi.size()) == 0) {
//				std::cout << "\nNo right nbrs...";
//				f = 0;
//				break;
//			}
//
//			int *coll = new int[int(labels_in_roi.size())];
//			for (int k = 0; k < int(labels_in_roi.size()); k++) {
//				int loc = labels_in_roi[k];
//				Point p[4];
//				p[0] = ileft;
//				p[1] = iright;
//				p[2] = Point(Lines[3][loc], Lines[2][loc]);
//				p[3] = Point(Lines[5][loc], Lines[4][loc]);
//				coll[k] = corrCoeff(p);
//			}
//			int minCollinearity = coll[0], minlbl = labels_in_roi[0];
//			for (int k = 1; k<int(labels_in_roi.size()); k++) {
//				if (coll[k] < minCollinearity) {
//					minCollinearity = coll[k];
//					minlbl = labels_in_roi[k];
//				}
//			}
//			std::cout << " --- next label = " << Lines[0][minlbl];
//			if (flag[minlbl][1] == -1 && flag[minlbl][0] == -1) { //no left nbr to this...
//				seq.push_back(Lines[0][minlbl]);
//				flag[i][3] = minlbl;
//				flag[i][4] = minCollinearity;
//				flag[minlbl][1] = i;
//				flag[minlbl][2] = minCollinearity;
//			}
//			else if (flag[minlbl][1] == -1 && flag[minlbl][0] > -1) { //no left nbr to this...
//				seq[0] = flag[minlbl][0];
//				seq.push_back(Lines[0][minlbl]);
//				for (int k = int(seq.size()) - 1; k > 1; k--)
//					seq[k] = seq[k - 1];
//				seq[1] = Lines[0][minlbl];
//				flag[i][3] = minlbl;
//				flag[i][4] = minCollinearity;
//				flag[minlbl][1] = i;
//				flag[minlbl][2] = minCollinearity;
//				f = 0;
//			}
//			else
//				f = 0;
//
//		} while (f);
//
//		if (seq[0]<int(lineSet.size())) //inserted in previous line...
//			lineSet[seq[0] - 1] = seq; //new line...
//		else
//			lineSet.push_back(seq);
//	}
//	for (int i = 0; i<int(lineSet.size()); i++) {
//		std::cout << "\nLine no. ";
//		for (int j = 0; j<int(lineSet[i].size()); j++) {
//			std::cout << lineSet[i][j] << " -> ";
//		}
//		std::cout << "end";
//	}
//}
//
//std::vector<std::vector<Point>> extendLines(std::vector<int> Lines[6], int S, int theta) {
//	//Sort the list...
//	int max0 = -1, max1 = -1;
//	for (int i = 0; i < S - 1; i++) {
//		for (int j = 0; j < S - i - 1; j++) {
//			if (Lines[1][j] < Lines[1][j + 1] || (Lines[1][j] == Lines[1][j + 1] && Lines[0][j] > Lines[0][j + 1])) {
//				int *temp = new int[6];
//				for (int k = 0; k < 6; k++) {
//					temp[k] = Lines[k][j];
//					Lines[k][j] = Lines[k][j + 1];
//					Lines[k][j + 1] = temp[k];
//				}
//			}
//		}
//	}
//	int **flag = new int*[S];  //line# : left label : collinearity coeff : right label : collinearity coeff
//	for (int i = 0; i < S; i++) {
//		flag[i] = new int[5];
//		for (int j = 0; j < 5; j++)
//			flag[i][j] = -1;
//	}
//
//	std::vector<std::vector<Point>> lineSet;
//	int sz = S, count = 0;
//	for (int i = 0; i < S; i++) {
//		Point left = Point(Lines[3][i], Lines[2][i]);
//		Point right = Point(Lines[5][i], Lines[4][i]);
//
//		//extend the line and find distance from every other label...
//		int theta = atan2(right.y - left.y, right.x - left.x) * 180 / CV_PI;
//		double **distance = new double*[sz];
//		for (int k = 0; k < sz; k++) {
//			distance[k] = new double[2];
//			distance[k][0] = 0; //label location
//			distance[k][1] = -1; //distance
//		}
//		for (int k = 0; k < S; k++) {
//			if (i == k) {
//				distance[k][0] = k;
//				continue;
//			}
//			if (flag[k][0] != -1) {
//				distance[k][0] = k;
//				distance[k][1] = 999999;
//				continue;
//			}
//
//			double dist_left, dist_right;
//			double m = tan(theta*CV_PI / 180);
//			dist_left = (abs(m*(Lines[3][i] - Lines[3][k]) - (Lines[2][i] - Lines[2][k]))) / sqrt(1 + m*m);
//			dist_right = (abs(m*(Lines[5][i] - Lines[5][k]) - (Lines[4][i] - Lines[4][k]))) / sqrt(1 + m*m);
//			distance[k][0] = k;
//			distance[k][1] = dist_left <= dist_right ? dist_left : dist_right;
//		}
//
//		//arrange distance in increasing order...
//		for (int k = 0; k < sz - 1; k++) {
//			for (int j = 0; j < sz - k - 1; j++) {
//				if (distance[j][1] > distance[j + 1][1]) {
//					double *temp = new double[2];
//					temp = distance[j];
//					distance[j] = distance[j + 1];
//					distance[j + 1] = temp;
//				}
//			}
//		}
//		std::cout << "\n\nDistance of line through label " << Lines[0][i] << " : ";
//		double diff_left, diff_right;
//		for (int k = 1; k < sz; k++) {
//			std::cout << "\nfrom label " << Lines[0][int(distance[k][0])] << " = " << distance[k][1];
//		}
//
//		//assign line nos....
//		int this_line = 1;
//		for (int k = 1; k < sz - 2; k++) {
//			/*diff_left = abs(distance[k][1] - distance[k - 1][1]);
//			diff_right = abs(distance[k + 1][1] - distance[k][1]);
//			*/
//			if (distance[k][1] > 70) {
//				this_line = k - 1;
//				break;
//			}
//		}
//		if (this_line >= 2) {
//			count++;
//			for (int k = 0; k <= this_line; k++) {
//				flag[int(distance[k][0])][0] = count;
//			}
//		}
//	}
//
//	//arrange increasing order of line no...
//	for (int i = 0; i < sz - 1; i++) {
//		for (int j = 0; j < sz - i - 1; j++) {
//			if (flag[j][0]> flag[j + 1][0] || (flag[j][0] == flag[j + 1][0] && Lines[0][j]>Lines[0][j + 1])) {
//				int *temp1 = new int[6];
//				for (int k = 0; k < 6; k++) {
//					temp1[k] = Lines[k][j];
//					Lines[k][j] = Lines[k][j + 1];
//					Lines[k][j + 1] = temp1[k];
//				}
//
//				int *temp2 = new int[5];
//				temp2 = flag[j];
//				flag[j] = flag[j + 1];
//				flag[j + 1] = temp2;
//			}
//		}
//	}
//	std::cout << "\nNo.of lines = " << count << "\n\nAssigned line nos.::";
//	std::vector<Point> seq;
//	for (int i = 0; i < S; i++) {
//		//std::cout << "\nLabel " << Lines[0][i] << " at " << Point(Lines[3][i], Lines[2][i]) << " is in line no. " << flag[i][0];
//		if (flag[i][0] < 1) continue;
//		if (i == 0 || (i > 0 && flag[i][0] != flag[i - 1][0])) {
//			if (int(seq.size()) > 0) {
//				lineSet.push_back(seq);
//				seq.clear();
//			}
//			//new line...
//			std::cout << "\n starting line " << flag[i][0];
//			seq.push_back(Point(flag[i][0], 0));
//		}
//		seq.push_back(Point(Lines[3][i], Lines[2][i]));
//		seq.push_back(Point(Lines[5][i], Lines[4][i]));
//	}
//	lineSet.push_back(seq);
//
//	for (int i = 0; i < int(lineSet.size()); i++) {
//		std::cout << "\n\nLine " << lineSet[i][0].x;
//		for (int k = 1; k < int(lineSet[i].size()) - 1; k += 2) {
//			for (int l = 1; l < int(lineSet[i].size()) - k - 2; l += 2) {
//				if (lineSet[i][l].x > lineSet[i][l + 2].x) {
//					Point temp;
//					temp = lineSet[i][l];
//					lineSet[i][l] = lineSet[i][l + 2];
//					lineSet[i][l + 2] = temp;
//
//					temp = lineSet[i][l + 1];
//					lineSet[i][l + 1] = lineSet[i][l + 3];
//					lineSet[i][l + 3] = temp;
//				}
//			}
//		}
//		for (int j = 1; j<int(lineSet[i].size()); j++) {
//
//			std::cout << "-->" << lineSet[i][j];
//		}
//	}
//	return(lineSet);
//}
//
//void initiatePath(int **visited, int **thin_roi, int h, int w, Point startVertex) {
//	Point adjVertex;
//	int window = 5;
//	int **nbd = new int*[window];
//	for (int i = 0; i < window; i++) {
//		nbd[i] = new int[window];
//		for (int j = 0; j < window; j++) {
//			if (startVertex.y - window + i < 0 || startVertex.y - window + i >= h || startVertex.x - window + j < 0 || startVertex.x - window + j >= w)
//				nbd[i][j] = 0;
//			else
//				nbd[i][j] = (255 - thin_roi[startVertex.y - window + i][startVertex.x - window + j]) / 255;
//			if (visited[startVertex.y - window + i][startVertex.x - window + j] == -2) //got a terminal point....
//				adjVertex = Point(startVertex.x - window + j, startVertex.y - window + i);
//		}
//	}
//	int sumN = 0, sumS = 0, sumE = 0, sumW = 0;
//	for (int k = 1; k < window - 1; k++) {
//		sumN += nbd[0][k];
//		sumS += nbd[window - 1][k];
//		sumW += nbd[k][0];
//		sumE += nbd[k][window - 1];
//	}
//}
//
//void closeCorners(std::vector<Point> &corners, std::vector<std::vector<Point>> &edgeLabels, int **roi_label) {
//	std::vector<Point> newCorners;
//	int *flag = new int[int(corners.size())];
//	for (int i = 0; i<int(corners.size()); i++)
//		flag[i] = 0;
//	for (int i = 0; i<int(corners.size()) - 1; i++) {
//		for (int j = i + 1; j<int(corners.size()); j++) {
//			Point ci = corners[i];
//			Point cj = corners[j];
//			if (abs(ci.x - cj.x) <= 5 && abs(ci.y - cj.y) <= 5) {
//				Point mean(int(ci.x + cj.x) / 2, int(ci.y + cj.y) / 2);
//				newCorners.push_back(mean);
//				std::cout << "\n" << ci << " and " << cj << " are too close! \t New corner at mean = " << mean;
//				flag[i] = 1;
//				flag[j] = 1;
//				std::vector<Point> newEdge;
//				newEdge.push_back(mean);
//				for (int k = 1; k<int(edgeLabels[i].size()); k++) {
//					if (roi_label[edgeLabels[i][k].y][edgeLabels[i][k].x] != 0) {
//						newEdge.push_back(edgeLabels[i][k]);
//					}
//				}
//				for (int k = 1; k<int(edgeLabels[j].size()); k++) {
//					if (roi_label[edgeLabels[j][k].y][edgeLabels[j][k].x] != 0) {
//						newEdge.push_back(edgeLabels[j][k]);
//					}
//				}
//				edgeLabels.push_back(newEdge);
//
//			}
//		}
//	}
//	for (int i = 0; i<int(newCorners.size()); i++) {
//		corners.push_back(newCorners[i]);
//	}
//	for (int i = 0; i<int(corners.size()); i++) {
//		if (flag[i] == 1) {
//			for (int k = i + 1; k<int(corners.size()); k++)
//				flag[k - 1] = flag[k];
//			corners.erase(corners.begin() + i);
//			edgeLabels.erase(edgeLabels.begin() + i);
//			i--;
//		}
//	}
//}
//
//void closeCorners(std::vector<Point> &corners, std::vector<Point> &terminals, std::vector<std::vector<Point>> &edgeLabels, int **roi_label) {
//	std::vector<Point> newCorners;
//	int *flag = new int[int(corners.size())];
//	for (int i = 0; i<int(corners.size()); i++)
//		flag[i] = 0;
//	for (int i = 0; i<int(corners.size()) - 1; i++) {
//		if (flag[i] == 1) continue;
//		std::vector<int> nbrs;
//		Point ci = corners[i];
//		nbrs.push_back(i);
//		for (int j = i + 1; j<int(corners.size()); j++) {
//			Point cj = corners[j];
//			if (abs(ci.x - cj.x) <= 7 && abs(ci.y - cj.y) <= 7) {
//				std::cout << "\nFound a nbr!" << i << " and " << j;
//				flag[j] = 1;
//				nbrs.push_back(j);
//			}
//		}
//		if (int(nbrs.size()) > 1) {
//			flag[i] = 1;
//			std::cout << "\nThe points ";
//			for (int j = 1; j<int(nbrs.size()); j++) {
//				std::cout << corners[nbrs[j]] << " , ";
//			}
//			std::cout << "\n are very close to " << ci;
//			Point mean = ci;
//			for (int j = 1; j<int(nbrs.size()); j++) {
//				Point cj = corners[nbrs[j]];
//				mean += cj;
//
//			}
//			mean = mean / int(nbrs.size());
//			newCorners.push_back(mean);
//			std::vector<Point> newEdge;
//			newEdge.push_back(mean);
//			for (int j = 0; j<int(nbrs.size()); j++) {
//				for (int k = 1; k<int(edgeLabels[nbrs[j]].size()); k++) {
//					if (roi_label[edgeLabels[nbrs[j]][k].y][edgeLabels[nbrs[j]][k].x] != 0) {
//						newEdge.push_back(edgeLabels[nbrs[j]][k]);
//					}
//				}
//			}
//			edgeLabels.push_back(newEdge);
//		}
//	}
//	for (int i = 0; i<int(newCorners.size()); i++) {
//		corners.push_back(newCorners[i]);
//	}
//	for (int i = 0; i<int(corners.size()); i++) {
//		if (flag[i] == 1) {
//			for (int k = i + 1; k<int(corners.size()); k++)
//				flag[k - 1] = flag[k];
//			corners.erase(corners.begin() + i);
//			edgeLabels.erase(edgeLabels.begin() + i);
//			i--;
//		}
//	}
//}
//
//std::vector<std::vector<int>> findTurns(int **roi_label, int &roi_labelNum, int h, int w, int **incidence, std::vector<int> &edges, std::vector<Point> &vertex, int E, int V) {
//	std::vector<std::vector<int>> IncMatrix;
//	for (int i = 0; i < V; i++) {
//		std::vector<int> vi_row;
//		for (int j = 0; j < E; j++) {
//			vi_row.push_back(incidence[i][j]);
//		}
//		IncMatrix.push_back(vi_row);
//	}
//
//	if (E > 1) {
//		for (int j = 0; j < E; j++) {
//			int lbl = edges[j];
//			std::cout << "\nChecking edge " << edges[j] << "...";
//			int vi = -1, vj = -1;
//			for (int i = 0; i < V; i++) {
//				if (incidence[i][j] && vi == -1)
//					vi = i;
//				else if (incidence[i][j] && vj == -1)
//					vj = i;
//			}
//			if (vj == -1) //connected to a loop...
//				continue;
//
//			int flag = 0;
//			double dist = 0, maxDist = -1;
//			Point farthestPoint, A = vertex[vi], B = vertex[vj];
//			double slope = tan(atan2(B.y - A.y, B.x - A.x));
//			for (int k = 0; k < h; k++) {
//				for (int l = 0; l < w; l++) {
//					if (roi_label[k][l] != lbl) continue;
//					dist = abs(slope*(l - A.x) - (k - A.y)) / sqrt(1 + slope * slope);
//					if (dist > maxDist) {
//						maxDist = dist;
//						farthestPoint = Point(l, k);
//					}
//				}
//			}
//			std::cout << "max distance = " << maxDist << " at " << farthestPoint;
//			if (maxDist > 25) {
//				flag = 1;
//				std::cout << "\nEdge " << edges[j] << " has a turn...at " << farthestPoint;
//				vertex.push_back(farthestPoint);
//				//V++;
//				//edges.push_back(lbl);
//				//new label...
//				double m1 = tan(atan2(farthestPoint.y - A.y, farthestPoint.x - A.x));
//				double m2 = tan(atan2(farthestPoint.y - B.y, farthestPoint.x - B.x));
//				double d1, d2;
//				for (int k = 0; k < h; k++) {
//					for (int l = 0; l < w; l++) {
//						if (roi_label[k][l] != lbl) continue;
//						d1 = abs(m1*(l - A.x) - (k - A.y)) / sqrt(1 + m1 * m1);
//						d2 = abs(m2*(l - B.x) - (k - B.y)) / sqrt(1 + m2 * m2);
//						if (d1 > d2) {
//							roi_label[k][l] = roi_labelNum + 1;
//						}
//					}
//				}
//				edges.push_back(roi_labelNum + 1);
//				roi_labelNum++;
//				//E++;
//				IncMatrix[vj][j] = 0;
//				for (int i = 0; i < V; i++) {
//					if (i == vj)
//						IncMatrix[i].push_back(1);
//					else
//						IncMatrix[i].push_back(0);
//				}
//				std::vector<int> new_row;
//				for (int jj = 0; jj < E; jj++) {
//					if (jj == j)
//						new_row.push_back(1);
//					else
//						new_row.push_back(0);
//				}
//				new_row.push_back(1);
//				IncMatrix.push_back(new_row);
//			}
//		}
//	}
//	else if (E == 1) {
//		std::cout << "\nChecking edge 1...";
//		int vi = 0, vj = 1;
//		double dist = 0, maxDist = -1;
//		Point farthestPoint, A = vertex[vi], B = vertex[vj];
//		double slope = tan(atan2(B.y - A.y, B.x - A.x));
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				if (roi_label[k][l] != 0) continue;
//				dist = abs(slope*(l - A.x) - (k - A.y)) / sqrt(1 + slope * slope);
//				if (dist > maxDist) {
//					maxDist = dist;
//					farthestPoint = Point(l, k);
//				}
//			}
//		}
//		std::cout << "max distance = " << maxDist << " at " << farthestPoint;
//		if (maxDist > 25) {
//			std::cout << "\nEdge 1 has a turn...at " << farthestPoint;
//			vertex.push_back(farthestPoint);
//			//V++;
//			//edges.push_back(2);
//			//new label...
//			double m1 = tan(atan2(farthestPoint.y - A.y, farthestPoint.x - A.x));
//			double m2 = tan(atan2(farthestPoint.y - B.y, farthestPoint.x - B.x));
//			double d1, d2;
//			for (int k = 0; k < h; k++) {
//				for (int l = 0; l < w; l++) {
//					if (roi_label[k][l] != 0) continue;
//					d1 = abs(m1*(l - A.x) - (k - A.y)) / sqrt(1 + m1 * m1);
//					d2 = abs(m2*(l - B.x) - (k - B.y)) / sqrt(1 + m2 * m2);
//					if (d1 > d2) {
//						roi_label[k][l] = 2;
//					}
//				}
//			}
//			edges.push_back(2);
//			roi_labelNum++;
//			//E++;
//			IncMatrix[vj][0] = 0;
//			IncMatrix[vi].push_back(0);
//			IncMatrix[vj].push_back(1);
//			std::vector<int> new_row;
//			new_row.push_back(1);
//			new_row.push_back(1);
//			IncMatrix.push_back(new_row);
//
//
//		}
//	}
//
//	std::cout << "\nNew vertex set : ";
//	for (int i = 0; i<int(vertex.size()); i++)
//		std::cout << vertex[i] << "    ";
//
//	std::cout << "\nNew edge set : ";
//	for (int i = 0; i<int(edges.size()); i++)
//		std::cout << edges[i] << "    ";
//	return(IncMatrix);
//}
//
//int countLabel(int **roi_label, int lbl, int h, int w) {
//	int count = 0;
//	for (int i = 0; i < h; i++) {
//		for (int j = 0; j < w; j++)
//			if (roi_label[i][j] == lbl)
//				count++;
//	}
//	return(count);
//}
//
//void removeParallelEdges(int **roi_label, int roi_labelNum, int h, int w, int **incidence, std::vector<int> &edges, std::vector<Point> vertex, int &E, int V) {
//	for (int i = 0; i < E - 1; i++) {
//		for (int k = i + 1; k < E; k++) {
//			int sum = 0;
//			for (int j = 0; j < V; j++) {
//				sum += abs(incidence[j][i] - incidence[j][k]);
//			}
//			int del_edge;
//			if (sum == 0) { //columns are same => parallel edges...
//				if (countLabel(roi_label, edges[i], h, w) > countLabel(roi_label, edges[k], h, w)) {
//					std::cout << "\nDelete e = " << edges[i];
//					del_edge = i;
//				}
//				else {
//					std::cout << "\nDelete e = " << edges[k];
//					del_edge = k;
//				}
//				//deleting column = del_edge...
//				for (int kk = del_edge + 1; kk < E; kk++) {
//					for (int j = 0; j < V; j++)
//						incidence[j][kk - 1] = incidence[j][kk];
//					edges[kk - 1] = edges[kk];
//				}
//				E--;
//				i--;
//			}
//		}
//	}
//}
//
//void findAdjacency(int **thin_roi, int h, int w, std::vector<Point> corners, std::vector<Point> terminals, std::vector<std::vector<int>> &adjacency, int &V, int &E, std::vector<Point> &vertex) {
//	int **roi_copy = thin_roi;
//	std::vector<std::vector<Point>> edgeLabels;
//	for (int i = 0; i<int(corners.size()); i++) {
//		Point ref = corners[i];
//		std::vector<Point> edges;
//		edges.push_back(ref);
//		for (int k = ref.y - 2; k <= ref.y + 2; k++) {
//			for (int l = ref.x - 2; l <= ref.x + 2; l++) {
//				roi_copy[k][l] = 255;
//			}
//		}
//		for (int k = -3; k <= 3; k++) {
//			if (roi_copy[ref.y + k][ref.x - 3] == 0)
//				edges.push_back(Point(ref.x - 3, ref.y + k));
//			if (roi_copy[ref.y + k][ref.x + 3] == 0)
//				edges.push_back(Point(ref.x + 3, ref.y + k));
//		}
//		for (int k = -2; k <= 2; k++) {
//			if (roi_copy[ref.y - 3][ref.x + k] == 0)
//				edges.push_back(Point(ref.x + k, ref.y - 3));
//			if (roi_copy[ref.y + 3][ref.x + k] == 0)
//				edges.push_back(Point(ref.x + k, ref.y + 3));
//		}
//		edgeLabels.push_back(edges);
//		std::cout << "\nSize=" << int(edges.size());
//		edges.clear();
//	}
//
//
//	//int V, E;
//	//Point *vertex;
//	std::vector<int> edges;
//
//	//Edge Set...
//	std::vector<std::vector<int>> edgeIncidence;
//	int roi_labelNum = 0;
//	int **roi_box;
//	int **roi_label = new int*[h];
//	for (int i = 0; i < h; i++) {
//		roi_label[i] = new int[w];
//		for (int j = 0; j < w; j++)
//			roi_label[i][j] = 0;
//	}
//	if (int(corners.size()) > 0) {  //If graph has more than 2 vertices...
//		roi_labelNum = labelling(roi_copy, roi_label, h, w);
//		roi_box = new int*[roi_labelNum];
//		for (int i = 0; i < roi_labelNum; i++)
//			roi_box[i] = new int[6];
//		boundingRectangles(roi_box, roi_label, h, w, roi_labelNum);
//		std::vector<Point> boundaryPoints;
//		int count = 0;
//		//std::cout << "\nBoundary Points : ";
//		for (int i = 1; i < roi_labelNum; i++) {
//			if (roi_box[i][2] == -1) continue;
//			count++;
//			for (int k = roi_box[i][4]; k <= roi_box[i][5]; k++) {
//				if (roi_copy[k][roi_box[i][2]] == 0) {
//					boundaryPoints.push_back(Point(roi_box[i][2], k));
//					//std::cout << Point(roi_box[i][2], k) << " , ";
//				}
//			}
//		}
//		//std::cout << "\nNo.of edges = " << count;
//		//closeCorners(corners, edgeLabels, roi_label);
//		closeCorners(corners, terminals, edgeLabels, roi_label);
//
//		//Vertex Set...
//		V = int(terminals.size()) + int(edgeLabels.size());
//		for (int k = 0; k < int(terminals.size()); k++) {
//			vertex.push_back(terminals[k]);
//		}
//		for (int k = 0; k < int(edgeLabels.size()); k++) {
//			vertex.push_back(edgeLabels[k][0]);
//		}
//
//		std::cout << "\nList of vertex: ";
//		for (int k = 0; k < V; k++)
//			std::cout << vertex[k] << " \t ";
//
//		for (int i = 0; i<int(terminals.size()); i++) {
//			std::vector<int>  pendantEdge;
//			pendantEdge.push_back(roi_label[terminals[i].y][terminals[i].x]);
//			edgeIncidence.push_back(pendantEdge);
//		}
//		for (int i = 0; i<int(terminals.size()); i++) {
//			std::cout << "\nList of labels for terminal at " << vertex[i] << " : ";
//			std::cout << edgeIncidence[i][0];
//		}
//
//		for (int i = 0; i<int(corners.size()); i++) {
//			std::vector<int>  cornerEdge;
//			for (int j = 1; j<int(edgeLabels[i].size()); j++) {
//				cornerEdge.push_back(roi_label[edgeLabels[i][j].y][edgeLabels[i][j].x]);
//			}
//			edgeIncidence.push_back(cornerEdge);
//		}
//		for (int i = 0; i<int(corners.size()); i++) {
//			std::cout << "\nList of labels for corner at " << vertex[i + int(terminals.size())] << " : ";
//			for (int j = 0; j<int(edgeIncidence[i + int(terminals.size())].size()); j++)
//				std::cout << edgeIncidence[i + int(terminals.size())][j] << " , ";
//		}
//
//		for (int i = 0; i<int(edgeLabels.size()); i++) {
//			for (int j = 1; j<int(edgeLabels[i].size()); j++) {
//				int flag = 0;
//				for (int k = 0; k<int(edges.size()); k++) {
//					if (roi_label[edgeLabels[i][j].y][edgeLabels[i][j].x] == edges[k]) {
//						flag = 1;
//						break;
//					}
//				}
//				if (flag == 0)
//					edges.push_back(roi_label[edgeLabels[i][j].y][edgeLabels[i][j].x]);
//			}
//		}
//		E = int(edges.size());
//		std::cout << "\nNo. of edges = " << int(edges.size());
//		std::cout << "\nList of edges: ";
//		for (int i = 0; i<int(edges.size()); i++)
//			std::cout << edges[i] << " \t ";
//	}
//	else {  //if graph has only 2 vertices...
//		E = 1;
//		V = 2;
//		vertex.push_back(terminals[0]);
//		vertex.push_back(terminals[1]);
//		edges.push_back(1);
//		edges.push_back(1);
//		std::cout << "\nNo.of edges = 1";
//	}
//
//	std::cout << "\nSize of incidence matrix = " << V << "x" << E;
//
//	int **incidence = new int*[V];
//	for (int i = 0; i < V; i++) {
//		incidence[i] = new int[E];
//		for (int j = 0; j < E; j++)
//			incidence[i][j] = 0;
//	}
//
//	if (E == 1) {
//		incidence[0][0] = 1;
//		incidence[1][0] = 1;
//	}
//	else {
//		for (int i = 0; i < E; i++) { //gives the edge/column no...
//			for (int k = 0; k<int(edgeIncidence.size()); k++) { //gives the vertex/row no...
//				for (int l = 0; l<int(edgeIncidence[k].size()); l++) {
//					if (edgeIncidence[k][l] == edges[i]) {
//						incidence[k][i] = 1;
//						break;
//					}
//				}
//			}
//		}
//	}
//	std::cout << "\nIncidence Matrix:\n\t";
//	for (int j = 0; j < E; j++)
//		std::cout << "\t" << edges[j];
//
//	for (int i = 0; i < V; i++) {
//		std::cout << "\n" << vertex[i];
//		for (int j = 0; j < E; j++)
//			std::cout << "\t" << incidence[i][j];
//	}
//
//	//edges[j] -> has label no. of each edge (int)
//	//vertex[i] -> has vertex pixel location (Point)
//	//incidence[i][j] -> has incidence relation (int-0/1)
//
//	//remove parallel edges...
//	removeParallelEdges(roi_label, roi_labelNum, h, w, incidence, edges, vertex, E, V);
//
//	std::cout << "\nIncidence Matrix after removing parallel edges, if any :\n\t";
//	for (int j = 0; j < E; j++)
//		std::cout << "\t" << edges[j];
//
//	for (int i = 0; i < V; i++) {
//		std::cout << "\n" << vertex[i];
//		for (int j = 0; j < E; j++)
//			std::cout << "\t" << incidence[i][j];
//	}
//
//	//check for turning points...
//	std::vector<std::vector<int>> IncMatrix;
//	if (E > 1)
//		IncMatrix = findTurns(roi_label, roi_labelNum, h, w, incidence, edges, vertex, E, V);
//	if (E == 1)
//		IncMatrix = findTurns(roi_copy, roi_labelNum, h, w, incidence, edges, vertex, E, V);
//	if (int(IncMatrix.size()) > V) {
//		std::cout << "\n****************New Incidence Matrix:\n\t";
//		for (int i = 0; i < int(IncMatrix.size()); i++) {
//			std::cout << "\n";
//			for (int j = 0; j < int(IncMatrix[i].size()); j++)
//				std::cout << "\t" << IncMatrix[i][j];
//		}
//	}
//
//	int VV = int(IncMatrix.size());
//	int EE = int(IncMatrix[0].size());
//
//	//repeat...
//	if (EE > 1) {
//		int **new_incidence = new int*[VV];
//		for (int i = 0; i < VV; i++) {
//			new_incidence[i] = new int[EE];
//			for (int j = 0; j < EE; j++)
//				new_incidence[i][j] = IncMatrix[i][j];
//		}
//		IncMatrix = findTurns(roi_label, roi_labelNum, h, w, new_incidence, edges, vertex, EE, VV);
//		if (int(IncMatrix.size()) > VV) {
//			std::cout << "\n****************New New Incidence Matrix:\n\t";
//			for (int i = 0; i < int(IncMatrix.size()); i++) {
//				std::cout << "\n";
//				for (int j = 0; j < int(IncMatrix[i].size()); j++)
//					std::cout << "\t" << IncMatrix[i][j];
//			}
//		}
//
//		VV = int(IncMatrix.size());
//		EE = int(IncMatrix[0].size());
//	}
//
//	for (int i = 0; i < VV; i++) {
//		std::vector<int> entry;
//		for (int j = 0; j < VV; j++) {
//			entry.push_back(0);
//		}
//		adjacency.push_back(entry);
//	}
//
//	for (int j = 0; j < EE; j++) {
//		int v1 = -1, v2 = -1;
//		for (int i = 0; i < VV; i++) {
//			if (v1 == -1 && IncMatrix[i][j] == 1)
//				v1 = i;
//			else if (v1 != -1 && v2 == -1 && IncMatrix[i][j] == 1)
//				v2 = i;
//		}
//		if (v1 == -1 || v2 == -1)
//			continue;
//		else {
//			adjacency[v1][v2] = 1;
//			adjacency[v2][v1] = 1;
//		}
//	}
//	V = VV;
//	E = EE;
//}
//
//void DFS(std::vector< pair<int, double> > *graph, int src, double prev_len, double *max_len, vector <bool> &visited)
//{
//	// Mark the src node visited
//	visited[src] = 1;
//
//	// curr_len is for distance from src point to its adjacent point
//	double curr_len = 0;
//
//	// Adjacent is pair type which stores destination point and distance
//	pair < int, double > adjacent;
//
//	// Traverse all adjacent
//	for (int i = 0; i<int(graph[src].size()); i++)
//	{
//		// Adjacent element
//		adjacent = graph[src][i];
//
//		// If point is not visited
//		if (!visited[adjacent.first])
//		{
//			// Total distance from src point to its adjacent
//			curr_len = prev_len + adjacent.second;
//
//			// Call DFS for adjacent city
//			DFS(graph, adjacent.first, curr_len, max_len, visited);
//		}
//
//		// If total cable length till now greater than previous length then update it
//		if ((*max_len) < curr_len)
//			*max_len = curr_len;
//
//		// make curr_len = 0 for next adjacent
//		curr_len = 0;
//	}
//}
//
//// n is number of cities or nodes in graph
//// cable_lines is total cable_lines among the cities
//// or edges in graph
//double longestPath(std::vector<pair<int, double> > *graph, int V)
//{
//	// maximum length of cable among the connected cities
//	double max_len = DBL_MIN;
//
//	// call DFS for each city to find maximum length of cable
//	for (int i = 1; i < V; i++)
//	{
//		// initialize visited array with 0
//		vector< bool > visited(V, false);
//
//		// Call DFS for src vertex i
//		DFS(graph, i, 0, &max_len, visited);
//	}
//
//	return max_len;
//}
//
//void dfs(std::vector<pair<int, double>> *graph, int vi, int count, double dist, vector<bool> &visited, int &q, double &m) {
//	visited[vi] = true;
//	for (int vj = 0; vj<int(graph[vi].size()); vj++) {
//		if (!visited[graph[vi][vj].first]) {
//			count++;
//			dist += graph[vi][vj].second;
//			dfs(graph, graph[vi][vj].first, count, dist, visited, q, m);
//			count--;
//			dist -= graph[vi][vj].second;
//		}
//	}
//	if (dist > m) {
//		m = dist;
//		q = vi;
//	}
//}
//
//void majorOrientation(std::vector<std::vector<int>> &adjacency, std::vector<Point> vertex, int V, int E) {
//	double theta_min = 0.0;
//	double numerator = 0, denominator = 0;
//	int **vec;
//	int *deg;
//	double *lengths;
//	vec = new int*[E];
//	for (int i = 0; i < E; i++) {
//		vec[i] = new int[2];
//	}
//	deg = new int[E];
//	lengths = new double[E];
//	int e = 0;
//	for (int i = 0; i < V; i++) {
//		for (int j = i + 1; j < V; j++) {
//			if (adjacency[i][j] == 1) {
//				Point p1, p2;
//				p1 = vertex[i].x <= vertex[j].x ? vertex[i] : vertex[j];
//				p2 = vertex[i].x > vertex[j].x ? vertex[i] : vertex[j];
//				int angle = atan2(p2.y - p1.y, p2.x - p1.x) * 180 / CV_PI;
//				/*if (angle < 0)
//				angle += 180;*/
//				numerator += norm(p2 - p1)*sin(angle*CV_PI / 180);
//				denominator += norm(p2 - p1)*cos(angle*CV_PI / 180);
//				//std::cout << "\norientation of vector from " << p1 << " to " << p2 << " = " << angle << " degrees and length = " << norm(p2 - p1);
//				vec[e][0] = vertex[i].x <= vertex[j].x ? i : j;
//				vec[e][1] = vertex[i].x > vertex[j].x ? i : j;
//				deg[e] = angle;
//				lengths[e] = norm(p2 - p1);
//				e++;
//			}
//		}
//	}
//	//sort in descending order of length...
//	for (int i = 0; i < E - 1; i++) {
//		for (int j = 0; j < E - i - 1; j++) {
//			if (lengths[j] < lengths[j + 1]) {
//				double temp = lengths[j];
//				lengths[j] = lengths[j + 1];
//				lengths[j + 1] = temp;
//
//				int temp_angle = deg[j];
//				deg[j] = deg[j + 1];
//				deg[j + 1] = temp_angle;
//
//				int *temp_vector = new int[2];
//				temp_vector = vec[j];
//				vec[j] = vec[j + 1];
//				vec[j + 1] = temp_vector;
//			}
//		}
//	}
//
//	double major_length = 0, minor_length = 0;
//	int major_n = 0, minor_n = 0;
//	int major_angle = 0, minor_angle = 0;
//	std::cout << "\nThe orientations in descending order of length...";
//	for (int i = 0; i < E; i++) {
//		std::cout << "\norientation of vector from " << vertex[vec[i][0]] << " to " << vertex[vec[i][1]] << " = " << deg[i] << " degrees and length = " << lengths[i];
//		std::cout << "\t deviated by " << abs(deg[0] - deg[i]) << " angles";
//		if (abs(deg[0] - deg[i]) <= 30) {
//			std::cout << "\nThis is a major edge...";
//			major_length += lengths[i];
//			major_angle += deg[i];
//			major_n++;
//		}
//		else if (abs(deg[0] - deg[i]) <= 50) {
//			std::cout << "\nThis is an ambiguous edge...";
//			//minor_length += lengths[i];
//			//minor_angle += deg[i];
//			//minor_n++;
//			adjacency[vec[i][0]][vec[i][1]] = 2;
//			adjacency[vec[i][1]][vec[i][0]] = 2;
//		}
//		else {
//			std::cout << "\nThis is a minor edge...";
//			minor_length += lengths[i];
//			minor_angle += deg[i];
//			minor_n++;
//			adjacency[vec[i][0]][vec[i][1]] = -1;
//			adjacency[vec[i][1]][vec[i][0]] = -1;
//		}
//	}
//	if (major_length > minor_length)
//		std::cout << "\nThe major orientation is " << major_angle / major_n << " degrees and the minor orientation is " << minor_angle / major_n;
//	else {
//		std::cout << "\nPossible orientations are " << major_angle / major_n << " degrees and " << minor_angle / major_n << " degrees";
//		for (int i = 0; i < E; i++) {
//			adjacency[vec[i][0]][vec[i][1]] = 2;
//			adjacency[vec[i][1]][vec[i][0]] = 2;
//		}
//	}
//}
//
//void trace(Mat img) {
//	Mat box_im, thin_im;
//	int **gaussimg, **gaussthin, **label, **box;
//
//	gaussimg = new int*[r];
//	gaussthin = new int*[r];
//	label = new int*[r];
//	for (int i = 0; i < r; i++) {
//		gaussimg[i] = new int[c];
//		gaussthin[i] = new int[c];
//		label[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			gaussimg[i][j] = 255;
//			gaussthin[i][j] = 255;
//			label[i][j] = 0;
//		}
//	}
//	mat2arr(img, gaussimg);
//	box_im.create(r, c, CV_8UC1);
//	cvtColor(img, box_im, CV_GRAY2BGR);
//	thin_im.create(r, c, CV_8UC1);
//	arr2mat(gaussthin, thin_im);
//	cvtColor(thin_im, thin_im, CV_GRAY2BGR);
//	//cvtColor(img, thin_im, CV_GRAY2BGR);
//
//	labelNum = labelling(gaussimg, label, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//		std::cout << "\nLabel " << i << " -> " << Point(box[i][2], box[i][4]);
//
//		/*	line(thin_im, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//		line(thin_im, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//		line(thin_im, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//		line(thin_im, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);*/
//	}
//	imshow("GaussBox.tif", box_im);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/GaussBox.tif", box_im);
//
//	//Find trace lines for each label...
//	std::vector<int> Lines[6];
//	int theta = 137;
//
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//
//		int **roi, **thin_roi;
//		int istart = box[i][4];
//		int jstart = box[i][2];
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//
//		roi = new int*[h];
//		thin_roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			thin_roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				thin_roi[k][l] = 255;
//				if (label[k + istart][l + jstart] == i) {
//					roi[k][l] = 0;
//				}
//				else
//					roi[k][l] = 1;
//			}
//		}
//
//		//Thin roi...
//		thinCC(roi, h, w, thin_roi);
//
//		for (int k = istart; k < h + istart; k++) {
//			for (int l = jstart; l < w + jstart; l++) {
//				if (thin_roi[k - istart][l - jstart] == 0) {
//					gaussthin[k][l] = 0;
//					thin_im.at<Vec3b>(k, l) = Vec3b({ 255,0,255 });
//				}
//			}
//		}
//
//		//Determine the terminal points...
//		std::vector<Point> terminals;
//		int sz;
//		sz = findTerminals(thin_roi, h, w, terminals);
//		std::cout << "\n Label " << i << " has " << sz << " terminal points";
//
//		for (int k = 0; k < sz; k++) {
//			gaussthin[terminals[k].y + istart][terminals[k].x + jstart] = 2;
//		}
//
//		//Determine the corner points...
//		std::vector<Point> corners;
//		int sz_c;
//		sz_c = findCorners(thin_roi, h, w, corners);
//		std::cout << "\n Label " << i << " has " << sz_c << " corner points";
//
//		for (int k = 0; k < sz_c; k++) {
//			gaussthin[corners[k].y + istart][corners[k].x + jstart] = 3;
//		}
//
//		int V, E;
//		std::vector<Point> vertex;
//		std::vector<std::vector<int>> adjacency;
//
//		findAdjacency(thin_roi, h, w, corners, terminals, adjacency, V, E, vertex);
//		std::cout << "\nE = " << E;
//
//		int e = 0;
//		std::cout << "\nAdjacency Matrix:\n\t";
//		for (int i = 0; i < V; i++)
//			std::cout << "   " << vertex[i];
//		for (int i = 0; i < V; i++) {
//			std::cout << "\n" << vertex[i];
//			for (int j = 0; j < V; j++) {
//				std::cout << "\t" << adjacency[i][j] << "   ";
//				if (adjacency[i][j] != 0)
//					e++;
//			}
//		}
//		E = e / 2;
//		std::cout << "\nNo. of edges = " << E;
//		//find orientations of each edge vector...
//		majorOrientation(adjacency, vertex, V, E);
//
//
//		std::cout << "\nAdjacency Matrix:\n\t";
//		for (int i = 0; i < V; i++)
//			std::cout << "   " << vertex[i];
//		for (int i = 0; i < V; i++) {
//			std::cout << "\n" << vertex[i];
//			for (int j = 0; j < V; j++) {
//				std::cout << "\t" << adjacency[i][j] << "   ";
//				if (adjacency[i][j] == 1) {
//					line(thin_im, vertex[i] + Point(jstart, istart), vertex[j] + Point(jstart, istart), Scalar(0, 0, 0), 2.5, 8);
//					circle(thin_im, vertex[i] + Point(jstart, istart), 5, Scalar(0, 0, 255), 2);
//					circle(thin_im, vertex[j] + Point(jstart, istart), 5, Scalar(0, 0, 255), 2);
//				}
//				else if (adjacency[i][j] == -1) {
//					line(thin_im, vertex[i] + Point(jstart, istart), vertex[j] + Point(jstart, istart), Scalar(175, 175, 175), 2.5, 8);
//					circle(thin_im, vertex[i] + Point(jstart, istart), 5, Scalar(0, 0, 255), 2);
//					circle(thin_im, vertex[j] + Point(jstart, istart), 5, Scalar(0, 0, 255), 2);
//				}
//				else if (adjacency[i][j] == 2) {
//					line(thin_im, vertex[i] + Point(jstart, istart), vertex[j] + Point(jstart, istart), Scalar(0, 255, 120), 2.5, 8);
//					circle(thin_im, vertex[i] + Point(jstart, istart), 5, Scalar(0, 0, 255), 2);
//					circle(thin_im, vertex[j] + Point(jstart, istart), 5, Scalar(0, 0, 255), 2);
//				}
//			}
//		}
//
//		//imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/GaussThin.tif", thin_im);
//
//
//		//std::vector<pair<int, double>> *graph = new std::vector<pair<int, double>>[V];
//		//for (int i = 0; i < V; i++) {
//		//	std::vector<pair<int, double>> adj;
//		//	for (int j = 0; j < V; j++) {
//		//		if (adjacency[i][j] == 1) {
//		//			adj.push_back(make_pair(j, norm(vertex[i] - vertex[j])));
//		//		}
//		//	}
//		//	graph[i] = adj;
//		//}
//
//		//std::cout << "\nAdjacency List:";
//		//for (int i = 0; i < V; i++) {
//		//	std::cout << "\nVertex " << vertex[i] << " -> ";
//		//	for (int j = 0; j < int(graph[i].size()); j++) {
//		//		std::cout << vertex[graph[i][j].first] << "  ,  " << graph[i][j].second << "\t ;";
//		//	}
//		//}
//
//		////std::cout << "Maximum length of path = " << longestPath(graph, V);
//
//		//std::vector<bool> visited;
//		//for (int i = 0; i < V; i++)
//		//	visited.push_back(false);
//		//int q = 0;
//		//double m = -1;
//		//dfs(graph, 0, 0, 0, visited, q, m);
//		//std::cout << "\nThe longest path from " << vertex[0] + Point(jstart, istart) << " is of length " << m << " to the vertex " << vertex[q] + Point(jstart, istart);
//		////reset...
//		//m = -1;
//		//int old_q = q;
//		//for (int i = 0; i < V; i++)
//		//	visited[i] = false;
//		//dfs(graph, old_q, 0, 0, visited, q, m);
//		//std::cout << "\nThe longest path from " << vertex[old_q] + Point(jstart, istart) << " is of length " << m << " to the vertex " << vertex[q] + Point(jstart, istart);
//
//		/*int **visited = new int*[h];
//		for (int k = 0; k < h; k++) {
//		visited[k] = new int[w];
//		for (int l = 0; l < w; l++) {
//		if (gaussthin[k + istart][l + istart] > 1)
//		visited[k][l] = -2;
//		else if (gaussthin[k + istart][l + istart] == 0)
//		visited[k][l] = -1;
//		else
//		visited[k][l] = 0;
//		}
//		}*/
//
//
//		//delete...
//		for (int k = 0; k < h; k++) {
//			delete[] roi[k];
//			delete[] thin_roi[k];
//		}
//		delete[] roi;
//		delete[] thin_roi;
//
//	}
//
//	//	arr2mat(gaussthin, thin_im);
//	//	cvtColor(thin_im, thin_im, CV_GRAY2BGR);
//	/*for (int i = 0; i < r; i++) {
//	for (int j = 0; j < c; j++) {
//	if (gaussthin[i][j] == 2) {
//	thin_im.at<Vec3b>(i, j) = Vec3b({ 0,0,255 });
//	circle(thin_im, Point(j, i), 5, Scalar(0, 0, 255), 2);
//	}
//	if (gaussthin[i][j] == 3) {
//	thin_im.at<Vec3b>(i, j) = Vec3b({ 255,0,0 });
//	circle(thin_im, Point(j, i), 5, Scalar(255, 0, 0), 2);
//	}
//	}
//	}*/
//
//	imshow("gaussthin.tif", thin_im);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/GaussThin.tif", thin_im);
//
//
//}
//
//void trace1(Mat img) {
//	Mat box_im, thin_im;
//	int **gaussimg, **gaussthin, **label, **box;
//
//	gaussimg = new int*[r];
//	gaussthin = new int*[r];
//	label = new int*[r];
//	for (int i = 0; i < r; i++) {
//		gaussimg[i] = new int[c];
//		gaussthin[i] = new int[c];
//		label[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			gaussimg[i][j] = 255;
//			gaussthin[i][j] = 255;
//			label[i][j] = 0;
//		}
//	}
//	mat2arr(img, gaussimg);
//	box_im.create(r, c, CV_8UC1);
//	cvtColor(img, box_im, CV_GRAY2BGR);
//	thin_im.create(r, c, CV_8UC1);
//	arr2mat(gaussthin, thin_im);
//	cvtColor(thin_im, thin_im, CV_GRAY2BGR);
//	//cvtColor(img, thin_im, CV_GRAY2BGR);
//
//	labelNum = labelling(gaussimg, label, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	/*for (int i = 0; i < labelNum; i++) {
//	if (box[i][2] == -1) continue;
//	line(box_im, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	line(box_im, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	line(box_im, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//	line(box_im, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	}*/
//	imshow("GaussBox.tif", box_im);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/GaussBox.tif", box_im);
//
//	//Find trace lines for each label...
//	std::vector<int> Lines[6];
//	int theta = 142;
//
//	for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//
//		int **roi, **thin_roi;
//		int istart = box[i][4];
//		int jstart = box[i][2];
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//
//		roi = new int*[h];
//		thin_roi = new int*[h];
//		for (int k = 0; k < h; k++) {
//			roi[k] = new int[w];
//			thin_roi[k] = new int[w];
//			for (int l = 0; l < w; l++) {
//				thin_roi[k][l] = 255;
//				if (label[k + istart][l + jstart] == i) {
//					roi[k][l] = 0;
//				}
//				else
//					roi[k][l] = 1;
//			}
//		}
//		thinCC(roi, h, w, thin_roi);
//
//		for (int k = istart; k < h + istart; k++) {
//			for (int l = jstart; l < w + jstart; l++) {
//				if (thin_roi[k - istart][l - jstart] == 0) {
//					gaussthin[k][l] = 0;
//					thin_im.at<Vec3b>(k, l) = Vec3b({ 255,0,255 });
//				}
//			}
//		}
//
//		//Determine the terminal points...
//		std::vector<Point> terminals;
//		int sz;
//		sz = findTerminals(thin_roi, h, w, terminals);
//		std::cout << "\n Label " << i << " has " << sz << " terminal points";
//
//		std::vector<Point> corners;
//		int sz_c;
//		sz_c = findCorners(thin_roi, h, w, corners);
//		std::cout << "\n Label " << i << " has " << sz_c << " corner points";
//
//		for (int k = 0; k < sz_c; k++) {
//			gaussthin[corners[k].y + istart][corners[k].x + jstart] = 3;
//		}
//
//		for (int k = 0; k < sz; k++) {
//			gaussthin[terminals[k].y + istart][terminals[k].x + jstart] = 2;
//			Point p1, p2;
//			p1 = terminals[k];
//			double slope = tan(theta * CV_PI / 180);
//			p2.x = w - 1;
//			p2.y = p1.y + slope*(p2.x - p1.x);
//			if (p2.y < 0) {
//				p2.y = 0;
//				p2.x = p1.x - p1.y / slope;
//			}
//			//line(box_im, p1 + Point(jstart, istart), p2 + Point(jstart, istart), Scalar(255, 0, 255), 2);
//			std::vector<Point> sequence;
//			//sequence.push_back(terminals[k]);
//			//int sz = traceLine(roi, h, w, terminals[k], theta, sequence);
//			int sz = longestRay(roi, h, w, terminals[k], theta, sequence);
//			if (sz < 10) continue;
//
//			Lines[0].push_back(i);
//			Lines[1].push_back(sz);
//			if ((terminals[k].x + jstart < sequence[sz - 1].x + jstart) || (terminals[k].x + jstart == sequence[sz - 1].x + jstart && terminals[k].y + istart < sequence[sz - 1].y + istart)) {
//				Lines[2].push_back(terminals[k].y + istart);
//				Lines[3].push_back(terminals[k].x + jstart);
//				/*Lines[4].push_back(end_point.y + istart);
//				Lines[5].push_back(end_point.x + jstart);*/
//				Lines[4].push_back(sequence[sz - 1].y + istart);
//				Lines[5].push_back(sequence[sz - 1].x + jstart);
//			}
//			else if ((terminals[k].x + jstart > sequence[sz - 1].x + jstart) || (terminals[k].x + jstart == sequence[sz - 1].x + jstart && terminals[k].y + istart > sequence[sz - 1].y + istart)) {
//				Lines[2].push_back(sequence[sz - 1].y + istart);
//				Lines[3].push_back(sequence[sz - 1].x + jstart);
//				/*Lines[2].push_back(end_point.y + istart);
//				Lines[3].push_back(end_point.x + jstart);*/
//				Lines[4].push_back(terminals[k].y + istart);
//				Lines[5].push_back(terminals[k].x + jstart);
//			}
//
//			std::cout << "\nStarts at " << terminals[k] + Point(jstart, istart) << " Size = " << sz;
//			for (int t = 0; t < sz - 1; t++) {
//				//	std::cout << "\n" << sequence[t] + Point(jstart, istart);
//				line(box_im, sequence[t] + Point(jstart, istart), sequence[t + 1] + Point(jstart, istart), Scalar(0, 0, 255), 3);
//				box_im.at<Vec3b>(sequence[t].y + istart, sequence[t].x + jstart) = Vec3b({ 255,0,0 });
//			}
//			//std::cout << "\n" << sequence[sz - 1] + Point(jstart, istart);
//			box_im.at<Vec3b>(sequence[sz - 1].y + istart, sequence[sz - 1].x + jstart) = Vec3b({ 255,0,0 });
//			//line(box_im, sequence[0] + Point(jstart, istart), end_point + Point(jstart, istart), Scalar(255, 0, 0), 3);
//		}
//
//		//delete...
//		for (int k = 0; k < h; k++) {
//			delete[] roi[k];
//			delete[] thin_roi[k];
//		}
//		delete[] roi;
//		delete[] thin_roi;
//
//	}
//	int S = int(Lines[0].size());
//	//std::vector<std::vector<Point>> lineSet = extendLines(Lines, S, theta);
//	/*std::cout << "\nSorted List of the labels and there sequence lines:";
//	for (int i = 0; i < S; i++) {
//	std::cout << "\nLabel " << Lines[0][i] << " has " << Lines[1][i] << " steps from (" << Lines[2][i] << " , " << Lines[3][i] << ") to (" << Lines[4][i] << " , " << Lines[5][i] << ")";
//	}*/
//
//	//std::cout << "\n\nblack:white ratios->";
//	//for (int i = 0; i<int(lineSet.size()); i++) {
//	//	double white_length = 0, black_length = 0;
//	//	for (int j = 1; j<int(lineSet[i].size()) - 1; j++) {
//	//		if (j % 2 == 1) { //intra-label segment => black...
//	//			black_length += norm(lineSet[i][j + 1] - lineSet[i][j]);
//	//		}
//	//		else { //inter-label segment => white...
//	//			double slope = atan2(lineSet[i][j + 1].y - lineSet[i][j].y, lineSet[i][j + 1].x - lineSet[i][j].x) * 180 / CV_PI;
//	//			int start = 0;
//	//			Point prev = lineSet[i][j];
//	//			for (int x = lineSet[i][j].x; x <= lineSet[i][j + 1].x; x++) {
//	//				int y = int(lineSet[i][j].y + tan(slope*CV_PI / 180)*(x - lineSet[i][j].x)); //y = y_0 + m*(x - x_0)
//	//				if (gaussimg[y][x] == 0 && start == 0) { //new black region...
//	//					start = 1;
//	//					white_length += norm(Point(x, y) - prev);
//	//					prev = Point(x, y);
//	//				}
//	//				if (gaussimg[y][x] != 0 && start == 1) { //new white region...
//	//					start = 0;
//	//					black_length += norm(Point(x, y) - prev);
//	//					prev = Point(x, y);
//	//				}
//	//			}
//	//			//white_length += norm(lineSet[i][j + 1] - lineSet[i][j]);
//	//		}
//	//		//line(box_im, lineSet[i][j], lineSet[i][j + 1], Scalar(255, 0, 0), 3);
//	//	}
//	//	std::cout << "\nLine " << i + 1 << " : white space = " << white_length << " , text space = " << black_length;
//	//	Vec3b color;
//	//	if (white_length < 0.6*black_length)
//	//		color = Vec3b({ 255,120,0 });
//	//	else
//	//		color = Vec3b({ 0,120,255 });
//	//	for (int j = 1; j<int(lineSet[i].size()) - 1; j++) {
//	//		line(box_im, lineSet[i][j], lineSet[i][j + 1], color, 3);
//	//	}
//
//	//}
//
//	imshow("gausstrace.tif", box_im);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/gausstrace.tif", box_im);
//
//	//	arr2mat(gaussthin, thin_im);
//	//	cvtColor(thin_im, thin_im, CV_GRAY2BGR);
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (gaussthin[i][j] == 2) {
//				thin_im.at<Vec3b>(i, j) = Vec3b({ 0,0,255 });
//				circle(thin_im, Point(j, i), 5, Scalar(0, 0, 255), 2);
//			}
//			if (gaussthin[i][j] == 3) {
//				thin_im.at<Vec3b>(i, j) = Vec3b({ 255,0,0 });
//				circle(thin_im, Point(j, i), 5, Scalar(255, 0, 0), 2);
//			}
//		}
//	}
//
//	imshow("gaussthin.tif", thin_im);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/GaussThin.tif", thin_im);
//
//
//}
//
//
//void invCC(Mat img, int flag = 0) {
//	int **orgimg, **label, **box, **binimg;
//	Mat colorimg, bin_img, thin_im, near4img;
//
//	orgimg = new int*[r];
//	binimg = new int*[r];
//	label = new int*[r];
//	for (int i = 0; i < r; i++) {
//		orgimg[i] = new int[c];
//		binimg[i] = new int[c];
//		label[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			orgimg[i][j] = 0;
//			binimg[i][j] = 0;
//			label[i][j] = 0;
//		}
//	}
//	mat2arr(img, orgimg);
//
//	//----------------------------------------------------------------------------------------------------------
//	//Binarize if needed...
//	if (flag) {
//		NICK(orgimg, r, c, binimg);
//		labelNum = labelling(binimg, label, r, c);
//
//		//Bound the components by rectangular & diagonal boxes
//		box = new int*[labelNum];
//		for (int i = 0; i < labelNum; i++)
//			box[i] = new int[6];
//		boundingRectangles(box, label, r, c, labelNum);
//
//		for (int i = 1; i < labelNum; i++) {
//			if (box[i][2] == -1) continue;
//			int h = box[i][5] - box[i][4] + 1;
//			int w = box[i][3] - box[i][2] + 1;
//			int count = 0;
//			for (int k = box[i][4]; k <= box[i][5]; k++) {
//				for (int l = box[i][2]; l <= box[i][3]; l++) {
//					if (binimg[k][l] == 0) count++;
//				}
//			}
//			if (h*w < 30 || count < 100) {
//				for (int k = box[i][4]; k <= box[i][5]; k++) {
//					for (int l = box[i][2]; l <= box[i][3]; l++) {
//						if (label[k][l] == i) {
//							label[k][l] = 0;
//							binimg[k][l] = 255;
//						}
//					}
//				}
//			}
//		}
//
//		bin_img.create(r, c, CV_8UC1);
//		arr2mat(binimg, bin_img);
//		GaussianBlur(bin_img, bin_img, Size(3, 3), 0, 0);
//		mat2arr(bin_img, orgimg);
//		NICK(orgimg, r, c, binimg);
//		orgimg = binimg;
//		arr2mat(binimg, bin_img);
//		imshow("binaryImg.tif", bin_img);
//		imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/binImg.tif", bin_img);
//		//imwrite("Data/ICDAR/Data218/binImg.tif", bin_img);
//
//		for (int i = 0; i < labelNum; i++)
//			delete[] box[i];
//		delete[] box;
//		for (int i = 0; i < r; i++) {
//			for (int j = 0; j < c; j++)
//				label[i][j] = 0;
//		}
//		labelNum = 0;
//	}
//
//	//----------------------------------------------------------------------------------------------------------
//	//Fill holes...
//	erode(orgimg);
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			orgimg[i][j] = 255 - orgimg[i][j];
//		}
//	}
//	//dilate(orgimg);
//
//	labelNum = labelling(orgimg, label, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	//fillHoles(orgimg, label, box);
//
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			orgimg[i][j] = 255 - orgimg[i][j];
//		}
//	}
//
//	//---------------------------------------------------------------------------------------------------------
//	//re-label components...
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			label[i][j] = 0;
//		}
//	}
//	for (int i = 0; i < labelNum; i++)
//		delete[] box[i];
//	delete[] box;
//	labelNum = 0;
//	std::cout << "\nRe-labelling...";
//	labelNum = labelling(orgimg, label, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	//dilate(orgimg);
//	arr2mat(orgimg, img);
//	colorimg.create(r, c, CV_8UC1);
//	colorimg = img;
//
//	imshow("invcc.tif", colorimg);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/invcc.tif", colorimg);
//
//	cvtColor(img, colorimg, CV_GRAY2BGR);
//	near4img.create(r, c, CV_8UC1);
//	cvtColor(img, near4img, CV_GRAY2BGR);
//	/*for (int i = 0; i < labelNum; i++) {
//	if (box[i][2] == -1) continue;
//	line(colorimg, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	line(colorimg, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	line(colorimg, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(0, 0, 255), 1.5, 8);
//	line(colorimg, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(0, 0, 255), 1.5, 8);
//	}*/
//
//	//imwrite("Data/ICDAR/Data218/invcc.tif", img);
//
//	//---------------------------------------------------------------------------------------------------------
//	//Exclude small components...
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//		int count = 0;
//		for (int k = box[i][4]; k <= box[i][5]; k++) {
//			for (int l = box[i][2]; l <= box[i][3]; l++) {
//				if (orgimg[k][l] == 0) count++;
//			}
//		}
//		if (h*w < 30 || count < 200) {
//			for (int k = 0; k < 6; k++)
//				box[i][k] = -1;
//		}
//	}
//
//	//Find medians and neighbour network...
//	std::vector<int> nodes;
//	std::vector<Point> centroids;
//	for (int i = 1; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//		int istart = box[i][4];
//		int jstart = box[i][2];
//		int h = box[i][5] - box[i][4] + 1;
//		int w = box[i][3] - box[i][2] + 1;
//
//		std::vector<int> X, Y;
//		for (int k = 0; k < h; k++) {
//			for (int l = 0; l < w; l++) {
//				if (label[istart + k][jstart + l] == i)
//					Y.push_back(istart + k);
//			}
//		}
//		for (int l = 0; l < w; l++) {
//			for (int k = 0; k < h; k++) {
//				if (label[istart + k][jstart + l] == i)
//					X.push_back(jstart + l);
//			}
//		}
//		int sz = int(X.size());
//		Point med = Point(int(median(X, sz)), int(median(Y, sz)));
//		nodes.push_back(i);
//		centroids.push_back(med);
//
//		circle(colorimg, med, 5, Scalar(0, 0, 255), -2);
//	}
//
//	int LN = int(nodes.size());
//	std::cout << "\nNo. of labels = " << LN;
//
//	for (int i = 0; i < LN; i++) {
//		std::cout << "\n Label " << nodes[i] << " -> " << centroids[i];
//	}
//
//	double **distance = new double*[LN];
//	for (int i = 0; i < LN; i++) {
//		distance[i] = new double[LN];
//		for (int j = 0; j < LN; j++) {
//			distance[i][j] = 0;
//		}
//	}
//	int countLines = 0;
//	for (int i = 0; i < LN; i++) {
//		int I = nodes[i];
//		std::cout << "\nJoining with label " << I << "--->";
//		for (int j = i + 1; j < LN; j++) {
//			int J = nodes[j];
//			Point p1 = centroids[i];
//			Point p2 = centroids[j];
//			int flag = 1;
//
//			if (abs(p1.y - p2.y) > 400 || abs(p1.x - p2.x) > 400) continue;
//
//			if (p1.x == p2.x) {
//				int imin = p1.y < p2.y ? p1.y : p2.y;
//				int imax = p1.y > p2.y ? p1.y : p2.y;
//				for (int k = imin; k <= imax; k++) {
//					if (orgimg[k][p1.x] == 0 && (label[k][p1.x] != nodes[i] && label[k][p1.x] != nodes[j])) {
//						flag = 0;
//						break;
//					}
//				}
//			}
//			else if (p1.y == p2.y) {
//				int jmin = p1.x < p2.x ? p1.x : p2.x;
//				int jmax = p1.x > p2.x ? p1.x : p2.x;
//				for (int l = jmin; l <= jmax; l++) {
//					if (orgimg[p1.y][l] == 0 && (label[p1.y][l] != nodes[i] && label[p1.y][l] != nodes[j])) {
//						flag = 0;
//						break;
//					}
//				}
//			}
//			else {
//				double slope = double(p1.y - p2.y) / double(p1.x - p2.x);
//				int jmin = p1.x <= p2.x ? p1.x : p2.x;
//				int jmax = p1.x >= p2.x ? p1.x : p2.x;
//				for (int l = jmin; l <= jmax; l++) {
//					int k = int(p1.y + slope*(l - p1.x));
//					if (orgimg[k][l] == 0 && (label[k][l] != nodes[i] && label[k][l] != nodes[j])) {
//						flag = 0;
//						break;
//					}
//				}
//				if (flag) {
//					int imin = p1.y <= p2.y ? p1.y : p2.y;
//					int imax = p1.y >= p2.y ? p1.y : p2.y;
//					for (int k = imin; k <= imax; k++) {
//						int l = int(p1.x + (k - p1.y) / slope);
//						if (orgimg[k][l] == 0 && (label[k][l] != nodes[i] && label[k][l] != nodes[j])) {
//							flag = 0;
//							break;
//						}
//					}
//				}
//			}
//			if (flag) {
//				line(colorimg, centroids[i], centroids[j], Scalar(120, 0, 255), 2);
//				countLines++;
//				int **roi;
//				int istart = box[I][4] <= box[J][4] ? box[I][4] : box[J][4];
//				int istop = box[I][5] >= box[J][5] ? box[I][5] : box[J][5];
//				int jstart = box[I][2] <= box[J][2] ? box[I][2] : box[J][2];
//				int jstop = box[I][3] >= box[J][3] ? box[I][3] : box[J][3];
//				int h = istop - istart + 1;
//				int w = jstop - jstart + 1;
//				roi = new int*[h];
//				for (int k = 0; k < h; k++) {
//					roi[k] = new int[w];
//					for (int l = 0; l < w; l++) {
//						if (label[k + istart][l + jstart] == I)
//							roi[k][l] = I;
//						else if (label[k + istart][l + jstart] == J)
//							roi[k][l] = J;
//						else
//							roi[k][l] = 0;
//					}
//				}
//				double dist = HDistance(roi, h, w, box[I][4] - istart, box[J][4] - istart, box[I][5] - istart, box[J][5] - istart, box[I][2] - jstart, box[J][2] - jstart, box[I][3] - jstart, box[J][3] - jstart, I, J);
//				std::cout << "\nDistance between label " << I << " and label " << J << " = " << dist;
//				distance[i][j] = dist;
//				distance[j][i] = dist;
//			}
//		}
//	}
//	std::cout << "\nNo. of lines = " << countLines;
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/AllNetwork.tif", colorimg);
//
//	//std::cout << "\nDistance matrix-->";
//
//	int **adjMatrix = new int*[LN];
//	for (int i = 0; i < LN; i++) {
//		adjMatrix[i] = new int[LN];
//		for (int j = 0; j < LN; j++) {
//			adjMatrix[i][j] = 0;
//		}
//	}
//	for (int i = 0; i < LN; i++) {
//		int j1 = -1, j2 = -1, j3 = -1, j4 = -1;
//		double min1 = 99999, min2 = 99999, min3 = 99999, min4 = 99999;
//		for (int j = 0; j < LN; j++) {
//			if (distance[i][j] == 0) continue;
//			if (distance[i][j] < min1) {
//				min1 = distance[i][j];
//				j1 = j;
//			}
//		}
//		for (int j = 0; j < LN; j++) {
//			if (distance[i][j] == 0) continue;
//			if (distance[i][j] < min2 && distance[i][j] > min1) {
//				min2 = distance[i][j];
//				j2 = j;
//			}
//		}
//		for (int j = 0; j < LN; j++) {
//			if (distance[i][j] == 0) continue;
//			if (distance[i][j] < min3 && distance[i][j] > min2) {
//				min3 = distance[i][j];
//				j3 = j;
//			}
//		}
//		for (int j = 0; j < LN; j++) {
//			if (distance[i][j] == 0) continue;
//			if (distance[i][j] < min4 && distance[i][j] > min3) {
//				min4 = distance[i][j];
//				j4 = j;
//			}
//		}
//		if (j1 != -1)
//			adjMatrix[i][j1] = 1;
//		if (j2 != -1)
//			adjMatrix[i][j2] = 1;
//		if (j3 != -1)
//			adjMatrix[i][j3] = 1;
//		if (j4 != -1)
//			adjMatrix[i][j4] = 1;
//		for (int j = 0; j < LN; j++) {
//			if ((j != j1 && j != j2 && j != j3 && j != j4))
//				distance[i][j] = 0;
//
//			//std::cout << distance[i][j] << "   ";
//		}
//		//std::cout << "\n";
//	}
//
//	int  *angles = new int[180];
//	for (int i = 0; i < 180; i++)
//		angles[i] = 0;
//	for (int i = 0; i < LN; i++) {
//		for (int j = i + 1; j < LN; j++) {
//			if (adjMatrix[i][j] == 0 && adjMatrix[j][i] == 0) continue;
//			double slope = atan2(centroids[i].y - centroids[j].y, centroids[i].x - centroids[j].x) * 180 / CV_PI;
//			int theta = floor(slope);
//			if (theta < 0) theta = theta + 180;
//			if (theta == 180) theta = 0;
//			angles[theta]++;
//			//std::cout << "\nlabel " << nodes[i] << " to label " << nodes[j] << " -> " << theta << " degrees";
//		}
//	}
//	std::cout << "\nNo. of lines in each bin->";
//	for (int i = 0; i < 180; i++)
//		std::cout << "\n" << i << " degrees : " << angles[i];
//	std::cout << "\nAdjacency matrix-->";
//	for (int i = 0; i < LN; i++) {
//		for (int j = 0; j < LN; j++) {
//			//std::cout << adjMatrix[i][j] << "   ";
//			if (adjMatrix[i][j] == 1) {
//				int I = nodes[i], J = nodes[j];
//				if (i < j)
//					arrowedLine(near4img, centroids[i], centroids[j], Scalar(0, 120, 255), 2, 8, 0, 0.1);
//				if (i > j)
//					arrowedLine(near4img, centroids[i], centroids[j], Scalar(255, 120, 0), 2, 8, 0, 0.1);
//			}
//
//		}
//		//std::cout << "\n";
//	}
//
//	//---------------------------------------------------------------------------------------------------------
//	//Find inclination of components...
//	//int **thin_img = new int*[r];
//	//for (int i = 0; i < r; i++) {
//	//	thin_img[i] = new int[c];
//	//	for (int j = 0; j < c; j++)
//	//		thin_img[i][j] = 255;
//	//}
//	//thin_im.create(r, c, CV_8UC1);
//	//arr2mat(thin_img, thin_im);
//	//cvtColor(thin_im, thin_im, CV_GRAY2BGR);
//	//for (int i = 1; i < labelNum; i++) {
//	//	if (box[i][2] == -1) continue;
//	//	int **roi, **thin_roi;
//	//	int istart = box[i][4];
//	//	int jstart = box[i][2];
//	//	int h = box[i][5] - box[i][4] + 1;
//	//	int w = box[i][3] - box[i][2] + 1;
//
//	//	std::vector<int> X, Y;
//	//	roi = new int*[h];
//	//	thin_roi = new int*[h];
//	//	for (int k = 0; k < h; k++) {
//	//		roi[k] = new int[w];
//	//		thin_roi[k] = new int[w];
//	//		for (int l = 0; l < w; l++) {
//	//			thin_roi[k][l] = 255;
//	//			if (label[k + istart][l + jstart] == i) {
//	//				roi[k][l] = 0;
//	//			}
//	//			else
//	//				roi[k][l] = 1;
//	//		}
//	//	}
//
//	//	findBoundary(roi, h, w);
//	//	/*thinCC(roi, h, w, thin_roi);*/
//
//	//	/*for (int k = istart; k < h + istart; k++) {
//	//		for (int l = jstart; l < w + jstart; l++) {
//	//			if (roi[k - istart][l - jstart] == 2) {
//	//				thin_im.at<Vec3b>(k, l) = Vec3b({ 0,0,0 });
//	//			}
//	//			else
//	//				thin_roi[k - istart][l - jstart] = 1;
//	//		}
//	//	}*/
//	//	for (int k = istart; k < h + istart; k++) {
//	//		for (int l = jstart; l < w + jstart; l++) {
//	//			if (roi[k - istart][l - jstart] == 2) 
//	//				thin_im.at<Vec3b>(k, l) = Vec3b({ 0,0,0 });
//	//			else if (roi[k - istart][l - jstart] == 1)
//	//				thin_im.at<Vec3b>(k, l) = Vec3b({ 255,255,255 });
//	//		}
//	//	}
//	//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/boundary.tif", thin_im);
//
//	//	//collect Y co-ordinates...
//	//	for (int k = 0; k < h; k++) {
//	//		for (int l = 0; l < w; l++) {
//	//			if (roi[k][l] == 0)
//	//				Y.push_back(k);
//	//		}
//	//	}
//	//	//collect X co-ordinates...
//	//	for (int l = 0; l < w; l++) {
//	//		for (int k = 0; k < h; k++) {
//	//			if (roi[k][l] == 0)
//	//				X.push_back(l);
//	//		}
//	//	}
//	//	std::cout << "\n***Label = " << i;
//	//	int inclineStdDev, inclineExtent, inclineCrossing;
//	//	//projectAll(roi, h, w, inclineStdDev, inclineExtent, inclineCrossing);
//	//	projectAllCross(roi, h, w, inclineCrossing);
//	//	//projectAllStdDev(roi, h, w, inclineStdDev);
//	//	Point p, q;
//	//	int sz = int(X.size());
//	//	std::cout << "\nsize = " << sz;
//	//	Point med = Point(int(mean(X, sz)), int(mean(Y, sz)));
//	//	std::cout << "\nMean = " << med;
//	//	//thin_im.at<Vec3b>(med.y + istart, med.x + jstart) = Vec3b({ 255,0,255 });
//
//	//	for (int theta = 0; theta <= 175; theta += 5) {
//	//		Vec3b color;
//	//		/*if (theta == inclineStdDev && theta != inclineExtent && theta != inclineCrossing)
//	//			color = { 255,0,255 };
//	//		else if (theta != inclineStdDev && theta == inclineExtent && theta != inclineCrossing)
//	//			color = { 255,0,0 };
//	//		else if (theta != inclineStdDev && theta != inclineExtent && theta == inclineCrossing)
//	//			color = { 0,120,255 };
//	//		else if ((theta == inclineStdDev && theta == inclineExtent && theta != inclineCrossing) || (theta == inclineStdDev && theta != inclineExtent && theta == inclineCrossing) || (theta != inclineStdDev && theta == inclineExtent && theta == inclineCrossing))
//	//			color = { 120,0,0 };
//	//		else if (theta == inclineStdDev && theta == inclineExtent && theta == inclineCrossing)
//	//			color = { 0,0,255 };	*/	
//	//		if (theta == inclineCrossing)
//	//			color = { 255,0,255 };
//	//		/*else if (theta == inclineStdDev)
//	//			color = { 0,120,255 };*/
//	//		else
//	//			continue;
//
//	//		if (theta == 0) {
//	//			p = Point(0, med.y);
//	//			q = Point(w - 1, med.y);
//	//		}
//	//		else if (theta == 90) {
//	//			p = Point(med.x, 0);
//	//			q = Point(med.x, h - 1);
//	//		}
//	//		else if (theta > 90) {
//	//			double slope = tan(theta*CV_PI / 180);
//	//			double intercept = double(med.y) - double(slope*med.x);
//	//			if (intercept < h) {
//	//				p.x = 0;
//	//				p.y = int(intercept);
//	//			}
//	//			else {
//	//				p.y = h - 1;
//	//				p.x = int((p.y - intercept) / slope);
//	//			}
//	//			if (-intercept / slope < w) {
//	//				q.x = int(-intercept / slope);
//	//				q.y = 0;
//	//			}
//	//			else {
//	//				q.x = w - 1;
//	//				q.y = int(q.x * slope + intercept);
//	//			}
//	//		}
//	//		else if (theta < 90) {
//	//			double slope = tan(theta*CV_PI / 180);
//	//			double intercept = double(med.y) - double(slope*med.x);
//	//			if (intercept >= 0) {
//	//				p.x = 0;
//	//				p.y = int(intercept);
//	//			}
//	//			else {
//	//				p.x = int(-intercept / slope);
//	//				p.y = 0;
//	//			}
//	//			if ((w - 1)*slope + intercept < h) {
//	//				q.x = w - 1;
//	//				q.y = (w - 1)*slope + intercept;
//	//			}
//	//			else {
//	//				q.y = h - 1;
//	//				q.x = int((q.y - intercept) / slope);
//	//			}
//	//		}
//	//		std::cout << "\nDraw line through " << p + Point(jstart, istart) << " , " << q + Point(jstart, istart);
//
//	//		line(colorimg, p + Point(jstart, istart), q + Point(jstart, istart), color, 3, 8);
//	//	}
//	//	
//	//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/inclinationsCross.tif", colorimg);
//	//	imwrite("Data/ICDAR/Data218/inclinations1.tif", colorimg);
//
//	//	//delete...
//	//	for (int k = 0; k < h; k++) {
//	//		delete[] roi[k];
//	//		delete[] thin_roi[k];
//	//	}
//	//	delete[] roi;
//	//	delete[] thin_roi;
//	//}
//
//
//
//	//imshow("invcc.tif", near4img);
//	//imwrite("Data/ICDAR/Data218/inclinations1.tif", colorimg);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/Near4Network.tif", near4img);
//
//	/*imshow("boundary.tif", thin_im);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan005/boundary.tif", thin_im);*/
//
//	//MultiOIC(orgimg, img);
//
//	//StripVer(orgimg, img);
//	//StripHor(orgimg, img);
//	//StripDiag(orgimg, img);
//	//StripOffDiag(orgimg, img);
//
//	//Delete...
//	for (int i = 0; i < r; i++) {
//		delete[] orgimg[i];
//	}
//	delete[] orgimg;
//}
//
///*********************************************************************************************************************************************/
//
//
//
//int main(int argc, char** argv) {
//
//	//-----------------------------Read image-----------------------
//	//Mat scanimg = imread("Data/ICDAR/Data218/scan218.tif", IMREAD_GRAYSCALE);
//	//Mat scanimg = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/scan005.tif", IMREAD_GRAYSCALE);
//	Mat scanimg = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/rotate-50org.tif", IMREAD_GRAYSCALE);
//	//Mat orgimg = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/invccorg.tif", IMREAD_GRAYSCALE);
//
//	//Mat img0 = imread("Data/ICDAR/Data218/Gauss0.tif", IMREAD_GRAYSCALE);
//	//Mat img0 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/Gauss0.tif", IMREAD_GRAYSCALE);
//	//Mat img-30 = imread("Data/ICDAR/Data218/Gauss-30.tif", IMREAD_GRAYSCALE);
//	//Mat img-30 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/Gauss-30.tif", IMREAD_GRAYSCALE);
//	//Mat img135 = imread("Data/ICDAR/Data218/Gauss135.tif", IMREAD_GRAYSCALE);
//	//Mat img135 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan005/Gauss137.tif", IMREAD_GRAYSCALE);
//	//Mat scancol;
//
//	imshow("Image", scanimg);
//	std::cout << "Image dimensions:" << scanimg.rows << "x" << scanimg.cols << endl;
//
//	r = scanimg.rows;
//	c = scanimg.cols;
//
//	/*scancol.create(r, c, CV_8UC1);
//	cvtColor(scanimg, scancol, CV_GRAY2BGR);*/
//	//invCC(scanimg, 1);
//	invCCProj(scanimg, 1);
//
//	/*imshow("Image", orgimg);
//	std::cout << "Image dimensions:" << orgimg.rows << "x" << orgimg.cols << endl;
//
//	r = orgimg.rows;
//	c = orgimg.cols;
//	regionMask2(orgimg);*/
//
//	/*scancol.create(r, c, CV_8UC1);
//	cvtColor(img135, scancol, CV_GRAY2BGR);
//	trace(img135);*/
//
//	std::cout << "\nDone!\n";
//	waitKey(0);
//	system("pause");
//	return 0;
//
//}
