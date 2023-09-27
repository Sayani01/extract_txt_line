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
//	std::cout << "\n" << arr[i];*/
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
///***********************************************************************************************************
//											PROJECTION PROFILES
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
//	if (abs(peaks[i - 1][1] - peaks[i][0]) < 30) {
//	peaks[i - 1][1] = peaks[i][1];
//	peaks[i - 1][2] = peaks[i - 1][2] > peaks[i][2] ? peaks[i - 1][2] : peaks[i][2];
//	peaks.erase(peaks.begin() + i);
//	i--;
//	}
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
//	hist1[i] = new int[2];
//	for (int j = 0; j < 2; j++)
//	hist1[i][j] = 0;
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
//	int start = -1, stop = -1;
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
//	std::vector<std::vector<int>> pv;
//
//	if (start < 0 || stop < 0)
//		return(pv);
//
//	int max = hist[0];
//	int *seq = new int[sz]; //1 => possible peak , -1 => possible valley
//	for (int i = 1; i < sz; i++) {
//		max = hist[i] > max ? hist[i] : max;
//	}
//	/*double avg = 0;
//	for (int i = start; i <= stop; i++) {
//	avg += hist[i];
//	}
//	avg /= (stop - start + 1);*/
//
//	double thresh = 0.3*max;
//
//	for (int i = 0; i < sz; i++) {
//		seq[i] = hist[i] >= thresh ? 1 : -1;
//	}
//
//
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
//			continue;*/
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
//			continue;*/
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
//	for (int j = c - 1; j >= 0; j--) {
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
//	for (int k = firstcol; k < c; k += int(N / 2)) {
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
//		Y[i] = 0;
//		whiteLines(hist, r, Y);
//		*/
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
//			if (peaks[t][2] > N && peaks[t][2] > 2.5*abs(peaks[t][1] - peaks[t][0])) {
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
//				for (int x = k; x <= k + w; x++) {
//				reg[y][x] = 0;
//				}
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
//				if ((vleft[1]) / double(p[1]) < 0.05 && (vright[1]) / double(p[1]) < 0.05 && p[1] > 7.0*double(2 * thickness)) {
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
//			mark[pv[t][1]] = 1;
//			line(regions, Point(k, pv[t][1]), Point(k + w, pv[t][1]), Scalar(255, 0, 0), 8, 8);
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
//			line(colorBinary, Point(k, i), Point(k + N, i), Scalar(255, 0, 0), 2.5, 8);*/
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
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/VerticalStripsProj15.tif", colorBinary);
//
//	reg_mask.create(r, c, CV_8UC1);
//	erode(reg);
//	arr2mat(reg, reg_mask);
//	cvtColor(reg_mask, reg_mask, CV_GRAY2BGR);
//
//	imshow("regions.tif", regions);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/VerticalStripsRegions15.tif", regions);
//
//	imshow("Mask_regions.tif", reg_mask);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/VerticalStripsRegionsMask15.tif", reg_mask);
//	//imwrite("Data/ICDAR/Data218/VerticalStrips.tif", colorBinary);
//}
//
//int findLongestWord(Mat bin_img, int **mask, int **label, int labelNum, int **box, int r, int c) {
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
//		if (h > w)
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
//					mask[box[i][4] + k][box[i][2] + l] = 0;
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
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/binImgProj.tif", bin_img);*/
//
//	return(int(1.25 * (box[max_lbl][3] - box[max_lbl][2] + 1)));
//	//return(max_hist);
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
////void invCCProj(Mat img, int flag = 0) {
////	int **orgimg, **label, **box, **binimg;
////	Mat colorimg, bin_img, thin_im, near4img;
////
////	orgimg = new int*[r];
////	binimg = new int*[r];
////	label = new int*[r];
////	for (int i = 0; i < r; i++) {
////		orgimg[i] = new int[c];
////		binimg[i] = new int[c];
////		label[i] = new int[c];
////		for (int j = 0; j < c; j++) {
////			orgimg[i][j] = 0;
////			binimg[i][j] = 0;
////			label[i][j] = 0;
////		}
////	}
////	mat2arr(img, orgimg);
////
////	//----------------------------------------------------------------------------------------------------------
////	//Binarize if needed...
////	if (flag) {
////		NICK(orgimg, r, c, binimg);
////		labelNum = labelling(binimg, label, r, c);
////
////		//Bound the components by rectangular & diagonal boxes
////		box = new int*[labelNum];
////		for (int i = 0; i < labelNum; i++)
////			box[i] = new int[6];
////		boundingRectangles(box, label, r, c, labelNum);
////
////		for (int i = 1; i < labelNum; i++) {
////			if (box[i][2] == -1) continue;
////			int h = box[i][5] - box[i][4] + 1;
////			int w = box[i][3] - box[i][2] + 1;
////			int count = 0;
////			for (int k = box[i][4]; k <= box[i][5]; k++) {
////				for (int l = box[i][2]; l <= box[i][3]; l++) {
////					if (binimg[k][l] == 0) count++;
////				}
////			}
////			if (h*w < 30 || count < 100) {
////				for (int k = box[i][4]; k <= box[i][5]; k++) {
////					for (int l = box[i][2]; l <= box[i][3]; l++) {
////						if (label[k][l] == i) {
////							label[k][l] = 0;
////							binimg[k][l] = 255;
////						}
////					}
////				}
////				box[i][2] = -1;
////				box[i][3] = -1;
////				box[i][4] = -1;
////				box[i][5] = -1;
////			}
////		}
////
////		bin_img.create(r, c, CV_8UC1);
////		arr2mat(binimg, bin_img);
////		GaussianBlur(bin_img, bin_img, Size(3, 3), 0, 0);
////		mat2arr(bin_img, orgimg);
////		NICK(orgimg, r, c, binimg);
////		orgimg = binimg;
////		arr2mat(binimg, bin_img);
////		imshow("binaryImg.tif", bin_img);
////		imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/binImg15.tif", bin_img);
////
////		//----------------------------------------------------------------------------------------------------------
////		for (int i = 0; i < labelNum; i++)
////			delete[] box[i];
////		delete[] box;
////		for (int i = 0; i < r; i++) {
////			for (int j = 0; j < c; j++)
////				label[i][j] = 0;
////		}
////		labelNum = 0;
////	}
////
////	//Fill holes...
////	erode(orgimg);
////	for (int i = 0; i < r; i++) {
////		for (int j = 0; j < c; j++) {
////			orgimg[i][j] = 255 - orgimg[i][j];
////		}
////	}
////	//dilate(orgimg);
////
////	labelNum = labelling(orgimg, label, r, c);
////
////	//Bound the components by rectangular & diagonal boxes
////	box = new int*[labelNum];
////	for (int i = 0; i < labelNum; i++)
////		box[i] = new int[6];
////	boundingRectangles(box, label, r, c, labelNum);
////
////	//fillHoles(orgimg, label, box);
////
////	for (int i = 0; i < r; i++) {
////		for (int j = 0; j < c; j++) {
////			orgimg[i][j] = 255 - orgimg[i][j];
////		}
////	}
////
////	//---------------------------------------------------------------------------------------------------------
////	//re-label components...
////	for (int i = 0; i < r; i++) {
////		for (int j = 0; j < c; j++) {
////			label[i][j] = 0;
////		}
////	}
////	for (int i = 0; i < labelNum; i++)
////		delete[] box[i];
////	delete[] box;
////	labelNum = 0;
////	std::cout << "\nRe-labelling...";
////	labelNum = labelling(orgimg, label, r, c);
////
////	//Bound the components by rectangular & diagonal boxes
////	box = new int*[labelNum];
////	for (int i = 0; i < labelNum; i++)
////		box[i] = new int[6];
////	boundingRectangles(box, label, r, c, labelNum);
////
////	//dilate(orgimg);
////	arr2mat(orgimg, img);
////	colorimg.create(r, c, CV_8UC1);
////	colorimg = img;
////
////	imshow("invcc.tif", colorimg);
////	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/invcc15.tif", colorimg);
////
////	cvtColor(img, colorimg, CV_GRAY2BGR);
////	near4img.create(r, c, CV_8UC1);
////	cvtColor(img, near4img, CV_GRAY2BGR);
////
////	//imwrite("Data/ICDAR/Data218/invcc.tif", img);
////
////	//---------------------------------------------------------------------------------------------------------
////	//Exclude small components...
////	for (int i = 1; i < labelNum; i++) {
////		if (box[i][2] == -1) continue;
////		int h = box[i][5] - box[i][4] + 1;
////		int w = box[i][3] - box[i][2] + 1;
////		int count = 0;
////		for (int k = box[i][4]; k <= box[i][5]; k++) {
////			for (int l = box[i][2]; l <= box[i][3]; l++) {
////				if (orgimg[k][l] == 0) count++;
////			}
////		}
////		if (h*w < 30 || count < 200) {
////			for (int k = 0; k < 6; k++)
////				box[i][k] = -1;
////		}
////	}
////
////	int N = findLongestWord(colorimg, label, labelNum, box, r, c);
////
////	for (int i = 0; i < labelNum; i++) {
////		if (box[i][2] == -1) continue;
////		line(colorimg, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(255, 0, 0), 2, 8);
////		line(colorimg, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(255, 0, 0), 2, 8);
////		line(colorimg, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(255, 0, 0), 2, 8);
////		line(colorimg, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(255, 0, 0), 2, 8);
////	}
////	imshow("binaryImgProj.tif", colorimg);
////	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/binImgProj15.tif", colorimg);
////	
////	StripVer(orgimg, img, N);
////
////	//Delete...
////	for (int i = 0; i < r; i++) {
////		delete[] orgimg[i];
////	}
////	delete[] orgimg;
////}
//
//int validAngles(Mat img, Mat rot, Mat mask_img) {
//	int **orgimg, **label, **box, **binimg, labelNum = 0, **mask;
//	Mat colorimg, bin_img;
//
//	orgimg = new int*[r];
//	binimg = new int*[r];
//	label = new int*[r];
//	mask = new int*[r];
//	for (int i = 0; i < r; i++) {
//		orgimg[i] = new int[c];
//		binimg[i] = new int[c];
//		label[i] = new int[c];
//		mask[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			orgimg[i][j] = 0;
//			binimg[i][j] = 0;
//			label[i][j] = 0;
//			mask[i][j] = 255;
//		}
//	}
//	mat2arr(img, orgimg);
//
//	labelNum = labelling(orgimg, label, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	box = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box[i] = new int[6];
//	boundingRectangles(box, label, r, c, labelNum);
//
//	//---------------------------------------------------------------------------------------------------------
//	////Color original image...
//	//colorimg.create(r, c, CV_8UC1);
//	//colorimg = img;
//	//cvtColor(img, colorimg, CV_GRAY2BGR);
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
//	int N = findLongestWord(rot, mask, label, labelNum, box, r, c);	
//	arr2mat(mask, mask_img);
//
//	/*for (int i = 0; i < labelNum; i++) {
//		if (box[i][2] == -1) continue;
//		line(rot, Point(box[i][2], box[i][4]), Point(box[i][2], box[i][5]), Scalar(255, 0, 0), 2, 8);
//		line(rot, Point(box[i][2], box[i][5]), Point(box[i][3], box[i][5]), Scalar(255, 0, 0), 2, 8);
//		line(rot, Point(box[i][2], box[i][4]), Point(box[i][3], box[i][4]), Scalar(255, 0, 0), 2, 8);
//		line(rot, Point(box[i][3], box[i][4]), Point(box[i][3], box[i][5]), Scalar(255, 0, 0), 2, 8);
//	}*/
//	//rot = mask_img;
//	
//
//	//Delete...
//	for (int i = 0; i < r; i++) {
//		delete[] orgimg[i];
//	}
//	delete[] orgimg;
//	return(N);
//}
//
//
////count the #text pixels in a label...
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
////count the #overlapping text pixels in 2 labels...
//int countCommonTextPixel(Mat orgimg, int **lbl1, int l1, int **lbl2, int l2, int *dim1, int *dim2) {
//	int count = 0;
//	Point limit_l, limit_u;
//	limit_l.x = dim1[0] > dim2[0] ? dim1[0] : dim2[0]; 
//	limit_l.y = dim1[2] > dim2[2] ? dim1[2] : dim2[2];
//	limit_u.x = dim1[1] < dim2[1] ? dim1[1] : dim2[1];
//	limit_u.y = dim1[3] < dim2[3] ? dim1[3] : dim2[3];
//	for (int i = limit_l.y; i <= limit_u.y; i++) {
//		for (int j = limit_l.x; j <= limit_u.x; j++) {
//			if (orgimg.at<uchar>(i, j) == 0 && lbl1[i][j] == l1 && lbl2[i][j] == l2)
//				count++;
//		}
//	}
//	return(count);
//}
//
//
//void regionMask1(Mat orgimg) {
//	Mat S1 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/-40Mask.tif", IMREAD_GRAYSCALE);
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
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasks.tif", allMask);
//	imshow("allMasks", allMaskCol);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasksColored.tif", allMaskCol);
//
//}
//
//void regionMask2(Mat orgimg) {
//	Mat S1 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/-45Mask.tif", IMREAD_GRAYSCALE);
//	Mat S2 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/0Mask.tif", IMREAD_GRAYSCALE);
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
//			if ((p1 && !p2) || (!p1 && p2)) {
//				continue;
//			}
//
//			//s1 and s2 overlapped...
//			if (p1 && p2) {
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
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasks.tif", allMask);
//	imshow("allMasks", allMaskCol);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasksColored.tif", allMaskCol);
//
//}
//
//void regionMask3(Mat orgimg) {
//	Mat S1 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/0Mask.tif", IMREAD_GRAYSCALE);
//	Mat S2 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/-90Mask.tif", IMREAD_GRAYSCALE);
//	Mat S3 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/-25Mask.tif", IMREAD_GRAYSCALE);
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
//			if ((p1 && !p2 && !p3) || (!p1 && p2 && !p3) || (!p1 && !p2 && p3)) {
//				continue;
//			}
//
//			//s1 and s2 overlapped...
//			if (p1 && p2 && !p3) {
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
//			else if (p1 && !p2 && p3) {
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
//			else if (!p1 && p2 && p3) {
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
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasks.tif", allMask);
//	imshow("allMasks", allMaskCol);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasksColored.tif", allMaskCol);
//
//}
//
//void regionMask4(Mat orgimg) {
//	Mat S1 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/0Mask.tif", IMREAD_GRAYSCALE);
//	Mat S2 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/15Mask.tif", IMREAD_GRAYSCALE);
//	Mat S3 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/-30Mask.tif", IMREAD_GRAYSCALE);
//	Mat S4 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/-45Mask.tif", IMREAD_GRAYSCALE);
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
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasks.tif", allMask);
//	imshow("allMasks", allMaskCol);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasksColored.tif", allMaskCol);
//
//}
//
//void regionMask1_1(Mat orgimg) {
//	Mat S1 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/-40Mask.tif", IMREAD_GRAYSCALE);
//
//	int **s1;
//	int **org;
//	int **labelorg;
//	int labelNum = 0;
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
//	labelorg = new int*[r];
//	for (int i = 0; i < r; i++) {
//		labelorg[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			labelorg[i][j] = 0;
//		}
//	}
//	labelNum = labelling(org, labelorg, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	int **box1 = new int*[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		box1[i] = new int[6];
//	boundingRectangles(box1, labelorg, r, c, labelNum);
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
//	int *label_color = new int[labelNum];
//	for (int i = 0; i < labelNum; i++)
//		label_color[i] = 0;
//
//	arr2mat(all_mask, allMask);
//	cvtColor(allMask, allMaskCol, CV_GRAY2BGR);
//	Vec3b color1 = { 255, 0, 255 };
//	Vec3b color2 = { 0, 153, 0 };
//	Vec3b color3 = { 255, 0, 0 };
//	Vec3b color4 = { 0, 120, 255 };
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (all_mask[i][j] == 150 && s1[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color1;
//				if (label_color[labelorg[i][j]] == 0)
//					label_color[labelorg[i][j]] = 1;
//			}
//		}
//	}
//
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (all_mask[i][j] == 150 && s1[i][j] != 0 && label_color[labelorg[i][j]] != 0) {
//				if (label_color[labelorg[i][j]] == 1)
//					allMaskCol.at<Vec3b>(i, j) = color1;
//			}
//
//		}
//	}
//
//	imshow("allMasks", allMask);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasks.tif", allMask);
//	imshow("allMasks", allMaskCol);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasksColored.tif", allMaskCol);
//
//}
//
//void regionMask2_1(Mat orgimg) {
//	Mat S1 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/15Mask.tif", IMREAD_GRAYSCALE);
//	Mat S2 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/-25Mask.tif", IMREAD_GRAYSCALE);
//
//	int **s1;
//	int **s2;
//	int **org;
//
//	int **label1;
//	int **label2;
//	int **labelorg;
//
//	int labelNum1 = 0;
//	int labelNum2 = 0;
//	int labelNumorg = 0;
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
//	//labelling org...
//	labelorg = new int*[r];
//	for (int i = 0; i < r; i++) {
//		labelorg[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			labelorg[i][j] = 0;
//		}
//	}
//	labelNumorg = labelling(org, labelorg, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	int **boxorg = new int*[labelNumorg];
//	for (int i = 0; i < labelNumorg; i++)
//		boxorg[i] = new int[6];
//	boundingRectangles(boxorg, labelorg, r, c, labelNumorg);
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
//	std::vector<pair<int, int>> comparedLabels;
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
//			if ((p1 && !p2) || (!p1 && p2)) {
//				continue;
//			}
//
//			int l1 = label1[i][j];
//			int l2 = label2[i][j];
//
//			int flag = 0;
//			for (int k = 0; k<int(comparedLabels.size()); k++) {
//				if ((comparedLabels[k].first == l1 && comparedLabels[k].second == l2) || (comparedLabels[k].first == l2 && comparedLabels[k].second == l1)) {
//					flag = 1;
//					break;
//				}
//			}
//			if (flag)
//				continue;
//
//			pair<int, int> new_pair(l1, l2);
//			comparedLabels.push_back(new_pair);
//
//			//s1 and s2 overlapped...
//			if (p1 && p2) {
//				double diag1 = (box1[l1][5] - box1[l1][4])*(box1[l1][5] - box1[l1][4]) + (box1[l1][3] - box1[l1][2])*(box1[l1][3] - box1[l1][2]);
//				double diag2 = (box2[l2][5] - box2[l2][4])*(box2[l2][5] - box2[l2][4]) + (box2[l2][3] - box2[l2][2])*(box2[l2][3] - box2[l2][2]);
//				std::cout << "\nDiag1 = " << sqrt(diag1) << " , Diag2 = " << sqrt(diag2);
//				
//				int *dim1 = new int[4];
//				int *dim2 = new int[4];
//				for (int k = 0; k < 4; k++) {
//					dim1[k] = box1[l1][k + 2];
//					dim2[k] = box2[l2][k + 2];
//				}
//				int c1 = countTextPixel(orgimg, label1, l1, (box1[l1][5] - box1[l1][4]), (box1[l1][3] - box1[l1][2]), Point(box1[l1][2], box1[l1][4]));
//				int c2 = countTextPixel(orgimg, label2, l2, (box2[l2][5] - box2[l2][4]), (box2[l2][3] - box2[l2][2]), Point(box2[l2][2], box2[l2][4]));
//				int c3 = countCommonTextPixel(orgimg, label1, l1, label2, l2, dim1, dim2);
//				std::cout << "\ncount1 = " << c1 << " , count2 = " << c2 << " , count3(overlap) = " << c3;
//				if (c1 < c2 && c3 / double(c1) > 0.3) {
//					std::cout << "\nRemove cc from s1";
//					for (int k = box1[l1][4]; k <= box1[l1][5]; k++) {
//						for (int l = box1[l1][2]; l <= box1[l1][3]; l++)
//							if (k < r && l < c && label1[k][l] == l1)
//								s1[k][l] = 255;
//					}
//				}
//				else if (c2 < c1 && c3 / double(c2) > 0.3) {
//					std::cout << "\nRemove cc from s2";
//					for (int k = box2[l2][4]; k <= box2[l2][5]; k++) {
//						for (int l = box2[l2][2]; l <= box2[l2][3]; l++)
//							if (k < r && l < c && label2[k][l] == l2)
//								s2[k][l] = 255;
//					}
//				}
//				else
//					std::cout << "\ns1 and s2 overlap...";
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
//	int *label_color = new int[labelNumorg];
//	for (int i = 0; i < labelNumorg; i++)
//		label_color[i] = 0;
//
//	arr2mat(all_mask, allMask);
//	cvtColor(allMask, allMaskCol, CV_GRAY2BGR);
//	Vec3b color0 = { 150, 150, 150 };
//	Vec3b color1 = { 255, 0, 255 };
//	Vec3b color2 = { 0, 153, 0 };
//	Vec3b color3 = { 255, 0, 0 };
//	Vec3b color4 = { 0, 120, 255 };
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (all_mask[i][j] == 150 && s1[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color1;
//				if (label_color[labelorg[i][j]] == 0)
//					label_color[labelorg[i][j]] = 1;
//				else if (label_color[labelorg[i][j]] != 1) {
//					label_color[labelorg[i][j]] = -1;
//				}
//			}
//			else if (all_mask[i][j] == 150 && s2[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color2;
//				if (label_color[labelorg[i][j]] == 0)
//					label_color[labelorg[i][j]] = 2;
//				else if (label_color[labelorg[i][j]] != 2) {
//					label_color[labelorg[i][j]] = -1;
//				}
//			}
//		}
//	}
//
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (label_color[labelorg[i][j]] == -1)
//				allMaskCol.at<Vec3b>(i, j) = color0;
//			if (label_color[labelorg[i][j]] == 1)
//				allMaskCol.at<Vec3b>(i, j) = color1;
//			if (label_color[labelorg[i][j]] == 2)
//				allMaskCol.at<Vec3b>(i, j) = color2;
//		}
//	}
//
//	imshow("allMasks", allMask);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasks.tif", allMask);
//	imshow("allMasks", allMaskCol);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasksColored.tif", allMaskCol);
//
//}
//
//void regionMask3_1(Mat orgimg) {
//	Mat S1 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/50Mask.tif", IMREAD_GRAYSCALE);
//	Mat S2 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/-10Mask.tif", IMREAD_GRAYSCALE);
//	Mat S3 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/-40Mask.tif", IMREAD_GRAYSCALE);
//
//	int **s1;
//	int **s2;
//	int **s3;
//	int **org;
//
//	int **label1;
//	int **label2;
//	int **label3;
//	int **labelorg;
//
//	int labelNum1 = 0;
//	int labelNum2 = 0;
//	int labelNum3 = 0;
//	int labelNumorg = 0;
//
//	Mat allMask, allMaskCol;
//	allMask.create(r, c, CV_8UC1);
//	allMaskCol.create(r, c, CV_8UC1);
//	int **all_mask;
//
//	//label original...
//	org = new int*[r];
//	for (int i = 0; i < r; i++) {
//		org[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			org[i][j] = 0;
//		}
//	}
//	mat2arr(orgimg, org);
//
//	labelorg = new int*[r];
//	for (int i = 0; i < r; i++) {
//		labelorg[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			labelorg[i][j] = 0;
//		}
//	}
//	labelNumorg = labelling(org, labelorg, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	int **boxorg = new int*[labelNumorg];
//	for (int i = 0; i < labelNumorg; i++)
//		boxorg[i] = new int[6];
//	boundingRectangles(boxorg, labelorg, r, c, labelNumorg);
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
//	std::vector<pair<int, int>> comparedLabels12;
//	std::vector<pair<int, int>> comparedLabels23;
//	std::vector<pair<int, int>> comparedLabels13;
//
//	//compare overlapped regions...
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			int p1 = (255 - s1[i][j]) / 255;
//			int p2 = (255 - s2[i][j]) / 255;
//			int p3 = (255 - s3[i][j]) / 255;
//
//			int l1 = label1[i][j];
//			int l2 = label2[i][j];
//			int l3 = label3[i][j];
//
//			//no component...
//			if (p1 + p2 + p3 == 0)
//				continue;
//			//no overlapping component...
//			if ((p1 && !p2 && !p3) || (!p1 && p2 && !p3) || (!p1 && !p2 && p3)) {
//				continue;
//			}
//
//			//s1 and s2 overlapped...
//			if (p1 && p2 && !p3) {
//				int flag = 0;
//				for (int k = 0; k<int(comparedLabels12.size()); k++) {
//					if ((comparedLabels12[k].first == l1 && comparedLabels12[k].second == l2) || (comparedLabels12[k].first == l2 && comparedLabels12[k].second == l1)) {
//						flag = 1;
//						break;
//					}
//				}
//				if (flag)
//					continue;
//
//				pair<int, int> new_pair(l1, l2);
//				comparedLabels12.push_back(new_pair);
//
//				double diag1 = (box1[l1][5] - box1[l1][4])*(box1[l1][5] - box1[l1][4]) + (box1[l1][3] - box1[l1][2])*(box1[l1][3] - box1[l1][2]);
//				double diag2 = (box2[l2][5] - box2[l2][4])*(box2[l2][5] - box2[l2][4]) + (box2[l2][3] - box2[l2][2])*(box2[l2][3] - box2[l2][2]);
//				
//				int *dim1 = new int[4];
//				int *dim2 = new int[4];
//				for (int k = 0; k < 4; k++) {
//					dim1[k] = box1[l1][k + 2];
//					dim2[k] = box2[l2][k + 2];
//				}
//
//				int c1 = countTextPixel(orgimg, label1, l1, (box1[l1][5] - box1[l1][4]), (box1[l1][3] - box1[l1][2]), Point(box1[l1][2], box1[l1][4]));
//				int c2 = countTextPixel(orgimg, label2, l2, (box2[l2][5] - box2[l2][4]), (box2[l2][3] - box2[l2][2]), Point(box2[l2][2], box2[l2][4]));
//				int c3 = countCommonTextPixel(orgimg, label1, l1, label2, l2, dim1, dim2);
//				std::cout << "\ncount1 = " << c1 << " , count2 = " << c2 << " , count3(overlap) = " << c3;
//				if (c1 < c2 && c3 / double(c1) > 0.3) {
//					std::cout << "\nRemove cc from s1";
//					for (int k = box1[l1][4]; k <= box1[l1][5]; k++) {
//						for (int l = box1[l1][2]; l <= box1[l1][3]; l++)
//							if (k < r && l < c && label1[k][l] == l1)
//								s1[k][l] = 255;
//					}
//				}
//				else if (c2 < c1 && c3 / double(c2) > 0.3) {
//					std::cout << "\nRemove cc from s2";
//					for (int k = box2[l2][4]; k <= box2[l2][5]; k++) {
//						for (int l = box2[l2][2]; l <= box2[l2][3]; l++)
//							if (k < r && l < c && label2[k][l] == l2)
//								s2[k][l] = 255;
//					}
//				}
//				else
//					std::cout << "\ns1 and s2 overlap...";
//			}
//			
//			//s1 and s3 overlapped...
//			else if (p1 && !p2 && p3) {
//				int flag = 0;
//				for (int k = 0; k<int(comparedLabels13.size()); k++) {
//					if ((comparedLabels13[k].first == l1 && comparedLabels13[k].second == l3) || (comparedLabels13[k].first == l3 && comparedLabels13[k].second == l1)) {
//						flag = 1;
//						break;
//					}
//				}
//				if (flag)
//					continue;
//
//				pair<int, int> new_pair(l1, l3);
//				comparedLabels13.push_back(new_pair);
//
//				double diag1 = (box1[l1][5] - box1[l1][4])*(box1[l1][5] - box1[l1][4]) + (box1[l1][3] - box1[l1][2])*(box1[l1][3] - box1[l1][2]);
//				double diag3 = (box3[l3][5] - box3[l3][4])*(box3[l3][5] - box3[l3][4]) + (box3[l3][3] - box3[l3][2])*(box3[l3][3] - box3[l3][2]);
//				std::cout << "\nDiag1 = " << sqrt(diag1) << " , Diag3 = " << sqrt(diag3);
//				
//				int *dim1 = new int[4];
//				int *dim3 = new int[4];
//				for (int k = 0; k < 4; k++) {
//					dim1[k] = box1[l1][k + 2];
//					dim3[k] = box3[l3][k + 2];
//				}
//
//				int c1 = countTextPixel(orgimg, label1, l1, (box1[l1][5] - box1[l1][4]), (box1[l1][3] - box1[l1][2]), Point(box1[l1][2], box1[l1][4]));
//				int c2 = countTextPixel(orgimg, label3, l3, (box3[l3][5] - box3[l3][4]), (box3[l3][3] - box3[l3][2]), Point(box3[l3][2], box3[l3][4]));
//				int c3 = countCommonTextPixel(orgimg, label1, l1, label3, l3, dim1, dim3);
//				std::cout << "\ncount1 = " << c1 << " , count2 = " << c2 << " , count3(overlap) = " << c3;
//				if (c1 < c2 && c3 / double(c1) > 0.3) {
//					std::cout << "\nRemove cc from s1";
//					for (int k = box1[l1][4]; k <= box1[l1][5]; k++) {
//						for (int l = box1[l1][2]; l <= box1[l1][3]; l++)
//							if (k < r && l < c && label1[k][l] == l1)
//								s1[k][l] = 255;
//					}
//				}
//				else if (c2 < c1 && c3 / double(c2) > 0.3) {
//					std::cout << "\nRemove cc from s3";
//					for (int k = box3[l3][4]; k <= box3[l3][5]; k++) {
//						for (int l = box3[l3][2]; l <= box3[l3][3]; l++)
//							if (k < r && l < c && label3[k][l] == l3)
//								s3[k][l] = 255;
//					}
//				}
//				else
//					std::cout << "\ns1 and s3 overlap...";
//			}
//			//s2 and s3 overlapped...
//			else if (!p1 && p2 && p3) {
//				int flag = 0;
//				for (int k = 0; k<int(comparedLabels23.size()); k++) {
//					if ((comparedLabels23[k].first == l3 && comparedLabels23[k].second == l2) || (comparedLabels23[k].first == l2 && comparedLabels23[k].second == l3)) {
//						flag = 1;
//						break;
//					}
//				}
//				if (flag)
//					continue;
//
//				pair<int, int> new_pair(l2, l3);
//				comparedLabels23.push_back(new_pair);
//
//				double diag2 = (box2[l2][5] - box2[l2][4])*(box2[l2][5] - box2[l2][4]) + (box2[l2][3] - box2[l2][2])*(box2[l2][3] - box2[l2][2]);
//				double diag3 = (box3[l3][5] - box3[l3][4])*(box3[l3][5] - box3[l3][4]) + (box3[l3][3] - box3[l3][2])*(box3[l3][3] - box3[l3][2]);
//				std::cout << "\nDiag2 = " << sqrt(diag2) << " , Diag3 = " << sqrt(diag3);
//				
//				int *dim2 = new int[4];
//				int *dim3 = new int[4];
//				for (int k = 0; k < 4; k++) {
//					dim2[k] = box2[l2][k + 2];
//					dim3[k] = box3[l3][k + 2];
//				}
//
//				int c1 = countTextPixel(orgimg, label2, l2, (box2[l2][5] - box2[l2][4]), (box2[l2][3] - box2[l2][2]), Point(box2[l2][2], box2[l2][4]));
//				int c2 = countTextPixel(orgimg, label3, l3, (box3[l3][5] - box3[l3][4]), (box3[l3][3] - box3[l3][2]), Point(box3[l3][2], box3[l3][4]));
//				int c3 = countCommonTextPixel(orgimg, label2, l2, label3, l3, dim2, dim3);
//				std::cout << "\ncount1 = " << c1 << " , count2 = " << c2 << " , count3(overlap) = " << c3;
//				if (c1 < c2 && c3 / double(c1) > 0.3) {
//					std::cout << "\nRemove cc from s2";
//					for (int k = box2[l2][4]; k <= box2[l2][5]; k++) {
//						for (int l = box2[l2][2]; l <= box2[l2][3]; l++)
//							if (k < r && l < c && label2[k][l] == l2)
//								s2[k][l] = 255;
//					}
//				}
//				else if (c2 < c1 && c3 / double(c2) > 0.3) {
//					std::cout << "\nRemove cc from s3";
//					for (int k = box3[l3][4]; k <= box3[l3][5]; k++) {
//						for (int l = box3[l3][2]; l <= box3[l3][3]; l++)
//							if (k < r && l < c && label3[k][l] == l3)
//								s3[k][l] = 255;
//					}
//				}
//				else
//					std::cout << "\ns2 and s3 overlap...";
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
//	int *label_color = new int[labelNumorg];
//	for (int i = 0; i < labelNumorg; i++)
//		label_color[i] = 0;
//
//	arr2mat(all_mask, allMask);
//	cvtColor(allMask, allMaskCol, CV_GRAY2BGR);
//	Vec3b color1 = { 255, 0, 255 };
//	Vec3b color2 = { 0, 153, 0 };
//	Vec3b color3 = { 255, 0, 0 };
//	Vec3b color4 = { 0, 120, 255 };
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (all_mask[i][j] == 150 && s1[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color1;
//				if (label_color[labelorg[i][j]] == 0)
//					label_color[labelorg[i][j]] = 1;
//			}
//			else if (all_mask[i][j] == 150 && s2[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color2;
//				if (label_color[labelorg[i][j]] == 0)
//					label_color[labelorg[i][j]] = 2;
//			}
//			else if (all_mask[i][j] == 150 && s3[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color3;
//				if (label_color[labelorg[i][j]] == 0)
//					label_color[labelorg[i][j]] = 3;
//			}
//		}
//	}
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (all_mask[i][j] == 150 && s1[i][j] != 0 && s2[i][j] != 0 && s3[i][j] != 0 && label_color[labelorg[i][j]] != 0) {
//				if(label_color[labelorg[i][j]] == 1)
//					allMaskCol.at<Vec3b>(i, j) = color1;
//				if (label_color[labelorg[i][j]] == 2)
//					allMaskCol.at<Vec3b>(i, j) = color2;
//				if (label_color[labelorg[i][j]] == 3)
//					allMaskCol.at<Vec3b>(i, j) = color3;
//			}
//
//		}
//	}
//
//	imshow("allMasks", allMask);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasks.tif", allMask);
//	imshow("allMasks", allMaskCol);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasksColored.tif", allMaskCol);
//
//}
//
//void regionMask4_1(Mat orgimg) {
//	Mat S1 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/0Mask.tif", IMREAD_GRAYSCALE);
//	Mat S2 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/30Mask.tif", IMREAD_GRAYSCALE);
//	Mat S3 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/90Mask.tif", IMREAD_GRAYSCALE);
//	Mat S4 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/-15Mask.tif", IMREAD_GRAYSCALE);
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
//	int **labelorg;
//
//	int labelNum1 = 0;
//	int labelNum2 = 0;
//	int labelNum3 = 0;
//	int labelNum4 = 0;
//	int labelNumorg = 0;
//
//	Mat allMask, allMaskCol;
//	allMask.create(r, c, CV_8UC1);
//	allMaskCol.create(r, c, CV_8UC1);
//	int **all_mask;
//
//	//label org...
//	org = new int*[r];
//	for (int i = 0; i < r; i++) {
//		org[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			org[i][j] = 0;
//		}
//	}
//	mat2arr(orgimg, org);
//
//	labelorg = new int*[r];
//	for (int i = 0; i < r; i++) {
//		labelorg[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			labelorg[i][j] = 0;
//		}
//	}
//	labelNumorg = labelling(org, labelorg, r, c);
//
//	//Bound the components by rectangular & diagonal boxes
//	int **boxorg = new int*[labelNumorg];
//	for (int i = 0; i < labelNumorg; i++)
//		boxorg[i] = new int[6];
//	boundingRectangles(boxorg, labelorg, r, c, labelNumorg);
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
//	int *label_color = new int[labelNumorg];
//	for (int i = 0; i < labelNumorg; i++)
//		label_color[i] = 0;
//
//	arr2mat(all_mask, allMask);
//	cvtColor(allMask, allMaskCol, CV_GRAY2BGR);
//	Vec3b color1 = { 255, 0, 255 };
//	Vec3b color2 = { 0, 153, 0 };
//	Vec3b color3 = { 255, 0, 0 };
//	Vec3b color4 = { 0, 120, 255 };
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (all_mask[i][j] == 150 && s1[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color1;
//				if (label_color[labelorg[i][j]] == 0)
//					label_color[labelorg[i][j]] = 1;
//			}
//			else if (all_mask[i][j] == 150 && s2[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color2;
//				if (label_color[labelorg[i][j]] == 0)
//					label_color[labelorg[i][j]] = 2;
//			}
//			else if (all_mask[i][j] == 150 && s3[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color3;
//				if (label_color[labelorg[i][j]] == 0)
//					label_color[labelorg[i][j]] = 3;
//			}
//			else if (all_mask[i][j] == 150 && s4[i][j] == 0) {
//				allMaskCol.at<Vec3b>(i, j) = color4;
//				if (label_color[labelorg[i][j]] == 0)
//					label_color[labelorg[i][j]] = 4;
//			}
//		}
//	}
//
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			if (all_mask[i][j] == 150 && s1[i][j] != 0 && s2[i][j] != 0 && s3[i][j] != 0 && label_color[labelorg[i][j]] != 0) {
//				if (label_color[labelorg[i][j]] == 1)
//					allMaskCol.at<Vec3b>(i, j) = color1;
//				if (label_color[labelorg[i][j]] == 2)
//					allMaskCol.at<Vec3b>(i, j) = color2;
//				if (label_color[labelorg[i][j]] == 3)
//					allMaskCol.at<Vec3b>(i, j) = color3;
//				if (label_color[labelorg[i][j]] == 4)
//					allMaskCol.at<Vec3b>(i, j) = color4;
//			}
//
//		}
//	}
//
//	imshow("allMasks", allMask);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasks.tif", allMask);
//	imshow("allMasks", allMaskCol);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/SuperimposeMasksColored.tif", allMaskCol);
//
//}
//
//
//int main(int argc, char** argv) {
//
//	//-----------------------------Read image-----------------------
//	//Mat scanimg = imread("Data/ICDAR/Data218/scan218.tif", IMREAD_GRAYSCALE);
//	//Mat scanimg = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/scan008.tif", IMREAD_GRAYSCALE);
//	//Mat scanimg = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/rotate15org.tif", IMREAD_GRAYSCALE);
//	//Mat orgimg = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/invccorg.tif", IMREAD_GRAYSCALE);
//
//	//Mat img0 = imread("Data/ICDAR/Data218/Gauss0.tif", IMREAD_GRAYSCALE);
//	//Mat img0 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/Gauss0.tif", IMREAD_GRAYSCALE);
//	//Mat img-30 = imread("Data/ICDAR/Data218/Gauss-30.tif", IMREAD_GRAYSCALE);
//	//Mat img-30 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/Gauss-30.tif", IMREAD_GRAYSCALE);
//	//Mat img135 = imread("Data/ICDAR/Data218/Gauss135.tif", IMREAD_GRAYSCALE);
//	//Mat img135 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/Gauss137.tif", IMREAD_GRAYSCALE);
//	//Mat scancol;
//
//	//imshow("Image", scanimg);
//	//std::cout << "Image dimensions:" << scanimg.rows << "x" << scanimg.cols << endl;
//
//	//r = scanimg.rows;
//	//c = scanimg.cols;
//
//	///*scancol.create(r, c, CV_8UC1);
//	//cvtColor(scanimg, scancol, CV_GRAY2BGR);*/
//	////invCC(scanimg, 1);
//	//invCCProj(scanimg, 1);
//
//	/*imshow("Image", orgimg);
//	std::cout << "Image dimensions:" << orgimg.rows << "x" << orgimg.cols << endl;
//
//	r = orgimg.rows;
//	c = orgimg.cols;
//	regionMask3_1(orgimg);*/
//
//	/*scancol.create(r, c, CV_8UC1);
//	cvtColor(img135, scancol, CV_GRAY2BGR);
//	trace(img135);*/
//
//	//-----------------------------------------------------------------------------------
//	int *longestWords = new int[18];
//
//	std::cout << "\nReading image1...";
//	Mat scan1 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate0org.tif", IMREAD_GRAYSCALE);
//	r = scan1.rows;
//	c = scan1.cols;
//	Mat rot1, mask_img1;
//	rot1.create(r, c, CV_8UC1);
//	cvtColor(scan1, rot1, CV_GRAY2BGR);
//	mask_img1.create(r, c, CV_8UC1);
//	longestWords[0] = validAngles(scan1, rot1, mask_img1);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot0.tif", rot1);
//
//	std::cout << "\nReading image2...";
//	Mat scan2 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate10org.tif", IMREAD_GRAYSCALE);
//	r = scan2.rows;
//	c = scan2.cols;
//	Mat rot2;
//	rot2.create(r, c, CV_8UC1);
//	cvtColor(scan2, rot2, CV_GRAY2BGR);
//	Mat mask_img2;
//	mask_img2.create(r, c, CV_8UC1);
//	longestWords[1] = validAngles(scan2, rot2, mask_img2);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot10.tif", rot2);
//
//	std::cout << "\nReading image3...";
//	Mat scan3 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate20org.tif", IMREAD_GRAYSCALE);
//	r = scan3.rows;
//	c = scan3.cols;
//	Mat rot3;
//	rot3.create(r, c, CV_8UC1);
//	cvtColor(scan3, rot3, CV_GRAY2BGR);
//	Mat mask_img3;
//	mask_img3.create(r, c, CV_8UC1);
//	longestWords[2] = validAngles(scan3, rot3, mask_img3);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot20.tif", rot3);
//
//	std::cout << "\nReading image4...";
//	Mat scan4 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate30org.tif", IMREAD_GRAYSCALE);
//	r = scan4.rows;
//	c = scan4.cols;
//	Mat rot4;
//	rot4.create(r, c, CV_8UC1);
//	cvtColor(scan4, rot4, CV_GRAY2BGR);
//	Mat mask_img4;
//	mask_img4.create(r, c, CV_8UC1);
//	longestWords[3] = validAngles(scan4, rot4, mask_img4);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot30.tif", rot4);
//
//	std::cout << "\nReading image5...";
//	Mat scan5 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate40org.tif", IMREAD_GRAYSCALE);
//	r = scan5.rows;
//	c = scan5.cols;
//	Mat rot5;
//	rot5.create(r, c, CV_8UC1);
//	cvtColor(scan5, rot5, CV_GRAY2BGR);
//	Mat mask_img5;
//	mask_img5.create(r, c, CV_8UC1);
//	longestWords[4] = validAngles(scan5, rot5, mask_img5);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot40.tif", rot5);
//
//	std::cout << "\nReading image6...";
//	Mat scan6 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate50org.tif", IMREAD_GRAYSCALE);
//	r = scan6.rows;
//	c = scan6.cols;
//	Mat rot6;
//	rot6.create(r, c, CV_8UC1);
//	cvtColor(scan6, rot6, CV_GRAY2BGR);
//	Mat mask_img6;
//	mask_img6.create(r, c, CV_8UC1);
//	longestWords[5] = validAngles(scan6, rot6, mask_img6);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot50.tif", rot6);
//
//	std::cout << "\nReading image7...";
//	Mat scan7 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate60org.tif", IMREAD_GRAYSCALE);
//	r = scan7.rows;
//	c = scan7.cols;
//	Mat rot7;
//	rot7.create(r, c, CV_8UC1);
//	cvtColor(scan7, rot7, CV_GRAY2BGR);
//	Mat mask_img7;
//	mask_img7.create(r, c, CV_8UC1);
//	longestWords[6] = validAngles(scan7, rot7, mask_img7);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot60.tif", rot7);
//
//	std::cout << "\nReading image8...";
//	Mat scan8 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate70org.tif", IMREAD_GRAYSCALE);
//	r = scan8.rows;
//	c = scan8.cols;
//	Mat rot8;
//	rot8.create(r, c, CV_8UC1);
//	cvtColor(scan8, rot8, CV_GRAY2BGR);
//	Mat mask_img8;
//	mask_img8.create(r, c, CV_8UC1);
//	longestWords[7] = validAngles(scan8, rot8, mask_img8);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot70.tif", rot8);
//
//	std::cout << "\nReading image9...";
//	Mat scan9 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate80org.tif", IMREAD_GRAYSCALE);
//	r = scan9.rows;
//	c = scan9.cols;
//	Mat rot9;
//	rot9.create(r, c, CV_8UC1);
//	cvtColor(scan9, rot9, CV_GRAY2BGR);
//	Mat mask_img9;
//	mask_img9.create(r, c, CV_8UC1);
//	longestWords[8] = validAngles(scan9, rot9, mask_img9);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot80.tif", rot9);
//
//	std::cout << "\nReading image10...";
//	Mat scan10 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate90org.tif", IMREAD_GRAYSCALE);
//	r = scan10.rows;
//	c = scan10.cols;
//	Mat rot10;
//	rot10.create(r, c, CV_8UC1);
//	cvtColor(scan10, rot10, CV_GRAY2BGR);
//	Mat mask_img10;
//	mask_img10.create(r, c, CV_8UC1);
//	longestWords[9] = validAngles(scan10, rot10, mask_img10);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot90.tif", rot10);
//
//	std::cout << "\nReading image11...";
//	Mat scan11 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate-10org.tif", IMREAD_GRAYSCALE);
//	r = scan11.rows;
//	c = scan11.cols;
//	Mat rot11;
//	rot11.create(r, c, CV_8UC1);
//	cvtColor(scan11, rot11, CV_GRAY2BGR);
//	Mat mask_img11;
//	mask_img11.create(r, c, CV_8UC1);
//	longestWords[10] = validAngles(scan11, rot11, mask_img11);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot-10.tif", rot11);
//
//	std::cout << "\nReading image12...";
//	Mat scan12 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate-20org.tif", IMREAD_GRAYSCALE);
//	r = scan12.rows;
//	c = scan12.cols;
//	Mat rot12;
//	rot12.create(r, c, CV_8UC1);
//	cvtColor(scan12, rot12, CV_GRAY2BGR);
//	Mat mask_img12;
//	mask_img12.create(r, c, CV_8UC1);
//	longestWords[11] = validAngles(scan12, rot12, mask_img12);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot-20.tif", rot12);
//
//	std::cout << "\nReading image13...";
//	Mat scan13 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate-30org.tif", IMREAD_GRAYSCALE);
//	r = scan13.rows;
//	c = scan13.cols;
//	Mat rot13;
//	rot13.create(r, c, CV_8UC1);
//	cvtColor(scan13, rot13, CV_GRAY2BGR);
//	Mat mask_img13;
//	mask_img13.create(r, c, CV_8UC1);
//	longestWords[12] = validAngles(scan13, rot13, mask_img13);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot-30.tif", rot13);
//
//	std::cout << "\nReading image14...";
//	Mat scan14 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate-40org.tif", IMREAD_GRAYSCALE);
//	r = scan14.rows;
//	c = scan14.cols;
//	Mat rot14;
//	rot14.create(r, c, CV_8UC1);
//	cvtColor(scan14, rot14, CV_GRAY2BGR);
//	Mat mask_img14;
//	mask_img14.create(r, c, CV_8UC1);
//	longestWords[13] = validAngles(scan14, rot14, mask_img14);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot-40.tif", rot14);
//
//	std::cout << "\nReading image15...";
//	Mat scan15 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate-50org.tif", IMREAD_GRAYSCALE);
//	r = scan15.rows;
//	c = scan15.cols;
//	Mat rot15;
//	rot15.create(r, c, CV_8UC1);
//	cvtColor(scan15, rot15, CV_GRAY2BGR);
//	Mat mask_img15;
//	mask_img15.create(r, c, CV_8UC1);
//	longestWords[14] = validAngles(scan15, rot15, mask_img15);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot-50.tif", rot15);
//
//	std::cout << "\nReading image16...";
//	Mat scan16 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate-60org.tif", IMREAD_GRAYSCALE);
//	r = scan16.rows;
//	c = scan16.cols;
//	Mat rot16;
//	rot16.create(r, c, CV_8UC1);
//	cvtColor(scan16, rot16, CV_GRAY2BGR);
//	Mat mask_img16;
//	mask_img16.create(r, c, CV_8UC1);
//	longestWords[15] = validAngles(scan16, rot16, mask_img16);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot-60.tif", rot16);
//
//	std::cout << "\nReading image17...";
//	Mat scan17 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate-70org.tif", IMREAD_GRAYSCALE);
//	r = scan17.rows;
//	c = scan17.cols;
//	Mat rot17;
//	rot17.create(r, c, CV_8UC1);
//	cvtColor(scan17, rot17, CV_GRAY2BGR);
//	Mat mask_img17;
//	mask_img17.create(r, c, CV_8UC1);
//	longestWords[16] = validAngles(scan17, rot17, mask_img17);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot-70.tif", rot17);
//
//	std::cout << "\nReading image18...";
//	Mat scan18 = imread("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/rotate-80org.tif", IMREAD_GRAYSCALE);
//	r = scan18.rows;
//	c = scan18.cols;
//	Mat rot18;
//	rot18.create(r, c, CV_8UC1);
//	cvtColor(scan18, rot18, CV_GRAY2BGR);
//	Mat mask_img18;
//	mask_img18.create(r, c, CV_8UC1);
//	longestWords[17] = validAngles(scan18, rot18, mask_img18);
//	imwrite("Data/MyScans/Scanned Data/Scanned Data/Scan008/AllRotations/binImgProjRot-80.tif", rot18);
//
//	for (int i = 9; i >= 0; i--)
//		std::cout << "\nRotate by " << i*10 << " degrees = " << longestWords[i];
//	for (int i = 10; i < 18; i++)
//		std::cout << "\nRotate by " << 90 - i*10 << " degrees = " << longestWords[i];
//
//	std::cout << "\nDone!\n";
//	waitKey(0);
//	system("pause");
//	return 0;
//
//}
