//#include <stdio.h>
//#include<conio.h>
//#include<iostream>
//#include "opencv2/imgproc/imgproc.hpp"
//#include "opencv2/highgui/highgui.hpp"
//#include "NICK.h"
//#include "Otsu.h"
//#include "RegLine.h"
//#include "ConnectedComponent.h"
//#include "BoundingBox.h"
//
//using namespace std;
//using namespace cv;
//
//int labelNum = 0, r, c;
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
//void arr2mat(int **arr, Mat A, int row, int col) {
//	for (int i = 0; i < row; i++) {
//		for (int j = 0; j < col; j++) {
//			A.at<uchar>(i, j) = arr[i][j];
//		}
//	}
//}
//void invImage(int **arr, int **inv_arr) {
//	for (int i = 0; i < r; i++) {
//		for (int j = 0; j < c; j++) {
//			inv_arr[i][j] = 255 - arr[i][j];
//		}
//	}
//}
//
///*----------------------------------------------------------------*/
///*                    Invert and fill holes                       */
///*----------------------------------------------------------------*/
//void fill(Mat img) {
//
//	Mat img1, bin_img, cc_img, main_img, inv_img, fill_img;
//	int **image, **bin_im1, **label, **box, **mainCC, *smallCC, **inv_im, *heights;
//
//	//----------------------Blur & Binarize image-----------------------
//	image = new int*[r];
//	bin_im1 = new int*[r];
//	inv_im = new int*[r];
//	for (int i = 0; i < r; i++) {
//		image[i] = new int[c];
//		bin_im1[i] = new int[c];
//		inv_im[i] = new int[c];
//		for (int j = 0; j < c; j++) {
//			image[i][j] = 0;
//			bin_im1[i][j] = 0;
//			inv_im[i][j] = 0;
//		}
//	}
//	mat2arr(img, image);
//	NICK(image, r, c, bin_im1);
//	bin_img.create(r, c, CV_8UC1);
//	arr2mat(bin_im1, bin_img);
//	imshow("binary", bin_img);
//	//imwrite("Data/Data29/OriginalBinary.tif", bin_img);
//
//	//----------------------Invert image-----------------------
//	invImage(bin_im1, inv_im);
//	inv_img.create(r, c, CV_8UC1);
//	arr2mat(inv_im, inv_img);
//	imshow("inverted", inv_img);
//
//	//--------------------------Fill holes-------------------------
//	fill_img = inv_img.clone();
//
//	vector<vector<Point>> contours;
//
//	findContours(fill_img, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
//	drawContours(fill_img, contours, -1, Scalar(255), CV_FILLED);
//
//	imshow("filled",fill_img);
//	mat2arr(fill_img, inv_im);
//	invImage(inv_im, bin_im1);
//	arr2mat(bin_im1, img);
//	imshow("filled", img);
//
//}
//
///*----------------------------------------------------------------*/
///*                        Main function                           */
///*----------------------------------------------------------------*/
//int main(int argc, char** argv) {
//
//	//-----------------------------Read image-----------------------
//	Mat img = imread("Data/Data29/scan29.jpg", IMREAD_GRAYSCALE);
//	imshow("Image", img);
//	std::cout << "Image dimensions:" << img.rows << "x" << img.cols << endl;
//
//	GaussianBlur(img, img, Size(3, 3), 0, 0);
//	r = img.rows;
//	c = img.cols;
//
//	fill(img);
//
//	std::cout << "\nDone!\n";
//	waitKey(0);
//	system("pause");
//	return 0;
//
//}