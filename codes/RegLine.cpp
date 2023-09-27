#include <stdio.h>
#include<conio.h>
#include<iostream>
#include "opencv2/imgproc/imgproc.hpp"
//#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "ConnectedComponent.h"

using namespace std;
using namespace cv;

double Mean(int *X, int N) {
	double avg = 0;
	for (int i = 0; i < N; i++) {
		avg += (double)X[i];
	}
	avg /= N;
	return(avg);
}

void regressionLine(int *X, int *Y, int N, double& slopeYX, double& slopeXY, double& interceptyx, double& interceptxy) {
	int meanX, meanY, meanXY, meanX2, meanY2, *XY, *X2, *Y2;
	double SumofSq = 0;
	XY = new int[N];
	X2 = new int[N];
	Y2 = new int[N];
	for (int i = 0; i < N; i++) {
		XY[i] = X[i] * Y[i];
		X2[i] = X[i] * X[i];
		Y2[i] = Y[i] * Y[i];
	}
	meanX = Mean(X, N);
	meanY = Mean(Y, N);
	meanXY = Mean(XY, N);
	meanX2 = Mean(X2, N);
	meanY2 = Mean(Y2, N);
	std::cout << "\nX' = " << meanX << " , Y' = " << meanY << " , XY = " << meanXY << " , X2 = " << meanX2;

	slopeYX = (meanX*meanY - meanXY) / (double)(meanX*meanX - meanX2);
	interceptyx = meanY - slopeYX*meanX;
	slopeXY = (meanX*meanY - meanXY) / (double)(meanY*meanY - meanY2);
	interceptxy = meanX - meanY * slopeXY;
	//std::cout << "\nSlopeYX = " << slopeYX << "\nSlopeXY = " << slopeXY;
}