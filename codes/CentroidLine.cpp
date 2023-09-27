#include <stdio.h>
#include<conio.h>
#include<iostream>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "VarDeclare.h"
#include "CentroidLine.h"

using namespace std;
using namespace cv;
namespace cv {
	using std::vector;
};

Point LocalCentroidDiag(int **img, int **label, int lbl, Point N, Point S, Point W, Point E) {
	Point cntr = Point(0, 0);
	int count = 0;
	/*if (N.y < 0) N.y = 0;
	if (S.y >= r) S.y = r - 1;*/

	int xmin = N.x, xmax = N.x;
	for (int i = N.y; i <= S.y; i++) {
		if (i < 0 || i >= r) {
			if (i <= W.y) xmin--;
			else xmin++;

			if (i <= E.y) xmax++;
			else xmax--;
		}
		else {
			if (xmin < 0) xmin = 0;
			if (xmax >= c) xmax = c - 1;
			for (int j = xmin; j <= xmax; j++) {
				if (j < 0 || j >= c) continue;
				if (label[i][j] == lbl) {
					cntr.x += j;
					cntr.y += i;
					count++;
				}
			}
			if (i <= W.y) xmin--;
			else xmin++;

			if (i <= E.y) xmax++;
			else xmax--;
		}
	}
	cntr.x = int(cntr.x / count);
	cntr.y = int(cntr.y / count);

	std::cout << "\nCentroids: (" << cntr.x << " , " << cntr.y << ")";
	return(cntr);
}

Point LocalCentroidUpright(int **img, int **label, int lbl, Point TL, Point TR, Point DR, Point DL) {
	Point cntr = Point(0, 0);
	int count = 0;

	for (int i = TL.y; i <= DL.y; i++) {
		for (int j = TL.x; j <= TR.x; j++) {
			if (label[i][j] == lbl) {
				cntr.x += j;
				cntr.y += i;
				count++;
			}
		}
	}
	cntr.x = int(cntr.x / count);
	cntr.y = int(cntr.y / count);

	return(cntr);
}

void Centroidline(int **img, int **label, int lbl, double *diagbox, int *box, Point ends[4][2]) {
	Point N, S, E, W, TL, TR, DL, DR;
	int factor = 12;

	//Draw horizontal line...
	TL = Point(box[2], box[4]);
	TR = Point(box[2] + factor, box[4]);
	DL = Point(box[2], box[5]);
	DR = Point(box[2] + factor, box[5]);
	ends[0][0] = LocalCentroidUpright(img, label, lbl, TL, TR, DL, DR);

	TL = Point(box[3] - factor, box[4]);
	TR = Point(box[3], box[4]);
	DL = Point(box[3] - factor, box[5]);
	DR = Point(box[3], box[5]);
	ends[0][1] = LocalCentroidUpright(img, label, lbl, TL, TR, DL, DR);

	//Draw vertical line...
	TL = Point(box[2], box[4]);
	TR = Point(box[3], box[4]);
	DL = Point(box[2], box[4] + factor);
	DR = Point(box[3], box[4] + factor);
	ends[1][0] = LocalCentroidUpright(img, label, lbl, TL, TR, DL, DR);

	TL = Point(box[2], box[5] - factor);
	TR = Point(box[3], box[5] - factor);
	DL = Point(box[2], box[5]);
	DR = Point(box[3], box[5]);
	ends[1][1] = LocalCentroidUpright(img, label, lbl, TL, TR, DL, DR);

	//Draw diag line...
	N = Point(diagbox[1], diagbox[0]);
	E = Point(diagbox[1] + factor, diagbox[0] + factor);
	W = Point(diagbox[7], diagbox[6]);
	S = Point(diagbox[7] + factor, diagbox[6] + factor);
	ends[2][0] = LocalCentroidDiag(img, label, lbl, N, S, W, E);

	N = Point(diagbox[5] - factor, diagbox[4] - factor);
	E = Point(diagbox[5], diagbox[4]);
	W = Point(diagbox[3] - factor, diagbox[2] - factor);
	S = Point(diagbox[3], diagbox[2]);
	ends[2][1] = LocalCentroidDiag(img, label, lbl, N, S, W, E);

	//Draw off-diag line...
	N = Point(diagbox[7] + factor, diagbox[6] - factor);
	E = Point(diagbox[3] + factor, diagbox[2] - factor);
	W = Point(diagbox[7], diagbox[6]);
	S = Point(diagbox[3], diagbox[2]);
	ends[3][0] = LocalCentroidDiag(img, label, lbl, N, S, W, E);

	N = Point(diagbox[1], diagbox[0]);
	E = Point(diagbox[5], diagbox[4]);
	W = Point(diagbox[1] - factor, diagbox[0] + factor);
	S = Point(diagbox[5] - factor, diagbox[4] + factor);
	ends[3][1] = LocalCentroidDiag(img, label, lbl, N, S, W, E);
}
