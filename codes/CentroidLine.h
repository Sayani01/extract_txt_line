#ifndef _CENTROIDLINE_H_
#define _CENTROIDLINE_H_

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
using namespace cv;
namespace cv {
	using std::vector;
};

extern Point LocalCentroidDiag(int **img, int **label, int lbl, Point N, Point S, Point W, Point E);

extern Point LocalCentroidUpright(int **img, int **label, int lbl, Point TL, Point TR, Point DR, Point DL);

extern void Centroidline(int **img, int **label, int lbl, double *diagbox, int *box, Point ends[4][2]);

#endif