#ifndef _OIC_H_
#define _OIC_H_

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
using namespace cv;
namespace cv {
	using std::vector;
};

extern void MakeOIC(int **img, vector<Point> &sequence, int g);

int ObjectOnBoundary(int **img, int ilimit[2], int jlimit[2]);

int vtype(int **img, Point q, int g);

Point start(int **img, Point q, int g, int &ts);

extern void MakeMultiOIC(Mat colimg, int **img, int R, int C, vector<Point> &sequence, int g);

#endif
