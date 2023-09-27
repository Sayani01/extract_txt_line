#ifndef _NEIGHBOURCC_H_
#define _NEIGHBOURCC_H_

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
using namespace cv;
namespace cv {
	using std::vector;
};

extern void nbr8Label(int **img, int **label, int **box, double **boundingBox, double **Data, int **change, Point **ends, int iter);

extern void allNbrLabels(int **img, int **label, int **box, double **boundingBox, double **Data, int **change, int iter);

extern void nearestLabel(int **img, int **label, int **box, double **Data, int **change);

extern void connectNbrs(int **img, int **label, int **box, double **boundingBox, double **Data, Point **ends, Mat pointers);

extern void near2axis(int **img, int **label, int **box, double **boundingBox, double **Data, Point **ends, Mat pointers);

extern void HausdorffNbrs(int **img, int **label, int **box, double **boundingBox, double **Data, Point **ends, Mat pointers);

extern void connectHNbrs(int **img, int **label, int **box, double **boundingBox, double **Data, Point **ends, double **connect, Mat pointers);

extern void connectVNbrs(int **img, int **label, int **box, double **boundingBox, double **Data, Point **ends, double **connect, Mat pointers);

extern void connecDNbrs(int **img, int **label, int **box, double **boundingBox, double **Data, Point **ends, double **connect, Mat pointers);
#endif
