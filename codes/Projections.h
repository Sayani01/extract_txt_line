#ifndef _PROJECTIONS_H_
#define _PROJECTIONS_H_

extern double stddev(int *x, int sz);
extern double stddev(double *x, int sz);
extern double project0(int **roi, int h, int w, int *hist);
extern double project90(int **roi, int h, int w, int *hist);
extern double project45(int **img, int **label, double *dim, int lbl, int *hist);
extern double project135(int **img, int **label, double *dim, int lbl, int *hist);
extern void projectAllCross(int **img, int h, int w, int &inclineCrossing);
extern void projectAllStdDev(int **img, int h, int w, int &inclineStdDev);
extern void projectAll(int **img, int h, int w, int &inclineStdDev, int &inclineExtent, int &inclineCrossing);
extern void trimCC(int **img, int **label, int **box, double **boundingBox);
#endif
