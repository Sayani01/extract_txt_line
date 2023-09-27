#ifndef _REGLINE_H_
#define _REGLINE_H_

extern double Mean(int *X, int N);
extern void regressionLine(int *X, int *Y, int N, double& slopeYX, double& slopeXY, double& interceptyx, double& interceptxy);
#endif
