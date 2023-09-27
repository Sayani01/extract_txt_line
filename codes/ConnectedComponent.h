#ifndef _CONNECTEDCOMPONENT_H_
#define _CONNECTEDCOMPONENT_H_

extern int labelling(int **img, int **label, int rows, int cols);

extern void boundingRectangles(int **box, int **label, int row, int col, int lbl);
#endif
