#ifndef _FPTA_H_
#define _FPTA_H_

extern int thinning(int **src, int r, int c, int **dest, int iter);
extern void thinCC(int **roi, int row, int col, int **thin_roi);

#endif
