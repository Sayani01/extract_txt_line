#ifndef _GHPTA_H_
#define _GHPTA_H_

extern int GHthinning(int **src, int r, int c, int **dest, int iter);
extern void GHthinCC(int **roi, int row, int col, int **thin_roi);

#endif

