#ifndef _KWK_H_
#define _KWK_H_

extern int KWKthinning(int **src, int r, int c, int **dest, int iter);
extern void KWKthinCC(int **roi, int row, int col, int **thin_roi);

#endif