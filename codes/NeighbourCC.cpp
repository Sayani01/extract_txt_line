#include <stdio.h>
#include<conio.h>
#include<iostream>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "VarDeclare.h"
#include "NeighbourCC.h"

using namespace std;
using namespace cv;
namespace cv {
	using std::vector;
};


double avgDensity(int **img, int **label, int **box, double **Data) {
	int *CCpixels = new int[labelNum];
	//Count the pixels of each label...
	for (int i = 0; i < labelNum; i++)
		CCpixels[i] = 0;
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			int lbl = label[i][j];
			CCpixels[lbl]++;
		}
	}
	double avgdensity = 0, n = 0;

	for (int i = 0; i < labelNum; i++) {
		if (box[i][2] == -1 || int(Data[i][3]) == 0) continue;
		int clr[9] = { 0 }; //count the #pixels of each neighboring color...
		double density = 0;
		int i_start = box[i][4], j_start = box[i][2];

		int hmin = i_start - (int)(2.5 * AH);
		if (hmin < 0) hmin = 0;
		int hmax = box[i][5] + (int)(2.5 * AH);
		if (hmax >= r) hmax = r - 1;

		int wmin = j_start - (int)(2.5 * AH);
		if (wmin < 0) wmin = 0;
		int wmax = box[i][3] + (int)(2.5 * AH);
		if (wmax >= c) wmax = c - 1;

		for (int k = hmin; k <= hmax; k++) {
			for (int l = wmin; l <= wmax; l++) {
				if (label[k][l] != i && img[k][l] == 0)
					density++;
			}
		}
		density /= ((hmax - hmin + 1)*(wmax - wmin + 1) - (box[i][5] - box[i][4] + 1)*(box[i][3] - box[i][2] + 1));
		avgdensity += density;
		n++;
	}
	avgdensity /= n;
	return(avgdensity);
}

int isMoreCut(int **img, int**label, int **box, double **Data, int oldlabel, int newlabel) {
	int countoldcut = 0, countnewcut = 0;
	int prev = 255, curr = 255;
	
	int **oldroi, **newroi;
	int h_old = box[oldlabel][5] - box[oldlabel][4] + 1;
	int w_old = box[oldlabel][3] - box[oldlabel][2] + 1;
	int h_new = box[newlabel][5] - box[newlabel][4] + 1;
	int w_new = box[newlabel][3] - box[newlabel][2] + 1;
	oldroi = new int*[h_old];
	for (int i = 0; i < h_old; i++) {
		oldroi[i] = new int[w_old];
		for (int j = 0; j < w_old; j++)
			if (label[box[oldlabel][4] + i][box[oldlabel][2] + j] == oldlabel)
				oldroi[i][j] = 0;
			else
				oldroi[i][j] = 255;
	}

	std::vector<Point> trackold, tracknew;
	double oldslope, newslope;
	int centroidx = Data[oldlabel][0] - box[oldlabel][2], centroidy = Data[oldlabel][1] - box[oldlabel][4];
	std::cout << "\n*****find tracks";
	//Track line with old slope...
	if (Data[oldlabel][2] == 90) {
		for(int i=0; i < h_old; i++)
			trackold.push_back(Point(centroidx, i));
	}
	else if (Data[oldlabel][2] == 0) {
		for (int j = 0; j < w_old; j++)
			trackold.push_back(Point(j, centroidy));
	}
	else {
		oldslope = tan(Data[oldlabel][2] * CV_PI / 180);
		int startx = 0, stopx = w_old - 1;
		int starty = int(centroidy + oldslope * (startx - centroidx));
		if (starty >= h_old) {
			starty = h_old - 1;
			startx = int(centroidx + (starty - centroidy) / oldslope);
		}
		if (starty < 0) {
			starty = 0;
			startx = int(centroidx - (centroidy / oldslope));
		}
		int stopy = int(centroidy + oldslope * (stopx - centroidx));
		if (stopy >= h_old) {
			stopy = h_old - 1;
			stopx = int(centroidx + (stopy - centroidy) / oldslope);
		}
		if (stopy < 0) {
			stopy = 0;
			stopx = int(centroidx - (centroidy / oldslope));
		}
		int i = starty, j = startx;
		do {
			trackold.push_back(Point(j, i));
			if (startx < stopx) j++;
			else j--;
			i = int(centroidy + oldslope * (j - centroidx));
		} while (j != stopx);
		trackold.push_back(Point(stopx, stopy));
	}
	std::cout << "\n*****got old track";
	//Track line with new slope...
	if (Data[newlabel][2] == 90) {
		for (int i = 0; i < h_old; i++)
			tracknew.push_back(Point(centroidx, i));
	}
	else if (Data[newlabel][2] == 0) {
		for (int j = 0; j < w_old; j++)
			tracknew.push_back(Point(j, centroidy));
	}
	else {
		newslope = tan(Data[newlabel][2] * CV_PI / 180);
		int startx = 0, stopx = w_old - 1;
		int starty = int(centroidy + newslope * (startx - centroidx));
		if (starty >= h_old) {
			starty = h_old - 1;
			startx = int(centroidx + (starty - centroidy) / newslope);
		}
		if (starty < 0) {
			starty = 0;
			startx = int(centroidx - (centroidy / newslope));
		}
		int stopy = int(centroidy + newslope * (stopx - centroidx));
		if (stopy >= h_old) {
			stopy = h_old - 1;
			stopx = int(centroidx + (stopy - centroidy) / newslope);
		}
		if (stopy < 0) {
			stopy = 0;
			stopx = int(centroidx - (centroidy / newslope));
		}
		int i = starty, j = startx;
		do {
			tracknew.push_back(Point(j, i));
			if (startx < stopx) j++;
			else j--;
			i = int(centroidy + newslope * (j - centroidx));
		} while (j != stopx);
		tracknew.push_back(Point(stopx, stopy));
	}
	std::cout << "\n*****got new track";
	//Count 1 -> 0 transitions in old line...
	prev = oldroi[trackold[0].y][trackold[0].x];
	for (int k = 1; k < int(trackold.size()); k++){
		curr = oldroi[trackold[k].y][trackold[k].x];
		if (prev != 0 && curr == 0) { //Bingo!...
			countoldcut++;
		}
		prev = curr;
	}
	std::cout << "\nNo. of cuts with old slope = " << countoldcut;

	//Count 1 -> 0 transitions in new line...
	prev = oldroi[tracknew[0].y][tracknew[0].x];
	for (int k = 1; k < int(tracknew.size()); k++) {
		curr = oldroi[tracknew[k].y][tracknew[k].x];
		if (prev != 0 && curr == 0) { //Bingo!...
			countnewcut++;
		}
		prev = curr;
	}
	std::cout << "\nNo. of cuts with new slope = " << countnewcut;

	//Delete memory...
	for (int i = 0; i < h_old; i++)
		delete[] oldroi[i];
	delete[] oldroi;

	if (countnewcut < countoldcut) return(0);
	else return(1);
}

double span(int **roi, int h, int w, double *data, int *box, Point ends[2], double angle) {
	Point first, last;
	int xmin = ends[0].x - box[2];
	int ymin = ends[0].y - box[4];
	int xmax = ends[1].x - box[2];
	int ymax = ends[1].y - box[4];
	int x0 = data[0] - box[2];
	int y0 = data[1] - box[4];
	int x, y;
	if (angle == 0 || angle == 180) {
		y = y0;
		//search first black point on horizontal line...
		for (int i = 0; i < w; i++) {
			if (roi[y][i] == 0) {
				first.x = i;
				first.y = y;
				break;
			}
		}
		//search last black point on horizontal line...
		for (int i = w - 1; i >= 0; i--) {
			if (roi[y][i] == 0) {
				first.x = i;
				first.y = y;
				break;
			}
		}
	}
	else if (angle == 90) {
		x = x0;
		for (int i = 0; i < h; i++) {
			if (roi[i][x] == 0) {
				first.x = x;
				first.y = i;
			}
		}
		for (int i = h - 1; i >= 0; i--) {
			if (roi[i][x] == 0) {
				first.x = x;
				first.y = i;
			}
		}
	}
	else {
		double slope = tan(angle * CV_PI / 180);
		//search first black point on line...
		for (int i = xmin + 1; i < xmax; i++) {
			y = y0 + slope*(i - x0);
			if (roi[y][i] == 0) {
				first.x = i;
				first.y = y;
				break;
			}
		}
		//std::cout << "\nhere " << first.x << " , " << first.y;
		//Search last black point on line...
		for (int i = xmax - 1; i > xmin; i--) {
			y = y0 + slope*(i - x0);
			if (roi[y][i] == 0) {
				last.x = i;
				last.y = y;
				break;
			}
		}
		//std::cout << "\nhere " << last.x << " , " << last.y;
	}
	
	return(norm(last - first));
}

void nbr8Label(int **img, int **label, int **box, double **boundingBox, double **Data, int **change, Point **ends, int iter) {
	int *CCpixels = new int[labelNum];
	//Count the pixels of each label...
	for (int i = 0; i < labelNum; i++)
		CCpixels[i] = 0;
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			int lbl = label[i][j];
			CCpixels[lbl]++;
		}
	}

	int *clr = new int[9]; //count the #pixels of each neighboring color...

						   //--------------------------Re-cluster based on nearest available neighbor-----------------------
	for (int i = 0; i < labelNum; i++) {
		if (box[i][2] == -1 || int(Data[i][5]) == 0) continue;

		int N_label = 0, E_label = 0, S_label = 0, W_label = 0, NE_label = 0, SE_label = 0, SW_label = 0, NW_label = 0, maxclr = 0, clr_count = 0;
		double N_dist = 999, E_dist = 999, S_dist = 999, W_dist = 999, NE_dist = 999, SE_dist = 999, SW_dist = 999, NW_dist = 999, mindist = 999;
		int newcolor = 0, distLabel = 0;
		int nbrlabels[8] = { 0 };

		//North nbr...
		int k = box[i][4] - 1;
		while (k >= 0 && box[i][4] - k < 4 * AH) {
			for (int l = box[i][2]; l <= box[i][3]; l++) {
				if (img[k][l] == 0 && Data[label[k][l]][5] != 0) {
					N_label = label[k][l];
					N_dist = box[i][4] - k;
					break;
				}
			}
			if (N_label != 0)
				break;
			k--;
		}
		if (N_dist < mindist && N_label != 0) {
			mindist = N_dist;
			distLabel = N_label;
		}
		//-------------------------------------------------------------------------------------------------------------------

		//South nbr...
		k = box[i][5] + 1;
		while (k < r && k - box[i][5] < 4 * AH) {
			for (int l = box[i][2]; l <= box[i][3]; l++) {
				if (img[k][l] == 0 && Data[label[k][l]][5] != 0) {
					S_label = label[k][l];
					S_dist = k - box[i][5];
					break;
				}
			}
			if (S_label != 0)
				break;
			k++;
		}
		if (S_dist < mindist && S_label != 0) {
			mindist = S_dist;
			distLabel = S_label;
		}
		//-------------------------------------------------------------------------------------------------------------------

		//West nbr...
		int l = box[i][2] - 1;
		while (l >= 0 && box[i][2] - l < 4 * AH) {
			for (int k = box[i][4]; k <= box[i][5]; k++) {
				if (img[k][l] == 0 && Data[label[k][l]][5] != 0) {
					W_label = label[k][l];
					W_dist = box[i][2] - l;
					break;
				}
			}
			if (W_label != 0)
				break;
			l--;
		}
		if (W_dist < mindist && W_label != 0) {
			mindist = W_dist;
			distLabel = W_label;
		}
		//-------------------------------------------------------------------------------------------------------------------

		//East nbr...
		l = box[i][3] + 1;
		while (l < c && l - box[i][3] < 4 * AH) {
			for (int k = box[i][4]; k <= box[i][5]; k++) {
				if (img[k][l] == 0 && Data[label[k][l]][5] != 0) {
					E_label = label[k][l];
					E_dist = l - box[i][3];
					break;
				}
			}
			if (E_label != 0)
				break;
			l++;
		}
		if (E_dist < mindist && E_label != 0) {
			mindist = E_dist;
			distLabel = E_label;
		}
		//-------------------------------------------------------------------------------------------------------------------

		//NE nbr...
		int hmin, hmax, wmin, wmax, h, w;
		double dist = 0;
		hmin = box[i][4] - 4 * AH;
		if (hmin < 0) hmin = 0;
		hmax = box[i][4] - 1;
		if (hmax < 0) hmax = 0;
		wmin = box[i][3] + 1;
		if (wmin >= c) wmin = c - 1;
		wmax = box[i][3] + 4 * AH;
		if (wmax >= c) wmax = c - 1;
		h = hmax - hmin + 1, w = wmax - wmin + 1;
		int **roi = new int*[h], **lblroi = new int*[h];
		for (int k = 0; k < h; k++) {
			roi[k] = new int[w];
			lblroi[k] = new int[w];
			for (int l = 0; l < w; l++) {
				roi[k][l] = img[hmin + k][wmin + l];
				lblroi[k][l] = label[hmin + k][wmin + l];
			}
		}
		//std::cout <<"****" << h << " , " << w <<"    ";
		int flag = 0;
		if (h <= w) {
			int b;
			for (b = 1 - h; b <= 0; b++) {
				int prod = 1;
				for (int y = -b; y < h; y++) {
					prod *= (roi[y][y + b] / 255);
					if (roi[y][y + b] == 0 && Data[lblroi[y][y + b]][5] != 0) {
						//std::cout << "yes!" << y + b << "," << y << "**" <<roi[y][y + b];
						dist = sqrt((hmax - 1 - y)*(hmax - 1 - y) + (y + b)*(y + b));
						if (NE_dist > dist) {
							NE_label = lblroi[y][y + b];
							NE_dist = dist;
						}
					}
				}
				if (prod == 0) {
					//std::cout << "***" << NE_label;
					flag = 1;
					break;
				}
			}
			if (flag != 1) {
				for (b = 1; b <= w - h; b++) {
					int prod = 1;
					for (int y = 0; y < h; y++) {
						prod *= (roi[y][y + b] / 255);
						if (roi[y][y + b] == 0 && Data[lblroi[y][y + b]][5] != 0) {
							//std::cout << "yes!" << y + b << "," << y << "**" << roi[y][y + b];
							dist = sqrt((hmax - 1 - y)*(hmax - 1 - y) + (y + b)*(y + b));
							if (NE_dist > dist) {
								NE_label = lblroi[y][y + b];
								NE_dist = dist;
							}
						}
					}
					if (prod == 0) {
						//std::cout << "\nHit black!";
						flag = 1;
						break;
					}
				}
			}
			if (flag != 1) {
				for (b = w - h + 1; b < w; b++) {
					int prod = 1;
					for (int y = 0; y < w - b; y++) {
						prod *= (roi[y][y + b] / 255);
						if (roi[y][y + b] == 0 && Data[lblroi[y][y + b]][5] != 0) {
							//std::cout << "yes!" << y + b << "," << y << "**" << roi[y][y + b];
							dist = sqrt((hmax - 1 - y)*(hmax - 1 - y) + (y + b)*(y + b));
							if (NE_dist > dist) {
								NE_label = lblroi[y][y + b];
								NE_dist = dist;
							}
						}
					}
					if (prod == 0) {
						//std::cout << "\nHit black!";
						flag = 1;
						break;
					}
				}
			}
			//std::cout << "\n h <= w : b = " << b;
		}
		else {
			int b;
			for (b = 1 - h; b <= w - h; b++) {
				int prod = 1;
				for (int x = 0; x < h + b; x++) {
					prod *= (roi[x - b][x] / 255);
					if (roi[x - b][x] == 0 && Data[lblroi[x - b][x]][5] != 0) {
						//std::cout << "yes!" << x << "," << x - c << "**" << roi[x - b][x];
						dist = sqrt((hmax - 1 - x + b)*(hmax - 1 - x + b) + x*x);
						if (NE_dist > dist) {
							NE_label = lblroi[x - b][x];
							NE_dist = dist;
						}
					}
				}
				if (prod == 0) {
					//std::cout << "\nHit black!";
					flag = 1;
					break;
				}
			}
			if (flag != 1) {
				for (b = w - h + 1; b <= 0; b++) {
					int prod = 1;
					for (int x = 0; x < w; x++) {
						prod *= (roi[x - b][x] / 255);
						if (roi[x - b][x] == 0 && Data[lblroi[x - b][x]][5] != 0) {
							//std::cout << "yes!" << x << "," << x - b << "**" << roi[x - b][x];
							dist = sqrt((hmax - 1 - x + b)*(hmax - 1 - x + b) + x*x);
							if (NE_dist > dist) {
								NE_label = lblroi[x - b][x];
								NE_dist = dist;
							}
						}
					}
					if (prod == 0) {
						//std::cout << "\nHit black!";
						flag = 1;
						break;
					}
				}
			}
			if (flag != 1) {
				for (b = 1; b < w; b++) {
					int prod = 1;
					for (int x = b; x < w; x++) {
						prod *= (roi[x - b][x] / 255);
						if (roi[x - b][x] == 0 && Data[lblroi[x - b][x]][5] != 0) {
							//std::cout << "yes!" << x << "," << x - c << "**" << roi[x - b][x];
							dist = sqrt((hmax - 1 - x + b)*(hmax - 1 - x + b) + x*x);
							if (NE_dist > dist) {
								NE_label = lblroi[x - b][x];
								NE_dist = dist;
							}
						}
					}
					if (prod == 0) {
						//std::cout << "\nHit black!";
						flag = 1;
						break;
					}
				}
			}
			//std::cout << "\n h > w : b = " << b;
		}
		if (NE_label == distLabel) {
			mindist = mindist > NE_dist ? NE_dist : mindist;
		}
		if (NE_dist < mindist && NE_dist > 0) {
			mindist = NE_dist;
			distLabel = NE_label;
		}
		//Delete...
		for (int k = 0; k < h; k++) {
			delete[] roi[k];
			delete[] lblroi[k];
		}
		delete[] roi;
		delete[] lblroi;
		//-------------------------------------------------------------------------------------------------------------------

		//NW nbr...
		dist = 0;
		hmin = box[i][4] - 4 * AH;
		if (hmin < 0) hmin = 0;
		hmax = box[i][4] - 1;
		if (hmax < 0) hmax = 0;
		wmin = box[i][2] - 4 * AH;
		if (wmin < 0) wmin = 0;
		wmax = box[i][2] - 1;
		if (wmax < 0) wmax = 0;
		h = hmax - hmin + 1, w = wmax - wmin + 1;
		int **roinw = new int*[h], **lblroinw = new int*[h];
		for (int k = 0; k < h; k++) {
			roinw[k] = new int[w];
			lblroinw[k] = new int[w];
			for (int l = 0; l < w; l++) {
				roinw[k][l] = img[hmin + k][wmin + l];
				lblroinw[k][l] = label[hmin + k][wmin + l];
			}
		}
		//std::cout << "****" << h << " , " << w << "    ";
		flag = 0;
		if (h <= w) {
			int b;
			for (b = h + w - 2; b >= w - 1; b--) {
				int prod = 1;
				for (int y = b - w + 1; y < h; y++) {
					prod *= (roinw[y][b - y] / 255);
					if (roinw[y][b - y] == 0 && Data[lblroinw[y][b - y]][5] != 0) {
						//std::cout << "yes!" << b - y << "," << y << "**" << roinw[y][b - y];
						dist = sqrt((hmax - 1 - y)*(hmax - 1 - y) + (wmax - 1 - b + y)*(wmax - 1 - b + y));
						if (NW_dist > dist) {
							NW_label = lblroinw[y][b - y];
							NW_dist = dist;
						}
					}
				}
				if (prod == 0) {
					//std::cout << "***" << NW_label;
					flag = 1;
					break;
				}
			}
			if (flag != 1) {
				for (b = w - 2; b >= h - 1; b--) {
					int prod = 1;
					for (int y = 0; y < h; y++) {
						prod *= (roinw[y][b - y] / 255);
						if (roinw[y][b - y] == 0 && Data[lblroinw[y][b - y]][5] != 0) {
							//std::cout << "yes!" << b - y << "," << y << "**" << roinw[y][b - y];
							dist = sqrt((hmax - 1 - y)*(hmax - 1 - y) + (wmax - 1 - b + y)*(wmax - 1 - b + y));
							if (NW_dist > dist) {
								NW_label = lblroinw[y][b - y];
								NW_dist = dist;
							}
						}
					}
					if (prod == 0) {
						//std::cout << "\nHit black!";
						flag = 1;
						break;
					}
				}
			}
			if (flag != 1) {
				for (b = h - 2; b >= 0; b--) {
					int prod = 1;
					for (int y = 0; y <= b; y++) {
						prod *= (roinw[y][b - y] / 255);
						if (roinw[y][b - y] == 0 && Data[lblroinw[y][b - y]][5] != 0) {
							//std::cout << "yes!" << b - y << "," << y << "**" << roinw[y][b - y];
							dist = sqrt((hmax - 1 - y)*(hmax - 1 - y) + (wmax - 1 - b + y)*(wmax - 1 - b + y));
							if (NW_dist > dist) {
								NW_label = lblroinw[y][b - y];
								NW_dist = dist;
							}
						}
					}
					if (prod == 0) {
						//std::cout << "\nHit black!";
						flag = 1;
						break;
					}
				}
			}
			//std::cout << "\n h <= w : b = " << b;
		}
		else {
			int b;
			for (b = h + w - 2; b >= h - 1; b--) {
				int prod = 1;
				for (int y = b - w + 1; y < h; y++) {
					prod *= (roinw[y][b - y] / 255);
					if (roinw[y][b - y] == 0 && Data[lblroinw[y][b - y]][5] != 0) {
						//std::cout << "yes!" << b - y << "," << y << "**" << roinw[y][b - y];
						dist = sqrt((hmax - 1 - y)*(hmax - 1 - y) + (wmax - 1 - b + y)*(wmax - 1 - b + y));
						if (NW_dist > dist) {
							NW_label = lblroinw[y][b - y];
							NW_dist = dist;
						}
					}
				}
				if (prod == 0) {
					//std::cout << "***" << NW_label;
					flag = 1;
					break;
				}
			}
			if (flag != 1) {
				for (b = h - 2; b >= w - 1; b--) {
					int prod = 1;
					for (int x = 0; x < w; x++) {
						prod *= (roinw[b - x][x] / 255);
						if (roinw[b - x][x] == 0 && Data[lblroinw[b - x][x]][5] != 0) {
							//std::cout << "yes!" << x << "," << b - x << "**" << roinw[b - x][x];
							dist = sqrt((hmax - 1 - b + x)*(hmax - 1 - b + x) + (wmax - 1 - x)*(wmax - 1 - x));
							if (NW_dist > dist) {
								NW_label = lblroinw[b - x][x];
								NW_dist = dist;
							}
						}
					}
					if (prod == 0) {
						//std::cout << "\nHit black!";
						flag = 1;
						break;
					}
				}
			}
			if (flag != 1) {
				for (b = w - 2; b >= 0; b--) {
					int prod = 1;
					for (int y = 0; y <= b; y++) {
						prod *= (roinw[y][b - y] / 255);
						if (roinw[y][b - y] == 0 && Data[lblroinw[y][b - y]][5] != 0) {
							//std::cout << "yes!" << b - y << "," << y << "**" << roinw[y][b - y];
							dist = sqrt((hmax - 1 - y)*(hmax - 1 - y) + (wmax - 1 - b + y)*(wmax - 1 - b + y));
							if (NW_dist > dist) {
								NW_label = lblroinw[y][b - y];
								NW_dist = dist;
							}
						}
					}
					if (prod == 0) {
						//std::cout << "\nHit black!";
						flag = 1;
						break;
					}
				}
			}
			//std::cout << "\n h > w : b = " << b;
		}
		if (NW_label == distLabel) {
			mindist = mindist > NW_dist ? NW_dist : mindist;
		}
		if (NW_dist < mindist && NW_dist > 0) {
			mindist = NW_dist;
			distLabel = NW_label;
		}

		//Delete...
		for (int k = 0; k < h; k++) {
			delete[] roinw[k];
			delete[] lblroinw[k];
		}
		delete[] roinw;
		delete[] lblroinw;
		//-------------------------------------------------------------------------------------------------------------------

		//SE nbr...
		dist = 0;
		hmin = box[i][5] + 1;
		if (hmin >= r) hmin = r - 1;
		hmax = box[i][5] + 4 * AH;
		if (hmax >= r) hmax = r - 1;
		wmin = box[i][3] + 1;
		if (wmin >= c) wmin = c - 1;
		wmax = box[i][3] + 4 * AH;
		if (wmax >= c) wmax = c - 1;
		h = hmax - hmin + 1, w = wmax - wmin + 1;
		int **roise = new int*[h], **lblroise = new int*[h];
		for (int k = 0; k < h; k++) {
			roise[k] = new int[w];
			lblroise[k] = new int[w];
			for (int l = 0; l < w; l++) {
				roise[k][l] = img[hmin + k][wmin + l];
				lblroise[k][l] = label[hmin + k][wmin + l];
			}
		}
		//std::cout <<"****" << h << " , " << w <<"    ";
		flag = 0;
		if (h <= w) {
			int b;
			for (b = 0; b <= h - 1; b++) {
				int prod = 1;
				for (int y = 0; y <= b; y++) {
					prod *= (roise[y][b - y] / 255);
					if (roise[y][b - y] == 0 && Data[lblroise[y][b - y]][5] != 0) {
						//std::cout << "yes!" << y + b << "," << y << "**" <<roi[y][y + b];
						dist = sqrt(y*y + (b - y)*(b - y));
						if (SE_dist > dist) {
							SE_label = lblroise[y][b - y];
							SE_dist = dist;
						}
					}
				}
				if (prod == 0) {
					//std::cout << "***" << SE_label;
					flag = 1;
					break;
				}
			}
			if (flag != 1) {
				for (b = h; b <= w - 1; b++) {
					int prod = 1;
					for (int y = 0; y < h; y++) {
						prod *= (roise[y][b - y] / 255);
						if (roise[y][b - y] == 0 && Data[lblroise[y][b - y]][5] != 0) {
							//std::cout << "yes!" << y + b << "," << y << "**" << roi[y][y + b];
							dist = sqrt(y*y + (b - y)*(b - y));
							if (SE_dist > dist) {
								SE_label = lblroise[y][b - y];
								SE_dist = dist;
							}
						}
					}
					if (prod == 0) {
						//std::cout << "\nHit black!";
						flag = 1;
						break;
					}
				}
			}
			if (flag != 1) {
				for (b = w ; b <= w + h - 2; b++) {
					int prod = 1;
					for (int y = b - w + 1; y < h; y++) {
						prod *= (roise[y][b - y] / 255);
						if (roise[y][b - y] == 0 && Data[lblroise[y][b - y]][5] != 0) {
							//std::cout << "yes!" << y + b << "," << y << "**" << roi[y][y + b];
							dist = sqrt(y*y + (b - y)*(b - y));
							if (SE_dist > dist) {
								SE_label = lblroise[y][b - y];
								SE_dist = dist;
							}
						}
					}
					if (prod == 0) {
						//std::cout << "\nHit black!";
						flag = 1;
						break;
					}
				}
			}
			//std::cout << "\n h <= w : b = " << b;
		}
		else {
			int b;
			for (b = 0; b <= w - 1; b++) {
				int prod = 1;
				for (int x = 0; x <= b; x++) {
					prod *= (roise[b - x][x] / 255);
					if (roise[b - x][x] == 0 && Data[lblroise[b - x][x]][5] != 0) {
						//std::cout << "yes!" << x << "," << x - c << "**" << roi[x - b][x];
						dist = sqrt((b - x)*(b - x) + x*x);
						if (SE_dist > dist) {
							SE_label = lblroise[b - x][x];
							SE_dist = dist;
						}
					}
				}
				if (prod == 0) {
					//std::cout << "\nHit black!";
					flag = 1;
					break;
				}
			}
			if (flag != 1) {
				for (b = w; b <= h - 1; b++) {
					int prod = 1;
					for (int x = 0; x < w ; x++) {
						prod *= (roise[b - x][x] / 255);
						if (roise[b - x][x] == 0 && Data[lblroise[b - x][x]][5] != 0) {
							//std::cout << "yes!" << x << "," << x - c << "**" << roi[x - b][x];
							dist = sqrt((b - x)*(b - x) + x*x);
							if (SE_dist > dist) {
								SE_label = lblroise[b - x][x];
								SE_dist = dist;
							}
						}
					}
					if (prod == 0) {
						//std::cout << "\nHit black!";
						flag = 1;
						break;
					}
				}
			}
			if (flag != 1) {
				for (b = h; b <= w + h - 2; b++) {
					int prod = 1;
					for (int x = b - h + 1; x < w; x++) {
						prod *= (roise[b - x][x] / 255);
						if (roise[b - x][x] == 0 && Data[lblroise[b - x][x]][5] != 0) {
							//std::cout << "yes!" << x << "," << x - c << "**" << roi[x - b][x];
							dist = sqrt((b - x)*(b - x) + x*x);
							if (SE_dist > dist) {
								SE_label = lblroise[b - x][x];
								SE_dist = dist;
							}
						}
					}
					if (prod == 0) {
						//std::cout << "\nHit black!";
						flag = 1;
						break;
					}
				}
			}
			//std::cout << "\n h > w : b = " << b;
		}
		if (SE_label == distLabel) {
			mindist = mindist > SE_dist ? SE_dist : mindist;
		}
		if (SE_dist < mindist && SE_dist > 0) {
			mindist = SE_dist;
			distLabel = SE_label;
		}
		//Delete...
		for (int k = 0; k < h; k++) {
			delete[] roise[k];
			delete[] lblroise[k];
		}
		delete[] roise;
		delete[] lblroise;
		//-------------------------------------------------------------------------------------------------------------------

		int currclr = int(Data[i][3]), nbrclr = int(Data[distLabel][3]);
		nbrlabels[1] = N_label;
		nbrlabels[2] = S_label;
		nbrlabels[3] = E_label;
		nbrlabels[4] = W_label;
		nbrlabels[5] = NE_label;
		nbrlabels[6] = NW_label;
		nbrlabels[7] = SE_label;
		/*std::cout << "\nNbrs for label " << i << " : ";
		for (int k = 1; k <= 6; k++) {
		std::cout << nbrlabels[k] << "   ";
		}
		std::cout << "\nDistances: " << N_dist << "," << S_dist << "," << E_dist << "," << W_dist << "," << NE_dist << "," << NW_dist;
		std::cout << "\nNearest label " << distLabel << " color = " << nbrclr ;*/

		clr[int(Data[nbrlabels[1]][3])] += CCpixels[nbrlabels[1]];
		for (int k = 2; k <= 7; k++) {
			int flag = 0;
			for (int l = 1; l < k; l++) {
				if (nbrlabels[l] == nbrlabels[k]) {
					//std::cout << "---same labels----";
					flag = 1;
					break;
				}
			}
			if (flag == 0) clr[int(Data[nbrlabels[k]][3])] += CCpixels[nbrlabels[k]];
		}

		int max = clr[1], maxcolor = 1;
		for (int k = 2; k < 9; k++) {
			if (max < clr[k]) {
				max = clr[k];
				maxcolor = k;
			}
		}
		double ht, wd;
		if (boundingBox[i][0] == 1) { //upright
			ht = box[i][5] - box[i][4] + 1;
			wd = box[i][3] - box[i][2] + 1;
		}
		if (boundingBox[i][0] == 2) { //diagonal
			ht = sqrt((boundingBox[i][5] - boundingBox[i][1])*(boundingBox[i][5] - boundingBox[i][1]) + (boundingBox[i][6] - boundingBox[i][2])*(boundingBox[i][6] - boundingBox[i][2]));
			wd = sqrt((boundingBox[i][3] - boundingBox[i][5])*(boundingBox[i][3] - boundingBox[i][5]) + (boundingBox[i][4] - boundingBox[i][6])*(boundingBox[i][4] - boundingBox[i][6]));
		}
		//------------------Re-color CC based on nearest neighbor & its gradient------------
		if (nbrclr != 0 && nbrclr != currclr) {
			if (iter <= 4 && mindist < 1.2*AH) {
				//must change case...
				if ((boundingBox[i][0] == 2 && (Data[i][3] == 1 || Data[i][3] == 3))
					|| (boundingBox[i][0] == 1 && ht < 0.5*wd && (Data[i][3] == 3 || Data[i][3] == 6 || Data[i][3] == 7))
					|| (boundingBox[i][0] == 1 && ht*0.5 > wd && (Data[i][3] == 1 || Data[i][3] == 5 || Data[i][3] == 8))) {  //diagonal dominant...
					change[i][0] = 1;
					change[i][1] = nbrclr;
					//change[i][2] = Data[distLabel][2];
					Data[i][3] = nbrclr;
					Data[i][2] = Data[distLabel][2];
				}
				
				else if (CCpixels[i] > CCpixels[distLabel]) {
					if (boundingBox[distLabel][0] == 2 && (Data[i][3] == 1 || Data[i][3] == 3 ) ) {  //diagonal dominant...
						continue;
					}
					if (boundingBox[distLabel][0] == 1 && ht < 0.5*wd && (Data[i][3] == 3 || Data[i][3] == 6 || Data[i][3] == 7)) {  //horizontal dominant...
						continue;
					}
					if (boundingBox[distLabel][0] == 1 && ht*0.5 > wd && (Data[i][3] == 1 || Data[i][3] == 5 || Data[i][3] == 8)) {  //vertical dominant...
						continue;
					}
					if (isMoreCut(img, label, box, Data, distLabel, i) ) {
						std::cout << "\nChanging near label with current slope!...";
						change[distLabel][0] = 1;
						change[distLabel][1] = currclr;
						//change[distLabel][2] = Data[i][2];
						Data[distLabel][3] = currclr;
						Data[distLabel][2] = Data[i][2];
						if (Data[distLabel][5] == 0) std::cout << "\nReplacing off CC!!";
					}
				}
				else {
					if (boundingBox[i][0] == 2 && (Data[distLabel][3] == 1 || Data[distLabel][3] == 3)) {  //diagonal dominant...
						continue;
					}
					if (boundingBox[i][0] == 1 && ht < 0.5*wd && (Data[distLabel][3] == 3 || Data[distLabel][3] == 6 || Data[distLabel][3] == 7)) {  //horizontal dominant...
						continue;
					}
					if (boundingBox[i][0] == 1 && ht > 0.5*wd && (Data[distLabel][3] == 1 || Data[distLabel][3] == 5 || Data[distLabel][3] == 8)) {  //vertical dominant...
						continue;
					}
					if (isMoreCut(img, label, box, Data, i, distLabel)) {
						std::cout << "\nChanging current label with near slope!...";
						change[i][0] = 1;
						change[i][1] = nbrclr;
						//change[i][2] = Data[distLabel][2];
						Data[i][3] = nbrclr;
						Data[i][2] = Data[distLabel][2];
						if (Data[distLabel][5] == 0) std::cout << "\nReplacing with off CC!!";
					}					
				}
			}
			else if (iter > 4) {
				double grad[8] = { 0 }, mingrad = 9999;
				int gradLabel = 0;
				for (int k = 1; k < 8; k++) {
					if (k % 2 == 0) { //S , W , NW...
						if (ends[i][0].x - ends[nbrlabels[k]][1].x != 0)
							grad[k] = (ends[i][0].y - ends[nbrlabels[k]][1].y) / (ends[i][0].x - ends[nbrlabels[k]][1].x);
						else
							grad[k] = 99999;
					}
					else { //N , E , NE, SE...
						if (ends[i][1].x - ends[nbrlabels[k]][0].x != 0)
							grad[k] = (ends[i][1].y - ends[nbrlabels[k]][0].y) / (ends[i][1].x - ends[nbrlabels[k]][0].x);
						else
							grad[k] = 99999;
					}
					if (grad[k] < mingrad) {
						mingrad = grad[k];
						gradLabel = nbrlabels[k];
					}
				}

				currclr = int(Data[i][3]), nbrclr = int(Data[gradLabel][3]);
				if (nbrclr != 0 && nbrclr != 0) {
					if (CCpixels[i] > CCpixels[gradLabel]) {
						change[gradLabel][0] = 1;
						change[gradLabel][1] = currclr;
						Data[gradLabel][3] = currclr;
						Data[gradLabel][2] = Data[i][2];
					}
					else {
						change[i][0] = 1;
						change[i][1] = nbrclr;
						Data[i][3] = nbrclr;
						Data[i][2] = Data[gradLabel][2];
					}

				}
			}
		}

	}
	delete[] clr;
}

void allNbrLabels(int **img, int **label, int **box, double **boundingBox, double **Data, int **change, int iter) {
	int *CCpixels = new int[labelNum];
	//Count the pixels of each label...
	for (int i = 0; i < labelNum; i++)
		CCpixels[i] = 0;
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			int lbl = label[i][j];
			CCpixels[lbl]++;
		}
	}

	int kk = 2.5, rho_k = 22;
	double avgdensity = avgDensity(img, label, box, Data);

	std::cout << "\nAllowed density = " << rho_k * avgdensity;

	for (int i = 0; i < labelNum; i++) {
		if (box[i][2] == -1 || int(Data[i][5]) == 0) continue;
		int clr[9] = { 0 }; //count the #pixels of each neighboring color...
		int totalTexts = 0, area;
		double density = 0, ht = 0;

		int i_start = box[i][4], j_start = box[i][2];
		int h = box[i][5] - box[i][4] + 1;
		int w = box[i][3] - box[i][2] + 1;

		//if (boundingBox[i][0] == 1) //upright
		//	ht = h;
		//if (boundingBox[i][0] == 2) //diagonal
		//	ht = sqrt((boundingBox[i][5] - boundingBox[i][1])*(boundingBox[i][5] - boundingBox[i][1]) + (boundingBox[i][6] - boundingBox[i][2])*(boundingBox[i][6] - boundingBox[i][2]));


		//if (ht > 2*AH || ht < 0.2 * AH) {
		//	std::cout << "\nLabel " << i << " excluded!";
		//	continue; //Exclude small pixels...
		//}

		//if (h >= 1.25*AH && iter > 12) {
		//	std::cout << "\nLabel " << i << " excluded!";
		//	continue; //Exclude larger pixels...
		//}

		int hmin = i_start - (int)(kk * AH);
		if (hmin < 0) hmin = 0;
		int hmax = box[i][5] + (int)(kk * AH);
		if (hmax >= r) hmax = r - 1;

		int wmin = j_start - (int)(kk * AH);
		if (wmin < 0) wmin = 0;
		int wmax = box[i][3] + (int)(kk * AH);
		if (wmax >= c) wmax = c - 1;

		area = (hmax - hmin + 1)*(wmax - wmin + 1) - h*w;
		std::vector<int> labels;

		//All nbrs on left...
		for (int k = hmin; k <= hmax; k++) {
			for (int l = wmin; l < j_start; l++) {
				int lbl = label[k][l], flag = 1;
				if (lbl == 0 || Data[lbl][5] == 0) continue;
				totalTexts++;
				for (int it = 0; it < int(labels.size()); it++) {
					if (lbl == labels[it]) {
						flag = 0;
						break;
					}
				}
				if (flag == 1) {
					labels.push_back(lbl);
				}
			}
		}

		//All nbrs on right...
		for (int k = hmin; k <= hmax; k++) {
			for (int l = box[i][3] + 1; l <= wmax; l++) {
				int lbl = label[k][l], flag = 1;
				if (lbl == 0 || Data[lbl][5] == 0) continue;
				totalTexts++;
				for (int it = 0; it < int(labels.size()); it++) {
					if (lbl == labels[it]) {
						flag = 0;
						break;
					}
				}
				if (flag == 1) {
					labels.push_back(lbl);
				}
			}
		}

		//All nbrs on top...
		for (int k = hmin; k <= i_start - 1; k++) {
			for (int l = box[i][2]; l < w; l++) {
				int lbl = label[k][l], flag = 1;
				if (lbl == 0 || Data[lbl][5] == 0) continue;
				totalTexts++;
				for (int it = 0; it < int(labels.size()); it++) {
					if (lbl == labels[it]) {
						flag = 0;
						break;
					}
				}
				if (flag == 1) {
					labels.push_back(lbl);
				}
			}
		}

		//All nbrs on bottom...
		for (int k = box[i][5] + 1; k <= hmax; k++) {
			for (int l = box[i][2]; l < w; l++) {
				int lbl = label[k][l], flag = 1;
				if (lbl == 0 || Data[lbl][5] == 0) continue;
				totalTexts++;
				for (int it = 0; it < int(labels.size()); it++) {
					if (lbl == labels[it]) {
						flag = 0;
						break;
					}
				}
				if (flag == 1) {
					labels.push_back(lbl);
				}
			}
		}

		/*std::cout << "\nLabel : " << i << " labels ";
		for (int it = 0; it<int(labels.size()); it++) {
		std::cout << labels[it] << "   ";
		}*/

		if (labels.size() == 0) {
			//std::cout << "No nbrs!";
			continue;
		}

		//Count the pixels of each color....
		int maxpixels = -1, maxclr = 0, max2pixels = -1, max2clr = 0;
		clr[int(Data[labels[0]][3])] += CCpixels[labels[0]];
		for (int k = 1; k < int(labels.size()); k++) {
			int lblclr = int(Data[labels[k]][3]);
			int flag = 0;
			for (int l = 0; l < k; l++) {
				if (labels[l] == labels[k]) {
					flag = 1;
					break;
				}
			}
			if(flag==0)
				clr[lblclr] += CCpixels[labels[k]];
		}

		for (int k = 1; k < 9; k++) {
			if (clr[k] == 0) continue;
			if (maxpixels < clr[k]) {
				maxpixels = clr[k];
				maxclr = k;
			}
		}
		for (int k = 1; k < 9; k++) {
			if (clr[k] == 0) continue;
			if (clr[k] == maxpixels) continue;
			if (max2pixels < clr[k]) {
				max2pixels = clr[k];
				max2clr = k;
			}
		}

		//Replace by max nbr color....
		//density = 100 * totalTexts / double(area);
		//std::cout << "\n Label " << i <<"  %density of texts = " << density;
		if (maxclr != int(Data[i][3])) {
			//std::cout << "  maxpixel = " << maxpixels << "  CCpixel = " << CCpixels[i];
			//if (density <= rho_k * avgdensity) continue;

			if (int(Data[i][3]) == 1) {  //0deg...
				if (maxclr == 5 && Data[i][2] == 5) {
					change[i][0] = 1;
					change[i][1] = maxclr;
					change[i][2] = 10;
				}
				if (maxclr == 8 && Data[i][2] == 175) {
					change[i][0] = 1;
					change[i][1] = maxclr;
					change[i][2] = 170;
				}
			}
			else if (int(Data[i][3]) == 2) {  //45deg...
				if (maxclr == 5 && Data[i][2] == 35) {
					change[i][0] = 1;
					change[i][1] = maxclr;
					change[i][2] = 30;
				}
				if (maxclr == 6 && Data[i][2] == 55) {
					change[i][0] = 1;
					change[i][1] = maxclr;
					change[i][2] = 60;
				}
			}
			else if (int(Data[i][3]) == 3) {  //90deg...
				if (maxclr == 6 && Data[i][2] == 85) {
					change[i][0] = 1;
					change[i][1] = maxclr;
					change[i][2] = 80;
				}
				if (maxclr == 7 && Data[i][2] == 95) {
					change[i][0] = 1;
					change[i][1] = maxclr;
					change[i][2] = 100;
				}
			}
			else if (int(Data[i][3]) == 4) {  //135deg...
				if (maxclr == 7 && Data[i][2] == 125) {
					change[i][0] = 1;
					change[i][1] = maxclr;
					change[i][2] = 120;
				}
				if (maxclr == 8 && Data[i][2] == 145) {
					change[i][0] = 1;
					change[i][1] = maxclr;
					change[i][2] = 150;
				}
			}
			/*else {
				change[i][0] = 1;
				change[i][1] = maxclr;
			}*/
		}
	}

	//--------------------Re-color CC's----------------------
	for (int i = 0; i < labelNum; i++) {
		if (Data[i][4] == 0) continue;
		//std::cout << "\nLabel " << i << " old = " << Data[i][3] << " new = " << change[i][1];
		if (change[i][0] == 1) {
			Data[i][3] = change[i][1];
			Data[i][2] = change[i][2];
		}
	}
}

void nearestLabel(int **img, int **label, int **box, double **Data, int **change) {
	int *CCpixels = new int[labelNum];
	for (int i = 0; i < labelNum; i++)
		CCpixels[i] = 0;
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			int lbl = label[i][j];
			CCpixels[lbl]++;
		}
	}

	int *clr = new int[9];

	//--------------------------Re-cluster based on nearest available neighbor----------------------- 
	for (int i = 0; i < labelNum; i++) {
		if (box[i][2] == -1 || Data[i][3] == 0) continue;

		double slope_i = Data[i][2];
		int color = Data[i][3];
		int top_label = 0, down_label = 0, right_label = 0, left_label = 0, maxPixel = 0, maxLabel = 0;
		double top_dist = 0, down_dist = 0, right_dist = 0, left_dist = 0, minDist = 0;
		int distLabel = 0;
		int i_start = box[i][4], j_start = box[i][2];
		int h = box[i][5] - box[i][4] + 1;
		int w = box[i][3] - box[i][2] + 1;
		for (int k = 0; k < 8; k++)
			clr[k] = 0;

		//Top nbr...
		int k = box[i][4] - 1;
		while (k >= 0 && box[i][4] - k < 4 * AH) {
			for (int l = box[i][2]; l <= box[i][3]; l++) {
				if (img[k][l] == 0) {
					top_label = label[k][l];
					top_dist = box[i][4] - k;
					break;
				}
			}
			if (top_label != 0) {
				//std::cout << "\nTop : " << top_label << " " << CCpixels[top_label];
				clr[int(Data[top_label][3])] += CCpixels[top_label];
				break;
			}
			k--;
		}
		minDist = top_dist;
		distLabel = top_label;

		//Down nbr...
		k = box[i][5] + 1;
		while (k < r && k - box[i][5] < 4 * AH) {
			for (int l = box[i][2]; l <= box[i][3]; l++) {
				if (img[k][l] == 0) {
					down_label = label[k][l];
					down_dist = k - box[i][5];
					break;
				}
			}
			if (down_label != 0) {
				//std::cout << "\nDown : " << down_label << " " << CCpixels[down_label];
				clr[int(Data[down_label][3])] += CCpixels[down_label];
				break;
			}
			k++;
		}
		if (down_dist < minDist && down_label != 0) {
			minDist = down_dist;
			distLabel = down_label;
		}

		//Left nbr...
		int l = box[i][2] - 1;
		while (l >= 0 && box[i][2] - l < 4 * AH) {
			for (int k = box[i][4]; k <= box[i][5]; k++) {
				if (img[k][l] == 0) {
					left_label = label[k][l];
					left_dist = box[i][2] - l;
					break;
				}
			}
			if (left_label != 0) {
				//std::cout << "\nLeft : " << left_label << " " << CCpixels[left_label];
				clr[int(Data[left_label][3])] += CCpixels[left_label];
				break;
			}
			l--;
		}
		if (left_dist < minDist && left_label != 0) {
			minDist = left_dist;
			distLabel = left_label;
		}

		//Right nbr...
		l = box[i][3] + 1;
		while (l < c && l - box[i][3] < 4 * AH) {
			for (int k = box[i][4]; k <= box[i][5]; k++) {
				if (img[k][l] == 0) {
					right_label = label[k][l];
					right_dist = l - box[i][3];
					break;
				}
			}
			if (right_label != 0) {
				//std::cout << "\nRight : " << right_label << " " << CCpixels[right_label];
				clr[int(Data[right_label][3])] += CCpixels[right_label];
				break;
			}
			l++;
		}
		if (right_dist < minDist && right_label != 0) {
			minDist = right_dist;
			distLabel = right_label;
		}

		int max = clr[1], max_clr = 1;
		for (int k = 2; k < 9; k++) {
			if (max < clr[k]) {
				max = clr[k];
				max_clr = k;
			}
		}
		//std::cout << "\nMaximum colour = " << max_clr;
		/*std::cout << "\nLabel " << i << " : Old = " << color << " New = " << Data[distLabel][3];
		if (color != Data[distLabel][3] && Data[distLabel][3] != 0) {
		std::cout << "\nnearest color = " << Data[distLabel][3] << " Label = " << distLabel;
		change[i][0] = 1;
		change[i][1] = Data[distLabel][3];
		}	*/
		//std::cout << "\nLabel " << i << " : Old = " << color << " New = " << max_clr;
		if (color != max_clr) {
			//std::cout << "\nnearest color = " << max_clr;
			change[i][0] = 1;
			change[i][1] = max_clr;
		}

	}

	//--------------------Re-color CC's----------------------
	for (int i = 0; i < labelNum; i++) {
		if (Data[i][4] == 0) continue;
		if (change[i][0] == 1)
			Data[i][3] = change[i][1];
	}

	std::cout << "\nhere!";
	delete[] clr;
}

void connectNbrs(int **img, int **label, int **box, double **boundingBox, double **Data, Point **ends, Mat pointers) {
	//Set up the connection list...
	int **connections = new int*[labelNum];
	for (int i = 0; i < labelNum; i++) {
		connections[i] = new int[2];
		for (int j = 0; j < 2; j++) {
			connections[i][0] = i;
			if (box[i][2] == -1 || Data[i][5] == 0)
				connections[i][1] = -1;
			else
				connections[i][1] = 0;
		}
	}

	//Connect CCs pairwise...
	int factor = 2;
	for (int i = 1; i < labelNum; i++) {
		if (box[i][2] == -1 || Data[i][5] == 0 || boundingBox[i][0] == 2) continue; //exclude the removed CCs and non horizontal ones...

		//Select the ROI...
		int **roi;
		int h = box[i][5] - box[i][4] + 1;
		int w = box[i][3] - box[i][2] + 1;
		roi = new int*[h];
		for (int k = 0; k < h; k++) {
			roi[k] = new int[w];
			for (int l = 0; l < w; l++) {
				if (img[box[i][4] + k][box[i][2] + l] == 175 && label[box[i][4] + k][box[i][2] + l] == i)
					roi[k][l] = 0;
				else
					roi[k][l] = 255;
			}
		}
		
		//Define neighborhood...
		int imin, imax, jmin, jmax;
		if (boundingBox[i][0] == 1) { //Upright...
			if (h < w) { //horizontal...
				imin = Data[i][1] - 2.5*AH;
				imax = Data[i][1] + 2.5*AH;
				jmin = Data[i][0];
				jmax = Data[i][0] + 5 * AH;
			}
			if (w < h) { //vertical...
				imin = Data[i][1];
				imax = Data[i][1] + 5 * AH;
				jmin = Data[i][0] - 2.5*AH;
				jmax = Data[i][0] + 2.5*AH;
			}
		}
		if (boundingBox[i][0] == 2) {  //diagonal...
			double ht = sqrt((boundingBox[i][5] - boundingBox[i][1])*(boundingBox[i][5] - boundingBox[i][1]) + (boundingBox[i][6] - boundingBox[i][2])*(boundingBox[i][6] - boundingBox[i][2]));
			double wd = sqrt((boundingBox[i][3] - boundingBox[i][5])*(boundingBox[i][3] - boundingBox[i][5]) + (boundingBox[i][4] - boundingBox[i][6])*(boundingBox[i][4] - boundingBox[i][6]));

			if (ht < wd) {  //nearly horizontal...
				imin = Data[i][1] - 2.5*AH;
				imax = Data[i][1] + 2.5*AH;
				jmin = Data[i][0];
				jmax = Data[i][0] + 5 * AH;
			}
			if (wd < ht) { //nearly vertical...
				imin = Data[i][1];
				imax = Data[i][1] + 5 * AH;
				jmin = Data[i][0] - 2.5*AH;
				jmax = Data[i][0] + 2.5*AH;
			}
		}

		if (imin < 0) imin = 0;
		if (imax >= r) imax = r - 1;
		if (jmin < 0) jmin = 0;
		if (jmax >= c) jmax = c - 1;
		
		//List the distinct neighboring labels...
		std::vector<int> nbrs;
		for (int k = imin; k <= imax; k++) {
			for (int l = jmin; l <= jmax; l++) {
				int nbrlbl = label[k][l];
				if (img[k][l] == 255 || nbrlbl == i || Data[nbrlbl][5] == 0) continue;
				int newnbr = 1;
				for (int count = 0; count < int(nbrs.size()); count++) {
					if (nbrs[count] == nbrlbl) {
						newnbr = 0;
						break;
					}
				}
				if (newnbr) { //Check if start point is in the nbrhood...
					if (ends[nbrlbl][0].x >= jmin && ends[nbrlbl][0].x <= jmax && ends[nbrlbl][0].y >= imin && ends[nbrlbl][0].y <= imax)
						nbrs.push_back(nbrlbl);
				}					
			}
		}
		int nbrno = int(nbrs.size());
		std::cout << "\nLabel " << i << " has " << nbrno << " nbrs...";

		if (nbrno < 2) continue;

		//Find the nearest CC...
		double mindist = 9999, dist, minlbl = 0;
		for (int k = 0; k < nbrno; k++) {
			dist = sqrt((ends[nbrs[k]][0].x - ends[i][1].x)*(ends[nbrs[k]][0].x - ends[i][1].x) + (ends[nbrs[k]][0].y - ends[i][1].y)*(ends[nbrs[k]][0].y - ends[i][1].y));
			/*if (ends[nbrs[k]][0].x - ends[i][1].x == 0) {
				dist = 0;
			}
			else {
				dist = (ends[nbrs[k]][0].y - ends[i][1].y) * (ends[nbrs[k]][0].x - ends[i][1].x);
			}*/
			if (dist < mindist) {
				mindist = dist;
				minlbl = nbrs[k];
			}
		}
		connections[i][1] = minlbl;

		//Delete memory...
		for (int k = 0; k < h; k++)
			delete[] roi[k];
		delete[] roi;
	}

	//Draw the connections...
	for (int i = 1; i < labelNum; i++) {
		if(connections[i][1]>0)
			arrowedLine(pointers, ends[connections[i][0]][1], ends[connections[i][1]][0], Scalar(0, 150, 255), 2, 8, 0, 0.05);
	}

	for (int i = 0; i < labelNum; i++)
		delete[] connections[i];
	delete[] connections;
}

//Consider the weighted average of angles of nbrs based on the distance of their centroids from the axis line... 
void near2axis(int **img, int **label, int **box, double **boundingBox, double **Data, Point **ends, Mat pointers) {
	int factor = 2.5;
	for (int i = 1; i < labelNum; i++) {
		if (Data[i][5] == 0) continue; //exclude the removed CCs...

		double slope = tan(Data[i][2] * CV_PI / 180);
		int imin, imax, jmin, jmax;
		jmin = box[i][2] - factor * AH;
		jmax = box[i][3] + factor * AH;
		imin = box[i][4] - factor * AH;
		imax = box[i][5] + factor * AH;

		if (imin < 0) imin = 0;
		if (imax >= r) imax = r - 1;
		if (jmin < 0) jmin = 0;
		if (jmax >= c) jmax = c - 1;

		//List the distinct neighboring labels...
		std::vector<int> nbrs;
		for (int k = imin; k <= imax; k++) {
			for (int l = jmin; l <= jmax; l++) {
				if (img[k][l] == 255 || label[k][l] == i || Data[label[k][l]][5] == 0) continue;
				int newnbr = 1;
				for (int count = 0; count<int(nbrs.size()); count++) {
					if (nbrs[count] == label[k][l]) {
						newnbr = 0;
						break;
					}
				}
				if (newnbr)
					nbrs.push_back(label[k][l]);
			}
		}
		int nbrno = int(nbrs.size());
		std::cout << "\nLabel " << i << " has " << nbrno << " nbrs...";
		if (nbrno == 0) continue;

		double *dist = new double[nbrno];
		double *weights = new double[nbrno];
		double maxdist = 0;
		for (int k = 0; k < nbrno; k++) {
			dist[k] = abs(slope*(Data[i][0] - Data[nbrs[k]][0]) - (Data[i][1] - Data[nbrs[k]][1])) / sqrt(1 + slope*slope);
			if (dist[k] > maxdist) maxdist = dist[k];
		}

		double newangle = 0, totalwt = 0;
		double newslope;
		//std:cout << "\nnbr-dist   wt      angle:";
		for (int k = 0; k < nbrno; k++) {
			weights[k] = 1 - (dist[k] / maxdist);
			totalwt += weights[k];
			newangle += weights[k] * Data[nbrs[k]][2];
			//std::cout << "\n" << dist[k] << "   " << weights[k] << "    " << Data[nbrs[k]][2];
		}
		newangle /= totalwt;
		newangle = int(newangle) - (int(newangle) % 5);
		newslope = tan(newangle * CV_PI / 180);
		std::cout << "\nLabel " << i << " : Old angle = " << Data[i][2] << " -> new angle = " << newangle;

		//Draw new axis line...
		Point endpts[2];
		if (newangle == 0 || newangle == 180) {
			endpts[0] = Point(box[i][2], Data[i][1]);
			endpts[1] = Point(box[i][3], Data[i][1]);
		}
		else if (newangle == 90) {
			endpts[0] = Point(Data[i][0], box[i][4]);
			endpts[1] = Point(Data[i][0], box[i][5]);
		}
		else {
			endpts[0].x = box[i][2];
			endpts[0].y = int(Data[i][1] + newslope*(box[i][2] - Data[i][0]));
			if (endpts[0].y < box[i][4]) {
				endpts[0].x = int(Data[i][0] + (box[i][4] - Data[i][1]) / newslope);
				endpts[0].y = box[i][4];
			}
			else if (endpts[0].y > box[i][5]) {
				endpts[0].x = int(Data[i][0] + (box[i][5] - Data[i][1]) / newslope);
				endpts[0].y = box[i][5];
			}

			endpts[1].x = box[i][3];
			endpts[1].y = int(Data[i][1] + newslope*(box[i][3] - Data[i][0]));
			if (endpts[1].y < box[i][4]) {
				endpts[1].x = int(Data[i][0] + (box[i][4] - Data[i][1]) / newslope);
				endpts[1].y = box[i][4];
			}
			else if (endpts[1].y > box[i][5]) {
				endpts[1].x = int(Data[i][0] + (box[i][5] - Data[i][1]) / newslope);
				endpts[1].y = box[i][5];
			}

		}


		//arrowedLine(pointers, endpts[0], endpts[1], Scalar(255, 0, 75), 1.75, 8, 0, 0.05);
		//Compare span...
		double orgspan, newspan=0;
		int **roi;
		int h = box[i][5] - box[i][4] + 1;
		int w = box[i][3] - box[i][2] + 1;
		roi = new int*[h];
		for (int k = 0; k < h; k++) {
			roi[k] = new int[w];
			for (int l = 0; l < w; l++) {
				if (img[box[i][4] + k][box[i][2] + l] == 175 && label[box[i][4] + k][box[i][2] + l] == i && Data[i][5] != 0)
					roi[k][l] = 0;
				else
					roi[k][l] = 255;
			}
		}
		double *data = Data[i];
		int *box_i = box[i];
		orgspan = span(roi, h, w, data, box_i, ends[i], Data[i][2]);
		newspan = span(roi, h, w, data, box_i, endpts, newangle);
		std::cout << "\nOld span = " << orgspan << " new span = " << newspan;

		if (orgspan < 1.0*newspan ) {
			//arrowedLine(pointers, endpts[0], endpts[1], Scalar(255, 0, 75), 2, 8, 0, 0.05);
			Data[i][2] = newangle;
			Data[i][4] = newspan;
			ends[i][0] = endpts[0];
			ends[i][1] = endpts[1];
		}
		else
			//arrowedLine(pointers, ends[i][0], ends[i][1], Scalar(255, 0, 175), 2, 8, 0, 0.05);

		//Delete memory...
		for (int k = 0; k < h; k++)
			delete[] roi[k];
		delete[] roi;
		
		delete[] dist;

	}
}

double Hausdorff(int **A, int **B, int ra, int ca, int rb, int cb, int *boxA, int *boxB) {
	double min, max = 0, d_AB, d_BA, dist;

	//d(A,B)...
	for (int i = 0; i < ra; i++) {
		for (int j = 0; j < ca; j++) {
			if (A[i][j] == 255) continue;
			min = 9999;
			for (int k = 0; k < rb; k++) {
				for (int l = 0; l < cb; l++) {
					if (B[k][l] == 255) continue;
					dist = sqrt(((boxA[4] + i) - (boxB[4] + k))*((boxA[4] + i) - (boxB[4] + k)) + ((boxA[2] + j) - (boxB[2] + l))*((boxA[2] + j) - (boxB[2] + l)));
					if (dist < min) min = dist;
				}
			}
			if (min > max) max = min;
		}
	}
	d_AB = max;

	//d(B,A)...
	max = 0;
	for (int i = 0; i < rb; i++) {
		for (int j = 0; j < cb; j++) {
			if (B[i][j] == 255) continue;
			min = 9999;
			for (int k = 0; k < ra; k++) {
				for (int l = 0; l < ca; l++) {
					if (A[k][l] == 255) continue;
					dist = sqrt(((boxB[4] + i) - (boxA[4] + k))*((boxB[4] + i) - (boxA[4] + k)) + ((boxB[2] + j) - (boxA[2] + l))*((boxB[2] + j) - (boxA[2] + l)));
					if (dist < min) min = dist;
				}
			}
			if (min > max) max = min;
		}
	}
	d_BA = max;

	if (d_AB >= d_BA) return(d_AB);
	if (d_AB < d_BA) return(d_BA);
}

//Consider the weighted average of angles of nbrs based on the Hausdorff distance between CCs...
void HausdorffNbrs(int **img, int **label, int **box, double **boundingBox, double **Data, Point **ends, Mat pointers) {
	int factor = 3.5;
	for (int i = 1; i < labelNum; i++) {
		if (box[i][2] == -1 || Data[i][5] == 0) continue; //exclude the removed CCs...

		//Select the ROI...
		int **roi;
		int h = box[i][5] - box[i][4] + 1;
		int w = box[i][3] - box[i][2] + 1;
		roi = new int*[h];
		for (int k = 0; k < h; k++) {
			roi[k] = new int[w];
			for (int l = 0; l < w; l++) {
				if (img[box[i][4] + k][box[i][2] + l] == 175 && label[box[i][4] + k][box[i][2] + l] == i )
					roi[k][l] = 0;
				else
					roi[k][l] = 255;
			}
		}

		int imin, imax, jmin, jmax;
		jmin = box[i][2] - factor * AH;
		jmax = box[i][3] + factor * AH;
		imin = box[i][4] - factor * AH;
		imax = box[i][5] + factor * AH;

		if (imin < 0) imin = 0;
		if (imax >= r) imax = r - 1;
		if (jmin < 0) jmin = 0;
		if (jmax >= c) jmax = c - 1;

		//List the distinct neighboring labels...
		std::vector<int> nbrs;
		for (int k = imin; k <= imax; k++) {
			for (int l = jmin; l <= jmax; l++) {
				if (img[k][l] == 255 || label[k][l] == i || Data[label[k][l]][5] == 0) continue;
				int newnbr = 1;
				for (int count = 0; count<int(nbrs.size()); count++) {
					if (nbrs[count] == label[k][l]) {
						newnbr = 0;
						break;
					}
				}
				if (newnbr)
					nbrs.push_back(label[k][l]);
			}
		}
		int nbrno = int(nbrs.size());
		//std::cout << "\nLabel " << i << " has " << nbrno << " nbrs...";

		//Find the Hausdorff distance between ROI and nbrs...
		double *dist = new double[nbrno];
		double maxdist = 0;
		double totalwt = 0, newangle = 0, newslope;

		if (nbrno == 0) continue;

		if (nbrno == 1) {
			int **roinbr;
			int hnbr = box[nbrs[0]][5] - box[nbrs[0]][4] + 1;
			int wnbr = box[nbrs[0]][3] - box[nbrs[0]][2] + 1;
			roinbr = new int*[hnbr];
			for (int k = 0; k < hnbr; k++) {
				roinbr[k] = new int[wnbr];
				for (int l = 0; l < wnbr; l++) {
					if (img[box[nbrs[0]][4] + k][box[nbrs[0]][2] + l] == 175 && label[box[nbrs[0]][4] + k][box[nbrs[0]][2] + l] == nbrs[0] && Data[nbrs[0]][5] != 0)
						roinbr[k][l] = 0;
					else
						roinbr[k][l] = 255;
				}
			}
			//Find distance...
			dist[0] = Hausdorff(roi, roinbr, h, w, hnbr, wnbr, box[i], box[nbrs[0]]);

			if (dist[0] < 100) {
				newangle = Data[nbrs[0]][2];
				newslope = tan(newangle * CV_PI / 180);
			}
			else
				continue;

			//Delete...
			for (int k = 0; k < hnbr; k++)
				delete[] roinbr[k];
			delete[] roinbr;
		}

		if (nbrno >= 2) {
			for (int kk = 0; kk < nbrno; kk++) {
				//Define nbrROI...
				int **roinbr;
				int hnbr = box[nbrs[kk]][5] - box[nbrs[kk]][4] + 1;
				int wnbr = box[nbrs[kk]][3] - box[nbrs[kk]][2] + 1;
				roinbr = new int*[hnbr];
				for (int k = 0; k < hnbr; k++) {
					roinbr[k] = new int[wnbr];
					for (int l = 0; l < wnbr; l++) {
						if (img[box[nbrs[kk]][4] + k][box[nbrs[kk]][2] + l] == 175 && label[box[nbrs[kk]][4] + k][box[nbrs[kk]][2] + l] == nbrs[kk] && Data[nbrs[kk]][5] != 0)
							roinbr[k][l] = 0;
						else
							roinbr[k][l] = 255;
					}
				}

				//Find distance...
				dist[kk] = Hausdorff(roi, roinbr, h, w, hnbr, wnbr, box[i], box[nbrs[kk]]);
				if (dist[kk] > maxdist) maxdist = dist[kk];

				//std::cout << "\nDistance between label " << i << " and label " << nbrs[kk] << " = " << dist[kk];

				//Delete...
				for (int k = 0; k < hnbr; k++)
					delete[] roinbr[k];
				delete[] roinbr;
			}

			//Find the weighted average angle of the nbrs...
			double *weights = new double[nbrno];
			//std::cout << "\nWeights :";
			for (int k = 0; k < nbrno; k++) {
				weights[k] = 1 - (dist[k] / maxdist);
				totalwt += weights[k];
				newangle += weights[k] * Data[nbrs[k]][2];
				//std::cout << "\nnbr " << nbrs[k] << ":" << Data[nbrs[k]][0] << " , " << Data[nbrs[k]][1] << " , " << Data[nbrs[k]][2] << " , " << Data[nbrs[k]][3] << " , " << Data[nbrs[k]][5];
				//std::cout << "\n" << weights[k] << " , " << Data[nbrs[k]][2] << "  //  ";
			}
			//std::cout << " Total = " << totalwt << " newangle = " << newangle;
			newangle /= totalwt;
			//newangle = int(newangle) - (int(newangle) % 5);
			newslope = tan(newangle * CV_PI / 180);
		}
		
		//Draw new axis line...
		Point endpts[2];
		if (newangle == 0 || newangle == 180) {
			endpts[0] = Point(box[i][2], Data[i][1]);
			endpts[1] = Point(box[i][3], Data[i][1]);
		}
		else if (newangle == 90) {
			endpts[0] = Point(Data[i][0], box[i][4]);
			endpts[1] = Point(Data[i][0], box[i][5]);
		}
		else {
			endpts[0].x = box[i][2];
			endpts[0].y = int(Data[i][1] + newslope*(box[i][2] - Data[i][0]));
			if (endpts[0].y < box[i][4]) {
				endpts[0].x = int(Data[i][0] + (box[i][4] - Data[i][1]) / newslope);
				endpts[0].y = box[i][4];
			}
			else if (endpts[0].y > box[i][5]) {
				endpts[0].x = int(Data[i][0] + (box[i][5] - Data[i][1]) / newslope);
				endpts[0].y = box[i][5];
			}

			endpts[1].x = box[i][3];
			endpts[1].y = int(Data[i][1] + newslope*(box[i][3] - Data[i][0]));
			if (endpts[1].y < box[i][4]) {
				endpts[1].x = int(Data[i][0] + (box[i][4] - Data[i][1]) / newslope);
				endpts[1].y = box[i][4];
			}
			else if (endpts[1].y > box[i][5]) {
				endpts[1].x = int(Data[i][0] + (box[i][5] - Data[i][1]) / newslope);
				endpts[1].y = box[i][5];
			}

		}
		//arrowedLine(pointers, endpts[0], endpts[1], Scalar(255, 0, 255), 2, 8, 0, 0.05);

		double orgspan = span(roi, h, w, Data[i], box[i], ends[i], Data[i][2]);
		double newspan = span(roi, h, w, Data[i], box[i], endpts, newangle);
		std::cout << "\nLabel " << i << ": Old span = " << orgspan << " new span = " << newspan;

		if (nbrno >= 2 &&  newspan > 1.0*orgspan) {
			//std::cout << "\nOld angle = " << Data[i][2] << " -> new angle = " << newangle;
			Data[i][2] = newangle;
			Data[i][4] = newspan;
			ends[i][0] = endpts[0];
			ends[i][1] = endpts[1];
		}
		else if (newspan > 1.0*orgspan) {
			//arrowedLine(pointers, endpts[0], endpts[1], Scalar(255, 0, 255), 2, 8, 0, 0.05);
			std::cout << "\nOld angle = " << Data[i][2] << " -> new angle = " << newangle;
			Data[i][2] = newangle;
			Data[i][4] = newspan;
			ends[i][0] = endpts[0];
			ends[i][1] = endpts[1];
		}
		else {
			std::cout << "\nOld angle = " << Data[i][2] << " -> new angle = " << Data[i][2] << "...No change";
		}

		//Delete memory...
		for (int k = 0; k < h; k++)
			delete[] roi[k];
		delete[] roi;

	}
}

//Connect the horizontal boxes based on proximity of overlapping CCs...
void connectHNbrs(int **img, int **label, int **box, double **boundingBox, double **Data, Point **ends, double **connect, Mat pointers) {
	//Set up the connection list...
	double **connections = new double*[labelNum];
	for (int i = 0; i < labelNum; i++) {
		connections[i] = new double[3];
		for (int j = 0; j < 3; j++) {
			connections[i][0] = i;
			if (box[i][2] == -1 || Data[i][5] == 0)
				connections[i][1] = -1;
			else
				connections[i][1] = 0;
			connections[i][2] = -1;
		}
	}

	//Connect CCs pairwise...
	int factor = 2;
	for (int i = 1; i < labelNum; i++) {
		if (box[i][2] == -1 || Data[i][5] == 0 || Data[i][3] != 1) continue; //exclude the removed CCs and non horizontal ones...

		//Select the ROI...
		int **roi;
		int h = box[i][5] - box[i][4] + 1;
		int w = box[i][3] - box[i][2] + 1;
		roi = new int*[h];
		for (int k = 0; k < h; k++) {
			roi[k] = new int[w];
			for (int l = 0; l < w; l++) {
				if (img[box[i][4] + k][box[i][2] + l] == 175 && label[box[i][4] + k][box[i][2] + l] == i)
					roi[k][l] = 0;
				else
					roi[k][l] = 255;
			}
		}

		//Define neighborhood...
		int imin, imax, jmin, jmax;
		imin = box[i][4] - 0.25*AH;
		imax = box[i][5] + 0.25*AH;
		jmin = box[i][2];
		jmax = c - 1;
		if (imin < 0) imin = 0;
		if (imax >= r) imax = r - 1;

		std::vector<int> nonhor;
		int countnonhor = 0, counthor = 0;
		for (int p = imin; p <= imax; p++) {
			for (int q = jmin; q <= jmax; q++) {
				if (img[p][q] != 255) { //inside the scope of roi...
					int lbl = label[p][q], flag = 0;
					for (int k = 0; k < int(nonhor.size()); k++) {
						if (nonhor[k] == lbl) {
							flag = 1;
							break;
						}
					}
					if (flag == 0 && (box[lbl][2] != -1 && Data[lbl][5] != 0 && Data[lbl][3] != 1 && lbl != i)) { //if valid component is of other orientation...
						nonhor.push_back(lbl);
						countnonhor++;
					}
					else if (flag == 0 && (box[lbl][2] != -1 && Data[lbl][5] != 0 && Data[lbl][3] == 1 && lbl != i)) {
						nonhor.push_back(lbl);
						counthor++;
					}
				}
			}
		}
		std::cout << "\nNo. of non-vertical components in range = " << countnonhor;
		std::cout << "\nNo. of vertical components in range = " << counthor;
		if (counthor <= countnonhor) continue;


		//List the distinct neighboring labels...
		std::vector<int> nbrs;
		for (int k = imin; k <= imax; k++) {
			for (int l = jmin; l <= jmax; l++) {
				int nbrlbl = label[k][l];
				if (img[k][l] == 255 || nbrlbl == i || Data[nbrlbl][5] == 0) continue;
				int newnbr = 1;
				for (int count = 0; count < int(nbrs.size()); count++) {
					if (nbrs[count] == nbrlbl) {
						newnbr = 0;
						break;
					}
				}
				if (newnbr) { //Check if nbr is horizontal...
					if (Data[nbrlbl][3] == 1)
						nbrs.push_back(nbrlbl);
				}
			}
		}
		int nbrno = int(nbrs.size());
		std::cout << "\nLabel " << i << " has " << nbrno << " nbrs...";

		//Find the nearest CC with at least 30% overlap...
		double mindist = 9999, dist;
		int minlbl = 0;
		double *prop = new double[nbrno];
		for (int k = 0; k < nbrno; k++) {
			int lbl = nbrs[k];
			int total = 0, section = 0;
			for (int y = box[lbl][4]; y <= box[lbl][5]; y++) {
				for (int x = box[lbl][2]; x <= box[lbl][3]; x++) {
					if (label[y][x] != lbl) continue;
					total++;
					if (y >= imin && y <= imax && x >= jmin && x <= jmax) section++;
					prop[k] = (double)section / total;
				}
			}
		}
		for (int k = 0; k < nbrno; k++) {
			if (prop[k] < 0.25) {
				for (int l = k; l < nbrno - 1; l++) {
					nbrs[l] = nbrs[l + 1];
					prop[l] = prop[l + 1];
				}
				k--;
				nbrno--;
				continue;
			}
		}
		for (int k = 0; k < nbrno; k++) {
			dist = sqrt((ends[nbrs[k]][0].x - ends[i][1].x)*(ends[nbrs[k]][0].x - ends[i][1].x) + (ends[nbrs[k]][0].y - ends[i][1].y)*(ends[nbrs[k]][0].y - ends[i][1].y));
			if (dist < mindist) {
				mindist = dist;
				minlbl = nbrs[k];
			}
		}

		//Check for previous connection...
		int tie = 0, tie_nbr = 0;
		for (int k = 0; k < labelNum; k++) {
			if (connections[k][1] == minlbl) {
				tie = 1;
				tie_nbr = k;
				break;
			}
		}
		if (tie) {
			if (connections[tie_nbr][2] > mindist) {
				connections[i][1] = minlbl;
				connections[i][2] = mindist;
				//Reset previous connection...
				connections[tie_nbr][1] = 0;
				connections[tie_nbr][2] = -1;
			}
		}
		if (!tie) {
			connections[i][1] = minlbl;
			connections[i][2] = mindist;
		}
			

		//Delete memory...
		for (int k = 0; k < h; k++)
			delete[] roi[k];
		delete[] roi;
	}

	//Draw the connections...
	for (int i = 1; i < labelNum; i++) {
		if (box[i][2] == -1 || Data[i][5] == 0 || Data[i][3] != 1) continue;
		//line(pointers, ends[i][0], ends[i][1], Scalar(0, 0, 150), 2, 8, 0);
		if (connections[i][1] > 0) {
			//arrowedLine(pointers, ends[int(connections[i][0])][1], ends[int(connections[i][1])][0], Scalar(0, 150, 10), 2, 8, 0, 0.05);
			connect[i][1] = connections[i][1];
			connect[i][2] = connections[i][2];
		}
	}

	for (int i = 0; i < labelNum; i++)
		delete[] connections[i];
	delete[] connections;
}

//Connect the vertical boxes based on proximity of overlapping CCs...
void connectVNbrs(int **img, int **label, int **box, double **boundingBox, double **Data, Point **ends, double **connect, Mat pointers) {
	//Set up the connection list...
	int **connections = new int*[labelNum];
	for (int i = 0; i < labelNum; i++) {
		connections[i] = new int[3];
		for (int j = 0; j < 3; j++) {
			connections[i][0] = i;
			if (box[i][2] == -1 || Data[i][5] == 0)
				connections[i][1] = -1;
			else
				connections[i][1] = 0;
			connections[i][2] = -1;
		}
	}

	//Connect CCs pairwise...
	int factor = 2;
	for (int i = 1; i < labelNum; i++) {
		if (box[i][2] == -1 || Data[i][5] == 0 || Data[i][3] != 2) continue; //exclude the removed CCs and non vertical ones...

		//Select the ROI...
		int **roi;
		int h = box[i][5] - box[i][4] + 1;
		int w = box[i][3] - box[i][2] + 1;
		roi = new int*[h];
		for (int k = 0; k < h; k++) {
			roi[k] = new int[w];
			for (int l = 0; l < w; l++) {
				if (img[box[i][4] + k][box[i][2] + l] == 175 && label[box[i][4] + k][box[i][2] + l] == i)
					roi[k][l] = 0;
				else
					roi[k][l] = 255;
			}
		}

		//Define neighborhood...
		int imin, imax, jmin, jmax;
		imin = box[i][5];
		imax = r - 1;
		jmin = box[i][2] - 0.25*AH;
		if (jmin < 0) jmin = 0;
		jmax = box[i][3] + 0.25*AH;
		if (jmax >= c) jmax = c - 1;

		std::vector<int> nonvert;
		int countnonvert = 0, countvert = 0;
		for (int p = imin; p <= imax; p++) {
			for (int q = jmin; q <= jmax; q++) {
				if (img[p][q] != 255) { //inside the scope of roi...
					int lbl = label[p][q], flag = 0;
					for (int k = 0; k < int(nonvert.size()); k++) {
						if (nonvert[k] == lbl) {
							flag = 1;
							break;
						}
					}
					if (flag == 0 && (box[lbl][2] != -1 && Data[lbl][5] != 0 && Data[lbl][3] != 2 && lbl != i)) { //if valid component is of other orientation...
						nonvert.push_back(lbl);
						countnonvert++;
					}
					else if (flag == 0 && (box[lbl][2] != -1 && Data[lbl][5] != 0 && Data[lbl][3] == 2 && lbl != i)) {
						nonvert.push_back(lbl);
						countvert++;
					}
				}
			}
		}
		std::cout << "\nNo. of non-vertical components in range = " << countnonvert;
		std::cout << "\nNo. of vertical components in range = " << countvert;
		if (countvert <= countnonvert) continue;


		//List the distinct neighboring labels...
		std::vector<int> nbrs;
		for (int k = imin; k <= imax; k++) {
			for (int l = jmin; l <= jmax; l++) {
				int nbrlbl = label[k][l];
				if (img[k][l] == 255 || nbrlbl == i || Data[nbrlbl][5] == 0) continue;
				int newnbr = 1;
				for (int count = 0; count < int(nbrs.size()); count++) {
					if (nbrs[count] == nbrlbl) {
						newnbr = 0;
						break;
					}
				}
				if (newnbr) { //Check if nbr is vertical...
					if (Data[nbrlbl][3] == 2)
						nbrs.push_back(nbrlbl);
				}
			}
		}
		int nbrno = int(nbrs.size());

		//Find the nearest CC...
		double mindist = 9999, dist;
		int minlbl = 0, tie = 0, tie_nbr = 0;
		double *prop = new double[nbrno];
		for (int k = 0; k < nbrno; k++) {
			int lbl = nbrs[k];
			int total = 0, section = 0;
			for (int y = box[lbl][4]; y <= box[lbl][5]; y++) {
				for (int x = box[lbl][2]; x <= box[lbl][3]; x++) {
					if (label[y][x] != lbl) continue;
					total++;
					if (y >= imin && y <= imax && x >= jmin && x <= jmax) section++;
					prop[k] = (double)section / total;
				}
			}
		}
		for (int k = 0; k < nbrno; k++) {
			if (prop[k] < 0.3) {
				for (int l = k; l < nbrno - 1; l++) {
					nbrs[l] = nbrs[l + 1];
					prop[l] = prop[l + 1];
				}
				k--;
				nbrno--;
				continue;
			}
		}
		for (int k = 0; k < nbrno; k++) {
			if (Data[nbrs[k]][2] <= 90 && Data[i][2] <= 90) {
				dist = sqrt((ends[nbrs[k]][0].x - ends[i][1].x)*(ends[nbrs[k]][0].x - ends[i][1].x) + (ends[nbrs[k]][0].y - ends[i][1].y)*(ends[nbrs[k]][0].y - ends[i][1].y));
			}
			else if (Data[nbrs[k]][2] > 90 && Data[i][2] <= 90) {
				dist = sqrt((ends[nbrs[k]][1].x - ends[i][1].x)*(ends[nbrs[k]][1].x - ends[i][1].x) + (ends[nbrs[k]][1].y - ends[i][1].y)*(ends[nbrs[k]][1].y - ends[i][1].y));
			}
			else if (Data[nbrs[k]][2] <= 90 && Data[i][2] > 90) {
				dist = sqrt((ends[nbrs[k]][0].x - ends[i][0].x)*(ends[nbrs[k]][0].x - ends[i][0].x) + (ends[nbrs[k]][0].y - ends[i][0].y)*(ends[nbrs[k]][0].y - ends[i][0].y));
			}
			else if (Data[nbrs[k]][2] > 90 && Data[i][2] > 90) {
				dist = sqrt((ends[nbrs[k]][1].x - ends[i][0].x)*(ends[nbrs[k]][1].x - ends[i][0].x) + (ends[nbrs[k]][1].y - ends[i][0].y)*(ends[nbrs[k]][1].y - ends[i][0].y));
			}

			if (dist < mindist) {
				mindist = dist;
				minlbl = nbrs[k];
			}
		}
		
		//Check for previous connection...
		for (int k = 0; k < labelNum; k++) {
			if (connections[k][1] == minlbl) {
				tie = 1;
				tie_nbr = k;
				break;
			}
		}
		if (tie) {
			if (connections[tie_nbr][2] > mindist) {
				connections[i][1] = minlbl;
				connections[i][2] = mindist;
				//Reset previous connection...
				connections[tie_nbr][1] = 0;
				connections[tie_nbr][2] = -1;
			}
		}
		if (!tie) {
			connections[i][1] = minlbl;
			connections[i][2] = mindist;
		}

		//Delete memory...
		for (int k = 0; k < h; k++)
			delete[] roi[k];
		delete[] roi;
	}

	//Draw the connections...
	for (int i = 1; i < labelNum; i++) {
		if (box[i][2] == -1 || Data[i][5] == 0 || Data[i][3] != 2) continue;
		//line(pointers, ends[i][0], ends[i][1], Scalar(0, 0, 150), 2, 8, 0);
		if (connections[i][1] > 0) {
			/*if (Data[connections[i][0]][2] <= 90 && Data[connections[i][1]][2] <= 90)
				arrowedLine(pointers, ends[connections[i][0]][1], ends[connections[i][1]][0], Scalar(0, 150, 10), 2, 8, 0, 0.05);
			else if (Data[connections[i][0]][2] <= 90 && Data[connections[i][1]][2] > 90)
				arrowedLine(pointers, ends[connections[i][0]][1], ends[connections[i][1]][1], Scalar(0, 150, 10), 2, 8, 0, 0.05);
			else if (Data[connections[i][0]][2] > 90 && Data[connections[i][1]][2] <= 90)
				arrowedLine(pointers, ends[connections[i][0]][0], ends[connections[i][1]][0], Scalar(0, 150, 10), 2, 8, 0, 0.05);
			else if (Data[connections[i][0]][2] > 90 && Data[connections[i][1]][2] > 90)
				arrowedLine(pointers, ends[connections[i][0]][0], ends[connections[i][1]][1], Scalar(0, 150, 10), 2, 8, 0, 0.05);*/
			connect[i][1] = connections[i][1];
			connect[i][2] = connections[i][2];
		}
	}

	for (int i = 0; i < labelNum; i++)
		delete[] connections[i];
	delete[] connections;
}

//Centroid of end strip...
Point  LocalCentroid(int **img, int **label, int lbl, Point N, Point S, Point W, Point E) {
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

	std::cout << "\nCentroids: (" << cntr.x << " , " << cntr.y << ")" ;
	return(cntr);
}

//Connect the diagonal boxes based on proximity of overlapping CCs...
void connecDNbrs(int **img, int **label, int **box, double **boundingBox, double **Data, Point **ends, double **connect, Mat pointers) {
	//Set up the connection list...
	int **connections = new int*[labelNum];
	for (int i = 0; i < labelNum; i++) {
		connections[i] = new int[4];
		for (int j = 0; j < 4; j++) {
			connections[i][0] = i;
			if (box[i][2] == -1 || Data[i][5] == 0)
				connections[i][1] = -1;
			else
				connections[i][1] = 0;
			connections[i][2] = -1;
			connections[i][2] = 0;
		}
	}

	//Connect CCs pairwise...
	int factor = 15;
	for (int i = 1; i < labelNum; i++) {
		if (box[i][2] == -1 || Data[i][5] == 0 || Data[i][3] != 4) continue; //exclude the removed CCs and non diagonal ones...

		Point centroidL, centroidR, N, S, W, E, vec_org;
		N = Point(boundingBox[i][2], boundingBox[i][1]);
		E = Point(boundingBox[i][6], boundingBox[i][5]);
		W = Point(boundingBox[i][2] - factor, boundingBox[i][1] + factor);
		S = Point(boundingBox[i][6] - factor, boundingBox[i][5] + factor);
		centroidR = LocalCentroid(img, label, i, N, S, W, E);

		N = Point(boundingBox[i][8] + factor, boundingBox[i][7] - factor);
		E = Point(boundingBox[i][4] + factor, boundingBox[i][3] - factor);
		W = Point(boundingBox[i][8], boundingBox[i][7]);
		S = Point(boundingBox[i][4], boundingBox[i][3]);
		centroidL = LocalCentroid(img, label, i, N, S, W, E);

		std::cout << "\nCentroids: (" << centroidL.x << " , " << centroidL.y << ") ; (" << centroidR.x << " , " << centroidR.y << ")";
		line(pointers, centroidL, centroidR, Scalar(0, 0, 150), 2, 8, 0);

		vec_org = centroidL - centroidR;
		int angle_org = atan2(vec_org.y, vec_org.x) * 180 / CV_PI;
		double m = tan(angle_org*CV_PI / 180);

		//Define neighborhood...
		int imin, imax, jmin, jmax;
		imin = 0;
		imax = ends[i][0].y;
		jmin = ends[i][0].x;
		jmax = -1 * int( (centroidL.y / m ) - centroidL.x );

		if (jmax >= c) jmax = c - 1;

		std::cout << "\n" << imin << "," << imax << "," << jmin << "," << jmax;
		
		//List the distinct neighboring labels...
		double c1 = (centroidL.y - m*centroidL.x) - 1.55 * AH;
		double c2 = (centroidL.y - m*centroidL.x) + 1.55 * AH;
		std::cout << "\nparallel lines at " << c1 << " , " << c2 << " at angle = " << angle_org << " , slope = " << m;

		std::vector<int> nondiag;
		int countnondiag = 0, countdiag = 0;
		for (int p = imin; p <= imax; p++) {
			for (int q = jmin; q <= jmax; q++) {
				if (p < 0 || p >= r || q < 0 || q >= c)
					continue;
				if ( img[p][q] != 255) { //inside the scope of roi...
					int lbl = label[p][q] , flag = 0;
					for (int k = 0; k < int(nondiag.size()); k++) {
						if (nondiag[k] == lbl) {
							flag = 1;
							break;
						}
					}
					if (flag == 0 && (box[lbl][2] != -1 && Data[lbl][5] != 0 && Data[lbl][3] != 4)) { //if valid component is of other orientation...
						nondiag.push_back(lbl);
						countnondiag++;
					}
					else if(flag == 0 && (box[lbl][2] != -1 && Data[lbl][5] != 0 && Data[lbl][3] == 4)){
						nondiag.push_back(lbl);
						countdiag++;
					}
				}
			}
		}
		std::cout << "\nNo. of non-diagonal components in range = " << countnondiag;
		std::cout << "\nNo. of diagonal components in range = " << countdiag;
		if (countdiag <= countnondiag) continue;

		std::vector<int> nbrs;
		for (int k = 1; k < labelNum; k++) {
			if (box[k][2] == -1 || Data[k][5] == 0 || Data[k][3] != 4 || k == i) continue;

			//else if (ends[k][0].x < jmin || ends[k][0].x > jmax || ends[k][0].y < imin || ends[k][0].y > imax) continue;  //Axis not in nbd region...
			//else if (ends[k][0].y - m*ends[k][0].x <= imax - 3.5*AH || ends[k][0].y - m*ends[k][0].x >= imax + 3.5*AH) continue;  // axis not inside bdd region...
			N = Point(boundingBox[k][8] + factor, boundingBox[k][7] - factor);
			E = Point(boundingBox[k][4] + factor, boundingBox[k][3] - factor);
			W = Point(boundingBox[k][8], boundingBox[k][7]);
			S = Point(boundingBox[k][4], boundingBox[k][3]);
			Point centerL = LocalCentroid(img, label, k, N, S, W, E);


			if (centerL.x < jmin || centerL.x > jmax || centerL.y < imin || centerL.y > imax) continue;
			else if (centerL.y - m*centerL.x <= c1 || centerL.y - m*centerL.x >= c2) continue;  // axis not inside bdd region...
			//else if (ends[k][0].y - m*ends[k][0].x <= c1 || ends[k][0].y - m*ends[k][0].x >= c2) continue;  // axis not inside bdd region...
			//else if (Data[k][0] < jmin || Data[k][0] > jmax || Data[k][1] < imin || Data[k][1] > imax) continue; //Axis not in nbd region... 
			else {                               //Include only the diagonal nbrs in range...
				nbrs.push_back(k);
				std::cout << "  ;  " << nbrs[int(nbrs.size())-1];
			}
		}

		int nbrno = int(nbrs.size());
		std::cout << "\nLabel " << i << " has " << nbrno << " nbrs...";
		if (nbrno == 0) continue;		

		//Find the nbrs with 3 least gradients of centroids...
		int minlbl = 0;
		double DIST;
		double *mingrad = new double[3];
		double *mindist = new double[3];
		int *gradlbl = new int[3];
		double *nbrSlope = new double[nbrno];
		double *angle_change = new double[nbrno];
		double *angleLL = new double[nbrno];
		double *dist = new double[nbrno];

		for (int t = 0; t < 3; t++) {
			mingrad[t] = 9999;
			mindist[t] = 9999;
			gradlbl[t] = 0;
			angleLL[t] = 0;
			//nbrSlope[t] = 0;
		}

		std::cout << "\nLabel " << i << " has change in angles = ";
		for (int nb = 0; nb < nbrno; nb++) {
			Point centerL, centerR;
			int nbr = nbrs[nb];

			N = Point(boundingBox[nbr][2], boundingBox[nbr][1]);
			E = Point(boundingBox[nbr][6], boundingBox[nbr][5]);
			W = Point(boundingBox[nbr][2] - factor, boundingBox[nbr][1] + factor);
			S = Point(boundingBox[nbr][6] - factor, boundingBox[nbr][5] + factor);
			centerR = LocalCentroid(img, label, nbr, N, S, W, E);

			N = Point(boundingBox[nbr][8] + factor, boundingBox[nbr][7] - factor);
			E = Point(boundingBox[nbr][4] + factor, boundingBox[nbr][3] - factor);
			W = Point(boundingBox[nbr][8], boundingBox[nbr][7]);
			S = Point(boundingBox[nbr][4], boundingBox[nbr][3]);
			centerL = LocalCentroid(img, label, nbr, N, S, W, E);

			Point vecL, vecR, vec_nbr;
			vecL = centerL - centroidL;
			vecR = centerL - centroidR;
			vec_nbr = centerR - centerL;

			/*int angleL = 180 - atan2(vecL.y, vecL.x) * 180 / CV_PI;
			int angleR = 180 - atan2(vecR.y, vecR.x) * 180 / CV_PI;
			int angle_nbr = 180 - atan2(vec_nbr.y, vec_nbr.x) * 180 / CV_PI;

			nbrSlope[nb] = (abs(angle_org - angleR) + abs(angle_nbr - angleR) + abs(angle_org - angleL))*norm(vecR);
			std::cout <<  nbrSlope[nb] << " , ";
			dist = abs(angle_org - angleL)+norm(vecR);*/

			double angleLR = acos((vecL.dot(vecR)) / (norm(vecR)*norm(vecL))) * 180 / CV_PI;
			double angleRnbr = acos((vecR.dot(vec_nbr)) / (norm(vecR)*norm(vec_nbr))) * 180 / CV_PI;

			/*if (angleLR > 180) angleLR = 360 - angleLR;
			if (angleRnbr > 180) angleRnbr = 360 - angleRnbr;
*/
			angle_change[nb] = abs((angleLR + angleRnbr));
			dist[nb] = norm(vecR);
			angleLL[nb] = acos((vecL.dot(vec_nbr)) / (norm(vecL)*norm(vec_nbr))) * 180 / CV_PI;;
			std::cout << angle_change[nb] << " + " << dist[nb] << ";";
		}

		double maxdst = 0;
		for (int nb = 0; nb < nbrno; nb++) {
			if (dist[nb] > maxdst) {
				maxdst = dist[nb];
			}
		}

		for(int nb = 0; nb < nbrno; nb++){
			nbrSlope[nb] = 0.75*(angle_change[nb] / 180) + 0.25*(dist[nb] / maxdst);
			//nbrSlope[nb] = angle_change[nb];
			if (mingrad[0] > nbrSlope[nb]) {
				mingrad[0] = nbrSlope[nb];
				gradlbl[0] = nbrs[nb];
				mindist[0] = angleLL[nb];
			}
			else if (mingrad[1] > nbrSlope[nb]) {
				mingrad[1] = nbrSlope[nb];
				gradlbl[1] = nbrs[nb];
				mindist[1] = (angle_change[nb]) * dist[nb] ;
			}
			else if (mingrad[2] > nbrSlope[nb]) {
				mingrad[2] = nbrSlope[nb];
				gradlbl[2] = nbrs[nb];
				mindist[2] = (angle_change[nb]) * dist[nb] ;
			}
		}
		std::cout << "\nlabels:" << gradlbl[0] << "," << gradlbl[1] << "," << gradlbl[2];
		std::cout << "\ngrads:" << mingrad[0] << "," << mingrad[1] << "," << mingrad[2];
		std::cout << "\ndist:" << mindist[0] << "," << mindist[1] << "," << mindist[2];

		////Find the one of the 3 which is nearest...	
		//if (nbrno == 1) {
		//	gradlbl[1] = 0;
		//	gradlbl[2] = 0;
		//}
		//if (nbrno == 2) {
		//	if (mindist[1] < mindist[0]) {
		//		double t = mingrad[0];
		//		mingrad[0] = mingrad[1];
		//		mingrad[1] = t;

		//		int tmp = gradlbl[0];
		//		gradlbl[0] = gradlbl[1];
		//		gradlbl[1] = tmp;

		//		t = mindist[0];
		//		mindist[0] = mindist[1];
		//		mindist[1] = t;
		//	}
		//}
		//else if (nbrno > 2) {
		//	for (int k = 2; k > 0; k--) {
		//		for (int j = 0; j < k; j++) {
		//			if (mindist[j + 1] < mindist[j]) {
		//				/*double t1 = mingrad[k];
		//				mingrad[k] = mingrad[k + 1];
		//				mingrad[k + 1] = t1;*/

		//				int tmp = gradlbl[j];
		//				gradlbl[j] = gradlbl[j + 1];
		//				gradlbl[j + 1] = tmp;

		//				double t2 = mindist[j];
		//				mindist[j] = mindist[j + 1];
		//				mindist[j + 1] = t2;
		//			}
		//		}
		//	}
		//}
		
		//std::cout << "\nnewlabels:" << gradlbl[0] << "," << gradlbl[1] << "," << gradlbl[2];
		////std::cout << "\nnewgrads:" << mingrad[0] << "," << mingrad[1] << "," << mingrad[2];
		//std::cout << "\nnewdist:" << mindist[0] << "," << mindist[1] << "," << mindist[2];


		//Check for previous connection...
		int count = 0 , flag = 0;
		while (count < 1 && flag == 0) {
			if (gradlbl[count] == 0) {
				break;
			}
			int tie = 0, tie_nbr = 0;
			for (int k = 0; k < labelNum; k++) {
				if (connections[k][1] == gradlbl[count]) {
					tie = 1;
					tie_nbr = k;
					break;
				}
			}
			if (tie) {
				if (connections[tie_nbr][3]  > dist[gradlbl[count]]) {
					connections[i][1] = gradlbl[count];
					connections[i][2] = mindist[count];
					connections[i][3] = dist[gradlbl[count]];
					flag = 1;
					//Reset previous connection...
					connections[tie_nbr][1] = 0;
					connections[tie_nbr][2] = -1;
					connections[tie_nbr][3] = 0;
				}
				//else if (abs(connections[tie_nbr][2] - mindist[count]) <= 5 && 0.5*connections[tie_nbr][3] > dist[gradlbl[count]]) {
				//	connections[i][1] = gradlbl[count];
				//	connections[i][2] = mindist[count];
				//	connections[i][3] = dist[gradlbl[count]];
				//	flag = 1;
				//	//Reset previous connection...
				//	connections[tie_nbr][1] = 0;
				//	connections[tie_nbr][2] = -1;
				//	connections[tie_nbr][2] = 0;
				//}
				else {
					count++;
					continue;
				}
			}
			if (!tie) {
				connections[i][1] = gradlbl[count];
				connections[i][2] = mindist[count];
				connections[i][3] = dist[gradlbl[count]];
				flag = 1;
			}
		}	


		//int tie = 0, tie_nbr = 0;
		//for (int k = 0; k < labelNum; k++) {
		//	if (connections[k][1] == minlbl) {
		//		tie = 1;
		//		tie_nbr = k;
		//		break;
		//	}
		//}
		//if (tie) {
		//	if (connections[tie_nbr][2] > DIST) {
		//		connections[i][1] = minlbl;
		//		connections[i][2] = DIST;
		//		//Reset previous connection...
		//		connections[tie_nbr][1] = 0;
		//		connections[tie_nbr][2] = -1;
		//	}
		//}
		//if (!tie) {
		//	connections[i][1] = minlbl;
		//	connections[i][2] = DIST;
		//}

		//Delete memory...
		delete[] nbrSlope;
		delete[] mingrad;
		delete[] mindist;
		delete[] gradlbl;
		delete[] angle_change;
		delete[] dist;
		
	}

	//Draw the connections...
	for (int i = 1; i < labelNum; i++) {
		if (connections[i][1] > 0) {
			//arrowedLine(pointers, ends[connections[i][0]][1], ends[connections[i][1]][0], Scalar(0, 150, 10), 2, 8, 0, 0.05);
			connect[i][1] = connections[i][1];
			connect[i][2] = connections[i][3];
		}
	}

	for (int i = 0; i < labelNum; i++)
		delete[] connections[i];
	delete[] connections;
}


//Find the 4 CCs nearest to the axis line...
/*double *perpdist = new double[nbrno];
double **near4CC;
near4CC = new double*[4];
for (int k = 0; k < 4; k++) {
near4CC[k] = new double[3];
near4CC[k][0] = 0;
near4CC[k][1] = 9999;
near4CC[k][2] = 0;
}
for (int k = 0; k < nbrno; k++) {
Point p1 = ends[i][0];
Point p2 = ends[i][1];
Point p0 = ends[nbrs[k]][0];
perpdist[k] = abs((p2.y - p1.y)*p0.x - (p2.x - p1.x)*p0.y + p2.x*p1.y - p1.x*p2.y) / norm(p2 - p1);

if (perpdist[k] < near4CC[0][1]) {
near4CC[0][1] = perpdist[k];
near4CC[0][0] = nbrs[k];
}
else if (perpdist[k] < near4CC[1][1]) {
near4CC[1][1] = perpdist[k];
near4CC[1][0] = nbrs[k];
}
else if (perpdist[k] < near4CC[2][1]) {
near4CC[2][1] = perpdist[k];
near4CC[2][0] = nbrs[k];
}
else if (perpdist[k] < near4CC[3][1]) {
near4CC[3][1] = perpdist[k];
near4CC[3][0] = nbrs[k];
}
}
std::cout << "\nNearest nbrs of Label " << i << " : " << near4CC[0][0] << " , " << near4CC[1][0] << " , " << near4CC[2][0] << " , " << near4CC[3][0];
//Find the closest out of the 4...
double mindist = 9999;
int minlbl = 0;
for (int k = 0; k < 4; k++) {
if (near4CC[k][0] == 0) continue;
int near = near4CC[k][0];
near4CC[k][2] = norm(ends[near][0] - ends[i][1]);
if (near4CC[k][2] < mindist) {
mindist = near4CC[k][2];
minlbl = near;
}
}*/

//Sort wrt col.3....
/*double *temp = new double[3];
for (int k = 3; k >= 0; k--) {
for (int j = 0; j < k; j++) {
if (near4CC[j][2] > near4CC[j + 1][2]) {
temp = near4CC[j];
near4CC[j] = near4CC[j + 1];
near4CC[j + 1] = temp;
}
}
}	*/

////Check for previous connection...
//int count = 0;
//while (count <= 3) {
//	if (near4CC[count][2] == -1) {
//		count++;
//		continue;
//	}
//	int tie = 0, tie_nbr = 0;
//	for (int k = 0; k < labelNum; k++) {
//		if (connections[k][1] == int(near4CC[count][0])) {
//			tie = 1;
//			tie_nbr = k;
//			break;
//		}
//	}
//	if (tie) {
//		if (connections[tie_nbr][2] > near4CC[count][2]) {
//			connections[i][1] = near4CC[count][0];
//			connections[i][2] = near4CC[count][2];
//			//Reset previous connection...
//			connections[tie_nbr][1] = 0;
//			connections[tie_nbr][2] = -1;
//		}
//		else {
//			count++;
//			continue;
//		}
//	}
//	if (!tie) {
//		connections[i][1] = near4CC[count][0];
//		connections[i][2] = near4CC[count][2];
//	}//}	

//Find the centroid of small text-nbd around the end point...
//Point centroidR = Point(0, 0);
//int count = 0;
//for (int k = 0; k <= 9; k++) {
//	for (int l = box[i][4]; l <= box[i][5]; l++) {
//		if (img[l][box[i][3] - k] == 175) {
//			centroidR.x += box[i][3] - k;
//			centroidR.y += l;
//			count++;
//		}
//	}
//}
//centroidR.x = int(centroidR.x / count);
//centroidR.y = int(centroidR.y / count);

//Point centroidL = Point(0, 0);
//count = 0;
//for (int k = 0; k <= 9; k++) {
//	for (int l = box[i][4]; l <= box[i][5]; l++) {
//		if (img[l][box[i][2] + k] == 175) {
//			centroidL.x += box[i][2] + k;
//			centroidL.y += l;
//			count++;
//		}
//	}
//}
//centroidL.x = int(centroidL.x / count);
//centroidL.y = int(centroidL.y / count);

//Find the minimum gradient of centroid and the centroids of the nbrs...
/*double *nbrSlope = new double[nbrno];*/
//std::cout << "\nLabel " << i << " has change in angles = ";
//for (int nb = 0; nb < nbrno; nb++) {
//	Point centreL = Point(0, 0);
//	int count = 0;
//	int nbr = nbrs[nb];
//	for (int k = 0; k <= 9; k++) {
//		for (int l = box[nbr][4]; l <= box[nbr][5]; l++) {
//			if (img[l][box[nbr][2] + k] == 175) {
//				centreL.x += box[nbr][2] + k;
//				centreL.y += l;
//				count++;
//			}
//		}
//	}
//	centreL.x = int(centreL.x / count);
//	centreL.y = int(centreL.y / count);

//	Point centreR = Point(0, 0);
//	count = 0;
//	for (int k = 0; k <= 9; k++) {
//		for (int l = box[nbr][4]; l <= box[nbr][5]; l++) {
//			if (img[l][box[nbr][3] - k] == 175) {
//				centreR.x += box[nbr][3] - k;
//				centreR.y += l;
//				count++;
//			}
//		}
//	}
//	centreR.x = int(centreR.x / count);
//	centreR.y = int(centreR.y / count);

//	/*centroid = Point(Data[i][0], Data[i][1]);
//	centre = Point(Data[nbr][0], Data[nbr][1]);*/
//	Point vec_org, vec, vec_nbr;
//	vec_org = centroidR - centroidL;
//	vec = centreL - centroidR;
//	vec_nbr = centreR - centreL;
//	nbrSlope[nb] = 1 - (abs(vec.dot(vec_org) / (norm(vec_org)*norm(vec))))*(abs(vec.dot(vec_nbr) / (norm(vec_nbr)*norm(vec))));
//	std::cout <<  nbrSlope[nb] << " , ";
//}

