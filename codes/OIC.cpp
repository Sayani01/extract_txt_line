#include <stdio.h>
#include<conio.h>
#include<iostream>
#include <vector>
#include <string>
#include "opencv2/imgproc/imgproc.hpp"
//#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "VarDeclare.h"

using namespace std;
using namespace cv;
namespace cv {
	using std::vector;
};

int ObjectOnBoundary(int **img, int ilimit[2], int jlimit[2]) {
	int flag = 0;
	for (int i = ilimit[0]; i <= ilimit[1]; i++) {
		if (img[i][jlimit[0]] == 0) { //found on left boundary
			flag = 1;
			return(1);
		}
	}
	if (flag == 0) {
		for (int j = jlimit[0]; j <= jlimit[1]; j++) {
			if (img[ilimit[0]][j] == 0) { //found on top boundary
				flag = 1;
				return(1);
			}
		}
	}
	if (flag == 0) {
		for (int i = ilimit[0]; i <= ilimit[1]; i++) {
			if (img[i][jlimit[1]] == 0) { //found on right boundary
				flag = 1;
				return(1);
			}
		}
	}
	if (flag == 0) {
		for (int j = jlimit[0]; j <= jlimit[1]; j++) {
			if (img[ilimit[1]][j] == 0) { //found on bottom boundary
				flag = 1;
				return(1);
			}
		}
	}
	return(0);
}

int vtype(int **img, Point q, int g) {
	int U[5] = { 0 };
	int ilimit[2], jlimit[2];
	int m = 0, rr = 0, t = 0;

	//Find object on boundary of U1...
	ilimit[0] = q.y - g ;
	ilimit[1] = q.y;
	jlimit[0] = q.x;
	jlimit[1] = q.x + g;
	U[1] = ObjectOnBoundary(img, ilimit, jlimit);

	//Find object on boundary of U2...
	ilimit[0] = q.y - g;
	ilimit[1] = q.y;
	jlimit[0] = q.x - g;
	jlimit[1] = q.x;
	U[2] = ObjectOnBoundary(img, ilimit, jlimit);

	//Find object on boundary of U3...
	ilimit[0] = q.y;
	ilimit[1] = q.y + g;
	jlimit[0] = q.x - g;
	jlimit[1] = q.x;
	U[3] = ObjectOnBoundary(img, ilimit, jlimit);
	
	//Find object on boundary of U4...
	ilimit[0] = q.y;
	ilimit[1] = q.y + g;
	jlimit[0] = q.x;
	jlimit[1] = q.x + g;
	U[4] = ObjectOnBoundary(img, ilimit, jlimit);

	//std::cout << "\n U1 = " << U[1] << "\n U2 = " << U[2] << "\n U3 = " << U[3] << "\n U4 = " << U[4];
	int flag = 0;
	for (int k = 1; k <= 4; k++) {
		if (U[k] == 1) {
			flag = 1;
			m++;
			rr += k;
		}
	}
	if (flag == 0) {
		std::cout << "\nHERE!";
		std::cout << "\nSTOP!";
	}
	if (m == 2 && (rr == 4 || rr == 6))
		t = -2;
	else if (m == 0 || m == 4)
		t = 0;
	else
		t = 2 - m;

	std::cout << " ts = " << t;
	return(t);	
}

Point start(int **img, Point q, int g, int &ts) {
	int is = (ceil(q.x / double(g)) - 1)*g;
	int js = (ceil(q.y / double(g)) - 1)*g;
	if (is == q.x) is -= g;
	if (js == q.y) js -= g;
	std::cout << "\n First black = " << js << " , " << is;
	ts = vtype(img, Point(is, js), g);
	return(Point(is, js));
}

void MakeOIC(int **img, vector<Point> &sequence, int g) {
	int ts = 0;
	Point qs = start(img, sequence[0], g, ts);
	int d = (2 + ts) % 4;
	int t = ts;
	Point q = qs;
	sequence.pop_back();
	do {
		std::cout << "\nt = " << t << " , d = " << d;
		if (abs(t) == 1)
			sequence.push_back(q);
		if (d == 0) {
			q = Point(q.x + g, q.y); //right
		}
		else if (d == 1) {
			q = Point(q.x, q.y - g); //up
		}
		else if (d == 2) {
			q = Point(q.x - g, q.y);  //left
		}
		else if (d == 3) {
			q = Point(q.x, q.y + g );  //down
		}
		t = vtype(img, q, g);
		if (t == -2)
			t = -1;
		if (d + t < 0) {
			d = d + 4;
		}
		d = (d + t) % 4;
		std::cout << "\nnext = " << q.y << " , " << q.x;
	} while (q != qs);
}

void MakeOIC(int **img, int **visited, int **T, vector<Point> &sequence, int g, int d, int ts, Point qs) {
	int t = ts;
	qs = Point(qs.x*g, qs.y*g);
	Point q = qs;
	sequence.push_back(q);
	if (t == 1) d = 3;
	if (t == -1) d = 0;
	std::cout << "\nstart = " << q.y << " , " << q.x;
	do {
		std::cout << "\nt = " << t << " , d = " << d;
		if (abs(t) == 1) {
			sequence.push_back(q);
			visited[q.y / g][q.x / g] = 1;
		}
		if (d == 0) {
			q = Point(q.x + g, q.y); //right
		}
		else if (d == 1) {
			q = Point(q.x, q.y - g); //up
		}
		else if (d == 2) {
			q = Point(q.x - g, q.y);  //left
		}
		else if (d == 3) {
			q = Point(q.x, q.y + g);  //down
		}
		t = T[q.y / g][q.x / g];
		if (t == -2)
			t = -1;
		if (d + t < 0) {
			d = d + 4;
		}
		d = (d + t) % 4;
		std::cout << "\nnext = " << q.y << " , " << q.x;
	} while (q != qs);
}


void MakeMultiOIC(Mat colimg, int **img, int R, int C, vector<Point> &sequence, int g) {
	int h = R/g, w = C/g;
	int **visited = new int*[h];
	int **t = new int*[h];
	for (int i = 0; i < h; i++) {
		visited[i] = new int[w];
		t[i] = new int[w];
		for (int j = 0; j < w; j++) {
			visited[i][j] = 0;
			t[i][j] = 0;
		}
	}
	Vec3b color1 = { 255,0,0 };
	Vec3b color2 = { 0,255,0 };
	for (int i = 1; i < h - 1; i++) {
		for (int j = 1; j < w - 1; j++) {
			Point q = Point(j*g, i*g);
			t[i][j] = vtype(img, q, g);
			if(t[i][j] == 1)
				colimg.at<Vec3b>(i*g, j*g) = color1;
			if (t[i][j] == -1)
				colimg.at<Vec3b>(i*g, j*g) = color2;
		}
	}
	imwrite("Data/ICDAR/Data218/MultiOIC.tif", colimg);
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			if (visited[i][j] == 0) {
				Point q = Point(j, i);
				if (abs(t[i][j]) == 1) {
					std::cout << "\nNEW OIC...";
					std::vector<Point> seq;
					MakeOIC(img, visited, t, seq, g, 3 , t[i][j], q);
					for (int k = 0; k<int(seq.size()) - 1; k++) {
						line(colimg, seq[k], seq[k + 1], Scalar(255, 0, 255), 1, 8);
						sequence.push_back(seq[k]);
					}		
					line(colimg, seq[int(seq.size()) - 1], seq[0], Scalar(255, 0, 255), 1, 8);
					sequence.push_back(seq[int(seq.size()) - 1]);
					sequence.push_back(seq[0]);
				}
				imwrite("Data/ICDAR/Data218/MultiOIC.tif", colimg);
			}
		}
	}
}
