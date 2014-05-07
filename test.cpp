#include "matrix.h"

#include <stdio.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <map>

using namespace std;

template<typename T> void out(const Matrix<T>& mat) {
	int R = mat.row(), C = mat.column();
	fprintf(stderr, "Matrix\n");
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			fprintf(stderr, "%-12.8lf%c", mat[i][j], "\t\n"[j == C - 1]);
		}
	}
}

double a[4][4] = {
	0.685694, -0.0392685, 1.27368, 0.516904,
	-0.0392685, 0.188004, -0.321713, -0.117981,
	1.27368, -0.321713, 3.11318, 1.29639,
	0.516904, -0.117981, 1.29639, 0.582414
};

int main() {
	time_t t1 = clock();
	Matrix<double> mat = Matrix<double>(4, 4);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			mat[i][j] = a[i][j];
		}
	}
	Matrix<double> U, B, V;
	mat.svd_jacobi(U, B, V);
	out(mat), printf("\n");
	out(U), printf("\n");
	out(B), printf("\n");
	out(V), printf("\n");
	cerr << clock() - t1 << endl;
	return 0;
}
