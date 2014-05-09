#include "matrix.h"
#include "datapoint.h"

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

int main() {
	return 0;
}
