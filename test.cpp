#include "matrix.h"
#include "datapoint.h"
#include "files.h"
#include "trie.h"

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
	Files fs;
	fs.clear();
	fs.getFiles("..\\20_newsgroups");
	for (int i = 0; i < (int)fs.size(); i++) {
		cout << fs[i] << endl;
	}
	cout << fs.size() << endl;
	return 0;
}
