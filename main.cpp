#include "kmeans.h"
#include "datapoint.h"
#include "matrix.h"
#include "tfidf.h"
#include "trie.h"

#include "files.h"
#include "constant.h"

#include <stdio.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <map>
#include <set>

using namespace std;

void out(const Matrix& mat) {
	fprintf(stderr, "mat:\n");
	int R = mat.row(), C = mat.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			fprintf(stderr, "%.9lf%c", mat[i][j], " \n"[j == C - 1]);
		}
	}
}

map<string, double> words;
char buffers[1 << 10];

extern "C" vector< DataPoint > get_datapoints(const char dic[][50], const int num[], int n) {
	set<int> vis;
	DataPoint point;
	vector< DataPoint > points;
	for (int i = 0; i < n; i++) {
		Files fs;
		fs.getFiles(dic[i]);
		vis.clear();
		for (int j = 0; j < num[i]; j++) {
			int x = rand() % (int)fs.size();
			if (vis.find(x) != vis.end()) {
				j--;
				continue;
			}
			vis.insert(x);
			point = DataPoint();
			int k = 0;
			for (map<string, double>::iterator iter = words.begin(); iter != words.end(); iter++) {
				strcpy(buffers, (iter->first).c_str());
				point.push_back(iter->second * document_word_num(fs[x].c_str(), buffers));
				k++;
			}
			point.normalize();
			points.push_back(point);
		}
	}
	return points;
}

extern "C" void gao(const char dic[][50], const int num[], int knum) {
	Kmeans K;
	//for (int i = 0; i < CNT; i++) {
	for (int i = 0; i < 1; i++) {
		fprintf(stderr, "Now selecting and read files... This will take about 1 minute.\n");
		vector< DataPoint > a2 = get_datapoints(dic, num, knum);
		fprintf(stderr, "Read files ends, start K-means... This will take about 5 seconds.\n");
		K.K_Means(a2, knum, 20);
		printf("Case #%d:\n", i + 1);
		fprintf(stderr, "Case #%d:\n", i + 1);
		printf("Km = %.9lf\n", K.JK);
		fprintf(stderr, "Km = %.9lf\n", K.JK);
		fprintf(stderr, "K-means ends, start svd... This will take about 1 second or 2 minutes.\n");
		Matrix m2 = DataPoint::get_matrix(a2);
		m2 = m2.column_center();
		a2 = DataPoint::get_datapoints(m2);
		double y = 0.0;
		for (int j = 0; j < (int)a2.size(); j++) {
			y += a2[j] * a2[j];
		}
		m2 = m2.cov2();
		Matrix u, s, v;
		m2.svn(s, u, v);
		fprintf(stderr, "ny^2 = %.9lf\n", y);
		for (int j = 0; j < knum - 1; j++) {
			y -= s[j][j];
		}
		printf("P2 = %.9lf\n", y);
		fprintf(stderr, "P2 = %.9lf\n", y);
		fprintf(stderr, "Now svd ends, start output the matrix with stdout.\n");
		Matrix q = Matrix((int)a2.size(), knum);
		for (int r = 0; r < q.row(); r++) {
			for (int c = 0; c < knum - 1; c++) {
				q[r][c] = v[r][c];
			}
			q[r][knum - 1] = 1.0 / sqrt(double(a2.size()));
		}
		Matrix C = q * q.transpose();
		Matrix P = Matrix(C.row(), C.column());
		for (int r = 0; r < C.row(); r++) {
			for (int c = 0; c < C.column(); c++) {
				if (r == c) {
					P[r][c] = 1;
				} else if (C[r][c] < 0.0) {
					P[r][c] = 0;
				} else if (C[r][c] / (sqrt(C[r][r]) * sqrt(C[c][c])) >= 0.5) {
					P[r][c] = 1;
				} else {
					P[r][c] = 0;
				}
			}
		}
		for (int r = 0; r < P.row(); r++) {
			for (int c = 0; c < P.column(); c++) {
				if (P[r][c] <= 0.1) {
					printf(" ");
				} else {
					printf("X");
				}
			}
			printf("\n");
		}
		fprintf(stderr, "Output the matrix with stdout ends, ");
		if (i == 1 - 1) {
			fprintf(stderr, "all the loops ends.\n");
		} else {
			fprintf(stderr, "now turn to the next loop.\n");
		}
	}
}

int main() {
	/*****************test start*******************
	double test[4][4] = {
		0.685694,        -0.0392685,      1.27368,     0.516904,        
		-0.0392685,      0.188004,        -0.321713,       -0.117981,       
		1.27368,     -0.321713,       3.11318,     1.29639,     
		0.516904,        -0.117981,       1.29639,     0.582414        
	};
	Matrix aa = Matrix(4, 4);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			aa[i][j] = test[i][j];
		}
	}
	Matrix U, S, V;
	aa.svn(S, U, V);
	out(aa), out(S), out(U), out(V);
	return 0;
	*******************test end******************/
	srand(time(NULL));
	FILE * _1000words;
	_1000words = fopen("1000words.txt", "r");
	if (_1000words == NULL) {
		get_1000_words();
		_1000words = fopen("1000words.txt", "r");
	}
	words.clear();
	double idf_val;
	int idf_num;
	while (fscanf(_1000words, "%s%lf%d", buffers, &idf_val, &idf_num) != EOF) {
		string tmp = buffers;
		words[tmp] = idf_val;
	}
	fclose(_1000words);

	printf("A2:\n"); gao(A2, BALANCE2, 2); printf("\n");
	printf("B2:\n"); gao(B2, BALANCE2, 2); printf("\n");
	printf("A5 Balance:\n"); gao(A5, BALANCE5, 5); printf("\n");
	printf("A5 Unbalance\n"); gao(A5, UNBALANCE5, 5); printf("\n");
	printf("B5 Balance:\n"); gao(B5, BALANCE5, 5); printf("\n");
	printf("B5 Unbalance\n"); gao(B5, UNBALANCE5, 5); printf("\n");

	return 0;
}
