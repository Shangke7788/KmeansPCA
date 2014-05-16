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

extern "C" void gao(const char dic[][50], const int num[], int knum, const char filename[]) {
	Kmeans K;
	string F = filename;
	string fname = "result\\" + F + ".txt";
	string picname = "result\\" + F + "\\";
	FILE *ftxt, *fpic;
	ftxt = fopen(fname.c_str(), "w");
	fprintf(stderr, "%s:\n", filename);
	for (int i = 0; i < CNT; i++) {
		fprintf(stderr, "Now selecting and read files... This will take about 1 minute.\n");
		vector< DataPoint > a2 = get_datapoints(dic, num, knum);
		fprintf(stderr, "Read files ends, start K-means... This will take about 5 seconds.\n");
		K.K_Means(a2, knum, 20);
		fprintf(ftxt, "Case #%d:\n", i + 1);
		fprintf(stderr, "Case #%d:\n", i + 1);
		fprintf(ftxt, "Km = %.9lf\n", K.JK);
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
		fprintf(ftxt, "P%d = %.9lf\n", knum, y);
		fprintf(stderr, "P%d = %.9lf\n", knum, y);
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
				if (C[r][c] < 0) {
					P[r][c] = 0;
				} else {
					P[r][c] = C[r][c] / (sqrt(C[r][r]) * sqrt(C[c][c]));
				}
			}
		}
		sprintf(buffers, "%s%d.txt", picname.c_str(), i + 1);
		fpic = fopen(buffers, "w");
		for (int r = 0; r < P.row(); r++) {
			for (int c = 0; c < P.column(); c++) {
				fprintf(fpic, "%.9lf%c", P[r][c], " \n"[c == P.column() - 1]);
			}
		}
		fclose(fpic);
		fprintf(stderr, "Output the matrix with stdout ends, ");
		if (i == CNT - 1) {
			fprintf(stderr, "all the loops ends.\n");
		} else {
			fprintf(stderr, "now turn to the next loop.\n");
		}
	}
	fclose(ftxt);
	fprintf(stderr, "\n");
}

int main() {
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

	gao(A2, BALANCE2, 2, "A2");
	gao(B2, BALANCE2, 2, "B2");
	gao(A5, BALANCE5, 5, "A5_Balance");
	gao(A5, UNBALANCE5, 5, "A5_Unbalance");
	gao(B5, BALANCE5, 5, "B5_Balance");
	gao(B5, UNBALANCE5, 5, "B5_Unbalance");

	return 0;
}
