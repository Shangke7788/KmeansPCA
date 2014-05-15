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

	Kmeans K;
	int knum = 2;
	for (int i = 0; i < 10; i++) {
		fprintf(stderr, "Now selecting and read files... This will take about 1 minute.\n");
		vector< DataPoint > a2 = get_datapoints(A2, BALANCE2, 2);
		fprintf(stderr, "Read files ends, start K-means... This will take about 5 seconds.\n");
		K.K_Means(a2, knum, 20);
		printf("Case #%d:\n", i + 1);
		printf("JK = %.9lf\n", K.JK);
		fprintf(stderr, "K-means ends, start svd...\n");
	//	Matrix m2 = DataPoint::get_matrix(a2);
	//	m2 = m2.cov();
	//	Matrix u, s, v;
	//	m2.svd_jacobi(u, s, v);
	}

	return 0;
}
