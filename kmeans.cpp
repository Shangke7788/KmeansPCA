#include "kmeans.h"
#include "datapoint.h"

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <set>
#include <algorithm>
#include <string.h>
#include <math.h>
#include <time.h>

using namespace std;

// 获取point元素最近的一个中心下标
int Kmeans::getLable(vector< DataPoint<double> > &means, DataPoint<double> &point) {
	double minDis = INF;
	int res = -1;
	for (int i = 0; i < kSize; i++) {
		double tmp = means[i].dis(point);
		if (minDis > tmp) {
			minDis = tmp, res = i;
		}
	}
	return res;
}

// 得到当前的标准差之和E
double Kmeans::getVal(vector< DataPoint<double> > &means, vector< DataPoint<double> > cluster[]) {
	double res = 0.0;
	for (int i = 0; i < kSize; i++) {
		for (int j = 0; j < (int)cluster[i].size(); j++) {
			res += means[i].dis(cluster[i][j]);
		}
	}
	return res;
}

// 各个集合得到新的质心，迭代的过程
void Kmeans::getMeans(vector< DataPoint<double> >& means, vector< DataPoint<double> > cluster[]) {
	means.clear();
	for (int i = 0; i < kSize; i++) {
		DataPoint<double> tmp = DataPoint<double>();
		for (int j = 0; j < (int)cluster[i].size(); j++) {
			tmp = tmp + cluster[i][j];
		}
		means.push_back(tmp / (double)cluster[i].size());
	}
}

// 执行K-means算法核心
void Kmeans::K_Means(vector< DataPoint<double> >& points, int k, int Tim = 20) {
	int srandnum = 0;
	JK = 0.0, minVal = INF, srand(srandnum = time(NULL));
	kSize = k, numSize = points.size();
	if (numSize) {
		dimSize = points[0].size();
	} else {
		return;
	}
	set< int > vis;
	vector< DataPoint<double> > cluster[MAX_K];
	vector< int > subscript[MAX_K];
	vector< DataPoint<double> > means;
	while (Tim--) {
		vis.clear(), means.clear();
		for (int i = 0; i < kSize; i++) {
			// 随机选取起始点
			int choose = (long long)rand() * numSize / RAND_MAX;
			if (vis.find(choose) != vis.end()) {
				i--;
				continue;
			}
			vis.insert(choose);
			means.push_back(points[choose]);
		}
		fill(cluster, cluster + kSize, vector< DataPoint<double> >());
		for (int i = 0; i < numSize; i++) {
			int lable = getLable(means, points[i]);
			cluster[lable].push_back(points[i]);
			subscript[lable].push_back(i);
		}

		double oldVal = INF, newVal = getVal(means, cluster);

		int tim = 0;
		// 迭代加深，直到前后两次相差在EPS以内
		while (fabs(oldVal - newVal) > EPS) {
			getMeans(means, cluster);
			fill(cluster, cluster + kSize, vector< DataPoint<double> >());
			for (int i = 0; i < numSize; i++) {
				int lable = getLable(means, points[i]);
				cluster[lable].push_back(points[i]);
				subscript[lable].push_back(i);
			}
			oldVal = newVal, newVal = getVal(means, cluster);
		}

		if (newVal < minVal) {
			minVal = newVal;
			for (int i = 0; i < kSize; i++) {
				this->cluster[i] = cluster[i];
				this->subscript[i] = subscript[i];
			}
			this->means = means;
		}
	}

	JK = 0.0;
	for (int i = 0; i < numSize; i++) {
		JK += points[i] * points[i];
	}
	for (int i = 0; i < kSize; i++) {
		int clusternum = this->cluster[i].size();
		double tmp = 0.0;
		for (int j = 0; j < clusternum; j++) {
			for (int k = 0; k < clusternum; k++) {
				tmp += points[this->subscript[i][j]] * points[this->subscript[i][k]];
			}
		}
		JK -= tmp / clusternum;
	}
}
