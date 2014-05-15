#ifndef _SHANGKE_KMEANS_H_
#define _SHANGKE_KMEANS_H_

#include "datapoint.h"
#include <set>
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif

using namespace std;

class Kmeans {
	private:
		const static int MAX_K = 1 << 10;
		const static double INF = 1e200;
		const static double EPS = 1.0e-8;

	public:
		int numSize, kSize, dimSize;
		vector< DataPoint > cluster[MAX_K];
		vector< int > subscript[MAX_K];
		vector< DataPoint > means;
		double minVal, JK;

		int getLable(vector< DataPoint >& means, DataPoint& point);
		double getVal(vector< DataPoint >& means, vector< DataPoint > cluster[]);
		void getMeans(vector< DataPoint >& means, vector< DataPoint > cluster[]);
		void K_Means(vector< DataPoint >& points, int k, int Tim);
};

#ifdef __cplusplus
}
#endif

#endif
