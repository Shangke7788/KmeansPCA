#ifndef _SHANGKE_KMEANS_H_
#define _SHANGKE_KMEANS_H_

#include "datapoint.h"
#include <set>
#include <vector>

using namespace std;

class Kmeans {
	private:
		const static int MAX_K = 1 << 10;
		const static double INF = 1e200;
		const static double EPS = 1.0e-8;

	public:
		int numSize, kSize, dimSize;
		vector< DataPoint<double> > cluster[MAX_K];
		vector< int > subscript[MAX_K];
		vector< DataPoint<double> > means;
		double minVal, JK;

		int getLable(vector< DataPoint<double> >& means, DataPoint<double>& point);
		double getVal(vector< DataPoint<double> >& means, vector< DataPoint<double> > cluster[]);
		void getMeans(vector< DataPoint<double> >& means, vector< DataPoint<double> > cluster[]);
		void K_Means(vector< DataPoint<double> >& points, int k, int Tim);
};

#endif
