#include "datapoint.h"
#include "matrix.h"

#include <math.h>
#include <algorithm>
#include <vector>

using namespace std;

DataPoint::DataPoint(): vector<double>() {
}

DataPoint::DataPoint(int dim): vector<double>(dim, double(0)) {
}

DataPoint::DataPoint(int dim, double value): vector<double>(dim, value) {
}

DataPoint::DataPoint(const DataPoint & o): vector<double>(o) {
}

DataPoint::DataPoint(const vector<double> & o): vector<double>(o) {
}

DataPoint DataPoint::operator + (const DataPoint& o) {
	int D = this->dim();
	DataPoint ans = DataPoint(D);
	for (int i = 0; i < D; i++) {
		ans[i] = (*this)[i] + o[i];
	}
	return ans;
}

DataPoint DataPoint::operator - (const DataPoint& o) {
	int D = this->dim();
	DataPoint ans = DataPoint(D);
	for (int i = 0; i < D; i++) {
		ans[i] = (*this)[i] - o[i];
	}
	return ans;
}

double DataPoint::operator * (const DataPoint& o) {
	int D = this->dim();
	double ans = double(0);
	for (int i = 0; i < D; i++) {
		ans += (*this)[i] * o[i];
	}
	return ans;
}

DataPoint DataPoint::operator + (const double o) {
	int D = this->dim();
	DataPoint ans = DataPoint(D);
	for (int i = 0; i < D; i++) {
		ans[i] = (*this)[i] + o;
	}
	return ans;
}

DataPoint DataPoint::operator - (const double o) {
	int D = this->dim();
	DataPoint ans = DataPoint(D);
	for (int i = 0; i < D; i++) {
		ans[i] = (*this)[i] - o;
	}
	return ans;
}

DataPoint DataPoint::operator * (const double o) {
	int D = this->dim();
	DataPoint ans = DataPoint(D);
	for (int i = 0; i < D; i++) {
		ans[i] = (*this)[i] * o;
	}
	return ans;
}

DataPoint DataPoint::operator / (const double o) {
	int D = this->dim();
	DataPoint ans = DataPoint(D);
	for (int i = 0; i < D; i++) {
		ans[i] = (*this)[i] / o;
	}
	return ans;
}

DataPoint& DataPoint::operator = (const DataPoint& o) {
	this->clear();
	int D = o.dim();
	for (int i = 0; i < D; i++) {
		this->push_back(o[i]);
	}
	return *this;
}

DataPoint& DataPoint::operator += (const DataPoint& o) {
	return *this = *this + o;
}

DataPoint& DataPoint::operator -= (const DataPoint& o) {
	return *this = *this - o;
}

DataPoint& DataPoint::operator += (const double o) {
	return *this = *this + o;
}

DataPoint& DataPoint::operator -= (const double o) {
	return *this = *this - o;
}

DataPoint& DataPoint::operator *= (const double o) {
	return *this = (*this) * o;
}

DataPoint& DataPoint::operator /= (const double o) {
	return *this = (*this) / o;
}

int DataPoint::dim() const {
	return int(this->size());
}

double DataPoint::form1() const {
	double ans = double(0);
	int D = this->dim();
	for (int i = 0; i < D; i++) {
		ans += (double)fabs(double((*this)[i]));
	}
	return ans;
}

double DataPoint::form2() const {
	int D = this->dim();
	double ans = 0.0;
	for (int i = 0; i < D; i++) {
		ans += (*this)[i] * (*this)[i];
	}
	return (double)sqrt(ans);
}

double DataPoint::forminf() const {
	double ans = double(0);
	int D = this->dim();
	for (int i = 0; i < D; i++) {
		ans = max(ans, (double)fabs(double((*this)[i])));
	}
	return ans;
}

double DataPoint::dis(const DataPoint& o) const {
	int D = this->dim();
	double ans = 0.0;
	for (int i = 0; i < D; i++) {
		ans += ((*this)[i] - o[i]) * ((*this)[i] - o[i]);
	}
	return (double)sqrt(ans);
}

void DataPoint::normalize() {
	double f = form2();
	int n = this->size();
	for (int i = 0; i < n; i++) {
		(*this)[i] /= f;
	}
}

DataPoint DataPoint::get_mean(const vector< DataPoint >& o) {
	int N = (int)o.size();
	if (N == 0) {
		return DataPoint();
	}
	int D = o[0].dim();
	DataPoint ans = DataPoint(D);
	for (int i = 0; i < N; i++) {
		ans += o[i];
	}
	return ans;
}

Matrix DataPoint::get_matrix(const vector< DataPoint >& o) {
	vector< vector<double> > ans;
	vector<double> tmp;
	for (int i = 0; i < (int)o.size(); i++) {
		tmp.clear();
		for (int j = 0; j < o[0].dim(); j++) {
			tmp.push_back(o[i][j]);
		}
		ans.push_back(tmp);
	}
	return Matrix(ans);
}

vector< DataPoint > DataPoint::get_datapoints(const Matrix& o) {
	vector< DataPoint > ans;
	DataPoint tmp;
	for (int i = 0; i < o.row(); i++) {
		tmp.clear();
		for (int j = 0; j < o.column(); j++) {
			tmp.push_back(o[i][j]);
		}
		ans.push_back(tmp);
	}
	return ans;
}
