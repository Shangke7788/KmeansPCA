#include "datapoint.h"
#include "matrix.h"

#include <math.h>
#include <algorithm>
#include <vector>

using namespace std;

template<typename T> DataPoint<T>::DataPoint(): vector<T>() {
}

template<typename T> DataPoint<T>::DataPoint(int dim): vector<T>(dim, T(0)) {
}

template<typename T> DataPoint<T>::DataPoint(int dim, T value): vector<T>(dim, value) {
}

template<typename T> DataPoint<T>::DataPoint(const DataPoint<T> & o): vector<T>(o) {
}

template<typename T> DataPoint<T>::DataPoint(const vector<T> & o): vector<T>(o) {
}

template<typename T> DataPoint<T> DataPoint<T>::operator + (const DataPoint<T>& o) {
	int D = this->dim();
	DataPoint<T> ans = DataPoint<T>(D);
	for (int i = 0; i < D; i++) {
		ans[i] = (*this)[i] + o[i];
	}
	return ans;
}

template<typename T> DataPoint<T> DataPoint<T>::operator - (const DataPoint<T>& o) {
	int D = this->dim();
	DataPoint<T> ans = DataPoint<T>(D);
	for (int i = 0; i < D; i++) {
		ans[i] = (*this)[i] - o[i];
	}
	return ans;
}

template<typename T> T DataPoint<T>::operator * (const DataPoint<T>& o) {
	int D = this->dim();
	T ans = T(0);
	for (int i = 0; i < D; i++) {
		ans += (*this)[i] * o[i];
	}
	return ans;
}

template<typename T> DataPoint<T> DataPoint<T>::operator + (const T o) {
	int D = this->dim();
	DataPoint<T> ans = DataPoint<T>(D);
	for (int i = 0; i < D; i++) {
		ans[i] = (*this)[i] + o;
	}
	return ans;
}

template<typename T> DataPoint<T> DataPoint<T>::operator - (const T o) {
	int D = this->dim();
	DataPoint<T> ans = DataPoint<T>(D);
	for (int i = 0; i < D; i++) {
		ans[i] = (*this)[i] - o;
	}
	return ans;
}

template<typename T> DataPoint<T> DataPoint<T>::operator * (const T o) {
	int D = this->dim();
	DataPoint<T> ans = DataPoint<T>(D);
	for (int i = 0; i < D; i++) {
		ans[i] = (*this)[i] * o;
	}
	return ans;
}

template<typename T> DataPoint<T> DataPoint<T>::operator / (const T o) {
	int D = this->dim();
	DataPoint<T> ans = DataPoint<T>(D);
	for (int i = 0; i < D; i++) {
		ans[i] = (*this)[i] / o;
	}
	return ans;
}

template<typename T> DataPoint<T>& DataPoint<T>::operator = (const DataPoint<T>& o) {
	this->clear();
	int D = o.dim();
	for (int i = 0; i < D; i++) {
		this->push_back(o[i]);
	}
	return *this;
}

template<typename T> DataPoint<T>& DataPoint<T>::operator += (const DataPoint<T>& o) {
	return *this = *this + o;
}

template<typename T> DataPoint<T>& DataPoint<T>::operator -= (const DataPoint<T>& o) {
	return *this = *this - o;
}

template<typename T> DataPoint<T>& DataPoint<T>::operator += (const T o) {
	return *this = *this + o;
}

template<typename T> DataPoint<T>& DataPoint<T>::operator -= (const T o) {
	return *this = *this - o;
}

template<typename T> DataPoint<T>& DataPoint<T>::operator *= (const T o) {
	return *this = (*this) * o;
}

template<typename T> DataPoint<T>& DataPoint<T>::operator /= (const T o) {
	return *this = (*this) / o;
}

template<typename T> int DataPoint<T>::dim() const {
	return int(this->size());
}

template<typename T> T DataPoint<T>::form1() const {
	T ans = T(0);
	int D = this->dim();
	for (int i = 0; i < D; i++) {
		ans += (T)fabs(double((*this)[i]));
	}
	return ans;
}

template<typename T> T DataPoint<T>::form2() const {
	return (T)sqrt(double((*this) * (*this)));
}

template<typename T> T DataPoint<T>::forminf() const {
	T ans = T(0);
	int D = this->dim();
	for (int i = 0; i < D; i++) {
		ans = max(ans, (T)fabs(double((*this)[i])));
	}
	return ans;
}

template<typename T> T DataPoint<T>::dis(const DataPoint<T>& o) const {
	return (T)sqrt(double((*this) * o));
}

template<typename T> DataPoint<T> DataPoint<T>::get_mean(const vector< DataPoint<T> >& o) {
	int N = o->size();
	if (N == 0) {
		return DataPoint<T>();
	}
	int D = o[0]->dim();
	DataPoint<T> ans = DataPoint<T>(D);
	for (int i = 0; i < N; i++) {
		ans += o[i];
	}
	return ans;
}

template<typename T> Matrix<T> DataPoint<T>::get_matrix(const vector< DataPoint<T> >& o) {
	vector< vector<T> > ans = o;
	return Matrix<T>(ans);
}

template<typename T> vector< DataPoint<T> > DataPoint<T>::get_datapoints(const Matrix<T>& o) {
	vector< vector<T> > ans = o;
	return vector< DataPoint<T> >(o);
}
