#include "matrix.h"

#include <vector>
#include <algorithm>
#include <math.h>
#include <map>

using namespace std;

template<typename T> bool Matrix<T>::isZero(const T x) {
	double y = (double)x;
	return ((y > EPS) - (y < -EPS)) == 0;
}

template<typename T> T Matrix<T>::sign(const T x) {
	double y = (double)x;
	return T((y > EPS) - (y < -EPS));
}

template<typename T> T Matrix<T>::sign_without0(const T x) {
	if (x < T(0)) {
		return T(-1);
	} else {
		return T(1);
	}
}

template<typename T> Matrix<T>::Matrix(): vector< vector<T> >() {
}

template<typename T> Matrix<T>::Matrix(int row, int column): vector< vector<T> >(row, vector<T>(column, T(0))) {
}

template<typename T> Matrix<T>::Matrix(int row): vector< vector<T> >(row, vector<T>(row, T(0))) {
	for (int i = 0; i < row; i++) {
		(*this)[i][i] = T(1);
	}
}

template<typename T> Matrix<T>::Matrix(int row, int column, T value): vector< vector<T> >(row, vector<T>(column, value)) {
}

template<typename T> Matrix<T>::Matrix(const vector< vector<T> > &mat): vector< vector<T> >(mat) {
}

template<typename T> Matrix<T>::Matrix(const Matrix<T> &mat): vector< vector<T> >(mat) {
}

template<typename T> Matrix<T> Matrix<T>::operator + (const Matrix<T> & o) {
	Matrix<T> ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] + o[i][j];
		}
	}
	return ans;
}

template<typename T> Matrix<T> Matrix<T>::operator - (const Matrix<T> & o) {
	Matrix<T> ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] - o[i][j];
		}
	}
	return ans;
}

template<typename T> Matrix<T> Matrix<T>::operator * (const Matrix<T> & o) {
	Matrix<T> ans = Matrix(this->row(), o.column());
	int R = ans.row(), M = this->column(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			for (int k = 0; k < M; k++) {
				ans[i][j] += (*this)[i][k] * o[k][j];
			}
		}
	}
	return ans;
}

template<typename T> Matrix<T> Matrix<T>::operator + (const T o) {
	Matrix<T> ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] + o;
		}
	}
	return ans;
}

template<typename T> Matrix<T> Matrix<T>::operator - (const T o) {
	Matrix<T> ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] - o;
		}
	}
	return ans;
}

template<typename T> Matrix<T> Matrix<T>::operator * (const T o) {
	Matrix<T> ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] * o;
		}
	}
	return ans;
}

template<typename T> Matrix<T> Matrix<T>::operator / (const T o) {
	Matrix<T> ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] / o;
		}
	}
	return ans;
}

template<typename T> Matrix<T>& Matrix<T>::operator = (const Matrix<T> & o) {
	if (this == &o) {
		return *this;
	}
	int R = this->row(), C = this->column();
	for (int i = 0; i < R; i++) {
		(*this)[i].clear();
	}
	this->clear();
	for (int i = 0; i < o.row(); i++) {
		this->push_back(vector<T>(o[i]));
	}
	return *this;
}

template<typename T> Matrix<T>& Matrix<T>::operator += (const Matrix<T> & o) {
	Matrix<T> ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] + o[i][j];
		}
	}
	return (*this) = ans;
}

template<typename T> Matrix<T>& Matrix<T>::operator -= (const Matrix<T> & o) {
	Matrix<T> ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] - o[i][j];
		}
	}
	return (*this) = ans;
}

template<typename T> Matrix<T>& Matrix<T>::operator *= (const Matrix<T> & o) {
	Matrix<T> ans = Matrix(this->row(), o.column());
	int R = ans.row(), M = this->column(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			for (int k = 0; k < M; k++) {
				ans[i][j] += (*this)[i][k] * o[k][j];
			}
		}
	}
	return (*this) = ans;
}

template<typename T> Matrix<T>& Matrix<T>::operator += (const T o) {
	Matrix<T> ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] + o;
		}
	}
	return (*this) = ans;
}

template<typename T> Matrix<T>& Matrix<T>::operator -= (const T o) {
	Matrix<T> ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] - o;
		}
	}
	return (*this) = ans;
}

template<typename T> Matrix<T>& Matrix<T>::operator *= (const T o) {
	Matrix<T> ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] * o;
		}
	}
	return (*this) = ans;
}

template<typename T> Matrix<T>& Matrix<T>::operator /= (const T o) {
	Matrix<T> ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] / o;
		}
	}
	return (*this) = ans;
}

template<typename T> int Matrix<T>::row() const {
	return (int)(this->size());
}

template<typename T> int Matrix<T>::column() const {
	if (this->row() == 0) {
		return 0;
	} else {
		return (int)(this->begin()->size());
	}
}

template<typename T> T Matrix<T>::get_value() const {
	int R = this->row(), C = this->column();
	if (R != C || R == 0) {
		return T(0);
	}
	Matrix<T> mat = *this;
	T ans = T(1);
	int n = R;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (!isZero(mat[j][i])) {
				T t = mat[i][i] / mat[j][i];
				for (int k = i; k < n; k++) {
					mat[i][k] -= t * mat[j][k];
				}
				for (int k = i; k < n; k++) {
					swap(mat[i][k], mat[j][k]);
				}
				ans = -ans;
			}
		}
		if (isZero(mat[i][i])) {
			return T(0);
		} else {
			ans = ans * mat[i][i];
		}
	}
	return ans;
}

template<typename T> int Matrix<T>::rank_det() const {
	Matrix mat = *this;
	int d = 0, R = mat.row(), C = mat.column();
	int ni = 0;
	for (int i = 0; i < R; i++) {
		for (int j = i + 1; j < R; j++) {
			if (!isZero(mat[j][ni])) {
				T t = mat[i][ni] / mat[j][ni];
				for (int k = ni; k < C; k++) {
					mat[i][k] -= t * mat[j][k];
				}
				for (int k = ni; k < C; k++) {
					swap(mat[i][k], mat[j][k]);
				}
			}
		}
		if (ni < R && isZero(mat[i][ni])) {
			ni++, i--;
			continue;
		}
		if (ni >= R) break;
		d++, ni++;
		if (ni >= R) break;
	}
	return d;
}

template<typename T> vector<T> Matrix<T>::get_row(int k) const {
	if (k >= this->row()) {
		return vector<T>();
	}
	vector<T> ans = (*this)[k];
	return ans;
}

template<typename T> vector<T> Matrix<T>::get_column(int k) const {
	if (k >= this->column()) {
		return vector<T>();
	}
	int R = this->row();
	vector<T> ans = vector<T>(R, T(0));
	for (int i = 0; i < R; i++) {
		ans[i] = (*this)[i][k];
	}
	return ans;
}

template<typename T> void Matrix<T>::set_row(const vector<T>& o, int k) {
	if (k >= this->row()) {
		return;
	}
	int C = this->column();
	for (int i = 0; i < C; i++) {
		(*this)[k][i] = o[i];
	}
}

template<typename T> void Matrix<T>::set_column(const vector<T>& o, int k) {
	if (k >= this->column()) {
		return;
	}
	int R = this->row();
	for (int i = 0; i < R; i++) {
		(*this)[i][k] = o[i];
	}
}

template<typename T> Matrix<T> Matrix<T>::transpose() const {
	Matrix<T> ans = Matrix(this->column(), this->row());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[j][i];
		}
	}
	return ans;
}

template<typename T> Matrix<T> Matrix<T>::column_center() const {
	Matrix<T> ans = Matrix<T>(this->row(), this->column());
	vector<T> ave = vector<T>(this->column(), T(0));
	int R = this->row(), C = this->column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ave[j] += (*this)[i][j];
		}
	}
	for (int i = 0; i < C; i++) {
		ave[i] /= R;
	}
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] - ave[j];
		}
	}
	return ans;
}

template<typename T> Matrix<T> Matrix<T>::row_center() const {
	Matrix<T> ans = Matrix<T>(this->row(), this->column());
	vector<T> ave = vector<T>(this->column(), T(0));
	int R = this->row(), C = this->column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ave[i] += (*this)[i][j];
		}
	}
	for (int i = 0; i < R; i++) {
		ave[i] /= C;
	}
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] - ave[i];
		}
	}
	return ans;
}

template<typename T> Matrix<T> Matrix<T>::center() const {
	Matrix<T> ans = Matrix<T>(this->row(), this->column());
	vector<T> r_ave = vector<T>(this->column(), T(0));
	vector<T> c_ave = vector<T>(this->row(), T(0));
	T ave = T(0);
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			r_ave[j] += (*this)[i][j];
			c_ave[i] += (*this)[i][j];
			ave += (*this)[i][j];
		}
	}
	for (int i = 0; i < C; i++) {
		c_ave[i] /= R;
	}
	for (int i = 0; i < R; i++) {
		r_ave[i] /= C;
	}
	ave /= R * C;
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] - r_ave[j] - c_ave[i] + ave;
		}
	}
	return ans;
}

template<typename T> Matrix<T> Matrix<T>::cov() const {
	Matrix<T> cen = this->column_center();
	Matrix<T> ans = Matrix(this->column(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < C; i++) {
		for (int j = 0; j < C; j++) {
			for (int k = 0; k < R; k++) {
				ans[i][j] += cen[i][k] * cen[j][k];
			}
			// ans[i][j] /= R - 1;
		}
	}
	return ans;
}

template<typename T> Matrix<T> Matrix<T>::inverse() const {
	int R = this->row(), C = this->column();
	if (R != C || R == 0) {
		return Matrix<T>();
	}
	Matrix<T> mat = Matrix<T>(R, C + C);
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			mat[i][j] = (*this)[i][j];
			mat[i][j + C] = 0;
		}
		mat[i][i + C] = 1;
	}
	for (int i = 0; i < R; i++) {
		int s = i;
		for (int j = i + 1; j < C; j++) {
			if (fabs((double)mat[j][i]) > fabs((double)mat[s][i])) {
				s = j;
			}
		}
		if (isZero(mat[s][i])) {
			return Matrix<T>();
		}
		T sav = mat[s][i];
		for (int j = i; j < (C << 1); j++) {
			swap(mat[i][j], mat[s][j]);
			mat[i][j] /= sav;
		}
		for (int j = i + 1; j < R; j++) {
			sav = -mat[j][i];
			for (int k = i; k < (C << 1); k++) {
				mat[j][k] += sav * mat[i][k];
			}
		}
	}
	for (int i = R - 1; i >= 0; i--) {
		for (int j = i - 1; j >= 0; j--) {
			T sav = -mat[j][i];
			for (int k = i; k < (C << 1); k++) {
				mat[j][k] += sav * mat[i][k];
			}
		}
	}
	Matrix<T> ans = Matrix<T>(R, C);
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = mat[i][j + C];
		}
	}
	return ans;
}

template<typename T> void Matrix<T>::svd_qr(Matrix<T>& U, Matrix<T>& B, Matrix<T>& V) const {
	Matrix A = *this;
	bool hastran = false;
	int m = this->row(), n = this->column();
	if (m < n) {
		A = A.transpose();
		hastran = true;
		swap(m, n);
	}
	// 转化为 主对角线+副对角线 矩阵
	vector<T> delta = vector<T>(n, T(0));
	vector<T> gama = vector<T>(n, T(0));
	U = Matrix<T>(m), V = Matrix<T>(n);
	for (int k = 0; k < n; k++) {
		T cc = T(0);
		for (int i = 0; i < m - k; i++) {
			cc += A[i][0] * A[i][0];
		}
		cc = (T)sqrt(double(cc));
		delta[k] = cc;
		vector<T> v1 = A.get_column(0);
		vector<T> e1 = vector<T>(m - k, T(0));
		e1[0] = cc;
		Matrix<T> P1 = householder(v1, e1);
		Matrix<T> tmp = Matrix<T>(m - k, n - k - 1);
		int R = tmp.row(), C = tmp.column();
		for (int i = 0; i < R; i++) {
			for (int j = 0; j < C; j++) {
				tmp[i][j] = A[i][j + 1];
			}
		}
		A = P1 * tmp;
		tmp = Matrix<T>(m);
		R = tmp.row(), C = tmp.column();
		for (int i = k; i < R; i++) {
			for (int j = k; j < C; j++) {
				tmp[i][j] = P1[i - k][j - k];
			}
		}
		U *= tmp;

		if (k == n - 1) {
			break;
		}

		cc = T(0);
		for (int i = 0; i < n - k - 1; i++) {
			cc += A[0][i] * A[0][i];
		}
		cc = (T)sqrt(double(cc));
		gama[k + 1] = cc;
		vector<T> u1 = A.get_row(0);
		e1 = vector<T>(n - k - 1, T(0));
		e1[0] = cc;
		Matrix<T> H1 = householder(u1, e1);
		tmp = Matrix<T>(m - k - 1, n - k - 1);
		R = tmp.row(), C = tmp.column();
		for (int i = 0; i < R; i++) {
			for (int j = 0; j < C; j++) {
				tmp[i][j] = A[i + 1][j];
			}
		}
		A = tmp * H1;

		tmp = Matrix<T>(n);
		R = tmp.row(), C = tmp.column();
		for (int i = k + 1; i < R; i++) {
			for (int j = k + 1; j < C; j++) {
				tmp[i][j] = H1[i - k - 1][j - k - 1];
			}
		}
		V *= tmp;
	}

	bool flag = true;
	while (flag) {
		flag = false;
		gama[0] = T(0);
		int p = n - 1, q = 1;
		for (int i = 1; i < n; i++) {
			if (fabs(double(gama[i])) <= EPS * (fabs(double(delta[i - 1])) + fabs(double(delta[i])))) {
				gama[i] = T(0);
			} else {
				flag = true;
				p = min(p, i), q = max(q, i);
			}
		}
		if (!flag) {
			break;
		}
		p--;
		T max_val = T(0);
		for (int i = 0; i < n; i++) {
			max_val = max((T)fabs(double(delta[i])), max_val);
			max_val = max((T)fabs(double(gama[i])), max_val);
		}
		bool has3 = false;
		int l = 0, L = 0;
		T x, y;
		for (int i = p; i < q; i++) {
			if ((T)fabs(double(delta[i])) <= max_val * (T)EPS) {
				delta[i] = 0, x = gama[i + 1], y = delta[i + 1];
				gama[i + 1] = 0, l = 1, L = i;
				has3 = true;
				break;
			}
		}
		if (has3) {
			while (true) {
				T c, s, r;
				givens(y, x, c, s, r);
				delta[L + l] = r;
				Matrix<T> G = Matrix<T>(m);
				G[L][L] = c, G[L][L + l] = -s;
				G[L + l][L] = s, G[L + l][L + l] = c;
				U = U * G;
				if (l < q - L) {
					x = s * gama[L + l + 1];
					gama[L + l + 1] = c * gama[L + l + 1];
					y = delta[L + l + 1];
					l = l + 1;
				} else {
					break;
				}
			}
		} else {
			Matrix<T> P, Q;
			vector<T> d = vector<T>(q - p + 1, T(0));
			vector<T> g = vector<T>(q - p + 1, T(0));
			for (int i = p; i <= q; i++) {
				d[i - p] = delta[i];
				if (i != p) {
					g[i - p] = gama[i];
				}
			}
			qr(d, g, P, Q);
			for (int i = p; i <= q; i++) {
				delta[i] = d[i - p];
				if (i != p) {
					gama[i] = g[i - p];
				}
			}
			Matrix<T> tmp = Matrix<T>(m);
			for (int i = p; i <= q; i++) {
				for (int j = p; j <= q; j++) {
					tmp[i][j] = P[i - p][j - p];
				}
			}
			U *= tmp;
			tmp = Matrix<T>(n);
			for (int i = p; i <= q; i++) {
				for (int j = q; j <= q; j++) {
					tmp[i][j] = Q[i - p][j - p];
				}
			}
			V *= tmp;
		}
	}

	B = Matrix<T>(m, n);
	for (int i = 0; i < n; i++) {
		B[i][i] = delta[i];
		if (i == n - 1) {
			break;
		}
		B[i][i + 1] = gama[i + 1];
	}
	if (hastran) {
		swap(U, V), U = U.transpose(), V = V.transpose(), B = B.transpose();
	}
}

template<typename T> void Matrix<T>::svd_jacobi(Matrix<T>& U, Matrix<T>& S, Matrix<T>& V) const {
	Matrix A = *this;
	bool hastran = false;
	int R = this->row(), C = this->column();
	if (R > C) {
		A = A.transpose();
		hastran = true;
		swap(R, C);
	}
	U = Matrix<T>(R, R), V = Matrix<T>(C, C);
	vector<T> E = vector<T>(C, 0);

	Matrix<T> B = A.transpose() * A;
	Matrix<T> J = Matrix<T>(C);
	vector<T> ss = vector<T>(C, 0);
	jacobi(B, C, ss, J);
	for (int i = 0; i < (int)ss.size(); i++) {
		ss[i] = (T)sqrt(double(ss[i]));
	}

	multimap<T, int> eigen;
	for (int i = 0; i < (int)ss.size(); i++) {
		eigen.insert(make_pair(ss[i], i));
	}
	typename multimap<T, int>::iterator iter = --eigen.end();
	int num_eig = 0;
	for (int i = 0; i < C; i++, iter--) {
		int index = iter->second;
		E[i] = ss[index];
		if (E[i] > EPS) {
			num_eig++;
		}
		for (int j = 0; j < C; j++) {
			V[j][i] = J[j][index];
		}
	}
	for (int i = 0; i < num_eig; i++) {
		Matrix<T> vi = vector< vector<T> >(1, V.get_column(i));
		T sigma = E[i];
		Matrix<T> ui = A * vi.transpose();
		for (int j = 0; j < R; j++) {
			U[j][i] = ui[j][0] / sigma;
		}
	}

	//U矩阵的后(rows-none_zero)列就不计算了，采用默认值0。因为这后几列对应的奇异值为0,在做数据压缩时用不到这几列

	S = Matrix<T>(R, C);
	for (int i = 0; i < R; i++) {
		S[i][i] = E[i];
	}
	if (hastran) {
		S = S.transpose(), swap(U, V), U = U.transpose(), V = V.transpose();
	}
}

template<typename T> Matrix<T> Matrix<T>::householder(const vector<T> &x, const vector<T> &y) {
	int n = x.size();
	if (n == 0) {
		return Matrix<T>();
	}
	T C = T(0);
	for (int i = 0; i < n; i++) {
		C += (x[i] - y[i]) * (x[i] - y[i]);
	}
	C = (T)sqrt(double(C));
	Matrix<T> H = Matrix<T>(n);
	if (isZero(C)) {
		return H;
	}
	Matrix<T> B = Matrix<T>(n, 1);
	for (int i = 0; i < n; i++) {
		B[i][0] = (x[i] - y[i]) / C;
	}
	Matrix<T> BT = B.transpose();
	H -= B * BT * T(2);
	return H;
}

template<typename T> void Matrix<T>::givens(const T x, const T y, T& c, T& s, T& r) {
	if (isZero(y)) {
		c = T(1), s = T(0), r = T(0);
	} else if (fabs(double(y)) > fabs(double(x))) {
		T t = -x / y;
		s = (T)sqrt(double(T(1) + t * t));
		r = -y * s;
		s = T(1) / s;
		c = s * t;
	} else {
		T t = -y / x;
		c = (T)sqrt(double(T(1) + t * t));
		r = x * c;
		c = T(1) / c;
		s = c * t;
	}
}

template<typename T> void Matrix<T>::update(const T c, const T s, vector<T>& v1, vector<T> &v2) {
	int n = v1.size();
	for (int i = 0; i < n; i++) {
		T t = v1[i];
		v1[i] = c * t - s * v2[i];
		v2[i] = s * t + c * v2[i];
	}
}

template<typename T> void Matrix<T>::qr(vector<T>& delta, vector<T>& gama, Matrix<T>& P, Matrix<T>& Q) {
	int n = delta.size();
	T d = ((delta[n - 2] * delta[n - 2] + gama[n - 2] * gama[n - 2]) - 
			delta[n - 1] * delta[n - 1] + gama[n - 1] * gama[n - 1]) / T(2);
	T mu = (delta[n - 1] * delta[n - 1] + gama[n - 1] * gama[n - 1]) + d - sign(d)
		* (T)sqrt(double(d * d + delta[n - 2] * delta[n - 2] *gama[n - 1] * gama[n - 1]));
	T x = delta[0] * delta[0] - mu, y = delta[0] * gama[1], k = 0;
	Q = Matrix<T>(n), P = Matrix<T>(n);
	while (true) {
		T c, s, r;
		givens(x, y, c, s, r);
		Matrix<T> cs = Matrix<T>(2, 2);
		cs[0][0] = c, cs[0][1] = s;
		cs[1][0] = -s, cs[1][1] = c;
		Matrix<T> tmp = Matrix<T>(2, 2);
		tmp[0][0] = delta[k], tmp[0][1] = gama[k + 1];
		tmp[1][0] = 0, tmp[1][1] = delta[k + 1];
		tmp *= cs;
		x = tmp[0][0], gama[k + 1] = tmp[0][1];
		y = tmp[1][0], delta[k + 1] = tmp[1][1];

		vector<T> k0 = Q.get_column(k);
		vector<T> k1 = Q.get_column(k + 1);
		update(c, s, k0, k1);
		Q.set_column(k0, k);
		Q.set_column(k1, k + 1);

		if (k > 0) {
			gama[k] = r;
		}
		givens(x, y, c, s, r);
		delta[k] = r;

		k0 = P.get_column(k);
		k1 = P.get_column(k + 1);
		update(c, s, k0, k1);
		P.set_column(k0, k);
		P.set_column(k1, k + 1);

		if (k < n - 2) {
			tmp = Matrix<T>(2, 2);
			tmp[0][0] = gama[k + 1], tmp[0][1] = T(0);
			tmp[1][0] = delta[k + 1], tmp[1][1] = gama[k + 2];
			cs = Matrix<T>(2, 2);
			cs[0][0] = c, cs[0][1] = s;
			cs[1][0] = -s, cs[1][1] = c;
			tmp = cs.transpose() * tmp;
			x = tmp[0][0], y = tmp[0][1];
			delta[k + 1] = tmp[1][0], gama[k + 2] = tmp[1][1];
			k = k + 1;
		} else {
			tmp = Matrix<T>(2, 1);
			tmp[0][0] = gama[n - 1], tmp[1][0] = delta[n - 1];
			cs = Matrix<T>(2, 2);
			cs[0][0] = c, cs[0][1] = s;
			cs[1][0] = -s, cs[1][1] = c;
			tmp = cs.transpose() * tmp;
			break;
		}
	}
}

template<typename T> void Matrix<T>::jacobi(Matrix<T>& mat, int size, vector<T>& E, Matrix<T>& J) {
	int iteration = ITERATION;
	while (iteration--) {
		bool pass = true;
		for (int i = 0; i < size; i++) {
			for (int j = i + 1; j < size; j++) {
				rotate(mat, i, j, pass, J);
			}
		}
		if (pass) break;
	}
	for (int i = 0; i < size; i++) {
		E[i] = mat[i][i];
		if (E[i] < EPS) {
			E[i] = 0.0;
		}
	}
}

template<typename T> void Matrix<T>::rotate(Matrix<T>& mat, int i, int j, bool& pass, Matrix<T>& J) {
	T ele = mat[i][j];
	if (fabs(double(ele)) < EPS) {
		return;
	}
	pass = false;
	T ele1 = mat[i][i];
	T ele2 = mat[j][j];
	int size = mat.row();
	T tao = (ele1 - ele2) / (ele * T(2));
	T Tan = sign_without0(tao) / ((T)fabs(double(tao)) + (T)sqrt(double(T(1) + tao * tao)));
	T Cos = T(1) / (T)sqrt(double(T(1) + Tan * Tan));
	T Sin = Cos * Tan;
	Matrix<T> G = Matrix<T>(size);
	G[i][i] = Cos, G[i][j] = -Sin;
	G[j][i] = Sin, G[j][j] = Cos;
	mat = G.transpose() * mat * G;
	J *= G;
}
