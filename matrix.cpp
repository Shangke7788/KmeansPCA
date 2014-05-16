#include "matrix.h"

#include <vector>
#include <algorithm>
#include <math.h>
#include <map>
#include <stdio.h>

using namespace std;

bool Matrix::isZero(const double x) {
	double y = (double)x;
	return ((y > EPS) - (y < -EPS)) == 0;
}

double Matrix::sign(const double x) {
	double y = (double)x;
	return double((y > EPS) - (y < -EPS));
}

double Matrix::sign_without0(const double x) {
	if (x < double(0)) {
		return double(-1);
	} else {
		return double(1);
	}
}

Matrix::Matrix(): vector< vector<double> >() {
}

Matrix::Matrix(int row, int column): vector< vector<double> >(row, vector<double>(column, double(0))) {
}

Matrix::Matrix(int row): vector< vector<double> >(row, vector<double>(row, double(0))) {
	for (int i = 0; i < row; i++) {
		(*this)[i][i] = double(1);
	}
}

Matrix::Matrix(int row, int column, double value): vector< vector<double> >(row, vector<double>(column, value)) {
}

Matrix::Matrix(const vector< vector<double> > &mat): vector< vector<double> >(mat) {
}

Matrix::Matrix(const Matrix &mat): vector< vector<double> >(mat) {
}

Matrix Matrix::operator + (const Matrix & o) {
	Matrix ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] + o[i][j];
		}
	}
	return ans;
}

Matrix Matrix::operator - (const Matrix & o) {
	Matrix ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] - o[i][j];
		}
	}
	return ans;
}

Matrix Matrix::operator * (const Matrix & o) {
	Matrix ans = Matrix(this->row(), o.column());
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

Matrix Matrix::operator + (const double o) {
	Matrix ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] + o;
		}
	}
	return ans;
}

Matrix Matrix::operator - (const double o) {
	Matrix ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] - o;
		}
	}
	return ans;
}

Matrix Matrix::operator * (const double o) {
	Matrix ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] * o;
		}
	}
	return ans;
}

Matrix Matrix::operator / (const double o) {
	Matrix ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] / o;
		}
	}
	return ans;
}

Matrix& Matrix::operator = (const Matrix & o) {
	if (this == &o) {
		return *this;
	}
	int R = this->row(), C = this->column();
	for (int i = 0; i < R; i++) {
		(*this)[i].clear();
	}
	this->clear();
	for (int i = 0; i < o.row(); i++) {
		this->push_back(vector<double>(o[i]));
	}
	return *this;
}

Matrix& Matrix::operator += (const Matrix & o) {
	Matrix ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] + o[i][j];
		}
	}
	return (*this) = ans;
}

Matrix& Matrix::operator -= (const Matrix & o) {
	Matrix ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] - o[i][j];
		}
	}
	return (*this) = ans;
}

Matrix& Matrix::operator *= (const Matrix & o) {
	Matrix ans = Matrix(this->row(), o.column());
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

Matrix& Matrix::operator += (const double o) {
	Matrix ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] + o;
		}
	}
	return (*this) = ans;
}

Matrix& Matrix::operator -= (const double o) {
	Matrix ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] - o;
		}
	}
	return (*this) = ans;
}

Matrix& Matrix::operator *= (const double o) {
	Matrix ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] * o;
		}
	}
	return (*this) = ans;
}

Matrix& Matrix::operator /= (const double o) {
	Matrix ans = Matrix(this->row(), this->column());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[i][j] / o;
		}
	}
	return (*this) = ans;
}

int Matrix::row() const {
	return (int)(this->size());
}

int Matrix::column() const {
	if (this->row() == 0) {
		return 0;
	} else {
		return (int)(this->begin()->size());
	}
}

double Matrix::get_value() const {
	int R = this->row(), C = this->column();
	if (R != C || R == 0) {
		return double(0);
	}
	Matrix mat = *this;
	double ans = double(1);
	int n = R;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (!isZero(mat[j][i])) {
				double t = mat[i][i] / mat[j][i];
				for (int k = i; k < n; k++) {
					mat[i][k] -= t * mat[j][k];
				}
				for (int k = i; k < n; k++) {
					std::swap(mat[i][k], mat[j][k]);
				}
				ans = -ans;
			}
		}
		if (isZero(mat[i][i])) {
			return double(0);
		} else {
			ans = ans * mat[i][i];
		}
	}
	return ans;
}

int Matrix::rank_det() const {
	Matrix mat = *this;
	int d = 0, R = mat.row(), C = mat.column();
	int ni = 0;
	for (int i = 0; i < R; i++) {
		for (int j = i + 1; j < R; j++) {
			if (!isZero(mat[j][ni])) {
				double t = mat[i][ni] / mat[j][ni];
				for (int k = ni; k < C; k++) {
					mat[i][k] -= t * mat[j][k];
				}
				for (int k = ni; k < C; k++) {
					std::swap(mat[i][k], mat[j][k]);
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

vector<double> Matrix::get_row(int k) const {
	if (k >= this->row()) {
		return vector<double>();
	}
	vector<double> ans = (*this)[k];
	return ans;
}

vector<double> Matrix::get_column(int k) const {
	if (k >= this->column()) {
		return vector<double>();
	}
	int R = this->row();
	vector<double> ans = vector<double>(R, double(0));
	for (int i = 0; i < R; i++) {
		ans[i] = (*this)[i][k];
	}
	return ans;
}

void Matrix::set_row(const vector<double>& o, int k) {
	if (k >= this->row()) {
		return;
	}
	int C = this->column();
	for (int i = 0; i < C; i++) {
		(*this)[k][i] = o[i];
	}
}

void Matrix::set_column(const vector<double>& o, int k) {
	if (k >= this->column()) {
		return;
	}
	int R = this->row();
	for (int i = 0; i < R; i++) {
		(*this)[i][k] = o[i];
	}
}

Matrix Matrix::transpose() const {
	Matrix ans = Matrix(this->column(), this->row());
	int R = ans.row(), C = ans.column();
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = (*this)[j][i];
		}
	}
	return ans;
}

Matrix Matrix::column_center() const {
	Matrix ans = Matrix(this->row(), this->column());
	vector<double> ave = vector<double>(this->column(), double(0));
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

Matrix Matrix::row_center() const {
	Matrix ans = Matrix(this->row(), this->column());
	vector<double> ave = vector<double>(this->column(), double(0));
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

Matrix Matrix::center() const {
	Matrix ans = Matrix(this->row(), this->column());
	vector<double> r_ave = vector<double>(this->column(), double(0));
	vector<double> c_ave = vector<double>(this->row(), double(0));
	double ave = double(0);
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

Matrix Matrix::cov() const {
	Matrix cen = *this;
	Matrix c2 = cen.transpose();
	Matrix ans = c2 * cen;
	return ans;
}

Matrix Matrix::cov2() const {
	Matrix cen = *this;
	Matrix c2 = cen.transpose();
	Matrix ans = cen * c2;
	return ans;
}

Matrix Matrix::inverse() const {
	int R = this->row(), C = this->column();
	if (R != C || R == 0) {
		return Matrix();
	}
	Matrix mat = Matrix(R, C + C);
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
			return Matrix();
		}
		double sav = mat[s][i];
		for (int j = i; j < (C << 1); j++) {
			std::swap(mat[i][j], mat[s][j]);
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
			double sav = -mat[j][i];
			for (int k = i; k < (C << 1); k++) {
				mat[j][k] += sav * mat[i][k];
			}
		}
	}
	Matrix ans = Matrix(R, C);
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			ans[i][j] = mat[i][j + C];
		}
	}
	return ans;
}

void Matrix::svd_qr(Matrix& U, Matrix& B, Matrix& V) const {
	Matrix A = *this;
	bool hastran = false;
	int m = this->row(), n = this->column();
	if (m < n) {
		A = A.transpose();
		hastran = true;
		std::swap(m, n);
	}
	// 转化为 主对角线+副对角线 矩阵
	vector<double> delta = vector<double>(n, double(0));
	vector<double> gama = vector<double>(n, double(0));
	U = Matrix(m), V = Matrix(n);
	for (int k = 0; k < n; k++) {
		double cc = double(0);
		for (int i = 0; i < m - k; i++) {
			cc += A[i][0] * A[i][0];
		}
		cc = (double)sqrt(double(cc));
		delta[k] = cc;
		vector<double> v1 = A.get_column(0);
		vector<double> e1 = vector<double>(m - k, double(0));
		e1[0] = cc;
		Matrix P1 = householder(v1, e1);
		Matrix tmp = Matrix(m - k, n - k - 1);
		int R = tmp.row(), C = tmp.column();
		for (int i = 0; i < R; i++) {
			for (int j = 0; j < C; j++) {
				tmp[i][j] = A[i][j + 1];
			}
		}
		A = P1 * tmp;
		tmp = Matrix(m);
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

		cc = double(0);
		for (int i = 0; i < n - k - 1; i++) {
			cc += A[0][i] * A[0][i];
		}
		cc = (double)sqrt(double(cc));
		gama[k + 1] = cc;
		vector<double> u1 = A.get_row(0);
		e1 = vector<double>(n - k - 1, double(0));
		e1[0] = cc;
		Matrix H1 = householder(u1, e1);
		tmp = Matrix(m - k - 1, n - k - 1);
		R = tmp.row(), C = tmp.column();
		for (int i = 0; i < R; i++) {
			for (int j = 0; j < C; j++) {
				tmp[i][j] = A[i + 1][j];
			}
		}
		A = tmp * H1;

		tmp = Matrix(n);
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
		gama[0] = double(0);
		int p = n - 1, q = 1;
		for (int i = 1; i < n; i++) {
			if (fabs(double(gama[i])) <= EPS * (fabs(double(delta[i - 1])) + fabs(double(delta[i])))) {
				gama[i] = double(0);
			} else {
				flag = true;
				p = min(p, i), q = max(q, i);
			}
		}
		if (!flag) {
			break;
		}
		p--;
		double max_val = double(0);
		for (int i = 0; i < n; i++) {
			max_val = max((double)fabs(double(delta[i])), max_val);
			max_val = max((double)fabs(double(gama[i])), max_val);
		}
		bool has3 = false;
		int l = 0, L = 0;
		double x, y;
		for (int i = p; i < q; i++) {
			if ((double)fabs(double(delta[i])) <= max_val * (double)EPS) {
				delta[i] = 0, x = gama[i + 1], y = delta[i + 1];
				gama[i + 1] = 0, l = 1, L = i;
				has3 = true;
				break;
			}
		}
		if (has3) {
			while (true) {
				double c, s, r;
				givens(y, x, c, s, r);
				delta[L + l] = r;
				Matrix G = Matrix(m);
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
			Matrix P, Q;
			vector<double> d = vector<double>(q - p + 1, double(0));
			vector<double> g = vector<double>(q - p + 1, double(0));
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
			Matrix tmp = Matrix(m);
			for (int i = p; i <= q; i++) {
				for (int j = p; j <= q; j++) {
					tmp[i][j] = P[i - p][j - p];
				}
			}
			U *= tmp;
			tmp = Matrix(n);
			for (int i = p; i <= q; i++) {
				for (int j = q; j <= q; j++) {
					tmp[i][j] = Q[i - p][j - p];
				}
			}
			V *= tmp;
		}
	}

	B = Matrix(m, n);
	for (int i = 0; i < n; i++) {
		B[i][i] = delta[i];
		if (i == n - 1) {
			break;
		}
		B[i][i + 1] = gama[i + 1];
	}
	if (hastran) {
		std::swap(U, V), U = U.transpose(), V = V.transpose(), B = B.transpose();
	}
}

void Matrix::svd_jacobi(Matrix& U, Matrix& S, Matrix& V) const {
	Matrix A = *this;
	bool hastran = false;
	int R = this->row(), C = this->column();
	if (R > C) {
		A = A.transpose();
		hastran = true;
		std::swap(R, C);
	}
	U = Matrix(R, R), V = Matrix(C, C);
	vector<double> E = vector<double>(C, 0);

	Matrix B = A.transpose() * A;
	Matrix J = Matrix(C);
	vector<double> ss = vector<double>(C, 0);
	jacobi(B, C, ss, J);
	for (int i = 0; i < (int)ss.size(); i++) {
		ss[i] = (double)sqrt(double(ss[i]));
	}

	multimap<double, int> eigen;
	for (int i = 0; i < (int)ss.size(); i++) {
		eigen.insert(make_pair(ss[i], i));
	}
	typename multimap<double, int>::iterator iter = --eigen.end();
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
		Matrix vi = vector< vector<double> >(1, V.get_column(i));
		double sigma = E[i];
		Matrix ui = A * vi.transpose();
		for (int j = 0; j < R; j++) {
			U[j][i] = ui[j][0] / sigma;
		}
	}

	//U矩阵的后(rows-none_zero)列就不计算了，采用默认值0。因为这后几列对应的奇异值为0,在做数据压缩时用不到这几列

	S = Matrix(R, C);
	for (int i = 0; i < R; i++) {
		S[i][i] = E[i];
	}
	if (hastran) {
		S = S.transpose(), std::swap(U, V), U = U.transpose(), V = V.transpose();
	}
}

Matrix Matrix::householder(const vector<double> &x, const vector<double> &y) {
	int n = x.size();
	if (n == 0) {
		return Matrix();
	}
	double C = double(0);
	for (int i = 0; i < n; i++) {
		C += (x[i] - y[i]) * (x[i] - y[i]);
	}
	C = (double)sqrt(double(C));
	Matrix H = Matrix(n);
	if (isZero(C)) {
		return H;
	}
	Matrix B = Matrix(n, 1);
	for (int i = 0; i < n; i++) {
		B[i][0] = (x[i] - y[i]) / C;
	}
	Matrix BT = B.transpose();
	H -= B * BT * double(2);
	return H;
}

void Matrix::givens(const double x, const double y, double& c, double& s, double& r) {
	if (isZero(y)) {
		c = double(1), s = double(0), r = double(0);
	} else if (fabs(double(y)) > fabs(double(x))) {
		double t = -x / y;
		s = (double)sqrt(double(double(1) + t * t));
		r = -y * s;
		s = double(1) / s;
		c = s * t;
	} else {
		double t = -y / x;
		c = (double)sqrt(double(double(1) + t * t));
		r = x * c;
		c = double(1) / c;
		s = c * t;
	}
}

void Matrix::update(const double c, const double s, vector<double>& v1, vector<double> &v2) {
	int n = v1.size();
	for (int i = 0; i < n; i++) {
		double t = v1[i];
		v1[i] = c * t - s * v2[i];
		v2[i] = s * t + c * v2[i];
	}
}

void Matrix::qr(vector<double>& delta, vector<double>& gama, Matrix& P, Matrix& Q) {
	int n = delta.size();
	double d = ((delta[n - 2] * delta[n - 2] + gama[n - 2] * gama[n - 2]) - 
			delta[n - 1] * delta[n - 1] + gama[n - 1] * gama[n - 1]) / double(2);
	double mu = (delta[n - 1] * delta[n - 1] + gama[n - 1] * gama[n - 1]) + d - sign(d)
		* (double)sqrt(double(d * d + delta[n - 2] * delta[n - 2] *gama[n - 1] * gama[n - 1]));
	double x = delta[0] * delta[0] - mu, y = delta[0] * gama[1], k = 0;
	Q = Matrix(n), P = Matrix(n);
	while (true) {
		double c, s, r;
		givens(x, y, c, s, r);
		Matrix cs = Matrix(2, 2);
		cs[0][0] = c, cs[0][1] = s;
		cs[1][0] = -s, cs[1][1] = c;
		Matrix tmp = Matrix(2, 2);
		tmp[0][0] = delta[k], tmp[0][1] = gama[k + 1];
		tmp[1][0] = 0, tmp[1][1] = delta[k + 1];
		tmp *= cs;
		x = tmp[0][0], gama[k + 1] = tmp[0][1];
		y = tmp[1][0], delta[k + 1] = tmp[1][1];

		vector<double> k0 = Q.get_column(k);
		vector<double> k1 = Q.get_column(k + 1);
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
			tmp = Matrix(2, 2);
			tmp[0][0] = gama[k + 1], tmp[0][1] = double(0);
			tmp[1][0] = delta[k + 1], tmp[1][1] = gama[k + 2];
			cs = Matrix(2, 2);
			cs[0][0] = c, cs[0][1] = s;
			cs[1][0] = -s, cs[1][1] = c;
			tmp = cs.transpose() * tmp;
			x = tmp[0][0], y = tmp[0][1];
			delta[k + 1] = tmp[1][0], gama[k + 2] = tmp[1][1];
			k = k + 1;
		} else {
			tmp = Matrix(2, 1);
			tmp[0][0] = gama[n - 1], tmp[1][0] = delta[n - 1];
			cs = Matrix(2, 2);
			cs[0][0] = c, cs[0][1] = s;
			cs[1][0] = -s, cs[1][1] = c;
			tmp = cs.transpose() * tmp;
			break;
		}
	}
}

void Matrix::jacobi(Matrix& mat, int size, vector<double>& E, Matrix& J) {
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

void Matrix::rotate(Matrix& mat, int i, int j, bool& pass, Matrix& J) {
	double ele = mat[i][j];
	if (fabs(double(ele)) < EPS) {
		return;
	}
	pass = false;
	double ele1 = mat[i][i];
	double ele2 = mat[j][j];
	int size = mat.row();
	double tao = (ele1 - ele2) / (ele * double(2));
	double Tan = sign_without0(tao) / ((double)fabs(double(tao)) + (double)sqrt(double(double(1) + tao * tao)));
	double Cos = double(1) / (double)sqrt(double(double(1) + Tan * Tan));
	double Sin = Cos * Tan;
	Matrix G = Matrix(size);
	G[i][i] = Cos, G[i][j] = -Sin;
	G[j][i] = Sin, G[j][j] = Cos;
	mat = G.transpose() * mat * G;
	J *= G;
}

int Matrix::svn(Matrix& S, Matrix& U, Matrix& V) const {
	Matrix A = *this;
	int R = this->row(), C = this->column();
	bool hastran = false;
	if (R > C) {
		A = A.transpose(), hastran = true;
		std::swap(R, C);
	}
	U = Matrix(R, R), V = Matrix(C), S = Matrix(R, C);

	A.hestens_jacobi(V);

	vector<double> E(C, 0.0);
	int none_zero = 0;
	for (int i = 0; i < C; i++) {
		vector<double> s1 = A.get_column(i);
		double norm = 0.0;
		for (int j = 0; j < R; j++) {
			norm += s1[j] * s1[j];
		}
		norm = sqrt(norm);
		if (norm > EPS) {
			none_zero++;
		}
		E[i] = norm;
	}

	/* *
	 * U矩阵的后(R - none_zero)列以及V的后(C - none_zero)列就不计算了，采用默认值0
	 * 对于奇异值分解A=U*Sigma*V^T，我们只需要U的前r列，V^T的前r行(即V的前r列)，就可以恢复A了。r是A的秩。
	 */
	for (int r = 0; r < R; r++) {
		S[r][r] = E[r];
		for (int c = 0; c < none_zero; c++) {
			U[r][c] = A[r][c] / E[c];
		}
	}

	if (hastran) {
		S = S.transpose(), std::swap(U, V), U = U.transpose(), V = V.transpose();
	}

	return none_zero;
}

void Matrix::hestens_jacobi(Matrix& V) {
	int R = this->row(), C = this->column();

	int iteration = ITERATION;
	while (iteration--) {
		bool pass = true;
		for (int i = 0; i < C; i++) {
			for (int j = i + 1; j < C; j++) {
				this->orthogonal(i, j, pass, V);
			}
		}
		if (pass) {
			break;
		}
	}
}

void Matrix::orthogonal(int i, int j, bool& pass, Matrix& V) {
	int R = this->row(), C = this->column();
	vector<double> ci = this->get_column(i);
	vector<double> cj = this->get_column(j);
	double ele = 0.0;
	for (int r = 0; r < R; r++) {
		ele += ci[r] * cj[r];
	}
	if (fabs(ele) < EPS) {
		return;
	}

	pass = false;
	double ele1 = 0.0;
	double ele2 = 0.0;
	for (int r = 0; r < R; r++) {
		ele1 += ci[r] * ci[r];
		ele2 += cj[r] * cj[r];
	}

	if (ele1 < ele2) {
		for (int r = 0; r < R; r++) {
			(*this)[r][i] = cj[r];
			(*this)[r][j] = ci[r];
		}
		for (int r = 0; r < C; r++) {
			double tmp = V[r][i];
			V[r][i] = V[r][j];
			V[r][j] = tmp;
		}
	}

	double tao = (ele1 - ele2) / (2 * ele);
	double Tan = sign_without0(tao) / (fabs(tao) + sqrt(1 + tao * tao));
	double Cos = 1.0 / sqrt(1 + Tan * Tan);
	double Sin = Cos * Tan;

	for (int r = 0; r < R; r++) {
		double var1 = (*this)[r][i] * Cos + (*this)[r][j] * Sin;
		double var2 = (*this)[r][j] * Cos - (*this)[r][i] * Sin;
		(*this)[r][i] = var1;
		(*this)[r][j] = var2;
	}

	for (int c = 0; c < C; c++) {
		double var1 = V[c][i] * Cos + V[c][j] * Sin;
		double var2 = V[c][j] * Cos - V[c][i] * Sin;
		V[c][i] = var1;
		V[c][j] = var2;
	}
}
