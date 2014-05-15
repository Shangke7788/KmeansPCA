#ifndef _SHANGKE_MAdoubleRIX_H_
#define _SHANGKE_MAdoubleRIX_H_

#include <vector>
#include <algorithm>

#ifdef __cplusplus
extern "C" {
#endif

using namespace std;

class Matrix: public vector< vector<double> > {
	private:
		const static double EPS = 1.0e-8; // 最小单元
		const static int ITERATION = 30;  // 迭代次数的上限
		static bool isZero(const double x);
		static double sign(const double x);
		static double sign_without0(const double x);
	public:
		/* *
		 * 构造器
		 */
		Matrix();
		Matrix(int row, int column);
		Matrix(int row, int column, double value);
		Matrix(int row);
		Matrix(const vector< vector<double> > &mat);
		Matrix(const Matrix &mat);

		/* *
		 * 运算符重载
		 */
		// 矩阵 op 矩阵
		Matrix operator + (const Matrix & o);
		Matrix operator - (const Matrix & o);
		Matrix operator * (const Matrix & o);

		// 矩阵 op 一般类型
		Matrix operator + (const double o);
		Matrix operator - (const double o);
		Matrix operator * (const double o);
		Matrix operator / (const double o);

		// 矩阵 op 矩阵 并保存
		Matrix &operator = (const Matrix & o);
		Matrix &operator += (const Matrix & o);
		Matrix &operator -= (const Matrix & o);
		Matrix &operator *= (const Matrix & o);

		// 矩阵 op 一般类型 并保存
		Matrix &operator += (const double o);
		Matrix &operator -= (const double o);
		Matrix &operator *= (const double o);
		Matrix &operator /= (const double o);

		// 返回矩阵的行
		int row() const;

		// 返回矩阵的列
		int column() const;

		// 返回矩阵行列式的值
		double get_value() const;

		// 返回矩阵的秩
		int rank_det() const;

		// 得到第k行的所有元素
		vector<double> get_row(int k) const;

		// 得到第k列的所有元素
		vector<double> get_column(int k) const;

		// 覆盖第k行的元素
		void set_row(const vector<double>& o, int k);

		// 覆盖第k列的元素
		void set_column(const vector<double>& o, int k);

		// 返回转置的矩阵
		Matrix transpose() const;

		// 转化成列中心矩阵
		Matrix column_center() const;

		// 转化成行中心矩阵
		Matrix row_center() const;

		// 转化成完全中心矩阵
		Matrix center() const;

		// 返回该矩阵的协方差矩阵
		Matrix cov() const;

		// 返回矩阵的逆矩阵
		Matrix inverse() const;

		// 使用传统Qr算法进行奇异值分解
		void svd_qr(Matrix& U, Matrix& B, Matrix& V) const;

		// 使用双边jacobi旋转算法进行奇异值分解
		void svd_jacobi(Matrix &U, Matrix& S, Matrix& V) const;

		/* *
		 * 静态函数代码
		 */
		// 将单位向量x转化成单位向量y的Householder变换
		static Matrix householder(const vector<double>& x, const vector<double>& y);

		// givens转换，给定x, y 得到c, s, r
		static void givens(const double x, const double y, double& c, double& s, double& r);

		// 将v1和v2转化成c * v1 - s * v2和s * v1 + c * v2
		static void update(const double c, const double s, vector<double>& v1, vector<double>& v2);

		// 传统QR迭代转移
		static void qr(vector<double>& delta, vector<double>& gama, Matrix& P, Matrix& Q);

		// 进行jacobi(雅可比)迭代转移
		static void jacobi(Matrix& mat, int size, vector<double>& E, Matrix& J);

		// 进行jacobi迭代中的旋转操作
		static void rotate(Matrix& mat, int i, int j, bool& pass, Matrix& J);
};

#ifdef __cplusplus
}
#endif

#endif
