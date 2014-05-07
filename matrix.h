#ifndef _SHANGKE_MATRIX_H_
#define _SHANGKE_MATRIX_H_

#include <vector>
#include <algorithm>

using namespace std;

template<typename T> class Matrix: public vector< vector<T> > {
	private:
		const static double EPS = 1.0e-8; // 最小单元
		const static int ITERATION = 30;  // 迭代次数的上限
		static bool isZero(const T x);
		static T sign(const T x);
		static T sign_without0(const T x);
	public:
		/* *
		 * 构造器
		 */
		Matrix();
		Matrix(int row, int column);
		Matrix(int row, int column, T value);
		Matrix(int row);
		Matrix(const vector< vector<T> > &mat);
		Matrix(const Matrix &mat);

		/* *
		 * 运算符重载
		 */
		// 矩阵 op 矩阵
		Matrix<T> operator + (const Matrix<T> & o);
		Matrix<T> operator - (const Matrix<T> & o);
		Matrix<T> operator * (const Matrix<T> & o);

		// 矩阵 op 一般类型
		Matrix<T> operator + (const T o);
		Matrix<T> operator - (const T o);
		Matrix<T> operator * (const T o);
		Matrix<T> operator / (const T o);

		// 矩阵 op 矩阵 并保存
		Matrix<T> &operator = (const Matrix<T> & o);
		Matrix<T> &operator += (const Matrix<T> & o);
		Matrix<T> &operator -= (const Matrix<T> & o);
		Matrix<T> &operator *= (const Matrix<T> & o);

		// 矩阵 op 一般类型 并保存
		Matrix<T> &operator += (const T o);
		Matrix<T> &operator -= (const T o);
		Matrix<T> &operator *= (const T o);
		Matrix<T> &operator /= (const T o);

		// 返回矩阵的行
		int row() const;

		// 返回矩阵的列
		int column() const;

		// 返回矩阵行列式的值
		T get_value() const;

		// 返回矩阵的秩
		int rank_det() const;

		// 得到第k行的所有元素
		vector<T> get_row(int k) const;

		// 得到第k列的所有元素
		vector<T> get_column(int k) const;

		// 覆盖第k行的元素
		void set_row(const vector<T>& o, int k);

		// 覆盖第k列的元素
		void set_column(const vector<T>& o, int k);

		// 返回转置的矩阵
		Matrix<T> transpose() const;

		// 转化成列中心矩阵
		Matrix<T> column_center() const;

		// 转化成行中心矩阵
		Matrix<T> row_center() const;

		// 转化成完全中心矩阵
		Matrix<T> center() const;

		// 返回该矩阵的协方差矩阵
		Matrix<T> cov() const;

		// 返回矩阵的逆矩阵
		Matrix<T> inverse() const;

		// 使用传统Qr算法进行奇异值分解
		void svd_qr(Matrix<T>& U, Matrix<T>& B, Matrix<T>& V) const;

		// 使用双边jacobi旋转算法进行奇异值分解
		void svd_jacobi(Matrix<T> &U, Matrix<T>& S, Matrix<T>& V) const;

		/* *
		 * 静态函数代码
		 */
		// 将单位向量x转化成单位向量y的Householder变换
		static Matrix<T> householder(const vector<T>& x, const vector<T>& y);

		// givens转换，给定x, y 得到c, s, r
		static void givens(const T x, const T y, T& c, T& s, T& r);

		// 将v1和v2转化成c * v1 - s * v2和s * v1 + c * v2
		static void update(const T c, const T s, vector<T>& v1, vector<T>& v2);

		// 传统QR迭代转移
		static void qr(vector<T>& delta, vector<T>& gama, Matrix<T>& P, Matrix<T>& Q);

		// 进行jacobi(雅可比)迭代转移
		static void jacobi(Matrix<T>& mat, int size, vector<T>& E, Matrix<T>& J);

		// 进行jacobi迭代中的旋转操作
		static void rotate(Matrix<T>& mat, int i, int j, bool& pass, Matrix<T>& J);
};

#endif
