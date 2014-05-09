#ifndef _SHANGKE_DATAPOINT_H_
#define _SHANGKE_DATAPOINT_H_

#include "matrix.h"
#include <vector>

using namespace std;

template<typename T> class DataPoint: public vector<T> {
	public:
		/* *
		 * 构造器
		 */
		DataPoint();
		DataPoint(int dim);
		DataPoint(int dim, T value);
		DataPoint(const DataPoint<T> & o);
		DataPoint(const vector<T> & o);

		/* *
		 * 运算符重载
		 */
		// 元素 op 元素
		DataPoint<T> operator + (const DataPoint<T> & o);
		DataPoint<T> operator - (const DataPoint<T> & o);
		T operator * (const DataPoint<T> & o);

		// 元素 op 一般类型
		DataPoint<T> operator + (const T o);
		DataPoint<T> operator - (const T o);
		DataPoint<T> operator * (const T o);
		DataPoint<T> operator / (const T o);

		// 元素 op 元素 并保存
		DataPoint<T> &operator = (const DataPoint<T> & o);
		DataPoint<T> &operator += (const DataPoint<T> & o);
		DataPoint<T> &operator -= (const DataPoint<T> & o);

		// 元素 op 一般类型 并保存
		DataPoint<T> &operator += (const T o);
		DataPoint<T> &operator -= (const T o);
		DataPoint<T> &operator *= (const T o);
		DataPoint<T> &operator /= (const T o);

		// 返回元素的维数
		int dim() const;

		// 返回元素x的||x||1
		T form1() const;

		// 返回元素x的||x||2
		T form2() const;

		// 返回元素x的||x||inf
		T forminf() const;

		// 返回两个元素之间的欧式距离
		T dis(const DataPoint<T> & o) const;

		/* *
		 * 静态函数代码
		 */
		// 得到一个元素集合的重心
		static DataPoint<T> get_mean(const vector< DataPoint<T> > & o);

		// 将一个元素集合映射成矩阵
		static Matrix<T> get_matrix(const vector< DataPoint<T> > & o);

		// 将一个矩阵映射成数据集合
		static vector< DataPoint<T> > get_datapoints(const Matrix<T> & o);
};

#endif
