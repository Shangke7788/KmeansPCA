#ifndef _SHANGKE_DATAPOINT_H_
#define _SHANGKE_DATAPOINT_H_

#include <vector>
#include "matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

using namespace std;

class DataPoint: public vector<double> {
	public:
		/* *
		 * 构造器
		 */
		DataPoint();
		DataPoint(int dim);
		DataPoint(int dim, double value);
		DataPoint(const DataPoint & o);
		DataPoint(const vector<double> & o);

		/* *
		 * 运算符重载
		 */
		// 元素 op 元素
		DataPoint operator + (const DataPoint & o);
		DataPoint operator - (const DataPoint & o);
		double operator * (const DataPoint & o);

		// 元素 op 一般类型
		DataPoint operator + (const double o);
		DataPoint operator - (const double o);
		DataPoint operator * (const double o);
		DataPoint operator / (const double o);

		// 元素 op 元素 并保存
		DataPoint &operator = (const DataPoint & o);
		DataPoint &operator += (const DataPoint & o);
		DataPoint &operator -= (const DataPoint & o);

		// 元素 op 一般类型 并保存
		DataPoint &operator += (const double o);
		DataPoint &operator -= (const double o);
		DataPoint &operator *= (const double o);
		DataPoint &operator /= (const double o);

		// 返回元素的维数
		int dim() const;

		// 返回元素x的||x||1
		double form1() const;

		// 返回元素x的||x||2
		double form2() const;

		// 返回元素x的||x||inf
		double forminf() const;

		// 返回两个元素之间的欧式距离
		double dis(const DataPoint & o) const;

		// 返回两个元素之间的欧式距离的平方
		double dis2(const DataPoint & o) const;

		// 归1化
		void normalize();

		/* *
		 * 静态函数代码
		 */
		// 得到一个元素集合的重心
		static DataPoint get_mean(const vector< DataPoint > & o);

		// 将一个元素集合映射成矩阵
		static Matrix get_matrix(const vector< DataPoint > & o);

		// 将一个矩阵映射成数据集合
		static vector< DataPoint > get_datapoints(const Matrix & o);
};

#ifdef __cplusplus
}
#endif

#endif
