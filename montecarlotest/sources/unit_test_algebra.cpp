/*
 * Copyright (C) 2013  Abel Walga
 *
 *  unit_test_algebra.cpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  unit_test_algebra.cpp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 22 août 2013
 *      Author: Abel Walga
 */

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#define BOOST_TEST_MODULE ThreadSafeContainerRegressionTest
#endif // STAND_ALONE
#include "algebra/matrix.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

using namespace std;

BOOST_AUTO_TEST_SUITE(AlgebraPckTest)

typedef boost::mpl::list<double> template_type;

template<typename T>
T value(int i) {
	return T();
}

float f_test_value[] = { 0.0f, 1.0f, 2.0f, 3.0f };
template<>
inline float value<float>(int i) {
	return f_test_value[i];
}
double d_test_value[] = { 0.0, 1.0, 2.0, 3.0 };
template<>
inline double value<double>(int i) {
	return d_test_value[i];
}
long double ld_test_value[] = { 0.0, 1.0, 2.0, 3.0 };
template<>
inline long double value<long double>(int i) {
	return ld_test_value[i];
}

	// 1. minimal vector regression test
BOOST_AUTO_TEST_CASE_TEMPLATE(VectorRegressionTestCase, T, template_type){

	using namespace linear_algebra;
	typedef Vector<T, RowVectorType> TestVector;
	typedef Vector<T, ColumnVectorType> TestVectorTrans;
	TestVector vec_1(3), vec_2(3);
	vec_1(0, value<T>(0));vec_2(0, value<T>(0));
	vec_1(1, value<T>(1));vec_2(1, value<T>(1));
	vec_1(2, value<T>(2));vec_2(2, value<T>(2));

	// check size
	BOOST_CHECK_EQUAL(vec_1.size(), 3);

	// check vector addition
	TestVector sum = vec_1+vec_2;
	BOOST_CHECK_EQUAL(sum[1], value<T>(1)+value<T>(1));
	// check addition scalar
	sum = vec_1+value<T>(1);
	BOOST_CHECK_EQUAL(sum[2], value<T>(2)+value<T>(1));

	// check substraction
	TestVector sub = vec_1-vec_2;
	BOOST_CHECK_EQUAL(sub[1], value<T>(2)-value<T>(2));

	// check scaling
	TestVector scl = vec_1*value<T>(1);
	BOOST_CHECK_EQUAL(scl[2], value<T>(2)*value<T>(1));

	// check vector multiplication;
	TestVector mul = vec_1*vec_2;
	BOOST_CHECK_EQUAL(mul[0], value<T>(0)*value<T>(0));
	BOOST_CHECK_EQUAL(mul[1], value<T>(1)*value<T>(1));
	BOOST_CHECK_EQUAL(mul[2], value<T>(2)*value<T>(2));

	// check transpose
	TestVectorTrans t_vec_1 = vec_1.transpose();
	BOOST_CHECK_EQUAL(t_vec_1[0], value<T>(0));
	BOOST_CHECK_EQUAL(t_vec_1[1], value<T>(1));
	BOOST_CHECK_EQUAL(t_vec_1[2], value<T>(2));

	// check transpose in place
	TestVectorTrans t_vec_2 = vec_1.transposeInPlace();
	BOOST_CHECK_EQUAL(t_vec_2[0], value<T>(0));
	BOOST_CHECK_EQUAL(t_vec_2[1], value<T>(1));
	BOOST_CHECK_EQUAL(t_vec_2[2], value<T>(2));

	// check inner product
	auto res = vec_2.dot(t_vec_2);
	BOOST_CHECK_EQUAL(res, 5);
	// check outer product
}
	// 2. vector test
BOOST_AUTO_TEST_CASE_TEMPLATE(VectorTestCase, T, template_type){

}
	// 1. minimal matrix regression test
BOOST_AUTO_TEST_CASE_TEMPLATE(MatrixRegressionTestCase, T, template_type){
	using namespace linear_algebra;
	typedef Matrix<T, Rectangular> TestMatrix;

	TestMatrix m_1(2,2);
	m_1(0,0,value<T>(0));m_1(0,1,value<T>(1));
	m_1(1,0,value<T>(2));m_1(1,1,value<T>(3));

	//
	BOOST_CHECK_EQUAL(m_1.size_i(), 2);
	BOOST_CHECK_EQUAL(m_1.size_j(), 2);
	BOOST_CHECK_EQUAL(m_1(0,0), 0);
	BOOST_CHECK_EQUAL(m_1(1,1), 3);
}
// 2. matrix test
BOOST_AUTO_TEST_CASE_TEMPLATE(MatrixTestCase, T, template_type){

}
BOOST_AUTO_TEST_SUITE_END()

