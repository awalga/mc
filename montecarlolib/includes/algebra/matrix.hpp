/*
 * Copyright (C) 2013  Abel Walga
 *
 *  matrix.hpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  matrix.hpp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 20 août 2013
 *      Author: Abel Walga
 */
#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <gsl/gsl_vector.h>					//gsl_vector
#include <gsl/gsl_matrix.h>					//gsl_matrix

#include <cassert>							// std::assert
#include <memory>							// std::unique_ptr

namespace linear_algebra {

/**************************************************************
 * 		Declaration of vector and matrixvectors.
 **************************************************************/
template<typename ScalarType, typename VectorType>
class Vector;
template<typename ScalarType, typename MatrixType>
class Matrix;
/***************************************************
 * 		Generic Vector type traits declaration
 ***************************************************/
struct RowVectorType {
};
struct ColumnVectorType {
};

// vector traits
template<typename VectorType>
struct vector_traits {
};

// row vector traits
template<>
struct vector_traits<RowVectorType> {
	typedef ColumnVectorType TransposeVectorType;
};

// column vector traits
template<>
struct vector_traits<ColumnVectorType> {
	typedef RowVectorType TransposeVectorType;
};
/***************************************************
 * 		Generic vector declaration
 ***************************************************/
/**
 * <p>Scalar row and column vector declaration.</p>
 *
 * <p>Minimal operation set includes index based getter/setter and vector
 * 	  operations addition, subtraction, scaling, product (element wise multiplication), inner product (dot product),
 * 	  outer product, transpose.
 * </p>
 *
 */
template<typename ScalarType, typename VectorType>
class Vector {

public:
	// vector type
	typedef typename vector_traits<VectorType>::TransposeVectorType TransposeType;

	/**
	 * <p>Constructor</p>
	 *
	 * </p>Construct a zero dimension vector</p>
	 */
	Vector();

	/**
	 * <p>Constructor</p>
	 *
	 * @param size vector dimension /=0
	 */
	Vector(size_t size);

	/**
	 * <p>Copy constructor</p>
	 *
	 * @param &v non zero dimension vector
	 */
	Vector(const Vector<ScalarType, VectorType> &v);

	/**
	 * <p>Move constructor</p>
	 *
	 * @param &&v non zero dimension vector
	 */
	Vector(Vector<ScalarType, VectorType> &&v) noexcept;

	/**
	 * <p>expensive copy assignment operator</p>
	 * <p>use move assignment when copy not needed</p>
	 *
	 * <p>Both vectors must be of same dimension or <code>this</code> is a zero dimension vector</p>
	 *
	 * @param &v non zero dimension vector
	 */
	Vector<ScalarType, VectorType>& operator=(
			const Vector<ScalarType, VectorType> &v);

	/**
	 * <p>move assignment operator</p>
	 *
	 * <p>Both vectors must be of same dimension or <code>this</code> is a zero dimension vector</p>
	 * <p>Self assignment through a moved rvalue reference is a no-operation</p>
	 *
	 * @param &&v non zero dimension vector.
	 */
	Vector<ScalarType, VectorType>& operator=(
			Vector<ScalarType, VectorType> && v);

	/***********************************************
	 * 		Vector accessors
	 ***********************************************/
	/**
	 * <p>get ièth elements : throw gsl error if i not in indexe</p>
	 *
	 * @param i index
	 */
	ScalarType& operator[](size_t i);

	/**
	 * <p>get i_th elements : throw gsl error if i not in indexe</p>
	 *
	 * @param i index
	 */
	const ScalarType& operator[](size_t i) const;

	/**
	 * <p>set i-th elements : throw gsl error if i not in indexe</p>
	 *
	 * @param i index
	 * @param v stored value
	 */
	void operator()(size_t i, const ScalarType& v);

	/**
	 * <p>set i-th elements : throw gsl error if i not in indexe</p>
	 *
	 * @param i index
	 * @param v stored value
	 */
	void operator()(size_t i, ScalarType&& v);

	/**
	 * <p>Vector size</p>
	 */
	size_t size();

	/***********************************************
	 * 		Vector operation
	 ***********************************************/
	/**
	 * <p>add two vector of same dimension</p>
	 *
	 * @param v
	 */
	Vector<ScalarType, VectorType>& operator+=(
			const Vector<ScalarType, VectorType>& v);

	/**
	 * <p> add two vector of same dimension</p>
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param right
	 */
	const Vector<ScalarType, VectorType> operator+(
			const Vector<ScalarType, VectorType>& right) const;

	/**
	 * <p>shift vector: transformation affine</p>
	 *
	 * @param shift
	 */
	Vector<ScalarType, VectorType>& operator+=(const ScalarType& shift);

	/**
	 * <p>shift vector: transformation affine</p>
	 *
	 * @param shift
	 */
	Vector<ScalarType, VectorType>& operator+=(ScalarType&& shift);

	/**
	 * <p> add two vector of same dimension</p>
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param shift
	 */
	const Vector<ScalarType, VectorType> operator+(
			const ScalarType& shift) const;
	/**
	 * <p> add two vector of same dimension</p>
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param shift
	 */
	const Vector<ScalarType, VectorType> operator+(ScalarType&& shift) const;
	/**
	 * <p>subtract two vector of same dimension</p>
	 *
	 * @param v
	 */
	Vector<ScalarType, VectorType>& operator-=(
			const Vector<ScalarType, VectorType>& v);

	/**
	 * <p> substract two vector of same dimension</p>
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param right
	 */
	const Vector<ScalarType, VectorType> operator-(
			const Vector<ScalarType, VectorType> &right) const;

	/**
	 * <p>scale by 1/N</p>
	 *
	 * @param N
	 */
	Vector<ScalarType, VectorType>& operator/=(unsigned long N);
	/**
	 * <p>scale by 1/N</p>
	 *
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param N scaling factor
	 */
	const Vector<ScalarType, VectorType> operator/(unsigned long N) const;

	/**
	 * <p>scale by coeff</p>
	 *
	 * @param coeff scaling factor
	 */
	Vector<ScalarType, VectorType>& operator*=(const ScalarType& coeff);

	/**
	 * <p>scale by coeff</p>
	 *
	 * @param coeff scaling factor
	 */
	Vector<ScalarType, VectorType>& operator*=(ScalarType&& coeff);

	/**
	 * <p>scale by coeff</p>
	 *
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param coeff scaling factor
	 */
	const Vector<ScalarType, VectorType> operator*(
			const ScalarType& coeff) const;
	/**
	 * <p>scale by coeff</p>
	 *
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param coeff scaling factor
	 */
	const Vector<ScalarType, VectorType> operator*(ScalarType&& coeff) const;

	/**
	 * <p>vector product : element wise product</p>
	 *
	 * @param right
	 */
	Vector<ScalarType, VectorType>& operator*=(
			const Vector<ScalarType, VectorType>& right);

	/**
	 * <p>element wise product</p>
	 *
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param right
	 */
	const Vector<ScalarType, VectorType> operator*(
			const Vector<ScalarType, VectorType> &right) const;

	/**
	 * <p>Transpose <code>this<code>.</p>
	 *
	 */
	const Vector<ScalarType, TransposeType> transpose() const;

	/**
	 * <p>Fast transpose as does not make a copy but <code>this</code> is discarded</p>
	 *
	 */
	const Vector<ScalarType, TransposeType> transposeInPlace() const;

	/**
	 * <p>Dot product</p>
	 *
	 * @param v columnVectorType
	 */
	ScalarType dot(const Vector<ScalarType, TransposeType> &v) const;
};

/***************************************************
 * 		Scalar Vector definition based on gsl lib.
 ***************************************************/
/**
 * <p><code>ScalarType</code> can be any type handled by <code>gsl_vector<code></p>
 * <p><t> complex type : complex_long_double, complex_double, complex_float</t></p>
 * <p><t> real type    : long_double, double, float</t></p>
 * <p><t> long		   : unsigned_long, long</t></p>
 * <p><t> int		   : unsigned_int,  int</t></p>
 * <p><t> char         : unsigned_char, char</t></p>
 */

/*****************************************************
 * 		Real vector definition based on gsl_vector
 *****************************************************/
template<typename MType>
class Matrix<double, MType> ;
// gsl based vector specialization
template<typename VectorType>
class Vector<double, VectorType> {
protected:

public:

	// transposed vector type
	typedef typename vector_traits<VectorType>::TransposeVectorType TransposeType;

	// define transposed friend
	friend class Vector<double, TransposeType> ;

	// define higher order tensor as friend
	template<typename SType, typename MType>
	friend class Matrix;

	/**
	 * <p>Constructor</p>
	 *
	 * </p>Construct a zero dimension vector</p>
	 */
	Vector() noexcept:
			vector(nullptr, GarbageVector()) {
	} // end default constructor

	/**
	 * <p>Constructor</p>
	 *
	 * @param size vector dimension /=0
	 */
	Vector(size_t size) :
			vector(gsl_vector_alloc(size), GarbageVector()) {
	} // end constructor

	/**
	 * <p>Copy constructor</p>
	 *
	 * @param &v non zero dimension vector
	 */
	Vector(const Vector<double, VectorType> &v) :
			vector(gsl_vector_alloc(v.vector->size), GarbageVector()) {
		// deep copy
		gsl_vector_memcpy(vector.get(), v.vector.get());
	} // end constructor

	/**
	 * <p>Move constructor</p>
	 *
	 * @param &&v non zero dimension vector
	 */
	Vector(Vector<double, VectorType> &&v) noexcept :
	vector(std::move(v.vector)) {
	} // end constructor

	/**
	 * <p>expensive copy assignment operator</p>
	 * <p>use move assignment when copy not needed</p>
	 *
	 * <p>Both vectors must be of same dimension or <code>this</code> is a zero dimension vector</p>
	 *
	 * @param &v non zero dimension vector
	 */
	Vector<double, VectorType>& operator=(const Vector<double, VectorType> &v) {
		// protect against self assignment temporary objects
		gsl_vector * tmpvec = gsl_vector_alloc(v.vector->size);
		// raise exception : GSL error if this->dimension differs from v.dimension.
		gsl_vector_memcpy(tmpvec, v.vector.get());

		// clean properly self and assign new pointer
		vector.reset(tmpvec);

		// return this
		return *this;
	} // end copy assignement

	/**
	 * <p>move assignment operator</p>
	 *
	 * <p>Both vectors must be of same dimension or <code>this</code> is a zero dimension vector</p>
	 * <p>Self assignment through a moved rvalue reference is a no-operation</p>
	 *
	 * @param &&v non zero dimension vector.
	 */
	Vector<double, VectorType>& operator=(Vector<double, VectorType>&& v) noexcept{
		// debug not same dimension
		if(vector) assert(vector->size==v.vector->size);

		// clean local and take other
		vector = std::move(v.vector);

		// return this
		return *this;
	} // end move assignement

	/***********************************************
	 * 		Vector accessors
	 ***********************************************/
	/**
	 * <p>get ièth elements : throw gsl error if i not in indexe</p>
	 *
	 * @param i index
	 */
	double& operator[](size_t i) {
		return *gsl_vector_ptr(vector.get(), i);
	} // end

	/**
	 * <p>get i_th elements : throw gsl error if i not in indexe</p>
	 *
	 * @param i index
	 */
	const double& operator[](size_t i) const {
		return *gsl_vector_const_ptr(vector.get(), i);
	} //

	/**
	 * <p>set i-th elements : throw gsl error if i not in indexe</p>
	 *
	 * @param i index
	 * @param v stored value
	 */
	void operator()(size_t i, const double& v) {
		gsl_vector_set(vector.get(), i, v);
	} // end

	/**
	 * <p>set i-th elements : throw gsl error if i not in indexe</p>
	 *
	 * @param i index
	 * @param v stored value
	 */
	void operator()(size_t i, double&& v) {
		gsl_vector_set(vector.get(), i, std::move(v));
	} // end

	/**
	 * <p>Vector size</p>
	 * <p>0 if null vector</p>
	 */
	constexpr size_t size() const noexcept{
		if(vector) return vector->size;
		else return 0;
	}
	/***********************************************
	 * 		Vector operation
	 ***********************************************/
	/**
	 * <p>add two vector of same dimension</p>
	 *
	 * @param v
	 */
	Vector<double, VectorType>& operator+=(const Vector<double, VectorType>& v) {
		gsl_vector_add(vector.get(), v.vector.get());
		return *this;
	} // end

	/**
	 * <p> add two vector of same dimension</p>
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param right
	 */
	const Vector<double, VectorType> operator+(const Vector<double, VectorType>& right) const {
		Vector<double, VectorType> res = *this;
		res += right;
		return res;
	} // end

	/**
	 * <p>shift vector: transformation affine</p>
	 *
	 * @param shift
	 */
	Vector<double, VectorType>& operator+=(const double& shift) {
		gsl_vector_add_constant(vector.get(), shift);
		return *this;
	} // end

	/**
	 * <p>shift vector: transformation affine</p>
	 *
	 * @param shift
	 */
	Vector<double, VectorType>& operator+=(double&& shift) {
		gsl_vector_add_constant(vector.get(), std::move(shift));
		return *this;
	} // end

	/**
	 * <p> add two vector of same dimension</p>
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param shift
	 */
	const Vector<double, VectorType> operator+(const double& shift) const {
		Vector<double, VectorType> res = *this;
		res += shift;
		return res;
	} // end
	/**
	 * <p> add two vector of same dimension</p>
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param shift
	 */
	const Vector<double, VectorType> operator+(double&& shift) const {
		Vector<double, VectorType> res = *this;
		res += std::move(shift);
		return res;
	} // end
	/**
	 * <p>subtract two vector of same dimension</p>
	 *
	 * @param v
	 */
	Vector<double, VectorType>& operator-=(const Vector<double, VectorType>& v) {
		gsl_vector_sub(vector.get(), v.vector.get());
		return *this;
	} // end

	/**
	 * <p> substract two vector of same dimension</p>
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param right
	 */
	const Vector<double, VectorType> operator-(const Vector<double, VectorType> &right) const {
		Vector<double, VectorType> res = (*this);
		res -= right;
		return res;
	} // end

	/**
	 * <p>scale by 1/N</p>
	 *
	 * @param N
	 */
	Vector<double, VectorType>& operator/=(unsigned long N) {
		gsl_vector_scale(vector.get(), (1.0)/N);
		return *this;
	} // end

	/**
	 * <p>scale by 1/N</p>
	 *
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param N scaling factor
	 */
	const Vector<double, VectorType> operator/(unsigned long N) const {
		Vector<double, VectorType> res = *this;
		res /= N;
		return res;
	} // end

	/**
	 * <p>scale by coeff</p>
	 *
	 * @param coeff scaling factor
	 */
	Vector<double, VectorType>& operator*=(const double& coeff) {
		gsl_vector_scale(vector.get(), coeff);
		return *this;
	} // end
	/**
	 * <p>scale by coeff</p>
	 *
	 * @param coeff scaling factor
	 */
	Vector<double, VectorType>& operator*=(double&& coeff) {
		gsl_vector_scale(vector.get(), std::move(coeff));
		return *this;
	} // end

	/**
	 * <p>scale by coeff</p>
	 *
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param coeff scaling factor
	 */
	const Vector<double, VectorType> operator*(const double& coeff) const {
		Vector<double, VectorType> res = *this;
		res /= coeff;
		return res;
	} // end
	/**
	 * <p>scale by coeff</p>
	 *
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param coeff scaling factor
	 */
	const Vector<double, VectorType> operator*(double&& coeff) const {
		Vector<double, VectorType> res = *this;
		res /= std::move(coeff);
		return res;
	} // end

	/**
	 * <p>vector product : element wise product</p>
	 *
	 * @param right
	 */
	Vector<double, VectorType>& operator*=(const Vector<double, VectorType>& right) {
		gsl_vector_mul(vector.get(), right.vector.get());
		return *this;
	} // end

	/**
	 * <p>element wise product</p>
	 *
	 * <p> performance penalty because involves a copy</p>
	 *
	 * @param right
	 */
	const Vector<double, VectorType> operator*(const Vector<double, VectorType> &right) const {
		Vector<double, VectorType> res = (*this);
		res *= right;
		return res;
	} // end

	/**
	 * <p>Transpose <code>this<code>.</p>
	 *
	 */
	const Vector<double, TransposeType> transpose() const {
		Vector<double, TransposeType> res(vector->size);
		gsl_vector_memcpy(res.vector.get(), vector.get()); // copy this data
		return res;
	} // end transpose

	/**
	 * <p>Fast transpose as does not make a copy but <code>this</code> is discarded</p>
	 *
	 */
	const Vector<double, TransposeType> transposeInPlace() {
		Vector<double, TransposeType> res;
		res.vector.reset(vector.release());
		return res;
	} // end transposeInPlace
	/*
	 * <Dot> product</p>
	 *
	 */
	double dot(const linear_algebra::Vector<double, TransposeType> &v) const;
	/***********************************************
	 * 		End Vector operation
	 ***********************************************/

	/**
	 * <p>Destructor</p>
	 */
	virtual ~Vector() {
	}
private:
	// custom deleter for a double <code>gsl_vector</code>
	struct GarbageVector {
		void operator ()(gsl_vector * pvector) const {
			if (pvector != nullptr) {
				gsl_vector_free(pvector);
				pvector = nullptr;
			}
		}
	};

	// gsl based back end vector for double ranged vectors
	std::unique_ptr<gsl_vector, GarbageVector> vector;
};
// define double column/row vectors
typedef Vector<double, ColumnVectorType> RealColumnVector;
typedef Vector<double, RowVectorType> RealRowVector;

/****************************************************
 * 		Generic Matrix type traits declaration
 ***************************************************/
template<typename MatrixType>
struct matrix_traits {
};

/*
 * <p> rectangular matrix most basic type matrix</p>
 */
struct Rectangular {
};
/**
 * <p> lower triangular matrix</p>
 */
struct LowerTriangular {
};
/**
 * <p> upper triangular matrix</p>
 */
struct UpperTriangular {
};
/*
 * <p> diagonal matrix</p>
 */
struct Diagonal {
};

/**
 * <p> Matrix type traits</p>
 */
template<>
struct matrix_traits<Rectangular> {
	typedef Rectangular TransposeMatrixType;
};
template<>
struct matrix_traits<LowerTriangular> {
	typedef UpperTriangular TransposeMatrixType;
};
template<>
struct matrix_traits<UpperTriangular> {
	typedef LowerTriangular TransposeMatrixType;
};
template<>
struct matrix_traits<Diagonal> {
	typedef Diagonal TransposeMatrixType;
};

/***************************************************
 * 		Generic Matrix declaration
 ***************************************************/
template<typename ScalarType, typename MatrixType>
class Matrix {
	// define matrix traits
	typedef typename matrix_traits<MatrixType>::TransposeMatrixType TransposeType;

	/**
	 * <p>Default Constructor</p>
	 */
	Matrix() noexcept;

	/**
	 * <p>Size based constructor</p>
	 * @param i row_size
	 * @param k column_size
	 */
	Matrix(size_t i, size_t j);

	/**
	 * <p>Copy constructor</p>
	 *
	 * @param matrix to copy from
	 */
	Matrix(const Matrix<ScalarType, MatrixType> &);

	/**
	 * <p>Move constructor</p>
	 *
	 * @param m matrix to take ownership of
	 */
	Matrix(Matrix<ScalarType, MatrixType> &&) noexcept;

	/**
	 * <p>Copy assignement</p>
	 *
	 * @param matrix to copy
	 */
	Matrix<ScalarType, MatrixType>& operator=(
			const Matrix<ScalarType, MatrixType>&);

	/**
	 * <p>Move assignement</p>
	 *	@param matrix to take ownership of
	 */
	Matrix<ScalarType, MatrixType>& operator=(
			Matrix<ScalarType, MatrixType> &&am) noexcept;

	/***********************************************
	 * 		Vector accessors
	 ***********************************************/

	/**
	 * <p>get M(i,j) elements : throw gsl error if i,j not in indexe/p>
	 *
	 * @param i row_index
	 * @param j col_index
	 */
	ScalarType& operator()(size_t i, size_t j);

	/**
	 * <p>get M(i,j) elements : throw gsl error if i,j not in index</p>
	 *
	 * @param i row_index
	 * @param j col_index
	 */
	const ScalarType& operator()(size_t i, size_t j) const;

	/**
	 * <p>set M(i,j) elements : throw gsl error if i,j not in index</p>
	 *
	 * @param i row_index
	 * @param j col_index
	 * @param v value to set
	 */
	void operator()(size_t i, size_t j, const ScalarType& v);

	/**
	 * <p>set M(i,j) elements : throw gsl error if i,j not in index</p>
	 *
	 * @param i row_index
	 * @param j col_index
	 * @param v value to set
	 */
	void operator()(size_t i, size_t j, ScalarType&& v);

	/**
	 * <p>Matrix first dimension size e.g nbr of rows</p>
	 * <p>0 if null matrix</p>
	 */
	constexpr size_t size_i() const noexcept;
	/**
	 * <p>Matrix second dimension size e.g nbr of columns</p>
	 * <p>0 if null matrix</p>
	 */
	constexpr size_t size_j() const noexcept;

	// destructor
	virtual ~Matrix();
};
/***************************************************
 * 		Real Matrix definition based on gsl_matrix
 ***************************************************/
/**
 * <p>Major row matrix implementation</p>
 */
template<typename MatrixType>
class Matrix<double, MatrixType> {
public:

	/**
	 * <p>Default constructor does not initialize</p>
	 *
	 */
	Matrix() noexcept:
			matrix(nullptr, GarbageMatrix()) {
	}
	/**
	 * <p>Square matrix constructor</p>
	 *
	 * <p>Allocate for a i x i matrix</p>
	 *
	 * @param i size
	 */
	Matrix(size_t i);
	/**
	 * <p>Size based constructor</p>
	 *
	 * <p>Allocate for a i x j matrix</p>
	 *
	 * @param i row_size
	 * @param j column_size
	 */
	Matrix(size_t i, size_t j) :
			matrix(gsl_matrix_alloc(i, j), GarbageMatrix()) {
	}

	/**
	 * <p>Deep Copy constructor</p>
	 *
	 * @param m matrix to copy from
	 */
	Matrix(const Matrix<double, MatrixType> &m) :
			matrix(gsl_matrix_alloc(m.matrix->size1, m.matrix->size2),
					GarbageMatrix()) {
		// deep copy
		gsl_matrix_memcpy(matrix.get(), m.matrix.get());
	}

	/**
	 * <p>Move constructor</p>
	 *
	 * @param m matrix to take ownership of
	 */
	Matrix(Matrix<double, MatrixType> &&m) noexcept: matrix(std::move(m.matrix)) {}

	/**
	 * <p>Copy assignement</p>
	 *
	 * <p>Make a deep copy of a matrix of a same dimension or <code>this</code> is a null dimension matrix</p>
	 *
	 * @param m matrix to copy
	 */
	Matrix<double, MatrixType>& operator=(const Matrix<double, MatrixType> &m) {
		// protect self assignement
		gsl_matrix * tmpmatrix = nullptr;
		if(matrix) {
			tmpmatrix = gsl_matrix_alloc(matrix->size1, matrix->size2);
		} else {
			tmpmatrix = gsl_matrix_alloc(m.matrix->size1, m.matrix->size2);
		}
		gsl_matrix_memcpy(tmpmatrix, m.matrix.get());

		// clean properly self and assign new data
		matrix.reset(tmpmatrix);

		// return owned copy
		return *this;
	}

	/**
	 * <p>Move assignement</p>
	 *
	 * <p>Take ownership of other matrix of same dimension</p>
	 *
	 * @param m matrix to take ownership of
	 */
	Matrix<double, MatrixType>& operator=(Matrix<double, MatrixType> &&m) noexcept{
		// assert size debug time only
		if(matrix) {
			assert(matrix->size1==m.matrix->size1);
			assert(matrix->size2==m.matrix->size2);
		}
		// clean properly self and assign new data
		matrix = std::move(m.matrix);

		// return owned copy
		return *this;
	}

	/***********************************************
	 * 		Matrix accessors
	 ***********************************************/
	/**
	 * <p>get M(i,j) elements : throw gsl error if i,j not in index</p>
	 *
	 * @param i row_index
	 * @param j col_index
	 */
	double& operator()(size_t i, size_t j) {
		return *gsl_matrix_ptr(matrix.get(), i, j);
	} // end

	/**
	 * <p>get M(i,j) elements : throw gsl error if i,j not in index</p>
	 *
	 * @param i row_index
	 * @param j col_index
	 */
	const double& operator()(size_t i, size_t j) const {
		return *gsl_matrix_const_ptr(matrix.get(), i, j);
	} //end

	/**
	 * <p>set M(i,j) elements : throw gsl error if i,j not in index</p>
	 *
	 * @param i row_index
	 * @param j col_index
	 * @param v value to set
	 */
	void operator()(size_t i, size_t j, const double& v) {
		gsl_matrix_set(matrix.get(), i, j, v);
	} // end

	/**
	 * <p>set M(i,j) elements : throw gsl error if i,j not in index</p>
	 *
	 * @param i row_index
	 * @param j col_index
	 * @param v value to set
	 */
	void operator()(size_t i, size_t j, double&& v) {
		gsl_matrix_set(matrix.get(), i, j, std::move(v));
	} // end

	/**
	 * <p>Matrix first dimension size e.g nbr of rows</p>
	 * <p>0 if null matrix</p>
	 */
	size_t size_i() const noexcept {
		if(matrix) return matrix->size1;
		else return 0;
	}
	/**
	 * <p>Matrix second dimension size e.g nbr of columns</p>
	 * <p>0 if null matrix</p>
	 */
	size_t size_j() const noexcept {
		if(matrix) return matrix->size2;
		else return 0;
	}
	/***********************************************
	 * 		Matrix operation
	 ***********************************************/
	virtual ~Matrix() {}
private:
	struct GarbageMatrix {
		void operator()(gsl_matrix *pmatrix) const {
			if (pmatrix != nullptr) {
				gsl_matrix_free(pmatrix);
				pmatrix = nullptr;
			}
		}
	};
	std::unique_ptr<gsl_matrix, GarbageMatrix> matrix;
};

} // end namespace linear_algebra
#endif /* MATRIX_HPP_ */
