/*
 * Copyright (C) 2013  Abel Walga
 *
 *  matrix.cpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  matrix.cpp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 25 août 2013
 *      Author: Abel Walga
 */

#include <gsl/gsl_vector.h>						//gsl_vector
#include <gsl/gsl_blas.h>						//gsl_blas_x
using namespace std;

#include "../includes/algebra/matrix.hpp"

namespace linear_algebra {
// define dot product between rowVector and columnVector
template<>
double linear_algebra::Vector<double, linear_algebra::RowVectorType>::dot(
		const linear_algebra::Vector<double, linear_algebra::ColumnVectorType> &v) const {
	double res;
	gsl_blas_ddot(vector.get(), v.vector.get(), &res);
	return res;
} // end dot
}

