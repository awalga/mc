/*
 * Copyright (C) 2013  Abel Walga
 *
 *  random_object.hpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  random_object.hpp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 14 août 2013
 *      Author: Abel Walga
 */

#ifndef RANDOM_OBJECT_HPP_
#define RANDOM_OBJECT_HPP_

#include "range_type.hpp"		// common type definition

namespace montecarlo {

/**
 * <p>
 * 		Define a random object based on <code>Sample</code> which gives a sequence of variable of range type <code>RangeType</code>
 *  	and scalar type <code>ScalarType</code>
 * </p>
 *
 * <p>
 *		The random objects defines the monte carlo expectation, variance for any <code>RangeType</code>
 *		and covariance, correlation and convariance when range type is greater than 1.
 *
 *		When <code>RangeType</code> has a dimension greater than 1 then it should define at least the subscript operator: <code>const RangeType& operator()[] const;</code>
 * </p>
 */
template<typename Sample, typename RangeType, typename ScalarType>
class RandomObject {

	// sample generator
	Sample sample;

public:

	/**
	 * <p>Constructor</p>
	 */
	template<typename ...Args>
	RandomObject(Args&&... args) :
			sample(Sample(std::forward<Args>(args)...)) {
	}

	/**
	 * <p>Monte Carlo expectation computed from sample of size N</p>
	 * <p><a name="mean"><b>Mean:<b></a> Computed as \f[E(X,N)=sum(Xi)/N\f]
	 * </p>
	 *
	 * @param N sample size
	 */
	RangeType expectation(unsigned long N) noexcept{
		// initialization
		RangeType x = sample(); // sample
		for (unsigned long n = 1; n < N; n++) {
			x += sample();
		}
		x /= N;
		return x;
	} // end mean

	/**
	 * <p>Monte Carlo variance computed from sample of size N</p>
	 * <p><a name="Variance"><b>Variance:</b></a> Computed as \f[Var(X,N)=E(X^2,N)-E(X,N)^2\f]
	 * </p>
	 *
	 * @param N sample size
	 */
	RangeType variance(unsigned long N) noexcept {
		// initialization
		RangeType x = sample(), mean = x, var = x * x;
		for (unsigned long n = 1; n < N; n++) {
			x = sample(); // use same sample to compute mean and variance
			mean += x;
			var += x * x;
		}
		mean /= N;
		var /= N;
		return var - mean * mean;
	} // end variance

	/**
	 * <p>Monte Carlo covariance of components i, j computed from sample size N</p>
	 * <p><a name="Covariance"><b>Covariance:</b></a> Computed as \f[Cov(X_i,X_j,N)=E(X_i*X_j,N)-E(X_i,N)E(X_j,N)\f]</p>
	 * <p>The type <code>ScalarType</code> must support operators
	 * <code>+=,-=,*=(const ScalarType&amp; u), /=(int N)</code>.</p>
	 *
	 * @param i first component index
	 * @param j second component index
	 * @param N sample size
	 */
	ScalarType covariance(int i, int j, unsigned long N) noexcept{
		// initialization
		RangeType x = sample(); // sample
		ScalarType xi = x[i], xj = x[j], meani = xi, meanj = xj, covarij = xi
				* xj;
		// sample N-1
		for (unsigned long n = 1; n < N; n++) {
			x = sample();
			xi = x[i];
			xj = x[j];
			meani += xi;
			meanj += xj;
			covarij += xi * xj;
		}

		meani /= N;
		meanj /= N;
		covarij /= N;
		return covarij - meani * meanj;
	} // end covariance

	/**
	 * <p>Monte Carlo correlation of components i, j computed from sample size N</p>
	 * <p><a name="Correlation squared"><b>rho^2:</b></a>
	 * Computed as \f[rho_squared(X_i,X_j,N)=(E(X_i*X_j,N)-E(X_i,N)*E(X_j,N))^2/(E(Xi^2,N)-E(Xi,N)^2)*(E(Xj^2,N)-E(Xj,N)^2)
	 * </p>
	 *
	 * @param i first component index
	 * @param j second component index
	 * @param N sample size
	 */
	ScalarType correlation(int i, int j, int N) noexcept{
		// initialization
		RangeType x = sample();
		ScalarType sum_i = x[i], sum_j = x[j], sum_i_j = x[i] * x[j], sum_i_i =
				x[i] * x[i], sum_j_j = x[j] * x[j];

		// sample N-1
		for (unsigned long n = 1; n < N; n++) {
			x = sample();
			sum_i += x[i];
			sum_j += x[j];
			sum_i_j += x[i] * x[j];
			sum_i_i += x[i] * x[i];
			sum_j_j += x[j] * x[j];
		}
		ScalarType num_rho = (N * sum_i_j - sum_i * sum_j);
		ScalarType den_rho = (N * sum_i_i - sum_i * sum_i)
				* (N * sum_j_j - sum_j * sum_j);
		return num_rho * num_rho / den_rho;
	} // end correlation

	/**
	 * <p>Monte Carlo covariance matrix computed from sample size N<p>
	 * <p><a name="CovarianceMatrix"><b>CovarianceMatrix:</b></a>
	 * Computed as \f[C(i,j)=Cov(X_i,X_j)\f]
	 * </p>
	 *
	 * @param N sample size
	 * @param matrix will holds covariance matrix result
	 */
	template<typename UTRMatrix>
	void covarianceMatrix(unsigned long N, UTRMatrix &cov){
		// initialization sample
		RangeType x = sample();
		ScalarType mean_x_i[sample.size()];

		for (size_t i = 0; i < sample.size(); i++) {
			mean_x_i[i] = x[i];
			for (size_t j = 0; j < sample.size(); j++) {
				cov(i, j) = x[i] * x[j];
			}
		}
		// sample N-1
		for (unsigned long n = 0; n < N; n++) {
			x = sample();
			for (size_t i = 0; i < sample.size(); i++) {
				mean_x_i[i] += x[i];
				for (size_t j = 0; j < sample.size(); j++) {
					cov(i, j) += x[i] * x[j];
				}
			}
		}
		// average
		for (size_t i = 0; i < sample.size(); i++) {
			mean_x_i[i] /= N;
			for (size_t j = 0; j < sample.size(); j++) {
				cov(i, j) = cov(i, j) / N - mean_x_i[i] * mean_x_i[j] / (N * N);
			}
		}
	} // end covariance matrix
};

} // end namespace montecarlo

#endif /* RANDOM_OBJECT_HPP_ */
