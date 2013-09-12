/*
 * Copyright (C) 2013  Abel Walga
 *
 *  random_variable.hpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  random_variable.hpp is distributed in the hope that it will be useful,
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

#ifndef RANDOM_VARIABLE_HPP_
#define RANDOM_VARIABLE_HPP_

#include "range_type.hpp"						// util::dimension, util::Real, util::PositiveInteger
#include "random_object.hpp"					// montecarlo::RandomObject
#include "algebra/matrix.hpp"					// linear_algebra::RealRowVector
#include <gsl/gsl_rng.h>						// gsl_rng, gsl_rng_type, gsl_rng_default, gsl_rng_alloc,gsl_rng_env_setup
#include <gsl/gsl_randist.h>					// gsl_rng_uniform, gsl_ran_gaussian
#include <vector>								// std:vector, std::array
#include <memory>								// std::unique_ptr
namespace randomvariable {

struct GarbageGenerator {
	void operator()(gsl_rng * prng) {
		if (prng) {
			gsl_rng_free(prng);
		}
		prng = nullptr;
	}
};
/**
 * Some random sequence
 */

// define an empirical random sequence
class EmpiricalRandomSequence {
private:
	// underlying data
	std::vector<double> data_set;

	// random number generator
	size_t n; // range
	std::unique_ptr<gsl_rng, GarbageGenerator> rng;

public:
	// sequence characteristics
	typedef double RangeType;
	typedef double ScalarType;

	/**
	 * <p>Constructor</p>
	 * <p>Default constructor builds a deterministic sequence of size one and element 1</p>
	 *
	 * @param const data& empirical data
	 */
	EmpiricalRandomSequence(const std::vector<double>& data = { 0 });
	/**
	 * <p>Constructor/</p>
	 *
	 * @param data&& empirical sequence
	 */
	EmpiricalRandomSequence(std::vector<double>&&);

	/**
	 * <p>Copy constructor</p>
	 */
	EmpiricalRandomSequence(const EmpiricalRandomSequence&);
	/**
	 * <p>Move constructor</p>
	 */
	EmpiricalRandomSequence(EmpiricalRandomSequence&&) noexcept = default;

	/**
	 * <p>Get a sample.</p>
	 */
	RangeType operator()() const noexcept {
		return data_set[gsl_rng_uniform_int(rng.get(), n)];
	}

	/**
	 * <p>Destructor</p>
	 */
	virtual ~EmpiricalRandomSequence() {
	}

};
// end of empirical_random_sequence

typedef montecarlo::RandomObject<EmpiricalRandomSequence,
		EmpiricalRandomSequence::RangeType, EmpiricalRandomSequence::ScalarType> EmpiricalRandomVariable;

// define a normal distribution sequence
class BernoulliSequence {
private:
	// normal distribution
	double p;
	// random number generator
	std::unique_ptr<gsl_rng, GarbageGenerator> rng;

public:
	// sequence characteristics
	typedef unsigned long RangeType;
	typedef unsigned long ScalarType;

	/**
	 * <p>Constructor</p>
	 * <p>Default constructor is the binomial(0.5)</p>
	 *
	 * @param data mean
	 * @param stddev standard deviation
	 */
	explicit BernoulliSequence(double _p = 0.5);
	/**
	 * <p>Copy constructor</p>
	 */
	BernoulliSequence(const BernoulliSequence&);

	/**
	 * <p>Move constructor</p>
	 */
	BernoulliSequence(BernoulliSequence&&) noexcept = default;

	/**
	 * <p>Get a sample.</p>
	 */
	RangeType operator()() const noexcept {
		return gsl_ran_bernoulli(rng.get(), p);
	}

	/**
	 * <p>Destructor</p>
	 */
	virtual ~BernoulliSequence() {
	}
};
// end of BernoulliSequence

typedef montecarlo::RandomObject<BernoulliSequence,
		BernoulliSequence::RangeType, BernoulliSequence::ScalarType> BernoulliRandomVariable;

// define a normal distribution sequence
class BinomialSequence {
private:
	// normal distribution
	double p;
	size_t trial;
	// random number generator
	std::unique_ptr<gsl_rng, GarbageGenerator> rng;

public:
	// sequence characteristics
	typedef unsigned long RangeType;
	typedef unsigned long ScalarType;

	/**
	 * <p>Constructor</p>
	 * <p>Default constructor is the binomial(0.5,1)</p>
	 *
	 * @param data mean
	 * @param stddev standard deviation
	 */
	BinomialSequence(double _p = 0.5, size_t _trial = 1);

	/**
	 * <p>Copy constructor</p>
	 */
	BinomialSequence(const BinomialSequence &);

	/**
	 * <p>Move constructor</p>
	 */
	BinomialSequence(BinomialSequence&&) noexcept = default;

	/**
	 * <p>Get a sample.</p>
	 */
	RangeType operator()() const noexcept {
		return gsl_ran_binomial(rng.get(), p, trial);
	}

	/**
	 * <p>Destructor</p>
	 */
	virtual ~BinomialSequence() {
	}
};
// end of BinomialSequence

typedef montecarlo::RandomObject<BinomialSequence, BinomialSequence::RangeType,
		BinomialSequence::ScalarType> BinomialRandomVariable;

// define a normal distribution sequence
class PoissonSequence {
private:
	// normal distribution
	double mu;
	// random number generator
	std::unique_ptr<gsl_rng, GarbageGenerator> rng;

public:
	// sequence characteristics
	typedef unsigned long RangeType;
	typedef unsigned long ScalarType;

	/**
	 * <p>Constructor</p>
	 * <p>Default constructor is the Poisson(1)</p>
	 *
	 * @param data mean
	 * @param stddev standard deviation
	 */
	explicit PoissonSequence(double _mu = 1.0);

	/**
	 * <p>Copy constructor</p>
	 */
	PoissonSequence(const PoissonSequence&);

	/**
	 * <p>Move constructor</p>
	 */
	PoissonSequence(PoissonSequence&&) noexcept = default;

	/**
	 * <p>Get a sample.</p>
	 */
	RangeType operator()() const noexcept {
		return gsl_ran_poisson(rng.get(), mu);
	}

	/**
	 * <p>Destructor</p>
	 */
	virtual ~PoissonSequence() {
	}
};
// end of PoissonSequence
typedef montecarlo::RandomObject<PoissonSequence, PoissonSequence::RangeType,
		PoissonSequence::ScalarType> PoissonRandomVariable;

// define a centered normal distribution sequence
class NormalSequence {
private:
	// normal distribution
	double sigma;
	// random number generator
	std::unique_ptr<gsl_rng, GarbageGenerator> rng;

public:
	// sequence characteristics
	typedef double RangeType;
	typedef double ScalarType;

	/**
	 * <p>Constructor</p>
	 * <p>Default constructor is the standard normal(0,1.0)</p>
	 *
	 * @param data mean
	 * @param stddev standard deviation
	 */
	NormalSequence(double _stddev = 1.0);

	/**
	 * <p>Copy constructor</p>
	 */
	NormalSequence(const NormalSequence &other);

	/**
	 * <p>Move constructor</p>
	 */
	NormalSequence(NormalSequence &&other) noexcept = default;

	/**
	 * <p>Get a sample.</p>
	 */
	RangeType operator()() const noexcept {
		return gsl_ran_gaussian_ziggurat(rng.get(), sigma);
	}

	/**
	 * <p>Destructor</p>
	 */
	virtual ~NormalSequence() {
	}
};
// end of NormalSequence

typedef montecarlo::RandomObject<NormalSequence, NormalSequence::RangeType,
		NormalSequence::ScalarType> STNRealRandomVariable;

// define a sequence of uncorrelated gaussian random variable vector.
class GaussianRealVectorSequence {

public:
	typedef linear_algebra::RealRowVector RangeType;
	typedef double ScalarType;

	/**
	 * <p>
	 * 		Default constructor
	 * </p>
	 *
	 * @param sigma
	 */
	GaussianRealVectorSequence(const std::vector<double> &_sigma = { 1.0 });

	/**
	 * <p>
	 * 		Copy constructor.
	 * </p>
	 *
	 * @param copy const GaussianRealVectorSequence&
	 */
	GaussianRealVectorSequence(const GaussianRealVectorSequence&);

	/**
	 * <p>
	 * 		Move constructor.
	 * 	</p>
	 *
	 * @param move GaussianRealVectorSequence&&
	 */
	GaussianRealVectorSequence(GaussianRealVectorSequence&&) noexcept = default;

	/**
	 * <p>
	 * 		Generate a random vector
	 * 	</p>
	 *
	 * @return RangeType
	 */
	RangeType operator()() const noexcept {
		linear_algebra::RealRowVector randVec(sigma.size());
		for (size_t i = 0; i < sigma.size(); i++) {
			randVec(i, gsl_ran_gaussian_ziggurat(rng.get(), sigma[i]));
		}
		return randVec;
	}

	/**
	 * <p>
	 * 		Get vector size.
	 * </p>
	 *
	 * @return size_t
	 */
	size_t size() const noexcept {
		return sigma.size();
	}

	~GaussianRealVectorSequence() {
	}
private:
	std::unique_ptr<gsl_rng, GarbageGenerator> rng;
	std::vector<double> sigma;
};
typedef montecarlo::RandomObject<GaussianRealVectorSequence,
		GaussianRealVectorSequence::RangeType,
		GaussianRealVectorSequence::ScalarType> GaussianRandomVector;
// end gaussian real vector sequence
}// end namespace randomvariable
#endif /* RANDOM_VARIABLE_HPP_ */
