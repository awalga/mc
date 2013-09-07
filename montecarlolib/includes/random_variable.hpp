/*
 * Copyright (C) 2013  Abel Walga
 *
 *  ${file_name} is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ${file_name} is distributed in the hope that it will be useful,
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

#include <vector>								// std::vector
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
	std::vector<util::Real> data_set;

	// random number generator
	util::dimension n; // range
	std::unique_ptr<gsl_rng, GarbageGenerator> rng;

public:
	// sequence characteristics
	typedef util::Real RangeType;
	typedef util::Real ScalarType;

	/**
	 * <p>Constructor</p>
	 * <p>Default constructor builds a deterministic sequence of size one and element 1</p>
	 *
	 * @param const data& empirical data
	 */
	EmpiricalRandomSequence(const std::vector<util::Real>& data = { 0 });
	/**
	 * <p>Constructor/</p>
	 *
	 * @param data&& empirical sequence
	 */
	EmpiricalRandomSequence(std::vector<util::Real> && data);

	/**
	 * <p>Copy constructor</p>
	 */
	EmpiricalRandomSequence(const EmpiricalRandomSequence&);

	/**
	 * <p>Move constructor</p>
	 */
	EmpiricalRandomSequence(EmpiricalRandomSequence&&);

	/**
	 * <p>Get a sample.</p>
	 */
	RangeType operator()() const {
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
	typedef util::PositiveInteger RangeType;
	typedef util::PositiveInteger ScalarType;

	/**
	 * <p>Constructor</p>
	 * <p>Default constructor is the binomial(0.5)</p>
	 *
	 * @param data mean
	 * @param stddev standard deviation
	 */
	explicit BernoulliSequence(double p = 0.5);

	/**
	 * <p>Copy constructor</p>
	 */
	BernoulliSequence(const BernoulliSequence&);

	/**
	 * <p>Move constructor</p>
	 */
	BernoulliSequence(BernoulliSequence&&);
	/**
	 * <p>Get a sample.</p>
	 */
	RangeType operator()() const {
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
	typedef util::PositiveInteger RangeType;
	typedef util::PositiveInteger ScalarType;

	/**
	 * <p>Constructor</p>
	 * <p>Default constructor is the binomial(0.5,1)</p>
	 *
	 * @param data mean
	 * @param stddev standard deviation
	 */
	BinomialSequence(double p = 0.5, size_t _trial = 1);

	/**
	 * <p>Copy constructor</p>
	 */
	BinomialSequence(const BinomialSequence&);

	/**
	 * <p>Move constructor</p>
	 */
	BinomialSequence(BinomialSequence&&);
	/**
	 * <p>Get a sample.</p>
	 */
	RangeType operator()() const {
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
	typedef util::PositiveInteger RangeType;
	typedef util::PositiveInteger ScalarType;

	/**
	 * <p>Constructor</p>
	 * <p>Default constructor is the Poisson(1)</p>
	 *
	 * @param data mean
	 * @param stddev standard deviation
	 */
	explicit PoissonSequence(double mu = 1.0);

	/**
	 * <p>Copy constructor</p>
	 */
	PoissonSequence(const PoissonSequence&);

	/**
	 * <p>Move constructor</p>
	 */
	PoissonSequence(PoissonSequence&&);
	/**
	 * <p>Get a sample.</p>
	 */
	RangeType operator()() const {
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

// define a normal distribution sequence
class STNRealSequence {
private:
	// normal distribution
	double mu, sigma;
	// random number generator
	std::unique_ptr<gsl_rng, GarbageGenerator> rng;

public:
	// sequence characteristics
	typedef util::Real RangeType;
	typedef util::Real ScalarType;

	/**
	 * <p>Constructor</p>
	 * <p>Default constructor is the standard normal(0,1.0)</p>
	 *
	 * @param data mean
	 * @param stddev standard deviation
	 */
	STNRealSequence(double mean = 0.0, double stddev = 1.0);

	/**
	 * <p>Copy constructor</p>
	 */
	STNRealSequence(const STNRealSequence&);

	/**
	 * <p>Move constructor</p>
	 */
	STNRealSequence(STNRealSequence&&);
	/**
	 * <p>Get a sample.</p>
	 */
	RangeType operator()() const {
		return gsl_ran_gaussian_ziggurat(rng.get(), sigma) + mu;
	}

	/**
	 * <p>Destructor</p>
	 */
	virtual ~STNRealSequence() {
	}
};
// end of STNRealSequence

typedef montecarlo::RandomObject<STNRealSequence, STNRealSequence::RangeType,
		STNRealSequence::ScalarType> STNRealRandomVariable;

// define a real iid normally distributed vector
class IIDSTNRealVectorSequence {
private:
	// rand var dimensions
	util::dimension range_size;
	util::dimension range_dim;

	// normal distribution
	double mu, sigma;

	// random number generator
	std::unique_ptr<gsl_rng, GarbageGenerator> rng;

public:
	// sequence characteristics
	typedef linear_algebra::RealRowVector RangeType;
	typedef util::Real ScalarType;

	/**
	 * <p>Constructor</p>
	 * <p>Default constructor is a sequence of size 1 of a standard normal(0,1.0)</p<
	 */
	IIDSTNRealVectorSequence(util::dimension n = 1, double commonmean = 0.0,
			double stddev = 1.0);

	/**
	 * <p>Copy constructor</p>
	 */
	IIDSTNRealVectorSequence(const IIDSTNRealVectorSequence&);

	/**
	 * <p>Move constructor</p>
	 */
	IIDSTNRealVectorSequence(IIDSTNRealVectorSequence&&);

	/**
	 * <p>Get a sample</p>
	 */
	RangeType operator()() const {
		linear_algebra::RealRowVector randVec(range_size);
		for (size_t i = 0; i < range_size; i++) {
			randVec(i, gsl_ran_gaussian_ziggurat(rng.get(), sigma) + mu);
		}
		return randVec;
	}

	util::dimension size() const {
		return range_size;
	}

	util::dimension dim() const {
		return range_dim;
	}
	/**
	 * <p>Destructor</p>
	 */
	virtual ~IIDSTNRealVectorSequence() {
	}
};

typedef montecarlo::RandomObject<IIDSTNRealVectorSequence,
		IIDSTNRealVectorSequence::RangeType,
		IIDSTNRealVectorSequence::ScalarType> IIDSTNRealRandomVectorVariable;

} // end namespace randomvariable

#endif /* RANDOM_VARIABLE_HPP_ */
