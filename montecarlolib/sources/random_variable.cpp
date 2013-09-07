/*
 * Copyright (C) 2013  Abel Walga
 *
 *  random_variable.cpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  random_variable.cpp is distributed in the hope that it will be useful,
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

using namespace std;

#include "../includes/range_type.hpp"
#include "../includes/random_variable.hpp"

using namespace util;
using namespace randomvariable;

// EmpiricalRandomSequence(data&) constructor
EmpiricalRandomSequence::EmpiricalRandomSequence(const vector<Real>& data) :
		data_set(data), n(data.size() - 1) {
	// get generator type from env var
	gsl_rng_env_setup();
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(
			gsl_rng_alloc(gsl_rng_default), GarbageGenerator());
} // end EmpiricalRandomSequence constructor

// EmpiricalRandomSequence(data&&) constructor
EmpiricalRandomSequence::EmpiricalRandomSequence(vector<Real> && data) :
		data_set(std::move(data)), n(data_set.size() - 1) {
	// get generator type from env var
	gsl_rng_env_setup();
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(
			gsl_rng_alloc(gsl_rng_default), GarbageGenerator());
} // end EmpiricalRandomSequence constructor

// EmpiricalRandomSequence copy constructor
EmpiricalRandomSequence::EmpiricalRandomSequence(
		const EmpiricalRandomSequence &other) :
		data_set(other.data_set), n(other.n) {
	gsl_rng * copy = gsl_rng_clone(other.rng.get());
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(copy, GarbageGenerator());
} // end EmpiricalRandomSequence copy constructor

// EmpiricalRandomSequence move constructor
EmpiricalRandomSequence::EmpiricalRandomSequence(
		EmpiricalRandomSequence &&other) :
		data_set(std::move(other.data_set)), n(std::move(other.n)) {
	rng = std::move(other.rng);
} // end EmpiricalRandomSequence move constructor

// BernoulliSequence(p) constructor
BernoulliSequence::BernoulliSequence(double _p) :
		p(_p) {
	// get generator type from env var
	gsl_rng_env_setup();
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(
			gsl_rng_alloc(gsl_rng_default), GarbageGenerator());
} // end BernoulliSequence constructor

// end BernoulliSequence copy constructor
BernoulliSequence::BernoulliSequence(const BernoulliSequence &other) :
		p(other.p) {
	gsl_rng * copy = gsl_rng_clone(other.rng.get());
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(copy, GarbageGenerator());
} // end BernoulliSequence copy constructor

// BernoulliSequence move constructor
BernoulliSequence::BernoulliSequence(BernoulliSequence &&other) :
		p(std::move(other.p)) {
	rng = std::move(other.rng);
} // end BernoulliSequence move constructor

// BinomialSequence(p, trial) constructor
BinomialSequence::BinomialSequence(double _p, size_t _trial) :
		p(_p), trial(_trial) {
	// get generator type from env var
	gsl_rng_env_setup();
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(
			gsl_rng_alloc(gsl_rng_default), GarbageGenerator());
} // end BinomialSequence(p, trial) constructor

// end BinomialSequence copy constructor
BinomialSequence::BinomialSequence(const BinomialSequence &other) :
		p(other.p), trial(other.trial) {
	gsl_rng * copy = gsl_rng_clone(other.rng.get());
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(copy, GarbageGenerator());
} // end BinomialSequence copy constructor

// BinomialSequence move constructor
BinomialSequence::BinomialSequence(BinomialSequence &&other) :
		p(std::move(other.p)), trial(std::move(other.trial)) {
	rng = std::move(other.rng);
} // end BinomialSequence move constructor

// PoissonSequence(mu) constructor
PoissonSequence::PoissonSequence(double _mu) :
		mu(_mu) {
	// get generator type from env var
	gsl_rng_env_setup();
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(
			gsl_rng_alloc(gsl_rng_default), GarbageGenerator());
} // end PoissonSequence(mu) constructor

// end PoissonSequence copy constructor
PoissonSequence::PoissonSequence(const PoissonSequence &other) :
		mu(other.mu) {
	gsl_rng * copy = gsl_rng_clone(other.rng.get());
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(copy, GarbageGenerator());
} // end PoissonSequence copy constructor

// PoissonSequence move constructor
PoissonSequence::PoissonSequence(PoissonSequence &&other) :
		mu(std::move(other.mu)) {
	rng = std::move(other.rng);
} // end PoissonSequence move constructor

// STNRealSequence(mu, sigma) constructor
STNRealSequence::STNRealSequence(double mean, double stddev) :
		mu(mean), sigma(stddev) {
	// get generator type from env var
	gsl_rng_env_setup();
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(
			gsl_rng_alloc(gsl_rng_default), GarbageGenerator());
} // end STNRealSequence(mu, sigma) constructor

// end STNRealSequence copy constructor
STNRealSequence::STNRealSequence(const STNRealSequence &other) :
		mu(other.mu), sigma(other.sigma) {
	gsl_rng * copy = gsl_rng_clone(other.rng.get());
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(copy, GarbageGenerator());
} // end STNRealSequence copy constructor

// STNRealSequence move constructor
STNRealSequence::STNRealSequence(STNRealSequence &&other) :
		mu(std::move(other.mu)), sigma(std::move(other.sigma)) {
	rng = std::move(other.rng);
} // end STNRealSequence move constructor

// IIDSTNRealVectorSequence(n, mu, sigma)  constructor
IIDSTNRealVectorSequence::IIDSTNRealVectorSequence(size_t n, double mean,
		double stddev) :
		range_size(n), range_dim(1), mu(mean), sigma(stddev) {
	// initialize generators
	gsl_rng_env_setup();
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(
			gsl_rng_alloc(gsl_rng_default), GarbageGenerator());
} // end IIDSTNRealVectorSequence constructor

// IIDSTNRealVectorSequence copy constructor
IIDSTNRealVectorSequence::IIDSTNRealVectorSequence(
		const IIDSTNRealVectorSequence &other) :
		range_size(other.range_size), range_dim(other.range_dim), mu(other.mu), sigma(
				other.sigma) {
	gsl_rng * copy = gsl_rng_clone(other.rng.get());
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(copy, GarbageGenerator());
} // end IIDSTNRealVectorSequence copy constructor

// IIDSTNRealVectorSequence move constructor
IIDSTNRealVectorSequence::IIDSTNRealVectorSequence(
		IIDSTNRealVectorSequence &&other) :
		range_size(std::move(other.range_size)), range_dim(
				std::move(other.range_dim)), mu(std::move(other.mu)), sigma(
				std::move(other.sigma)) {
	rng = std::move(other.rng);
} // end IIDSTNRealVectorSequence move constructor
