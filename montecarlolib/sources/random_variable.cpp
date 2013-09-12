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
#include <initializer_list>				// std::initializer_list
using namespace std;

#include "../includes/range_type.hpp"
#include "../includes/random_variable.hpp"

namespace randomvariable { //

// EmpiricalRandomSequence(data&) constructor
EmpiricalRandomSequence::EmpiricalRandomSequence(const vector<double>& data) :
		data_set(data), n(data.size() - 1) {
	// get generator type from env var
	gsl_rng_env_setup();
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(
			gsl_rng_alloc(gsl_rng_default), GarbageGenerator());
} // end EmpiricalRandomSequence constructor

// EmpiricalRandomSequence(data&&) constructor
EmpiricalRandomSequence::EmpiricalRandomSequence(vector<double> && data) :
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

// NormalSequence(0, sigma) constructor
NormalSequence::NormalSequence(double stddev) :
		sigma(stddev) {
	// get generator type from env var
	gsl_rng_env_setup();
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(
			gsl_rng_alloc(gsl_rng_default), GarbageGenerator());
} // end NormalSequence(0, sigma) constructor

// end NormalSequence copy constructor
NormalSequence::NormalSequence(const NormalSequence &other) :
		sigma(other.sigma) {
	gsl_rng * copy = gsl_rng_clone(other.rng.get());
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(copy, GarbageGenerator());
} // end NormalSequence copy constructor

// Gaussian vector
randomvariable::GaussianRealVectorSequence::GaussianRealVectorSequence(
		const std::vector<double>& _sigma) :
		sigma(_sigma) {
	// get generator type from env var
	gsl_rng_env_setup();
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(
			gsl_rng_alloc(gsl_rng_default), GarbageGenerator());

} // end Gaussian vector

// copy gaussian vector sequence
randomvariable::GaussianRealVectorSequence::GaussianRealVectorSequence(
		const GaussianRealVectorSequence &other) :
		sigma(other.sigma) {
	gsl_rng * copy = gsl_rng_clone(other.rng.get());
	rng = std::unique_ptr<gsl_rng, GarbageGenerator>(copy, GarbageGenerator());
} // end copy Gaussian vector

} // end namespace randomvariable
