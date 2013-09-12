/*
 * Copyright (C) 2013  Abel Walga
 *
 *  stochastic_process.hpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  stochastic_process.hpp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 26 août 2013
 *      Author: Abel Walga
 */
#ifndef STOCHASTIC_PROCESS_HPP_
#define STOCHASTIC_PROCESS_HPP_

#include "random_object.hpp"

#include <vector>					// std::vector
#include <array>					// std::array

namespace randomprocess {

/*****************************************************************************************************************
 * 						Stochastic process, path sample generator, path function, stopping time definition
 *****************************************************************************************************************/
// forward declaration of a path sample generator used as the sequence generator for a conditioned random variable
template<typename StochasticDriver, typename Functional>
class PathSequence;

/******************************************************************************************************************
 *
 ******************************************************************************************************************/

/**
 * <p>
 * 	  Define a stochastic process indexed by a discrete time set t_0, t_1, ..., t_n, for a given state space <code>X_t_Type</code>
 * 	  and stochastic driver <code>StochasticDriver</code> : X(t,N) = (X_0, X_1, X_2, X_3, ...)
 * </p>
 *
 * <p>
 *	  The stochastic process yields a random variable H = f(X(t)) functional of a conditioned path X(t)=(X_0 ,..., X_i, ..., X_i+k, X_i+k+1,...,X_t) where:
 *
 * 			i.	i>=0 (X_0,..., X_i) is the initial deterministic path,
 * 			ii.	(X_i+1, ...,X_i+k)\ k>=1  or {}\ k=0  the random conditioning path,
 * 			iii. (X_i+k+1, ..., X_t) the random path sample conditioned on (X_0,...,X_i+k)
 * 			iV.  T>t>=i+k>=0 the stopping time - 1 unit of time.
 *	</p>
 *	<p>
 *	 Example: a sample path of stopping time T with a deterministic initial condition X_0:
 * 	  		i.  conditioned at time k = 0 is given by:
 * 	  		X(t) = (X_0, ..., X_T) where X_0 is deterministic for any sample of the conditioned random variable.
 * 	  		ii. conditioned at time k = 3 is given by:
 * 	  		X(t) = (X_0, X_1, X_2, X_3,..., X_T) where X_0, X_1, X_2, X_3 are deterministic for any sample of the conditioned random path samples
 * 	</p>
 */
template<typename StochasticDriver, typename Functional>
class StochasticProcess {

public:
	/**
	 * \def X_i_Generator the stochastic driver used to advance the path at each time step.
	 */
	typedef StochasticDriver X_t_Generator;
	/**
	 * \def SequenceGenerator used to sample a path up to stopping time.
	 */
	typedef PathSequence<X_t_Generator, Functional> PathSampleGenerator;

	// type def path range : PathRangeType is X_i structure and if dimension of X_i =1 then PathRangeType and PathScalarType are of same range.
	/**
	 * \def X_i RangeType.
	 */
	typedef typename X_t_Generator::X_t_Type PathElementRangeType;
	/**
	 * \def X_i ScalarType same as RangeType if X_i is one dimensional.
	 */
	typedef typename X_t_Generator::X_t_ScalarType PathElementScalarType;

	/**
	 * \def H=f(X(t)) RangeType.
	 */
	typedef typename PathSampleGenerator::RangeType RandomVarRangeType;

	/**
	 * \def H=f(X(t)) ScalarType same as RangeType if H is one dimensional.
	 */
	typedef typename PathSampleGenerator::ScalarType RandomVarScalarType;

	/**
	 * <p>Arguments forwarding constructor</p>
	 *
	 * @param args
	 */
	template<typename ...Args>
	StochasticProcess(Args&&... args) :
			generator(X_t_Generator(std::forward<Args>(args)...)) {
	}

	/**
	 * <p>Constructor</p>
	 *
	 * @param _generator
	 */
	StochasticProcess(const X_t_Generator& _generator) :
			generator(_generator) {
	}

	/**
	 * <p>get X(t) the random variable conditioned at a given time <code>conditionedAt</code></p>
	 *
	 * @param conditionedAt
	 * @return RandomObject conditioned at conditionedAt
	 */
	montecarlo::RandomObject<PathSampleGenerator, RandomVarRangeType, RandomVarScalarType> operator()(
			size_t conditionedAt) {
		// return random variable conditioned at time t.
		return montecarlo::RandomObject<PathSampleGenerator, RandomVarRangeType,
				RandomVarScalarType>(
				generator.getConditionedAtDriver(conditionedAt)); // at i+k
	}

	/**
	 *
	 */
	/**
	 * <p>get X(0) the random variable conditioned on the initial deterministic condition</p>
	 *
	 * @return  RandomObject conditioned on the deterministic initial condition (X_0, ..., X_i)
	 */
	montecarlo::RandomObject<PathSampleGenerator, RandomVarRangeType, RandomVarScalarType> operator()() {
		return (*this)(0);
	} // end
private:
	X_t_Generator generator;
};

// path generator definition
/**
 * <p>
 *	  Random path sample generator.
 * </p>
 * <p>
 * 	  The path sequence generates a stochastic path (X_0,...,X_t) driven by the <code>StochasticDriver</code> where t<T is the stopping time-1
 * 	  and evaluates the <code>Functional</code> of this path.
 * </p>
 */
template<typename StochasticDriver, typename Functional>
class PathSequence {

public:

	// path builder range type
	typedef typename Functional::RangeType RangeType;
	typedef typename Functional::ScalarType ScalarType;
	typedef typename StochasticDriver::X_t_Type PathRangeType;

	/**
	 * <p>Constructor</p>
	 *
	 * @param _generator pathGenerator
	 */
	PathSequence(const StochasticDriver &_generator) :
			generator(_generator), function(Functional()) {
	}
	/**
	 * <p>Get path end point.</p>
	 */
	RangeType operator()() {
		std::vector<PathRangeType> h;
		// get path
		generator.conditionAt(0); // at i
		while (!generator.isStoppingTime()) {
			h.push_back(generator());
		}
		// return H
		return function(std::move(h), generator.getCurrentTime());
	} // end
private:
	StochasticDriver generator;
	Functional function;
};

// specialization when the functional returns the last element of the path
struct StoppingTimeFunctionalPath {
};

template<typename StochasticDriver>
class PathSequence<StochasticDriver, StoppingTimeFunctionalPath> {

public:

	// path builder range type
	typedef typename StochasticDriver::X_t_Type RangeType;
	typedef typename StochasticDriver::X_t_Type ScalarType;

	/**
	 * <p>Constructor</p>
	 *
	 * @param _generator
	 */
	PathSequence(const StochasticDriver &_generator) :
			generator(_generator) {
	}

	/**
	 * <p>Get path end point.</p>
	 */
	RangeType operator()() {
		RangeType value = generator.conditionAt(0); // at i
		while (!generator.isStoppingTime()) {
			value = generator();
		}
		return value;
	} // end
private:
	StochasticDriver generator;
};

/*****************************************************************************************************************
 * 							path driver interface
 *****************************************************************************************************************/
/**
 * <p>
 * 		A stochastic path driver is used to generate the next value for a family of X(t): (X_0 ,..., X_i, ..., X_i+k, X_i+k+1, ..., X_t) where T>t>=k>=i>=0
 * 		where (X_0, ..., X_i) is a deterministic conditioning path prior to simulation, (X_i+1, ...,X_i+k) is a random conditioning path
 * 		and (X_i+k+1, ..., ,X_t-1) is a sample path where T>t-1 \ T being the stopping time. i+1 is the offSet of the generator,
 * 		k is the conditioning time. @see StochasticProcess::operator(const size_t&), T is the stopping time for a given sample path.
 * </p>
 */
/**
 * <p>
 * 		Define a constant time (h = t_i_1 -t_i = t_j_1 - t_j) n-th order stochastic driver :
 *
 * 						X_t_1= f(X_t, ..., X_t-n)  t_1 = t+h.
 * </p>
 * <p>
 * 		Sequence must define the return value type <code>X_t_Type</code> and <code>X_t_ScalarType</code>, and the following operators:
 *
 * 		<p>
 * 		i. <code>X_t_Type Sequence::operator(const X_t_Type &X_t_1, const X_t_Type &X_t, const size_t &t_2, const double &h)</code>
 * 			   where X_t_1 and X_t are the values at time t_1 and t_2 respectively, t_2 current computing time, h time step.
 * 		</p>
 * 		<p>
 * 		ii.<code>int Sequence::operator(const size_t& t_2, const X_t_Type &X_t_1, const X_t_Type &X_t, const X_t_Type &X_1, const X_t_Type &X_0) const</code>
 * 			   where t_2 is current time, X_t_1, X_t previous time value, and X_0 initial condition
 * 		</p>
 * 	</p>
 */
template<typename Sequence, size_t Order>
class StochasticDriver {
public:
	/**
	 * \def X_t_Type path coordinate type.
	 */
	typedef typename Sequence::RangeType X_t_Type;
	/**
	 * \def X_t_ScalarType path scalar type if X_t_Type is a scalar then X_t_Type and X_t_ScalarType are of same type.
	 */
	typedef typename Sequence::ScalarType X_t_ScalarType; // if range type dimension is greater than 1.

	/**
	 * <p>Constructor</p>
	 *
	 * @param _seq Sequence
	 * @param _X_i i initial conditions
	 * @param T size_t simulation horizon
	 * @param _dt double time step
	 */
	StochasticDriver(Sequence _seq, std::array<X_t_Type, Order> _X_i, size_t T,
			double _dt) :
			seq(_seq), h(_dt), stoppingTime(T), currentTime(Order), X_i(_X_i), X_t(
					_X_i) {
	}

	/**
	 * <p>generate path up to current time</p>
	 *
	 * @return RangeType pathValue
	 */
	X_t_Type operator()() {
		X_t[currentTime % Order] = seq(X_t, currentTime, h);
		currentTime++;
		return X_t[(currentTime - 1) % Order];
	}
	/**
	 * <p>condition this path generator to t=(i_1 + k)</p>
	 *
	 * @param const size_t& k
	 * @param RangeType X_i_k
	 */

	X_t_Type conditionAt(const size_t &k) {
		size_t conditionJustBefore = Order + k; // i+1+k
		size_t conditioningTime = Order; // i+1
		// k=0;
		X_t = X_i; // copy same as std::copy
		for (; conditioningTime < conditionJustBefore; conditioningTime++) { // k>=1 && i+1 < i+k+1
			X_t[conditioningTime % Order] = seq(X_t, conditioningTime, h);
		}
		currentTime = conditionJustBefore;
		return X_t[(currentTime - 1) % Order]; //(i+k)%(i+1)
	}
	/**
	 * <p>stopping time condition</p>
	 *
	 * @param const size_t& any time t
	 */
	int isStoppingTime() const {
		return (currentTime >= stoppingTime || seq(currentTime, X_t, X_i));
	}

	/**
	 * <p>get deterministic/initial condition path time + 1</p>
	 *
	 * @return size_t offSet
	 */
	size_t getOffSet() const {
		return Order;
	}

	/**
	 * <p>get path time + 1</p>
	 *
	 * @return size_t t
	 */
	size_t getCurrentTime() const {
		return currentTime;
	}

	/**
	 * <p>get a stochastic driver where the initial condition is the conditioning path of <code>this</code></p>
	 *
	 * @param k const size_t& conditioning time
	 * @return driver StochasticDriver<Sequence, Order>
	 */
	StochasticDriver<Sequence, Order> getConditionedAtDriver(const size_t &k) {
		this->conditionAt(k);
		std::array<X_t_Type, Order> init;
		// reorder time
		for (size_t j = Order, i = currentTime - 1; j > 0; j--, i--) {
			init[j-1] = X_t[i % Order];
		}
		return StochasticDriver<Sequence, Order>(seq, init, stoppingTime, h);
	}
private:
	/**
	 *  <p>hands of god if stochastic or your will if deterministic</p>
	 */
	Sequence seq;

	/**
	 * <p>time increment t_1 = t_0 + h
	 */
	double h;

	/**
	 * for a given  X(t)=(X_0 ,..., X_i, ..., X_i+k, X_i+k+1, ..., X_t,...,X_T) i>=order-1
	 *
	 * offSet = order, stoppingTime>=0, conditionningTime = order-1+k and currentTime>= order+k
	 */
	const size_t stoppingTime; // offset = Order
	size_t currentTime;

	/**
	 * <p>last n=order coordinates</p>
	 */
	const std::array<X_t_Type, Order> X_i;
	/**
	 * <p>last n=order coordinates</p>
	 */
	std::array<X_t_Type, Order> X_t;
};
/*****************************************************************************************************************
 * 							some path drivers
 *****************************************************************************************************************/
/**
 * \def constant for a first order stochastic generator
 */
const size_t FIRST_ORDER = 1;
/**
 * \def constant for a second order stochastic generator
 */
const size_t SECOND_ORDER = 2;

/**
 * <p>
 * 		Define a constant time (h = t_i_1 -t_i = t_j_1 - t_j) first order stochastic path driver:
 *
 * 						X_t_1 = f(X_t) t_1 = t+h.
 * 	</p>
 * 	<p>
 * 		Sequence must define the return value type <code>X_t_Type</code> and <code>X_t_Type</code>, and the following operators:
 *
 * 		<p>
 * 		i. <code>X_t_Type Sequence::operator(const X_t_Type &X_t, const size_t &t_1, const double &h)</code>
 * 			  where X_t is the current value, t_1 current computing time, h time step.
 * 		</p>
 * 		<p>
 * 		ii.<code>int Sequence::operator(const size_t& t_1, const X_t_Type &X_t, const X_t_Type &X_0) const</code>
 * 			  where t_1 is current time, X_t previous time value, and X_0 initial condition
 * 		</p>
 * 	</p>
 */
template<typename Sequence>
class StochasticDriver<Sequence, FIRST_ORDER> {
public:

	/**
	 * \def X_t_Type path coordinate type.
	 */
	typedef typename Sequence::RangeType X_t_Type;
	/**
	 * \def X_t_ScalarType path scalar type if X_t_Type is a scalar then X_t_Type and X_t_ScalarType are of same type.
	 */
	typedef typename Sequence::ScalarType X_t_ScalarType; // if range type dimension is greater than 1.

	/**
	 * <p>Constructor</p>
	 *
	 * @param Sequence _seq stochastic driver
	 * @param RangeType _X_i initial condtion
	 * @param size_t T simulation horizon
	 * @param double _dt time increment
	 */
	StochasticDriver(Sequence _seq, X_t_Type _X_i, size_t T, double _dt) :
			seq(_seq), h(_dt), offSet(1), stoppingTime(T), currentTime(1), X_i(
					_X_i), X_t(_X_i) {
	}

	/**
	 * <p>generate path up to current time t</p>
	 *
	 * @return RangeType pathValue
	 */
	X_t_Type operator()() {
		X_t = seq(X_t, currentTime, h);
		currentTime++;
		return X_t;
	}

	/**
	 * <p>condition this path generator to t=(i + k)</p>
	 *
	 * @param const size_t& k
	 * @param RangeType X_i_k
	 */
	X_t_Type conditionAt(const size_t& k) {
		size_t conditionJustBefore = offSet + k; // i+k+1
		size_t conditioningTime = offSet; // i+1
		X_t = X_i; // k=0
		for (; conditioningTime < conditionJustBefore; conditioningTime++) { // k>=1 && i+k < i+k+1
			X_t = seq(X_t, conditioningTime, h);
		}
		currentTime = conditionJustBefore;
		return X_t;
	}
	/**
	 * <p>stopping time condition</p>
	 *
	 * @param const size_t& any time t
	 */
	int isStoppingTime() const {
		return (currentTime >= stoppingTime || seq(currentTime, X_t, X_i));
	}

	/**
	 * <p>get deterministic/initial condition path time + 1</p>
	 *
	 * @return size_t offSet
	 */
	size_t getOffSet() const {
		return offSet;
	}

	/**
	 * <p>get path time + 1. potentially a stopping time if stopping time reached then return the stopping time</p>
	 *
	 * @return size_t t
	 */
	size_t getCurrentTime() const {
		return currentTime;
	}

	/**
	 * <p>get a stochastic driver where the initial condition is the conditioning path of <code>this</code></p>
	 *
	 * @param k const size_t& conditioning time
	 * @return driver StochasticDriver<Sequence, Order>
	 */
	StochasticDriver<Sequence, 1> getConditionedAtDriver(const size_t &k) {
		this->conditionAt(k);
		return StochasticDriver<Sequence, 1>(seq, X_t, stoppingTime, h);
	}
private:
	/**
	 *  <p>hands of god if stochastic or your will if deterministic</p>
	 */
	Sequence seq;
	/**
	 * <p>time increment t_1 = t_0 + h
	 */
	double h;
	/**
	 * for a given  X(t)=(X_0 ,..., X_i, ..., X_i+k, X_i+k+1, ..., X_t,...,X_T) i>=0
	 *
	 * invariant : offSet = i+1, stoppingTime>=0, conditionningTime = i+k and currentTime>=i+1+k
	 */
	const size_t offSet, stoppingTime;
	size_t currentTime;

	/**
	 * <p>for a random walk we only need the last step.</p>
	 */
	const X_t_Type X_i;
	X_t_Type X_t;
};

/**
 * \def First order path driver where the amplitude and the type of walk is determined by <code>Sequence</code>.
 */
template<typename Sequence> using RandomWalkGenerator = StochasticDriver<Sequence,FIRST_ORDER>;
// end

/**
 * <p>
 * 		Define a constant time (h = t_i_1 -t_i = t_j_1 - t_j) second order stochastic path driver :
 *
 * 						X_t_2 = f(X_t_1, X_t) t_1 = t+h, t_2 = t+2*h.
 * </p>
 * <p>
 * 		Sequence must define the return value type <code>X_t_Type</code> and <code>X_t_ScalarType</code>, and the following operators:
 *
 * 		<p>
 * 		i. <code>RangeType Sequence::operator(const X_t_Type &X_t_1, const X_t_Type &X_t, const size_t &t_2, const double &h)</code>
 * 			  where X_t_1 and X_t are the values at time t_1 and t_2 respectively, t_2 current computing time, h time step.
 * 		</p>
 * 		<p>
 * 		ii.<code>int Sequence::operator(const size_t& t_2, const X_t_Type &X_t_1, const X_t_Type &X_t, const X_t_Type &X_1, const X_t_Type &X_0) const</code>
 * 			  where t_2 is current time, X_t_1, X_t previous time value, and X_0 initial condition.
 * 		</p>
 * 	</p>
 */
template<typename Sequence>
class StochasticDriver<Sequence, SECOND_ORDER> {
public:
	/**
	 * \def X_t_Type path coordinate type.
	 */
	typedef typename Sequence::RangeType X_t_Type;
	/**
	 * \def X_t_ScalarType path scalar type if X_t_Type is a scalar then X_t_Type and X_t_ScalarType are of same type.
	 */
	typedef typename Sequence::ScalarType X_t_ScalarType; // if range type dimension is greater than 1.

	/**
	 * <p>Constructor</p>
	 * @param _seq source of stochastic
	 * @param _X_i 1st initial condition
	 * @param _X_i_1 2nd initial condition
	 * @param T simulation time
	 * @param _dt time step
	 */
	StochasticDriver(Sequence _seq, X_t_Type _X_i, X_t_Type _X_i_1, size_t T,
			double _dt) :
			seq(_seq), h(_dt), offSet(2), stoppingTime(T), currentTime(2), X_i(
					_X_i), X_i_1(_X_i_1), X_t(_X_i), X_t_1(_X_i_1) {
	}
	/**
	 * <p>generate path up to current time t_2</p>
	 *
	 * @return RangeType pathValue
	 */
	X_t_Type operator()() {
		auto tmp = X_t_1;
		X_t_1 = seq(X_t_1, X_t, currentTime, h);
		X_t = std::move(tmp);
		currentTime++;
		return X_t_1;
	}

	/**
	 * <p>condition this path generator to t=(i_1 + k)</p>
	 *
	 * @param const size_t& k
	 * @param RangeType X_i_k
	 */
	X_t_Type conditionAt(const size_t& k) {
		size_t conditionJustBefore = offSet + k; // i+k+2
		size_t conditioningTime = offSet; // i+2
		X_t = X_i;
		X_t_1 = X_i_1; // k=0
		X_t_Type tmp;
		for (; conditioningTime < conditionJustBefore; conditioningTime++) { // k>=1 && i+2 < i+k+2
			tmp = X_t_1;
			X_t_1 = seq(X_t_1, X_t, conditioningTime, h); // one step forward for t_1
			std::swap(X_t, tmp); // one step forward for t
		}
		currentTime = conditionJustBefore;
		return X_t_1;
	}
	/**
	 * <p>stopping time condition</p>
	 *
	 * @param const size_t& any time t
	 */
	int isStoppingTime() const {
		return (currentTime >= stoppingTime
				|| seq(currentTime, X_t_1, X_t, X_i_1, X_i));
	}

	/**
	 * <p>get deterministic/initial condition path time + 1</p>
	 *
	 * @return size_t offSet
	 */
	size_t getOffSet() const {
		return offSet;
	}

	/**
	 * <p>get path time + 1</p>
	 *
	 * @return size_t t
	 */
	size_t getCurrentTime() const {
		return currentTime;
	}

	/**
	 * <p>get a stochastic driver where the initial condition is the conditioning path of <code>this</code></p>
	 *
	 * @param k const size_t& conditioning time
	 * @return driver StochasticDriver<Sequence, Order>
	 */
	StochasticDriver<Sequence, 2> getConditionedAtDriver(const size_t &k) {
		this->conditionAt(k);
		return StochasticDriver(seq, X_t, X_t_1, stoppingTime, h);
	}
private:

	/**
	 *  <p>hands of god if stochastic or your will if deterministic</p>
	 */
	Sequence seq;
	/**
	 * <p>time increment t_1 = t_0 + h
	 */
	double h;
	/**
	 * for a given  X(t)=(X_0 ,..., X_i, ..., X_i+k, X_i+k+1, ..., X_t,...,X_T) i>0
	 *
	 * invariant : offSet = i+2, stoppingTime>=0, conditionningTime = i+1+k and currentTime>=i+k+2
	 */
	const size_t offSet, stoppingTime;
	size_t currentTime;

	/**
	 * <p>for a random walk we only need the last step.</p>
	 */
	const X_t_Type X_i, X_i_1; // X_0, X_1
	X_t_Type X_t, X_t_1;
};

/**
 * \def Second order path driver where X_t_2 = f(X_t1, X_t) and f is defined by <code>Sequence</code>.
 */
template<typename Sequence> using SecondOrderStochasticDriver = StochasticDriver<Sequence, SECOND_ORDER>;

} // end montecarlo

#endif /* STOCHASTIC_PROCESS_HPP_ */
