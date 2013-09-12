/*
 * Copyright (C) 2013  Abel Walga
 *
 *  unit_test_stochastic_process.cpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  unit_test_stochastic_process.cpp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 29 août 2013
 *      Author: Abel Walga
 */

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#define BOOST_TEST_MODULE ThreadSafeContainerRegressionTest
#endif // STAND_ALONE

#include <stochastic_process.hpp>
#include <boost/test/unit_test.hpp>

#include <array>

using namespace std;
using namespace randomprocess;

BOOST_AUTO_TEST_SUITE(StochasticProcessTest)

// first order stochastic X_t_1 = f(X_t)
struct DeterministicWalk {
	typedef long RangeType;
	typedef long ScalarType;

	RangeType operator()() {
		return 1;
	}
	// time step function keep it up
	RangeType operator()(const RangeType &X_t, const size_t &t_1,
			const double &h) {
		return X_t + 1;
	}
	// stopping time keep walking till the cliff
	int operator()(const size_t &t, const RangeType &X_t,
			const RangeType &X_0) const {
		return 0;
	}
};
// second order stochastic X_t_2 = f(X_t_1, X_t)
struct SequenceOrdreDeuxHomogene {
	typedef long RangeType;
	typedef long ScalarType;

	RangeType operator()() {
		return 1;
	}
	// time step function keep it up
	RangeType operator()(const RangeType &X_t_1, const RangeType &X_t,
			const size_t &t_1, const double &h) {
		return X_t_1 + X_t;
	}
	// stopping time keep walking till the cliff
	int operator()(const size_t &t, const RangeType &X_t_1,
			const RangeType &X_t, const RangeType &X_1,
			const RangeType &X_0) const {
		return 0;
	}
};
// third order stochastic X_t_3 = f(X_t_2, X_t_1, X_t)
struct SequenceOrdreTroisHomogene {
	typedef long RangeType;
	typedef long ScalarType;

	RangeType operator()() {
		return 1;
	}
	// time step function keep it up
	RangeType operator()(const std::array<RangeType, 3> &X_t,
			const size_t &t_1, const double &h) {
		return X_t[0] + X_t[1] + X_t[2];
	}
	// stopping time keep walking till the cliff
	int operator()(const size_t &t, const std::array<RangeType, 3> &X_t,
			const std::array<RangeType, 3> &X_0) const {
		return 0;
	}
};
// one dimension walk

typedef RandomWalkGenerator<DeterministicWalk> DeterministicWalkGenerator;
typedef SecondOrderStochasticDriver<SequenceOrdreDeuxHomogene> DeterministicSequenceOrdreDeux;
typedef StochasticDriver<SequenceOrdreTroisHomogene, 3> DeterministicSequenceOrdreTrois;

//// path functional test
template<long Step>
struct HeavySide {

	typedef long RangeType;
	typedef long ScalarType;

	long operator()(const vector<RangeType> &h,
			const size_t &stoppingTime) {
		return Step <= h[h.size() - 1] ? 1 : 0;
	}
};

// first order sequence generator test
BOOST_AUTO_TEST_CASE(RandomWalkGenerator) {
	// test driver
	DeterministicWalkGenerator walk(DeterministicWalk(), 0, 10, 1.0);
	long X_0 = walk.conditionAt(0);
	long X_T = X_0;
	BOOST_CHECK_EQUAL(X_0, 0); // first element is 0;
	BOOST_CHECK_EQUAL(walk.getCurrentTime(), 1); // current time = offSet = 1;
	BOOST_CHECK_EQUAL(walk.getOffSet(), walk.getCurrentTime()); // current time = offSet = 1;
	while (!walk.isStoppingTime()) {
		X_T = walk();
	}
	BOOST_CHECK_EQUAL(walk.getCurrentTime(), 10); // current time = stopping time
	BOOST_CHECK_EQUAL(walk.getOffSet(), 1); // current time = offSet = 1;
	BOOST_CHECK_EQUAL(X_T, 9); // last element is 9;

	// same resut expected on re conditioned
	X_0 = walk.conditionAt(0);
	BOOST_CHECK_EQUAL(X_0, 0); // frist element is 0;
	BOOST_CHECK_EQUAL(walk.getCurrentTime(), 1); // current time = offSet = 1;
	BOOST_CHECK_EQUAL(walk.getOffSet(), walk.getCurrentTime()); // current time = offSet = 1;
	while (!walk.isStoppingTime()) {
		X_T = walk();
	}
	BOOST_CHECK_EQUAL(walk.getCurrentTime(), 10); // current time = stopping time
	BOOST_CHECK_EQUAL(walk.getOffSet(), 1); // current time = offSet = 1;
	BOOST_CHECK_EQUAL(X_T, 9); // last element is 9;

	// test stochastic process
	StochasticProcess<DeterministicWalkGenerator, StoppingTimeFunctionalPath> deterministicWalkProcess(
			DeterministicWalk(), 0, 10, 1.0);
	auto var = deterministicWalkProcess();
	auto mean = var.expectation(1000);
	BOOST_CHECK_EQUAL(mean, 9);

	// test functional stochastic process
	StochasticProcess<DeterministicWalkGenerator, HeavySide<9> > processBiggerThanMe(
			DeterministicWalk(), 0, 10, 1.0);
	BOOST_CHECK_EQUAL(processBiggerThanMe().expectation(1000), 1);

	StochasticProcess<DeterministicWalkGenerator, HeavySide<10> > processShorterThanMe(
			DeterministicWalk(), 0, 10, 1.0);
	BOOST_CHECK_EQUAL(processShorterThanMe().expectation(1000), 0);
}
// second order sequence generator test
BOOST_AUTO_TEST_CASE(SecondOrderSequenceGenerator) {
	DeterministicSequenceOrdreDeux sequence(SequenceOrdreDeuxHomogene(), 0, 1,
			10, 1.0);
	long X_1 = sequence.conditionAt(0);
	long X_T = X_1;
	BOOST_CHECK_EQUAL(X_1, 1); // first element is 1;
	BOOST_CHECK_EQUAL(sequence.getCurrentTime(), 2); // current time = offSet = 2;
	BOOST_CHECK_EQUAL(sequence.getOffSet(), sequence.getCurrentTime()); // current time = offSet = 2;
	while (!sequence.isStoppingTime()) {
		X_T = sequence();
	}
	BOOST_CHECK_EQUAL(sequence.getCurrentTime(), 10); // current time = stopping time
	BOOST_CHECK_EQUAL(sequence.getOffSet(), 2); // current time = offSet = 2;
	BOOST_CHECK_EQUAL(X_T, 34); // last element is 34;

	// same resut expected on re conditioned
	X_1 = sequence.conditionAt(0);
	BOOST_CHECK_EQUAL(X_1, 1); // frist element is 0;
	BOOST_CHECK_EQUAL(sequence.getCurrentTime(), 2); // current time = offSet = 2;
	BOOST_CHECK_EQUAL(sequence.getOffSet(), sequence.getCurrentTime()); // current time = offSet = 2;
	while (!sequence.isStoppingTime()) {
		X_T = sequence();
	}
	BOOST_CHECK_EQUAL(sequence.getCurrentTime(), 10); // current time = stopping time
	BOOST_CHECK_EQUAL(sequence.getOffSet(), 2); // current time = offSet = 2;
	BOOST_CHECK_EQUAL(X_T, 34); // last element is 34;
	// test stochastic process
	StochasticProcess<DeterministicSequenceOrdreDeux, StoppingTimeFunctionalPath> process(
			SequenceOrdreDeuxHomogene(), 0, 1, 10, 1.0);
	auto var = process();
	auto mean = var.expectation(1000);
	BOOST_CHECK_EQUAL(mean, 34);

	// test functional stochastic process
	StochasticProcess<DeterministicSequenceOrdreDeux, HeavySide<34> > processBiggerThanMe(
			SequenceOrdreDeuxHomogene(), 0, 1, 10, 1.0);
	BOOST_CHECK_EQUAL(processBiggerThanMe().expectation(1000), 1);

	StochasticProcess<DeterministicSequenceOrdreDeux, HeavySide<35> > processShorterThanMe(
			SequenceOrdreDeuxHomogene(), 0, 1, 10, 1.0);
	BOOST_CHECK_EQUAL(processShorterThanMe().expectation(1000), 0);
}
// third order sequence generator test
BOOST_AUTO_TEST_CASE(ThirddOrderSequenceGenerator) {
	DeterministicSequenceOrdreTrois sequence(SequenceOrdreTroisHomogene(), { 0,
			1, 2 }, 10, 1.0);
	long X_2 = sequence.conditionAt(0);
	long X_T = X_2;
	BOOST_CHECK_EQUAL(X_2, 2); // first element is 2;
	BOOST_CHECK_EQUAL(sequence.getCurrentTime(), 3); // current time = offSet = 3;
	BOOST_CHECK_EQUAL(sequence.getOffSet(), sequence.getCurrentTime()); // current time = offSet = 3;
	while (!sequence.isStoppingTime()) {
		X_T = sequence();
	}
	BOOST_CHECK_EQUAL(sequence.getCurrentTime(), 10); // current time = stopping time
	BOOST_CHECK_EQUAL(sequence.getOffSet(), 3); // current time = offSet = 3;
	BOOST_CHECK_EQUAL(X_T, 125); // last element is 125;

	// same resut expected on re conditioned
	X_2 = sequence.conditionAt(0);
	BOOST_CHECK_EQUAL(X_2, 2); // first element is 2;
	BOOST_CHECK_EQUAL(sequence.getCurrentTime(), 3); // current time = offSet = 3;
	BOOST_CHECK_EQUAL(sequence.getOffSet(), sequence.getCurrentTime()); // current time = offSet = 3;
	while (!sequence.isStoppingTime()) {
		X_T = sequence();
	}
	BOOST_CHECK_EQUAL(sequence.getCurrentTime(), 10); // current time = stopping time
	BOOST_CHECK_EQUAL(sequence.getOffSet(), 3); // current time = offSet = 3;
	BOOST_CHECK_EQUAL(X_T, 125); // last element is 125;

	// test stochastic process
	std::array<long, 3> init { 0, 1, 2 };
	StochasticProcess<DeterministicSequenceOrdreTrois,
			StoppingTimeFunctionalPath> process(SequenceOrdreTroisHomogene(),
			init, 10, 1.0);
	auto var = process();
	auto mean = var.expectation(1000);
	BOOST_CHECK_EQUAL(mean, 125);

	// test functional stochastic process
	StochasticProcess<DeterministicSequenceOrdreTrois, HeavySide<125> > processBiggerThanMe(
			SequenceOrdreTroisHomogene(), init, 10, 1.0);
	BOOST_CHECK_EQUAL(processBiggerThanMe().expectation(1000), 1);

	StochasticProcess<DeterministicSequenceOrdreTrois, HeavySide<126> > processShorterThanMe(
			SequenceOrdreTroisHomogene(), init, 10, 1.0);
	BOOST_CHECK_EQUAL(processShorterThanMe().expectation(1000), 0);
}
BOOST_AUTO_TEST_SUITE_END()

