/*
 * Copyright (C) 2013  Abel Walga
 *
 *  montecarlo.cpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  montecarlo.cpp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 31 août 2013
 *      Author: Abel Walga
 */

#include <range_type.hpp>
#include <random_variable.hpp>
#include <stochastic_process.hpp>

#include <iostream>
#include <vector>
#include <random>
#include <functional>

using namespace std;

using namespace randomvariable;
using namespace randomprocess;
using namespace montecarlo;
using namespace linear_algebra;

class STN {
public:
	typedef double RangeType;
	typedef double ScalarType;
	typedef normal_distribution<RangeType> Dist;

	default_random_engine generator;
	Dist::param_type resultType;
	Dist stn;

	function<double(void)> stn_gen = std::bind(stn, generator);

	double operator()() {
		return stn_gen();
	}
};
void test_random_object() {
	STN stn;
	std::vector<double> data;
	for (int i = 0; i < 200; i++) {
		data.push_back(stn());
	}
	EmpiricalRandomVariable empVar(std::move(data));
	cout << "Empirical Distribution mean" << endl;
	cout << empVar.expectation(50000) << endl;
	cout << "Empirical Distribtion variance" << endl;
	cout << empVar.variance(50000) << endl;

	STNRealRandomVariable stnVar(1.0);
	cout << "STN Distribution mean" << endl;
	cout << stnVar.expectation(100) << endl;
	cout << "STN Distribution mean" << endl;
	cout << stnVar.variance(100) << endl;

	GaussianRandomVector iidvector(std::vector<double>( { 10.0, 0.0, 1.0 }));
	auto vecMean = iidvector.expectation(10000);
	cout << "Mean vector" << endl;
	for (int i = 0; i < 10; i++) {
		cout << vecMean[i] << endl;
	}
	cout << "Variance vector" << endl;
	auto vecVar = iidvector.variance(10000);
	for (int i = 0; i < 10; i++) {
		cout << vecVar[i] << endl;
	}
	cout << "Covariance i,j" << endl;
	cout << iidvector.covariance(0, 1, 10000) << endl;
	cout << "Correlation i,j squared" << endl;
	cout << iidvector.correlation(0, 1, 10000) << endl;

	typedef Matrix<double, UpperTriangular> UTRMatrix;
	UTRMatrix cov(10, 10);
	iidvector.covarianceMatrix(10000, cov);
	cout << "Coveriance Matrix" << endl;
	for (size_t i = 0; i < cov.size_i(); i++) {
		for (size_t j = 0; j < cov.size_j(); j++) {
			cout << cov(i, j) << "\t";
		}
		cout << endl;
	}
	cout << endl;
}
// first order stochastic X_t_1 = f(X_t)
struct DeterministicWalk {
	typedef long RangeType;
	typedef long ScalarType;

	inline RangeType operator()() {
		return 1;
	}
	// time step function keep it up
	inline RangeType operator()(const RangeType &X_t, const size_t &t_1,
			const double &h) {
		return X_t + 1;
	}
	// stopping time keep walking till the cliff
	inline int operator()(const size_t &t, const RangeType &X_t,
			const RangeType &X_0) const {
		return 0;
	}
};
// second order stochastic X_t_2 = f(X_t_1, X_t)
struct SequenceOrdreDeuxHomogene {
	typedef long RangeType;
	typedef long ScalarType;

	inline RangeType operator()() {
		return 1;
	}
	// time step function keep it up
	inline RangeType operator()(const RangeType &X_t_1, const RangeType &X_t,
			const size_t &t_1, const double &h) {
		return X_t_1 + X_t;
	}
	// stopping time keep walking till the cliff
	inline int operator()(const size_t &t, const RangeType &X_t_1,
			const RangeType &X_t, const RangeType &X_1,
			const RangeType &X_0) const {
		return 0;
	}
};
// third order stochastic X_t_3 = f(X_t_2, X_t_1, X_t)
struct SequenceOrdreTroisHomogene {
	typedef long RangeType;
	typedef long ScalarType;

	inline RangeType operator()() {
		return 1;
	}
	// time step function keep it up
	inline RangeType operator()(const array<RangeType, 3> &X_t,
			const size_t &t_1, const double &h) {
		return X_t[0] + X_t[1] + X_t[2];
	}
	// stopping time keep walking till the cliff
	inline int operator()(const size_t &t, const array<RangeType, 3> &X_t,
			const array<RangeType, 3> &X_0) const {
		return 0;
	}
};
// one dimension walk

typedef RandomWalkGenerator<DeterministicWalk> DeterministicWalkGenerator;
typedef SecondOrderStochasticDriver<SequenceOrdreDeuxHomogene> DeterministicSequenceOrdreDeux;
typedef StochasticDriver<SequenceOrdreTroisHomogene, 3> DeterministicSequenceOrdreTrois;

#include <array>
// path functional test
template<long Step>
struct HeavySide {

	typedef long RangeType;
	typedef long ScalarType;

	inline long operator()(const vector<RangeType> &h,
			const size_t &stoppingTime) {
		return Step <= h[stoppingTime - 1] ? 1 : 0;
	}
};

void testProcessOrderDeux() {
	// test stochastic process
	array<long, 3> init { 0, 1, 2 };
	StochasticProcess<DeterministicSequenceOrdreTrois,
			StoppingTimeFunctionalPath> process(SequenceOrdreTroisHomogene(),
			init, 10, 1.0);
	auto var = process();
	auto mean = var.expectation(1000);
	cout << mean << endl;
}

void testGaussianRandomVector() {
	GaussianRandomVector vector1(std::vector<double>( { 1.0 }));
	GaussianRandomVector vector2(std::move(vector1));
	auto meanVec = vector2.expectation(5000);
	for (size_t i = 0; i < meanVec.size(); i++) {
		cout << meanVec[i] << endl;
	}
	auto variance = vector2.variance(5000);
	for (size_t i = 0; i < variance.size(); i++) {
		cout << variance[i] << endl;
	}
}
// random walk
const double X_0 = 1.0;
const double dt = 1.0;
const size_t T = 100;

//  one dimension random walk X_t_1 = f(X_t)
struct OneDimensionRandomWalkGenerator {
	// rvar range type
	typedef double RangeType;
	typedef double ScalarType;
	// rng
	NormalSequence stn;

	OneDimensionRandomWalkGenerator() :
			stn(std::sqrt(dt)) {
	}
	// time step function keep it up
	RangeType operator()(const RangeType &X_t, const size_t &t_1,
			const double &h) {
		return X_t + stn();// X_t+sqrt(dt)*STN();
	}
	// stopping time keep walking till the cliff
	int operator()(const size_t &t, const RangeType &X_t,
			const RangeType &X_0) const {
		return 0;
	}
};
// one dimension random walk functional defined for a simulation time of 100 and starting point 1.0
struct OneDimensionRandomWalkFunctional {

	typedef RealRowVector RangeType;
	typedef double ScalarType;

	size_t simulationTime;
	double meanValue;

	OneDimensionRandomWalkFunctional(): simulationTime(T), meanValue(X_0){}

	/**
	 * @param path
	 * @param stoppingTime
	 * @return randomvalue
	 */
	RangeType operator()(const vector<double> &path, const size_t &stoppingTime) {
		// TODO check path size is not bigger than simul time.
		RangeType h(simulationTime-1);
		for(size_t i = 0; i<path.size();i++){
			h[i] = path[i];
		}
		for(size_t i = path.size(); i<simulationTime-1; i++){
			h[i] = meanValue;
		}
		return h;
	}
};
typedef StochasticDriver<OneDimensionRandomWalkGenerator, 1> RandomWalkDriver;
typedef StochasticProcess<RandomWalkDriver, OneDimensionRandomWalkFunctional> RandomWalkStochasticProcess;
int main() {
	RandomWalkStochasticProcess randomWalkProcess(OneDimensionRandomWalkGenerator(), X_0, T, dt);
	auto randomWalkVar = randomWalkProcess();
	auto meanVector = randomWalkVar.expectation(1000);
	for(size_t i = 0; i<meanVector.size(); i++){
		cout << i << ":" << meanVector[i] << endl;
	}
}
