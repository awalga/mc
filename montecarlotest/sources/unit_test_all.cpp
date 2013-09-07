/*
 * Copyright (C) 2013  Abel Walga
 *
 *  unit_test_all.cpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  unit_test_all.cpp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 22 août 2013
 *      Author: Abel Walga
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN AllMonteCarloLibTest

#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;

struct AllMonteCarloLibTest {
	AllMonteCarloLibTest() {
		unit_test_log.set_threshold_level(
				boost::unit_test::log_level::log_successful_tests);
	}
	~AllMonteCarloLibTest() {
	}
};

BOOST_GLOBAL_FIXTURE(AllMonteCarloLibTest)

