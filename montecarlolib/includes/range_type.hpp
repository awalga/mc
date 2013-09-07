/*
 * Copyright (C) 2013  Abel Walga
 *
 *  range_type.hpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  range_type.hpp is distributed in the hope that it will be useful,
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

#ifndef RANGE_TYPE_HPP_
#define RANGE_TYPE_HPP_

#include <iostream>
#include <cstddef>

namespace util {

/**
 * <p>Define dimension type as size_t which is dependant on the plateform</p>
 */
typedef size_t dimension;

/**
 * <p>Define real value as double</p>
 */
typedef double Real;

/**
 * <p>Define integer value as int</p>
 */
typedef long Integer;

/**
 * <p>Define positive integer as unsigned int</p>
 */
typedef unsigned int PositiveInteger;

} // util
#endif /* RANGE_TYPE_HPP_ */
