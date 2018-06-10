/***************************************************************************
*   Copyright (C) 2007 by Pablo Diaz-Gutierrez   *
*   pablo@ics.uci.edu   *
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU Library General Public License as       *
*   published by the Free Software Foundation; either version 2 of the    *
*   License, or (at your option) any later version.                       *
*                                                                         *
*   This program is distributed in the hope that it will be useful,       *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*   GNU General Public License for more details.                          *
*                                                                         *
*   You should have received a copy of the GNU Library General Public     *
*   License along with this program; if not, write to the                 *
*   Free Software Foundation, Inc.,                                       *
*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
***************************************************************************/

#include <math.h>
#include <cassert>
#include "util.h"

double rangerand(double min, double max, long steps)
{
	assert(max > min);
	return min + ((rand() % steps) * (max - min)) / steps;
}

double gaussianrand(double min, double max, long steps)
{
	assert(max > min);
	double halfinterval = (max - min) / 2;
	double center = min + halfinterval;
	double factor = rangerand(0, 2 * M_PI, steps);
	double value = cos(factor)*sin(factor);

	return center + value*halfinterval;
}