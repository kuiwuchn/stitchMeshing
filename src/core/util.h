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

#ifndef __UTIL_H__
#define __UTIL_H__

#include <exception>
#include <string>

#define M_PI 3.14159265358979323846264338327950288

class HE_exception : public std::exception
{
public:
	/**
	*    Constructor.
	* @param msg Message to be shown when what() is invoked. It should inform about the cause of the exception.
	* @see what
	*/
	HE_exception(const char* msg) : _msg(msg) {}

	/**
	*    Constructor.
	* @param msg Message to be shown when what() is invoked. It should inform about the cause of the exception.
	* @param msg2 Message to be shown when what() is invoked. It should complement the information in msg.
	* @see what
	*/
	HE_exception(const char* msg, const char* msg2) : _msg((std::string(msg) + ": ") + msg2) {}

	/**
	* Destructor.
	*/
	virtual ~HE_exception() throw() {}

	/**
	*    Returns a string with information about the exception.
	* @return The string with information about the exception.
	*/
	virtual const char* what() const throw() { return _msg.c_str(); }

protected:
	std::string _msg;
};

/**
* Ranged uniform random number generator. It produces a random number in the [min,max] interval, and has the granularity indicated by steps.
* @param min Lower bound of the interval.
* @param max Upper bound of the interval.
* @param steps Granularity of the interval. Number of distinct values that can be returned.
* @return A value in the interval.
*/
double gaussianrand(double min, double max, long steps = 10000);

/**
*  Gaussian-alike (bell shaped) random number generator. It produces a random number in the [min,max] interval, with a
* distribution that is similar to that of a Gaussian function. Note, however, that the actual distribution is NOT Gaussian.
* @param min Lower bound of the interval.
* @param max Upper bound of the interval.
* @param steps Granularity of the interval. Number of distinct values that can be returned.
* @return A value in the interval.
*/
double rangerand(double min, double max, long steps = 10000);

#endif // __UTIL_H__