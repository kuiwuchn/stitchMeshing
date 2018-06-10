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

#ifdef _WIN32
#pragma warning(disable:4786)
#endif

#include <cassert>
#include <iostream>
#include "HE_Vertex.h"
#include "HE_HalfEdge.h"
#include "HE_Face.h"
#include "util.h"

cyPoint3f HE_Vertex::normal() const
{
	if (!_edge)
		return cyPoint3f(0, 0, 0);
	const_edge_circulator e = begin();
	const_edge_circulator sentinel = e;

	const double TINY = 0.0000001;
	cyPoint3f N(float(rangerand(-TINY, TINY)), float(rangerand(-TINY, TINY)), float(rangerand(-TINY, TINY)));
	do
	{
		assert(this != (*e)->dst());
		cyPoint3f v((*e)->tangent());
		++e;

		// Do not use the normal of faces which are actually holes
		if ((*e)->face()->hole())
			continue;

		cyPoint3f v2((*e)->tangent());
		float a = 1.0f;// v.angle(v2);
		if (a > 0.0001)
		{
			cyPoint3f n = v2 ^ v;
			N += a*n;
		}
	} while (e != sentinel);

	N.Normalize();
	return N;
}

int HE_Vertex::valence() const
{
	if (!_edge)
		return 0;
	const_edge_circulator e = begin();
	const_edge_circulator sentinel = e;

	int n(0);
	do
	{
		n++;
		++e;
	} while (e != sentinel);

	return n;
}


const HE_HalfEdge* HE_Vertex::edgeAdjacentTo(const HE_Face* f) const
{
	if (!_edge)
		return 0;
	const_edge_circulator e = begin();
	const_edge_circulator sentinel = e;

	do
	{
		if (f == (*e)->face())
			return *e;
		++e;
	} while (e != sentinel);

	return 0;
}


HE_HalfEdge* HE_Vertex::edgeTo(int v)
{
	if (!_edge)
		return 0;
	edge_circulator e = begin();
	edge_circulator sentinel = e;

	do
	{
		if (v == (*e)->dst()->index())
			return *e;
		++e;
	} while (e != sentinel);

	return 0;
}

const HE_HalfEdge* HE_Vertex::edgeTo(int v) const
{
	if (!_edge)
		return 0;
	const_edge_circulator e = begin();
	const_edge_circulator sentinel = e;

	do
	{
		if (v == (*e)->dst()->index())
			return *e;
		++e;
	} while (e != sentinel);

	return 0;
}


void HE_Vertex::edge_circulator::operator++()
{
	assert(_here->check());
	assert(_here->twin()->dst() == _base);
	_here = _here->twin()->next();
}

void HE_Vertex::const_edge_circulator::operator++()
{
	assert(_here->check());
	assert(_here->twin()->dst() == _base);
	_here = _here->twin()->next();
}

void HE_Vertex::edge_circulator::operator--()
{
	assert(_here->check());
	assert(_here->prev()->dst() == _base);
	_here = _here->prev()->twin();
}

void HE_Vertex::const_edge_circulator::operator--()
{
	assert(_here->check());
	assert(_here->prev()->dst() == _base);
	_here = _here->prev()->twin();
}


bool HE_Vertex::isBoundary() const
{
	if (!_edge)
		return true;
	const_edge_circulator e = begin();
	const_edge_circulator sentinel = e;

	do
	{
		if ((*e)->isBoundary())
			return true;
		++e;
	} while (e != sentinel);

	return false;
}

bool HE_Vertex::adjacent(const HE_Vertex* v2) const
{
	if (!_edge)
		return false;
	HE_Vertex::const_edge_circulator eit = begin();
	HE_Vertex::const_edge_circulator sentinel = eit;
	do
	{
		if ((*eit)->dst() == v2)
			return true;
		++eit;
	} while (eit != sentinel);
	return false;
}


bool HE_Vertex::check(bool holesFilled) const
{
	assert(_index >= 0);
	if (_edge && holesFilled)
	{
		HE_Vertex::const_edge_circulator eit = begin();
		HE_Vertex::const_edge_circulator sentinel = eit;
		do
		{
			const HE_HalfEdge* const heSentinel = *eit;
			const HE_HalfEdge* he = heSentinel;
			assert(heSentinel->twin()->dst() == this);
			do
			{
				assert(he->dst());
				assert(he->face());
				assert(he->next());
				assert(he->twin());
				assert(he->twin()->twin());
				assert(he == he->twin()->twin());
				he = he->next();
			} while (he != heSentinel);

			++eit;
		} while (eit != sentinel);
	}

	return true;
}

