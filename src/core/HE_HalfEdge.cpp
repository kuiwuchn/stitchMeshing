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

#include <cstdio>
#include <cassert>
#include <set>
#include "HE_HalfEdge.h"
#include "HE_Face.h"
#include "HE_Vertex.h"

HE_HalfEdge* HE_HalfEdge::prev()
{
	HE_HalfEdge* he = _next;
	assert(he);
	while (he && he->next() != this)
	{
		assert(he->next());
		he = he->next();
		assert(he);
	}

	return he;
}

const HE_HalfEdge* HE_HalfEdge::prev() const
{
	HE_HalfEdge* he = _next;
	assert(he);
	while (he && he->next() != this)
	{
		assert(he->next());
		he = he->next();
		assert(he);
	}

	return he;
}

cyPoint3f HE_HalfEdge::normal() const
{
	return (src()->normal() + dst()->normal()) / 2;
}


cyPoint3f HE_HalfEdge::tangent() const
{
	return dst()->position() - src()->position();
}

cyPoint3f HE_HalfEdge::midpoint() const
{
	return (src()->position() + dst()->position()) / 2;
}

bool HE_HalfEdge::isBoundary() const
{
	return face()->hole() || twin()->face()->hole();
}

bool HE_HalfEdge::check(bool holesFilled) const
{
	assert(prev()->next() == this);
	if (holesFilled)
	{
		assert(twin());
		assert(twin() != this);
		assert(twin()->twin() == this);
		assert(prev()->dst() == twin()->dst());
		assert(prev()->dst() == src());
		assert(prev()->dst() == twin()->next()->src());
		assert(src() == twin()->dst());
		assert(src() == twin()->next()->src());
		assert(src() == twin()->next()->src());
		assert(src() != dst());
		//	 assert(next()->dst() != twin()->next()->dst());
	}
	assert(next() != 0);
	assert(face() != 0);
	assert(dst() != 0);
	assert(!(*dst() == *prev()->dst()));

	// Make sure the half edge is correctly linked in a face, and that the vertices are unique
	{
		std::set<const HE_HalfEdge*> visited;
		std::set<int> iVisited;
		const HE_HalfEdge* e = this;
		const HE_HalfEdge* sentinel = e;
		do
		{
			assert(visited.end() == visited.find(e));
			visited.insert(e);
			assert(iVisited.end() == iVisited.find(e->dst()->index()));
			iVisited.insert(e->dst()->index());
			e = e->next();
		} while (e != sentinel);
	}

	assert(prev()->next() == this);

	return true;
}
