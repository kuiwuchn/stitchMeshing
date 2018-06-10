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

#include <cassert>
#include <cmath>
#include "HE_Face.h"
#include "HE_Vertex.h"
#include "HE_HalfEdge.h"
#include "util.h"

cyPoint3f HE_Face::normal() const
{
	//static Vector3Df prev(1,0,0); // $$$ Hack to prevent problems with degenerate normals...
	const double TINY = 0.0000001;
	cyPoint3f N(float(rangerand(-TINY, TINY)), float(rangerand(-TINY, TINY)), float(rangerand(-TINY, TINY)));
#if 0 // We can still compute the normal of a hole... Let the application decide its use.
	if (hole())
		return N;
#endif

	const HE_HalfEdge* sentinel = edge();
	const HE_HalfEdge* he = sentinel;
	assert(he);
	const HE_HalfEdge* he2 = he->next();
	assert(he2);
	do
	{
		cyPoint3f v(-he->tangent());
		cyPoint3f v2(he2->tangent());
		float a = 1.0f;// v.angle(v2);
		cyPoint3f n = v2 ^ v;
		N += a*n;
		he = he2;
		he2 = he2->next();
	} while (he != sentinel);

	N.Normalize();
#if 1
	//assert(!isnan(N.x()));
	return N;
#else
	if (!isnan(N[0]) && !isnan(N[1]) && !isnan(N[2]))
	{
		prevN = N;
		return N;
	}
	else
	{
		assert(!isnan(prevN.x()));
		return prevN;
	}
#endif
}

cyPoint3f HE_Face::centroid() const
{
	const_edge_circulator e = begin();
	assert(*e);
	const_edge_circulator sentinel = e;
	int n = 0;
	cyPoint3f C(0, 0, 0);
	do
	{
		n++;
		C += (*e)->dst()->position();
		++e;
	} while (e != sentinel);

	return C / float(n);
}

int HE_Face::size() const
{
	HE_Face::const_edge_circulator e = begin();
	assert(*e);
	HE_Face::const_edge_circulator sentinel = e;
	int n(0);
	do
	{
		n++;
		++e;
	} while (sentinel != e);

	return n;
}

int HE_Face::NumHalfEdge() const
{
	HE_Face::const_edge_circulator e = begin();
	assert(*e);
	int i = 0;
	HE_Face::const_edge_circulator sentinel = e;
	do
	{
		++e;
		++i;
	} while (sentinel != e);

	return i;
}

HE_EDGETYPE HE_Face::HalfEdgeType(const HE_HalfEdge* he) const
{
	HE_Face::const_edge_circulator e = begin();
	assert(*e);
	int i = 0;
	HE_Face::const_edge_circulator sentinel = e;
	do
	{
		if (*e == he) {
			if (TestFlag(VEF_FLAG_QUAD))
			{
				switch (i)
				{
				case 0: return HE_TOP;
				case 1: case 3: return HE_SIDE;
				case 2: return HE_BTM;
				default:
					break;
				}
			}
			else if (TestFlag(VEF_FLAG_INC))
			{
				switch (i)
				{
				case 0: case 4: return HE_TOP;
				case 1: case 3: return HE_SIDE;
				case 2: return HE_BTM;
				default:
					break;
				}
			}
			else if (TestFlag(VEF_FLAG_DEC))
			{
				switch (i)
				{
				case 0: return HE_TOP;
				case 1: case 4: return HE_SIDE;
				case 2: case 3: return HE_BTM;
				default:
					break;
				}
			}
			else if (TestFlag(VEF_FLAG_QUAD_SHORTROW_L))
			{
				switch (i)
				{
				case 0: return HE_TOP;
				case 2: case 3: return HE_SIDE;
				case 1: return HE_BTM;
				default:
					break;
				}
			}
			else if (TestFlag(VEF_FLAG_QUAD_SHORTROW_R))
			{
				switch (i)
				{
				case 1: return HE_TOP;
				case 2: case 3: return HE_SIDE;
				case 0: return HE_BTM;
				default:
					break;
				}
			}
		}
		++e;
		++i;
	} while (sentinel != e);

	return HE_NONE;
}

int HE_Face::HalfEdgeIdx(const HE_HalfEdge* he) const
{
	HE_Face::const_edge_circulator e = begin();
	assert(*e);
	int i = 0;
	HE_Face::const_edge_circulator sentinel = e;
	do
	{
		if (*e == he)
			return i;
		++e;
		++i;
	} while (sentinel != e);

	return -1;
}

bool HE_Face::contains(const HE_HalfEdge* he) const
{
	HE_Face::const_edge_circulator e = begin();
	assert(*e);
	HE_Face::const_edge_circulator sentinel = e;
	do
	{
		if (*e == he)
			return true;
		++e;
	} while (sentinel != e);

	return false;
}

bool HE_Face::contains(const HE_Vertex* v) const
{
	HE_Face::const_edge_circulator e = begin();
	assert(*e);
	HE_Face::const_edge_circulator sentinel = e;
	do
	{
		if ((*e)->dst() == v)
			return true;
		++e;
	} while (sentinel != e);

	return false;
}

const HE_HalfEdge* HE_Face::fromVertex(int v) const
{
	HE_Face::const_edge_circulator e = begin();
	assert(*e);
	HE_Face::const_edge_circulator sentinel = e;
	do
	{
		if ((*e)->dst()->index() == v)
			return (*e)->next();
		++e;
	} while (sentinel != e);

	return 0;
}

const HE_HalfEdge* HE_Face::toVertex(int v) const
{
	HE_Face::const_edge_circulator e = begin();
	assert(*e);
	HE_Face::const_edge_circulator sentinel = e;
	do
	{
		if ((*e)->dst()->index() == v)
			return *e;
		++e;
	} while (sentinel != e);

	return 0;
}

const HE_HalfEdge* HE_Face::fromVertex(HE_Vertex* v) const
{
	return fromVertex(v->index());
}

const HE_HalfEdge* HE_Face::toVertex(HE_Vertex* v) const
{
	return toVertex(v->index());
}


bool HE_Face::adjacent(const HE_Face* f2) const
{
	HE_Face::const_edge_circulator eit = begin();
	assert(*eit);
	HE_Face::const_edge_circulator sentinel = eit;
	do
	{
		if ((*eit)->twin()->face() == f2)
			return true;
		++eit;
	} while (eit != sentinel);
	return false;
}

bool HE_Face::degenerate() const
{
	HE_Face::const_edge_circulator eit = begin();
	assert(*eit);
	HE_Face::const_edge_circulator sentinel = eit;
	do
	{
		const HE_HalfEdge* he = *eit;
		if (he->src() == he->dst())
			return true;
		++eit;
	} while (eit != sentinel);
	return false;
}
