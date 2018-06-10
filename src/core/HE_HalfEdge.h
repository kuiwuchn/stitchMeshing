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

#ifndef __HALFEDGE_H__
#define __HALFEDGE_H__

#include <ostream>
#include "cyPoint.h"

#define EFLAG_HORIZONTAL		(1<<16)
#define EFLAG_VERTICAL			(1<<17)

class HE_Vertex;
class HE_Face;
class HE_Polyhedron;

/**
* Class that represents the half-edges of the mesh.
* Each edge [a,b] of the mesh is represented by two HalfEdge objects, one from 'a' to 'b' and one the other way around.
* Each HalfEdge is also associated to the polygon (or hole) that is adjacent to that edge on a given side. By convention,
* the associated face for a HalfEdge is that which the HalfEdge traverses in countercloskwise order. For example, face [a,b,c]
* is associated to HalfEdge objects [a,b], [b,c] and [c,a]. The symmetric HalfEdge objects are associated to other faces.
*
* The four differentiating pieces of data a HalfEdge contains are its twin() (or symmetric) HalfEdge, the next() HalfEdge in
* the associated Face, the face() itself to which the HalfEdge is associated, and the Vertex dst() where the HalfEdge ends.
*/
class HE_HalfEdge
{
	friend class HE_Polyhedron;
public:
	/// Constructor.
	explicit HE_HalfEdge(HE_Vertex* v, HE_Face* f, HE_HalfEdge* t, HE_HalfEdge* n, int idx = -1)
		: _dst(v), _face(f), _twin(t), _next(n), _index(idx), _flags(0) {}
	/// Constructor.
	explicit HE_HalfEdge(HE_Vertex* v, HE_Face* f, int idx = -1)
		: _dst(v), _face(f), _twin(0), _next(0), _index(idx), _flags(0) {}

	// Set methods
	void next(HE_HalfEdge* n) { _next = n; } ///< Sets the next HalfEdge in the Face
	void twin(HE_HalfEdge* p) { _twin = p; } ///< Sets the symmetric HalfEdge to this one.
	void face(HE_Face* f) { _face = f; }  ///< Sets the associated Face object.
	void dst(HE_Vertex* v) { _dst = v; }  ///< Sets the destination Vertex object.

	// Get methods
	HE_HalfEdge* next() { return _next; }             ///< Returns the next HalfEdge in the Face.
	const HE_HalfEdge* next() const { return _next; } ///< Returns the next HalfEdge in the Face. (Const method)
	HE_HalfEdge* prev();                              ///< Returns the previous HalfEdge in the Face.
	const HE_HalfEdge* prev() const;                  ///< Returns the previous HalfEdge in the Face. (Const method)
	HE_HalfEdge* twin() { return _twin; }             ///< Returns the twin (symmetric) HalfEdge.
	const HE_HalfEdge* twin() const { return _twin; } ///< Returns the twin (symmetric) HalfEdge. (Const method)
	HE_Face* face() { return _face; }                 ///< Returns the associated Face object.
	const HE_Face* face() const { return _face; }     ///< Returns the associated Face object. (Const method)
	HE_Vertex* dst() { return _dst; }                 ///< Returns the destination Vertex of the HalfEdge.
	const HE_Vertex* dst() const { return _dst; }     ///< Returns the destination Vertex of the HalfEdge. (Const method)
	int index() const { return _index; }           ///< Returns the index of the HalfEdge in the Polyhedron.
	void index(int i) { _index = i; }                ///< Modifies the index of the HalfEdge (dangerous!).

													
	HE_Vertex* src() { return _twin->dst(); }			/// Returns the source Vertex of the HalfEdge. (DANGEROUS if _twin is not set!)
	/// Returns the source Vertex of the HalfEdge. (Const method) (DANGEROUS if _twin is not set!)
	const HE_Vertex* src() const { return _twin->dst(); }

	// Simple computations
	cyPoint3f normal() const;   ///< Computes the normal of the HalfEdge (average of the two adjacent faces).
	cyPoint3f tangent() const;  ///< Returns the tangent direction of the HalfEdge.
	cyPoint3f midpoint() const;  ///< Returns the middle point of the HalfEdge.
	float squaredLength() const { return float(tangent().LengthSquared()); } ///< Squared length of the HalfEdge.
	float length() const { return float(tangent().Length()); }  ///< Length of the segment.
	bool isBoundary() const;  ///< Tells if the half-edge is in a boundary of the mesh.

	bool degenerate() const { return src() == _dst; } ///< Is this a degenerate edge?
	bool check(bool holesFilled = true) const;   ///< Checks the local connectivity of the HalfEdge

	// flags
	void SetFlag(unsigned int flag) { _flags |= flag; }
	void ClearFlag(unsigned int flag) { _flags &= ~flag; }
	void ChangeFlag(unsigned int flag, bool set) { if (set) SetFlag(flag); else ClearFlag(flag); }
	bool TestFlag(unsigned int flag) const { return (_flags & flag) > 0; }
	void ClearFlags() { _flags = 0; }
	unsigned int GetFlags() const { return _flags; }

protected:
	HE_Vertex*			_dst;
	HE_Face*			_face;
	HE_HalfEdge*		_twin;
	HE_HalfEdge*		_next;
	int				_index;
	unsigned int	_flags;
};

#endif // __HALFEDGE_H__
