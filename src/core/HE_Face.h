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

#ifndef __FACE_H__
#define __FACE_H__

#include "cyPoint.h"
#include "HE_HalfEdge.h"

class HE_Vertex;
class HE_Polyhedron;

#define VEF_FLAG_NONE				(1<<0)
#define VEF_FLAG_QUAD				(1<<1)
#define VEF_FLAG_SHORTROW_L			(1<<2)
#define VEF_FLAG_SHORTROW_R			(1<<3)
#define VEF_FLAG_INC				(1<<4)
#define VEF_FLAG_DEC				(1<<5)

#define VEF_FLAG_QUAD_SHORTROW_L	(1<<6)
#define VEF_FLAG_QUAD_SHORTROW_R	(1<<7)
#define VEF_FLAG_QUAD_INC			(1<<8)
#define VEF_FLAG_QUAD_DEC			(1<<9)

enum HE_EDGETYPE
{
	HE_NONE = 0,
	HE_TOP,
	HE_BTM,
	HE_SIDE,
};

/**
* Class that represents the faces of the mesh. It contains a point to an arbitrary HalfEdge in its boundary,
* which can be used to traverse the whole face.
*/
class HE_Face
{
	friend class HE_Polyhedron;

public:
	class const_edge_circulator;
	class edge_circulator;

	/// Circulator that traverses the HalfEdge(s) that compound the boundary of the Race.
	class edge_circulator
	{
		friend class const_edge_circulator;
	public:
		edge_circulator(const HE_Face* f, HE_HalfEdge* h) : _base(f), _here(h) {}

		HE_HalfEdge* operator*() { return _here; }
		void operator++() { _here = _here->next(); }
		void operator++(int) { ++(*this); }
		void operator--() { _here = _here->prev(); }
		void operator--(int) { --(*this); }
		bool operator==(const edge_circulator& e) const { return _base == e._base && _here == e._here; }
		bool operator!=(const edge_circulator& e) const { return !(e == *this); }
	protected:
		const HE_Face* _base;
		HE_HalfEdge* _here;
	};
	/**
	*    Reference to an arbitrary HalfEdge of the Race, used for traversal.
	* @return An edge_circulator that refers to an arbitrary HalfEdge of the Race.
	*/
	edge_circulator begin() { return edge_circulator(this, _edge); }

	/// Constant circulator that traverses the HalfEdge(s) that compound the boundary of the Race.
	class const_edge_circulator
	{
		friend class edge_circulator;
	public:
		const_edge_circulator(const HE_Face* f, const HE_HalfEdge* h) : _base(f), _here(h) {}
		const_edge_circulator(const edge_circulator& ec) : _base(ec._base), _here(ec._here) {}

		const HE_HalfEdge* operator*() { return _here; }
		void operator++() { _here = _here->next(); }
		void operator++(int) { ++(*this); }
		void operator--() { _here = _here->prev(); }
		void operator--(int) { --(*this); }
		bool operator==(const const_edge_circulator& e) const { return _base == e._base && _here == e._here; }
		bool operator!=(const const_edge_circulator& e) const { return !(e == *this); }
	protected:
		const HE_Face* _base;
		const HE_HalfEdge* _here;
	};
	/**
	*    Reference to an arbitrary const HalfEdge of the Race, used for traversal.
	* @return A const_edge_circulator that refers to an arbitrary const HalfEdge of the Race.
	*/
	const_edge_circulator begin() const { return const_edge_circulator(this, _edge); }

	/**
	*    Constructor.
	* @param e A HalfEdge that is to be associated to the Race.
	* @param hole Indicates whether the Race is actually a hole. By default, false.
	*/
	explicit HE_Face(HE_HalfEdge* e, bool hole = false, int idx = -1) : _edge(e), _hole(hole), _index(idx), _flags(0) {}

	// Set methods
	void hole(bool h) { _hole = h; }       ///< Sets the 'hole' flag of the Race to true of false.
	void edge(HE_HalfEdge* e) { _edge = e; }  ///< Sets the referenc HalfEdge of the Race.

	// Get methods
	bool hole() const { return _hole; }            ///< Tells if the Race is actually a hole.
	HE_HalfEdge* edge() { return _edge; }             ///< Returns an arbitrary HalfEdge on the Race's boundary.
	const HE_HalfEdge* edge() const { return _edge; } ///< Returns an arbitrary const HalfEdge on the Race's boundary.
	const HE_HalfEdge* fromVertex(int v) const;       ///< Returns a const HalfEdge on the Race that starts at vertex with index 'v'
	const HE_HalfEdge* fromVertex(HE_Vertex* v) const;   ///< Returns a const HalfEdge on the Race that starts at Vertex 'v'
	const HE_HalfEdge* toVertex(int v) const;         ///< Returns a const HalfEdge on the Race that ends at vertex with index 'v'
	const HE_HalfEdge* toVertex(HE_Vertex* v) const;     ///< Returns a const HalfEdge on the Race that ends at Vertex 'v'

	// Simple computations
	int size() const;                        ///< Number of HalfEdge's in the Race.
	int index() const { return _index; }     ///< Index of the face in the Polyhedron.
	void index(int i) { _index = i; }          ///< Modifies the index of the Face (dangerous!).
	bool contains(const HE_HalfEdge* he) const; ///< Tells if this Face contains HalfEdge he.
	bool contains(const HE_Vertex* v) const;    ///< Tells if this Face contains Vertex v.
	bool adjacent(const HE_Face* f2) const;     ///< Tells if this Face is adjacent to f2.
	bool degenerate() const;                 ///< Tells is the Face contains degenerate features like repeated vertices.
	cyPoint3f normal() const;      ///< Computes the normal of the Face.
	cyPoint3f centroid() const;    ///< Computes the centroid of the Face.

	int HalfEdgeIdx(const HE_HalfEdge* he) const;
	int NumHalfEdge() const;

	HE_EDGETYPE HalfEdgeType(const HE_HalfEdge* he) const;

	// flags
	void SetFlag(unsigned int flag) { _flags |= flag; }
	void ClearFlag(unsigned int flag) { _flags &= ~flag; }
	void ChangeFlag(unsigned int flag, bool set) { if (set) SetFlag(flag); else ClearFlag(flag); }
	bool TestFlag(unsigned int flag) const { return (_flags & flag) > 0; }
	void ClearFlags() { _flags = 0; }
	unsigned int GetFlags() const { return _flags; }

protected:
	HE_HalfEdge*		_edge;
	bool			_hole;
	int				_index;
	unsigned int	_flags;
};

#endif // __FACE_H__
