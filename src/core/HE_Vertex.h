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

#ifndef __VERTEX_H__
#define __VERTEX_H__

#include "cyPoint.h"

class HE_HalfEdge;
class HE_Face;
class HE_Polyhedron;

/// Class that represents the vertices of the mesh. It basic topologic and geometric information.
class HE_Vertex
{
	friend class HE_Polyhedron;
public:
	/// Circulator for the edges departing from this Vertex.
	class edge_circulator;
	class const_edge_circulator;
	class edge_circulator
	{
		friend class const_edge_circulator;
	public:
		edge_circulator(const HE_Vertex* v, HE_HalfEdge* h) : _base(v), _here(h) {}

		const HE_Vertex* base() const { return _base; };
		HE_HalfEdge* operator*() { return _here; }
		void operator++();
		void operator++(int) { ++(*this); }
		void operator--();
		void operator--(int) { --(*this); }
		bool operator==(const edge_circulator& e) const { return _base == e._base && _here == e._here; }
		bool operator!=(const edge_circulator& e) const { return !(e == *this); }
	protected:
		const HE_Vertex* _base;
		HE_HalfEdge* _here;
	};

	/// Const circulator for the edges departing from this Vertex.
	class const_edge_circulator
	{
		friend class edge_circulator;
	public:
		const_edge_circulator(const HE_Vertex* v, const HE_HalfEdge* h) : _base(v), _here(h) {}
		const_edge_circulator(const edge_circulator& ec) : _base(ec.base()), _here(ec._here) {}

		const HE_Vertex* base() const { return _base; };
		const HE_HalfEdge* operator*() const { return _here; }
		void operator++();
		void operator++(int) { ++(*this); }
		void operator--();
		void operator--(int) { --(*this); }
		bool operator==(const const_edge_circulator& e) const { return _base == e._base && _here == e._here; }
		bool operator!=(const const_edge_circulator& e) const { return !(e == *this); }
	protected:
		const HE_Vertex* _base;
		const HE_HalfEdge* _here;
	};

	/// Circulator referring to an arbitrary departing HalfEdge.
	edge_circulator begin() { return edge_circulator(this, _edge); }
	/// Const circulator referring to an arbitrary departing HalfEdge.
	const_edge_circulator begin() const { return const_edge_circulator(this, _edge); }

	explicit HE_Vertex(int i = -1) : _edge(0), _index(i) {}    ///< Constructor
	explicit HE_Vertex(const cyPoint3f& pos, int i = -1) : _pos(pos), _edge(0), _index(i) {}    ///< Constructor
	explicit HE_Vertex(const cyPoint3f& pos, HE_HalfEdge* e, int i = -1) : _pos(pos), _edge(e), _index(i) {}  ///< Constructor

																										// Set methods
	void index(int i) { _index = i; } ///< Sets the index of the Vertex in the mesh.
	void position(const cyPoint3f& pos) { _pos = pos; } ///< Sets the position of the Vertex.
	void edge(HE_HalfEdge* e) { _edge = e; }  ///< Associates an arbitrary HalfEdge to the Vertex.

										   // Get methods
	int index() const { return _index; }  ///< Returns the index of the Vertex object.
	const cyPoint3f& position() const { return _pos; } ///< Returns the position of the Vertex.
	HE_HalfEdge* edge() { return _edge; }              ///< Returns an arbitrary HalfEdge associated to the Vertex.
	const HE_HalfEdge* edge() const { return _edge; }  ///< Returns an arbitrary const HalfEdge associated to the Vertex.
	bool adjacent(const HE_Vertex* V2) const;

	// Simple computations
	cyPoint3f normal() const;    ///< Computes the normal of the Vertex.
	int valence() const;    ///< Counts the number of HalfEdge objects departing from the Vertex.
	int degree() const { return valence(); }    ///< Counts the number of HalfEdge objects departing from the Vertex.
	bool isBoundary() const;

	bool operator==(const HE_Vertex& v) const { return _index == v._index; } /// Are these vertices the same one?
																		  /// Returns the HalfEdge that departs from this Vertex and ends in Vertex 'v'.
	HE_HalfEdge* edgeTo(const HE_Vertex* v) { return edgeTo(v->index()); }
	/// Returns the HalfEdge that departs from this Vertex and ends in Vertex 'v'.
	const HE_HalfEdge* edgeTo(const HE_Vertex* v) const { return edgeTo(v->index()); }

	/// Returns the HalfEdge that departs from this Vertex and ends in Vertex with index 'v'.
	HE_HalfEdge* edgeTo(int v);
	/// Returns the HalfEdge that departs from this Vertex and ends in Vertex with index 'v'.
	const HE_HalfEdge* edgeTo(int v) const;

	/// Returns the HalfEdge that departs from this Vertex and is adjacent to Face 'f'.
	const HE_HalfEdge* edgeAdjacentTo(const HE_Face* f) const;

	bool check(bool holesFilled = true) const; ///< Sanity check.

											   /// Returns the cyPoint3f between the positions of two vertices.
	cyPoint3f operator-(const HE_Vertex& V) const { return V._pos - _pos; }

protected:
	cyPoint3f _pos;
	HE_HalfEdge* _edge;
	int _index;
};

#endif // __VERTEX_H__