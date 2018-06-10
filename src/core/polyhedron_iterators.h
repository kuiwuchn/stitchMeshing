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

// This file is to be included in polyhedron.h, just to keep it simple and short.
// Here we declare (and implement most of) the iterators that traverse faces,
// edges and vertices of a Polyhedron.

class const_face_iterator;
class const_edge_iterator;
class const_vertex_iterator;

/// Face iterator for the Polyhedron class. It traverses the faces of the Polyhedron, one at a time.
class face_iterator : public std::iterator<std::forward_iterator_tag, HE_Face*>
{
public:
	friend class const_face_iterator;
	face_iterator(const HE_Polyhedron* p, unsigned f) : _base(p), _here(f) {}

	HE_Face* operator*() { return _base->_faces[_here]; }
	void operator++() { ++_here; }
	void operator++(int) { ++_here; }
	bool operator==(const face_iterator& fi) const { return _here == fi._here && _base == fi._base; }
	bool operator!=(const face_iterator& fi) const { return !(*this == fi); }

protected:
	const HE_Polyhedron* _base;
	unsigned _here;
};

/// Const Face iterator for the Polyhedron class. It traverses the faces of the Polyhedron, one at a time.
class const_face_iterator : public std::iterator<std::forward_iterator_tag, const HE_Face*>
{
public:
	const_face_iterator(const HE_Polyhedron* p, unsigned f) : _base(p), _here(f) {}
	const_face_iterator(const face_iterator& f) : _base(f._base), _here(f._here) {}

	const HE_Face* operator*() const { return _base->_faces[_here]; }
	void operator++() { ++_here; }
	void operator++(int) { ++_here; }
	bool operator==(const const_face_iterator& fi) const { return _here == fi._here && _base == fi._base; }
	bool operator!=(const const_face_iterator& fi) const { return !(*this == fi); }

protected:
	const HE_Polyhedron* _base;
	unsigned _here;
};

face_iterator fBegin() { return face_iterator(this, 0); }        ///< Iterator for the first Face of the Polyhedron
face_iterator fEnd() { return face_iterator(this, numFaces()); } ///< Iterator past the last Face of the Polyhedron
																 /// Const Iterator for the first Face of the Polyhedron
const_face_iterator fBegin() const { return const_face_iterator(this, 0); }
/// Const Iterator past the last Face of the Polyhedron
const_face_iterator fEnd() const { return const_face_iterator(this, numFaces()); }



/// HalfEdge iterator for the Polyhedron class. It traverses the hald-edges of the Polyhedron, one at a time.
class edge_iterator
{
public:
	friend class const_edge_iterator;
	edge_iterator(const HE_Polyhedron* p, unsigned e) : _base(p), _here(e) {}

	HE_HalfEdge* operator*() { return _base->_halfEdges[_here]; }
	void operator++() { ++_here; }
	void operator++(int) { ++_here; }
	bool operator==(const edge_iterator& ei) const { return _here == ei._here && _base == ei._base; }
	bool operator!=(const edge_iterator& ei) const { return !(*this == ei); }

protected:
	const HE_Polyhedron* _base;
	unsigned _here;
};

/// Const HalfEdge iterator for the Polyhedron class. It traverses the hald-edges of the Polyhedron, one at a time.
class const_edge_iterator
{
public:
	const_edge_iterator(const HE_Polyhedron* p, unsigned e) : _base(p), _here(e) {}
	const_edge_iterator(const edge_iterator& e) : _base(e._base), _here(e._here) {}

	const HE_HalfEdge* operator*() const { return _base->_halfEdges[_here]; }
	void operator++() { ++_here; }
	void operator++(int) { ++_here; }
	bool operator==(const const_edge_iterator& ei) const { return _here == ei._here && _base == ei._base; }
	bool operator!=(const const_edge_iterator& ei) const { return !(*this == ei); }

protected:
	const HE_Polyhedron* _base;
	unsigned _here;
};
edge_iterator eBegin() { return edge_iterator(this, 0); }            ///< Iterator for the first HalfEdge of the Polyhedron
edge_iterator eEnd() { return edge_iterator(this, numHalfEdges()); } ///< Iterator past the last HalfEdge of the Polyhedron

																	 /// Const Iterator for the first HalfEdge of the Polyhedron
const_edge_iterator eBegin() const { return const_edge_iterator(this, 0); }
/// Const Iterator past the last HalfEdge of the Polyhedron
const_edge_iterator eEnd() const { return const_edge_iterator(this, numHalfEdges()); }


/// Vertex iterator for the Polyhedron class. It traverses the vertices of the Polyhedron, one at a time.
class vertex_iterator
{
public:
	friend class const_vertex_iterator;
	vertex_iterator(const HE_Polyhedron* p, unsigned v) : _base(p), _here(v) {}

	HE_Vertex* operator*() { return _base->_vertices[_here]; }
	void operator++() { ++_here; }
	void operator++(int) { ++_here; }
	bool operator==(const vertex_iterator& vi) const { return _here == vi._here && _base == vi._base; }
	bool operator!=(const vertex_iterator& vi) const { return !(*this == vi); }

protected:
	const HE_Polyhedron* _base;
	unsigned _here;
};

/// Const Vertex iterator for the Polyhedron class. It traverses the vertices of the Polyhedron, one at a time.
class const_vertex_iterator
{
public:
	const_vertex_iterator(const HE_Polyhedron* p, unsigned v) : _base(p), _here(v) {}
	const_vertex_iterator(const vertex_iterator& v) : _base(v._base), _here(v._here) {}

	const HE_Vertex* operator*() { return _base->_vertices[_here]; }
	void operator++() { ++_here; }
	void operator++(int) { ++_here; }
	bool operator==(const const_vertex_iterator& vi) const { return _here == vi._here && _base == vi._base; }
	bool operator!=(const const_vertex_iterator& vi) const { return !(*this == vi); }

protected:
	const HE_Polyhedron* _base;
	unsigned _here;
};

vertex_iterator vBegin() { return vertex_iterator(this, 0); }           ///< Iterator for the first Vertex of the Polyhedron
vertex_iterator vEnd() { return vertex_iterator(this, numVertices()); } ///< Iterator past the last Vertex of the Polyhedron

																		/// Const Iterator for the first Vertex of the Polyhedron
const_vertex_iterator vBegin() const { return const_vertex_iterator(this, 0); }

/// Const Iterator past the last Vertex of the Polyhedron
const_vertex_iterator vEnd() const { return const_vertex_iterator(this, numVertices()); }

