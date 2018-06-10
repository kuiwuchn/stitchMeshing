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

#include <algorithm>
#include <stack>
#include <iterator>
#include <float.h>
#include <assert.h>
#include <iostream>
#include "HE_Polyhedron.h"
#include "util.h"

// Duplicates Face* data from 'faces', adding it to 'p'.
void duplicateFaces(HE_Polyhedron& p, const std::set<const HE_Face*>& faces)
{
	// First, copy the vertices used in 'faces'.
	std::map<const HE_Vertex*, HE_Vertex*> vertMap;
	for (std::set<const HE_Face*>::const_iterator it = faces.begin(); it != faces.end(); ++it)
	{
		const HE_Face* f = *it;
		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;
		do
		{
			const HE_HalfEdge* he = *e;
			const HE_Vertex* V = he->dst();
			HE_Vertex* vHere;
			if (vertMap.end() == vertMap.find(V))
				vertMap[V] = p.addVertex(V->position());
			vHere = vertMap[V];
			++e;
		} while (e != sentinel);
	}

	// Now, copy the faces themselves
	std::vector<int> corners;
	for (std::set<const HE_Face*>::const_iterator it = faces.begin(); it != faces.end(); ++it)
	{
		const HE_Face* f = *it;
		corners.clear();

		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;
		do
		{
			const HE_HalfEdge* he = *e;
			const HE_Vertex* V = he->dst();
			assert(vertMap.end() != vertMap.find(V));
			HE_Vertex* vHere = vertMap[V];
			corners.push_back(vHere->index());
			++e;
		} while (e != sentinel);

		p.addFace(corners);
	}
}


HE_Polyhedron::HE_Polyhedron(const HE_Polyhedron& p)
{
	std::set<const HE_Face*> faces;
	for (const_face_iterator fit = p.fBegin(); fit != p.fEnd(); ++fit)
	{
		const HE_Face* f = *fit;
		if (f->hole())
			continue;
		faces.insert(f);
	}
	duplicateFaces(*this, faces);
	finalize();
}


HE_Polyhedron::HE_Polyhedron(const char* filename) : _type(FT_OFF)
{
	if (!load(filename))
		throw HE_exception("Cannot load input file", filename);
}

HE_Polyhedron::HE_Polyhedron(const MatrixXf &pPos, std::vector<std::vector<uint32_t>> &pIndices) : _type(FT_NONE)
{
	std::vector<std::vector<int> > faces;

	//std::cout << "v " << pPos.cols() << std::endl;
	_vertices.resize((uint32_t)pPos.cols());
	for (unsigned v = 0; v < pPos.cols(); v++)
		_vertices[v] = new HE_Vertex(cyPoint3f(pPos(0, v), pPos(1, v), pPos(2, v)), v);

	//std::cout << "f " << pIndices.size() << std::endl;
	faces.resize((uint32_t)pIndices.size()); 

	for (int i = 0; i < pIndices.size(); i++) {
		if (pIndices[i].size() == 5) std::cout << "ERROR: there should not be pentagon!\n";
		for (int j = 0; j < pIndices[i].size(); j++)
			faces[i].push_back(pIndices[i][j]);
	}

	loadFaces(faces);
	finalize();
}

HE_Polyhedron::HE_Polyhedron()
	: _type(FT_NONE)
{
}

HE_Polyhedron::HE_Polyhedron(const std::vector<HE_Vertex>& verts, const std::vector<std::vector<int> >& faces) : _type(FT_NONE)
{
	_vertices.resize(verts.size());
	for (unsigned v = 0; v<verts.size(); v++)
		_vertices[v] = new HE_Vertex(verts[v].position(), v);
	loadFaces(faces);
	finalize();
}

HE_Polyhedron::HE_Polyhedron(const std::set<const HE_Face*>& faces) : _type(FT_NONE)
{
	duplicateFaces(*this, faces);
	finalize();
}

HE_Polyhedron::~HE_Polyhedron()
{
	for (face_iterator fit = fBegin(); fit != fEnd(); ++fit)
		delete *fit;
	_faces.clear();

	for (edge_iterator eit = eBegin(); eit != eEnd(); ++eit)
		delete *eit;
	_halfEdges.clear();

	for (vertex_iterator vit = vBegin(); vit != vEnd(); ++vit)
		delete *vit;
	_vertices.clear();
}

HE_Vertex* HE_Polyhedron::addVertex(const cyPoint3f& v)
{
	HE_Vertex* V = new HE_Vertex(v, _vertices.size());
	_vertices.push_back(V);
	return V;
}

void HE_Polyhedron::clearData()
{
	_faces.clear();
	_holes.clear();
	_halfEdges.clear();
	_vertices.clear();
}

bool HE_Polyhedron::load(const char* filename)
{
	filetype_e format = fileFormatFromFileName(filename);

	clearData();
	bool success = false;
	switch (format)
	{
	case FT_PLY_ASCII:
		success = readPly(filename);
		type(FT_PLY_ASCII);
		break;

	case FT_OBJ:
		std::cout << "Loading .obj file\n";
		success = readObj(filename);
		type(FT_OBJ);
		break;

	case FT_NONE:
		throw HE_exception("Unknown input file format", filename);
		type(FT_NONE);
		break;
	}

	//  packVertices();
	return success;
}


struct IndexedVertex
{
	typedef float value_type;

	IndexedVertex(HE_Vertex* V) : _V(V) {}
	HE_Vertex* _V;
	inline value_type operator[] (size_t const k) const { return _V->position()[k]; }
};
inline bool operator== (const IndexedVertex& a, const IndexedVertex& b) { return a._V->index() == b._V->index(); }

void insertAdjacentFaces(HE_Vertex* V, std::set<HE_Face*>& faces)
{
	HE_Vertex::edge_circulator eit = V->begin();
	const HE_Vertex::edge_circulator sentinel = eit;
	do
	{
		HE_HalfEdge* he = *eit;
		assert(he);
		HE_Face* f = he->face();
		assert(f);
		if (faces.end() == faces.find(f))
			faces.insert(f);
		++eit;
	} while (eit != sentinel);
}


// Adds all the faces that are adjacent to any vertices in 'verts'
void insertAdjacentFaces(HE_Polyhedron* poly, std::set<HE_Face*>& faces, const std::set<const HE_Vertex*>& verts)
{
	for (HE_Polyhedron::edge_iterator it = poly->eBegin(); it != poly->eEnd(); ++it)
	{
		HE_HalfEdge* he = *it;
		const HE_Vertex* V = he->dst();
		if (verts.end() != verts.find(V))
		{
			HE_Face* f = he->face();
			assert(f);
			if (faces.end() == faces.find(f))
				faces.insert(f);
		}
	}
}


void HE_Polyhedron::finalize()
{
	assert(check(false));

	// First, find the connections that need finalization. During construction, this is all of them.
	std::map<std::pair<int, int>, HE_HalfEdge*> connections;
	std::set<const HE_Vertex*> problemVerts;
	int n(0);
	for (edge_iterator it = eBegin(); it != eEnd(); ++it)
	{
		HE_HalfEdge* he = *it;
		assert(he->check(false));
		assert(!he->twin());
		HE_Vertex* Va = he->prev()->dst();
		HE_Vertex* Vb = he->dst();
		assert(Va->check(false));
		assert(Vb->check(false));
		int a = Va->index();
		int b = Vb->index();
		assert(a != b);
		std::pair<int, int> idx(a, b);
		if (connections.end() != connections.find(idx))
		{
#if 1
			std::cerr << ++n << "th non-manifold edge in finalize(), endpoints [" << a << ',' << b << "]\n";
			if (problemVerts.end() == problemVerts.find(Va))
				problemVerts.insert(Va);
			if (problemVerts.end() == problemVerts.find(Vb))
				problemVerts.insert(Vb);
			continue;
#else
			throw HE_exception("Repeated half-edge: The half-edge data structure only supports manifolds.");
#endif
		}
		connections[idx] = he;

		std::pair<int, int> idx2(b, a);
		std::map<std::pair<int, int>, HE_HalfEdge*>::iterator cit = connections.find(idx2);
		if (connections.end() != cit)
		{
			HE_HalfEdge* he2 = cit->second;
			assert(he != he2);
			he->twin(he2);
			he2->twin(he);
			assert(he->prev()->dst() == he2->dst());
			assert(he2->prev()->dst() == he->dst());
			assert(he->dst() == he2->src());
			assert(he2->dst() == he->src());
			assert(he->twin()->dst() == he->prev()->dst());
			assert(he2->twin()->dst() == he2->prev()->dst());
			//         connections.erase ( cit );
		}
	}
	assert(connections.size() <= (unsigned)numHalfEdges());

	assert(check(false));
	closeSurface();
	assert(check(true));

	std::cout << "F: " << numFaces() << " E: " << numHalfEdges() << " V: " << numVertices() << std::endl;
}

/**
* Adds the necessary hole-tagged faces to complete the manifold.
* @param connections Precomputed set of halfedges without a twin.
*/
void closeSurface(HE_Polyhedron* poly, std::map<std::pair<int, int>, HE_HalfEdge*>& connections)
{
	std::cerr << "closeSurface with " << connections.size() << " unmatched half-edges\n";
	// Add half-edges around the detected boundaries, each connected set being assigned to a face/hole
	while (!connections.empty())
	{
		std::map<std::pair<int, int>, HE_HalfEdge*>::iterator pit = connections.begin();
		std::pair<int, int> pts = pit->first;
		HE_HalfEdge* e = pit->second;

		//cerr << "processing edge <" << pts.first << ',' << pts.second << ">\n";
		const int start = pts.second; // pts.second is the dst() of a surface half-edge
		int v = start;
		std::stack<int> s;
		s.push(pts.second);
		s.push(pts.first);
		//   cerr << "v1 = " << pts.second << endl;
		//   cerr << "v0 = " << pts.first << endl;
		std::set<int> visited;
		visited.insert(pts.first);
		visited.insert(pts.second);
		assert(visited.size() == s.size());
		do
		{
			//      cerr << "'q' has " << s.size() << " elements\n";

			// Set 'e' to the next HalfEdge in the hole, crossing through pinched vertices.
			e = e->prev();
			std::set<const HE_HalfEdge*> seen;
			while (e->twin() && seen.end() == seen.find(e))
			{
				seen.insert(e);
				e = e->twin()->prev();
			}
			assert(e->face());
			if (seen.end() != seen.find(e))
				break;

			// If 'v' has been visited, pop back the stack until we find it again.
			//      v = e->dst()->index();
			v = e->prev()->dst()->index();
			//      cerr << "v = " << v << endl;
			if (visited.end() != visited.find(v))
			{
				//        cerr << "Found a repetition when 'q' has " << s.size() << " elements\n";
				// Construct a face with the required vertices
				std::vector<int> corners;
				corners.push_back(v);
				assert(s.top() != v);
				do
				{
					const int prev = s.top();
					corners.push_back(prev);
					//           assert(prev != v);
					s.pop();
					assert(visited.end() != visited.find(prev));
					visited.erase(prev);
					assert(visited.size() == s.size());

					//           cerr << "When 'q' has " << s.size() << " elements, (unreverted) corners contains: ";
					//           copy (corners.begin(), corners.end(), ostream_iterator<int>(cerr, " "));
					//           cerr << endl;
				} while (s.top() != v);
				assert(s.top() == v);
				s.pop();

				reverse(corners.begin(), corners.end());
				//        cerr << "The reverted corners are: ";
				//        copy (corners.begin(), corners.end(), ostream_iterator<int>(cerr, " "));
				//        cerr << endl;

				HE_Face* f = poly->addFace(corners);
				f->hole(true);

				// Connect the edges of the new face to their twins
				HE_HalfEdge* he = f->edge();
				const HE_HalfEdge* const sentinel = he;
				do
				{
					HE_HalfEdge* aux = he->prev();
					const int a = aux->dst()->index();
					const int b = he->dst()->index();
					assert(a != b);
					//           cerr << "Seeking " << b << ',' << a << endl;

					std::pair<int, int> ridx(b, a);
					assert(connections.end() != connections.find(ridx));
					HE_HalfEdge* t = connections[ridx];
					connections.erase(ridx);
					assert(t);
					assert(t != he);
					assert(t->dst() != he->dst());
					t->twin(he);
					he->twin(t);
					//           cerr << he->prev()->dst()->index() << ',' << he->dst()->index() << " paired to "
					//              << t->prev()->dst()->index() << ',' << t->dst()->index() << endl;

					assert(he->src()->index() == a);

					he = aux;
				} while (he != sentinel);
			}

			visited.insert(v);
			s.push(v);
			assert(visited.size() == s.size());
		} while (v != start);
		//   cerr << "In the end, 'q' has " << s.size() << " elements, with " << s.top() << " on top\n";
	}
}


void findUnmatchedHalfEdges(HE_Polyhedron* poly, std::map<std::pair<int, int>, HE_HalfEdge*>& connections)
{
	for (HE_Polyhedron::edge_iterator it = poly->eBegin(); it != poly->eEnd(); ++it)
	{
		HE_HalfEdge* he = *it;
		if (!he->twin())
		{
			HE_Vertex* Va = he->prev()->dst();
			HE_Vertex* Vb = he->dst();
			int a = Va->index();
			int b = Vb->index();
			assert(a != b);
			std::pair<int, int> idx(a, b);
			if (connections.end() != connections.find(idx))
				throw HE_exception("Repeated half-edge after processing: The half-edge data structure only supports manifolds.");

			assert(connections.end() == connections.find(std::pair<int, int>(b, a)));
			connections[idx] = he;
		}
		else
		{
			assert(he->twin()->twin() == he);
		}
	}
}


/// Adds the necessary hole-tagged faces to complete the manifold.
void HE_Polyhedron::closeSurface()
{
	// Find edges without a twin
	std::map<std::pair<int, int>, HE_HalfEdge*> connections;
	findUnmatchedHalfEdges(this, connections);
	if (!connections.empty())
		::closeSurface(this, connections);
}


HE_HalfEdge* HE_Polyhedron::edge(int from, int to)
{
	HE_Vertex* V = vertex(from);
	if (V)
		return V->edgeTo(to);
	return 0;
}

const HE_HalfEdge* HE_Polyhedron::edge(int from, int to) const
{
	const HE_Vertex* V = vertex(from);
	if (V)
		return V->edgeTo(to);
	return 0;
}


void HE_Polyhedron::boundingBox(cyPoint3f& m, cyPoint3f& M) const
{
	m = cyPoint3f(FLT_MAX, FLT_MAX, FLT_MAX);
	M = -m;

	// Find the bounding min and max values for each coordinate, each point
	for (const_vertex_iterator vit = vBegin(); vit != vEnd(); ++vit)
	{
		const cyPoint3f& v = (*vit)->position();
		for (int i = 0; i<3; i++)
		{
			if (v[i] > M[i])
				M[i] = v[i];
			if (v[i] < m[i])
				m[i] = v[i];
		}
	}
}


void HE_Polyhedron::scaleAndCenter(float radius)
{
	cyPoint3f m, M;
	boundingBox(m, M);
	cyPoint3f ctr = (m + M) / 2;
	float diff = 0;
	for (int i = 0; i<3; i++)
	{
		assert(M[i] >= m[i]);
		if (diff < M[i] - m[i])
			diff = M[i] - m[i];
	}

	const float FACTOR = radius / diff;
	for (vertex_iterator vit = vBegin(); vit != vEnd(); ++vit)
	{
		HE_Vertex* V = *vit;
		cyPoint3f v = V->position() - ctr;
		v *= FACTOR;
		V->position(v);
	}

}


void HE_Polyhedron::invertNormals()
{
	// Invert the order of the edges around the faces
	std::map<const HE_HalfEdge*, HE_HalfEdge*> next;
	std::map<const HE_HalfEdge*, HE_Vertex*> dst;
	for (edge_iterator eit = eBegin(); eit != eEnd(); ++eit)
	{
		HE_HalfEdge* e = *eit;

		assert(next.end() == next.find(e->next()));
		next[e->next()] = e;

		assert(dst.end() == dst.find(e));
		dst[e] = e->src();
	}

	// Apply such changes
	for (edge_iterator eit = eBegin(); eit != eEnd(); ++eit)
	{
		HE_HalfEdge* e = *eit;
		e->next(next[e]);

		HE_Vertex* v = e->dst();
		e->dst(dst[e]);
		v->edge(e);
	}

	assert(check());
}


bool faceIsHole(const HE_Face* f) { return f->hole(); }

int HE_Polyhedron::numHoles() const
{
	return count_if(fBegin(), fEnd(), faceIsHole);
}


bool HE_Polyhedron::check(bool holesFilled) const
{
#if 0
	assert(*vBegin());
	assert(*eBegin());
	assert(*fBegin());
	for (Polyhedron::const_vertex_iterator it = vBegin(); it != vEnd(); ++it)
	{
		assert((*it)->edge() != 0);
	}
#endif
	for (unsigned i = 0; i<_vertices.size(); i++)
	{
		const HE_Vertex* V = _vertices[i];
		assert(V->index() == (int)i);
		assert(V->check(holesFilled));

		//     cerr << i << " -> index " << _vertices[i]->index() << endl;
	}

	for (HE_Polyhedron::const_edge_iterator it = eBegin(); it != eEnd(); ++it)
	{
		const HE_HalfEdge* he = *it;
		assert(_halfEdges[he->index()] == *it);
		assert(he->check(holesFilled));
		if (holesFilled)
		{
			// No adjacent holes, please.
			assert(!he->face()->hole() || !he->twin()->face()->hole());
		}
	}

	for (HE_Polyhedron::const_face_iterator it = fBegin(); it != fEnd(); ++it)
	{
		std::set<const HE_HalfEdge*> visited;
		const HE_Face* f = *it;
		assert(f == _faces[f->index()]);
		HE_Face::const_edge_circulator eit = f->begin();
		HE_Face::const_edge_circulator sentinel = eit;
		do
		{
			assert((*eit)->face() == f);
			assert(visited.end() == visited.find(*eit));
			visited.insert(*eit);
			++eit;
		} while (eit != sentinel);
	}

	return true;
}

/// Tells if the Polyhedron contains a given Vertex
bool HE_Polyhedron::contains(const HE_Vertex* V) const
{
	if (V->index() >= numVertices())
		return false;
	return vertex(V->index()) == V;
}

/// Tells if the Polyhedron contains a given Face
bool HE_Polyhedron::contains(const HE_Face* F) const
{
	if (F->index() >= numFaces())
		return false;
	return face(F->index()) == F;
}

/// Tells if the Polyhedron contains a given HalfEdge
bool HE_Polyhedron::contains(const HE_HalfEdge* he) const
{
	if (he->index() >= numHalfEdges())
		return false;
	return halfedge(he->index()) == he;
}

HE_Face* HE_Polyhedron::_face(const std::vector<int>& corners) const
{
	assert(corners.size() >= 2);
	const HE_HalfEdge* he = edge(corners[0], corners[1]);
	if (!he)
		return 0;
	return (HE_Face*)he->face();
}

void HE_Polyhedron::setEdgeFlags(std::vector<std::vector<bool>> edgeFlags)
{
	for (int fi = 0; fi < numFaces(); fi++)
	{
		HE_Face* f = face(fi);
		if (f->hole()) continue;
		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;
		int ei = 0;
		do
		{
			HE_HalfEdge* he = halfedge((*e)->index());
			if (edgeFlags[fi][ei])
				he->SetFlag(EFLAG_HORIZONTAL);
			else
				he->SetFlag(EFLAG_VERTICAL);
			++e;
			++ei;
		} while (e != sentinel);
	}
}