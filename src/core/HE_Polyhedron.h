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

/** @file
* Declaration of classes Polyhedron and Polyhedron_GL, enumeration filetype_e and other functions.
*/

#ifndef __POLYHEDRON_H__
#define __POLYHEDRON_H__

#include <cassert>
#include <vector>
#include <map>
#include <set>
#include <fstream>

#include "HE_HalfEdge.h"
#include "HE_Face.h"
#include "HE_Vertex.h"

#include "../common.h"

/** enum that identifies the different types of mesh file formats supported. Only a few are actually supported.
*/
enum filetype_e
{
	FT_PLY_ASCII, ///< PLY ASCII format.
	FT_PLY_BIN, ///< PLY binary format (currently unsupported).
	FT_OFF, ///< Object File Format format.
	FT_OBJ, ///< OBJ format by Alias (currently unsupported).
	FT_TRI, ///< Simple non-indexed triangle format. Troublesome, because of the lack of indices (everything is boundary).
	FT_NONE ///< Not identified format. Indicates an error situation.
};

#ifndef INVALID
# define INVALID -1
#endif
#define ISVALID(v) (v>=0)

/** Class that holds the representation of a manifold mesh in a half-edge data structure.
*  It contains a set of Vertex, HalfEdge and Face objects that together represent the mesh.
* I/O functionality, get/set methods and other simple manipulation routines are directly implemented in the class.
* More advanced algorithms are implemented through external functions.
* Note that HE_exception objects may be thrown if something goes wrong inside this class.
*/
class HE_Polyhedron
{
public:
	/// Default constructor. Creates an empty Polyhedron.
	HE_Polyhedron();

	/// Copy constructor. Ignores isolated vertices in 'p'.
	HE_Polyhedron(const HE_Polyhedron& p);

	/**
	*    Constructor for Polyhedron class.
	* @param filename File from where the mesh is to be read. Must be one of the supported formats.
	* @return A valid Polyhedron object initialized with the data in the file.
	*/
	explicit HE_Polyhedron(const char* filename);

	explicit HE_Polyhedron(const MatrixXf &pPos, std::vector<std::vector<uint32_t>> &pIndices);

	/**
	*    Constructor that takes an array of vertices and an array of faces.
	* @param verts Set of Vertex objects to be copied onto the Polyhedron.
	* @param faces Set of faces to be loaded onto the Polyhedron.
	* @return A valid Polyhedron initialized with the indicated data.
	*/
	explicit HE_Polyhedron(const std::vector<HE_Vertex>& verts, const std::vector<std::vector<int> >& faces);

	/**
	* Constructor that duplicates the faces from some other Polyhedron, adjacent vertices too.
	* @param faces Set of faces to be duplicated.
	* @return A valid Polyhedron initialized with the indicated data.
	*/
	explicit HE_Polyhedron(const std::set<const HE_Face*>& faces);

	/**
	* Destructor. Deletes all the used memory.
	*/
	~HE_Polyhedron();

	/**
	*    Loads a mesh from a file.
	* Clears the current contents of the Polyhedron and loads the mesh described in the indicated file.
	* @param filename File with the mesh to be loaded.
	* @return True in case of success; false otherwise.
	*/
	bool load(const char* filename);

	/**
	*    Loads a mesh from a PLY ASCII file.
	* Clears the current contents of the Polyhedron and loads the mesh described in the indicated file, which should be of format PLY ASCII.
	* @param filename File with the mesh to be loaded.
	* @return True in case of success; false otherwise.
	*/
	bool readPly(const char* filename);

	/**
	*    Loads a mesh from a OBJ ASCII file.
	* Clears the current contents of the Polyhedron and loads the mesh described in the indicated file, which should be of format PLY ASCII.
	* @param filename File with the mesh to be loaded.
	* @return True in case of success; false otherwise.
	*/
	bool readObj(const char* filename);

	/**
	* Clears all the data in the mesh (vertices, faces...).
	*/
	void clearData();

	/**
	* Adds a new vertex to the Polyhedron. Its adjacency information is empty.
	* @param x X coordinate of the position of the vertex.
	* @param y Y coordinate of the position of the vertex.
	* @param z Z coordinate of the position of the vertex.
	* @return A pointer to the newly created Vertex.
	*/
	HE_Vertex* addVertex(float x, float y, float z) { return addVertex(cyPoint3f(x, y, z)); }

	/**
	* Adds a new vertex to the Polyhedron. Its adjacency information is empty.
	* @param v Position of the vertex to be added.
	* @return A pointer to the newly created Vertex.
	*/
	HE_Vertex* addVertex(const cyPoint3f& v);


	/**
	* Adds a new face to the Polyhedron. Its adjacency information (esp. some _twin links) may be empty.
	* @param corners
	* @return A pointer to the newly created Face.
	*/
	HE_Face* addFace(const std::vector<int>& corners);

	/**
	* Adds a new triangular face to the Polyhedron. Its adjacency information (esp. some _twin links) may be empty.
	* @param a First corner of the triangle.
	* @param b Second corner of the triangle.
	* @param c Third corner of the triangle.
	* @return The newly added Face.
	*/
	HE_Face* addFace(int a, int b, int c);

	/// Bulk-adding vertices
	void loadVertices(const std::vector<float>& verts);

	/// Bulk-adding faces
	void loadFaces(const std::vector<std::vector<int> >& faces);

	/**
	*    Tells the number of vertices in the Polyhedron.
	* @return Number of vertices in the Polyhedron.
	*/
	int numVertices() const { return (int)_vertices.size(); }

	/**
	*    Tells the number of faces in the Polyhedron.
	* @return Number of faces in the Polyhedron.
	*/
	int numFaces() const { return (int)_faces.size(); }

	/**
	*    Tells the number of half-edges in the Polyhedron.
	* @return Number of half-edges in the Polyhedron.
	*/
	int numHalfEdges() const { return (int)_halfEdges.size(); }

	/**
	*    Tells the number of holes in the Polyhedron.
	*  Some faces are tagged as being holes, thus not representing any surface, but a gap in it. This method counts the number of such faces.
	* @return Number of holes in the Polyhedron.
	*/
	int numHoles() const;

	/**
	* Tells if two vertices are adjacent.
	* @param v1 First Vertex
	* @param v2 Second Vertex
	* @return True if they're adjacent; false otherwise.
	*/
	bool adjacent(const HE_Vertex* v1, const HE_Vertex* v2) const { return v1->adjacent(v2); }

	/**
	* Tells if two faces are adjacent.
	* @param f1 First Face
	* @param f2 Second Face
	* @return True if they're adjacent; false otherwise.
	*/
	bool adjacent(const HE_Face* f1, const HE_Face* f2) const { return f1->adjacent(f2); }

	/**
	*    Tells the format of the mesh.
	* @return A filetype_e value corresponding to the current format of the mesh.
	* This is useful for knowing the type of file from where the mesh was read, and to properly save it.
	*/
	filetype_e type() const { return _type; }

	/**
	*    Sets the format of the mesh. This is only relevant if the mesh is going to be saved to file.
	* @param ft filetype_e value that indicates the desired format.
	*/
	void type(filetype_e ft) { _type = ft; }

	/**
	* Finds the HalfEdge that connects two vertices.
	* These vertices are indicated by their indices, and should be within range.
	* @param from Index of the vertex from which the searched HalfEdge departs.
	* @param to Index of the vertex where the searched HalfEdge ends.
	* @return The desired HalfEdge, or 0 if it does not exist.
	*/
	HE_HalfEdge* edge(int from, int to);
	HE_HalfEdge* edge(const HE_Vertex* vFrom, const HE_Vertex* vTo) { return edge(vFrom->index(), vTo->index()); }
	const HE_HalfEdge* edge(int from, int to) const;
	const HE_HalfEdge* edge(const HE_Vertex* vFrom, const HE_Vertex* vTo) const { return edge(vFrom->index(), vTo->index()); }

	/// Returns the indicated Vertex.
	HE_Vertex* vertex(int v) { return (v >= 0 && v< (int)_vertices.size()) ? _vertices[v] : 0; }

	/// Returns the indicated Vertex.
	const HE_Vertex* vertex(int v) const { return (v >= 0 && v< (int)_vertices.size()) ? _vertices[v] : 0; }

	/// Returns the indicated Face.
	HE_Face* face(int f) { return (f >= 0 && f< (int)_faces.size()) ? _faces[f] : 0; }

	/// Returns the indicated Face.
	const HE_Face* face(int f) const { return (f >= 0 && f< (int)_faces.size()) ? _faces[f] : 0; }

	/// Returns the indicated Face.
	const HE_Face* face(const std::vector<int>& corners) const { return _face(corners); }

	/// Returns the indicated Face.
	HE_Face* face(const std::vector<int>& corners) { return _face(corners); }

	/// Returns the indicated HalfEdge.
	HE_HalfEdge* halfedge(int h) { return (h >= 0 && h< (int)_halfEdges.size()) ? _halfEdges[h] : 0; }

	/// Returns the indicated HalfEdge.
	const HE_HalfEdge* halfedge(int h) const { return (h >= 0 && h< (int)_halfEdges.size()) ? _halfEdges[h] : 0; }

	const std::vector<HE_Vertex*>& vertices() const { return _vertices; } ///<Gets the array of vertices being used.

	void invertNormals();												/// Inverts all the normals of the mesh, by swapping the direction of the HalfEdge objects.
	
	void setEdgeFlags(std::vector<std::vector<bool>> edgeFlags);

	/**
	*    Centers the mesh around the origin, and scales it to fit in a bounding box of the indicated size.
	* @param side Side of the desired bounding box after calling this method. Optional parameter.
	*/
	void scaleAndCenter(float side = 10);

	/**
	*    Computes the (axis aligned) bounding box of the mesh.
	* @param m Coordinates of the lower corner of the bounding box.
	* @param M Coordinates of the upper corner of the bounding box.
	*/
	void boundingBox(cyPoint3f& m, cyPoint3f& M) const;

	/// Computes the diagonal of the bounding box of the Polyhedron.
	float diagonal() const { cyPoint3f m, M; boundingBox(m, M); return float((m - M).Length()); }

	/**
	* Makes the _twin connections for all halfedges, and closes the holes.
	*/
	void finalize();

	// Iterators/circulators
#include "polyhedron_iterators.h"

	/**
	* Verifies that the data structure is internally consistent.
	* @param holesFilled Have the holes been filled yet?
	* @return True if everything was fine; False otherwise.
	*/
	bool check(bool holesFilled = true) const;

	bool contains(const HE_Vertex* V) const;
	bool contains(const HE_Face* F) const;
	bool contains(const HE_HalfEdge* he) const;

protected:
	void closeSurface();
	bool readPlyHeader(std::ifstream& in, int& nv, int& nf, short order[11]);
	bool readPlyData(std::ifstream& in, int& nv, int& nf, const short order[11], std::vector<float>& verts, std::vector<std::vector<int> >& faces);
	HE_Face* _face(const std::vector<int>& corners) const;

	std::vector<HE_Face*> _faces;
	std::set<HE_Face*> _holes;
	std::vector<HE_HalfEdge*> _halfEdges;
	std::vector<HE_Vertex*> _vertices;
	filetype_e _type;
};

/**
*    Opens the indicated file and checks magic strings and other information in order to determine the file format.
* @param filename File to be checked.
* @return An indicator of the file's format, or FT_NONE, if no supported format is identified.
*/
filetype_e fileFormat(const char* filename);

/**
* Determines the file format given its name.
* @param filename File name to be checked.
* @return An indicator of the file's format, or FT_NONE, if no supported format is identified.
*/
filetype_e fileFormatFromFileName(const char* filename);

#endif // __POLYHEDRON_H__
