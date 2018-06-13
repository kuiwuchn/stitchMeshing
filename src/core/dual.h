/*********************************************************************************
*   Copyright (C) 2018 by Kui Wu                                                 *
*   kwu@cs.utah.edu                                                              *
*                                                                                *
*   This program is free software: you can redistribute it and/or modify         *
*   it under the terms of the GNU Lesser General Public License as published by  *
*   the Free Software Foundation, either version 3 of the License, or	         *
*   (at your option) any later version.									         *
*   																	         *
*   This program is distributed in the hope that it will be useful,		         *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of		         *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the		         *
*   GNU Lesser General Public License for more details.						     *
*   																	         *
*   You should have received a copy of the GNU Lesser General Public License	 *
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.        *
*********************************************************************************/

#ifndef __Single_Strip_Mesh_Simplification__Dual__
#define __Single_Strip_Mesh_Simplification__Dual__

#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
#include <string> 

#include "cyPoint.h"
#include "HE_Polyhedron.h"

#include "gurobi_c++.h"

class HE_Polyhedron;

struct DualNode {
	/* A vertex in the dual graph, */
	/* which is actually a face in the mesh. */

	/* indices into the initial vertex array */
	std::vector<int> _v;			/// this ordering is important for culling calculation
	std::vector<int> _e;			/// This ordering is not particularly important
	cyPoint3f _fn;					/// The face normal
	cyPoint3f _fc;					/// The face center
	bool _labeled;
	int _groupId;					/// group id

	/* Convenience methods */
	cyPoint3f normal() const { return _fn; }
	int vertex(int i) const { return _v[i]; };
	int edge(int i) const { return _e[i]; };
};

enum EDGE_TYPE
{
	EG_NONE = 0,
	EG_TOP,
	EG_BTM,
};

struct GroupNode
{
	std::vector<int>		_topGroupIds;		/// top group ids
	std::vector<int>		_topGroupCounts;
	std::vector<EDGE_TYPE>	_topEdgeType;

	std::vector<int>		_btmGroupIds;		/// bottom group ids
	std::vector<int>		_btmGroupCounts;
	std::vector<EDGE_TYPE>	_btmEdgeType;

	bool					_needFlip;
};

struct GroupEdge
{
	int			_fi[2];
	EDGE_TYPE	_eType[2];
	int			_weight;
};

struct DualEdge {
	/* An edge in the dual graph */
	/* Ordering unimportant */
	/* These are the nodes of the graph (faces in the mesh) */
	int _idx;
	int _fi[2];
	float _wt;		/// The weight
	int _hei[2];		/// The index of twin halfedges
};

struct avertex {
	cyPoint3f _coords;	/// 3D coordinates of the vertex
	cyPoint3f _norm;	/// The normals of the vertex
	cyPoint3f coords() const { return _coords; } 
	cyPoint3f normal() const { return _norm; }
};

struct virtualQuad
{
	int _fi[2];
	int _he[4];
	float _wt;
};

enum MISMATCH_TYPE
{
	MM_NONE = 0,
	MM_WALE_COURSE,
	MM_WALE,
};

/* The Dual graph */
class DualGraph {
public:
	std::vector<MISMATCH_TYPE>		_mismatchTypes;
	std::vector<int>				_idxOffsets;
	std::vector<std::vector<bool>>	_gurobiResultPool;
	int								_nSolutions;
	std::vector<std::vector<int>>	_groupFaceIdx;
	std::vector<bool>				_isLoop;
	std::vector<GroupNode>			_groupnodes;
	std::vector<GroupEdge>			_groupedges;

	std::vector<bool>				_isBadFaces;
	std::vector<int>				_badAngleIdx;

	std::vector<std::tuple<int, int, int>>	_trianglePairs;
	std::vector<virtualQuad>		_virtualQuads;

	friend class HE_Polyhedron;

	DualGraph(HE_Polyhedron* ply);
	~DualGraph();

	void constructDualGraph(HE_Polyhedron* ply);

	const avertex& vertex(int i) const	{ return _avertices[i]; }
	const DualNode& node(int i) const	{ return _dualnodes[i]; }
	const DualEdge& edge(int i) const	{ return _dualedges[i]; }
	float weight(int i) const			{ return _dualedges[i]._wt; }

	int numVertices() const { return int(_avertices.size()); }
	int numNodes() const	{ return int(_dualnodes.size()); }
	int numEdges() const	{ return int(_dualedges.size()); }

	void gurobiSolver(std::string pFilename);
	void loadGurobiResult(std::string pFilename);

	void gurobiSolverHori();

	void cutUVMismatchQuad(std::vector<HE_Vertex>& pVerts, std::vector<std::vector<int> >& pFaces, std::vector<std::vector<bool>>& edgeFlags);

	void findUVMismatch();
	void findWaleMismatch();
	void findFaceGroups();

	void findFaceGroups_v2();

	void findLoops();

	void subdivide(std::vector<HE_Vertex>& pVerts, std::vector<std::vector<int> >& pFaces);

	void cut(std::vector<HE_Vertex>& pVerts, std::vector<std::vector<int> >& pFaces);

	void fitFaceEdgesOrder();

	void labelQuadFace();
	void labelTriangleFace();

	void findTopBtmGroupId();

	void createGroupEdges();

	void waleMismatchSolver();

	void setSubGroupFace(DualGraph* pDual);

	void setGroupEdgeFlag();

	void reorderSubdivisionGroupFaceIdx(DualGraph* pDual);

	void findBadQuads();

	void mergeTriangles(std::vector<HE_Vertex>& pVerts, std::vector<std::vector<int> >& pFaces);

	void mergeTriangles(std::vector<HE_Vertex>& pVerts, std::vector<std::vector<int> >& pFaces, std::vector<std::vector<bool>>& edgeFlags);

	void removeBadVertex(std::vector<HE_Vertex>& pVerts, std::vector<std::vector<int> >& pFaces);

	void flipTriangles(std::vector<HE_Vertex>& pVerts, std::vector<std::vector<int> >& pFaces, std::vector<std::vector<bool>>& edgeFlags);

	void breakConnectIncOrDec();

	void HVCheck();

	void BadVertexCheck();

	void fixUVMismatch();

	bool CheckBadVertex();

	void flipTriangleUV();
	void flipEdgeUV();

	void fixTriangleUV(int fi, int i);

	void exportQuadDominantMesh(const char* filename);

	void exportUVMesh(const char* filename);

	void exportTextureMesh(const char* filename);

protected:
	bool addEdgeIndexToDualnode(int i, int f);

	struct {
		bool operator()(DualEdge a, DualEdge b) const { return a._wt < b._wt; }
	} customLess;

	std::vector<avertex>	_avertices;
	std::vector<DualNode>	_dualnodes;
	std::vector<DualEdge>	_dualedges;

	HE_Polyhedron *_poly;
};

#endif 