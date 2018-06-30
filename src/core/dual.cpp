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

#include "Dual.h"

DualGraph::DualGraph(HE_Polyhedron* ply) {
    _poly = ply;
    /* Construct the dual graph from ply */
    constructDualGraph(ply);
}


DualGraph::~DualGraph() {
    _avertices.clear();
    _dualnodes.clear();
    _dualedges.clear();
}

/* Use the connectivity information already built in ply to create the dual graph */
void DualGraph::constructDualGraph(HE_Polyhedron* ply) {
    /* Store the vertex information locally in a vector */
    for(int i=0; i< ply->numVertices(); i++)
    {
        avertex v;
        v._coords = ply->vertex(i)->position();        
        v._norm = ply->vertex(i)->normal();       
        _avertices.push_back(v);
    }
    
    /* Iterate over the nodes */
    for(int i=0; i<ply->numFaces(); i++)
    {
        DualNode n;
		const HE_Face* f = ply->face(i);
		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;
		do
		{
			const HE_HalfEdge* he = *e;
			n._v.push_back(he->dst()->index());
			++e;
		} while (e != sentinel);
        
		n._e.resize(n._v.size()); 
		for (int j = 0; j < (int)n._v.size(); j++) n._e[j] = -1;

        n._fn = ply->face(i)->normal();
		n._fc = ply->face(i)->centroid();
		n._labeled = false;
        _dualnodes.push_back(n);
    }
    
    int n = 0;
	cyPoint3f n1, n2;
    for(int i=0; i< ply->numHalfEdges(); i++)
    {
        DualEdge e;
        
        if (ply->halfedge(i)->face()->index() < ply->halfedge(i)->twin()->face()->index()) {
            e._fi[0] = ply->halfedge(i)->face()->index();
            e._fi[1] = ply->halfedge(i)->twin()->face()->index();       
            e._hei[0] = ply->halfedge(i)->index();
            e._hei[1] = ply->halfedge(i)->twin()->index();

			n1 = ply->halfedge(i)->face()->normal();
			n2 = ply->halfedge(i)->twin()->face()->normal();

			e._wt = 0.25f - 0.25f * n1.Dot(n2);
			e._idx = n;

            _dualedges.push_back(e);
            
            addEdgeIndexToDualnode(n, ply->halfedge(i)->face()->index());
            addEdgeIndexToDualnode(n, ply->halfedge(i)->twin()->face()->index());

			n++;
        }
    }

	for (int i = 0; i < numEdges(); i++)
		_mismatchTypes.push_back(MM_NONE);
}

bool DualGraph::addEdgeIndexToDualnode(int i, int fi) {
	for (int j = 0; j < int(_dualnodes[fi]._e.size()); j++)
		if (_dualnodes[fi]._e[j] == -1) { _dualnodes[fi]._e[j] = i; return true; }
	return false;
}

void DualGraph::cutUVMismatchQuad(std::vector<HE_Vertex>& pVerts, std::vector<std::vector<int> >& pFaces, std::vector<std::vector<bool>>& pEdgeFlags)
{
	int nv = _poly->numVertices();
	for (int vi = 0; vi < nv; vi++)
	{
		HE_Vertex v = HE_Vertex(_poly->vertex(vi)->position(), vi);
		pVerts.push_back(v);
	}

	for (int fi = 0; fi < _poly->numFaces(); fi++)
	{
		HE_Face* f = _poly->face(fi);
		if (f->hole()) continue;
		std::vector<int> cornerList;
		std::vector<bool> edgeFlag;
		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;

		bool isBadFace = false;
		bool typeA = false;
		int badAngleIdx;
		if (f->NumHalfEdge() == 4)
		{
			do
			{
				if ((*e)->TestFlag(EFLAG_HORIZONTAL) && (*e)->twin()->TestFlag(EFLAG_VERTICAL))
				{
					typeA = true;	isBadFace = true;	break;
				}
				else if ((*e)->TestFlag(EFLAG_VERTICAL) && (*e)->twin()->TestFlag(EFLAG_HORIZONTAL) && (*e)->twin()->face()->NumHalfEdge() == 3)
				{
					typeA = false;	isBadFace = true;	break;
				}
				++e;
			} while (e != sentinel);

			if (isBadFace)
			{
				e = f->begin();
				e--;
				sentinel = e;
				HE_Face::const_edge_circulator nxtE = e;
				nxtE--;
				cyPoint3f t0, t1;
				int i = 0;
				float worstAngle = 0;
				float dotValue;
				do
				{
					t0 = (*e)->tangent();		t0.Normalize();
					t1 = (*nxtE)->tangent();	t1.Normalize();
					dotValue = fabsf(t0.Dot(t1));
	
					if (dotValue > worstAngle)
					{
						worstAngle = dotValue;
						badAngleIdx = i;
					}					

					++e; ++nxtE; ++i;
				} while (e != sentinel);
			}
		}

		e = f->begin();
		sentinel = e;
		if (!isBadFace)
		{
			do
			{
				cornerList.push_back((*e)->dst()->index());
				edgeFlag.push_back((*e)->TestFlag(EFLAG_HORIZONTAL));
				++e;
			} while (e != sentinel);

			pFaces.push_back(cornerList);
			pEdgeFlags.push_back(edgeFlag);
		}
		else
		{
			int idx = badAngleIdx;
			if (idx == 0 || idx == 2)
			{
				cornerList.push_back((*e)->dst()->index());
				edgeFlag.push_back(typeA);
				e++;
				cornerList.push_back((*e)->dst()->index());
				edgeFlag.push_back((*e)->twin()->TestFlag(EFLAG_HORIZONTAL));
				e++;
				cornerList.push_back((*e)->dst()->index());
				edgeFlag.push_back((*e)->twin()->TestFlag(EFLAG_HORIZONTAL));
				pFaces.push_back(cornerList);
				pEdgeFlags.push_back(edgeFlag);

				cornerList.clear();
				edgeFlag.clear();

				cornerList.push_back((*e)->dst()->index());
				edgeFlag.push_back(typeA);
				e++;
				cornerList.push_back((*e)->dst()->index());
				edgeFlag.push_back((*e)->twin()->TestFlag(EFLAG_HORIZONTAL));
				e++;
				cornerList.push_back((*e)->dst()->index());
				edgeFlag.push_back((*e)->twin()->TestFlag(EFLAG_HORIZONTAL));
				pFaces.push_back(cornerList);
				pEdgeFlags.push_back(edgeFlag);
			}
			else if (idx == 1 || idx == 3)
			{
				e++;
				cornerList.push_back((*e)->dst()->index());
				edgeFlag.push_back(typeA);
				e++;
				cornerList.push_back((*e)->dst()->index());
				edgeFlag.push_back((*e)->twin()->TestFlag(EFLAG_HORIZONTAL));
				e++;
				cornerList.push_back((*e)->dst()->index());
				edgeFlag.push_back((*e)->twin()->TestFlag(EFLAG_HORIZONTAL));
				pFaces.push_back(cornerList);
				pEdgeFlags.push_back(edgeFlag);

				cornerList.clear();
				edgeFlag.clear();
				cornerList.push_back((*e)->dst()->index());
				edgeFlag.push_back(typeA);
				e++;
				cornerList.push_back((*e)->dst()->index());
				edgeFlag.push_back((*e)->twin()->TestFlag(EFLAG_HORIZONTAL));
				e++;
				cornerList.push_back((*e)->dst()->index());
				edgeFlag.push_back((*e)->twin()->TestFlag(EFLAG_HORIZONTAL));
				pFaces.push_back(cornerList);
				pEdgeFlags.push_back(edgeFlag);
			}
		}
	}
}

void DualGraph::findUVMismatch()
{
	int uvMismatchCount = 0;
	HE_HalfEdge* he1; HE_HalfEdge* he2;
	for (int i = 0; i < numEdges(); i++)
	{
		he1 = _poly->halfedge(_dualedges[i]._hei[0]);
		he2 = _poly->halfedge(_dualedges[i]._hei[1]);

		if (he1->TestFlag(EFLAG_VERTICAL) && he2->TestFlag(EFLAG_HORIZONTAL))
		{
			uvMismatchCount++;
			_mismatchTypes[i] = MM_WALE_COURSE;
		}

		if (he2->TestFlag(EFLAG_VERTICAL) && he1->TestFlag(EFLAG_HORIZONTAL))
		{
			uvMismatchCount++;
			_mismatchTypes[i] = MM_WALE_COURSE;
		}
	}
	std::cout << "UV mismatch number " << uvMismatchCount << std::endl;
}

bool vqSort(const virtualQuad &a, const virtualQuad &b)
{
	return a._wt < b._wt;
}

void DualGraph::removeBadVertex(std::vector<HE_Vertex>& pVerts, std::vector<std::vector<int> >& pFaces)
{
	std::vector<bool> removedFaceIdx;
	for (int fi = 0; fi < _poly->numFaces(); fi++) removedFaceIdx.push_back(false);

	std::vector<int> mappingFromOriToNew;
	int i = 0;
	for (int vi = 0; vi < _poly->numVertices(); vi++)
	{
		HE_Vertex* v = _poly->vertex(vi);
		if (v->degree() <= 3)
		{
			HE_Vertex::const_edge_circulator e = v->begin();
			HE_Vertex::const_edge_circulator sentinel = e;

			do
			{
				std::cout << (*e)->dst()->index() << " ";
				removedFaceIdx[(*e)->face()->index()] = true;
				++e;
			} while (e != sentinel);
			std::cout << std::endl;
			mappingFromOriToNew.push_back(-1);
		}
		else
			mappingFromOriToNew.push_back(i++);
	}

	int nv = _poly->numVertices();
	for (int vi = 0; vi < nv; vi++)
	{
		HE_Vertex v = HE_Vertex(_poly->vertex(vi)->position(), mappingFromOriToNew[vi]);
		if (_poly->vertex(vi)->degree() > 3)
		{
			pVerts.push_back(v);
		}
	}

	for (int fi = 0; fi < _poly->numFaces(); fi++)
	{
		HE_Face* f = _poly->face(fi);
		std::vector<int> cornerList;
		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;
		if (!removedFaceIdx[fi])
		{
			do
			{
				cornerList.push_back(mappingFromOriToNew[(*e)->dst()->index()]);
				++e;
			} while (e != sentinel);
			pFaces.push_back(cornerList);
		}
	}

	for (int vi = 0; vi < _poly->numVertices(); vi++)
	{
		std::vector<int> cornerList;
		HE_Vertex* v = _poly->vertex(vi);
		if (v->degree() <= 3)
		{
			HE_Vertex::const_edge_circulator e = v->begin();
			HE_Vertex::const_edge_circulator sentinel = e;

			do
			{
				//std::cout << mappingFromOriToNew[(*e)->dst()->index()] << " ";
				cornerList.push_back(mappingFromOriToNew[(*e)->dst()->index()]);

				if ((*e)->face()->NumHalfEdge() == 4) cornerList.push_back(mappingFromOriToNew[(*e)->next()->dst()->index()]);

				e--;
			} while (e != sentinel);
			pFaces.push_back(cornerList);
		}
	}

	//for (int ci = 0; ci < (int)cornerLists.size(); ci++)
	//{
	//	std::vector<int> cornerList;
	//	for (int cii = 0; cii < (int)cornerLists[ci].size(); cii++)
	//	{
	//		cornerList.push_back(mappingFromOriToNew[cornerLists[ci][cii]]);
	//	}
	//	pFaces.push_back(cornerList);
	//}
}

void DualGraph::mergeTriangles(std::vector<HE_Vertex>& pVerts, std::vector<std::vector<int> >& pFaces, std::vector<std::vector<bool>>& pEdgeFlags)
{
	/* find triangles pair */
	for (int fi = 0; fi < (int)_dualnodes.size(); fi++)
	{
		HE_Face* f = _poly->face(fi);
		if (f->NumHalfEdge() == 3)
		{
			int added = false;
#if 0
			if (f->TestFlag(VEF_FLAG_SHORTROW_L) || f->TestFlag(VEF_FLAG_SHORTROW_R))
			{
				HE_HalfEdge* he = f->edge();
				while(he->TestFlag(EFLAG_HORIZONTAL))
				{
					he = he->next();
				}
				
				if (he->twin()->face()->TestFlag(VEF_FLAG_INC) || he->twin()->face()->TestFlag(VEF_FLAG_DEC))
				{
					if (he->twin()->next()->TestFlag(EFLAG_VERTICAL) && 
						(he->next()->twin()->face()->TestFlag(VEF_FLAG_SHORTROW_L) || he->next()->twin()->face()->TestFlag(VEF_FLAG_SHORTROW_R)))
					{
						if (he->TestFlag(EFLAG_VERTICAL) && he->prev()->TestFlag(EFLAG_VERTICAL))
						{
						}
						else if (he->next()->twin()->next()->TestFlag(EFLAG_VERTICAL) && he->next()->twin()->prev()->TestFlag(EFLAG_VERTICAL))
						{
						}
						else if (he->dst()->degree() <= 4 || he->src()->degree() <= 4)
						{
						}
						else
						{
							he->next()->ClearFlags();
							he->next()->twin()->ClearFlags();
							he->next()->SetFlag(EFLAG_VERTICAL);
							he->next()->twin()->SetFlag(EFLAG_VERTICAL);

							_trianglePairs.push_back(std::tuple<int, int, int>(fi, he->twin()->face()->index(), he->index()));

							added = true;
						}
					}
					else if (he->twin()->prev()->TestFlag(EFLAG_VERTICAL) &&
						(he->prev()->twin()->face()->TestFlag(VEF_FLAG_SHORTROW_L) || he->prev()->twin()->face()->TestFlag(VEF_FLAG_SHORTROW_R)))
					{
						if (he->TestFlag(EFLAG_VERTICAL) && he->next()->TestFlag(EFLAG_VERTICAL))
						{
						}
						else if (he->prev()->twin()->next()->TestFlag(EFLAG_VERTICAL) && he->prev()->twin()->prev()->TestFlag(EFLAG_VERTICAL))
						{
						}
						else if (he->dst()->degree() <= 4 || he->src()->degree() <= 4)
						{
						}
						else
						{
							he->prev()->ClearFlags();
							he->prev()->twin()->ClearFlags();
							he->prev()->SetFlag(EFLAG_VERTICAL);
							he->prev()->twin()->SetFlag(EFLAG_VERTICAL);

							_trianglePairs.push_back(std::tuple<int, int, int>(fi, he->twin()->face()->index(), he->index()));

							added = true;
						}
					}
				}
			}
#endif
			if (!added)
			{
				HE_Face::const_edge_circulator e = f->begin();
				HE_Face::const_edge_circulator sentinel = e;
				do
				{
					const HE_Face* twinF = (*e)->twin()->face();
					if (twinF->NumHalfEdge() == 3 && twinF->index() < fi)
					{
						if ((*e)->next()->TestFlag(EFLAG_HORIZONTAL) && (*e)->prev()->TestFlag(EFLAG_VERTICAL) &&
							(*e)->twin()->next()->TestFlag(EFLAG_HORIZONTAL) && (*e)->twin()->prev()->TestFlag(EFLAG_VERTICAL))
						{
							if ((*e)->dst()->degree() > 4 && (*e)->src()->degree() > 4)
								_trianglePairs.push_back(std::tuple<int, int, int>(fi, twinF->index(), (*e)->index()));
						}
						else if ((*e)->next()->TestFlag(EFLAG_VERTICAL) && (*e)->prev()->TestFlag(EFLAG_HORIZONTAL) &&
							(*e)->twin()->next()->TestFlag(EFLAG_VERTICAL) && (*e)->twin()->prev()->TestFlag(EFLAG_HORIZONTAL))
						{
							if ((*e)->dst()->degree() > 4 && (*e)->src()->degree() > 4)
								_trianglePairs.push_back(std::tuple<int, int, int>(fi, twinF->index(), (*e)->index()));
						}
					}
					++e;
				} while (e != sentinel);
			}
		}
	}
#if 1
	/* create virtual quads */
	_virtualQuads.clear();
	HE_HalfEdge* tmpHe0;
	HE_HalfEdge* tmpHe1;
	virtualQuad vq;
	cyPoint3f t0, t1;
	float dotValue;
	for (int ti = 0; ti < (int)_trianglePairs.size(); ti++)
	{
		HE_HalfEdge* he0 = _poly->halfedge(std::get<2>(_trianglePairs[ti]));
		HE_HalfEdge* he1 = he0->twin();

		vq._fi[0] = he0->face()->index();
		vq._fi[1] = he1->face()->index();
		vq._he[0] = he0->next()->index();
		vq._he[1] = he0->next()->next()->index();
		vq._he[2] = he1->next()->index();
		vq._he[3] = he1->next()->next()->index();

		float maxDot = 0;
		for (int i = 0; i < 4; i++)
		{
			tmpHe0 = _poly->halfedge(vq._he[i]);
			if (i == 3)
				tmpHe1 = _poly->halfedge(vq._he[0]);
			else
				tmpHe1 = _poly->halfedge(vq._he[i + 1]);
			t0 = tmpHe0->tangent();	t0.Normalize();
			t1 = tmpHe1->tangent();	t1.Normalize();

			dotValue = fabsf(t0.Dot(t1));

			maxDot = fmaxf(maxDot, dotValue);
		}

		//if (maxDot < 0.3f)
		{
			vq._wt = maxDot;
			_virtualQuads.push_back(vq);
		}
	}

	/* sort virtual quads */
	std::sort(_virtualQuads.begin(), _virtualQuads.end(), vqSort);

	/* find merging virtual quads */
	std::vector<bool> isTriAvaible;
	std::vector<int> isTriMerged;
	for (int fi = 0; fi < (int)_dualnodes.size(); fi++)
	{
		isTriAvaible.push_back(true);
		isTriMerged.push_back(-1);
	}

	virtualQuad* vqp;
	std::vector<bool> markers;
	for (int vqi = 0; vqi < (int)_virtualQuads.size(); vqi++)
	{
		vqp = &(_virtualQuads[vqi]);
		if (isTriAvaible[vqp->_fi[0]] && isTriAvaible[vqp->_fi[1]])
		{
			isTriAvaible[vqp->_fi[0]] = false;
			isTriAvaible[vqp->_fi[1]] = false;
			markers.push_back(true);
		}
		else
		{
			markers.push_back(false);
		}
	}

	// assign vertices list
	int nv = _poly->numVertices();
	for (int vi = 0; vi < nv; vi++)
	{
		HE_Vertex v = HE_Vertex(_poly->vertex(vi)->position(), vi);
		pVerts.push_back(v);
	}

	// face idx
	for (int fi = 0; fi < _poly->numFaces(); fi++)
	{
		HE_Face* f = _poly->face(fi);

		if (f->hole()) continue;

		std::vector<int> cornerList;
		std::vector<bool> edgeFlag;
		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;
		if (isTriAvaible[fi])
		{
			do
			{
				cornerList.push_back((*e)->dst()->index());
				edgeFlag.push_back((*e)->TestFlag(EFLAG_HORIZONTAL));
				++e;
			} while (e != sentinel);
			pFaces.push_back(cornerList);
			pEdgeFlags.push_back(edgeFlag);

			//if (count == 0)
			//{
			//	for (int i = 0; i < edgeFlag.size(); i++)
			//	{
			//		std::cout << edgeFlag[i] << " ";
			//	}
			//	std::cout << std::endl;
			//}
		}
	}

	for (int vqi = 0; vqi < (int)_virtualQuads.size(); vqi++)
	{
		std::vector<int> cornerList;
		std::vector<bool> edgeFlag;
		if (markers[vqi])
		{
			vqp = &(_virtualQuads[vqi]);
			for (int i = 0; i < 4; i++)
			{
				cornerList.push_back(_poly->halfedge(vqp->_he[i])->dst()->index());
				edgeFlag.push_back(_poly->halfedge(vqp->_he[i])->TestFlag(EFLAG_HORIZONTAL));
				/*std::cout << _poly->halfedge(vqp->_he[i])->src()->index() << " "
				<< _poly->halfedge(vqp->_he[i])->dst()->index() << ", ";*/
				//std::cout << edgeFlag[i] << " ";
			}
			//std::cout << std::endl;
			//std::cout << std::endl;
			pFaces.push_back(cornerList);
			pEdgeFlags.push_back(edgeFlag);
		}
	}
#endif
}

void DualGraph::flipTriangles(std::vector<HE_Vertex>& pVerts, std::vector<std::vector<int> >& pFaces, std::vector<std::vector<bool>>& pEdgeFlags)
{
	int nv = _poly->numVertices();
	for (int vi = 0; vi < nv; vi++)
	{
		HE_Vertex v = HE_Vertex(_poly->vertex(vi)->position(), vi);
		pVerts.push_back(v);
	}

	_trianglePairs.clear();
	for (int fi = 0; fi < (int)_dualnodes.size(); fi++)
	{
		HE_Face* f = _poly->face(fi);
		if (f->NumHalfEdge() == 3)
		{
			int Hcount = 0, Vcount = 0;
			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;
			do
			{
				if ((*e)->TestFlag(EFLAG_HORIZONTAL)) Hcount++;
				if ((*e)->TestFlag(EFLAG_VERTICAL)) Vcount++;
				e++;
			} while (e != sentinel);

			if (Hcount == 0 || Hcount > 2 || Vcount == 0 || Vcount > 2)
			{
				e = f->begin();
				do
				{
					if ((*e)->twin()->face()->NumHalfEdge() == 3)
					{
						_trianglePairs.push_back(std::tuple<int, int, int>(fi, (*e)->twin()->face()->index(), (*e)->index()));
						break;
					}
					e++;
				} while (e != sentinel);
			}
		}
	}

	std::vector<bool> needRemovedMarker;
	for (int fi = 0; fi < numNodes(); fi++) needRemovedMarker.push_back(false);

	for (int i = 0; i < _trianglePairs.size(); i++)
	{
		std::cout << std::get<0>(_trianglePairs[i]) << " " << std::get<1>(_trianglePairs[i]) << std::endl;
		needRemovedMarker[std::get<0>(_trianglePairs[i])] = true;
		needRemovedMarker[std::get<1>(_trianglePairs[i])] = true;
	}

	for (int fi = 0; fi < _poly->numFaces(); fi++)
	{
		HE_Face* f = _poly->face(fi);
		if (f->hole() || needRemovedMarker[fi]) continue;
		std::vector<int> cornerList;
		std::vector<bool> edgeFlag;
		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;

		do
		{
			cornerList.push_back((*e)->dst()->index());
			edgeFlag.push_back((*e)->TestFlag(EFLAG_HORIZONTAL));
			++e;
		} while (e != sentinel);

		pFaces.push_back(cornerList);
		pEdgeFlags.push_back(edgeFlag);
	}

	int fi0, fi1, e0i;
	HE_HalfEdge* he;
	for (int i = 0; i < _trianglePairs.size(); i++)
	{
		fi0 = std::get<0>(_trianglePairs[i]);
		fi1 = std::get<1>(_trianglePairs[i]);
		e0i = std::get<2>(_trianglePairs[i]);
		he = _poly->halfedge(e0i);

		std::vector<int> cornerList;
		std::vector<bool> edgeFlag;
		cornerList.push_back(he->next()->dst()->index());
		edgeFlag.push_back(he->next()->twin()->TestFlag(EFLAG_HORIZONTAL));
		cornerList.push_back(he->twin()->next()->dst()->index());
		edgeFlag.push_back(true);
		cornerList.push_back(he->dst()->index());
		edgeFlag.push_back(he->twin()->prev()->TestFlag(EFLAG_HORIZONTAL));
		pFaces.push_back(cornerList);
		pEdgeFlags.push_back(edgeFlag);

		cornerList.clear();
		edgeFlag.clear();

		cornerList.push_back(he->twin()->next()->dst()->index());
		edgeFlag.push_back(he->twin()->next()->twin()->TestFlag(EFLAG_HORIZONTAL));
		cornerList.push_back(he->next()->dst()->index());
		edgeFlag.push_back(true);
		cornerList.push_back(he->src()->index());
		edgeFlag.push_back(he->prev()->TestFlag(EFLAG_HORIZONTAL));
		pFaces.push_back(cornerList);
		pEdgeFlags.push_back(edgeFlag);
	}
}

void DualGraph::mergeTriangles(std::vector<HE_Vertex>& pVerts, std::vector<std::vector<int> >& pFaces)
{
	/* find triangles pair */
	for (int fi = 0; fi < (int)_dualnodes.size(); fi++)
	{
		HE_Face* f = _poly->face(fi);
		if (f->NumHalfEdge() == 3)
		{
			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;
			do
			{
				const HE_Face* twinF = (*e)->twin()->face();
				if (twinF->NumHalfEdge() == 3 && twinF->index() < fi)
				{
					_trianglePairs.push_back(std::tuple<int, int, int>(fi, twinF->index(), (*e)->index()));
				}
				++e; 
			} while (e != sentinel);
		}
	}

	/* create virtual quads */
	_virtualQuads.clear();
	HE_HalfEdge* tmpHe0;
	HE_HalfEdge* tmpHe1;
	virtualQuad vq;
	cyPoint3f t0, t1;
	float dotValue;
	for (int ti = 0; ti < (int)_trianglePairs.size(); ti++)
	{
		HE_HalfEdge* he0 = _poly->halfedge(std::get<2>(_trianglePairs[ti]));
		HE_HalfEdge* he1 = he0->twin();

		vq._fi[0] = he0->face()->index();
		vq._fi[1] = he1->face()->index();
		vq._he[0] = he0->next()->index();
		vq._he[1] = he0->next()->next()->index();
		vq._he[2] = he1->next()->index();
		vq._he[3] = he1->next()->next()->index();

		float maxDot = 0;
		for (int i = 0; i < 4; i++)
		{
			tmpHe0 = _poly->halfedge(vq._he[i]);
			if (i == 3)
				tmpHe1 = _poly->halfedge(vq._he[0]);
			else 
				tmpHe1 = _poly->halfedge(vq._he[i+1]);
			t0 = tmpHe0->tangent();	t0.Normalize();
			t1 = tmpHe1->tangent();	t1.Normalize();

			dotValue = fabsf(t0.Dot(t1));

			maxDot = fmaxf(maxDot, dotValue);
		}

		if (maxDot < 0.1f) {
			vq._wt = maxDot;
			_virtualQuads.push_back(vq);
		}
	}

	/* sort virtual quads */
	std::sort(_virtualQuads.begin(), _virtualQuads.end(), vqSort);

	/* find merging virtual quads */
	std::vector<bool> isTriAvaible;
	std::vector<int> isTriMerged;
	for (int fi = 0; fi < (int)_dualnodes.size(); fi++)
	{
		isTriAvaible.push_back(true);
		isTriMerged.push_back(-1);
	}
	
	virtualQuad* vqp;
	std::vector<bool> markers;
	for (int vqi = 0; vqi < (int)_virtualQuads.size(); vqi++)
	{
		vqp = &(_virtualQuads[vqi]);
		if (isTriAvaible[vqp->_fi[0]] && isTriAvaible[vqp->_fi[1]])
		{
			isTriAvaible[vqp->_fi[0]] = false;
			isTriAvaible[vqp->_fi[1]] = false;
			markers.push_back(true);
		}
		else
		{
			markers.push_back(false);
		}
	}

	// assign vertices list
	int nv = _poly->numVertices();
	for (int vi = 0; vi < nv; vi++)
	{
		HE_Vertex v = HE_Vertex(_poly->vertex(vi)->position(), vi);
		pVerts.push_back(v);
	}

	// face idx
	for (int fi = 0; fi < _poly->numFaces(); fi++)
	{
		HE_Face* f = _poly->face(fi);
		if (f->hole()) continue;
		std::vector<int> cornerList;
		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;
		if (isTriAvaible[fi])
		{
			do
			{
				cornerList.push_back((*e)->dst()->index());
				++e;
			} while (e != sentinel);

			pFaces.push_back(cornerList);
		}
	}

	for (int vqi = 0; vqi < (int)_virtualQuads.size(); vqi++)
	{
		std::vector<int> cornerList;
		if (markers[vqi])
		{
			vqp = &(_virtualQuads[vqi]);
			for (int i = 0; i < 4; i++)
			{
				cornerList.push_back(_poly->halfedge(vqp->_he[i])->dst()->index());
				/*std::cout << _poly->halfedge(vqp->_he[i])->src()->index() << " "
					<< _poly->halfedge(vqp->_he[i])->dst()->index() << ", ";*/
			}
			//std::cout << std::endl;
			pFaces.push_back(cornerList);
		}
	}
}



void DualGraph::findBadQuads()
{
	_isBadFaces.clear();
	_badAngleIdx.clear();
	cyPoint3f t0, t1;
	for (int fi = 0; fi < (int)_dualnodes.size(); fi++)
	{
		HE_Face* f = _poly->face(fi);

		if (f->NumHalfEdge() == 3 || f->hole()) 
		{
			_isBadFaces.push_back(false);
			_badAngleIdx.push_back(-1);
			continue;
		}

		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;

		bool determined = false;
		do 
		{
			if ((*e)->twin()->face()->hole())
			{
				_isBadFaces.push_back(false);
				_badAngleIdx.push_back(-1);
				determined = true;
				break;
			}
			e++;
		} while (e != sentinel);

		if (determined) continue;

		bool isBad = false;
		int badIdx = -1;
		e = f->begin();
		int i = 1;
		do
		{
			if ((*e)->dst()->degree() == 3)
			{
				isBad = true;
				badIdx = i - 1;
				break;
			}
			++e; ++i;
		} while (e != sentinel);
#if 0
		if (isBad == false) {
			e = f->begin();
			e--;
			sentinel = e;
			HE_Face::const_edge_circulator nxtE = e;
			nxtE--;
			isBad = false;
			badIdx = -1;
			i = 0;
			float worstAngle = 0;
			float dotValue;
			do
			{
				t0 = (*e)->tangent();		t0.Normalize();
				t1 = (*nxtE)->tangent();	t1.Normalize();
				dotValue = fabsf(t0.Dot(t1));
				if (dotValue > 0.8f)
				{
					if (dotValue > worstAngle)
					{
						worstAngle = dotValue;
						badIdx = i;
					}
					isBad = true;
				}

				++e; ++nxtE; ++i;
			} while (e != sentinel);
		}
#endif
		_isBadFaces.push_back(isBad);
		_badAngleIdx.push_back(badIdx);
	}
}

void DualGraph::findWaleMismatch()
{
	HE_HalfEdge* he1; HE_HalfEdge* he2;
	//int he1_idx, he2_idx;
	for (int i = 0; i < numEdges(); i++)
	{
		he1 = _poly->halfedge(_dualedges[i]._hei[0]);
		he2 = _poly->halfedge(_dualedges[i]._hei[1]);

		if (he1->TestFlag(EFLAG_HORIZONTAL) && he2->TestFlag(EFLAG_HORIZONTAL))
		{
			HE_EDGETYPE he1_type = he1->face()->HalfEdgeType(he1);
			HE_EDGETYPE he2_type = he2->face()->HalfEdgeType(he2);

			if (he1_type == HE_BTM && he2_type == HE_BTM)
			{
				_mismatchTypes[i] = MM_WALE;
			}
			else if (he2_type == HE_TOP && he1_type == HE_TOP)
			{
				_mismatchTypes[i] = MM_WALE;
			}
			//else std::cout << "here\n";
		}
	}
}

void DualGraph::findFaceGroups_v2()
{
	int solutionIdx = 0;
	std::vector<bool> faceMarkers;

	for (int fi = 0; fi < _poly->numFaces(); fi++) faceMarkers.push_back(false);

	for (int fi = 0; fi < _poly->numFaces(); fi++)
	{
		if (faceMarkers[fi]) continue;

		faceMarkers[fi] = true;

		std::vector<int> faceList;
		std::vector<int> tmpFaceList;

		faceList.push_back(fi);
		HE_Face* f = _poly->face(fi);

		// get next face
		int nxtFi = -1;
		int outdegree = f->NumHalfEdge();
		HE_Face::const_edge_circulator e = f->begin();
		for (int i = 0; i < outdegree; i++, e++)
		{
			if ((*e)->TestFlag(EFLAG_VERTICAL) && (*e)->twin()->TestFlag(EFLAG_VERTICAL)) {
				const HE_Face* nf = (*e)->twin()->face();
				if (!faceMarkers[nf->index()]) { nxtFi = nf->index(); break; }
			}
		}

		while (nxtFi != -1)
		{
			int nnxtFi = -1;
			faceList.push_back(nxtFi);
			HE_Face* nf = _poly->face(nxtFi);
			faceMarkers[nxtFi] = true;
			int noutdegree = nf->NumHalfEdge();
			HE_Face::const_edge_circulator e = nf->begin();
			for (int i = 0; i < noutdegree; i++, e++)
			{
				if ((*e)->TestFlag(EFLAG_VERTICAL) && (*e)->twin()->TestFlag(EFLAG_VERTICAL)) {
					const HE_Face* nnf = (*e)->twin()->face();
					if (!faceMarkers[nnf->index()]) { nnxtFi = nnf->index(); break; }
				}
			}
			nxtFi = nnxtFi;
		}

		/* another direction */
		nxtFi = -1;
		e = f->begin();
		for (int i = 0; i < outdegree; i++, e++)
		{
			if ((*e)->TestFlag(EFLAG_VERTICAL) && (*e)->twin()->TestFlag(EFLAG_VERTICAL)) {
				const HE_Face* nf = (*e)->twin()->face();
				if (!faceMarkers[nf->index()]) { nxtFi = nf->index(); break; }
			}
		}

		while (nxtFi != -1)
		{
			int nnxtFi = -1;
			tmpFaceList.push_back(nxtFi);
			HE_Face* nf = _poly->face(nxtFi);
			faceMarkers[nxtFi] = true;
			int noutdegree = nf->NumHalfEdge();
			HE_Face::const_edge_circulator e = nf->begin();
			for (int i = 0; i < noutdegree; i++, e++)
			{
				if ((*e)->TestFlag(EFLAG_VERTICAL) && (*e)->twin()->TestFlag(EFLAG_VERTICAL)) {
					const HE_Face* nnf = (*e)->twin()->face();
					if (!faceMarkers[nnf->index()]) { nnxtFi = nnf->index(); break; }
				}
			}

			nxtFi = nnxtFi;
		}

		/* connect two face lists together*/
		reverse(faceList.begin(), faceList.end());
		for (int i = 0; i < (int)tmpFaceList.size(); i++) faceList.push_back(tmpFaceList[i]);

		_groupFaceIdx.push_back(faceList);
	}

	/* group id */
	for (int fj = 0; fj < (int)_groupFaceIdx.size(); fj++)
	{
		for (int fi = 0; fi < (int)_groupFaceIdx[fj].size(); fi++)
		{
			_dualnodes[_groupFaceIdx[fj][fi]]._groupId = fj;
		}
	}
}

void DualGraph::findFaceGroups()
{
	for (int ni = 0; ni < (int)_dualnodes.size(); ni++)
		_dualnodes[ni]._groupId = -1;

	_groupFaceIdx.clear();

	int solutionIdx = 0;
	std::vector<bool> faceMarkers;

	for (int fi = 0; fi < _poly->numFaces(); fi++) faceMarkers.push_back(false);

	for (int fi = 0; fi < _poly->numFaces(); fi++) 
	{
		if (faceMarkers[fi]) continue;

		faceMarkers[fi] = true;

		std::vector<int> faceList;
		std::vector<int> tmpFaceList;

		faceList.push_back(fi);
		HE_Face* f = _poly->face(fi);

		if (f->hole()) continue;

		// get next face
		int nxtFi = -1;
		int outdegree = f->NumHalfEdge();
		HE_Face::const_edge_circulator e = f->begin();
		for (int i = 0; i < outdegree; i++, e++)
		{
			if ((*e)->TestFlag(EFLAG_VERTICAL) && (*e)->twin()->TestFlag(EFLAG_VERTICAL)) {
				const HE_Face* nf = (*e)->twin()->face();
				if (!faceMarkers[nf->index()]) { nxtFi = nf->index(); break; }			
			}
		}

		while (nxtFi != -1)
		{
			int nnxtFi = -1;
			faceList.push_back(nxtFi);
			HE_Face* nf = _poly->face(nxtFi);
			faceMarkers[nxtFi] = true;
			int noutdegree = nf->NumHalfEdge();
			HE_Face::const_edge_circulator e = nf->begin();
			for (int i = 0; i < noutdegree; i++, e++)
			{
				if ((*e)->TestFlag(EFLAG_VERTICAL) && (*e)->twin()->TestFlag(EFLAG_VERTICAL)) {
					const HE_Face* nnf = (*e)->twin()->face();
					if (!faceMarkers[nnf->index()]) { nnxtFi = nnf->index(); break; }
				}
			}
			nxtFi = nnxtFi;
		}
			
		/* another direction */	
		nxtFi = -1;
		e = f->begin();
		for (int i = 0; i < outdegree; i++, e++)
		{
			if ((*e)->TestFlag(EFLAG_VERTICAL) && (*e)->twin()->TestFlag(EFLAG_VERTICAL)) {
				const HE_Face* nf = (*e)->twin()->face();
				if (!faceMarkers[nf->index()]) { nxtFi = nf->index(); break; }
			}
		}

		while (nxtFi != -1)
		{
			int nnxtFi = -1;
			tmpFaceList.push_back(nxtFi);
			HE_Face* nf = _poly->face(nxtFi);
			faceMarkers[nxtFi] = true;
			int noutdegree = nf->NumHalfEdge();
			HE_Face::const_edge_circulator e = nf->begin();
			for (int i = 0; i < noutdegree; i++, e++)
			{
				if ((*e)->TestFlag(EFLAG_VERTICAL) && (*e)->twin()->TestFlag(EFLAG_VERTICAL)) {
					const HE_Face* nnf = (*e)->twin()->face();
					if (!faceMarkers[nnf->index()]) { nnxtFi = nnf->index(); break; }
				}
			}

			nxtFi = nnxtFi;
		}

		/* connect two face lists together*/
		reverse(faceList.begin(),faceList.end());
		for (int i = 0; i < (int)tmpFaceList.size(); i++) faceList.push_back(tmpFaceList[i]);	

		_groupFaceIdx.push_back(faceList);
	}
	
	/* group id */
	for (int fj = 0; fj < (int)_groupFaceIdx.size(); fj++)
	{
		for (int fi = 0; fi < (int)_groupFaceIdx[fj].size(); fi++)
		{
			_dualnodes[_groupFaceIdx[fj][fi]]._groupId = fj;
		}
	}
}

void DualGraph::loadGurobiResult(std::string pFilename, bool pFlip)
{
	int numVars = 0;

	std::cout << "Face number: " << _poly->numFaces() << std::endl;
	std::cout << "Halfedge number: " << _poly->numHalfEdges() << std::endl;

	for (int fi = 0; fi < _poly->numFaces(); fi++)
	{
		HE_Face* f = _poly->face(fi);

		if (f->hole()) {
			continue;
		}

		_idxOffsets.push_back(numVars);
		numVars += f->NumHalfEdge();
	}

	/* load gurobi result */
	std::string uv_config = pFilename + "_uv_config.txt";
	std::ifstream uvfile;
	uvfile.open(uv_config);

	std::vector<bool> gurobiResult;
	gurobiResult.clear();
	int fidx = 1;
	for (int i = 0; i < numVars; i++)
	{
		int tmp;
		uvfile >> tmp;
		gurobiResult.push_back(bool(tmp == 1));

		if (i == (_idxOffsets[fidx] - 1)) {
			fidx++;
		}
	}

	/* set halfedge flag */
	int fi = 0;
	for (HE_Polyhedron::const_face_iterator it = _poly->fBegin(); it != _poly->fEnd(); ++it, fi++)
	{
		const HE_Face* f = *it;
		if (f->hole())	continue;

		int offset = _idxOffsets[fi];
		int outdegree = f->NumHalfEdge();
		HE_Face::const_edge_circulator e = f->begin();
		for (int i = 0; i < outdegree; i++, e++)
		{
			HE_HalfEdge* he = _poly->halfedge((*e)->index());
			if (pFlip)
			{
				if (gurobiResult[offset + i])
					he->SetFlag(EFLAG_HORIZONTAL);
				else
					he->SetFlag(EFLAG_VERTICAL);
			}
			else 
			{
				if (!gurobiResult[offset + i])
					he->SetFlag(EFLAG_HORIZONTAL);
				else
					he->SetFlag(EFLAG_VERTICAL);
			}
		}
	}

	_gurobiResultPool.push_back(gurobiResult);

	uvfile.close();
}

void DualGraph::gurobiSolverHori()
{
#ifdef USE_GUROBI
	try {
		GRBEnv* env = 0;
		env = new GRBEnv();

		env->set(GRB_IntParam_PoolSolutions, 1024);
		//env->set(GRB_DoubleParam_PoolGap, 1e-12f);

		GRBModel model = GRBModel(*env);

		std::vector<int> hePolyMappingHori;
		std::vector<int> heHoriMappingPoly;

		int numVars = 0;
		for (int ei = 0; ei < _poly->numHalfEdges(); ei++)
		{
			if (_poly->halfedge(ei)->TestFlag(EFLAG_HORIZONTAL))
			{
				hePolyMappingHori.push_back(numVars);
				heHoriMappingPoly.push_back(ei);
				numVars++;
			}
			else
			{
				hePolyMappingHori.push_back(-1);
			}
		}
		std::cout << "Vars number: " << numVars << std::endl;

		// Create variables
		GRBVar* grbVars = (GRBVar*)malloc(sizeof(GRBVar) * numVars);
		for (int i = 0; i < numVars; i++)
		{
			std::string x = 'x' + std::to_string(i);
			grbVars[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, x);
		}

		// Add constraint: 
		std::string c;
		int ci = 0;
		for (int fi = 0; fi < _poly->numFaces(); fi++)
		{
			const HE_Face* f = _poly->face(fi);
			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;
			int horiCount = 0;
			std::vector<int> horiIdx;
			do
			{
				const HE_HalfEdge* he = *e;

				if (he->TestFlag(EFLAG_HORIZONTAL))
				{
					horiIdx.push_back(he->index());
					horiCount++;
				}

				++e;
			} while (e != sentinel);

			if (horiCount == 2)
			{
				c = 'c' + std::to_string(ci) + '0';
				model.addConstr(grbVars[hePolyMappingHori[horiIdx[0]]] + grbVars[hePolyMappingHori[horiIdx[1]]], GRB_EQUAL, 1, c);
				ci++;
			}
		}

		GRBQuadExpr obje;
		GRBVar* tmpVar_0;
		GRBVar* tmpVar_1;
		
		for (int ei = 0; ei < _poly->numHalfEdges(); ei++)
		{
			HE_HalfEdge* he = _poly->halfedge(ei);
			if (he->TestFlag(EFLAG_HORIZONTAL) && he->twin()->TestFlag(EFLAG_HORIZONTAL))
			{
				tmpVar_0 = &grbVars[hePolyMappingHori[he->index()]];
				tmpVar_1 = &grbVars[hePolyMappingHori[he->twin()->index()]];
				obje += (1.0f - (*tmpVar_0 - *tmpVar_1) * (*tmpVar_0 - *tmpVar_1));
			}
		}

		model.setObjective(obje, GRB_MINIMIZE);

		// Optimize model		
		model.optimize();

		if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
		{
			for (int fi = 0; fi < _poly->numFaces(); fi++)
			{
				HE_Face* f = _poly->face(fi);

				if (f->TestFlag(VEF_FLAG_QUAD))
				{
					if (grbVars[hePolyMappingHori[f->edge()->index()]].get(GRB_DoubleAttr_X) > 0.5f)
						f->edge(f->edge()->next()->next());
				}

				//HE_Face::const_edge_circulator e = f->begin();
				//HE_Face::const_edge_circulator sentinel = e;
				//int horiCount = 0;
				//std::vector<int> horiIdx;
				//do
				//{
				//	const HE_HalfEdge* he = *e;

				//	if (he->TestFlag(EFLAG_HORIZONTAL))
				//	{
				//		horiIdx.push_back(he->index());
				//		horiCount++;
				//	}

				//	++e;
				//} while (e != sentinel);

				//if (horiCount == 2)
				//{
				//	std::cout 
				//		<< grbVars[hePolyMappingHori[horiIdx[0]]].get(GRB_DoubleAttr_X) << " "
				//		<< grbVars[hePolyMappingHori[horiIdx[1]]].get(GRB_DoubleAttr_X) << " | ";
				//}
			}

			//for (int i = 0; i < (int)_groupnodes.size(); i++)
			//{
			//	if (_groupnodes[i]._needFlip)
			//	{
			//		for (int fi = 0; fi < (int)_groupFaceIdx[i].size(); fi++)
			//		{
			//			HE_Face* f = _poly->face(_groupFaceIdx[i][fi]);
			//			if (f->TestFlag(VEF_FLAG_QUAD))
			//			{
			//				f->edge(f->edge()->next());
			//				f->edge(f->edge()->next());
			//			}
			//			else if (f->TestFlag(VEF_FLAG_INC))
			//			{
			//				f->ClearFlag(VEF_FLAG_INC);
			//				f->SetFlag(VEF_FLAG_DEC);
			//			}
			//			else if (f->TestFlag(VEF_FLAG_DEC))
			//			{
			//				f->ClearFlag(VEF_FLAG_DEC);
			//				f->SetFlag(VEF_FLAG_INC);
			//			}
			//			else if (f->TestFlag(VEF_FLAG_SHORTROW_L))
			//			{
			//				f->ClearFlag(VEF_FLAG_SHORTROW_L);
			//				f->edge(f->edge()->next());
			//				f->SetFlag(VEF_FLAG_SHORTROW_R);
			//			}
			//			else if (f->TestFlag(VEF_FLAG_SHORTROW_R))
			//			{
			//				f->ClearFlag(VEF_FLAG_SHORTROW_R);
			//				f->edge(f->edge()->prev());
			//				f->SetFlag(VEF_FLAG_SHORTROW_L);
			//			}
			//		}
			//	}
			//}

			//for (int i = 0; i < (int)_groupnodes.size(); i++)
			//{
			//	if (_groupnodes[i]._needFlip)
			//	{
			//		reverse(_groupFaceIdx[i].begin(), _groupFaceIdx[i].end());
			//	}
			//}
		}
		else
		{
			std::cout << "No solution" << std::endl;
		}
	}
	catch (GRBException e) {
		std::cout << "Error code = " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
	}
	catch (...) {
		std::cout << "Exception during optimization" << std::endl;
	}
#endif
}

void DualGraph::gurobiSolver(std::string pFilename)
{
#ifdef USE_GUROBI
	try {
		GRBEnv* env = 0;
		env = new GRBEnv();

		GRBModel model = GRBModel(*env);

		int numVars = 0;

		std::cout << "Face number: " << _poly->numFaces() << std::endl;
		std::cout << "Halfedge number: " << _poly->numHalfEdges() << std::endl;

 		for (int fi = 0; fi < _poly->numFaces(); fi++)
		{
			HE_Face* f = _poly->face(fi);

			if (f->hole())  continue;

			_idxOffsets.push_back(numVars);
			numVars += f->NumHalfEdge();
		}

		// Create variables
		GRBVar* grbVars = (GRBVar*)malloc(sizeof(GRBVar) * numVars);

		for (int i = 0; i < numVars; i++)
		{
			std::string x = 'x' + std::to_string(i);
			grbVars[i] = model.addVar(0.0, 1.0, 0.5, GRB_BINARY, x);
		}

		// Add constraint: 
		// Fi0 = Fi2
		// Fi1 = Fi3
		// Fi0 + Fi1 = 1

		std::string c;
		for (int fi = 0; fi < _poly->numFaces(); fi++)
		{
			HE_Face* f = _poly->face(fi);

			if (f->hole()) continue;

			int offset = _idxOffsets[fi];
			int outdegree = f->NumHalfEdge();
			if (outdegree == 3) {
				c = 'c' + std::to_string(fi) + '0';
				model.addConstr(grbVars[offset + 0] + grbVars[offset + 1] + grbVars[offset + 2], GRB_GREATER_EQUAL, 1, c);

				c = 'c' + std::to_string(fi) + '1';
				model.addConstr(grbVars[offset + 0] + grbVars[offset + 1] + grbVars[offset + 2], GRB_LESS_EQUAL, 2, c);
			}
			else if (outdegree == 4) {
				c = 'c' + std::to_string(fi) + '0';
				model.addConstr(grbVars[offset + 0], GRB_EQUAL, grbVars[offset + 2], c);

				c = 'c' + std::to_string(fi) + '1';
				model.addConstr(grbVars[offset + 1], GRB_EQUAL, grbVars[offset + 3], c);

				c = 'c' + std::to_string(fi) + '2';
				model.addConstr(grbVars[offset + 0] + grbVars[offset + 1], GRB_EQUAL, 1.0f, c);

				c = 'c' + std::to_string(fi) + '3';
				model.addConstr(grbVars[offset + 0] + grbVars[offset + 1] + grbVars[offset + 2] + grbVars[offset + 3], GRB_EQUAL, 2.0f, c);
			}
			else {
				std::cout << "ERROR: shouldn't be here!\n";
			}
		}
	
		GRBQuadExpr obje;
		float weight;
		for (int i = 0; i < numEdges(); i++)
		{
			DualEdge tmpde = _dualedges[i];
			int f0 = tmpde._fi[0], f1 = tmpde._fi[1];
			HE_HalfEdge* he0 = _poly->halfedge(tmpde._hei[0]);
			HE_HalfEdge* he1 = _poly->halfedge(tmpde._hei[1]);
			int idx0 = _poly->face(f0)->HalfEdgeIdx(he0);
			int idx1 = _poly->face(f1)->HalfEdgeIdx(he1);

			if (_poly->face(f0)->hole() || _poly->face(f1)->hole()) continue;

			weight = 1.0f;
			obje += weight * (grbVars[_idxOffsets[f0] + idx0] - grbVars[_idxOffsets[f1] + idx1]) * (grbVars[_idxOffsets[f0] + idx0] - grbVars[_idxOffsets[f1] + idx1]);
		}

		std::cout << "variables number: " << numVars << std::endl;
#if 1
		model.setObjective(obje, GRB_MINIMIZE);

		// Optimize model		
		model.optimize();

		if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
		{
			_nSolutions = model.get(GRB_IntAttr_SolCount);
			std::cout << "Number of solutions found: " << _nSolutions << std::endl;

			_nSolutions = 1;

			for (int e = 0; e < _nSolutions; e++) {

				std::ofstream uvfile;
				std::string uv_config = pFilename + "_uv_config.txt";
				uvfile.open(uv_config);
				
				model.set(GRB_IntParam_SolutionNumber, e);

				std::vector<bool> gurobiResult;
				gurobiResult.clear();
				int fidx = 1;
				for (int i = 0; i < numVars; i++)
				{
					gurobiResult.push_back(grbVars[i].get(GRB_DoubleAttr_X) > 0.5f);

					if (grbVars[i].get(GRB_DoubleAttr_X) > 0.5f)
						uvfile << int(1) << " ";
					else 
						uvfile << int(0) << " ";

					if (i == (_idxOffsets[fidx] - 1)) 
						fidx++;
				}

				uvfile.close();

				_gurobiResultPool.push_back(gurobiResult);
			}
		}
		else
		{
			std::cout << "No solution" << std::endl;
		}

#endif
	}
	catch (GRBException e) {
		std::cout << "Error code = " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
	}
	catch (...) {
		std::cout << "Exception during optimization" << std::endl;
	}
#endif
}

void DualGraph::cut(std::vector<HE_Vertex>& pVerts, std::vector<std::vector<int> >& pFaces)
{
	// assign vertices list
	int nv = _poly->numVertices();
	for (int vi = 0; vi < nv; vi++)
	{
		HE_Vertex v = HE_Vertex(_poly->vertex(vi)->position(), vi);
		pVerts.push_back(v);
	}

	// face idx
	for (int fi = 0; fi < _poly->numFaces(); fi++)
	{
		HE_Face* f = _poly->face(fi);

		if (f->hole()) continue;

		std::vector<int> cornerList;
		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;
		if (!_isBadFaces[fi])
		{
			do
			{
				cornerList.push_back((*e)->dst()->index());
				++e;
			} while (e != sentinel);

			pFaces.push_back(cornerList);
		}
		else
		{
			int idx = _badAngleIdx[fi];
			if (idx == 0 || idx == 2)
			{
				cornerList.push_back((*e)->dst()->index());
				e++;
				cornerList.push_back((*e)->dst()->index());
				e++;
				cornerList.push_back((*e)->dst()->index());
				pFaces.push_back(cornerList);

				if ((cornerList[0] == cornerList[1]) || (cornerList[0] == cornerList[2]) || (cornerList[1] == cornerList[2]))
					std::cout << cornerList[0] << " " << cornerList[1] << " " << cornerList[2] << std::endl;

				cornerList.clear();
				cornerList.push_back((*e)->dst()->index());
				e++;
				cornerList.push_back((*e)->dst()->index());
				e++;
				cornerList.push_back((*e)->dst()->index());
				pFaces.push_back(cornerList);

				if ((cornerList[0] == cornerList[1]) || (cornerList[0] == cornerList[2]) || (cornerList[1] == cornerList[2]))
					std::cout << cornerList[0] << " " << cornerList[1] << " " << cornerList[2] << std::endl;
			}
			else if (idx == 1 || idx == 3)
			{
				e++;
				cornerList.push_back((*e)->dst()->index());
				e++;
				cornerList.push_back((*e)->dst()->index());
				e++;
				cornerList.push_back((*e)->dst()->index());
				pFaces.push_back(cornerList);

				if ((cornerList[0] == cornerList[1]) || (cornerList[0] == cornerList[2]) || (cornerList[1] == cornerList[2]))
					std::cout << cornerList[0] << " " << cornerList[1] << " " << cornerList[2] << std::endl;

				cornerList.clear();
				cornerList.push_back((*e)->dst()->index());
				e++;
				cornerList.push_back((*e)->dst()->index());
				e++;
				cornerList.push_back((*e)->dst()->index());
				pFaces.push_back(cornerList);

				if ((cornerList[0] == cornerList[1]) || (cornerList[0] == cornerList[2]) || (cornerList[1] == cornerList[2]))
					std::cout << cornerList[0] << " " << cornerList[1] << " " << cornerList[2] << std::endl;
			}		
		}
	}
}

void DualGraph::subdivide(std::vector<HE_Vertex>& pVerts, std::vector<std::vector<int> >& pFaces)
{
	// assign vertices list
	int nv = _poly->numVertices();
	for (int vi = 0; vi < nv; vi++)
	{
		HE_Vertex v = HE_Vertex(_poly->vertex(vi)->position(), vi);
		pVerts.push_back(v);
	}

	std::vector<int> edgeIdxMapping;
	for (int i = 0; i < _poly->numHalfEdges(); i++) edgeIdxMapping.push_back(-1);

	cyPoint3f edgeCenter;
	int hei0, hei1;
	HE_Vertex* v0; HE_Vertex* v1;
	int ne = numEdges();
	for (int ei = 0; ei < ne; ei++ )
	{
		hei0 = _dualedges[ei]._hei[0];
		hei1 = _dualedges[ei]._hei[1];
		edgeIdxMapping[hei0] = ei;
		edgeIdxMapping[hei1] = ei;
		v0 = _poly->halfedge(hei0)->dst();
		v1 = _poly->halfedge(hei1)->dst();
		edgeCenter = (v0->position() + v1->position()) * 0.5f;
		HE_Vertex vv = HE_Vertex(edgeCenter, nv + ei);
		pVerts.push_back(vv);
	}

	cyPoint3f faceCenter;
	for (int fi = 0; fi < numNodes(); fi++ )
	{
		faceCenter = _dualnodes[fi]._fc;
		HE_Vertex vv = HE_Vertex(faceCenter, nv + ne + fi);
		pVerts.push_back(vv);
	}

	// assign faces list
	HE_Face* f;
	for (int fi = 0; fi < numNodes(); fi++)
	{
		f = _poly->face(fi);
		if (f->hole()) continue;

		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;

		int i = 0;
		do
		{
			int v0 = nv + edgeIdxMapping[(*e)->index()];
			int v1 = (*e)->dst()->index();
			int v2 = nv + edgeIdxMapping[(*e)->next()->index()];
			int v3 = nv + ne + fi;

			std::vector<int> viList;
			if (i == 0)
			{
				viList.push_back(v1); viList.push_back(v2); viList.push_back(v3); viList.push_back(v0);
			}
			else if (i == 1)
			{
				viList.push_back(v0); viList.push_back(v1); viList.push_back(v2); viList.push_back(v3);
			}
			else if (i == 2)
			{	
				viList.push_back(v3); viList.push_back(v0); viList.push_back(v1); viList.push_back(v2);
			}
			else if (i == 3)
			{	
				viList.push_back(v2); viList.push_back(v3); viList.push_back(v0); viList.push_back(v1);
			}
 
			pFaces.push_back(viList);

			++e; ++i;
		} while (e != sentinel);
	}
}

void DualGraph::fitFaceEdgesOrder()
{
	/* rotate face edge to make sure the first halfedge is horizontal */
	for (int fj = 0; fj < (int)_groupFaceIdx.size(); fj++)
	{
		for (int fi = 0; fi < (int)_groupFaceIdx[fj].size(); fi++)
		{
			HE_Face* f = _poly->face(_groupFaceIdx[fj][fi]);
			if (f->hole()) continue;

			while (f->edge()->TestFlag(EFLAG_VERTICAL)) f->edge(f->edge()->next());

			if (f->NumHalfEdge() == 4) f->SetFlag(VEF_FLAG_QUAD);
		}
	}

	/*
	 *		   H(0) Top
	 *		 -------
	 *	V(1) |  /\ | V(3)
	 *       |  |  |     
	 *		 -------
	 *		   H(2) Bottom
	 */

	/* flip each quad faces */
	for (int fj = 0; fj < (int)_groupFaceIdx.size(); fj++)
	{
		for (int fi = 0; fi < (int)_groupFaceIdx[fj].size() - 1; fi++)
		{
			HE_Face* f0 = _poly->face(_groupFaceIdx[fj][fi]);
			HE_Face* f1 = _poly->face(_groupFaceIdx[fj][fi + 1]);
			if (f0->hole() || f1->hole()) std::cout << "error!\n";

			int tmpFi0 = _groupFaceIdx[fj][fi];
			int tmpFi1 = _groupFaceIdx[fj][fi + 1];

			/* find halfedge of f0, the twin of which belongs to f1 */
			HE_Face::const_edge_circulator e = f0->begin();
			while ((*e)->twin()->face()->index() != tmpFi1 || (*e)->TestFlag(EFLAG_HORIZONTAL))  e++;
			
			/* rotate f0 */
			if (fi == 0) {
				if (f0->NumHalfEdge() == 4) {
					while (f0->HalfEdgeIdx((*e)) != 3)
					{
						f0->edge(f0->edge()->next());
					}
				}
				else {
					while (f0->HalfEdgeIdx((*e)) != 2)
					{
						f0->edge(f0->edge()->next());
					}
				}
			}

			/* rotate f1 */
			const HE_HalfEdge* tmpHe = (*e)->twin();
			while (f1->HalfEdgeIdx(tmpHe) != 1)
			{
				f1->edge(f1->edge()->next());
			}
		}
	}
}

void DualGraph::labelQuadFace()
{
	for (int fi = 0; fi < _poly->numFaces(); fi++)
	{
		HE_Face* f = _poly->face(fi);
		if (f->hole()) continue;

		while (f->edge()->TestFlag(EFLAG_VERTICAL)) f->edge(f->edge()->next());

		if (f->NumHalfEdge() == 4) f->SetFlag(VEF_FLAG_QUAD);
	}
}

void DualGraph::labelTriangleFace()
{
	for (int fi = 0; fi < _poly->numFaces(); fi++)
	{
		if (!_poly->face(fi)->TestFlag(VEF_FLAG_QUAD))
			_poly->face(fi)->ClearFlags();
	}

	for (int fj = 0; fj < (int)_groupFaceIdx.size(); fj++)
	{
		for (int fi = 0; fi < (int)_groupFaceIdx[fj].size(); fi++)
		{
			HE_Face* f = _poly->face(_groupFaceIdx[fj][fi]);
			if (!f->TestFlag(VEF_FLAG_QUAD)) {
				int numHoriEdge = 0;
				int numVertEdge = 0;
				
				HE_Face::const_edge_circulator e = f->begin();
				HE_Face::const_edge_circulator sentinel = e;
				bool isL = false;
				int i = 0;
				do
				{
					if ((*e)->TestFlag(EFLAG_HORIZONTAL))	numHoriEdge++;
					if ((*e)->TestFlag(EFLAG_VERTICAL))		numVertEdge++;

					if (i == 2 && (*e)->TestFlag(EFLAG_VERTICAL)) isL = true;

					++e; ++i;
				} while (e != sentinel);

				if (numHoriEdge == 2 && numVertEdge == 1) {
					if (isL)	f->SetFlag(VEF_FLAG_SHORTROW_L);
					else		f->SetFlag(VEF_FLAG_SHORTROW_R);
				}

				if (numHoriEdge == 1 && numVertEdge == 2) {

					while (f->edge()->TestFlag(EFLAG_VERTICAL)) f->edge(f->edge()->next());

					HE_HalfEdge* he = f->edge()->next();
					HE_Face* nf = he->twin()->face();

					if (nf->NumHalfEdge() == 4)
					{
						if (nf->HalfEdgeIdx(he->twin()) == 3)
							f->SetFlag(VEF_FLAG_INC);
						else if (nf->HalfEdgeIdx(he->twin()) == 1)
							f->SetFlag(VEF_FLAG_DEC);
					}
					else if (nf->NumHalfEdge() == 3)
					{
						if (fi == 0)
						{
							HE_Face* prev_f = _poly->face(_groupFaceIdx[fj][fi + 1]);

							HE_Face::const_edge_circulator e = f->begin();
							HE_Face::const_edge_circulator sentinel = e;
							while (true)
							{
								if ((*e)->twin()->face()->index() == prev_f->index())
								{
									break;
								}
								e++;
							}

							if (f->HalfEdgeIdx((*e)) == 2)
								f->SetFlag(VEF_FLAG_INC);
							else if (f->HalfEdgeIdx((*e)) == 1)
								f->SetFlag(VEF_FLAG_DEC);
						}
						else {
							HE_Face* prev_f = _poly->face(_groupFaceIdx[fj][fi - 1]);

							HE_Face::const_edge_circulator e = f->begin();
							HE_Face::const_edge_circulator sentinel = e;
							while (true)
							{
								if ((*e)->twin()->face()->index() == prev_f->index())
								{
									break;
								}
								e++;
							}

							if (f->HalfEdgeIdx((*e)) == 1)
								f->SetFlag(VEF_FLAG_INC);
							else if (f->HalfEdgeIdx((*e)) == 2)
								f->SetFlag(VEF_FLAG_DEC);
							
						}
					}
				}
			}
		}
	}
}

void DualGraph::breakConnectIncOrDec()
{
	for (int gi = 0; gi < (int)_groupFaceIdx.size(); gi++)
	{
		for (int i = 0; i < (int)_groupFaceIdx[gi].size(); i++) 
		{
	
			int f1i = _groupFaceIdx[gi][i];
			int f2i = (i == (_groupFaceIdx[gi].size() - 1)) ? _groupFaceIdx[gi][0] : _groupFaceIdx[gi][i+1];

			HE_Face* f1 = _poly->face(f1i);
			HE_Face* f2 = _poly->face(f2i);

			if (f1->hole() || f2->hole()) continue;

			if (f1->TestFlag(VEF_FLAG_INC) && f2->TestFlag(VEF_FLAG_INC))// && !f0->TestFlag(VEF_FLAG_DEC))
			{
				HE_Face::const_edge_circulator e = f1->begin();
				HE_Face::const_edge_circulator sentinel = e;
				while (true)
				{
					if ((*e)->twin()->face()->index() == f2i)
					{
						break;
					}
					e++;
				}

				HE_HalfEdge* he = _poly->halfedge((*e)->index());

				if ((he->next()->TestFlag(EFLAG_HORIZONTAL) &&
					he->prev()->TestFlag(EFLAG_HORIZONTAL) ) || (
					he->twin()->next()->TestFlag(EFLAG_HORIZONTAL) &&
					he->twin()->prev()->TestFlag(EFLAG_HORIZONTAL)))
				{
					;
				}
				else
				{
					he->ClearFlags();
					he->twin()->ClearFlags();
					he->SetFlag(EFLAG_HORIZONTAL);
					he->twin()->SetFlag(EFLAG_HORIZONTAL);
				}
			}

			else if (f1->TestFlag(VEF_FLAG_DEC) && f2->TestFlag(VEF_FLAG_DEC))// && !f0->TestFlag(VEF_FLAG_INC))
			{
				HE_Face::const_edge_circulator e = f1->begin();
				HE_Face::const_edge_circulator sentinel = e;
				while (true)
				{
					if ((*e)->twin()->face()->index() == f2i)
					{
						break;
					}
					e++;
				} 

				HE_HalfEdge* he = _poly->halfedge((*e)->index());

				if ((he->next()->TestFlag(EFLAG_HORIZONTAL) && he->prev()->TestFlag(EFLAG_HORIZONTAL) )||
					(he->twin()->next()->TestFlag(EFLAG_HORIZONTAL) && he->twin()->prev()->TestFlag(EFLAG_HORIZONTAL)))
				{
					;
				}
				else
				{
					he->ClearFlags();
					he->twin()->ClearFlags();
					he->SetFlag(EFLAG_HORIZONTAL);
					he->twin()->SetFlag(EFLAG_HORIZONTAL);
				}
			}
		}
	}
}

void DualGraph::findTopBtmGroupId()
{
	for (int gi = 0; gi < (int)_groupFaceIdx.size(); gi++) 
	{
		GroupNode gn;
		_groupnodes.push_back(gn);
	}

	/* find top groups id */
	for (int gi = 0; gi < (int)_groupFaceIdx.size(); gi++)
	{
		std::vector<int> topGroupIdsCount, btmGroupIdsCount;
		std::vector<EDGE_TYPE> topGroupEdgeTypes, btmGroupEdgeTypes;
		for (int i = 0; i < (int)_groupFaceIdx.size(); i++) {
			topGroupIdsCount.push_back(0); btmGroupIdsCount.push_back(0);
			topGroupEdgeTypes.push_back(EG_NONE);
			btmGroupEdgeTypes.push_back(EG_NONE);
		}

		for (int fi = 0; fi < (int)_groupFaceIdx[gi].size(); fi++)
		{
			HE_Face* f = _poly->face(_groupFaceIdx[gi][fi]);

			/* top halfedge */
			if (f->TestFlag(VEF_FLAG_QUAD) && f->edge()->twin()->face()->TestFlag(VEF_FLAG_QUAD))
			{
				if (f->edge()->TestFlag(EFLAG_HORIZONTAL) && f->edge()->twin()->TestFlag(EFLAG_HORIZONTAL))
				{
					topGroupIdsCount[_dualnodes[f->edge()->twin()->face()->index()]._groupId]++;
					if (f->edge()->twin()->face()->HalfEdgeIdx(f->edge()->twin()) == 0)
						topGroupEdgeTypes[_dualnodes[f->edge()->twin()->face()->index()]._groupId] = EG_TOP;
					else if (f->edge()->twin()->face()->HalfEdgeIdx(f->edge()->twin()) == 2)
						topGroupEdgeTypes[_dualnodes[f->edge()->twin()->face()->index()]._groupId] = EG_BTM;
					else std::cout << "Error! - findTopGroup\n";
				}
			}
			else if (f->TestFlag(VEF_FLAG_QUAD) && f->edge()->twin()->face()->TestFlag(VEF_FLAG_DEC))
			{
				topGroupIdsCount[_dualnodes[f->edge()->twin()->face()->index()]._groupId]++;
				topGroupEdgeTypes[_dualnodes[f->edge()->twin()->face()->index()]._groupId] = EG_BTM;
			}
			else if (f->TestFlag(VEF_FLAG_QUAD) && f->edge()->twin()->face()->TestFlag(VEF_FLAG_INC))
			{
				topGroupIdsCount[_dualnodes[f->edge()->twin()->face()->index()]._groupId]++;
				topGroupEdgeTypes[_dualnodes[f->edge()->twin()->face()->index()]._groupId] = EG_TOP;
			}
			else if (f->TestFlag(VEF_FLAG_INC) && f->edge()->twin()->face()->TestFlag(VEF_FLAG_QUAD))
			{
				topGroupIdsCount[_dualnodes[f->edge()->twin()->face()->index()]._groupId]++;
				if (f->edge()->twin()->face()->HalfEdgeIdx(f->edge()->twin()) == 0)
					topGroupEdgeTypes[_dualnodes[f->edge()->twin()->face()->index()]._groupId] = EG_TOP;
				else if (f->edge()->twin()->face()->HalfEdgeIdx(f->edge()->twin()) == 2)
					topGroupEdgeTypes[_dualnodes[f->edge()->twin()->face()->index()]._groupId] = EG_BTM;
				else std::cout << "Error! - findTopGroup\n";
			} 
			else if (f->TestFlag(VEF_FLAG_QUAD) && f->edge()->twin()->face()->TestFlag(VEF_FLAG_SHORTROW_L))
			{
				if (f->edge()->TestFlag(EFLAG_HORIZONTAL) && f->edge()->twin()->TestFlag(EFLAG_HORIZONTAL))
				{
					topGroupIdsCount[_dualnodes[f->edge()->twin()->face()->index()]._groupId]++;
					if (f->edge()->twin()->face()->HalfEdgeIdx(f->edge()->twin()) == 0)
						topGroupEdgeTypes[_dualnodes[f->edge()->twin()->face()->index()]._groupId] = EG_TOP;
					else if (f->edge()->twin()->face()->HalfEdgeIdx(f->edge()->twin()) == 1)
						topGroupEdgeTypes[_dualnodes[f->edge()->twin()->face()->index()]._groupId] = EG_BTM;
					else std::cout << "Error! - findTopGroup\n";
				}
			}
			else if (f->TestFlag(VEF_FLAG_QUAD) && f->edge()->twin()->face()->TestFlag(VEF_FLAG_SHORTROW_R))
			{
				if (f->edge()->TestFlag(EFLAG_HORIZONTAL) && f->edge()->twin()->TestFlag(EFLAG_HORIZONTAL))
				{
					topGroupIdsCount[_dualnodes[f->edge()->twin()->face()->index()]._groupId]++;
					if (f->edge()->twin()->face()->HalfEdgeIdx(f->edge()->twin()) == 0)
						topGroupEdgeTypes[_dualnodes[f->edge()->twin()->face()->index()]._groupId] = EG_TOP;
					else if (f->edge()->twin()->face()->HalfEdgeIdx(f->edge()->twin()) == 2)
						topGroupEdgeTypes[_dualnodes[f->edge()->twin()->face()->index()]._groupId] = EG_BTM;
					else std::cout << "Error! - findTopGroup\n";
				}
			}
			else if ((f->TestFlag(VEF_FLAG_SHORTROW_R) || f->TestFlag(VEF_FLAG_SHORTROW_L)) && f->edge()->twin()->face()->TestFlag(VEF_FLAG_QUAD))
			{
				if (f->edge()->TestFlag(EFLAG_HORIZONTAL) && f->edge()->twin()->TestFlag(EFLAG_HORIZONTAL))
				{
					topGroupIdsCount[_dualnodes[f->edge()->twin()->face()->index()]._groupId]++;
					if (f->edge()->twin()->face()->HalfEdgeIdx(f->edge()->twin()) == 0)
						topGroupEdgeTypes[_dualnodes[f->edge()->twin()->face()->index()]._groupId] = EG_TOP;
					else if (f->edge()->twin()->face()->HalfEdgeIdx(f->edge()->twin()) == 2)
						topGroupEdgeTypes[_dualnodes[f->edge()->twin()->face()->index()]._groupId] = EG_BTM;
					else std::cout << "Error! - findTopGroup\n";
				}
			}

			/* bottom halfedge */
			HE_HalfEdge* he = f->edge()->next()->next();
			if (f->TestFlag(VEF_FLAG_QUAD) && he->twin()->face()->TestFlag(VEF_FLAG_QUAD))
			{
				if (he->TestFlag(EFLAG_HORIZONTAL) && he->twin()->TestFlag(EFLAG_HORIZONTAL))
				{
					btmGroupIdsCount[_dualnodes[he->twin()->face()->index()]._groupId]++;

					if (he->twin()->face()->HalfEdgeIdx(he->twin()) == 0)
						btmGroupEdgeTypes[_dualnodes[he->twin()->face()->index()]._groupId] = EG_TOP;
					else if (he->twin()->face()->HalfEdgeIdx(he->twin()) == 2)
						btmGroupEdgeTypes[_dualnodes[he->twin()->face()->index()]._groupId] = EG_BTM;
					else std::cout << "Error! - findTopGroup\n";
				}
			}
			else if (f->TestFlag(VEF_FLAG_QUAD) && he->twin()->face()->TestFlag(VEF_FLAG_DEC))
			{
				btmGroupIdsCount[_dualnodes[f->edge()->twin()->face()->index()]._groupId]++;
				btmGroupEdgeTypes[_dualnodes[f->edge()->twin()->face()->index()]._groupId] = EG_BTM;
			}
			else if (f->TestFlag(VEF_FLAG_QUAD) && he->twin()->face()->TestFlag(VEF_FLAG_INC))
			{
				btmGroupIdsCount[_dualnodes[f->edge()->twin()->face()->index()]._groupId]++;
				btmGroupEdgeTypes[_dualnodes[f->edge()->twin()->face()->index()]._groupId] = EG_TOP;
			}
			else if (f->TestFlag(VEF_FLAG_QUAD) && he->twin()->face()->TestFlag(VEF_FLAG_SHORTROW_L))
			{
				if (he->TestFlag(EFLAG_HORIZONTAL) && he->twin()->TestFlag(EFLAG_HORIZONTAL))
				{
					btmGroupIdsCount[_dualnodes[he->twin()->face()->index()]._groupId]++;

					if (he->twin()->face()->HalfEdgeIdx(he->twin()) == 0)
						btmGroupEdgeTypes[_dualnodes[he->twin()->face()->index()]._groupId] = EG_TOP;
					else if (he->twin()->face()->HalfEdgeIdx(he->twin()) == 1)
						btmGroupEdgeTypes[_dualnodes[he->twin()->face()->index()]._groupId] = EG_BTM;
					else std::cout << "Error! - findTopGroup\n";
				}
			}
			else if (f->TestFlag(VEF_FLAG_QUAD) && he->twin()->face()->TestFlag(VEF_FLAG_SHORTROW_R))
			{
				if (he->TestFlag(EFLAG_HORIZONTAL) && he->twin()->TestFlag(EFLAG_HORIZONTAL))
				{
					btmGroupIdsCount[_dualnodes[he->twin()->face()->index()]._groupId]++;

					if (he->twin()->face()->HalfEdgeIdx(he->twin()) == 0)
						btmGroupEdgeTypes[_dualnodes[he->twin()->face()->index()]._groupId] = EG_TOP;
					else if (he->twin()->face()->HalfEdgeIdx(he->twin()) == 2)
						btmGroupEdgeTypes[_dualnodes[he->twin()->face()->index()]._groupId] = EG_BTM;
					else std::cout << "Error! - findTopGroup\n";
				}
			}
			else if (f->TestFlag(VEF_FLAG_DEC) && f->edge()->twin()->face()->TestFlag(VEF_FLAG_QUAD))
			{
				he = f->edge();
				btmGroupIdsCount[_dualnodes[he->twin()->face()->index()]._groupId]++;

				if (he->twin()->face()->HalfEdgeIdx(he->twin()) == 0)
					btmGroupEdgeTypes[_dualnodes[he->twin()->face()->index()]._groupId] = EG_TOP;
				else if (he->twin()->face()->HalfEdgeIdx(he->twin()) == 2)
					btmGroupEdgeTypes[_dualnodes[he->twin()->face()->index()]._groupId] = EG_BTM;
				else std::cout << "Error! - findTopGroup\n";			
			}
		}

		for (int i = 0; i < (int)topGroupIdsCount.size(); i++)
		{
			if (topGroupIdsCount[i] != 0 && i != gi)
			{
				_groupnodes[gi]._topGroupIds.push_back(i);
				_groupnodes[gi]._topGroupCounts.push_back(topGroupIdsCount[i]);
				_groupnodes[gi]._topEdgeType.push_back(topGroupEdgeTypes[i]);
			}
		}

		for (int i = 0; i < (int)btmGroupIdsCount.size(); i++)
		{
			if (btmGroupIdsCount[i] != 0 && i != gi)
			{
				_groupnodes[gi]._btmGroupIds.push_back(i);
				_groupnodes[gi]._btmGroupCounts.push_back(btmGroupIdsCount[i]);
				_groupnodes[gi]._btmEdgeType.push_back(btmGroupEdgeTypes[i]);
			}
		}
	}

	createGroupEdges();
}

void DualGraph::createGroupEdges()
{
	for (int gi = 0; gi < (int)_groupFaceIdx.size(); gi++)
	{
		for (int i = 0; i < (int)_groupnodes[gi]._topGroupIds.size(); i++)
		{
			int topId = _groupnodes[gi]._topGroupIds[i];

			GroupEdge ge;
			ge._fi[0] = gi;
			ge._fi[1] = topId;
			ge._eType[0] = EG_TOP;
			ge._weight = _groupnodes[gi]._topGroupCounts[i];

			if (_groupnodes[gi]._topEdgeType[i] == EG_TOP)
			{
				for (int j = 0; j < (int)_groupnodes[topId]._topGroupIds.size(); j++)
				{
					if (_groupnodes[topId]._topGroupIds[j] == gi) ge._eType[1] = EG_TOP;
				}
			} 
			else if (_groupnodes[gi]._topEdgeType[i] == EG_BTM)
			{
				for (int j = 0; j < (int)_groupnodes[topId]._btmGroupIds.size(); j++)
				{
					if (_groupnodes[topId]._btmGroupIds[j] == gi) ge._eType[1] = EG_BTM;
				}
			}

			_groupedges.push_back(ge);
		}

		for (int i = 0; i < (int)_groupnodes[gi]._btmGroupIds.size(); i++)
		{
			int btmId = _groupnodes[gi]._btmGroupIds[i];

			GroupEdge ge;
			ge._fi[0] = gi;
			ge._fi[1] = btmId;
			ge._eType[0] = EG_BTM;
			ge._weight = _groupnodes[gi]._btmGroupCounts[i];

			if (_groupnodes[gi]._btmEdgeType[i] == EG_TOP)
			{
				for (int j = 0; j < (int)_groupnodes[btmId]._topGroupIds.size(); j++)
				{
					if (_groupnodes[btmId]._topGroupIds[j] == gi) ge._eType[1] = EG_TOP;					
				}
			}

			if (_groupnodes[gi]._btmEdgeType[i] == EG_BTM)
			{
				for (int j = 0; j < (int)_groupnodes[btmId]._btmGroupIds.size(); j++)
				{
					if (_groupnodes[btmId]._btmGroupIds[j] == gi) ge._eType[1] = EG_BTM;
				}
			}
			_groupedges.push_back(ge);
		}
	}
}

void DualGraph::waleMismatchSolver()
{
#ifdef USE_GUROBI
	try {
		GRBEnv* env = 0;
		env = new GRBEnv();

		env->set(GRB_IntParam_PoolSolutions, 1024);
		//env->set(GRB_DoubleParam_PoolGap, 1e-12f);

		GRBModel model = GRBModel(*env);

		int numVars = 2 * (int)_groupnodes.size();
		std::cout << "Vars number: " << numVars << std::endl;
		// Create variables
		GRBVar* grbVars = (GRBVar*)malloc(sizeof(GRBVar) * numVars);
		for (int i = 0; i < numVars; i++)
		{
			std::string x = 'x' + std::to_string(i);
			grbVars[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, x);
		}

		// Add constraint: 
		std::string c;
		for (int gi = 0; gi < (int)_groupnodes.size(); gi++)
		{
			c = 'c' + std::to_string(gi) + '0';
			model.addConstr(grbVars[gi * 2 + 0] + grbVars[gi * 2 + 1], GRB_EQUAL, 1, c);
		}

		GRBQuadExpr obje;
		float weight;
		GRBVar* tmpVar_0;
		GRBVar* tmpVar_1;
		for (int i = 0; i < (int)_groupedges.size(); i++)
		{
			if (_groupedges[i]._eType[0] == EG_TOP) tmpVar_0 = &grbVars[_groupedges[i]._fi[0] * 2];
			else  tmpVar_0 = &grbVars[_groupedges[i]._fi[0] * 2 + 1];

			if (_groupedges[i]._eType[1] == EG_TOP) tmpVar_1 = &grbVars[_groupedges[i]._fi[1] * 2];
			else  tmpVar_1 = &grbVars[_groupedges[i]._fi[1] * 2 + 1];

			weight = (float)_groupedges[i]._weight;
			//obje += weight * (grbVars[_idxOffsets[f0] + idx0] - grbVars[_idxOffsets[f1] + idx1]) * (grbVars[_idxOffsets[f0] + idx0] - grbVars[_idxOffsets[f1] + idx1]);

			obje += weight * (1.0f - (*tmpVar_0 - *tmpVar_1) * (*tmpVar_0 - *tmpVar_1));
		}

		model.setObjective(obje, GRB_MINIMIZE);

		// Optimize model		
		model.optimize();

		//model.set(GRB SOLUTION_LIMIT, 1);

		if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
		{
			for (int i = 0; i < numVars; i++)
			{
				std::cout << (grbVars[i].get(GRB_DoubleAttr_X) > 0.5f) << " ";

				if (i % 2 == 0) {
					if ((grbVars[i].get(GRB_DoubleAttr_X) < 0.5f))
						_groupnodes[i / 2]._needFlip = true;
					else 
						_groupnodes[i / 2]._needFlip = false;

	/*				if (i == 22) {
						if ((grbVars[i].get(GRB_DoubleAttr_X) > 0.5f))
							_groupnodes[i / 2]._needFlip = false;
						else
							_groupnodes[i / 2]._needFlip = true;
					}*/
				}
				if (i%2 == 1) 
					std::cout << "| ";
			}
			std::cout << std::endl;

			for (int i = 0; i < (int)_groupnodes.size(); i++)
			{
				if (_groupnodes[i]._needFlip) 
				{
					for (int fi = 0; fi < (int)_groupFaceIdx[i].size(); fi++)
					{
						HE_Face* f = _poly->face(_groupFaceIdx[i][fi]);
						if (f->TestFlag(VEF_FLAG_QUAD))
						{
							f->edge(f->edge()->next());
							f->edge(f->edge()->next());
						}
						else if (f->TestFlag(VEF_FLAG_INC)) 
						{ 
							f->ClearFlag(VEF_FLAG_INC);
							f->SetFlag(VEF_FLAG_DEC);
						}
						else if (f->TestFlag(VEF_FLAG_DEC)) 
						{ 
							f->ClearFlag(VEF_FLAG_DEC);
							f->SetFlag(VEF_FLAG_INC); 
						}
						else if (f->TestFlag(VEF_FLAG_SHORTROW_L)) 
						{ 
							f->ClearFlag(VEF_FLAG_SHORTROW_L);
							f->edge(f->edge()->next());
							f->SetFlag(VEF_FLAG_SHORTROW_R); 
						}
						else if (f->TestFlag(VEF_FLAG_SHORTROW_R)) 
						{ 
							f->ClearFlag(VEF_FLAG_SHORTROW_R);
							f->edge(f->edge()->prev());
							f->SetFlag(VEF_FLAG_SHORTROW_L); 
						}
					}
				}
			}
		
			for (int i = 0; i < (int)_groupnodes.size(); i++)
			{
				if (_groupnodes[i]._needFlip)
				{
					reverse(_groupFaceIdx[i].begin(), _groupFaceIdx[i].end());
				}
			}
		}
		else
		{
			std::cout << "No solution" << std::endl;
		}
	}
	catch (GRBException e) {
		std::cout << "Error code = " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
	}
	catch (...) {
		std::cout << "Exception during optimization" << std::endl;
	}
#endif
}

void DualGraph::setSubGroupFace(DualGraph* pDual)
{
	/* compute face idx offset */
	std::vector<int> faceIdxOffset;
	int accuFaceCount = 0;
	for (int i = 0; i < pDual->_poly->numFaces(); i++)
	{
		faceIdxOffset.push_back(accuFaceCount);
		accuFaceCount += pDual->_poly->face(i)->NumHalfEdge();
	}

	/* set subdivision face group idx */
	/* set face flag */
	std::vector<std::vector<int>>* coarseGroupFaceIdx = &pDual->_groupFaceIdx;
	for (int i = 0; i < (int)coarseGroupFaceIdx->size(); i++)
	{
		std::vector<int> faceList;
		for (int j = 0; j < (int)(*coarseGroupFaceIdx)[i].size(); j++)
		{
			HE_Face* f = pDual->_poly->face((*coarseGroupFaceIdx)[i][j]);
			int numFace = f->NumHalfEdge();

			for (int fi = 0; fi < numFace; fi++)
			{
				int coarseFi = faceIdxOffset[(*coarseGroupFaceIdx)[i][j]] + fi;
				faceList.push_back(coarseFi);
				if (f->TestFlag(VEF_FLAG_QUAD))
				{
					_poly->face(coarseFi)->SetFlag(VEF_FLAG_QUAD);
				}
				else if (f->TestFlag(VEF_FLAG_INC))
				{
					if (fi == 0 || fi == 2) _poly->face(coarseFi)->SetFlag(VEF_FLAG_QUAD);
					else _poly->face(coarseFi)->SetFlag(VEF_FLAG_QUAD_INC);

					if (fi == 2) _poly->face(coarseFi)->edge(_poly->face(coarseFi)->edge()->prev());
					else if (fi == 1) _poly->face(coarseFi)->edge(_poly->face(coarseFi)->edge()->prev());
				}
				else if (f->TestFlag(VEF_FLAG_DEC))
				{
					if (fi == 0 || fi == 2) _poly->face(coarseFi)->SetFlag(VEF_FLAG_QUAD);
					else _poly->face(coarseFi)->SetFlag(VEF_FLAG_QUAD_DEC);

					if (fi == 0) _poly->face(coarseFi)->edge(_poly->face(coarseFi)->edge()->next()->next());
					else if (fi == 2) _poly->face(coarseFi)->edge(_poly->face(coarseFi)->edge()->next());
					else if (fi == 1) _poly->face(coarseFi)->edge(_poly->face(coarseFi)->edge()->prev());
				}
				else if (f->TestFlag(VEF_FLAG_SHORTROW_L))
				{
					if (fi == 1 || fi == 2)
					{
						_poly->face(coarseFi)->SetFlag(VEF_FLAG_QUAD);
						_poly->face(coarseFi)->edge(_poly->face(coarseFi)->edge()->prev());
					}
					else _poly->face(coarseFi)->SetFlag(VEF_FLAG_QUAD_SHORTROW_L);
				}
				else if (f->TestFlag(VEF_FLAG_SHORTROW_R))
				{
					if (fi == 0 || fi == 1)
					{
						_poly->face(coarseFi)->SetFlag(VEF_FLAG_QUAD);
					}
					else
					{
						_poly->face(coarseFi)->SetFlag(VEF_FLAG_QUAD_SHORTROW_R);
						_poly->face(coarseFi)->edge(_poly->face(coarseFi)->edge()->next()->next());
					}
				}
			}
		}
		_groupFaceIdx.push_back(faceList);
	}

	//for (int i = 0; i < _groupFaceIdx.size(); i++)
	//{
	//	for (int j = 0; j < _groupFaceIdx[i].size(); j++)
	//	{
	//		std::cout << _groupFaceIdx[i][j] << " ";
	//	}
	//	std::cout << std::endl;
	//}
}

void DualGraph::setGroupEdgeFlag()
{
	for (int fi = 0; fi < numNodes(); fi++)
	{
		HE_Face* f = _poly->face(fi);
		if (f->TestFlag(VEF_FLAG_QUAD))
		{
			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;
			int i = 0;
			do
			{
				HE_HalfEdge* he = _poly->halfedge((*e)->index());

				if (i == 0 || i == 2) he->SetFlag(EFLAG_HORIZONTAL);
				else if (i == 1 || i == 3) he->SetFlag(EFLAG_VERTICAL);

				++e; ++i;
			} while (e != sentinel);
		}
		else if (	f->TestFlag(VEF_FLAG_QUAD_INC)	|| f->TestFlag(VEF_FLAG_QUAD_SHORTROW_L) ||
					f->TestFlag(VEF_FLAG_QUAD_DEC)	|| f->TestFlag(VEF_FLAG_QUAD_SHORTROW_R) )
		{
			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;
			int i = 0;
			do
			{
				HE_HalfEdge* he = _poly->halfedge((*e)->index());

				if (i == 0 || i == 1) he->SetFlag(EFLAG_HORIZONTAL);
				else if (i == 2 || i == 3) he->SetFlag(EFLAG_VERTICAL);

				++e; ++i;
			} while (e != sentinel);
		}
	}
}

void DualGraph::reorderSubdivisionGroupFaceIdx(DualGraph* pDual)
{
	/* compute face idx offset */
	std::vector<int> faceIdxOffset;
	int accuFaceCount = 0;
	for (int i = 0; i < pDual->_poly->numFaces(); i++)
	{
		faceIdxOffset.push_back(accuFaceCount);
		accuFaceCount += pDual->_poly->face(i)->NumHalfEdge();
	}

	/* set subdivision face group idx */
	std::vector<std::vector<int>>* coarseGroupFaceIdx = &pDual->_groupFaceIdx;
	_groupFaceIdx.clear();
	for (int i = 0; i < (int)coarseGroupFaceIdx->size(); i++)
	{
		std::vector<int> faceList;
		for (int j = 0; j < (int)(*coarseGroupFaceIdx)[i].size(); j++)
		{
			HE_Face* f = pDual->_poly->face((*coarseGroupFaceIdx)[i][j]);
			int numFace = f->NumHalfEdge();

			for (int fi = 0; fi < numFace; fi++)
			{
				int coarseFi = faceIdxOffset[(*coarseGroupFaceIdx)[i][j]] + fi;
				if (f->TestFlag(VEF_FLAG_QUAD))
				{
					if (fi == 1 || fi == 2) faceList.push_back(coarseFi);
				}
				else if (f->TestFlag(VEF_FLAG_INC))
				{
					if (fi == 1) faceList.push_back(coarseFi);
				}
				else if (f->TestFlag(VEF_FLAG_DEC))
				{
					if (fi == 0) faceList.push_back(coarseFi);
					if (fi == 2) faceList.insert(faceList.end() - 1, coarseFi);
				}
				else if (f->TestFlag(VEF_FLAG_SHORTROW_L))
				{
					if (fi == 0 || fi == 1) faceList.push_back(coarseFi);
				}
				else if (f->TestFlag(VEF_FLAG_SHORTROW_R))
				{
					if (fi == 1 || fi == 2) faceList.push_back(coarseFi);
				}
				else
					std::cout << "Error: unknown face flag - reorderSubdivisionGroupFaceIdx\n";
			}
		}

		if (pDual->_isLoop[i]) { _groupFaceIdx.push_back(faceList); faceList.clear(); }
	
		for (int j = (int)(*coarseGroupFaceIdx)[i].size() - 1; j >= 0 ; j--)
		{
			HE_Face* f = pDual->_poly->face((*coarseGroupFaceIdx)[i][j]);
			int numFace = f->NumHalfEdge();

			for (int fi = numFace - 1; fi >= 0; fi--)
			{
				int coarseFi = faceIdxOffset[(*coarseGroupFaceIdx)[i][j]] + fi;
				if (f->TestFlag(VEF_FLAG_QUAD))
				{
					if (fi == 0 || fi == 3) faceList.push_back(coarseFi);
				}
				else if (f->TestFlag(VEF_FLAG_INC))
				{
					if (fi == 0 || fi == 2) faceList.push_back(coarseFi);
				}
				else if (f->TestFlag(VEF_FLAG_DEC))
				{
					if (fi == 1) faceList.push_back(coarseFi);
				}
				else if (f->TestFlag(VEF_FLAG_SHORTROW_L))
				{
					if (fi == 2) faceList.push_back(coarseFi);
				}
				else if (f->TestFlag(VEF_FLAG_SHORTROW_R))
				{
					if (fi == 0) faceList.push_back(coarseFi);
				}
				else
					std::cout << "Error: unknown face flag - reorderSubdivisionGroupFaceIdx\n";
			}
		}

		_groupFaceIdx.push_back(faceList);
	}
}

void DualGraph::findLoops()
{
	for (int gi = 0; gi < (int)_groupFaceIdx.size(); gi++)
	{
		bool isloop = false;

		int startFi = _groupFaceIdx[gi][0];
		int endFi = _groupFaceIdx[gi][_groupFaceIdx[gi].size() - 1];
		HE_Face* startF = _poly->face(startFi);
		HE_Face* endF = _poly->face(endFi);
		//if (startF->TestFlag(VEF_FLAG_QUAD) && endF->TestFlag(VEF_FLAG_QUAD))
		{
			HE_Face::const_edge_circulator e = startF->begin();
			int numHalfedge = startF->NumHalfEdge();
			for (int i = 0; i < numHalfedge; i++, e++)
			{
				if ((*e)->twin()->face()->index() == endF->index())
				{
					if ((*e)->TestFlag(EFLAG_VERTICAL) && (*e)->twin()->TestFlag(EFLAG_VERTICAL))
					{
						isloop = true;
					}
				}
			}
		}

		_isLoop.push_back(isloop);
	}
}

void DualGraph::fixUVMismatch()
{
	HE_HalfEdge* he1; HE_HalfEdge* he2;
	for (int i = 0; i < numEdges(); i++)
	{
		he1 = _poly->halfedge(_dualedges[i]._hei[0]);
		he2 = _poly->halfedge(_dualedges[i]._hei[1]);

		if (he1->TestFlag(EFLAG_VERTICAL) && he2->TestFlag(EFLAG_HORIZONTAL))
		{
			he1->ClearFlags(); he1->SetFlag(EFLAG_HORIZONTAL);
		}
		else if (he2->TestFlag(EFLAG_VERTICAL) && he1->TestFlag(EFLAG_HORIZONTAL))
		{
			he2->ClearFlags(); he2->SetFlag(EFLAG_HORIZONTAL);
		}
	}
}

bool DualGraph::CheckBadVertex()
{
	for (int vi = 0; vi < _poly->numVertices(); vi++)
	{
		HE_Vertex* v = _poly->vertex(vi);
		if (v->degree() <= 3)
		{
			return true;
		}
	}
	return false;
}

void DualGraph::flipEdgeUV()
{
	for (int fi = 0; fi < _poly->numFaces(); fi++)
	{
		HE_Face* f = _poly->face(fi);

		int HoriCount = 0;
		int VertCount = 0;
		HE_Face::const_edge_circulator e = f->begin();
		for (int ei = 0; ei < f->NumHalfEdge(); ei++)
		{
			if ((*e)->TestFlag(EFLAG_HORIZONTAL)) HoriCount++;
			if ((*e)->TestFlag(EFLAG_VERTICAL)) VertCount++;
			e++;
		}

		if (HoriCount == 0 || HoriCount == 3 || VertCount == 0 || VertCount == 3)
		{
			HE_HalfEdge* he = f->edge();
			for (int ei = 0; ei < 3; ei++)
			{
				if (he->twin()->face()->NumHalfEdge() == 3 &&
					(he->twin()->next()->twin()->face()->NumHalfEdge() == 3 || he->twin()->prev()->twin()->face()->NumHalfEdge() == 3))
				{
					he->ClearFlags();
					he->twin()->ClearFlags();
					if (he->next()->TestFlag(EFLAG_HORIZONTAL))
					{
						he->SetFlag(EFLAG_VERTICAL);
						he->twin()->SetFlag(EFLAG_VERTICAL);
					}
					else {
						he->SetFlag(EFLAG_HORIZONTAL);
						he->twin()->SetFlag(EFLAG_HORIZONTAL);
					}

					fixTriangleUV(he->twin()->face()->index(), he->twin()->face()->HalfEdgeIdx(he->twin()));
					break;
				}
				he = he->next();
			}
		}
	}
}

void DualGraph::flipTriangleUV()
{
	for (int fi = 0; fi < _poly->numFaces(); fi++)
	{
		HE_Face* f = _poly->face(fi);
		if (f->hole()) continue;

		bool hasQuadNeighbor = false;
		if (f->NumHalfEdge() == 3)
		{
			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;
			bool neighborIsQuad[3];
			int i = 0;
			do 
			{
				if ((*e)->twin()->face()->NumHalfEdge() == 4)
				{
					hasQuadNeighbor = true;
					neighborIsQuad[i] = true;
				}
				else
				{
					neighborIsQuad[i] = false;
				}
				e++; i++;
			} while (e != sentinel);

			if (hasQuadNeighbor)
			{
				HE_HalfEdge* he = f->edge();
				for (int ei = 0; ei < 3; ei++)
				{
					if (!neighborIsQuad[ei])
					{
						if (he->twin()->next()->TestFlag(EFLAG_HORIZONTAL) != he->next()->TestFlag(EFLAG_HORIZONTAL) &&
							he->twin()->prev()->TestFlag(EFLAG_HORIZONTAL) == he->prev()->TestFlag(EFLAG_HORIZONTAL))
						{
							if (he->twin()->next()->twin()->face()->NumHalfEdge() == 3 &&
								(he->twin()->next()->twin()->next()->twin()->face()->NumHalfEdge() == 3 ||
									he->twin()->next()->twin()->prev()->twin()->face()->NumHalfEdge() == 3))
							{
								he->twin()->next()->ClearFlags();
								he->twin()->next()->twin()->ClearFlags();

								if (he->next()->TestFlag(EFLAG_HORIZONTAL))
								{
									he->twin()->next()->SetFlag(EFLAG_HORIZONTAL);
									he->twin()->next()->twin()->SetFlag(EFLAG_HORIZONTAL);
								}
								else if (he->next()->TestFlag(EFLAG_VERTICAL))
								{
									he->twin()->next()->SetFlag(EFLAG_VERTICAL);
									he->twin()->next()->twin()->SetFlag(EFLAG_VERTICAL);
								}

								fixTriangleUV(he->twin()->next()->twin()->face()->index(),
									he->twin()->next()->twin()->face()->HalfEdgeIdx(he->twin()->next()->twin()));
							}
						}
						else if (he->twin()->prev()->TestFlag(EFLAG_HORIZONTAL) != he->prev()->TestFlag(EFLAG_HORIZONTAL) &&
							he->twin()->next()->TestFlag(EFLAG_HORIZONTAL) == he->next()->TestFlag(EFLAG_HORIZONTAL))
						{
							if (he->twin()->prev()->twin()->face()->NumHalfEdge() == 3 && 
								(he->twin()->prev()->twin()->next()->twin()->face()->NumHalfEdge() == 3 || 
									he->twin()->prev()->twin()->prev()->twin()->face()->NumHalfEdge() == 3))
							{
								he->twin()->prev()->ClearFlags();
								he->twin()->prev()->twin()->ClearFlags();

								if (he->prev()->TestFlag(EFLAG_HORIZONTAL))
								{
									he->twin()->prev()->SetFlag(EFLAG_HORIZONTAL);
									he->twin()->prev()->twin()->SetFlag(EFLAG_HORIZONTAL);
								}
								else if (he->prev()->TestFlag(EFLAG_VERTICAL))
								{
									he->twin()->prev()->SetFlag(EFLAG_VERTICAL);
									he->twin()->prev()->twin()->SetFlag(EFLAG_VERTICAL);
								}

								fixTriangleUV(he->twin()->prev()->twin()->face()->index(),
									he->twin()->prev()->twin()->face()->HalfEdgeIdx(he->twin()->prev()->twin()));
							}
						}
						else if (he->twin()->prev()->TestFlag(EFLAG_HORIZONTAL) != he->prev()->TestFlag(EFLAG_HORIZONTAL) &&
							he->twin()->next()->TestFlag(EFLAG_HORIZONTAL) != he->next()->TestFlag(EFLAG_HORIZONTAL))
						{
							if (he->twin()->prev()->twin()->face()->NumHalfEdge() == 3 &&
								(he->twin()->prev()->twin()->next()->twin()->face()->NumHalfEdge() == 3 ||
									he->twin()->prev()->twin()->prev()->twin()->face()->NumHalfEdge() == 3) &&
								he->twin()->next()->twin()->face()->NumHalfEdge() == 3 &&
								(he->twin()->next()->twin()->next()->twin()->face()->NumHalfEdge() == 3 ||
									he->twin()->next()->twin()->prev()->twin()->face()->NumHalfEdge() == 3)
								)
							{
								he->twin()->prev()->ClearFlags();
								he->twin()->prev()->twin()->ClearFlags();

								if (he->prev()->TestFlag(EFLAG_HORIZONTAL))
								{
									he->twin()->prev()->SetFlag(EFLAG_HORIZONTAL);
									he->twin()->prev()->twin()->SetFlag(EFLAG_HORIZONTAL);
								}
								else if (he->prev()->TestFlag(EFLAG_VERTICAL))
								{
									he->twin()->prev()->SetFlag(EFLAG_VERTICAL);
									he->twin()->prev()->twin()->SetFlag(EFLAG_VERTICAL);
								}

								fixTriangleUV(he->twin()->prev()->twin()->face()->index(),
									he->twin()->prev()->twin()->face()->HalfEdgeIdx(he->twin()->prev()->twin()));

								he->twin()->next()->ClearFlags();
								he->twin()->next()->twin()->ClearFlags();

								if (he->next()->TestFlag(EFLAG_HORIZONTAL))
								{
									he->twin()->next()->SetFlag(EFLAG_HORIZONTAL);
									he->twin()->next()->twin()->SetFlag(EFLAG_HORIZONTAL);
								}
								else if (he->next()->TestFlag(EFLAG_VERTICAL))
								{
									he->twin()->next()->SetFlag(EFLAG_VERTICAL);
									he->twin()->next()->twin()->SetFlag(EFLAG_VERTICAL);
								}

								fixTriangleUV(he->twin()->next()->twin()->face()->index(),
									he->twin()->next()->twin()->face()->HalfEdgeIdx(he->twin()->next()->twin()));
							}
						}
						//break;
					}
					he = he->next();
				}
			}
		}
	}
}

void DualGraph::fixTriangleUV(int fi, int i)
{
	HE_Face* f = _poly->face(fi);
	int HoriCount = 0;
	int VertCount = 0;
	HE_Face::const_edge_circulator e = f->begin();
	for (int ei = 0; ei < f->NumHalfEdge(); ei++)
	{
		if ((*e)->TestFlag(EFLAG_HORIZONTAL)) HoriCount++;
		if ((*e)->TestFlag(EFLAG_VERTICAL)) VertCount++;
		e++;
	}

	if (HoriCount == 0 || HoriCount == 3 || VertCount == 0 || VertCount == 3)
	{
		HE_Face::edge_circulator ee = f->begin();
		for (int ei = 0; ei < f->NumHalfEdge(); ei++)
		{
			if (ei == i) { ee++;  continue; }
			
			if ((*ee)->twin()->face()->NumHalfEdge() == 3 && 
				((*ee)->twin()->next()->twin()->face()->NumHalfEdge() == 3 || (*ee)->twin()->prev()->twin()->face()->NumHalfEdge() == 3))
			{
				if ((*ee)->TestFlag(EFLAG_HORIZONTAL))
				{
					(*ee)->ClearFlags();
					(*ee)->SetFlag(EFLAG_VERTICAL);
					(*ee)->twin()->ClearFlags();
					(*ee)->twin()->SetFlag(EFLAG_VERTICAL);
				} 
				else if ((*ee)->TestFlag(EFLAG_VERTICAL))
				{
					(*ee)->ClearFlags();
					(*ee)->SetFlag(EFLAG_HORIZONTAL);
					(*ee)->twin()->ClearFlags();
					(*ee)->twin()->SetFlag(EFLAG_HORIZONTAL);
				}

				fixTriangleUV((*ee)->twin()->face()->index(), (*ee)->twin()->face()->HalfEdgeIdx((*ee)->twin()));
				break;
			}

			ee++;
		}
	}
}