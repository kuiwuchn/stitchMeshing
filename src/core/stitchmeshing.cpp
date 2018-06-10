#include "../hierarchy.h"

void MultiResolutionHierarchy::convert2Poly()
{
	std::cout << "Convert 2 poly...";
	mPoly = new HE_Polyhedron(mV_tag, F_tag);

	mDual = new DualGraph(mPoly);

	//mDual->BadVertexCheck();
	std::cout << "Done\n";
}

void MultiResolutionHierarchy::stitchMeshing()
{
	//////////////////////////////////////////////////////////////////////////
	std::cout << "------------ cut bad quad ------------\n";

	mDual->findBadQuads();
	std::vector<HE_Vertex> cutVerts;
	std::vector<std::vector<int>> cutFaces;
	mDual->cut(cutVerts, cutFaces);

	delete mPoly;
	delete mDual;

	mPoly = new HE_Polyhedron(cutVerts, cutFaces);
	mDual = new DualGraph(mPoly);

	//////////////////////////////////////////////////////////////////////////
#if 1
	std::cout << "------------ remove bad vertex ------------\n";
	std::vector<HE_Vertex> mergeVerts;
	std::vector<std::vector<int>> mergeFaces;
	mDual->mergeTriangles(mergeVerts, mergeFaces);

	delete mPoly;
	delete mDual;

	mPoly = new HE_Polyhedron(mergeVerts, mergeFaces);
	mDual = new DualGraph(mPoly);
#endif

	string filename = "gurobi_result.txt";

	std::cout << "------------ UV minimization ------------\n";
	mDual->gurobiSolver(filename);
	mDual->loadGurobiResult(filename);
	mDual->findUVMismatch();

	//////////////////////////////////////////////////////////////////////////
	std::cout << "------------ cut UV mismatch quad ------------\n";
#if 1
	std::vector<HE_Vertex> cut2ndVerts;
	std::vector<std::vector<int>> cut2ndFaces;
	std::vector<std::vector<bool>> cut2ndEdgeFlags;
	mDual->cutUVMismatchQuad(cut2ndVerts, cut2ndFaces, cut2ndEdgeFlags);

	delete mPoly;
	delete mDual;

	mPoly = new HE_Polyhedron(cut2ndVerts, cut2ndFaces);
	mDual = new DualGraph(mPoly);
	mPoly->setEdgeFlags(cut2ndEdgeFlags);

	mDual->fixUVMismatch();
	
	mDual->findFaceGroups();
	mDual->fitFaceEdgesOrder();
	mDual->labelTriangleFace();

#if 1
	std::vector<HE_Vertex> flipVerts;
	std::vector<std::vector<int>> flipFaces;
	std::vector<std::vector<bool>> flipEdgeFlags;
	mDual->flipTriangles(flipVerts, flipFaces, flipEdgeFlags);

	delete mPoly;
	delete mDual;

	mPoly = new HE_Polyhedron(flipVerts, flipFaces);
	mDual = new DualGraph(mPoly);
	mPoly->setEdgeFlags(flipEdgeFlags);
#endif

	mDual->HVCheck();
	mDual->BadVertexCheck();
	mDual->findUVMismatch();
#endif

	//////////////////////////////////////////////////////////////////////////
	std::cout << "------------ break INC/DEC ------------\n";
#if 1
	mDual->breakConnectIncOrDec();
	mDual->findFaceGroups();
	mDual->fitFaceEdgesOrder();
	mDual->labelTriangleFace();
	mDual->HVCheck();
	mDual->BadVertexCheck();
#endif
	////////////////////////////////////////////////////////////////////////
	std::cout << "------------ flip triangle UV ------------\n";
#if 1
	mDual->flipTriangleUV();
	mDual->HVCheck();
#endif
	////////////////////////////////////////////////////////////////////////
	std::cout << "------------ flip Edge UV ------------\n";
#if 1
	mDual->flipEdgeUV();
	mDual->HVCheck();
	mDual->findFaceGroups();
	mDual->fitFaceEdgesOrder();
	mDual->labelTriangleFace();
	mDual->BadVertexCheck();
#endif
	////////////////////////////////////////////////////////////////////////
	std::cout << "------------ merge triangles ------------\n";
#if 1
	std::vector<HE_Vertex> mergeVerts_2;
	std::vector<std::vector<int>> mergeFaces_2;
	std::vector<std::vector<bool>> mergeEdgeFlags_2;
	mDual->mergeTriangles(mergeVerts_2, mergeFaces_2, mergeEdgeFlags_2);

	delete mPoly;
	delete mDual;

	mPoly = new HE_Polyhedron(mergeVerts_2, mergeFaces_2);
	mDual = new DualGraph(mPoly);
	mPoly->setEdgeFlags(mergeEdgeFlags_2);
	
	mDual->findFaceGroups();
	mDual->fitFaceEdgesOrder();
	mDual->labelTriangleFace();
	
	mDual->HVCheck();
	mDual->BadVertexCheck();
#endif

	////////////////////////////////////////////////////////////////////////
	std::cout << "------------ break INC/DEC again ------------\n";
#if 1	
	mDual->breakConnectIncOrDec();
#endif

	std::cout << "------------------------------------------ labeling\n";

#if 1
	mDual->findFaceGroups();
	mDual->fitFaceEdgesOrder();
	mDual->labelTriangleFace();
	mDual->HVCheck();
	mDual->BadVertexCheck();
#endif

	//std::string outputQuadFilename = filename + "_uv_final_mesh.obj";
	//gDual->exportUVMesh(outputQuadFilename.c_str());

	////////////////////////////////////////////////////////////////////////
	std::cout << "------------ wale mismatch minimization ------------\n";
#if 1
	mDual->findTopBtmGroupId();
	mDual->waleMismatchSolver();
	mDual->findLoops();
	mDual->findUVMismatch();
	mDual->findWaleMismatch();
	mDual->HVCheck();
	mDual->BadVertexCheck();
#endif

	std::cout << "------------------------------------------ alignment\n";

	//////////////////////////////////////////////////////////////////////////
	std::cout << "------------ subdivision ------------\n";
#if 1
	std::vector<HE_Vertex> subVerts;
	std::vector<std::vector<int>> subFaces;

	mDual->subdivide(subVerts, subFaces);
	mSubPoly = new HE_Polyhedron(subVerts, subFaces);
	mSubDual = new DualGraph(mSubPoly);
	mSubDual->setSubGroupFace(mDual);
	mSubDual->setGroupEdgeFlag();
	mSubDual->findUVMismatch();
	mSubDual->findWaleMismatch();
	mSubDual->reorderSubdivisionGroupFaceIdx(mDual);
	mSubDual->findLoops();
#endif

	//////////////////////////////////////////////////////////////////////////
	std::cout << "------------ merge triangles to INC/DEC ------------\n";
#if 1
	removeQuadDecInc();
#endif
}

void MultiResolutionHierarchy::convert2Rend()
{
	std::cout << "Convert 2 render buffer...";

#if 1
	mV_StMesh_rend.resize(3, mCleanPoly->numVertices());
	for (int vi = 0; vi < mCleanPoly->numVertices(); vi++)
	{
		mV_StMesh_rend(0, vi) = mCleanPoly->vertex(vi)->position().x;
		mV_StMesh_rend(1, vi) = mCleanPoly->vertex(vi)->position().y;
		mV_StMesh_rend(2, vi) = mCleanPoly->vertex(vi)->position().z;
	}

	int triNum = 0, quadNum = 0, penNum = 0;
	for (int fi = 0; fi < mCleanPoly->numFaces(); fi++) {
		HE_Face* f = mCleanPoly->face(fi);
		if (f->hole()) continue;

		if (f->NumHalfEdge() == 3) triNum++;
		else if (f->NumHalfEdge() == 4) quadNum++;
		else if (f->NumHalfEdge() == 5) penNum++;
	}
	

	int triTotalNum = triNum + 2 * quadNum + 3 * penNum;

	std::cout << triTotalNum << std::endl;

	mF_StMesh_rend.resize(3, triTotalNum);

	std::cout << mCleanPoly->numFaces() << std::endl;

	int c = 0;
	for (int fi = 0; fi < mCleanPoly->numFaces(); fi++)
	{
		HE_Face* f = mCleanPoly->face(fi);
		if (f->hole()) continue;

		std::vector<int> viList; viList.clear();
		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;
		do
		{
			viList.push_back((*e)->dst()->index());
			++e;
		} while (e != sentinel);

		if (f->NumHalfEdge() == 3)
		{
			mF_StMesh_rend(0, c) = viList[0];
			mF_StMesh_rend(1, c) = viList[1];
			mF_StMesh_rend(2, c++) = viList[2];
		}
		else if (f->NumHalfEdge() == 4)
		{
			mF_StMesh_rend(0, c) = viList[0];
			mF_StMesh_rend(1, c) = viList[1];
			mF_StMesh_rend(2, c++) = viList[3];

			mF_StMesh_rend(0, c) = viList[1];
			mF_StMesh_rend(1, c) = viList[2];
			mF_StMesh_rend(2, c++) = viList[3];
		}
		else if (f->NumHalfEdge() == 5)
		{
			mF_StMesh_rend(0, c) = viList[0];
			mF_StMesh_rend(1, c) = viList[1];
			mF_StMesh_rend(2, c++) = viList[4];

			mF_StMesh_rend(0, c) = viList[1];
			mF_StMesh_rend(1, c) = viList[2];
			mF_StMesh_rend(2, c++) = viList[4];

			mF_StMesh_rend(0, c) = viList[2];
			mF_StMesh_rend(1, c) = viList[3];
			mF_StMesh_rend(2, c++) = viList[4];
		}
	}
	std::cout << c << std::endl;
#else
	mV_StMesh_rend = mV_tag;

	int triNum = 0, quadNum = 0, penNum = 0;
	for (int i = 0; i < F_tag.size(); i++) {
		if (F_tag[i].size() == 3) triNum++;
		else if (F_tag[i].size() == 4) quadNum++;
		else if (F_tag[i].size() == 5) penNum++;
	}

	int triTotalNum = triNum + 2 * quadNum + 3 * penNum;

	mF_StMesh_rend.resize(3, triTotalNum);

	int c = 0;
	for (int i = 0; i < F_tag.size(); i++)
	{
		if (F_tag[i].size() == 3)
		{
			mF_StMesh_rend(0, c) = F_tag[i][0];
			mF_StMesh_rend(1, c) = F_tag[i][1];
			mF_StMesh_rend(2, c++) = F_tag[i][2];
		}
		else if (F_tag[i].size() == 4)
		{
			mF_StMesh_rend(0, c) = F_tag[i][0];
			mF_StMesh_rend(1, c) = F_tag[i][1];
			mF_StMesh_rend(2, c++) = F_tag[i][3];

			mF_StMesh_rend(0, c) = F_tag[i][1];
			mF_StMesh_rend(1, c) = F_tag[i][2];
			mF_StMesh_rend(2, c++) = F_tag[i][3];
		}
		else if (F_tag[i].size() == 5)
		{
			mF_StMesh_rend(0, c) = F_tag[i][0];
			mF_StMesh_rend(1, c) = F_tag[i][1];
			mF_StMesh_rend(2, c++) = F_tag[i][4];

			mF_StMesh_rend(0, c) = F_tag[i][1];
			mF_StMesh_rend(1, c) = F_tag[i][2];
			mF_StMesh_rend(2, c++) = F_tag[i][4];

			mF_StMesh_rend(0, c) = F_tag[i][2];
			mF_StMesh_rend(1, c) = F_tag[i][3];
			mF_StMesh_rend(2, c++) = F_tag[i][4];
		}
	}
#endif
	std::cout << "Done\n";
}

void MultiResolutionHierarchy::removeQuadDecInc()
{
	/* create face idx mapping */
	std::vector<int> mappingSubToClean;
	std::vector<int> mappingCleanToSub;
	int accuFaceCount = 0;
	for (int fi = 0; fi < mSubPoly->numFaces(); fi++)
	{
		HE_Face* f = mSubPoly->face(fi);
		if (f->hole()) continue;
		if (f->TestFlag(VEF_FLAG_QUAD_INC) || f->TestFlag(VEF_FLAG_QUAD_DEC))
		{
			mappingSubToClean.push_back(-1);
			continue;
		}
		else
		{
			mappingSubToClean.push_back(accuFaceCount++);
			mappingCleanToSub.push_back(fi);
		}
	}

	/* set vertices positions */
	std::vector<HE_Vertex> cleanVerts;
	for (int vi = 0; vi < mSubPoly->numVertices(); vi++)
	{
		HE_Vertex v = HE_Vertex(mSubPoly->vertex(vi)->position(), vi);
		cleanVerts.push_back(v);
	}

	/* set faces index */
	std::vector<std::vector<int>> cleanFaces;
	for (int fi = 0; fi < mSubPoly->numFaces(); fi++)
	{
		HE_Face* f = mSubPoly->face(fi);
		if (f->hole()) continue;
		if (f->TestFlag(VEF_FLAG_QUAD_INC) || f->TestFlag(VEF_FLAG_QUAD_DEC))
		{
			continue;
		}
		else
		{
			std::vector<int> cornerList;
			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;

			do
			{
				cornerList.push_back((*e)->dst()->index());
				++e;
			} while (e != sentinel);

			cleanFaces.push_back(cornerList);
		}
	}

	/* split quad inc and dec */
	for (int fi = 0; fi < mSubPoly->numFaces(); fi++)
	{
		HE_Face* f = mSubPoly->face(fi);
		if (f->hole()) continue;
		//std::cout << fi << std::endl;
		if (f->TestFlag(VEF_FLAG_QUAD_INC))
		{
			int VertIdx = f->edge()->dst()->index();
			int f_left_idx = f->edge()->next()->next()->twin()->face()->index();
			int f_right_idx = f->edge()->prev()->twin()->face()->index();
			int mapped_f_left_idx = mappingSubToClean[f_left_idx];
			int mapped_f_right_idx = mappingSubToClean[f_right_idx];

			//if (mapped_f_left_idx < 0 || mapped_f_right_idx < 0) continue;

			cleanFaces[mapped_f_left_idx].insert(cleanFaces[mapped_f_left_idx].begin() + 3, VertIdx);
			cleanFaces[mapped_f_right_idx].insert(cleanFaces[mapped_f_right_idx].begin() + 1, VertIdx);
			std::rotate(cleanFaces[mapped_f_right_idx].begin(), cleanFaces[mapped_f_right_idx].begin() + 1, cleanFaces[mapped_f_right_idx].end());
		}
		else if (f->TestFlag(VEF_FLAG_QUAD_DEC))
		{
			int VertIdx = f->edge()->dst()->index();
			int f_right_idx = f->edge()->next()->next()->twin()->face()->index();
			int f_left_idx = f->edge()->prev()->twin()->face()->index();
			int mapped_f_left_idx = mappingSubToClean[f_left_idx];
			int mapped_f_right_idx = mappingSubToClean[f_right_idx];

			//if (mapped_f_left_idx < 0 || mapped_f_right_idx < 0) continue;
			//std::cout << mapped_f_left_idx << " " << mapped_f_right_idx << std::endl;
			cleanFaces[mapped_f_left_idx].insert(cleanFaces[mapped_f_left_idx].begin() + 3, VertIdx);
			cleanFaces[mapped_f_right_idx].insert(cleanFaces[mapped_f_right_idx].begin() + 1, VertIdx);
		}
		//else 
		//	std::cout << "ERROR! expected face flag " << fi << std::endl;
	}

	mCleanPoly = new HE_Polyhedron(cleanVerts, cleanFaces);
	mCleanDual = new DualGraph(mCleanPoly);

	/* set face flag */
	for (int fi = 0; fi < mSubPoly->numFaces(); fi++)
	{
		HE_Face* f = mSubPoly->face(fi);
		if (f->hole()) continue;
		if (f->TestFlag(VEF_FLAG_QUAD))
		{
			HE_Face* cleanF = mCleanPoly->face(mappingSubToClean[fi]);
			cleanF->SetFlag(VEF_FLAG_QUAD);
		}
		else if (f->TestFlag(VEF_FLAG_QUAD_SHORTROW_L))
		{
			HE_Face* cleanF = mCleanPoly->face(mappingSubToClean[fi]);
			cleanF->SetFlag(VEF_FLAG_QUAD_SHORTROW_L);
		}
		else if (f->TestFlag(VEF_FLAG_QUAD_SHORTROW_R))
		{
			HE_Face* cleanF = mCleanPoly->face(mappingSubToClean[fi]);
			cleanF->SetFlag(VEF_FLAG_QUAD_SHORTROW_R);
		}
	}

	for (int fi = 0; fi < mSubPoly->numFaces(); fi++)
	{
		HE_Face* f = mSubPoly->face(fi);
		if (f->hole()) continue;
		if (f->TestFlag(VEF_FLAG_QUAD_INC))
		{
			int f_left_idx = f->edge()->next()->next()->twin()->face()->index();
			int f_right_idx = f->edge()->prev()->twin()->face()->index();
			int mapped_f_left_idx = mappingSubToClean[f_left_idx];
			int mapped_f_right_idx = mappingSubToClean[f_right_idx];

			mCleanPoly->face(mapped_f_left_idx)->ClearFlag(VEF_FLAG_QUAD);
			mCleanPoly->face(mapped_f_left_idx)->SetFlag(VEF_FLAG_INC);
			mCleanPoly->face(mapped_f_right_idx)->ClearFlag(VEF_FLAG_QUAD);
			mCleanPoly->face(mapped_f_right_idx)->SetFlag(VEF_FLAG_INC);
		}
		else if (f->TestFlag(VEF_FLAG_QUAD_DEC))
		{
			int f_right_idx = f->edge()->next()->next()->twin()->face()->index();
			int f_left_idx = f->edge()->prev()->twin()->face()->index();
			int mapped_f_left_idx = mappingSubToClean[f_left_idx];
			int mapped_f_right_idx = mappingSubToClean[f_right_idx];

			mCleanPoly->face(mapped_f_left_idx)->ClearFlag(VEF_FLAG_QUAD);
			mCleanPoly->face(mapped_f_left_idx)->SetFlag(VEF_FLAG_DEC);
			mCleanPoly->face(mapped_f_right_idx)->ClearFlag(VEF_FLAG_QUAD);
			mCleanPoly->face(mapped_f_right_idx)->SetFlag(VEF_FLAG_DEC);
		}
	}

	/* set edge flag */
	for (int fi = 0; fi < mCleanPoly->numFaces(); fi++)
	{
		HE_Face* f = mCleanPoly->face(fi);
		if (f->hole()) continue;
		if (f->TestFlag(VEF_FLAG_QUAD))
		{
			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;
			int i = 0;
			do
			{
				HE_HalfEdge* he = mCleanPoly->halfedge((*e)->index());

				if (i == 0 || i == 2) he->SetFlag(EFLAG_HORIZONTAL);
				else if (i == 1 || i == 3) he->SetFlag(EFLAG_VERTICAL);

				++e; ++i;
			} while (e != sentinel);
		}
		else if (f->TestFlag(VEF_FLAG_QUAD_SHORTROW_L) || f->TestFlag(VEF_FLAG_QUAD_SHORTROW_R))
		{
			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;
			int i = 0;
			do
			{
				HE_HalfEdge* he = mCleanPoly->halfedge((*e)->index());

				if (i == 0 || i == 1) he->SetFlag(EFLAG_HORIZONTAL);
				else if (i == 2 || i == 3) he->SetFlag(EFLAG_VERTICAL);

				++e; ++i;
			} while (e != sentinel);
		}
		else if (f->TestFlag(VEF_FLAG_DEC))
		{
			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;
			int i = 0;
			do
			{
				HE_HalfEdge* he = mCleanPoly->halfedge((*e)->index());

				if (i == 0 || i == 2 || i == 3) he->SetFlag(EFLAG_HORIZONTAL);
				else if (i == 1 || i == 4) he->SetFlag(EFLAG_VERTICAL);

				++e; ++i;
			} while (e != sentinel);
		}
		else if (f->TestFlag(VEF_FLAG_INC))
		{
			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;
			int i = 0;
			do
			{
				HE_HalfEdge* he = mCleanPoly->halfedge((*e)->index());

				if (i == 0 || i == 2 || i == 4) he->SetFlag(EFLAG_HORIZONTAL);
				else if (i == 1 || i == 3) he->SetFlag(EFLAG_VERTICAL);

				++e; ++i;
			} while (e != sentinel);
		}
		else
		{
			std::cout << "ERROR! unknown face flag " << fi << std::endl;
		}
	}

	/* find mismatch */
	mCleanDual->findUVMismatch();
	mCleanDual->findWaleMismatch();

	/* set group face id */
	for (int j = 0; j < (int)mSubDual->_groupFaceIdx.size(); j++)
	{
		std::vector<int> groupIdx;
		for (int i = 0; i < (int)mSubDual->_groupFaceIdx[j].size(); i++)
		{
			if (mappingSubToClean[mSubDual->_groupFaceIdx[j][i]] > -1)
			{
				groupIdx.push_back(mappingSubToClean[mSubDual->_groupFaceIdx[j][i]]);
			}
		}
		mCleanDual->_groupFaceIdx.push_back(groupIdx);
	}

	/* find loops */
	mCleanDual->findLoops();
}

void MultiResolutionHierarchy::exportResult()
{
	std::string objFilename = "result_mesh.obj";
	std::ofstream exportFile;
	exportFile.open(objFilename.c_str());

	for (int vi = 0; vi < mCleanPoly->numVertices(); vi++)
	{
		exportFile << "v " << mCleanPoly->vertex(vi)->position().x << " " << mCleanPoly->vertex(vi)->position().y << " " << mCleanPoly->vertex(vi)->position().z << std::endl;
	}

	for (int fi = 0; fi < mCleanPoly->numFaces(); fi++)
	{
		HE_Face* f = mCleanPoly->face(fi);

		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;

		exportFile << "f ";
		do
		{
			exportFile << (*e)->dst()->index() + 1 << " ";
			++e;
		} while (e != sentinel);

		exportFile << std::endl;
	}
	exportFile.close();
}