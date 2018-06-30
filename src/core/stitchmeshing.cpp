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

#include "../hierarchy.h"

// https://stackoverflow.com/questions/2328339/how-to-generate-n-different-colors-for-any-natural-number-n
int indexcolors[128] = {
	0xEEEEEE, 0xFFFF00, 0x1CE6FF, 0xFF34FF, 0xFF4A46, 0x008941, 0x006FA6, 0xA30059,
	0xFFDBE5, 0x7A4900, 0x0000A6, 0x63FFAC, 0xB79762, 0x004D43, 0x8FB0FF, 0x997D87,
	0x5A0007, 0x809693, 0xFEFFE6, 0x1B4400, 0x4FC601, 0x3B5DFF, 0x4A3B53, 0xFF2F80,
	0x61615A, 0xBA0900, 0x6B7900, 0x00C2A0, 0xFFAA92, 0xFF90C9, 0xB903AA, 0xD16100,
	0xDDEFFF, 0x000035, 0x7B4F4B, 0xA1C299, 0x300018, 0x0AA6D8, 0x013349, 0x00846F,
	0x372101, 0xFFB500, 0xC2FFED, 0xA079BF, 0xCC0744, 0xC0B9B2, 0xC2FF99, 0x001E09,
	0x00489C, 0x6F0062, 0x0CBD66, 0xEEC3FF, 0x456D75, 0xB77B68, 0x7A87A1, 0x788D66,
	0x885578, 0xFAD09F, 0xFF8A9A, 0xD157A0, 0xBEC459, 0x456648, 0x0086ED, 0x886F4C,

	0x34362D, 0xB4A8BD, 0x00A6AA, 0x452C2C, 0x636375, 0xA3C8C9, 0xFF913F, 0x938A81,
	0x575329, 0x00FECF, 0xB05B6F, 0x8CD0FF, 0x3B9700, 0x04F757, 0xC8A1A1, 0x1E6E00,
	0x7900D7, 0xA77500, 0x6367A9, 0xA05837, 0x6B002C, 0x772600, 0xD790FF, 0x9B9700,
	0x549E79, 0xFFF69F, 0x201625, 0x72418F, 0xBC23FF, 0x99ADC0, 0x3A2465, 0x922329,
	0x5B4534, 0xFDE8DC, 0x404E55, 0x0089A3, 0xCB7E98, 0xA4E804, 0x324E72, 0x6A3A4C,
	0x83AB58, 0x001C1E, 0xD1F7CE, 0x004B28, 0xC8D0F6, 0xA3A489, 0x806C66, 0x222800,
	0xBF5650, 0xE83000, 0x66796D, 0xDA007C, 0xFF1A59, 0x8ADBB4, 0x1E0200, 0x5B4E51,
	0xC895C5, 0x320033, 0xFF6832, 0x66E1D3, 0xCFCDAC, 0xD0AC94, 0x7ED379, 0x012C58
};

void MultiResolutionHierarchy::convert2Poly()
{
	std::cout << "Convert 2 poly...";
	mPoly = new HE_Polyhedron(mV_tag, F_tag);

	mDual = new DualGraph(mPoly);

	std::cout << "Done\n";
}

void MultiResolutionHierarchy::labelMesh(bool pFlip)
{
	//////////////////////////////////////////////////////////////////////////
	//std::cout << "------------ cut bad quad ------------\n";

	//mDual->findBadQuads();
	//std::vector<HE_Vertex> cutVerts;
	//std::vector<std::vector<int>> cutFaces;
	//mDual->cut(cutVerts, cutFaces);

	//delete mPoly;
	//delete mDual;

	//mPoly = new HE_Polyhedron(cutVerts, cutFaces);
	//mDual = new DualGraph(mPoly);

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
	mDual->loadGurobiResult(filename, pFlip);
	mDual->findUVMismatch();

	//////////////////////////////////////////////////////////////////////////
	std::cout << "------------ cut UV mismatch quad ------------\n";
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

	std::vector<HE_Vertex> flipVerts;
	std::vector<std::vector<int>> flipFaces;
	std::vector<std::vector<bool>> flipEdgeFlags;
	mDual->flipTriangles(flipVerts, flipFaces, flipEdgeFlags);

	delete mPoly;
	delete mDual;

	mPoly = new HE_Polyhedron(flipVerts, flipFaces);
	mDual = new DualGraph(mPoly);
	mPoly->setEdgeFlags(flipEdgeFlags);

	mDual->findUVMismatch();

	//////////////////////////////////////////////////////////////////////////
	std::cout << "------------ break INC/DEC ------------\n";
	mDual->breakConnectIncOrDec();
	mDual->findFaceGroups();
	mDual->fitFaceEdgesOrder();
	mDual->labelTriangleFace();

	////////////////////////////////////////////////////////////////////////
	std::cout << "------------ flip triangle UV ------------\n";
	mDual->flipTriangleUV();

	////////////////////////////////////////////////////////////////////////
	std::cout << "------------ flip Edge UV ------------\n";
	mDual->flipEdgeUV();
	mDual->findFaceGroups();
	mDual->fitFaceEdgesOrder();
	mDual->labelTriangleFace();

	////////////////////////////////////////////////////////////////////////
	std::cout << "------------ merge triangles ------------\n";
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

	////////////////////////////////////////////////////////////////////////
	std::cout << "------------ break INC/DEC again ------------\n";
	mDual->breakConnectIncOrDec();
	std::cout << "------------------------------------------ labeling\n";
}

void MultiResolutionHierarchy::alignMesh()
{
	mDual->findFaceGroups();
	mDual->fitFaceEdgesOrder();
	mDual->labelTriangleFace();

	////////////////////////////////////////////////////////////////////////
	std::cout << "------------ wale mismatch minimization ------------\n";
	mDual->findTopBtmGroupId();
	mDual->waleMismatchSolver();
	mDual->findLoops();
	mDual->findUVMismatch();
	mDual->findWaleMismatch();

	std::cout << "------------------------------------------ alignment\n";
}

void MultiResolutionHierarchy::stitchMeshing()
{
	//////////////////////////////////////////////////////////////////////////
	std::cout << "------------ subdivision ------------\n";
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

	//////////////////////////////////////////////////////////////////////////
	std::cout << "------------ merge triangles to INC/DEC ------------\n";
	removeQuadDecInc();
}

cyPoint3f colorConverter(int hexValue)
{
	cyPoint3f rgbColor;
	rgbColor.x = ((hexValue >> 16) & 0xFF) / 255.0f;  // Extract the RR byte
	rgbColor.y = ((hexValue >> 8) & 0xFF) / 255.0f;   // Extract the GG byte
	rgbColor.z = ((hexValue) & 0xFF) / 255.0f;        // Extract the BB byte

	return rgbColor;
}

void MultiResolutionHierarchy::convertLabelMesh2Rend()
{
	int triNum = 0, quadNum = 0, penNum = 0;
	for (int fi = 0; fi < mPoly->numFaces(); fi++) {
		HE_Face* f = mPoly->face(fi);
		if (f->hole()) continue;

		if (f->NumHalfEdge() == 3) triNum++;
		else if (f->NumHalfEdge() == 4) quadNum++;
		else if (f->NumHalfEdge() == 5) penNum++;
	}

	int totalVNum = triNum * 3 + quadNum * 4 + penNum * 5;

	mV_LbMesh_rend.resize(3, totalVNum);
	mT_LbMesh_rend.resize(2, totalVNum);

	int vtCount = 0;
	for (int fi = 0; fi < mPoly->numFaces(); fi++)
	{
		HE_Face* f = mPoly->face(fi);
		if (f->hole()) continue;

		if (f->NumHalfEdge() == 4)
		{
			while (!f->edge()->TestFlag(EFLAG_VERTICAL))
			{
				f->edge(f->edge()->next());
			}

			mV_LbMesh_rend(0, vtCount + 0) = f->edge()->src()->position().x;
			mV_LbMesh_rend(1, vtCount + 0) = f->edge()->src()->position().y;
			mV_LbMesh_rend(2, vtCount + 0) = f->edge()->src()->position().z;
			mV_LbMesh_rend(0, vtCount + 1) = f->edge()->next()->src()->position().x;
			mV_LbMesh_rend(1, vtCount + 1) = f->edge()->next()->src()->position().y;
			mV_LbMesh_rend(2, vtCount + 1) = f->edge()->next()->src()->position().z;
			mV_LbMesh_rend(0, vtCount + 2) = f->edge()->next()->next()->src()->position().x;
			mV_LbMesh_rend(1, vtCount + 2) = f->edge()->next()->next()->src()->position().y;
			mV_LbMesh_rend(2, vtCount + 2) = f->edge()->next()->next()->src()->position().z;
			mV_LbMesh_rend(0, vtCount + 3) = f->edge()->next()->next()->next()->src()->position().x;
			mV_LbMesh_rend(1, vtCount + 3) = f->edge()->next()->next()->next()->src()->position().y;
			mV_LbMesh_rend(2, vtCount + 3) = f->edge()->next()->next()->next()->src()->position().z;

			mT_LbMesh_rend(0, vtCount + 1) = 0;
			mT_LbMesh_rend(1, vtCount + 1) = 0;
			mT_LbMesh_rend(0, vtCount + 2) = 0.5;
			mT_LbMesh_rend(1, vtCount + 2) = 0;
			mT_LbMesh_rend(0, vtCount + 3) = 0.5;
			mT_LbMesh_rend(1, vtCount + 3) = 1;
			mT_LbMesh_rend(0, vtCount + 0) = 0;
			mT_LbMesh_rend(1, vtCount + 0) = 1;
			vtCount += 4;
		}
		else if (f->NumHalfEdge() == 3)
		{
			int hCount = 0, vCount = 0;
			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;
			do
			{
				if ((*e)->TestFlag(EFLAG_VERTICAL)) vCount++;
				if ((*e)->TestFlag(EFLAG_HORIZONTAL)) hCount++;
				++e;
			} while (e != sentinel);

			if (hCount == 2 && vCount == 1)
			{
				while (!f->edge()->TestFlag(EFLAG_VERTICAL))
				{
					f->edge(f->edge()->next());
				}

				mV_LbMesh_rend(0, vtCount + 0) = f->edge()->src()->position().x;
				mV_LbMesh_rend(1, vtCount + 0) = f->edge()->src()->position().y;
				mV_LbMesh_rend(2, vtCount + 0) = f->edge()->src()->position().z;
				mV_LbMesh_rend(0, vtCount + 1) = f->edge()->next()->src()->position().x;
				mV_LbMesh_rend(1, vtCount + 1) = f->edge()->next()->src()->position().y;
				mV_LbMesh_rend(2, vtCount + 1) = f->edge()->next()->src()->position().z;
				mV_LbMesh_rend(0, vtCount + 2) = f->edge()->next()->next()->src()->position().x;
				mV_LbMesh_rend(1, vtCount + 2) = f->edge()->next()->next()->src()->position().y;
				mV_LbMesh_rend(2, vtCount + 2) = f->edge()->next()->next()->src()->position().z;

				mT_LbMesh_rend(0, vtCount + 1) = 0.5;
				mT_LbMesh_rend(1, vtCount + 1) = 1;
				mT_LbMesh_rend(0, vtCount + 2) = 1;
				mT_LbMesh_rend(1, vtCount + 2) = 1;
				mT_LbMesh_rend(0, vtCount + 0) = 1;
				mT_LbMesh_rend(1, vtCount + 0) = 0;
			}
			else if (hCount == 1 && vCount == 2)
			{
				while (!f->edge()->TestFlag(EFLAG_HORIZONTAL))
				{
					f->edge(f->edge()->next());
				}

				mV_LbMesh_rend(0, vtCount + 0) = f->edge()->src()->position().x;
				mV_LbMesh_rend(1, vtCount + 0) = f->edge()->src()->position().y;
				mV_LbMesh_rend(2, vtCount + 0) = f->edge()->src()->position().z;
				mV_LbMesh_rend(0, vtCount + 1) = f->edge()->next()->src()->position().x;
				mV_LbMesh_rend(1, vtCount + 1) = f->edge()->next()->src()->position().y;
				mV_LbMesh_rend(2, vtCount + 1) = f->edge()->next()->src()->position().z;
				mV_LbMesh_rend(0, vtCount + 2) = f->edge()->next()->next()->src()->position().x;
				mV_LbMesh_rend(1, vtCount + 2) = f->edge()->next()->next()->src()->position().y;
				mV_LbMesh_rend(2, vtCount + 2) = f->edge()->next()->next()->src()->position().z;

				mT_LbMesh_rend(0, vtCount + 1) = 1;
				mT_LbMesh_rend(1, vtCount + 1) = 0;
				mT_LbMesh_rend(0, vtCount + 2) = 0.5;
				mT_LbMesh_rend(1, vtCount + 2) = 0;
				mT_LbMesh_rend(0, vtCount + 0) = 0.5;
				mT_LbMesh_rend(1, vtCount + 0) = 1;
			}
			else std::cout << "ERROR: export uv mesh\n";
			vtCount += 3;
		}
		else std::cout << "ERROR: export uv mesh\n";
	}

	int triTotalNum = triNum + 2 * quadNum + 3 * penNum;

	mF_LbMesh_rend.resize(3, triTotalNum);

	int c = 0;
	int vi = 0;
	for (int fi = 0; fi < mPoly->numFaces(); fi++)
	{
		const HE_Face* f = mPoly->face(fi);

		if (f->hole()) continue;

		std::vector<int> viList; viList.clear();
		HE_Face::const_edge_circulator e = f->begin();
		HE_Face::const_edge_circulator sentinel = e;
		do
		{
			viList.push_back(vi);
			vi++;
			++e;
		} while (e != sentinel);

		if (f->NumHalfEdge() == 3)
		{
			mF_LbMesh_rend(0, c) = viList[0];
			mF_LbMesh_rend(1, c) = viList[1];
			mF_LbMesh_rend(2, c++) = viList[2];
		}
		else if (f->NumHalfEdge() == 4)
		{
			mF_LbMesh_rend(0, c) = viList[0];
			mF_LbMesh_rend(1, c) = viList[1];
			mF_LbMesh_rend(2, c++) = viList[3];
			   
			mF_LbMesh_rend(0, c) = viList[1];
			mF_LbMesh_rend(1, c) = viList[2];
			mF_LbMesh_rend(2, c++) = viList[3];
		}
		else if (f->NumHalfEdge() == 5)
		{
			mF_LbMesh_rend(0, c) = viList[0];
			mF_LbMesh_rend(1, c) = viList[1];
			mF_LbMesh_rend(2, c++) = viList[4];
			   
			mF_LbMesh_rend(0, c) = viList[1];
			mF_LbMesh_rend(1, c) = viList[2];
			mF_LbMesh_rend(2, c++) = viList[4];
			   
			mF_LbMesh_rend(0, c) = viList[2];
			mF_LbMesh_rend(1, c) = viList[3];
			mF_LbMesh_rend(2, c++) = viList[4];
		}
	}
}

void MultiResolutionHierarchy::convertAlignMesh2Rend()
{
	int triNum = 0, quadNum = 0, penNum = 0;
	for (int fi = 0; fi < mPoly->numFaces(); fi++) {
		HE_Face* f = mPoly->face(fi);
		if (f->hole()) continue;

		if (f->NumHalfEdge() == 3) triNum++;
		else if (f->NumHalfEdge() == 4) quadNum++;
		else if (f->NumHalfEdge() == 5) penNum++;
	}

	int totalVNum = triNum * 3 + quadNum * 4 + penNum * 5;

	mV_AlMesh_rend.resize(3, totalVNum);
	mC_AlMesh_rend.resize(3, totalVNum);
	mT_AlMesh_rend.resize(2, totalVNum);
 
	int vi = 0;
	for (int j = 0; j < (int)mDual->_groupFaceIdx.size(); j++)
	{
		cyPoint3f c = colorConverter(indexcolors[j % 128]);
		for (int i = 0; i < (int)mDual->_groupFaceIdx[j].size(); i++)
		{
			const HE_Face* f = mPoly->face(mDual->_groupFaceIdx[j][i]);
			if (f->hole())
				continue;

			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;

			int idx = 0;
			if (f->TestFlag(VEF_FLAG_QUAD))
			{
				do
				{
					const HE_Vertex* v = (*e)->dst();
					switch (idx)
					{
					case 0:
						mT_AlMesh_rend(0, vi) = 1;
						mT_AlMesh_rend(1, vi) = 0;
						break;
					case 1:
						mT_AlMesh_rend(0, vi) = 1;
						mT_AlMesh_rend(1, vi) = 1;
						break;
					case 2:
						mT_AlMesh_rend(0, vi) = 0;
						mT_AlMesh_rend(1, vi) = 1;
						break;
					case 3:
						mT_AlMesh_rend(0, vi) = 0;
						mT_AlMesh_rend(1, vi) = 0;
						break;
					default: break;
					}
					mV_AlMesh_rend(0, vi) = v->position().x;
					mV_AlMesh_rend(1, vi) = v->position().y;
					mV_AlMesh_rend(2, vi) = v->position().z;

					mC_AlMesh_rend(0, vi) = c.x;
					mC_AlMesh_rend(1, vi) = c.y;
					mC_AlMesh_rend(2, vi++) = c.z;
					++e; ++idx;
				} while (e != sentinel);
			}
			else if (f->TestFlag(VEF_FLAG_INC))
			{
				do
				{
					const HE_Vertex* v = (*e)->dst();
					switch (idx)
					{
					case 0:
						mT_AlMesh_rend(0, vi) = 1;
						mT_AlMesh_rend(1, vi) = 0;
						break;
					case 1:
						mT_AlMesh_rend(0, vi) = 1;
						mT_AlMesh_rend(1, vi) = 1;
						break;
					case 2:
						mT_AlMesh_rend(0, vi) = 0;
						mT_AlMesh_rend(1, vi) = 1;
						break;
					case 3:
						mT_AlMesh_rend(0, vi) = 0;
						mT_AlMesh_rend(1, vi) = 0;
						break;
					case 4:
						mT_AlMesh_rend(0, vi) = 0.5;
						mT_AlMesh_rend(1, vi) = 0;
						break;
					default: break;
					}
					mV_AlMesh_rend(0, vi) = v->position().x;
					mV_AlMesh_rend(1, vi) = v->position().y;
					mV_AlMesh_rend(2, vi) = v->position().z;
					mC_AlMesh_rend(0, vi) = c.x;
					mC_AlMesh_rend(1, vi) = c.y;
					mC_AlMesh_rend(2, vi++) = c.z;
					++e; ++idx;
				} while (e != sentinel);
			}
			else if (f->TestFlag(VEF_FLAG_DEC))
			{
				do
				{
					const HE_Vertex* v = (*e)->dst();
					switch (idx)
					{
					case 0:
						mT_AlMesh_rend(0, vi) = 1;
						mT_AlMesh_rend(1, vi) = 0;
						break;
					case 1:
						mT_AlMesh_rend(0, vi) = 1;
						mT_AlMesh_rend(1, vi) = 1;
						break;
					case 2:
						mT_AlMesh_rend(0, vi) = 0.5;
						mT_AlMesh_rend(1, vi) = 1;
						break;
					case 3:
						mT_AlMesh_rend(0, vi) = 0;
						mT_AlMesh_rend(1, vi) = 1;
						break;
					case 4:
						mT_AlMesh_rend(0, vi) = 0;
						mT_AlMesh_rend(1, vi) = 0;
						break;
					default: break;
					}
					mV_AlMesh_rend(0, vi) = v->position().x;
					mV_AlMesh_rend(1, vi) = v->position().y;
					mV_AlMesh_rend(2, vi) = v->position().z;
					mC_AlMesh_rend(0, vi) = c.x;
					mC_AlMesh_rend(1, vi) = c.y;
					mC_AlMesh_rend(2, vi++) = c.z;
					++e; ++idx;
				} while (e != sentinel);
			}
			else if (f->TestFlag(VEF_FLAG_QUAD_SHORTROW_L) ||
				f->TestFlag(VEF_FLAG_QUAD_SHORTROW_R))
			{
				do
				{
					const HE_Vertex* v = (*e)->dst();
					mT_AlMesh_rend(0, vi) = 0;
					mT_AlMesh_rend(1, vi) = 0;
					mV_AlMesh_rend(0, vi) = v->position().x;
					mV_AlMesh_rend(1, vi) = v->position().y;
					mV_AlMesh_rend(2, vi) = v->position().z;
					mC_AlMesh_rend(0, vi) = c.x;
					mC_AlMesh_rend(1, vi) = c.y;
					mC_AlMesh_rend(2, vi++) = c.z;
					++e;
				} while (e != sentinel);
			}
			else {
				do
				{
					const HE_Vertex* v = (*e)->dst();
					mT_AlMesh_rend(0, vi) = 0;
					mT_AlMesh_rend(1, vi) = 0;
					mV_AlMesh_rend(0, vi) = v->position().x;
					mV_AlMesh_rend(1, vi) = v->position().y;
					mV_AlMesh_rend(2, vi) = v->position().z;
					mC_AlMesh_rend(0, vi) = c.x;
					mC_AlMesh_rend(1, vi) = c.y;
					mC_AlMesh_rend(2, vi++) = c.z;
					++e;
				} while (e != sentinel);
			}
		}
	}

	int triTotalNum = triNum + 2 * quadNum + 3 * penNum;

	mF_AlMesh_rend.resize(3, triTotalNum);

	int c = 0;
	vi = 0;
	for (int j = 0; j < (int)mDual->_groupFaceIdx.size(); j++)
	{
		for (int i = 0; i < (int)mDual->_groupFaceIdx[j].size(); i++)
		{
			const HE_Face* f = mPoly->face(mDual->_groupFaceIdx[j][i]);

			if (f->hole()) continue;

			std::vector<int> viList; viList.clear();
			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;
			do
			{
				viList.push_back(vi);
				vi++;
				++e;
			} while (e != sentinel);

			if (f->NumHalfEdge() == 3)
			{
				mF_AlMesh_rend(0, c) = viList[0];
				mF_AlMesh_rend(1, c) = viList[1];
				mF_AlMesh_rend(2, c++) = viList[2];
			}
			else if (f->NumHalfEdge() == 4)
			{
				mF_AlMesh_rend(0, c) = viList[0];
				mF_AlMesh_rend(1, c) = viList[1];
				mF_AlMesh_rend(2, c++) = viList[3];

				mF_AlMesh_rend(0, c) = viList[1];
				mF_AlMesh_rend(1, c) = viList[2];
				mF_AlMesh_rend(2, c++) = viList[3];
			}
			else if (f->NumHalfEdge() == 5)
			{
				mF_AlMesh_rend(0, c) = viList[0];
				mF_AlMesh_rend(1, c) = viList[1];
				mF_AlMesh_rend(2, c++) = viList[4];
				   
				mF_AlMesh_rend(0, c) = viList[1];
				mF_AlMesh_rend(1, c) = viList[2];
				mF_AlMesh_rend(2, c++) = viList[4];
				   
				mF_AlMesh_rend(0, c) = viList[2];
				mF_AlMesh_rend(1, c) = viList[3];
				mF_AlMesh_rend(2, c++) = viList[4];
			}
		}
	}

	mE_AlMesh_rend.resize(6, mPoly->numHalfEdges() * 2);

	for (int ei = 0; ei < mPoly->numHalfEdges(); ei++)
	{
		HE_HalfEdge* he = mPoly->halfedge(ei);
		mE_AlMesh_rend(0, ei * 2 + 0) = he->src()->position().x;
		mE_AlMesh_rend(1, ei * 2 + 0) = he->src()->position().y;
		mE_AlMesh_rend(2, ei * 2 + 0) = he->src()->position().z;
		mE_AlMesh_rend(3, ei * 2 + 0) = 0;
		mE_AlMesh_rend(4, ei * 2 + 0) = 0;
		mE_AlMesh_rend(5, ei * 2 + 0) = 0;
		   
		mE_AlMesh_rend(0, ei * 2 + 1) = he->dst()->position().x;
		mE_AlMesh_rend(1, ei * 2 + 1) = he->dst()->position().y;
		mE_AlMesh_rend(2, ei * 2 + 1) = he->dst()->position().z;
		mE_AlMesh_rend(3, ei * 2 + 1) = 0;
		mE_AlMesh_rend(4, ei * 2 + 1) = 0;
		mE_AlMesh_rend(5, ei * 2 + 1) = 0;
	}
}

void MultiResolutionHierarchy::convertStitchMesh2Rend()
{
	std::cout << "Convert 2 render buffer...";

#if 0
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

	mF_StMesh_rend.resize(3, triTotalNum);

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
#else
	int triNum = 0, quadNum = 0, penNum = 0;
	for (int fi = 0; fi < mCleanPoly->numFaces(); fi++) {
		HE_Face* f = mCleanPoly->face(fi);
		if (f->hole()) continue;

		if (f->NumHalfEdge() == 3) triNum++;
		else if (f->NumHalfEdge() == 4) quadNum++;
		else if (f->NumHalfEdge() == 5) penNum++;
	}

	int totalVNum = triNum * 3 + quadNum * 4 + penNum * 5;

	mV_StMesh_rend.resize(3, totalVNum);
	mC_StMesh_rend.resize(3, totalVNum);
	mT_StMesh_rend.resize(2, totalVNum);
	//for (int vi = 0; vi < mCleanPoly->numVertices(); vi++)
	//{
	//	mV_StMesh_rend(0, vi) = mCleanPoly->vertex(vi)->position().x;
	//	mV_StMesh_rend(1, vi) = mCleanPoly->vertex(vi)->position().y;
	//	mV_StMesh_rend(2, vi) = mCleanPoly->vertex(vi)->position().z;
	//}
	int vi = 0;
	for (int j = 0; j < (int)mCleanDual->_groupFaceIdx.size(); j++)
	{
		cyPoint3f c = colorConverter(indexcolors[j % 128]);
		for (int i = 0; i < (int)mCleanDual->_groupFaceIdx[j].size(); i++)
		{
			const HE_Face* f = mCleanPoly->face(mCleanDual->_groupFaceIdx[j][i]);
			if (f->hole())
				continue;

			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;

			int idx = 0;
			if (f->TestFlag(VEF_FLAG_QUAD))
			{
				do
				{
					const HE_Vertex* v = (*e)->dst();
					switch (idx)
					{
					case 0:
						mT_StMesh_rend(0, vi) = 1;
						mT_StMesh_rend(1, vi) = 0;
						break;
					case 1:
						mT_StMesh_rend(0, vi) = 1;
						mT_StMesh_rend(1, vi) = 1;
						break;
					case 2:
						mT_StMesh_rend(0, vi) = 0;
						mT_StMesh_rend(1, vi) = 1;
						break;
					case 3:
						mT_StMesh_rend(0, vi) = 0;
						mT_StMesh_rend(1, vi) = 0;
						break;
					default: break;
					}
					mV_StMesh_rend(0, vi) = v->position().x;
					mV_StMesh_rend(1, vi) = v->position().y;
					mV_StMesh_rend(2, vi) = v->position().z;

					mC_StMesh_rend(0, vi) = c.x;
					mC_StMesh_rend(1, vi) = c.y;
					mC_StMesh_rend(2, vi++) = c.z;
					++e; ++idx;
				} while (e != sentinel);
			}
			else if (f->TestFlag(VEF_FLAG_INC))
			{
				do
				{
					const HE_Vertex* v = (*e)->dst();
					switch (idx)
					{
					case 0:
						mT_StMesh_rend(0, vi) = 1;
						mT_StMesh_rend(1, vi) = 0;
						break;
					case 1:
						mT_StMesh_rend(0, vi) = 1;
						mT_StMesh_rend(1, vi) = 1;
						break;
					case 2:
						mT_StMesh_rend(0, vi) = 0;
						mT_StMesh_rend(1, vi) = 1;
						break;
					case 3:
						mT_StMesh_rend(0, vi) = 0;
						mT_StMesh_rend(1, vi) = 0;
						break;
					case 4:
						mT_StMesh_rend(0, vi) = 0.5;
						mT_StMesh_rend(1, vi) = 0;
						break;
					default: break;
					}
					mV_StMesh_rend(0, vi) = v->position().x;
					mV_StMesh_rend(1, vi) = v->position().y;
					mV_StMesh_rend(2, vi) = v->position().z;
					mC_StMesh_rend(0, vi) = c.x;
					mC_StMesh_rend(1, vi) = c.y;
					mC_StMesh_rend(2, vi++) = c.z;
					++e; ++idx;
				} while (e != sentinel);
			}
			else if (f->TestFlag(VEF_FLAG_DEC))
			{
				do
				{
					const HE_Vertex* v = (*e)->dst();
					switch (idx)
					{
					case 0:
						mT_StMesh_rend(0, vi) = 1;
						mT_StMesh_rend(1, vi) = 0;
						break;
					case 1:
						mT_StMesh_rend(0, vi) = 1;
						mT_StMesh_rend(1, vi) = 1;
						break;
					case 2:
						mT_StMesh_rend(0, vi) = 0.5;
						mT_StMesh_rend(1, vi) = 1;
						break;
					case 3:
						mT_StMesh_rend(0, vi) = 0;
						mT_StMesh_rend(1, vi) = 1;
						break;
					case 4:
						mT_StMesh_rend(0, vi) = 0;
						mT_StMesh_rend(1, vi) = 0;
						break;
					default: break;
					}
					mV_StMesh_rend(0, vi) = v->position().x;
					mV_StMesh_rend(1, vi) = v->position().y;
					mV_StMesh_rend(2, vi) = v->position().z;
					mC_StMesh_rend(0, vi) = c.x;
					mC_StMesh_rend(1, vi) = c.y;
					mC_StMesh_rend(2, vi++) = c.z;
					++e; ++idx;
				} while (e != sentinel);
			}
			else if (f->TestFlag(VEF_FLAG_QUAD_SHORTROW_L) ||
				f->TestFlag(VEF_FLAG_QUAD_SHORTROW_R))
			{
				do
				{
					const HE_Vertex* v = (*e)->dst();
					mT_StMesh_rend(0, vi) = 0;
					mT_StMesh_rend(1, vi) = 0;
					mV_StMesh_rend(0, vi) = v->position().x;
					mV_StMesh_rend(1, vi) = v->position().y;
					mV_StMesh_rend(2, vi) = v->position().z;
					mC_StMesh_rend(0, vi) = c.x;
					mC_StMesh_rend(1, vi) = c.y;
					mC_StMesh_rend(2, vi++) = c.z;
					++e;
				} while (e != sentinel);
			}
			else {
				std::cout << "ERROR!\n";
			}
		}
	}

	int triTotalNum = triNum + 2 * quadNum + 3 * penNum;

	mF_StMesh_rend.resize(3, triTotalNum);

	int c = 0;
	vi = 0;
	for (int j = 0; j < (int)mCleanDual->_groupFaceIdx.size(); j++)
	{
		for (int i = 0; i < (int)mCleanDual->_groupFaceIdx[j].size(); i++)
		{
			const HE_Face* f = mCleanPoly->face(mCleanDual->_groupFaceIdx[j][i]);

			if (f->hole()) continue;

			std::vector<int> viList; viList.clear();
			HE_Face::const_edge_circulator e = f->begin();
			HE_Face::const_edge_circulator sentinel = e;
			do
			{
				viList.push_back(vi);
				vi++;
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
	}
#endif
	mE_StMesh_rend.resize(6, mCleanPoly->numHalfEdges() * 2);

	for (int ei = 0; ei < mCleanPoly->numHalfEdges(); ei++)
	{
		HE_HalfEdge* he = mCleanPoly->halfedge(ei);
		mE_StMesh_rend(0, ei * 2 + 0) = he->src()->position().x;
		mE_StMesh_rend(1, ei * 2 + 0) = he->src()->position().y;
		mE_StMesh_rend(2, ei * 2 + 0) = he->src()->position().z;
		mE_StMesh_rend(3, ei * 2 + 0) = 0;
		mE_StMesh_rend(4, ei * 2 + 0) = 0;
		mE_StMesh_rend(5, ei * 2 + 0) = 0;

		mE_StMesh_rend(0, ei * 2 + 1) = he->dst()->position().x;
		mE_StMesh_rend(1, ei * 2 + 1) = he->dst()->position().y;
		mE_StMesh_rend(2, ei * 2 + 1) = he->dst()->position().z;
		mE_StMesh_rend(3, ei * 2 + 1) = 0;
		mE_StMesh_rend(4, ei * 2 + 1) = 0;
		mE_StMesh_rend(5, ei * 2 + 1) = 0;
	}

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

void MultiResolutionHierarchy::exportResult(char * path)
{
	std::fstream exportFile(path, std::ios::out);

	for (int vi = 0; vi < mCleanPoly->numVertices(); vi++)
	{
		exportFile << "v " << mCleanPoly->vertex(vi)->position().x << " " << mCleanPoly->vertex(vi)->position().y << " " << mCleanPoly->vertex(vi)->position().z << std::endl;
	}

	for (int fi = 0; fi < mCleanPoly->numFaces(); fi++)
	{
		HE_Face* f = mCleanPoly->face(fi);
		if (f->hole()) continue;

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