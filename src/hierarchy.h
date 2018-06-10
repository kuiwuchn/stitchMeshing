#pragma once

#include "common.h"
#include "aabb.h"
#include "lock.h"
#include "adjacency.h"
#include "meshio.h"
#include <nanogui/serializer/core.h>
#include <nanogui/serializer/sparse.h>
#include <unordered_map>
#include <algorithm> 
#include <queue>
#include "global_types.h"
#include <set>

#include "core/HE_Polyhedron.h"

using nanogui::Serializer;
using namespace std;
class BVH;

struct MeshStats {
	AABB mAABB;
	Vector3f mWeightedCenter;
	double mAverageEdgeLength;
	double mMaximumEdgeLength;
	double mSurfaceArea;

	MeshStats() :
		mWeightedCenter(Vector3f::Zero()),
		mAverageEdgeLength(0.0f),
		mMaximumEdgeLength(0.0f),
		mSurfaceArea(0.0f) { }
};
class MultiResolutionHierarchy {
public:
    MultiResolutionHierarchy();

    bool load(const std::string &filename);
	MeshStats compute_mesh_stats(const MatrixXu &F_, const MatrixXf &V_, bool deterministic = false);

	//protected:
	void build();
	void construct_tEs_tFEs(MatrixXu & F, std::vector<std::vector<uint32_t>> &mtFes, std::vector<tuple_E> &mtEs);
	void construct_tEs_tFEs(std::vector<std::vector<uint32_t>> &F, std::vector<std::vector<uint32_t>> &mtFes, std::vector<tuple_E> &mtEs);
	void orient_polygon_mesh(MatrixXf &HV, vector<vector<uint32_t>> &HF, vector<vector<uint32_t>> &HFE, vector<tuple_E> &Es);
	void orient_polygon_mesh(MatrixXf &HV, vector<vector<uint32_t>> &HF);
	void laplacian_smoothing(MatrixXf &V_, int smooth_iterations);

    bool tetMesh() const { return mT.cols() > 0; }
    bool triMesh() const { return mT.cols() == 0; }

    MatrixXf &V(uint32_t i = 0) { return mV[i]; }
    const MatrixXf &V(uint32_t i = 0) const { return mV[i]; }

    MatrixXf &N(uint32_t i = 0) { return mN[i]; }
    const MatrixXf &N(uint32_t i = 0) const { return mN[i]; }

    MatrixXf &Q(uint32_t i = 0) { return mQ[i]; }
    const MatrixXf &Q(uint32_t i = 0) const { return mQ[i]; }

    MatrixXf &O(uint32_t i = 0) { return mO[i]; }
    const MatrixXf &O(uint32_t i = 0) const { return mO[i]; }

    MatrixXf &C(uint32_t i = 0) { return mC[i]; }
    const MatrixXf &C(uint32_t i = 0) const { return mC[i]; }

    MatrixXu &F() { return mF; }
    const MatrixXu &F() const { return mF; }

    MatrixXu &T() { return mT; }
    const MatrixXu &T() const { return mT; }

    SMatrix &L(uint32_t i = 0) { return mL[i]; }
    const SMatrix &L(uint32_t i) const { return mL[i]; }

    const BVH *bvh() const { return mBVH; }

    void smoothOrientationsTet(uint32_t l, bool alignment, bool randomization);
    void smoothOrientationsTri(uint32_t l, bool alignment, bool randomization, bool extrinsic);
    void detectOrientationSingularitiesTri();
    void detectOrientationSingularitiesTet();
    void prolongOrientations(int level);

    void smoothPositionsTet(uint32_t l, bool alignment, bool randomization);
    void smoothPositionsTri(uint32_t l, bool alignment, bool randomization, bool extrinsic);
    
    void detectPositionSingularitiesTri();
    void detectPositionSingularitiesTet();
    void prolongPositions(int level);

//init edge tagging
	void init_edge_tagging2D();
	void init_edge_tagging3D();
//mesh extraction
		bool tagging_collapseTri(bool triangle_Switch);

	bool meshExtraction2D();
	void edge_tagging2D(vector<uint32_t> &ledges);
	void swap_data2D();
	bool split_long_edge2D(vector<uint32_t> &ledges);
	bool split_face2D(bool red_edge);
	bool split_pentagon();
	bool remove_doublets2D();
	Float compute_cost_edge2D_angle(int32_t v0, int32_t v1, vector<uint32_t> &vs, MatrixXf &V_);


	bool meshExtraction3D();
	void construct_Es_Fs_Polyhedral();
	void orient_hybrid_mesh(MatrixXf &HV, vector<vector<uint32_t>> &HF, vector<vector<uint32_t>> &HP, vector<vector<bool>> &HPF_flag);
	void swap_data3D();

		bool edge_tagging3D(vector<uint32_t> &ledges);
		void tagging_collapseTet();
		bool split_long_edge3D(vector<uint32_t> &ledges);
		bool split_face3D(bool red_edge);
		bool split_polyhedral3D();
		void candidate_loops(vector<uint32_t> &fs, vector<vector<uint32_t>> &nfes, vector<vector<uint32_t>> &nfvs, vector<pair<double, uint32_t>> &fs_rank);
		bool recursive_ring(vector<uint32_t> &rvs, vector<uint32_t> &pvs, vector<uint32_t> &pes, uint32_t v0, uint32_t ve[2]);
		int32_t share_the_same_fs(vector<uint32_t> &es, const uint32_t eid);
		Float compute_cost_face3D(vector<uint32_t> &vs, uint32_t rv, bool single);
		Float compute_cost_edge3D(uint32_t v0, uint32_t v1);
		void tagging_singularities_T_nodes(MatrixXf &V_tagging, vector<tuple_E> &E_tagging, vector<vector<uint32_t>> &F_tagging);

	void composit_edges_colors(MatrixXf &Result_Vs, std::vector<tuple_E> &Es_to_render, MatrixXf &Result_edges);
	void composit_edges_centernodes_triangles(std::vector<std::vector<uint32_t>> &Actual_Fs, MatrixXf &nodes, MatrixXf &Result_edges, MatrixXf &center_nodes, MatrixXu &Triangles);
	
    const MatrixXf &orientationSingularities() const { return mOrientationSingularities; }
    const MatrixXf &positionSingularities() const { return mPositionSingularities; }

    AABB aabb() const { return mAABB; }

    size_t vertexCount() const { return mV.size() > 0 ? mV[0].cols() : 0; }
    size_t faceCount() const { return mF.cols(); }
    size_t tetCount() const { return mT.cols(); }

	void set_tet_elen_ratio(Float ratio) { tElen_ratio = ratio; tet_elen = ratio * ms.mAverageEdgeLength; };
	Float tet_elen_ratio() { return tElen_ratio; };
    Float averageEdgeLength() const { return mAverageEdgeLength; }
	Float scale() const { return ratio_scale; }
	void setScale(Float scale) { 
		ratio_scale = scale; 
		mScale = diagonalLen * scale; 
		mInvScale = 1.f / mScale;
		tet_elen = tElen_ratio * ratio_scale * diagonalLen * 0.3;
	}

    ordered_lock &mutex() const { return mMutex; }

    int levels() const { return mL.size(); }

// stitch meshing
	void convert2Poly();

public:
	//for both 2D & 3D 
	int laplacian_iteration = 0;
    std::vector<MatrixXf> mV;
    std::vector<MatrixXf> mN;
    std::vector<MatrixXf> mQ;
    std::vector<MatrixXf> mO;
    std::vector<MatrixXf> mC;
    std::vector<SMatrix> mL;
    std::vector<SMatrix> mP;
    MatrixXu mF;
    MatrixXu mT;
	bool Q_FROM_FILE = false;
	MatrixXf QoF;

	std::vector<std::vector<uint32_t>> nFes;
	std::vector<tuple_E> nEs;
	vector<vector<bool>> nV_boundary_flag;
	std::vector<std::vector<uint32_t>> nV_nes;


	vector<vector<uint32_t>> vnfs;
	Float quadricW = 1;

    MatrixXf mOrientationSingularities;
    MatrixXf mPositionSingularities;
    MatrixXf mNF, mCF;
    uint32_t mOrientationIterations;
    uint32_t mPositionIterations;
    AABB mAABB;
    Float mAverageEdgeLength;
    mutable ordered_lock mMutex;
    Float mScale, mInvScale;
	Float diagonalLen;
	Float ratio_scale;
	Float tet_elen, tElen_ratio;

    BVH *mBVH;

    MatrixXf mQ_combed;
    MatrixXi mO_combed;

public:
	std::string outpath;
	statistics sta;
	MeshStats ms;
	bool Two_rosy_flag=true;
	bool re_color;
	bool doublets, triangles, decomposes, splitting, Qquadric;
	bool o_flag;
	MatrixXf mV_final;
	std::vector<tuple_E> Es_reddash_left;

	MatrixXf mV_tag;
	std::vector<std::vector<uint32_t>> FEs_tag, F_tag, P_tag;
	std::vector<bool> Hex_flag;
	std::vector<std::vector<bool>> PF_flag;

	std::vector<int> F_tag_type;
	std::vector<std::vector<uint32_t>> F_final, P_final;
	std::vector<Vector4f> ECs;
	//for t-mesh vertex tag
	vector<int> V_flag;
	//for rendering
	MatrixXf mO_center, mV_tag_rend;
	MatrixXf mV_final_rend;
	MatrixXf E_rend, E_O_rend, E_I_rend;
	MatrixXf E_rend_o, E_O_rend_o, E_I_rend_o;
	MatrixXf E_tag_rend, E_tag_left_rend;
	MatrixXf E_final_rend;
	MatrixXu F_tag_rend;
	MatrixXu F_final_rend;

	//////////////////////////////////////////////////////////////////////////
	// stitch meshing 

	HE_Polyhedron* mPoly;
};
