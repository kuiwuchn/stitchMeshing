#pragma once
#include <cstdlib>
#include <vector>
#include <Eigen/Dense>
//using namespace Eigen;
using namespace std;

/*typedefs*/
#if defined(SINGLE_PRECISION)
typedef float Float;
#else
typedef double Float;
#endif

#define Interior_RegularE 4
#define Boundary_RegularE 2

#define Precision 1.e-7
#define Precision_Pro 1.0e-5
#define Jacobian_Bound 1.0e-3

#define PAI 3.1415926


enum Feature_V_Type {
	INTERIOR = -4,
	CORNER,
	LINE,
	REGULAR
};

const int tetra_table[4][4] =
{
	{ 0,2,1,3 },
	{ 1,0,2,3 },
	{ 2,1,0,3 },
	{ 3,0,1,2 },
};

const int hex_face_table[6][4] =
{
	{ 0,1,2,3 },
	{ 4,5,6,7 },
	{ 0,1,5,4 },
	{ 0,4,7,3 },
	{ 3,2,6,7 },
	{ 1,5,6,2 },
};
const int hex_tetra_table[8][4] =
{
	{ 0,3,4,1 },
	{ 1,0,5,2 },
	{ 2,1,6,3 },
	{ 3,2,7,0 },
	{ 4,7,5,0 },
	{ 5,4,6,1 },
	{ 6,5,7,2 },
	{ 7,6,4,3 },
};
const double hex_ref_shape[8][3] =
{
	{ 0, 0, 0 },
	{ 1, 0, 0 },
	{ 1, 1, 0 },
	{ 0, 1, 0 },
	{ 0, 0, 1 },
	{ 1, 0, 1 },
	{ 1, 1, 1 },
	{ 0, 1, 1 },
};

//-------------------------------------------------------------------
//---For Hybrid mesh-------------------------------------------------
struct Hybrid_V
{
	uint32_t id, svid, fvid;
	vector<Float> v;
	vector<uint32_t> neighbor_vs;
	vector<uint32_t> neighbor_es;
	vector<uint32_t> neighbor_fs;
	vector<uint32_t> neighbor_hs;

	bool boundary;
};
struct Hybrid_E
{
	uint32_t id;
	vector<uint32_t> vs;
	vector<uint32_t> neighbor_fs;
	vector<uint32_t> neighbor_hs;
	
	bool boundary;

	bool hex_edge = false;
};
struct Hybrid_F
{
	uint32_t id;
	vector<uint32_t> vs;
	vector<uint32_t> es;
	vector<uint32_t> neighbor_hs;
	bool boundary;
};

struct Hybrid
{
	uint32_t id;
	vector<uint32_t> vs;
	vector<uint32_t> es;
	vector<uint32_t> fs;
	bool boundary;
	bool hex = false;
};
//-------------------------------------------------------------------
//-------------------------------------------------------------------
struct Singular_V
{
	uint32_t id, hid;
	bool boundary;

	vector<uint32_t> neighbor_svs;//singular vs
	vector<uint32_t> neighbor_ses;//singular es
	bool fake;

	uint32_t which_singularity;
	uint32_t which_singularity_type;
};
struct Singular_E
{
	uint32_t id;
	std::vector<uint32_t> vs;//singular v
	vector<uint32_t> es_link;//hex e
	vector<uint32_t> vs_link;//hex v
	bool boundary;

	vector<uint32_t> neighbor_ses;//singular es
	bool circle;//0 for normal, 2 for circle, 
};
struct Frame_V
{
	uint32_t id, hid, svid;
	uint32_t what_type;//1 singular, 2, extra-nordinary, 3, regular
				  //0 non extra-ordinary node and at surface,1 extra-ordinary node and at surface;
				  //2 non extra-ordinary node and inside volume,3 extra-ordinary node and inside volume;
	vector<uint32_t>  neighbor_fvs;
	vector<uint32_t>  neighbor_fes;
	vector<uint32_t>  neighbor_ffs;
	vector<uint32_t>  neighbor_fhs;
	bool boundary;
};
struct Frame_E
{
	uint32_t id;
	bool singular = false;
	std::vector<uint32_t> vs;
	bool boundary;//0 interior, 1 boundary, 2 first and last slice boundary
	vector<uint32_t>  vs_link;//v of hex_v
	vector<uint32_t>  es_link;//e of hex_e
	vector<uint32_t>  neighbor_fes;
	vector<uint32_t>  neighbor_ffs;
	vector<uint32_t>  neighbor_fhs;
};
struct Frame_F
{
	uint32_t id;
	bool boundary;
	uint32_t F_location;//for rendering

	vector<uint32_t> vs;
	vector<uint32_t> es;
	vector<uint32_t>  fvs_net;
	vector<uint32_t>  ffs_net;
	vector<uint32_t>  neighbor_ffs;
	vector<uint32_t>  neighbor_fhs;
	uint32_t Color_ID;
};
struct Frame_H
{
	uint32_t id;

	std::vector<uint32_t> vs;
	std::vector<uint32_t> es;
	std::vector<uint32_t> fs;
	vector<vector<vector<uint32_t> >> vs_net;
	vector<uint32_t>  fs_net;
	vector<uint32_t>  hs_net;
	vector<uint32_t>  neighbor_fhs;//neighboring cube	
	uint32_t Color_ID;
};

struct Singularity
{
	vector<Singular_V> SVs;
	vector<Singular_E> SEs;
};
struct Frame
{
	vector<Frame_V> FVs;
	vector<Frame_E> FEs;
	vector<Frame_F> FFs;
	vector<Frame_H> FHs;
};
enum Mesh_type {
	Tri = 0,
	Qua,
	Tet,
	Hyb,
	Hex
};
struct Mesh_Topology
{
	bool euler_problem;
	bool manifoldness_problem;

	int genus;
	int surface_euler;
	int volume_euler;
	bool surface_manifoldness;
	bool volume_manifoldness;

	bool frame_euler_problem;
	bool frame_manifoldness_problem;

	int frame_genus;
	int frame_surface_euler;
	int frame_volume_euler;
	bool frame_surface_manifoldness;
	bool frame_volume_manifoldness;
};

struct Mesh_Quality
{
	double min_Jacobian;
	double ave_Jacobian;
	Eigen::VectorXd V_Js;
	Eigen::VectorXd H_Js;
	Eigen::VectorXi Self_In;//self intersected 1, othersize 0
	Eigen::VectorXd H_Vols;//polyhedral volume
};
struct Mesh
{
	short type;//Mesh_type
	Eigen::MatrixXd V;
	vector<Hybrid_V> Vs;
	vector<Hybrid_E> Es;
	vector<Hybrid_F> Fs;
	vector<Hybrid> Hs;
};