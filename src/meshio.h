#pragma once
#include "common.h"
#include "global_types.h"

extern void loadTetMesh(const std::string &prefix, MatrixXf &V, MatrixXu &F, MatrixXu &T);
extern void loadTriMesh(const std::string &filename, MatrixXf &V, MatrixXu &F);
extern void load_obj(const std::string &filename, MatrixXu &F, MatrixXf &V);
extern void load_off(const std::string &filename, std::vector<std::vector<uint32_t>> &F, MatrixXf &V);
extern void loadTetMesh_VTK(const std::string &prefix, MatrixXf &V, MatrixXu &F, MatrixXu &T);
extern void loadTetMesh_MESH(const std::string &prefix, MatrixXf &V, MatrixXu &F, MatrixXu &T);

extern void write_obj(const std::string &filename, const MatrixXu &F,
	const MatrixXf &V,
	const MatrixXf &N = MatrixXf(),
	const MatrixXf &Nf = MatrixXf(),
	const MatrixXf &UV = MatrixXf(),
	const MatrixXf &C = MatrixXf());
extern void write_surface_mesh_VTK(MatrixXf &V, std::vector<std::vector<uint32_t>> &F, char * path);
extern void write_surface_mesh_OFF(MatrixXf &V, std::vector<std::vector<uint32_t>> &F, char * path);
extern void write_surface_mesh_OBJ(MatrixXf &V, std::vector<std::vector<uint32_t>> &F, char * path);
extern void write_surface_mesh_OBJ(MatrixXf &V, MatrixXu &F, char * path);

extern void write_volume_mesh_VTK(MatrixXf &V, std::vector<std::vector<uint32_t>> &T, char * path);
extern void write_volume_mesh_VTK(MatrixXf &V, std::vector<tuple_E> &E, std::vector<std::vector<uint32_t>> &F, std::vector<int> &F_type, char * path);
extern void write_volume_mesh_HYBRID(MatrixXf &V, std::vector<std::vector<uint32_t>> &F, std::vector<std::vector<uint32_t>> &P, std::vector<bool> &P_flag, std::vector<std::vector<bool>> &PF_flag, char * path);
extern void write_volume_mesh_MESH(MatrixXf &V, std::vector<std::vector<uint32_t>> &T, char * path);

extern void write_edge_coloring_TXT(std::vector<MatrixXf> &E_Colors, char * path);

extern void write_singularities_SING(MatrixXf &Singularity, char * path);
extern void write_Vertex_Types_TXT(std::vector<int> &V_types, char * path);

extern void write_statistics_TXT(statistics &sta, char * path);


extern void load_HYBRID_mesh(Mesh &mesh, string path);

