#pragma once
#include <queue>
#include<algorithm>
#include<vector>
#include <iterator>
#include<set>
#include "global_types.h"

class base_complex
{
public:
	base_complex() {};
	void singularity_structure(Singularity &si, Mesh &mesh);
	void base_complex_extraction(Singularity &si, Frame &frame, Mesh &mesh);
	void base_complex_node_edge_extraction(Singularity &si, Frame &frame, Mesh &mesh);
	void node_on_circular_singularity(Singularity &si, Frame &frame, Mesh &mesh, vector<uint32_t> &Nodes);
	void base_complex_face_extraction(Singularity &si, Frame &frame, Mesh &mesh);
	void base_complex_cuboid_extraction(Singularity &si, Frame &frame, Mesh &mesh);

	void singularity_base_complex(Singularity &si, Frame &frame, Mesh &mesh);

	void assign_color(Frame &frame);
	~base_complex() {};
};

