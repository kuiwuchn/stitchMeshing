#include "meshio.h"
#include <fstream>
#include <rply.h>
#include <unordered_map>
#include "timer.h"
#include <iomanip>

template <typename ParseHeader, typename ParseLine>
void loadTextFile(const std::string &filename, ParseHeader parseHeader,
	ParseLine parseLine) {
	std::ifstream is(filename);
	if (is.fail())
		throw std::runtime_error("Unable to open tetgen file \"" + filename + "\"!");

	std::string line_str;
	int line = 0, actual_line = 0;

	try {
		while (std::getline(is, line_str)) {
			++actual_line;
			std::istringstream iss(line_str);

			if (line_str.empty() || line_str[0] == '#')
				continue;

			if (++line == 1)
				parseHeader(iss);
			else
				parseLine(iss);
		}
	}
	catch (const std::exception &e) {
		throw std::runtime_error(
			"Unable to parse tetgen file \"" + filename + ".node\" (line "
			+ std::to_string(actual_line) + "): " + e.what());
	}
}

void loadTetMesh(const std::string &prefix, MatrixXf &V, MatrixXu &F, MatrixXu &T) {
	uint32_t vertexCount, faceCount, tetCount, dimension, nodesPerTet, boundaryInfo;

	Timer<> timer;
	timer.beginStage("Loading tetrahedral mesh \"" + prefix + ".{node/face/ele}\"");

	loadTextFile(
		prefix + ".node",
		[&](std::istringstream &is) {
		if (!(is >> vertexCount) || !(is >> dimension))
			throw std::runtime_error("Unable to parse header!");
		if (dimension != 3)
			throw std::runtime_error("Invalid dimension!");
		V.resize(3, vertexCount);
	},
		[&](std::istringstream &is) {
		uint32_t index;
		Float x, y, z;
		if (!(is >> index >> x >> y >> z))
			throw std::runtime_error("Unable to parse vertex data!");
		V.col(index) << x, y, z;
	}
	);

	loadTextFile(
		prefix + ".ele",
		[&](std::istringstream &is) {
		if (!(is >> tetCount) || !(is >> nodesPerTet))
			throw std::runtime_error("Unable to parse header!");
		if (nodesPerTet != 4)
			throw std::runtime_error("Only files with 4 nodes per tetrahedron are supported!");
		T.resize(4, tetCount);
	},
		[&](std::istringstream &is) {
		uint32_t index, i0, i1, i2, i3;
		if (!(is >> index >> i0 >> i1 >> i2 >> i3))
			throw std::runtime_error("Unable to parse element data!");
		T.col(index) << i0, i1, i2, i3;
	}
	);

	uint32_t actualFaceCount = 0;
	loadTextFile(
		prefix + ".face",
		[&](std::istringstream &is) {
		if (!(is >> faceCount) || !(is >> boundaryInfo))
			throw std::runtime_error("Unable to parse header!");
		F.resize(3, faceCount);
	},
		[&](std::istringstream &is) {
		uint32_t index, i0, i1, i2;
		if (!(is >> index >> i0 >> i1 >> i2))
			throw std::runtime_error("Unable to parse face data!");
		if (boundaryInfo) {
			int b;
			if (!(is >> b))
				throw std::runtime_error("Unable to parse face data!");
			if (b == 0)
				return;
		}
		F.col(actualFaceCount++) << i0, i1, i2;
	}
	);
	F.row(0).swap(F.row(1));
	F.conservativeResize(3, actualFaceCount);

	timer.endStage("V=" + std::to_string(V.cols()) + ", F = " +
		std::to_string(F.cols()) + ", T = " +
		std::to_string(T.cols()));
}

void loadTriMesh(const std::string &filename, MatrixXf &V, MatrixXu &F) {
	auto message_cb = [](p_ply ply, const char *msg) { cerr << "rply: " << msg << endl; };

	Timer<> timer;
	timer.beginStage("Loading triangle mesh \"" + filename + "\"");

	p_ply ply = ply_open(filename.c_str(), message_cb, 0, nullptr);
	if (!ply)
		throw std::runtime_error("Unable to open PLY file \"" + filename + "\"!");

	if (!ply_read_header(ply)) {
		ply_close(ply);
		throw std::runtime_error("Unable to open PLY header of \"" + filename + "\"!");
	}

	p_ply_element element = nullptr;
	uint32_t vertexCount = 0, faceCount = 0;

	/* Inspect the structure of the PLY file, load number of faces if avaliable */
	while ((element = ply_get_next_element(ply, element)) != nullptr) {
		const char *name;
		long nInstances;

		ply_get_element_info(element, &name, &nInstances);
		if (!strcmp(name, "vertex"))
			vertexCount = (uint32_t)nInstances;
		else if (!strcmp(name, "face"))
			faceCount = (uint32_t)nInstances;
	}

	if (vertexCount == 0 || faceCount == 0)
		throw std::runtime_error("PLY file \"" + filename + "\" is invalid! No face/vertex/elements found!");

	F.resize(3, faceCount);
	V.resize(3, vertexCount);

	struct VertexCallbackData { MatrixXf &V; };
	struct FaceCallbackData { MatrixXu &F; };

	auto rply_vertex_cb = [](p_ply_argument argument) -> int {
		VertexCallbackData *data; long index, coord;
		ply_get_argument_user_data(argument, (void **)&data, &coord);
		ply_get_argument_element(argument, nullptr, &index);
		data->V(coord, index) = (Float)ply_get_argument_value(argument);
		return 1;
	};

	auto rply_index_cb = [](p_ply_argument argument) -> int {
		FaceCallbackData *data;
		long length, value_index, index;
		ply_get_argument_property(argument, nullptr, &length, &value_index);

		if (length != 3)
			throw std::runtime_error("Only triangle faces are supported!");

		ply_get_argument_user_data(argument, (void **)&data, nullptr);
		ply_get_argument_element(argument, nullptr, &index);

		if (value_index >= 0)
			data->F(value_index, index) = (uint32_t)ply_get_argument_value(argument);

		return 1;
	};

	VertexCallbackData vcbData{ V };
	FaceCallbackData fcbData{ F };

	if (!ply_set_read_cb(ply, "vertex", "x", rply_vertex_cb, &vcbData, 0) ||
		!ply_set_read_cb(ply, "vertex", "y", rply_vertex_cb, &vcbData, 1) ||
		!ply_set_read_cb(ply, "vertex", "z", rply_vertex_cb, &vcbData, 2)) {
		ply_close(ply);
		throw std::runtime_error("PLY file \"" + filename + "\" does not contain vertex position data!");
	}

	if (!ply_set_read_cb(ply, "face", "vertex_indices", rply_index_cb, &fcbData, 0)) {
		ply_close(ply);
		throw std::runtime_error("PLY file \"" + filename + "\" does not contain vertex indices!");
	}

	if (!ply_read(ply)) {
		ply_close(ply);
		throw std::runtime_error("Error while loading PLY data from \"" + filename + "\"!");
	}

	ply_close(ply);
	timer.endStage("V=" + std::to_string(V.cols()) + ", F = " + std::to_string(F.cols()));
}

void load_obj(const std::string &filename, MatrixXu &F, MatrixXf &V) {
	/// Vertex indices used by the OBJ format
	struct obj_vertex {
		uint32_t p = (uint32_t)-1;
		uint32_t n = (uint32_t)-1;
		uint32_t uv = (uint32_t)-1;

		inline obj_vertex() { }

		inline obj_vertex(const std::string &string) {
			std::vector<std::string> tokens = str_tokenize(string, '/', true);

			if (tokens.size() < 1 || tokens.size() > 3)
				throw std::runtime_error("Invalid vertex data: \"" + string + "\"");

			p = str_to_uint32_t(tokens[0]);

#if 0
			if (tokens.size() >= 2 && !tokens[1].empty())
				uv = str_to_uint32_t(tokens[1]);

			if (tokens.size() >= 3 && !tokens[2].empty())
				n = str_to_uint32_t(tokens[2]);
#endif
		}

		inline bool operator==(const obj_vertex &v) const {
			return v.p == p && v.n == n && v.uv == uv;
		}
	};

	/// Hash function for obj_vertex
	struct obj_vertexHash : std::unary_function<obj_vertex, size_t> {
		std::size_t operator()(const obj_vertex &v) const {
			size_t hash = std::hash<uint32_t>()(v.p);
			hash = hash * 37 + std::hash<uint32_t>()(v.uv);
			hash = hash * 37 + std::hash<uint32_t>()(v.n);
			return hash;
		}
	};

	typedef std::unordered_map<obj_vertex, uint32_t, obj_vertexHash> VertexMap;

	size_t last_slash_idx = filename.rfind('.');
	if (filename.substr(last_slash_idx) != ".OBJ" && filename.substr(last_slash_idx) != ".obj")
		throw std::runtime_error("Unable to open OBJ file \"" + filename + "\"!");

	std::ifstream is(filename);
	if (is.fail())
		throw std::runtime_error("Unable to open OBJ file \"" + filename + "\"!");
	cout << "Loading \"" << filename << "\" .. ";
	cout.flush();
	Timer<> timer;

	std::vector<Vector3f>   positions;
	//std::vector<Vector2f>   texcoords;
	//std::vector<Vector3f>   normals;
	std::vector<uint32_t>   indices;
	std::vector<obj_vertex> vertices;
	VertexMap vertexMap;

	std::string line_str;
	while (std::getline(is, line_str)) {
		std::istringstream line(line_str);

		std::string prefix;
		line >> prefix;

		if (prefix == "v") {
			Vector3f p;
			line >> p.x() >> p.y() >> p.z();
			positions.push_back(p);
		}
		else if (prefix == "vt") {
			/*
			Vector2f tc;
			line >> tc.x() >> tc.y();
			texcoords.push_back(tc);
			*/
		}
		else if (prefix == "vn") {
			/*
			Vector3f n;
			line >> n.x() >> n.y() >> n.z();
			normals.push_back(n);
			*/
		}
		else if (prefix == "f") {
			std::string v1, v2, v3, v4;
			line >> v1 >> v2 >> v3 >> v4;
			obj_vertex tri[6];
			int nVertices = 3;

			tri[0] = obj_vertex(v1);
			tri[1] = obj_vertex(v2);
			tri[2] = obj_vertex(v3);

			if (!v4.empty()) {
				/* This is a quad, split into two triangles */
				tri[3] = obj_vertex(v4);
				tri[4] = tri[0];
				tri[5] = tri[2];
				nVertices = 6;
			}
			/* Convert to an indexed vertex list */
			for (int i = 0; i<nVertices; ++i) {
				const obj_vertex &v = tri[i];
				VertexMap::const_iterator it = vertexMap.find(v);
				if (it == vertexMap.end()) {
					vertexMap[v] = (uint32_t)vertices.size();
					indices.push_back((uint32_t)vertices.size());
					vertices.push_back(v);
				}
				else {
					indices.push_back(it->second);
				}
			}
		}
	}
	F.resize(3, indices.size() / 3);
	memcpy(F.data(), indices.data(), sizeof(uint32_t)*indices.size());
	V.resize(3, vertices.size());
	for (uint32_t i = 0; i<vertices.size(); ++i)
		V.col(i) = positions.at(vertices[i].p - 1);
}

void load_off(const std::string &filename, std::vector < std::vector<uint32_t>> &F, MatrixXf &V) {
	FILE *ff = fopen(filename.data(), "rt");
	char s[1024], sread[1024], sread2[1024];
	int vnum, hnum, r;	float x, y, z;
	if (fscanf(ff, "%s", &sread) != 1 || (strcmp(sread, "OFF") != 0))
		throw std::runtime_error("cannot find head of OFF!");

	fscanf(ff, "%d %d %d", &vnum, &hnum, &r);
	V.resize(3, vnum);
	V.setZero();
	for (int i = 0; i<vnum; i++)
	{
		int framid;
		fscanf(ff, "%f %f %f", &x, &y, &z);

		V(0, i) = x;
		V(1, i) = y;
		V(2, i) = z;
	}
	F.clear();
	F.resize(hnum);
	for (int i = 0; i<hnum; i++) {
		int nw;
		fscanf(ff, "%d", &nw);
		F[i].resize(nw);
		for (int j = 0; j<nw; j++) {
			fscanf(ff, "%d", &(F[i][j]));
			F[i][j]--;
		}
	}

	fclose(ff);
}
void loadTetMesh_VTK(const std::string &filename, MatrixXf &V, MatrixXu &F, MatrixXu &T)
{
	char file[300];
	//std::ifstream ff(filename+".vtk", std::ios::in);
	std::ifstream ff(filename, std::ios::in);
	char s[1024], sread[1024], sread2[1024];
	int vnum, hnum;	float x, y, z; int is_boundary;

	bool find = false; int lines = 0;
	while (!find)
	{
		ff.getline(s, 1023);
		if (sscanf(s, "%s %d %s", &sread, &vnum, &sread2) == 3 && (strcmp(sread, "POINTS") == 0))
			find = true;
		if (++lines>10)
			throw std::runtime_error("cannot find head of VTK!");
	}
	V.resize(3, vnum);
	V.setZero();
	for (int i = 0; i<vnum; i++)
	{
		ff.getline(s, 1023);
		int framid;
		sscanf(s, "%f %f %f", &x, &y, &z);

		V(0, i) = x;
		V(1, i) = y;
		V(2, i) = z;
	}

	find = false;
	while (!find)
	{
		int temp_int;
		ff.getline(s, 1023);
		if (sscanf(s, "%s %d %d", &sread, &hnum, &temp_int) == 3 && (strcmp(sread, "CELLS") == 0))
			find = true;
	}
	std::vector<std::vector<int>> T_, F_;
	T_.reserve(hnum); F_.reserve(hnum);
	int a, b, c, d;
	while (ff.getline(s, 1023))
	{
		std::vector<int> vs;
		if (sscanf(s, "%d %d %d %d", &vnum, &a, &b, &c) == 4 && vnum == 3)
		{
			vs.push_back(a);
			vs.push_back(b);
			vs.push_back(c);
			F_.push_back(vs);
		}
		else if (sscanf(s, "%d %d %d %d %d", &vnum, &a, &b, &c, &d) == 5 && vnum == 4)
		{
			vs.push_back(a);
			vs.push_back(b);
			vs.push_back(c);
			vs.push_back(d);
			T_.push_back(vs);
		}
	}
	F.resize(3, F_.size()); T.resize(4, T_.size());
	for (int i = 0; i<F_.size(); i++)
	{
		F(0, i) = F_[i][0];
		F(1, i) = F_[i][1];
		F(2, i) = F_[i][2];
	}
	for (int i = 0; i<T_.size(); i++)
	{
		T(0, i) = T_[i][0];
		T(1, i) = T_[i][1];
		T(2, i) = T_[i][2];
		T(3, i) = T_[i][3];
	}
	ff.close();
}
void loadTetMesh_MESH(const std::string &prefix, MatrixXf &V, MatrixXu &F, MatrixXu &T)
{
	//Vs.clear(); Hexs.clear();

	//std::fstream f(fname, std::ios::in);
	//char s[1024], sread[1024];
	//int vnum, tnum;	float x, y, z;

	//int find = false;
	//while (!find)
	//{
	//	f.getline(s, 1023);
	//	if (sscanf(s, "%s", &sread) == 1 && (strcmp(sread, "Vertices") == 0))
	//		find = true;
	//}
	//f.getline(s, 1023);
	//sscanf(s, "%d", &vnum);
	//for (int i = 0; i<vnum; i++)
	//{
	//	f.getline(s, 1023);
	//	int framid;
	//	sscanf(s, "%f %f %f %d", &x, &y, &z, &framid);

	//	Hex_V v;
	//	v.index = i;
	//	v.where_location = -1;
	//	v.fixed = false;
	//	v.slice_id = framid;
	//	v.v[0] = x; v.v[1] = y; v.v[2] = z;
	//	Vs.push_back(v);
	//}
	//find = false;
	//while (!find)
	//{
	//	int temp_int;
	//	f.getline(s, 1023);
	//	if (sscanf(s, "%s", &sread) == 1 && (strcmp(sread, "Hexahedra") == 0))
	//		find = true;
	//}
	//f.getline(s, 1023);
	//sscanf(s, "%d", &tnum);
	//int hid = 0;
	//for (int i = 0; i<tnum; i++)
	//{
	//	f.getline(s, 1023);
	//	int a, b, c, d, e, f, g, m;
	//	//sscanf(s,"%d %d %d %d %d %d %d %d %d",&vnum,&a,&b,&c,&d,&e,&f,&g,&m);
	//	sscanf(s, "%d %d %d %d %d %d %d %d %d", &a, &b, &c, &d, &e, &f, &g, &m, &vnum);
	//	Hex h;

	//	a--; b--; c--; d--; e--; f--; g--; m--;

	//	h.V_Ids[0] = a;
	//	h.V_Ids[1] = b;
	//	h.V_Ids[2] = c;
	//	h.V_Ids[3] = d;
	//	h.V_Ids[4] = e;
	//	h.V_Ids[5] = f;
	//	h.V_Ids[6] = g;
	//	h.V_Ids[7] = m;
	//	// 		h.v_ids.push_back(vnum-1);h.v_ids.push_back((int)x-1);h.v_ids.push_back((int)y-1);h.v_ids.push_back((int)z-1);
	//	// 		h.v_ids.push_back(a-1);h.v_ids.push_back(b-1);h.v_ids.push_back(c-1);h.v_ids.push_back(d-1);

	//	//  		bool all=true;//deal with deformed torus base-complex
	//	//  		for(int j = 0; j < 8; j++)
	//	//  			if(!(h.V_Ids[j] >= 0 && h.V_Ids[j] <= 961))
	//	//  				all = false;
	//	//  		if(all)
	//	//  			continue;

	//	h.index = hid++;
	//	h.frame_component = -1;
	//	Hexs.push_back(h);

	//	Vs[a].neighbor_Hs.push_back(h.index);
	//	Vs[b].neighbor_Hs.push_back(h.index);
	//	Vs[c].neighbor_Hs.push_back(h.index);
	//	Vs[d].neighbor_Hs.push_back(h.index);
	//	Vs[e].neighbor_Hs.push_back(h.index);
	//	Vs[f].neighbor_Hs.push_back(h.index);
	//	Vs[g].neighbor_Hs.push_back(h.index);
	//	Vs[m].neighbor_Hs.push_back(h.index);
	//}

	//f.close();
}

void write_obj(const std::string &filename, const MatrixXu &F,
	const MatrixXf &V, const MatrixXf &N, const MatrixXf &Nf,
	const MatrixXf &UV, const MatrixXf &C) {
	Timer<> timer;
	cout << "Writing \"" << filename << "\" (V=" << V.cols()
		<< ", F=" << F.cols() << ") .. ";
	cout.flush();
	std::ofstream os(filename);
	if (os.fail())
		throw std::runtime_error("Unable to open OBJ file \"" + filename + "\"!");
	if (N.size() > 0 && Nf.size() > 0)
		throw std::runtime_error("Please specify either face or vertex normals but not both!");

	for (uint32_t i = 0; i<V.cols(); ++i)
		os << "v " << V(0, i) << " " << V(1, i) << " " << V(2, i) << endl;

	for (uint32_t i = 0; i<N.cols(); ++i)
		os << "vn " << N(0, i) << " " << N(1, i) << " " << N(2, i) << endl;

	for (uint32_t i = 0; i<Nf.cols(); ++i)
		os << "vn " << Nf(0, i) << " " << Nf(1, i) << " " << Nf(2, i) << endl;

	for (uint32_t i = 0; i<UV.cols(); ++i)
		os << "vt " << UV(0, i) << " " << UV(1, i) << endl;

	/* Check for irregular faces */
	std::map<uint32_t, std::pair<uint32_t, std::map<uint32_t, uint32_t>>> irregular;
	size_t nIrregular = 0;

	for (uint32_t f = 0; f<F.cols(); ++f) {
		if (F.rows() == 4) {
			if (F(2, f) == F(3, f)) {
				nIrregular++;
				auto &value = irregular[F(2, f)];
				value.first = f;
				value.second[F(0, f)] = F(1, f);
				continue;
			}
		}
		os << "f ";
		for (uint32_t j = 0; j<F.rows(); ++j) {
			uint32_t idx = F(j, f);
			idx += 1;
			os << idx;
			if (Nf.size() > 0)
				idx = f + 1;
			os << "//" << idx << " ";
		}
		os << endl;
	}

	for (auto item : irregular) {
		auto face = item.second;
		uint32_t v = face.second.begin()->first, first = v;
		os << "f ";
		while (true) {
			uint32_t idx = v + 1;
			os << idx;
			if (Nf.size() > 0)
				idx = face.first + 1;
			os << "//" << idx << " ";

			v = face.second[v];
			if (v == first)
				break;
		}
		os << endl;
	}
}

void write_surface_mesh_VTK(MatrixXf &V, std::vector<std::vector<uint32_t>> &F, char * path)
{
	std::fstream f(path, std::ios::out);

	f << "# vtk DataFile Version 2.0" << std::endl << "mesh vtk data - converted from .off" << std::endl;
	f << "ASCII" << std::endl;
	f << "DATASET POLYDATA" << std::endl;

	f << "POINTS " << V.cols() << " double" << std::endl;
	for (int i = 0; i<V.cols(); i++)
		f << V(0, i) << " " << V(1, i) << " " << V(2, i) << std::endl;

	uint32_t cell = 0;
	for (int i = 0; i < F.size(); i++)
		cell += F[i].size() + 1;
	f << "POLYGONS " << F.size() << " " << cell << std::endl;

	for (int i = 0; i < F.size(); i++)
	{
		f << " " << F[i].size() << " ";
		for (int j = 0; j < F[i].size(); j++)
			f << F[i][j] << " ";
		f << std::endl;
	}
	f << "CELL_TYPES " << F.size() << std::endl;
	//f << "SCALARS cell_scalars int 1" << std::endl;
	//f << "LOOKUP_TABLE default" << std::endl;
	//for (int i = 0; i < F.size(); i++)
	//	f << i << std::endl;

	for (int i = 0; i < F.size(); i++)
	{
		if (F[i].size() == 3)
			f << 5 << std::endl;
		else if (F[i].size() == 4)
			f << 9 << std::endl;
		else if (F[i].size() >= 5)
			f << 7 << std::endl;
	}
	f.close();
}
void write_surface_mesh_OFF(MatrixXf &V, std::vector<std::vector<uint32_t>> &F, char * path)
{
	std::fstream f(path, std::ios::out);

	f << "OFF" << endl;
	f << V.cols() << " " << F.size() << " " << 0 << endl;
	for (int i = 0; i<V.cols(); i++)
		f << V(0, i) << " " << V(1, i) << " " << V(2, i) << std::endl;

	for (int i = 0; i < F.size(); i++) {
		f << F[i].size() << " ";
		for (int j = 0; j < F[i].size(); j++)
			f << F[i][j] << " ";
		f << endl;
	}
	f.close();
}
void write_surface_mesh_OBJ(MatrixXf &V, std::vector<std::vector<uint32_t>> &F, char * path)
{
	std::fstream f(path, std::ios::out);
	for (int i = 0; i<V.cols(); i++)
		f << "v " << V(0, i) << " " << V(1, i) << " " << V(2, i) << std::endl;

	for (int i = 0; i < F.size(); i++) {
		f << "f";
		for (int j = 0; j < F[i].size(); j++)
			f << " " << F[i][j] + 1;
		f << endl;
	}
	f.close();

}
void write_surface_mesh_OBJ(MatrixXf &V, MatrixXu &F, char * path)
{
	std::fstream f(path, std::ios::out);
	for (int i = 0; i<V.cols(); i++)
		f << "v " << V(0, i) << " " << V(1, i) << " " << V(2, i) << std::endl;

	for (int i = 0; i < F.cols(); i++) {
		f << "f";
		for (int j = 0; j < F.rows(); j++)
			f << " " << F(j, i) + 1;
		f << endl;
	}
	f.close();
}

void write_volume_mesh_VTK(MatrixXf &V, std::vector<std::vector<uint32_t>> &T, char * path)
{
	std::fstream f(path, std::ios::out);

	f << "# vtk DataFile Version 2.0" << std::endl << "mesh vtk data - converted from .off" << std::endl;
	f << "ASCII" << std::endl;
	f << "DATASET UNSTRUCTURED_GRID" << std::endl;

	f << "POINTS " << V.cols() << " double" << std::endl;
	for (int i = 0; i<V.cols(); i++)
		f << V(0, i) << " " << V(1, i) << " " << V(2, i) << std::endl;

	uint32_t cell = 0;
	for (int i = 0; i < T.size(); i++)
		cell += T[i].size() + 1;
	f << "CELLS " << T.size() << " " << cell << std::endl;

	for (int i = 0; i < T.size(); i++)
	{
		f << " " << T[i].size() << " ";
		for (int j = 0; j < T[i].size(); j++)
			f << T[i][j] << " ";
		f << std::endl;
	}
	f << "CELL_TYPES " << T.size() << std::endl;

	for (int i = 0; i < T.size(); i++)
	{
		if (T[i].size() == 4)
			f << 10 << std::endl;
		else if (T[i].size() == 6)
			f << 13 << std::endl;
		else if (T[i].size() >= 8)
			f << 12 << std::endl;
	}
	f.close();
}
void write_volume_mesh_VTK(MatrixXf &V, std::vector<tuple_E> &E, std::vector<std::vector<uint32_t>> &F, std::vector<int> &F_type, char * path)
{
	std::fstream f(path, std::ios::out);

	f << "# vtk DataFile Version 2.0" << std::endl << "hex-dominant mesh vtk data" << std::endl;
	f << "ASCII" << std::endl;
	f << "DATASET UNSTRUCTURED_GRID" << std::endl;
	//f << "DATASET POLYDATA" << std::endl;

	f << "POINTS " << V.cols() << " double" << std::endl;
	for (int i = 0; i<V.cols(); i++)
		f << V(0, i) << " " << V(1, i) << " " << V(2, i) << std::endl;

	//f << "LINES " << E.size() << " " << 3 * E.size() << endl;
	//for (int i = 0; i<E.size(); i++)
	//	f << 2 << " "<< std::get<0>(E[i]) << " "<< std::get<1>(E[i])<< endl;

	int polygon_num = 0;
	for (int i = 0; i<F.size(); i++)
		polygon_num += F[i].size();

	//f << "POLYGONS " << F.size() << " " << polygon_num + F.size() << endl;
	f << "CELLS " << F.size() << " " << polygon_num + F.size() << endl;

	for (int i = 0; i < F.size(); i++)
	{
		f << " " << F[i].size() << " ";
		for (int j = 0; j < F[i].size(); j++)
			f << F[i][j] << " ";
		f << std::endl;
	}

	f << "CELL_TYPES " << F.size() << std::endl;
	for (int i = 0; i < F.size(); i++)
		f << 7 << std::endl;

	std::vector<short> V_color(V.cols(), -1);
	for (int i = 0; i < E.size(); i++) {
		V_color[std::get<0>(E[i])] = std::get<4>(E[i]);
		V_color[std::get<1>(E[i])] = std::get<4>(E[i]);
	}



	f << "POINT_DATA " << V.cols() << std::endl;
	f << "SCALARS E_Scalars int" << std::endl;
	f << "LOOKUP_TABLE E_Table" << std::endl;
	for (int i = 0; i<V.cols(); i++)
		f << V_color[i] << std::endl;

	//f << "CELL_DATA " << E.size()+F_type.size() << std::endl;
	//f << "SCALARS F_Scalars int" << std::endl;
	//f << "LOOKUP_TABLE F_Table" << std::endl;
	//for (int i = 0; i<E.size(); i++)
	//	f << std::get<4>(E[i]) << std::endl;
	//for (int i = 0; i<F_type.size(); i++)
	//	f << F_type[i] << std::endl;
	f << "CELL_DATA " << F_type.size() << std::endl;
	f << "SCALARS F_Scalars float" << std::endl;
	f << "LOOKUP_TABLE F_Table" << std::endl;

	for (int i = 0; i<F_type.size(); i++)
		f << F_type[i] << std::endl;

	f.close();
}
void write_volume_mesh_HYBRID(MatrixXf &V, std::vector<std::vector<uint32_t>> &F, std::vector<std::vector<uint32_t>> &P, std::vector<bool> &P_flag, std::vector<std::vector<bool>> &PF_flag, char * path)
{
	std::fstream f(path, std::ios::out);

	f << V.cols() << " " << F.size() << " " << 3 * P.size() << std::endl;
	for (int i = 0; i<V.cols(); i++)
		f << V(0, i) << " " << V(1, i) << " " << V(2, i) << std::endl;

	for (auto fvs : F) {
		f << fvs.size() << " ";
		for (auto vid : fvs)
			f << vid << " ";
		f << std::endl;
	}

	for (uint32_t i = 0; i < P.size(); i++) {
		f << P[i].size() << " ";
		for (auto fid : P[i])
			f << fid << " ";
		f << std::endl;
		f << PF_flag[i].size() << " ";
		for (auto f_flag : PF_flag[i])
			f << f_flag << " ";
		f << std::endl;
	}

	for (uint32_t i = 0; i < P_flag.size(); i++) {
		f << P_flag[i] << std::endl;
	}

	f.close();
}
void write_volume_mesh_MESH(MatrixXf &V, std::vector<std::vector<uint32_t>> &T, char * path)
{
	std::fstream f_out_meshlab_(path, std::ios::out);

	f_out_meshlab_ << "MeshVersionFormatted 1" << std::endl;
	f_out_meshlab_ << "Dimension 3" << std::endl;
	f_out_meshlab_ << "Vertices" << " " << V.cols() << std::endl;
	for (int i = 0; i<V.cols(); i++)
		f_out_meshlab_ << V(0, i) << " " << V(1, i) << " " << V(2, i) << " " << 0 << std::endl;

	f_out_meshlab_ << "Tetrahedra" << " " << T.size() << std::endl;

	for (int i = 0; i<T.size(); i++)
	{
		f_out_meshlab_ << T[i][0] + 1 << " " << T[i][1] + 1 << " " << T[i][2] + 1 << " " << T[i][3] + 1 << " " << 0 << std::endl;
	}
	f_out_meshlab_ << "End";
	f_out_meshlab_.close();
}

void write_edge_coloring_TXT(std::vector<MatrixXf> &E_Colors, char * path) {
	std::fstream f(path, std::ios::out);
	for (uint32_t i = 0; i < E_Colors.size(); i++)
		f << E_Colors[i].cols() << " ";
	f << endl;

	for (uint32_t i = 0; i < E_Colors.size(); i++) {
		for (uint32_t j = 0; j < E_Colors[i].cols(); j++) {
			for (uint32_t k = 0; k < E_Colors[i].rows(); k++)
				//f << std::setprecision(10) << E_Colors[i](k, j) << " ";
				f << std::fixed << E_Colors[i](k, j) << " ";
			f << endl;
		}
	}
	f.close();
}
void write_singularities_SING(MatrixXf &Singularity, char * path) {
	std::fstream ff(path, std::ios::out);
	ff << Singularity.cols() / 2 << std::endl;
	for (uint32_t i = 0; i < Singularity.cols() / 2; i++) {
		for (uint32_t j = 0; j<3; j++)
			ff << std::fixed << Singularity(j, 2 * i) << " ";
		for (uint32_t j = 0; j < 3; j++)
			ff << std::fixed << Singularity(j, 2 * i + 1) << " ";
		ff << endl;
	}
	ff.close();
}
void write_Vertex_Types_TXT(std::vector<int> &V_types, char * path) {
	std::fstream ff(path, std::ios::out);
	ff << V_types.size() << std::endl;
	for (auto type : V_types) {
		ff << std::fixed << type << endl;
	}
	ff.close();
}
void write_statistics_TXT(statistics &sta, char * path) {
	std::fstream f(path, std::ios::out);

	f << 1 + 1 + 1 + sta.timings.size() + sta.polyhedral_ratios.size() << std::endl;
	f << sta.hex_ratio << "\t hex - ratio" << std::endl;
	f << sta.tN << " " << sta.tetN << "\t #triangle #tet" << std::endl;
	f << sta.hN << " " << sta.pN << "\t #hex #total" << std::endl;
	//for(auto timing:sta.timings)
	if (sta.timings.size()) {
		f << sta.timings[0] / 1000 << "\t\t " << "timing: data pre-processing" << std::endl;
		f << sta.timings[1] / 1000 << "\t\t " << "timing: rosy optimization" << std::endl;
		f << sta.timings[2] / 1000 << "\t\t " << "timing: posy optimization" << std::endl;
		f << sta.timings[3] / 1000 << "\t\t " << "timing: mesh extraction" << std::endl;
	}

	for (uint32_t i = 0; i < sta.polyhedral_nums.size(); i++) {
		if (sta.polyhedral_nums[i]>0)
			f << sta.polyhedral_ratios[i] << "\t " << i << "th polyhedral " << sta.polyhedral_nums[i] << std::endl;
	}
	f.close();
}

void load_HYBRID_mesh(Mesh &mesh, string path) {
	string filename = path;
	FILE *f = fopen(filename.data(), "rt");
	if (!f) return;
	/*check(f,"MeshVersionFormatted",1);
	check(f,"Dimension",3);*/

	int nv, np, nh;
	fscanf(f, "%d %d %d", &nv, &np, &nh);
	nh /= 3; // hack, bug in files?

	mesh.V.resize(3, nv);
	for (int i = 0; i<nv; i++) {
		Float x, y, z;
		fscanf(f, "%f %f %f", &x, &y, &z);
		mesh.V(0, i) = x;
		mesh.V(1, i) = y;
		mesh.V(2, i) = z;
	}
	mesh.Fs.resize(np);
	for (int i = 0; i<np; i++) {
		Hybrid_F &p = mesh.Fs[i];
		p.id = i;
		int nw;

		fscanf(f, "%d", &nw);
		p.vs.resize(nw);
		for (int j = 0; j<nw; j++) {
			fscanf(f, "%d", &(p.vs[j]));
		}
	}
	mesh.Hs.resize(nh);
	for (int i = 0; i<nh; i++) {
		Hybrid &h = mesh.Hs[i];
		h.id = i;

		int nf;
		fscanf(f, "%d", &nf);
		h.fs.resize(nf);

		for (int j = 0; j<nf; j++) {
			fscanf(f, "%d", &(h.fs[j]));
		}

		int tmp; fscanf(f, "%d", &tmp);
		for (int j = 0; j<nf; j++) {
			int s;
			fscanf(f, "%d", &s);
		}
	}
	for (int i = 0; i<nh; i++) {
		int tmp;
		fscanf(f, "%d", &(mesh.Hs[i].hex));
	}

	fclose(f);
}