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

#include <string>
#include <iterator>
#include <assert.h>
#include <iostream>
#include "HE_Polyhedron.h"
#include "util.h"

/// Here we implement the methods that read and save data from and to files,
/// in different formats. Detection of the file format is also done automagically.

const int LEN = 128;

filetype_e fileFormatFromFileName(const char* filename)
{
	std::string name(filename);
	for (unsigned i = 0; i<name.size(); i++)
		name[i] = tolower(name[i]);

	if (name.find(".obj") < name.size())
		return FT_OBJ;
	else if (name.find(".tri") < name.size())
		return FT_TRI;
	if (name.find(".ply") < name.size())
		return FT_PLY_ASCII; // We prefer ASCII, but this is not sure...
	else if (name.find(".off") < name.size())
		return FT_OFF;

	return FT_NONE;
}


filetype_e fileFormat(const char* filename)
{
	std::string name(filename);
	for (unsigned i = 0; i<name.size(); i++)
		name[i] = tolower(name[i]);

	// Can we detect by just looking at the filename?
	if (name.find(".obj") < name.size())
		return FT_OBJ;
	else if (name.find(".tri") < name.size())
		return FT_TRI;

	std::ifstream in(filename);
	if (!in)
		throw HE_exception("Cannot open file", filename);

	// Check the header's magic string
	std::string s;
	in >> s;
	if (s == "OFF")
	{
		assert(name.find(".off") < name.size());
		in.close();
		return FT_OFF;
	}

	if (s == "ply")
	{
		assert(name.find(".ply") < name.size());

		in >> s;
		assert(s == "format");
		in >> s;
		in.close();
		if (s == "ascii")
			return FT_PLY_ASCII;
		else
			return FT_PLY_BIN;
	}

	in.close();
	return FT_NONE;
}

bool HE_Polyhedron::readObj(const char* filename)
{
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		std::cout << "ERROR: Cannot open file " << filename << std::endl;
		return false;
	}

	class Buffer
	{
		char data[1024];
		int readLine;
	public:
		int ReadLine(FILE *fp)
		{
			char c = fgetc(fp);
			while (!feof(fp)) {
				while (isspace(c) && (!feof(fp) || c != '\0')) c = fgetc(fp);	// skip empty space
				if (c == '#') while (!feof(fp) && c != '\n' && c != '\r' && c != '\0') c = fgetc(fp);	// skip comment line
				else break;
			}
			int i = 0;
			bool inspace = false;
			while (i < 1024 - 1) {
				if (feof(fp) || c == '\n' || c == '\r' || c == '\0') break;
				if (isspace(c)) {	// only use a single space as the space character
					inspace = true;
				}
				else {
					if (inspace) data[i++] = ' ';
					inspace = false;
					data[i++] = c;
				}
				c = fgetc(fp);
			}
			data[i] = '\0';
			readLine = i;
			return i;
		}
		char& operator[](int i) { return data[i]; }
		void ReadVertex(cyPoint3f &v) const { v.Zero(); sscanf(data + 2, "%f %f %f", &v.x, &v.y, &v.z); }
		void ReadFloat3(float f[3]) const { f[2] = f[1] = f[0] = 0; int n = sscanf(data + 2, "%f %f %f", &f[0], &f[1], &f[2]); if (n == 1) f[2] = f[1] = f[0]; }
		void ReadFloat(float *f) const { sscanf(data + 2, "%f", f); }
		void ReadInt(int *i, int start) const { sscanf(data + start, "%d", i); }
		bool IsCommand(const char *cmd) const {
			int i = 0;
			while (cmd[i] != '\0') {
				if (cmd[i] != data[i]) return false;
				i++;
			}
			return (data[i] == '\0' || data[i] == ' ');
		}
		const char* Data(int start = 0) { return data + start; }
	};
	Buffer buffer;

	std::vector<cyPoint3f>	_v;		// vertices
	std::vector<std::vector<int> > faces;
	while (int rb = buffer.ReadLine(fp)) {
		if (buffer.IsCommand("v")) {
			cyPoint3f vertex;
			buffer.ReadVertex(vertex);
			_v.push_back(vertex);
		}
		else if (buffer.IsCommand("f")) {
			int facevert = -1;
			bool inspace = true;
			int type = 0;
			unsigned int index;
			std::vector<int> face;
			for (int i = 2; i < rb; i++) {
				if (buffer[i] == ' ') inspace = true;
				else {
					if (inspace) {
						inspace = false;
						type = 0;
						index = 0;
						switch (facevert) {
						case -1:
						case 0:
						default:
							face.push_back(0);
							facevert++;
							break;
						}
					}
					if (buffer[i] == '/') { type++; index = 0; }
					if (buffer[i] >= '0' && buffer[i] <= '9') {
						index = index * 10 + (buffer[i] - '0');
						switch (type) {case 0: face[facevert] = index - 1; break;}
					}
				}
			}
			faces.push_back(face);
		}
		if (feof(fp)) break;
	}

	fclose(fp);

	if (faces.size() == 0) return true; // No faces found

	std::vector<float> verts;
	for (int i = 0; i < (int)_v.size(); i++) { verts.push_back(_v[i].x); verts.push_back(_v[i].y); verts.push_back(_v[i].z); }

	loadVertices(verts);
	loadFaces(faces);
	finalize();

	return true;
}

bool HE_Polyhedron::readPly(const char* filename)
{
	std::ifstream in(filename);
	if (!in)
	{
		throw HE_exception("Cannot open file", filename);
		return false;
	}

	int nv, nf;
	short order[11];		// order of x,y,z, nx,ny,nz, red,green,blue, tu,tv vertex properties
	bool result = readPlyHeader(in, nv, nf, order);
	assert(result);

	std::vector<float> verts;
	std::vector<std::vector<int> > faces;
	result = readPlyData(in, nv, nf, order, verts, faces);
	assert(result);
	loadVertices(verts);
	loadFaces(faces);
	finalize();

	return true;
}

bool HE_Polyhedron::readPlyHeader(std::ifstream& in, int& nv, int& nf, short order[11])
{
	char buf[LEN], type[LEN], c[LEN];
	int i;

	// read ply file header
	in.getline(buf, LEN);
	assert(strncmp(buf, "ply", 3) == 0);

	in.getline(buf, LEN);
	if (strncmp(buf, "format ascii", 12) != 0)
	{
		std::cerr << "Error: Input file is not in ASCII format. Line read: " << buf << std::endl;
		return false;
	}

	in.getline(buf, LEN);
	while (strncmp(buf, "comment", 7) == 0)
		in.getline(buf, LEN);

	// read number of vertices
	if (strncmp(buf, "element vertex", 14) == 0)
	{
		sscanf(buf, "element vertex %d\n", &nv);
	}
	else
	{
		std::cerr << "Error: number of vertices expected\n";
		return false;
	}

	for (i = 0; i<11; ++i)
		order[i] = -1;

	// read vertex properties order
	i = 0;
	in.getline(buf, LEN);
	while (strncmp(buf, "property", 8) == 0)
	{
		sscanf(buf, "property %s %s\n", type, c);
		if (strncmp(c, "x", 1) == 0)
			order[0] = i;
		else if (strncmp(c, "y", 1) == 0)
			order[1] = i;
		else if (strncmp(c, "z", 1) == 0)
			order[2] = i;

		else if (strncmp(c, "nx", 2) == 0)
			order[3] = i;
		else if (strncmp(c, "ny", 2) == 0)
			order[4] = i;
		else if (strncmp(c, "nz", 2) == 0)
			order[5] = i;

		else if (strncmp(c, "red", 3) == 0)
			order[6] = i;
		else if (strncmp(c, "green", 5) == 0)
			order[7] = i;
		else if (strncmp(c, "blue", 4) == 0)
			order[8] = i;

		else if (strncmp(c, "tu", 2) == 0)
			order[9] = i;
		else if (strncmp(c, "tv", 2) == 0)
			order[10] = i;

		i++;
		in.getline(buf, LEN);
	}
	//nproperties = i;

	for (i = 0; i<3; i++)
	{
		if (order[i] < 0)
		{
			std::cerr << "Error: not enough vertex coordinate fields (nx, ny, nz)\n";
			return false;
		}
	}

	// number of faces and face properties
	if (strncmp(buf, "element face", 12) == 0)
		sscanf(buf, "element face %d\n", &nf);
	else
	{
		std::cerr << "Error: number of faces expected\n";
		return false;
	}

	in.getline(buf, LEN);
	if (strncmp(buf, "property list", 13) != 0)
	{
		std::cerr << "Error: property list expected\n";
		return false;
	}

	// Skip everything up to the end of the header
	in.getline(buf, LEN);
	while (strncmp(buf, "end_header", 10) != 0)
		in.getline(buf, LEN);

	return true;
}

bool HE_Polyhedron::readPlyData(std::ifstream& in, int& nv, int& nf, const short order[11], std::vector<float>& verts, std::vector<std::vector<int> >& faces)
{
	/////////////////////////////////////// Read vertices ///////////////////////////////////////
	char buf[LEN];
	int i, k, v;
	float values[32];
	verts.clear();
	for (i = 0; i < nv; i++)
	{
		in.getline(buf, LEN);
		sscanf(buf, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &values[0], &values[1], &values[2], &values[3],
			&values[4], &values[5], &values[6], &values[7], &values[8], &values[9], &values[10], &values[11],
			&values[12], &values[13], &values[14], &values[15]);

		verts.push_back(values[order[0]]);
		verts.push_back(values[order[1]]);
		verts.push_back(values[order[2]]);
#if 0
		if (hasnormal)
			for (j = 0; j < 3; j++)
				normals[i][j] = values[order[3 + j]];
		if (hascolor)
			for (j = 0; j < 3; j++)
				colors[i][j] = (unsigned char)values[order[6 + j]];
		if (hastexture)
			for (j = 0; j < 2; j++)
				texcoords[i][j] = values[order[9 + j]];
#endif
	}


	/////////////////////////////////////// Read faces ///////////////////////////////////////
	// read in face connectivity
	for (i = 0; i < nf; i++)
	{
		in >> k;
		std::vector<int> vi(k);
		for (int j = 0; j<k; j++)
		{
			in >> v;
			vi[j] = v;
		}

		faces.push_back(vi);
	}

	return true;
}


void HE_Polyhedron::loadVertices(const std::vector<float>& verts)
{
	assert(0 == verts.size() % 3);

	int nv = verts.size() / 3;
	_vertices.resize(nv);
	for (int i = 0; i<nv; ++i)
	{
		cyPoint3f v(verts[3 * i], verts[3 * i + 1], verts[3 * i + 2]);
		_vertices[i] = new HE_Vertex(v, i);
		assert(!_vertices[i]->edge());
	}
}

void HE_Polyhedron::loadFaces(const std::vector<std::vector<int> >& faces)
{
	for (unsigned i = 0; i<faces.size(); ++i)
		addFace(faces[i]);
}


HE_Face* HE_Polyhedron::addFace(int a, int b, int c)
{
	std::vector<int> corners(3);
	corners[0] = a;
	corners[1] = b;
	corners[2] = c;
	return addFace(corners);
}

HE_Face* HE_Polyhedron::addFace(const std::vector<int>& corners)
{
#if 0
	cerr << "Adding face #" << numFaces() << ": ";
	copy(corners.begin(), corners.end(), ostream_iterator<int>(cerr, ", "));
	cerr << endl;
#endif

	HE_Face* f = new HE_Face(0);
	f->_index = _faces.size();
	_faces.push_back(f);

	// Create one HalfEdge per consecutive pair of corners
	std::vector<HE_HalfEdge*> edges;
	const unsigned nc = corners.size();
	HE_HalfEdge* he = 0;
	for (unsigned i = 0; i<nc; ++i)
	{
		int c1 = corners[i];
		int c2 = corners[(i + 1) % nc];
		he = new HE_HalfEdge(_vertices[c2], f, _halfEdges.size());
		assert(!he->twin());
		assert(he->face());
		_halfEdges.push_back(he);
		_vertices[c1]->edge(he);
		//     cerr << c1 << " -> edge: " << he->index() << ", which connects " << c1 << ',' << c2 << endl;
		edges.push_back(he);
	}
	f->edge(he);

	// Connect each HalfEdge to its next in the face
	assert(nc == edges.size());
	const unsigned& ne = nc;
	for (unsigned e = 0; e<ne; e++)
		edges[e]->next(edges[(e + 1) % ne]);

	return f;
}
