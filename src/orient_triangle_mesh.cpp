#include "orient_triangle_mesh.h"

float	uctet(vector<float> a, vector<float> b, vector<float> c, vector<float> d)
{
	float res = 0;

	vector<float> x, y, z;

	for (int i = 0; i<3; i++)
	{
		x.push_back(b[i] - a[i]);
		y.push_back(c[i] - a[i]);
		z.push_back(d[i] - a[i]);
	}

	res = -((x[0] * y[1] * z[2] + x[1] * y[2] * z[0] + x[2] * y[0] * z[1]) - (x[2] * y[1] * z[0] + x[1] * y[0] * z[2] + x[0] * y[2] * z[1]));

	return res;
}
void orient_triangle_mesh_index(MatrixXf &Vs, MatrixXu &Ts)
{
	vector<vector<float>> p; vector<vector<int>> f;
	for (int i = 0; i<Vs.cols(); i++)
	{
		vector<float> pv;
		pv.push_back(Vs(0, i));
		pv.push_back(Vs(1, i));
		pv.push_back(Vs(2, i));
		p.push_back(pv);
	}
	for (int i = 0; i<Ts.cols(); i++)
	{
		vector<int> fv;
		fv.push_back(Ts(0, i));
		fv.push_back(Ts(1, i));
		fv.push_back(Ts(2, i));
		f.push_back(fv);
	}

	std::map<std::set<int>, vector<int>> edge_2_neb_tri;
	std::set<vector<int>> direct_edges;
	std::set<std::set<int>> no_direct_edges;
	vector<bool> whe_tri_in;

	vector<vector<int>>	nf;

	nf = f;

	for (int i = 0; i<f.size(); i++)
	{
		std::set<int> xa, xb, xc;

		xa.insert(f[i][0]);
		xa.insert(f[i][1]);

		xb.insert(f[i][1]);
		xb.insert(f[i][2]);

		xc.insert(f[i][2]);
		xc.insert(f[i][0]);

		whe_tri_in.push_back(false);

		vector<int> ct;

		ct.push_back(i);

		if (edge_2_neb_tri.find(xa) == edge_2_neb_tri.end())
		{
			edge_2_neb_tri.insert(std::pair<std::set<int>, vector<int>>(xa, ct));
		}
		else
		{
			edge_2_neb_tri[xa].push_back(i);
		}

		if (edge_2_neb_tri.find(xb) == edge_2_neb_tri.end())
		{
			edge_2_neb_tri.insert(std::pair<std::set<int>, vector<int>>(xb, ct));
		}
		else
		{
			edge_2_neb_tri[xb].push_back(i);
		}

		if (edge_2_neb_tri.find(xc) == edge_2_neb_tri.end())
		{
			edge_2_neb_tri.insert(std::pair<std::set<int>, vector<int>>(xc, ct));
		}
		else
		{
			edge_2_neb_tri[xc].push_back(i);
		}
	}

	std::set<int> xa, xb, xc;
	vector<int> ya, yb, yc;

	xa.insert(f[0][0]);
	xa.insert(f[0][1]);

	xb.insert(f[0][1]);
	xb.insert(f[0][2]);

	xc.insert(f[0][2]);
	xc.insert(f[0][0]);

	ya.push_back(f[0][0]);
	ya.push_back(f[0][1]);

	yb.push_back(f[0][1]);
	yb.push_back(f[0][2]);

	yc.push_back(f[0][2]);
	yc.push_back(f[0][0]);

	no_direct_edges.insert(xa);
	no_direct_edges.insert(xb);
	no_direct_edges.insert(xc);

	direct_edges.insert(ya);
	direct_edges.insert(yb);
	direct_edges.insert(yc);

	whe_tri_in[0] = true;

	std::queue<int> queue_loop;

	for (int i = 0; i<2; i++)
	{
		if (!whe_tri_in[edge_2_neb_tri[xa][i]])
		{
			queue_loop.push(edge_2_neb_tri[xa][i]);

			whe_tri_in[edge_2_neb_tri[xa][i]] = true;
		}

		if (!whe_tri_in[edge_2_neb_tri[xb][i]])
		{
			queue_loop.push(edge_2_neb_tri[xb][i]);

			whe_tri_in[edge_2_neb_tri[xb][i]] = true;
		}

		if (!whe_tri_in[edge_2_neb_tri[xc][i]])
		{
			queue_loop.push(edge_2_neb_tri[xc][i]);

			whe_tri_in[edge_2_neb_tri[xc][i]] = true;
		}
	}

	while (!queue_loop.empty())
	{
		xa.clear();
		xb.clear();
		xc.clear();

		ya.clear();
		yb.clear();
		yc.clear();

		int c;

		c = queue_loop.front();

		xa.insert(f[c][0]);
		xa.insert(f[c][1]);

		xb.insert(f[c][1]);
		xb.insert(f[c][2]);

		xc.insert(f[c][2]);
		xc.insert(f[c][0]);

		ya.push_back(f[c][0]);
		ya.push_back(f[c][1]);

		yb.push_back(f[c][1]);
		yb.push_back(f[c][2]);

		yc.push_back(f[c][2]);
		yc.push_back(f[c][0]);

		int cnt, ct;

		cnt = 0;
		ct = 0;

		if (no_direct_edges.find(xa) != no_direct_edges.end())
		{
			cnt++;
		}

		if (no_direct_edges.find(xb) != no_direct_edges.end())
		{
			cnt++;
		}

		if (no_direct_edges.find(xc) != no_direct_edges.end())
		{
			cnt++;
		}

		if (direct_edges.find(ya) != direct_edges.end())
		{
			ct++;
		}

		if (direct_edges.find(yb) != direct_edges.end())
		{
			ct++;
		}

		if (direct_edges.find(yc) != direct_edges.end())
		{
			ct++;
		}

		if (cnt == 0)
		{
			std::cout << "Error in triangle direction solving!" << std::endl;
			exit(0);
		}

		if (ct != 0)
		{
			ya.clear();
			yb.clear();
			yc.clear();

			ya.push_back(f[c][1]);
			ya.push_back(f[c][0]);

			yb.push_back(f[c][2]);
			yb.push_back(f[c][1]);

			yc.push_back(f[c][0]);
			yc.push_back(f[c][2]);

			nf[c][0] = f[c][1];
			nf[c][1] = f[c][0];

		}
		else
		{
			// do nothing
		}

		no_direct_edges.insert(xa);
		no_direct_edges.insert(xb);
		no_direct_edges.insert(xc);

		direct_edges.insert(ya);
		direct_edges.insert(yb);
		direct_edges.insert(yc);

		for (int i = 0; i<2; i++)
		{
			if (!whe_tri_in[edge_2_neb_tri[xa][i]])
			{
				queue_loop.push(edge_2_neb_tri[xa][i]);

				whe_tri_in[edge_2_neb_tri[xa][i]] = true;
			}

			if (!whe_tri_in[edge_2_neb_tri[xb][i]])
			{
				queue_loop.push(edge_2_neb_tri[xb][i]);

				whe_tri_in[edge_2_neb_tri[xb][i]] = true;
			}

			if (!whe_tri_in[edge_2_neb_tri[xc][i]])
			{
				queue_loop.push(edge_2_neb_tri[xc][i]);

				whe_tri_in[edge_2_neb_tri[xc][i]] = true;
			}
		}

		queue_loop.pop();
	}

	float res = 0;

	vector<float> ori;

	ori.push_back(0);
	ori.push_back(0);
	ori.push_back(0);

	for (int i = 0; i<nf.size(); i++)
	{
		res += uctet(ori, p[nf[i][0]], p[nf[i][1]], p[nf[i][2]]);
	}

	if (res > 0)
	{
		int tmi;

		for (int i = 0; i<nf.size(); i++)
		{
			tmi = nf[i][0];
			nf[i][0] = nf[i][1];
			nf[i][1] = tmi;
		}
	}

	for (int i = 0; i<Ts.cols(); i++)
	{
		Ts(0, i) = nf[i][0];
		Ts(1, i) = nf[i][1];
		Ts(2, i) = nf[i][2];
	}
}