#include "base_complex.h"


void base_complex::singularity_structure(Singularity &si, Mesh &mesh){
	si.SVs.clear(); si.SEs.clear();

	uint32_t INVALID_V = mesh.Vs.size(), INVALID_E = mesh.Vs.size();
	std::vector<uint32_t> V_flag(mesh.Vs.size(), INVALID_V), E_flag(mesh.Es.size(),false);

	for (auto &v : mesh.Vs) { v.fvid = -1; v.svid = -1; }

	uint32_t SV_count = 0, SE_count = 0;
	vector<Singular_E> circle_ses;//circular singular edges

	for (int i = 0; i < mesh.Es.size(); i++) {
		if (E_flag[i]) continue;
		if ((!mesh.Es[i].boundary && mesh.Es[i].neighbor_hs.size() == Interior_RegularE) ||
			(mesh.Es[i].boundary && mesh.Es[i].neighbor_hs.size() == Boundary_RegularE))
			continue;//non-singular e-->continue;
		
		std::function<bool(uint32_t, uint32_t, uint32_t &)> singular_proceed = [&](uint32_t vid, uint32_t eid, uint32_t &neid)->bool {
			uint32_t num1 = 0, num2 = 0;
			for (uint32_t j = 0; j<mesh.Vs[vid].neighbor_es.size(); j++){
				uint32_t cur_e = mesh.Vs[vid].neighbor_es[j];
				if (cur_e == eid) continue;
				if ((!mesh.Es[cur_e].boundary && mesh.Es[cur_e].neighbor_hs.size() != Interior_RegularE) ||
					(mesh.Es[cur_e].boundary && mesh.Es[cur_e].neighbor_hs.size() != Boundary_RegularE))
					num1++;

				if (mesh.Es[cur_e].boundary == mesh.Es[eid].boundary&&mesh.Es[cur_e].neighbor_hs.size() == mesh.Es[eid].neighbor_hs.size()){
					num2++; neid = cur_e;
				}
			}
			if (num1 == 1 && num2 == 1)
				return true;
			return false;
		};
		uint32_t v_left = mesh.Es[i].vs[0], v_right = mesh.Es[i].vs[1];
		uint32_t sv_left, sv_right;
		std::vector<uint32_t> vs_left, vs_right, es_left, es_right;
		
		bool is_circle = false;
		//left 
		es_left.push_back(i); vs_left.push_back(v_left);
		uint32_t cur_e = i, next_e = INVALID_E;
		while (singular_proceed(v_left, cur_e, next_e)) {
			cur_e = next_e;
			if (cur_e == i){ is_circle = true; break; }
			es_left.push_back(next_e);
			if (mesh.Es[cur_e].vs[0] == v_left) v_left = mesh.Es[cur_e].vs[1]; else v_left = mesh.Es[cur_e].vs[0];
			vs_left.push_back(v_left);
		}

		if (is_circle){
			Singular_E se; se.id = si.SEs.size();
			se.circle = true;
			se.es_link = es_left; se.vs_link = vs_left;
			for(uint32_t j=0;j<se.es_link.size();j++) E_flag[se.es_link[j]] = true;
			se.boundary = mesh.Es[i].boundary;
			si.SEs.push_back(se);
			circle_ses.push_back(se);
			continue;
		}

		sv_left = V_flag[v_left];
		if (V_flag[v_left] == INVALID_V) {
			sv_left = si.SVs.size(); V_flag[v_left] = sv_left;
			Singular_V sv; sv.fake = false;
			sv.id = sv_left;
			sv.hid = v_left;
			sv.boundary = mesh.Vs[v_left].boundary;
			si.SVs.push_back(sv);
		}
		//right
		vs_right.push_back(v_right);
		cur_e = i, next_e = INVALID_E;
		while (singular_proceed(v_right, cur_e, next_e)) {
			cur_e = next_e;
			if (mesh.Es[cur_e].vs[0] == v_right) v_right = mesh.Es[cur_e].vs[1]; else v_right = mesh.Es[cur_e].vs[0];
			vs_right.push_back(v_right);
			es_right.push_back(next_e);
		}
		sv_right = V_flag[v_right];
		if (V_flag[v_right] == INVALID_V) {
			sv_right = si.SVs.size(); V_flag[v_right] = sv_right;
			Singular_V sv; sv.fake = false;
			sv.id = sv_right;
			sv.hid = v_right;
			sv.boundary = mesh.Vs[v_right].boundary;
			si.SVs.push_back(sv);
		}
		//se
		Singular_E se; se.id = si.SEs.size(); se.circle = false;
		se.vs.resize(2); se.vs[0] = sv_left; se.vs[1] = sv_right;
		std::reverse(vs_left.begin(), vs_left.end());
		se.vs_link = vs_left; se.vs_link.insert(se.vs_link.end(), vs_right.begin(), vs_right.end());
		std::reverse(es_left.begin(), es_left.end());
		se.es_link = es_left; se.es_link.insert(se.es_link.end(), es_right.begin(), es_right.end());
		for (uint32_t j = 0; j<se.es_link.size(); j++) E_flag[se.es_link[j]] = true;
		se.boundary = mesh.Es[i].boundary;
		si.SEs.push_back(se);
	}

	for (auto sv : si.SVs) mesh.Vs[sv.hid].svid = sv.id;
}
void base_complex::base_complex_extraction(Singularity &si, Frame &frame, Mesh &mesh) {
	frame.FVs.clear(); frame.FEs.clear(); frame.FFs.clear(); frame.FHs.clear();

	base_complex_node_edge_extraction(si, frame,mesh);
	base_complex_face_extraction(si, frame, mesh);
	base_complex_cuboid_extraction(si, frame, mesh);

	singularity_base_complex(si, frame, mesh);
	assign_color(frame);
}
enum V_tag {
	R = -5,//regular
	S,//singular node
	E,//extra-ordinary: on singular edge
	Bn,//base-complex node
	Be//on base-complex edge
};
void base_complex::base_complex_node_edge_extraction(Singularity &si, Frame &frame, Mesh &mesh) {
	uint32_t NON_SE = mesh.Vs.size() + 1, MULTIPLE_SE = mesh.Vs.size(), INVALIDE_FE = mesh.Es.size();

	frame.FVs.clear(); frame.FEs.clear();

	std::vector<V_tag> v_tag(mesh.Vs.size(), V_tag::R);
	std::vector<uint32_t> v_neibor_se(mesh.Vs.size(), NON_SE), v_on_fe(mesh.Vs.size(), INVALIDE_FE);
	std::vector<bool> e_flag(mesh.Es.size(), false);

	for (uint32_t i = 0; i < si.SEs.size(); i++)
		for (uint32_t j = 0; j < si.SEs[i].vs_link.size(); j++) {
			uint32_t vid = si.SEs[i].vs_link[j]; v_tag[vid] = V_tag::E;
			if (v_neibor_se[vid] != NON_SE) v_neibor_se[vid] = MULTIPLE_SE;
			else v_neibor_se[vid] = i;
		}
	//nodes
	vector<uint32_t> nodes(si.SVs.size());
	for (uint32_t i = 0; i < si.SVs.size(); i++) nodes[i] = si.SVs[i].hid;
	//nodes on circular singularities
	vector<uint32_t> Snodes;
	node_on_circular_singularity(si, frame, mesh, Snodes);
	if (Snodes.size()) nodes.insert(nodes.end(), Snodes.begin(), Snodes.end());
	//node_pool
	std::queue<uint32_t> node_pool;
	for (uint32_t i = 0; i < nodes.size(); i++) {
		Frame_V fv;
		fv.id = frame.FVs.size(); fv.hid = nodes[i];
		v_tag[fv.hid] = V_tag::Bn;
		mesh.Vs[fv.hid].fvid = fv.id; frame.FVs.push_back(fv);
		node_pool.push(fv.id);
	}
	//return;
	//nodes --> edges
	while (!node_pool.empty()) {
		uint32_t fvid = node_pool.front(); node_pool.pop();
		uint32_t vid = frame.FVs[fvid].hid;
		for (uint32_t i = 0; i < mesh.Vs[vid].neighbor_es.size(); i++) {
			uint32_t eid = mesh.Vs[vid].neighbor_es[i];
			if (e_flag[eid]) continue; e_flag[eid] = true;

			uint32_t vid_next = mesh.Es[eid].vs[0]; if (vid_next == vid) vid_next = mesh.Es[eid].vs[1];

			Frame_E fe; fe.id = frame.FEs.size();
			fe.vs.push_back(fvid);
			fe.vs_link.push_back(vid); fe.vs_link.push_back(vid_next);
			fe.es_link.push_back(eid);
			//from fvid, trace along eid direction...
			uint32_t vid_pre = vid;
			while (true) {
				if (v_tag[vid_next] == V_tag::Bn) {//hit another Bc node
					fe.vs.push_back(mesh.Vs[vid_next].fvid);
					frame.FEs.push_back(fe);
					break;
				}
				else if (v_tag[vid_next] == V_tag::Be) {//hit a vertex on Be
														//new node
					Frame_V fv; fv.id = frame.FVs.size(); fv.hid = vid_next;
					v_tag[fv.hid] = V_tag::Bn;
					mesh.Vs[fv.hid].fvid = fv.id; frame.FVs.push_back(fv);
					node_pool.push(fv.id);
					//edge end
					fe.vs.push_back(fv.id); frame.FEs.push_back(fe);
					//split feid --> feid, feid1
					uint32_t feid = v_on_fe[vid_next];
					uint32_t feid_v1 = frame.FEs[feid].vs[1]; frame.FEs[feid].vs[1] = fv.id;

					Frame_E fe_new; fe_new.id = frame.FEs.size();
					fe_new.vs.push_back(fv.id); fe_new.vs.push_back(feid_v1);

					std::vector<uint32_t> &vs_link = frame.FEs[feid].vs_link, &es_link = frame.FEs[feid].es_link;
					uint32_t vid_next_pos = std::find(vs_link.begin(), vs_link.end(), vid_next) - vs_link.begin();
					fe_new.vs_link.insert(fe_new.vs_link.end(), vs_link.begin() + vid_next_pos, vs_link.end());
					fe_new.es_link.insert(fe_new.es_link.end(), es_link.begin() + vid_next_pos, es_link.end());
					vs_link.erase(vs_link.begin() + vid_next_pos + 1, vs_link.end());
					es_link.erase(es_link.begin() + vid_next_pos, es_link.end());

					for (uint32_t j = 0; j < fe_new.vs_link.size(); j++) v_on_fe[fe_new.vs_link[j]] = fe_new.id;
					frame.FEs.push_back(fe_new);
					break;
				}
				else if (v_tag[vid_next] == V_tag::E && (v_neibor_se[vid_pre] != MULTIPLE_SE && v_neibor_se[vid_next] != v_neibor_se[vid_pre])) {//hit a vertex on Se
																										//new node
					Frame_V fv; fv.id = frame.FVs.size(); fv.hid = vid_next;
					v_tag[fv.hid] = V_tag::Bn;
					mesh.Vs[fv.hid].fvid = fv.id; frame.FVs.push_back(fv);
					node_pool.push(fv.id);
					//edge end
					fe.vs.push_back(fv.id); frame.FEs.push_back(fe);
					break;
				}

				v_tag[vid_next] = V_tag::Be;
				v_on_fe[vid_next] = fe.id;
				bool find_next = false;
				for (uint32_t j = 0; j < mesh.Vs[vid_next].neighbor_es.size(); j++) {
					uint32_t neid = mesh.Vs[vid_next].neighbor_es[j];
					uint32_t nvid = mesh.Es[neid].vs[0]; if (nvid == vid_next) nvid = mesh.Es[neid].vs[1];
					if (nvid == vid_pre) continue;

					bool found = true;
					for (uint32_t k = 0; k < mesh.Es[neid].neighbor_hs.size(); k++) {
						uint32_t nhid = mesh.Es[neid].neighbor_hs[k];
						if (std::find(mesh.Hs[nhid].vs.begin(), mesh.Hs[nhid].vs.end(), vid_pre) != mesh.Hs[nhid].vs.end()){
							found = false; break;
						}
					}
					if (found) {
						vid_pre = vid_next; vid_next = nvid;
						fe.vs_link.push_back(vid_next);
						fe.es_link.push_back(neid);
						find_next = true;
						e_flag[neid] = true;
						break;
					}
				}
				if (!find_next) {
					//new node
					Frame_V fv; fv.id = frame.FVs.size(); fv.hid = vid_next;
					v_tag[fv.hid] = V_tag::Bn;
					mesh.Vs[fv.hid].fvid = fv.id; frame.FVs.push_back(fv);
					node_pool.push(fv.id);
					//edge end
					fe.vs.push_back(fv.id); frame.FEs.push_back(fe);
					break;
				}
			}
		}
	}

	for (uint32_t i = 0; i < frame.FEs.size(); i++) {
		uint32_t v0 = frame.FEs[i].vs[0], v1 = frame.FEs[i].vs[1];
		frame.FVs[v0].neighbor_fvs.push_back(v1);
		frame.FVs[v1].neighbor_fvs.push_back(v0);
		frame.FVs[v0].neighbor_fes.push_back(i);
		frame.FVs[v1].neighbor_fes.push_back(i);
	}

	for (auto fv : frame.FVs)mesh.Vs[fv.hid].fvid = fv.id;

	//char path[300] = "../datasets/orig.hex_frame.vtk";
	//h_io io;
	//io.write_Frame_VTK(frame,mesh,path);
}
void base_complex::node_on_circular_singularity(Singularity &si, Frame &frame, Mesh &mesh, vector<uint32_t> &Nodes) {
	//nodes on circular singularities
	std::vector<int> v_flag(mesh.Vs.size(), V_tag::R);
	std::vector<uint32_t> circulars;
	for (uint32_t i = 0; i < si.SEs.size(); i++) {
		if (si.SEs[i].circle) circulars.push_back(i);
		for (uint32_t j = 0; j < si.SEs[i].vs_link.size(); j++)
			if (si.SEs[i].circle) v_flag[si.SEs[i].vs_link[j]] = V_tag::E;
			else v_flag[si.SEs[i].vs_link[j]] = i;
	}
	for (uint32_t i = 0; i < si.SVs.size(); i++) v_flag[si.SVs[i].hid] = V_tag::S;

	vector<vector<uint32_t>> Node_cans(circulars.size());
	//intersection with singular nodes
	function<bool()>intersect_singularVs = [&]() -> bool {
		bool new_S = false;

		for (uint32_t i = 0; i < circulars.size(); i++) {
			
			if (Node_cans[i].size() >= 3)continue;

			auto &vs_link = si.SEs[circulars[i]].vs_link;
			auto &es_link = si.SEs[circulars[i]].es_link;
			for (uint32_t j = 0; j < vs_link.size(); j++) {
				auto vid = vs_link[j];
				if (v_flag[vid] == V_tag::S) continue;

				for (auto eid : mesh.Vs[vid].neighbor_es) {
					if (eid == es_link[j] || eid == es_link[(j + 1) % es_link.size()]) continue;

					uint32_t vid_pre = vid, vid_aft = mesh.Es[eid].vs[0]; if (vid_aft == vid_pre) vid_aft = mesh.Es[eid].vs[1];
					bool found = false;

					while (true) {
						if (v_flag[vid_aft] == V_tag::S) {
							Node_cans[i].push_back(vid); v_flag[vid] = V_tag::S;
							found = true; new_S = true; break;
						}

						bool find_next = false;
						for (auto nvid : mesh.Vs[vid_aft].neighbor_vs) {
							if (nvid == vid_pre || nvid == vid || v_flag[nvid] == V_tag::E || v_flag[nvid] >= 0) continue;
							std::sort(mesh.Vs[nvid].neighbor_fs.begin(), mesh.Vs[nvid].neighbor_fs.end());
							std::sort(mesh.Vs[vid_pre].neighbor_fs.begin(), mesh.Vs[vid_pre].neighbor_fs.end());
							vector<uint32_t> common_fs;
							std::set_intersection(mesh.Vs[vid_pre].neighbor_fs.begin(), mesh.Vs[vid_pre].neighbor_fs.end(),
								mesh.Vs[nvid].neighbor_fs.begin(), mesh.Vs[nvid].neighbor_fs.end(), std::back_inserter(common_fs));
							if (!common_fs.size()) {
								vid_pre = vid_aft; vid_aft = nvid; find_next = true; break;
							}
						}
						if (!find_next) break;
					}
					if (found) break;
				}
			}
		}
		return new_S;
	};
	function<bool()>intersect_singularEVs = [&]() -> bool {
		bool new_S = false;

		for (uint32_t i = 0; i < circulars.size(); i++) {

			if (Node_cans[i].size() >= 3)continue;

			auto &vs_link = si.SEs[circulars[i]].vs_link;
			auto &es_link = si.SEs[circulars[i]].es_link;

			vector<pair<uint32_t, uint32_t>> candidates; vector<int32_t> se_vote(si.SEs.size(), 0);

			for (uint32_t j = 0; j < vs_link.size(); j++) {
				auto vid = vs_link[j];
				if (v_flag[vid] == V_tag::S) continue;

				for (auto eid : mesh.Vs[vid].neighbor_es) {
					if (eid == es_link[j] || eid == es_link[(j + 1) % es_link.size()]) continue;

					uint32_t vid_pre = vid, vid_aft = mesh.Es[eid].vs[0]; if (vid_aft == vid_pre) vid_aft = mesh.Es[eid].vs[1];
					bool found = false;

					while (true) {
						if (v_flag[vid_aft] >= 0) {
							candidates.push_back(make_pair(vid, v_flag[vid_aft]));
							se_vote[v_flag[vid_aft]]++;
							break;
						}

						bool find_next = false;
						for (auto nvid : mesh.Vs[vid_aft].neighbor_vs) {
							if (nvid == vid_pre || nvid == vid || v_flag[nvid] == V_tag::E) continue;
							std::sort(mesh.Vs[nvid].neighbor_fs.begin(), mesh.Vs[nvid].neighbor_fs.end());
							std::sort(mesh.Vs[vid_pre].neighbor_fs.begin(), mesh.Vs[vid_pre].neighbor_fs.end());
							vector<uint32_t> common_fs;
							std::set_intersection(mesh.Vs[vid_pre].neighbor_fs.begin(), mesh.Vs[vid_pre].neighbor_fs.end(),
								mesh.Vs[nvid].neighbor_fs.begin(), mesh.Vs[nvid].neighbor_fs.end(), std::back_inserter(common_fs));
							if (!common_fs.size()) {
								vid_pre = vid_aft; vid_aft = nvid; find_next = true; break;
							}
						}
						if (!find_next) break;
					}
					if (found) break;
				}
			}
			for (uint32_t j = 0; j < se_vote.size();j++) if (se_vote[j] > 0 && se_vote[j] <= 2) {
				for (auto a_pair : candidates) if (a_pair.second == j) {
					Node_cans[i].push_back(a_pair.first); v_flag[a_pair.first] = V_tag::S;
					new_S = true; 
				}
			}
			if (new_S) {
				std::sort(Node_cans[i].begin(), Node_cans[i].end());
				Node_cans[i].erase(std::unique(Node_cans[i].begin(), Node_cans[i].end()), Node_cans[i].end());
				return true;
			}
		}
		return false;
	};
	function<bool()>intersect_parallelVs = [&]() -> bool {
		bool new_S = false;

		for (uint32_t i = 0; i < circulars.size(); i++) {

			if (Node_cans[i].size() >= 3)continue;

			auto &vs_link = si.SEs[circulars[i]].vs_link;
			auto &es_link = si.SEs[circulars[i]].es_link;

			uint32_t len = vs_link.size() / 3;
			vector<uint32_t> candidates(3);
			candidates[0]=vs_link[0];
			candidates[1]=vs_link[len];
			candidates[2]=vs_link[2 * len];

			for (uint32_t j = 0; j < 3; j++) {
				if (v_flag[candidates[j]] == V_tag::S) continue;
				Node_cans[i].push_back(candidates[j]); 
				v_flag[candidates[j]] = V_tag::S;
				new_S = true; 
				if(Node_cans[i].size() >= 3) break;
			}
			return new_S;
		}
		return false;
	};

	bool still_circular = true;
	while (still_circular) {
		while(intersect_singularVs());
		while (intersect_singularEVs()) {
			while (intersect_singularVs());
		}
		while (intersect_parallelVs()) {
			while (intersect_singularVs());
			while (intersect_singularVs());
			while (intersect_singularEVs()) {
				while (intersect_singularVs());
			}
		}
		still_circular = false;
		for (auto can : Node_cans)if (can.size() <= 2) {
			still_circular = true; break;
		}
	}
	Nodes.clear();
	for (auto cans : Node_cans) Nodes.insert(Nodes.end(), cans.begin(), cans.end());
	std::sort(Nodes.begin(), Nodes.end());
	Nodes.erase(std::unique(Nodes.begin(), Nodes.end()), Nodes.end());

	//std::vector<int> v_flag(mesh.Vs.size(), V_tag::R);
	//std::vector<uint32_t> circulars;
	//for (uint32_t i = 0; i < si.SVs.size(); i++) v_flag[si.SVs[i].hid] = V_tag::S;
	//for (uint32_t i = 0; i < si.SEs.size(); i++) {
	//	if (si.SEs[i].circle) circulars.push_back(i);
	//	for (uint32_t j = 0; j < si.SEs[i].vs_link.size(); j++)
	//		if (si.SEs[i].circle) v_flag[si.SEs[i].vs_link[j]] = i;
	//		else v_flag[si.SEs[i].vs_link[j]] = V_tag::E;
	//}
	//std::vector<bool> Vfinds(circulars.size(), false);
	//for (uint32_t i = 0; i < circulars.size(); i++) {
	//	for (uint32_t j = 0; j < si.SEs[circulars[i]].vs_link.size(); j++) {
	//		uint32_t vid = si.SEs[circulars[i]].vs_link[j];
	//		for (uint32_t k = 0; k < mesh.Vs[vid].neighbor_es.size(); k++) {
	//			uint32_t eid = mesh.Vs[vid].neighbor_es[k];
	//			if (eid == si.SEs[circulars[i]].es_link[j] || eid == si.SEs[circulars[i]].es_link[(j + 1) % si.SEs[circulars[i]].es_link.size()])
	//				continue;
	//			uint32_t vid_pre = vid, vid_aft = mesh.Es[eid].vs[0];
	//			if (vid_aft == vid_pre) vid_aft = mesh.Es[eid].vs[1];
	//			bool found = false;
	//			while (true) {
	//				if (v_flag[vid_aft] == V_tag::S || v_flag[vid_aft] == V_tag::E) {
	//					nodes.push_back(vid);  found = true;
	//					break;
	//				}
	//				bool find_next = false;
	//				for (uint32_t m = 0; m < mesh.Vs[vid_aft].neighbor_vs.size(); m++) {
	//					uint32_t nvid = mesh.Vs[vid_aft].neighbor_vs[m];
	//					if (nvid == vid_pre || nvid == vid) continue;
	//					std::sort(mesh.Vs[nvid].neighbor_fs.begin(), mesh.Vs[nvid].neighbor_fs.end());
	//					std::sort(mesh.Vs[vid_pre].neighbor_fs.begin(), mesh.Vs[vid_pre].neighbor_fs.end());
	//					std::vector<uint32_t> common_fs;
	//					std::set_intersection(mesh.Vs[vid_pre].neighbor_fs.begin(), mesh.Vs[vid_pre].neighbor_fs.end(),
	//						mesh.Vs[nvid].neighbor_fs.begin(), mesh.Vs[nvid].neighbor_fs.end(), std::back_inserter(common_fs));
	//					if (!common_fs.size()) {
	//						vid_pre = vid_aft; vid_aft = nvid; find_next = true; break;
	//					}
	//				}
	//				if (!find_next) break;
	//			}
	//			if (found) { Vfinds[i] = true; break; }
	//		}
	//	}
	//	if (Vfinds[i]) continue;
	//	uint32_t circle_len = si.SEs[circulars[i]].vs_link.size() / 3;
	//	nodes.push_back(si.SEs[circulars[i]].vs_link[0]);
	//	nodes.push_back(si.SEs[circulars[i]].vs_link[circle_len]);
	//	nodes.push_back(si.SEs[circulars[i]].vs_link[2 * circle_len]);
	//}
	//if (circulars.size()) {
	//	std::sort(nodes.begin(), nodes.end()); nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
	//}
}
void base_complex::base_complex_face_extraction(Singularity &si, Frame &frame, Mesh &mesh)
{
	frame.FFs.clear();
	//faces from edges.
	uint32_t INVALID_E = frame.FEs.size(), INVALID_F = mesh.Fs.size();
	std::vector<uint32_t> e_tag(mesh.Es.size(), INVALID_E);
	for (uint32_t i = 0; i < frame.FEs.size(); i++) 
		for (uint32_t j = 0; j < frame.FEs[i].es_link.size(); j++) e_tag[frame.FEs[i].es_link[j]] = i;
	std::vector<uint32_t> f_flag(mesh.Fs.size(), INVALID_F);
	std::vector<bool> fe_flag(frame.FEs.size(), false);
	
	for (uint32_t i = 0; i < frame.FEs.size(); i++) {
		uint32_t eid = frame.FEs[i].es_link[0];
		for (uint32_t j = 0; j < mesh.Es[eid].neighbor_fs.size(); j++) {
			uint32_t fid = mesh.Es[eid].neighbor_fs[j];
			if (f_flag[fid] != INVALID_F) continue;

			Frame_F ff; ff.id = frame.FFs.size();
			ff.boundary = mesh.Fs[fid].boundary;
			ff.es.push_back(i); fe_flag[i] = true;

			std::queue<uint32_t> f_pool; f_pool.push(fid);
			while (!f_pool.empty()) {
				fid = f_pool.front(); f_pool.pop();
				if (f_flag[fid] != INVALID_F) continue;
				f_flag[fid] = ff.id;
				ff.ffs_net.push_back(fid);

				for (uint32_t k = 0; k < 4; k++) {
					uint32_t feid = mesh.Fs[fid].es[k];
					if (e_tag[feid] != INVALID_E) { 
						if (!fe_flag[e_tag[feid]]) {
							fe_flag[e_tag[feid]] = true;
							ff.es.push_back(e_tag[feid]);
						}
						continue; 
					}
					for (uint32_t m = 0; m < mesh.Es[feid].neighbor_fs.size(); m++) {
						uint32_t enfid = mesh.Es[feid].neighbor_fs[m];
						if (f_flag[enfid] != INVALID_F) continue;
						bool pass = true;
						for (uint32_t n = 0; n < mesh.Fs[enfid].neighbor_hs.size(); n++) {
							uint32_t fnhid = mesh.Fs[enfid].neighbor_hs[n];
							for (uint32_t p = 0; p < 6; p++) {
								uint32_t cur_fid = mesh.Hs[fnhid].fs[p];
								if (f_flag[cur_fid] == ff.id) { pass = false; break; }
							}
							if (!pass) break;
						}
						if (pass) f_pool.push(enfid);
					}
				}
			}
			frame.FFs.push_back(ff);
			for (uint32_t k = 0; k < ff.es.size(); k++) fe_flag[ff.es[k]] = false;
		}
	}
	//re-order es
	for (uint32_t i = 0; i < frame.FFs.size(); i++) {

		uint32_t e0 = frame.FFs[i].es[0], e1 = -1, e2 = -1, e3 = -1;
		uint32_t v0 = frame.FEs[e0].vs[0], v1 = frame.FEs[e0].vs[1], v2 = -1, v3 = -1;
		//e1, v2
		for (uint32_t j = 1; j < 4; j++) {
			uint32_t eid = frame.FFs[i].es[j], v0_= frame.FEs[eid].vs[0], v1_ = frame.FEs[eid].vs[1];
			if (v0_ == v1 || v1_ == v1) {
				if (v0_ == v1) v2 = v1_; else v2 = v0_;
				e1 = eid; break;
			}
		}
		//e3, v3
		for (uint32_t j = 1; j < 4; j++) {
			uint32_t eid = frame.FFs[i].es[j], v0_ = frame.FEs[eid].vs[0], v1_ = frame.FEs[eid].vs[1];
			if (v0_ == v0 || v1_ == v0) {
				if (v0_ == v0) v3 = v1_; else v3 = v0_;
				e3 = eid; break;
			}
		}


		for (uint32_t j = 1; j < 4; j++) 
			if (frame.FFs[i].es[j] != e1 && frame.FFs[i].es[j] != e3) 
				e2 = frame.FFs[i].es[j];
		frame.FFs[i].es[0] = e0;
		frame.FFs[i].es[1] = e1;
		frame.FFs[i].es[2] = e2;
		frame.FFs[i].es[3] = e3;

		for (uint32_t j = 0; j < 4; j++) frame.FEs[frame.FFs[i].es[j]].neighbor_ffs.push_back(i);

		frame.FFs[i].vs.resize(4);
		frame.FFs[i].vs[0] = v0;
		frame.FFs[i].vs[1] = v1;
		frame.FFs[i].vs[2] = v2;
		frame.FFs[i].vs[3] = v3;

		for (uint32_t j = 0; j < 4; j++) frame.FVs[frame.FFs[i].vs[j]].neighbor_ffs.push_back(i);
	}
}
void base_complex::base_complex_cuboid_extraction(Singularity &si, Frame &frame, Mesh &mesh)
{
	frame.FHs.clear();

	uint32_t INVALID_F = mesh.Fs.size(), INVALID_H=mesh.Hs.size();
	std::vector<uint32_t> f_flag(mesh.Fs.size(), INVALID_F), h_flag(mesh.Hs.size(), INVALID_H);
	std::vector<bool> ff_flag(frame.FFs.size(),false);

	for (uint32_t i = 0; i < frame.FFs.size(); i++) for (uint32_t j = 0; j < frame.FFs[i].ffs_net.size(); j++)
		f_flag[frame.FFs[i].ffs_net[j]] = i;

	while (true) {
		Frame_H fh; fh.id = frame.FHs.size(); fh.Color_ID = -1;
		uint32_t start_h = INVALID_H;
		for (uint32_t i = 0; i < h_flag.size(); i++) if (h_flag[i] == INVALID_H) {start_h = i; break;}
		if (start_h == INVALID_H) break;

		std::queue<uint32_t> h_pool; h_pool.push(start_h);
		while (!h_pool.empty()) {
			start_h = h_pool.front(); h_pool.pop();
			if (h_flag[start_h] != INVALID_H) continue;
			h_flag[start_h] = fh.id;
			fh.hs_net.push_back(start_h);
			for (uint32_t i = 0; i < 6; i++) {
				uint32_t fid = mesh.Hs[start_h].fs[i];
				if (f_flag[fid] != INVALID_F) {
					if (!ff_flag[f_flag[fid]]) {
						ff_flag[f_flag[fid]] = true;
						fh.fs.push_back(f_flag[fid]);
					}
					continue;
				}
				for (uint32_t j = 0; j < mesh.Fs[fid].neighbor_hs.size(); j++) {
					uint32_t hid = mesh.Fs[fid].neighbor_hs[j];
					if (h_flag[hid] != INVALID_H) continue;
					h_pool.push(hid);
				}
			}
		}
		frame.FHs.push_back(fh);
		for (uint32_t k = 0; k < fh.fs.size(); k++) ff_flag[fh.fs[k]] = false;
	}

	for (uint32_t i = 0; i < frame.FHs.size(); i++) {
		frame.FHs[i].es.reserve(12);
		for (uint32_t j = 0; j < frame.FHs[i].fs.size(); j++) {
			uint32_t fid = frame.FHs[i].fs[j];
			for (uint32_t k = 0; k < 4; k++) frame.FHs[i].es.push_back(frame.FFs[fid].es[k]);
		}
		std::sort(frame.FHs[i].es.begin(), frame.FHs[i].es.end());
		frame.FHs[i].es.erase(std::unique(frame.FHs[i].es.begin(), frame.FHs[i].es.end()), frame.FHs[i].es.end());
		frame.FHs[i].vs= frame.FFs[frame.FHs[i].fs[0]].vs;
		std::vector<uint32_t> vs = frame.FFs[frame.FHs[i].fs[0]].vs;
		std::sort(vs.begin(), vs.end());
		short cors_f = 1;
		for (uint32_t j = 1; j < 6; j++) {
			std::vector<uint32_t> vsj = frame.FFs[frame.FHs[i].fs[j]].vs;
			std::sort(vsj.begin(), vsj.end());
			std::vector<uint32_t> common_vs;
			std::set_intersection(vs.begin(), vs.end(), vsj.begin(), vsj.end(), std::back_inserter(common_vs));
			if (common_vs.size())continue; else { cors_f = j; break; }
		}
		vs = frame.FFs[frame.FHs[i].fs[cors_f]].vs;
		for (uint32_t j = 0; j < 4; j++) {
			std::vector<uint32_t> nvs = frame.FVs[frame.FHs[i].vs[j]].neighbor_fvs;
			for (uint32_t k = 0; k < nvs.size(); k++)
				if (std::find(vs.begin(), vs.end(), nvs[k]) != vs.end()){
					frame.FHs[i].vs.push_back(nvs[k]); break;
				}
		}

		for (uint32_t j = 0; j < 8; j++)frame.FVs[frame.FHs[i].vs[j]].neighbor_fhs.push_back(i);
		for (uint32_t j = 0; j < 12; j++)frame.FEs[frame.FHs[i].es[j]].neighbor_fhs.push_back(i);
		for (uint32_t j = 0; j < 6; j++)frame.FFs[frame.FHs[i].fs[j]].neighbor_fhs.push_back(i);
	}

	for (auto &v : frame.FVs) v.boundary = false;
	for (auto &e : frame.FEs) e.boundary = false;
	for (auto &f : frame.FFs) {
		f.boundary = false;
		if (f.neighbor_fhs.size() == 1) {
			f.boundary = true;
			for (auto vid : f.vs) frame.FVs[vid].boundary = true;
			for (auto eid : f.es) frame.FEs[eid].boundary = true;
		}
	}
}
void base_complex::singularity_base_complex(Singularity &si, Frame &frame, Mesh &mesh) {

	for (auto &fv : frame.FVs) {
		if (mesh.Vs[fv.hid].svid == -1) fv.svid = -1;
		else fv.svid = mesh.Vs[fv.hid].svid;
	}
	for (auto &fe : frame.FEs) {
		if ((!fe.boundary && fe.neighbor_fhs.size() == Interior_RegularE) ||
			(fe.boundary && fe.neighbor_fhs.size() == Boundary_RegularE))
			fe.singular = false;
		else
			fe.singular = true;
	}
}

void base_complex::assign_color(Frame &frame)
{
	int Color_id = 0; int Color_Total = 8;
	for (uint32_t i = 0; i<frame.FHs.size(); i++){
		if (frame.FHs[i].Color_ID == -1){
			vector<int> ids;
			for (auto fid:frame.FHs[i].fs){
				for (auto hid:frame.FFs[fid].neighbor_fhs){
					if (hid != i&&frame.FHs[hid].Color_ID != -1)
					{
						ids.push_back(frame.FHs[hid].Color_ID);
					}
				}
			}
			while (true){
				//Ind_C+=2;
				if (find(ids.begin(), ids.end(), Color_id) != ids.end()){
					Color_id++;
					Color_id = Color_id%Color_Total;
				}
				else
					break;
			}
			frame.FHs[i].Color_ID = Color_id;
		}
	}
}
