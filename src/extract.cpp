#include "hierarchy.h"
#include "positions.h"
//2D&3D========================================================================================================//
std::priority_queue<tuple_E, std::vector<tuple_E>, LessThan> Es_reddash;
std::vector<uint32_t> V_map;
std::vector<std::vector<uint32_t>> Reverse_V_map;

std::vector<tuple_E> mEs;
std::vector<tuple_F> mFs;
std::vector<std::vector<uint32_t>> mFs2D;
std::vector<std::vector<uint32_t>> FEs;

vector<vector<uint32_t>> V_pvs, V_pes, E_pfs;

std::vector<short> mV_flag;
std::vector<bool> mE_flag, mF_flag;
//re-coloring
MatrixXf mQ_copy2D, mO_copy2D, mN_copy2D, mV_copy2D;
MatrixXf newQ2D, newN2D, newV2D;
bool non_manifold = false;
//2D===========================================================================================================//
void construct_Es_FEs()
{
	FEs.clear(); mEs.clear(); V_pvs.clear(); E_pfs.clear();

	std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, int>> temp;
	temp.reserve(mFs2D.size() * 3);
	FEs.resize(mFs2D.size());
	for (uint32_t f = 0; f < mFs2D.size(); ++f) {
		for (uint32_t e = 0; e < mFs2D[f].size(); ++e) {
			uint32_t v0 = mFs2D[f][e], v1 = mFs2D[f][(e + 1) % mFs2D[f].size()];
			if (v0 > v1) std::swap(v0, v1);
			temp.push_back(std::make_tuple(v0, v1, f, e, Edge_tag::B));
		}
		std::vector<uint32_t> fes(mFs2D[f].size());
		FEs[f] = fes;
	}
	std::sort(temp.begin(), temp.end());
	mEs.reserve(temp.size() / 2);
	int E_num = -1;
	for (uint32_t i = 0; i < temp.size(); ++i) {
		if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
			std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
			E_num++;
			mEs.push_back(std::make_tuple(std::get<0>(temp[i]), std::get<1>(temp[i]), true, 0, std::get<4>(temp[i]), E_num, -1, 0));
		}
		else if (i != 0 && (std::get<0>(temp[i]) == std::get<0>(temp[i - 1]) &&
			std::get<1>(temp[i]) == std::get<1>(temp[i - 1])))
			std::get<2>(mEs[E_num]) = false;

		FEs[std::get<2>(temp[i])][std::get<3>(temp[i])] = E_num;
	}

	if (mV_copy2D.cols()) {
		V_pvs.resize(mV_copy2D.cols());
		for (uint32_t i = 0; i < mEs.size(); i++) {
			uint32_t v0 = get<0>(mEs[i]), v1 = get<1>(mEs[i]);
			V_pvs[v0].push_back(v1);
			V_pvs[v1].push_back(v0);
		}
	}
	E_pfs.resize(mEs.size());	
	for (uint32_t i = 0; i < FEs.size(); i++) for (auto eid : FEs[i]) { E_pfs[eid].push_back(i); }
}
bool simple_polygon(std::vector<std::vector<uint32_t>> &fvs, std::vector<std::vector<uint32_t>> &fes, std::vector<uint32_t> &pvs, 
	std::vector<uint32_t> &pes, std::vector<uint32_t> &vs_disgard, std::vector<uint32_t> &es_disgard)
{
	es_disgard.clear();//es_disgard is in the interior
	for (int i = 0; i < fes.size(); i++) {
		for (int j = 0; j < fes[i].size(); j++) 
			if (mE_flag[fes[i][j]]) {
				es_disgard.push_back(fes[i][j]); 
				mE_flag[fes[i][j]] = false;
			}
			else mE_flag[fes[i][j]] = true;
	}
	short which_polygon = 0;
	for (int i = 0; i < fes.size(); i++) {
		for (int j = 0; j < fes[i].size(); j++)
			if (mE_flag[fes[i][j]]) { 
				if (!pes.size()) 
					which_polygon = i;

				pes.push_back(fes[i][j]); mE_flag[fes[i][j]] = false;
			}
	}
	//test nvs for each v
	for (uint32_t i = 0; i < pes.size(); ++i) {
		uint32_t v0 = V_map[std::get<0>(mEs[pes[i]])], v1 = V_map[std::get<1>(mEs[pes[i]])];
		mV_flag[v0]++; mV_flag[v1]++;
		if (mV_flag[v0] > 2 || mV_flag[v1] > 2) {
			for (uint32_t j = 0; j < pes.size(); ++j) {
				uint32_t v0_ = V_map[std::get<0>(mEs[pes[j]])], v1_ = V_map[std::get<1>(mEs[pes[j]])];
				mV_flag[v0_] = mV_flag[v1_] = 0;
			}
			return false;
		}
	}
	for (uint32_t i = 0; i < pes.size(); ++i) {
		uint32_t v0 = V_map[std::get<0>(mEs[pes[i]])], v1 = V_map[std::get<1>(mEs[pes[i]])];
		mV_flag[v0] = mV_flag[v1] = 0;
	}
	//extract the polygon	
	pvs.clear();
	pvs.reserve(pes.size());
	std::vector<bool> e_flag(pes.size(), false);
	pvs.push_back(V_map[std::get<0>(mEs[pes[0]])]);
	pvs.push_back(V_map[std::get<1>(mEs[pes[0]])]);
	e_flag[0] = true;
	uint32_t start_v = pvs[1];
	for (uint32_t i = 2; i < pes.size(); i++) {
		for (uint32_t j = 1; j < pes.size(); j++) {
			if (!e_flag[j]) {
				if (V_map[std::get<0>(mEs[pes[j]])] == start_v) {
					e_flag[j] = true;
					pvs.push_back(V_map[std::get<1>(mEs[pes[j]])]);
					start_v = V_map[std::get<1>(mEs[pes[j]])];
					break;
				}
				else if (V_map[std::get<1>(mEs[pes[j]])] == start_v) {
					e_flag[j] = true;
					pvs.push_back(V_map[std::get<0>(mEs[pes[j]])]);
					start_v = V_map[std::get<0>(mEs[pes[j]])];
					break;
				}
			}
		}
	}
	if (pvs.size() != pes.size())
		return false;
	//judge direction
	bool correct = true;
	for (int i = 0; i < fvs[which_polygon].size(); i++)
		if (fvs[which_polygon][i] == V_map[std::get<0>(mEs[pes[0]])])
			if (fvs[which_polygon][(i + 1) % fvs[which_polygon].size()] != V_map[std::get<1>(mEs[pes[0]])])
				correct = false;
	if (!correct)
		std::reverse(pvs.begin(), pvs.end());
	//vs_disgard
	for (uint32_t i = 0; i < fvs.size(); ++i)
		for (uint32_t j = 0; j < fvs[i].size(); ++j)
			mV_flag[fvs[i][j]] = 1; 
	for (uint32_t i = 0; i < pvs.size(); ++i) mV_flag[pvs[i]] = 0;
	for (uint32_t i = 0; i < fvs.size(); ++i)
		for (uint32_t j = 0; j < fvs[i].size(); ++j)
			if (mV_flag[fvs[i][j]] == 1) {
				vs_disgard.push_back(fvs[i][j]); mV_flag[fvs[i][j]] = 0;
			}
	return true;
}
void reindex_2D(MatrixXf &Vs, std::vector<std::vector<uint32_t>> &F_Vs)
{
	std::vector<int> V_flag(Vs.size(), -1);
	for (uint32_t i = 0; i < F_Vs.size(); i++)
		for (uint32_t j = 0; j < F_Vs[i].size(); j++)
			V_flag[F_Vs[i][j]] = 0;

	uint32_t v_num = 0;
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1)
			V_flag[i] = v_num++;
	MatrixXf mV_local_(3, v_num);
	v_num = 0;
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1)
			mV_local_.col(v_num++) = Vs.col(i);
	std::vector<std::vector<uint32_t>> mFs_local_(F_Vs.size());
	for (uint32_t i = 0; i < F_Vs.size(); i++)
		for (uint32_t j = 0; j < F_Vs[i].size(); j++)
			mFs_local_[i].push_back(V_flag[F_Vs[i][j]]);

	mV_local_.swap(Vs);
	mFs_local_.swap(F_Vs);
}
void topology_check_2D(std::vector<std::vector<uint32_t>> &F_Vs, std::vector<std::vector<uint32_t>> &F_Es, int &genus, bool &manifoldness)
{
	std::set<uint32_t> Vs_;
	for (uint32_t i = 0; i < F_Vs.size(); i++) {
		for (uint32_t j = 0; j < F_Vs[i].size(); j++)
			Vs_.insert(F_Vs[i][j]);
	}
	int v_num = Vs_.size();

	std::set<uint32_t> Es_; uint32_t max_eid = 0;
	for (uint32_t i = 0; i < F_Es.size(); i++) {
		for (uint32_t j = 0; j < F_Es[i].size(); j++) {
			Es_.insert(F_Es[i][j]);
			if (max_eid < F_Es[i][j])
				max_eid = F_Es[i][j];
		}
	}
	int e_num = Es_.size();

	std::vector<std::vector<uint32_t>> E_nts(max_eid+1);
	for (uint32_t i = 0; i < F_Es.size(); i++)
		for (uint32_t j = 0; j < F_Es[i].size(); j++)
			E_nts[F_Es[i][j]].push_back(i);

	genus = -(v_num + F_Vs.size() - e_num - 2) / 2;

	manifoldness = true;
	for (uint32_t i = 0; i < E_nts.size(); i++)
		if (E_nts[i].size() > 2) manifoldness = false;
}

void MultiResolutionHierarchy::init_edge_tagging2D() {
	mFs2D.clear();
	mFs2D.resize(mF.cols());
	for (uint32_t i = 0; i < mF.cols(); ++i) for (int j = 0; j < 3; ++j) mFs2D[i].push_back(mF(j, i));
	construct_Es_FEs();

	for (uint32_t i = 0; i < mEs.size(); ++i) {
		uint32_t v0 = std::get<0>(mEs[i]), v1 = std::get<1>(mEs[i]);
		Vector3f q_cur = mQ[0].col(v0);;
		Vector3f q_next = applyRotationKeep(q_cur, mN[0].col(v0), mQ[0].col(v1), mN[0].col(v1));
		std::pair<int, Float> a_pair= assignColorWeighted2D(mO[0].col(v0), q_cur, mN[0].col(v0), mV[0].col(v0),
			mO[0].col(v1), q_next, mN[0].col(v1), mV[0].col(v1),
			mScale, mInvScale);
		std::get<4>(mEs[i]) = a_pair.first;
		std::get<3>(mEs[i]) = a_pair.second;

		const Vector3f o0 = mO[0].col(v0), o1 = mO[0].col(v1);
		Float energy = (o0 - o1).norm();
	}

	E_rend.resize(6, mEs.size() * 2);
	composit_edges_colors(mV[0], mEs, E_rend);

	std::vector<std::vector<uint32_t>> PV_npes(mV[0].cols());
	for (uint32_t i = 0; i < mEs.size(); i++) {
		uint32_t v0 = std::get<0>(mEs[i]), v1 = std::get<1>(mEs[i]);
		PV_npes[v0].push_back(i); PV_npes[v1].push_back(i);
	}
	std::vector<bool> V_flag_(mV[0].cols(), false);
	mO_center.resize(3, mV[0].cols());
	while (true) {
		std::vector<uint32_t> v_pool, v_set;
		for (uint32_t i = 0; i < V_flag_.size(); i++) if (!V_flag_[i]) { v_pool.push_back(i); break; }
		if (!v_pool.size()) break;
		v_set = v_pool;
		V_flag_[v_pool[0]] = true;
		while (v_pool.size()) {
			std::vector<uint32_t> v_pool_sudo;
			for (uint32_t j = 0; j < v_pool.size(); j++)
				for (uint32_t k = 0; k < PV_npes[v_pool[j]].size(); k++) {
					uint32_t eid = PV_npes[v_pool[j]][k];
					uint32_t v0 = std::get<0>(mEs[eid]), v1 = std::get<1>(mEs[eid]);
					if (std::get<4>(mEs[eid]) == Edge_tag::R) {
						if (!V_flag_[v0]) v_pool_sudo.push_back(v0);
						if (!V_flag_[v1]) v_pool_sudo.push_back(v1);
					}
				}
			v_pool.clear();
			if (v_pool_sudo.size()) {
				v_pool.clear();
				for (uint32_t j = 0; j < v_pool_sudo.size(); j++) if (!V_flag_[v_pool_sudo[j]]) {
					v_pool.push_back(v_pool_sudo[j]); V_flag_[v_pool_sudo[j]] = true;
				}
				v_set.insert(v_set.end(), v_pool.begin(), v_pool.end());
			}
		}
		Vector3f center; center.setZero();
		for (uint32_t j = 0; j < v_set.size(); j++)
			center += mO[0].col(v_set[j]);
		center /= v_set.size();
		for (uint32_t j = 0; j < v_set.size(); j++)
			mO_center.col(v_set[j]) = center;
	}

	E_O_rend.resize(6, mEs.size() * 2);
	composit_edges_colors(mO_center, mEs, E_O_rend);

	E_I_rend = E_rend;
}

bool MultiResolutionHierarchy::tagging_collapseTri(bool triangle_Switch)
{
	uint32_t INVALID_V = mV_copy2D.cols(), INVALID_E = mEs.size();
	std::vector<uint32_t> E_map(mEs.size());
	std::vector<std::vector<uint32_t>> Reverse_E_map(mEs.size());
	std::vector<std::vector<uint32_t>> Reverse_V_nts(mV_copy2D.cols());
	std::vector<std::vector<uint32_t>> Reverse_E_nts(mEs.size());
	std::vector<bool> mV_B_flag(mV_copy2D.cols(), false);//boundary flag											 
	std::vector<int32_t> F_map(mFs2D.size(), -1);
	//V_map, Reverse_V_map, E_map, Reverse_E_map
	V_map.resize(mV_copy2D.cols()); Reverse_V_map.clear(); Reverse_V_map.resize(mV_copy2D.cols());
	for (uint32_t i = 0; i < mV_copy2D.cols(); i++) { V_map[i] = i; Reverse_V_map[i].push_back(i); }
	for (uint32_t i = 0; i < mEs.size(); i++) { E_map[i] = i; Reverse_E_map[i].push_back(i); }
	//Reverse_V_nts, Reverse_E_nts
	for (uint32_t i = 0; i < mFs2D.size(); i++)
		for (uint32_t j = 0; j < mFs2D[i].size(); ++j) { Reverse_V_nts[mFs2D[i][j]].push_back(i); Reverse_E_nts[FEs[i][j]].push_back(i); }
	//mV_B_flag
	for (auto e : mEs) if (std::get<2>(e) == 1) { mV_B_flag[std::get<1>(e)] = mV_B_flag[std::get<0>(e)] = true; }
	//mV_flag, mE_flag, mF_flag
	mV_flag.resize(mV_copy2D.size()); std::fill(mV_flag.begin(), mV_flag.end(), 0);
	mE_flag.resize(mEs.size()); std::fill(mE_flag.begin(), mE_flag.end(), false);
	mF_flag.resize(mFs2D.size()); std::fill(mF_flag.begin(), mF_flag.end(), true);

	//re-coloring
	std::vector<uint32_t> E_TimeStamp(mEs.size(), 0);
	std::vector<std::vector<uint32_t>> pV_npes(mV_copy2D.cols());
	for (auto e : mEs) {
		pV_npes[std::get<0>(e)].push_back(std::get<5>(e));
		pV_npes[std::get<1>(e)].push_back(std::get<5>(e));
	}
	//triangle
	uint32_t Remove_Tri_Num = 1; bool entered_triangle = false;
	uint32_t Degenerate_triangle = 1;
	//
	uint32_t iteration = 0, once = false;
	bool topology = true; uint32_t Es_reddash_N = Es_reddash.size();
	while ((!Es_reddash.empty() && topology)||(Remove_Tri_Num && triangle_Switch) ) {
		iteration++;
		topology = false;
		uint32_t size_pool = Es_reddash.size();
		std::vector<tuple_E> Es_nonmanifold;
		for (uint32_t i = 0; i < size_pool; i++) {
			tuple_E e = Es_reddash.top();
			Es_reddash.pop();
			uint32_t e_map = E_map[std::get<5>(e)];
			uint32_t v0 = std::get<0>(e), v1 = std::get<1>(e);
			uint32_t v0_map = V_map[v0], v1_map = V_map[v1];

			if (e_map == INVALID_E) { std::get<4>(mEs[std::get<5>(e)]) = std::get<4>(e); continue; }

			if (std::get<7>(e) < E_TimeStamp[e_map]) continue;

			if (!Reverse_E_nts[e_map].size()) { std::get<4>(mEs[e_map]) = std::get<4>(e); continue; }
			if (std::get<4>(mEs[e_map]) != Edge_tag::B) continue;
			if (v0_map == v1_map) { std::cout << "somewhere is wrong" << endl; system("PAUSE"); }

			

			if (std::get<4>(e) == Edge_tag::D) { //edge fuse to form a simple polygon
				if (Reverse_E_nts[e_map].size() == 1) continue;
				std::vector<std::vector<uint32_t>> fvs, fes;
				for (auto fid : Reverse_E_nts[e_map]) {
					fvs.push_back(mFs2D[fid]); fes.push_back(FEs[fid]);
				}
				std::vector<uint32_t> pvs, pes, vs_disgard, es_disgard;
				if (!simple_polygon(fvs, fes, pvs, pes, vs_disgard, es_disgard)) {

					std::get<3>(e) *= 1.2; Es_nonmanifold.push_back(e); continue;
				}
				std::get<4>(mEs[e_map]) = Edge_tag::D;

				uint32_t f0 = Reverse_E_nts[e_map][0], f1 = Reverse_E_nts[e_map][1];
				mFs2D[f0] = pvs; FEs[f0] = pes; mF_flag[f1] = false; mFs2D[f1].clear(); FEs[f1].clear();
				//pvs, pes
				for (auto eid : pes) std::replace(Reverse_E_nts[eid].begin(), Reverse_E_nts[eid].end(), f1, f0);
				for (auto vid : pvs) {
					std::replace(Reverse_V_nts[vid].begin(), Reverse_V_nts[vid].end(), f1, f0);
					std::sort(Reverse_V_nts[vid].begin(), Reverse_V_nts[vid].end());
					Reverse_V_nts[vid].erase(std::unique(Reverse_V_nts[vid].begin(), Reverse_V_nts[vid].end()), Reverse_V_nts[vid].end());
				}
				//es_disgard, vs_disgard
				for (auto eid : es_disgard) {
					Reverse_E_nts[eid].clear();
					for (auto reid : Reverse_E_map[eid]) E_map[reid] = INVALID_E;
					Reverse_E_map[eid].clear();
				}
				for (auto vid : vs_disgard) {
					Reverse_V_nts[vid].clear();
					for (auto rvid : Reverse_V_map[vid]) V_map[rvid] = INVALID_V;
					Reverse_V_map[vid].clear();
				}
			}
			else if (std::get<4>(e) == Edge_tag::R) {
				if (mV_B_flag[v0_map] && mV_B_flag[v1_map] && !std::get<2>(mEs[e_map])) {
					std::get<3>(e) *= 1.2; Es_nonmanifold.push_back(e);
					continue;
				}
				std::sort(Reverse_V_nts[v0_map].begin(), Reverse_V_nts[v0_map].end());
				std::sort(Reverse_V_nts[v1_map].begin(), Reverse_V_nts[v1_map].end());
				std::vector<uint32_t> common_ts;
				std::set_intersection(Reverse_V_nts[v0_map].begin(), Reverse_V_nts[v0_map].end(),
					Reverse_V_nts[v1_map].begin(), Reverse_V_nts[v1_map].end(), std::back_inserter(common_ts));

				if (common_ts.size() != Reverse_E_nts[e_map].size()) {
					std::get<3>(e) *= 1.2;
					Es_nonmanifold.push_back(e);
					continue;
				}

				int genus_pre = -1, genus_aft = -1; bool manifoldness_aft = true;

				std::function<void(std::vector<uint32_t> &)> check_pre = [&](std::vector<uint32_t> &ts) -> void {
					int v_num = 0, e_num = 0;
					//v_num
					for (auto fid : ts) for (auto vid : mFs2D[fid]) mV_flag[vid] = true;
					for (auto fid : ts) for (auto vid : mFs2D[fid]) if (mV_flag[vid]) { v_num++; mV_flag[vid] = false; }
					//e_num
					for (auto fid : ts) for (auto eid : FEs[fid]) mE_flag[eid] = true;
					for (auto fid : ts) for (auto eid : FEs[fid]) if (mE_flag[eid]) { e_num++; mE_flag[eid] = false; }

					genus_pre = v_num + ts.size() - e_num;
				};
				std::function<bool(std::vector < std::vector<uint32_t>> &, std::vector < std::vector<uint32_t>> &)> check_aft =
					[&](std::vector < std::vector<uint32_t>> &fvs, std::vector < std::vector<uint32_t>> &fes) {
					//each polygon should be simple
					for (int j = 0; j < fvs.size(); j++) {
						std::vector<std::vector<uint32_t>> fvs_local(1), fes_local(1);
						fvs_local[0] = fvs[j]; fes_local[0] = fes[j];
						std::vector<uint32_t> pvs, pes, vs_disgard, es_disgard;
						if (!simple_polygon(fvs_local, fes_local, pvs, pes, vs_disgard, es_disgard))
							return false;
					}
					//genus
					int v_num = 0, e_num = 0;
					//v_num
					for (auto vs : fvs) for (auto vid : vs) mV_flag[vid] = true;
					for (auto vs : fvs) for (auto vid : vs) if (mV_flag[vid]) { v_num++; mV_flag[vid] = false; }

					//e_num
					std::vector<std::tuple<uint32_t, uint32_t, uint32_t, int>> temp;
					temp.reserve(fvs.size() * 3);
					std::vector < std::vector<uint32_t>> FEs_temp(fvs.size());
					for (uint32_t f = 0; f < fvs.size(); ++f) {
						for (uint32_t e = 0; e < fvs[f].size(); ++e) {
							uint32_t v0_ = fvs[f][e], v1_ = fvs[f][(e + 1) % fvs[f].size()];
							if (v0_ > v1_) std::swap(v0_, v1_);
							temp.push_back(std::make_tuple(v0_, v1_, f, e));
						}
						std::vector<uint32_t> es(fvs[f].size());
						FEs_temp[f] = es;
					}
					std::sort(temp.begin(), temp.end());
					e_num = -1;
					for (uint32_t j = 0; j < temp.size(); ++j) {
						if (j == 0 || (j != 0 && (std::get<0>(temp[j]) != std::get<0>(temp[j - 1]) ||
							std::get<1>(temp[j]) != std::get<1>(temp[j - 1])))) {
							e_num++;
						}
						FEs_temp[std::get<2>(temp[j])][std::get<3>(temp[j])] = e_num;
					}
					e_num++;
					genus_aft = v_num + fvs.size() - e_num;
					if (genus_pre != genus_aft) return false;
					//manifoldness
					std::vector<std::vector<uint32_t>> E_nts_temp(e_num);
					for (uint32_t j = 0; j < FEs_temp.size(); j++)
						for (uint32_t k = 0; k < FEs_temp[j].size(); k++)
							E_nts_temp[FEs_temp[j][k]].push_back(j);

					for (uint32_t j = 0; j < E_nts_temp.size(); j++)
						if (E_nts_temp[j].size() > 2) {
							manifoldness_aft = false; return false;
						}
					//two same polygons is not allowed
					for (auto &vs : fvs) std::sort(vs.begin(), vs.end());
					std::sort(fvs.begin(), fvs.end());
					for (uint32_t j = 1; j < fvs.size(); ++j)
						if (fvs[j - 1].size() == fvs[j].size() && std::equal(fvs[j - 1].begin(), fvs[j - 1].end(), fvs[j].begin()))
							return false;
					return true;
				};
				std::function<bool()> topology_check = [&]() -> bool {
					std::vector<uint32_t> ts = Reverse_V_nts[v0_map];
					ts.insert(ts.end(), Reverse_V_nts[v1_map].begin(), Reverse_V_nts[v1_map].end());
					std::sort(ts.begin(), ts.end());
					ts.erase(std::unique(ts.begin(), ts.end()), ts.end());
					check_pre(ts);

					std::vector<bool> ts_flag(ts.size(), true);
					std::vector < std::vector<uint32_t>> mFs2D_local(ts.size()), FEs_local(ts.size()), mFs2D_temp, FEs_temp;
					for (int j = 0; j < ts.size(); j++) { mFs2D_local[j] = mFs2D[ts[j]]; FEs_local[j] = FEs[ts[j]]; F_map[ts[j]] = j; }

					for (auto vid : Reverse_V_map[v1_map]) V_map[vid] = v0_map;
					for (auto fid : common_ts) {
						mFs2D_local[F_map[fid]].erase(std::remove(mFs2D_local[F_map[fid]].begin(), mFs2D_local[F_map[fid]].end(), v1_map), mFs2D_local[F_map[fid]].end());
						FEs_local[F_map[fid]].erase(std::remove(FEs_local[F_map[fid]].begin(), FEs_local[F_map[fid]].end(), e_map), FEs_local[F_map[fid]].end());

						if (mFs2D_local[F_map[fid]].size() == 2) ts_flag[F_map[fid]] = false;
					}

					mFs2D_temp.reserve(ts_flag.size()); FEs_temp.reserve(ts_flag.size());
					for (int j = 0; j < ts_flag.size(); j++) {
						if (ts_flag[j]) {
							std::vector<uint32_t> vs_;
							for (int k = 0; k < mFs2D_local[j].size(); k++)
								vs_.push_back(V_map[mFs2D_local[j][k]]);
							mFs2D_temp.push_back(vs_);
							FEs_temp.push_back(FEs_local[j]);
						}
						F_map[ts[j]] = -1;
					}
					mFs2D_temp.swap(mFs2D_local); FEs_temp.swap(FEs_local);
					bool pass_check = check_aft(mFs2D_local, FEs_local);

					for (auto vid : Reverse_V_map[v1_map]) V_map[vid] = v1_map;
					if (pass_check) return true;
					return false;
				};

				if (!topology_check()) {
					std::get<3>(e) *= 1.2;
					Es_nonmanifold.push_back(e);
					continue;
				}

				std::get<4>(mEs[e_map]) = Edge_tag::R;
				//merge v1_map and v0_map -> v0_map
				for (auto vid : Reverse_V_map[v1_map]) V_map[vid] = v0_map;

				Reverse_V_map[v0_map].insert(Reverse_V_map[v0_map].end(), Reverse_V_map[v1_map].begin(), Reverse_V_map[v1_map].end());
				std::vector<uint32_t>().swap(Reverse_V_map[v1_map]);
				//boundaryness of v0_map
				if (!mV_B_flag[v0_map]) for (auto vid : Reverse_V_map[v0_map]) if (mV_B_flag[vid]) { mV_B_flag[v0_map] = true; break; }
				//boundaryness of es
				for (auto fid : common_ts) {
					mFs2D[fid].erase(std::remove(mFs2D[fid].begin(), mFs2D[fid].end(), v1_map), mFs2D[fid].end());
					FEs[fid].erase(std::remove(FEs[fid].begin(), FEs[fid].end(), e_map), FEs[fid].end());

					if (mFs2D[fid].size() == 2) {
						mF_flag[fid] = false;
						//Reverse_V_nts, Reverse_E_nts
						for (short k = 0; k < 2; k++) {
							uint32_t vid = mFs2D[fid][k];
							Reverse_V_nts[vid].erase(std::remove(Reverse_V_nts[vid].begin(), Reverse_V_nts[vid].end(), fid), Reverse_V_nts[vid].end());
							uint32_t eid = FEs[fid][k];
							Reverse_E_nts[eid].erase(std::remove(Reverse_E_nts[eid].begin(), Reverse_E_nts[eid].end(), fid), Reverse_E_nts[eid].end());
						}
						//Reverse_E_nts
						uint32_t eid0 = FEs[fid][0], eid1 = FEs[fid][1];
						Reverse_E_nts[eid0].insert(Reverse_E_nts[eid0].end(), Reverse_E_nts[eid1].begin(), Reverse_E_nts[eid1].end());
						for (auto eid_temp : Reverse_E_nts[eid1])
							std::replace(FEs[eid_temp].begin(), FEs[eid_temp].end(), eid1, eid0);
						Reverse_E_nts[eid1].clear();
						//Reverse_E_map
						Reverse_E_map[eid0].insert(Reverse_E_map[eid0].end(), Reverse_E_map[eid1].begin(), Reverse_E_map[eid1].end());
						if (std::get<2>(mEs[eid0]) || std::get<2>(mEs[eid1]))
							for (auto eid : Reverse_E_map[eid0]) std::get<2>(mEs[eid]) = true;
						for (auto eid : Reverse_E_map[eid1]) E_map[eid] = eid0;
						Reverse_E_map[eid1].clear();
					}
				}
				//e_map 
				Reverse_E_nts[e_map].clear();
				for (auto eid : Reverse_E_map[e_map]) E_map[eid] = INVALID_E;
				Reverse_E_map[e_map].clear();
				//nts
				Reverse_V_nts[v0_map].insert(Reverse_V_nts[v0_map].end(), Reverse_V_nts[v1_map].begin(), Reverse_V_nts[v1_map].end());
				std::sort(Reverse_V_nts[v0_map].begin(), Reverse_V_nts[v0_map].end());
				Reverse_V_nts[v0_map].erase(std::unique(Reverse_V_nts[v0_map].begin(), Reverse_V_nts[v0_map].end()), Reverse_V_nts[v0_map].end());

				for (auto fid : Reverse_V_nts[v1_map]) std::replace(mFs2D[fid].begin(), mFs2D[fid].end(), v1_map, v0_map);
				std::vector<uint32_t>().swap(Reverse_V_nts[v1_map]);
				std::vector<uint32_t> nts_; nts_.reserve(Reverse_V_nts[v0_map].size());
				for (auto fid : Reverse_V_nts[v0_map]) if (mF_flag[fid]) nts_.push_back(fid); nts_.swap(Reverse_V_nts[v0_map]);

				//update V
				
				Vector3f newv, newn, posy; posy.setZero(); newv.setZero(); newn.setZero();
				Vector3f rosy = Vector3f::Zero(), q; q = mQ_copy2D.col(Reverse_V_map[v0_map][0]);
				for (auto rvid : Reverse_V_map[v0_map]) {
					newv += mV_copy2D.col(rvid);
					newn += mN_copy2D.col(rvid);
					posy += mO_copy2D.col(rvid);
					rosy += applyRotationKeep(q, mN_copy2D.col(Reverse_V_map[v0_map][0]), mQ_copy2D.col(rvid), mN_copy2D.col(rvid));
				}
				posy /= Reverse_V_map[v0_map].size();
				newv /= Reverse_V_map[v0_map].size();
				newV2D.col(v0_map) = newv;
				newN2D.col(v0_map) = newn.normalized();
				newQ2D.col(v0_map) = rosy.normalized();
				mV_tag.col(v0_map) = posy;

				pV_npes[v0_map].insert(pV_npes[v0_map].end(), pV_npes[v1_map].begin(), pV_npes[v1_map].end());
				std::vector<uint32_t>().swap(pV_npes[v1_map]);
				std::vector<uint32_t> nes_temp;
				for (auto eid : pV_npes[v0_map]) {
					if (E_map[eid] == INVALID_E) continue;
					nes_temp.push_back(E_map[eid]);
				}
				std::sort(nes_temp.begin(), nes_temp.end()); nes_temp.erase(std::unique(nes_temp.begin(), nes_temp.end()), nes_temp.end());
				pV_npes[v0_map] = nes_temp;
				//re-coloring
				if (re_color) {
					for (auto eid : nes_temp) {
						tuple_E e_;
						uint32_t v0_ = V_map[std::get<0>(mEs[eid])], v1_ = V_map[std::get<1>(mEs[eid])];
						const Vector3f o0 = mV_tag.col(v0_), o1 = mV_tag.col(v1_);
						Float energy = (o0 - o1).norm();
						std::get<3>(e_) = energy;
						std::get<0>(e_) = v0_; std::get<1>(e_) = v1_; std::get<2>(e_) = std::get<2>(mEs[eid]);
						std::get<5>(e_) = eid;
						
						//majority voting for eid color 
						uint32_t latest_time = 0;
						Vector3f q0, q1, n0, n1; std::vector<uint32_t> votes(4, 0);
						for (auto eo : Reverse_E_map[eid]) {
							uint32_t v0_o = std::get<0>(mEs[eo]), v1_o = std::get<1>(mEs[eo]), v0_om = V_map[v0_o], v1_om = V_map[v1_o];
							q0 = newQ2D.col(v0_o), q1 = newQ2D.col(v1_o), n0 = newN2D.col(v0_o), n1 = newN2D.col(v1_o);
							Vector3f q_next = applyRotationKeep(q0, n0, q1, n1);
							std::pair<int, Float> a_pair = assignColorWeighted2D(mV_tag.col(v0_om), q0, n0, newV2D.col(v0_o),
								mV_tag.col(v1_om), q_next, n1, newV2D.col(v1_o),
								mScale, mInvScale);
							votes[a_pair.first]++;
							std::get<3>(e_) = a_pair.second;

							E_TimeStamp[eo]++;
							if (latest_time < E_TimeStamp[eo])latest_time = E_TimeStamp[eo];
						}
						E_TimeStamp[eid] = latest_time;
						uint32_t pos = 0, num = votes[pos];
						for (uint32_t k = 1; k < 4; k++) if (votes[k] > num) { num = votes[k]; pos = k; }
						std::get<4>(e_) = pos;
						std::get<7>(e_) = E_TimeStamp[eid];
						
						if (std::get<4>(e_) == Edge_tag::B) continue;
						
						Es_nonmanifold.push_back(e_);
					}
				}
			}

			topology = true;
			break;
		}
		for (uint32_t i = 0; i < Es_nonmanifold.size(); i++)
			Es_reddash.push(Es_nonmanifold[i]);
		
		if ((!topology || Es_reddash.empty()) && triangle_Switch) {
			std::cout << "entering triangles removal procedure" << endl;
			Degenerate_triangle = false;
			Remove_Tri_Num = 1;
			while (Remove_Tri_Num) {

				for (uint32_t i = 0; i < Reverse_E_nts.size(); i++) {
					if (!Reverse_E_map[i].size()) continue;
					if ((!get<2>(mEs[i]) && Reverse_E_nts[i].size() != 2) ||
						(get<2>(mEs[i]) && Reverse_E_nts[i].size() != 1)) cout << "Non-manifold before" << endl;
				}

				Remove_Tri_Num = 0;
				std::vector<uint32_t>mF_flag_temp(mFs2D.size(), true);

				std::function<void(std::vector<uint32_t> &, Float &)> angle_2es =
					[&](std::vector<uint32_t> & es, Float & angle) -> void {
					Vector3f vec0 = mV_tag.col(V_map[std::get<0>(mEs[es[0]])]) - mV_tag.col(V_map[std::get<1>(mEs[es[0]])]);
					Vector3f vec1 = mV_tag.col(V_map[std::get<0>(mEs[es[1]])]) - mV_tag.col(V_map[std::get<1>(mEs[es[1]])]);
					if (V_map[std::get<0>(mEs[es[0]])] != V_map[std::get<0>(mEs[es[1]])])
						vec0 = -vec0;
					vec0.normalize(); vec1.normalize();
					Float dot_ = vec0.dot(vec1);
					angle = std::acos(dot_);
				};
				std::function<bool(std::tuple<double, uint32_t> &, std::tuple<double, uint32_t> &)> greater = [&](
					std::tuple<double, uint32_t> &a, std::tuple<double, uint32_t> &b)->bool {
					return std::get<0>(a)> std::get<0>(b);
				};

				for (uint32_t i = 0; i < FEs.size(); i++) {
					if (mF_flag_temp[i] && FEs[i].size() == 3) {
						Float cost = std::numeric_limits<Float>::infinity();
						std::vector<std::tuple<double, uint32_t>> candidates;

						bool all_reverse = true;
						for (uint32_t j = 0; j < 3; j++) {
							uint32_t eid = FEs[i][j];
							Vector3f normal0, normal1;
							if (Reverse_E_nts[eid].size() == 1) continue;
							else if (Reverse_E_nts[eid].size() == 2) {
								uint32_t f0 = Reverse_E_nts[eid][0], f1 = Reverse_E_nts[eid][1];
								Vector3f vec0, vec1;
								vec0 = (mV_tag.col(mFs2D[f0][0]) - mV_tag.col(mFs2D[f0][1])).normalized();
								vec1 = (mV_tag.col(mFs2D[f0][0]) - mV_tag.col(mFs2D[f0][2])).normalized();
								normal0 = (vec0.cross(vec1)).normalized();
								vec0 = (mV_tag.col(mFs2D[f1][0]) - mV_tag.col(mFs2D[f1][1])).normalized();
								vec1 = (mV_tag.col(mFs2D[f1][0]) - mV_tag.col(mFs2D[f1][2])).normalized();
								normal1 = (vec0.cross(vec1)).normalized();
							}
							else
								;

							Float n01 = normal0.dot(normal1);
							if (n01 > 0) all_reverse = false;
						}
						for (uint32_t j = 0; j < 3; j++) {
							uint32_t eid = FEs[i][j];
							Vector3f normal0, normal1;
							if (Reverse_E_nts[eid].size() == 1) continue;
							else if (Reverse_E_nts[eid].size() == 2) {
								uint32_t f0 = Reverse_E_nts[eid][0], f1 = Reverse_E_nts[eid][1];
								Vector3f vec0, vec1;
								vec0 = (mV_tag.col(mFs2D[f0][0]) - mV_tag.col(mFs2D[f0][1])).normalized();
								vec1 = (mV_tag.col(mFs2D[f0][0]) - mV_tag.col(mFs2D[f0][2])).normalized();
								normal0 = (vec0.cross(vec1)).normalized();
								vec0 = (mV_tag.col(mFs2D[f1][0]) - mV_tag.col(mFs2D[f1][1])).normalized();
								vec1 = (mV_tag.col(mFs2D[f1][0]) - mV_tag.col(mFs2D[f1][2])).normalized();
								normal1 = (vec0.cross(vec1)).normalized();
							}
							else std::cout << "non-manifold" << endl;

							std::vector<uint32_t> es_2(2);
							Float angle;
							es_2[0] = FEs[i][(j + 1) % 3];
							es_2[1] = FEs[i][(j + 2) % 3];
							angle_2es(es_2, angle);
							Float n01 = normal0.dot(normal1);
							if (n01 < 0 && !all_reverse) continue;
							double cost_ = std::abs(std::abs(angle - PAI / 2) - PAI / 2) * n01;
							candidates.push_back(std::make_tuple(cost_, FEs[i][j]));
						}
						std::sort(candidates.begin(), candidates.end(), greater);

						for (auto e_cur : candidates) {
							uint32_t id = std::get<1>(e_cur);
							std::vector<std::vector<uint32_t>> fvs, fes;
							for (auto fid : Reverse_E_nts[id]) {
								fvs.push_back(mFs2D[fid]); fes.push_back(FEs[fid]);
							}
							std::vector<uint32_t> pvs, pes, vs_disgard, es_disgard;
							if (!simple_polygon(fvs, fes, pvs, pes, vs_disgard, es_disgard)) continue;

							//two same polygons is not allowed
							vector<uint32_t> nfs, nfs_;
							for (auto eid : pes)nfs.insert(nfs.end(), Reverse_E_nts[eid].begin(), Reverse_E_nts[eid].end());
							sort(nfs.begin(), nfs.end()); nfs.erase(unique(nfs.begin(), nfs.end()), nfs.end());
							for (auto fid : nfs) {
								bool find_f = false;
								for (auto fid_ : Reverse_E_nts[id]) if (fid_ == fid)find_f = true;
								if (!find_f) nfs_.push_back(fid);
							}
							bool two_poly_same = false;
							vector<uint32_t> pvs_ = pvs; sort(pvs_.begin(), pvs_.end());
							for (auto fid : nfs_) {
								vector<uint32_t> vs_ = mFs2D[fid];
								std::sort(vs_.begin(), vs_.end());

								if (pvs_.size() == vs_.size() && std::equal(pvs_.begin(), pvs_.end(), vs_.begin())){
									two_poly_same = true;
									break;
								}
							}
							if (two_poly_same) continue;

							std::get<4>(mEs[id]) = Edge_tag::D;

							uint32_t f0 = Reverse_E_nts[id][0], f1 = Reverse_E_nts[id][1];
							mFs2D[f0] = pvs; FEs[f0] = pes; mF_flag[f1] = false; mFs2D[f1].clear(); FEs[f1].clear();
							//pvs, pes
							for (auto eid : pes) std::replace(Reverse_E_nts[eid].begin(), Reverse_E_nts[eid].end(), f1, f0);
							for (auto vid : pvs) {
								std::replace(Reverse_V_nts[vid].begin(), Reverse_V_nts[vid].end(), f1, f0);
								std::sort(Reverse_V_nts[vid].begin(), Reverse_V_nts[vid].end());
								Reverse_V_nts[vid].erase(std::unique(Reverse_V_nts[vid].begin(), Reverse_V_nts[vid].end()), Reverse_V_nts[vid].end());
							}
							//es_disgard, vs_disgard
							for (auto eid : es_disgard) {
								Reverse_E_nts[eid].clear();
								for (auto reid : Reverse_E_map[eid]) E_map[reid] = INVALID_E;
								Reverse_E_map[eid].clear();
							}
							for (auto vid : vs_disgard) {
								Reverse_V_nts[vid].clear();
								for (auto rvid : Reverse_V_map[vid]) V_map[rvid] = INVALID_V;
								Reverse_V_map[vid].clear();
							}

							mF_flag_temp[f0] = mF_flag_temp[f1] = false;
							Remove_Tri_Num++;
							break;
						}
					}
				}

				if (Remove_Tri_Num > 0) {
					std::cout << "removed triangles # " << Remove_Tri_Num << endl;
					Degenerate_triangle = true;
					entered_triangle = true;
				}
			}
			if (Degenerate_triangle) {
				once = false; topology = true;
				for (uint32_t i = 0; i < FEs.size(); i++) {
					if (FEs[i].size()) {
						tuple_E e = mEs[FEs[i][0]];
						uint32_t v0 = V_map[std::get<0>(e)], v1 = V_map[std::get<1>(e)];

						std::get<7>(e) = 0;
						std::get<4>(e) = Edge_tag::B;
						std::get<5>(e) = FEs[i][0];
						Es_reddash.push(e);
					}
				}
			}
		}
	}
	for (uint32_t i = 0; i < Reverse_E_map.size(); i++) {
		if (Reverse_E_map[i].size())
			for (uint32_t j = 0; j < Reverse_E_map[i].size(); j++)
				std::get<4>(mEs[Reverse_E_map[i][j]]) = std::get<4>(mEs[i]);
	}

	std::cout << Es_reddash.size() << " remaining" << endl;

	return entered_triangle;
}

bool MultiResolutionHierarchy::meshExtraction2D() {
	mV_copy2D = mV[0];  mO_copy2D = mO[0]; mQ_copy2D = mQ[0]; mN_copy2D = mN[0];

	mFs2D.clear();
	mFs2D.resize(mF.cols());
	for (uint32_t i = 0; i < mF.cols(); ++i) for (int j = 0; j < 3; ++j) mFs2D[i].push_back(mF(j, i));
	construct_Es_FEs();
	
	int genus_pre = -1, genus_aft = -1; bool manifoldness_pre = true, manifoldness_aft = true;
	topology_check_2D(mFs2D, FEs, genus_pre, manifoldness_pre);
	std::cout << "Input genus: " << genus_pre << "		manifoldness: " << manifoldness_pre << endl;
	if (!manifoldness_pre)
		return false;
	//=============START MESH EXTRACTION=============//
	vector<uint32_t> ledges;
	edge_tagging2D(ledges);

	mV_tag = mO[0]; newQ2D = mQ[0]; newN2D = mN[0]; newV2D = mV[0];
	//edge split
	if (ledges.size() && splitting) {
		split_long_edge2D(ledges);
	}

	std::vector<tuple_E> omEs = mEs;
	MatrixXf mV_tago;

	std::cout << "Compute target edge tagging" << endl;
	Es_reddash = std::priority_queue<tuple_E, std::vector<tuple_E>, LessThan>();
	for (uint32_t i = 0; i < mEs.size(); ++i) {
		switch (std::get<4>(mEs[i]))
		{
		case Edge_tag::R: Es_reddash.push(mEs[i]); break;
		case Edge_tag::B: break;
		case Edge_tag::D: Es_reddash.push(mEs[i]);  break;
		default: throw std::runtime_error("stuck at a invalid color!");
		}

		std::get<4>(mEs[i]) = Edge_tag::B;
	}

	cout << "loop starting for collapsing & splitting" << endl;
	int32_t f_num = 0, times = 0, equal_times = 0; bool triangle_Switch = false; 
	while (true) {
		cout << "triangle_Switch " << triangle_Switch << endl;
		tagging_collapseTri(triangle_Switch);
		
		swap_data2D();

		if (non_manifold) break;
		
		if (f_num == F_tag.size()) {
			equal_times++; 
			if (equal_times == 2) {
				if (doublets) {

					uint32_t Degenerate_edge = 1, Triangle = 1;
					while (Degenerate_edge || Triangle) {
						while (remove_doublets2D()) {
							mF_flag.resize(F_tag.size()); fill(mF_flag.begin(), mF_flag.end(), true);
							mFs2D = F_tag;
							swap_data2D();
						}
						Degenerate_edge = 0;
						Triangle = 0;

						if (triangles) {
							edge_tagging2D(ledges);
							while (Es_reddash.size()) Es_reddash.pop();
							for (uint32_t i = 0; i < mEs.size(); ++i) {
								switch (get<4>(mEs[i]))
								{
								case Edge_tag::R: Es_reddash.push(mEs[i]);  break;// break; //
								case Edge_tag::B: break;
								case Edge_tag::D: Es_reddash.push(mEs[i]);  break; //break; 

								default: throw std::runtime_error("stuck at a invalid color!");
								}
								get<4>(mEs[i]) = Edge_tag::B;
							}
							bool entered = tagging_collapseTri(true);

							swap_data2D();

							if (entered) {
								Degenerate_edge = 1;
								Triangle = 1;
							}
						}
					}
				}

				if (decomposes) {
					while (split_face2D(false));
					F_tag = mFs2D;
					FEs_tag = FEs;
				}

				break;
			}
		}
		else f_num = F_tag.size();
		
		times++;
		if (times > 10) {
			if (doublets) {

				uint32_t Degenerate_edge = 1, Triangle = 1;
				while (Degenerate_edge || Triangle) {
					while (remove_doublets2D()) {
						mF_flag.resize(F_tag.size()); fill(mF_flag.begin(), mF_flag.end(), true);
						mFs2D = F_tag;
						swap_data2D();
					}

					Degenerate_edge = 0;
					Triangle = 0;

					if (triangles) {
						edge_tagging2D(ledges);
						while (Es_reddash.size()) Es_reddash.pop();
						for (uint32_t i = 0; i < mEs.size(); ++i) {
							switch (get<4>(mEs[i]))
							{
							case Edge_tag::R: Es_reddash.push(mEs[i]);  break;// break; //
							case Edge_tag::B: break;
							case Edge_tag::D: Es_reddash.push(mEs[i]);  break; //break; 

							default: throw std::runtime_error("stuck at a invalid color!");
							}
							get<4>(mEs[i]) = Edge_tag::B;
						}
						bool entered = tagging_collapseTri(true);

						swap_data2D();

						if (entered) {
							Degenerate_edge = 1;
							Triangle = 1;
						}
					}
				}
			}

			if (decomposes) {
				while (split_face2D(false));
				F_tag = mFs2D;
				FEs_tag = FEs;
			}
			break;
		}
		edge_tagging2D(ledges);
		
		triangle_Switch = false;
		if (ledges.size() && splitting) {
			split_long_edge2D(ledges);
		}
		else if (triangles) triangle_Switch = true;
			
		if (splitting) {
			while (split_face2D(true));
		}

		while (Es_reddash.size()) Es_reddash.pop();
		for (uint32_t i = 0; i < mEs.size(); ++i) {
			switch (get<4>(mEs[i]))
			{
			case Edge_tag::R: Es_reddash.push(mEs[i]);  break;// break; //
			case Edge_tag::B: break;
			case Edge_tag::D: Es_reddash.push(mEs[i]);  break; //break; 
											
			default: throw std::runtime_error("stuck at a invalid color!");
			}
			get<4>(mEs[i]) = Edge_tag::B;
		}
	}

	//split_pentagon();
	F_tag = mFs2D;
	FEs_tag = FEs;

	if (doublets) {
		while (remove_doublets2D()) {
			mF_flag.resize(F_tag.size()); fill(mF_flag.begin(), mF_flag.end(), true);
			mFs2D = F_tag;
			swap_data2D();
		}
	}

	topology_check_2D(F_tag, FEs_tag, genus_aft, manifoldness_aft);
	std::cout << "Output genus: " << genus_aft << "		manifoldness: " << manifoldness_aft << endl;

	if (genus_pre != genus_aft || !manifoldness_aft)
		std::cout << "TOPOLOGY ISSUE ---->>>>>>>>DOUBLE CHECK!" << endl;

	reindex_2D(mV_tag, F_tag);
	composit_edges_centernodes_triangles(F_tag, mV_tag, E_final_rend, mV_final_rend, F_final_rend);

	ECs.clear();
	ECs.resize(F_tag.size(), Eigen::Vector4f::Zero());
	for (uint32_t i = 0; i < F_tag.size(); i++) {
		vector<uint32_t> vs = F_tag[i];
		for (auto vid : vs) {
			Vector4f v;
			v[0] = mV_tag(0,vid);
			v[1] = mV_tag(1, vid);
			v[2] = mV_tag(2, vid);
			v[3] = 1;
			ECs[i] += v;
		}
		ECs[i] /= vs.size();
	}

	//mesh element statistics
	uint32_t largest_polygon = 0;
	for (auto vs : F_tag) if (vs.size() > largest_polygon) largest_polygon = vs.size();
	uint32_t quad_num = 0; std::vector<uint32_t> Q_type(largest_polygon + 1, 0);
	for (auto vs : F_tag) {
		Q_type[vs.size()]++;
		if (vs.size() == 4)  quad_num++;
	}
	std::cout << "Total polygon: " << F_tag.size() << "  Quad_num: " << quad_num << "  Ratio: " << Float(quad_num) / F_tag.size() << endl;
	for (uint32_t i = 0; i < Q_type.size(); i++)
		if (Q_type[i]) std::cout << i << "_thgon: " << Q_type[i] << endl;
	
	return true;
}

void MultiResolutionHierarchy::edge_tagging2D(vector<uint32_t> &ledges) {

	ledges.clear();
	V_pes.clear(); V_pes.resize(mO_copy2D.cols());
	for (uint32_t i = 0; i < mEs.size(); i++) {
		uint32_t v0 = get<0>(mEs[i]), v1 = get<1>(mEs[i]);
		V_pes[v0].push_back(i);
		V_pes[v1].push_back(i);
	}
	for (auto &pes : V_pes) sort(pes.begin(), pes.end());
	
	for (auto &e : mEs) {
		uint32_t v0 = std::get<0>(e), v1 = std::get<1>(e);
		Vector3f q0 = mQ_copy2D.col(v0), q1 = applyRotationKeep(q0, mN_copy2D.col(v0), mQ_copy2D.col(v1), mN_copy2D.col(v1));
		tuple<short, Float, Vector2f> a_posy = posy2D_completeInfo(mO_copy2D.col(v0), q0, mN_copy2D.col(v0), mV_copy2D.col(v0),
			mO_copy2D.col(v1), q1, mN_copy2D.col(v1), mV_copy2D.col(v1), mScale, mInvScale);

		std::get<3>(e) = std::get<1>(a_posy);
		std::get<4>(e) = std::get<0>(a_posy);
	}

	for (uint32_t i = 0; i < mO_copy2D.cols(); i++) {
		vector<uint32_t> &vs = V_pvs[i];
		//orient rosy
		vector<Vector3f> qs(vs.size() + 1); Vector3f q0 = mQ_copy2D.col(i);
		qs[0] = q0;
		for (uint32_t j = 0; j < vs.size(); j++) qs[j + 1] = applyRotationKeep(q0,mN_copy2D.col(i), mQ_copy2D.col(vs[j]), mN_copy2D.col(vs[j]));

		//find long edges
		for (int32_t j = 0; j < vs.size(); j++) {
			for (int32_t k = j + 1; k < vs.size(); k++) {
				int32_t v0s[2], v1s[2], pos0[2], pos1[2];
				v0s[0] = i; v0s[1] = vs[j];
				v1s[0] = i; v1s[1] = vs[k];

				pos0[0] = 0;
				pos0[1] = find(vs.begin(), vs.end(), v0s[1]) - vs.begin() + 1;
				pos1[0] = 0;
				pos1[1] = find(vs.begin(), vs.end(), v1s[1]) - vs.begin() + 1;

				Vector3f qs0[2], qs1[2];
				qs0[0] = qs[pos0[0]];
				qs0[1] = qs[pos0[1]];
				qs1[0] = qs[pos1[0]];
				qs1[1] = qs[pos1[1]];


				tuple<int, Float, Vector2f> a_posy0 = posy2D_completeInfo(mO_copy2D.col(v0s[0]), qs0[0], mN_copy2D.col(v0s[0]), mV_copy2D.col(v0s[0]),
					mO_copy2D.col(v0s[1]), qs0[1], mN_copy2D.col(v0s[1]), mV_copy2D.col(v0s[1]), mScale, mInvScale);

				if (std::get<0>(a_posy0) != Edge_tag::B) {
					continue;
				}

				Vector2f len0 = std::get<2>(a_posy0);
				short direction0 = -1;
				for (uint32_t j = 0; j < 2; j++) if (len0[j] > 0.5) direction0 = j;

				tuple<int, Float, Vector2f> a_posy1 = posy2D_completeInfo(mO_copy2D.col(v1s[0]), qs1[0], mN_copy2D.col(v1s[0]), mV_copy2D.col(v1s[0]),
					mO_copy2D.col(v1s[1]), qs1[1], mN_copy2D.col(v1s[1]), mV_copy2D.col(v1s[1]), mScale, mInvScale);

				if (std::get<0>(a_posy1) != Edge_tag::B) {
					continue;
				}
				Vector2f len1 = get<2>(a_posy1);
				short direction1 = -1;
				for (uint32_t j = 0; j < 2; j++) if (len1[j] > 0.5) direction1 = j;

				if (direction0 == direction1) {
					if (std::round(len0[direction0] / len1[direction1]) >= 2) {
						uint32_t which_v = 1;
						//compute middle point info
						Vector3f qn = (qs0[0] + qs0[1]).normalized();
						Vector3f gn = (mO_copy2D.col(v0s[0]) + mO_copy2D.col(v0s[1])) * 0.5;
						Vector3f nn = (mN_copy2D.col(v0s[0]) + mN_copy2D.col(v0s[1])).normalized();
						Vector3f vn = (mV_copy2D.col(v0s[0]) + mV_copy2D.col(v0s[1])) * 0.5;

						tuple<int, Float, Vector2f> a_posy = posy2D_completeInfo(mO_copy2D.col(v1s[which_v]), qs1[which_v], mN_copy2D.col(v1s[which_v]), mV_copy2D.col(v1s[which_v]),
							gn, qn, nn, vn, mScale, mInvScale);

						if (std::get<0>(a_posy) != Edge_tag::R) continue;

						vector<uint32_t> sharede;
						set_intersection(V_pes[v0s[0]].begin(), V_pes[v0s[0]].end(), V_pes[v0s[1]].begin(), V_pes[v0s[1]].end(), back_inserter(sharede));
						if (!sharede.size()) {
							cout << "error" << endl; 
							system("PAUSE");
						}
						ledges.push_back(sharede[0]);
					}
					else if (std::round(len1[direction1] / len0[direction0]) >= 2) {
						uint32_t which_v = 1;
						//compute middle point info
						Vector3f qn = (qs1[0] + qs1[1]).normalized();
						Vector3f gn = (mO_copy2D.col(v1s[0]) + mO_copy2D.col(v1s[1])) * 0.5;
						Vector3f nn = (mN_copy2D.col(v1s[0]) + mN_copy2D.col(v1s[1])).normalized();
						Vector3f vn = (mV_copy2D.col(v1s[0]) + mV_copy2D.col(v1s[1])) * 0.5;

						tuple<int, Float, Vector2f> a_posy = posy2D_completeInfo(mO_copy2D.col(v0s[which_v]), qs0[which_v], mN_copy2D.col(v0s[which_v]), mV_copy2D.col(v0s[which_v]),
							gn, qn, nn, vn, mScale, mInvScale);

						if (std::get<0>(a_posy) != Edge_tag::R) continue;

						vector<uint32_t> sharede;
						set_intersection(V_pes[v1s[0]].begin(), V_pes[v1s[0]].end(), V_pes[v1s[1]].begin(), V_pes[v1s[1]].end(), back_inserter(sharede));

						if (!sharede.size()) {
							cout << "error" << endl; 
							system("PAUSE");
						}
						ledges.push_back(sharede[0]);
					}
				}
			}
		}
	}

	sort(ledges.begin(), ledges.end()); ledges.erase(unique(ledges.begin(), ledges.end()), ledges.end());
}
void reindex_2D(MatrixXf &HV, MatrixXf &HQ, MatrixXf &HN, MatrixXf &HO, std::vector<std::vector<uint32_t>> &HFv) {
	//re-index V
	std::vector<int32_t> V_flag(HV.size(), -1);
	for (auto fvs : HFv) for (auto vid : fvs) V_flag[vid] = 0;
	uint32_t v_num = 0;
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) V_flag[i] = v_num++;
	MatrixXf mV_local_(3, v_num), mQ_local_(3, v_num), mN_local_(3, v_num), mO_local_(3, v_num);
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) {
			mV_local_.col(V_flag[i]) = HV.col(i);
			mQ_local_.col(V_flag[i]) = HQ.col(i);
			mN_local_.col(V_flag[i]) = HN.col(i);
			mO_local_.col(V_flag[i]) = HO.col(i);
		}
	for (auto &fvs : HFv) for (uint32_t j = 0; j < fvs.size(); j++)
		fvs[j] = V_flag[fvs[j]];
	mV_local_.swap(HV);
	mQ_local_.swap(HQ);
	mN_local_.swap(HN);
	mO_local_.swap(HO);
}
void MultiResolutionHierarchy::swap_data2D(){
	F_tag.clear(); 
	for (uint32_t i = 0; i < mFs2D.size();i++) if (mF_flag[i]) F_tag.push_back(mFs2D[i]);
	reindex_2D(newV2D, newQ2D, newN2D, mV_tag, F_tag);
	
	mV_copy2D = newV2D;  mO_copy2D = mV_tag; mQ_copy2D = newQ2D; mN_copy2D = newN2D;

	mFs2D = F_tag;
	construct_Es_FEs();
	FEs_tag = FEs;
}
bool MultiResolutionHierarchy::split_long_edge2D(vector<uint32_t> &ledges) {

	uint32_t Nvo = mO_copy2D.cols(), Nv = Nvo + ledges.size();
	uint32_t Ne = mEs.size(), Nf = mFs2D.size();
	mV_copy2D.conservativeResize(3, Nv);
	mN_copy2D.conservativeResize(3, Nv);
	mQ_copy2D.conservativeResize(3, Nv);
	mO_copy2D.conservativeResize(3, Nv);
	V_pvs.resize(Nv);
	for (auto e_pos : ledges) {

		tuple_E e = mEs[e_pos]; uint32_t eid = e_pos;

		uint32_t v0 = std::get<0>(e), v1 = std::get<1>(e);
		Vector3f q0 = mQ_copy2D.col(v0), q1 = applyRotationKeep(q0, mN_copy2D.col(v0), mQ_copy2D.col(v1), mN_copy2D.col(v1));
		//compute middle point info
		Vector3f qn = (q0 + q1).normalized();
		Vector3f gn = (mO_copy2D.col(v0) + mO_copy2D.col(v1)) * 0.5;
		Vector3f nn = (mN_copy2D.col(v0) + mN_copy2D.col(v1)).normalized();
		Vector3f vm = (mV_copy2D.col(v0) + mV_copy2D.col(v1)) * 0.5;

		uint32_t vn = Nvo++;
		uint32_t en = Ne++;

		mV_copy2D.col(vn) = vm;
		mN_copy2D.col(vn) = nn;
		mQ_copy2D.col(vn) = qn;
		mO_copy2D.col(vn) = gn;

		replace(V_pvs[v0].begin(), V_pvs[v0].end(), v1, vn);
		replace(V_pvs[v1].begin(), V_pvs[v1].end(), v0, vn);
		V_pvs[vn].push_back(v0);
		V_pvs[vn].push_back(v1);

		tuple_E e0, e1; std::tuple<short, Float, Vector2f> a_posy;
		get<0>(e0) = v0;
		get<1>(e0) = vn;
		get<2>(e0) = get<2>(e);
		a_posy = posy2D_completeInfo(mO_copy2D.col(v0), q0, mN_copy2D.col(v0), mV_copy2D.col(v0), gn, qn, nn, vm, mScale, mInvScale);
		get<3>(e0) = get<1>(a_posy);
		get<4>(e0) = get<0>(a_posy);
		get<5>(e0) = eid;
		get<6>(e0) = 0;
		get<7>(e0) = 0;

		get<0>(e1) = v1;
		get<1>(e1) = vn;
		get<2>(e1) = get<2>(e);
		a_posy = posy2D_completeInfo(mO_copy2D.col(v1), q1, mN_copy2D.col(v1), mV_copy2D.col(v1), gn, qn, nn, vm, mScale, mInvScale);
		get<4>(e1) = get<0>(a_posy);
		get<3>(e1) = get<1>(a_posy);
		get<5>(e1) = en;
		get<6>(e1) = 0;
		get<7>(e1) = 0;


		if (Ne > mEs.size()) {
			mEs.resize(std::max(Ne, (uint32_t)mEs.size() * 2));
			E_pfs.resize(mEs.size());
		}
		mEs[eid] = e0; mEs[en] = e1;
		E_pfs[en] = E_pfs[eid];

		vector<uint32_t> fids, vnpos_, vnpos, vposs;
		for (uint32_t j = 0; j < E_pfs[eid].size(); j++) {
			auto fid = E_pfs[eid][j];
			//update mpFvs, mpFes
			vector<uint32_t> fvs;
			for (uint32_t k = 0; k < mFs2D[fid].size(); k++) {
				auto v_cur = mFs2D[fid][k], v_aft = mFs2D[fid][(k + 1) % mFs2D[fid].size()];
				fvs.push_back(v_cur);
				if ((v_cur == v0 && v_aft == v1) || (v_cur == v1 && v_aft == v0)) {
					fvs.push_back(vn);
					vnpos_.push_back(k + 1);
				}
			}
			mFs2D[fid] = fvs;
			FEs[fid].push_back(en);

			//collect v
			for (uint32_t k = 0; k < fvs.size(); k++) {
				auto vid = fvs[k];
				if (vid == vn || vid == v0 || vid == v1) continue;
				Vector3f q_test = applyRotationKeep(qn, nn, mQ_copy2D.col(vid), mN_copy2D.col(vid));
				a_posy = posy2D_completeInfo(mO_copy2D.col(vid), q_test, mN_copy2D.col(vid), mV_copy2D.col(vid), gn, qn, nn, vm, mScale, mInvScale);

				//if (vn == 52360 && vid == 45673) {
				//	cout << "52360  45673 color " << get<0>(a_posy) << endl;
				//}

				if (get<0>(a_posy) != Edge_tag::R) continue;
				vnpos.push_back(vnpos_[j]); fids.push_back(fid); vposs.push_back(k);
				break;
			}
		}
		if (!fids.size()) continue;
		//update es
		for (uint32_t j = 0; j < fids.size(); j++) {
			auto fid = fids[j];
			std::vector<uint32_t> &es = FEs[fid], es_temp, &vs = mFs2D[fid];
			auto vid = vs[vposs[j]];
			//orient es direction
			for (uint32_t k = 0; k < vs.size(); k++) {
				for (auto e : es) {
					if ((get<0>(mEs[e]) == vs[k] || get<0>(mEs[e]) == vs[(k + 1) % vs.size()]) &&
						(get<1>(mEs[e]) == vs[k] || get<1>(mEs[e]) == vs[(k + 1) % vs.size()])) {
						es_temp.push_back(e); break;
					}
				}
			}
			es = es_temp;

			vector<vector<uint32_t>> nes2(2), nvs2(2);
			int32_t start = vnpos[j], end = vposs[j];

			if (find(V_pvs[vid].begin(), V_pvs[vid].end(), vs[end]) != V_pvs[vid].end()) continue;

			en = Ne++;
			uint32_t fn = Nf++;

			tuple_E new_e;
			get<0>(new_e) = vs[start];
			get<1>(new_e) = vs[end];
			get<2>(new_e) = false;
			Vector3f q_test = applyRotationKeep(qn, nn, mQ_copy2D.col(vid), mN_copy2D.col(vid));
			a_posy = posy2D_completeInfo(mO_copy2D.col(vid), q_test, mN_copy2D.col(vid), mV_copy2D.col(vid), gn, qn, nn, vm, mScale, mInvScale);
			get<4>(new_e) = get<0>(a_posy);
			get<3>(new_e) = get<1>(a_posy);
			get<5>(new_e) = en;
			get<6>(new_e) = 0;
			get<7>(new_e) = 0;

			if (Ne > mEs.size()) {
				mEs.resize(std::max(Ne, (uint32_t)mEs.size() * 2));
				E_pfs.resize(mEs.size());
			}
			mEs[en] = new_e;


			for (uint32_t m = 0; m < 2; m++) {
				if (m == 1) std::swap(start, end);
				int32_t length = (end - start + vs.size() + 1) % vs.size();
				std::vector<uint32_t> nes(length), nvs(length);
				for (uint32_t k = 0; k < length; k++) {
					nvs[k] = vs[(start + k) % vs.size()];
					if (k + 1 == length) nes[k] = std::get<5>(new_e);
					else nes[k] = es[(start + k) % es.size()];
				}
				nes2[m] = nes; nvs2[m] = nvs;
			}	

			V_pvs[vs[start]].push_back(vs[end]);
			V_pvs[vs[end]].push_back(vs[start]);

			E_pfs[en].push_back(fid);
			E_pfs[en].push_back(fn);

			for (auto eid_ : nes2[1])
				if (eid_ != en && find(E_pfs[eid_].begin(), E_pfs[eid_].end(), fid) != E_pfs[eid_].end())
					replace(E_pfs[eid_].begin(), E_pfs[eid_].end(), fid, fn);

			if (Nf > mFs2D.size()) {
				mFs2D.resize(std::max(Nf, (uint32_t)mFs2D.size() * 2));
				FEs.resize(std::max(Nf, (uint32_t)FEs.size() * 2));
			}
			mFs2D[fid] = nvs2[0];
			mFs2D[fn] = nvs2[1];
			FEs[fid] = nes2[0];
			FEs[fn] = nes2[1];
		}
	}

	mEs.resize(Ne);
	E_pfs.resize(Ne);
	mFs2D.resize(Nf);
	FEs.resize(Nf);
	
	mV_tag = mO_copy2D;
	newQ2D = mQ_copy2D;
	newN2D = mN_copy2D;
	newV2D = mV_copy2D;
	return true;
}
bool MultiResolutionHierarchy::split_face2D(bool red_edge) {
	int32_t split_num = 0;
	uint32_t Ne = mEs.size(), Nf = mFs2D.size();

	for (uint32_t i = 0; i < mFs2D.size(); i++) {

		if (mFs2D[i].size() < 6) continue;

		std::vector<uint32_t> &es = FEs[i], es_temp, &vs = mFs2D[i];
		//orient es direction
		for (uint32_t k = 0; k < vs.size(); k++) {
			for (auto e : es) {
				if ((get<0>(mEs[e]) == vs[k] || get<0>(mEs[e]) == vs[(k + 1) % vs.size()]) &&
					(get<1>(mEs[e]) == vs[k] || get<1>(mEs[e]) == vs[(k + 1) % vs.size()])) {
					es_temp.push_back(e); break;
				}
			}
		}
		es = es_temp;
		int32_t start = -1, end = -1;
		vector<tuple<double, uint32_t, uint32_t>> es_rank;
		tuple<short, Float, Vector2f> a_posy;

		for (int32_t j = 0; j < vs.size(); j++) {
			int32_t pre = (j - 1 + vs.size()) % vs.size();
			for (int32_t k = j + 2; k < vs.size(); k++) {
				if (pre == k) continue;//share an original edge

				uint32_t v0 = vs[j], v1 = vs[k];
				Vector3f q_test = applyRotationKeep(mQ_copy2D.col(v0), mN_copy2D.col(v0), mQ_copy2D.col(v1), mN_copy2D.col(v1));
				a_posy = posy2D_completeInfo(mO_copy2D.col(v0), mQ_copy2D.col(v0), mN_copy2D.col(v0), mV_copy2D.col(v0), 
					mO_copy2D.col(v1), q_test, mN_copy2D.col(v1), mV_copy2D.col(v1), mScale, mInvScale);

				if (find(V_pvs[vs[j]].begin(), V_pvs[vs[j]].end(), vs[k]) != V_pvs[vs[j]].end()) {
					continue;
				}

				if (red_edge) {
					if (get<0>(a_posy) != Edge_tag::R) continue;
					start = j;
					end = k;
					break;
				}
				else if (!((k - j + vs.size() + 1) % vs.size() < 4 || (j - k + vs.size() + 1) % vs.size() < 4)) {
					Float cost = compute_cost_edge2D_angle(j, k, vs, mO_copy2D);
					es_rank.push_back(std::make_tuple(cost, j, k));
				}
			}
		}
		if (start == -1) {
			if (!es_rank.size()) continue;
			sort(es_rank.begin(), es_rank.end());
			start = get<1>(es_rank[0]);
			end = get<2>(es_rank[0]);
		}

		vector<vector<uint32_t>> nes2(2), nvs2(2);
		uint32_t en = Ne++;
		uint32_t fn = Nf++;

		tuple_E new_e;
		get<0>(new_e) = vs[start];
		get<1>(new_e) = vs[end];
		get<2>(new_e) = false;

		Vector3f q_test = applyRotationKeep(mQ_copy2D.col(vs[start]), mN_copy2D.col(vs[start]), mQ_copy2D.col(vs[end]), mN_copy2D.col(vs[end]));
		a_posy = posy2D_completeInfo(mO_copy2D.col(vs[start]), mQ_copy2D.col(vs[start]), mN_copy2D.col(vs[start]), mV_copy2D.col(vs[start]),
			mO_copy2D.col(vs[end]), q_test, mN_copy2D.col(vs[end]), mV_copy2D.col(vs[end]), mScale, mInvScale);
		
		get<4>(new_e) = get<0>(a_posy);
		get<3>(new_e) = get<1>(a_posy);
		get<5>(new_e) = en;
		get<6>(new_e) = 0;
		get<7>(new_e) = 0;

		if (Ne > mEs.size()) {
			mEs.resize(std::max(Ne, (uint32_t)mEs.size() * 2));
			E_pfs.resize(mEs.size());
		}
		mEs[en] = new_e;
		
		for (uint32_t m = 0; m < 2; m++) {
			if (m == 1) std::swap(start, end);
			int32_t length = (end - start + vs.size() + 1) % vs.size();
			std::vector<uint32_t> nes(length), nvs(length);
			for (uint32_t k = 0; k < length; k++) {
				nvs[k] = vs[(start + k) % vs.size()];
				if (k + 1 == length) nes[k] = std::get<5>(new_e);
				else nes[k] = es[(start + k) % es.size()];
			}
			nes2[m] = nes; nvs2[m] = nvs;
		}
		
		V_pvs[vs[start]].push_back(vs[end]);
		V_pvs[vs[end]].push_back(vs[start]);

		E_pfs[en].push_back(i);
		E_pfs[en].push_back(fn);

		for (auto eid_ : nes2[1])
			if (eid_ != en && find(E_pfs[eid_].begin(), E_pfs[eid_].end(), i) != E_pfs[eid_].end())
				replace(E_pfs[eid_].begin(), E_pfs[eid_].end(), i, fn);

		if (Nf > mFs2D.size()) {
			mFs2D.resize(std::max(Nf, (uint32_t)mFs2D.size() * 2));
			FEs.resize(std::max(Nf, (uint32_t)FEs.size() * 2));
		}
		mFs2D[i] = nvs2[0];
		mFs2D[fn] = nvs2[1];
		FEs[i] = nes2[0];
		FEs[fn] = nes2[1];

		split_num++;
	}
	mEs.resize(Ne);
	E_pfs.resize(Ne);
	mFs2D.resize(Nf);
	FEs.resize(Nf);

	return split_num;
}
bool MultiResolutionHierarchy::split_pentagon() {
	int32_t split_num = 0;
	uint32_t Ne = mEs.size(), Nf = mFs2D.size();

	std::function<void(std::vector<uint32_t> &, Float &)> angle_2es =
		[&](std::vector<uint32_t> & vs, Float & angle) -> void {
		Vector3f vec0 = mO_copy2D.col(vs[0]) - mV_tag.col(vs[1]);
		Vector3f vec1 = mO_copy2D.col(vs[2]) - mV_tag.col(vs[1]);
		vec0.normalize(); vec1.normalize();
		Float dot_ = vec0.dot(vec1);
		angle = std::acos(dot_);
	};


	for (uint32_t i = 0; i < mFs2D.size(); i++) {

		if (mFs2D[i].size() != 5) continue;

		std::vector<uint32_t> &es = FEs[i], es_temp, &vs = mFs2D[i];
		int32_t start = -1, end = -1;

		vector<tuple<Float, uint32_t>> vs_rank(vs.size());
		vector<Float> angles(vs.size(),0);
		for (uint32_t k = 0; k < vs.size(); k++) {
			vector<uint32_t> v3; 
			v3.push_back(vs[(k - 1 + vs.size()) % vs.size()]);
			v3.push_back(vs[k]);
			v3.push_back(vs[(k+1)%vs.size()]);
			get<1>(vs_rank[k]) = k;
			angle_2es(v3, get<0>(vs_rank[k]));
		}
		sort(vs_rank.begin(), vs_rank.end());
		start = get<1>(vs_rank[vs_rank.size()-1]);

		uint32_t vpre = vs[(start - 1 + vs.size()) % vs.size()];
		uint32_t vaft = vs[(start + 1) % vs.size()];
		Float elenpre = (mO_copy2D.col(vpre) - mO_copy2D.col(vs[start])).norm();
		Float elenaft = (mO_copy2D.col(vaft) - mO_copy2D.col(vs[start])).norm();

		uint32_t vprepre = vs[(start - 2 + vs.size()) % vs.size()];
		uint32_t vaftaft = vs[(start + 2) % vs.size()];
		Float eleno = (mO_copy2D.col(vprepre) - mO_copy2D.col(vaftaft)).norm();
		
		if (elenpre / eleno > elenaft / eleno) end = (start + 2 + vs.size()) % vs.size();
		else end = (start - 2 + vs.size()) % vs.size();

		vector<vector<uint32_t>> nes2(2), nvs2(2);
		uint32_t en = Ne++;
		uint32_t fn = Nf++;

		tuple_E new_e;
		get<0>(new_e) = vs[start];
		get<1>(new_e) = vs[end];
		get<2>(new_e) = false;

		Vector3f q_test = applyRotationKeep(mQ_copy2D.col(vs[start]), mN_copy2D.col(vs[start]), mQ_copy2D.col(vs[end]), mN_copy2D.col(vs[end]));
		tuple<short, Float, Vector2f> a_posy = posy2D_completeInfo(mO_copy2D.col(vs[start]), mQ_copy2D.col(vs[start]), mN_copy2D.col(vs[start]), mV_copy2D.col(vs[start]),
			mO_copy2D.col(vs[end]), q_test, mN_copy2D.col(vs[end]), mV_copy2D.col(vs[end]), mScale, mInvScale);

		get<4>(new_e) = get<0>(a_posy);
		get<3>(new_e) = get<1>(a_posy);
		get<5>(new_e) = en;
		get<6>(new_e) = 0;
		get<7>(new_e) = 0;

		if (Ne > mEs.size()) {
			mEs.resize(std::max(Ne, (uint32_t)mEs.size() * 2));
			E_pfs.resize(mEs.size());
		}
		mEs[en] = new_e;

		for (uint32_t m = 0; m < 2; m++) {
			if (m == 1) std::swap(start, end);
			int32_t length = (end - start + vs.size() + 1) % vs.size();
			std::vector<uint32_t> nes(length), nvs(length);
			for (uint32_t k = 0; k < length; k++) {
				nvs[k] = vs[(start + k) % vs.size()];
				if (k + 1 == length) nes[k] = std::get<5>(new_e);
				else nes[k] = es[(start + k) % es.size()];
			}
			nes2[m] = nes; nvs2[m] = nvs;
		}

		V_pvs[vs[start]].push_back(vs[end]);
		V_pvs[vs[end]].push_back(vs[start]);

		E_pfs[en].push_back(i);
		E_pfs[en].push_back(fn);

		for (auto eid_ : nes2[1])
			if (eid_ != en && find(E_pfs[eid_].begin(), E_pfs[eid_].end(), i) != E_pfs[eid_].end())
				replace(E_pfs[eid_].begin(), E_pfs[eid_].end(), i, fn);

		if (Nf > mFs2D.size()) {
			mFs2D.resize(std::max(Nf, (uint32_t)mFs2D.size() * 2));
			FEs.resize(std::max(Nf, (uint32_t)FEs.size() * 2));
		}
		mFs2D[i] = nvs2[0];
		mFs2D[fn] = nvs2[1];
		FEs[i] = nes2[0];
		FEs[fn] = nes2[1];

		split_num++;
	}
	mEs.resize(Ne);
	E_pfs.resize(Ne);
	mFs2D.resize(Nf);
	FEs.resize(Nf);

	return split_num;

}
bool MultiResolutionHierarchy::remove_doublets2D(){
	uint32_t n_dashed = 0;

	if (F_tag.size() == 3) return 0;//hack, avoid degeneracy of the simplest shape

	mV_flag.resize(mV_tag.size()); std::fill(mV_flag.begin(), mV_flag.end(), false);
	mE_flag.resize(mEs.size()); std::fill(mE_flag.begin(), mE_flag.end(), false);
	mF_flag.resize(F_tag.size()); std::fill(mF_flag.begin(), mF_flag.end(), true);
	std::vector<std::vector<uint32_t>> Vs_nes(mV_tag.cols()), Es_nfs(mEs.size());
	std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> Es_tuples;
	std::vector<std::vector<uint32_t>> Es_sets;

	//Es_nfs, mE_flag
	for (uint32_t i = 0; i < FEs_tag.size(); i++)
		for (uint32_t j = 0; j < FEs_tag[i].size(); j++) {
			Es_nfs[FEs_tag[i][j]].push_back(i);
			mE_flag[FEs_tag[i][j]] = true;
		}
	for (uint32_t i = 0; i < mE_flag.size(); i++)
		if (mE_flag[i]) {
			uint32_t v0 = std::get<0>(mEs[i]), v1 = std::get<1>(mEs[i]);
			Vs_nes[v0].push_back(i);
			Vs_nes[v1].push_back(i);
		}
	//Es_tuples, Es_sets
	Es_tuples.reserve(Es_nfs.size() / 2);
	for (uint32_t i = 0; i < Es_nfs.size(); i++) {
		if (Es_nfs[i].size() == 2) {
			if (Es_nfs[i][0] > Es_nfs[i][1]) std::swap(Es_nfs[i][0], Es_nfs[i][1]);
			Es_tuples.push_back(std::make_tuple(Es_nfs[i][0], Es_nfs[i][1], i));
		}
	}
	std::sort(Es_tuples.begin(), Es_tuples.end());
	std::vector<uint32_t> es;
	for (uint32_t i = 0; i < Es_tuples.size(); i++) {
		if (i == 0 || (i != 0 && (std::get<0>(Es_tuples[i]) != std::get<0>(Es_tuples[i - 1]) || std::get<1>(Es_tuples[i]) != std::get<1>(Es_tuples[i - 1])))) {
			if (i != 0 && es.size() > 1) Es_sets.push_back(es);
			es.clear(); es.reserve(1); es.push_back(std::get<2>(Es_tuples[i]));
		}
		else es.push_back(std::get<2>(Es_tuples[i]));

		if (i + 1 == Es_tuples.size() && es.size() > 1)
			Es_sets.push_back(es);
	}
	//reset mE_flag
	std::fill(mE_flag.begin(), mE_flag.end(), false);
	//remove Es_sets
	for (uint32_t i = 0; i < Es_sets.size(); i++) {
		//find source, destination: vs_vd
		std::vector<uint32_t> vs_vd; vs_vd.reserve(2);
		for (uint32_t j = 0; j < Es_sets[i].size(); j++) {
			tuple_E e = mEs[Es_sets[i][j]];
			uint32_t v0 = std::get<0>(e), v1 = std::get<1>(e);
			if (!mV_flag[v0]) mV_flag[v0] = true; else if (mV_flag[v0]) mV_flag[v0] = false;
			if (!mV_flag[v1]) mV_flag[v1] = true; else if (mV_flag[v1]) mV_flag[v1] = false;
		}
		short ind = 0;
		for (uint32_t j = 0; j < Es_sets[i].size(); j++) {
			tuple_E e = mEs[Es_sets[i][j]];
			uint32_t v0 = std::get<0>(e), v1 = std::get<1>(e);
			if (mV_flag[v0]) { vs_vd.push_back(v0); mV_flag[v0] = false; }
			if (mV_flag[v1]) { vs_vd.push_back(v1); mV_flag[v1] = false; }
		}
		if (vs_vd.size() > 2) { continue; std::cout << "two face have more than one shared segments (happens when coarse resolution at cylinder-like shapes)!" << endl; }
		//common_es
		std::sort(Vs_nes[vs_vd[0]].begin(), Vs_nes[vs_vd[0]].end());
		std::sort(Vs_nes[vs_vd[1]].begin(), Vs_nes[vs_vd[1]].end());
		std::vector<uint32_t> common_es;
		std::set_intersection(Vs_nes[vs_vd[0]].begin(), Vs_nes[vs_vd[0]].end(), Vs_nes[vs_vd[1]].begin(), Vs_nes[vs_vd[1]].end(), std::back_inserter(common_es));
		//f0, f1
		uint32_t f0 = Es_nfs[Es_sets[i][0]][0], f1 = Es_nfs[Es_sets[i][0]][1];
		//if exist edge && belong to one of the f
		if (common_es.size()) {
			uint32_t fid = -1, fid_expand = -1;
			if (std::find(FEs_tag[f0].begin(), FEs_tag[f0].end(), common_es[0]) != FEs_tag[f0].end())
			{
				fid = f0; fid_expand = f1;
			}
			else if (std::find(FEs_tag[f1].begin(), FEs_tag[f1].end(), common_es[0]) != FEs_tag[f1].end())
			{
				fid = f1; fid_expand = f0;
			}
			else { continue; }

			mF_flag[fid] = false;
			//remove vs from fid_expand
			for (uint32_t j = 0; j < Es_sets[i].size(); j++) {
				tuple_E e = mEs[Es_sets[i][j]];
				uint32_t v0 = std::get<0>(e), v1 = std::get<1>(e);
				mV_flag[v0] = true; mV_flag[v1] = true;
			}
			std::vector<uint32_t> vs_; vs_.reserve(F_tag[fid_expand].size());
			for (uint32_t j = 0; j < F_tag[fid_expand].size(); j++)
				if (!mV_flag[F_tag[fid_expand][j]] || (mV_flag[F_tag[fid_expand][j]] && (F_tag[fid_expand][j] == vs_vd[0] || F_tag[fid_expand][j] == vs_vd[1])))
					vs_.push_back(F_tag[fid_expand][j]);
			vs_.swap(F_tag[fid_expand]);
			for (uint32_t j = 0; j < Es_sets[i].size(); j++) {
				tuple_E e = mEs[Es_sets[i][j]];
				uint32_t v0 = std::get<0>(e), v1 = std::get<1>(e);
				mV_flag[v0] = false; mV_flag[v1] = false;
			}
			//remove es from fid_expand
			for (uint32_t j = 0; j < Es_sets[i].size(); j++)
				mE_flag[Es_sets[i][j]] = true;
			std::vector<uint32_t> es_; es_.reserve(FEs_tag[fid_expand].size());
			for (uint32_t j = 0; j < FEs_tag[fid_expand].size(); j++)
				if (!mE_flag[FEs_tag[fid_expand][j]])
					es_.push_back(FEs_tag[fid_expand][j]);
			es_.swap(FEs_tag[fid_expand]);
			FEs_tag[fid_expand].push_back(common_es[0]);
			//update v_nes
			for (short j = 0; j < 2; j++) {
				es_.clear(); es_.reserve(Vs_nes[vs_vd[j]].size());
				for (int k = 0; k < Vs_nes[vs_vd[j]].size(); k++)
					if (!mE_flag[Vs_nes[vs_vd[j]][k]]) es_.push_back(Vs_nes[vs_vd[j]][k]);
				es_.swap(Vs_nes[vs_vd[j]]);
			}
			for (uint32_t j = 0; j < Es_sets[i].size(); j++)
				mE_flag[Es_sets[i][j]] = false;
		}
		else {
			//remove vs from f0, f1
			for (uint32_t j = 0; j < Es_sets[i].size(); j++) {
				tuple_E e = mEs[Es_sets[i][j]];
				uint32_t v0 = std::get<0>(e), v1 = std::get<1>(e);
				mV_flag[v0] = true; mV_flag[v1] = true;
			}
			std::vector<uint32_t> vs_; vs_.reserve(F_tag[f0].size());
			for (uint32_t j = 0; j < F_tag[f0].size(); j++)
				if (!mV_flag[F_tag[f0][j]] || (mV_flag[F_tag[f0][j]] && (F_tag[f0][j] == vs_vd[0] || F_tag[f0][j] == vs_vd[1])))
					vs_.push_back(F_tag[f0][j]);
			vs_.swap(F_tag[f0]);
			vs_.clear(); vs_.reserve(F_tag[f1].size());
			for (uint32_t j = 0; j < F_tag[f1].size(); j++)
				if (!mV_flag[F_tag[f1][j]] || (mV_flag[F_tag[f1][j]] && (F_tag[f1][j] == vs_vd[0] || F_tag[f1][j] == vs_vd[1])))
					vs_.push_back(F_tag[f1][j]);
			vs_.swap(F_tag[f1]);
			for (uint32_t j = 0; j < Es_sets[i].size(); j++) {
				tuple_E e = mEs[Es_sets[i][j]];
				uint32_t v0 = std::get<0>(e), v1 = std::get<1>(e);
				mV_flag[v0] = false; mV_flag[v1] = false;
			}

			for (uint32_t j = 0; j < Es_sets[i].size(); j++)
				mE_flag[Es_sets[i][j]] = true;
			//remove es from f0, f1
			std::vector<uint32_t> es_; es_.reserve(FEs_tag[f0].size());
			for (uint32_t j = 0; j < FEs_tag[f0].size(); j++)
				if (!mE_flag[FEs_tag[f0][j]])
					es_.push_back(FEs_tag[f0][j]);
			es_.swap(FEs_tag[f0]);
			es_.clear(); es_.reserve(FEs_tag[f1].size());
			for (uint32_t j = 0; j < FEs_tag[f1].size(); j++)
				if (!mE_flag[FEs_tag[f1][j]])
					es_.push_back(FEs_tag[f1][j]);
			es_.swap(FEs_tag[f1]);
			//add an updated e to f0 f1
			std::get<0>(mEs[Es_sets[i][0]]) = vs_vd[0]; std::get<1>(mEs[Es_sets[i][0]]) = vs_vd[1];
			FEs_tag[f0].push_back(Es_sets[i][0]); FEs_tag[f1].push_back(Es_sets[i][0]);
			//update v_nes
			for (short j = 0; j < 2; j++) {
				es_.clear(); es_.reserve(Vs_nes[vs_vd[j]].size());
				for (int k = 0; k < Vs_nes[vs_vd[j]].size(); k++)
					if (!mE_flag[Vs_nes[vs_vd[j]][k]]) es_.push_back(Vs_nes[vs_vd[j]][k]);
				es_.swap(Vs_nes[vs_vd[j]]);
				Vs_nes[vs_vd[j]].push_back(Es_sets[i][0]);
			}
			for (uint32_t j = 0; j < Es_sets[i].size(); j++)
				mE_flag[Es_sets[i][j]] = false;
		}
		n_dashed++;
	}
	std::vector<std::vector<uint32_t>> F_tag_, FEs_tag_;
	F_tag_.reserve(F_tag.size()); FEs_tag_.reserve(F_tag.size());
	for (uint32_t i = 0; i < F_tag.size(); i++)
		if (mF_flag[i]) {
			F_tag_.push_back(F_tag[i]);
			FEs_tag_.push_back(FEs_tag[i]);
		}
	F_tag_.swap(F_tag); FEs_tag_.swap(FEs_tag);

	return n_dashed;
}
Float MultiResolutionHierarchy::compute_cost_edge2D_angle(int32_t v0, int32_t v1, vector<uint32_t> &vs, MatrixXf &V_) {
	int32_t v0_pre = (v0 - 1 + vs.size()) % vs.size(), v0_aft = (v0 + 1) % vs.size();

	int32_t v1_pre = (v1 - 1 + vs.size()) % vs.size(), v1_aft = (v1 + 1) % vs.size();

	MatrixXf bes_vec(3, 4); Vector3f e_vec = (V_.col(vs[v0]) - V_.col(vs[v1])).normalized();
	bes_vec.col(0) = (V_.col(vs[v0]) - V_.col(vs[v0_pre])).normalized();
	bes_vec.col(1) = (V_.col(vs[v0]) - V_.col(vs[v0_aft])).normalized();

	bes_vec.col(2) = (V_.col(vs[v1]) - V_.col(vs[v1_pre])).normalized();
	bes_vec.col(3) = (V_.col(vs[v1]) - V_.col(vs[v1_aft])).normalized();

	std::vector<double> angles(4);
	for (uint32_t k = 0; k < 4; k++) {
		if (k > 1) e_vec *= -1;
		angles[k] = std::acos(bes_vec.col(k).dot(e_vec));
		angles[k] = std::abs(angles[k] - PAI / 2);
	}
	std::sort(angles.begin(), angles.end(), std::greater<double>());

	return angles[0];
}
void MultiResolutionHierarchy::tagging_singularities_T_nodes(MatrixXf &V_tagging, vector<tuple_E> &E_tagging, vector<vector<uint32_t>> &F_tagging) {
	enum  V_type
	{
		regular = 0,
		singular,
		t_node,
		boundary,
		s_t_node,
	};

	vector<std::vector<uint32_t>> Vs_nes(V_tagging.cols()), Vs_nfs(V_tagging.cols());
	vector<bool> mV_B_flag(V_tagging.cols(), false);
	V_flag.resize(V_tagging.cols());
	fill(V_flag.begin(), V_flag.end(), V_type::regular);//boundary flag

	for (auto e : E_tagging) if (std::get<2>(e) == 1) { mV_B_flag[std::get<1>(e)] = mV_B_flag[std::get<0>(e)] = true; }
	for (uint32_t i = 0; i < E_tagging.size(); i++) {
		uint32_t v0 = std::get<0>(E_tagging[i]), v1 = std::get<1>(E_tagging[i]);
		Vs_nes[v0].push_back(i);
		Vs_nes[v1].push_back(i);
	}
	for (uint32_t f = 0; f < F_tagging.size(); ++f) for (auto vid : F_tagging[f]) Vs_nfs[vid].push_back(f);

	vector<int32_t> V_tags(V_tagging.cols(), 0);//0-regular, 1-singular, 2-t_node, 3 -boundary, and 4 both singular & t_node
	for (uint32_t i = 0; i < Vs_nfs.size(); i++) {
		if (mV_B_flag[i]) {
			V_flag[i] = V_type::boundary;
			continue;
		}

		if (Vs_nes[i].size() != 4) {
			V_flag[i] = V_type::singular;
			continue;
		}
	}
	for (uint32_t i = 0; i < F_tagging.size(); i++) {
		if (F_tagging[i].size() == 5) {
			vector<int32_t> t_candidates;
			for (uint32_t j = 0; j < F_tagging[i].size(); j++) {
				//if (Vs_nes[F_tag[i][j]].size() == 3 && !mV_B_flag[F_tag[i][j]]) t_candidates.push_back(j);
				if (mV_B_flag[F_tagging[i][j]]) {
					if (Vs_nes[F_tagging[i][j]].size() == 3 || Vs_nes[F_tagging[i][j]].size() == 2) t_candidates.push_back(j);
				}
				else {
					if (Vs_nes[F_tagging[i][j]].size() == 3) t_candidates.push_back(j);
				}
			}
			if (t_candidates.size()) {

				vector<tuple<double, uint32_t>> vs_rank;
				for (auto v_id : t_candidates) {
					int32_t v0_pre = (v_id - 1 + F_tagging[i].size()) % F_tagging[i].size(), v0_aft = (v_id + 1) % F_tagging[i].size();

					MatrixXf bes_vec(3, 2);
					bes_vec.col(0) = (V_tagging.col(F_tagging[i][v_id]) - V_tagging.col(F_tagging[i][v0_pre])).normalized();
					bes_vec.col(1) = (V_tagging.col(F_tagging[i][v_id]) - V_tagging.col(F_tagging[i][v0_aft])).normalized();

					Float dot_ = bes_vec.col(0).dot(bes_vec.col(1));
					Float angle = std::acos(dot_);
					Float cost = std::abs(angle - PAI);// (std::abs(n.sum()));
					vs_rank.push_back(std::make_tuple(cost, F_tagging[i][v_id]));
				}
				sort(vs_rank.begin(), vs_rank.end());
				
				V_flag[get<1>(vs_rank[0])] = V_type::t_node;
			}
			else {
				cout << "singular & t-node pentagon " << i << endl;
			}
		}
	}
}
void MultiResolutionHierarchy::composit_edges_colors(MatrixXf &Result_Vs, std::vector<tuple_E> &Es_to_render, MatrixXf &Result_edges)
{
	Result_edges.resize(6, 2 * Es_to_render.size());
	//for rendering edges
	for (uint32_t i = 0; i < Es_to_render.size(); ++i) {
		Vector3f color;
		if (std::get<4>(Es_to_render[i]) == Edge_tag::R)
			color = Vector3f(1, 0, 0);
		else if (std::get<4>(Es_to_render[i]) == Edge_tag::B)
			color = Vector3f(0, 0, 1);
		else if (std::get<4>(Es_to_render[i]) == Edge_tag::D)
			color = Vector3f(0, 1, 0);
		else if (std::get<4>(Es_to_render[i]) == Edge_tag::H)
		color = Vector3f(1, 1, 1);

		uint32_t i0 = std::get<0>(Es_to_render[i]), i1 = std::get<1>(Es_to_render[i]);

		Result_edges.col(i * 2 + 0) << Result_Vs.col(i0), color;
		Result_edges.col(i * 2 + 1) << Result_Vs.col(i1), color;
	}
}
void MultiResolutionHierarchy::composit_edges_centernodes_triangles(std::vector<std::vector<uint32_t>> &Actual_Fs, MatrixXf &nodes, MatrixXf &Result_edges, MatrixXf &center_nodes, MatrixXu &Triangles)
{
	//for rendering edges
	int num_e = 0;
	for (uint32_t f = 0; f < Actual_Fs.size(); ++f) {
		num_e += Actual_Fs[f].size();
	}
	Result_edges.resize(6, num_e * 2);
	Result_edges.setZero();
	num_e = 0;
	for (uint32_t f = 0; f < Actual_Fs.size(); ++f) {
		for (int k = 0; k < Actual_Fs[f].size(); ++k) {
			Vector3f color = Vector3f(0, 0, 1);
			uint32_t i0 = Actual_Fs[f][k], i1 = Actual_Fs[f][(k + 1) % Actual_Fs[f].size()];

			Result_edges.col((num_e) * 2 + 0) << nodes.col(i0), color;
			Result_edges.col((num_e++) * 2 + 1) << nodes.col(i1), color;
		}
	}
	//for rendering faces
	uint32_t col_ = 0;
	for (int i = 0; i < Actual_Fs.size(); i++)
		col_ += Actual_Fs[i].size();
	Triangles.resize(3, col_);
	Triangles.setZero();
	center_nodes.resize(3, nodes.cols() + Actual_Fs.size());
	center_nodes.setZero();
	center_nodes.block(0, 0, 3, nodes.cols()) = nodes;

	col_ = 0;
	for (int i = 0; i < Actual_Fs.size(); i++)
	{
		int f_size = Actual_Fs[i].size();
		Vector3f center;
		center.setZero();
		for (int j = 0; j < f_size; j++)
			center += nodes.col(Actual_Fs[i][j]);
		center /= Actual_Fs[i].size();
		center_nodes.col(nodes.cols() + i) = center;

		for (int j = 0; j < f_size; j++)
		{
			Triangles(0, col_) = nodes.cols() + i;
			Triangles(1, col_) = Actual_Fs[i][j];
			Triangles(2, col_++) = Actual_Fs[i][(j + 1) % f_size];
		}
	}
}