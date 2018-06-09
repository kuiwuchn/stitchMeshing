#include "hierarchy.h"
#include "positions.h"
#include "timer.h"
#include"quadratic.h"
//3D===========================================================================================================//
std::vector<std::vector<uint32_t>> mTs;

std::vector<tuple_E> mpEs;
std::vector<std::vector<uint32_t>> mpFvs, mpFes, mpPs;
std::vector<bool>mpF_boundary_flag; std::vector<std::vector<uint32_t>> PV_npfs, PE_npfs, PF_npps;

std::vector<short> mpV_flag, mpE_flag, mpF_flag, mpP_flag;

std::priority_queue<tuple_E, std::vector<tuple_E>, LessThan> Es_red;
std::vector<uint32_t> pV_map;
std::vector<std::vector<uint32_t>> Reverse_pV_map;
std::vector<uint32_t> pE_map;
std::vector<std::vector<uint32_t>> Reverse_pE_map;
std::vector<uint32_t> pF_map;
std::vector<std::vector<uint32_t>> Reverse_pF_map;


std::vector<std::vector<uint32_t>> PV_npvs;
std::vector<std::vector<uint32_t>> PV_npes_sudo;
std::vector<std::vector<uint32_t>> PV_npfs_sudo;
std::vector<std::vector<uint32_t>> PE_npfs_sudo;
std::vector<std::vector<uint32_t>> PF_npps_sudo;

//re-coloring
MatrixXf mQ_copy, mO_copy, mN_copy;
vector<Quadric> Quadric_copy, newQu3D;
MatrixXf newQ, newN3D, newV3D, newC3D;

//timeing
long long topo_check_time = 0, decomposition_time = 0, total_time = 0;
//===========================================================================================================//

void construct_Es_TetEs_Fs_TetFs_FEs()
{
	mpEs.clear(); mpFvs.clear(); mpFes.clear(); mpPs.clear(); mpF_boundary_flag.clear();
	//mpFvs, mpPs, mpF_boundary_flag
	std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>> tempF;
	tempF.reserve(mTs.size() * 4);
	mpPs.resize(mTs.size());
	for (uint32_t t = 0; t < mTs.size(); ++t) {
		for (uint32_t f = 0; f < 4; ++f) {
			uint32_t v0 = mTs[t][tet_faces[f][0]], v1 = mTs[t][tet_faces[f][1]], v2 = mTs[t][tet_faces[f][2]];
			if (v0 > v1) std::swap(v0, v1);
			if (v1 > v2) std::swap(v2, v1);
			if (v0 > v1) std::swap(v0, v1);
			tempF.push_back(std::make_tuple(v0, v1, v2, t, f));
		}
		std::vector<uint32_t> fs(4);
		mpPs[t] = fs;
	}
	std::sort(tempF.begin(), tempF.end());
	mpFvs.reserve(tempF.size() / 3); mpF_boundary_flag.reserve(tempF.size() / 3);
	int F_num = -1;
	std::vector<uint32_t> fi(3);
	for (uint32_t i = 0; i < tempF.size(); ++i) {
		if (i == 0 || (i != 0 &&
			(std::get<0>(tempF[i]) != std::get<0>(tempF[i - 1]) ||
				std::get<1>(tempF[i]) != std::get<1>(tempF[i - 1]) ||
				std::get<2>(tempF[i]) != std::get<2>(tempF[i - 1])))) {
			F_num++;
			fi[0] = std::get<0>(tempF[i]); fi[1] = std::get<1>(tempF[i]); fi[2] = std::get<2>(tempF[i]);
			mpFvs.push_back(fi); mpF_boundary_flag.push_back(true);
		}
		else if (i != 0 && (std::get<0>(tempF[i]) == std::get<0>(tempF[i - 1]) &&
			std::get<1>(tempF[i]) == std::get<1>(tempF[i - 1]) &&
			std::get<2>(tempF[i]) == std::get<2>(tempF[i - 1])))
			mpF_boundary_flag[F_num] = false;

		mpPs[std::get<3>(tempF[i])][std::get<4>(tempF[i])] = F_num;
	}
	//mpFes, mEs
	std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> temp;
	temp.reserve(mpFvs.size() * 3);
	mpFes.resize(mpFvs.size()); std::vector<uint32_t> fes(3);
	for (uint32_t i = 0; i < mpFvs.size(); ++i) {
		for (uint32_t e = 0; e < 3; ++e) {
			uint32_t v0 = mpFvs[i][e], v1 = mpFvs[i][(e + 1) % 3];
			if (v0 > v1) std::swap(v0, v1);
			temp.push_back(std::make_tuple(v0, v1, i, e));
		}
		mpFes[i] = fes;
	}
	std::sort(temp.begin(), temp.end());
	mpEs.reserve(temp.size() / 3);
	int E_num = -1;
	for (uint32_t i = 0; i < temp.size(); ++i) {
		if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
			std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
			E_num++;
			mpEs.push_back(std::make_tuple(std::get<0>(temp[i]), std::get<1>(temp[i]), false, 0, Edge_tag::B, E_num, -1, 0));
		}
		mpFes[std::get<2>(temp[i])][std::get<3>(temp[i])] = E_num;
	}
	//
	for (uint32_t i = 0; i < mpFes.size(); ++i)
			if (mpF_boundary_flag[i]) for (uint32_t e = 0; e < 3; ++e) std::get<2>(mpEs[mpFes[i][e]]) = true;


	PE_npfs.clear(); PE_npfs.resize(mpEs.size());
	PF_npps.clear(); PF_npps.resize(mpFvs.size());
	PV_npvs.clear(); PV_npvs.resize(mO_copy.cols());
	for (uint32_t i = 0; i < mpEs.size(); i++) {
		uint32_t v0 = get<0>(mpEs[i]), v1 = get<1>(mpEs[i]);
		PV_npvs[v0].push_back(v1);
		PV_npvs[v1].push_back(v0);
	}
	for (uint32_t i = 0; i < mpFes.size(); i++) for (auto eid:mpFes[i]) { PE_npfs[eid].push_back(i); }
	for (uint32_t i = 0; i < mpPs.size(); i++) for (auto fid:mpPs[i]) PF_npps[fid].push_back(i);
}
void MultiResolutionHierarchy::construct_Es_Fs_Polyhedral()
{
	mpEs.clear(); mpFes.clear();
	//mpFes, mEs
	std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> temp;
	mpFes.resize(mpFvs.size()); std::vector<uint32_t> fes;
	for (uint32_t i = 0; i < mpFvs.size(); ++i) {
		for (uint32_t e = 0; e < mpFvs[i].size(); ++e) {
			uint32_t v0 = mpFvs[i][e], v1 = mpFvs[i][(e + 1) % mpFvs[i].size()];
			if (v0 > v1) std::swap(v0, v1);
			temp.push_back(std::make_tuple(v0, v1, i, e));
		}
		fes.resize(mpFvs[i].size());
		mpFes[i] = fes;
	}
	std::sort(temp.begin(), temp.end());
	mpEs.reserve(temp.size() / 3);
	int32_t E_num = -1;
	for (uint32_t i = 0; i < temp.size(); ++i) {
		if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
			std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
			E_num++;
			mpEs.push_back(std::make_tuple(std::get<0>(temp[i]), std::get<1>(temp[i]), false, 0, Edge_tag::B, E_num, -1, 0));
		}
		mpFes[std::get<2>(temp[i])][std::get<3>(temp[i])] = E_num;
	}
	//
	mpF_boundary_flag.clear(); mpF_boundary_flag.resize(mpFes.size(),false);
	for (auto fs : mpPs)for (auto fid : fs) if(mpF_boundary_flag[fid]) mpF_boundary_flag[fid] = false; else mpF_boundary_flag[fid] = true;

	for (uint32_t i = 0; i < mpFes.size(); ++i)
		if (mpF_boundary_flag[i]) for (uint32_t e = 0; e < mpFes[i].size(); ++e) std::get<2>(mpEs[mpFes[i][e]]) = true;

	PE_npfs.clear(); PE_npfs.resize(mpEs.size());
	PF_npps.clear(); PF_npps.resize(mpFvs.size());
	PV_npvs.clear(); PV_npvs.resize(mO_copy.cols());
	for (uint32_t i = 0; i < mpEs.size(); i++) {
		uint32_t v0 = get<0>(mpEs[i]), v1 = get<1>(mpEs[i]);
		PV_npvs[v0].push_back(v1);
		PV_npvs[v1].push_back(v0);
	}
	for (uint32_t i = 0; i < mpFes.size(); i++) for (auto eid : mpFes[i]) { PE_npfs[eid].push_back(i); }
	for (uint32_t i = 0; i < mpPs.size(); i++) for (auto fid : mpPs[i]) PF_npps[fid].push_back(i);

}
bool simple_polygon_3D(std::vector<std::vector<uint32_t>> &fvs, std::vector<std::vector<uint32_t>> &fes, std::vector<uint32_t> &pvs,
	std::vector<uint32_t> &pes, std::vector<uint32_t> &vs_disgard, std::vector<uint32_t> &es_disgard, bool v_involve)
{
	es_disgard.clear();//es_disgard is in the interior
	for (int i = 0; i < fes.size(); i++) {
		for (int j = 0; j < fes[i].size(); j++)
			if (mpE_flag[fes[i][j]]) {
				es_disgard.push_back(fes[i][j]);
				mpE_flag[fes[i][j]] = false;
			}
			else mpE_flag[fes[i][j]] = true;
	}
	short which_polygon = 0;
	for (int i = 0; i < fes.size(); i++) {
		for (int j = 0; j < fes[i].size(); j++)
			if (mpE_flag[fes[i][j]]) {
				if (!pes.size())
					which_polygon = i;

				pes.push_back(fes[i][j]); mpE_flag[fes[i][j]] = false;
			}
	}
	//test nvs for each v
	for (uint32_t i = 0; i < pes.size(); ++i) {
		uint32_t v0 = pV_map[std::get<0>(mpEs[pes[i]])], v1 = pV_map[std::get<1>(mpEs[pes[i]])];
		mpV_flag[v0]++; mpV_flag[v1]++;
		if (mpV_flag[v0] > 2 || mpV_flag[v1] > 2) {
			for (uint32_t j = 0; j < pes.size(); ++j) {
				uint32_t v0_ = pV_map[std::get<0>(mpEs[pes[j]])], v1_ = pV_map[std::get<1>(mpEs[pes[j]])];
				mpV_flag[v0_] = mpV_flag[v1_] = 0;
			}
			return false;
		}
	}
	for (uint32_t i = 0; i < pes.size(); ++i) {
		uint32_t v0 = pV_map[std::get<0>(mpEs[pes[i]])], v1 = pV_map[std::get<1>(mpEs[pes[i]])];
		mpV_flag[v0] = mpV_flag[v1] = 0;
	}
	//extract the polygon	
	if (!pes.size()) return false;
	pvs.clear();
	pvs.reserve(pes.size());
	std::vector<bool> e_flag(pes.size(), false);
	pvs.push_back(pV_map[std::get<0>(mpEs[pes[0]])]);
	pvs.push_back(pV_map[std::get<1>(mpEs[pes[0]])]);
	e_flag[0] = true;
	uint32_t start_v = pvs[1];
	for (uint32_t i = 2; i < pes.size(); i++) {
		for (uint32_t j = 1; j < pes.size(); j++) {
			if (!e_flag[j]) {
				if (pV_map[std::get<0>(mpEs[pes[j]])] == start_v) {
					e_flag[j] = true;
					pvs.push_back(pV_map[std::get<1>(mpEs[pes[j]])]);
					start_v = pV_map[std::get<1>(mpEs[pes[j]])];
					break;
				}
				else if (pV_map[std::get<1>(mpEs[pes[j]])] == start_v) {
					e_flag[j] = true;
					pvs.push_back(pV_map[std::get<0>(mpEs[pes[j]])]);
					start_v = pV_map[std::get<0>(mpEs[pes[j]])];
					break;
				}
			}
		}
	}
	if (pvs.size() != pes.size())
		return false;
	if (!v_involve) return true;
	//judge direction
	bool correct = true;
	for (int i = 0; i < fvs[which_polygon].size(); i++)
		if (fvs[which_polygon][i] == pV_map[std::get<0>(mpEs[pes[0]])])
			if (fvs[which_polygon][(i + 1) % fvs[which_polygon].size()] != pV_map[std::get<1>(mpEs[pes[0]])])
				correct = false;
	if (!correct)
		std::reverse(pvs.begin(), pvs.end());
	//vs_disgard
	for (uint32_t i = 0; i < fvs.size(); ++i)
		for (uint32_t j = 0; j < fvs[i].size(); ++j)
			mpV_flag[fvs[i][j]] = 1;
	for (uint32_t i = 0; i < pvs.size(); ++i) mpV_flag[pvs[i]] = 0;
	for (uint32_t i = 0; i < fvs.size(); ++i)
		for (uint32_t j = 0; j < fvs[i].size(); ++j)
			if (mpV_flag[fvs[i][j]] == 1) {
				vs_disgard.push_back(fvs[i][j]); mpV_flag[fvs[i][j]] = 0;
			}
	return true;
}
bool simple_polygon_3D_v2(vector<vector<uint32_t>> &fvs, vector<vector<uint32_t>> &fes, vector<uint32_t> &pvs,
	vector<uint32_t> &pes, vector<uint32_t> &vs_disgard, vector<uint32_t> &es_disgard, bool v_involve){
	es_disgard.clear();//es_disgard is in the interior
	for (int i = 0; i < fes.size(); i++) {
		for (int j = 0; j < fes[i].size(); j++)
			if (mpE_flag[fes[i][j]]) {
				es_disgard.push_back(fes[i][j]);
				mpE_flag[fes[i][j]] = false;
			}
			else mpE_flag[fes[i][j]] = true;
	}
	short which_polygon = 0;
	for (int i = 0; i < fes.size(); i++) {
		for (int j = 0; j < fes[i].size(); j++)
			if (mpE_flag[fes[i][j]]) {
				if (!pes.size())
					which_polygon = i;

				pes.push_back(fes[i][j]); mpE_flag[fes[i][j]] = false;
			}
	}
	//test nvs for each v
	for (uint32_t i = 0; i < pes.size(); ++i) {
		uint32_t v0 = get<0>(mpEs[pes[i]]), v1 = get<1>(mpEs[pes[i]]);
		mpV_flag[v0]++; mpV_flag[v1]++;
		if (mpV_flag[v0] > 2 || mpV_flag[v1] > 2) {
			for (uint32_t j = 0; j < pes.size(); ++j) {
				uint32_t v0_ = get<0>(mpEs[pes[j]]), v1_ = get<1>(mpEs[pes[j]]);
				mpV_flag[v0_] = mpV_flag[v1_] = 0;
			}
			return false;
		}
	}
	for (uint32_t i = 0; i < pes.size(); ++i) {
		uint32_t v0 = get<0>(mpEs[pes[i]]), v1 = get<1>(mpEs[pes[i]]);
		mpV_flag[v0] = mpV_flag[v1] = 0;
	}
	//extract the polygon	
	if (!pes.size()) return false;
	pvs.clear();
	pvs.reserve(pes.size());
	std::vector<bool> e_flag(pes.size(), false);
	pvs.push_back(get<0>(mpEs[pes[0]]));
	pvs.push_back(get<1>(mpEs[pes[0]]));
	e_flag[0] = true;
	uint32_t start_v = pvs[1];
	for (uint32_t i = 2; i < pes.size(); i++) {
		for (uint32_t j = 1; j < pes.size(); j++) {
			if (!e_flag[j]) {
				if (get<0>(mpEs[pes[j]]) == start_v) {
					e_flag[j] = true;
					pvs.push_back(get<1>(mpEs[pes[j]]));
					start_v = get<1>(mpEs[pes[j]]);
					break;
				}
				else if (get<1>(mpEs[pes[j]]) == start_v) {
					e_flag[j] = true;
					pvs.push_back(get<0>(mpEs[pes[j]]));
					start_v = get<0>(mpEs[pes[j]]);
					break;
				}
			}
		}
	}
	if (pvs.size() != pes.size())
		return false;
	if (!v_involve) return true;
	//judge direction
	bool correct = true;
	for (int i = 0; i < fvs[which_polygon].size(); i++)
		if (fvs[which_polygon][i] == get<0>(mpEs[pes[0]]))
			if (fvs[which_polygon][(i + 1) % fvs[which_polygon].size()] != get<1>(mpEs[pes[0]]))
				correct = false;
	if (!correct)
		std::reverse(pvs.begin(), pvs.end());
	//vs_disgard
	for (uint32_t i = 0; i < fvs.size(); ++i)
		for (uint32_t j = 0; j < fvs[i].size(); ++j)
			mpV_flag[fvs[i][j]] = 1;
	for (uint32_t i = 0; i < pvs.size(); ++i) mpV_flag[pvs[i]] = 0;
	for (uint32_t i = 0; i < fvs.size(); ++i)
		for (uint32_t j = 0; j < fvs[i].size(); ++j)
			if (mpV_flag[fvs[i][j]] == 1) {
				vs_disgard.push_back(fvs[i][j]); mpV_flag[fvs[i][j]] = 0;
			}
	return true;
}
bool simple_polyhedral(std::vector<std::vector<uint32_t>> &pfs, std::vector<uint32_t> &pf,
	std::vector<uint32_t> &vs_disgard, std::vector<uint32_t> &es_disgard, std::vector<uint32_t> &fs_disgard)
{
	fs_disgard.clear();//fs_disgard is in the interior
	for (int i = 0; i < pfs.size(); i++) {
		for (int j = 0; j < pfs[i].size(); j++)
			if (mpF_flag[pfs[i][j]]) {
				fs_disgard.push_back(pfs[i][j]);
				mpF_flag[pfs[i][j]] = false;
			}
			else mpF_flag[pfs[i][j]] = true;
	}
	for (int i = 0; i < pfs.size(); i++) {
		for (int j = 0; j < pfs[i].size(); j++)
			if (mpF_flag[pfs[i][j]]) {
				pf.push_back(pfs[i][j]); mpF_flag[pfs[i][j]] = false;
			}
	}
	//test each e whether non-manifold
	bool non_simple = false;
	for (uint32_t i = 0; i < pf.size(); ++i) {
		for (uint32_t j = 0; j < mpFes[pf[i]].size(); ++j) {
			mpE_flag[mpFes[pf[i]][j]]++;
			if (mpE_flag[mpFes[pf[i]][j]] > 2) non_simple = true;
		}
		if (non_simple) {
			for (uint32_t k = 0; k < pf.size(); ++k) for (uint32_t j = 0; j < mpFes[pf[k]].size(); ++j) mpE_flag[mpFes[pf[k]][j]] = false;
			return false;
		}
	}
	for (uint32_t k = 0; k < pf.size(); ++k) for (uint32_t j = 0; j < mpFes[pf[k]].size(); ++j) if (mpE_flag[mpFes[pf[k]][j]] != 2) non_simple = true;
	for (uint32_t i = 0; i < pf.size(); ++i) for (uint32_t j = 0; j < mpFes[pf[i]].size(); ++j) mpE_flag[mpFes[pf[i]][j]] = false;
	if (non_simple) return false;
	//test each v whether non-manifold
	std::vector<uint32_t> vs_set;
	for (uint32_t i = 0; i < pf.size(); ++i) {
		for (uint32_t j = 0; j <mpFes[pf[i]].size(); ++j) {
			uint32_t v0 = pV_map[std::get<0>(mpEs[mpFes[pf[i]][j]])];
			uint32_t v1 = pV_map[std::get<1>(mpEs[mpFes[pf[i]][j]])];
			PV_npfs_sudo[v0].push_back(i);
			PV_npfs_sudo[v1].push_back(i);
			if (!mpV_flag[v0]) { vs_set.push_back(v0); mpV_flag[v0] = true; }
			if (!mpV_flag[v1]) { vs_set.push_back(v1); mpV_flag[v1] = true; }
		}
	}
	for (uint32_t i = 0; i < pf.size(); ++i) {
		for (uint32_t j = 0; j <mpFes[pf[i]].size(); ++j) {
			uint32_t v0 = pV_map[std::get<0>(mpEs[mpFes[pf[i]][j]])];
			uint32_t v1 = pV_map[std::get<1>(mpEs[mpFes[pf[i]][j]])];
			mpV_flag[v0] = mpV_flag[v1] = false;
		}
	}
	for (uint32_t m = 0; m<vs_set.size(); m++) {
		std::sort(PV_npfs_sudo[vs_set[m]].begin(), PV_npfs_sudo[vs_set[m]].end());
		PV_npfs_sudo[vs_set[m]].erase(std::unique(PV_npfs_sudo[vs_set[m]].begin(), PV_npfs_sudo[vs_set[m]].end()), PV_npfs_sudo[vs_set[m]].end());
		if (PV_npfs_sudo[vs_set[m]].size() == pf.size()) continue;

		std::vector<std::vector<uint32_t>> fes(PV_npfs_sudo[vs_set[m]].size()), fvs;
		for (uint32_t k = 0; k < fes.size(); k++) fes[k] = mpFes[pf[PV_npfs_sudo[vs_set[m]][k]]];
		std::vector<uint32_t> fv, fe, vs_dis, es_dis;
		if (!simple_polygon_3D(fvs, fes, fv, fe, vs_dis, es_dis, false))
		{
			non_simple = true; break;
		}
	}
	for (uint32_t m = 0; m<vs_set.size(); m++)
		PV_npfs_sudo[vs_set[m]].clear();
	if (non_simple) return false;
	//test single layer of polyhedral, non non-manifold v and multi-layers of polyhedral
	std::vector<uint32_t> pf_temp, pf_sudo; pf_sudo.reserve(pf.size()); pf_temp.reserve(pf.size());
	for (uint32_t i = 0; i < pf.size(); ++i) mpF_flag[pf[i]] = true;
	pf_temp.push_back(pf[0]); pf_sudo = pf_temp; mpF_flag[pf_temp[0]] = false;
	while (pf_temp.size()) {
		std::vector<uint32_t> pf_;
		for (int i = 0; i < pf_temp.size(); i++)
			for (int j = 0; j < mpFes[pf_temp[i]].size(); j++) {
				uint32_t eid = mpFes[pf_temp[i]][j];
				for (int k = 0; k<PE_npfs[eid].size(); k++)
					if (mpF_flag[PE_npfs[eid][k]]) {
						pf_.push_back(PE_npfs[eid][k]);
						mpF_flag[PE_npfs[eid][k]] = false;
					}
			}
		if (pf_.size()) {
			pf_temp = pf_;
			pf_sudo.insert(pf_sudo.end(), pf_temp.begin(), pf_temp.end());
		}
		else break;
	}

	if (pf_sudo.size() != pf.size()) {
		for (uint32_t i = 0; i < pf.size(); ++i) mpF_flag[pf[i]] = false;
		return false;
	}
	//es_disgard, vs_disgard
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			for (uint32_t k = 0; k < mpFes[pfs[i][j]].size(); k++)
				mpE_flag[mpFes[pfs[i][j]][k]] = true;
			for (uint32_t k = 0; k < mpFvs[pfs[i][j]].size(); k++)
				mpV_flag[mpFvs[pfs[i][j]][k]] = true;
		}
	}
	for (uint32_t j = 0; j < pf.size(); ++j) {
		for (uint32_t k = 0; k < mpFes[pf[j]].size(); k++)
			mpE_flag[mpFes[pf[j]][k]] = false;
		for (uint32_t k = 0; k < mpFvs[pf[j]].size(); k++)
			mpV_flag[mpFvs[pf[j]][k]] = false;
	}
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			for (uint32_t k = 0; k < mpFes[pfs[i][j]].size(); k++)
				if (mpE_flag[mpFes[pfs[i][j]][k]]) {
					mpE_flag[mpFes[pfs[i][j]][k]] = false;
					es_disgard.push_back(mpFes[pfs[i][j]][k]);
				}
			for (uint32_t k = 0; k < mpFvs[pfs[i][j]].size(); k++)
				if (mpV_flag[mpFvs[pfs[i][j]][k]]) {
					mpV_flag[mpFvs[pfs[i][j]][k]] = false;
					vs_disgard.push_back(mpFvs[pfs[i][j]][k]);
				}
		}
	}

	return true;
}
bool simple_polyhedral_v2(std::vector<std::vector<uint32_t>> &pfs)
{
	//test each e whether non-manifold
	bool non_simple = false;
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			mpE_flag[pfs[i][j]]++;
			if (mpE_flag[pfs[i][j]] > 2) non_simple = true;
		}
		if (non_simple) {
			for (uint32_t k = 0; k < pfs.size(); ++k) for (uint32_t j = 0; j < pfs[k].size(); ++j) mpE_flag[pfs[k][j]] = false;
			return false;
		}
	}
	for (uint32_t k = 0; k < pfs.size(); ++k) 
		for (uint32_t j = 0; j < pfs[k].size(); ++j)
			if (mpE_flag[pfs[k][j]] != 2) non_simple = true;
	for (uint32_t k = 0; k < pfs.size(); ++k)
		for (uint32_t j = 0; j < pfs[k].size(); ++j)
			mpE_flag[pfs[k][j]] = false;

	if (non_simple) return false;
	//test each v whether non-manifold
	std::vector<uint32_t> vs_set;
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			uint32_t v0 = pV_map[std::get<0>(mpEs[pfs[i][j]])];
			uint32_t v1 = pV_map[std::get<1>(mpEs[pfs[i][j]])];
			PV_npfs_sudo[v0].push_back(i);
			PV_npfs_sudo[v1].push_back(i);
			if (!mpV_flag[v0]) { vs_set.push_back(v0); mpV_flag[v0] = true; }
			if (!mpV_flag[v1]) { vs_set.push_back(v1); mpV_flag[v1] = true; }
		}
	}
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			uint32_t v0 = pV_map[std::get<0>(mpEs[pfs[i][j]])];
			uint32_t v1 = pV_map[std::get<1>(mpEs[pfs[i][j]])];
			mpV_flag[v0] = mpV_flag[v1] = false;
		}
	}
	for (uint32_t m = 0; m<vs_set.size(); m++){
		std::sort(PV_npfs_sudo[vs_set[m]].begin(), PV_npfs_sudo[vs_set[m]].end()); 
		PV_npfs_sudo[vs_set[m]].erase(std::unique(PV_npfs_sudo[vs_set[m]].begin(), PV_npfs_sudo[vs_set[m]].end()), PV_npfs_sudo[vs_set[m]].end());
		if (PV_npfs_sudo[vs_set[m]].size() == pfs.size()) continue;

		std::vector<std::vector<uint32_t>> fes(PV_npfs_sudo[vs_set[m]].size()), fvs;
		for (uint32_t k = 0; k < fes.size(); k++) fes[k] = pfs[PV_npfs_sudo[vs_set[m]][k]];	
		std::vector<uint32_t> fv, fe, vs_dis, es_dis;
		if (!simple_polygon_3D(fvs, fes, fv, fe, vs_dis, es_dis, false))
		{
			non_simple = true; break;
		}
	}
	for (uint32_t m = 0; m<vs_set.size(); m++)
		PV_npfs_sudo[vs_set[m]].clear();
	if (non_simple) return false;

	return true;
}
bool simple_polyhedral_v3(vector<vector<uint32_t>> &pfs)
{
	//test each e whether non-manifold
	bool non_simple = false;
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			mpE_flag[pfs[i][j]]++;
			if (mpE_flag[pfs[i][j]] > 2) non_simple = true;
		}
		if (non_simple) {
			for (uint32_t k = 0; k < pfs.size(); ++k) for (uint32_t j = 0; j < pfs[k].size(); ++j) mpE_flag[pfs[k][j]] = false;
			return false;
		}
	}
	for (uint32_t k = 0; k < pfs.size(); ++k)
		for (uint32_t j = 0; j < pfs[k].size(); ++j)
			if (mpE_flag[pfs[k][j]] != 2) non_simple = true;
	for (uint32_t k = 0; k < pfs.size(); ++k)
		for (uint32_t j = 0; j < pfs[k].size(); ++j)
			mpE_flag[pfs[k][j]] = false;

	if (non_simple) return false;
	//test each v whether non-manifold
	std::vector<uint32_t> vs_set;
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			uint32_t v0 = get<0>(mpEs[pfs[i][j]]);
			uint32_t v1 = get<1>(mpEs[pfs[i][j]]);
			PV_npfs_sudo[v0].push_back(i);
			PV_npfs_sudo[v1].push_back(i);
			if (!mpV_flag[v0]) { vs_set.push_back(v0); mpV_flag[v0] = true; }
			if (!mpV_flag[v1]) { vs_set.push_back(v1); mpV_flag[v1] = true; }
		}
	}
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			uint32_t v0 = get<0>(mpEs[pfs[i][j]]);
			uint32_t v1 = get<1>(mpEs[pfs[i][j]]);
			mpV_flag[v0] = mpV_flag[v1] = false;
		}
	}
	for (uint32_t m = 0; m<vs_set.size(); m++) {
		std::sort(PV_npfs_sudo[vs_set[m]].begin(), PV_npfs_sudo[vs_set[m]].end());
		PV_npfs_sudo[vs_set[m]].erase(std::unique(PV_npfs_sudo[vs_set[m]].begin(), PV_npfs_sudo[vs_set[m]].end()), PV_npfs_sudo[vs_set[m]].end());
		if (PV_npfs_sudo[vs_set[m]].size() == pfs.size()) continue;

		std::vector<std::vector<uint32_t>> fes(PV_npfs_sudo[vs_set[m]].size()), fvs;
		for (uint32_t k = 0; k < fes.size(); k++) fes[k] = pfs[PV_npfs_sudo[vs_set[m]][k]];
		std::vector<uint32_t> fv, fe, vs_dis, es_dis;
		if (!simple_polygon_3D_v2(fvs, fes, fv, fe, vs_dis, es_dis, false))
		{
			non_simple = true; break;
		}
	}
	for (uint32_t m = 0; m<vs_set.size(); m++)
		PV_npfs_sudo[vs_set[m]].clear();
	if (non_simple) return false;

	return true;
}
void cut_a_polyhedral(std::vector<uint32_t> &ps, std::vector<std::vector<uint32_t>> &pfs, std::vector<uint32_t> &e_circle, 
	std::vector<uint32_t> &ps0, std::vector<uint32_t> &ps1)
{
	for (uint32_t i = 0; i < pfs.size(); ++i) for (uint32_t j = 0; j < pfs[i].size(); ++j) PE_npfs_sudo[pfs[i][j]].push_back(i);
	for (uint32_t i = 0; i < ps.size(); ++i) mpF_flag[i] = true;

	std::vector<uint32_t> pf_temp; pf_temp.push_back(0); ps0 = pf_temp; mpF_flag[0] = false;
	while (pf_temp.size()) {
		std::vector<uint32_t> pf_;
		for (int i = 0; i < pf_temp.size(); i++)
			for (int j = 0; j < pfs[pf_temp[i]].size(); j++) {
				uint32_t eid = pfs[pf_temp[i]][j];
				if (std::find(e_circle.begin(), e_circle.end(), eid) != e_circle.end()) continue;
				for (int k = 0; k<PE_npfs_sudo[eid].size(); k++)
					if (mpF_flag[PE_npfs_sudo[eid][k]]) {
						pf_.push_back(PE_npfs_sudo[eid][k]);
						mpF_flag[PE_npfs_sudo[eid][k]] = false;
					}
			}
		if (pf_.size()) {
			pf_temp = pf_;
			ps0.insert(ps0.end(), pf_temp.begin(), pf_temp.end());
		}
		else break;
	}
	for (uint32_t i = 0; i < ps.size(); ++i) mpF_flag[i] = false;
	for (uint32_t i = 0; i < pfs.size(); ++i) for (uint32_t j = 0; j < pfs[i].size(); ++j) if(PE_npfs_sudo[pfs[i][j]].size()) PE_npfs_sudo[pfs[i][j]].clear();
	for (uint32_t i = 0; i < ps0.size(); ++i) ps0[i] = ps[ps0[i]];
	std::set<uint32_t> s_model(ps.begin(), ps.end());
	std::set<uint32_t> s_pattern(ps0.begin(), ps0.end());
	std::set_difference(s_model.begin(), s_model.end(), s_pattern.begin(), s_pattern.end(), std::back_inserter(ps1));
}
void reindex_3D(MatrixXf &HV, std::vector<std::vector<uint32_t>> &HFv, std::vector<std::vector<uint32_t>> &HPf){
	//re-index F
	std::vector<int32_t> F_flag(HFv.size(), -1);
	for (auto pfs : HPf) for (auto fid : pfs) F_flag[fid] = 0;
	uint32_t f_num = 0;
	for (uint32_t i = 0; i < F_flag.size(); i++)
		if (F_flag[i] != -1) F_flag[i] = f_num++;
	std::vector<std::vector<uint32_t>> mFs_local_(f_num);
	for (uint32_t i = 0; i < F_flag.size(); i++)
		if (F_flag[i] != -1) mFs_local_[F_flag[i]] = HFv[i];
	mFs_local_.swap(HFv);
	for (auto &pfs : HPf)for (uint32_t i = 0; i < pfs.size(); i++)pfs[i] = F_flag[pfs[i]];
	//re-index V
	std::vector<int32_t> V_flag(HV.size(), -1);
	for (auto fvs:HFv) for (auto vid:fvs) V_flag[vid] = 0;
	uint32_t v_num = 0;
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) V_flag[i] = v_num++;
	MatrixXf mV_local_(3, v_num);
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) mV_local_.col(V_flag[i]) = HV.col(i);
	for (auto &fvs:HFv) for (uint32_t j = 0; j < fvs.size(); j++)
			fvs[j] = V_flag[fvs[j]];
	mV_local_.swap(HV);
}
void reindex_3D(MatrixXf &HV, MatrixXf &HQ, std::vector<std::vector<uint32_t>> &HFv, std::vector<std::vector<uint32_t>> &HPf) {
	//re-index F
	std::vector<int32_t> F_flag(HFv.size(), -1);
	for (auto pfs : HPf) for (auto fid : pfs) F_flag[fid] = 0;
	uint32_t f_num = 0;
	for (uint32_t i = 0; i < F_flag.size(); i++)
		if (F_flag[i] != -1) F_flag[i] = f_num++;
	std::vector<std::vector<uint32_t>> mFs_local_(f_num);
	for (uint32_t i = 0; i < F_flag.size(); i++)
		if (F_flag[i] != -1) mFs_local_[F_flag[i]] = HFv[i];
	mFs_local_.swap(HFv);
	for (auto &pfs : HPf)for (uint32_t i = 0; i < pfs.size(); i++)pfs[i] = F_flag[pfs[i]];
	//re-index V
	std::vector<int32_t> V_flag(HV.size(), -1);
	for (auto fvs : HFv) for (auto vid : fvs) V_flag[vid] = 0;
	uint32_t v_num = 0;
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) V_flag[i] = v_num++;
	MatrixXf mV_local_(3, v_num), mQ_local_(4, v_num);
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) {
			mV_local_.col(V_flag[i]) = HV.col(i);
			mQ_local_.col(V_flag[i]) = HQ.col(i);
		}
	for (auto &fvs : HFv) for (uint32_t j = 0; j < fvs.size(); j++)
		fvs[j] = V_flag[fvs[j]];
	mV_local_.swap(HV);
	mQ_local_.swap(HQ);
}
void reindex_3D(MatrixXf &HV, MatrixXf &HQ, vector<Quadric> &HQU, std::vector<std::vector<uint32_t>> &HFv, std::vector<std::vector<uint32_t>> &HPf) {
	//re-index F
	std::vector<int32_t> F_flag(HFv.size(), -1);
	for (auto pfs : HPf) for (auto fid : pfs) F_flag[fid] = 0;
	uint32_t f_num = 0;
	for (uint32_t i = 0; i < F_flag.size(); i++)
		if (F_flag[i] != -1) F_flag[i] = f_num++;
	std::vector<std::vector<uint32_t>> mFs_local_(f_num);
	for (uint32_t i = 0; i < F_flag.size(); i++)
		if (F_flag[i] != -1) mFs_local_[F_flag[i]] = HFv[i];
	mFs_local_.swap(HFv);
	for (auto &pfs : HPf)for (uint32_t i = 0; i < pfs.size(); i++)pfs[i] = F_flag[pfs[i]];
	//re-index V
	std::vector<int32_t> V_flag(HV.size(), -1);
	for (auto fvs : HFv) for (auto vid : fvs) V_flag[vid] = 0;
	uint32_t v_num = 0;
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) V_flag[i] = v_num++;
	MatrixXf mV_local_(3, v_num), mQ_local_(4, v_num);
	vector<Quadric> HQU_local_(v_num);
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) {
			mV_local_.col(V_flag[i]) = HV.col(i);
			mQ_local_.col(V_flag[i]) = HQ.col(i);
			HQU_local_[V_flag[i]] = HQU[i];
		}
	for (auto &fvs : HFv) for (uint32_t j = 0; j < fvs.size(); j++)
		fvs[j] = V_flag[fvs[j]];
	mV_local_.swap(HV);
	mQ_local_.swap(HQ);
	HQU_local_.swap(HQU);
}
void MultiResolutionHierarchy::orient_hybrid_mesh(MatrixXf &HV, vector<vector<uint32_t>> &HF, vector<vector<uint32_t>> &HP, vector<vector<bool>> &HPF_flag) {

	mpF_boundary_flag.clear(); mpF_boundary_flag.resize(HF.size()); std::fill(mpF_boundary_flag.begin(), mpF_boundary_flag.end(), false);
	for (auto pfs : HP)for (auto fid : pfs) if (mpF_boundary_flag[fid]) mpF_boundary_flag[fid] = false; else mpF_boundary_flag[fid] = true;
	//mpFes, mEs
	mpEs.clear(); mpFes.clear();
	std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> temp;
	temp.reserve(HF.size() * 3);
	mpFes.resize(HF.size()); std::vector<uint32_t> fes;
	for (uint32_t i = 0; i < HF.size(); ++i) {
		for (uint32_t e = 0; e < HF[i].size(); ++e) {
			uint32_t v0 = HF[i][e], v1 = HF[i][(e + 1) % HF[i].size()];
			if (v0 > v1) std::swap(v0, v1);
			temp.push_back(std::make_tuple(v0, v1, i, e));
		}
		fes.resize(HF[i].size());
		mpFes[i] = fes;
	}
	std::sort(temp.begin(), temp.end());
	mpEs.reserve(temp.size() / 3);
	int32_t E_num = -1;
	for (uint32_t i = 0; i < temp.size(); ++i) {
		if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
			std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
			E_num++;
			mpEs.push_back(std::make_tuple(std::get<0>(temp[i]), std::get<1>(temp[i]), false, 0, Edge_tag::B, E_num, -1, 0));
		}
		mpFes[std::get<2>(temp[i])][std::get<3>(temp[i])] = E_num;
	}
	//PE_npfs_sudo
	PE_npfs_sudo.clear(); PE_npfs_sudo.resize(mpEs.size());
	for (uint32_t i = 0; i < mpFes.size();i++)for (auto eid : mpFes[i]) PE_npfs_sudo[eid].push_back(i);
	//PF_npps_sudo
	PF_npps_sudo.clear(); PF_npps_sudo.resize(HF.size());
	for (uint32_t i = 0; i < HP.size(); i++)for (auto fid : HP[i])PF_npps_sudo[fid].push_back(i);
	//orient surface
	mpF_flag.clear(); mpF_flag.resize(mpF_boundary_flag.size()); std::fill(mpF_flag.begin(), mpF_flag.end(), true);

	uint32_t start_f = 0;
	for (uint32_t i = 0; i < mpF_boundary_flag.size(); i++)if (mpF_boundary_flag[i]) {
		start_f = i; break;
	}
	mpF_flag[start_f] = false;
	std::queue<uint32_t> pf_temp; pf_temp.push(start_f);
	while (!pf_temp.empty()) {
		uint32_t fid = pf_temp.front(); pf_temp.pop();
		for (auto eid : mpFes[fid]) for (auto nfid : PE_npfs_sudo[eid]) {
			if (!mpF_boundary_flag[nfid] || !mpF_flag[nfid]) continue;
			uint32_t v0 = std::get<0>(mpEs[eid]), v1 = std::get<1>(mpEs[eid]);
			int32_t v0_pos = std::find(HF[fid].begin(), HF[fid].end(), v0) - HF[fid].begin();
			int32_t v1_pos = std::find(HF[fid].begin(), HF[fid].end(), v1) - HF[fid].begin();

			if ((v0_pos + 1) % HF[fid].size() != v1_pos) swap(v0, v1);

			int32_t v0_pos_ = std::find(HF[nfid].begin(), HF[nfid].end(), v0) - HF[nfid].begin();
			int32_t v1_pos_ = std::find(HF[nfid].begin(), HF[nfid].end(), v1) - HF[nfid].begin();

			if ((v0_pos_ + 1) % HF[nfid].size() == v1_pos_) std::reverse(HF[nfid].begin(), HF[nfid].end());

			pf_temp.push(nfid); mpF_flag[nfid] = false;
		}
	}

	Float res = 0;
	Vector3f ori; ori.setZero();
	for (uint32_t i = 0; i < HF.size();i++) {
		if (!mpF_boundary_flag[i]) continue;
		auto &fvs = HF[i];
		Vector3f center; center.setZero(); for (auto vid:fvs) center += HV.col(vid); center /= fvs.size();

		for (uint32_t j = 0; j < fvs.size(); j++) {
			Vector3f x = HV.col(fvs[j]) - ori, y = HV.col(fvs[(j+1)%fvs.size()]) - ori, z = center - ori;
			res += -((x[0] * y[1] * z[2] + x[1] * y[2] * z[0] + x[2] * y[0] * z[1]) - (x[2] * y[1] * z[0] + x[1] * y[0] * z[2] + x[0] * y[2] * z[1]));
		}
	}
	if (res > 0){
		for (uint32_t i = 0; i < HF.size(); i++) if (mpF_boundary_flag[i]) std::reverse(HF[i].begin(),HF[i].end());
	}
	//orient polyhedral
	std::vector<short> F_visit(mpF_boundary_flag.size(), 0);//0 un-visited, 1 visited once, 2 visited twice
	for (uint32_t i = 0; i < mpF_boundary_flag.size(); i++)if (mpF_boundary_flag[i]) { F_visit[i]++; mpF_flag[i] = true; }
	std::vector<bool> F_state(mpF_boundary_flag.size(), false);//false is the reverse direction, true is the same direction
	std::vector<bool> P_visit(HP.size(), false);
	HPF_flag.resize(HP.size());
	while (true) {
		std::vector<uint32_t> candidates;
		for (uint32_t i = 0; i < F_visit.size(); i++)if (F_visit[i] == 1)candidates.push_back(i);
		if (!candidates.size()) break;
		for (auto ca : candidates) {
			if (F_visit[ca] == 2) continue;
			uint32_t pid = PF_npps_sudo[ca][0];
			if (P_visit[pid]) if (PF_npps_sudo[ca].size() == 2) pid = PF_npps_sudo[ca][1];
			if (P_visit[pid]) {
				if (F_visit[ca] == 1) { 
					for (auto p : PF_npps_sudo[ca])
						;
					F_visit[ca] = 2; 
				}
				continue;
			}

			auto &fs = HP[pid];
			for (auto fid : fs)mpF_flag[fid] = false;
			uint32_t start_f = ca;
			mpF_flag[start_f] = true; F_visit[ca]++; if (F_state[ca]) F_state[ca] = false; else F_state[ca] = true;

			std::queue<uint32_t> pf_temp; pf_temp.push(start_f);
			while (!pf_temp.empty()) {
				uint32_t fid = pf_temp.front(); pf_temp.pop();
				for (auto eid : mpFes[fid]) for (auto nfid : PE_npfs_sudo[eid]) {
					if (mpF_flag[nfid]) continue;
					uint32_t v0 = std::get<0>(mpEs[eid]), v1 = std::get<1>(mpEs[eid]);
					int32_t v0_pos = std::find(HF[fid].begin(), HF[fid].end(), v0) - HF[fid].begin();
					int32_t v1_pos = std::find(HF[fid].begin(), HF[fid].end(), v1) - HF[fid].begin();

					if ((v0_pos + 1) % HF[fid].size() != v1_pos) std::swap(v0, v1);

					int32_t v0_pos_ = std::find(HF[nfid].begin(), HF[nfid].end(), v0) - HF[nfid].begin();
					int32_t v1_pos_ = std::find(HF[nfid].begin(), HF[nfid].end(), v1) - HF[nfid].begin();

					if (F_state[fid]) {
						if ((v0_pos_ + 1) % HF[nfid].size() == v1_pos_) F_state[nfid] = false;
						else F_state[nfid] = true;
					}
					else if (!F_state[fid]) {
						if((v0_pos_ + 1) % HF[nfid].size() == v1_pos_) F_state[nfid] = true;
						else F_state[nfid] = false;
					} 

					F_visit[nfid]++;

					pf_temp.push(nfid); mpF_flag[nfid] = true;
				}
			}
			P_visit[pid] = true;
			for (auto fid : fs) {
				if (!mpF_flag[fid]) {
					;
				}
				HPF_flag[pid].push_back(F_state[fid]);
			}
		}
	}
	for (auto &visit : F_visit) {
		if (visit != 2) { 
		} else visit = 0;
	}
	std::fill(F_state.begin(), F_state.end(), false);
	for (uint32_t i = 0; i < HP.size(); i++) for (uint32_t j = 0; j < HP[i].size(); j++) {
		switch (F_visit[HP[i][j]]) {
		case 0: F_visit[HP[i][j]]++; F_state[HP[i][j]] = HPF_flag[i][j]; break;
		case 1: F_visit[HP[i][j]]++; 
			break;
		case 2: {
		}
		}
	}
}
void MultiResolutionHierarchy::swap_data3D() {
	F_tag.clear(); P_tag.clear();
	for (auto p : mpPs) if (p.size()) P_tag.push_back(p);
	F_tag = mpFvs;
	if (Qquadric) {
		reindex_3D(mV_tag, newQ, newQu3D, F_tag, P_tag);
		mO_copy = mV_tag; mQ_copy = newQ; Quadric_copy = newQu3D;
	}
	else {
		reindex_3D(mV_tag, newQ, F_tag, P_tag);
		mO_copy = mV_tag; mQ_copy = newQ;
	}
	mpFvs = F_tag; mpPs = P_tag;
	construct_Es_Fs_Polyhedral();
}

bool MultiResolutionHierarchy::meshExtraction3D() {
	mQ_copy = mQ[0]; mO_copy = mO[0]; mN_copy = mN[0]; 
	Quadric_copy.clear(); Quadric_copy.resize(mQ_copy.cols());

	for (uint32_t i = 0; i < Quadric_copy.size(); i++) {
		Quadric_copy[i].initByPointAndNormal(mO_copy.col(i), mN_copy.col(i), quadricW);
	}

	mTs.clear();
	mTs.resize(mT.cols());
	for (uint32_t i = 0; i < mT.cols(); ++i)
		for (uint32_t j = 0; j < 4; ++j)
			mTs[i].push_back(mT(j, i));
	construct_Es_TetEs_Fs_TetFs_FEs();

//=============START MESH EXTRACTION=============//
	vector<uint32_t> ledges;
	edge_tagging3D(ledges);

	mV_tag = mO[0]; newQ = mQ[0]; newN3D = mN[0]; newQu3D = Quadric_copy;
	//vertex insert.
	if (ledges.size() && splitting) {
		split_long_edge3D(ledges);
	}
	for (uint32_t i = 0; i < mpEs.size(); ++i) {
		switch (get<4>(mpEs[i]))
		{
		case Edge_tag::R: Es_red.push(mpEs[i]);  break;
		case Edge_tag::B: break;
		case Edge_tag::D: Es_red.push(mpEs[i]);  break; 
		case Edge_tag::H: Es_red.push(mpEs[i]); break; 
		default: throw std::runtime_error("stuck at a invalid color!");
		}
		get<4>(mpEs[i]) = Edge_tag::B;
	}

	//loop
	int32_t h_num = 0, times =0, equal_times = 0;
	while (true) {
		//edge-collapse & fusion
		tagging_collapseTet();

		swap_data3D();
		
		if (times > 10) break;

		Timer<> time_decompose;
		//face insert
		if (splitting) {
			while (split_polyhedral3D()) {
				split_face3D(false);
			}
		}
		
		edge_tagging3D(ledges);
		//vertex insert
		if (ledges.size() && splitting) {
			split_long_edge3D(ledges);
			split_polyhedral3D();
		}
		//edge insert
		if (splitting) {
			while (split_face3D(true)) {
				split_polyhedral3D();
			}
		}

		if (splitting) {
			while (split_face3D(false)) {
				split_polyhedral3D();
			}
		}

		while (Es_red.size()) Es_red.pop();
		for (uint32_t i = 0; i < mpEs.size(); ++i) {
			switch (std::get<4>(mpEs[i]))
			{
			case Edge_tag::R: Es_red.push(mpEs[i]);  break;
			case Edge_tag::B: break;
			case Edge_tag::D: Es_red.push(mpEs[i]);  break;
			case Edge_tag::H: Es_red.push(mpEs[i]); break;
			default: throw std::runtime_error("stuck at a invalid color!");
			}
			std::get<4>(mpEs[i]) = Edge_tag::B;
		}

		if (h_num == P_tag.size()) {
			equal_times++; if (equal_times == 2) {
				swap_data3D();
				break;
			}
		}
		else h_num = P_tag.size();
		times++;
	}


	uint32_t largest_polyhdral = 0;
	for (uint32_t i = 0; i < P_tag.size(); i++)
		if (P_tag[i].size() > largest_polyhdral) largest_polyhdral = P_tag[i].size();

	uint32_t hex_num = 0; std::vector<uint32_t> H_type(largest_polyhdral + 1, 0), F_type(F_tag.size(), 0);
	Hex_flag.clear(); Hex_flag.resize(P_tag.size(), false);
	ECs.clear();
	ECs.resize(P_tag.size(), Eigen::Vector4f::Zero());
	for (uint32_t i = 0; i < P_tag.size(); i++) {
		H_type[P_tag[i].size()]++;
		//find hex
		if (P_tag[i].size() == 6) {
			bool hex = true;
			for (uint32_t j = 0; j < 6; j++)
				if (F_tag[P_tag[i][j]].size() != 4)
					hex = false;
			if (hex) {
				for (uint32_t j = 0; j < 6; j++)
					F_type[P_tag[i][j]] = 1;
				hex_num++;
				Hex_flag[i] = true;
			}
		}
		vector<uint32_t> vs;
		for (auto f : P_tag[i])vs.insert(vs.end(), F_tag[f].begin(), F_tag[f].end());
		sort(vs.begin(), vs.end()); vs.erase(std::unique(vs.begin(), vs.end()), vs.end());
		for (auto vid : vs) {
			Vector4f v;
			v[0] = mV_tag(0, vid);
			v[1] = mV_tag(1, vid);
			v[2] = mV_tag(2, vid);
			v[3] = 1;

			ECs[i] += v;
		}
		ECs[i] /= vs.size();
	}
	//statistics
	sta.hex_ratio = Float(hex_num) / P_tag.size();
	sta.hN = hex_num;
	sta.pN = P_tag.size();
	sta.polyhedral_nums = H_type;
	if (sta.pN > 0) {
		for (uint32_t i = 0; i < H_type.size(); i++) {
			sta.polyhedral_ratios.push_back(Float(H_type[i]) / sta.pN);
		}
	}
	//output
	PF_flag.clear();
	orient_hybrid_mesh(mV_tag, F_tag, P_tag, PF_flag);

	E_final_rend.setZero();
	E_final_rend.resize(6, 2 * mpEs.size());
	composit_edges_colors(mV_tag, mpEs, E_final_rend);
	composit_edges_centernodes_triangles(F_tag, mV_tag, E_final_rend, mV_final_rend, F_final_rend);



	cout << "done with extraction!" << endl;
}

bool MultiResolutionHierarchy::edge_tagging3D(vector<uint32_t> &ledges) {
	
	bool hyperlong_edge = false; ledges.clear();
	for (auto &e : mpEs) {
		uint32_t v0 = std::get<0>(e), v1 = std::get<1>(e);
		Quaternion q0 = mQ_copy.col(v0), q1 = Quaternion::applyRotation(mQ_copy.col(v1), q0);
		std::tuple<short, Float, Vector3f> a_posy = posy3D_completeInfo(mO_copy.col(v0), q0, mO_copy.col(v1), q1, mScale, mInvScale);

		Vector3f len = std::get<2>(a_posy);
		for (uint32_t j = 0; j < 3; j++) if (std::round(len[j]) > 1) { 
			hyperlong_edge = true; 
			if(get<0>(a_posy)==Edge_tag::B)
				ledges.push_back(std::get<5>(e)); 
			break; 
		}

		std::get<4>(e) = std::get<0>(a_posy);
		std::get<3>(e) = std::get<1>(a_posy);
	}
	PV_npes_sudo.clear(); PV_npes_sudo.resize(mO_copy.cols());
	for (uint32_t i = 0; i < mpEs.size(); i++) {
		uint32_t v0 = get<0>(mpEs[i]), v1 = get<1>(mpEs[i]);
		PV_npes_sudo[v0].push_back(i);
		PV_npes_sudo[v1].push_back(i);
	}
	for (auto &pes : PV_npes_sudo) sort(pes.begin(), pes.end());
	for (uint32_t i = 0; i < mO_copy.cols(); i++) {
		vector<uint32_t> vs = PV_npvs[i];
		//orient rosy
		std::vector<Quaternion> qs(vs.size() + 1); Quaternion q0 = mQ_copy.col(i);
		qs[0] = q0;
		for (uint32_t j = 0; j < vs.size(); j++) qs[j + 1] = Quaternion::applyRotation(mQ_copy.col(vs[j]), q0);
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

				Quaternion qs0[2], qs1[2];
				qs0[0] = qs[pos0[0]];
				qs0[1] = qs[pos0[1]];
				qs1[0] = qs[pos1[0]];
				qs1[1] = qs[pos1[1]];
				
				tuple<short, Float, Vector3f> a_posy0 = posy3D_completeInfo(mO_copy.col(v0s[0]), qs0[0], mO_copy.col(v0s[1]), qs0[1], mScale, mInvScale);
				if (std::get<0>(a_posy0) != Edge_tag::B) {
					continue;
				}

				Vector3f len0 = std::get<2>(a_posy0);
				short direction0 = -1;
				for (uint32_t j = 0; j < 3; j++) if (len0[j] > 0.5) direction0 = j;

				tuple<short, Float, Vector3f> a_posy1 = posy3D_completeInfo(mO_copy.col(v1s[0]), qs1[0], mO_copy.col(v1s[1]), qs1[1], mScale, mInvScale);
				if (std::get<0>(a_posy1) != Edge_tag::B) {
					continue;
				}
				Vector3f len1 = get<2>(a_posy1);
				short direction1 = -1;
				for (uint32_t j = 0; j < 3; j++) if (len1[j] > 0.5) direction1 = j;

				if (direction0 == direction1) {
					if (std::round(len0[direction0] / len1[direction1]) >= 2) {
						uint32_t which_v = 1;
						//compute middle point info
						Quaternion qn = (qs0[0] + qs0[1]).normalized();
						Vector3f gn = (mO_copy.col(v0s[0]) + mO_copy.col(v0s[1])) * 0.5;
						tuple<short, Float, Vector3f> a_posy = posy3D_completeInfo(mO_copy.col(v1s[which_v]), qs1[which_v], gn, qn, mScale, mInvScale);

						if (std::get<0>(a_posy) != Edge_tag::R) continue;

						vector<uint32_t> sharede;
						set_intersection(PV_npes_sudo[v0s[0]].begin(), PV_npes_sudo[v0s[0]].end(), PV_npes_sudo[v0s[1]].begin(), PV_npes_sudo[v0s[1]].end(), back_inserter(sharede));
						if (!sharede.size()) {
							system("PAUSE");
						}
						ledges.push_back(sharede[0]);
					}
					else if (std::round(len1[direction1] / len0[direction0]) >= 2) {
						uint32_t which_v = 1;
						//compute middle point info
						Quaternion qn = (qs1[0] + qs1[1]).normalized();
						Vector3f gn = (mO_copy.col(v1s[0]) + mO_copy.col(v1s[1])) * 0.5;
						tuple<short, Float, Vector3f> a_posy = posy3D_completeInfo(mO_copy.col(v0s[which_v]), qs0[which_v], gn, qn, mScale, mInvScale);
						if (std::get<0>(a_posy) != Edge_tag::R) continue;

						vector<uint32_t> sharede;
						set_intersection(PV_npes_sudo[v1s[0]].begin(), PV_npes_sudo[v1s[0]].end(), PV_npes_sudo[v1s[1]].begin(), PV_npes_sudo[v1s[1]].end(), back_inserter(sharede));
						if (!sharede.size()) {
							system("PAUSE");
						}
						ledges.push_back(sharede[0]);
					}
				}
			}
		}
	}
	for (uint32_t i = 0; i < mpFes.size(); i++) {
		auto &fes = mpFes[i];
		auto &fvs = mpFvs[i];
		std::vector<uint32_t> es_temp;
		//orient es direction
		for (uint32_t j = 0; j < fvs.size(); j++) {
			for (auto e : fes) {
				if ((get<0>(mpEs[e]) == fvs[j] || get<0>(mpEs[e]) == fvs[(j + 1) % fvs.size()]) &&
					(get<1>(mpEs[e]) == fvs[j] || get<1>(mpEs[e]) == fvs[(j + 1) % fvs.size()])) {
					es_temp.push_back(e); break;
				}
			}
		}
		fes = es_temp;
		//orient rosy
		std::vector<Quaternion> qs(fvs.size()); Quaternion q0 = mQ_copy.col(fvs[0]);
		qs[0] = q0;
		for (uint32_t j = 1; j < fvs.size(); j++) qs[j] = Quaternion::applyRotation(mQ_copy.col(fvs[j]), q0);
		//find long edges
		for (uint32_t j = 0; j < fes.size(); j++) {
			uint32_t e0 = fes[j], e1 = fes[(j + 1) % fes.size()];

			int32_t v0s[2], v1s[2], pos0[2], pos1[2];
			v0s[0] = get<0>(mpEs[e0]); v0s[1] = get<1>(mpEs[e0]);
			v1s[0] = get<0>(mpEs[e1]); v1s[1] = get<1>(mpEs[e1]);

			pos0[0] = find(fvs.begin(), fvs.end(), v0s[0]) - fvs.begin();
			pos0[1] = find(fvs.begin(), fvs.end(), v0s[1]) - fvs.begin();
			pos1[0] = find(fvs.begin(), fvs.end(), v1s[0]) - fvs.begin();
			pos1[1] = find(fvs.begin(), fvs.end(), v1s[1]) - fvs.begin();

			Quaternion qs0[2], qs1[2];
			qs0[0] = qs[pos0[0]];
			qs0[1] = qs[pos0[1]];
			qs1[0] = qs[pos1[0]];
			qs1[1] = qs[pos1[1]];
			tuple<short, Float, Vector3f> a_posy0 = posy3D_completeInfo(mO_copy.col(v0s[0]), qs0[0], mO_copy.col(v0s[1]), qs0[1], mScale, mInvScale);
			if (std::get<0>(a_posy0) != Edge_tag::B) continue;

			Vector3f len0 = std::get<2>(a_posy0);
			short direction0 = -1;
			for (uint32_t j = 0; j < 3; j++) if (len0[j] > 0.5) direction0 = j;

			tuple<short, Float, Vector3f> a_posy1 = posy3D_completeInfo(mO_copy.col(v1s[0]), qs1[0], mO_copy.col(v1s[1]), qs1[1], mScale, mInvScale);
			if (std::get<0>(a_posy1) != Edge_tag::B) continue;

			Vector3f len1 = get<2>(a_posy1);
			short direction1 = -1;
			for (uint32_t j = 0; j < 3; j++) if (len1[j] > 0.5) direction1 = j;

			if (direction0 == direction1) {
				if (std::round(len0[direction0] / len1[direction1]) >= 2) {
					uint32_t which_v = 0;
					if (v1s[1] != v0s[0] && v1s[1] != v0s[1]) which_v = 1;
					//compute middle point info
					Quaternion qn = (qs0[0] + qs0[1]).normalized();
					Vector3f gn = (mO_copy.col(v0s[0]) + mO_copy.col(v0s[1])) * 0.5;
					tuple<short, Float, Vector3f> a_posy = posy3D_completeInfo(mO_copy.col(v1s[which_v]), qs1[which_v], gn, qn, mScale, mInvScale);
					if (std::get<0>(a_posy) != Edge_tag::R) continue;

					ledges.push_back(e0);
				}
				else if (std::round(len1[direction1] / len0[direction0]) >= 2) {
					uint32_t which_v = 0;
					if (v0s[1] != v1s[0] && v0s[1] != v1s[1]) which_v = 1;
					//compute middle point info
					Quaternion qn = (qs1[0] + qs1[1]).normalized();
					Vector3f gn = (mO_copy.col(v1s[0]) + mO_copy.col(v1s[1])) * 0.5;
					tuple<short, Float, Vector3f> a_posy = posy3D_completeInfo(mO_copy.col(v0s[which_v]), qs0[which_v], gn, qn, mScale, mInvScale);
					if (std::get<0>(a_posy) != Edge_tag::R) continue;

					ledges.push_back(e1);
				}
			}
		}
	}

	sort(ledges.begin(), ledges.end()); ledges.erase(unique(ledges.begin(), ledges.end()), ledges.end());
	return hyperlong_edge;
}
void MultiResolutionHierarchy::init_edge_tagging3D()
{
	mTs.clear(); mO_copy = mO[0];
	mTs.resize(mT.cols());
	for (uint32_t i = 0; i < mT.cols(); ++i)
		for (uint32_t j = 0; j < 4; ++j)
			mTs[i].push_back(mT(j, i));
	construct_Es_TetEs_Fs_TetFs_FEs();

	std::cout << "Compute target edge tagging" << endl;
	for (uint32_t i = 0; i < mpEs.size(); ++i) {
		uint32_t v0 = std::get<0>(mpEs[i]), v1 = std::get<1>(mpEs[i]);
		Quaternion q_next = Quaternion::applyRotation(mQ[0].col(v1), mQ[0].col(v0));

		 std::pair<int, Float> a_pair = assignColorWeighted3D(mO[0].col(v0), mQ[0].col(v0), mO[0].col(v1), q_next, mScale, mInvScale);
		 std::get<4>(mpEs[i]) = a_pair.first;
		 std::get<3>(mpEs[i]) = a_pair.second;

		const Vector3f o0 = mO[0].col(v0), o1 = mO[0].col(v1);
		Float energy = (o0 - o1).norm();
	}
	for (uint32_t i = 0; i < mpFes.size(); ++i) {
		short n_R = 0;
		for (uint32_t j = 0; j < 3; j++) if (std::get<4>(mpEs[mpFes[i][j]]) == Edge_tag::R) n_R++;
		if (n_R == 2)
			for (uint32_t j = 0; j < 3; j++) std::get<4>(mpEs[mpFes[i][j]]) = Edge_tag::R;
	}

	E_rend.resize(6, mpEs.size() * 2);
	composit_edges_colors(mV[0], mpEs, E_rend);

	std::vector<std::vector<uint32_t>> PV_npes(mV[0].cols());
	for (uint32_t i = 0; i < mpEs.size(); i++) {
		uint32_t v0 = std::get<0>(mpEs[i]), v1 = std::get<1>(mpEs[i]);
		PV_npes[v0].push_back(i); PV_npes[v1].push_back(i);
	}
	std::vector<bool> V_flag_(mV[0].cols(),false);
	vector<vector<uint32_t>> VSets;
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
					uint32_t v0 = std::get<0>(mpEs[eid]), v1 = std::get<1>(mpEs[eid]);
					if (std::get<4>(mpEs[eid]) == Edge_tag::R) {
						if (!V_flag_[v0]) v_pool_sudo.push_back(v0);
						if (!V_flag_[v1]) v_pool_sudo.push_back(v1);
					}
				}
			v_pool.clear();
			if (v_pool_sudo.size()) {
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
		VSets.push_back(v_set);
	}

	E_O_rend.resize(6, mpEs.size() * 2);
	composit_edges_colors(mO_center, mpEs, E_O_rend);

	E_I_rend = E_rend;

	vector<int> colors_num(4, 0);
	for (uint32_t i = 0; i < mpEs.size(); ++i) {
		switch (get<4>(mpEs[i]))
		{
		case Edge_tag::R: colors_num[0]++;  break;
		case Edge_tag::B: colors_num[1]++;  break;
		case Edge_tag::D: colors_num[2]++;  break;
		case Edge_tag::H: colors_num[3]++; break;
		default: throw std::runtime_error("stuck at a invalid color!");
		}
	}
	vector<MatrixXf> E_colrs(4); 
	for (uint32_t i = 0; i < 4; i++) {
		E_colrs[i].resize(12, colors_num[i]);
		colors_num[i] = 0;
	}	
	for (uint32_t i = 0; i < mpEs.size(); ++i) {
		uint32_t id0 = get<0>(mpEs[i]), id1 = get<1>(mpEs[i]);
		Vector3f v0=mV[0].col(id0), v1 = mV[0].col(id1), v0_ = mO_center.col(id0), v1_ = mO_center.col(id1);
		switch (get<4>(mpEs[i]))
		{
		case Edge_tag::R: E_colrs[0].col(colors_num[0]++) << v0, v1, v0_, v1_;  break;
		case Edge_tag::B: E_colrs[1].col(colors_num[1]++) << v0, v1, v0_, v1_;  break;
		case Edge_tag::D: E_colrs[2].col(colors_num[2]++) << v0, v1, v0_, v1_;  break;
		case Edge_tag::H: E_colrs[3].col(colors_num[3]++) << v0, v1, v0_, v1_;  break;
		default: throw std::runtime_error("stuck at a invalid color!");
		}
	}
}
void MultiResolutionHierarchy::tagging_collapseTet()
{
	uint32_t INVALID_V = mV_tag.cols(), INVALID_E = mpEs.size(), INVALID_F = mpFvs.size(), INVALID_P = mpPs.size();
	pV_map.clear(); pV_map.resize(mV_tag.cols());
	Reverse_pV_map.clear(); Reverse_pV_map.resize(mV_tag.cols());
	pE_map.clear();pE_map.resize(mpEs.size());
	Reverse_pE_map.clear(); Reverse_pE_map.resize(mpEs.size());
	pF_map.clear(); pF_map.resize(mpFvs.size());
	Reverse_pF_map.clear(); Reverse_pF_map.resize(mpFvs.size());

	PV_npfs.clear(); PV_npfs.resize(mV_tag.cols());
	PV_npfs_sudo.clear(); PV_npfs_sudo.resize(mV_tag.cols());
	PE_npfs.clear(); PE_npfs.resize(mpEs.size());
	PE_npfs_sudo.clear(); PE_npfs_sudo.resize(mpEs.size());
	PF_npps.clear(); PF_npps.resize(mpFvs.size());
	PF_npps_sudo.clear(); PF_npps_sudo.resize(mpFvs.size());

	std::vector<bool> mV_B_flag(mV_tag.cols(), false);//boundary flag	
	std::vector<uint32_t> F_mapping(mpFvs.size(), INVALID_F), P_mapping(mpPs.size(), INVALID_P);//for local use in topology check
//V_map, Reverse_V_map, E_map, Reverse_E_map, pF_map, Reverse_pF_map
	for (uint32_t i = 0; i < mV_tag.cols(); i++) { pV_map[i] = i; Reverse_pV_map[i].push_back(i); }
	for (uint32_t i = 0; i < mpEs.size(); i++) { pE_map[i] = i; Reverse_pE_map[i].push_back(i); }
	for (uint32_t i = 0; i < mpFvs.size(); i++) { pF_map[i] = i; Reverse_pF_map[i].push_back(i); }
	//PV_npfs, PE_npfs, PF_npps
	for (uint32_t i = 0; i < mpFvs.size(); i++) for (uint32_t j = 0; j < mpFvs[i].size(); ++j) { PV_npfs[mpFvs[i][j]].push_back(i); PE_npfs[mpFes[i][j]].push_back(i); }
	for (uint32_t i = 0; i < mpPs.size(); i++) for (uint32_t j = 0; j < mpPs[i].size(); ++j) PF_npps[mpPs[i][j]].push_back(i);
	//mV_B_flag
	for (uint32_t i = 0; i < mpEs.size(); i++) if (std::get<2>(mpEs[i])) { mV_B_flag[std::get<1>(mpEs[i])] = mV_B_flag[std::get<0>(mpEs[i])] = true; }
	//mV_flag, mE_flag, mF_flag, mpP_flag
	mpV_flag.resize(mV_tag.cols()); std::fill(mpV_flag.begin(), mpV_flag.end(), false);
	mpE_flag.resize(mpEs.size()); std::fill(mpE_flag.begin(), mpE_flag.end(), false);
	mpF_flag.resize(mpFvs.size()); std::fill(mpF_flag.begin(), mpF_flag.end(), false);
	mpP_flag.resize(mpPs.size()); std::fill(mpP_flag.begin(), mpP_flag.end(), true);
	
	//re-coloring
	std::vector<uint32_t> E_TimeStamp(mpEs.size(), 0);
	std::vector<std::vector<uint32_t>> pV_npes(mV_tag.cols());
	for (auto e : mpEs) { 
		pV_npes[std::get<0>(e)].push_back(std::get<5>(e)); 
		pV_npes[std::get<1>(e)].push_back(std::get<5>(e));
	}
	//point doublets
	std::vector<bool> V_doublets_Flag(mV_tag.cols(), false);
	std::vector<uint32_t> v0_mapR, v1_mapR; uint32_t Remove_Doublets_Num = 0;
	//
	uint32_t iteration = 0, once = false; 
	bool topology = true; uint32_t Es_reddash_N = Es_red.size(); uint32_t left_NUM = -1, left_NUM2 = -1;
	while (!Es_red.empty() && topology) {

		std::function<bool(std::vector<std::vector<uint32_t>> &, std::vector<uint32_t> &, std::vector<uint32_t> &)> merge_polyhedra = [&](std::vector<std::vector<uint32_t>> &pfs, std::vector<uint32_t> &ps, std::vector<uint32_t> &pf) -> bool {
			//check simplicity of the polyhedral
			std::vector<uint32_t>vs_disgard(0), es_disgard(0), fs_disgard(0);
			if (!simple_polyhedral(pfs, pf, vs_disgard, es_disgard, fs_disgard))
				return false;
			//perform connectivity editing
			mpPs[ps[0]] = pf;
			for (uint32_t i = 1; i < ps.size(); i++) {
				std::vector<uint32_t>().swap(mpPs[ps[i]]);
				mpP_flag[ps[i]] = false;
			}

			for (uint32_t i = 0; i < pf.size(); i++) {
				if (std::find(PF_npps[pf[i]].begin(), PF_npps[pf[i]].end(), ps[0]) == PF_npps[pf[i]].end())
					PF_npps[pf[i]].push_back(ps[0]);
				std::vector<uint32_t> nps; nps.reserve(PF_npps[pf[i]].size());
				for (uint32_t j = 0; j < PF_npps[pf[i]].size(); j++)
					if (mpP_flag[PF_npps[pf[i]][j]]) nps.push_back(PF_npps[pf[i]][j]);
				nps.swap(PF_npps[pf[i]]);
			}

			for (uint32_t i = 0; i < es_disgard.size(); i++) {
				pE_map[es_disgard[i]] = INVALID_E; std::vector<uint32_t>().swap(Reverse_pE_map[es_disgard[i]]); std::vector<uint32_t>().swap(PE_npfs[es_disgard[i]]);
			}
			for (uint32_t i = 0; i < vs_disgard.size(); i++) {
				pV_map[vs_disgard[i]] = INVALID_V; std::vector<uint32_t>().swap(Reverse_pV_map[vs_disgard[i]]); std::vector<uint32_t>().swap(PV_npfs[vs_disgard[i]]);
			}
			for (uint32_t i = 0; i < fs_disgard.size(); i++) mpF_flag[fs_disgard[i]] = true;
			for (uint32_t j = 0; j < pf.size(); ++j) {
				for (uint32_t k = 0; k < mpFes[pf[j]].size(); k++)
					mpE_flag[mpFes[pf[j]][k]] = true;
				for (uint32_t k = 0; k < mpFvs[pf[j]].size(); k++)
					mpV_flag[mpFvs[pf[j]][k]] = true;
			}

			for (uint32_t j = 0; j < pf.size(); ++j) {
				for (uint32_t k = 0; k < mpFes[pf[j]].size(); k++)
					if (mpE_flag[mpFes[pf[j]][k]]) {
						uint32_t eid = mpFes[pf[j]][k];
						std::vector<uint32_t> nfs; nfs.reserve(PE_npfs[eid].size());
						for (uint32_t m = 0; m < PE_npfs[eid].size(); m++) if (!mpF_flag[PE_npfs[eid][m]]) nfs.push_back(PE_npfs[eid][m]);
						nfs.swap(PE_npfs[eid]); mpE_flag[mpFes[pf[j]][k]] = false;
					}
				for (uint32_t k = 0; k < mpFvs[pf[j]].size(); k++)
					if (mpV_flag[mpFvs[pf[j]][k]]) {
						uint32_t vid = mpFvs[pf[j]][k];
						std::vector<uint32_t> nfs; nfs.reserve(PV_npfs[vid].size());
						for (uint32_t m = 0; m < PV_npfs[vid].size(); m++) if (!mpF_flag[PV_npfs[vid][m]]) nfs.push_back(PV_npfs[vid][m]);
						nfs.swap(PV_npfs[vid]); mpV_flag[mpFvs[pf[j]][k]] = false;
					}
			}

			for (uint32_t i = 0; i < fs_disgard.size(); i++) {
				std::vector<uint32_t>().swap(mpFvs[fs_disgard[i]]);
				std::vector<uint32_t>().swap(mpFes[fs_disgard[i]]);
				std::vector<uint32_t>().swap(Reverse_pF_map[fs_disgard[i]]);
				std::vector<uint32_t>().swap(PF_npps[fs_disgard[i]]);
				pF_map[fs_disgard[i]] = INVALID_F;
				mpF_flag[fs_disgard[i]] = false;
			}
			return true;
		};
		std::function<void(std::vector<uint32_t> &)> degenerate_polyhedra = [&](std::vector<uint32_t> &PPs_ring) -> void {
			//degenerate polyhedra
			for (uint32_t j = 0; j < PPs_ring.size(); j++) {
				uint32_t p = PPs_ring[j];
				if (mpPs[p].size() == 2) {//merge pf0 and pf1 -> pf0
					uint32_t pf0 = mpPs[p][0], pf1 = mpPs[p][1];
					//update mPPs
					mpP_flag[p] = false; std::vector<uint32_t>().swap(mpPs[p]);
					//update Reverse_pF_map, boundaryness
					Reverse_pF_map[pf0].insert(Reverse_pF_map[pf0].end(), Reverse_pF_map[pf1].begin(), Reverse_pF_map[pf1].end());
					for (uint32_t m = 0; m < Reverse_pF_map[pf1].size(); m++) pF_map[Reverse_pF_map[pf1][m]] = pf0;
					if (mpF_boundary_flag[pf0] || mpF_boundary_flag[pf1])
						for (uint32_t m = 0; m < Reverse_pF_map[pf0].size(); m++) mpF_boundary_flag[Reverse_pF_map[pf0][m]] = true;
					//mpPs
					for (uint32_t k = 0; k < PF_npps[pf1].size(); k++) {
						uint32_t po = PF_npps[pf1][k];
						if (PF_npps[pf0].size() && std::find(mpPs[po].begin(), mpPs[po].end(), pf0) == mpPs[po].end())
							std::replace(mpPs[po].begin(), mpPs[po].end(), pf1, pf0);
						else mpPs[po].erase(std::remove(mpPs[po].begin(), mpPs[po].end(), pf1), mpPs[po].end());
					}
					//PF_npps
					PF_npps[pf0].insert(PF_npps[pf0].end(), PF_npps[pf1].begin(), PF_npps[pf1].end());
					std::vector<uint32_t> nps_temp; nps_temp.reserve(2);
					for (uint32_t pp = 0; pp < PF_npps[pf0].size(); pp++) if (mpP_flag[PF_npps[pf0][pp]]) nps_temp.push_back(PF_npps[pf0][pp]);
					std::sort(nps_temp.begin(), nps_temp.end()); nps_temp.erase(std::unique(nps_temp.begin(), nps_temp.end()), nps_temp.end());
					nps_temp.swap(PF_npps[pf0]);
					//PE_npfs
					for (uint32_t k = 0; k < mpFes[pf1].size(); k++)
						PE_npfs[mpFes[pf1][k]].erase(std::remove(PE_npfs[mpFes[pf1][k]].begin(), PE_npfs[mpFes[pf1][k]].end(), pf1), PE_npfs[mpFes[pf1][k]].end());
					//change mPV_npfs
					for (uint32_t v = 0; v < mpFvs[pf1].size(); v++) {
						uint32_t vid = mpFvs[pf1][v];
						PV_npfs[vid].erase(std::remove(PV_npfs[vid].begin(), PV_npfs[vid].end(), pf1), PV_npfs[vid].end());
					}
					//delete pf1
					std::vector<uint32_t>().swap(mpFes[pf1]);
					std::vector<uint32_t>().swap(mpFvs[pf1]);
					std::vector<uint32_t>().swap(PF_npps[pf1]);
					std::vector<uint32_t>().swap(Reverse_pF_map[pf1]);
				}
				else if (mpPs[p].size() == 3) {//one face pf0 is composed by another two: pf1 and pf2 //find pf0;
					uint32_t pf0 = mpPs[p][0];
					uint32_t pf12[2]; pf12[0] = mpPs[p][1]; pf12[1] = mpPs[p][2];

					if (mpFes[pf0].size() < mpFes[pf12[0]].size()) std::swap(pf0, pf12[0]);
					if (mpFes[pf0].size() < mpFes[pf12[1]].size()) std::swap(pf0, pf12[1]);
					if (mpFvs[pf0].size() + 2 != (mpFvs[pf12[0]].size() + mpFvs[pf12[1]].size())) {
						continue; 
					}
					if (PF_npps[pf0].size() == 2 && PF_npps[pf12[0]].size() == 2) {
						if (PF_npps[pf0][0] > PF_npps[pf0][1]) std::swap(PF_npps[pf0][0], PF_npps[pf0][1]);
						if (PF_npps[pf12[0]][0] > PF_npps[pf12[0]][1]) std::swap(PF_npps[pf12[0]][0], PF_npps[pf12[0]][1]);
						if (PF_npps[pf0][0] == PF_npps[pf12[0]][0] && PF_npps[pf0][1] == PF_npps[pf12[0]][1]) std::swap(pf12[1], pf0);
					}
					if (PF_npps[pf0].size() == 2 && PF_npps[pf12[1]].size() == 2) {
						if (PF_npps[pf0][0] > PF_npps[pf0][1]) std::swap(PF_npps[pf0][0], PF_npps[pf0][1]);
						if (PF_npps[pf12[1]][0] > PF_npps[pf12[1]][1]) std::swap(PF_npps[pf12[1]][0], PF_npps[pf12[1]][1]);
						if (PF_npps[pf0][0] == PF_npps[pf12[1]][0] && PF_npps[pf0][1] == PF_npps[pf12[1]][1]) std::swap(pf12[0], pf0);
					}

					mpP_flag[p] = false; std::vector<uint32_t>().swap(mpPs[p]);
					for (uint32_t j = 0; j<2; j++)
						if (mpF_boundary_flag[pf0] && mpF_boundary_flag[pf12[j]]) {
							//PE_npfs
							for (uint32_t k = 0; k < mpFes[pf12[j]].size(); k++)
								PE_npfs[mpFes[pf12[j]][k]].erase(std::remove(PE_npfs[mpFes[pf12[j]][k]].begin(), PE_npfs[mpFes[pf12[j]][k]].end(), pf12[j]), PE_npfs[mpFes[pf12[j]][k]].end());
							//PV_npfs
							for (uint32_t v = 0; v < mpFvs[pf12[j]].size(); v++) {
								uint32_t vid = mpFvs[pf12[j]][v];
								PV_npfs[vid].erase(std::remove(PV_npfs[vid].begin(), PV_npfs[vid].end(), pf12[j]), PV_npfs[vid].end());
							}
							//delete pf12[j]
							std::vector<uint32_t>().swap(mpFes[pf12[j]]);
							std::vector<uint32_t>().swap(mpFvs[pf12[j]]);
							std::vector<uint32_t>().swap(PF_npps[pf12[j]]);
							pF_map[pf12[j]] = INVALID_F;
						}

					//update Reverse_pF_map, boundaryness
					if (mpF_boundary_flag[pf0] || mpF_boundary_flag[pf12[0]])
						for (uint32_t m = 0; m < Reverse_pF_map[pf12[0]].size(); m++) {
							uint32_t fff = Reverse_pF_map[pf12[0]][m];
							mpF_boundary_flag[fff] = true;
							for (uint32_t n = 0; n < mpFes[fff].size(); n++)
								std::get<2>(mpEs[mpFes[fff][n]]) = true;
						}
					if (mpF_boundary_flag[pf0] || mpF_boundary_flag[pf12[1]])
						for (uint32_t m = 0; m < Reverse_pF_map[pf12[1]].size(); m++) {
							uint32_t fff = Reverse_pF_map[pf12[1]][m];
							mpF_boundary_flag[fff] = true;
							for (uint32_t n = 0; n < mpFes[fff].size(); n++)
								std::get<2>(mpEs[mpFes[fff][n]]) = true;
						}
					//mpPs
					for (uint32_t pp = 0; pp < PF_npps[pf0].size(); pp++) {
						uint32_t po = PF_npps[pf0][pp];
						if (mpP_flag[po]) {
							if (PF_npps[pf12[0]].size() && std::find(mpPs[po].begin(), mpPs[po].end(), pf12[0]) == mpPs[po].end())
								std::replace(mpPs[po].begin(), mpPs[po].end(), pf0, pf12[0]);
							else mpPs[po].erase(std::remove(mpPs[po].begin(), mpPs[po].end(), pf0), mpPs[po].end());
							if (PF_npps[pf12[1]].size() && std::find(mpPs[po].begin(), mpPs[po].end(), pf12[1]) == mpPs[po].end())
								mpPs[po].push_back(pf12[1]);
						}
					}
					//PF_npps
					for (uint32_t k = 0; k < 2; k++) {
						PF_npps[pf12[k]].insert(PF_npps[pf12[k]].end(), PF_npps[pf0].begin(), PF_npps[pf0].end());
						std::vector<uint32_t> nps_temp; nps_temp.reserve(2);
						for (uint32_t pp = 0; pp < PF_npps[pf12[k]].size(); pp++)
							if (mpP_flag[PF_npps[pf12[k]][pp]]) nps_temp.push_back(PF_npps[pf12[k]][pp]);
						std::sort(nps_temp.begin(), nps_temp.end()); nps_temp.erase(std::unique(nps_temp.begin(), nps_temp.end()), nps_temp.end());
						nps_temp.swap(PF_npps[pf12[k]]);
					}
					//PE_npfs
					for (uint32_t k = 0; k < mpFes[pf0].size(); k++)
						PE_npfs[mpFes[pf0][k]].erase(std::remove(PE_npfs[mpFes[pf0][k]].begin(), PE_npfs[mpFes[pf0][k]].end(), pf0), PE_npfs[mpFes[pf0][k]].end());
					//PV_npfs
					for (uint32_t v = 0; v < mpFvs[pf0].size(); v++) {
						uint32_t vid = mpFvs[pf0][v];
						PV_npfs[vid].erase(std::remove(PV_npfs[vid].begin(), PV_npfs[vid].end(), pf0), PV_npfs[vid].end());
					}
					//delete pf0
					std::vector<uint32_t>().swap(mpFes[pf0]);
					std::vector<uint32_t>().swap(mpFvs[pf0]);
					std::vector<uint32_t>().swap(PF_npps[pf0]);
					pF_map[pf0] = INVALID_F;
				}
			}
		};
		std::function<bool(std::vector<uint32_t> &, std::vector<uint32_t> &)> edge_fuse_polygons = [&](std::vector<uint32_t> &pf_test, std::vector<uint32_t> & PPs_ring_set)->bool {
			std::sort(pf_test.begin(), pf_test.end());
			pf_test.erase(std::unique(pf_test.begin(), pf_test.end()), pf_test.end());
			std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> fs_tuples;
			for (uint32_t j = 0; j < pf_test.size(); j++) {
				if (PF_npps[pf_test[j]].size() == 2) {
					uint32_t p0 = PF_npps[pf_test[j]][0], p1 = PF_npps[pf_test[j]][1];
					if (p0 > p1) std::swap(p0, p1);
					fs_tuples.push_back(std::make_tuple(p0, p1, pf_test[j]));
				}
				else fs_tuples.push_back(std::make_tuple(PF_npps[pf_test[j]][0], INVALID_P, pf_test[j]));
			}
			std::sort(fs_tuples.begin(), fs_tuples.end());
			std::vector<std::vector<uint32_t>> f_sets, f_sets_sudo; std::vector<uint32_t> f_set;
			for (uint32_t j = 0; j < fs_tuples.size(); j++) {
				if (j == 0 || (j != 0 && (std::get<0>(fs_tuples[j]) != std::get<0>(fs_tuples[j - 1]) || std::get<1>(fs_tuples[j]) != std::get<1>(fs_tuples[j - 1])))) {
					if (j != 0 && f_set.size() >= 2) f_sets.push_back(f_set);
					f_set.clear(); f_set.reserve(2); f_set.push_back(std::get<2>(fs_tuples[j]));
				}
				else f_set.push_back(std::get<2>(fs_tuples[j]));
				if (j + 1 == fs_tuples.size() && f_set.size() >= 2) f_sets.push_back(f_set);
			}
			//decompose into pairwise faces
			for (uint32_t j = 0; j < f_sets.size(); j++) {
				f_set = f_sets[j];
				while (f_set.size() > 1) {
					uint32_t fid0 = f_set[0], fid1 = INVALID_F;
					f_set.erase(std::remove(f_set.begin(), f_set.end(), fid0), f_set.end());
					std::sort(f_set.begin(), f_set.end());

					for (uint32_t k = 0; k < mpFes[fid0].size(); k++) {
						uint32_t eid = mpFes[fid0][k];
						if (std::get<4>(mpEs[eid]) == Edge_tag::D)
						{
							std::sort(PE_npfs[eid].begin(), PE_npfs[eid].end());
							std::vector<uint32_t> common_f;
							std::set_intersection(PE_npfs[eid].begin(), PE_npfs[eid].end(), f_set.begin(), f_set.end(), std::back_inserter(common_f));
							if (common_f.size()) { fid1 = common_f[0]; break; }
							else continue;
						}
					}
					if (once && fid1 == INVALID_F) {
						for (uint32_t k = 0; k < mpFes[fid0].size(); k++) {
							uint32_t eid = mpFes[fid0][k];
							if ((mpF_boundary_flag[fid0] && std::get<4>(mpEs[eid]) == Edge_tag::D) || !mpF_boundary_flag[fid0])
							{
								std::sort(PE_npfs[eid].begin(), PE_npfs[eid].end());
								std::vector<uint32_t> common_f;
								std::set_intersection(PE_npfs[eid].begin(), PE_npfs[eid].end(), f_set.begin(), f_set.end(), std::back_inserter(common_f));
								if (common_f.size()) { fid1 = common_f[0]; break; }
								else continue;
							}
						}
					}
					if (fid1 != INVALID_F) {
						std::vector<uint32_t> a_pair;
						a_pair.push_back(fid0); a_pair.push_back(fid1);
						f_sets_sudo.push_back(a_pair);
						f_set.erase(std::remove(f_set.begin(), f_set.end(), fid1), f_set.end());
					}
				}
			}
			f_sets_sudo.swap(f_sets);

			std::vector<uint32_t> disgarded_fs; disgarded_fs.reserve(f_sets.size());
			for (uint32_t j = 0; j < f_sets.size(); j++) {
				//find shared e
				std::sort(mpFes[f_sets[j][0]].begin(), mpFes[f_sets[j][0]].end());
				std::sort(mpFes[f_sets[j][1]].begin(), mpFes[f_sets[j][1]].end());
				std::vector<uint32_t> common_es;
				std::set_intersection(mpFes[f_sets[j][0]].begin(), mpFes[f_sets[j][0]].end(), mpFes[f_sets[j][1]].begin(), mpFes[f_sets[j][1]].end(), std::back_inserter(common_es));
				if (common_es.size() != 1) continue;
				if (once) {
					if (mpF_boundary_flag[f_sets[j][0]] && std::get<4>(mpEs[common_es[0]]) != Edge_tag::D) continue;
				}
				else {
					if (std::get<4>(mpEs[common_es[0]]) != Edge_tag::D)
						continue;
				}
				//if e == D judge is_simple_polygon
				std::vector<std::vector<uint32_t>> fv_sets(2), fe_sets(2);
				std::vector<uint32_t> poly_vs, poly_es, vs_disgard, es_disgard;
				fv_sets[0] = mpFvs[f_sets[j][0]]; fv_sets[1] = mpFvs[f_sets[j][1]];
				fe_sets[0] = mpFes[f_sets[j][0]]; fe_sets[1] = mpFes[f_sets[j][1]];
				if (!simple_polygon_3D(fv_sets, fe_sets, poly_vs, poly_es, vs_disgard, es_disgard, true))
					continue;
				//connectivity editing
				uint32_t f0 = f_sets[j][0], f1 = f_sets[j][1]; disgarded_fs.push_back(f1);

				PPs_ring_set.insert(PPs_ring_set.end(), PF_npps[f0].begin(), PF_npps[f0].end());

				mpFes[f0] = poly_es; mpFvs[f0] = poly_vs;

				for (uint32_t k = 0; k < PF_npps[f0].size(); k++) {
					uint32_t pid = PF_npps[f0][k];
					mpPs[pid].erase(std::remove(mpPs[pid].begin(), mpPs[pid].end(), f1), mpPs[pid].end());
				}
				std::vector<uint32_t>().swap(mpFvs[f1]); std::vector<uint32_t>().swap(mpFes[f1]);
				std::vector<uint32_t>().swap(PF_npps[f1]); std::vector<uint32_t>().swap(Reverse_pF_map[f1]);

				for (uint32_t k = 0; k < es_disgard.size(); k++) { pE_map[es_disgard[k]] = INVALID_E; std::vector<uint32_t>().swap(Reverse_pE_map[es_disgard[k]]); std::vector<uint32_t>().swap(PE_npfs[es_disgard[k]]); }
				for (uint32_t k = 0; k < vs_disgard.size(); k++) {
					pV_map[vs_disgard[k]] = INVALID_V; std::vector<uint32_t>().swap(Reverse_pV_map[vs_disgard[k]]); std::vector<uint32_t>().swap(PV_npfs[vs_disgard[k]]);
				}
				mpF_flag[f1] = true;
				for (uint32_t k = 0; k < poly_es.size(); k++) {
					uint32_t eid = poly_es[k];
					if (std::find(PE_npfs[eid].begin(), PE_npfs[eid].end(), f0) == PE_npfs[eid].end())
						PE_npfs[eid].push_back(f0);
					std::vector<uint32_t> nfs; nfs.reserve(PE_npfs[eid].size());
					for (uint32_t m = 0; m < PE_npfs[eid].size(); m++) if (!mpF_flag[PE_npfs[eid][m]]) nfs.push_back(PE_npfs[eid][m]);
					nfs.swap(PE_npfs[eid]);
				}
				for (uint32_t k = 0; k < poly_vs.size(); k++) {
					uint32_t vid = poly_vs[k];
					if (std::find(PV_npfs[vid].begin(), PV_npfs[vid].end(), f0) == PV_npfs[vid].end())
						PV_npfs[vid].push_back(f0);
					std::vector<uint32_t> nfs; nfs.reserve(PV_npfs[vid].size());
					for (uint32_t m = 0; m < PV_npfs[vid].size(); m++) if (!mpF_flag[PV_npfs[vid][m]]) nfs.push_back(PV_npfs[vid][m]);
					nfs.swap(PV_npfs[vid]);
				}
				mpF_flag[f1] = false;
			}

			if (!disgarded_fs.size()) return false;
			for (uint32_t j = 0; j < disgarded_fs.size(); j++) mpF_flag[disgarded_fs[j]] = true;
			std::vector<uint32_t> fs_total_temp; fs_total_temp.reserve(pf_test.size());
			for (uint32_t j = 0; j < pf_test.size(); j++) 
				if (!mpF_flag[pf_test[j]]) { fs_total_temp.push_back(pf_test[j]); }
			fs_total_temp.swap(pf_test);
			for (uint32_t j = 0; j < disgarded_fs.size(); j++) mpF_flag[disgarded_fs[j]] = false;

			return true;
		};
		std::function<void(uint32_t &, uint32_t &)> face_swap = [&](uint32_t &Reid, uint32_t &Heid) {
			uint32_t v0_ = pV_map[std::get<0>(mpEs[Reid])], v1_ = pV_map[std::get<1>(mpEs[Reid])];
			std::vector<uint32_t> fs_total = PV_npfs[v0_], ps_total; fs_total.insert(fs_total.end(), PV_npfs[v1_].begin(), PV_npfs[v1_].end());

			for (uint32_t j = 0; j < fs_total.size(); j++) ps_total.insert(ps_total.end(), PF_npps[fs_total[j]].begin(), PF_npps[fs_total[j]].end());
			std::sort(ps_total.begin(), ps_total.end()); ps_total.erase(std::unique(ps_total.begin(), ps_total.end()), ps_total.end());
			fs_total.clear(); for (uint32_t j = 0; j < ps_total.size(); j++) fs_total.insert(fs_total.end(), mpPs[ps_total[j]].begin(), mpPs[ps_total[j]].end());
			std::sort(fs_total.begin(), fs_total.end()); fs_total.erase(std::unique(fs_total.begin(), fs_total.end()), fs_total.end());
			
			std::vector<uint32_t> candidate_es;
			for (uint32_t j = 0; j < fs_total.size(); j++) for (uint32_t k = 0; k < mpFes[fs_total[j]].size(); k++) {
				uint32_t eid_c = mpFes[fs_total[j]][k];
				if (std::get<2>(mpEs[eid_c])) continue;//should not be on boundary
				if (PE_npfs[eid_c].size() != 3) continue;//only has three neighbor fs && ps.
				uint32_t ev0 = pV_map[std::get<0>(mpEs[eid_c])], ev1 = pV_map[std::get<1>(mpEs[eid_c])];
				if (ev0 == v0_ || ev0 == v1_ || ev1 == v0_ || ev1 == v1_) continue;// should not touch Reid
				candidate_es.push_back(eid_c);
			}
			std::sort(candidate_es.begin(), candidate_es.end()); candidate_es.erase(std::unique(candidate_es.begin(), candidate_es.end()), candidate_es.end());

			Heid = INVALID_E;
			for (uint32_t j = 0; j < candidate_es.size(); j++) {
				
				std::vector<uint32_t> ps; ps.reserve(PE_npfs[candidate_es[j]].size());
				for (uint32_t k = 0; k< PE_npfs[candidate_es[j]].size(); k++)
					ps.insert(ps.end(), PF_npps[PE_npfs[candidate_es[j]][k]].begin(), PF_npps[PE_npfs[candidate_es[j]][k]].end());
				std::sort(ps.begin(), ps.end()); ps.erase(std::unique(ps.begin(), ps.end()), ps.end());	
				std::vector<std::vector<uint32_t>> pfs(ps.size());
				for (uint32_t k = 0; k < ps.size(); k++) pfs[k] = mpPs[ps[k]];
				std::vector<uint32_t> pf, vs_disgard(0), es_disgard(0), fs_disgard(0);
				if (!simple_polyhedral(pfs, pf, vs_disgard, es_disgard, fs_disgard)) continue;
				if (pf.size() != 6 || vs_disgard.size() != 0 || !(es_disgard.size() == 1 && es_disgard[0] == candidate_es[j]) || fs_disgard.size() != 3)
					continue;
				bool on_p_boundary = false;
				for (uint32_t k = 0; k < pf.size(); k++) for (uint32_t m = 0; m < mpFes[pf[k]].size(); m++) if (Reid == mpFes[pf[k]][m]) on_p_boundary = true;
				if (on_p_boundary) { Heid = candidate_es[j]; break; }
			}
		};
		iteration++;
		topology = false;
		int size_pool = Es_red.size();
		std::vector<tuple_E> Es_nonmanifold;
		for (uint32_t i = 0; i < size_pool; i++) {
			tuple_E e = Es_red.top(); Es_red.pop();
			uint32_t e_map = pE_map[std::get<5>(e)];

			if (e_map == INVALID_E) { 
				continue; }

			if (std::get<7>(e) < E_TimeStamp[e_map]) continue;

			if (!PE_npfs[e_map].size()) { 
				continue; }
			if (std::get<4>(e) == std::get<4>(mpEs[e_map])) { 
				continue; }
			if (std::get<4>(mpEs[e_map]) != Edge_tag::B) {
				continue;
			}

			uint32_t v0 = std::get<0>(e), v1 = std::get<1>(e);
			uint32_t v0_map = pV_map[v0], v1_map = pV_map[v1];
			if (v0_map == v1_map) { std::cout << "somewhere is wrong" << endl; system("PAUSE"); }

			std::vector<uint32_t> pf_test(0);
			if (std::get<4>(e) == Edge_tag::H || std::get<4>(e) == Edge_tag::D) {
				if (std::get<4>(e) == Edge_tag::H) {
					std::get<4>(mpEs[e_map]) = Edge_tag::H;
					//all the surrounding faces should be glued together!
					std::vector<uint32_t> ps; ps.reserve(PE_npfs[e_map].size());
					std::vector<std::vector<uint32_t>> pfs;
					for (uint32_t j = 0; j < PE_npfs[e_map].size(); j++)
						ps.insert(ps.end(), PF_npps[PE_npfs[e_map][j]].begin(), PF_npps[PE_npfs[e_map][j]].end());
					std::sort(ps.begin(), ps.end());
					ps.erase(std::unique(ps.begin(), ps.end()), ps.end());
					pfs.resize(ps.size());
					for (uint32_t j = 0; j < ps.size(); j++) pfs[j]= mpPs[ps[j]];
					std::vector<uint32_t> pf;
					if(merge_polyhedra(pfs,ps,pf))
						pf_test = pf;
				}
				else if (std::get<4>(e) == Edge_tag::D) {
					std::get<4>(mpEs[e_map]) = Edge_tag::D;

					std::vector<uint32_t> npfs_temp = PE_npfs[e_map];
					for (uint32_t k = 0; k < npfs_temp.size(); k++) {
						uint32_t fid = npfs_temp[k];
						if (pF_map[fid] == INVALID_F) { continue; }
						if (mpF_boundary_flag[fid]) continue;//f on boundary
						std::vector<short> sides;
						for (uint32_t j = 0; j < mpFes[fid].size(); j++) {
							uint32_t eid = mpFes[fid][j];
							if (std::get<4>(mpEs[eid]) == Edge_tag::D) sides.push_back(std::get<6>(mpEs[eid]));
						}
						if (sides.size() < 2) continue;
						short which_side = sides[0]; bool pass = false;
						for (uint32_t m = 1; m < sides.size();m++) if (which_side != sides[m]) { pass = true; break; }
						std::vector<std::vector<uint32_t>> pfs(2);
						for (uint32_t j = 0; j < PF_npps[fid].size(); j++)
							pfs[j] = mpPs[PF_npps[fid][j]];
						//check simplicity of the polyhedral
						std::vector<uint32_t> pf, ps = PF_npps[fid];
						if (merge_polyhedra(pfs, ps, pf)) {
							if (pf_test.size()) {
								pf_test.insert(pf_test.end(), pf.begin(), pf.end());
								std::sort(pf_test.begin(), pf_test.end());
								pf_test.erase(std::unique(pf_test.begin(), pf_test.end()), pf_test.end());
								std::vector<uint32_t> pf_test_sudo; pf_test_sudo.reserve(pf_test.size());
								for (uint32_t j = 0; j < pf_test.size(); j++) if (pF_map[pf_test[j]] != INVALID_F) pf_test_sudo.push_back(pf_test[j]);
								pf_test_sudo.swap(pf_test);
							}else pf_test.insert(pf_test.end(), pf.begin(), pf.end());
						}
					}
					pf_test.insert(pf_test.end(), PE_npfs[e_map].begin(), PE_npfs[e_map].end());
				}
			}
			else if (std::get<4>(e) == Edge_tag::R) {
				if (mV_B_flag[v0_map] && mV_B_flag[v1_map] && !std::get<2>(mpEs[e_map])) {
					std::get<3>(e) *= 2;
					Es_nonmanifold.push_back(e);
					continue;
				}
				std::sort(PV_npfs[v0_map].begin(), PV_npfs[v0_map].end());
				std::sort(PV_npfs[v1_map].begin(), PV_npfs[v1_map].end());
				std::vector<uint32_t> common_fs;
				std::set_intersection(PV_npfs[v0_map].begin(), PV_npfs[v0_map].end(), PV_npfs[v1_map].begin(), PV_npfs[v1_map].end(), std::back_inserter(common_fs));
				
				if (common_fs.size() != PE_npfs[e_map].size()){
					std::get<3>(e) *= 2;
					Es_nonmanifold.push_back(e);
					continue;
				}
				int genus_pre = -1, genus_aft = -1; bool manifoldness_aft = true;

				std::function<bool(std::vector<std::vector<uint32_t>> &, std::vector<std::vector<uint32_t>> &)> check_pre = [&](
					std::vector<std::vector<uint32_t>> &pFes_local, std::vector<std::vector<uint32_t>> &pPPs_local) -> bool {
					for (uint32_t i = 0; i < pFes_local.size(); i++) {
						for (uint32_t j = 0; j < pFes_local[i].size(); j++) {
							mpE_flag[pFes_local[i][j]] = true;
							uint32_t v0 = pV_map[std::get<0>(mpEs[pFes_local[i][j]])], v1 = pV_map[std::get<1>(mpEs[pFes_local[i][j]])];
							mpV_flag[v0] = mpV_flag[v1] = true;
						}
					}
					uint32_t e_num = 0, v_num = 0;
					for (uint32_t i = 0; i < pFes_local.size(); i++) {
						for (uint32_t j = 0; j < pFes_local[i].size(); j++) {
							if (mpE_flag[pFes_local[i][j]]) { e_num++; mpE_flag[pFes_local[i][j]] = false;}
							uint32_t v0 = pV_map[std::get<0>(mpEs[pFes_local[i][j]])], v1 = pV_map[std::get<1>(mpEs[pFes_local[i][j]])];
							if (mpV_flag[v0]) { v_num++; mpV_flag[v0] = false; }
							if (mpV_flag[v1]) { v_num++; mpV_flag[v1] = false; }
						}
					}
					genus_pre = (v_num + pFes_local.size() - e_num - pPPs_local.size() - 2);
					return true;
				};
				
				std::function<bool(std::vector<std::vector<uint32_t>> &, std::vector<bool> &, std::vector<std::vector<uint32_t>> &)> check_aft = [&](
					std::vector<std::vector<uint32_t>> &pFes_local, std::vector<bool> &pF_boundaryness, std::vector<std::vector<uint32_t>> &pPPs_local) -> bool {
					for (uint32_t i = 0; i < pFes_local.size(); i++) {
						for (uint32_t j = 0; j < pFes_local[i].size(); j++) {
							mpE_flag[pFes_local[i][j]] = true;
							uint32_t v0 = pV_map[std::get<0>(mpEs[pFes_local[i][j]])], v1 = pV_map[std::get<1>(mpEs[pFes_local[i][j]])];
							mpV_flag[v0] = mpV_flag[v1] = true;
						}
					}
					uint32_t e_num = 0, v_num = 0; std::vector<tuple_E> es_local; es_local.reserve(pFes_local.size());

					for (uint32_t i = 0; i < pFes_local.size(); i++) {
						for (uint32_t j = 0; j < pFes_local[i].size(); j++) {
							if (mpE_flag[pFes_local[i][j]]) { 
								e_num++; tuple_E e = mpEs[pFes_local[i][j]];
								std::get<0>(e) = pV_map[std::get<0>(e)];
								std::get<1>(e) = pV_map[std::get<1>(e)]; 
								if (std::get<0>(e) > std::get<1>(e)) std::swap(std::get<0>(e), std::get<1>(e));
								es_local.push_back(e);
								mpE_flag[pFes_local[i][j]] = false;
							}
							uint32_t v0 = pV_map[std::get<0>(mpEs[pFes_local[i][j]])], v1 = pV_map[std::get<1>(mpEs[pFes_local[i][j]])];
							if (mpV_flag[v0]) { v_num++; mpV_flag[v0] = false; }
							if (mpV_flag[v1]) { v_num++; mpV_flag[v1] = false; }
						}
					}
					genus_aft = (v_num + pFes_local.size() - e_num - pPPs_local.size() - 2);
					if (genus_pre != genus_aft) return false;

					//non-manifold face?
					std::vector<std::vector<uint32_t>> f_nts(pFes_local.size());
					for (uint32_t i = 0; i < pPPs_local.size(); i++)
						for (uint32_t j = 0; j < pPPs_local[i].size(); j++)
							f_nts[pPPs_local[i][j]].push_back(i);
					for (uint32_t i = 0; i < f_nts.size(); i++)
						if ((pF_boundaryness[i] && f_nts[i].size() != 1) || (!pF_boundaryness[i] && f_nts[i].size() > 2)) manifoldness_aft = false;
					if(!manifoldness_aft) return false;

					//non-manifold boundary edge?
					for (uint32_t i = 0; i < pFes_local.size(); i++) for (uint32_t j = 0; j < pFes_local[i].size(); j++) PE_npfs_sudo[pFes_local[i][j]].push_back(i);
					for (uint32_t i = 0; i < pFes_local.size(); i++) for (uint32_t j = 0; j < pFes_local[i].size(); j++) {
						if (!std::get<2>(mpEs[pFes_local[i][j]])) continue;
						uint32_t bf_num = 0;
						for (uint32_t k = 0; k < PE_npfs_sudo[pFes_local[i][j]].size(); k++)
							if (pF_boundaryness[PE_npfs_sudo[pFes_local[i][j]][k]]) bf_num++;
						if (bf_num > 2) { manifoldness_aft = false; break; }
					}
					for (uint32_t i = 0; i < pFes_local.size(); i++) for (uint32_t j = 0; j < pFes_local[i].size(); j++) PE_npfs_sudo[pFes_local[i][j]].clear();
					if (!manifoldness_aft) return false;

					//check duplicated es
					std::sort(es_local.begin(), es_local.end());
					for (uint32_t i = 0; i < es_local.size(); ++i)
						if (i != 0 && (std::get<0>(es_local[i]) == std::get<0>(es_local[i - 1]) && std::get<1>(es_local[i]) == std::get<1>(es_local[i - 1])))
							return false;

					for (uint32_t i = 0; i < pPPs_local.size(); ++i) {
						std::vector<std::vector<uint32_t>> pfes(pPPs_local[i].size());
						for (uint32_t j = 0; j < pPPs_local[i].size(); ++j)
							pfes[j] = pFes_local[pPPs_local[i][j]];
						if (!simple_polyhedral_v2(pfes))
							return false;
					}
					//check the simplicity of polygons.
					for (uint32_t i = 0; i < pFes_local.size(); i++) {
						std::vector<std::vector<uint32_t>> fv_sets(1), fe_sets(1);
						std::vector<uint32_t> poly_vs, poly_es, vs_disgard, es_disgard;
						fe_sets[0] = pFes_local[i];
						if (!simple_polygon_3D(fv_sets, fe_sets, poly_vs, poly_es, vs_disgard, es_disgard, false))
							return false;
					}
					//check duplicated fs
					std::vector<std::vector<uint32_t>> pFvs_local(pFes_local.size());
					for (uint32_t i = 0; i < pFes_local.size(); i++) {
						for (uint32_t j = 0; j < pFes_local[i].size(); j++) {
							uint32_t v0 = pV_map[std::get<0>(mpEs[pFes_local[i][j]])], v1 = pV_map[std::get<1>(mpEs[pFes_local[i][j]])];
							mpV_flag[v0] = mpV_flag[v1] = true;
						}
						std::vector<uint32_t> fvs;
						for (uint32_t j = 0; j < pFes_local[i].size(); j++) {
							uint32_t v0 = pV_map[std::get<0>(mpEs[pFes_local[i][j]])], v1 = pV_map[std::get<1>(mpEs[pFes_local[i][j]])];
							if (mpV_flag[v0]) { fvs.push_back(v0); mpV_flag[v0] = false; }
							if (mpV_flag[v1]) { fvs.push_back(v1); mpV_flag[v1] = false; }
						}
						std::sort(fvs.begin(), fvs.end());
						pFvs_local[i] = fvs;
					}
					std::sort(pFvs_local.begin(), pFvs_local.end());
					for (uint32_t j = 1; j < pFvs_local.size(); ++j)
						if (pFvs_local[j - 1].size() == pFvs_local[j].size() && std::equal(pFvs_local[j - 1].begin(), pFvs_local[j - 1].end(), pFvs_local[j].begin()))
							return false;
					return true;
				};

				std::function<bool()> check_topology = [&]() -> bool {
					std::vector<uint32_t> fs_total = PV_npfs[v0_map], ps_total;
					std::vector<std::vector<uint32_t>> pFes_local, pPPs_local;
					std::vector<bool> pF_boundaryness_local;

					fs_total.insert(fs_total.end(), PV_npfs[v1_map].begin(), PV_npfs[v1_map].end());

					for (uint32_t j = 0; j < fs_total.size(); j++) ps_total.insert(ps_total.end(), PF_npps[fs_total[j]].begin(), PF_npps[fs_total[j]].end());
					std::sort(ps_total.begin(), ps_total.end()); ps_total.erase(std::unique(ps_total.begin(), ps_total.end()), ps_total.end());

					fs_total.clear();
					for (uint32_t j = 0; j < ps_total.size(); j++) {
						fs_total.insert(fs_total.end(), mpPs[ps_total[j]].begin(), mpPs[ps_total[j]].end());
						P_mapping[ps_total[j]] = j;
					}
					std::sort(fs_total.begin(), fs_total.end()); fs_total.erase(std::unique(fs_total.begin(), fs_total.end()), fs_total.end());
					for (uint32_t j = 0; j < fs_total.size(); j++) for (uint32_t k = 0; k < PF_npps[fs_total[j]].size(); k++)  
						if(P_mapping[PF_npps[fs_total[j]][k]] != INVALID_P)
							PF_npps_sudo[fs_total[j]].push_back(PF_npps[fs_total[j]][k]);

					pFes_local.resize(fs_total.size()); pF_boundaryness_local.resize(fs_total.size());
					for (uint32_t j = 0; j < fs_total.size(); j++) {
						pFes_local[j] = mpFes[fs_total[j]]; F_mapping[fs_total[j]] = j; pF_boundaryness_local[j] = mpF_boundary_flag[fs_total[j]];
					}
					pPPs_local.resize(ps_total.size());
					for (uint32_t j = 0; j < ps_total.size(); j++)
						for (uint32_t k = 0; k < mpPs[ps_total[j]].size(); k++)
							pPPs_local[j].push_back(F_mapping[mpPs[ps_total[j]][k]]);

					//genus before
					check_pre(pFes_local, pPPs_local);
					//pseudo removal process
					bool pass_check = true;
					std::vector<bool> F_flag_local(pFes_local.size(), true), P_flag_local(pPPs_local.size(), true);
					//pseudo update mappings merge v1_map and v0_map -> v0_map
					for (uint32_t j = 0; j < Reverse_pV_map[v1_map].size(); j++)
						pV_map[Reverse_pV_map[v1_map][j]] = v0_map;

					std::vector<uint32_t> PPs_ring, newIds;
					for (uint32_t j = 0; j < PE_npfs[e_map].size(); j++) PPs_ring.insert(PPs_ring.end(), PF_npps[PE_npfs[e_map][j]].begin(), PF_npps[PE_npfs[e_map][j]].end());
					std::sort(PPs_ring.begin(), PPs_ring.end()); PPs_ring.erase(std::unique(PPs_ring.begin(), PPs_ring.end()), PPs_ring.end());
					
					//////////////////////////////////////////////////////////////////////////////
					for (uint32_t j = 0; j < PE_npfs[e_map].size(); j++) {
						uint32_t fid = PE_npfs[e_map][j];
						pFes_local[F_mapping[fid]].erase(std::remove(pFes_local[F_mapping[fid]].begin(), pFes_local[F_mapping[fid]].end(), e_map), pFes_local[F_mapping[fid]].end());
					}
					//ffs_set
					std::vector<uint32_t> ffs_set; ffs_set.reserve(PPs_ring.size() * 3);
					for (uint32_t j = 0; j < PPs_ring.size(); j++) for (uint32_t k = 0; k < mpPs[PPs_ring[j]].size(); k++) mpF_flag[mpPs[PPs_ring[j]][k]] = true;
					for (uint32_t j = 0; j < PPs_ring.size(); j++) for (uint32_t k = 0; k < mpPs[PPs_ring[j]].size(); k++)
						if (mpF_flag[mpPs[PPs_ring[j]][k]]) { ffs_set.push_back(mpPs[PPs_ring[j]][k]); mpF_flag[mpPs[PPs_ring[j]][k]] = false; }
					//ees_tuples
					std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> ees_tuples; ees_tuples.reserve(ffs_set.size() * 3);
					for (uint32_t j = 0; j < ffs_set.size(); j++) for (uint32_t k = 0; k < pFes_local[F_mapping[ffs_set[j]]].size(); k++) mpE_flag[pFes_local[F_mapping[ffs_set[j]]][k]] = true;
					for (uint32_t j = 0; j < ffs_set.size(); j++) for (uint32_t k = 0; k < pFes_local[F_mapping[ffs_set[j]]].size(); k++)
						if (mpE_flag[pFes_local[F_mapping[ffs_set[j]]][k]]) {
							uint32_t eeid = pFes_local[F_mapping[ffs_set[j]]][k];
							uint32_t evid0 = pV_map[std::get<0>(mpEs[eeid])], evid1 = pV_map[std::get<1>(mpEs[eeid])];
							if (evid0 > evid1) std::swap(evid0, evid1);
							ees_tuples.push_back(std::make_tuple(evid0, evid1, eeid));
							mpE_flag[eeid] = false;
						}
					std::sort(ees_tuples.begin(), ees_tuples.end());
					//ees_sets: collapsing pairs
					std::vector<std::vector<uint32_t>> ees_sets; std::vector<uint32_t> es_set;
					for (uint32_t j = 0; j < ees_tuples.size(); j++) {
						if (j == 0 || (j != 0 && (std::get<0>(ees_tuples[j]) != std::get<0>(ees_tuples[j - 1]) || std::get<1>(ees_tuples[j]) != std::get<1>(ees_tuples[j - 1])))) {
							if (j != 0 && es_set.size() == 2) ees_sets.push_back(es_set);
							es_set.clear(); es_set.reserve(2); es_set.push_back(std::get<2>(ees_tuples[j]));
						}
						else es_set.push_back(std::get<2>(ees_tuples[j]));
						if (j + 1 == ees_tuples.size() && es_set.size() == 2) ees_sets.push_back(es_set);
					}
					for (uint32_t j = 0; j < ees_sets.size(); j++) {
						//judge if they share a common face
						uint32_t eeid0 = ees_sets[j][0], eeid1 = ees_sets[j][1];
						std::sort(PE_npfs[eeid0].begin(), PE_npfs[eeid0].end()); std::sort(PE_npfs[eeid1].begin(), PE_npfs[eeid1].end());
						std::vector<uint32_t> ecom_fs;
						std::set_intersection(PE_npfs[eeid0].begin(), PE_npfs[eeid0].end(), PE_npfs[eeid1].begin(), PE_npfs[eeid1].end(), std::back_inserter(ecom_fs));
						if (ecom_fs.size() > 1) { 
							pass_check = false; break; }
						else if (ecom_fs.size()) {
							uint32_t fid = ecom_fs[0]; if (pFes_local[F_mapping[fid]].size() != 2) { pass_check = false; break; }
							//merge e0 and e1 -> e0; remove e1 from pFes_local 
							for (uint32_t k = 0; k < PE_npfs[eeid1].size(); k++) {
								uint32_t nfid = PE_npfs[eeid1][k]; std::replace(pFes_local[F_mapping[nfid]].begin(), pFes_local[F_mapping[nfid]].end(), eeid1, eeid0);
							}
							//update pPPs_local
							for (uint32_t k = 0; k < PF_npps_sudo[fid].size(); k++) {
								uint32_t pp = P_mapping[PF_npps_sudo[fid][k]];
								pPPs_local[pp].erase(std::remove(pPPs_local[pp].begin(), pPPs_local[pp].end(), F_mapping[fid]), pPPs_local[pp].end());
							}

							F_flag_local[F_mapping[fid]] = false;
							std::vector<uint32_t>().swap(PF_npps_sudo[fid]);
						}
						else {
							//shared common polyhedral
							std::vector<uint32_t> pps0, pps1, ecom_ps;
							for (uint32_t k = 0; k < PE_npfs[eeid0].size(); k++) pps0.insert(pps0.end(), PF_npps_sudo[PE_npfs[eeid0][k]].begin(), PF_npps_sudo[PE_npfs[eeid0][k]].end());
							std::sort(pps0.begin(), pps0.end()); pps0.erase(std::unique(pps0.begin(), pps0.end()), pps0.end());
							for (uint32_t k = 0; k < PE_npfs[eeid1].size(); k++) pps1.insert(pps1.end(), PF_npps_sudo[PE_npfs[eeid1][k]].begin(), PF_npps_sudo[PE_npfs[eeid1][k]].end());
							std::sort(pps1.begin(), pps1.end()); pps1.erase(std::unique(pps1.begin(), pps1.end()), pps1.end());
							std::set_intersection(pps0.begin(), pps0.end(), pps1.begin(), pps1.end(), std::back_inserter(ecom_ps));
							if (ecom_ps.size() != 1) { 
								pass_check = false; break; 
							}
							//cut the polyhedral into two
							uint32_t p = P_mapping[ecom_ps[0]];
							std::vector<std::vector<uint32_t>> pfs(pPPs_local[p].size()); std::vector<uint32_t> e_circle(3), ps0, ps1;
							for (uint32_t k = 0; k < pPPs_local[p].size(); k++) pfs[k] = pFes_local[pPPs_local[p][k]];
							e_circle[0] = e_map; e_circle[1] = eeid0; e_circle[2] = eeid1;
							cut_a_polyhedral(pPPs_local[p], pfs, e_circle, ps0, ps1);
							//assign an id to the new polyhedral
							uint32_t newId = -1;
							for (uint32_t k = 0; k < mpP_flag.size(); k++) 
								if (!mpP_flag[k]) { newId = k; mpP_flag[k] = true; P_mapping[k] = pPPs_local.size(); newIds.push_back(k); ps_total.push_back(k); break; }
							if (newId == -1) { pass_check = false; 
							break; }
							//update neighborhood info
							P_flag_local.push_back(true); pPPs_local[p] = ps0;  pPPs_local.push_back(ps1);
							for (uint32_t k = 0; k < ps1.size(); k++) std::replace(PF_npps_sudo[fs_total[ps1[k]]].begin(), PF_npps_sudo[fs_total[ps1[k]]].end(), ecom_ps[0], newId);
							//merge e0 and e1 -> e0; remove e1 from pFes_local 
							for (uint32_t k = 0; k < PE_npfs[eeid1].size(); k++) {
								uint32_t nfid = PE_npfs[eeid1][k]; std::replace(pFes_local[F_mapping[nfid]].begin(), pFes_local[F_mapping[nfid]].end(), eeid1, eeid0);
							}

							PPs_ring.push_back(newId);
						}
					}
					//////////////////////////////////////////////////////////////////////////////
					if (pass_check) {
						//pseudo delete singlets
						for (uint32_t j = 0; j < PPs_ring.size(); j++) {
							uint32_t p = P_mapping[PPs_ring[j]];
							if (pPPs_local[p].size() == 2) {//merge pf0 and pf1 -> pf0
								uint32_t pf0 = pPPs_local[p][0], pf1 = pPPs_local[p][1];
								P_flag_local[p] = false; F_flag_local[pf1] = false;
								if (pF_boundaryness_local[pf1]) pF_boundaryness_local[pf0] = true;
								//mpPs
								for (uint32_t k = 0; k < PF_npps_sudo[fs_total[pf1]].size(); k++) {
									uint32_t po = PF_npps_sudo[fs_total[pf1]][k];
									if (F_flag_local[pf0] && std::find(pPPs_local[P_mapping[po]].begin(), pPPs_local[P_mapping[po]].end(), pf0) == pPPs_local[P_mapping[po]].end())
										std::replace(pPPs_local[P_mapping[po]].begin(), pPPs_local[P_mapping[po]].end(), pf1, pf0);
									else pPPs_local[P_mapping[po]].erase(std::remove(pPPs_local[P_mapping[po]].begin(), pPPs_local[P_mapping[po]].end(), pf1), pPPs_local[P_mapping[po]].end());
								}
								//PF_npps
								PF_npps_sudo[fs_total[pf0]].insert(PF_npps_sudo[fs_total[pf0]].end(), PF_npps_sudo[fs_total[pf1]].begin(), PF_npps_sudo[fs_total[pf1]].end());
								std::vector<uint32_t> nps_temp; nps_temp.reserve(2);
								for (uint32_t pp = 0; pp < PF_npps_sudo[fs_total[pf0]].size(); pp++) {
									uint32_t po = PF_npps_sudo[fs_total[pf0]][pp];
									if (P_flag_local[P_mapping[po]]) nps_temp.push_back(po);
								}
								std::sort(nps_temp.begin(), nps_temp.end()); nps_temp.erase(std::unique(nps_temp.begin(), nps_temp.end()), nps_temp.end());
								nps_temp.swap(PF_npps_sudo[fs_total[pf0]]);
								std::vector<uint32_t>().swap(PF_npps_sudo[fs_total[pf1]]);
							}
							else if (pPPs_local[p].size() == 3) {//one face pf0 is composed by another two: pf1 and pf2
								uint32_t pf0 = pPPs_local[p][0];
								uint32_t pf12[2]; pf12[0] = pPPs_local[p][1]; pf12[1] = pPPs_local[p][2];
								if (pFes_local[pf0].size() < pFes_local[pf12[0]].size()) std::swap(pf0, pf12[0]);
								if (pFes_local[pf0].size() < pFes_local[pf12[1]].size()) std::swap(pf0, pf12[1]);

								if (pFes_local[pf0].size() + 2 != (pFes_local[pf12[0]].size() + pFes_local[pf12[1]].size())) {
									continue;
								}

								if (PF_npps_sudo[fs_total[pf0]].size() == 2 && PF_npps_sudo[fs_total[pf12[0]]].size() == 2) {
									if (PF_npps_sudo[fs_total[pf0]][0] > PF_npps_sudo[fs_total[pf0]][1]) std::swap(PF_npps_sudo[fs_total[pf0]][0], PF_npps_sudo[fs_total[pf0]][1]);
									if (PF_npps_sudo[fs_total[pf12[0]]][0] > PF_npps_sudo[fs_total[pf12[0]]][1]) std::swap(PF_npps_sudo[fs_total[pf12[0]]][0], PF_npps_sudo[fs_total[pf12[0]]][1]);
									
									if (PF_npps_sudo[fs_total[pf0]][0] == PF_npps_sudo[fs_total[pf12[0]]][0] && PF_npps_sudo[fs_total[pf0]][1] == PF_npps_sudo[fs_total[pf12[0]]][1]) std::swap(pf12[1], pf0);
								}
								if (PF_npps_sudo[fs_total[pf0]].size() == 2 && PF_npps_sudo[fs_total[pf12[1]]].size() == 2) {
									if (PF_npps_sudo[fs_total[pf0]][0] > PF_npps_sudo[fs_total[pf0]][1]) std::swap(PF_npps_sudo[fs_total[pf0]][0], PF_npps_sudo[fs_total[pf0]][1]);
									if (PF_npps_sudo[fs_total[pf12[1]]][0] > PF_npps_sudo[fs_total[pf12[1]]][1]) std::swap(PF_npps_sudo[fs_total[pf12[1]]][0], PF_npps_sudo[fs_total[pf12[1]]][1]);
									
									if (PF_npps_sudo[fs_total[pf0]][0] == PF_npps_sudo[fs_total[pf12[1]]][0] && PF_npps_sudo[fs_total[pf0]][1] == PF_npps_sudo[fs_total[pf12[1]]][1]) std::swap(pf12[0], pf0);
								}

								P_flag_local[p] = false; F_flag_local[pf0] = false;

								if (pF_boundaryness_local[pf0] && pF_boundaryness_local[pf12[0]]) { F_flag_local[pf12[0]] = false; }
								if (pF_boundaryness_local[pf0] && pF_boundaryness_local[pf12[1]]) { F_flag_local[pf12[1]] = false; }

								if (pF_boundaryness_local[pf0]) { pF_boundaryness_local[pf12[0]] = pF_boundaryness_local[pf12[1]] = true; }

								for (uint32_t pp = 0; pp < PF_npps_sudo[fs_total[pf0]].size(); pp++) {
									uint32_t po = PF_npps_sudo[fs_total[pf0]][pp];
									{
										if(F_flag_local[pf12[0]] && std::find(pPPs_local[P_mapping[po]].begin(), pPPs_local[P_mapping[po]].end(), pf12[0]) == pPPs_local[P_mapping[po]].end())
											std::replace(pPPs_local[P_mapping[po]].begin(), pPPs_local[P_mapping[po]].end(), pf0, pf12[0]); 
										else pPPs_local[P_mapping[po]].erase(std::remove(pPPs_local[P_mapping[po]].begin(), pPPs_local[P_mapping[po]].end(), pf0), pPPs_local[P_mapping[po]].end());
										if(F_flag_local[pf12[1]] && std::find(pPPs_local[P_mapping[po]].begin(), pPPs_local[P_mapping[po]].end(), pf12[1]) == pPPs_local[P_mapping[po]].end()) 
											pPPs_local[P_mapping[po]].push_back(pf12[1]);
									}
								}
								//PF_npps
								for (uint32_t k = 0; k < 2; k++) {
									PF_npps_sudo[fs_total[pf12[k]]].insert(PF_npps_sudo[fs_total[pf12[k]]].end(), PF_npps_sudo[fs_total[pf0]].begin(), PF_npps_sudo[fs_total[pf0]].end());
									std::vector<uint32_t> nps_temp; nps_temp.reserve(2);
									for (uint32_t pp = 0; pp < PF_npps_sudo[fs_total[pf12[k]]].size(); pp++) {
										uint32_t po = PF_npps_sudo[fs_total[pf12[k]]][pp];
										if (P_flag_local[P_mapping[po]]) nps_temp.push_back(po);
									}
									std::sort(nps_temp.begin(), nps_temp.end()); nps_temp.erase(std::unique(nps_temp.begin(), nps_temp.end()), nps_temp.end());
									nps_temp.swap(PF_npps_sudo[fs_total[pf12[k]]]);
								}
								std::vector<uint32_t>().swap(PF_npps_sudo[fs_total[pf0]]);
							}
						}
						//pseudo new PFes, PPs
						std::vector<uint32_t> F_mapping_(pFes_local.size(), INVALID_F);
						std::vector<std::vector<uint32_t>> pFes_local_, pPPs_local_;
						pFes_local_.reserve(pFes_local.size()); pPPs_local_.reserve(pPPs_local.size());
						std::vector<bool> pF_boundaryness_local_; pF_boundaryness_local_.reserve(pF_boundaryness_local.size());
						for (uint32_t j = 0; j < F_flag_local.size(); j++)
							if (F_flag_local[j]) {
								pFes_local_.push_back(pFes_local[j]);
								F_mapping_[j] = pFes_local_.size() - 1;
								pF_boundaryness_local_.push_back(pF_boundaryness_local[j]);
							}
						for (uint32_t j = 0; j < P_flag_local.size(); j++)
							if (P_flag_local[j]) {
								std::vector<uint32_t> temp_fs;
								for (uint32_t k = 0; k < pPPs_local[j].size(); k++) temp_fs.push_back(F_mapping_[pPPs_local[j][k]]);
								pPPs_local_.push_back(temp_fs);
							}

						pass_check = check_aft(pFes_local_, pF_boundaryness_local_, pPPs_local_);
					}
					//update back
					for (uint32_t j = 0; j < fs_total.size(); j++) { F_mapping[fs_total[j]] = INVALID_F; PF_npps_sudo[fs_total[j]].clear(); }
					for (uint32_t j = 0; j < ps_total.size(); j++) P_mapping[ps_total[j]] = INVALID_P;
					for (uint32_t j = 0; j < newIds.size(); j++) mpP_flag[newIds[j]] = false;
					for (uint32_t j = 0; j < Reverse_pV_map[v1_map].size(); j++)
						pV_map[Reverse_pV_map[v1_map][j]] = v1_map;
					//compare
					if (!pass_check)
						return false;
					return true;
				};
				
				if (!check_topology()) {
					std::get<3>(e) *= 2;
					Es_nonmanifold.push_back(e);
					continue;
				}

				std::get<4>(mpEs[e_map]) = Edge_tag::R;
				for (uint32_t j = 0; j < Reverse_pE_map[e_map].size(); j++) pE_map[Reverse_pE_map[e_map][j]] = INVALID_E;

				if (doublets) {
					v0_mapR = Reverse_pV_map[v0_map];
					v1_mapR = Reverse_pV_map[v1_map];
				}

				topology = true;
				//update mappings merge v1_map and v0_map -> v0_map
				for (uint32_t j = 0; j < Reverse_pV_map[v1_map].size(); j++) pV_map[Reverse_pV_map[v1_map][j]] = v0_map;
				Reverse_pV_map[v0_map].insert(Reverse_pV_map[v0_map].end(), Reverse_pV_map[v1_map].begin(), Reverse_pV_map[v1_map].end());
				if (mV_B_flag[v1_map])  mV_B_flag[v0_map] = true;
				std::vector<uint32_t>().swap(Reverse_pV_map[v1_map]);

				PV_npfs[v0_map].insert(PV_npfs[v0_map].end(), PV_npfs[v1_map].begin(), PV_npfs[v1_map].end());
				std::sort(PV_npfs[v0_map].begin(), PV_npfs[v0_map].end()); PV_npfs[v0_map].erase(std::unique(PV_npfs[v0_map].begin(), PV_npfs[v0_map].end()), PV_npfs[v0_map].end());
				std::vector<uint32_t>().swap(PV_npfs[v1_map]);
				for (uint32_t j = 0; j < PV_npfs[v0_map].size(); j++) {
					uint32_t fid = PV_npfs[v0_map][j];
					if (std::find(mpFvs[fid].begin(), mpFvs[fid].end(), v0_map) != mpFvs[fid].end())
						mpFvs[fid].erase(std::remove(mpFvs[fid].begin(), mpFvs[fid].end(), v1_map), mpFvs[fid].end());
					else std::replace(mpFvs[fid].begin(), mpFvs[fid].end(), v1_map, v0_map);
				}

				std::vector<uint32_t> PPs_ring_set;
				for (uint32_t j = 0; j < PE_npfs[e_map].size(); j++) PPs_ring_set.insert(PPs_ring_set.end(), PF_npps[PE_npfs[e_map][j]].begin(), PF_npps[PE_npfs[e_map][j]].end());
				std::sort(PPs_ring_set.begin(), PPs_ring_set.end()); PPs_ring_set.erase(std::unique(PPs_ring_set.begin(), PPs_ring_set.end()), PPs_ring_set.end());
//////////////////////////////////////////////////////////////////////////////
				for (uint32_t j = 0; j < PE_npfs[e_map].size(); j++) { uint32_t fid = PE_npfs[e_map][j]; mpFes[fid].erase(std::remove(mpFes[fid].begin(), mpFes[fid].end(), e_map), mpFes[fid].end()); }
				//ffs_set
				std::vector<uint32_t> ffs_set; ffs_set.reserve(PPs_ring_set.size()*3);
				for (uint32_t j = 0; j < PPs_ring_set.size(); j++) for (uint32_t k = 0; k < mpPs[PPs_ring_set[j]].size(); k++) mpF_flag[mpPs[PPs_ring_set[j]][k]] = true;
				for (uint32_t j = 0; j < PPs_ring_set.size(); j++) for (uint32_t k = 0; k < mpPs[PPs_ring_set[j]].size(); k++) 
					if (mpF_flag[mpPs[PPs_ring_set[j]][k]]) { ffs_set.push_back(mpPs[PPs_ring_set[j]][k]); mpF_flag[mpPs[PPs_ring_set[j]][k]] = false; }
				//ees_tuples
				std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> ees_tuples; ees_tuples.reserve(ffs_set.size() * 3);
				for (uint32_t j = 0; j < ffs_set.size(); j++) for (uint32_t k = 0; k < mpFes[ffs_set[j]].size(); k++) mpE_flag[mpFes[ffs_set[j]][k]] = true;
				for (uint32_t j = 0; j < ffs_set.size(); j++) for (uint32_t k = 0; k < mpFes[ffs_set[j]].size(); k++) 
					if (mpE_flag[mpFes[ffs_set[j]][k]]) {
						uint32_t eeid = mpFes[ffs_set[j]][k]; mpE_flag[eeid] = false;
						uint32_t evid0 = pV_map[std::get<0>(mpEs[eeid])], evid1 = pV_map[std::get<1>(mpEs[eeid])];
						if (evid0 > evid1) std::swap(evid0, evid1);
						ees_tuples.push_back(std::make_tuple(evid0, evid1, eeid));
					}
				std::sort(ees_tuples.begin(), ees_tuples.end());
				//ees_sets: collapsing pairs
				std::vector<std::vector<uint32_t>> ees_sets; std::vector<uint32_t> es_set;
				for (uint32_t j = 0; j < ees_tuples.size(); j++) {
					if (j == 0 || (j != 0 && (std::get<0>(ees_tuples[j]) != std::get<0>(ees_tuples[j - 1]) || std::get<1>(ees_tuples[j]) != std::get<1>(ees_tuples[j - 1])))) {
						if (j != 0 && es_set.size() == 2) ees_sets.push_back(es_set);
						es_set.clear(); es_set.reserve(2); es_set.push_back(std::get<2>(ees_tuples[j]));
					}
					else es_set.push_back(std::get<2>(ees_tuples[j]));
					if (j + 1 == ees_tuples.size() && es_set.size() == 2) ees_sets.push_back(es_set);
				}
				for (uint32_t j = 0; j < ees_sets.size(); j++) {
				//judge if they share a common face
					uint32_t eeid0 = ees_sets[j][0], eeid1 = ees_sets[j][1];
					std::sort(PE_npfs[eeid0].begin(), PE_npfs[eeid0].end()); std::sort(PE_npfs[eeid1].begin(), PE_npfs[eeid1].end());
					std::vector<uint32_t> ecom_fs;
					std::set_intersection(PE_npfs[eeid0].begin(), PE_npfs[eeid0].end(), PE_npfs[eeid1].begin(), PE_npfs[eeid1].end(), std::back_inserter(ecom_fs));
					if (ecom_fs.size() > 1) { std::cout << "ERROR here!" << endl; system("PAUSE"); }
					else if (ecom_fs.size()) {
						uint32_t fid = ecom_fs[0]; if (mpFes[fid].size() != 2) { std::cout << "ERROR here!" << endl; system("PAUSE"); }
						//update mPPs
						for (uint32_t k = 0; k < PF_npps[fid].size(); k++) {
							uint32_t pp = PF_npps[fid][k]; mpPs[pp].erase(std::remove(mpPs[pp].begin(), mpPs[pp].end(), fid), mpPs[pp].end());
						}
						//update PV_npfs
						for (uint32_t k = 0; k < mpFvs[fid].size(); k++)
							PV_npfs[mpFvs[fid][k]].erase(std::remove(PV_npfs[mpFvs[fid][k]].begin(), PV_npfs[mpFvs[fid][k]].end(), fid), PV_npfs[mpFvs[fid][k]].end());
						//update Reverse_pE_map, boundaryness
						uint32_t e0 = mpFes[fid][0], e1 = mpFes[fid][1];
						Reverse_pE_map[e0].insert(Reverse_pE_map[e0].end(), Reverse_pE_map[e1].begin(), Reverse_pE_map[e1].end());
						for (uint32_t m = 0; m < Reverse_pE_map[e1].size(); m++) pE_map[Reverse_pE_map[e1][m]] = e0;

						if (std::get<2>(mpEs[e0]) || std::get<2>(mpEs[e1]))
							for (uint32_t m = 0; m < Reverse_pE_map[e0].size(); m++) std::get<2>(mpEs[Reverse_pE_map[e0][m]]) = true;
						//mPFes
						for (uint32_t k = 0; k < PE_npfs[e1].size(); k++) {
							uint32_t nfid = PE_npfs[e1][k];
							if (std::find(mpFes[nfid].begin(), mpFes[nfid].end(), e0) != mpFes[nfid].end())
								mpFes[nfid].erase(std::remove(mpFes[nfid].begin(), mpFes[nfid].end(), e1), mpFes[nfid].end());
							else std::replace(mpFes[nfid].begin(), mpFes[nfid].end(), e1, e0);
						}
						//PE_npfs
						PE_npfs[e0].insert(PE_npfs[e0].end(), PE_npfs[e1].begin(), PE_npfs[e1].end());
						std::vector<uint32_t> temp_fs; temp_fs.reserve(PE_npfs[e0].size());
						for (uint32_t k = 0; k < PE_npfs[e0].size(); k++)
							if (PE_npfs[e0][k] != fid) temp_fs.push_back(PE_npfs[e0][k]);
						std::sort(temp_fs.begin(), temp_fs.end()); temp_fs.erase(std::unique(temp_fs.begin(), temp_fs.end()), temp_fs.end());
						temp_fs.swap(PE_npfs[e0]);

						//delete e1
						std::vector<uint32_t>().swap(Reverse_pE_map[e1]);
						std::vector<uint32_t>().swap(PE_npfs[e1]);
						//delete fid
						std::vector<uint32_t>().swap(mpFvs[fid]);
						std::vector<uint32_t>().swap(mpFes[fid]);
						std::vector<uint32_t>().swap(PF_npps[fid]);
						pF_map[fid] = INVALID_F;
					}
					else {
						//shared common polyhedral
						std::vector<uint32_t> pps0, pps1, ecom_ps;
						for (uint32_t k = 0; k < PE_npfs[eeid0].size(); k++) pps0.insert(pps0.end(), PF_npps[PE_npfs[eeid0][k]].begin(), PF_npps[PE_npfs[eeid0][k]].end());
						std::sort(pps0.begin(), pps0.end()); pps0.erase(std::unique(pps0.begin(), pps0.end()), pps0.end());
						for (uint32_t k = 0; k < PE_npfs[eeid1].size(); k++) pps1.insert(pps1.end(), PF_npps[PE_npfs[eeid1][k]].begin(), PF_npps[PE_npfs[eeid1][k]].end());
						std::sort(pps1.begin(), pps1.end()); pps1.erase(std::unique(pps1.begin(), pps1.end()), pps1.end());
						std::set_intersection(pps0.begin(), pps0.end(), pps1.begin(), pps1.end(), std::back_inserter(ecom_ps));
						if (ecom_ps.size() != 1) { std::cout << "ERROR here!" << endl; system("PAUSE"); }
						//cut the polyhedral into two
						std::vector<std::vector<uint32_t>> pfs(mpPs[ecom_ps[0]].size()); std::vector<uint32_t> e_circle(3), ps0, ps1;
						for (uint32_t k = 0; k < mpPs[ecom_ps[0]].size(); k++) pfs[k] = mpFes[mpPs[ecom_ps[0]][k]];
						e_circle[0] = e_map; e_circle[1] = eeid0; e_circle[2] = eeid1;
						cut_a_polyhedral(mpPs[ecom_ps[0]], pfs, e_circle, ps0, ps1);
						//assign an id to the new polyhedral
						uint32_t newId = -1;
						for (uint32_t k = 0; k < mpP_flag.size(); k++) if (!mpP_flag[k]) {newId = k; mpP_flag[newId] = true; break; }
						//update neighborhood info
						mpPs[ecom_ps[0]] = ps0;  mpPs[newId] = ps1;
						for (uint32_t k = 0; k < ps1.size(); k++) std::replace(PF_npps[ps1[k]].begin(), PF_npps[ps1[k]].end(), ecom_ps[0], newId);
						PPs_ring_set.push_back(newId);

						//update Reverse_pE_map, boundaryness
						Reverse_pE_map[eeid0].insert(Reverse_pE_map[eeid0].end(), Reverse_pE_map[eeid1].begin(), Reverse_pE_map[eeid1].end());
						for (uint32_t m = 0; m < Reverse_pE_map[eeid1].size(); m++) pE_map[Reverse_pE_map[eeid1][m]] = eeid0;
						if (std::get<2>(mpEs[eeid0]) || std::get<2>(mpEs[eeid1]))
							for (uint32_t m = 0; m < Reverse_pE_map[eeid0].size(); m++) std::get<2>(mpEs[Reverse_pE_map[eeid0][m]]) = true;
						//mPFes
						for (uint32_t k = 0; k < PE_npfs[eeid1].size(); k++) {
							uint32_t nfid = PE_npfs[eeid1][k];
							if (std::find(mpFes[nfid].begin(), mpFes[nfid].end(), eeid0) != mpFes[nfid].end())
								mpFes[nfid].erase(std::remove(mpFes[nfid].begin(), mpFes[nfid].end(), eeid1), mpFes[nfid].end());
							else std::replace(mpFes[nfid].begin(), mpFes[nfid].end(), eeid1, eeid0);
						}
						//PE_npfs
						PE_npfs[eeid0].insert(PE_npfs[eeid0].end(), PE_npfs[eeid1].begin(), PE_npfs[eeid1].end());
						std::sort(PE_npfs[eeid0].begin(), PE_npfs[eeid0].end()); 
						PE_npfs[eeid0].erase(std::unique(PE_npfs[eeid0].begin(), PE_npfs[eeid0].end()), PE_npfs[eeid0].end());
						//delete eeid
						std::vector<uint32_t>().swap(Reverse_pE_map[eeid1]);
						std::vector<uint32_t>().swap(PE_npfs[eeid1]);
					}
				}
//////////////////////////////////////////////////////////////////////////////
				degenerate_polyhedra(PPs_ring_set);

				PPs_ring_set.clear();
				for (uint32_t j = 0; j < PV_npfs[v0_map].size(); j++) PPs_ring_set.insert(PPs_ring_set.end(), PF_npps[PV_npfs[v0_map][j]].begin(), PF_npps[PV_npfs[v0_map][j]].end());
				std::sort(PPs_ring_set.begin(), PPs_ring_set.end()); PPs_ring_set.erase(std::unique(PPs_ring_set.begin(), PPs_ring_set.end()), PPs_ring_set.end());

				for (uint32_t j = 0; j < PPs_ring_set.size(); j++)
					if(mpP_flag[PPs_ring_set[j]]) pf_test.insert(pf_test.end(), mpPs[PPs_ring_set[j]].begin(), mpPs[PPs_ring_set[j]].end());
				std::sort(pf_test.begin(), pf_test.end()); pf_test.erase(std::unique(pf_test.begin(), pf_test.end()), pf_test.end());
			}

			topology = true;
			//all faces of the new polyhedral, check if merge or not, and their simplicity
			if (pf_test.size()) {
				std::vector<uint32_t> PPs_ring_set;

				edge_fuse_polygons(pf_test, PPs_ring_set);
				//remove degenerate polyhedrals.
				if (PPs_ring_set.size()) {
					std::sort(PPs_ring_set.begin(), PPs_ring_set.end()); PPs_ring_set.erase(std::unique(PPs_ring_set.begin(), PPs_ring_set.end()), PPs_ring_set.end());
					degenerate_polyhedra(PPs_ring_set);
				}
			}


			if (std::get<4>(e) == Edge_tag::R) {
				//update V

				Vector3f posy; posy.setZero();
				Quaternion rosy = Quaternion::Zero(), q; q = mQ_copy.col(Reverse_pV_map[v0_map][0]);
				for (auto rvid : Reverse_pV_map[v0_map]) {
					posy += mO_copy.col(rvid);

					rosy = (rosy + Quaternion::applyRotation(mQ_copy.col(rvid), q)).normalized();
				}
				posy /= Reverse_pV_map[v0_map].size();

				mV_tag.col(v0_map) = posy;
				newQ.col(v0_map) = rosy.normalized();

				pV_npes[v0_map].insert(pV_npes[v0_map].end(), pV_npes[v1_map].begin(), pV_npes[v1_map].end());
				std::vector<uint32_t>().swap(pV_npes[v1_map]);
				std::vector<uint32_t> nes_temp;
				for (auto eid : pV_npes[v0_map]) {
					if (pE_map[eid] == INVALID_E) continue;
					nes_temp.push_back(pE_map[eid]);
				}
				std::sort(nes_temp.begin(), nes_temp.end()); nes_temp.erase(std::unique(nes_temp.begin(), nes_temp.end()), nes_temp.end());
				pV_npes[v0_map] = nes_temp;

				if (!re_color) break;

				for (auto eid : nes_temp) {

					tuple_E e_;
					uint32_t v0_ = pV_map[std::get<0>(mpEs[eid])], v1_ = pV_map[std::get<1>(mpEs[eid])];
					const Vector3f o0 = mV_tag.col(v0_), o1 = mV_tag.col(v1_);
					Float energy = (o0 - o1).norm();
					std::get<0>(e_) = v0_; std::get<1>(e_) = v1_; std::get<2>(e_) = std::get<2>(mpEs[eid]);
					std::get<5>(e_) = eid; std::get<7>(e_) = ++E_TimeStamp[eid];
					Quaternion q0, q1, n0, n1; std::vector<uint32_t> votes(4, 0);
					for (auto eo : Reverse_pE_map[eid]) {
						uint32_t v0_o = std::get<0>(mpEs[eo]), v1_o = std::get<1>(mpEs[eo]), v0_om = pV_map[v0_o], v1_om = pV_map[v1_o];
						q0 = newQ.col(v0_o), q1 = newQ.col(v1_o);
						Quaternion q_next = Quaternion::applyRotation(q1, q0);
						std::pair<int, Float> a_pair = assignColorWeighted3D(mV_tag.col(v0_om), q0, mV_tag.col(v1_om), q_next, mScale, mInvScale);
						std::get<3>(e_) = a_pair.second;

						votes[std::get<0>(a_pair)]++; std::get<6>(e_) = std::get<1>(a_pair);
					}
					uint32_t pos = 0, num = votes[pos];
					for (uint32_t k = 1; k < 4; k++) if (votes[k] > num) { num = votes[k]; pos = k; }
					std::get<4>(e_) = pos;

					if (std::get<4>(e_) == Edge_tag::B) continue;
					Es_nonmanifold.push_back(e_);
				}
			}

			break;
		}

		for (uint32_t i = 0; i < Es_nonmanifold.size(); i++)
			Es_red.push(Es_nonmanifold[i]);

		if ((!topology || Es_red.empty())) {
			topology = true;
			once = true;
			std::vector<uint32_t> PPs_ring_set;

			std::vector<uint32_t> all_fs; all_fs.reserve(mpFes.size());
			for (uint32_t i = 0; i < mpFes.size(); ++i) if (mpFes[i].size()) all_fs.push_back(i);
			
			while (edge_fuse_polygons(all_fs, PPs_ring_set));
			std::sort(PPs_ring_set.begin(), PPs_ring_set.end()); PPs_ring_set.erase(std::unique(PPs_ring_set.begin(), PPs_ring_set.end()), PPs_ring_set.end());
			degenerate_polyhedra(PPs_ring_set); if (!PPs_ring_set.size()) once=false;

			std::fill(E_TimeStamp.begin(), E_TimeStamp.end(), 0);
			while (Es_red.size()) {
				tuple_E e = Es_red.top(); Es_red.pop();
				std::get<4>(mpEs[std::get<5>(e)]) = std::get<4>(e);
				E_TimeStamp[std::get<5>(e)]++;
				std::get<7>(mpEs[std::get<5>(e)]) = E_TimeStamp[std::get<5>(e)];
			}
			tuple_E ei;
			for (uint32_t i = 0; i < mpFes.size(); ++i) {
				if (mpFes[i].size()) {
					for (uint32_t j = 0; j < mpFes[i].size(); ++j) {
						if (mpE_flag[mpFes[i][j]]) continue;
						if (std::get<4>(mpEs[mpFes[i][j]]) == Edge_tag::R) {
							Es_red.push(mpEs[mpFes[i][j]]); mpE_flag[mpFes[i][j]] = true;
							std::get<4>(mpEs[mpFes[i][j]]) = Edge_tag::B;
						}
						if (std::get<4>(mpEs[mpFes[i][j]]) == Edge_tag::D || std::get<4>(mpEs[mpFes[i][j]]) == Edge_tag::H) {
							Es_red.push(mpEs[mpFes[i][j]]); mpE_flag[mpFes[i][j]] = true;
							std::get<4>(mpEs[mpFes[i][j]]) = Edge_tag::B;
						}
					}
				}
			}
			std::fill(mpE_flag.begin(), mpE_flag.end(), false);

			if (!once) {
				std::priority_queue<tuple_E, std::vector<tuple_E>, LessThan> Es_red_sudo;
				while (Es_red.size())
				{
					tuple_E e = Es_red.top(); Es_red.pop();
					Es_red_sudo.push(e);
					if (std::get<4>(e) == Edge_tag::R) {
						uint32_t heid; face_swap(std::get<5>(e), heid);
						if (heid != INVALID_E) {
							e = mpEs[heid]; std::get<4>(e) = Edge_tag::H;
							Es_red_sudo.push(e);
							once = true;
						}
					}
				}
				Es_red_sudo.swap(Es_red);
			}
			
			uint32_t cur_num = 0;
			for (auto fes : mpFes) if (fes.size()) cur_num++;
			if (left_NUM == cur_num) {
				once = false;
			}
			else left_NUM = cur_num;
			//
			if (!once) {
				topology = false; 
				left_NUM = 0;
				break;
			}
		}
	}
	std::set<tuple_E> lefted_es;
	while (Es_red.size())
	{
		tuple_E e = Es_red.top(); Es_red.pop();
		std::get<0>(e) = pV_map[std::get<0>(e)];
		std::get<1>(e) = pV_map[std::get<1>(e)];
		if (std::get<0>(e) > std::get<1>(e)) std::swap(std::get<0>(e), std::get<1>(e));
		std::get<4>(mpEs[std::get<5>(e)]) = std::get<4>(e);
		lefted_es.insert(e);
	}
}

bool MultiResolutionHierarchy::split_long_edge3D(vector<uint32_t> &ledges) {

	uint32_t Nvo = mO_copy.cols(), Nv= Nvo + ledges.size();
	uint32_t Ne = mpEs.size(), Nf = mpFvs.size();
	mQ_copy.conservativeResize(4, Nv);
	mO_copy.conservativeResize(3, Nv);
	if(Qquadric)
		Quadric_copy.resize(Nv);

	vector<short>().swap(mpV_flag); mpV_flag.resize(mO_copy.cols(), false);
	vector<short>().swap(mpE_flag); mpE_flag.resize(mpEs.size(), false);
	PV_npfs_sudo.clear(); PV_npfs_sudo.resize(Nv);
	PV_npvs.resize(Nv);

	for (auto e_pos : ledges) {

		tuple_E e = mpEs[e_pos]; uint32_t eid = get<5>(e);

		uint32_t v0 = std::get<0>(e), v1 = std::get<1>(e);
		Quaternion q0 = mQ_copy.col(v0), q1 = Quaternion::applyRotation(mQ_copy.col(v1), q0);
		//compute middle point info
		Quaternion qn = (q0 + q1).normalized();
		Vector3f gn = (mO_copy.col(v0) + mO_copy.col(v1)) * 0.5;
		uint32_t vn = Nvo++;
		uint32_t en = Ne++;

		if (Qquadric) {
			Quadric_copy[vn] = Quadric_copy[v0] + Quadric_copy[v1];
			Quadric_copy[vn].getMinimum(gn);
		}

		mQ_copy.col(vn) = qn;
		mO_copy.col(vn) = gn;

		replace(PV_npvs[v0].begin(), PV_npvs[v0].end(), v1, vn);
		replace(PV_npvs[v1].begin(), PV_npvs[v1].end(), v0, vn);
		PV_npvs[vn].push_back(v0);
		PV_npvs[vn].push_back(v1);

		tuple_E e0, e1; std::tuple<short, Float, Vector3f> a_posy;
		get<0>(e0) = v0;
		get<1>(e0) = vn;
		get<2>(e0) = get<2>(e);
		a_posy = posy3D_completeInfo(mO_copy.col(v0), q0, gn, qn, mScale, mInvScale);
		get<3>(e0) = get<1>(a_posy);
		get<4>(e0) = get<0>(a_posy);
		get<5>(e0) = eid;
		get<6>(e0) = 0;
		get<7>(e0) = 0;

		get<0>(e1) = v1;
		get<1>(e1) = vn;
		get<2>(e1) = get<2>(e);
		a_posy = posy3D_completeInfo(mO_copy.col(v1), q1, gn, qn, mScale, mInvScale);
		get<4>(e1) = get<0>(a_posy);
		get<3>(e1) = get<1>(a_posy);
		get<5>(e1) = en;
		get<6>(e1) = 0;
		get<7>(e1) = 0;

		if (Ne > mpEs.size()) {
			mpEs.resize(std::max(Ne, (uint32_t)mpEs.size() * 2));
			PE_npfs.resize(mpEs.size());
			mpV_flag.resize(mpEs.size(), false);
			mpE_flag.resize(mpEs.size(), false);
		}
		mpEs[eid] = e0; mpEs[en] = e1;
		PE_npfs[en] = PE_npfs[eid];

		vector<uint32_t> fids, vnpos, vposs; vector<vector<uint32_t>> vsss, vsss_;
		for (uint32_t j = 0; j < PE_npfs[eid].size(); j++) {
			auto fid = PE_npfs[eid][j];
			//update mpFvs, mpFes
			vector<uint32_t> fvs;
			for (uint32_t k = 0; k < mpFvs[fid].size(); k++) {
				auto v_cur = mpFvs[fid][k], v_aft = mpFvs[fid][(k + 1) % mpFvs[fid].size()];
				fvs.push_back(v_cur);
				if ((v_cur == v0 && v_aft == v1) || (v_cur == v1 && v_aft == v0)){
					fvs.push_back(vn);
					vnpos.push_back(k+1);
				}
			}
			vsss.push_back(mpFvs[fid]);
			mpFvs[fid] = fvs;
			mpFes[fid].push_back(en);
			vsss_.push_back(mpFvs[fid]);

			//collect v
			for (uint32_t k = 0; k < mpFvs[fid].size(); k++) {
				auto vid = mpFvs[fid][k];
				if (vid == vn) continue;
				if (vid == v0 || vid == v1) continue;
				Quaternion q_test = Quaternion::applyRotation(mQ_copy.col(vid), qn);
				a_posy = posy3D_completeInfo(gn, qn, mO_copy.col(vid), q_test, mScale, mInvScale);

				fids.push_back(fid); vposs.push_back(k);
				break;
			}
		}
		if (!fids.size()) continue;
		//update es
		for (uint32_t j = 0; j < fids.size(); j++) {
			auto fid = fids[j];
			std::vector<uint32_t> &es = mpFes[fid], es_temp, &vs = mpFvs[fid];
			auto vid = vs[vposs[j]];
			//orient es direction
			for (uint32_t k = 0; k < vs.size(); k++) {
				for (auto e : es) {
					if ((get<0>(mpEs[e]) == vs[k] || get<0>(mpEs[e]) == vs[(k + 1) % vs.size()]) &&
						(get<1>(mpEs[e]) == vs[k] || get<1>(mpEs[e]) == vs[(k + 1) % vs.size()])) {
						es_temp.push_back(e); break;
					}
				}
			}
			es = es_temp;

			vector<vector<uint32_t>> nes2(2), nvs2(2);
			int32_t start = vnpos[j], end = vposs[j];
			
			if (find(PV_npvs[vid].begin(), PV_npvs[vid].end(), vs[end]) != PV_npvs[vid].end()) continue;

			en = Ne++;
			uint32_t fn = Nf++;

			tuple_E new_e;
			get<0>(new_e) = vs[start];
			get<1>(new_e) = vs[end];
			get<2>(new_e) = mpF_boundary_flag[fid];
			Quaternion q_test = Quaternion::applyRotation(mQ_copy.col(vid), qn);
			a_posy = posy3D_completeInfo(gn, qn, mO_copy.col(vid), q_test, mScale, mInvScale);
			get<4>(new_e) = get<0>(a_posy);
			get<3>(new_e) = get<1>(a_posy);
			get<5>(new_e) = en;
			get<6>(new_e) = 0;
			get<7>(new_e) = 0;
			
			if (Ne > mpEs.size()) {
				mpEs.resize(std::max(Ne, (uint32_t)mpEs.size() * 2));
				PE_npfs.resize(mpEs.size());
				mpE_flag.resize(mpEs.size(), false);
			}
			mpEs[en] = new_e;


			for (uint32_t m = 0; m < 2; m++) {
				if (m == 1) std::swap(start, end);
				int32_t length = (end - start + vs.size() + 1) % vs.size();
				std::vector<uint32_t> nes(length), nvs(length);
				for (uint32_t k = 0; k < length; k++) {
					nvs[k] = vs[(start + k) % vs.size()];
					if (k + 1 == length) nes[k] = std::get<5>(new_e);
					else nes[k] = es[(start + k) % es.size()];
				}

				vector<vector<uint32_t>> fvs_simple_test, fes_simple_test;
				fvs_simple_test.push_back(nvs);
				fes_simple_test.push_back(nes);
				vector<uint32_t> pvs, pes, vs_disgard, es_disgard;
				if (!simple_polygon_3D_v2(fvs_simple_test, fes_simple_test, pvs, pes, vs_disgard, es_disgard, false))
					;

				nes2[m] = nes; nvs2[m] = nvs;
			}
			bool this_one = true;
			//simple polyhedral
			for (auto pid : PF_npps[fid]) {
				std::vector<std::vector<uint32_t>> pfes;
				for (uint32_t n = 0; n < mpPs[pid].size(); n++) {
					if(mpPs[pid][n] != fid) pfes.push_back(mpFes[mpPs[pid][n]]);
					else {
						pfes.push_back(nes2[0]);
						pfes.push_back(nes2[1]);
					}
				}
				if (!simple_polyhedral_v3(pfes)) { this_one = false; break; }
			}

			if (!this_one) {
				tuple_E e_ran;
				mpEs[en] = e_ran;
				Ne--;
				Nf--;
				continue;
			}
			PV_npvs[vs[start]].push_back(vs[end]);
			PV_npvs[vs[end]].push_back(vs[start]);


			PE_npfs[en].push_back(fid);
			PE_npfs[en].push_back(fn);

			for (auto eid_ : nes2[1])
				if (eid_ != en && find(PE_npfs[eid_].begin(), PE_npfs[eid_].end(), fid) != PE_npfs[eid_].end())
					replace(PE_npfs[eid_].begin(), PE_npfs[eid_].end(), fid, fn);

			for (auto pid : PF_npps[fid]) mpPs[pid].push_back(fn);

			if (Nf > mpFvs.size()) {
				mpFvs.resize(std::max(Nf, (uint32_t)mpFvs.size() * 2));
				mpFes.resize(std::max(Nf, (uint32_t)mpFes.size() * 2));
				PF_npps.resize(mpFvs.size());
				mpF_boundary_flag.resize(mpFvs.size());
			}
			mpFvs[fid] = nvs2[0];
			mpFvs[fn] = nvs2[1];
			mpFes[fid] = nes2[0];
			mpFes[fn] = nes2[1];
			PF_npps[fn] = PF_npps[fid];
			mpF_boundary_flag[fn] = mpF_boundary_flag[fid];
		}
	}

	mpEs.resize(Ne);
	PE_npfs.resize(Ne);
	mpFvs.resize(Nf);
	mpFes.resize(Nf);
	PF_npps.resize(Nf);
	mpF_boundary_flag.resize(Nf);

	mV_tag = mO_copy;
	newQ = mQ_copy;
	
	if (Qquadric)
		newQu3D = Quadric_copy;

	return true;
}
bool MultiResolutionHierarchy::split_face3D(bool red_edge) {
	int32_t split_num = 0;
	uint32_t Ne = mpEs.size(), Nf = mpFvs.size();

	vector<short>().swap(mpV_flag); mpV_flag.resize(mO_copy.cols(), false);
	vector<short>().swap(mpE_flag); mpE_flag.resize(mpEs.size(), false);
	PV_npfs_sudo.clear(); PV_npfs_sudo.resize(mO_copy.cols());

	for (uint32_t i = 0; i < mpFvs.size(); i++) {
		
		std::vector<uint32_t> &es = mpFes[i], es_temp, &vs = mpFvs[i];
		//orient es direction
		for (uint32_t j = 0; j < vs.size(); j++) {
			for (auto e : es) {
				if ((get<0>(mpEs[e]) == vs[j] || get<0>(mpEs[e]) == vs[(j + 1) % vs.size()]) &&
					(get<1>(mpEs[e]) == vs[j] || get<1>(mpEs[e]) == vs[(j + 1) % vs.size()])) {
					es_temp.push_back(e); break;
				}
			}
		}
		es = es_temp;
		int32_t start = -1, end = -1;
		vector<tuple<double, uint32_t, uint32_t>> es_rank;
		tuple<short, Float, Vector3f> a_posy;
		for (int32_t j = 0; j < vs.size(); j++) {
			int32_t pre = (j - 1 + vs.size()) % vs.size();
			for (int32_t k = j + 2; k < vs.size(); k++) {
				if (pre == k) continue;//share an original edge
				
				uint32_t v0 = vs[j], v1 = vs[k];
				Quaternion q0 = mQ_copy.col(v0), q1 = Quaternion::applyRotation(mQ_copy.col(v1), q0);
				a_posy = posy3D_completeInfo(mO_copy.col(v0), q0, mO_copy.col(v1), q1, mScale, mInvScale);


				if (red_edge) {
					if (get<0>(a_posy) != Edge_tag::R) continue;
					start = j;
					end = k;
					break;
				}else if (!((k - j + vs.size() + 1) % vs.size() < 4 || (j - k + vs.size() + 1) % vs.size() < 4)) {
					Float cost = compute_cost_edge3D(vs[j], vs[k]);
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

		if (find(PV_npvs[vs[start]].begin(), PV_npvs[vs[start]].end(), vs[end]) != PV_npvs[vs[start]].end()) continue;


		vector<vector<uint32_t>> nes2(2), nvs2(2);
		uint32_t en = Ne++;
		uint32_t fn = Nf++;

		tuple_E new_e;
		get<0>(new_e) = vs[start];
		get<1>(new_e) = vs[end];
		get<2>(new_e) = mpF_boundary_flag[i];
		get<4>(new_e) = get<0>(a_posy);
		get<3>(new_e) = get<1>(a_posy);
		get<5>(new_e) = en;
		get<6>(new_e) = 0;
		get<7>(new_e) = 0;


		if (Ne > mpEs.size()) {
			mpEs.resize(std::max(Ne, (uint32_t)mpEs.size() * 2));
			PE_npfs.resize(mpEs.size());
			mpE_flag.resize(mpEs.size(), false);
		}
		mpEs[en] = new_e;

		for (uint32_t m = 0; m < 2; m++) {
			if (m == 1) std::swap(start, end);
			int32_t length = (end - start + vs.size() + 1) % vs.size();
			std::vector<uint32_t> nes(length), nvs(length);
			for (uint32_t k = 0; k < length; k++) {
				nvs[k] = vs[(start + k) % vs.size()];
				if (k + 1 == length) nes[k] = std::get<5>(new_e);
				else nes[k] = es[(start + k) % es.size()];
			}

			vector<vector<uint32_t>> fvs_simple_test, fes_simple_test;
			fvs_simple_test.push_back(nvs);
			fes_simple_test.push_back(nes);
			vector<uint32_t> pvs, pes, vs_disgard, es_disgard;
			if (!simple_polygon_3D_v2(fvs_simple_test, fes_simple_test, pvs, pes, vs_disgard, es_disgard, false))
				;

			nes2[m] = nes; nvs2[m] = nvs;
		}

		bool this_one = true;
		//simple polyhedral
		for (auto pid : PF_npps[i]) {
			std::vector<std::vector<uint32_t>> pfes;
			for (uint32_t n = 0; n < mpPs[pid].size(); n++) {
				if (mpPs[pid][n] != i) pfes.push_back(mpFes[mpPs[pid][n]]);
				else {
					pfes.push_back(nes2[0]);
					pfes.push_back(nes2[1]);
				}
			}
			if (!simple_polyhedral_v3(pfes)) { this_one = false; break; }
		}

		if (!this_one) {
			tuple_E e_ran;
			mpEs[en] = e_ran;
			Ne--;
			Nf--;
			continue;
		}
		PV_npvs[vs[start]].push_back(vs[end]);
		PV_npvs[vs[end]].push_back(vs[start]);

		PE_npfs[en].push_back(i);
		PE_npfs[en].push_back(fn);

		for (auto eid_ : nes2[1])
			if (en!=eid_ && find(PE_npfs[eid_].begin(), PE_npfs[eid_].end(), i) != PE_npfs[eid_].end())
				replace(PE_npfs[eid_].begin(), PE_npfs[eid_].end(), i, fn);

		for (auto pid : PF_npps[i]) mpPs[pid].push_back(fn);
		
		if (Nf > mpFvs.size()) {
			mpFvs.resize(std::max(Nf, (uint32_t)mpFvs.size() * 2));
			mpFes.resize(std::max(Nf, (uint32_t)mpFes.size() * 2));
			PF_npps.resize(mpFvs.size());
			mpF_boundary_flag.resize(mpFvs.size());
		}
		mpFvs[i] = nvs2[0];
		mpFvs[fn] = nvs2[1];
		mpFes[i] = nes2[0];
		mpFes[fn] = nes2[1];
		PF_npps[fn] = PF_npps[i];
		mpF_boundary_flag[fn] = mpF_boundary_flag[i];

		split_num++;
	}
	mpEs.resize(Ne);
	PE_npfs.resize(Ne);
	mpFvs.resize(Nf);
	mpFes.resize(Nf);
	PF_npps.resize(Nf);
	mpF_boundary_flag.resize(Nf);

	return split_num;
}
bool MultiResolutionHierarchy::split_polyhedral3D() {
	int32_t split_num = 0;

	PV_npes_sudo.clear(); PV_npes_sudo.resize(mO_copy.cols());
	PV_npfs_sudo.clear(); PV_npfs_sudo.resize(mO_copy.cols());
	PE_npfs_sudo.clear(); PE_npfs_sudo.resize(mpEs.size());

	vector<short>().swap(mpV_flag); mpV_flag.resize(mO_copy.size(), false);
	vector<short>().swap(mpE_flag); mpE_flag.resize(mpEs.size(), false);
	vector<short>().swap(mpF_flag); mpF_flag.resize(mpFvs.size(), false);

	uint32_t Nf = mpFvs.size(), Np = mpPs.size();
	for (uint32_t i = 0; i < mpPs.size(); i++) {
		vector<uint32_t> &fs = mpPs[i], es, vs;
		//red vertex pairs detection
		for (auto fid : fs) {
			es.insert(es.end(), mpFes[fid].begin(), mpFes[fid].end());
			for (auto vid : mpFvs[fid]) PV_npfs_sudo[vid].push_back(fid);
		}
		sort(es.begin(), es.end()); es.erase(unique(es.begin(), es.end()), es.end());
		for (auto eid : es) {
			uint32_t v0 = std::get<0>(mpEs[eid]), v1 = std::get<1>(mpEs[eid]);
			PV_npes_sudo[v0].push_back(eid);
			PV_npes_sudo[v1].push_back(eid);

			vs.push_back(v0); vs.push_back(v1);
		}
		sort(vs.begin(), vs.end()); vs.erase(unique(vs.begin(), vs.end()), vs.end());

		for (auto vid : vs) std::sort(PV_npes_sudo[vid].begin(), PV_npes_sudo[vid].end());
		for (auto vid : vs) std::sort(PV_npfs_sudo[vid].begin(), PV_npfs_sudo[vid].end());

		vector<pair<uint32_t, uint32_t>> red_pairs;
		for (uint32_t j = 0; j < vs.size(); j++)for (uint32_t k = j + 1; k < vs.size(); k++) {

			uint32_t v0 = vs[j], v1 = vs[k];
			Quaternion q0 = mQ_copy.col(v0), q1 = Quaternion::applyRotation(mQ_copy.col(v1), q0);
			tuple<short, Float, Vector3f> a_posy = posy3D_completeInfo(mO_copy.col(v0), q0, mO_copy.col(v1), q1, mScale, mInvScale);

			vector<uint32_t> sharedes, sharedfs;
			set_intersection(PV_npes_sudo[vs[j]].begin(), PV_npes_sudo[vs[j]].end(), PV_npes_sudo[vs[k]].begin(), PV_npes_sudo[vs[k]].end(), back_inserter(sharedes));
			if (sharedes.size()) continue;
			set_intersection(PV_npfs_sudo[vs[j]].begin(), PV_npfs_sudo[vs[j]].end(), PV_npfs_sudo[vs[k]].begin(), PV_npfs_sudo[vs[k]].end(), back_inserter(sharedfs));
			if (sharedfs.size()) {
				continue;
			}

			red_pairs.push_back(make_pair(vs[j], vs[k]));
		}

		for (auto vid : vs) PV_npes_sudo[vid].clear();
		for (auto vid : vs) PV_npfs_sudo[vid].clear();

		if (!red_pairs.size())continue;
		//candidates
		vector<pair<double, uint32_t>> fs_rank;
		vector<vector<uint32_t>> nfes, nfvs;
		candidate_loops(fs, nfes, nfvs, fs_rank);
		//find the best
		std::sort(fs_rank.begin(), fs_rank.end());
		for (auto f : fs_rank) {
			bool this_one = false;
			for (auto a_pair : red_pairs) {
				if (find(nfvs[f.second].begin(), nfvs[f.second].end(), a_pair.first) != nfvs[f.second].end() &&
					find(nfvs[f.second].begin(), nfvs[f.second].end(), a_pair.second) != nfvs[f.second].end())
					this_one = true;
			}
			if(!this_one) continue;
			//cut a polyhedral into two
			vector<vector<uint32_t>> nps(2);
			vector<vector<uint32_t>> pfes_temp;
			for (auto fid : fs) pfes_temp.push_back(mpFes[fid]);
			cut_a_polyhedral(fs, pfes_temp, nfes[f.second], nps[0], nps[1]);

			uint32_t fn = Nf++;
			uint32_t pn = Np++;

			for (uint32_t m = 0; m < 2; m++) nps[m].push_back(fn);

			if (Nf > mpFvs.size()) {
				mpFvs.resize(std::max(Nf, (uint32_t)mpFvs.size() * 2));
				mpFes.resize(std::max(Nf, (uint32_t)mpFes.size() * 2));
				PF_npps.resize(mpFvs.size());
				mpF_boundary_flag.resize(mpFvs.size());
				mpF_flag.resize(mpFvs.size(), false);
			}
			mpFvs[fn] = nfvs[f.second];
			mpFes[fn] = nfes[f.second];
			PF_npps[fn].push_back(i);
			PF_npps[fn].push_back(pn);
			mpF_boundary_flag[fn] = false;
			for (auto eid : nfes[f.second]) PE_npfs[eid].push_back(fn);

			split_num++;

			if (Np > mpPs.size()) {
				mpPs.resize(std::max(Nf, (uint32_t)mpPs.size() * 2));
			}
			
			mpPs[i] = nps[0];
			mpPs[pn] = nps[1];

			for (auto fid : mpPs[pn]) {
				if (fid != fn) replace(PF_npps[fid].begin(), PF_npps[fid].end(), i, pn);
			}

			break;
		}
	}
	mpFvs.resize(Nf);
	mpFes.resize(Nf);
	PF_npps.resize(Nf);
	mpF_boundary_flag.resize(Nf);
	mpF_flag.resize(Nf);
	mpPs.resize(Np);
	return split_num;
}
void MultiResolutionHierarchy::candidate_loops(vector<uint32_t> &fs, vector<vector<uint32_t>> &nfes, vector<vector<uint32_t>> &nfvs, vector<pair<double, uint32_t>> &fs_rank){
	//es, vs
	vector<uint32_t> es;
	for (auto fid : fs)es.insert(es.end(),mpFes[fid].begin(), mpFes[fid].end());
	sort(es.begin(), es.end()); es.erase(unique(es.begin(), es.end()), es.end());

	for (auto eid : es) {
		uint32_t v0 = get<0>(mpEs[eid]), v1 = get<1>(mpEs[eid]);
		PV_npes_sudo[v0].push_back(eid);
		PV_npes_sudo[v1].push_back(eid);
	}
	for (auto eid : es) sort(PE_npfs[eid].begin(), PE_npfs[eid].end());

	//edge pairs
	fs_rank.clear();
	for (uint32_t i = 0; i < es.size(); i++) {
		uint32_t v0[2], v1[2];
		v0[0] = get<0>(mpEs[es[i]]), v0[1] = get<1>(mpEs[es[i]]);
		for (uint32_t j = i + 1; j < es.size(); j++) {
			v1[0] = get<0>(mpEs[es[j]]); v1[1] = get<1>(mpEs[es[j]]);

			std::vector<uint32_t> fvs, fes, rvs;
			rvs.push_back(v0[0]); rvs.push_back(v0[1]);
			rvs.push_back(v1[0]); rvs.push_back(v1[1]);
			fvs.push_back(v0[0]); fes.push_back(es[i]);
			if (!recursive_ring(rvs, fvs, fes, v0[0], v1)) continue;
			std::reverse(fvs.begin(), fvs.end());
			std::reverse(fes.begin(), fes.end());
			fvs.push_back(v0[1]);
			if (!recursive_ring(rvs, fvs, fes, v0[1], v1)) continue;
			int32_t fid = share_the_same_fs(fes, es[j]);
			if (fid != -1) continue; //share boundary triangles/quads/pentagons

			fes.push_back(es[j]);

			vector<vector<uint32_t>> fvs_simple_test, fes_simple_test;
			fvs_simple_test.push_back(fvs);
			fes_simple_test.push_back(fes);
			vector<uint32_t> pvs, pes, vs_disgard, es_disgard;
			if (!simple_polygon_3D_v2(fvs_simple_test, fes_simple_test, pvs, pes, vs_disgard, es_disgard, false)) {
				continue;
			}
			Float cost = compute_cost_face3D(fvs, -1, false); ;
			fs_rank.push_back(std::make_pair(cost, nfes.size()));
			nfes.push_back(fes);
			nfvs.push_back(fvs);
		}
	}
	for (auto eid : es) {
		uint32_t v0 = get<0>(mpEs[eid]), v1 = get<1>(mpEs[eid]);
		PV_npes_sudo[v0].clear();
		PV_npes_sudo[v1].clear();
	}	
}
bool MultiResolutionHierarchy::recursive_ring(vector<uint32_t> &rvs, vector<uint32_t> &pvs, vector<uint32_t> &pes, uint32_t v0, uint32_t ve[2]) {
	if (v0 == ve[0] || v0 == ve[1]) return true;

	std::vector<std::tuple<Float, uint32_t, uint32_t>> potentials;
	for (auto eid : PV_npes_sudo[v0]) {

		if (find(pes.begin(), pes.end(), eid) != pes.end()) continue;

		int32_t fid = share_the_same_fs(pes, eid);
		if (fid != -1) continue; //share boundary triangles/quads/pentagons

		uint32_t vid = get<0>(mpEs[eid]);
		if (get<0>(mpEs[eid]) == v0) vid = get<1>(mpEs[eid]);

		Float cost_ = compute_cost_face3D(rvs, vid, true);
		potentials.push_back(make_tuple(cost_, eid, vid));
	}
	sort(potentials.begin(), potentials.end());
	if (!potentials.size()) return false;
	for (auto a_tuple : potentials) {
		pvs.push_back(get<2>(a_tuple));
		pes.push_back(get<1>(a_tuple));
		rvs.push_back(get<2>(a_tuple));

		if (!recursive_ring(rvs, pvs, pes, get<2>(a_tuple), ve)) {
			pvs.pop_back();
			pes.pop_back();
			rvs.pop_back();
			continue;
		}
		return true;
	}
	return false;
};
int32_t MultiResolutionHierarchy::share_the_same_fs(vector<uint32_t> &es, const uint32_t eid) {
	for (uint32_t i = 0; i < es.size(); i++) {
		std::vector<uint32_t> sharedf;
		std::set_intersection(PE_npfs[es[i]].begin(), PE_npfs[es[i]].end(), PE_npfs[eid].begin(), PE_npfs[eid].end(), back_inserter(sharedf));
		if (sharedf.size()) {
			if (mpFes[sharedf[0]].size() <= 4) return sharedf[0]; //share boundary triangles/quads
			uint32_t v0_e0 = get<0>(mpEs[eid]), v1_e0 = get<1>(mpEs[eid]),
				v0_ei = get<0>(mpEs[es[i]]), v1_ei = get<1>(mpEs[es[i]]);
			if (!(v0_e0 == v0_ei || v1_e0 == v0_ei || v0_e0 == v1_ei || v1_e0 == v1_ei)) return sharedf[0];//isolated two edges
																										   //not in the same direction
		}
	}
	return -1;
}
Float MultiResolutionHierarchy::compute_cost_face3D(vector<uint32_t> &vs, uint32_t rv, bool single) {
	vector<uint32_t> vs0; vector<Quaternion> qs(vs.size());
	Float min_cost = numeric_limits<Float>::max();
	vector<Vector3f> dirs;
	if (single) {
		vs0.push_back(rv);
		Quaternion q0 = mQ_copy.col(rv);
		for (uint32_t i = 0; i < vs.size(); i++) qs[i] = Quaternion::applyRotation(mQ_copy.col(vs[i]), q0);

		dirs.resize(vs.size());
		for (uint32_t i = 0; i < vs.size(); i++)
			dirs[i] = exact_3dir(mO_copy.col(rv), q0, mO_copy.col(vs[i]), qs[i], mScale, mInvScale);
	}
	else {
		vs0 = vs;
		Quaternion q0 = mQ_copy.col(vs[0]); qs[0] = q0;
		for (uint32_t i = 1; i < vs.size(); i++) qs[i] = Quaternion::applyRotation(mQ_copy.col(vs[i]), q0);

		vector<Vector3f> dirs;
		for (uint32_t i = 0; i<vs0.size(); i++)for (uint32_t j = i + 1; j<vs.size(); j++)
			dirs.push_back(exact_3dir(mO_copy.col(vs0[i]), qs[i], mO_copy.col(vs[j]), qs[j], mScale, mInvScale));
	}
	for (uint32_t i = 0; i < 3; i++) {
		Float cost_i = 0;
		for (auto dir : dirs) {
			dir[i] = std::abs(dir[i]);
			cost_i += dir[i];
		}
		if (cost_i < min_cost) min_cost = cost_i;
	}
	return min_cost;
}
Float MultiResolutionHierarchy::compute_cost_edge3D(uint32_t v0, uint32_t v1) {
	Float min_cost = std::numeric_limits<Float>::max();
	std::vector<Vector3f> dirs;
	Quaternion q0 = mQ_copy.col(v0), q1;
	q1 = Quaternion::applyRotation(mQ_copy.col(v1), q0);
	Vector3f dir = exact_3dir(mO_copy.col(v0), q0, mO_copy.col(v1), q1, mScale, mInvScale);

	for (uint32_t i = 0; i < 3; i++) {
		Float cost_i = 0;
		for (uint32_t j = 0; j < 3; j++) {
			if (i != j) cost_i += dir[j] * dir[j];
		}

		if (cost_i < min_cost) min_cost = cost_i;
	}
	return min_cost;
}
