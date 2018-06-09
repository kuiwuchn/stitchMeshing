#include "hierarchy.h"
#include "quat.h"
#include "timer.h"
#include "positions.h"

void MultiResolutionHierarchy::smoothOrientationsTri(uint32_t l, bool alignment, bool randomization, bool extrinsic) {
    const MatrixXf &N = mN[l];
    const SMatrix &L = mL[l];
    MatrixXf &Q = mQ[l];

    Timer<> timer;
    double error = 0;
    int nLinks = 0;
    MatrixXf Q_new(Q.rows(), Q.cols());
    tbb::spin_mutex mutex;

    tbb::parallel_for(
        tbb::blocked_range<uint32_t>(0u, (uint32_t) L.outerSize(), GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
            std::vector<std::pair<uint32_t, Float>> neighbors;
            double errorLocal = 0;
            int nLinksLocal = 0;
            for (uint32_t k = range.begin(); k != range.end(); ++k) {
                SMatrix::InnerIterator it(L, k);

                uint32_t i = it.row();
                Vector3f q_i = Vector3f::Zero();
                Vector3f n_i = N.col(i);

				if (nV_boundary_flag[l][i]) {
					Q_new.col(i) = Q.col(i);
					continue;
				}

                neighbors.clear();
                for (; it; ++it) {
                    uint32_t j = it.col();
                    if (i == j)
                        continue;
                    neighbors.push_back(std::make_pair(j, it.value()));
                }

                if (randomization && neighbors.size() > 0)
                    pcg32(mOrientationIterations, k)
                        .shuffle(neighbors.begin(), neighbors.end());

                for (auto n : neighbors) {
                    uint32_t j = n.first;
                    Float value = n.second;
                    Float dp;

                    Vector3f q_j = Q.col(j), n_j = N.col(j);
					if (Two_rosy_flag) {
						if (extrinsic) {
							q_j = compat_orientation_extrinsic_2((q_i == Vector3f::Zero()) ? Q.col(i) : q_i, n_i, q_j, n_j);
						}
						else {
							q_j = compat_orientation_intrinsic_2((q_i == Vector3f::Zero()) ? Q.col(i) : q_i, n_i, q_j, n_j);
						}
					}
					else {
						if (extrinsic) {
							q_j = applyRotationExtrinsic((q_i == Vector3f::Zero()) ? Q.col(i) : q_i, n_i, q_j, n_j);
							dp = Q.col(i).dot(applyRotation(Q.col(i), N.col(i), Q.col(j), N.col(j)));
						}
						else {
							q_j = applyRotation((q_i == Vector3f::Zero()) ? Q.col(i) : q_i, n_i, q_j, n_j);
							dp = Q.col(i).dot(applyRotation(Q.col(i), N.col(i), Q.col(j), N.col(j)));
						}
					}

                    errorLocal += std::abs(std::acos(std::min(dp, (Float) 1)));
                    ++nLinksLocal;

                    q_i += q_j * value;
                }

                if (q_i != Vector3f::Zero())
                    Q_new.col(i) = (q_i - n_i.dot(q_i) * n_i).normalized();
            }
            tbb::spin_mutex::scoped_lock guard(mutex);
            error += errorLocal;
            nLinks += nLinksLocal;
        }
    );

    mOrientationIterations++;
    Q = std::move(Q_new);
}

void MultiResolutionHierarchy::smoothOrientationsTet(uint32_t l, bool alignment, bool randomization) {
    const MatrixXf &N = mN[l];
    const SMatrix &L = mL[l];
    MatrixXf &Q = mQ[l];

    Timer<> timer;

    double error = 0;
    int nLinks = 0;
    MatrixXf Q_new(Q.rows(), Q.cols());
    tbb::spin_mutex mutex;

#if 1
    tbb::parallel_for(
        tbb::blocked_range<uint32_t>(0u, (uint32_t) L.outerSize(), GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
            std::vector<std::pair<uint32_t, Float>> neighbors;
            double errorLocal = 0;
            int nLinksLocal = 0;
            for (uint32_t k = range.begin(); k != range.end(); ++k) {
#endif 
				SMatrix::InnerIterator it(L, k);

                uint32_t i = it.row();
                Quaternion q_i = Quaternion::Zero();
                Vector3f n_i = N.col(i);

                neighbors.clear();
                for (; it; ++it) {
                    uint32_t j = it.col();
                    if (i == j)
                        continue;
                    neighbors.push_back(std::make_pair(j, it.value()));
                }

                if (randomization && neighbors.size() > 0)
                    pcg32(mOrientationIterations, k)
                        .shuffle(neighbors.begin(), neighbors.end());

                for (auto n : neighbors) {
                    uint32_t j = n.first;
                    Float value = n.second;

                    Quaternion q_j = Q.col(j);
                    q_j = Quaternion::applyRotation(
                        q_j, (q_i == Quaternion::Zero()) ? Q.col(i) : q_i);

                    Float dp = Q.col(i).dot(Quaternion::applyRotation(Q.col(j), Q.col(i)));
                    errorLocal += std::abs(std::acos(std::min(dp, (Float) 1)));
                    ++nLinksLocal;

                    q_i += q_j * value;

                    if (alignment && n_i != Vector3f::Zero()) {
                        Float magnitude = q_i.norm();
                        q_i = Quaternion(q_i / magnitude).align(n_i) * magnitude;
                    }
                }

                if (q_i != Quaternion::Zero())
                    Q_new.col(i) = q_i.normalized();
            }
            tbb::spin_mutex::scoped_lock guard(mutex);
            error += errorLocal;
            nLinks += nLinksLocal;
#if 1
        }
    );
#endif
    mOrientationIterations++;
    Q = std::move(Q_new);
}

void MultiResolutionHierarchy::prolongOrientations(int level) {
	const SMatrix &P = mP[level];
	const MatrixXf &N = mN[level];
	for (int k = 0; k < P.outerSize(); ++k) {
		SMatrix::InnerIterator it(P, k);
		for (; it; ++it) {
			if (!tetMesh()) {
				if (nV_boundary_flag[level][it.row()]) continue;
				Vector3f q_j = mQ[level + 1].col(it.col());
				Vector3f n_i = N.col(it.row());
				mQ[level].col(it.row()) = q_j - n_i * n_i.dot(q_j);
			}
			else {
				Quaternion q_j = mQ[level + 1].col(it.col());
				Vector3f n_i = N.col(it.row());
				if (n_i != Vector3f::Zero()) {
					Float magnitude = q_j.norm();
					q_j = Quaternion(q_j / magnitude).align(n_i) * magnitude;
				}
				mQ[level].col(it.row()) = q_j;
			}
		}
	}
}

void MultiResolutionHierarchy::detectOrientationSingularitiesTri() {
    uint32_t singularityCount = 0;
    auto const &V = mV[0], &Q = mQ[0], &N = mN[0], &NF = mNF;
    MatrixXf &S = mOrientationSingularities;
    tbb::spin_mutex mutex;

    Timer<> timer;
    timer.beginStage("Detecting orientation singularities");
    tbb::parallel_for(
        tbb::blocked_range<uint32_t>(0u, (uint32_t) mF.cols(), GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
            for (uint32_t f = range.begin(); f != range.end(); ++f) {
                Vector3f face_center = Vector3f::Zero();
                int index = 0;
                for (int i = 0; i < 3; ++i) {
                    int i0 = mF(i, f), i1 = mF((i + 1) % 3, f);
                    index += findRotation(Q.col(i0), N.col(i0),
                                          Q.col(i1), N.col(i1));
                    face_center += V.col(i0);
                }
                index = index % 4;

                if (index != 0) {
                    tbb::spin_mutex::scoped_lock guard(mutex);
                    face_center *= 1.f / 3.f;
                    Vector3f color = Vector3f::Constant(1.f);
                    if (index == 1)
                        color = Vector3f::UnitZ();
                    else if (index == 3)
                        color = Vector3f::UnitX();
                    if (singularityCount + 1 > S.cols())
                        S.conservativeResize(9, S.cols() * 2 + 1);
                    S.col(singularityCount++)
                        << face_center + NF.col(f) * mAverageEdgeLength / 3, NF.col(f), color;
                }
            }
        }
    );

    S.conservativeResize(9, singularityCount);
    timer.endStage("Found " + std::to_string(singularityCount) + " singular faces");
}

void MultiResolutionHierarchy::detectOrientationSingularitiesTet() {
    uint32_t singularityCount = 0;
    const uint8_t tet_faces[4][3] = { { 1, 0, 2 }, { 3, 2, 0 }, { 1, 2, 3 }, { 0, 1, 3 } };
    auto const &V = mV[0], &Q = mQ[0];
    MatrixXf &S = mOrientationSingularities;
    tbb::spin_mutex mutex;

    Timer<> timer;
    timer.beginStage("Detecting orientation singularities");

#if 0
    tbb::parallel_for(
        tbb::blocked_range<uint32_t>(0u, (uint32_t) mF.cols(), GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
            for (uint32_t f = range.begin(); f != range.end(); ++f) {
#endif

    for (uint32_t t = 0; t < mT.cols(); ++t) {
        for (auto face : tet_faces) {
            Vector3f face_center = Vector3f::Zero();

            AxisToAxisRot rot1, rot2, rot3;
            for (int i = 0; i < 3; ++i) {
                int i0 = mT(face[i], t), i1 = mT(face[(i + 1) % 3], t);
                rot1 *= AxisToAxisRot(Q.col(i0), Q.col(i1));
                face_center += V.col(i0);
            }

            for (int i = 0; i < 3; ++i) {
                int i0 = mT(face[(i+1)%3], t), i1 = mT(face[(i + 2) % 3], t);
                rot2 *= AxisToAxisRot(Q.col(i0), Q.col(i1));
            }

            for (int i = 0; i < 3; ++i) {
                int i0 = mT(face[(i+2)%3], t), i1 = mT(face[(i + 3) % 3], t);
                rot3 *= AxisToAxisRot(Q.col(i0), Q.col(i1));
            }


            if (!rot1.isIdentity()) {
                face_center *= 1.f / 3.f;

                Vector3f tet_center =
                    0.25f * (V.col(mT(0, t)) + V.col(mT(1, t)) +
                             V.col(mT(2, t)) + V.col(mT(3, t)));
                Quaternion q0 = Q.col(mT(face[0], t));
                Quaternion q1 = Q.col(mT(face[1], t));
                Quaternion q2 = Q.col(mT(face[2], t));

                Vector3f normal =
                    Vector3f(V.col(mT(face[1], t)) - V.col(mT(face[0], t))).cross(
                    Vector3f(V.col(mT(face[2], t)) - V.col(mT(face[0], t)))).normalized();

                Vector3f rot_axis1 = q0.toMatrix().col(rot1.preservedAxis()) * (rot1.index > 6 ? 1 : -1);
                Vector3f rot_axis2 = q1.toMatrix().col(rot2.preservedAxis()) * (rot2.index > 6 ? 1 : -1);
                Vector3f rot_axis3 = q2.toMatrix().col(rot3.preservedAxis()) * (rot3.index > 6 ? 1 : -1);

                Vector3f rot_axis = (rot_axis1 + rot_axis2 + rot_axis3).normalized();

                Vector3f color = rot_axis.dot(normal) > 0 ? Vector3f::UnitX() : Vector3f::UnitZ();

                if (color == Vector3f::UnitZ()) {
					;/* std::cout << endl;
                    std::cout 
                        << rot_axis1.dot(normal) << " vs "
                        << rot_axis2.dot(normal) << " vs "
                        << rot_axis3.dot(normal) << std::endl;*/
                }
				if (rot1.index < 4 || rot1.index > 9)
					;

                tbb::spin_mutex::scoped_lock guard(mutex);
                if (singularityCount + 1 > S.cols())
                    S.conservativeResize(6, S.cols() * 2 + 2);
                S.col(singularityCount++) << face_center, color;
                S.col(singularityCount++) << tet_center, color;
            }
        }
    }

    S.conservativeResize(6, singularityCount);
    timer.endStage("Found " + std::to_string(singularityCount) + " singular faces");
}
