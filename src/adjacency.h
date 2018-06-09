#pragma once

#include "common.h"
#include <map>

#define INVALID ((uint32_t) -1)

/**
 * \param C
 *    #F by 3/4 list of triangle/quad indices
 * \param[out] CC
 *    #T by 3/4 list indicating adjacent triangles/quads
 * \param[out] CCi
 *    #T by 3/4 list indicating the corresponding edge index in the adjacent triangle/quad
 */
inline void compute_tri_adjacency_2d(const MatrixXu &F, MatrixXu &FF, MatrixXu &FFi) {
    FF.setConstant(F.rows(), F.cols(), INVALID);
    FFi.setConstant(F.rows(), F.cols(), INVALID);

    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> temp;
    temp.reserve(F.size());
    for (uint32_t f = 0; f < F.cols(); ++f) {
        for (uint32_t e = 0; e < F.rows(); ++e) {
            uint32_t v0 = F(e, f), v1 = F((e + 1) % F.rows(), f);
            if (v0 > v1)
                std::swap(v0, v1);
            temp.push_back(std::make_tuple(v0, v1, f, e));
        }
    }

    std::sort(temp.begin(), temp.end());

    for (uint32_t i = 1; i < temp.size(); ++i) {
        auto e0 = temp[i-1], e1 = temp[i];

        if (std::get<0>(e0) == std::get<0>(e1) &&
            std::get<1>(e0) == std::get<1>(e1)) {
            FF (std::get<3>(e0), std::get<2>(e0)) = std::get<2>(e1);
            FFi(std::get<3>(e0), std::get<2>(e0)) = std::get<3>(e1);
            FF (std::get<3>(e1), std::get<2>(e1)) = std::get<2>(e0);
            FFi(std::get<3>(e1), std::get<2>(e1)) = std::get<3>(e0);
            ++i;
        }
    }
}

/**
 * Compute adjacency information in a tetrahedral mesh
 *
 * \param T
 *    #T by 4 list of tetrahedron indices
 * \param[out] TT
 *    #T by 4 list indicating adjacent tetrahedra
 * \param[out] TTi
 *    #T by 4 list indicating the corresponding face index in the adjacent tetrahedron
 * \param[out] TTe
 *    #T by 4 list indicating the edge ID corresponding to the first edge of each face in the adjacent tetrahedron
 */
inline void compute_tet_adjacency_3d(const MatrixXu &T, MatrixXu &TT, MatrixXu &TTi, MatrixXu &TTe) {
    TT.setConstant(T.rows(), T.cols(), INVALID);
    TTi.setConstant(T.rows(), T.cols(), INVALID);

    const uint8_t tet_faces[4][3] = {
        { 1, 0, 2 }, { 3, 2, 0 }, { 1, 2, 3 }, { 0, 1, 3 } };

    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>> temp;
    temp.reserve(T.size());
    for (uint32_t t = 0; t < T.cols(); ++t) {
        for (uint32_t f = 0; f < 4; ++f) {
            uint32_t v0 = T(tet_faces[f][0], t),
                     v1 = T(tet_faces[f][1], t),
                     v2 = T(tet_faces[f][2], t);
            if (v1 > v2) std::swap(v1, v2);
            if (v0 > v2) std::swap(v0, v2);
            if (v0 > v1) std::swap(v0, v1);
            temp.push_back(std::make_tuple(v0, v1, v2, t, f));
        }
    }

    std::sort(temp.begin(), temp.end());

    for (uint32_t i = 1; i < temp.size(); ++i) {
        auto e0 = temp[i-1], e1 = temp[i];
        if (std::get<0>(e0) == std::get<0>(e1) &&
            std::get<1>(e0) == std::get<1>(e1) &&
            std::get<2>(e0) == std::get<2>(e1)) {
            uint32_t t0 = std::get<3>(e0), f0 = std::get<4>(e0);
            uint32_t t1 = std::get<3>(e1), f1 = std::get<4>(e1);

            TT(f0, t0) = t1; TTi(f0, t0) = f1;
            TT(f1, t1) = t0; TTi(f1, t1) = f0;
            ++i;
        }
    }

    TTe.setConstant(4, T.cols(), INVALID);
    for (uint32_t t = 0; t < T.cols(); ++t) {
        for (uint32_t f = 0; f < 4; ++f) {
            uint32_t t2 = TT(f, t);
            uint32_t f2 = TTi(f, t);

            if (t2 == INVALID)
                continue;

            uint32_t f0v0 = T(tet_faces[f][0], t);
            uint32_t f0v1 = T(tet_faces[f][1], t);

            for (uint32_t i = 0; i < 3; ++i) {
                uint32_t f1v0 = T(tet_faces[f2][i], t2);
                uint32_t f1v1 = T(tet_faces[f2][(i+1)%3], t2);

                if (std::make_pair(f0v0, f0v1) ==
                    std::make_pair(f1v1, f1v0)) {
                    TTe(f, t) = i;
                    break;
                }
            }
        }
    }
}

inline void compute_vertex_to_face_2d(const MatrixXf &V, const MatrixXu &F, MatrixXu &VFi) {
    VFi.setConstant(2, V.cols(), INVALID);
    for (uint32_t f = 0; f < F.cols(); ++f)
        for (uint32_t i = 0; i < F.rows(); ++i)
            VFi.col(F(i, f)) << f, i;
}

inline void compute_edge_to_tet_3d(const MatrixXu &T, MatrixXu &E, MatrixXu &ETi) {
    const uint8_t tet_faces[4][3] = {
        { 1, 0, 2 }, { 3, 2, 0 }, { 1, 2, 3 }, { 0, 1, 3 } };

    typedef std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t> Tuple;
    std::vector<Tuple> temp;
    temp.reserve(T.cols() * 6);

    for (uint32_t t = 0; t < T.cols(); ++t) {
        for (uint32_t f = 0; f < 4; ++f) {
            for (uint32_t i = 0; i < 3; ++i) {
                uint32_t v0 = T(tet_faces[f][i], t);
                uint32_t v1 = T(tet_faces[f][(i+1)%3], t);
                if (v0 < v1)
                    temp.push_back(std::make_tuple(v0, v1, t, f, i));
            }
        }
    }

    std::sort(temp.begin(), temp.end());

    temp.erase(std::unique(temp.begin(), temp.end(),
                           [](const Tuple &a, const Tuple &b) {
                               return std::get<0>(a) == std::get<0>(b) &&
                                      std::get<1>(a) == std::get<1>(b);
                           }),
               temp.end());

    E.resize(2, temp.size());
    ETi.resize(3, temp.size());

    for (uint32_t i=0; i<temp.size(); ++i) {
        E.col(i) << std::get<0>(temp[i]), std::get<1>(temp[i]);
        ETi.col(i) << std::get<2>(temp[i]), std::get<3>(temp[i]), std::get<4>(temp[i]);
    }
}

class circulator_3d {
public:
    struct iterator {
        iterator(circulator_3d *circ, uint32_t t, uint32_t f, uint32_t e)
            : m_circ(circ), m_t(t), m_f(f), m_e(e), m_first(t), m_success(false) { }

        iterator(const iterator &it) = default;
        iterator &operator=(const iterator &) = default;

        iterator &operator++() {
            jump_to_opposite_tet();

            if (m_t == m_first || m_t == INVALID) {
                m_success = m_t == m_first;
                m_t = m_f = m_e = INVALID;
                return *this;
            }

            jump_to_opposite_face();

            return *this;
        }

        uint32_t operator*() const { return m_t; }
        uint32_t t() const { return m_t; }
        uint32_t f() const { return m_f; }
        uint32_t e() const { return m_e; }

        bool operator==(const iterator &it) const { return m_t == it.m_t && m_f == it.m_f && m_e == it.m_e; }
        bool operator!=(const iterator &it) const { return m_t != it.m_t || m_f != it.m_f || m_e != it.m_e; }

        bool success() { return m_success; }

        void jump_to_opposite_tet() {
            std::tie(m_t, m_f, m_e) =
                std::make_tuple(m_circ->m_TT(m_f, m_t), m_circ->m_TTi(m_f, m_t),
                                (m_circ->m_TTe(m_f, m_t) - m_e + 3) % 3);
        }

        void jump_to_opposite_face() {
            const uint8_t jump_tbl[4][3][2] = {{{3, 0}, {1, 1}, {2, 0}},
                                               {{2, 1}, {0, 1}, {3, 2}},
                                               {{0, 2}, {1, 0}, {3, 1}},
                                               {{0, 0}, {2, 2}, {1, 2}}};
            std::tie(m_f, m_e) =
                std::make_pair(jump_tbl[m_f][m_e][0], jump_tbl[m_f][m_e][1]);
        }

        circulator_3d *m_circ;
        uint32_t m_t, m_f, m_e;
        uint32_t m_first;
        bool m_success;
    };

public:
    circulator_3d(const MatrixXu &TT, const MatrixXu &TTi, const MatrixXu &TTe,
                  const MatrixXu &ETi, uint32_t e, bool ccw = true)
        : m_TT(TT), m_TTi(TTi), m_TTe(TTe),
          m_start(this, ETi(0, e), ETi(1, e), ETi(2, e)) {
        if (!ccw)
            m_start.jump_to_opposite_face();
    }

    iterator begin() { return m_start; }
    iterator end() { return iterator(this, INVALID, INVALID, INVALID); }
private:
    const MatrixXu &m_TT, &m_TTi, &m_TTe;
    iterator m_start;
};

class circulator_2d {
public:
    struct iterator {
        iterator() { }
        iterator(circulator_2d *circ, uint32_t f, uint32_t e, uint32_t v)
            : m_circ(circ), m_f(f), m_e(e), m_v(v), m_first(f), m_success(false) { }

        iterator(const iterator &it) = default;
        iterator &operator=(const iterator &) = default;

        iterator &operator++() {
            jump_to_opposite_face();
            if (m_f == INVALID || m_f == m_first) {
                m_success = m_f == m_first;
                m_f = m_e = m_v = INVALID;
                return *this;
            }
            jump_to_opposite_edge();
            return *this;
        }

        void jump_to_opposite_face() {
            std::tie(m_f, m_e) = std::make_pair(m_circ->m_FF(m_e, m_f), m_circ->m_FFi(m_e, m_f));
            m_v = m_circ->m_ccw ? m_e : (m_e + 1) % 3;
        }

        void jump_to_opposite_edge() {
            m_e = (5 - m_v - m_e) % 3;
        }

        uint32_t operator*() { return m_f; }
        uint32_t f() const { return m_f; }
        uint32_t e() const { return m_e; }
        uint32_t v() const { return m_v; }
        bool success() const { return m_success; }

        bool operator==(const iterator &it) const { return m_f == it.m_f && m_e == it.m_e && m_v == it.m_v; }
        bool operator!=(const iterator &it) const { return m_f != it.m_f || m_e != it.m_e || m_v != it.m_v; }

        operator bool() const { return m_f != INVALID && m_e != INVALID && m_v != INVALID; };

        circulator_2d *m_circ = nullptr;
        uint32_t m_f = INVALID;
        uint32_t m_e = INVALID;
        uint32_t m_v = INVALID;
        uint32_t m_first;
        bool m_success;
    };

public:
    circulator_2d(const MatrixXu &FF, const MatrixXu &FFi, const MatrixXu &VFi,
                  uint32_t v, bool ccw = true)
        : m_FF(FF), m_FFi(FFi),
          m_start(this, VFi(0, v), (VFi(1, v) + 2) % 3, VFi(1, v)), m_ccw(ccw) {
        if (!ccw)
            m_start.jump_to_opposite_edge();
    }

    iterator begin() { return m_start; }
    iterator end() { return iterator(this, INVALID, INVALID, INVALID); }

    const MatrixXu &m_FF, &m_FFi;
    iterator m_start;
    bool m_ccw;
};

inline void rewind_vertex_to_face_2d(const MatrixXu &FF, const MatrixXu &FFi,
                                  MatrixXu &VFi) {
    for (uint32_t v = 0; v < VFi.cols(); ++v) {
        circulator_2d circ(FF, FFi, VFi, v, false);

        auto it = circ.begin(), saved_it = it;
        for (; it != circ.end(); ++it)
            saved_it = it;

        if (!it.success() && saved_it != circ.begin())
            VFi.col(v) << saved_it.f(), saved_it.v();
    }
}

inline void rewind_edge_to_tet_3d(const MatrixXu &TT, const MatrixXu &TTi,
                                  const MatrixXu &TTe, const MatrixXu &E,
                                  MatrixXu &ETi) {
    for (uint32_t e = 0; e < E.cols(); ++e) {
        circulator_3d circ(TT, TTi, TTe, ETi, e, false);
        auto it = circ.begin(), saved_it = it;
        for (; it != circ.end(); ++it)
            saved_it = it;
        if (!it.success() && saved_it != circ.begin()) {
            saved_it.jump_to_opposite_face();
            ETi.col(e) << saved_it.t(), saved_it.f(), saved_it.e();
        }
    }
}
