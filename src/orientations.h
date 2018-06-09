#pragma once

#include "quat.h"

inline int findRotation(const Vector3f &q0, const Vector3f &n0, const Vector3f &_q1, const Vector3f &n1) {
    const Vector3f q1 = rotateVectorIntoPlane(_q1, n1, n0);
    const Vector3f t1 = n0.cross(q1);
    const Float dp0 = q1.dot(q0), dp1 = t1.dot(q0);

    if (std::abs(dp0) > std::abs(dp1))
        return dp0 > 0 ? 0 : 2;
    else
        return dp1 > 0 ? 1 : 3;
}

inline Vector3f applyRotation(const Vector3f &q0, const Vector3f &n0,
                              const Vector3f &_q1, const Vector3f &n1) {
    const Vector3f q1 = rotateVectorIntoPlane(_q1, n1, n0);
    const Vector3f t1 = n0.cross(q1);
    const Float dp0 = q1.dot(q0), dp1 = t1.dot(q0);

    if (std::abs(dp0) > std::abs(dp1))
        return q1 * std::copysign(1.f, dp0);
    else
        return t1 * std::copysign(1.f, dp1);
}

inline Vector3f applyRotationKeep(const Vector3f &q0, const Vector3f &n0,
                                  const Vector3f &_q1, const Vector3f &n1) {
    const Vector3f q1 = rotateVectorIntoPlane(_q1, n1, n0);
    const Vector3f t1 = n0.cross(q1);
    const Float dp0 = q1.dot(q0), dp1 = t1.dot(q0);

    if (std::abs(dp0) > std::abs(dp1))
        return _q1 * std::copysign(1.f, dp0);
    else
        return n1.cross(_q1) * std::copysign(1.f, dp1);
}

inline Vector3f applyRotationExtrinsic(const Vector3f &q0, const Vector3f &n0,
                                       const Vector3f &q1, const Vector3f &n1) {
    const Vector3f list[4] = { q1, n1.cross(q1), n0.cross(q1),
                               n0.cross(n1.cross(q1)) };

    Float best_score = -std::numeric_limits<Float>::infinity();
    int best = -1;

    for (int i = 0; i < 4; ++i) {
        Float score = std::abs(q0.dot(list[i]));
        if (score > best_score) {
            best = i;
            best_score = score;
        }
    }

    const Float dp = q0.dot(list[best]);
    return list[best] * std::copysign(1.f, dp);
}

// a succinct class for AxisToAxis rotations
struct AxisToAxisRot {
    typedef unsigned char Index;

    Index index; /* the only member: an index, 0 to 23.  */

    inline bool isIdentity() const { return index == 0; }

    inline bool preservesOneAxis() const { return index <= 9; }
    inline int preservedAxis() const { return (30 - index) % 3; }

    inline void operator*=(const AxisToAxisRot &b) {
        index = mult(index, b.index);
    }
    inline void invert() { index = inverse(index); }
    AxisToAxisRot inverse() { return AxisToAxisRot(inverse(index)); }

    inline Quaternion toQuat() { return toQuat(index); }

    AxisToAxisRot() : index(0) {}

    AxisToAxisRot(int i) : index(i) {}

    AxisToAxisRot(const Quaternion &q, const Quaternion &ref) {
        Float dp[4] = {
            ref.dot(Quaternion(q.x(), q.y(), q.z(), q.w())),
            ref.dot(Quaternion(q.y(), -q.x(), q.w(), -q.z())),
            ref.dot(Quaternion(-q.z(), q.w(), q.x(), -q.y())),
            ref.dot(Quaternion(q.w(), q.z(), -q.y(), -q.x())),
        };
        Float a[4] = { std::abs(dp[0]), std::abs(dp[1]), std::abs(dp[2]),
                       std::abs(dp[3]) };
        // find M, N, max indices of a[i]
        int M = 0, N = 1, m = 2, n = 3;
        if (a[M] < a[m])
            std::swap(M, m);
        if (a[N] < a[n])
            std::swap(N, n);
        if (a[M] < a[N]) {
            std::swap(M, N);
            m = n;
        }
        if (a[N] < a[m])
            N = m;

        const Float s = std::sqrt(.5f);
        const Float h = .5f;

        Float vA = a[M];
        Float vB = (a[M] + a[N]) * s;
        Float vC = (a[0] + a[1] + a[2] + a[3]) * h;

        #define signDiff(A, B) (std::signbit(dp[A]) != std::signbit(dp[B]))

        if (vA > vB && vA > vC) {
            index = M;
        } else if (vB > vC) {
            index = 3 + M + N + (M * N != 0) * 4 + signDiff(M, N) * 3;
        } else {
            index = (4 + 6 + 6) + signDiff(1, 0) + signDiff(2, 0) * 2 +
                    signDiff(3, 0) * 4;
        }

        #undef signDiff
    }

    static float oneOverQuatNorm(Index i) {
        return (i < 4) ? 1 : (i < 16) ? sqrt(0.5) : 0.5;
    }

    static const inline Quaternion toQuat(Index i) {
        return toIntegerQuat(i) * oneOverQuatNorm(i);
    }

    static const Quaternion toIntegerQuat(Index i) {
        #define M -1
        static const char table[24][4] = {
            /* 0*/ { 0, 0, 0, 1 }, /* class A */

            /* 1*/ { 0, 0, 1, 0 },
            /* 2*/ { 0, 1, 0, 0 },
            /* 3*/ { 1, 0, 0, 0 },

            /* 4*/ { 0, 0, 1, 1 }, /* class B */
            /* 5*/ { 0, 1, 0, 1 },
            /* 6*/ { 1, 0, 0, 1 },
            /* 7*/ { 0, 0, M, 1 },
            /* 8*/ { 0, M, 0, 1 },
            /* 9*/ { M, 0, 0, 1 },

            /*10*/ { 0, 1, 1, 0 }, /* first "monster" */
            /*11*/ { 1, 0, 1, 0 },
            /*12*/ { 1, 1, 0, 0 },
            /*13*/ { 0, M, 1, 0 },
            /*14*/ { M, 0, 1, 0 },
            /*15*/ { M, 1, 0, 0 },

            /*16*/ { 1, 1, 1, 1 }, /* class C */
            /*17*/ { 1, 1, M, 1 },
            /*18*/ { 1, M, 1, 1 },
            /*19*/ { 1, M, M, 1 },
            /*20*/ { M, 1, 1, 1 },
            /*21*/ { M, 1, M, 1 },
            /*22*/ { M, M, 1, 1 },
            /*23*/ { M, M, M, 1 },
        };
        #undef M
        return Quaternion(table[i][0], table[i][1], table[i][2], table[i][3]);
    }

    static Index inverse(Index a) {
        static const Index table[24] = {
            0,  1,  2,  3,                  // class A
            7,  8,  9,  4,  5,  6,          // class B1
            10, 11, 12, 13, 14, 15,         // class B2
            23, 22, 21, 20, 19, 18, 17, 16, // class C
        };
        // if (mult(a,table[a])!=0) cout<<"NOOOO! A"<<int(a)<<endl;
        // if (mult(table[a],a)!=0) cout<<"NOOOO! B"<<int(a)<<endl;
        return table[a];
    }

    static Index mult(Index a, Index b) {
        static const Index table[24][24] = {
            { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 },
            { 1, 0, 3, 2, 7, 14, 10, 4, 11, 13, 6, 8, 15, 9, 5, 12, 19, 20, 23, 16, 17, 22, 21, 18 },
            { 2, 3, 0, 1, 12, 8, 13, 15, 5, 10, 9, 14, 4, 6, 11, 7, 22, 18, 17, 21, 23, 19, 16, 20 },
            { 3, 2, 1, 0, 15, 11, 9, 12, 14, 6, 13, 5, 7, 10, 8, 4, 21, 23, 20, 22, 18, 16, 19, 17 },
            { 4, 7, 15, 12, 1, 20, 16, 0, 18, 22, 19, 23, 2, 21, 17, 3, 10, 5, 11, 6, 14, 9, 13, 8 },
            { 5, 11, 8, 14, 16, 2, 17, 21, 0, 20, 23, 3, 22, 18, 1, 19, 12, 13, 6, 7, 10, 15, 4, 9 },
            { 6, 13, 10, 9, 18, 16, 3, 17, 19, 0, 1, 21, 23, 2, 22, 20, 11, 12, 15, 14, 4, 5, 8, 7 },
            { 7, 4, 12, 15, 0, 17, 19, 1, 23, 21, 16, 18, 3, 22, 20, 2, 6, 14, 8, 10, 5, 13, 9, 11 },
            { 8, 14, 5, 11, 22, 0, 18, 19, 2, 23, 20, 1, 16, 17, 3, 21, 4, 6, 13, 15, 9, 7, 12, 10 },
            { 9, 10, 13, 6, 20, 21, 0, 23, 22, 3, 2, 16, 17, 1, 19, 18, 5, 7, 4, 8, 15, 11, 14, 12 },
            { 10, 9, 6, 13, 23, 19, 2, 20, 16, 1, 0, 22, 18, 3, 21, 17, 8, 15, 12, 5, 7, 14, 11, 4 },
            { 11, 5, 14, 8, 21, 1, 23, 16, 3, 18, 17, 0, 19, 20, 2, 22, 7, 10, 9, 12, 13, 4, 15, 6 },
            { 12, 15, 7, 4, 3, 23, 22, 2, 17, 16, 21, 20, 0, 19, 18, 1, 9, 8, 14, 13, 11, 10, 6, 5 },
            { 13, 6, 9, 10, 17, 22, 1, 18, 21, 2, 3, 19, 20, 0, 16, 23, 14, 4, 7, 11, 12, 8, 5, 15 },
            { 14, 8, 11, 5, 19, 3, 20, 22, 1, 17, 18, 2, 21, 23, 0, 16, 15, 9, 10, 4, 6, 12, 7, 13 },
            { 15, 12, 4, 7, 2, 18, 21, 3, 20, 19, 22, 17, 1, 16, 23, 0, 13, 11, 5, 9, 8, 6, 10, 14 },
            { 16, 21, 19, 22, 11, 10, 12, 5, 6, 4, 7, 9, 8, 15, 13, 14, 23, 2, 3, 17, 1, 20, 18, 0 },
            { 17, 18, 23, 20, 6, 12, 14, 13, 7, 5, 11, 15, 9, 8, 4, 10, 3, 22, 19, 1, 16, 2, 0, 21 },
            { 18, 17, 20, 23, 13, 4, 11, 6, 15, 8, 14, 7, 10, 5, 12, 9, 1, 16, 21, 3, 22, 0, 2, 19 },
            { 19, 22, 16, 21, 8, 6, 15, 14, 10, 7, 4, 13, 11, 12, 9, 5, 18, 3, 2, 20, 0, 17, 23, 1 },
            { 20, 23, 18, 17, 10, 15, 5, 9, 4, 14, 8, 12, 13, 11, 7, 6, 2, 21, 16, 0, 19, 3, 1, 22 },
            { 21, 16, 22, 19, 5, 13, 7, 11, 9, 15, 12, 6, 14, 4, 10, 8, 17, 1, 0, 23, 2, 18, 20, 3 },
            { 22, 19, 21, 16, 14, 9, 4, 8, 13, 12, 15, 10, 5, 7, 6, 11, 20, 0, 1, 18, 3, 23, 17, 2 },
            { 23, 20, 17, 18, 9, 7, 8, 10, 12, 11, 5, 4, 6, 14, 15, 13, 0, 19, 22, 2, 21, 1, 3, 16 },
        };
        return table[a][b];
    }

    /* only needed to construct mult table */
    static bool matchesUpToSign(Quaternion q0, Quaternion q1) {
        return (q0 - q1).squaredNorm() < 0.01 || (q0 + q1).squaredNorm() < 0.01;
    }

    /* only needed to construct mult table */
    static Index fromQuat(const Quaternion &q) {
        for (int i = 0; i < 24; i++) {
            if (matchesUpToSign(q, toQuat(i))) {
                return i;
            }
        }
        cout << "ERROR?!?" << endl;
        return 0;
    }

    /* only needed to construct mult table */
    static void printMultiplicationTable() {

        /*
        static char str[24][3] = {
            "A0","A1","A2","A3",
            "B0","B1","B2","B3","B4","B5","B6","B7","B8","B9","Ba","Bb",
            "C0","C1","C2","C3","C4","C5","C6","C7",
        };*/
        static char str[24][3] = {
            " 0", " 1", " 2", " 3", " 4", " 5", " 6", " 7",
            " 8", " 9", "10", "11", "12", "13", "14", "15",
            "16", "17", "18", "19", "20", "21", "22", "23",
        };
        cout << endl
             << endl
             << "{" << endl;
        for (int i = 0; i < 24; i++) {
            cout << "  {";
            for (int j = 0; j < 24; j++)
                cout << str[fromQuat(toQuat(i) * toQuat(j))] << ",";
            cout << "}," << endl;
        }
        cout << "}" << endl
             << endl;
    }
};
//2-Rosy
inline Vector3f compat_orientation_extrinsic_2(const Vector3f &q0, const Vector3f &n0,
	const Vector3f &q1, const Vector3f &n1) {
	return q1 * signum(q0.dot(q1));
}
inline Vector3f compat_orientation_intrinsic_2(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &_q1, const Vector3f &n1) {
	const Vector3f q1 = rotateVectorIntoPlane(_q1, n1, n0);
	return q1 * signum(q1.dot(q0));
}