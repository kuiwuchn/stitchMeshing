#pragma once

#include "orientations.h"

inline Vector3f middle_point(const Vector3f &p0, const Vector3f &n0, const Vector3f &p1, const Vector3f &n1) {
    Float n0p0 = n0.dot(p0), n0p1 = n0.dot(p1),
          n1p0 = n1.dot(p0), n1p1 = n1.dot(p1),
          n0n1 = n0.dot(n1),
          denom = 1.0f / (1.0f - n0n1*n0n1 + 1e-4f),
          lambda_0 = 2.0f*(n0p1 - n0p0 - n0n1*(n1p0 - n1p1))*denom,
          lambda_1 = 2.0f*(n1p0 - n1p1 - n0n1*(n0p1 - n0p0))*denom;
    return 0.5f * (p0 + p1) - 0.25f * (n0 * lambda_0 + n1 * lambda_1);
}

inline std::pair<Vector3f, Vector3i>
findClosestPair(const Vector3f &o0, const Quaternion &q0,
                const Vector3f &o1, const Quaternion &q1,
                Float scale, Float invScale) {
    Quaternion qn = (q0 + Quaternion::applyRotation(q1, q0)).normalized();
    Matrix3f M = qn.toMatrix();
    Vector3f rel = round(Vector3f(M.transpose() * ((o0 - o1) * invScale)));
    return std::make_pair(o1 + (M * rel) * scale, rel.cast<int>());
}
inline  Vector3f exact_3dir(const Vector3f &o0, const Quaternion &q0,
	const Vector3f &o1, const Quaternion &q1,
	Float scale, Float invScale) {
	Quaternion qn = (q0 + Quaternion::applyRotation(q1, q0)).normalized();
	Matrix3f M = qn.toMatrix();
	return M.transpose() * ((o0 - o1) * invScale);
	
}
inline std::tuple<short, Float, Vector3f> posy3D_completeInfo(const Vector3f &o0, const Quaternion &q0,
	const Vector3f &o1, const Quaternion &q1,
	Float scale, Float invScale) {
	Quaternion qn = (q0 + Quaternion::applyRotation(q1, q0)).normalized();
	Matrix3f M = qn.toMatrix();
	Vector3f vvv = (M.transpose() * ((o0 - o1) * invScale));

	// we want to count how many coords of vvv will be rounded to 1
	for (uint32_t i = 0; i < 3; i++) {
		vvv[i] = std::abs(vvv[i]);
	}
	short res = 0;  Float Weight = 0; std::vector<bool> bits(3);
	for (uint32_t i = 0; i < 3; i++) {
		if (vvv[i] > 0.5) { res++; bits[i] = true; }
		else bits[i] = false;
	}
	if (res == 0) { // it's a red! 
		for (int i = 0; i < 3; i++) Weight += vvv[i] * vvv[i];
	}
	else if (res == 1) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 1.0;
	}
	else if (res == 2) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 2.0;
	}
	else if (res == 3) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 3.0;
	}

	return std::make_tuple(res, Weight, vvv);
}

inline std::pair<short, short> assignColorBiased(const Vector3f &o0, const Quaternion &q0,
	const Vector3f &o1, const Quaternion &q1,
	Float scale, Float invScale) {
	Quaternion qn = (q0 + Quaternion::applyRotation(q1, q0)).normalized();
	Matrix3f M = qn.toMatrix();
	Vector3f vvv = (M.transpose() * ((o0 - o1) * invScale));

	// we want to count how many coords of vvv will be rounded to 1
	for (uint32_t i = 0; i < 3; i++) {
		vvv[i] = std::abs(vvv[i]);
	}
	short res = 0, ind; std::vector<bool> bits(3);
	for (uint32_t i = 0; i < 3; i++) {
		if (vvv[i] > 0.5) { res++; bits[i] = true; }
		else bits[i] = false;
	}
	if (res == 0) { // it's a red! 
		ind=-1;// for (int i = 0; i < 3; i++) if (vvv[i] > 0.2) res = 1;
	}
	else if (res == 1) {
		for (short i = 0; i < 3; i++) if (bits[i]) ind = i;
	}
	else if (res == 2) {
		for (short i = 0; i < 3; i++) if (!bits[i]) ind = i;
	}else if(res==3)
		ind = -1;
	return std::make_pair(res,ind);
}
inline std::pair<short, Float> assignColorWeighted3D(const Vector3f &o0, const Quaternion &q0,
	const Vector3f &o1, const Quaternion &q1,
	Float scale, Float invScale) {
	Quaternion qn = (q0 + Quaternion::applyRotation(q1, q0)).normalized();
	Matrix3f M = qn.toMatrix();
	Vector3f vvv = (M.transpose() * ((o0 - o1) * invScale));

	// we want to count how many coords of vvv will be rounded to 1
	for (uint32_t i = 0; i < 3; i++) {
		vvv[i] = std::abs(vvv[i]);
	}
	short res = 0;  Float Weight = 0; std::vector<bool> bits(3);
	for (uint32_t i = 0; i < 3; i++) {
		if (vvv[i] > 0.5) { res++; bits[i] = true; }
		else bits[i] = false;
	}
	if (res == 0) { // it's a red! 
		for (int i = 0; i < 3; i++) Weight += vvv[i] * vvv[i];
	}
	else if (res == 1) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 1.0;
	}
	else if (res == 2) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 2.0;
	}
	else if (res == 3) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 3.0;
	}

	return std::make_pair(res, Weight);
}


inline Vector3f
findClosestPairExtrinsic(const Vector3f &o0_, const Vector3f &q0, const Vector3f &n0, const Vector3f &p0,
                         const Vector3f &o1_, const Vector3f &q1, const Vector3f &n1, const Vector3f &p1,
                         Float scale, Float invScale) {
    typedef Eigen::Matrix<Float, 3, 2> Matrix;
    Matrix M0 = (Matrix() << q0, n0.cross(q0)).finished();
    Matrix M1 = (Matrix() << q1, n1.cross(q1)).finished();
    Vector3f middle = middle_point(p0, n0, p1, n1);
    Vector3f rel0_ =  M0 * floor(Vector2f(M0.transpose() * (middle - o0_) * invScale)) * scale;
    Vector3f rel1_ =  M1 * floor(Vector2f(M1.transpose() * (middle - o1_) * invScale)) * scale;

    Float best_score = std::numeric_limits<Float>::infinity();
    Vector3f o1 = o1_;

    for (int i=0; i<16; ++i) {
        Vector3f rel0 = rel0_ + M0 * Vector2f((i & (1 << 0)) >> 0,
                                              (i & (1 << 1)) >> 1) * scale;
        Vector3f rel1 = rel1_ + M1 * Vector2f((i & (1 << 2)) >> 2,
                                              (i & (1 << 3)) >> 3) * scale;

        Float score = (o0_ + rel0 - (o1_ + rel1)).squaredNorm();

        if (score < best_score) {
            best_score = score;
            o1 = o1_ + rel1 - rel0;
        }
    }

    return o1;
}

inline std::pair<Vector3f, Vector2i>
findClosestPair(const Vector3f &o0,  const Vector3f &q0,  const Vector3f &n0, const Vector3f &p0,
                const Vector3f &o1_, const Vector3f &q1_, const Vector3f &n1, const Vector3f &p1,
                Float scale, Float invScale) {
    typedef Eigen::Matrix<Float, 3, 2> Matrix;

    Vector3f qn = (q0 + applyRotation(q0, n0, q1_, n1)).normalized();
    Vector3f middle = middle_point(p0, n0, p1, n1);
    Vector3f o1 = rotateVectorIntoPlane(o1_ - middle, n1, n0) + middle;

    Matrix M = (Matrix() << qn, n0.cross(qn)).finished();
    Vector2f rel = round(Vector2f(M.transpose() * ((o0 - o1) * invScale)));
    return std::make_pair(o1 + (M * rel) * scale, rel.cast<int>());
}
inline std::pair<int, Float> assignColorWeighted2D(const Vector3f &o0, const Vector3f &q0, const Vector3f &n0, const Vector3f &p0,
	const Vector3f &o1_, const Vector3f &q1_, const Vector3f &n1, const Vector3f &p1,
	Float scale, Float invScale) {
	typedef Eigen::Matrix<Float, 3, 2> Matrix;

	Vector3f qn = (q0 + applyRotation(q0, n0, q1_, n1)).normalized();
	Vector3f middle = middle_point(p0, n0, p1, n1);
	Vector3f o1 = rotateVectorIntoPlane(o1_ - middle, n1, n0) + middle;

	Matrix M = (Matrix() << qn, n0.cross(qn)).finished();
	//Vector2f rel = round();
	Vector2f vv = Vector2f(M.transpose() * ((o0 - o1) * invScale));
	for (uint32_t i = 0; i < 2; i++) vv[i] = std::abs(vv[i]);

	short res = 0; Float Weight=0; std::vector<bool> bits(2);
	for (uint32_t i = 0; i < 2; i++) {
		if (vv[i] > 0.5) { res++; bits[i] = true; }
		else bits[i] = false;
	}
	if (res == 0) { // it's a red! 
		for (int i = 0; i < 2; i++) Weight += vv[i] * vv[i];
	}
	else if (res == 1) {
		for (int i = 0; i < 2; i++) {
			if (bits[i]) { Weight += (1 - vv[i]) * (1 - vv[i]); }
			else Weight += vv[i] * vv[i];
		}
		Weight += 1.0;
	}
	else if (res == 2) {
		for (int i = 0; i < 2; i++)
			Weight += (1 - vv[i]) * (1 - vv[i]);
		Weight += std::sqrt(2.0);
	}

	return std::make_pair(res, Weight);
}
inline std::tuple<int, Float, Vector2f> posy2D_completeInfo(const Vector3f &o0, const Vector3f &q0, const Vector3f &n0, const Vector3f &p0,
	const Vector3f &o1_, const Vector3f &q1_, const Vector3f &n1, const Vector3f &p1,
	Float scale, Float invScale) {
	typedef Eigen::Matrix<Float, 3, 2> Matrix;

	Vector3f qn = (q0 + applyRotation(q0, n0, q1_, n1)).normalized();
	Vector3f middle = middle_point(p0, n0, p1, n1);
	Vector3f o1 = rotateVectorIntoPlane(o1_ - middle, n1, n0) + middle;

	Matrix M = (Matrix() << qn, n0.cross(qn)).finished();
	//Vector2f rel = round();
	Vector2f vv = Vector2f(M.transpose() * ((o0 - o1) * invScale));
	for (uint32_t i = 0; i < 2; i++) vv[i] = std::abs(vv[i]);

	short res = 0; Float Weight = 0; std::vector<bool> bits(2);
	for (uint32_t i = 0; i < 2; i++) {
		if (vv[i] > 0.5) { res++; bits[i] = true; }
		else bits[i] = false;
	}
	if (res == 0) { // it's a red! 
		for (int i = 0; i < 2; i++) Weight += vv[i] * vv[i];
	}
	else if (res == 1) {
		for (int i = 0; i < 2; i++) {
			if (bits[i]) { Weight += (1 - vv[i]) * (1 - vv[i]); }
			else Weight += vv[i] * vv[i];
		}
		Weight += 1.0;
	}
	else if (res == 2) {
		for (int i = 0; i < 2; i++)
			Weight += (1 - vv[i]) * (1 - vv[i]);
		Weight += std::sqrt(2.0);
	}

	return std::make_tuple(res, Weight, vv);
}
inline std::tuple<int, Float, Vector2f> posy2D_completeInfo_global(const Vector3f &o0, const Vector3f &o1, Float scale, Float invScale) {
	typedef Eigen::Matrix<Float, 3, 2> Matrix;
	Vector2f vv = ((round(o0) - round(o1)) * 1).segment(0, 2);
	for (uint32_t i = 0; i < 2; i++) vv[i] = std::abs(vv[i]);

	short res = 0; Float Weight = 0; std::vector<bool> bits(2);
	for (uint32_t i = 0; i < 2; i++) {
		if (vv[i] > 0.5) { res++; bits[i] = true; }
		else bits[i] = false;
	}
	if (res == 0) { // it's a red! 
		for (int i = 0; i < 2; i++) Weight += vv[i] * vv[i];
	}
	else if (res == 1) {
		for (int i = 0; i < 2; i++) {
			if (bits[i]) { Weight += (1 - vv[i]) * (1 - vv[i]); }
			else Weight += vv[i] * vv[i];
		}
		Weight += 1.0;
	}
	else if (res == 2) {
		for (int i = 0; i < 2; i++)
			Weight += (1 - vv[i]) * (1 - vv[i]);
		Weight += std::sqrt(2.0);
	}

	return std::make_tuple(res, Weight, vv);
}
inline Vector3f findClosest(const Vector3f &o, const Quaternion &q, 
                            const Vector3f &ref, Float scale, Float invScale) {
    Matrix3f M = q.toMatrix();
    Vector3f rel = M.transpose() * ((ref - o) * invScale);
    return o + (M * round(rel)) * scale;
}

inline Vector3f findClosest(const Vector3f &o, const Vector3f &q, const Vector3f &n,
                            const Vector3f &ref, Float scale, Float invScale) {
    typedef Eigen::Matrix<Float, 3, 2> Matrix;
    Matrix M = (Matrix() << q, n.cross(q)).finished();
    Vector2f rel = M.transpose() * ((ref - o) * invScale);
    return o + (M * round(rel)) * scale;
}


//projection
void projectPointOnQuad(const vector<Vector3d>& quad_vs, vector<Vector3d> & vs_normals, const Vector3d& p, Vector2d& uv, Vector3d& interpolP, Vector3d& interpolN);
void projectPointOnTriangle(const vector<Vector3d>& tri_vs, const vector<Vector3d> & vs_normals, const Vector3d& p, Vector2d& uv, Vector3d& interpolP, Vector3d& interpolN);
template <typename T>
T bilinear(const T& v1, const T& v2, const T& v3, const T& v4, const Vector2d& uv);
template <typename T>
T barycentric(const T& v1, const T& v2, const T& v3, const Vector2d& uv);
template <typename T>
bool num_equal(const T& x, const T& y, const double &precision);