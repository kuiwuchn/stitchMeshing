#pragma once

#if defined(_WIN32)
    #define NOMINMAX
    #pragma warning(disable: 4244 4018 4100 4610 4510 4127 4512 4146 4267 4503 4800 4706)
#endif

#define EIGEN_DONT_PARALLELIZE
#define EIGEN_DEFAULT_IO_FORMAT \
    Eigen::IOFormat(7, Eigen::DontAlignCols, ", ", ";\n", "", "", "[", "]")
#define GRAIN_SIZE 1024

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <tbb/tbb.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <atomic>

/* Application precision -- can be set to single or double precision */
#if defined(SINGLE_PRECISION)
typedef float Float;
#else
typedef double Float;
#endif

/* Useful Eigen typedefs based on the current precision */
typedef Eigen::Matrix<int32_t, 2, 1>                            Vector2i;
typedef Eigen::Matrix<int32_t, 3, 1>                            Vector3i;
typedef Eigen::Matrix<int32_t, 4, 1>                            Vector4i;
typedef Eigen::Matrix<uint32_t, 2, 1>                           Vector2u;
typedef Eigen::Matrix<uint32_t, 3, 1>                           Vector3u;
typedef Eigen::Matrix<uint32_t, 4, 1>                           Vector4u;
typedef Eigen::Matrix<uint8_t, 4, 1>                            Vector4u8;
typedef Eigen::Matrix<Float, 2, 1>                              Vector2f;
typedef Eigen::Matrix<Float, 3, 1>                              Vector3f;
typedef Eigen::Matrix<double, 2, 1>                              Vector2d;
typedef Eigen::Matrix<double, 3, 1>                              Vector3d;
typedef Eigen::Matrix<Float, 4, 1>                              Vector4f;
typedef Eigen::Matrix<int32_t, Eigen::Dynamic, 1>               VectorXi;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, 1>              VectorXu;
typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, 1>               VectorXu8;
typedef Eigen::Matrix<bool, Eigen::Dynamic, 1>                  VectorXb;
typedef Eigen::Matrix<Float, Eigen::Dynamic, 1>                 VectorXf;
typedef Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic>  MatrixXi;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXu;
typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>  MatrixXu8;
typedef Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>    MatrixXf;
typedef Eigen::Matrix<Float, 2, 2>                              Matrix2f;
typedef Eigen::Matrix<Float, 3, 3>                              Matrix3f;
typedef Eigen::Matrix<int, 3, 3>                                Matrix3i;
typedef Eigen::Matrix<Float, 4, 4>                              Matrix4f;
typedef Eigen::SparseMatrix<Float, Eigen::RowMajor>             SMatrix;
typedef Eigen::Triplet<Float>                                   Triplet;

using Eigen::Quaternionf;

using std::cout;
using std::cerr;
using std::endl;

typedef std::tuple<uint32_t, uint32_t, bool, Float, short, uint32_t, short, uint32_t> tuple_E;//v0, v1, boundary, energy, color, edge index, xy/yz/xz plane, timestamp
typedef std::tuple<uint32_t, uint32_t, uint32_t, bool, uint32_t> tuple_F;//v0, v1, v2, boundary, face index; for tet face
struct LessThan { bool operator()(const tuple_E& lhs, const tuple_E& rhs) const { return (std::get<3>(lhs) > std::get<3>(rhs)); } };

struct statistics {
	double hex_ratio;//hexa number / total element num
	int tN, tetN, hN, pN;//hexa number, total element num
	std::vector<double> timings;//tetrahedralizing+hierarchy-tree construction, rosy optimizaiton, posy optimization, and mesh extraction. 
	std::vector<uint32_t> polyhedral_nums;
	std::vector<double> polyhedral_ratios;
};

#define PAI 3.1415926585
#define Precision_Pro 1.e-5
																	  //for volume
const uint8_t tet_faces[4][3] = { { 1, 0, 2 },{ 3, 2, 0 },{ 1, 2, 3 },{ 0, 1, 3 } };
const uint8_t tet_face_edge[4][3] = { { 0, 1, 3 },{ 5, 1, 2 },{ 3, 5, 4 },{ 0, 4, 2 } };//
const uint32_t tet_edges[6][2] = { { 0, 1 },{ 0, 2 },{ 0, 3 },{ 1, 2 },{ 1, 3 },{ 2, 3 } };


enum Edge_tag {
	R = 0,
	B,
	D,
	H
};

inline std::string str_tolower(std::string str) {
	std::transform(str.begin(), str.end(), str.begin(), ::tolower);
	return str;
}

inline uint32_t str_to_uint32_t(const std::string &str) {
    char *end_ptr = nullptr;
    uint32_t result = (uint32_t) strtoul(str.c_str(), &end_ptr, 10);
    if (*end_ptr != '\0')
        throw std::runtime_error("Could not parse unsigned integer \"" + str + "\"");
    return result;
}

inline uint32_t str_to_int32_t(const std::string &str) {
    char *end_ptr = nullptr;
    int32_t result = (int32_t) strtol(str.c_str(), &end_ptr, 10);
    if (*end_ptr != '\0')
        throw std::runtime_error("Could not parse signed integer \"" + str + "\"");
    return result;
}

inline Float str_to_float(const std::string &str) {
    char *end_ptr = nullptr;
    Float result = (Float) strtod(str.c_str(), &end_ptr);
    if (*end_ptr != '\0')
        throw std::runtime_error("Could not parse floating point value \"" + str + "\"");
    return result;
}

inline Vector3f hsv_to_rgb(Float h, Float s, Float v) {
    if (s == 0.f) { // achromatic (grey)
        return Vector3f::Constant(v);
    }
    h *= 6;
    int i = std::floor(h);
    Float f = h - i; // fractional part of h
    Float p = v * (1 - s);
    Float q = v * (1 - s * f);
    Float t = v * (1 - s * (1 - f));
    switch (i) {
        case 0:  return Vector3f(v, t, p); break;
        case 1:  return Vector3f(q, v, p); break;
        case 2:  return Vector3f(p, v, t); break;
        case 3:  return Vector3f(p, q, v); break;
        case 4:  return Vector3f(t, p, v); break;
        default: return Vector3f(v, p, q); break;
    }
}

inline std::string timeString(double time, bool precise = false) {
    if (std::isnan(time) || std::isinf(time))
        return "inf";

    std::string suffix = "ms";
    if (time > 1000) {
        time /= 1000; suffix = "s";
        if (time > 60) {
            time /= 60; suffix = "m";
            if (time > 60) {
                time /= 60; suffix = "h";
                if (time > 12) {
                    time /= 12; suffix = "d";
                }
            }
        }
    }

    std::ostringstream os;
    os << std::setprecision(precise ? 4 : 1)
       << std::fixed << time << suffix;

    return os.str();
}

inline std::string memString(size_t size, bool precise = false) {
    double value = (double) size;
    const char *suffixes[] = {
        "B", "KiB", "MiB", "GiB", "TiB", "PiB"
    };
    int suffix = 0;
    while (suffix < 5 && value > 1024.0f) {
        value /= 1024.0f; ++suffix;
    }

    std::ostringstream os;
    os << std::setprecision(suffix == 0 ? 0 : (precise ? 4 : 1))
       << std::fixed << value << " " << suffixes[suffix];

    return os.str();
}

inline bool atomicCompareAndExchange(volatile uint32_t *v, uint32_t newValue, uint32_t oldValue) {
#if defined(WIN32)
    return _InterlockedCompareExchange(
        reinterpret_cast<volatile long *>(v), (long) newValue, (long) oldValue) == (long) oldValue;
#else
    return __sync_bool_compare_and_swap(v, oldValue, newValue);
#endif
}

inline void coordinate_system(const Vector3f &a, Vector3f &b, Vector3f &c) {
    if (std::abs(a.x()) > std::abs(a.y())) {
        Float invLen = 1.0f / std::sqrt(a.x() * a.x() + a.z() * a.z());
        c = Vector3f(a.z() * invLen, 0.0f, -a.x() * invLen);
    } else {
        Float invLen = 1.0f / std::sqrt(a.y() * a.y() + a.z() * a.z());
        c = Vector3f(0.0f, a.z() * invLen, -a.y() * invLen);
    }
    b = c.cross(a);
}

template <typename T> T round(const T &value) {
    return value.unaryExpr([](Float f) { return std::round(f); });
}

template <typename T> T floor(const T &value) {
    return value.unaryExpr([](Float f) { return std::floor(f); });
}

template <typename T> T ceil(const T &value) {
    return value.unaryExpr([](Float f) { return std::ceil(f); });
}

inline Vector3f rotateVectorIntoPlane(Vector3f q, const Vector3f &source_normal, const Vector3f &target_normal) {
    const Float cosTheta = source_normal.dot(target_normal);
    if (cosTheta < 0.9999f) {
        Vector3f axis = source_normal.cross(target_normal);
        q = q * cosTheta + axis.cross(q) +
             axis * (axis.dot(q) * (1.0f - cosTheta) / axis.dot(axis));
    }
    return q;
}

inline std::vector<std::string> &str_tokenize(const std::string &s, char delim, std::vector<std::string> &elems, bool include_empty = false) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim))
		if (!item.empty() || include_empty)
			elems.push_back(item);
	return elems;
}

inline std::vector<std::string> str_tokenize(const std::string &s, char delim, bool include_empty) {
	std::vector<std::string> elems;
	str_tokenize(s, delim, elems, include_empty);
	return elems;
}
inline Float signum(Float value) {
	return std::copysign((Float)1, value);
}