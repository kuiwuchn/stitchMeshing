// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
///
/// \file		cyPoint.h 
/// \author		Cem Yuksel
/// \version	1.4
/// \date		September 3, 2012
///
/// \brief 2D, 3D and 4D point classes.
///
///
/// @copydoc cyPoint3f
///
/// cyPoint3f is 3D point class
///
/// You can use point classes with OpenGL as
///
/// glVertex3f( myvector.x, myvector.y, myvector.z );
///
/// or
///
/// glVertex3fv( &myvector.x );
///
//-------------------------------------------------------------------------------

#ifndef _CY_POINT_H_INCLUDED_
#define _CY_POINT_H_INCLUDED_

//-------------------------------------------------------------------------------

#include <math.h>

//-------------------------------------------------------------------------------

/// 2D point class

class cyPoint2f
{
	friend cyPoint2f operator+( const float v, const cyPoint2f &pt ) { return pt+v; }		///< Addition with a constant
	friend cyPoint2f operator-( const float v, const cyPoint2f &pt ) { return -(pt-v); }	///< Subtraction from a constant
	friend cyPoint2f operator*( const float v, const cyPoint2f &pt ) { return pt*v; }		///< Multiplication with a constant

public:

	float x, y;

	///@name Constructors
	cyPoint2f() {}
	cyPoint2f( float _x, float _y ) { x=_x; y=_y; }
	cyPoint2f( const float *pt ) { x=pt[0]; y=pt[1]; }
	cyPoint2f( const cyPoint2f &pt ) { x=pt.x; y=pt.y; }

	///@name Set & Get value functions
	cyPoint2f& Zero() { x=0; y=0; return *this; }						///< Sets x and y coordinates as zero
	cyPoint2f& Set( float _x, float _y ) { x=_x; y=_y; return *this; }	///< Sets x and y coordinates as given
	cyPoint2f& Set( const float *v ) { x=v[0]; y=v[1]; return *this; }	///< Sets x and y coordinates using the values in the given array
	void GetValue( float *v ) const { v[0]=x; v[1]=y; }					///< Puts x and y values into the array

	///@name Length and Normalize functions
	void		Normalize()		{ float s = (1.0f)/Length(); *this *= s; }
	cyPoint2f	GetNormalized() const { float s = (1.0f)/Length(); return *this * s; }
	float		LengthSquared()	const { return x*x + y*y; }
	float		Length()		const { return (float) sqrt(LengthSquared()); }

	///@name Limit functions
	void ClampMinMax( float minValue, float maxValue ) { ClampMin(minValue); ClampMax(maxValue); }
	void ClampMin( float n ) { if(x<n)x=n; if(y<n)y=n; }
	void ClampMax( float n ) { if(x>n)x=n; if(y>n)y=n; }

	///@name Unary operators
	cyPoint2f operator-() const { return cyPoint2f(-x,-y); } 
	cyPoint2f operator+() const { return *this; }

	///@name Binary operators
	cyPoint2f operator+( const cyPoint2f &pt ) const { return cyPoint2f(x+pt.x, y+pt.y); }
	cyPoint2f operator-( const cyPoint2f &pt ) const { return cyPoint2f(x-pt.x, y-pt.y); }
	cyPoint2f operator*( const cyPoint2f &pt ) const { return cyPoint2f(x*pt.x, y*pt.y); }
	cyPoint2f operator/( const cyPoint2f &pt ) const { return cyPoint2f(x/pt.x, y/pt.y); }
	cyPoint2f operator+(float n) const { return cyPoint2f(x+n, y+n); }
	cyPoint2f operator-(float n) const { return cyPoint2f(x-n, y-n); }
	cyPoint2f operator*(float n) const { return cyPoint2f(x*n, y*n); }
	cyPoint2f operator/(float n) const { return cyPoint2f(x/n, y/n); }

	///@name Assignment operators
	cyPoint2f& operator+=( const cyPoint2f &pt ) { x+=pt.x; y+=pt.y; return *this; }
	cyPoint2f& operator-=( const cyPoint2f &pt ) { x-=pt.x; y-=pt.y; return *this; }
	cyPoint2f& operator*=( const cyPoint2f &pt ) { x*=pt.x; y*=pt.y; return *this; }
	cyPoint2f& operator/=( const cyPoint2f &pt ) { x/=pt.x; y/=pt.y; return *this; }
	cyPoint2f& operator+=(float n) { x+=n; y+=n; return *this; }
	cyPoint2f& operator-=(float n) { x-=n; y-=n; return *this; }
	cyPoint2f& operator*=(float n) { x*=n; y*=n; return *this; }
	cyPoint2f& operator/=(float n) { x/=n; y/=n; return *this; }

	///@name Test operators
	int operator==( const cyPoint2f& pt ) const { return ( (pt.x==x) && (pt.y==y) ); }
	int operator!=( const cyPoint2f& pt ) const { return ( (pt.x!=x) || (pt.y!=y) ); }

	///@name Access operators
	float& operator[]( int i ) { return (&x)[i]; }
	float  operator[]( int i ) const { return (&x)[i]; }

	///@name Cross product and dot product
	float Cross	   ( const cyPoint2f &pt ) const { return x*pt.y-y*pt.x; }		///< Cross product
	float operator^( const cyPoint2f &pt ) const { return Cross(pt); }			///< Cross product operator
	float Dot	   ( const cyPoint2f &pt ) const { return x*pt.x + y*pt.y; }	///< Dot product
	float operator%( const cyPoint2f &pt ) const { return Dot(pt); }			///< Dot product operator

};




//-------------------------------------------------------------------------------


/// 3D point class

class cyPoint3f
{
	friend cyPoint3f operator+( const float v, const cyPoint3f &pt ) { return pt+v; }		///< Addition with a constant
	friend cyPoint3f operator-( const float v, const cyPoint3f &pt ) { return -(pt-v); }	///< Subtraction from a constant
	friend cyPoint3f operator*( const float v, const cyPoint3f &pt ) { return pt*v; }		///< Multiplication with a constant

public:

	float x, y, z;

	///@name Constructors
	cyPoint3f() { }
	cyPoint3f( float _x, float _y, float _z ) { x=_x; y=_y; z=_z; }
	cyPoint3f( const float *pt ) { x=pt[0]; y=pt[1]; z=pt[2]; }
	cyPoint3f( const cyPoint3f &pt ) { x=pt.x; y=pt.y; z=pt.z; }
	cyPoint3f( const cyPoint2f &pt ) { x=pt.x; y=pt.y; z=0.0f; }
	cyPoint3f( const cyPoint2f &pt, float _z ) { x=pt.x; y=pt.y; z=_z; }

	///@name Set & Get value functions
	cyPoint3f& Zero() { x=0; y=0; z=0; return *this; }									///< Sets x, y and z coordinates as zero
	cyPoint3f& Set( float _x, float _y, float _z ) { x=_x; y=_y; z=_z; return *this; }	///< Sets x, y and z coordinates as given
	cyPoint3f& Set( const float *v ) { x=v[0]; y=v[1]; z=v[2]; return *this; }			///< Sets x, y and z coordinates using the values in the given array
	void GetValue( float *v ) const { v[0]=x; v[1]=y; v[2]=z; }							///< Puts x, y and z values into the array

	///@name Length and Normalize functions
	void		Normalize()		{ float s = (1.0f)/Length(); *this *= s; }
	cyPoint3f	GetNormalized()	const { float s = (1.0f)/Length(); return *this * s; }
	float		LengthSquared()	const { return x*x + y*y + z*z; }
	float		Length()		const { return (float) sqrt(LengthSquared()); }

	///@name Limit functions
	void ClampMinMax( float min, float max ) { ClampMin(min); ClampMax(max); }
	void ClampMin( float n ) { if(x<n)x=n; if(y<n)y=n; if(z<n)z=n; }
	void ClampMax( float n ) { if(x>n)x=n; if(y>n)y=n; if(z>n)z=n; }

	///@name Unary operators
	cyPoint3f operator-() const { return cyPoint3f(-x,-y,-z); } 
	cyPoint3f operator+() const { return *this; }

	///@name Binary operators
	cyPoint3f operator+( const cyPoint3f &pt ) const { return cyPoint3f(x+pt.x, y+pt.y, z+pt.z); }
	cyPoint3f operator-( const cyPoint3f &pt ) const { return cyPoint3f(x-pt.x, y-pt.y, z-pt.z); }
	cyPoint3f operator*( const cyPoint3f &pt ) const { return cyPoint3f(x*pt.x, y*pt.y, z*pt.z); }
	cyPoint3f operator/( const cyPoint3f &pt ) const { return cyPoint3f(x/pt.x, y/pt.y, z/pt.z); }
	cyPoint3f operator+(float n) const { return cyPoint3f(x+n, y+n, z+n); }
	cyPoint3f operator-(float n) const { return cyPoint3f(x-n, y-n, z-n); }
	cyPoint3f operator*(float n) const { return cyPoint3f(x*n, y*n, z*n); }
	cyPoint3f operator/(float n) const { return cyPoint3f(x/n, y/n, z/n); }

	///@name Assignment operators
	cyPoint3f& operator+=( const cyPoint3f &pt ) { x+=pt.x; y+=pt.y; z+=pt.z; return *this; }
	cyPoint3f& operator-=( const cyPoint3f &pt ) { x-=pt.x; y-=pt.y; z-=pt.z; return *this; }
	cyPoint3f& operator*=( const cyPoint3f &pt ) { x*=pt.x; y*=pt.y; z*=pt.z; return *this; }
	cyPoint3f& operator/=( const cyPoint3f &pt ) { x/=pt.x; y/=pt.y; z/=pt.z; return *this; }
	cyPoint3f& operator+=(float n) { x+=n; y+=n; z+=n; return *this; }
	cyPoint3f& operator-=(float n) { x-=n; y-=n; z-=n; return *this; }
	cyPoint3f& operator*=(float n) { x*=n; y*=n; z*=n; return *this; }
	cyPoint3f& operator/=(float n) { x/=n; y/=n; z/=n; return *this; }

	///@name Test operators
	int operator==( const cyPoint3f& pt ) const { return ( (pt.x==x) && (pt.y==y) && (pt.z==z) ); }
	int operator!=( const cyPoint3f& pt ) const { return ( (pt.x!=x) || (pt.y!=y) || (pt.z!=z) ); }

	///@name Access operators
	float& operator[]( int i ) { return (&x)[i]; }
	float  operator[]( int i ) const { return (&x)[i]; }

	///@name Cross product and dot product
	cyPoint3f	Cross	 ( const cyPoint3f pt ) const { return cyPoint3f(y*pt.z-z*pt.y, z*pt.x-x*pt.z, x*pt.y-y*pt.x); }	///< Cross product
	cyPoint3f	operator^( const cyPoint3f pt ) const { return Cross(pt); }					///< Cross product
	float		Dot		 ( const cyPoint3f pt ) const { return x*pt.x + y*pt.y + z*pt.z; }	///< Dot product
	float		operator%( const cyPoint3f pt ) const { return Dot(pt); }					///< Dot product

	///@name Conversion Methods
	cyPoint2f	XY() const { return cyPoint2f(x,y); }
};

//-------------------------------------------------------------------------------


/// 4D point class
class cyPoint4f
{
	friend cyPoint4f operator+( const float v, const cyPoint4f &pt ) { return pt+v; }		///< Addition with a constant
	friend cyPoint4f operator-( const float v, const cyPoint4f &pt ) { return -(pt-v); }	///< Subtraction from a constant
	friend cyPoint4f operator*( const float v, const cyPoint4f &pt ) { return pt*v; }		///< Multiplication with a constant

public:

	float x, y, z, w;

	///@name Constructors
	cyPoint4f() { }
	cyPoint4f( float _x, float _y, float _z, float _w ) { x=_x; y=_y; z=_z; w=_w; }
	cyPoint4f( const float *pt ) { x=pt[0]; y=pt[1]; z=pt[2]; w=pt[3]; }
	cyPoint4f( const cyPoint4f &pt ) { x=pt.x; y=pt.y; z=pt.z; w=pt.w; }
	cyPoint4f( const cyPoint3f &pt ) { x=pt.x; y=pt.y; z=pt.z; w=1.0f; }
	cyPoint4f( const cyPoint3f &pt, float _w ) { x=pt.x; y=pt.y; z=pt.z; w=_w; }

	///@name Set & Get value functions
	cyPoint4f& Zero() { x=0; y=0; z=0; w=0; return *this; }												///< Sets x, y, z and w coordinates as zero
	cyPoint4f& Set( float _x, float _y, float _z, float _w ) { x=_x; y=_y; z=_z; w=_w; return *this; }	///< Sets x, y, z and w coordinates as given
	cyPoint4f& Set( const float *v ) { x=v[0]; y=v[1]; z=v[2]; w=v[3]; return *this; }					///< Sets x, y, z and w coordinates using the values in the given array
	void GetValue( float *v ) const { v[0]=x; v[1]=y; v[2]=z; v[3]=w; }									///< Puts x, y, z and w values into the array

	///@name Length and Normalize functions
	void		Normalize()		{ float s = (1.0f)/Length(); *this *= s; }
	cyPoint4f	GetNormalized() const { float s = (1.0f)/Length(); return *this * s; }
	float		LengthSquared()	const { return x*x + y*y + z*z + w*w; }
	float		Length()		const { return (float) sqrt(LengthSquared()); }

	///@name Limit functions
	void ClampMinMax( float min, float max ) { ClampMin(min); ClampMax(max); }
	void ClampMin( float n ) { if(x<n)x=n; if(y<n)y=n; if(z<n)z=n; if(w<n)w=n; }
	void ClampMax( float n ) { if(x>n)x=n; if(y>n)y=n; if(z>n)z=n; if(w>n)w=n; }

	///@name Unary operators
	cyPoint4f operator-() const { return cyPoint4f(-x,-y,-z,-w); } 
	cyPoint4f operator+() const { return *this; }

	///@name Binary operators
	cyPoint4f operator+( const cyPoint4f &pt ) const { return cyPoint4f(x+pt.x, y+pt.y, z+pt.z, w+pt.w); }
	cyPoint4f operator-( const cyPoint4f &pt ) const { return cyPoint4f(x-pt.x, y-pt.y, z-pt.z, w-pt.w); }
	cyPoint4f operator*( const cyPoint4f &pt ) const { return cyPoint4f(x*pt.x, y*pt.y, z*pt.z, w*pt.w); }
	cyPoint4f operator/( const cyPoint4f &pt ) const { return cyPoint4f(x/pt.x, y/pt.y, z/pt.z, w/pt.w); }
	cyPoint4f operator+(float n) const { return cyPoint4f(x+n, y+n, z+n, w+n); }
	cyPoint4f operator-(float n) const { return cyPoint4f(x-n, y-n, z-n, w-n); }
	cyPoint4f operator*(float n) const { return cyPoint4f(x*n, y*n, z*n, w*n); }
	cyPoint4f operator/(float n) const { return cyPoint4f(x/n, y/n, z/n, w/n); }

	///@name Assignment operators
	cyPoint4f& operator+=( const cyPoint4f &pt ) { x+=pt.x; y+=pt.y; z+=pt.z; w+=pt.w; return *this; }
	cyPoint4f& operator-=( const cyPoint4f &pt ) { x-=pt.x; y-=pt.y; z-=pt.z; w-=pt.w; return *this; }
	cyPoint4f& operator*=( const cyPoint4f &pt ) { x*=pt.x; y*=pt.y; z*=pt.z; w*=pt.w; return *this; }
	cyPoint4f& operator/=( const cyPoint4f &pt ) { x/=pt.x; y/=pt.y; z/=pt.z; w/=pt.w; return *this; }
	cyPoint4f& operator+=(float n) { x+=n; y+=n; z+=n; w+=n; return *this; }
	cyPoint4f& operator-=(float n) { x-=n; y-=n; z-=n; w-=n; return *this; }
	cyPoint4f& operator*=(float n) { x*=n; y*=n; z*=n; w*=n; return *this; }
	cyPoint4f& operator/=(float n) { x/=n; y/=n; z/=n; w/=n; return *this; }

	///@name Test operators
	int operator==( const cyPoint4f& pt ) const { return ( (pt.x==x) && (pt.y==y) && (pt.z==z) && (pt.w==w) ); }
	int operator!=( const cyPoint4f& pt ) const { return ( (pt.x!=x) || (pt.y!=y) || (pt.z!=z) || (pt.w!=w) ); }

	///@name Access operators
	float& operator[]( int i ) { return (&x)[i]; }
	float  operator[]( int i ) const { return (&x)[i]; }

	///@ Dot product
	float		Dot		 ( const cyPoint4f pt ) const { return x*pt.x + y*pt.y + z*pt.z + w*pt.w; }	///< Dot product
	float		operator%( const cyPoint4f pt ) const { return Dot(pt); }							///< Dot product

	///@name Conversion Methods
	cyPoint2f	XY() const { return cyPoint2f(x,y); }
	cyPoint3f	XYZ() const { return cyPoint3f(x,y,z); }
};

//-------------------------------------------------------------------------------


namespace cy {
	typedef cyPoint2f Point2f;
	typedef cyPoint3f Point3f;
	typedef cyPoint4f Point4f;
}


//-------------------------------------------------------------------------------

#endif

