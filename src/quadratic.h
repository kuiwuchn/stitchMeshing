#pragma once
#include "common.h"
typedef float Scalar;

/* 
   xT*A*x + b*x + c
*/
class Quadric
{
public:
    typedef Scalar ScalarType;
	ScalarType a[6];		// Matrix 3x3: a11 a12 a13 a22 a23 a33
	ScalarType b[3];		// b vector
	ScalarType c;			// scalar

	inline Quadric() { c = -1; }

    // if you call initByXXXX, the quadric will be the DISTANCE of a given x from that XXXX

    // plane normal n, point p, k = n*p
	void initByPlane( const Vector3f & n, Scalar k )
	{
		a[0] =  n[0]*n[0];	// a11
		a[1] =  n[1]*n[0];	// a12 (=a21)
		a[2] =  n[2]*n[0];	// a13 (=a31)
		a[3] =  n[1]*n[1];	// a22
		a[4] =  n[2]*n[1];	// a23 (=a32)
		a[5] =  n[2]*n[2];	// a33
		b[0] = (-2.0)*k*n[0];
		b[1] = (-2.0)*k*n[1];
		b[2] = (-2.0)*k*n[2];
		c    =  k*k;
	}
    
    // do we need initByLine ? not for now.    
	/*
  void initByLine( const LineType & r ) // Init dato un raggio
  {
    ScalarType K = (ScalarType)(r.Origin()*r.Direction());
    a[0] = (ScalarType)1.0-r.Direction()[0]*r.Direction()[0]; // a11
    a[1] = (ScalarType)-r.Direction()[0]*r.Direction()[1]; // a12 (=a21)
    a[2] = (ScalarType)-r.Direction()[0]*r.Direction()[2]; // a13 (=a31)
    a[3] = (ScalarType)1.0-r.Direction()[1]*r.Direction()[1]; // a22
    a[4] = (ScalarType)-r.Direction()[1]*r.Direction()[2]; // a23 (=a32)
    a[5] = (ScalarType)1.0-r.Direction()[2]*r.Direction()[2]; // a33
    b[0] = (ScalarType)2.0*(r.Direction()[0]*K - r.Origin()[0]);
    b[1] = (ScalarType)2.0*(r.Direction()[1]*K - r.Origin()[1]);
    b[2] = (ScalarType)2.0*(r.Direction()[2]*K - r.Origin()[2]);
    c = -K*K + (ScalarType)(r.Origin()*r.Origin());
  }
 */

  void initByPoint( const Vector3f & p ) {
        // a is the ide
		a[0] =  1;	// a11
		a[1] =  0;	// a12 (=a21)
		a[2] =  0;	// a13 (=a31)
		a[3] =  1;	// a22
		a[4] =  0;	// a23 (=a32)
		a[5] =  1;	// a33
		b[0] = (-2.0)*p[0];
		b[1] = (-2.0)*p[1];
		b[2] = (-2.0)*p[2];
		c    = p.dot(p); 
  }

    // inits the quadric to be the wieghted distance:
    // apply(x) = sqrardDistanceFromPoint * weight + squaredDistFromPlane
    // use small weight! 0.001.  ACtually you are using FLOATS so maybe not too little 0.001 is fine
  void initByPointAndNormal(const Vector3f & p, const Vector3f & n, Scalar weight ){
        initByPoint( p );
        (*this) *= weight;
        
        //Quadric b;
        //b.initByPlane( n, n.dot(p));
        //(*this) += b;
  }
    
  /* 1) assign a quadric to each vertex.
     2) when init: each vertex, initByPointAndNormal
        (maybe, debug: verify that apply(vertexPos) == 0 )
        For internal points: use 0 as normal. 
     3) when collapse v0 and v1 into v_new
        v_new.quadric = v0 + v1   // no wiughting acutally
        // put the vertex pos in the minimum of the quadric:
        v_new.quadric.getMinimum( v_new.pos );
    DONE
    */
    
    

  void setZero() {
		a[0] = 0;
		a[1] = 0;
		a[2] = 0;
		a[3] = 0;
		a[4] = 0;
		a[5] = 0;
		b[0] = 0;
		b[1] = 0;
		b[2] = 0;
		c    = 0;
    }

void operator = ( const Quadric & q )
	{
		a[0] = q.a[0];
		a[1] = q.a[1];
		a[2] = q.a[2];
		a[3] = q.a[3];
		a[4] = q.a[4];
		a[5] = q.a[5];
		b[0] = q.b[0];
		b[1] = q.b[1];
		b[2] = q.b[2];
		c    = q.c;
	}

  void operator += ( const Quadric & q )
	{
		a[0] += q.a[0];
		a[1] += q.a[1];
		a[2] += q.a[2];
		a[3] += q.a[3];
		a[4] += q.a[4];
		a[5] += q.a[5];
		b[0] += q.b[0];
		b[1] += q.b[1];
		b[2] += q.b[2];
		c    += q.c;
	}

    // just for debugging...
    Scalar apply( const Vector3f & p ) const	
	{
	  return (
        p[0]*p[0]*a[0] + 2*p[0]*p[1]*a[1] + 2*p[0]*p[2]*a[2] + p[0]*b[0]
		   	               +   p[1]*p[1]*a[3] + 2*p[1]*p[2]*a[4] + p[1]*b[1]
			                                  +   p[2]*p[2]*a[5] + p[2]*b[2]	+ c);
	}

bool getMinimum(Vector3f &result) const
{
  ScalarType c0=-b[0]/2;
  ScalarType c1=-b[1]/2;
  ScalarType c2=-b[2]/2;

  ScalarType t125 = a[4]*a[1];
  ScalarType t124 = a[4]*a[4]-a[3]*a[5];
  ScalarType t123 = -a[1]*a[5]+a[4]*a[2];
  ScalarType t122 = a[2]*a[3]-t125;
  ScalarType t121 = -a[2]*a[1]+a[0]*a[4];
  ScalarType t120 = a[2]*a[2];
  ScalarType t119 = a[1]*a[1];
  ScalarType t117 = 1.0/(-a[3]*t120+2*a[2]*t125-t119*a[5]-t124*a[0]);
  result[0] = -(t124*c0+t122*c2-t123*c1)*t117;
  result[1] = (t123*c0-t121*c2+(-t120+a[0]*a[5])*c1)*t117;
  result[2] = -(t122*c0+(t119-a[0]*a[3])*c2+t121*c1)*t117;
  return true;
}
    
  void operator *= ( const ScalarType & w )
	{
		a[0] *= w;
		a[1] *= w;
		a[2] *= w;
		a[3] *= w;
		a[4] *= w;
		a[5] *= w;
		b[0] *= w;
		b[1] *= w;
		b[2] *= w;
		c    *= w;
	}
    
  Quadric operator * ( const ScalarType & w )
	{
        Quadric res = *this;
        res*=w;
		return res;
	}

  Quadric operator + ( const Quadric & w )
	{
        Quadric res = *this;	
        res+=w;
		return res;
	}
};