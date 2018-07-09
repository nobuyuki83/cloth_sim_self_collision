//
//  vector3d.h
//  self_contact
//
//  Created by Nobuyuki Umetani on 11/8/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//

#ifndef vector3d_h
#define vector3d_h

#include <cassert>
#include <math.h>
#include <iostream>

#define NEARLY_ZERO 1.e-16


// rule about naming, the method starts "Set" change it self (not const)
class CVector3D;


//! @{
CVector3D operator+(const CVector3D&, const CVector3D&);
CVector3D operator-(const CVector3D&, const CVector3D&);
CVector3D operator*(double, const CVector3D&);
CVector3D operator*(const CVector3D&, double);
double operator*(const CVector3D&, const CVector3D& );
CVector3D operator^(const CVector3D&, const CVector3D&);
//! @}

//! 3 dimentional vector class
class CVector3D
{
public:
	CVector3D(double vx, double vy, double vz) : x(vx), y(vy), z(vz){}
	CVector3D(): x(0.0), y(0.0), z(0.0){}
	CVector3D(const CVector3D& rhs){
		x = rhs.x; y = rhs.y; z = rhs.z;
	}
	virtual ~CVector3D(){}
  
	void SetVector(double vx, double vy, double vz){ x = vx; y = vy; z = vz; }
  
	inline const CVector3D operator-() const{ return -1.0*(*this); }
	inline const CVector3D operator+() const{ return *this; }
	inline CVector3D& operator=(const CVector3D& rhs){
		if( this != &rhs ){ x = rhs.x; y = rhs.y; z = rhs.z; }
		return *this;
	}
	inline CVector3D& operator+=(const CVector3D& rhs){
		x += rhs.x; y += rhs.y; z += rhs.z;
		return *this;
	}
	inline CVector3D& operator-=(const CVector3D& rhs){
		x -= rhs.x; y -= rhs.y; z -= rhs.z;
		return *this;
	}
	inline CVector3D& operator*=(double d){
		x *= d; y *= d; z *= d;
		return *this;
	}
	inline CVector3D& operator/=(double d){
		if( fabs(d) < NEARLY_ZERO ){ return *this; }
		x /= d; y /= d; z /= d;
		return *this;
	}
  inline double operator[](int i) const{
    if( i == 0 ) return x;
    if( i == 1 ) return y;
    if( i == 2 ) return z;
    return 0;
  }
  inline double& operator[](int i){
    if( i == 0 ) return x;
    if( i == 1 ) return y;
    if( i == 2 ) return z;
    assert(0);
    return x;
  }
	inline CVector3D operator+(){ return *this; }
	inline CVector3D operator-(){ return CVector3D(-x,-y,-z); }
  
	friend bool operator==(const CVector3D&, const CVector3D&);
	friend bool operator!=(const CVector3D&, const CVector3D&);
  
	friend CVector3D Cross(const CVector3D&, const CVector3D&);
	friend double Dot(const CVector3D&, const CVector3D&);
  
	inline double Length()  const{ return sqrt( x*x+y*y+z*z ); }
	inline double DLength() const{ return x*x+y*y+z*z; }
  
  void SetNormalizedVector()
  {
    double invmag = 1.0/Length();
    x *= invmag;
    y *= invmag;
    z *= invmag;
  }
  
  void SetZero()
  {
    x = 0.0;
    y = 0.0;
    z = 0.0;
  }
public:
	double x;	//!< x axis coordinate
	double y;	//!< y axis coordinate
	double z;	//!< z axis coordinate
};

//! add
inline CVector3D operator + (const CVector3D& lhs, const CVector3D& rhs){
	CVector3D temp = lhs;
	temp += rhs;
	return temp;
}

//! subtract
inline CVector3D operator - (const CVector3D& lhs, const CVector3D& rhs){
	CVector3D temp = lhs;
	temp -= rhs;
	return temp;
}

//! scale
inline CVector3D operator * (double d, const CVector3D& rhs){
	CVector3D temp = rhs;
	temp *= d;
	return temp;
}

//! scale
inline CVector3D operator * (const CVector3D& vec, double d){
  CVector3D temp = vec;
  temp *= d;
  return temp;
}

//! divide by real number
inline CVector3D operator / (const CVector3D& vec, double d){
	CVector3D temp = vec;
	temp /= d;
	return temp;
}


//! mult
inline double operator * (const CVector3D& lhs, const CVector3D& rhs){
	return Dot(lhs,rhs);
}


//! mult
inline CVector3D operator ^ (const CVector3D& lhs, const CVector3D& rhs){
	return Cross(lhs,rhs);
}

inline std::ostream &operator<<(std::ostream &output, const CVector3D& v)
{
  output.setf(std::ios::scientific);
  output << v.x << " " << v.y << " " << v.z;
  return output;
  }
  
inline bool operator == (const CVector3D& lhs, const CVector3D& rhs){
  if( fabs(lhs.x - rhs.x) < NEARLY_ZERO
     && fabs(lhs.y - rhs.y) < NEARLY_ZERO
     && fabs(lhs.z - rhs.z) < NEARLY_ZERO )
    return true;
  else return false;
}
  
inline bool operator != (const CVector3D& lhs, const CVector3D& rhs){
  if( lhs == rhs )	return false;
  else return true;
}
  
  
  
//! length of vector
inline double Length(const CVector3D& point)
{
  return	sqrt( point.x*point.x + point.y*point.y + point.z*point.z );
}
  
inline double SquareDistance(const CVector3D& ipo0, const CVector3D& ipo1)
{
  return
    ( ipo1.x - ipo0.x )*( ipo1.x - ipo0.x )
  + ( ipo1.y - ipo0.y )*( ipo1.y - ipo0.y )
  + ( ipo1.z - ipo0.z )*( ipo1.z - ipo0.z );
 }
  
//! distance between two points
inline double Distance(const CVector3D& ipo0, const CVector3D& ipo1)
{
  return	sqrt( SquareDistance(ipo0,ipo1) );
}  
  
inline double Dot(const CVector3D &arg1, const CVector3D &arg2)
{
  return arg1.x*arg2.x + arg1.y*arg2.y + arg1.z*arg2.z;
}
 
inline CVector3D Cross(const CVector3D& arg1, const CVector3D& arg2)
{
  CVector3D temp;
  temp.x = arg1.y*arg2.z - arg1.z*arg2.y;
  temp.y = arg1.z*arg2.x - arg1.x*arg2.z;
  temp.z = arg1.x*arg2.y - arg1.y*arg2.x;
  return temp;
}
  
inline double ScalarTripleProduct(const CVector3D& a, const CVector3D& b, const CVector3D& c){
  return a.x*(b.y*c.z - b.z*c.y) + a.y*(b.z*c.x - b.x*c.z) + a.z*(b.x*c.y - b.y*c.x);
}
  
  
  
//! Volume of a tetrahedra
inline double TetVolume(const CVector3D& v1,
                        const CVector3D& v2,
                        const CVector3D& v3,
                        const CVector3D& v4 )
{
  return
  (   ( v2.x - v1.x )*( ( v3.y - v1.y )*( v4.z - v1.z ) - ( v4.y - v1.y )*( v3.z - v1.z ) )
   -	( v2.y - v1.y )*( ( v3.x - v1.x )*( v4.z - v1.z ) - ( v4.x - v1.x )*( v3.z - v1.z ) )
   +	( v2.z - v1.z )*( ( v3.x - v1.x )*( v4.y - v1.y ) - ( v4.x - v1.x )*( v3.y - v1.y ) )
   ) * 0.16666666666666666666666666666667;
}
  
  
  
inline void UnitNormal
(CVector3D& vnorm,
 const CVector3D& v1,
 const CVector3D& v2,
 const CVector3D& v3)
{
  vnorm.x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
  vnorm.y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
  vnorm.z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);
  const double dtmp1 = 1.0 / Length(vnorm);
  vnorm.x *= dtmp1;
  vnorm.y *= dtmp1;
  vnorm.z *= dtmp1;
}
  
  
//! Hight of a tetrahedra
inline double Height(const CVector3D& v1, const CVector3D& v2, const CVector3D& v3, const CVector3D& v4){
  // get normal vector
  double dtmp_x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
  double dtmp_y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
  double dtmp_z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);
    
  // normalize normal vector
  const double dtmp1 = 1.0 / sqrt( dtmp_x*dtmp_x + dtmp_y*dtmp_y + dtmp_z*dtmp_z );
  dtmp_x *= dtmp1;
  dtmp_y *= dtmp1;
  dtmp_z *= dtmp1;
    
  return (v4.x-v1.x)*dtmp_x+(v4.y-v1.y)*dtmp_y+(v4.z-v1.z)*dtmp_z;
}

  
inline CVector3D GetMinDist_LinePoint(const CVector3D& p, // point
                                      const CVector3D& s, // source
                                      const CVector3D& d) // direction
{
  assert( Dot(d,d) > 1.0e-20 );
  const CVector3D ps = s-p;
  double a = Dot(d,d);
  double b = Dot(d,s-p);
  double t = -b/a;
  return s+t*d;
}
  
inline CVector3D GetMinDist_LineSegPoint(const CVector3D& p, // point
                                         const CVector3D& s, // source
                                         const CVector3D& e) // direction
{
  CVector3D d = e-s;
  assert( Dot(d,d) > 1.0e-20 );
  const CVector3D ps = s-p;
  double a = Dot(d,d);
  double b = Dot(d,s-p);
  double t = -b/a;
  if( t < 0 ) t = 0;
  if( t > 1 ) t = 1;
  return s+t*d;
}

#endif
