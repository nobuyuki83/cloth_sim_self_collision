//
//  aabb.h
//  self_contact
//
//  Created by Nobuyuki Umetani on 11/8/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//

#ifndef aabb_h
#define aabb_h

#include <math.h>
#include <assert.h>

//! 3D bounding box class
class CAABB3D
{
public:
	CAABB3D(){
		x_min=0;	x_max=0;
		y_min=0;	y_max=0;
		z_min=0;	z_max=0;
		isnt_empty = false;
	}
	CAABB3D(double x_min0,double x_max0,  double y_min0,double y_max0,  double z_min0,double z_max0)
  : x_min(x_min0),x_max(x_max0),
  y_min(y_min0),y_max(y_max0),
  z_min(z_min0),z_max(z_max0)
	{
		assert( x_min <= x_max );
		assert( y_min <= y_max );
		assert( z_min <= z_max );
		isnt_empty = true;
	}
	CAABB3D( const CAABB3D& bb )
  : x_min(bb.x_min),x_max(bb.x_max),
  y_min(bb.y_min),y_max(bb.y_max),
  z_min(bb.z_min),z_max(bb.z_max),
  isnt_empty(bb.isnt_empty){}
  void SetCenterWidth(double cx, double cy, double cz,  double wx, double wy, double wz)
  {
    x_min = cx-wx*0.5; x_max = cx+wx*0.5;
    y_min = cy-wy*0.5; y_max = cy+wy*0.5;
    z_min = cz-wz*0.5; z_max = cz+wz*0.5;
  }
  void GetCenterWidth(double& cx, double& cy, double& cz,  double& wx, double& wy, double& wz)
  {
    cx = (x_max+x_min)*0.5;
    cy = (y_max+y_min)*0.5;
    cz = (z_max+z_min)*0.5;
    wx = (x_max-x_min);
    wy = (y_max-y_min);
    wz = (z_max-z_min);
  }
	CAABB3D& operator+=(const CAABB3D& bb)
	{
		if( !bb.isnt_empty ) return *this;
		if( !isnt_empty ){
			x_max = bb.x_max;	x_min = bb.x_min;
			y_max = bb.y_max;	y_min = bb.y_min;
			z_max = bb.z_max;	z_min = bb.z_min;
      this->isnt_empty = bb.isnt_empty;
			return *this;
		}
		x_max = ( x_max > bb.x_max ) ? x_max : bb.x_max;
		x_min = ( x_min < bb.x_min ) ? x_min : bb.x_min;
		y_max = ( y_max > bb.y_max ) ? y_max : bb.y_max;
		y_min = ( y_min < bb.y_min ) ? y_min : bb.y_min;
		z_max = ( z_max > bb.z_max ) ? z_max : bb.z_max;
		z_min = ( z_min < bb.z_min ) ? z_min : bb.z_min;
		return *this;
	}
  bool IsIntersect(const CAABB3D& bb) const
  {
    if( !isnt_empty ) return false;
    if( !bb.isnt_empty ) return false;
    if( x_max < bb.x_min ) return false;
    if( x_min > bb.x_max ) return false;
    if( y_max < bb.y_min ) return false;
    if( y_min > bb.y_max ) return false;
    if( z_max < bb.z_min ) return false;
    if( z_min > bb.z_max ) return false;
    return true;
  }
  CAABB3D& operator+=(const double v[3])
	{
		if( !isnt_empty ){
			x_max = v[0];	x_min = v[0];
			y_max = v[1];	y_min = v[1];
			z_max = v[2];	z_min = v[2];
      this->isnt_empty = true;
			return *this;
		}
		x_max = ( x_max > v[0] ) ? x_max : v[0];
		x_min = ( x_min < v[0] ) ? x_min : v[0];
		y_max = ( y_max > v[1] ) ? y_max : v[1];
		y_min = ( y_min < v[1] ) ? y_min : v[1];
		z_max = ( z_max > v[2] ) ? z_max : v[2];
		z_min = ( z_min < v[2] ) ? z_min : v[2];
		return *this;
	}
  bool AddPoint(double x, double y, double z, double eps){
    if( eps <= 0 ){ return false; }
    if( isnt_empty ){
      x_min = ( x_min < x-eps ) ? x_min : x-eps;
      y_min = ( y_min < y-eps ) ? y_min : y-eps;
      z_min = ( z_min < z-eps ) ? z_min : z-eps;
      x_max = ( x_max > x+eps ) ? x_max : x+eps;
      y_max = ( y_max > y+eps ) ? y_max : y+eps;
      z_max = ( z_max > z+eps ) ? z_max : z+eps;
    }
    else{
      isnt_empty = true;
      x_min = x-eps;
      y_min = y-eps;
      z_min = z-eps;
      x_max = x+eps;
      y_max = y+eps;
      z_max = z+eps;
    }
    return true;
  }
  double MinimumDistance(double x, double y, double z) const
  {
    double x0, y0, z0;
    if(      x < x_min ){ x0 = x_min; }
    else if( x < x_max ){ x0 = x;     }
    else{                 x0 = x_max; }
    if(      y < y_min ){ y0 = y_min; }
    else if( y < y_max ){ y0 = y;     }
    else{                 y0 = y_max; }
    if(      z < z_min ){ z0 = z_min; }
    else if( z < z_max ){ z0 = z;     }
    else{                 z0 = z_max; }
    return sqrt( (x0-x)*(x0-x) + (y0-y)*(y0-y) + (z0-z)*(z0-z) );
  }
  bool IsInside(double x, double y, double z) const
  {
   if( !isnt_empty ) return false; // âΩÇ‡Ç»Ç¢èÍçáÇÕèÌÇ…ãU
   if(   x >= x_min && x <= x_max
      && y >= y_min && y <= y_max
      && z >= z_min && z <= z_max ) return true;
   return false;
  }
public:
	double x_min,x_max,  y_min,y_max,  z_min,z_max;
	bool isnt_empty;	//!< false if there is nothing inside
};


#endif
