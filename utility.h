//
//  utility.h
//
//  Created by Nobuyuki Umetani on 11/7/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//

#ifndef utility_h
#define utility_h

#include <math.h>

// inner product of 3-vector a[3] and b[3]
// 内積を求める
inline double Dot3D(const double a[], const double b[]){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

// cross product of 3-vector v1[3] and v2[3]
// 外積を求める
inline void Cross3D(double r[3], const double v1[3], const double v2[3]){
  r[0] = v1[1]*v2[2] - v2[1]*v1[2];
  r[1] = v1[2]*v2[0] - v2[2]*v1[0];
  r[2] = v1[0]*v2[1] - v2[0]*v1[1];
}

// length of vector v[3]
// ベクトルの長さを求める
inline double Length3D(const double v[3]){
  return sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

// distance betweenr v1[3] and v2[3]
// ２つのベクトルの距離を求める
inline double Distance3D(const double v1[3], const double v2[3])
{
  const double t = (v2[0]-v1[0])*(v2[0]-v1[0]) + (v2[1]-v1[1])*(v2[1]-v1[1]) + (v2[2]-v1[2])*(v2[2]-v1[2]);
  return sqrt(t);
}

// area of a triangle
// 三角形の面積を求める
inline double TriArea3D
(const double v1[3], // position of a 1st vertex
 const double v2[3], // position of a 2nd vertex
 const double v3[3]) // position of a 3rd vertex
{
  double n[3];
  n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  return sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
}


// compute unitnomal and area of a triangle
// 三角形の単位法線ベクトルを求める
inline void  UnitNormalAreaTri3D
(double n[3], // unit normal
 double& a, // area
 ////
 const double v1[3], // position of a 1st vertex
 const double v2[3], // position of a 2nd vertex
 const double v3[3]  // position of a 3rd vertex
 )
{
  n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
  const double invlen = 0.5/a;
  n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}

/* ---------------------------------------------------------------------- */


// multiple quaternions (for camera control)
// 四元数の掛け算
inline void QuatMult(double r[], const double p[], const double q[])
{
  r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
  r[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
  r[2] = p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1];
  r[3] = p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0];
}

// copy quaternions (for camera control)
// 四元数のコピー
inline void QuatCopy(double r[], const double p[])
{
  r[0] = p[0];
  r[1] = p[1];
  r[2] = p[2];
  r[3] = p[3];
}

// rotate quaternion (for camera control)
// 四元数の回転
inline void QuatRot(double r[], const double q[])
{
  double x2 = q[1] * q[1] * 2.0;
  double y2 = q[2] * q[2] * 2.0;
  double z2 = q[3] * q[3] * 2.0;
  double xy = q[1] * q[2] * 2.0;
  double yz = q[2] * q[3] * 2.0;
  double zx = q[3] * q[1] * 2.0;
  double xw = q[1] * q[0] * 2.0;
  double yw = q[2] * q[0] * 2.0;
  double zw = q[3] * q[0] * 2.0;
  
  r[ 0] = 1.0 - y2 - z2;
  r[ 1] = xy + zw;
  r[ 2] = zx - yw;
  r[ 4] = xy - zw;
  r[ 5] = 1.0 - z2 - x2;
  r[ 6] = yz + xw;
  r[ 8] = zx + yw;
  r[ 9] = yz - xw;
  r[10] = 1.0 - x2 - y2;
  r[ 3] = r[ 7] = r[11] = r[12] = r[13] = r[14] = 0.0;
  r[15] = 1.0;
}


#endif
