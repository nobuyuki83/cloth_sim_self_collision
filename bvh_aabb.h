//
//  bvh_aabb.h
//  self_contact
//
//  Created by Nobuyuki Umetani on 11/8/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//

#ifndef bvh_aabb_h
#define bvh_aabb_h

#include <stack>
#include <vector>

#include "aabb.h"
#include "vector3d.h"

class CNodeBVH
{
public:
  int iroot;
  int ichild[2];
};

int MakeBVHTopology_TopDown
(const std::vector<int>& aTri,
 const std::vector<double>& aXYZ,
 std::vector<CNodeBVH>& aNodeBVH);

// BVHのBounding Boxを構築
void BuildBoundingBoxChild_Prx
(int ibvh,
 double delta,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 const std::vector<CNodeBVH>& aNodeBVH,
 std::vector<CAABB3D>& aBB);

void BuildBoundingBoxChild_CCD
(int ibvh,
 double dt,
 const std::vector<double>& aXYZ,
 const std::vector<double>& aUVW,
 const std::vector<int>& aTri,
 const std::vector<CNodeBVH>& aNodeBVH,
 std::vector<CAABB3D>& aBB);


double DistanceFaceVertex
(const CVector3D& p0, const CVector3D& p1, const CVector3D& p2,
 const CVector3D& p3,
 double& w0, double& w1);

double DistanceEdgeEdge
(const CVector3D& p0, const CVector3D& p1,
 const CVector3D& q0, const CVector3D& q1,
 double& ratio_p, double& ratio_q);

double FindCoplanerInterp
(const CVector3D& s0, const CVector3D& s1, const CVector3D& s2, const CVector3D& s3,
 const CVector3D& e0, const CVector3D& e1, const CVector3D& e2, const CVector3D& e3);

// smallest element of contact (vertex-face or edge-edge)
// 接触要素のクラス（点と面，辺と辺）
class CContactElement
{
public:
  CContactElement(bool is_fv,int j0, int j1, int j2, int j3)
  {
    this->is_fv = is_fv;
    if( is_fv ){
      this->is_fv = true;
      // 検索が高速になるように並び替える
      if(      j0 < j1 && j0 < j2 && j1 < j2 ){ ino0=j0;  ino1=j1;  ino2=j2;  ino3=j3; }
      else if( j0 < j1 && j0 < j2 && j2 < j1 ){ ino0=j0;  ino1=j2;  ino2=j1;  ino3=j3; }
      else if( j1 < j0 && j1 < j2 && j0 < j2 ){ ino0=j1;  ino1=j0;  ino2=j2;  ino3=j3; }
      else if( j1 < j0 && j1 < j2 && j2 < j0 ){ ino0=j1;  ino1=j2;  ino2=j0;  ino3=j3; }
      else if( j2 < j0 && j2 < j1 && j0 < j1 ){ ino0=j2;  ino1=j0;  ino2=j1;  ino3=j3; }
      else if( j2 < j0 && j2 < j1 && j1 < j0 ){ ino0=j2;  ino1=j1;  ino2=j0;  ino3=j3; }
      else { assert(0); }
    }
    else{
      this->is_fv = false;
      // 検索が高速になるように並び替える      
      if(      j0 < j1 && j0 < j2 && j0 < j3 && j2 < j3 ){ ino0=j0;  ino1=j1;  ino2=j2;  ino3=j3; }
      else if( j0 < j1 && j0 < j2 && j0 < j3 && j3 < j2 ){ ino0=j0;  ino1=j1;  ino2=j3;  ino3=j2; }
      else if( j1 < j0 && j1 < j2 && j1 < j3 && j2 < j3 ){ ino0=j1;  ino1=j0;  ino2=j2;  ino3=j3; }
      else if( j1 < j0 && j1 < j2 && j1 < j3 && j3 < j2 ){ ino0=j1;  ino1=j0;  ino2=j3;  ino3=j2; }
      else if( j2 < j0 && j2 < j1 && j2 < j3 && j0 < j1 ){ ino0=j2;  ino1=j3;  ino2=j0;  ino3=j1; }
      else if( j2 < j0 && j2 < j1 && j2 < j3 && j1 < j0 ){ ino0=j2;  ino1=j3;  ino2=j1;  ino3=j0; }
      else if( j3 < j0 && j3 < j1 && j3 < j2 && j0 < j1 ){ ino0=j3;  ino1=j2;  ino2=j0;  ino3=j1; }
      else if( j3 < j0 && j3 < j1 && j3 < j2 && j1 < j0 ){ ino0=j3;  ino1=j2;  ino2=j1;  ino3=j0; }
      else { assert(0); }
    }
  }
  bool operator < (const CContactElement& p2) const
  {
    if( ino0 != p2.ino0 ){ return ino0 < p2.ino0; }
    if( ino1 != p2.ino1 ){ return ino1 < p2.ino1; }
    if( ino2 != p2.ino2 ){ return ino2 < p2.ino2; }
    return ino3 < p2.ino3;
  }
  public:
    bool is_fv; // true: ee contact, false: vf contact, 真ならFV，偽ならEE
    int ino0, ino1, ino2, ino3; // four points in the contact, 接触要素のの中の４点
};



void GetContactElement_Proximity
(std::set<CContactElement>& aContactElem,
 ////
 double delta,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 int ibvh,
 const std::vector<CNodeBVH>& aBVH,
 const std::vector<CAABB3D>& aBB);
      
void GetContactElement_CCD
(std::set<CContactElement>& aContactElem,
 /////
 double dt,
 double delta,
 const std::vector<double>& aXYZ,
 const std::vector<double>& aUVW,
 const std::vector<int>& aTri,
 int ibvh,
 const std::vector<CNodeBVH>& aBVH,
 const std::vector<CAABB3D>& aBB);

#endif
