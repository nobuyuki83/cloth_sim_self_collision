//
//  bvh_aabb.cpp
//
//  Created by Nobuyuki Umetani on 11/10/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//

#include <stdio.h>
#include <set>

#include "bvh_aabb.h"
#include "jagged_array.h"
#include "vector3d.h"


static inline void AddPoint
(CAABB3D& bb, const CVector3D& p, double eps)
{
  bb.AddPoint(p.x, p.y, p.z, eps);
}

// BVHのBounding Boxを構築,三角形同士の近接計算用
void BuildBoundingBoxChild_Prx
(int ibvh,
 double delta,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 const std::vector<CNodeBVH>& aNodeBVH,
 std::vector<CAABB3D>& aBB)
{
  aBB.resize( aNodeBVH.size() );
  assert( ibvh < aNodeBVH.size() );
  int ichild0 = aNodeBVH[ibvh].ichild[0];
  int ichild1 = aNodeBVH[ibvh].ichild[1];
  if( ichild1 == -1 ){ // leaf、葉ノード
    int itri = ichild0;
    assert( itri < aTri.size() );
    int ino0 = aTri[itri*3+0];
    int ino1 = aTri[itri*3+1];
    int ino2 = aTri[itri*3+2];
    const CVector3D p0( aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2] );
    const CVector3D p1( aXYZ[ino1*3+0], aXYZ[ino1*3+1], aXYZ[ino1*3+2] );
    const CVector3D p2( aXYZ[ino2*3+0], aXYZ[ino2*3+1], aXYZ[ino2*3+2] );
    CAABB3D& bb = aBB[ibvh];
    bb.isnt_empty = false;
    bb.AddPoint(p0.x,p0.y,p0.z, delta*0.5);
    bb.AddPoint(p1.x,p1.y,p1.z, delta*0.5);
    bb.AddPoint(p2.x,p2.y,p2.z, delta*0.5);
    return;
  }
  // internal node,内部ノードは子ノードのBounding Volume
  assert( aNodeBVH[ichild0].iroot == ibvh );
  assert( aNodeBVH[ichild1].iroot == ibvh );
  BuildBoundingBoxChild_Prx(ichild0,delta, aXYZ,aTri,aNodeBVH,aBB);
  BuildBoundingBoxChild_Prx(ichild1,delta, aXYZ,aTri,aNodeBVH,aBB);
  CAABB3D& bb = aBB[ibvh];
  bb.isnt_empty = false;
  bb  = aBB[ichild0];
  bb += aBB[ichild1];
  return;
}

// BVHのBounding Boxを構築,三角形同士のCCD交差計算用
void BuildBoundingBoxChild_CCD
(int ibvh,
 double dt,
 const std::vector<double>& aXYZ,
 const std::vector<double>& aUVW,
 const std::vector<int>& aTri,
 const std::vector<CNodeBVH>& aNodeBVH,
 std::vector<CAABB3D>& aBB)
{
  double eps = 1.0e-10;
  //  std::cout << ibvh << " " << aNodeBVH_TreeTopo.size() << std::endl;
  assert( ibvh < aNodeBVH.size() );
  int ichild0 = aNodeBVH[ibvh].ichild[0];
  int ichild1 = aNodeBVH[ibvh].ichild[1];
  if( ichild1 == -1 ){ // leaf
    const int itri = ichild0;
    assert( itri < aTri.size() );
    const int ino0 = aTri[itri*3+0];
    const int ino1 = aTri[itri*3+1];
    const int ino2 = aTri[itri*3+2];
    const CVector3D p0( aXYZ[ino0*3+0],                   aXYZ[ino0*3+1],                   aXYZ[ino0*3+2] );
    const CVector3D p1( aXYZ[ino1*3+0],                   aXYZ[ino1*3+1],                   aXYZ[ino1*3+2] );
    const CVector3D p2( aXYZ[ino2*3+0],                   aXYZ[ino2*3+1],                   aXYZ[ino2*3+2] );
    const CVector3D q0( aXYZ[ino0*3+0]+dt*aUVW[ino0*3+0], aXYZ[ino0*3+1]+dt*aUVW[ino0*3+1], aXYZ[ino0*3+2]+dt*aUVW[ino0*3+2] );
    const CVector3D q1( aXYZ[ino1*3+0]+dt*aUVW[ino1*3+0], aXYZ[ino1*3+1]+dt*aUVW[ino1*3+1], aXYZ[ino1*3+2]+dt*aUVW[ino1*3+2] );
    const CVector3D q2( aXYZ[ino2*3+0]+dt*aUVW[ino2*3+0], aXYZ[ino2*3+1]+dt*aUVW[ino2*3+1], aXYZ[ino2*3+2]+dt*aUVW[ino2*3+2] );
    CAABB3D& bb = aBB[ibvh];
    bb.isnt_empty = false;
    AddPoint(bb,p0, eps);  AddPoint(bb,p1, eps);  AddPoint(bb,p2, eps);
    AddPoint(bb,q0, eps);  AddPoint(bb,q1, eps);  AddPoint(bb,q2, eps);
    return;
  }
  // internal node,内部ノードは子ノードのBounding Volume  
  assert( aNodeBVH[ichild0].iroot == ibvh );
  assert( aNodeBVH[ichild1].iroot == ibvh );
  BuildBoundingBoxChild_CCD(ichild0,dt, aXYZ,aUVW,aTri,aNodeBVH,aBB);
  BuildBoundingBoxChild_CCD(ichild1,dt, aXYZ,aUVW,aTri,aNodeBVH,aBB);
  CAABB3D& bb = aBB[ibvh];
  bb.isnt_empty = false;
  bb  = aBB[ichild0];
  bb += aBB[ichild1];
  return;
}

// 三角形を囲む三角形を構築する
static void MakeTriSurTri
(std::vector<int>& aTriSur,
 const int nnode,
 const std::vector<int>& aTri)
{
  const int nTri = (int)aTri.size()/3;
  const int e2n[3][2] = {{1,2},{2,0},{0,1}};
  aTriSur.clear();
  aTriSur.resize(nTri*3,-1);
  CJaggedArray aTSV;  aTSV.SetNodeToElem(aTri, nTri, 3, nnode);
  for(int itri=0;itri<nTri;itri++){
    for(int iedge=0;iedge<3;iedge++){
      int ino0 = aTri[itri*3+e2n[iedge][0]];
      int ino1 = aTri[itri*3+e2n[iedge][1]];
      for(int itrisur=aTSV.index[ino0];itrisur<aTSV.index[ino0+1];itrisur++){
        int jtri = aTSV.array[itrisur];
        if( jtri == itri ) continue;
        for(int jedge=0;jedge<3;jedge++){
          int jno0 = aTri[jtri*3+e2n[jedge][0]];
          int jno1 = aTri[jtri*3+e2n[jedge][1]];
          if( ino0 == jno1 && ino1 == jno0 ){
            aTriSur[itri*3+iedge] = jtri;
            break;
          }
        }
        if( aTriSur[itri*3+iedge] == jtri ) break;
      }
    }
  }
}


// the center of gravity of a triangle
CVector3D TriGravityCenter
(int itri,
 const std::vector<int>& aTri,
 const std::vector<double>& aXYZ)
{
  int ino0 = aTri[itri*3+0];
  int ino1 = aTri[itri*3+1];
  int ino2 = aTri[itri*3+2];
  double cgx = (aXYZ[ino0*3+0] + aXYZ[ino1*3+0] + aXYZ[ino2*3+0])/3.0;
  double cgy = (aXYZ[ino0*3+1] + aXYZ[ino1*3+1] + aXYZ[ino2*3+1])/3.0;
  double cgz = (aXYZ[ino0*3+2] + aXYZ[ino1*3+2] + aXYZ[ino2*3+2])/3.0;
  return CVector3D(cgx,cgy,cgz);
}

void DevideTriAryConnex
(int iroot_node,
 std::vector<int>& aTri2Node,
 std::vector<CNodeBVH>& aNodeBVH,
 ////
 const std::vector<int>& list,
 const std::vector<int>& aTri,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTriSur)
{
  assert( list.size() > 1 );
  double eps = 1.0e-10;
// リストに含まれる三角形の重心のAABBを作る  
  CAABB3D bb;
  for(int il=0;il<list.size();il++){
    int itri = list[il];
    assert( itri < aTri.size()/3 );
    assert( aTri2Node[itri] == iroot_node );
    CVector3D cg = TriGravityCenter(itri,aTri,aXYZ);
    AddPoint(bb,cg, eps);
  }
  double lenx = bb.x_max - bb.x_min;
  double leny = bb.y_max - bb.y_min;
  double lenz = bb.z_max - bb.z_min;
  CVector3D dir; // AABBの最も長い辺の向き
  if( lenx > leny && lenx > lenz ){ dir = CVector3D(1,0,0); }
  if( leny > lenz && leny > lenx ){ dir = CVector3D(0,1,0); }
  if( lenz > lenx && lenz > leny ){ dir = CVector3D(0,0,1); }
  CVector3D org((bb.x_min+bb.x_max)*0.5,  (bb.y_min+bb.y_max)*0.5,  (bb.z_min+bb.z_max)*0.5);
  int itri_ker = -1;
  for(int il=0;il<list.size();il++){
    int itri0 = list[il];
    const CVector3D& gc0 = TriGravityCenter(itri0,aTri,aXYZ);
    if( fabs(Dot(gc0-org,dir)) < 1.0e-10 ) continue;
    if( Dot(gc0-org,dir) < 0 ){ dir *= -1; }
    itri_ker = itri0;
    break;
  }
  if( itri_ker == -1 ){
    org = CVector3D(0,0,0);
    for(int il=0;il<list.size();il++){
      int itri0 = list[il];
      const CVector3D& gc0 = TriGravityCenter(itri0,aTri,aXYZ);
      org += gc0;
    }
    org = org/list.size();
    double mat[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
    for(int il=0;il<list.size();il++){
      int itri0 = list[il];
      const CVector3D& gc0 = TriGravityCenter(itri0,aTri,aXYZ);
      CVector3D vcg = gc0-org;
      for(int i=0;i<3;i++){ for(int j=0;j<3;j++){ mat[i][j] += vcg[i]*vcg[j]; } }
    }
    dir = CVector3D(1,1,1);
    for(int i=0;i<10;i++){
      double tmp[3] = {
        mat[0][0]*dir[0] + mat[0][1]*dir[1] + mat[0][2]*dir[2],
        mat[1][0]*dir[0] + mat[1][1]*dir[1] + mat[1][2]*dir[2],
        mat[2][0]*dir[0] + mat[2][1]*dir[1] + mat[2][2]*dir[2] };
      double len = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
      dir[0] = tmp[0]/len;
      dir[1] = tmp[1]/len;
      dir[2] = tmp[2]/len;
    }
    for(int il=0;il<list.size();il++){
      int itri0 = list[il];
      const CVector3D& gc0 = TriGravityCenter(itri0,aTri,aXYZ);
      if( fabs(Dot(gc0-org,dir)) < 1.0e-10 ) continue;
      if( Dot(gc0-org,dir) < 0 ){ dir *= -1; }
      itri_ker = itri0;
      break;
    }
  }
  int inode_ch0 = (int)aNodeBVH.size();
  int inode_ch1 = (int)aNodeBVH.size()+1;
  aNodeBVH.resize(aNodeBVH.size()+2);
  aNodeBVH[inode_ch0].iroot = iroot_node;
  aNodeBVH[inode_ch1].iroot = iroot_node;
  aNodeBVH[iroot_node].ichild[0] = inode_ch0;
  aNodeBVH[iroot_node].ichild[1] = inode_ch1;
  std::vector<int> list_ch0;
  { // 子ノード０に含まれる三角形を抽出（itri_kerと接続していて，dirベクトルの正方向にある三角形）
    aTri2Node[itri_ker] = inode_ch0;
    list_ch0.push_back(itri_ker);
    std::stack<int> stack;
    stack.push(itri_ker);
    while(!stack.empty()){
      int itri0 = stack.top();
      stack.pop();
      for(int iedge=0;iedge<3;iedge++){
        int jtri = aTriSur[itri0*3+iedge];
        if( jtri == -1 ) continue;
        if( aTri2Node[jtri] != iroot_node ) continue;
        assert( jtri < aTri.size()/3 );
        const CVector3D& gc1 = TriGravityCenter(jtri,aTri,aXYZ);
        if( Dot(gc1-org,dir) < 0 ) continue;
        stack.push(jtri);
        aTri2Node[jtri] = inode_ch0;
        list_ch0.push_back(jtri);
      }
    }
    assert( list_ch0.size() > 0 );
  }
  // 子ノード１に含まれる三角形を抽出(入力リストから子ノード0に含まれる物を除外)
  std::vector<int> list_ch1;
  for(int il=0;il<list.size();il++){
    int itri = list[il];
    if( aTri2Node[itri] == inode_ch0 ) continue;
    assert( aTri2Node[itri] == iroot_node );
    aTri2Node[itri] = inode_ch1;
    list_ch1.push_back(itri);
  }
  assert( list_ch1.size() > 0 );
  
  /////
  if( list_ch0.size() == 1 ){
    aNodeBVH[inode_ch0].ichild[0] = list_ch0[0];
    aNodeBVH[inode_ch0].ichild[1] = -1;
  }
  else{ // 子ノード0にある三角形を再度分割
    DevideTriAryConnex(inode_ch0,aTri2Node,aNodeBVH,
                       list_ch0,aTri,aXYZ,aTriSur);
  }
  list_ch0.clear();
  
  //////
  if( list_ch1.size() == 1 ){
    aNodeBVH[inode_ch1].ichild[0] = list_ch1[0];
    aNodeBVH[inode_ch1].ichild[1] = -1;
  }
  else{ // 子ノード1にある三角形を再度分割
    DevideTriAryConnex(inode_ch1,aTri2Node,aNodeBVH,
                       list_ch1,aTri,aXYZ,aTriSur);
  }
}



int MakeBVHTopology_TopDown
(const std::vector<int>& aTri,
 const std::vector<double>& aXYZ,
 std::vector<CNodeBVH>& aNodeBVH)
{
  std::vector<int> aTriSur;
  MakeTriSurTri(aTriSur,(int)aXYZ.size()/3,aTri);
  
  aNodeBVH.clear();
  const int ntri = (int)aTri.size()/3;
  std::vector<int> list(ntri);
  for(int itri=0;itri<ntri;itri++){ list[itri] = itri; }
  std::vector<int> aTri2Node;
  aTri2Node.resize(ntri,0);
  aNodeBVH.resize(1);
  aNodeBVH[0].iroot = -1;
  DevideTriAryConnex(0,aTri2Node,aNodeBVH,
                     list,aTri,aXYZ,aTriSur);
  return 0;
}



/* ---------------------------------------------------------------------------------- */


// VFの距離
double DistanceFaceVertex
(const CVector3D& p0, const CVector3D& p1, const CVector3D& p2,
 const CVector3D& p3,
 double& w0, double& w1)
{
  CVector3D v20 =p0-p2;
  CVector3D v21 =p1-p2;
  double t0 = Dot(v20,v20);
  double t1 = Dot(v21,v21);
  double t2 = Dot(v20,v21);
  double t3 = Dot(v20,p3-p2);
  double t4 = Dot(v21,p3-p2);
  double det = t0*t1-t2*t2;
  double invdet = 1.0/det;
  w0 = (+t1*t3-t2*t4)*invdet;
  w1 = (-t2*t3+t0*t4)*invdet;
  const double w2 = 1-w0-w1;
  CVector3D pw = w0*p0 + w1*p1 + w2*p2;
  return (pw-p3).Length();
}

//　EEの距離
double DistanceEdgeEdge
(const CVector3D& p0, const CVector3D& p1,
 const CVector3D& q0, const CVector3D& q1,
 double& ratio_p, double& ratio_q)
{
  const CVector3D& vp =p1-p0;
  const CVector3D& vq =q1-q0;
  if( Cross(vp,vq).Length() < 1.0e-10 ){ // handling parallel edge
    CVector3D pq0 = p0-q0;
    CVector3D nvp = vp; nvp.SetNormalizedVector();
    CVector3D vert = pq0 - Dot(pq0,nvp)*nvp;
    double dist = vert.Length();
    double lp0 = Dot(p0,nvp);
    double lp1 = Dot(p1,nvp);
    double lq0 = Dot(q0,nvp);
    double lq1 = Dot(q1,nvp);
    double p_min  = ( lp0 < lp1 ) ? lp0 : lp1;
    double p_max  = ( lp0 > lp1 ) ? lp0 : lp1;
    double q_min  = ( lq0 < lq1 ) ? lq0 : lq1;
    double q_max  = ( lq0 > lq1 ) ? lq0 : lq1;
    double lm;
    if(      p_max < q_min ){ lm = (p_max+q_min)*0.5; }
    else if( q_max < p_min ){ lm = (q_max+p_min)*0.5; }
    else if( p_max < q_max ){ lm = (p_max+q_min)*0.5; }
    else{                     lm = (q_max+p_min)*0.5; }
    ratio_p = (lm-lp0)/(lp1-lp0);
    ratio_q = (lm-lq0)/(lq1-lq0);
    return dist;
  }
  double t0 = Dot(vp,vp);
  double t1 = Dot(vq,vq);
  double t2 = Dot(vp,vq);
  double t3 = Dot(vp,q0-p0);
  double t4 = Dot(vq,q0-p0);
  double det = t0*t1-t2*t2;
  double invdet = 1.0/det;
  ratio_p = (+t1*t3-t2*t4)*invdet;
  ratio_q = (+t2*t3-t0*t4)*invdet;
  CVector3D pc = p0 + ratio_p*vp;
  CVector3D qc = q0 + ratio_q*vq;
  return (pc-qc).Length();
}



/////////////////////////////////////////////////////////////////////


// VFの距離が所定の距離以下にあるかどうか
bool IsContact_FV_Proximity
(int ino0, int ino1, int ino2, int ino3,
 const CVector3D& p0, const CVector3D& p1, const CVector3D& p2, const CVector3D& p3,
 const CAABB3D& bb,
 const double delta)
{
  if( ino3 == ino0 || ino3 == ino1 || ino3 == ino2 ){ return false; }
  if( !bb.IsInside(p3.x,p3.y,p3.z) ) return false;
  double height = fabs( Height(p0,p1,p2,p3) );
  if( height > delta ) return false;
  double w0,w1;
  const double dist = DistanceFaceVertex(p0,p1,p2,p3,w0,w1);
  const double w2 = 1-w0-w1;
  if( dist > delta ) return false;
  double mgn = ( Distance(p0, p1) + Distance(p1, p2) + Distance(p2, p3) ) / 3.0;
  mgn = 0;
  if( w0 < -mgn || w0 > 1+mgn ) return false;
  if( w1 < -mgn || w1 > 1+mgn ) return false;
  if( w2 < -mgn || w2 > 1+mgn ) return false;
  return true;
}

// EEの距離が所定の距離以下にあるかどうか
bool IsContact_EE_Proximity
(int ino0,        int ino1,        int jno0,        int jno1,
 const CVector3D& p0, const CVector3D& p1, const CVector3D& q0, const CVector3D& q1,
 const double delta)
{
  if( ino0 == jno0 || ino0 == jno1 || ino1 == jno0 || ino1 == jno1 ) return false;
  if( q0.x+delta < p0.x && q0.x+delta < p1.x && q1.x+delta < p0.x && q1.x+delta < p1.x ) return false;
  if( q0.x-delta > p0.x && q0.x-delta > p1.x && q1.x-delta > p0.x && q1.x-delta > p1.x ) return false;
  if( q0.y+delta < p0.y && q0.y+delta < p1.y && q1.y+delta < p0.y && q1.y+delta < p1.y ) return false;
  if( q0.y-delta > p0.y && q0.y-delta > p1.y && q1.y-delta > p0.y && q1.y-delta > p1.y ) return false;
  if( q0.z+delta < p0.z && q0.z+delta < p1.z && q1.z+delta < p0.z && q1.z+delta < p1.z ) return false;
  if( q0.z-delta > p0.z && q0.z-delta > p1.z && q1.z-delta > p0.z && q1.z-delta > p1.z ) return false;
  double ratio_p, ratio_q;
  double dist = DistanceEdgeEdge(p0, p1, q0, q1, ratio_p, ratio_q);
  if( dist > delta ) return false;
  if( ratio_p < 0 ) return false;
  if( ratio_p > 1 ) return false;
  if( ratio_q < 0 ) return false;
  if( ratio_q > 1 ) return false;
  const CVector3D& pm = (1-ratio_p)*p0 + ratio_p*p1;
  const CVector3D& qm = (1-ratio_q)*q0 + ratio_q*q1;
  if( (pm-qm).Length() > delta ) return false;
  return true;
}

// ２つのノードの下の階層ので接触する要素を抽出
void GetContactElement_Proximity
(std::set<CContactElement>& aContactElem,
 ////
 double delta,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 int ibvh0,
 int ibvh1,
 const std::vector<CNodeBVH>& aBVH,
 const std::vector<CAABB3D>& aBB)
{
  assert( ibvh0 < aBB.size() );
  assert( ibvh1 < aBB.size() );
  if( !aBB[ibvh0].IsIntersect(aBB[ibvh1]) ) return;
  const int ichild0_0 = aBVH[ibvh0].ichild[0];
  const int ichild0_1 = aBVH[ibvh0].ichild[1];
  const int ichild1_0 = aBVH[ibvh1].ichild[0];
  const int ichild1_1 = aBVH[ibvh1].ichild[1];
  const bool is_leaf0 = (ichild0_1 == -1);
  const bool is_leaf1 = (ichild1_1 == -1);
  if(      !is_leaf0 && !is_leaf1 ){
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0_0,ichild1_0, aBVH,aBB);
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0_1,ichild1_0, aBVH,aBB);
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0_0,ichild1_1, aBVH,aBB);
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0_1,ichild1_1, aBVH,aBB);
  }
  else if( !is_leaf0 &&  is_leaf1 ){
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0_0,ibvh1,aBVH,aBB);
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0_1,ibvh1,aBVH,aBB);
  }
  else if(  is_leaf0 && !is_leaf1 ){
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ibvh0,ichild1_0,aBVH,aBB);
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ibvh0,ichild1_1,aBVH,aBB);
  }
  else if(  is_leaf0 &&  is_leaf1 ){
    const int itri = ichild0_0;
    const int jtri = ichild1_0;
    const int in0 = aTri[itri*3+0];
    const int in1 = aTri[itri*3+1];
    const int in2 = aTri[itri*3+2];
    const int jn0 = aTri[jtri*3+0];
    const int jn1 = aTri[jtri*3+1];
    const int jn2 = aTri[jtri*3+2];
    const CVector3D p0(aXYZ[in0*3+0], aXYZ[in0*3+1], aXYZ[in0*3+2]);
    const CVector3D p1(aXYZ[in1*3+0], aXYZ[in1*3+1], aXYZ[in1*3+2]);
    const CVector3D p2(aXYZ[in2*3+0], aXYZ[in2*3+1], aXYZ[in2*3+2]);
    const CVector3D q0(aXYZ[jn0*3+0], aXYZ[jn0*3+1], aXYZ[jn0*3+2]);
    const CVector3D q1(aXYZ[jn1*3+0], aXYZ[jn1*3+1], aXYZ[jn1*3+2]);
    const CVector3D q2(aXYZ[jn2*3+0], aXYZ[jn2*3+1], aXYZ[jn2*3+2]);
    if( IsContact_FV_Proximity(   in0,in1,in2,jn0, p0,p1,p2,q0, aBB[ichild0_0], delta) ){
      aContactElem.insert( CContactElement(true,    in0,in1,in2,jn0) );
    }
    if( IsContact_FV_Proximity(   in0,in1,in2,jn1, p0,p1,p2,q1, aBB[ichild0_0], delta) ){
      aContactElem.insert( CContactElement(true,    in0,in1,in2,jn1) );
    }
    if( IsContact_FV_Proximity(   in0,in1,in2,jn2, p0,p1,p2,q2, aBB[ichild0_0], delta) ){
      aContactElem.insert( CContactElement(true,    in0,in1,in2,jn2) );
    }
    if( IsContact_FV_Proximity(   jn0,jn1,jn2,in0, q0,q1,q2,p0, aBB[ichild1_0], delta) ){
      aContactElem.insert( CContactElement(true,    jn0,jn1,jn2,in0) );
    }
    if( IsContact_FV_Proximity(   jn0,jn1,jn2,in1, q0,q1,q2,p1, aBB[ichild1_0], delta) ){
      aContactElem.insert( CContactElement(true,    jn0,jn1,jn2,in1) );
    }
    if( IsContact_FV_Proximity(   jn0,jn1,jn2,in2, q0,q1,q2,p2, aBB[ichild1_0], delta) ){
      aContactElem.insert( CContactElement(true,    jn0,jn1,jn2,in2) );
    }
    ////
    if( IsContact_EE_Proximity(      in0,in1,jn0,jn1, p0,p1,q0,q1, delta) ){
      aContactElem.insert( CContactElement(false,    in0,in1,jn0,jn1) );
    }
    if( IsContact_EE_Proximity(      in0,in1,jn1,jn2, p0,p1,q1,q2, delta) ){
      aContactElem.insert( CContactElement(false,    in0,in1,jn1,jn2) );
    }
    if( IsContact_EE_Proximity(      in0,in1,jn2,jn0, p0,p1,q2,q0, delta) ){
      aContactElem.insert( CContactElement(false,    in0,in1,jn2,jn0) );
    }
    if( IsContact_EE_Proximity(      in1,in2,jn0,jn1, p1,p2,q0,q1, delta) ){
      aContactElem.insert( CContactElement(false,    in1,in2,jn0,jn1) );
    }
    if( IsContact_EE_Proximity(      in1,in2,jn1,jn2, p1,p2,q1,q2, delta) ){
      aContactElem.insert( CContactElement(false,    in1,in2,jn1,jn2) );
    }
    if( IsContact_EE_Proximity(      in1,in2,jn2,jn0, p1,p2,q2,q0, delta) ){
      aContactElem.insert( CContactElement(false,    in1,in2,jn2,jn0) );
    }
    if( IsContact_EE_Proximity(      in2,in0,jn0,jn1, p2,p0,q0,q1, delta) ){
      aContactElem.insert( CContactElement(false,    in2,in0,jn0,jn1) );
    }
    if( IsContact_EE_Proximity(      in2,in0,jn1,jn2, p2,p0,q1,q2, delta) ){
      aContactElem.insert( CContactElement(false,    in2,in0,jn1,jn2) );
    }
    if( IsContact_EE_Proximity(      in2,in0,jn2,jn0, p2,p0,q2,q0, delta) ){
      aContactElem.insert( CContactElement(false,    in2,in0,jn2,jn0) );
    }
  }
}

// この木の階層の中で近くにある要素を抽出
void GetContactElement_Proximity
(std::set<CContactElement>& aContactElem,
 ////
 double delta,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 int ibvh,
 const std::vector<CNodeBVH>& aBVH,
 const std::vector<CAABB3D>& aBB)
{
  const int ichild0 = aBVH[ibvh].ichild[0];
  const int ichild1 = aBVH[ibvh].ichild[1];
  const bool is_leaf = (ichild1 == -1);
  if( is_leaf ) return;
  GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0,ichild1,aBVH,aBB);
  GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0,        aBVH,aBB);
  GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild1,        aBVH,aBB);
}


/* ------------------------------------------------------------------------------------- */

// 三次関数を評価する関数
inline double EvaluateCubic
(double r2, // 入力の値
 double k0, double k1, double k2, double k3) // 三次関数の係数
{
  return k0 + k1*r2 + k2*r2*r2 + k3*r2*r2*r2;
}

// 二分法における三次関数の根を探す範囲を狭める関数
static void BisectRangeCubicRoot
(int& icnt,                  // (in out)何回幅を狭めたかというカウンタ
 double& r0, double& r1, // (in,out)入力の範囲から階を探して、狭められた後の範囲を返す
 double v0, double v1, // (in)入力の範囲の両端における値
 double k0, double k1, double k2, double k3) // 三次関数の係数
{
  icnt--;
  if( icnt <= 0 ) return;
  double r2 = 0.5*(r0+r1); // r2はr0とr1の中点
  double v2 = EvaluateCubic(r2, k0,k1,k2,k3); // v2はr2における値
  if( v0*v2 < 0 ){ r1 = r2; } // r0とr2の間で符号が変化する
  else{            r0 = r2; } // r1とr2の間で符号が変化する
  BisectRangeCubicRoot(icnt,r0,r1,v0,v2,k0,k1,k2,k3); // r0とr1の間でさらに範囲を狭める
}

// 三次関数の根を探す関数
static double FindRootCubic
(double r0, double r1,
 double v0, double v1,
 double k0, double k1, double k2, double k3)
{
  int icnt=15; // １５回範囲を狭める
  BisectRangeCubicRoot(icnt, r0,r1, v0,v1, k0,k1,k2,k3);
  return 0.5*(r0+r1);
}

// ４つの点が同一平面上にならぶような補間係数を探す
double FindCoplanerInterp
(const CVector3D& s0, const CVector3D& s1, const CVector3D& s2, const CVector3D& s3,
 const CVector3D& e0, const CVector3D& e1, const CVector3D& e2, const CVector3D& e3)
{
  const CVector3D x1 = s1-s0;
  const CVector3D x2 = s2-s0;
  const CVector3D x3 = s3-s0;
  const CVector3D v1 = e1-e0-x1;
  const CVector3D v2 = e2-e0-x2;
  const CVector3D v3 = e3-e0-x3;
  // 三次関数の係数の計算
  const double k0 = ScalarTripleProduct(x3,x1,x2);
  const double k1 = ScalarTripleProduct(v3,x1,x2)+ScalarTripleProduct(x3,v1,x2)+ScalarTripleProduct(x3,x1,v2);
  const double k2 = ScalarTripleProduct(v3,v1,x2)+ScalarTripleProduct(v3,x1,v2)+ScalarTripleProduct(x3,v1,v2);
  const double k3 = ScalarTripleProduct(v3,v1,v2);
  double r0=-0.0;
  double r1=+1.0;
  const double f0 = EvaluateCubic(r0,k0,k1,k2,k3);
  const double f1 = EvaluateCubic(r1,k0,k1,k2,k3);
  double det = k2*k2-3*k1*k3;
  if( fabs(k3) < 1.0e-10 && fabs(k2) > 1.0e-10 ){ // quadric function、二次関数
    double r2 = -k1/(2*k2); // 極値をとるr
    const double f2 = EvaluateCubic(r2, k0,k1,k2,k3);
    if( r2 > 0 && r2 < 1 ){
      if(      f0*f2 < 0 ){
        return FindRootCubic(r0,r2, f0,f2, k0,k1,k2,k3);

      }
      else if( f2*f1 < 0 ){
        return FindRootCubic(r2,r1, f2,f1, k0,k1,k2,k3);
      }
    }
  }
  if( det > 0 && fabs(k3) > 1.0e-10 ) // cubic function with two extream value、三次関数で極値がある場合
  {
    double r3 = (-k2-sqrt(det))/(3*k3); // 極値をとる小さい方のr
    const double f3 = EvaluateCubic(r3, k0,k1,k2,k3);
    if( r3 > 0 && r3 < 1 ){
      if(      f0*f3 < 0 ){
        return FindRootCubic(r0,r3, f0,f3, k0,k1,k2,k3);
      }
      else if( f3*f1 < 0 ){
        return FindRootCubic(r3,r1, f3,f1, k0,k1,k2,k3);
      }
    }
    double r4 = (-k2+sqrt(det))/(3*k3); // 極値をとる大きい方のr
    const double f4 = EvaluateCubic(r4, k0,k1,k2,k3);
    if( r3 > 0 && r3 < 1 && r4 > 0 && r4 < 1 ){
      if( f3*f4 < 0 ){
        return FindRootCubic(r3,r4, f3,f4, k0,k1,k2,k3);
      }
    }
    if( r4 > 0 && r4 < 1 ){
      if(      f0*f4 < 0 ){
        return FindRootCubic(r0,r4, f0,f4, k0,k1,k2,k3);
      }
      else if( f4*f1 < 0 ){
        return FindRootCubic(r4,r1, f4,f1, k0,k1,k2,k3);
      }
    }
  }
  // monotonus function、0と１の間で短調増加関数
  if( f0*f1 > 0 ){ return -1; } // 根がない場合
  return FindRootCubic(r0,r1, f0,f1, k0,k1,k2,k3);
}


// CCDのFVで接触する要素を検出
bool IsContact_FV_CCD
(int ino0,        int ino1,        int ino2,        int ino3,
 const CVector3D& p0, const CVector3D& p1, const CVector3D& p2, const CVector3D& p3,
 const CVector3D& q0, const CVector3D& q1, const CVector3D& q2, const CVector3D& q3,
 const CAABB3D& bb)
{
  double eps = 1.0e-10;
  if( ino3 == ino0 || ino3 == ino1 || ino3 == ino2 ){ return false; }
  CAABB3D bbp;
  AddPoint(bbp,p3, eps);
  AddPoint(bbp,q3, eps);
  if( !bb.IsIntersect(bbp) ) return false;
  { // CSAT
    CVector3D n = Cross(p1-p0,p2-p0);
    double t0 = Dot(p0-p3,n);
    double t1 = Dot(q0-q3,n);
    double t2 = Dot(q1-q3,n);
    double t3 = Dot(q2-q3,n);
    if( t0*t1 > 0 && t0*t2 > 0 && t0*t3 > 0 ){ return false; }
  }
  double r0,r1;
  double dist = DistanceFaceVertex(p0, p1, p2, p3, r0,r1);
  {
    double vn0 = (p0-q0).Length();
    double vn1 = (p1-q1).Length();
    double vn2 = (p2-q2).Length();
    double vn3 = (p3-q3).Length();
    double vnt = ( vn0 > vn1 ) ? vn0 : vn1;
    vnt = ( vn2 > vnt ) ? vn2 : vnt;
    double max_app = (vnt+vn3);
    ////
    const double r2 = 1-r0-r1;
    if( dist > max_app ) return false;
    if( r0 < 0 || r0 > 1 || r1 < 0 || r1 > 1 || r2 < 0 || r2 > 1 ){
      double dist01 = (GetMinDist_LineSegPoint(p3, p0, p1)-p3).Length();
      double dist12 = (GetMinDist_LineSegPoint(p3, p1, p2)-p3).Length();
      double dist20 = (GetMinDist_LineSegPoint(p3, p2, p0)-p3).Length();
      if( dist01 > max_app && dist12 > max_app && dist20 > max_app ){ return false; }
    }
  }
  double t = FindCoplanerInterp(p0,p1,p2,p3, q0,q1,q2,q3);
  if( t < 0 || t > 1 ) return false;
  CVector3D p0m = (1-t)*p0 + t*q0;
  CVector3D p1m = (1-t)*p1 + t*q1;
  CVector3D p2m = (1-t)*p2 + t*q2;
  CVector3D p3m = (1-t)*p3 + t*q3;
  double w0, w1;
  DistanceFaceVertex(p0m, p1m, p2m, p3m, w0,w1);
  double w2 = 1-w0-w1;
  if( w0 < 0 || w0 > 1 ) return false;
  if( w1 < 0 || w1 > 1 ) return false;
  if( w2 < 0 || w2 > 1 ) return false;
  return true;
}


// CCDのEEで接触する要素を検出
bool IsContact_EE_CCD
(int ino0,         int ino1,         int jno0,         int jno1,
 const CVector3D& p0s, const CVector3D& p1s, const CVector3D& q0s, const CVector3D& q1s,
 const CVector3D& p0e, const CVector3D& p1e, const CVector3D& q0e, const CVector3D& q1e)
{
  double eps = 1.0e-10;
  if( ino0 == jno0 || ino0 == jno1 || ino1 == jno0 || ino1 == jno1 ) return false;
  CAABB3D bbq;
  AddPoint(bbq,q0s, eps);
  AddPoint(bbq,q1s, eps);
  AddPoint(bbq,q0e, eps);
  AddPoint(bbq,q1e, eps);
  
  CAABB3D bbp;
  AddPoint(bbp,p0s, eps);
  AddPoint(bbp,p1s, eps);
  AddPoint(bbp,p0e, eps);
  AddPoint(bbp,p1e, eps);
  if( !bbp.IsIntersect(bbq) ) return false;
  /*
   ////
   { // CSAT
   CVector3D n = Cross(p1s-p0s,q1s-q0s);
   double t0 = Dot(q0s-p0s,n);
   double t1 = Dot(q0e-p0e,n);
   double t2 = Dot(q1e-p0e,n);
   double t3 = Dot(q0e-p1e,n);
   double t4 = Dot(q1e-p1e,n);
   if( t0*t1 > 0 && t0*t2 > 0 && t0*t3 > 0 && t0*t4 > 0 ){ return false; }
   }
   double app_max;
   {
   double vnp0 = (p0e-p0s).Length();
   double vnp1 = (p1e-p1s).Length();
   double vnq0 = (q0e-q0s).Length();
   double vnq1 = (q1e-q1s).Length();
   double vnp = ( vnp0 > vnp1 ) ? vnp0 : vnp1;
   double vnq = ( vnq0 > vnq1 ) ? vnq0 : vnq1;
   app_max =(vnp+vnq);
   }
   ////
   double r0,r1;
   double dist = DistanceEdgeEdge(p0s, p1s, q0s, q1s, r0,r1);
   {
   if( dist > app_max ) return false;
   if( r0 < 0 || r0 > 1 || r1 < 0 || r1 > 1 ){
   CVector3D prj_p0s = GetMinDist_LineSegPoint(p0s, q0s, q1s);
   double dist0 = (prj_p0s-p0s).Length();
   CVector3D prj_p1s = GetMinDist_LineSegPoint(p1s, q0s, q1s);
   double dist1 = (prj_p1s-p1s).Length();
   CVector3D prj_q0s = GetMinDist_LineSegPoint(q0s, p0s, p1s);
   double dist2 = (prj_q0s-q0s).Length();
   CVector3D prj_q1s = GetMinDist_LineSegPoint(q1s, p0s, p1s);
   double dist3 = (prj_q1s-q1s).Length();
   if( dist0 > app_max && dist1 > app_max && dist2 > app_max && dist3 > app_max ) return false;
   }
   }
   */
  const double t = FindCoplanerInterp(p0s,p1s,q0s,q1s, p0e,p1e,q0e,q1e);
  if( t < 0 || t > 1 ) return false;
  CVector3D p0m = (1-t)*p0s + t*p0e;
  CVector3D p1m = (1-t)*p1s + t*p1e;
  CVector3D q0m = (1-t)*q0s + t*q0e;
  CVector3D q1m = (1-t)*q1s + t*q1e;
  double w0,w1;
  double dist = DistanceEdgeEdge(p0m, p1m, q0m, q1m, w0,w1);
  if( w0 < 0 || w0 > 1 ) return false;
  if( w1 < 0 || w1 > 1 ) return false;
  if( dist > 1.0e-2 ) return false;
  return true;
}


// CCDで接触する要素を検出
void GetContactElement_CCD
(std::set<CContactElement>& aContactElem,
 ////
 double dt,
 double delta,
 const std::vector<double>& aXYZ,
 const std::vector<double>& aUVW,
 const std::vector<int>& aTri,
 int ibvh0,
 int ibvh1,
 const std::vector<CNodeBVH>& aBVH,
 const std::vector<CAABB3D>& aBB)
{
  assert( ibvh0 < aBB.size() );
  assert( ibvh1 < aBB.size() );
  if( !aBB[ibvh0].IsIntersect(aBB[ibvh1]) ) return;
  const int ichild0_0 = aBVH[ibvh0].ichild[0];
  const int ichild0_1 = aBVH[ibvh0].ichild[1];
  const int ichild1_0 = aBVH[ibvh1].ichild[0];
  const int ichild1_1 = aBVH[ibvh1].ichild[1];
  const bool is_leaf0 = (ichild0_1 == -1);
  const bool is_leaf1 = (ichild1_1 == -1);
  if(      !is_leaf0 && !is_leaf1 ){
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0_0,ichild1_0, aBVH,aBB);
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0_1,ichild1_0, aBVH,aBB);
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0_0,ichild1_1, aBVH,aBB);
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0_1,ichild1_1, aBVH,aBB);
  }
  else if( !is_leaf0 &&  is_leaf1 ){
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0_0,ibvh1,     aBVH,aBB);
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0_1,ibvh1,     aBVH,aBB);
  }
  else if(  is_leaf0 && !is_leaf1 ){
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ibvh0,    ichild1_0, aBVH,aBB);
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ibvh0,    ichild1_1, aBVH,aBB);
  }
  else if(  is_leaf0 &&  is_leaf1 ){
    const int itri = ichild0_0;
    const int jtri = ichild1_0;
    int in0 = aTri[itri*3+0];
    int in1 = aTri[itri*3+1];
    int in2 = aTri[itri*3+2];
    int jn0 = aTri[jtri*3+0];
    int jn1 = aTri[jtri*3+1];
    int jn2 = aTri[jtri*3+2];
    const CVector3D p0s(aXYZ[in0*3+0],                  aXYZ[in0*3+1],                  aXYZ[in0*3+2]);
    const CVector3D p1s(aXYZ[in1*3+0],                  aXYZ[in1*3+1],                  aXYZ[in1*3+2]);
    const CVector3D p2s(aXYZ[in2*3+0],                  aXYZ[in2*3+1],                  aXYZ[in2*3+2]);
    const CVector3D q0s(aXYZ[jn0*3+0],                  aXYZ[jn0*3+1],                  aXYZ[jn0*3+2]);
    const CVector3D q1s(aXYZ[jn1*3+0],                  aXYZ[jn1*3+1],                  aXYZ[jn1*3+2]);
    const CVector3D q2s(aXYZ[jn2*3+0],                  aXYZ[jn2*3+1],                  aXYZ[jn2*3+2]);
    const CVector3D p0e(aXYZ[in0*3+0]+dt*aUVW[in0*3+0], aXYZ[in0*3+1]+dt*aUVW[in0*3+1], aXYZ[in0*3+2]+dt*aUVW[in0*3+2]);
    const CVector3D p1e(aXYZ[in1*3+0]+dt*aUVW[in1*3+0], aXYZ[in1*3+1]+dt*aUVW[in1*3+1], aXYZ[in1*3+2]+dt*aUVW[in1*3+2]);
    const CVector3D p2e(aXYZ[in2*3+0]+dt*aUVW[in2*3+0], aXYZ[in2*3+1]+dt*aUVW[in2*3+1], aXYZ[in2*3+2]+dt*aUVW[in2*3+2]);
    const CVector3D q0e(aXYZ[jn0*3+0]+dt*aUVW[jn0*3+0], aXYZ[jn0*3+1]+dt*aUVW[jn0*3+1], aXYZ[jn0*3+2]+dt*aUVW[jn0*3+2]);
    const CVector3D q1e(aXYZ[jn1*3+0]+dt*aUVW[jn1*3+0], aXYZ[jn1*3+1]+dt*aUVW[jn1*3+1], aXYZ[jn1*3+2]+dt*aUVW[jn1*3+2]);
    const CVector3D q2e(aXYZ[jn2*3+0]+dt*aUVW[jn2*3+0], aXYZ[jn2*3+1]+dt*aUVW[jn2*3+1], aXYZ[jn2*3+2]+dt*aUVW[jn2*3+2]);
    
    if( IsContact_FV_CCD(      in0,in1,in2,jn0, p0s,p1s,p2s,q0s, p0e,p1e,p2e,q0e, aBB[ibvh0]) ){
      aContactElem.insert( CContactElement(true, in0,in1,in2,jn0) );
    }
    if( IsContact_FV_CCD(      in0,in1,in2,jn1, p0s,p1s,p2s,q1s, p0e,p1e,p2e,q1e, aBB[ibvh0]) ){
      aContactElem.insert( CContactElement(true, in0,in1,in2,jn1) );
    }
    if( IsContact_FV_CCD(      in0,in1,in2,jn2, p0s,p1s,p2s,q2s, p0e,p1e,p2e,q2e, aBB[ibvh0]) ){
      aContactElem.insert( CContactElement(true, in0,in1,in2,jn2) );
    }
    if( IsContact_FV_CCD(      jn0,jn1,jn2,in0, q0s,q1s,q2s,p0s, q0e,q1e,q2e,p0e, aBB[ibvh1]) ){
      aContactElem.insert( CContactElement(true, jn0,jn1,jn2,in0) );
    }
    if( IsContact_FV_CCD(      jn0,jn1,jn2,in1, q0s,q1s,q2s,p1s, q0e,q1e,q2e,p1e, aBB[ibvh1]) ){
      aContactElem.insert( CContactElement(true, jn0,jn1,jn2,in1) );
    }
    if( IsContact_FV_CCD(      jn0,jn1,jn2,in2, q0s,q1s,q2s,p2s, q0e,q1e,q2e,p2e, aBB[ibvh1]) ){
      aContactElem.insert( CContactElement(true, jn0,jn1,jn2,in2) );
    }
    ////
    if( IsContact_EE_CCD(          in0,in1,jn0,jn1, p0s,p1s,q0s,q1s,  p0e,p1e,q0e,q1e) ){
      aContactElem.insert( CContactElement(false,  in0,in1,jn0,jn1) );
    }
    if( IsContact_EE_CCD(          in0,in1,jn1,jn2, p0s,p1s,q1s,q2s,  p0e,p1e,q1e,q2e) ){
      aContactElem.insert( CContactElement(false,  in0,in1,jn1,jn2) );
    }
    if( IsContact_EE_CCD(          in0,in1,jn2,jn0, p0s,p1s,q2s,q0s,  p0e,p1e,q2e,q0e) ){
      aContactElem.insert( CContactElement(false,  in0,in1,jn2,jn0) );
    }
    if( IsContact_EE_CCD(          in1,in2,jn0,jn1, p1s,p2s,q0s,q1s,  p1e,p2e,q0e,q1e) ){
      aContactElem.insert( CContactElement(false,  in1,in2,jn0,jn1) );
    }
    if( IsContact_EE_CCD(          in1,in2,jn1,jn2, p1s,p2s,q1s,q2s,  p1e,p2e,q1e,q2e) ){
      aContactElem.insert( CContactElement(false,  in1,in2,jn1,jn2) );
    }
    if( IsContact_EE_CCD(          in1,in2,jn2,jn0, p1s,p2s,q2s,q0s,  p1e,p2e,q2e,q0e) ){
      aContactElem.insert( CContactElement(false,  in1,in2,jn2,jn0) );
    }
    if( IsContact_EE_CCD(          in2,in0,jn0,jn1, p2s,p0s,q0s,q1s,  p2e,p0e,q0e,q1e) ){
      aContactElem.insert( CContactElement(false,  in2,in0,jn0,jn1) );
    }
    if( IsContact_EE_CCD(          in2,in0,jn1,jn2, p2s,p0s,q1s,q2s,  p2e,p0e,q1e,q2e) ){
      aContactElem.insert( CContactElement(false,  in2,in0,jn1,jn2) );
    }
    if( IsContact_EE_CCD(          in2,in0,jn2,jn0, p2s,p0s,q2s,q0s,  p2e,p0e,q2e,q0e) ){
      aContactElem.insert( CContactElement(false,  in2,in0,jn2,jn0) );
    }
  }
}

void GetContactElement_CCD
(std::set<CContactElement>& aContactElem,
 ////
 double dt,
 double delta,
 const std::vector<double>& aXYZ,
 const std::vector<double>& aUVW,
 const std::vector<int>& aTri,
 int ibvh,
 const std::vector<CNodeBVH>& aBVH,
 const std::vector<CAABB3D>& aBB)
{
  const int ichild0 = aBVH[ibvh].ichild[0];
  const int ichild1 = aBVH[ibvh].ichild[1];
  const bool is_leaf = (ichild1 == -1);
  if( is_leaf ) return;
  GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0,ichild1,aBVH,aBB);
  GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0,        aBVH,aBB);
  GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild1,        aBVH,aBB);
}

