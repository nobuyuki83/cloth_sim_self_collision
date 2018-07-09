//
//  cloth_internal_physics.h
//
//  布の内部物理を解くための関数群
//
//  Created by Nobuyuki Umetani on 11/7/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//

#ifndef cloth_internal_physics_h
#define cloth_internal_physics_h

#include "utility.h"

void MakePositiveDefinite_Sim22(const double s2[3],double s3[3])
{
  const double b = (s2[0]+s2[1])*0.5;
  const double d = (s2[0]-s2[1])*(s2[0]-s2[1])*0.25 + s2[2]*s2[2];
  const double e = sqrt(d);
  if( b-e > 1.0e-20 ){
    s3[0] = s2[0];
    s3[1] = s2[1];
    s3[2] = s2[2];
    return;
  }
  if( b+e < 0 ){
    s3[0] = 0;
    s3[1] = 0;
    s3[2] = 0;
    return;
  }
  const double l = b+e;
  double t0[2] = { s2[0]-l, s2[2]   };
  double t1[2] = { s2[2],   s2[1]-l };
  //  std::cout << t0[0]*t1[1]-t0[1]*t1[0] << std::endl;
  const double sqlen_t0 = t0[0]*t0[0]+t0[1]*t0[1];
  const double sqlen_t1 = t1[0]*t1[0]+t1[1]*t1[1];
  if( sqlen_t0 > sqlen_t1 ){
    if( sqlen_t0 < 1.0e-20 ){
      s3[0] = 0;
      s3[1] = 0;
      s3[2] = 0;
      return;
    }
    const double invlen_t0 = 1.0/sqrt(sqlen_t0);
    t0[0] *= invlen_t0;
    t0[1] *= invlen_t0;
    s3[0] = l*t0[0]*t0[0];
    s3[1] = l*t0[1]*t0[1];
    s3[2] = l*t0[0]*t0[1];
  }
  else{
    if( sqlen_t1 < 1.0e-20 ){
      s3[0] = 0;
      s3[1] = 0;
      s3[2] = 0;
      return;
    }
    const double invlen_t1 = 1.0/sqrt(sqlen_t1);
    t1[0] *= invlen_t1;
    t1[1] *= invlen_t1;
    s3[0] = l*t1[0]*t1[0];
    s3[1] = l*t1[1]*t1[1];
    s3[2] = l*t1[0]*t1[1];
  }
  return;
}


void WdWddW_CST
(double& W, // (out) energy，歪エネルギー
 double dW[3][3], // (out) 1st derivative of energy，歪エネルギーの一階微分
 double ddW[3][3][3][3], // (out) 2nd derivative of energy，歪エネルギーの二階微分
 ////
 const double C[3][3], // (in) undeformed triangle vertex positions，変形前の三角形の三頂点の位置
 const double c[3][3], // (in) deformed triangle vertex positions，変形後の三角形の山頂点の位置
 const double lambda, // (in) Lame's 1st parameter，ラメ第一定数
 const double myu     // (in) Lame's 2nd parameter，ラメ第二定数
 )
{
	double D[3][3] = { // undeformed edge vector，変形前の辺ベクトル
		{ C[1][0]-C[0][0], C[1][1]-C[0][1], C[1][2]-C[0][2] },
		{ C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] }, { 0,0,0 } };
  double Area;
  UnitNormalAreaTri3D(D[2], Area, C[0], C[1], C[2]);
	
	double G[3][3]; // inverse of D, Dの逆行列
	{
		Cross3D(G[0], D[1], D[2]);
		const double invtmp1 = 1.0/Dot3D(G[0],D[0]);
		G[0][0] *= invtmp1;	G[0][1] *= invtmp1;	G[0][2] *= invtmp1;
		////
		Cross3D(G[1], D[2], D[0]);
		const double invtmp2 = 1.0/Dot3D(G[1],D[1]);
		G[1][0] *= invtmp2;	G[1][1] *= invtmp2;	G[1][2] *= invtmp2;
	}
	
	const double d[2][3] = { // deformed edge vector, 変形後の辺ベクトル
		{ c[1][0]-c[0][0], c[1][1]-c[0][1], c[1][2]-c[0][2] },
		{ c[2][0]-c[0][0], c[2][1]-c[0][1], c[2][2]-c[0][2] } };
  
	const double dNdr[3][2] = { {-1.0, -1.0}, {+1.0, +0.0}, {+0.0, +1.0} };
  
  const double E2[3] = {  // green lagrange strain，グリーンラグランジュ歪(工学歪表記)
		0.5*( Dot3D(d[0],d[0]) - Dot3D(D[0],D[0]) ),
		0.5*( Dot3D(d[1],d[1]) - Dot3D(D[1],D[1]) ),
		1.0*( Dot3D(d[0],d[1]) - Dot3D(D[0],D[1]) ) };
  const double GuGu2[3] = { Dot3D(G[0],G[0]), Dot3D(G[1],G[1]), Dot3D(G[1],G[0]) };
  const double Cons2[3][3] = { // constitutive tensor, 構成則テンソル
    { lambda*GuGu2[0]*GuGu2[0] + 2*myu*(GuGu2[0]*GuGu2[0]),
      lambda*GuGu2[0]*GuGu2[1] + 2*myu*(GuGu2[2]*GuGu2[2]),
      lambda*GuGu2[0]*GuGu2[2] + 2*myu*(GuGu2[0]*GuGu2[2]) },
    { lambda*GuGu2[1]*GuGu2[0] + 2*myu*(GuGu2[2]*GuGu2[2]),
      lambda*GuGu2[1]*GuGu2[1] + 2*myu*(GuGu2[1]*GuGu2[1]),
      lambda*GuGu2[1]*GuGu2[2] + 2*myu*(GuGu2[2]*GuGu2[1]) },
    { lambda*GuGu2[2]*GuGu2[0] + 2*myu*(GuGu2[0]*GuGu2[2]),
      lambda*GuGu2[2]*GuGu2[1] + 2*myu*(GuGu2[2]*GuGu2[1]),
      lambda*GuGu2[2]*GuGu2[2] + 1*myu*(GuGu2[0]*GuGu2[1] + GuGu2[2]*GuGu2[2]) } };
  const double S2[3] = {  // 2nd Piola-Kirchhoff stress，第二ピオラ・キルヒホッフ応力
    Cons2[0][0]*E2[0] + Cons2[0][1]*E2[1] + Cons2[0][2]*E2[2],
    Cons2[1][0]*E2[0] + Cons2[1][1]*E2[1] + Cons2[1][2]*E2[2],
    Cons2[2][0]*E2[0] + Cons2[2][1]*E2[1] + Cons2[2][2]*E2[2] };
  
  // compute energy，歪エネルギー
  W = 0.5*Area*(E2[0]*S2[0] + E2[1]*S2[1] + E2[2]*S2[2]);
  
  // compute 1st derivative，エネルギーの一階微分の計算
  for(int ino=0;ino<3;ino++){
    for(int idim=0;idim<3;idim++){
      dW[ino][idim] = Area*
      (+S2[0]*d[0][idim]*dNdr[ino][0]
       +S2[2]*d[0][idim]*dNdr[ino][1]
       +S2[2]*d[1][idim]*dNdr[ino][0]
       +S2[1]*d[1][idim]*dNdr[ino][1]);
    }
  }
  
  double S3[3] = { S2[0], S2[1], S2[2] };
  MakePositiveDefinite_Sim22(S2,S3);
	
  // compute second derivative，エネルギーの二階微分の計算
	for(int ino=0;ino<3;ino++){
    for(int jno=0;jno<3;jno++){
      for(int idim=0;idim<3;idim++){
        for(int jdim=0;jdim<3;jdim++){
          double dtmp0 = 0;
          dtmp0 += d[0][idim]*dNdr[ino][0]*Cons2[0][0]*d[0][jdim]*dNdr[jno][0];
          dtmp0 += d[0][idim]*dNdr[ino][0]*Cons2[0][1]*d[1][jdim]*dNdr[jno][1];
          dtmp0 += d[0][idim]*dNdr[ino][0]*Cons2[0][2]*d[0][jdim]*dNdr[jno][1];
          dtmp0 += d[0][idim]*dNdr[ino][0]*Cons2[0][2]*d[1][jdim]*dNdr[jno][0];
          dtmp0 += d[1][idim]*dNdr[ino][1]*Cons2[1][0]*d[0][jdim]*dNdr[jno][0];
          dtmp0 += d[1][idim]*dNdr[ino][1]*Cons2[1][1]*d[1][jdim]*dNdr[jno][1];
          dtmp0 += d[1][idim]*dNdr[ino][1]*Cons2[1][2]*d[0][jdim]*dNdr[jno][1];
          dtmp0 += d[1][idim]*dNdr[ino][1]*Cons2[1][2]*d[1][jdim]*dNdr[jno][0];
          dtmp0 += d[0][idim]*dNdr[ino][1]*Cons2[2][0]*d[0][jdim]*dNdr[jno][0];
          dtmp0 += d[0][idim]*dNdr[ino][1]*Cons2[2][1]*d[1][jdim]*dNdr[jno][1];
          dtmp0 += d[0][idim]*dNdr[ino][1]*Cons2[2][2]*d[0][jdim]*dNdr[jno][1];
          dtmp0 += d[0][idim]*dNdr[ino][1]*Cons2[2][2]*d[1][jdim]*dNdr[jno][0];
          dtmp0 += d[1][idim]*dNdr[ino][0]*Cons2[2][0]*d[0][jdim]*dNdr[jno][0];
          dtmp0 += d[1][idim]*dNdr[ino][0]*Cons2[2][1]*d[1][jdim]*dNdr[jno][1];
          dtmp0 += d[1][idim]*dNdr[ino][0]*Cons2[2][2]*d[0][jdim]*dNdr[jno][1];
          dtmp0 += d[1][idim]*dNdr[ino][0]*Cons2[2][2]*d[1][jdim]*dNdr[jno][0];
          ddW[ino][jno][idim][jdim] = dtmp0*Area;
        }
      }
      const double dtmp1 = Area*
      (+S3[0]*dNdr[ino][0]*dNdr[jno][0]
       +S3[2]*dNdr[ino][0]*dNdr[jno][1]
       +S3[2]*dNdr[ino][1]*dNdr[jno][0]
       +S3[1]*dNdr[ino][1]*dNdr[jno][1]);
      ddW[ino][jno][0][0] += dtmp1;
      ddW[ino][jno][1][1] += dtmp1;
      ddW[ino][jno][2][2] += dtmp1;
    }
	}
}

// compute energy and its 1st and 2nd derivative for cloth bending
// 布の曲げエネルギーと，その節点の変位に対する１階と２階微分を求める
void WdWddW_Bend
(double& W,  // (out) energy，歪エネルギー
 double dW[4][3], // (out) 1st derivative of energy，歪エネルギーの一階微分
 double ddW[4][4][3][3], // (out) 2nd derivative of energy，歪エネルギーの二階微分
 ////
 const double C[4][3], // (in) undeformed triangle vertex positions，変形前の辺周りの四頂点の位置
 const double c[4][3], // (in) deformed triangle vertex positions，変形後の辺周りの四頂点の位置
 double stiff)
{
  const double A0 = TriArea3D(C[0],C[2],C[3]);
  const double A1 = TriArea3D(C[1],C[3],C[2]);
  const double L0 = Distance3D(C[2],C[3]);
  const double H0 = A0*2.0/L0;
  const double H1 = A1*2.0/L0;
  const double e23[3] = { C[3][0]-C[2][0], C[3][1]-C[2][1], C[3][2]-C[2][2] };
  const double e02[3] = { C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] };
  const double e03[3] = { C[3][0]-C[0][0], C[3][1]-C[0][1], C[3][2]-C[0][2] };
  const double e12[3] = { C[2][0]-C[1][0], C[2][1]-C[1][1], C[2][2]-C[1][2] };
  const double e13[3] = { C[3][0]-C[1][0], C[3][1]-C[1][1], C[3][2]-C[1][2] };
  double cot023, cot032;
  {
    const double r2 = -Dot3D(e02,e23);
    const double r3 = +Dot3D(e03,e23);
    cot023 = r2/H0;
    cot032 = r3/H0;
  }
  double cot123, cot132;
  {
    const double r2 = -Dot3D(e12,e23);
    const double r3 = +Dot3D(e13,e23);
    cot123 = r2/H1;
    cot132 = r3/H1;
  }
  const double tmp0 = stiff/((A0+A1)*L0*L0);
  const double K[4] = { -cot023-cot032, -cot123-cot132, cot032+cot132, cot023+cot123 };
  
  // compute 2nd derivative of energy，エネルギーの二階微分を求める
  for(int i=0;i<4*4*3*3;i++){ (&ddW[0][0][0][0])[i] = 0; }
  for(int ino=0;ino<4;ino++){
    for(int jno=0;jno<4;jno++){
      const double tmp = K[ino]*K[jno]*tmp0;
      ddW[ino][jno][0][0] = tmp;
      ddW[ino][jno][1][1] = tmp;
      ddW[ino][jno][2][2] = tmp;
    }
  }
  // compute 1st derivative of energy，エネルギーとその一階微分を求める
  for(int ino=0;ino<4;ino++){
    for(int idim=0;idim<3;idim++){
      dW[ino][idim] = 0;
      for(int jno=0;jno<4;jno++){
        for(int jdim=0;jdim<3;jdim++){
          dW[ino][idim] += ddW[ino][jno][idim][jdim]*c[jno][jdim];
        }
      }
      W += dW[ino][idim]*c[ino][idim];
    }
  }
}


// compute energy and its 1st and 2nd derivative for contact against object
// 接触エネルギーと，その節点の変位に対する１階と２階微分を求める
void WdWddW_Contact
(double& W,  // (out) energy，接触エネルギー
 double dW[3], // (out) 1st derivative of energy，接触エネルギーの一階微分
 double ddW[3][3], // (out) 2nd derivative of energy，接触エネルギーの二階微分
 ////
 const double c[3], // (in) deformed triangle vertex positions，変形後の三角形の山頂点の位置
 double stiff_contact,
 double contact_clearance,
 void (*penetrationDepth)(double& , double* , const double*)
 )
{
  double pd, n[3];
  penetrationDepth(pd,n,c);
  pd += contact_clearance;
  if( pd  < 0 ){
    W = 0;
    dW[0] = 0;  dW[1] = 0;  dW[2] = 0;
    ddW[0][0] = 0;  ddW[0][1] = 0;  ddW[0][2] = 0;
    ddW[1][0] = 0;  ddW[1][1] = 0;  ddW[1][2] = 0;
    ddW[2][0] = 0;  ddW[2][1] = 0;  ddW[2][2] = 0;
    return;
  }
  W = 0.5*stiff_contact*pd*pd;
  
  dW[0] = -stiff_contact*pd*n[0];
  dW[1] = -stiff_contact*pd*n[1];
  dW[2] = -stiff_contact*pd*n[2];
  
  ddW[0][0] = stiff_contact*n[0]*n[0];
  ddW[0][1] = stiff_contact*n[0]*n[1];
  ddW[0][2] = stiff_contact*n[0]*n[2];
  ddW[1][0] = stiff_contact*n[1]*n[0];
  ddW[1][1] = stiff_contact*n[1]*n[1];
  ddW[1][2] = stiff_contact*n[1]*n[2];
  ddW[2][0] = stiff_contact*n[2]*n[0];
  ddW[2][1] = stiff_contact*n[2]*n[1];
  ddW[2][2] = stiff_contact*n[2]*n[2];
}







#endif
