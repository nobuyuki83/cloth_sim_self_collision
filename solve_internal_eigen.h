//
//  solve_internal_eigen.h
//  self_contact
//
//  Created by Nobuyuki Umetani on 11/17/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//

#ifndef solve_internal_eigen_h
#define solve_internal_eigen_h

#include "Eigen/Dense"
#include "cloth_internal_physics.h"



// compute total energy and its first and second derivatives
// 全体のエネルギーとその，節点位置における一階微分，二階微分を計算
void AddWdWddW_Cloth
(double& W, // (out) energy，歪エネルギー
 Eigen::VectorXd& dW, // (out) first derivative of energy，歪エネルギーの一階微分
 Eigen::MatrixXd& ddW, // (out) second derivative of energy，歪エネルギーの二階微分
 ////
 const std::vector<double>& aXYZ, // (in) deformed vertex positions，現在の頂点の座標配列
 const std::vector<double>& aXYZ0, // (in) initial vertex positions，変形前の頂点の座標配列
 const std::vector<int>& aTri, // (in) triangle index，三角形の頂点インデックス配列
 const std::vector<int>& aQuad, // (in) index of 4 vertices required for bending，曲げ計算のための４頂点のインデックス配列
 double lambda, // (in) Lame's 1st parameter，ラメの第一定数
 double myu,  // (in) Lame's 2nd parameter　ラメの第二定数
 double stiff_bend // (in) bending stiffness，曲げ剛性
 )
{
  // marge element in-plane strain energy
  // 面内歪エネルギーを追加
  for(int itri=0;itri<aTri.size()/3;itri++){
    const int aIP[3] = { aTri[itri*3+0], aTri[itri*3+1], aTri[itri*3+2] };
    double C[3][3]; double c[3][3];
    for(int ino=0;ino<3;ino++){
      const int ip = aIP[ino];
      for(int i=0;i<3;i++){ C[ino][i] = aXYZ0[ip*3+i]; }
      for(int i=0;i<3;i++){ c[ino][i] = aXYZ [ip*3+i]; }
    }
    double e, de[3][3], dde[3][3][3][3];
    WdWddW_CST( e,de,dde, C,c, lambda,myu );
    W += e;  // marge energy
    // marge de
    for(int ino=0;ino<3;ino++){
      const int ip = aIP[ino];
      for(int i =0;i<3;i++){ dW(ip*3+i) += de[ino][i]; }
    }
    // marge dde
    for(int ino=0;ino<3;ino++){
      const int ip = aIP[ino];
      for(int jno=0;jno<3;jno++){
        const int jp = aIP[jno];
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            ddW(ip*3+i,jp*3+j) += dde[ino][jno][i][j];
          }
        }
      }
    }
  }
  
  // marge element bending energy
  // 曲げエネルギーを追加
  for(int iq=0;iq<aQuad.size()/4;iq++){
    const int aIP[4] = { aQuad[iq*4+0], aQuad[iq*4+1], aQuad[iq*4+2], aQuad[iq*4+3] };
    double C[4][3]; double c[4][3];
    for(int ino=0;ino<4;ino++){
      const int ip = aIP[ino];
      for(int i=0;i<3;i++){ C[ino][i] = aXYZ0[ip*3+i]; }
      for(int i=0;i<3;i++){ c[ino][i] = aXYZ [ip*3+i]; }
    }
    double e, de[4][3], dde[4][4][3][3];
    WdWddW_Bend( e,de,dde, C,c, stiff_bend );
    W += e;  // marge energy
    // marge de
    for(int ino=0;ino<4;ino++){
      const int ip = aIP[ino];
      for(int i =0;i<3;i++){ dW(ip*3+i) += de[ino][i]; }
    }
    // marge dde
    for(int ino=0;ino<4;ino++){
      const int ip = aIP[ino];
      for(int jno=0;jno<4;jno++){
        const int jp = aIP[jno];
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            ddW(ip*3+i,jp*3+j) += dde[ino][jno][i][j];
          }
        }
      }
    }
  }
}

// compute total energy and its first and second derivatives
// 全体のエネルギーとその，節点位置における一階微分，二階微分を計算
void AddWdWddW_Contact
(double& W, // (out) energy，歪エネルギー
 Eigen::VectorXd& dW, // (out) first derivative of energy，歪エネルギーの一階微分
 Eigen::MatrixXd& ddW, // (out) second derivative of energy，歪エネルギーの二階微分
 ////
 const std::vector<double>& aXYZ, // (in) deformed vertex positions，現在の頂点の座標配列
 double stiff_contact,
 double contact_clearance,
 void (*penetrationDepth)(double& , double* , const double*)
 )
{
  // marge point-wise external contact energy
  // 面内歪エネルギーを追加
  for(int ip=0;ip<aXYZ.size()/3;ip++){
    double c[3] = { aXYZ[ip*3+0], aXYZ[ip*3+1], aXYZ[ip*3+2] };
    double e, de[3], dde[3][3];
    WdWddW_Contact( e,de,dde, c, stiff_contact,contact_clearance, penetrationDepth );
    W += e;  // marge energy
    // marge de
    for(int i =0;i<3;i++){ dW(ip*3+i) += de[i]; }
    // marge dde
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        ddW(ip*3+i,ip*3+j) += dde[i][j];
      }
    }
  }
}

// compute total energy and its first and second derivatives
// 全体のエネルギーとその，節点位置における一階微分，二階微分を計算
void AddWdWddW_Gravity
(double& W, // (out) energy，歪エネルギー
 Eigen::VectorXd& dW, // (out) first derivative of energy，歪エネルギーの一階微分
 Eigen::MatrixXd& ddW, // (out) second derivative of energy，歪エネルギーの二階微分
 ////
 const std::vector<double>& aXYZ, // (in) deformed vertex positions，現在の頂点の座標配列
 const double gravity[3], // (in) gravitational accerelation，重力加速度
 double mass_point // (in) mass of a point，頂点の質量
 )
{
  // marge potential energy
  // 重力ポテンシャル・エネルギーを追加
  for(int ip=0;ip<aXYZ.size()/3;ip++){
    const double c[3] = {aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]};
    W -= mass_point*( c[0]*gravity[0] + c[1]*gravity[1] + c[2]*gravity[2] );
    dW(ip*3+0) -= mass_point*gravity[0];
    dW(ip*3+1) -= mass_point*gravity[1];
    dW(ip*3+2) -= mass_point*gravity[2];
  }
}

void StepTime_InternalDynamics
(
 std::vector<double>& aXYZ, // (in,out) deformed vertex positions，現在の頂点位置配列
 std::vector<double>& aUVW, // (in,out) deformed vertex velocity，現在の頂点速度配列
 ////
 const std::vector<double>& aXYZ0,// (in) initial vertex positions，変形前の頂点の座標配列
 const std::vector<int>& aBCFlag, // (in) boundary condition flag (0:free 1:fixed)，境界条件フラグの配列
 const std::vector<int>& aTri, // (in) triangle index，三角形の頂点インデックス配列
 const std::vector<int>& aQuad, // (in) index of 4 vertices required for bending，曲げ計算のための４頂点のインデックス配列
 const double dt, // (in) size of time step，時間ステップの大きさ
 double lambda, // (in) Lame's 1st parameter，ラメ第一定数
 double myu, // (in) Lame's 2nd parameter，ラメ第二定数
 double stiff_bend, // (in) bending stiffness 曲げ剛性
 const double gravity[3], // (in) gravitatinal accereration，重力加速度
 double mass_point, // (in) mass for a point，頂点あたりの質量
 double stiff_contact, // (in) stiffness of contact, 接触の硬さ
 double contact_clearance, // (in) clearance distance of contact, 接触距離
 void (*penetrationDepth)(double& , double* , const double*) // (in) 節点位置から接触物との距離を返す関数
 )
{
  const int np = (int)aXYZ.size()/3; // number of point，頂点数
  const int nDof = np*3; // degree of freedom，全自由度数
  // compute total energy and its first and second derivatives
  // 全体のエネルギーとその，節点位置における一階微分，二階微分を計算
  double W = 0; // 全体のエネルギー
  Eigen::VectorXd dW(nDof);        dW.setConstant(0.0); // 全体のエネルギーの節点での一階微分
  Eigen::MatrixXd ddW(nDof,nDof);  ddW.setConstant(0.0); // 全体のエネルギーの節点での二階微分
  AddWdWddW_Cloth(W,dW,ddW, // 布のエネルギーとその微分を追加
                  aXYZ,aXYZ0,
                  aTri,aQuad,
                  lambda,myu,stiff_bend);
  AddWdWddW_Contact(W,dW,ddW, // 接触のエネルギーとその微分を追加
                    aXYZ,
                    stiff_contact,contact_clearance,penetrationDepth);
  AddWdWddW_Gravity(W,dW,ddW, // 重力のエネルギーとその微分を追加
                    aXYZ,
                    gravity,mass_point);
  std::cout << "energy : " << W << std::endl;
  // compute coefficient matrix and left-hand-side vector
  // Back-ward Eular time integration
  Eigen::VectorXd b(nDof); // left-hand-size vector，右辺ベクトル
  Eigen::MatrixXd A(nDof,nDof); // coefficient matrix，係数行列
  b = -dW;
  for(int i=0;i<nDof;i++){
    b(i) += mass_point*aUVW[i]/dt;
  }
  A = ddW;
  for(int i=0;i<nDof;i++){
    A(i,i) += mass_point / (dt*dt);
  }
  // put fixed boundary condition
  for(int ip=0;ip<np;ip++){
    if( aBCFlag[ip] == 0 ) continue;
    for(int i=0;i<3;i++){  b(ip*3+i) = 0; }
    for(int jp=0;jp<np;jp++){
      for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
          A(ip*3+i,jp*3+j) = 0;
          A(jp*3+j,ip*3+i) = 0;
        }
      }
    }
    for(int i=0;i<3;i++){ A(ip*3+i,ip*3+i) = 1.0; }
  }
  Eigen::VectorXd x = A.partialPivLu().solve(b);   // solve linear system，連立一次方程式を解く
  for(int i=0;i<nDof;i++){ aXYZ[i] += x(i); }  // update position，頂点位置の更新
  for(int i=0;i<nDof;i++){ aUVW[i] = x(i)/dt; }  // update velocity，頂点の速度の更新
}


void UpdateIntermidiateVelocity
(std::vector<double>& aUVW, // (in,out) deformed vertex velocity，現在の頂点速度配列
 ////
 const std::vector<double>& aXYZ, // (in,out) deformed vertex positions，現在の頂点位置配列
 const std::vector<double>& aXYZ0,// (in) initial vertex positions，変形前の頂点の座標配列
 const std::vector<int>& aBCFlag, // (in) boundary condition flag (0:free 1:fixed)，境界条件フラグの配列
 const std::vector<int>& aTri, // (in) triangle index，三角形の頂点インデックス配列
 const std::vector<int>& aQuad, // (in) index of 4 vertices required for bending，曲げ計算のための４頂点のインデックス配列
 const double dt, // (in) size of time step，時間ステップの大きさ
 double lambda, // (in) Lame's 1st parameter，ラメ第一定数
 double myu, // (in) Lame's 2nd parameter，ラメ第二定数
 double stiff_bend, // (in) bending stiffness 曲げ剛性
 const double gravity[3], // (in) gravitatinal accereration，重力加速度
 double mass_point // (in) mass for a point，頂点あたりの質量
 )
{
  const int np = (int)aXYZ.size()/3; // number of point，頂点数
  const int nDof = np*3; // degree of freedom，全自由度数
  // compute total energy and its first and second derivatives
  // 全体のエネルギーとその，節点位置における一階微分，二階微分を計算
  double W = 0;
  Eigen::VectorXd dW(nDof);        dW.setConstant(0.0);
  Eigen::MatrixXd ddW(nDof,nDof);  ddW.setConstant(0.0);
  AddWdWddW_Cloth(W,dW,ddW,
                  aXYZ,aXYZ0,
                  aTri,aQuad,
                  lambda,myu,stiff_bend);
  AddWdWddW_Gravity(W,dW,ddW,
                    aXYZ,
                    gravity,mass_point);
  for(int ip=0;ip<np;ip++){
    if( aBCFlag[ip] == 0 ) continue;
    aUVW[ip*3+0] = 0;
    aUVW[ip*3+1] = 0;
    aUVW[ip*3+2] = 0;
  }
  for(unsigned int ip=0;ip<np;ip++){
    aUVW[ip*3+0] += -(0.5*dt/mass_point)*dW[ip*3+0];
    aUVW[ip*3+1] += -(0.5*dt/mass_point)*dW[ip*3+1];
    aUVW[ip*3+2] += -(0.5*dt/mass_point)*dW[ip*3+2];
  }
}



#endif
