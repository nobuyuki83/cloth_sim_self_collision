//
//  main.cpp
//
//  simple_cloth, 単純な布のシミュレーション
//
//  Created by Nobuyuki Umetani on 11/3/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//

// 視点の回転：左ボタンドラッグ
// 衝突オブジェクトの変更：スペースキー
// アニメーションの計算・停止: 'a'キー

#include <iostream>
#include <vector>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "../utility.h" // ベクトル，四元数の演算
#include "../solve_internal_eigen.h" // Eigenライブラリを使った布のソルバ
#include "../self_collision_cloth.h" // 自己衝突を解くライブラリ

/* ------------------------------------------------------------------------ */


// Setting problem here
void SetClothShape_Square
(std::vector<double>& aXYZ0, // (out) undeformed vertex positions，変形前の頂点の位置配列
 std::vector<int>& aBCFlag, // (out) boundary condition flag (0:free 1:fixed)，境界条件フラグ
 std::vector<int>& aTri, // (out) index of triangles，三角形の頂点インデックス
 std::vector<int>& aQuad, // (out) index of 4 vertices required for bending，曲げ計算のための４頂点の配列
 double& total_area, // (out) total area of cloth，布の面積
 ///
 double elem_length, // (in) number of division of the square cloth edge, 一辺の分割数
 double cloth_size_x,
 double cloth_size_z) // (in) size of square cloth，一辺の長さ
{
  // make vertex potision array 頂点位置配列を作る
  const int ndiv_x = cloth_size_x/elem_length; // size of an element
  const int ndiv_z = cloth_size_z/elem_length; // size of an element
  const int nxyz =(ndiv_x+1)*(ndiv_z+1); // number of points
  aXYZ0.reserve( nxyz*3 );
  for(int ix=0;ix<ndiv_x+1;ix++){
    for(int iz=0;iz<ndiv_z+1;iz++){
      aXYZ0.push_back( ix*elem_length );
      aXYZ0.push_back( 0.0 );                  
      aXYZ0.push_back( iz*elem_length );
    }
  }
  
  // make triangle index array, 三角形の頂点配列を作る
  const int ntri = ndiv_x*ndiv_z*2;
  aTri.reserve( ntri*3 );
  for(int ix=0;ix<ndiv_x;ix++){
    for(int iz=0;iz<ndiv_z;iz++){
      aTri.push_back(  ix   *(ndiv_z+1)+ iz    );
      aTri.push_back( (ix+1)*(ndiv_z+1)+ iz    );
      aTri.push_back(  ix   *(ndiv_z+1)+(iz+1) );
      ////
      aTri.push_back( (ix+1)*(ndiv_z+1)+(iz+1) );
      aTri.push_back(  ix   *(ndiv_z+1)+(iz+1) );
      aTri.push_back( (ix+1)*(ndiv_z+1)+ iz    );
    }
  }
  
  aQuad.reserve( ndiv_x*ndiv_z + ndiv_x*(ndiv_z-1) + (ndiv_x-1)*ndiv_z );
  for(int ix=0;ix<ndiv_x;ix++){
    for(int iz=0;iz<ndiv_z;iz++){
      aQuad.push_back( (ix+0)*(ndiv_z+1)+(iz+0) );
      aQuad.push_back( (ix+1)*(ndiv_z+1)+(iz+1) );
      aQuad.push_back( (ix+1)*(ndiv_z+1)+(iz+0) );
      aQuad.push_back( (ix+0)*(ndiv_z+1)+(iz+1) );
    }
  }
  for(int ix=0;ix<ndiv_x;ix++){
    for(int iz=0;iz<ndiv_z-1;iz++){
      aQuad.push_back( (ix+1)*(ndiv_z+1)+(iz+0) );
      aQuad.push_back( (ix+0)*(ndiv_z+1)+(iz+2) );
      aQuad.push_back( (ix+1)*(ndiv_z+1)+(iz+1) );
      aQuad.push_back( (ix+0)*(ndiv_z+1)+(iz+1) );
    }
  }
  for(int ix=0;ix<ndiv_x-1;ix++){
    for(int iz=0;iz<ndiv_z;iz++){
      aQuad.push_back( (ix+0)*(ndiv_z+1)+(iz+1) );
      aQuad.push_back( (ix+2)*(ndiv_z+1)+(iz+0) );
      aQuad.push_back( (ix+1)*(ndiv_z+1)+(iz+0) );
      aQuad.push_back( (ix+1)*(ndiv_z+1)+(iz+1) );
    }
  }
  
  // compute total area 全体の面積を計算
  total_area = cloth_size_x*cloth_size_z;
  
  // set fixed boundary condition 固定境界条件を指定する
  aBCFlag = std::vector<int>(nxyz,0);
  /*
  for(int ix=0;ix<ndiv+1;ix++){
    aBCFlag[ ix*(ndiv+1) + ndiv ] = 1;
  }
   */
}


// active if( imode_contact == 0 )
void penetrationDepth_Plane(double& pd, double* n, const double* p)
{
  n[0] = 0.0;  n[1] = 0.0;  n[2] = 1.0; // normal of the plane
  pd = -0.5 - p[2]; // penetration depth
};

// active if( imode_contact == 1 )
void penetrationDepth_Sphere(double& pd, double* n, const double* p)
{
  const double center[3] = { 0.5, 0.0, +2 };
  const double radius = 3.0;
  n[0] = p[0]-center[0];
  n[1] = p[1]-center[1];
  n[2] = p[2]-center[2];
  double len = Length3D(n);
  n[0] /= len;
  n[1] /= len;
  n[2] /= len;
  pd = radius-len; // penetration depth
  n[0] *= -1;
  n[1] *= -1;
  n[2] *= -1;
  pd *= -1;
};


/* ------------------------------------------------------------------------ */
// input parameter for simulation
double elem_length = 0.08;
double cloth_size_x = 0.2;
double cloth_size_z = 3.0;
std::vector<double> aXYZ0; // undeformed vertex positions，変形前の頂点の位置配列
std::vector<double> aXYZ; // deformed vertex positions，変形中の頂点の位置配列
std::vector<double> aUVW; // deformed vertex velocity，変形中の頂点の速度
std::vector<int> aBCFlag;  // boundary condition flag (0:free 1:fixed)，境界条件フラグ
std::vector<int> aTri;  // index of triangles，三角形の頂点インデックス
std::vector<int> aQuad; // index of 4 vertices required for bending，曲げ計算のための４頂点の配列
const double lambda = 0.0; // Lame's 1st parameter，ラメ第一定数
const double myu    =30.0; // Lame's 2nd parameter，ラメ第二定数
const double stiff_bend = 1.0e-2; // bending stiffness 曲げ剛性
const double areal_density = 0.1; // areal density of a cloth, 布の面密度
const double gravity[3] = {0,0.05,-10}; // gravitatinal accereration，重力加速度
double time_step_size = 0.015; // size of time step，時間ステップの大きさ
const double stiff_contact = 300;
const double contact_clearance = 0.01;
int imode_contact; // mode of contacting object
double mass_point; // mass for a point，頂点あたりの質量

// 自己接触のための関数群
int iroot_bvh; // BVH木構造のルートノードのインデックス
std::vector<CNodeBVH> aNodeBVH; // BVHのノードの配列
std::vector<CAABB3D> aBB_BVH; // AABBの配列，サイズはノードの数
CJaggedArray aEdge;

std::vector<double> aNormal; // deformed vertex noamals，変形中の頂点の法線(可視化用)

// data for camera
double view_height = 2.0;
bool is_animation;
double camera_scale;
double camera_qrot[4];
double camera_trans[3] = {0,0,0};
int imodifier;
int ibutton;
double mov_begin_x, mov_begin_y;
int imode_draw = 0;
/* ------------------------------------------------------------------------ */


void MakeNormal()
{ // make normal
  const int np = (int)aXYZ.size()/3;
  const int ntri = (int)aTri.size()/3;
  aNormal.assign(np*3,0);
  for(int itri=0;itri<ntri;itri++){
    const int ip0 = aTri[itri*3+0];
    const int ip1 = aTri[itri*3+1];
    const int ip2 = aTri[itri*3+2];
    double c[3][3] = {
      { aXYZ[ip0*3+0],aXYZ[ip0*3+1],aXYZ[ip0*3+2] },
      { aXYZ[ip1*3+0],aXYZ[ip1*3+1],aXYZ[ip1*3+2] },
      { aXYZ[ip2*3+0],aXYZ[ip2*3+1],aXYZ[ip2*3+2] } };
    double n[3],area; UnitNormalAreaTri3D(n, area, c[0], c[1], c[2]);
    aNormal[ip0*3+0] += n[0];  aNormal[ip0*3+1] += n[1];  aNormal[ip0*3+2] += n[2];
    aNormal[ip1*3+0] += n[0];  aNormal[ip1*3+1] += n[1];  aNormal[ip1*3+2] += n[2];
    aNormal[ip2*3+0] += n[0];  aNormal[ip2*3+1] += n[1];  aNormal[ip2*3+2] += n[2];
  }
  for(unsigned int ip=0;ip<np;ip++){
    double sqlen =
    + aNormal[ip*3+0]*aNormal[ip*3+0]
    + aNormal[ip*3+1]*aNormal[ip*3+1]
    + aNormal[ip*3+2]*aNormal[ip*3+2];
    double invlen = 1.0/sqrt(sqlen);
    aNormal[ip*3+0] *= invlen;
    aNormal[ip*3+1] *= invlen;
    aNormal[ip*3+2] *= invlen;
  }
}

void StepTime()
{
  
  void (*penetrationDepth)(double& , double* , const double*);
  if(      imode_contact == 0 ){ // contact with plane，床との衝突
    penetrationDepth = penetrationDepth_Plane;
  }
  if(      imode_contact == 1 ){ // contact with sphere，球との衝突
    penetrationDepth = penetrationDepth_Sphere;
  }
  
  std::vector<double> aXYZ1 = aXYZ;
  ::StepTime_InternalDynamics
  (aXYZ, aUVW,
   aXYZ0, aBCFlag,
   aTri, aQuad,
   time_step_size,
   lambda, myu, stiff_bend,
   gravity, mass_point,
   stiff_contact,contact_clearance,penetrationDepth);
  ////
  bool is_impulse_applied;
  GetIntermidiateVelocityContactResolved
  (aUVW,
   is_impulse_applied,
   time_step_size,
   contact_clearance,
   mass_point,
   stiff_contact,
   aXYZ1,
   aTri,
   aEdge,
   iroot_bvh,  aNodeBVH, aBB_BVH);
  if( is_impulse_applied ){
    std::cout << "update middle velocity" << std::endl;
    for(unsigned int ip=0;ip<aXYZ.size()/3;ip++){
      aXYZ[ip*3+0] =  aXYZ1[ip*3+0] + aUVW[ip*3+0]*time_step_size;
      aXYZ[ip*3+1] =  aXYZ1[ip*3+1] + aUVW[ip*3+1]*time_step_size;
      aXYZ[ip*3+2] =  aXYZ1[ip*3+2] + aUVW[ip*3+2]*time_step_size;
    }
    // now aXYZ is collision free
    ::UpdateIntermidiateVelocity
    (aUVW,
     aXYZ, aXYZ0, aBCFlag,
     aTri, aQuad,
     time_step_size,
     lambda, myu, stiff_bend,
     gravity, mass_point);
  }
  MakeNormal();
}

void myGlutDisplay(void)
{
  //	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
	::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);
  
	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );
  
	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
  ::glTranslated(camera_trans[0],camera_trans[1],camera_trans[2]);      
  {
    double R_view_3d[16];
    QuatRot(R_view_3d, camera_qrot);        
    ::glMultMatrixd(R_view_3d);
  }
  
  bool is_lighting = glIsEnabled(GL_LIGHTING);
  
  
  { // draw triangle
    ::glDisable(GL_LIGHTING);    
    ::glLineWidth(2);
    ::glColor3d(0,0,0);
    ::glBegin(GL_LINES);
    for(int itri=0;itri<aTri.size()/3;itri++){
      const int ip0 = aTri[itri*3+0];
      const int ip1 = aTri[itri*3+1];
      const int ip2 = aTri[itri*3+2];
      ::glVertex3d(aXYZ[ip0*3+0],aXYZ[ip0*3+1],aXYZ[ip0*3+2]);
      ::glVertex3d(aXYZ[ip1*3+0],aXYZ[ip1*3+1],aXYZ[ip1*3+2]);
      ::glVertex3d(aXYZ[ip1*3+0],aXYZ[ip1*3+1],aXYZ[ip1*3+2]);
      ::glVertex3d(aXYZ[ip2*3+0],aXYZ[ip2*3+1],aXYZ[ip2*3+2]);
      ::glVertex3d(aXYZ[ip2*3+0],aXYZ[ip2*3+1],aXYZ[ip2*3+2]);
      ::glVertex3d(aXYZ[ip0*3+0],aXYZ[ip0*3+1],aXYZ[ip0*3+2]);
    }
    ::glEnd();
  }  
  { // draw triangle face
    ::glEnable(GL_LIGHTING);
    float color[4] = {200.0/256.0, 200.0/256.0, 200.0/256.0,1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
    glShadeModel(GL_SMOOTH);
    ::glBegin(GL_TRIANGLES);
    for(int itri=0;itri<aTri.size()/3;itri++){
      const int ip0 = aTri[itri*3+0];
      const int ip1 = aTri[itri*3+1];
      const int ip2 = aTri[itri*3+2];
      double c[3][3] = {
        { aXYZ[ip0*3+0],aXYZ[ip0*3+1],aXYZ[ip0*3+2] },
        { aXYZ[ip1*3+0],aXYZ[ip1*3+1],aXYZ[ip1*3+2] },
        { aXYZ[ip2*3+0],aXYZ[ip2*3+1],aXYZ[ip2*3+2] } };
      ::glNormal3d(aNormal[ip0*3+0],aNormal[ip0*3+1],aNormal[ip0*3+2]);
      ::glVertex3dv(c[0]);
      ::glNormal3d(aNormal[ip1*3+0],aNormal[ip1*3+1],aNormal[ip1*3+2]);
      ::glVertex3dv(c[1]);
      ::glNormal3d(aNormal[ip2*3+0],aNormal[ip2*3+1],aNormal[ip2*3+2]);
      ::glVertex3dv(c[2]);
    }
    ::glEnd();
  }

  
  { // fixed boundary condition
    ::glDisable(GL_LIGHTING);    
    ::glPointSize(5);
    ::glColor3d(0,0,1);
    ::glBegin(GL_POINTS);
    for(int ip=0;ip<aXYZ.size()/3;ip++){
      if( aBCFlag[ip] == 0 ) continue;
      ::glVertex3d(aXYZ0[ip*3+0],aXYZ0[ip*3+1],aXYZ0[ip*3+2]);
    }
    ::glEnd();
  }
  
  if(      imode_contact == 0 ){     // draw floor
    ::glDisable(GL_LIGHTING);    
    ::glLineWidth(1);
    ::glColor3d(1,0,0);
    ::glBegin(GL_LINES);
    double grid_x_min = -10;
    double grid_x_max = +10;
    double grid_y_min = -10;
    double grid_y_max = +10;
    int ndiv_grid = 30;
    for(unsigned int ix=0;ix<ndiv_grid+1;ix++){
      double x0 = (grid_x_max-grid_x_min) / ndiv_grid * ix + grid_x_min;
      ::glVertex3d(x0,grid_y_min,-0.5);
      ::glVertex3d(x0,grid_y_max,-0.5);
    }
    for(unsigned int iz=0;iz<ndiv_grid+1;iz++){
      double z0 = (grid_y_max-grid_y_min) / ndiv_grid * iz + grid_y_min;
      ::glVertex3d(grid_x_min,z0,-0.5);
      ::glVertex3d(grid_x_max,z0,-0.5);
    }
    ::glEnd();        
  }
  else if( imode_contact == 1 ){
    ::glDisable(GL_LIGHTING);    
    ::glLineWidth(1);
    ::glColor3d(1,0,0);    
    ::glPushMatrix();
    ::glTranslated(0.5, 0.0, 2.0);
    ::glutWireSphere(3, 16, 16);
    ::glPopMatrix();
  }
  
  
  if( is_lighting ){ ::glEnable(GL_LIGHTING); }
  else{              ::glDisable(GL_LIGHTING); }  
  
  ::glutSwapBuffers();
}

void myGlutIdle(){
    
  if( is_animation ){ StepTime(); }
  
  ::glutPostRedisplay();
}


void myGlutResize(int w, int h)
{
  int win_w, win_h;
  if( w < 0 ){
    ::GLint viewport[4];
    ::glGetIntegerv(GL_VIEWPORT,viewport);
    win_w = viewport[2];
    win_h = viewport[3];
  }
  else{
    win_w = w;
    win_h = h;
    ::glViewport(0,0,w,h);
  }
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
  double asp = (double)win_w/win_h;
  ::glOrtho(-asp*view_height*camera_scale,
            +asp*view_height*camera_scale,
            -view_height*camera_scale,
            +view_height*camera_scale,
            -3*view_height,
            +3*view_height);
	::glutPostRedisplay();
}

void myGlutSpecial(int Key, int x, int y)
{
	switch(Key)
	{
    case GLUT_KEY_PAGE_UP:
      camera_scale *= 1.1;
      break;
    case GLUT_KEY_PAGE_DOWN:
      camera_scale *= 1.0/1.1;
      break;
    case GLUT_KEY_HOME :
      break;
    case GLUT_KEY_END :
      break;
    default:
      break;
  }
  myGlutResize(-1, -1);
	::glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
  const double asp = (double)win_w/win_h;
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
  double dx = mov_end_x - mov_begin_x;
  double dy = mov_end_y - mov_begin_y;  
  if( ibutton == GLUT_LEFT_BUTTON && imodifier == GLUT_ACTIVE_CTRL ){
    // rotate view
    double a = sqrt(dx * dx + dy * dy);
    double ar = a*0.5; // angle
    double dq[4] = { cos(ar), -dy*sin(ar)/a, dx*sin(ar)/a, 0.0 };
    if (a != 0.0) {
      double qtmp[4]; QuatMult(qtmp, dq, camera_qrot);
      double R_view_3d[16];
      QuatRot(R_view_3d, qtmp);
      QuatCopy(camera_qrot,qtmp);
    }
  }
  if( ibutton == GLUT_LEFT_BUTTON && imodifier == GLUT_ACTIVE_SHIFT ){
    double v[3] = {
      dx*view_height*camera_scale*asp,
      dy*view_height*camera_scale,
      0};
    camera_trans[0] += v[0];
    camera_trans[1] += v[1];
    camera_trans[2] += v[2];
  }
	mov_begin_x = mov_end_x;
	mov_begin_y = mov_end_y;
	::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
  imodifier = glutGetModifiers();
  ibutton = button;
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	mov_begin_x = (2.0*x-win_w)/win_w;
	mov_begin_y = (win_h-2.0*y)/win_h;
}

void myGlutKeyboard(unsigned char Key, int x, int y)
{
	switch(Key)
	{
    case 'q':
    case 'Q':
    case '\033':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
    case 'a':
      is_animation = !is_animation;
      break;
    case 'd': // change draw mode
      imode_draw++;
      if( imode_draw >= 2 ){
        imode_draw = 0;
      }
      break;
    case ' ':
      imode_contact++;
      aXYZ = aXYZ0;
      aUVW.assign(aUVW.size(),0.0);
      if( imode_contact >= 2 ){
        imode_contact = 0;
      }
      break;
  }
	::glutPostRedisplay();
}

void InitializeGL()
{
  
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);
  
  { // initialize light parameter
    GLfloat light0_Kd[]   = {0.9f, 0.3f, 0.3f, 1.0f};
    GLfloat light0_Pos[4] = {+0.5f, -0.5f, +1.0f, 0.0f};
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light0_Kd);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_Pos);
    
    GLfloat light1_Kd[]   = {0.3f, 0.3f, 0.9f, 1.0f};
    GLfloat light1_Pos[4] = {-0.5f, +0.5f, +1.0f, 0.0f};
    glLightfv(GL_LIGHT1, GL_DIFFUSE,  light1_Kd);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_Pos);
  }
  //////
  // initialize camera parameters
  ::glMatrixMode(GL_MODELVIEW);
  ::glLoadIdentity();
  camera_scale = 0.5;
  myGlutResize(-1,-1);
  { // initialize view rotation
    double a[3] = {1,0,0.3};
    double invla = 1.0/Length3D(a);
    a[0] *= invla;  a[1] *= invla;  a[2] *= invla;
    double theta = -0.5;
    camera_qrot[0]=cos(theta*0.5);
    camera_qrot[1]=a[0]*sin(theta*0.5);
    camera_qrot[2]=a[1]*sin(theta*0.5);
    camera_qrot[3]=a[2]*sin(theta*0.5);
  }
}


int main(int argc,char* argv[])
{
  { // initialze data
    double total_area;
    SetClothShape_Square(aXYZ0,aBCFlag,aTri,aQuad,total_area,
                         elem_length,cloth_size_x,cloth_size_z);
    mass_point = total_area*areal_density / (aXYZ0.size()/3.0);
    // initialize deformation
    aXYZ = aXYZ0;
    aUVW.assign(aXYZ.size(),0.0);
    MakeNormal();    
    ////
    iroot_bvh = MakeBVHTopology_TopDown(aTri,aXYZ,aNodeBVH);
    aEdge.SetEdgeOfElem(aTri,(int)aTri.size()/3,3,(int)aXYZ.size()/3,false);
  }
  
  
  glutInit(&argc, argv);
  
	// Initialize GLUT window 3D
  glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
 	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  glutCreateWindow("3D View");
	glutDisplayFunc(myGlutDisplay);
	glutIdleFunc(myGlutIdle);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
  
  ////////////////////////
  
  InitializeGL();
 
  glutMainLoop();
	return 0;
}


