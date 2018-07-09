
//
//  matrix square sparse.h
//
//  Created by Nobuyuki Umetani on 11/7/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//


#if defined(__VISUALC__)
#pragma warning( disable : 4786 ) 
#endif
#define for if(0); else for

#include <iostream>
#include <cassert>
#include <math.h>
#include <vector>

#include "matrix_square_sparse.h"


CMatrixSquareSparse::CMatrixSquareSparse()
{
	m_nblk = 0;
	m_len = 0;
  
	m_ncrs = 0;
	m_colInd = 0;
	m_rowPtr = 0;
  
	m_valCrs = 0;
	m_valDia = 0;
}

CMatrixSquareSparse::~CMatrixSquareSparse()
{
  if( m_colInd != 0 ){ delete[] m_colInd; m_colInd = 0; }
	if( m_rowPtr != 0 ){ delete[] m_rowPtr; m_rowPtr = 0; }
  
	if( m_valCrs != 0 ){ delete[] m_valCrs; m_valCrs = 0; }
	if( m_valDia != 0 ){ delete[] m_valDia; m_valDia = 0; }
}


void CMatrixSquareSparse::Initialize(int nblk, int len)
{
  this->m_nblk = nblk;
  this->m_len = len;
  const int blksize = m_len*m_len;
  m_valDia = new double [nblk*blksize];
  for(int i=0;i<nblk*blksize;i++){ m_valDia[i] = 0; }
  m_colInd = new int [nblk+1];
  for(int i=0;i<nblk+1;i++){ m_colInd[i] = 0; }
  ////  
	if( m_rowPtr != 0 ){ delete[] m_rowPtr; m_rowPtr = 0; }
	if( m_valCrs != 0 ){ delete[] m_valCrs; m_valCrs = 0; }  
}

void CMatrixSquareSparse::operator = (const CMatrixSquareSparse& m)
{
  this->m_nblk = m.m_nblk;
  this->m_len  = m.m_len;
  const int blksize = m_len*m_len;
  this->m_ncrs = m.m_ncrs;  
  if( m_colInd != 0 ){ delete[] m_colInd; m_colInd = 0; }
  if( m_rowPtr != 0 ){ delete[] m_rowPtr; m_rowPtr = 0; }
	if( m_valCrs != 0 ){ delete[] m_valCrs; m_valCrs = 0; }
  if( m_valDia != 0 ){ delete[] m_valDia; m_valDia = 0; }
  m_colInd = new int    [m_nblk+1];
  m_rowPtr = new int    [m_ncrs];
  m_valDia = new double [m_nblk*blksize];
  m_valCrs = new double [m_ncrs*blksize];
  for(int i=0;i<m_nblk+1;      i++){ m_colInd[i] = m.m_colInd[i]; }
  for(int i=0;i<m_ncrs;        i++){ m_rowPtr[i] = m.m_rowPtr[i]; }
  for(int i=0;i<m_nblk*blksize;i++){ m_valDia[i] = m.m_valDia[i]; }
  for(int i=0;i<m_ncrs*blksize;i++){ m_valCrs[i] = m.m_valCrs[i]; }
}


bool CMatrixSquareSparse::SetZero()
{
	for(int i=0;i<m_len*m_len*m_nblk;i++){ m_valDia[i] = 0.0; }
  for(int i=0;i<m_len*m_len*m_ncrs;i++){ m_valCrs[i] = 0.0; }
	return true;
}

bool CMatrixSquareSparse::Mearge
(int nblkel_col, const int* blkel_col,
 int nblkel_row, const int* blkel_row,
 int blksize, const double* emat,
 std::vector<int>& marge_buffer)
{
  assert( m_colInd != 0 );
	assert( m_valCrs != 0 );
	assert( m_valDia != 0 );

	assert( nblkel_col == nblkel_row );
  assert( blksize == m_len*m_len );

	const int* colind = m_colInd;
	const int* rowptr = m_rowPtr;
	double* vcrs = m_valCrs;
	double* vdia = m_valDia;

	for(int iblkel=0;iblkel<nblkel_col;iblkel++){
		const int iblk1 = blkel_col[iblkel];
		for(int jpsup=colind[iblk1];jpsup<colind[iblk1+1];jpsup++){
			assert( jpsup < m_ncrs );
			const int jblk1 = rowptr[jpsup];
			marge_buffer[jblk1] = jpsup;
		}
		for(int jblkel=0;jblkel<nblkel_row;jblkel++){
			if( iblkel == jblkel ){	// Marge Diagonal
				const double* pval_in = &emat[(iblkel*nblkel_row+iblkel)*blksize];
				double* pval_out = &vdia[iblk1*blksize];
				for(int i=0;i<blksize;i++){ pval_out[i] += pval_in[i]; }
			}
			else{	// Marge Non-Diagonal
				const int jblk1 = blkel_row[jblkel];
				assert( jblk1 < m_nblk );
				if( marge_buffer[jblk1] == -1 ) continue;
        assert( marge_buffer[jblk1] >= 0 && marge_buffer[jblk1] < (int)m_ncrs );
				const int jpsup1 = marge_buffer[jblk1];
				assert( m_rowPtr[jpsup1] == jblk1 );
				const double* pval_in = &emat[(iblkel*nblkel_row+jblkel)*blksize];
				double* pval_out = &vcrs[jpsup1*blksize];
				for(int i=0;i<blksize;i++){ pval_out[i] += pval_in[i]; }
			}
		}
		for(int jpsup=colind[iblk1];jpsup<colind[iblk1+1];jpsup++){
			assert( jpsup < m_ncrs );
			const int jblk1 = rowptr[jpsup];
			marge_buffer[jblk1] = -1;
		}
	}
	return true;
}

void CMatrixSquareSparse::SetPattern
(const std::vector<int>& colind,
 const std::vector<int>& rowptr)
{
//  CMatrixSquareSparse::DeletePattern();

  assert( m_valDia != 0);
  assert( m_colInd != 0 );
	assert( m_ncrs == 0 );
  assert( m_rowPtr == 0 );
  
  assert( colind.size() == m_nblk+1 );
  for(int iblk=0;iblk<m_nblk+1;iblk++){
    m_colInd[iblk] = colind[iblk];
  }
  m_ncrs = colind[m_nblk];
  m_rowPtr = new int [m_ncrs];
  for(int icrs=0;icrs<m_ncrs;icrs++){
    m_rowPtr[icrs] = rowptr[icrs];
  }
  
  const int blksize = m_len*m_len;
  m_valCrs = new double [m_ncrs*blksize];
}

// Calc Matrix Vector Product
// {y} = alpha*[A]{x} + beta*{y}
void CMatrixSquareSparse::MatVec
(double alpha,
 const std::vector<double>& x,
 double beta,
 std::vector<double>& y) const
{
	const int blksize = m_len*m_len;

	if( m_len == 1 ){
		const double* vcrs  = m_valCrs;
		const double* vdia = m_valDia;
		const int* colind = m_colInd;
		const int* rowptr = m_rowPtr;
		////////////////
		for(int iblk=0;iblk<m_nblk;iblk++){
			double& vy = y[iblk];
			vy *= beta;
			const int colind0 = colind[iblk];
			const int colind1 = colind[iblk+1];
			for(int icrs=colind0;icrs<colind1;icrs++){
				assert( icrs < m_ncrs );
				const int jblk0 = rowptr[icrs];
				assert( jblk0 < m_nblk );
				vy += alpha * vcrs[icrs] * x[jblk0];
			}
			vy += alpha * vdia[iblk] * x[iblk];
		}
	}
	else if( m_len == 2 ){
		const double* vcrs  = m_valCrs;
		const double* vdia = m_valDia;
		const int* colind = m_colInd;
		const int* rowptr = m_rowPtr;
    const int nblk = m_nblk;    
		////////////////
		for(int iblk=0;iblk<nblk;iblk++){
			y[iblk*2+0] *= beta;
			y[iblk*2+1] *= beta;
			const int icrs0 = colind[iblk];
			const int icrs1 = colind[iblk+1];
			for(int icrs=icrs0;icrs<icrs1;icrs++){
				assert( icrs < m_ncrs );
				const int jblk0 = rowptr[icrs];
				assert( jblk0 < nblk );
				y[iblk*2+0] += alpha * ( vcrs[icrs*4  ]*x[jblk0*2+0] + vcrs[icrs*4+1]*x[jblk0*2+1] );
				y[iblk*2+1] += alpha * ( vcrs[icrs*4+2]*x[jblk0*2+0] + vcrs[icrs*4+3]*x[jblk0*2+1] );
			}
			y[iblk*2+0] += alpha * ( vdia[iblk*4  ]*x[iblk*2+0] + vdia[iblk*4+1]*x[iblk*2+1] );
			y[iblk+2*1] += alpha * ( vdia[iblk*4+2]*x[iblk*2+0] + vdia[iblk*4+3]*x[iblk*2+1] );
		}
	}
	else if( m_len == 3 ){
		const double* vcrs  = m_valCrs;
		const double* vdia = m_valDia;
		const int* colind = m_colInd;
		const int* rowptr = m_rowPtr;
		////////////////
		for(int iblk=0;iblk<m_nblk;iblk++){
			y[iblk*3+0] *= beta;
			y[iblk*3+1] *= beta;
			y[iblk*3+2] *= beta;
			const int icrs0 = colind[iblk];
			const int icrs1 = colind[iblk+1];
			for(int icrs=icrs0;icrs<icrs1;icrs++){
				assert( icrs < m_ncrs );
				const int jblk0 = rowptr[icrs];
				assert( jblk0 < m_nblk );
				y[iblk*3+0] += alpha * ( vcrs[icrs*9  ]*x[jblk0*3+0] + vcrs[icrs*9+1]*x[jblk0*3+1] + vcrs[icrs*9+2]*x[jblk0*3+2] );
				y[iblk*3+1] += alpha * ( vcrs[icrs*9+3]*x[jblk0*3+0] + vcrs[icrs*9+4]*x[jblk0*3+1] + vcrs[icrs*9+5]*x[jblk0*3+2] );
				y[iblk*3+2] += alpha * ( vcrs[icrs*9+6]*x[jblk0*3+0] + vcrs[icrs*9+7]*x[jblk0*3+1] + vcrs[icrs*9+8]*x[jblk0*3+2] );
			}
			y[iblk*3+0] += alpha * ( vdia[iblk*9  ]*x[iblk*3+0] + vdia[iblk*9+1]*x[iblk*3+1] + vdia[iblk*9+2]*x[iblk*3+2] );
			y[iblk*3+1] += alpha * ( vdia[iblk*9+3]*x[iblk*3+0] + vdia[iblk*9+4]*x[iblk*3+1] + vdia[iblk*9+5]*x[iblk*3+2] );
			y[iblk*3+2] += alpha * ( vdia[iblk*9+6]*x[iblk*3+0] + vdia[iblk*9+7]*x[iblk*3+1] + vdia[iblk*9+8]*x[iblk*3+2] );
		}
	}
	else{
		const double* vcrs  = m_valCrs;
		const double* vdia = m_valDia;
		const int* colind = m_colInd;
		const int* rowptr = m_rowPtr;
		////////////////
		for(int iblk=0;iblk<m_nblk;iblk++){
			for(int idof=0;idof<m_len;idof++){ y[iblk*m_len+idof] *= beta; }
			const int colind0 = colind[iblk];
			const int colind1 = colind[iblk+1];
			for(int icrs=colind0;icrs<colind1;icrs++){
				assert( icrs < m_ncrs );
				const int jblk0 = rowptr[icrs];
				assert( jblk0 < m_nblk );
				for(int idof=0;idof<m_len;idof++){
				for(int jdof=0;jdof<m_len;jdof++){
					y[iblk*m_len+idof] += alpha * vcrs[icrs*blksize+idof*m_len+jdof] * x[jblk0*m_len+jdof];
				}
				}
			}
			for(int idof=0;idof<m_len;idof++){
			for(int jdof=0;jdof<m_len;jdof++){
				y[iblk*m_len+idof] += alpha * vdia[iblk*blksize+idof*m_len+jdof] * x[iblk*m_len+jdof];
			}
			}
		}
	}
}

void CMatrixSquareSparse::SetBoundaryCondition
(const std::vector<int>& bc_flag)
{  
	const int blksize = m_len*m_len;
	
	for(int iblk=0;iblk<m_nblk;iblk++){
    if( bc_flag[iblk] == 0 ) continue;
		for(int idof=0;idof<m_len;idof++){
			for(int jdof=0;jdof<m_len;jdof++){
        m_valDia[iblk*blksize+idof*m_len+jdof] = 0.0;
      }
			m_valDia[iblk*blksize+idof*m_len+idof] = 1.0;
    }
    for(int icrs=m_colInd[iblk];icrs<m_colInd[iblk+1];icrs++){
      for(int idof=0;idof<m_len;idof++){
        for(int jdof=0;jdof<m_len;jdof++){
          m_valCrs[icrs*blksize+idof*m_len+jdof] = 0.0;
        }
			}
    }
  }
  for(int icrs=0;icrs<m_ncrs;icrs++){
		const int jblk1 = m_rowPtr[icrs];
    if( bc_flag[jblk1] == 0 ) continue;
		for(int idof=0;idof<m_len;idof++){
			for(int jdof=0;jdof<m_len;jdof++){
        m_valCrs[icrs*blksize+idof*m_len+jdof] = 0.0;
      }
		}
	}
}








//////////////////////////////////////////////////////////////////////////

double InnerProduct
(std::vector<double>& r_vec,
 std::vector<double>& u_vec)
{
  const int n = (int)r_vec.size();
  assert( u_vec.size() == n );
  double r = 0.0;
  for(int i=0;i<n;i++){
    r += r_vec[i]*u_vec[i];
  }
  return r;
}

// {y} = {y} + a * {x}
void AXPY(double a,
          const std::vector<double>& x,
          std::vector<double>& y)
{
  const int n = (int)x.size();
  assert( y.size() == n );
  for(int i=0;i<n;i++){
    y[i] += a*x[i];
  }
}

void Solve_CG
(double& conv_ratio,
 int& iteration,
 const CMatrixSquareSparse& mat,
 std::vector<double>& r_vec,
 std::vector<double>& x_vec)
{
	const double conv_ratio_tol = conv_ratio;
	const int mx_iter = iteration;
  
	const int nblk = mat.m_nblk;
  const int len = mat.m_len;
  assert(r_vec.size() == nblk*len);
  const int ndof = nblk*len;
  
  // {x} = 0
  x_vec.resize(ndof);
  for(int i=0;i<ndof;i++){ x_vec[i] = 0; }
  
  double sqnorm_res = InnerProduct(r_vec,r_vec);
  if( sqnorm_res < 1.0e-30 ){
    conv_ratio = 0.0;
    iteration = 0;
    return;
  }
	double inv_sqnorm_res_ini = 1.0 / sqnorm_res;
  
  std::vector<double> Ap_vec(ndof);
  
	// Set Initial Serch Direction
	// {p} = {r}
  std::vector<double>  p_vec = r_vec;
  
	iteration = mx_iter;
	for(int iitr=1;iitr<mx_iter;iitr++){
    
		double alpha;
		{	// alpha = (r,r) / (p,Ap)
			mat.MatVec(1.0,p_vec,0.0,Ap_vec);
			const double pAp = InnerProduct(p_vec,Ap_vec);
			alpha = sqnorm_res / pAp;
		}
    
		// update x
		// {x} = +alpha*{ p} + {x}
		AXPY(alpha,p_vec,x_vec);
    
		// {r} = -alpha*{Ap} + {r}
		AXPY(-alpha,Ap_vec,r_vec);
    
		double sqnorm_res_new = InnerProduct(r_vec,r_vec);
    // Converge Judgement
		if( sqnorm_res_new * inv_sqnorm_res_ini < conv_ratio_tol*conv_ratio_tol ){
      conv_ratio = sqrt( sqnorm_res * inv_sqnorm_res_ini );
      iteration = iitr;
      return;
		}
    
		// beta = (r1,r1) / (r0,r0)
		const double beta = sqnorm_res_new / sqnorm_res;
		sqnorm_res = sqnorm_res_new;
    
    
		// {p} = {r} + beta*{p}
    for(int i=0;i<ndof;i++){ p_vec[i] = r_vec[i] + beta*p_vec[i]; }
	}
  
	return;
}

