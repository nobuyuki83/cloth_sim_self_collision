//
//  matrix_sparse.h
//
//  Created by Nobuyuki Umetani on 11/7/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//


#if !defined(MATRIX_SQUARE_SPARSE_H)
#define MATRIX_SQUARE_SPARSE_H

#include <vector>

class CMatrixSquareSparse
{
public:
	CMatrixSquareSparse(int nblk, int len);
  
	CMatrixSquareSparse();
	virtual ~CMatrixSquareSparse();
  
  void Initialize(int nblk, int len);
  void operator = (const CMatrixSquareSparse& m);
  void SetPattern(const std::vector<int>& colind, const std::vector<int>& rowptr);

	bool SetZero();
	bool Mearge(int nblkel_col, const int* blkel_col,
              int nblkel_row, const int* blkel_row,
              int blksize, const double* emat,
              std::vector<int>& m_marge_tmp_buffer);
  // Calc Matrix Vector Product
  // {y} = alpha * [A]{x} + beta * {y}  
	void MatVec(double alpha,
              const std::vector<double>& x,
              double beta,
              std::vector<double>& y) const;
  void SetBoundaryCondition(const std::vector<int>& bc_flag);
public:
	int m_nblk;
	int m_len;
  
	int  m_ncrs;
	int* m_colInd;
	int* m_rowPtr;
  
	double* m_valCrs;
	double* m_valDia;
};

double InnerProduct
(std::vector<double>& r_vec,
 std::vector<double>& u_vec);

void AXPY(double a,
          const std::vector<double>& x,
          std::vector<double>& y);

void Solve_CG
(double& conv_ratio,
 int& iteration,
 const CMatrixSquareSparse& mat,
 std::vector<double>& r_vec,
 std::vector<double>& u_vec);


#endif // MATDIA_CRS_H
