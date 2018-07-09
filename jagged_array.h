//
//  jagged array.h
//
//  Created by Nobuyuki Umetani on 11/7/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//


#if !defined(JAGGED_ARRAY_H)
#define JAGGED_ARRAY_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <vector>
#include <assert.h>
#include <iostream>

//! Jagged array
class CJaggedArray
{
public:
	CJaggedArray(){}

  void InitializeSize(int size){
    index.clear();
    index.resize(size+1,0);
    array.clear();
  }
  int Size() const{
    if( index.size() == 0 ){ return 0; }
    return (int)index.size()-1;
  }

	//! initialize with transpose pattern
	void SetTranspose(int size, const CJaggedArray& crs){	
    index.clear();
		index.resize( size+1, 0 );
		for(int icrs=0;icrs<crs.array.size();icrs++){
			const int jno = crs.array[icrs];
			if( jno < size ){ index[jno+1]++; }
		}
		for(int j=0;j<size;j++){ index[j+1] += index[j]; }
		array.resize(index[size]);
		for(int i=0;i<crs.Size();i++){
		for(int icrs=crs.index[i];icrs<crs.index[i+1];icrs++){
			const int jno = crs.array[icrs];
			const int jcrs = index[jno];
			array[jcrs] = i;
			index[jno]++;
		}
		}
		for(int k=size-1;k>=0;k--){ index[k+1] = index[k]; }
		index[0] = 0;
	}
public:
  //! initialize as a dense matrix
	void Fill(int ncol, int nrow){
		index.resize(ncol+1);
		index[0] = 0;
		array.resize(ncol*nrow);
		for(int i=0;i<ncol;i++){
			index[i+1] = (i+1)*nrow;
			for(int j=0;j<nrow;j++){
				array[i*nrow+j] = j;
			}
		}
	}
	//! sort that each sub array have incremental odering
	void Sort(){
    if( index.size() == 0 ) return;
    const int size = (int)index.size()-1;
		for(int ipoin=0;ipoin<size;ipoin++){
			const int is = index[ipoin  ];
			const int ie = index[ipoin+1];
			if( is == ie ) continue;
			assert( is < ie );
			int itmp;
			for(int i=is;i<ie-1;i++){
				for(int j=ie-1;j>i;j--){
					if( array[j] < array[j-1] ){
						itmp = array[j];
						array[j] = array[j-1];
						array[j-1] = itmp;
					}
				}
			}
		}
	}
	bool CheckValid() const {
    int size = (int)index.size()-1;
    if( index[size] > (int)array.size() ){ return false; }
		return true;
	}
  void SetNodeToElem(const std::vector<int>& elem,
                     int nelem,
                     int nnoel,
                     int nnode)
  {
    index.clear();
    index.resize(nnode+1,0);
    for(int ielem=0;ielem<nelem;ielem++){
      for(int inoel=0;inoel<nnoel;inoel++){
        int ino0 = elem[ielem*nnoel+inoel];
        index[ino0+1]++;
      }
    }
    for(int ino=0;ino<nnode;ino++){
      index[ino+1] += index[ino];
    }
    const int nelsup = index[nnode];
    array.resize(nelsup);
    for(int ielem=0;ielem<nelem;ielem++){
      for(int inoel=0;inoel<nnoel;inoel++){
        int ino0 = elem[ielem*nnoel+inoel];
        const int ind = index[ino0];
        array[ind] = ielem;
        index[ino0]++;
      }
    }
    for(int inode=nnode;inode>0;inode--){
      index[inode] = index[inode-1];
    }
    index[0] = 0;
  }
  void SetNodeToElem(const int* paElem,
                     int nelem,
                     int nnoel,
                     int nnode)
  {
    index.clear();
    index.resize(nnode+1,0);
    for(int ielem=0;ielem<nelem;ielem++){
      for(int inoel=0;inoel<nnoel;inoel++){
        int ino0 = paElem[ielem*nnoel+inoel];
        index[ino0+1]++;
      }
    }
    for(int ino=0;ino<nnode;ino++){
      index[ino+1] += index[ino];
    }
    const int nelsup = index[nnode];
    array.resize(nelsup);
    for(int ielem=0;ielem<nelem;ielem++){
      for(int inoel=0;inoel<nnoel;inoel++){
        int ino0 = paElem[ielem*nnoel+inoel];
        const int ind = index[ino0];
        array[ind] = ielem;
        index[ino0]++;
      }
    }
    for(int inode=nnode;inode>0;inode--){
      index[inode] = index[inode-1];
    }
    index[0] = 0;
  }
  void SetEdgeOfElem(const std::vector<int>& elem,
                     int nelem,
                     int nnoel,
                     int nnode,
                     bool is_dia)
  {
    CJaggedArray crs;
    crs.SetNodeToElem(elem,nelem,nnoel,nnode);
    std::vector<int> aflg;
    aflg.resize(nnode,-1);
    index.resize(nnode+1,0);
    for(int inode=0;inode<nnode;inode++){
      if( !is_dia ){ aflg[inode] = inode; }
      for(int icrs=crs.index[inode];icrs<crs.index[inode+1];icrs++){
        int jelem = crs.array[icrs];
        for(int jnoel=0;jnoel<nnoel;jnoel++){
          int jnode = elem[jelem*nnoel+jnoel];
          if( aflg[jnode] != inode ){
            aflg[jnode] = inode;            
            index[inode+1]++;
          }
        }
      }
    }
    for(int ino=0;ino<nnode;ino++){
      index[ino+1] += index[ino];
    }
    const int nelsup = index[nnode];
    array.resize(nelsup);
    for(int ino=0;ino<nnode;ino++){ aflg[ino] = -1; }
    for(int inode=0;inode<nnode;inode++){
      if( !is_dia ){ aflg[inode] = inode; }      
      for(int icrs=crs.index[inode];icrs<crs.index[inode+1];icrs++){
        int jelem = crs.array[icrs];
        for(int jnoel=0;jnoel<nnoel;jnoel++){
          int jnode = elem[jelem*nnoel+jnoel];
          if( aflg[jnode] != inode ){
            aflg[jnode] = inode;            
            const int ind = index[inode];
            array[ind] = jnode;
            index[inode]++;
          }
        }
      }
    }
    for(int inode=nnode;inode>0;inode--){
      index[inode] = index[inode-1];
    }
    index[0] = 0;
  }
  
  void SetEdgeOfElem(const int* paElem,
                     int nelem,
                     int nnoel,
                     int nnode,
                     bool is_dia)
  {
    CJaggedArray crs;
    crs.SetNodeToElem(paElem,nelem,nnoel,nnode);
    std::vector<int> aflg;
    aflg.resize(nnode,-1);
    index.resize(nnode+1,0);
    for(int inode=0;inode<nnode;inode++){
      if( !is_dia ){ aflg[inode] = inode; }
      for(int icrs=crs.index[inode];icrs<crs.index[inode+1];icrs++){
        int jelem = crs.array[icrs];
        for(int jnoel=0;jnoel<nnoel;jnoel++){
          int jnode = paElem[jelem*nnoel+jnoel];
          if( aflg[jnode] != inode ){
            aflg[jnode] = inode;
            index[inode+1]++;
          }
        }
      }
    }
    for(int ino=0;ino<nnode;ino++){
      index[ino+1] += index[ino];
    }
    const int nelsup = index[nnode];
    array.resize(nelsup);
    for(int ino=0;ino<nnode;ino++){ aflg[ino] = -1; }
    for(int inode=0;inode<nnode;inode++){
      if( !is_dia ){ aflg[inode] = inode; }
      for(int icrs=crs.index[inode];icrs<crs.index[inode+1];icrs++){
        int jelem = crs.array[icrs];
        for(int jnoel=0;jnoel<nnoel;jnoel++){
          int jnode = paElem[jelem*nnoel+jnoel];
          if( aflg[jnode] != inode ){
            aflg[jnode] = inode;
            const int ind = index[inode];
            array[ind] = jnode;
            index[inode]++;
          }
        }
      }
    }
    for(int inode=nnode;inode>0;inode--){
      index[inode] = index[inode-1];
    }
    index[0] = 0;
  }
public:
	std::vector<int> index;
	std::vector<int> array;
};

#endif
