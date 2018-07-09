//
//  self_collision_cloth.h
//
//  Created by Nobuyuki Umetani on 10/28/12.
//  Copyright (c) 2012 Nobuyuki Umetani. All rights reserved.
//

#ifndef contact_self_collision_cloth_h
#define contact_self_collision_cloth_h

#include <vector>
#include <set>

#include "jagged_array.h"
#include "aabb.h"
#include "bvh_aabb.h"



// 衝突が解消された中間速度を返す
void GetIntermidiateVelocityContactResolved
(std::vector<double>& aUVWm,
 bool& is_impulse_applied,
 ////
 double dt,
 double contact_clearance,
 double mass_point,
 double cloth_contact_stiffness,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aTri,
 const CJaggedArray& aEdge,
 int iroot_bvh,
 const std::vector<CNodeBVH>& aNodeBVH,
 std::vector<CAABB3D>& aBB);
    
#endif
