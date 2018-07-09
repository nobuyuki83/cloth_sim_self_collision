project(internal_cloth_sparse)

cmake_minimum_required(VERSION 2.8)
set( CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g" )

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
)

add_executable(${PROJECT_NAME}
  main.cpp
  ../cloth_internal_physics.h  
  ../ilu_sparse.cpp                   
  ../ilu_sparse.h
  ../jagged_array.h
  ../matrix_square_sparse.cpp         
  ../matrix_square_sparse.h
  ../solve_internal_sparse.h
  ../utility.h
  ../vector3d.h
)

target_link_libraries(${PROJECT_NAME} 
  ${GLUT_LIBRARY} 
  ${OPENGL_LIBRARY}
)
