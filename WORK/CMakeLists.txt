
cmake_minimum_required(VERSION 3.0)

project(external-library Fortran)

enable_language(Fortran)


# find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(MPI REQUIRED)


if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

if(Apple)
endif(Apple)
# intel imac
# include_directories(/usr/local/Cellar/hdf5-mpi/1.12.2/include)
# link_directories(/usr/local/Cellar/hdf5-mpi/1.12.2/lib)

# m1 mac
# include_directories(/opt/homebrew/Cellar/hdf5-mpi/1.13.0/include)
# link_directories(/opt/homebrew/Cellar/hdf5-mpi/1.13.0/lib)
include_directories(/opt/homebrew/Cellar/hdf5-mpi/1.14.1/include)
link_directories(/opt/homebrew/Cellar/hdf5-mpi/1.14.1/lib)

if(LINUX)
    # skywalker
    include_directories(/home/skywalker/spack/opt/spack/linux-ubuntu21.10-zen2/gcc-11.2.0/hdf5-1.12.2-mhxavduujpiuoqw7cat6nmckqvz35k7n/include)
    link_directories(/home/skywalker/spack/opt/spack/linux-ubuntu21.10-zen2/gcc-11.2.0/hdf5-1.12.2-mhxavduujpiuoqw7cat6nmckqvz35k7n/lib)

    # hufflepuff
    include_directories(/home/hashimoto/spack/opt/spack/linux-centos8-zen2/gcc-10.2.0/hdf5-1.10.7-ckxmqvh3sbside7v3wpm54fm4665otid/include)
    link_directories(/home/hashimoto/spack/opt/spack/linux-centos8-zen2/gcc-10.2.0/hdf5-1.10.7-ckxmqvh3sbside7v3wpm54fm4665otid/lib)
endif(LINUX)

set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch")


set(PROJECT_SOURCE_DIR "../SRC")

set(EXECUTABLE_OUTPUT_PATH "..")
set(EXECUTABLE MPI-OpenMPI)

add_executable(${EXECUTABLE}
  ${PROJECT_SOURCE_DIR}/cklib.f
  ${PROJECT_SOURCE_DIR}/thermo.f90
  ${PROJECT_SOURCE_DIR}/main.f90
  ${PROJECT_SOURCE_DIR}/flux.f90
  ${PROJECT_SOURCE_DIR}/params.f90
  ${PROJECT_SOURCE_DIR}/variables.f90
  ${PROJECT_SOURCE_DIR}/source.f90
  ${PROJECT_SOURCE_DIR}/mach.f
  ${PROJECT_SOURCE_DIR}/math.f
  ${PROJECT_SOURCE_DIR}/tranlib.f
  ${PROJECT_SOURCE_DIR}/trans.f90
  ${PROJECT_SOURCE_DIR}/xerror.f
  ${PROJECT_SOURCE_DIR}/odes.f90
  ${PROJECT_SOURCE_DIR}/muscl.f90
  ${PROJECT_SOURCE_DIR}/weno.f90
  ${PROJECT_SOURCE_DIR}/mfr.f90
  ${PROJECT_SOURCE_DIR}/radiation.f90
  ${PROJECT_SOURCE_DIR}/ignitionsource.f90
  ${PROJECT_SOURCE_DIR}/boundary.f90
  ${PROJECT_SOURCE_DIR}/flow_bound.f90
  ${PROJECT_SOURCE_DIR}/update.f90
  ${PROJECT_SOURCE_DIR}/timeintegration.f90
  ${PROJECT_SOURCE_DIR}/write_hdf5.f90
  ${PROJECT_SOURCE_DIR}/read_hdf5.f90
  ${PROJECT_SOURCE_DIR}/parallel.f90
  ${PROJECT_SOURCE_DIR}/eos.f90
  ${PROJECT_SOURCE_DIR}/block.f90
  ${PROJECT_SOURCE_DIR}/mkdir.f90
  )

target_link_libraries(${EXECUTABLE}
  PRIVATE
  MPI::MPI_Fortran
  hdf5_hl
  hdf5_fortran
  hdf5
  )

