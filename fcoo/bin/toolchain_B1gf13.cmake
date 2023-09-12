# File toolchain_B1gf13.cmake for gfortran-13
# override DEBUG and RELEASE Fortran flags:
  set(CMAKE_Fortran_FLAGS_RELEASE "-g -O3 -fno-second-underscore -ffree-line-length-none -fallow-argument-mismatch" CACHE STRING "")
# The last option is to accept some mismatch with MPI_SEND/RECV/WAITALL/STARTALL/BCAST
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0  -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow -frecursive -fbacktrace -march=native -fno-second-underscore -ffree-line-length-none -fallow-argument-mismatch -Wno-unused-label" CACHE STRING "")
# The last option here to avoid that some unused format lines cause errors with -Wall (-Wunused)

