# File toolchain_B1if19.cmake for ifort-19
# override DEBUG and RELEASE Fortran flags:
  set(CMAKE_Fortran_FLAGS_RELEASE "-g -O3 -no-fma -ip -traceback -i4 -real-size 32 -fp-model precise -assume byterecl -fno-alias -fno-fnalias -sox -xhost" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -no-fma -ip -traceback -i4 -real-size 32 -fp-model precise -assume byterecl -fno-alias -fno-fnalias -sox -debug all -warn all -check all -check noarg_temp_created -fp-stack-check -heap-arrays -fpe0 -xhost" CACHE STRING "")

