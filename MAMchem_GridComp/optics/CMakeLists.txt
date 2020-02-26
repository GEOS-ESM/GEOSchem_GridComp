set(miev_srcs 
   miev/ErrPack.f
   miev/MIEV0.F
   )

add_library(miev ${miev_srcs})
target_link_libraries(miev PUBLIC OpenMP::OpenMP_Fortran)

string(REPLACE " " ";" tmp ${FREAL8})
foreach(flag ${tmp})
   target_compile_options (miev PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach()

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include(UseF2Py)

add_f2py_module(optics_
   SOURCES mie.F90 
   DESTINATION lib/Python/optics
   ONLY scattering_lognormal
   LIBRARIES miev
   INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR}
   USE_OPENMP DOUBLE_PRECISION
   )
add_dependencies(optics_ miev)

install(
   FILES radiation.py gads.py mam7-optics.lut.py
   DESTINATION lib/Python/optics
   )
