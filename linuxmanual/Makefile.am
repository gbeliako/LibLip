lib_LIBRARIES = liblip.a
liblip_a_SOURCES = forest.cpp interpol.cpp liblip.cpp slipint.cpp forest.h interpol.h memblock.h liblipc.h slipint.h liblip.h
liblip_a_FLAGS = -version-info 2:0:0

INCLUDES = -Iglpk/ -Iqhull/ -Itnt/

CXXFLAGS = -Wno-deprecated

liblip_a_LIBADD= -lglpk -lqhull -lm

liblip_a_CC = g++

nobase_include_HEADERS = memblock.h forest.h interpol.h slipint.h liblip.h  liblipc.h 
