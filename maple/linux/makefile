# C++ compiler variables
CC = c++
CCC = cc
LL = c++

INCLUDE_PATH= -I/usr/include/ -I/usr/include/malloc/ -Itnt/  -Iglpk/  -I/opt/maple2017/extern/include/

CPPFLAGS= -O2 -Wno-deprecated -fPIC
LINKFLAGS= -O
LINK_LIBS= -lm 
OBJS=	forest.o \
	interpol.o \
	liblip.o \
	slipint.o \
	mwrap_dliblip.o

FLAGS =  -O -Wno-deprecated -fPIC -g -m64 
LDFLAGS=  -lc -lmaplec  -lrt -lmaple  -shared


C_FILES := $(wildcard glpk/*.c)
OBJ_FILES := $(addprefix obj/,$(notdir $(C_FILES:.c=.o)))

LIB_PATH = /opt/maple2017/bin.X86_64_LINUX/

all:	mapleliblip.so

mapleliblip.so: $(OBJS) $(OBJ_FILES)
		$(LL)  $(OBJS) $(OBJ_FILES)  -o $@ $(LDFLAGS) -L$(LIB_PATH) -v

forest.o:	forest.cpp
		$(CC) -c forest.cpp $(CPPFLAGS) $(INCLUDE_PATH)
		
interpol.o:	interpol.cpp
		$(CC) -c interpol.cpp $(CPPFLAGS) -I$(INCLUDE_PATH)
		
liblip.o:	liblip.cpp
		$(CC) -c liblip.cpp $(CPPFLAGS) -I$(INCLUDE_PATH)
		
slipint.o:	slipint.cpp
		$(CC) -c slipint.cpp $(CPPFLAGS) -I$(INCLUDE_PATH)
		
mwrap_dliblip.o:			 mwrap_dliblip.c
				$(CCC) -c  mwrap_dliblip.c $(FLAGS) -I$(INCLUDE_PATH)		

obj/%.o:	glpk/%.c
	$(CCC)  -c $(CPPFLAGS) -o $@ $<
   
   
.PHONY: clean all

clean:
	rm -f *.o
	rm -f *.so
	rm -f obj/*.o

	
