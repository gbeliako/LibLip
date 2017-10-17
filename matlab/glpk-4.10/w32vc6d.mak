## w32vc6d.mak (W32 multi-threaded dynamic library, VC++ 6.0) ##

##----------------------------------------------------------------------
## To build W32 multi-threaded dynamic library for GLPK with VC++ 6.0
## enter the subdirectory where the file you are reading is placed and
## type the command:
##
##    nmake.exe /f w32vc6d.mak
##
## The following three files will be created in the same subdirectory:
##
##    glpk410.dll, GLPK dynamic library,
##
##    glpk410.lib, GLPK import library, and
##
##    glpsol.exe, stand-alone LP/MIP solver.
##
## To check if GLPK has been successfully compiled type the command:
##
##    nmake.exe /f w32vc6d.mak check
##
## Header files needed to use the GLPK dynamic library are contained in
## the subdirectory 'include'.
##
## Written by Andrew Makhorin <mao@mai2.rcnet.ru>, <mao@gnu.org>.
##----------------------------------------------------------------------

CFLAGS = /nologo /W3 /MT

.c.obj:
        cl.exe $(CFLAGS) /Iinclude /Fo$*.obj /c $*.c

all: glpk410.dll glpsol.exe tspsol.exe

check: glpsol.exe
        glpsol.exe --mps examples\plan.mps

glpk410.dll: src\glpavl.obj \
        src\glpbfi.obj \
        src\glpdmp.obj \
        src\glphbm.obj \
        src\glpiet.obj \
        src\glpinv.obj \
        src\glpios1.obj \
        src\glpios2.obj \
        src\glpios3.obj \
        src\glpipm.obj \
        src\glpipp1.obj \
        src\glpipp2.obj \
        sysdep\w32\glplib1a.obj \
        sysdep\w32\glplib1b.obj \
        src\glplib2.obj \
        src\glplib3.obj \
        src\glplib4.obj \
        src\glplpp1.obj \
        src\glplpp2.obj \
        src\glplpx1.obj \
        src\glplpx2.obj \
        src\glplpx3.obj \
        src\glplpx4.obj \
        src\glplpx5.obj \
        src\glplpx6a.obj \
        src\glplpx6b.obj \
        src\glplpx6c.obj \
        src\glplpx6d.obj \
        src\glplpx7.obj \
        src\glplpx7a.obj \
        src\glplpx7b.obj \
        src\glplpx8a.obj \
        src\glplpx8b.obj \
        src\glplpx8c.obj \
        src\glplpx8d.obj \
        src\glplpx8e.obj \
        src\glpluf.obj \
        src\glpmat.obj \
        src\glpmip1.obj \
        src\glpmip2.obj \
        src\glpmpl1.obj \
        src\glpmpl2.obj \
        src\glpmpl3.obj \
        src\glpmpl4.obj \
        src\glpqmd.obj \
        src\glprng.obj \
        src\glpspx1.obj \
        src\glpspx2.obj \
        src\glpstr.obj \
        src\glptsp.obj
        cl.exe $(CFLAGS) /LD /Feglpk410.dll src\*.obj sysdep\w32\*.obj \
                w32vc6d.def

glpsol.exe: examples\glpsol.obj glpk410.dll
        cl.exe $(CFLAGS) /Feglpsol.exe examples\glpsol.obj glpk410.lib

tspsol.exe: examples\tspsol.obj glpk410.dll
        cl.exe $(CFLAGS) /Fetspsol.exe examples\tspsol.obj glpk410.lib

## eof ##
