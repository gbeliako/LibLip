mex -output mlliblip matlabliblip.c src/*.cpp glpk-4.10/src/*.c -Iglpk-4.10/include -Isrc CFLAGS="\$CFLAGS -fpermissive " -v
