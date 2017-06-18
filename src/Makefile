objects=main.o potts.o pottsSqure.o
all:$(objects)
	g++ -o potts.exe $(objects)
main.o:main.cpp
	g++ -c  main.cpp
potts.o:potts.cpp potts.h
	g++ -c potts.cpp
pottsSqure.o:pottsSqure.cpp pottsSqure.h
	g++ -c  pottsSqure.cpp