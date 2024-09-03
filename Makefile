### Compilers r###
CC=g++## for C++
FF=gfortran## for Fortran

###  Directories ###
SRCDIR=${PWD}/src
INCDIR=${PWD}/inc
LHAPDFDIR=$//home/yehudi/LHAPDF6/include
LHAPDFDIR2=$//home/yehudi/LHAPDF6/share/LHAPDF
LIBDIR=${PWD}/lib
BLDDIR=${PWD}/bld
GSLDIR=$//home/yehudi/gsl/include
BOOSTDIR=$//usr/local/boost_1_77_0

### Flags ###
CFLAG=-O3 -Wall -I${INCDIR} -I${LHAPDFDIR} -I${LHAPDFDIR2} -I${GSLDIR} -I${BOOSTDIR}
FFLAG=-O3 -I${INCDIR}

### Paths ###
VPATH=${SRCDIR}

### Source files ###
CFILES=$(wildcard ${SRCDIR}/*.cpp)
FFILES=$(wildcard ${SRCDIR}/*.f)

### Object files ###
COBJS=$(subst .cpp,.o,$(subst ${SRCDIR},${BLDDIR},${CFILES}))
FOBJS=$(subst .f,.o,$(subst ${SRCDIR},${BLDDIR},${FFILES}))

### Libraries ###
LIB=${LIBDIR}/libresum.a
GSLIB=-L/home/yehudi/gsl/lib -lgsl -lgslcblas
STLIB= -L/usr/local/gfortran/lib -lm -lgfortran 
LHAPDFLIB= -L/home/yehudi/LHAPDF6/lib -lLHAPDF

### Commands ###
all: RUN
lib: ${LIB}

RUN: main.cpp ${LIB}
	${CC} ${CFLAG} -o $@ main.cpp ${LIB} ${STLIB} ${GSLIB} ${LHAPDFLIB}

${LIB}:	${COBJS} ${FOBJS}
	ar -ruc $@ ${BLDDIR}/*.o

${BLDDIR}%.o: ${SRCDIR}%.f
	cd ${BLDDIR}; ${FF} -c ${FFLAG} $?;

${BLDDIR}%.o: ${SRCDIR}%.cpp
	cd ${BLDDIR}; ${CC} -c ${CFLAG} $<;
	
### Cleaning ###
clean:
	rm -f RUN ${LIBDIR}/*.a ${BLDDIR}/*.o *.log output/*; clear;

