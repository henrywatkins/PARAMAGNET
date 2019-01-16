#==============================================================
#
#        MAKEFILE FOR MAGNETISED COLLISIONAL TRANSPORT CODE
#
#        The .o files are compiled separately then linked in
#        the final make MAGNET command.
#
#
#
#===============================================================


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

CC= mpicc#${PETSC_DIR}/${PETSC_ARCH}/bin/mpicxx
CFLAGS= -fopenmp -O3
LDFLAGS= -fopenmp
OBJECTS= MAGNET.o Global.o TransportCoeff.o Output.o Initial.o 
INCLUDES= -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
LIBS= 
#COMPILE=mpicc

MAGNET: ${OBJECTS}
	${CC} ${LDFLAGS} ${OBJECTS}  ${PETSC_LIB} -o MAGNET
	rm *.o

MAGNET.o: MAGNET.cpp  Global.h TransportCoeff.h Output.h Initial.h
	${CC} -c ${CFLAGS} MAGNET.cpp ${INCLUDES}

Global.o: Global.cpp Global.h
	${CC} -c ${CFLAGS} Global.cpp

TransportCoeff.o: TransportCoeff.cpp TransportCoeff.h
	${CC} -c ${CFLAGS} TransportCoeff.cpp

Output.o: Output.cpp Output.h
	${CC} -c ${CFLAGS} Output.cpp

Initial.o: Initial.cpp Initial.h
	${CC} -c ${CFLAGS} Initial.cpp
