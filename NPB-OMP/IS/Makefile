SHELL=/bin/sh
BENCHMARK=is
BENCHMARKU=IS

include ../config/make.def

include ../sys/make.common

OBJS = is.o \
       ${COMMON}/c_print_results.o \
       ${COMMON}/c_timers.o \
       ${COMMON}/c_wtime.o


${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB}

is.o:             is.cpp  npbparams.hpp
	${CCOMPILE} is.cpp

clean:
	- rm -f *.o *~ mputil*
	- rm -f npbparams.hpp core
	- if [ -d rii_files ]; then rm -r rii_files; fi
