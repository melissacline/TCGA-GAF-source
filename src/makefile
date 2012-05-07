include ~/kent/src/inc/common.mk

KENT = /cluster/home/cline/kent
L += -lm ${MYSQLLIBS}
MYLIBDIR = ${KENT}/src/lib/${MACHTYPE}
MYLIBS =  ${MYLIBDIR}/jkhgap.a ${MYLIBDIR}/jkweb.a

CFLAGS += -I${KENT}/src/inc
CFLAGS += -I${KENT}/src/hg/inc
O = transcriptsToGenesAndExons.o

transcriptsToGenesAndExons: $O 
	${CC} ${COPT} -o transcriptsToGenesAndExons $O ${MYLIBS} $L

transcriptsToGenesAndExons.o: transcriptsToGenesAndExons.c
	${CC} ${COPT} ${CFLAGS} -o $@ -c $<

clean::
	rm -f ${O} transcriptsToGenesAndExons
