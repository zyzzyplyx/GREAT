BASE_DIR = /home/chuck/github-repos/GREAT
KENT_DIR = ${BASE_DIR}/kent/src
KENT_LIBS =  $(KENT_DIR)/lib/$(MACHTYPE)/jkweb.a

# summarize libs and includes needed, note that libs are order dependent
INCS = -I$(KENT_DIR)/inc
LIBS = $(KENT_LIBS) -lm

CC=gcc
CFLAGS= ${COPT} ${INCS} -c -Wall -static
LDFLAGS=

RDOBJECTS = createRegulatoryDomains.o regdom.o
POBJECTS = calculateBinomialP.o regdom.o
BETAPOBJECTS = calculateBetaP.o regdom.o

all: createRegulatoryDomains calculateBinomialP calculateBetaP

createRegulatoryDomains: $(RDOBJECTS)
	$(CC) $(LDFLAGS) ${COPT} -o $@ $(RDOBJECTS) ${LIBS}

calculateBinomialP: $(POBJECTS)
	$(CC) $(LDFLAGS) ${COPT} -o $@ $(POBJECTS) ${LIBS}
	
calculateBetaP: $(BETAPOBJECTS)
	$(CC) $(LDFLAGS) ${COPT} -o $@ $(BETAPOBJECTS) ${LIBS}

clean:
	rm -f *.o createRegulatoryDomains calculateBinomialP calculateBetaP

regdom.o:	regdom.h
createRegulatoryDomains.o:	regdom.h
calculateBinomialP.o:	regdom.h
calculateBetaP.o: regdom.h
