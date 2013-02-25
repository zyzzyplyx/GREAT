KENT_DIR = path/to/your/kent/src
KENT_LIBS =  $(KENT_DIR)/lib/$(MACHTYPE)/jkweb.a

# summarize libs and includes needed, note that libs are order dependent
INCS = -I$(KENT_DIR)/inc
LIBS = $(KENT_LIBS) -lm

CC=gcc
CFLAGS= ${COPT} ${INCS} -c -Wall -static
LDFLAGS=

RDOBJECTS = createRegulatoryDomains.o regdom.o
POBJECTS = calculateBinomialP.o regdom.o

all: createRegulatoryDomains calculateBinomialP

createRegulatoryDomains: $(RDOBJECTS)
	$(CC) $(LDFLAGS) ${COPT} -o $@ $(RDOBJECTS) ${LIBS}

calculateBinomialP: $(POBJECTS)
	$(CC) $(LDFLAGS) ${COPT} -o $@ $(POBJECTS) ${LIBS}

clean:
	rm -f *.o createRegulatoryDomains calculateBinomialP

regdom.o:	regdom.h
createRegulatoryDomains.o:	regdom.h
calculateBinomialP.o:	regdom.h

