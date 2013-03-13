import scipy
from GREATx import WeightedDart, TermDartTSSTriple

inFile = open('../data/SRFtoTerms.data')
outFile = open('../data/GiLocal.data','w')

lineObjects = []
for line in inFile:
    lineObjects.append(TermDartTSSTriple(line))
termIDs = list(set([lineObject.termID for lineObject in lineObjects])) ##unique set of term IDS
dartNames = list(set([lineObject.dartName for lineObject in lineObjects])) ##unique set of dart names

for termID in termIDs:
    if termID != 'UNKNOWN':
        termIDObjects = filter(lambda x: x.termID == termID,lineObjects)##All the lines associated with the term
        weightedDarts = []
        for dartName in dartNames:
            termIDObjectsForDartName =\
                    filter(lambda x: x.dartName == dartName, termIDObjects)
            if len(termIDObjectsForDartName) > 0:
                dartWeight = sum([x.weight for x in termIDObjectsForDartName])
                dartPosition = termIDObjectsForDartName[0].dartPosition
                chrName = termIDObjectsForDartName[0].chrName
                weightedDarts.append(WeightedDart(chrName=chrName,\
                                                  name=dartName,\
                                                  position=dartPosition,\
                                                  weight=dartWeight))

        X_bar = sum([x.weight for x in weightedDarts])/len(weightedDarts)
        S = scipy.sqrt(sum([x.weight**2/len(weightedDarts)\
                            for x in weightedDarts]) - X_bar**2)

        for weightedDartObjectI in weightedDarts: 
            numerator1 = 0.0
            numerator2 = 0.0
            denom1 = 0.0
            denom2 = 0.0
            for weightedDartObjectJ in weightedDarts:
                if weightedDartObjectJ.name == weightedDartObjectI.name:
                    continue
                if weightedDartObjectJ.chrName != weightedDartObjectI.chrName:
                    continue

                if (scipy.fabs(weightedDartObjectI.position - weightedDartObjectJ.position) > 46944323): # size of the smallest chromosome
                    spatialWeight = 0
                else:
                    spatialWeight = 1

                numerator1 += spatialWeight*weightedDartObjectJ.weight
                numerator2 += spatialWeight*X_bar
                denom1 += len(weightedDarts)*(spatialWeight**2)            
                denom2 += spatialWeight            

            denom2 = denom2**2
            numerator = numerator1 - numerator2
            if(len(weightedDarts) == 1): 
                ZScore = "Only one dart on " + weightedDartObjectI.chrName
            else:
                denominator = S*scipy.sqrt(\
                        (denom1 - denom2)/(len(weightedDarts) - 1))
                if denominator == 0:
                    ZScore = "0 denom"
                else:
                    ZScore = numerator/denominator
            outFile.write(\
                    "\t".join([termID,\
                               weightedDartObjectI.name,\
                               weightedDartObjectI.chrName,\
                               str(weightedDartObjectI.position),\
                               str(ZScore)])\
                    + "\n")

