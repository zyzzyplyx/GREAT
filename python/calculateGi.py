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
            relevantDarts = filter(lambda x: x.dartName == dartName,termIDObjects)
            if len(relevantDarts) > 0:
                dartWeight = sum([x.weight for x in relevantDarts])
                dartPosition = relevantDarts[0].dartPosition
                chrName = relevantDarts[0].chrName
                weightedDarts.append(WeightedDart(chrName = chrName, name = dartName, position = dartPosition, weight = dartWeight))

        X_bar = sum([x.weight for x in weightedDarts])/len(weightedDarts)
        S = scipy.sqrt( sum([x.weight**2/len(weightedDarts) for x in weightedDarts]) - X_bar**2 )
	for weightedDartObjecti in weightedDarts: 
            numerator1 = 0.0
            numerator2 = 0.0
            denom1 = 0.0
            denom2 = 0.0
            for weightedDartObjectj in weightedDarts:
                if(weightedDartObjectj.name == weightedDartObjecti.name): continue
                if weightedDartObjectj.chrName != weightedDartObjecti.chrName:
                    spatialWeight = 0.0
                else:
                    spatialWeight = float(scipy.fabs(weightedDartObjecti.position - weightedDartObjectj.position))
                    numerator1 = numerator1 + spatialWeight*weightedDartObjectj.weight
                    numerator2 = numerator2 + spatialWeight*X_bar
                    denom1 = denom1 + len(weightedDarts)*(spatialWeight**2)			
                    denom2 = denom2 + spatialWeight			
            denom2 = denom2**2
            numerator = numerator1 - numerator2
            if(len(weightedDarts) == 1): 
                Zscore = "Only one dart on " + weightedDartObjecti.chrName
            else:
                denominator = S*scipy.sqrt((denom1 -denom2)/(len(weightedDarts) - 1))
                if denominator == 0:
                    Zscore = "O denom"
                else:
                    Zscore = numerator/denominator
            outFile.write(termID + "\t" + weightedDartObjecti.name +"\t" +str(Zscore) + "\n")

            
                
		
            


