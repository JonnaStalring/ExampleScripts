import Orange
from AZutilities import evalUtilities
from AZutilities import dataUtilities
from AZutilities import AZOrangePredictor
from trainingMethods import AZorngRF
import orngTest
import orngStat
import orange
import string
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem
import numpy
from rdkit.Chem import AllChem



def getDist(queryFp, fps):

    distList = []
    for fp in fps:
        dist = DataStructs.FingerprintSimilarity(queryFp,fp) # Tanimoto
        #dist = DataStructs.DiceSimilarity(queryFp,fp) # Dice
        distList.append(dist)
    distList.sort()
    #medDist = numpy.median(distList[len(distList)-9:len(distList)])
    medDist = numpy.median(distList[len(distList)-4:len(distList)])
    return medDist



def getFps(data):

    molList = []
    for ex in data:
        mol = Chem.MolFromSmiles(ex["Smiles"].value)
        if mol:
            molList.append(mol)
        else:
            print ex["Smiles"].value
            print ex["Leonumber"].value
    fps = [FingerprintMols.FingerprintMol(x) for x in molList] # Topological
    #fps = [AllChem.GetMorganFingerprint(x, 2) for x in molList]
    #print "Length of data and fp ", len(data), len(fps)
    return fps

THRS = 0.75

model = AZorngRF.RFread("OI_RFmodel")
predictor = AZOrangePredictor.AZOrangePredictor("OI_RFmodel")

train = dataUtilities.DataTable("Assay1350AZOdesc.txt")

test = dataUtilities.DataTable("nonIIDextTestAZOdesc.txt")

# Calculate fingerprints for train and test sets
fps = getFps(train)
fpsTest = getFps(test)

for idx in range(len(test[0:2])):

     smiles = test[idx]["Smiles"].value
     predictor.getDescriptors(smiles)
     prediction, prob = predictor.predict(True)
     # Normalize prob
     prob = 200*abs(prob)
     if prob > 100:
         prob = 0.98
     #print "****************** Prediction with AZOrangePredictor ******************"
     print prediction
     print prob

     # Create an Orange data set to calculate the fp of the smiles
     features = [Orange.data.variable.String("Smiles")]
     domain = Orange.data.Domain(features)
     smiData = Orange.data.Table(domain)
     smiData.append([smiles])
     fpsSmiles = getFps(smiData) 
     distSmi = getDist(fpsSmiles[0], fps)
     print "Dist from smiles ", distSmi

     predList = model(test[idx], returnDFV = True)
     pred = predList[0].value
     prob = predList[1]
     # Normalize prob
     prob = 200*abs(prob)
     if prob > 100:
         prob = 0.98
     actual = test[idx]["BioActivity"].value
     if pred == actual:
         correct = True
     else:
         correct = False
     print "ID, pred, prob, actual, correct  ",test[idx]["ID"].value, pred, prob, actual, correct

     # Calculate the median of the topological distance to the 10 NN in the train set
     dist = getDist(fpsTest[idx], fps)
     print "Dist from test ", dist

     #if dist > THRS:
     #    print "ID, pred, prob, actual, correct  ",test[idx]["ID"].value, pred, prob, actual, correct
     #else:
     #    outAD = outAD + 1
























