import Orange
from AZutilities import evalUtilities
from AZutilities import dataUtilities
from trainingMethods import AZorngRF
from trainingMethods import AZorngCvBoost
from trainingMethods import AZorngCvSVM
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


def getCM(pred, actual, CM):
    if actual == "Inactive":
        if pred == "Inactive":
            CM[1][1] = CM[1][1] + 1
        else:
            CM[1][0] = CM[1][0] + 1
    elif actual == "Active":
        if pred == "Active":
            CM[0][0] = CM[0][0] + 1
        else:
            CM[0][1] = CM[0][1] + 1
    return CM


def predict(model, test, fps, fpsTest, label, thrs):

    # Predict test set using AD
    fid = open("predictions_"+label+".txt", "w")
    fid.write("Pred\tProb\tActual\tCorrect\n")
    fid.close()
    CM = [[0, 0], [0, 0]]
    outAD = 0
    for idx in range(len(test)):
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
    
        # Calculate the median of the topological distance to the 10 NN in the train set
        dist = getDist(fpsTest[idx], fps)
    
        if dist > thrs:
            CM = getCM(pred, actual, CM)
            print "pred, prob, actual, correct  ", pred, prob, actual, correct
            fid = open("predictions_"+label+".txt", "a")
            fid.write(pred+"\t"+str(prob)+"\t"+actual+"\t"+str(correct)+"\n")
            fid.close()
        else:
            outAD = outAD + 1
    
    print CM
    MCC = round(evalUtilities.calcMCC(CM),3)
    print "MCC of test set ", MCC
    print "Fraction of outAD in test set ", float(outAD)/len(test)
    return MCC, float(outAD)/len(test)
    


train = dataUtilities.DataTable("IIDsetAZOdesc.txt")
extTest = dataUtilities.DataTable("nonIIDextTestAZOdesc.txt")
randTest = dataUtilities.DataTable("IIDextTestAZOdesc.txt")
print "Train set ", len(train)
print "randTest set ", len(randTest)
print "extTest set ", len(extTest)

# Calculate fingerprints for train and test sets
fps = getFps(train)
fpsRandTest = getFps(randTest)
fpsExtTest = getFps(extTest)

# Deselect descriptors with no variance
descList = ["ID", "Smiles", "Conc", "Effect", "Conc_1",  "Effect_1", "ID_1", "origSmiles_1", "BioActivity_1" \
,"rdk.fr_dihydropyridine", "rdk.fr_nitroso", "rdk.fr_benzodiazepine", "rdk.fr_thiocyan", "rdk.VSA_EState4" ,"rdk.VSA_EState6" \
,"rdk.VSA_EState7" ,"rdk.VSA_EState1" ,"rdk.VSA_EState2" ,"rdk.VSA_EState3" ,"rdk.SlogP_VSA9" ,"rdk.SMR_VSA8" ,"rdk.fr_diazo" \
,"rdk.fr_prisulfonamd" ,"rdk.fr_isocyan" ,"rdk.fr_azide" ,"rdk.fr_isothiocyan"]
train = dataUtilities.attributeDeselectionData(train, descList)
print "Length domain ", len(train.domain)

#learner = AZorngCvSVM.CvSVMLearner(C=32, gamma=0.03125)
learner = AZorngRF.RFLearner()
#learner = AZorngRF.RFLearner(stratify = "Yes") # No effect
#learner = AZorngCvBoost.CvBoostLearner()
#learner.stratify = "Yes" # No effect
#learner.priors = {"Active":0.80, "Inactive":0.20}
model = learner(train)

thrs = 0.75
print "Predictions for random test set"
MCC_rand, outAD_rand = predict(model, randTest, fps, fpsRandTest, "extIID_0.75", thrs)
print "Predictions for external test set"
MCC_ext, outAD_ext = predict(model, extTest, fps, fpsExtTest, "extNonIID_0.75", thrs)









