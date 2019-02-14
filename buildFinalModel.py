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

train = dataUtilities.DataTable("Assay1350AZOdesc.txt")

descList = ["ID", "Smiles", "Conc", "Effect", "Conc_1",  "Effect_1", "ID_1", "origSmiles_1", "BioActivity_1" \
,"rdk.fr_dihydropyridine", "rdk.fr_nitroso", "rdk.fr_benzodiazepine", "rdk.fr_thiocyan", "rdk.VSA_EState4" ,"rdk.VSA_EState6" \
,"rdk.VSA_EState7" ,"rdk.VSA_EState1" ,"rdk.VSA_EState2" ,"rdk.VSA_EState3" ,"rdk.SlogP_VSA9" ,"rdk.SMR_VSA8" ,"rdk.fr_diazo" \
,"rdk.fr_prisulfonamd" ,"rdk.fr_isocyan" ,"rdk.fr_azide" ,"rdk.fr_isothiocyan"]
train = dataUtilities.attributeDeselectionData(train, descList)
print "Length domain ", len(train.domain)
train.save("finalTrainData.txt")

learner = AZorngRF.RFLearner()

# CV accuracy
res = orngTest.crossValidation([learner], train, strat=orange.MakeRandomIndices.StratifiedIfPossible, folds = 10)
CM = evalUtilities.ConfMat(res)[0]
CA = round(orngStat.CA(res)[0],3)
MCC = round(evalUtilities.calcMCC(CM),3)
# TH, FL, FH, TL
resList = [str(CM[0][0]), str(CM[0][1]), str(CM[1][0]), str(CM[1][1]), str(CA), str(MCC)]
wrtStr = string.join(resList, "\t")
print "CV results"
print wrtStr

# Build model
model = learner(train)

model.write("OI_RFmodel")

