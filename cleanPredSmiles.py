import sys
from molvs import standardize_smiles
from molvs import Standardizer
from molvs import tautomer
from molvs import charge
from molvs import fragment
from molvs import normalize
from rdkit import Chem

def getCleanSmiles(inSmiles):
    #smiles = standardize_smiles(inSmiles)
    smiles = inSmiles

    inMol = Chem.MolFromSmiles(smiles)

    t = tautomer.TautomerCanonicalizer()
    outMol = t.canonicalize(inMol)

    f = fragment.LargestFragmentChooser()
    outMol = f.choose(outMol)

    c = charge.Uncharger()
    outMol = c.uncharge(outMol)

    smiles = Chem.MolToSmiles(outMol)

    #print(inSmiles+"\t -> \t"+smiles)
    return smiles


if __name__ == "__main__":
    smiles = sys.argv[1]
    cleanSmiles = getCleanSmiles(smiles)
    print(cleanSmiles)
