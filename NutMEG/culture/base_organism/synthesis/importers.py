""" functions for importing data """

from BioMolecule import BioMolecule
import os
this_dir = os.path.dirname(__file__)


def copies_moles(filename='EColiSynth_Ishihama2008_dataset3.csv'):

    dataset = open(this_dir+'/../../../data/synthesis/'+filename, "r")

    CopyNo = []
    MolMass = []

    params_input = dataset.readlines()
    firstline=True # which line number of data we're on
    for line in params_input:
        if not line.startswith("#") and len(line) > 0:   # Kill comments and blanks, though there shouldn't be any anyway
            token = line.split("\t") #split via commas as its a csv file
            if firstline:
                firstline=False
                continue
            #elif str(token[31]) == 'no': # use to limit to cytosolic proteins (by removing robosomic ones.
            else:
				## now we have all of the relevant data
                CopyNo.append(float(token[8]))
                MolMass.append(float(token[10]))

    dataset.close()

    return CopyNo, MolMass


def AAreader(host_drymass, host_volume, filename='EColiSynth_AminoAcids.csv'):

    dataset = open(this_dir+'/../../../data/synthesis/'+filename, "r")

    AA = []

    params_input = dataset.readlines()
    firstline=True # which line number of data we're on
    for line in params_input:
        if not line.startswith("#") and len(line) > 0:   # Kill comments and blanks, though there shouldn't be any anyway
            token = line.split("\t") #split via tabs as that's how I saved it
            if firstline:
                firstline=False
                continue
            else:
                ## add in the AA name, its frequency, molecular weight, and concentration (g/mol)
                Bio = BioMolecule(str(token[0]), float(token[5])) #name and molecular weight (also, Gibbs data)
                Bio.frequency = float(token[6]) # frequency
                Bio.conc_mol_per_cell = float(token[3])*host_drymass*1000 # data is in mol/(g cell) so convert
                # for newconc, if implemented, use float[7] below
                Bio.conc_mol_per_l = Bio.conc_mol_per_cell/(host_volume*1000) # convert to mol/l (volume is in SI units)
                # print(str(token[0]), Bio.conc_mol_per_cell/(host_volume*1000), float(token[7]))
                AA.append(Bio)
                # the final 0 is for protein link counting

    dataset.close()


    return AA
