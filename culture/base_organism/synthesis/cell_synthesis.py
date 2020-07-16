import sys,os, math, statistics
sys.path.append(os.path.dirname(__file__))
import importers
from BioMolecule import BioMolecule

from NutMEG.util.sqlshortcuts import sqlshortcuts as sqlsc

import sqlite3

this_dir = os.path.dirname(__file__)
class cell_synthesis:

    @staticmethod
    def getreplistAA(AAList, listlen):
        """Get a representative list of Amino acids to build proteins
        which matches average cell composition."""
        replist = []
        for a in AAList:
            su =round(a.frequency*listlen) #vnumber of times to show up
            for i in range(0,su):
                a.counter=1
                replist.append(a)
        return replist



    # for calculating the quotient. Make sure concentrations are sent in mol/l
    @staticmethod
    def lnProteinReactionQuotient(PC, AAList):
        Q = math.log(PC)
        for a in AAList:
            # multiply by aa count each time
            if a.counter != 0 and a.conc_mol_per_l != 0.:
                # print(a.name, a.counter, math.log(a.conc_mol_per_l))
                Q -= a.counter*(math.log(a.conc_mol_per_l)) # again, in moles/cell --- g/mol /MM, then to mol/L...
        return Q




    @staticmethod
    def ProteinSynthCost(TK, host_drymass, host_volume):
        """Calculate and return the cost of protein synthesis per gram.
        Assumes amino acids already synthesised.

        Parameters
        -----------
        TK : float, int
            Temperature in Kelvin
        host_drymass : float
            dry mass of organism in kg/cell
        host_volume : float
            total (wet) volume of one organism.
        """
        # get temperature as an int to save unneccesary reruns
        TK =round(TK)

        # get arrays of number of copies and molecular mass for each protein
        # Ishihama (2008) has data for
        CopyNo, MolMass = importers.copies_moles()

        # extract data about 20 amino acids, get a list of BioMolecules
        # representing them.
        AAList = importers.AAreader(host_drymass, host_volume)
        #backbones and glycine groups
        AABB = BioMolecule('AABB', 74.05866, T=TK)
        PBB = BioMolecule('PBB', 56.04918, T=TK)
        GLY = BioMolecule('GLY', 57.05892, T=TK)
        Water = BioMolecule('H2O(l)', 18, T=TK) # not a BioMolecule but this
        #class will get us all the info we need.

        # average amino acid molecular mass.
        avgAAMM = 108.0197 # u


        # estimate the length of each protein in the dataset
        Length = []
        for m, c in zip(MolMass, CopyNo):
            Length.append(m/avgAAMM)


        sumlens = 0. # total length (number of amino acids)
        sumcops = 0. # total copies (number of proteins)
        summass = 0. # total protein mass
        for m, c, l in zip(MolMass, CopyNo, Length):
            sumlens += (c*l)
            sumcops += c
            summass += (c*m)
        meanlength = int(sumlens/sumcops) #Mean protein length in dataset
        PMass_cell = summass #total protein massin dataset in u, (For EColi, suggested to be c. 32% of total)



        # build the mean protein in the dataset.
        # built proportionally with available amino acids.
        # uses composition of E Coli

        ProteinGibbs = AABB.std_formation_gibbs #std fomration gibbs for the whole protein
        #start with the amino acid backbone, update as we add more

        ProteinMass = AABB.Mr # whole protein mass. update as we add more.

        ReactionStdGibbs = 0. #reaction std energy cost. update as we add more.

        #representative list of Amino acids to build protein which matches average cell composition.
        replist = cell_synthesis.getreplistAA(AAList, meanlength)


        # build the protein
        i=0 #counter
        for AA in replist:
        # while i<a:

            # AA = getRandomAA(AAList)

            # input it to the overall reaction calculation (AA is a reagent, so take away)
            ReactionStdGibbs -= AA.std_formation_gibbs


            #Glycine has special requirements...
            if AA.name == 'Glycine(aq)':
                # add GLY and no PBB
                ProteinGibbs += GLY.std_formation_gibbs
                ProteinMass += GLY.Mr
            else:
                if i!=0:
                    ProteinGibbs += PBB.std_formation_gibbs # add a protein backbone
                    ProteinMass += PBB.Mr # add the backbone mass
                    # add the formed water molecule to the overall reaction calculation
                    ReactionStdGibbs += Water.std_formation_gibbs
                ProteinGibbs += AA.std_formation_R # add the R group Gibbs
                ProteinMass += (AA.Mr - AABB.Mr) # add the R group mass


            i+=1

        # we now have the mass of this protein [u] and its standard free energy of formation [J/mol].

	    # now complete the Gibbs of the reaction sum(AA) -> protein + water.
        # we added the AA and water as we built the protein
        ReactionStdGibbs += ProteinGibbs

        #data covers c. 32% of cell protein so scale up to get the number of proteins per cell.
        ProteinConc_NoPerCell = meanlength*sumcops*3.

		# Note: the protein concentration could actually be less/more than this,
        #      as we assume the final conc of an individual cell.

        # convert to SI mol/l
        ProteinConc_MolPerL = ProteinConc_NoPerCell/(6.022e23 * host_volume*1000)

        # get the ln reaction quotient which is a damn pain
        lnQ = cell_synthesis.lnProteinReactionQuotient(ProteinConc_MolPerL, replist)

        # now get the Gibbs of reaction
        ReactionGibbs = ReactionStdGibbs+(8.31*(TK)*lnQ)

        TotalEReq = ReactionGibbs/ProteinMass
        # energy required to build 1 gram of the average protein.
        # dimensional anlaysis: J/mol * (1/u) = J/g
        # Assuming peptide costs scale with the whole cell,
        # this is the approx cost of building macromolecules from AA in 1 dry g of cells.

        return TotalEReq




    @staticmethod
    def a_e(T, Eh=-0.27):#-0.27):
        """ Calculate the activity of an electron
        Assumes neutral conditons for Eh"""
        return(math.exp(-96485.332*Eh/(8.314462618*T)))

    @staticmethod
    def AAsynthCost(TK):
        """Get the cost of synthesising all of the amino acids in one dry g of
        cells using the method from McCollom and Amend (2005).
        This is currently fixed to E Coli.

        #TODO: Export the list of lists to be in a csv file in /data/, then
        import it here. This will allow for other model organisms to be
        simulated.
        #TODO: Work in methods for oxic environments
        """
        TK = round(TK)
        AAinfo = {}
        # build a dictionary of the table in McCollom and Amend (2005)
        for n, a, g, HCO3, H, NH4, e, H2S, umol in zip(
          ['Alanine', 'Arginine', 'Asparginine','Aspartate', 'Cysteine',
            'Glutamate', 'Glutamine', 'Glycine', 'Histidine', 'Isoleucine',
            'Leucine', 'Lysine', 'Methionine', 'Phenylalanine', 'Proline',
            'Serine', 'Threonine', 'Tryptophan', 'Tyrosine', 'Valine'],
          [48.4,49.0,30.3,30.5,10.5,40.8,36.5,43.7,14.0,36.2,56.1,47.7,21.8,
            29.1,24.2,21.5,28.7,11.0,23.7,27.1],
          [191.7,248.2,153.1,169.8,128.2,294.6,281.1,76,237.6,537.1,546.4,
            503.3,544.1,774.8,377.2,101.4,208.9,850.0,714.9,426.4],
          [3,6,4,4,3,5,5,2,6,6,6,6,5,9,5,3,4,11,9,5],
          [14,25,14,14,12,21,21,7,23,35,35,33,26,48,26,12,9,55,46,28],
          [1,4,2,1,1,1,2,1,3,1,1,2,1,1,1,1,4,2,1,1],
          [12,22,12,12,10,18,18,6,20,30,30,28,22,40,22,10,9,46,38,24],
          [0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
          [543,281,229,229,87,278,250,582,90,276,428,326,146,176,210,205,241,
            54,131,402]):
            # a in mg/(gcells)
            # std g in kJ/mol formation
            # umol is per dry g
            AAinfo[n] = {'Amount':a, 'stdGr':-g*1000, 'HCO3':HCO3, 'H':H, 'NH4':NH4, 'e':e, 'H2S':H2S, 'umol':umol}


        SynthCost=0 # in J
        for key, v in AAinfo.items():
            # exponents denominator, e.g. log10(denominator of Q)
            # there's 10^-2 M HCO3, 1o^-7 H+ (pH 7), 10^-8 NH4 M, 10^-5 H2S M
            #Again, from McCollom and Amend (I think)
            Qdenom_exp = ((-2)*v['HCO3'])+((-7)*v['H'])+((-8)*v['NH4'])+((-5)*v['H2S'])
            l = cell_synthesis.a_e(TK)**v['e']+Qdenom_exp #electon activity
            # overall log10 reaction quootient
            log10Q = -12-(math.log10(cell_synthesis.a_e(TK)**v['e'])+Qdenom_exp)
            # get G but don't forget to convert log to ln
            G = v['stdGr']+(8.31*(TK)*log10Q*math.log(10))
            SynthCost += (G*v['umol']*1e-6) # J per dry g
        return SynthCost



    @staticmethod
    def get_ESynth_density(T, compute=None, AA=False):
        """Calculate the energy cost to synthesise one dry g of cells.

        If you want to compute ESynth density for a new organism,
        pass compute as a dictionary in the following format:
        {'dbfilename':None (or path/to_database),
        'host':base_organism}
        """
        AAsynth, Psynth = None, None

        if compute != None:
            # first see if a new db filename has been passed
            try:
                AAsynth, Psynth, = cell_synthesis.extract_AAP_from_db(T=T,
                  filename=compute['dbfilename'])
            except:
                AAsynth = cell_synthesis.AAsynthCost(T)
                Psynth  = cell_synthesis.ProteinSynthCost(
                  T, compute['host'].dry_mass, compute['host'].base_volume)
                if type(compute['dbfilename']) is type('str'):
                    # add results to the database
                    cell_synthesis.add_AAP_to_db(
                      T, AAsynth, Psynth, filename=compute['dbfilename'])
        else:
            # use E Coli
            try:
                AAsynth, Psynth, = cell_synthesis.extract_AAP_from_db(
                  T=round(T), filename='EColiSynth') #csv might be easier
            except:
                AAsynth = cell_synthesis.AAsynthCost(T)
                Psynth  = cell_synthesis.ProteinSynthCost(T, 2.8e-16, 1e-18)
                cell_synthesis.add_AAP_to_db(T, AAsynth, Psynth)
        if AA:
            # return costs for AA synthesis and protein synthesis together.
            # left this in so previous versions stay compatable
            return AAsynth + Psynth
        else:
            # return all costs as a dictionary, leaving room to add more.
            return {'Psynth' : Psynth, 'AAsynth':AAsynth}


    @staticmethod
    def add_AAP_to_db(T, AAsynth, Psynth, filename='EColiSynth'):
        """Add results to an SQL database to save having to recompute them.

        Parameters
        __________
        T : Temperature in K
        AAsynth : cost of amino acid synthesis in J per dry g cells
        Psynth : cost of 1 dry g protein synthesis from amino acids

        #TODO: Perhaps would be more intuitive as a csv rather than db.
        """
        db = sqlite3.connect(this_dir+'/../../../data/synthesis/'+filename)
        cursor = db.cursor()
        try:
            # attempt to add in an entry, with OrgID 'Tester'
            Paramsdict = {'Temp':['REAL', round(T)],
              'costAASynth':['REAL', AAsynth],
              'costPeptideSynth':['REAL', Psynth]}

            cursor.execute('CREATE TABLE IF NOT EXISTS ' +\
            'synthesis (Temp REAL, costAAsynth REAL, costPeptideSynth REAL)')

            command, vals = sqlsc.INSERTlst('synthesis', Paramsdict)
            cursor.execute(command, vals)
            db.commit()
        except sqlite3.Error or TypeError as e:
            raise e
        except:
            raise
        finally:
            db.close()


    @staticmethod
    def extract_AAP_from_db(T, filename='EColiSynth'):
        """Extract amino acid and protein synthesis from database"""

        db=sqlite3.connect(this_dir+'/../../../data/synthesis/'+filename)
        cursor = db.cursor()
        AA = 0.0
        P = 0.0
        try:
            cursor.execute('SELECT costAASynth FROM synthesis WHERE Temp = '+str(round(T)))
            AA = cursor.fetchone()[0]
            cursor.execute('SELECT costPeptideSynth FROM synthesis WHERE Temp = '+str(round(T)))
            P = cursor.fetchone()[0]
        except sqlite3.Error or TypeError as e:
            raise e
        finally:
            db.close()
        return AA, P


    # TODO LATER
    # distance this from E Coli
