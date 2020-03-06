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
        replist = []
        for a in AAList:
            su =round(a.frequency*listlen) #number of times to show up
            for i in range(0,su):
                a.counter=1
                replist.append(a)
        return replist



    # for calculating the quotient later. Make sure concentrations are sent in mol/l
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
    def ProteinSynthCost(TK, host_drymass, host_volume, cmfile=None, aafile=None):
        TK =round(TK)

        # get arrays of copynumber and molecular mass for each protein Ishihama has data for

        CopyNo, MolMass = importers.copies_moles()

        AAList = importers.AAreader(host_drymass, host_volume)
        #generateMCParams(AAList) # assign on the parameters we need for Monte Carlo
        AABB = BioMolecule('AABB', 74.05866, T=TK)
        PBB = BioMolecule('PBB', 56.04918, T=TK)
        GLY = BioMolecule('GLY', 57.05892, T=TK)
        Water = BioMolecule('H2O(l)', 18, T=TK) #lol

        avgAAMM = 108.0197 # u



        # produce a gaussian kernel density
        # create arrays which contain every single protein, one for mass one for length
        MolMassCount = []
        LengthCount = []
        Length = []
        # print(CopyNo)
        # USE THIS TO JUST GET THE LENGTH
        for m, c in zip(MolMass, CopyNo):
            Length.append(m/avgAAMM)


        # we need to keep count of how much dry weight of protein we have produced
        TotalPMassSim = 0. #total mass of protein simulated in u
        TotalESim = 0. #total energy requirement simulated in J
        """
        PMass_cell = host_proteinfraction*host_drymass/(1.660539e-27) #in u
        """

        sumlens = 0.
        sumcops = 0.
        summass = 0.
        for m, c, l in zip(MolMass, CopyNo, Length):
            sumlens += (c*l)
            sumcops += c
            summass += (c*m)
        meanlength = int(sumlens/sumcops) #Mean protein length in dataset
        PMass_cell = summass #total protein massin dataset in u, (For EColi, suggested to be c. 32% of total)



        # build the mean protein in the dataset.
        # built proportionally with available amino acids.
        # uses composition of E Coli

        ProteinGibbs = AABB.std_formation_gibbs #start with the amino acid backbone., will become std fomration gibbs for the whole protein Update as we add more
        ProteinMass = AABB.Mr # its mass. update as we add more.
        ReactionStdGibbs = 0. #std energy cost. update as we add more.
        replist = cell_synthesis.getreplistAA(AAList, meanlength) #representative list of Amino acids to build one with matches average cell composition.
        # build the protein
        i=0
        for AA in replist:
        # while i<a:

            # AA = getRandomAA(AAList)

            # input it to the overall reaction calculation (AA is a reagent)
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

	    # now complete the Gibbs of the reaction AA -> protein + water, we added the AA and water as we built the protein
        ReactionStdGibbs += ProteinGibbs

        # ProteinConc_NoPerCell = meanlength*len(LengthCount)*3.
        ProteinConc_NoPerCell = meanlength*sumcops*3. #data covers c. 32% of cell protein so scale up to get the number of proteins per cell.


		# the protein concentration could actually be less/more than this, as we assume the final conc of an individual cell.

        ProteinConc_MolPerL = ProteinConc_NoPerCell/(6.022e23 * host_volume*1000) # convert to SI mol/l

        # get the ln reaction quotient which is a damn pain
        lnQ = cell_synthesis.lnProteinReactionQuotient(ProteinConc_MolPerL, replist)


        # now get the Gibbs of reaction
        ReactionGibbs = ReactionStdGibbs+(8.31*(TK)*lnQ)

        #EperProtein += (ReactionGibbs/(6.022e23)) # J/protein

        #TotalEReq = EperProtein/(ProteinMass*1.660539e-24)

        TotalEReq = ReactionGibbs/ProteinMass
        #energy required to build 1 gram of the average protein. Assuming peptide costs scale with the whole cell, this is the approx cost of building macromolecules from AA in 1 dry g of cells.


        """USE THIS.
        As we just have the average protein, with mass ProteinMass. In one ProteinMass there are (1/(ProteinMass*1.660539e-24)) proteins in one g. (e.g. one dry gram). x this by the energy required to build this mean protein and you have the total energy to produce one dry gram of proteins. As we were always scaling up anyway, this is equivalent to 1 dry g of cells for our model. THAT WILL DO. AND ITS EASY TO EXPLAIN. The Msc project can iron out the big kinks, stop wasting time and rewrite this method to that.


        """

        # return as a density (per dry g)
        return TotalEReq




    """ Calculate the activity of an electron
    Assumes neutral conditons for Eh"""
    @staticmethod
    def a_e(T, Eh=-0.27):#-0.27):
        # print(math.exp(-96485.332*Eh/(8.314462618*T)))
        return(math.exp(-96485.332*Eh/(8.314462618*T)))

    @staticmethod
    def AAsynthCost(TK):
        """Get the cost of synthesising all of the amino acids in one dry g of cells using the method from McCollom and Amend (2005). This is currently fixed to E Coli.

        #TODO: Export the dictionary to be in a csv file in /data/, then import it here. This will allow for other model organisms to be simulated.
        #TODO: Work in methods for oxic environments
        """
        TK = round(TK)
        AAinfo = {}
        for n, a, g, HCO3, H, NH4, e, H2S, umol in zip(['Alanine', 'Arginine', 'Asparginine','Aspartate', 'Cysteine', 'Glutamate', 'Glutamine', 'Glycine', 'Histidine', 'Isoleucine', 'Leucine', 'Lysine', 'Methionine', 'Phenylalanine', 'Proline', 'Serine', 'Threonine', 'Tryptophan', 'Tyrosine', 'Valine'], [48.4,49.0,30.3,30.5,10.5,40.8,36.5,43.7,14.0,36.2,56.1,47.7,21.8,29.1,24.2,21.5,28.7,11.0,23.7,27.1], [191.7,248.2,153.1,169.8,128.2,294.6,281.1,76,237.6,537.1,546.4,503.3,544.1,774.8,377.2,101.4,208.9,850.0,714.9,426.4], [3,6,4,4,3,5,5,2,6,6,6,6,5,9,5,3,4,11,9,5],[14,25,14,14,12,21,21,7,23,35,35,33,26,48,26,12,9,55,46,28],[1,4,2,1,1,1,2,1,3,1,1,2,1,1,1,1,4,2,1,1],[12,22,12,12,10,18,18,6,20,30,30,28,22,40,22,10,9,46,38,24],[0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0], [543,281,229,229,87,278,250,582,90,276,428,326,146,176,210,205,241,54,131,402]):
            # a in mg/(gcells)
            # std g in kJ/mol formation
            # umol is per dry g
            AAinfo[n] = {'Amount':a, 'stdGr':-g*1000, 'HCO3':HCO3, 'H':H, 'NH4':NH4, 'e':e, 'H2S':H2S, 'umol':umol}


        SynthCost=0
        for key, v in AAinfo.items():
            Qdenom_exp = ((-2)*v['HCO3'])+((-7)*v['H'])+((-8)*v['NH4'])+((-5)*v['H2S'])
            l = cell_synthesis.a_e(TK)**v['e']+Qdenom_exp
            #Q = 1e-12 / 10**(math.log10(cell_synthesis.a_e(TK)**v['e'])+Qdenom_exp)
            log10Q = -12-(math.log10(cell_synthesis.a_e(TK)**v['e'])+Qdenom_exp)
            # print(10**(math.log10(a_e(273.15+25)**v['e'])+Qdenom_exp))
             # = (1e-12) / (((1e-2)**v['HCO3'])*((1e-7)**v['H'])*((1e-8)**v['NH4'])*((1e-5)**v['H2S'])*(a_e(273.15+25)**v['e']))
            G = v['stdGr']+(8.31*(TK)*log10Q*math.log(10))
            SynthCost += (G*v['umol']*1e-6)
        return SynthCost



    @staticmethod
    def get_ESynth_density(T, AA=True, compute=None):
        """Calculate the energy cost to synthesise one cell per dry g.

        If you want to compute ESynth density for a new organism, pass compute as a dictionary in the following format: {'dbfilename':None, 'host':base_organism}
        """
        if compute != None:
            # first see if a new db filename has been passed
            try:
                AAsynth, Psynth, = cell_synthesis.extract_AAP_from_db(T=T, filename=compute['dbfilename'])
            except:
                AAsynth = cell_synthesis.AAsynthCost(T)
                Psynth  = cell_synthesis.ProteinSynthCost(T, compute['host'].dry_mass, compute['host'].base_volume, host_proteinfraction=compute['host'].proteinfraction)
                Psynth = Psynth
                if type(compute['dbfilename']) is type('str'):
                    cell_synthesis.add_AAP_to_db(T, AAsynth, Psynth, filename=compute['dbfilename'])
        else:
            try:
                AAsynth, Psynth, = cell_synthesis.extract_AAP_from_db(T=round(T), filename='EColiSynth') #csv?
            except:
                AAsynth = cell_synthesis.AAsynthCost(T)
                Psynth  = cell_synthesis.ProteinSynthCost(T, 2.8e-16, 1e-18)
                cell_synthesis.add_AAP_to_db(T, AAsynth, Psynth)
                Psynth = Psynth
        if AA:
            # return costs for AA synthesis and protein synthesis
            return AAsynth + Psynth
        else:
            return Psynth


    @staticmethod
    def add_AAP_to_db(T, AAsynth, Psynth, filename='EColiSynth'):
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


    # get E_synth density function for use in base_organism
    # database for saved densities at certain temperatures. Do it for 0-100 E Coli?


    # TODO LATER
    # distance this from E Coli











    """ LEGACY CODE FROM OLD PROTEIN SYNTH METHOD
    @staticmethod
    def ProteinSynthCost_legacy(TK, host_drymass, host_volume, host_proteinfraction=0.55, cmfile=None, aafile=None):
        TK =round(TK)

        # get arrays of copynumber and molecular mass for each protein Ishihama has data for

        CopyNo, MolMass = importers.copies_moles()

        AAList = importers.AAreader(host_drymass, host_volume)
        #generateMCParams(AAList) # assign on the parameters we need for Monte Carlo
        AABB = BioMolecule('AABB', 74.05866, T=TK)
        PBB = BioMolecule('PBB', 56.04918, T=TK)
        GLY = BioMolecule('GLY', 57.05892, T=TK)
        Water = BioMolecule('H2O(l)', 18, T=TK) #lol

        avgAAMM = 108.0197 # u



        # produce a gaussian kernel density
        # create arrays which contain every single protein, one for mass one for length
        MolMassCount = []
        LengthCount = []
        Length = []
        # print(CopyNo)
        # USE THIS TO JUST GET THE LENGTH
        for m, c in zip(MolMass, CopyNo):
            if c>1.e5: # limit the really big copy counts (could be a good idea to throttle them instead?)
                c = 1.e5 #throttle down
            l = m/avgAAMM
            counter=1
            while counter<c:
                MolMassCount.append(m)
                LengthCount.append(l)
                counter+=1


        # get a gaussian kernel density
        # gkde_Ishihama = gaussian_kde(LengthCount, bw_method=0.017)

        # we need to keep count of how much dry weight of protein we have produced
        TotalPMassSim = 0. #total mass of protein simulated in u
        TotalESim = 0. #total energy requirement simulated in J

        PMass_cell = host_proteinfraction*host_drymass/(1.660539e-27) #in u


        # get a test sample of the desired size
        #sample = gkde_Ishihama.resample(samplesize)[0]
        # sample = [278]
        # round each length to a whole number
        #mean = sum(sample)/len(sample)

        meanlength = int(statistics.mean(LengthCount))


        #loop through each protein, and add a proportionate AA each time, given the composition of E Coli
        ProteinGibbs = AABB.std_formation_gibbs #the amino acid backbone
        ProteinMass = AABB.Mr # its mass
        ReactionStdGibbs = 0.
        replist = cell_synthesis.getreplistAA(AAList, meanlength)
        # build the protein
        i=0
        for AA in replist:
        # while i<a:

            # AA = getRandomAA(AAList)

            # input it to the overall reaction calculation
            ReactionStdGibbs -= AA.std_formation_gibbs


            #Glycine has special requirements...
            if AA.name == 'Glycine(aq)':
                # add GLY and no PBB
                ProteinGibbs += GLY.std_formation_gibbs
                ProteinMass += GLY.Mr
            else:
                if i!=0:
                    ProteinGibbs += PBB.std_formation_gibbs # add a protein backbone
                    ProteinMass += PBB.Mr # at the backbone mass
                    # add the formed water molecule to the overall reaction calculation
                    ReactionStdGibbs += Water.std_formation_gibbs
                ProteinGibbs += AA.std_formation_R # add the R group Gibbs
                ProteinMass += (AA.Mr - AABB.Mr) # add the R group mass


            i+=1

        # we now have the mass of this protein [u] and its free energy of formation [J/mol].

	    # now complete the Gibbs of the reaction AA -> protein + water, we added the AA and water as we built the protein
        ReactionStdGibbs += ProteinGibbs

        ProteinConc_NoPerCell = meanlength*len(LengthCount)*3.


        # ProteinConc_NoPerCell = gkde_Ishihama.evaluate([a])*len(LengthCount)*3. # no proteins / cell , supposedly data covers approx one third of cell

		# the protein concentration could actually be less than this, as we assume the final conc of the formed cell.

        ProteinConc_MolPerL = ProteinConc_NoPerCell/(6.022e23 * host_volume*1000) # convert to SI mol/l
        # print(ProteinConc_MolPerL)
        # get the ln reaction quotient which is a damn pain
        lnQ = cell_synthesis.lnProteinReactionQuotient(ProteinConc_MolPerL, replist)


        # now get the Gibbs of reaction

        ReactionGibbs = ReactionStdGibbs+(8.31*(TK)*lnQ)


        TotalPMassSim += ProteinMass
        TotalESim += (ReactionGibbs/(6.022e23)) # J/protein

        TotalEReq = TotalESim*(PMass_cell/TotalPMassSim)/host_proteinfraction #tot energy for proteins in cell

        # return as a density (per dry g)
        return TotalEReq/(host_drymass*1000)
    """
