"""
	THIS IS PYTHON 2.7 CODE FOR MAKING USE OF REAKTORO
	
	reads in gift.csv to calculate the equilibrium constant of a reaction,
	may be extended in the future for more calculations
"""



import sys
sys.path.append("/home/s1771909/rto/lib/python2.7/site-packages")
from reaktoro import *


def reader():
	params = open('gift.csv', "r")
	
	reaction = ""
	eT = []
	eP = []
	
	params_input = params.readlines()
	firstline=True # which line number of data we're on
	for line in params_input:
		if not line.startswith("#") and len(line) > 0:   # Kill comments and blanks, though there shouldn't be any anyway
			token = line.split(",") #split via commas as its a csv file
			if firstline: 
				# store the reaction name, this will only be printed on the first line
				reaction = token[0]
				firstline=False
			else:
				# now we have all of the relevant data
				eT.append(float(token[0]))
				eP.append(float(token[1]))
						
	params.close()
	return reaction, eT, eP



reaction, eT, eP = reader()

#use the SUPCRT07 database as the BL one is bonkers
db = Database("supcrt07.xml")
thermo = Thermo(db)

lnK = [thermo.lnEquilibriumConstant(T, P, reaction).val for T, P in zip(eT, eP)]
stdG = [-8.314472*T*lnEQ for T, lnEQ in zip(eT, lnK)]

export = open('reak.csv', 'w')
export.write(reaction + ",\n")
for i in range(0, len(lnK)):
	export.write(str(eT[i]) + "," + str(eP[i]) + "," + str(lnK[i]) + "," + str(stdG[i]) + ",\n")

print "Fe++ " + (str((thermo.standardPartialMolarGibbsEnergy(T=273, P=100000, species="Fe++")).val)+",\n")
print "O2(g) " + (str(thermo.standardPartialMolarGibbsEnergy(T=273, P=100000, species="O2(g)").val)+",\n")
print "H+" + (str(thermo.standardPartialMolarGibbsEnergy(T=273, P=100000, species="H+").val)+",\n")
print "Fe+++" + (str(thermo.standardPartialMolarGibbsEnergy(T=273, P=100000, species="Fe+++").val)+",\n")
print "H2O(l) " + (str(thermo.standardPartialMolarGibbsEnergy(T=273, P=100000, species="H2O(l)").val)+",\n")
print "G " + (str(stdG[0]))

export.close()


	
