import sqlite3

TPdb = sqlite3.connect('../data/TPdb')

cursor = TPdb.cursor()

#This is for creating the database
cursor.execute('CREATE TABLE reagents(name TEXT, T TEXT, ' + \
				'P TEXT, G TEXT, Hz TEXT, ' + \
				'H TEXT, S TEXT, Cp TEXT, Cv TEXT)') # all of these are standard values
TPdb.commit()
TPdb.close()
