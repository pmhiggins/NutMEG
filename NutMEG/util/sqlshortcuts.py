

class sqlshortcuts:
    """Shortcuts for common SQL commands we use throughout NutMEG"""

    @staticmethod
    def SELECTpreamble(tabname, SELECT):
        """Return the preamble to a SELECT Query

        tabname : name of tble to select from.
        SELECT : name of column from which to extract data.
        """
        return 'SELECT ' + SELECT + ' FROM ' + tabname + ' WHERE '

    @staticmethod
    def SELECTcolumns(tabname, dbdict, WHERE):
        """Return string for query to get all the values corresponding
        to the keys in dbdict

        dbdict : dictionary of data to select. Only keys (column names) are
        important for this method.
        WHERE : column name to use as identifier (e.g. ID).
        """
        First = True
        lst = 'SELECT '
        for k in sorted(dbdict.keys()):
            if First:
                First = False
                lst += k
                continue
            else:
                lst += ', '
                lst += k
        lst += ' FROM ' + tabname + ' WHERE ' + WHERE + ' = ?'
        return lst

    @staticmethod
    def INSERTlst(tabname, dbdict):
        """Return string for query to insert all keys and values of dbdict.

        tabname : name of table.
        dbdict : dictionary of data to insert. Keys should be column titles
        and values should be the values to insert.
        """
        First=True
        vals=[]
        lst = 'INSERT INTO ' + tabname + '(' + \
         sqlshortcuts.commas_string_keys(dbdict) + ')'
        lst += ' VALUES('
        First=True
        for k in sorted(dbdict.keys()):
            if dbdict[k][0]=='TEXT':
                vals.append(str(dbdict[k][1]))
            else:
                vals.append(dbdict[k][1])

            if First:
                First = False
                lst += '?'
                continue
            else:
                lst += ','
                lst += '?'
        lst += ')'

        return lst, tuple(vals)

    @staticmethod
    def commas_string_keys(dbdict):
        """Return keys of dbdict as one string separated by commas"""
        First=True
        lst= ''
        for k in sorted(dbdict.keys()):
            if First:
                First = False
                lst += k
                continue
            else:
                lst += ', '
                lst += k
        return lst

    @staticmethod
    def commas_string_keys_types(dbdict):
        """Return keys of dbdict and their SQL data type as one string
        separated by commas"""
        First=True
        lst= ''
        for k in sorted(dbdict.keys()):
            if First:
                First = False
                lst += k + ' ' + dbdict[k][0]
                continue
            else:
                lst += ', '
                lst += k + ' ' + dbdict[k][0]
        return lst

    @staticmethod
    def ANDlst(dbdict):
        """Return keys of dbdict in the sql AND x=? format. Also returns the
        corresponding vales as a tuple."""
        First = True
        vals=[]
        lst =''
        for k in sorted(dbdict):
            if dbdict[k][0]=='TEXT':
                vals.append(str(dbdict[k][1]))
            else:
                vals.append(dbdict[k][1])
            if First:
                First = False
                lst += k + ' = ?'
                continue
            else:
                lst += ' AND '
                lst += k + ' = ?'
        return lst, tuple(vals)
