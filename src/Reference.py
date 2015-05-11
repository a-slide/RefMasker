# -*- coding: utf-8 -*-

"""
@package    RefMasker
@brief      Helper class for RefMasker to represent References
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com> <adrien.leger@inserm.fr> <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""
# Standard library imports
from os import R_OK, access
#from pyDNA.FileUtils import is_readable_file

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Reference(object):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    REFERENCE_NAMES = []

    #~~~~~~~CLASS METHODS~~~~~~~#

    @ classmethod
    def ADD_TO_REFERENCE_NAMES(self, name):
        self.REFERENCE_NAMES.append(name)

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, name, fasta):
        # Create self variables
        self.name = name
        self.fasta = fasta
        self._test_values()

        ###### UNZIP FASTA IF compressed
        ###### LOAD FASTA IN A PYTHON OBJECT (Biopython?)

        # Define additional self variables
        # list of hits found
        # nb of bp masked

        # Add name to a class list
        self.ADD_TO_REFERENCE_NAMES(self.name)

    # Fundamental class methods str and repr
    def __str__(self):
        msg = "REFERENCE CLASS\n\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in self.__dict__.keys()]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        return (msg)

    def __repr__(self):
        #return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
        return self.__str__()

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _test_values(self):
        assert self.name not in self.REFERENCE_NAMES, "Reference name <{}> is duplicated".format(self.name)
        self._is_readable_file(self.fasta)

    def _is_readable_file (self, fp):
        """ Verify the readability of a file or list of file """
        if not access(fp, R_OK):
            raise IOError ("{} is not a valid file".format(fp))
