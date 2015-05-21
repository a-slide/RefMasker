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
from os import remove
from FileUtils import is_readable_file, is_gziped, gunzip_file
from collections import OrderedDict

# Specific Third party import
import pyfasta # install with pip

# Local imports
from Sequence import Sequence

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Reference(object):
    """ Represent a reference containing several Sequences """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    REFERENCE_NAMES = []

    #~~~~~~~CLASS METHODS~~~~~~~#

    @ classmethod
    def ADD_TO_REFERENCE_NAMES(self, name):
        self.REFERENCE_NAMES.append(name)

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, name, fasta):
        """
        Create a reference object extract fasta ref if needed and create a sequence object per
        sequences found in the fasta file
        """
        print ("Create {} object".format(name))
        # Create self variables
        self.name = name

        # Test values
        assert self.name not in self.REFERENCE_NAMES, "Reference name <{}> is duplicated".format(self.name)
        is_readable_file(fasta)

        # Ungzip Reference file if needed for Blast that can't manage compressed files
        if is_gziped(fasta):
            print ("\tUnzip file for blast")
            self.unziped_fasta = True
            self.fasta = gunzip_file (fasta)
        else:
            self.unziped_fasta = False
            self.fasta = fasta

        # Loading the fasta sequence in a pyfasta.Fasta (seq_record is a mapping and not a str)
        print ("\tParsing the file with pyfasta")
        self.seq_dict = OrderedDict ()
        fasta_record = pyfasta.Fasta(self.fasta, flatten_inplace=True)
        print ("\tFound {} sequences in {}".format(len(fasta_record) , self.name))

        for id, seq_record in fasta_record.items():
            name, sep, descr = id.partition(" ")
            assert name not in self.seq_dict, "2 sequences with the same name found in {}".format(self.name)
            self.seq_dict [name] = Sequence(name = name, seq_record = seq_record, descr = descr)

        # Add name to a class list
        self.ADD_TO_REFERENCE_NAMES(self.name)

    def __del__(self):
        """Destructor to remove the database and unziped fasta files if needed"""
        print ("Cleaning and De-initializing the Reference \"{}\"".format(self.name))
        remove (self.fasta+".flat")
        remove (self.fasta+".gdx")
        if self.unziped_fasta:
            remove (self.fasta)

    def __str__(self):
        msg = "REFERENCE CLASS\tParameters list\n"
        msg+= "  Name: {}\n".format(self.name)
        msg+= "  Fasta path: {}\n".format(self.fasta)
        msg+= "  Fasta unziped: {}\n".format(self.unziped_fasta)
        msg+= "  Sequences found\n".format(self.unziped_fasta)
        for s in self.seq_dict.values():
            msg+= "    Name: {}\tSeq: {}...\tNumber of hits: {}\n".format(
                s.name, s.seq_record[0:10], len(s.hit_list))
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def add_hit_list (self, hit_list):
        """Parse a list of BlastHit objects and attibute each of them to its matching Sequence"""

        for hit in hit_list:
            try:
                self.seq_dict[hit.s_id].add_hit(hit)

            except KeyError as E:
                print ("No sequence matching with the hit subject id")
