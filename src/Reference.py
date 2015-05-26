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
from os import path
from shutil import rmtree
from collections import OrderedDict
from gzip import open as gopen
from tempfile import mkdtemp

# Specific Third party import
import pyfasta # install with pip

# Local imports
from FileUtils import is_readable_file, is_gziped, gunzip, cp
from Sequence import Sequence

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Reference(object):
    """ Represent a reference fasta file containing several Sequences """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    REFERENCE_NAMES = []

    #~~~~~~~CLASS METHODS~~~~~~~#

    @ classmethod
    def ADD_TO_REFERENCE_NAMES(self, name):
        self.REFERENCE_NAMES.append(name)

    @ classmethod
    def RESET_REFERENCE_NAMES (self):
        self.REFERENCE_NAMES = []

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, name, fasta):
        """
        Create a reference object extract fasta ref if needed and create a sequence object per
        sequences found in the fasta file
        """
        print ("Create {} object".format(name))
        # Create self variables
        self.name = name
        self.temp_dir = mkdtemp()

        try:
            # Test values
            assert self.name not in self.REFERENCE_NAMES, "Reference name <{}> is duplicated".format(self.name)
            assert is_readable_file(fasta), "{} is not a valid file".format(fasta)

            # If gziped, ungzip the reference fasta file in the temporary folder. If not compress copy
            # in the temporary folder

            if is_gziped(fasta):
                print ("\tUnzip fasta file in a temporary directory")
                self.fasta = gunzip(fasta, self.temp_dir)
            else:
                print ("\tCopy fasta file in a temporary directory")
                self.fasta = cp(fasta, self.temp_dir)

            # Loading the fasta sequence in a pyfasta.Fasta (seq_record is a mapping and not a str)
            print ("\tParsing the file with pyfasta")
            self.seq_dict = OrderedDict ()
            fasta_record = pyfasta.Fasta(self.fasta, flatten_inplace=True)
            print ("\tFound {} sequences in {}".format(len(fasta_record) , self.name))

            for id, seq_record in fasta_record.items():
                seq_name, sep, descr = id.partition(" ")
                assert seq_name not in self.seq_dict, "Reference name <{}> is duplicated in <{}>".format(seq_name,self.name)
                self.seq_dict [seq_name] = Sequence(name=seq_name, seq_record=seq_record, descr=descr)

            # Add name to a class list
            self.ADD_TO_REFERENCE_NAMES(self.name)

        except Exception as E:
            self.clean()
            raise E

    # Enter and exit are defined to use the context manager "with"
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        """Destructor to remove the file created by Reference"""
        self.clean()

    # Typical string methods
    def __str__(self):
        msg = "REFERENCE CLASS\tParameters list\n"
        msg+= "  Name: {}\n".format(self.name)
        msg+= "  Temporary dir: {}\n".format(self.temp_dir)
        msg+= "  Fasta path: {}\n".format(self.fasta)
        msg+= "  Number of sequences: {}\n".format(self.n_seq)
        msg+= "  Number of hit(s) in sequences: {}\n".format(self.n_hit)
        for s in self.seq_dict.values():
            msg+= "    Name: {}\tSeq: {}...\tNumber of hits: {}\n".format(
                s.name, s.seq_record[0:10], len(s.hit_list))
        return (msg)

    def __repr__(self):
        return ("<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__))

    #~~~~~~~PROPERTIES AND MAGIC~~~~~~~#

    @property
    def n_hit(self ):
        """Count the number of hits in all the Sequences of the Reference"""
        n_hit = 0
        for sequence in self.seq_dict.values():
            n_hit += sequence.n_hit
        return n_hit

    @property
    def n_seq(self ):
        return len(self.seq_dict)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def add_hit_list (self, hit_list):
        """Parse a list of BlastHit objects and attibute each of them to its matching Sequence"""

        for hit in hit_list:
            try:
                self.seq_dict[hit.s_id].add_hit(hit)
            except KeyError as E:
                print ("No sequence matching with the hit subject id")

    def output_masked_reference (self):

        # Count the number of hit in all Sequence objects from the Reference
        if not self.n_hit:
            print ("No hit found in all sequence from the reference {}".format(self.name))
            return None

        # Else = write a new reference in the current folder
        with gopen ("{}_modified.fa.gz".format(self.name), "wb") as fasta:
            for seq in self.seq_dict.values():
                fasta.write("{}_{}\n{}\n".format(seq.name, seq.descr, seq.output_masked_sequence()))

    def clean (self):
        print ("Cleaning up temporary files for the reference \"{}\"".format(self.name))
        rmtree(self.temp_dir)
