# -*- coding: utf-8 -*-

"""
@package    Refeed
@brief      Helper class for Refeed to represent References
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
from Sequence import Sequence, Sequence_with_masker, Sequence_with_replacer

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Reference(object):
    """
    Represent a reference fasta file containing several sequences. Use with the context manager
    to remove temporary files generated during the the parsing and indexation of the fasta
    sequence. Alternatively, the clean method can be called at the end of the object usage
    """
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

    def __init__ (self, name, fasta, masking=True):
        """
        Create a reference object extract fasta ref if needed and create a sequence object per
        sequences found in the fasta file
        @param name     Name of the Reference
        @param fasta    Path to a fasta file (can be gzipped)
        @param masking  Will create a Sequence object with the ability to output a sequence where
        areas of the sequence overlapping hits will be masked.
        if True. If not the Sequence object will have the ability to replace the original sequence
        with the hit sequences
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

            for seq_name, seq_record in fasta_record.items():
                assert seq_name not in self.seq_dict, "Reference name <{}> is duplicated in <{}>".format(seq_name,self.name)

                # Define a Sequence object depending of the type of output required
                if masking:
                    self.seq_dict[seq_name] = Sequence_with_masker(name=seq_name, seq_record=seq_record)
                else:
                    self.seq_dict[seq_name] = Sequence_with_replacer(name=seq_name, seq_record=seq_record)

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
        """
        Parse a list of BlastHit objects and attibute each of them to its matching Sequence
        @param hit_list A list of BlastHit objects
        """

        for hit in hit_list:
            try:
                self.seq_dict[hit.s_id].add_hit(hit)
            except KeyError as E:
                print ("No sequence matching with the hit subject id")

    def output_reference (self, compress=True):
        """
        Output a reference corresponding to the original sequenced but masked with a masking
        character for bases overlapped by a BlastHit.
        @param compress Compress the output fasta file
        @param masking_char Character that will be used to mask hits
        """
        # Count the number of hit in all Sequence objects from the Reference
        if not self.n_hit:
            print ("No hit found in all sequence of the reference {}".format(self.name))
            fasta_path = None

        # Write a new compressed reference in the current folder
        elif compress:
            fasta_path = "{}_modified.fa.gz".format(self.name)
            with gopen (fasta_path, "wb") as fasta:
                for seq in self.seq_dict.values():
                    # Write the sequence in the fasta file
                    fasta.write(">{}\n{}\n".format(seq.name, seq.output_sequence()))

        # Write a new uncompressed reference in the current folder
        else:
            fasta_path = "{}_modified.fa".format(self.name)
            with open (fasta_path, "w") as fasta:
                for seq in self.seq_dict.values():
                    # Write the sequence in the fasta file
                    fasta.write(">{}\n{}\n".format(seq.name, seq.output_sequence()))

        return fasta_path

    def get_report (self):
        """Generate a report under the form of a list"""
        pass

    def clean (self):
        print ("Cleaning up temporary files for the reference \"{}\"".format(self.name))
        # Remove the temporary directory containing files generated during program execution
        rmtree(self.temp_dir)
        # Cleanup the self dictionary
        self.__dict__ = {}
