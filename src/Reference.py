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

    def __init__ (self, name, fasta, compress=True):
        """
        Create a reference object extract fasta ref if needed and create a sequence object per
        sequences found in the fasta file
        @param name     Name of the Reference
        @param fasta    Path to a fasta file (can be gzipped)
        @param compress Fasta output will be gzipped if True
        """
        print(("Create {} object".format(name)))
        # Create self variables
        self.name = name
        self.temp_dir = mkdtemp()
        self.compress = compress

        # Create a name for the fasta file to be generated
        self.modified_fasta = "{}_masked.fa{}".format(self.name, ".gz" if self.compress else "")

        try:
            # Test values
            assert self.name not in self.REFERENCE_NAMES, "Reference name <{}> is duplicated".format(self.name)
            assert is_readable_file(fasta), "{} is not a valid file".format(fasta)

            # If gziped, ungzip the reference fasta file in the temporary folder. If not compress
            # copy in the temporary folder

            if is_gziped(fasta):
                print (" * Unzip fasta file in a temporary directory")
                self.fasta = gunzip(fasta, self.temp_dir)
            else:
                print (" * Copy fasta file in a temporary directory")
                self.fasta = cp(fasta, self.temp_dir)

            # Loading the fasta sequence in a pyfasta.Fasta (seq_record is a mapping)
            print (" * Parsing the file with pyfasta")
            seq_dict = {}
            fasta_record = pyfasta.Fasta(self.fasta, flatten_inplace=True)
            print((" * Found {} sequences in {}".format (len (fasta_record), self.name)))

            for name, seq_record in list(fasta_record.items()):

                # Remove additional sequence descriptor in fasta header and create a Sequence object
                short_name = name.partition(" ")[0]
                assert short_name not in seq_dict, "Reference name <{}> is duplicated in <{}>".format(short_name,self.name)
                seq_dict[short_name] = Sequence(name=short_name, seq_record=seq_record)

            # Save to a name sorted ordered dict
            self.seq_dict = OrderedDict(sorted(list(seq_dict.items()), key=lambda x: x))

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
        for s in list(self.seq_dict.values()):
            msg+= "    Name: {}\tSeq: {}...\tNumber of hits: {}\n".format(
                s.name, s.seq_record[0:10], len(s.hit_list))
        return (msg)

    def __repr__(self):
        return ("<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__))

    #~~~~~~~PROPERTIES AND MAGIC~~~~~~~#

    @property
    def n_hit(self ):
        """Count the number of hits in all the Sequences of the Reference"""
        return sum([sequence.n_hit for sequence in list(self.seq_dict.values())])

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

    def output_reference (self):
        """
        Output a reference corresponding to the original sequenced but masked with a masking
        character for bases overlapped by a BlastHit.
        """
        # Count the number of hit in all Sequence objects from the Reference
        if not self.n_hit:
            return None

        # Write a new compressed reference in the current folder
        elif self.compress:
            with gopen (self.modified_fasta, "wb") as fasta:
                for seq in list(self.seq_dict.values()):
                    # Write the sequence in the fasta file
                    fasta.write(">{}\n{}\n".format(seq.name, seq.output_sequence()).encode())
            return self.modified_fasta

        # Write a new uncompressed reference in the current folder
        else:
            with open (self.modified_fasta, "wb") as fasta:
                for seq in list(self.seq_dict.values()):
                    # Write the sequence in the fasta file
                    fasta.write(">{}\n{}\n".format(seq.name, seq.output_sequence()).encode())
            return self.modified_fasta

    def get_report (self, full=False):
        """
        Generate a report under the form of an Ordered dictionary
        @param full If true a dict containing all self parameters will be returned
        """
        report = OrderedDict ()
        report["Reference Name"] = self.name
        report["Number of sequences"] = self.n_seq
        report["Number of hit(s)"] = self.n_hit

        # Include in report only if hit where found in the reference
        if self.n_hit:
            report["Number of base(s) modified"] = sum([seq.mod_bases for seq in list(self.seq_dict.values())])
            report["Modified fasta"] = self.modified_fasta
            report["Modified Sequences"] = OrderedDict ()
            for seq in list(self.seq_dict.values()):
                if seq.hit_list:
                    report["Modified Sequences"][seq.name] = seq.get_report(full=full)

        return report

    def clean (self):
        print((" * Cleaning up temporary files for the reference \"{}\"".format(self.name)))
        # Remove the temporary directory containing files generated during program execution
        rmtree(self.temp_dir)
        # Cleanup the self dictionary
        self.__dict__ = {}
