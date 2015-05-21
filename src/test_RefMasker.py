# -*- coding: utf-8 -*-

"""
@package    RefMasker
@brief      Test function to be used with python package RefMasker
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com> <adrien.leger@inserm.fr> <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# IMPORTS #########################################################################################

# Standard library packages import
import sys, os, string, tempfile
from random import randint as ri
from random import uniform as rf
from random import choice as rc

# Third party packages import
import pytest
import pyfasta

# Import the current working dir in the python path to allow local package imports
sys.path.append(os.getcwd())

# HELPER FUNCTIONS AND CLASSES ####################################################################

"""Generate a random string"""
def rs (length):
    return ''.join(rc(string.ascii_lowercase + string.digits) for _ in range(length))

"""Generate a random DNA string"""
def rDNA (length):
    return ''.join(rc(["A","T","C","G"]) for _ in range(length))

"""Helper class to generates a random fasta file on the fly"""
class rand_fasta (object):
    def __init__(self, len_seq=500, n_seq=1):
        self.seq_dict = {}
        self.fasta_path = tempfile.mkstemp()[1]
        with open (self.fasta_path, "w") as fp:
            for i in range(n_seq):
                seq = rDNA(len_seq)
                self.seq_dict["ref_{}".format(i)] = seq
                fp.write (">ref_{}\n{}\n".format(i, seq))

    def __enter__(self):
        print ("Create a random fasta file")
        return self

    def __exit__(self, type, value, traceback):
        print ("Destroy fasta files")
        os.remove(self.fasta_path)

# TESTS BLAST HIT #################################################################################

@pytest.yield_fixture
def pyfasta_record(len_seq=500, n_seq=1):
    with rand_fasta(len_seq, n_seq) as r:
        pyfasta_Fasta = pyfasta.Fasta(r.fasta_path, flatten_inplace=True)
        for id, seq_record in pyfasta_Fasta.items():
            yield seq_record

def test_Sequence():
    from Sequence import Sequence

    for seq_record in pyfasta_record(200, 10):
        sequence = Sequence(name = rs(10), seq_record = seq_record, descr = rs(10))
        print (sequence)
        assert sequence.length == 200



#def test_2():
    #from X import X
    #assert(Y)

#def test_3():
    #from X import X
    #assert(Y)
