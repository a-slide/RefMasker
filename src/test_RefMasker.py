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

# Standard library packages import
import sys, os, random, string

# Third party packages import
import pytest

# Import the current working dir in the python path to allow local package imports
sys.path.append(os.getcwd())

# Helper functions
def rand_string (length):
    return ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(length))

def rand_DNA (length):
    return ''.join(random.choice(["A","T","C","G"]) for _ in range(length))

#def rand_fasta_file (n_seq, len_seq):
#   create a tempory fasta file
#   for i in range (n_seq):
#          fp.write (">Ref_{}\n{}\n".format(i, rand_DNA(len_seq)))
#   return file_path

# Start tests
def test_Sequence():
    from Sequence import Sequence
    from pyfasta import Fasta
    fasta_file = rand_fasta_file(1, 200)
    seq_record = ""
    seq = Sequence (

#def test_2():
    #from X import X
    #assert(Y)

#def test_3():
    #from X import X
    #assert(Y)
