# -*- coding: utf-8 -*-

"""
@package    ReMasker
@brief      Test function to be used with python package ReMasker
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com> <adrien.leger@inserm.fr> <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# IMPORTS #########################################################################################

# Standard library packages import
import sys, string, filecmp
from os import getcwd, path, remove
from random import randint as ri
from random import uniform as rf
from random import choice as rc
from shutil import rmtree
from tempfile import mkdtemp
from gzip import open as gopen
from collections import OrderedDict

# Third party packages import
import pytest
import pyfasta

# Import the current working dir in the python path to allow local package imports
sys.path.append(getcwd())

# local package imports
from Sequence import Sequence
from Reference import Reference
from pyBlast.BlastHit import BlastHit
from pyBlast.Blastn import Blastn

# HELPER FUNCTIONS AND CLASSES ####################################################################

def rs (length):
    """Generate a random string"""
    return ''.join(rc(string.ascii_lowercase + string.digits) for _ in range(length))

def rDNA (length):
    """Generate a random DNA string"""
    return ''.join(rc(["A","T","C","G"]) for _ in range(length))

class rand_fasta (object):
    """Helper class to generates a random fasta file on the fly"""
    def __init__(self, len_seq=500, n_seq=1, gziped=False):
        self.seq_dict = {}
        self.temp_dir = mkdtemp()

        if gziped:
            self.fasta_path = path.join(self.temp_dir+"/random.fa.gz")
            with gopen (self.fasta_path, "w") as fp:
                for i in range(n_seq):
                    seq = rDNA(len_seq)
                    self.seq_dict["seq_{}".format(i)] = seq
                    fp.write (">seq_{}\n{}\n".format(i, seq))
        else:
            self.fasta_path = path.join(self.temp_dir+"/random.fa")
            with open (self.fasta_path, "w") as fp:
                for i in range(n_seq):
                    seq = rDNA(len_seq)
                    self.seq_dict["seq_{}".format(i)] = seq
                    fp.write (">seq_{}\n{}\n".format(i, seq))

    def __enter__(self):
        print ("\nCreate a random fasta file")
        return self

    def __exit__(self, type, value, traceback):
        print ("Destroy temporary files")
        rmtree(self.temp_dir)

class defined_fasta (object):
    """Helper class to generates fasta file from a defined dictionary of sequence"""
    def __init__(self, seq_dict={"seq_0":"ATCGATCG", "seq_1":"CGTATCGA",}):
        self.seq_dict = seq_dict
        self.temp_dir = mkdtemp()
        self.fasta_path = path.join(self.temp_dir+"/defined.fa")
        with open (self.fasta_path, "w") as fp:
            for name, seq in list(self.seq_dict.items()):
                fp.write (">{}\n{}\n".format(name, seq))

    def __enter__(self):
        print ("\nCreate a random fasta file")
        return self

    def __exit__(self, type, value, traceback):
        print ("Destroy temporary files")
        rmtree(self.temp_dir)

# FIXTURES FOR TESTS ##############################################################################

@pytest.yield_fixture
def yield_sequence(len_seq=500, n_seq=1):
    with rand_fasta(len_seq, n_seq) as r:
        print ("Create a pyfasta record")
        pyfasta_Fasta = pyfasta.Fasta(r.fasta_path, flatten_inplace=True)
        for name, seq_record in list(pyfasta_Fasta.items()):
            yield Sequence(name = name, seq_record = seq_record)

@pytest.yield_fixture
def yield_BlastHit(len_seq=100, n_hit=5, min_len_hit=5, max_len_hit=20, s_id="seq_0"):
    for _ in range (n_hit):
        length = ri(min_len_hit, max_len_hit)
        start = ri(1, len_seq-length)
        end = start+length
        seq = length*"N"
        yield(BlastHit(s_id=s_id, s_start=start, s_end=end, q_seq=seq))

@pytest.yield_fixture
def yield_sequence_and_BlastHit(len_seq=500, n_query=10, min_len_hit=20, max_len_hit=40):

    # Generate a random subject sequence
    with rand_fasta(len_seq=len_seq, n_seq=1) as subject:
        pyfasta_gen = pyfasta.Fasta(subject.fasta_path, flatten_inplace=True)
        seq_record = pyfasta_gen["seq_0"]
        sequence = Sequence(name = "seq_0", seq_record = seq_record)

        seq_dict = {}
        # Generate random hits from the subject
        for i in range(n_query):
            start = ri(1, len_seq-max_len_hit)
            end = start+ri(min_len_hit, max_len_hit)
            seq_dict["query_{}:{}-{}".format(i, start, end)] = str(seq_record[start:end])
        with defined_fasta(seq_dict) as query:

            # Create Blast DB and perform a blast to generate a list of hits
            with Blastn(subject.fasta_path) as blastn:
                hit_list = blastn(query.fasta_path, task="blastn", best_query_hit=True,)

            yield (sequence, hit_list)

@pytest.yield_fixture
def yield_reference(n_ref=1, len_seq=500, n_seq=1, gziped=False):
    for i in range(n_ref):
        with rand_fasta(len_seq, n_seq, gziped) as fasta:
            with Reference(name="ref_{}".format(i), fasta=fasta.fasta_path) as ref:
                yield ref

# TESTS SEQUENCE CLASS ############################################################################

def test_Sequence_create():
    """Test to verify sequence creation"""
    for sequence in yield_sequence(len_seq=200, n_seq=1):
        print ("Test Sequence object creation")
        assert len(sequence) == 200

@pytest.mark.parametrize("s_len, s_id, s_start, s_end, q_start, q_end, q_seq, s_start_res, s_end_res, q_start_res, q_end_res, q_seq_res", [
    pytest.mark.xfail((100, "seq_0", 90, 110, 0, 0, "" ,0, 0, 0, 0, "")),
    pytest.mark.xfail((100, "seq_1", 80, 100, 0, 0, "" ,0, 0, 0, 0, "")),
    (100, "seq_0", 80, 90, 20, 30, "ATCG" ,79, 90, 19, 30, "ATCG"),
    (100, "seq_0", 90, 80, 20, 30, "ATCG" ,79, 90, 19, 30, "CGAT"),
    (100, "seq_0", 80, 90, 30, 20, "ATCG" ,79, 90, 19, 30, "CGAT"),
    (100, "seq_0", 90, 80, 30, 20, "ATCG" ,79, 90, 19, 30, "ATCG")])

def test_Sequence_add_hit(s_len, s_id, s_start, s_end, q_start, q_end, q_seq, s_start_res, s_end_res, q_start_res, q_end_res, q_seq_res):
    """Test the addition of valid and invalid hit to a Sequence object"""
    for sequence in yield_sequence(len_seq=s_len, n_seq=1):

        hit = BlastHit(s_id=s_id, s_start=s_start, s_end=s_end, q_start=q_start, q_end=q_end, q_seq=q_seq)
        sequence.add_hit(hit)
        assert sequence.n_hit == 1
        assert sequence.hit_list[0].s_start == s_start_res
        assert sequence.hit_list[0].s_end == s_end_res
        assert sequence.hit_list[0].q_start == q_start_res
        assert sequence.hit_list[0].q_end == q_end_res
        #assert sequence.hit_list[0].q_seq == q_seq_res

@pytest.mark.parametrize("len_seq, n_hit", [(100, 1), (100, 5), (200, 10)])

def test_Sequence_output_sequence_1 (len_seq, n_hit):
    """
    Test the capacity of a Sequence object to generate a masked sequence from fake a list of fake
    BlastHit objects
    """
    for sequence in yield_sequence(len_seq=len_seq, n_seq=1):

        # If no hit the return seq should be the same than the original one
        assert sequence.output_sequence() == str(sequence.seq_record)

        hit_list =[]
        for hit in yield_BlastHit(len_seq=len_seq, n_hit=n_hit, min_len_hit=15, max_len_hit=25):
            hit_list.append(hit)
            sequence.add_hit(hit)

        # Verify that the correct number of hit were added to the Sequence
        assert sequence.n_hit == n_hit

        # With the list of hit the function should now return a modified sequence
        assert sequence.output_sequence() != str(sequence.seq_record)

        # For visual confirmation of proper masking
        print((sequence.seq_record))
        print((sequence.output_sequence()))
        for hit in hit_list:
            print((" "*hit.s_start+hit.q_seq))

def test_Sequence_output_sequence_2 ():
    """
    Test the capacity of a Sequence object to generate a masked sequence from a list of BlastHit
    objects generated by query sequences sampled from a random subject and aligned with Blastn.
    """
    for sequence, hit_list in yield_sequence_and_BlastHit(len_seq=100, n_query=5, min_len_hit=15, max_len_hit=25):
        for hit in hit_list:
            sequence.add_hit(hit)

        for i, j in list(sequence.get_report().items()):
            print(("{}\t{}".format(i,j)))

        masked_sequence = (sequence.output_sequence())

        # For visual confirmation of proper masking
        print((sequence.seq_record))
        print (masked_sequence)
        for hit in hit_list:
            print((" "*hit.s_start+hit.q_seq))

        # Verify than all hit coordinates on subject query are properly masked
        for hit in hit_list:
            for position, base in enumerate (masked_sequence):
                if hit.s_start <= position < hit.s_end:
                    assert base == 'N', "The base in position {} of {} should be masked".format(position,sequence.name)

# TESTS REFERENCE CLASS ############################################################################

@pytest.mark.parametrize("n_ref, len_seq, n_seq, gziped", [
    (1, 1000, 1, False),
    (1, 1000, 1, True),
    (2, 10000, 2, False),
    (2, 10000, 2, True)])

def test_Reference_create(n_ref, len_seq, n_seq, gziped):
    """Test to verify Reference creation"""
    # Generates different Reference object from different combinations
    for ref in yield_reference(n_ref, len_seq, n_seq, gziped):
        assert ref.n_seq == n_seq, "Not enough sequences generated"
    # Reset ref name to start a new round
    Reference.RESET_REFERENCE_NAMES()

@pytest.mark.parametrize("n_ref, len_seq, n_seq" , [(1, 1000, 1), (2, 10000, 2)])

def test_Reference_add_hit_list(n_ref, len_seq, n_seq):
    """Test the ability of add_hit_list to add a proper number of fake hits to each sequence in ref"""
    for ref in yield_reference(n_ref=n_ref, len_seq=len_seq, n_seq=n_seq):

        seq_dict = OrderedDict()
        hit_list = []
        for i in range (n_seq):
            n_hit = ri (0,10)
            s_id = "seq_{}".format(i)
            seq_dict[s_id] = n_hit
            hit_list.extend([hit for hit in yield_BlastHit(len_seq=len_seq, n_hit=n_hit, s_id=s_id)])

        ref.add_hit_list(hit_list)

        for name, n_hit in list(seq_dict.items()):
            assert ref.seq_dict[name].n_hit == n_hit

        #ref.output_masked_reference(compress=False)
    # Reset ref name to start a new round
    Reference.RESET_REFERENCE_NAMES()


def test_Reference_output_masked_reference():
    """Test the all Reference methods with predetermined ref and hits"""

    original_ref = OrderedDict()
    original_ref ["s0"] = "ATATTCGCTGCTCAGCTCGTCGATCGCGCACAGCTACGTCGCGCTGTCATGACGTTAGCATCGTACGCATATACTTCAGCTCTGACGATCAGCTACGACTCGCAGCTAGCCTACGACTAGCAGTCTCGTACCGTATGCGTCAGCTCAGCATCAT"
    original_ref ["s1"] = "CGCTAGCTAGCGATTATTCGCTAGCTCGTCGATCGTCGCCTAGATCGTCGATCAGTCAGTCGTCAGTGTGTAACACGTCAGTCAAGCTCAGCTAGCTCGTCGTATGACGCGCGCTATATCGTAGACGATACGCTCGACTAGCATGCATCAGCAT"

    expected_ref = OrderedDict()
    expected_ref ["s0"] = "NNNNNNNNNNNNNNNNNNNNNNNNNNCGCACAGCTACGTCGCGCTGTCATGACGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAGCCTACGACTAGCAGTCNNNNNNNNNNNNNNNNNNNNNNNNNNNAT"
    expected_ref ["s1"] = "CGCTAGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTAGATCGTCNNNNNNNNNNNNNNNNNNNNNNNNNACGTCAGTCAAGCTCAGCTAGCTCGTCGTATGACGCGCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNCATGCATCAGCAT"

    query = OrderedDict()
    query ["q0"] = "CGATCGACGAGCTGAGCAGCGAATAT"       # reverse match in s0
    query ["q1"] = "TAGCATCGTACGCATATACTTCAGCTCTGA"   # forward match in s0
    query ["q2"] = "CTCTGACGATCAGCTACCACTCGCAGC"      # forward match in s0 with 1 mismatch
    query ["q3"] = "TCGTACCGTATGCGTACAGCTCAGCATC"     # forward match in s0 with 2 mismatches and 1 insertion
    query ["q4"] = "CGATCGACGAGCTAGCGAATAATCGCTA"     # reverse match in s1
    query ["q5"] = "GATCAGTCAGTCGTCAGTGTGTAAC"        # forward match in s1
    query ["q6"] = "GATTATTCGCTAGTCGTCGATCGTCGC"      # forward match in s1 with 1 deletion
    query ["q7"] = "CTATATACGTAGACGATACGTCGACTAGG"    # forward match in s1 with 1 mismatches, 1 deletion and 1 insertion

    # Create a reference fasta file
    with defined_fasta(seq_dict=original_ref) as ref_fasta:
        # Create a query fasta file
        with defined_fasta(seq_dict=query) as query_fasta:
            # Blast reference vs query and return a list of hit found
            with Blastn(ref_fasta.fasta_path) as blastn:
                hit_list = blastn(query_fasta.fasta_path, task="blastn", best_query_hit=True,)

        # Create a Reference Object for the reference fasta file
        with Reference(name="ref0", fasta=ref_fasta.fasta_path, compress=False) as reference:

             # Should return None since no hit were added
            print ("Test output_masked_reference without hit added")
            assert not reference.output_reference ()

            # Add test
            reference.add_hit_list (hit_list)
            print ("Test output_masked_reference with hits added")
            masked_obtained = reference.output_reference ()

            with defined_fasta(seq_dict=expected_ref) as masked_expected:
                # files should be equal
                with open(masked_obtained, "r") as obtained:
                    with open(masked_expected.fasta_path, "r") as expected:
                        for l1, l2 in zip(obtained.readlines(), expected.readlines()):
                            assert l1 == l2

            remove (masked_obtained)

    Reference.RESET_REFERENCE_NAMES()
