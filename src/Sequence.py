# -*- coding: utf-8 -*-

"""
@package    RefMasker
@brief      Helper class for RefMasker to represent Sequences
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com> <adrien.leger@inserm.fr> <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Sequence(object):
    """ Represent a single sequence from a fasta file """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    # Class dict to compute the complement of a DNA sequence
    DNA_COMPLEMENT = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, name, seq_record, descr):
        """
        Create a Sequence object that will store
        """
        # Create self variables
        self.name = name
        self.seq_record = seq_record
        self.descr = descr
        self.length = len(self.seq_record)

        # Will be used later to store blast hits
        self.hit_list = []

    def __str__(self):
        msg = "SEQUENCE CLASS\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in self.__dict__.keys()]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __len__ (self):
        """Support for len method"""
        return len(self.seq_record)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def add_hit (self, hit):
        """Parse a list of BlastHit objects and attibute each of them to its matching Sequence"""

        # Hit cannot have borders outside the sequence size
        if hit.s_start > self.length or hit.s_end > self.length:
            raise ValueError, ("Invalid hit: Outside of sequence \"{}\" borders".format(self.name))

        # If the forward strand of the hit is aligned on the forward strand of the subject
        if hit.s_orient and hit.q_orient:
            self.hit_list.append(hit)

        # If only the reference strand is reverse = reverse the subject coordinates and make
        # the reverse complement of the query sequence
        elif not hit.s_orient:
            hit.s_start, hit.s_end = hit.s_end, hit.s_start
            hit.q_seq = self._reverse_complement(hit.q_seq)
            self.hit_list.append(hit)

        # If only the query strand is reverse = reverse the query coordinates and make
        # the reverse complement of the query sequence
        elif not hit.s_orient
            hit.q_start, hit.q_end = hit.q_end, hit.q_start
            hit.q_seq = self._reverse_complement(hit.q_seq)
            self.hit_list.append(hit)

        # If the reverse strand of the hit is aligned on the reverse strand of the subject
        # reverse the coordinates of the sequence but does not modify the query sequence
        else
            it.s_start, hit.s_end = hit.s_end, hit.s_start
            it.q_start, hit.q_end = hit.q_end, hit.q_start
            self.hit_list.append(hit)

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _DNA_reverse_comp (self, sequence):
    """
    Generate the reverese complementary sequence of a given DNA sequence
    @param sequence DNA sequence
    @return Reverse complement of the sequence
    """
    rc = ""
    for base in sequence[::-1]:
        # Return the complement of a base except if an illegal character is found (such as - or .)
        rc += self.DNA_COMPLEMENT.get(base, "")

    return rc
