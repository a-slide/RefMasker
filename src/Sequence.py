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

    def __init__ (self, name, seq_record, descr=""):
        """
        Create a Sequence object that will store
        """
        # Create self variables
        self.name = name
        self.seq_record = seq_record
        self.descr = descr
        self.seq_len = len(self.seq_record)

        # Will be used later to store blast hits
        self.hit_list = []
        self.mask_base = 0

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

    #~~~~~~~PROPERTIES AND MAGIC~~~~~~~#

    @property
    def n_hit(self ):
        return len(self.hit_list)

    def __len__ (self):
        """Support for len method"""
        return len(self.seq_record)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def add_hit (self, hit):
        """Add a hit to hit_list after verification and eventual modifications"""

        # Hit cannot have borders outside the sequence size
        if hit.s_start > self.seq_len or hit.s_end > self.seq_len:
            raise ValueError, ("Invalid hit: Outside of sequence borders")
        if hit.s_id != self.name:
            print hit.s_id
            print self.name
            raise ValueError, ("Invalid hit: hit subject name does not match Sequence name")

        if hit.s_orient == "+":
            # If the forward strand of the hit is aligned on the forward strand of the subject
            if hit.q_orient == "+":
                pass

            # If only the query strand is reverse = reverse the query coordinates and make
            # the reverse complement of the query sequence
            else:
                hit.q_start, hit.q_end = hit.q_end, hit.q_start
                hit.q_seq = self._DNA_reverse_comp(hit.q_seq)

        else:
            # If the reverse strand of the hit is aligned on the reverse strand of the subject
            # reverse the coordinates of the sequence but does not modify the query sequence
            if hit.q_orient == "-":
                hit.s_start, hit.s_end = hit.s_end, hit.s_start
                hit.q_start, hit.q_end = hit.q_end, hit.q_start

            # If only the reference strand is reverse = reverse the subject coordinates and make
            # the reverse complement of the query sequence
            else:
                hit.s_start, hit.s_end = hit.s_end, hit.s_start
                hit.q_seq = self._DNA_reverse_comp(hit.q_seq)

        # Finally append the modified hit to the list
        self.hit_list.append(hit)

    def output_masked_sequence (self):
        """
        Output a sequence corresponding to the original seq record sequence but masked with
        a masking character for bases overlapped by a BlastHit
        """

        ######################################### TODO = Count the number of masked bases

        if not self.hit_list:
            # No need to modify the sequence
            return str(self.seq_record)

        # init a STR buffer to create the masked sequence
        masked_seq = ""
        # Sort list by hit subject start position
        self.hit_list.sort(key=lambda x: x.s_start)

        # Head of the hit_list. If the first hit do not start at the beginning of the subject
        # keep the original sequence until the first hit and save start and end of the first hit
        first_hit = self.hit_list[0]
        start_mask, end_mask = first_hit.s_start, first_hit.s_end
        if first_hit.s_start > 0:
            masked_seq += str(self.seq_record[0:first_hit.s_start])

        for hit in self.hit_list[1:self.n_hit]:
            if hit.s_start <= end_mask+1:
                if hit.s_end > end_mask:
                    end_mask = hit.s_end

            else:
                # Write previous masked interval
                masked_seq += "N"*(end_mask-start_mask)
                # Write current sequence interval
                masked_seq += str(self.seq_record[end_mask:hit.s_start])
                # Update start_mask and end_mask borders
                start_mask, end_mask = hit.s_start, hit.s_end

        # Tail of the hit_list. Write previous masked interval and write the end of sequence if
        # there is something left to write
        masked_seq += "N"*(end_mask-start_mask)
        masked_seq += str(self.seq_record[end_mask:self.seq_len])

        return masked_seq

    #def output_replaced_sequence (self):
        #"""Add a hit to hit_list after verification and eventual modifications"""
        #if not self.hit_list:
            #return None # No need to modify the sequence

        ## init a STR buffer to create the masked sequence
        #Masked_seq = ""
        ## Sort list by hit subject start position
        #self.hit_list.sort(key=lambda x: x.s_start)

        ## CREATE META HITS... I DON'T KNOW HOW ... TAKE INTO ACCOUNT ILLEGAL CHARACTERS...


    def get_report (self):
        """Generate a report under the form of a list"""
        pass

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _DNA_reverse_comp (self, sequence):
        """
        Generate the reverese complementary sequence of a given DNA sequence
        @param sequence DNA sequence
        @return Reverse complement of the sequence
        """
        rc = ""
        for base in sequence[::-1]:
            # Return the complement of a base except if an illegal character is found such as "-"
            # in which case the char is returned)
            try:
                rc += self.DNA_COMPLEMENT[base]
            except KeyError:
                rc += base

        return rc
