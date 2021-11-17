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

# Standard library imports
from collections import OrderedDict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Sequence(object):
    """
    Represent a single sequence from a fasta file and can store a list of Blasthits found by Blastn
    Support masking of the original sequence with the *N* where blast hits were found
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, name, seq_record):
        """
        Create a Sequence object that will store
        """
        # Create self variables
        self.name = name
        self.seq_record = seq_record
        self.seq_len = len(self.seq_record)

        # Will be used later to store blast hits
        self.hit_list = []
        self.mod_bases = 0

    def __str__(self):
        msg = "SEQUENCE CLASS\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in list(self.__dict__.keys())]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PROPERTIES AND MAGIC~~~~~~~#

    @property
    def n_hit (self):
        return len(self.hit_list)

    def __len__ (self):
        """Support for len method"""
        return len(self.seq_record)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def add_hit (self, hit):
        """
        Add a hit to hit_list after verification and eventual modifications depending of the
        subject and query orientation
        """
        # Hit cannot have borders outside the sequence size
        if hit.s_start > self.seq_len or hit.s_end > self.seq_len:
            raise ValueError(("Invalid hit: Outside of sequence borders"))
        if hit.s_id != self.name:
            print(hit.s_id)
            print(self.name)
            raise ValueError(("Invalid hit: hit subject name does not match Sequence name"))

        # Reverse coordinates of the subject if the orientation of read is negative
        if hit.s_orient == "-":
            hit.s_start, hit.s_end = hit.s_end, hit.s_start

        # Reverse coordinates of the query if the orientation of read is negative
        if hit.q_orient == "-":
            hit.q_start, hit.q_end = hit.q_end, hit.q_start

        # Finally append the modified hit to the list
        self.hit_list.append(hit)


    # TODO : Create an option for soft masking
    def output_sequence (self):
        """
        Output a sequence corresponding to the original seq record sequence but masked with
        a masking character for bases overlapped by a BlastHit
        """
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
                # Write previous masked interval and update the counter of masked bases
                masked_seq += "N"*(end_mask-start_mask)
                self.mod_bases += end_mask-start_mask
                # Write current sequence interval
                masked_seq += str(self.seq_record[end_mask:hit.s_start])
                # Update start_mask and end_mask borders
                start_mask, end_mask = hit.s_start, hit.s_end

        # Tail of the list. Write previous masked interval and update the counter of masked bases
        masked_seq += "N"*(end_mask-start_mask)
        self.mod_bases += end_mask-start_mask
        # Write the end of sequence if there is something left to write
        masked_seq += str(self.seq_record[end_mask:self.seq_len])

        return masked_seq

    def get_report (self, full=False):
        """
        Generate a report under the form of an Ordered dictionary
        @param full If true a dict containing all self parameters will be returned
        """
        report = OrderedDict ()
        report["Number of hits"] = self.n_hit

        # Include in report only if hit where found in the sequence
        if self.n_hit:
            report["Number of base modified"] = self.mod_bases
            if full:
                report["Sequence length"] = self.seq_len
                report["Percent of modified bases"] = float(self.mod_bases)/self.seq_len*100.0

                # Details of blast hits
                report["Blast Hits"] = OrderedDict ()

                # Sort the hit list according to the name and start coordinate of the query
                hit_list = sorted( self.hit_list, key=lambda x: (x.q_id, x.q_start))
                for i, hit in enumerate(hit_list):
                    report["Blast Hits"]["Hit {:03d}".format(i+1)] = hit.get_report(full=full)

        return report
