#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
@package    RefMasker
@brief      Main file of the program
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2015
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

try:
    # Standard library imports
    import ConfigParser
    import optparse
    import sys
    from os import R_OK, access
    from time import time
    #from datetime import datetime

    # Local imports
    #from pyDNA.FileUtils import is_readable_file, rm_blank
    from Conf_file import write_example_conf
    from Reference import Reference
    from pyBlast.MakeBlastDB import MakeBlastDB
    from pyBlast.MakeBlastn import MakeBlastn

except ImportError as E:
    print (E)
    print ("Please verify your dependencies. See Readme for more informations\n")
    exit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class RefMasker(object):
    """
    Fastq file demultiplexer, handling double indexing, molecular indexing and filtering based
    on index quality
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "RefMasker 0.1"
    USAGE = "Usage: %prog -c Conf.txt [-i -h]"

    #~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init (self):
        """
        init class method for instantiation from command line. Parse arguments parse CL arguments
        """

        # Define parser usage, options
        optparser = optparse.OptionParser(usage = self.USAGE, version = self.VERSION)
        optparser.add_option('-c', dest="conf_file",
            help= "Path to the configuration file [Mandatory]")
        optparser.add_option('-i', dest="init_conf", action='store_true',
            help= "Generate an example configuration file and exit [Facultative]")

        # Parse arguments
        options, args = optparser.parse_args()

        return RefMasker(options.conf_file, options.init_conf)

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, conf_file=None, init_conf=None):
        """
        Initialization function, parse options from configuration file and verify their values.
        All self.variables are initialized explicitly in init.
        """

        # Create a example conf file if needed
        if init_conf:
            print("Create an example configuration file in the current folder")
            write_example_conf()
            sys.exit(0)

        print("Initialize Quade")
        # Parse the configuration file and verify the values of variables
        try:

            # Verify if conf file was given and is valid
            assert conf_file, "A path to the configuration file is mandatory"
            self._is_readable_file(conf_file)
            self.conf = conf_file

            # Define a configuration file parser object and load the configuration file
            cp = ConfigParser.RawConfigParser(allow_no_value=True)
            cp.read(self.conf)

            print("\tParse Blast options")
            # Blast parameters section
            self.blastn_opt = cp.get("Blast", "blastn_opt")
            self.blast_task = cp.get("Blast", "blast_task")
            self.best_per_query_seq = cp.get("Blast", "best_per_query_seq")
            self.evalue = cp.getfloat("Blast", "evalue")

            assert self.evalue > 0, "Authorized values for evalue: float > 0"

            self.mkblastdb_opt= cp.get("Blast", "mkblastdb_opt")
            self.blastn = cp.get("Blast", "blastn")
            self.mkblastdb = cp.get("Blast", "mkblastdb")

            print("\tParse output options")
            # Output parameters section

            self.repl_char = cp.get("Output", "repl_char")
            self.repl_with_query = cp.getboolean("Output", "repl_with_query")
            self.modif_seq_only = cp.getboolean("Output", "modif_seq_only")
            self.merge_ref = cp.getboolean("Output", "merge_ref")
            self.compress_output = cp.getboolean("Output", "compress_output")

            print("\tParse Reference sequences")
            # Iterate only on sections starting by "reference", create Reference objects
            # And store them in a list
            self.reference_list = []
            for reference in [i for i in cp.sections() if i.startswith("reference")]:
                # Create Reference objects
                self.reference_list.append (
                    Reference (
                        name = self._rm_blank(cp.get(reference, "name"), replace ='_'),
                        fasta = self._rm_blank(cp.get(reference, "fasta"), replace ='\ ')))

        # Handle the many possible errors occurring during conf file parsing or variable test
        except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as E:
            print ("Option or section missing. Report to the template configuration file\n" + E.message)
            sys.exit(1)

        except (ValueError, AssertionError) as E:
            print ("One of the value in the configuration file is not correct\n" + E.message)
            sys.exit(1)

        except (IOError) as E:
            print ("One of the file is incorrect or unreadable\n" + E.message)
            sys.exit(1)

    def __str__(self):
        msg = "RefMasker CLASS\n\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in self.__dict__.keys()]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """
        Mask references homologies iteratively, starting by the last reference which is masked by
        all the others then to the penultimate masked by all others except the last and and so
        forth until there is only 1 reference remaining
        """
        start_time = time()

        print ("Start to process files")
        # Iterate over index in Reference.instances staring by the last one until the 2nd one
        for i in range(len(self.reference_list)-1, 0, -1):
            subject = self.reference_list[i]
            query_list = self.reference_list[0:i]

            print ("Processing Reference {}".format(subject.name)
            hit_list = []

            # Create a MakeBlastDB object of nucl type and fasta input (default)
            dbmaker = MakeBlastDB (
                makeblastdb = self.mkblastdb,
                makeblastdb_opt = self.makeblastdb_opt

            # Create the blast database
            subject_db = dbmaker(ref_path = subject.fasta, db_path = subject.name)

            # Init a MakeBlast object of blastn type (default)
            blaster = MakeBlastn (
                blastn_opt = self.mkblastdb_opt,
                blastn = self.mkblast
                task=self.blast_task
                evalue = self.evalue
                best_per_query_seq = self.best_per_query_seq)

            # Blast each query in query list against the subject database
            for query in query_list:
                print ("Blast {} against {} database".format(query.name, subject.name)
                hit_list = blaster(query_path = query.fasta, db_path = self.makeblastdb_opt)
                ########## subject.add_hits(hit_list)

            # Remove DB files at the and of the blast against all queries
            dbmaker.remove_db_files()

        # Write a report
        print ("Generate_a csv report")
        print ("Done in {}s".format(round(time()-start_time, 3)))
        return(0)

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _is_readable_file (self, fp):
        """ Verify the readability of a file or list of file """
        if not access(fp, R_OK):
            raise IOError ("{} is not a valid file".format(fp))

    def _rm_blank (self, name, replace=""):
        """ Replace blank spaces in a name by a given character (default = remove)"""
        return replace.join(name.split())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    refmasker = RefMasker.class_init()
    refmasker()


## Standard library packages import
#from os import remove, path
#import gzip
#from time import time
#from sys import stdout

## Third party package import
#from Bio import SeqIO

## Local library packages import
#from pyDNA.Utilities import import_seq, file_basename, mkdir
#from Blast import Blastn

##~~~~~~~MAIN METHODS~~~~~~~#

#def _iterative_masker (self): #### TODO The fuction directly manipulate reference field= change that
    #"""
    #Mask references homologies iteratively, starting by the last reference which is masked by
    #all the others then to the penultimate masked by all others except the last and and so
    #forth until there is only 1 reference remaining
    #"""
    ## Iterate over index in Reference.instances staring by the last one until the 2nd one
    #for i in range(Reference.countInstances()-1, 0, -1):

        ## Extract subject and query_list from ref_list
        #subject = Reference.Instances[i]
        #query_list = Reference.Instances[0:i]
        #print ("\n# PROCESSING REFERENCE {} #\n".format(subject.name))

        ## Perform a blast of query list against subject
        #hit_list = Blastn.align (
            #query_list = [ref.ref_fasta for ref in query_list],
            #subject_fasta = subject.ref_fasta,
            #align_opt = self.blastn_opt,
            #db_opt = self.mkblastdb_opt,
            #db_outdir = self.db_dir,
            #db_outname = subject.name)

        ## Masking hits in suject fasta if hits in hit_list
        #subject.ref_fasta = mask (
            #subject_fasta= subject.ref_fasta,
            #hit_list = hit_list,
            #ref_outdir = self.ref_dir,
            #ref_outname = "masked_{}.fa".format(subject.name),
            #compress_ouput = False)


#def mask (  subject_fasta,
            #hit_list,
            #ref_outdir="./references/",
            #ref_outname="masked_ref.fa",
            #compress_ouput=True ):
    #"""
    #Import a reference fasta sequence, Mask positions indicated by hits from a hit_list and write
    #the modified fasta sequence in a new file.
    #@param subject_fasta Fasta sequence of the subject to edit (can be gzipped)
    #@param hit_list List of hit objects. Hits need at least 3 fields named s_id, s_start and s_end
    #coresponding to the name of the sequence matched, and the hit start/end (0 based).
    #@param ref_outdir Directory where the masked reference will be created
    #@param ref_outname Name of the masked reference
    #@param compress_ouput If true the output will be gzipped
    #@return A path to the modified sequence if the hit list was valid.
    #"""

    ## Test if object the first object of hit_list have the require s_id, s_start and s_end fields
    #try:
        #a = hit_list[0].s_id
        #a = hit_list[0].s_start
        #a = hit_list[0].s_end

    #except IndexError:
        #print ("No hit found, The subject fasta file will not be edited")
        #return subject_fasta
    #except AttributeError as E:
        #print ("The list provided does not contain suitable hit object, The subject fasta file will not be edited")
        #return subject_fasta

    ## Initialize output folder
    #mkdir(ref_outdir)

    ## Initialize input fasta file
    #if subject_fasta[-2:].lower() == "gz":
        #in_handle = gzip.open(subject_fasta, "r")
    #else:
        #in_handle = open(subject_fasta, "r")

    ## Initialize output fasta file
    #if compress_ouput:
        #ref_path = path.join (ref_outdir, ref_outname+".gz")
        #out_handle = gzip.open(ref_path, 'w')
    #else:
        #ref_path = path.join (ref_outdir, ref_outname)
        #out_handle = open(ref_path, 'w')

    ## Generate a list of ref that will need to be modified
    #id_list = {hit.s_id:0 for hit in hit_list}.keys()

    ## Iterate over record in the subject fasta file
    #print ("Masking hit positions and writting a new reference for {} ".format(ref_outname))
    #i=j=0
    #start_time = time()
    #for record in SeqIO.parse(in_handle, "fasta"):
        ## Progress Marker
        #stdout.write("*")
        #stdout.flush()

        ## Check if the record is in the list of record to modify
        #if record.id in id_list:
            #i+=1
            ##~print ("Hit found in {}. Editing the sequence".format(record.id))
            ## Casting Seq type to MutableSeq Type to allow string editing
            #record.seq = record.seq.tomutable()

            ## For each hit in the list of hit found
            #for hit in hit_list:
                #if record.id == hit.s_id:

                    ## For all position between start and end coordinates modify the base by N
                    #for position in range (hit.s_start, hit.s_end):
                        #record.seq[position]= 'n'
        #else:
            #j+=1
            ##~print ("No hit found in {}".format(record.id))

        ## Finally write the sequence modified or not
        #out_handle.write(record.format("fasta"))
    #print("")
    ## Report informations
    #print("{} sequence(s) from {} modified in {}s".format(i,ref_outname, round(time()-start_time),2))

    ## Close files and return the masked ref path
    #in_handle.close()
    #out_handle.close()
    #return ref_path
