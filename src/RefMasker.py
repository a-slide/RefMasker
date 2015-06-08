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
    from collections import OrderedDict
    from datetime import datetime

    # Third party import
    import pyfasta # mandatory for fasta reading in the reference class

    # Local imports
    from FileUtils import is_readable_file, rm_blank
    from Conf_file import write_example_conf
    from Reference import Reference
    from pyBlast.BlastHit import BlastHit
    from pyBlast.Blastn import Blastn

except ImportError as E:
    print (E)
    print ("Please verify your dependencies. See Readme for more informations\n")
    exit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class RefMasker(object):
    """
    Main class of RefMasker program. Can output a model configuration file OR parse a filed
    configuration file, load and prepare References, perform iterative blast starting from the
    last reference against all others, then the penultimate against References listed before, and
    so one until there is only 1 reference. For References in which at least one read was found
    positions overlapped by blast hits are hard masked with *N* and written in a new fasta file.
    Finally, CSV reports are generated and the temporary files are deleted.
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

        print("Initialize RefMasker")
        # Parse the configuration file and verify the values of variables
        try:

            # Verify if conf file was given and is valid
            assert conf_file, "A path to the configuration file is mandatory"
            assert is_readable_file(conf_file), "{} is not a valid file".format(conf_file)
            self.conf = conf_file

            # Define a configuration file parser object and load the configuration file
            cp = ConfigParser.RawConfigParser(allow_no_value=True)
            cp.read(self.conf)

            print(" * Parse output options")
            # Output parameters section
            self.summary_report = cp.getboolean("Output", "summary_report")
            self.detailed_report = cp.getboolean("Output", "detailed_report")
            self.compress_output = cp.getboolean("Output", "compress_output")

            print(" * Parse Blast options")
            # Blast parameters section
            self.blastn_exec = cp.get("Blast", "blastn_exec")
            self.makeblastdb_exec = cp.get("Blast", "makeblastdb_exec")
            self.blast_task = cp.get("Blast", "blast_task")
            self.best_query_hit = cp.getboolean("Blast", "best_query_hit")
            self.evalue = cp.getfloat("Blast", "evalue")
            assert self.evalue > 0, "Authorized values for evalue: float > 0"

            print(" * Parse Reference sequences")
            # Iterate only on sections starting by "reference", create Reference objects
            # And store them in a list
            self.reference_list = []
            for reference in [i for i in cp.sections() if i.startswith("reference")]:
                # Create Reference objects
                self.reference_list.append (
                    Reference (
                        name = rm_blank(cp.get(reference, "name"), replace ='_'),
                        fasta = rm_blank(cp.get(reference, "fasta"), replace ='\ '),
                        compress = self.compress_output))

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
        print ("\nStart to process files")
        # Iterate over index in Reference.instances staring by the last one until the 2nd one

        try:
            for i in range(len(self.reference_list)-1, 0, -1):
                subject = self.reference_list[i]
                query_list = self.reference_list[0:i]

                print ("\nProcessing Reference \"{}\"".format(subject.name))

                # Create a blast database for the current subject sequence
                with Blastn(ref_path=subject.fasta, makeblastdb_exec=self.makeblastdb_exec) as blastn:

                    # Blast each query file of the query list against the subject
                    for query in query_list:
                        print (" * Blast against \"{}\"".format(query.name))

                        # Save the list of hit in a local variable
                        hit_list = blastn (
                            query_path = query.fasta,
                            blastn_exec = self.blastn_exec,
                            task = self.blast_task,
                            evalue = self.evalue,
                            best_query_hit = self.best_query_hit)

                        # Add the hit of list found to the subject
                        if hit_list:
                            print("   * {} hit(s) found".format(len(hit_list)))
                            subject.add_hit_list(hit_list)

                        else:
                            print ("   * No hit found")

                # if hits were found output the new fasta file in the current folder
                if subject.n_hit:
                    print (" * Write a modified reference fasta file in the current directory")
                    subject.output_reference ()
                else:
                    print (" * Reference file unmodified")

            # Write reports if requested
            if self.summary_report:
                print ("\nGenerate a summary report")
                with open ("Summary_report.csv", "w") as report:
                    report.write ("Program {}\tDate {}\n\n".format(self.VERSION,str(datetime.today())))
                    for ref in self.reference_list:
                        report.write(self._dict_to_report(ref.get_report(full=False)))
                        report.write("\n")

            if self.detailed_report:
                with open ("Detailed_report.csv", "w") as report:
                    print ("\nGenerate a detailed report")
                    report.write ("Program {}\tDate {}\n\n".format(self.VERSION,str(datetime.today())))
                    for ref in self.reference_list:
                        report.write(self._dict_to_report(ref.get_report(full=True)))
                        report.write("\n")

        # Catch possible exceptions
        except Exception as E:
            print ("ERROR during execution of RefMasker")
            print (E.message)

        # Even in case of exception this block will  be executed to remove temporary files
        finally:
            print ("\nCleanup temporary files")
            for ref in self.reference_list:
                ref.clean()

            print ("\nDone in {}s".format(round(time()-start_time, 3)))
            return(0)

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _dict_to_report(self, d, tab=""):
        """
        Recursive function to return a text report from nested dict or OrderedDict objects
        """
        report = ""
        for name, value in d.items():
            if type(value) == OrderedDict or type(value) == dict:
                report += "{}{}\n".format(tab, name)
                report += self._dict_to_report(value, tab=tab+"\t")
            else:
                report += "{}{}\t{}\n".format(tab, name, value)
        return report

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    refmasker = RefMasker.class_init()
    refmasker()
