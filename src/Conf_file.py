# -*- coding: utf-8 -*-

"""
@package    Refeed
@brief      Contain the template of the empty configuration file for Refeed
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

def write_example_conf():

    with open ("Refeed_conf_file.txt", 'wb') as fp:
        fp.write ("""
###################################################################################################
#                                   REFEED CONFIGURATION FILE                                      #
###################################################################################################
# Values can by customized with users values, but the file template must remain unchanged,
# otherwise the program will not be able to load default values.
# File path should be indicated as absolute path preferably and should not contain blank spaces
# Values identified with '**' in the descriptor are not recommended to be modified

###################################################################################################
[Output]

# Define the output options

# Output a simple and summary csv report at the end of execution
summary_report = True

# Output a detailed report including all blast hit coordinates and evalue
detailed_report = True

# Replace homologies in the reference with the query sequence instead of masking with "N"
# (incompatible with replace_homologies)(BOOLEAN)
mask_homologies : True

# Replace homologies in the reference with the query sequence instead of masking with the
# replacement character (incompatible with mask_homologies) (BOOLEAN)
replace_homologies : False

# Gzip fasta output (BOOLEAN)
compress_output : True

###################################################################################################
[Blast]

# Options to parameter blastn from blast+ used to find homologies between references

# Path to the blastn executable (not required is added to system path) (STRING) **
blastn_exec :

# Path to the mkblastdb executable (not required is added to system path) (STRING) **
makeblastdb_exec :

# Output only the best hit per query sequence in the query fasta file. Better to false
# if masking homologies (BOOLEAN) **
best_query_hit : False

# E Value cuttoff to retain alignments (FLOAT)
evalue : 0.1

# type of blast algorithm to perform 'blastn', 'blastn-short', 'dc-megablast', 'megablast',
# 'rmblastn'. Default = dc-megablast
blast_task = dc-megablast

###################################################################################################
# REFERENCE DEFINITION

# These sections corresponds to the definition of References sequence objects. It is possible to
# include as many independent reference as required by duplicating a entire Ref section.
# The order of references is CRITICAL. Homologies between the last reference and all the others
# reference above, will be masked (or replaced) in the last reference. Then, the homologies between
# the penultimate sequence and all the others references above (ie all except the last) will be
# masked, and so on until there is only one reference remaining
# All references have to be unique and to be organized as follow :
#   * [referenceX] = Reference identifier section, where X is the reference id starting from 1 for
#     the first reference and incrementing by 1 for each additional reference
#   * name = A name to define the reference (STRING)
#   * fasta = A path to a fasta file that can contain multiple sequences, can be gziped (STRING)

[reference1]
name : ref1
fasta : ../dataset/ref1.fa.gz

[reference2]
name : ref2
fasta : ../dataset/ref2.fa.gz

[reference3]
name : ref3
fasta : ../dataset/ref3.fa.gz

[reference4]
name : ref4
fasta : ../dataset/ref4.fa.gz

[reference5]
name : ref5
fasta : ../dataset/ref5.fa.gz

""")
