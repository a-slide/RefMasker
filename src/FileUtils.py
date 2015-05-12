# -*- coding: utf-8 -*-

"""
@package
@brief
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2015
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Standard library imports
from os import access, R_OK, remove, path
from gzip import open as gopen

#~~~~~~~ PREDICATES ~~~~~~~#

def is_readable_file (fp):
    """ Verify the readability of a file or list of file """
    if not access(fp, R_OK):
        raise IOError ("{} is not a valid file".format(fp))

def is_gziped (fp):
    """ Return True if the file is Gziped else False """
    return fp[-2:].lower() == "gz"

#~~~~~~~ PATH MANIPULATION ~~~~~~~#

def file_basepath (fp):
    """Return the path of a file without the last extension """
    return fp.rpartition('.')[0]

def file_basename (fp):
    """Return the basename of a file without folder location and extension """
    return fp.rpartition('/')[2].partition('.')[0]

def file_extension (fp):
    """ Return The extension of a file in lowercase """
    return fp.rpartition(".")[2].lower()

def file_name (fp):
    """ Return The complete name of a file with the extension but without folder location """
    return fp.rpartition("/")[2]

def dir_name (fp):
    """ Return the complete path where is located the file without the file name """
    return fp.rpartition("/")[0].rpartition("/")[2]

def rm_blank (name, replace=""):
    """ Replace blank spaces in a name by a given character (default = remove)
    Blanks at extremities are always removed and nor replaced """
    return replace.join(name.split())

#~~~~~~~ FILE MANIPULATION ~~~~~~~#

def gunzip_file (in_path):
    """
    @param in_path Path of the input compressed file
    @param out_path Path of the output uncompressed file (facultative)
    @exception  OSError Can be raise by open
    """
    # Generate a automatic name without .gz extension
    out_path = file_basepath(in_path)

    try:
        # Try to initialize handle for
        with gopen(in_path, 'rb') as in_handle:
            with open(out_path, "wb") as out_handle:
            # Write input file in output file
                out_handle.write (in_handle.read())

        return path.abspath(out_path)

    except IOError as E:
        print(E)
        if path.isfile (out_path):
            try:
                remove (out_path)
            except OSError:
                print "Can't remove {}".format(out_path)
