###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ 		= "Joel Boyd"
__copyright__ 	= "Copyright 2017"
__credits__ 	= ["Joel Boyd"]
__license__ 	= "GPL3"
__maintainer__ 	= "Joel Boyd"
__email__ 		= "joel.boyd near uq.net.au"
__status__ 		= "Development"
__version__ 	= "0.0.1"

###############################################################################

# Imports
import logging

# Local imports 
from tpmgenerator import TPMGenerator 
from tester import Tester
from dirseq import DirSeq

###############################################################################

class Run:

	def _check_args(self, args):
		pass

	def main(self, args):
		self._check_args(args)
		logging.info('Running TPM normalisier')
		ds = DirSeq()
		ds.main(args.bam_files, args.gff_files, args.forward_reads_only, args.cutoff, args.null)

