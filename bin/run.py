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

# System imports
import logging
import os
import shutil
import time
import sys

# Local imports 
from tpmgenerator import TPMGenerator 
from tester import Tester
from dirseq import DirSeq

###############################################################################

debug={1:logging.CRITICAL,
       2:logging.ERROR,
       3:logging.WARNING,
       4:logging.INFO,
       5:logging.DEBUG}

###############################################################################

class Run:
	FEATURE_TYPE_CDS 	= b'CDS'
	COUNT_TYPE_COVERAGE = 'coverage'
	COUNT_TYPE_COUNT 	= 'count'

	NULL_TYPE_TWO_SIDED	= 'two-sided'
	NULL_TYPE_GREATER	= 'greater'
	NULL_TYPE_LESS 		= 'less'
	COUNT_TYPES 		= [COUNT_TYPE_COUNT,
	  			   		   COUNT_TYPE_COVERAGE]
	FEATURE_TYPES 		= [FEATURE_TYPE_CDS]
	NULL_TYPES 			= [NULL_TYPE_TWO_SIDED, NULL_TYPE_GREATER, NULL_TYPE_LESS]

	# Output files
	DIST_FILE			= 'detectm_quality.tsv'
	TPM_FILE 			= 'detectm_tpm_scores.tsv'

	def _logging_setup(self, args):
		'''
		Set up the logging output to stdout and a logging file in the output
		directory.
		
		Inputs
		------
		args    - object. Argparse object
		'''

		logger = logging.getLogger('')
		logger.setLevel(debug[args.verbosity])
		log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
									   datefmt="%Y-%m-%d %H:%M:%S %p")

		stream_logger = logging.StreamHandler(sys.stdout)
		stream_logger.setFormatter(log_format)
		stream_logger.setLevel(debug[args.verbosity])
		logger.addHandler(stream_logger)
		file_logger = logging.FileHandler(os.path.join(args.output_directory, args.log), 'a')
		file_logger.setFormatter(log_format)
		logger.addHandler(file_logger)

	def _check_general(self, args):
		'''
		Check general input and output options are valid.

		Parameters
		----------
		args    - object. Argparse object
		'''

		# Coverage option is deprecated
		if args.measure_type==self.COUNT_TYPE_COVERAGE:
			raise Exception("'coverage' measure type has been deprecated as of 0.0.3")
		# Check logging level is valid
		if args.verbosity not in range(1, 6):
			raise Exception("Logging verbosity must be a positive integer between 1 and 5.")
		
		# Set a working directory name, if none is specified
		if not args.output_directory:
			args.output_directory = '%s-detectm_output' % (time.strftime("%Y-%m-%d_%H-%M"))

		# Create working directory
		if(os.path.isdir(args.output_directory) or os.path.isfile(args.output_directory)):

			if args.force:

				if os.path.isdir(args.output_directory):
					shutil.rmtree(args.output_directory)

				else:
					os.remove(args.output_directory)

			else:
				raise Exception("File '%s' exists." % args.output_directory)
		os.mkdir(args.output_directory)  

		if(args.ignore_directions and 
		   args.measure_type!=self.COUNT_TYPE_COVERAGE):
			raise Exception("ignore_directions and count_type != coverage is currently unsupported")

	def main(self, args, command):

		self._check_general(args)
		self._logging_setup(args)
		
		logging.info("Running command: %s" % ' '.join(command))
		
		ds = DirSeq()
		tg = TPMGenerator()
		
		logging.info('Running DirSeq component')
		
		dirseq_output_lines = ds.main(args.bam,
									  args.gff,
									  args.forward_reads_only,
									  args.ignore_directions,
									  args.measure_type,
									  args.accepted_feature_types,
									  args.cutoff,
									  args.null,
									  os.path.join(args.output_directory, self.DIST_FILE))
		
		logging.info('Calculating TPM scores')

		tg.main(dirseq_output_lines, # Dirseq output
				args.rl,
                args.bam,
				os.path.join(args.output_directory, self.TPM_FILE))

		logging.info('Done')
