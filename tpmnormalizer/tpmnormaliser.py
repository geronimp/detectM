#!/usr/bin/env python
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
import argparse
import logging
import sys

# Local imports
from run import Run

###############################################################################

debug={1:logging.CRITICAL,
       2:logging.ERROR,
       3:logging.WARNING,
       4:logging.INFO,
       5:logging.DEBUG}

###############################################################################
      
class CustomHelpFormatter(argparse.HelpFormatter):
	
    def _split_lines(self, text, width):
        return text.splitlines()
    
description = u""" _______ _____  __  __                                    _ _               
|__   __|  __ \|  \/  |                                  | (_)              
   | |  | |__) | \  / |  _ __   ___  _ __ _ __ ___   __ _| |_ ___  ___ _ __ 
   | |  |  ___/| |\/| | | '_ \ / _ \| '__| '_ ` _ \ / _` | | / __|/ _ \ '__|
   | |  | |    | |  | | | | | | (_) | |  | | | | | | (_| | | \__ \  __/ |   
   |_|  |_|    |_|  |_| |_| |_|\___/|_|  |_| |_| |_|\__,_|_|_|___/\___|_|   
  ------------------------------------------------------------------------------------

  Authors: J. Boyd, B. Woodcroft

"""

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=description,
									 formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--bam_files', type=str, required=True, 
		help='Input alignment')
	parser.add_argument('--gff_files', type=str, required=True, 
		help='Output alignment')
	parser.add_argument('--forward_reads_only', action='store_true', default=False, 
		help='Only consider forward reads')
	parser.add_argument('--cutoff', type=float, default='0.05', 
		help='pvalue cutoff for the binomial test (Default = 0.05).')
	parser.add_argument('--null', type=str, default = 'two-sided', 
		help='Null hypothesis for the biomial test (two-sided, greater, or less; Default = two-sided)')
	parser.add_argument('--log', type=str, default=False,
		help='Output logging information to file')
	parser.add_argument('--verbosity', type=int, default=4,
		help='1 - 5, 1 being silent, 5 being noisy indeed. Default = 4')
	
	args = parser.parse_args()
	
	if args.log:
		if os.path.isfile(args.log):
			raise Exception("File %s exists" % args.log)
		logging.basicConfig(filename=args.log,
							level=debug[args.verbosity],
							format='%(asctime)s %(levelname)s: %(message)s',
							datefmt='%m/%d/%Y %I:%M:%S %p')
	else:
		logging.basicConfig(level=debug[args.verbosity],
							format='%(asctime)s %(levelname)s: %(message)s',
							datefmt='%m/%d/%Y %I:%M:%S %p')

	logging.info('Command: %s' % ' '.join(sys.argv))
	r = Run()
	r.main(args)
