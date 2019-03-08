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
from scipy import stats
import os
import logging

###############################################################################
class Tester:
	
	def binom(self, forward_count, reverse_count, cutoff, null):
		'''
		Apply a binomial distribution test to the forward and reverse
		read counts to determine whether there is directionality.

		Parameters
		----------
		forward_count	- Float. Number of forward reads mapping in the 
						  forward direction to a gene feature
		reverse_count	- Float. Number of forward reads mapping in the 
						  reverse direction to a gene feature
		cutoff			- Float. Significance cutoff for a binomial 
						  distribution test.
		null			- String. The null hypothesis of the test (i.e.
						  two-sided, greater, or less)

		Output
		------
		True if significant (p-value equal to or less than the cutoff), 
		otherwise false.
		'''
		counts = [forward_count, reverse_count]
		p_value = stats.binom_test( counts, alternative = null )
		if p_value <= cutoff:
			rg = str(max(counts) - min(counts))
			return p_value, rg
		else:
			return False
	def kolmogorov_smirnov(self, counts):
		'''
		description
		
		Inputs
		------
		
		Outputs
		-------
		
		'''
		if any(counts):
			d, pvalue = stats.kstest(counts, 'norm')
			return d, pvalue
		else:
			logging.warning('No results for Kolmogorov-Smirnov test')
			return 'NA', 'NA'
			
