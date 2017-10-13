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
__author__ 		= "Ben Woodcroft, Joel Boyd"
__copyright__ 	= "Copyright 2017"
__credits__ 	= ["Ben Woodcroft", "Joel Boyd"]
__license__ 	= "GPL3"
__maintainer__ 	= "Joel Boyd"
__email__ 		= "joel.boyd near uq.net.au"
__status__ 		= "Development"
__version__ 	= "0.0.1"
###############################################################################

# System imports
import tempfile
import os
import logging
import subprocess

from tester import Tester

###############################################################################

class DirSeq:
	'''
	This is a direct pythonic translation of Ben Woodcroft's DirSeq: 
		https://github.com/wwood/dirseq
	
	I've tried to be as faithful possible to the original, besides small 
	changes reflecting the differences between ruby and python. If there are 
	any bugs, they are probably the result of my lack of understanding of Ruby. 
	'''
	BAM_INDEX_SUFFIX = '.bai'
	STAR = '*'

	def _calculate_cov(self, covs, num_covs):
		return float(sum(covs)) / num_covs

	def _get_covs(self, cov_lines):
		feature_to_covs = {}
		previous_feature = None
		covs = []
		num_covs = 0
		for line in cov_lines:
			sline = line.split("\t")
			if sline[0]=='all': break
			feat = sline[8]
			if 'product' in sline[-2]:
				description = sline[-2].split(';')[-1].split('=')[1]
			else:
				description = 'None'

			if previous_feature:
				if feat != previous_feature:

					feature_to_covs[previous_feature] \
						= [sline[8][3:].split(';')[0], # Gene
						   sline[0], # Contig
						   sline[2], # Type
						   sline[3], # Start
						   sline[4], # Stop
						   sline[6], # Strand
						   str(self._calculate_cov(covs, num_covs)),
						   description]
			if len(sline) == 13:
				pass ### ~ TODO: Not implemented -hist
			elif len(sline) == 10:
				covs.append(int(sline[9]))
				num_covs += 1
			else:
				raise Exception("Unexpected bedtools output line: %s" % line)
			previous_feature=feat

		feature_to_covs[previous_feature] \
			= [sline[8][3:].split(';')[0], # Gene
						   sline[0], # Contig
						   sline[2], # Type
						   sline[3], # Start
						   sline[4], # Stop
						   sline[6], # Strand
						   str(self._calculate_cov(covs, num_covs)),
						   description]

		return feature_to_covs

	def _command_to_parsed(self, cmds):
		covs_initial = []
		for cmd in cmds:
			logging.info('Command: %s' % (cmd))
			covs_lines_initial = subprocess.check_output(cmd, shell=True).strip().split('\n')
			covs_initial.append(self._get_covs(covs_lines_initial))
		covs = covs_initial[0]
		if len(covs_initial) > 1:
			for key, coverage in covs_initial[1].items():
				coverage[6] = str(float(covs[key][6]) + float(coverage[6]))
				covs_initial[1][key] = coverage

		return covs

	def _write(self, covs_fwd, covs_rev, cutoff, null): 
		
		t=Tester()
		
		header = ['gene',
				  'contig',
				  'type',
				  'start',
				  'end',
				  'strand',
				  'forward_read_count',
				  'reverse_read_count',
				  'pvalue',
				  'normalized_read_count',
				  'annotation']
		print '\t'.join(header)

		for feature, forward_line in covs_fwd.iteritems():
			reverse_line = covs_rev[feature]
			forward_count = float(forward_line[6])
			reverse_count = float(reverse_line[6])
			result = t.binom(float(forward_count), float(reverse_count), cutoff, null)
			if result:
				pvalue, normalized_read_count = result
				if pvalue<=cutoff:
					output_line = forward_line[:6] + [str(forward_count), str(reverse_count), str(pvalue), str(normalized_read_count)] + [covs_fwd.values()[1][7]]
					print '\t'.join(output_line)

	def main(self, bam, gff, forward_reads_only, cutoff, null):
		
		nofastagff = tempfile.NamedTemporaryFile(suffix='.gff')
		cmd = "sed '/^##FASTA$/,$d' %s > %s" \
					% (gff, nofastagff.name)
		logging.info('Command: %s' % (cmd))
		subprocess.call(cmd, shell=True)


		bam_index = bam + self.BAM_INDEX_SUFFIX
		if not os.path.isfile(bam_index):
			raise Exception('Bam index file does not exist. Please index the BAM file')
		
		# Listing contigs in sorted order
		logging.info('Listing contigs in sorted order')
		bam_contigs = tempfile.NamedTemporaryFile(suffix='.tsv')
		cmd = "samtools idxstats %s | cut -f1,2 > %s" % (bam, bam_contigs.name)
		logging.info('Command: %s' % (cmd))
		subprocess.call(cmd, shell=True)

		logging.info('Finding featureless contigs')
		cmd = "grep -v '^#' %s | cut -f1 | sort | uniq | grep -vFw -f /dev/stdin %s | cut -f1" % (gff, bam_contigs.name)
		logging.info('Command: %s' % (cmd))
		featureless_contigs = [x for x in subprocess.check_output(cmd, shell=True).strip().split('\n')
				 			   if(x!=self.STAR and x!='')]
		logging.info('Found %i featureless contigs' % len(featureless_contigs))
		
		dummy_lines = []
		for featureless_contig in featureless_contigs:
			dummy_lines.append([featureless_contig,
								'dirseq',
								'misc_RNA',
								'1',
								'2',
								'.',
								'+',
								'0',
								"ID=%s_dummy_feature" % featureless_contig])
		
		sorted_gff_file = tempfile.NamedTemporaryFile(suffix='.gff')
		with tempfile.NamedTemporaryFile(suffix='.gff') as extra_features_file:
			for dummy_line in dummy_lines:
				extra_features_file.write('\t'.join(dummy_line) + '\n')
			extra_features_file.flush()
			
			cmd = "cat %s %s | bedtools sort -i /dev/stdin -faidx %s > %s" % (extra_features_file.name, nofastagff.name, bam_contigs.name, sorted_gff_file.name)
			logging.info('Command: %s' % (cmd))
			subprocess.call(cmd, shell=True)

		read1_flag = '-F128' #account for read1 in pair, as well as single reads mapping
		read2_flag = '-f128'

		cmdf1 = "samtools view -u %s %s | bedtools coverage -sorted -g %s -b /dev/stdin -a %s -s -counts" % (read1_flag, bam, bam_contigs.name, sorted_gff_file.name)
		cmdf2 = "samtools view -u %s %s | bedtools coverage -sorted -g %s -b /dev/stdin -a %s -s -counts" % (read2_flag, bam, bam_contigs.name, sorted_gff_file.name)
		cmdr1 = "samtools view -u %s %s | bedtools coverage -sorted -g %s -b /dev/stdin -a %s -S -counts" % (read1_flag, bam, bam_contigs.name, sorted_gff_file.name)
		cmdr2 = "samtools view -u %s %s | bedtools coverage -sorted -g %s -b /dev/stdin -a %s -S -counts" % (read2_flag, bam, bam_contigs.name, sorted_gff_file.name)
		
		if forward_reads_only:
			commands_fwd = [cmdf1]
			commands_rev = [cmdr1]
		else:
			commands_fwd = [cmdf1, cmdf2]
			commands_rev = [cmdr1, cmdr2]

		covs_fwd = self._command_to_parsed(commands_fwd)
		covs_rev = self._command_to_parsed(commands_rev)

		self._write(covs_fwd, covs_rev, cutoff, null)

		nofastagff.close()
		bam_contigs.close()
		sorted_gff_file.close()