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
__author__ = "Ben Woodcroft, Joel Boyd"
__copyright__ = "Copyright 2017"
__credits__ = ["Ben Woodcroft", "Joel Boyd"]
__license__ = "GPL3"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"
__version__ = "0.0.1"

import tempfile
import os
import logging
import subprocess
import csv
import io
from tester import Tester
from collections import defaultdict

class DirSeq:
    '''
    This is a direct pythonic translation of Ben Woodcroft's DirSeq: 
        https://github.com/wwood/dirseq
    
    I've tried to be as faithful possible to the original, besides small 
    changes due to the differences between ruby and python. If there are 
    any bugs, they are probably the result of my lack of understanding of Ruby. 
    '''
    COUNT_TYPE_COVERAGE = 'coverage'
    COUNT_TYPE_COUNT = 'count'
    BAM_INDEX_SUFFIX = '.bai'
    STAR = b'*'

    def _calculate_cov(self, covs, num_covs):
        return float(sum(covs)) / num_covs

    def _get_covs(self, cov_lines, accepted_feature_types):

        removed 			= set()
        reader 				= cov_lines.split(b'\n')
        keys 				= [x.split(b'\t')[8] for x in reader if not x.split(b'\t')[0]=="all"]
        total 				= float(len(keys))
        feature_to_covs 	= defaultdict.fromkeys(keys)

        for idx, line in enumerate(reader):
            sline = line.split(b'\t')
            logging.debug('Parsing output [%f%% complete]' % ((idx / total)*100))

            if sline[0]=='all': break

            feat = sline[8]
            feature_type = sline[2]

            if feature_type not in accepted_feature_types:

                if feat not in removed:
                    del feature_to_covs[feat]
                    removed.add(feat)
                
                logging.debug('Skipping feature as it is of type %s' % feature_type)
                continue
            if(b'product' in sline[8] or b'annotations' in sline[8]):
                description = sline[8].split(b';')[-1].split(b'=')[1]
            else:
                description = 'None'
            

            if len(sline) == 13:
                num = int(sline[10])
                cov = num*int(sline[9])
            
            elif len(sline) == 10:
                cov = int(sline[9])
            
            else:
                raise Exception("Unexpected bedtools output line: %s" % line)

            if(feature_to_covs[feat]!=None):
                feature_to_covs[feat][6].append(cov)
            
            else:
                feature_to_covs[feat] \
                    = [sline[8].split(b';')[0].split(b'=')[1], # Gene
                       sline[0], # Contig
                       sline[2], # Type
                       sline[3], # Start
                       sline[4], # Stop
                       sline[6], # Strand
                       [cov],
                       description]
        
        logging.info('Parsing output [%f%% complete]' % ((total / total)*100))
        
        for key in set(keys):
        
            if key in feature_to_covs:
                feature_to_covs[key][6] = self._calculate_cov(feature_to_covs[key][6], len(feature_to_covs[key][6]))
        
        return feature_to_covs

    def _command_to_parsed(self, cmds, accepted_feature_types):
        
        covs_initial = list()

        for cmd in cmds:
            logging.info('Command: %s' % (cmd))
            covs_lines_initial = subprocess.check_output(cmd, shell=True).strip()

            covs_initial.append(self._get_covs(covs_lines_initial, accepted_feature_types))

        pooled_covs = dict()

        if len(covs_initial) > 1:

            for key, entry in covs_initial[0].items():
                pooled_coverage = str(float(covs_initial[1][key][6]) + float(entry[6]))	

                pooled_covs[key] = [covs_initial[1][key][0],
                                    covs_initial[1][key][1],
                                    covs_initial[1][key][2],
                                    covs_initial[1][key][3],
                                    covs_initial[1][key][4],
                                    covs_initial[1][key][5],
                                    str(pooled_coverage),
                                    covs_initial[1][key][7]]
        else:
            pooled_covs = covs_initial[0]

        return pooled_covs

    def _compile_output(self, covs_fwd, covs_rev, cutoff, null, 
               ignore_directions, measure_type, distribution_output): 
        
        logging.info('Compiling results')
        t = Tester()
        total = float(len(covs_fwd))
        header = ['gene', 'contig', 'type', 'start', 'end', 'strand']
        directionality_list_all = list()
        directionality_list_sig = list()
        
        if ignore_directions:
            header.append('average_coverage')
        
        else:
        
            if measure_type == self.COUNT_TYPE_COVERAGE:
                header.append('forward_average_coverage')
                header.append('reverse_average_coverage')
        
            elif measure_type == self.COUNT_TYPE_COUNT:
                header.append('forward_read_count')
                header.append('reverse_read_count')
        
            header += ['pvalue', 'normalized_read_count', 'directionality']
        
        header.append('annotation')

        output_lines = list()
        output_lines.append(header)
        
        if ignore_directions:

            for feature, forward_line in covs_fwd.items():
                output_line = forward_line[:6] + [forward_line[6], feature]
                output_lines.append(output_line)
        else:
            for idx, (feature, forward_line) in enumerate(covs_fwd.items()):
                logging.debug('Preparing output [%f%% complete]' % ((idx/total)*100))
                reverse_line = covs_rev[feature]
                forward_count = float(forward_line[6])	
                reverse_count = float(reverse_line[6])
                count_sum = (forward_count + reverse_count)
            
                if count_sum>0:
                    directionality = forward_count / count_sum
                else:
                    directionality = 0
            
                directionality_list_all.append(directionality)
                result = t.binom(float(forward_count), float(reverse_count), cutoff, null)
                
                if result:
                    pvalue, normalized_read_count = result
            
                    if pvalue<=cutoff:
                        directionality_list_sig.append(directionality)
                        output_line = forward_line[:6] + [str(forward_count),
                                                          str(reverse_count),
                                                          str(pvalue),
                                                          str(normalized_read_count),
                                                          str(directionality)] + [forward_line[-1]]
                        output_lines.append(output_line)
            
            directionality_all_result \
                    = t.kolmogorov_smirnov(directionality_list_all)
            directionality_sig_result \
                    = t.kolmogorov_smirnov(directionality_list_sig)

            self._distribution_output(distribution_output, directionality_all_result, directionality_sig_result)
        
        return output_lines
    
    def _distribution_output(self, distribution_output, directionality_all_result, directionality_sig_result):
        
        with open(distribution_output, 'w') as out_io:
            out_io.write('\t'.join(['Type', 'D', 'pvalue']) + '\n')
            out_io.write('All distribution result\t' + '\t'.join([str(x) for x in directionality_all_result]) + '\n')
            out_io.write('Passed distribution result\t' + '\t'.join([str(x) for x in directionality_sig_result]) + '\n')

    def main(self, bam, gff, forward_reads_only, ignore_directions, 
             measure_type, accepted_feature_types, cutoff, null,
             distribution_output):
        '''
        Run direq
        
        Inputs
        ------
        
        Outputs
        -------
        
        '''

        # Removing FASTA component from GFF file
        logging.info('Removing FASTA component from GFF file')
        nofastagff = tempfile.NamedTemporaryFile(suffix='.gff')
        cmd = "sed '/^##FASTA$/,$d' %s > %s" \
                    % (gff, nofastagff.name)
        logging.info('Command: %s' % (cmd))
        subprocess.call(cmd, shell=True)

        # Check for the existance of a BAM index file
        bam_index = bam + self.BAM_INDEX_SUFFIX
        
        if not os.path.isfile(bam_index):
            raise Exception('Bam index file does not exist. Please index the BAM file')
        
        # Listing contigs in sorted order
        logging.info('Listing contigs in sorted order')
        bam_contigs = tempfile.NamedTemporaryFile(suffix='.tsv')
        cmd = "samtools idxstats %s | cut -f1,2 | grep -v ^* > %s" % (bam, bam_contigs.name)
        logging.info('Command: %s' % (cmd))
        subprocess.call(cmd, shell=True)
        
        # Identifying contigs in feature file but not in BAM file
        logging.info('Finding contigs with no coverage')
        bam_contig_set = set()
        
        for line in open(bam_contigs.name):
            bam_contig_set.add(line.split()[0])
        
        gff_contig_set = list()
        lengths = dict()
        
        for line in open(gff):
        
            if line.startswith('##sequence-region'):
                _, contig, _, length = line.strip().split()
                gff_contig_set.append(contig)
                lengths[contig] = length
        coverageless_contigs = [contig for contig in gff_contig_set 
                                if contig not in bam_contig_set]
        
        # Identifying contigs that have no features in the GFF file
        logging.info('Finding featureless contigs')
        cmd = "grep -v '^#' %s | cut -f1 | sort | uniq | grep -vFw -f /dev/stdin %s | cut -f1" % (gff, bam_contigs.name)
        logging.info('Command: %s' % (cmd))

        featureless_contigs = [x for x in subprocess.check_output(cmd, shell=True).strip().split(str.encode('\n'))
                               if(x!=self.STAR and x!=b'')]
        logging.info('Found %i featureless contigs' % len(featureless_contigs))
        
        # Filling in missing contigs
        with open(bam_contigs.name, "a") as bam_contigs_io:
            
            for contig in coverageless_contigs:
                bam_contigs_io.write('\t'.join([contig, lengths[contig]]) + '\n')
        
        # Create dummy entries for featureless contigs and create a "full" GFF file
        dummy_lines = list()

        for featureless_contig in featureless_contigs:
            dummy_lines.append([featureless_contig,
                            b'dirseq',
                            b'misc_RNA',
                            b'1',
                            b'2',
                            b'.',
                            b'+',
                            b'0',
                            b"ID=%s_dummy_feature" % featureless_contig])

        sorted_gff_file = tempfile.NamedTemporaryFile(suffix='.gff')
        
        with tempfile.NamedTemporaryFile(suffix='.gff') as extra_features_file:
        
            for dummy_line in dummy_lines:
                extra_features_file.write(b'\t'.join(dummy_line) + b'\n')
            extra_features_file.flush()
            
            cmd = "cat %s %s | bedtools sort -i /dev/stdin -faidx %s > %s" % (extra_features_file.name, nofastagff.name, bam_contigs.name, sorted_gff_file.name)
            logging.info('Command: %s' % (cmd))
            subprocess.call(cmd, shell=True)

        if ignore_directions:
            cmd = "bedtools coverage -b %s -a %s -hist" % (bam, sorted_gff_file.name)
            logging.info('Command: %s' % cmd)
            cov_lines_fwd = subprocess.check_output(cmd, shell=True).strip()
            logging.info('Parsing coverage profiles')
            covs_fwd = self._get_covs(cov_lines_fwd, accepted_feature_types)
            covs_rev = None

        else:
            read1_flag = '-F128' #account for read1 in pair, as well as single reads mapping
            read2_flag = '-f128'

            if measure_type==self.COUNT_TYPE_COUNT:
                bedtools_type_flag = '-counts'
        
            elif measure_type==self.COUNT_TYPE_COVERAGE:
                bedtools_type_flag = '-hist'
        
            else:
                raise Exception("Measure type not recognised: %s" % (measure_type))

            cmdf1 = f"samtools view -u {read1_flag} {bam} | bedtools coverage -sorted -g {bam_contigs.name} -b /dev/stdin -a {sorted_gff_file.name} -s {bedtools_type_flag}"
            cmdf2 = f"samtools view -u {read2_flag} {bam} | bedtools coverage -sorted -g {bam_contigs.name} -b /dev/stdin -a {sorted_gff_file.name} -s {bedtools_type_flag}"
            cmdr1 = f"samtools view -u {read1_flag} {bam} | bedtools coverage -sorted -g {bam_contigs.name} -b /dev/stdin -a {sorted_gff_file.name} -S {bedtools_type_flag}"
            cmdr2 = f"samtools view -u {read2_flag} {bam} | bedtools coverage -sorted -g {bam_contigs.name} -b /dev/stdin -a {sorted_gff_file.name} -S {bedtools_type_flag}"
            
            if forward_reads_only:
                commands_fwd = [cmdf1]
                commands_rev = [cmdr1]

            else:
                commands_fwd = [cmdf1, cmdr2]
                commands_rev = [cmdf2, cmdr1]

            covs_fwd = self._command_to_parsed(commands_fwd, accepted_feature_types)
            covs_rev = self._command_to_parsed(commands_rev, accepted_feature_types)

        output_lines = self._compile_output(covs_fwd, covs_rev, cutoff, null, 
                                            ignore_directions, measure_type, 
                                            distribution_output)

        nofastagff.close()
        bam_contigs.close()
        sorted_gff_file.close()

        return output_lines
