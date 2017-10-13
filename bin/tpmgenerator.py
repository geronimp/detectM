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

__author__      = "Joel Boyd"
__copyright__   = "Copyright 2017"
__credits__     = ["Joel Boyd"]
__license__     = "GPL3"
__maintainer__  = "Joel Boyd"
__email__       = "joel.boyd near uq.net.au"
__status__      = "Development"
__version__     = "0.0.1"

###############################################################################
#       rg x rl x 10^6
# TPM = --------------
#         flg x T
#
# rg:    reads mapped to gene g
# rl:    read length
# flg:   feature length
# T:     sum of rgxrl/flg for all genes

class TPMGenerator:

    def parse_read_lengths(self, lengths):
        d = {}
        for line in open(lengths):
            sample, rl = line.strip().split('\t')
            d[sample] = float(rl)
        return d

    def get_t(self, feature_dict, rl):
        
        T_f = 0
        T_r = 0
        T = 0

        forward_all = float(sum([x[0] for x in feature_dict.values()]))
        reverse_all = float(sum([x[1] for x in feature_dict.values()]))

        for gene, features in feature_dict.items():
            T_f += (features[0]*rl)/features[2]
            T_r += (features[1]*rl)/features[2]
            T += ((features[1] + features[0])*rl)/features[2]
            features.append( str(features[0] / forward_all) )
            features.append( str(features[1] / reverse_all) )
            feature_dict[gene] = features

        return T_f, T_r, T, feature_dict

    def main(self, dirseq_file, lengths, sample):

        read_lengths    = self.parse_read_lengths(lengths)
        rl              = read_lengths[sample]
        dirseq_io       = open(dirseq_file)
        headers         = dirseq_io.readline().split()
        feature_dict    = {}
        
        for idx, line in enumerate(dirseq_io):
            
            contig, type, start, end, strand, rg_f, rg_r, annotation \
                = line.strip().split('\t')
            genome \
                = '_'.join(contig.split('_')[:2])
            flg \
                = float(end) - float(start) 
            feature_dict[contig + '_' + str(idx)] \
                = [float(rg_f), float(rg_r), float(flg), 
                   annotation, start, end, contig, genome]
        
        T_f, T_r, T, feature_dict = self.get_t(feature_dict, rl)

        # Print header
        print '\t'.join(['gene_id',
                         'description',
                         'start_pos',
                         'end_pos',
                         'contig',
                         'genome',
                         'sample', 
                         'count_forward',
                         'count_reverse',
                         'perc_forward',
                         'perc_reverse',
                         'TPM_forward_specific',
                         'TPM_reverse_specific',
                         'TPM_forward',
                         'TPM_reverse',
                         'directionality',
                         'group'])

        
        for gene, features in feature_dict.items():
            tpm_f_specific = (features[0]*rl*1e6)/(features[2]*T_f)
            tpm_r_specific = (features[1]*rl*1e6)/(features[2]*T_r)
            tpm_f = (features[0]*rl*1e6)/(features[2]*T)
            tpm_r = (features[1]*rl*1e6)/(features[2]*T)
            if(tpm_r>0 or tpm_f>0):    
                directionality = (tpm_f/sum([tpm_f, tpm_r])) * 100

                for i in range(len(bin_list)):
                    if(directionality>bin_list[i] and directionality<=bin_list[i+1]):
                        group = str(bin_list[i]) +'-'+ str(bin_list[i+1])
                        break
                    elif directionality==0.0:
                        group='0-5'
                        break
                # Print output line
                print '\t'.join([gene,
                                 features[3],
                                 features[4],
                                 features[5],
                                 features[6],
                                 features[7],
                                 sample,
                                 str(features[0]),
                                 str(features[1]),
                                 features[8],
                                 features[9],
                                 str(tpm_f_specific),
                                 str(tpm_r_specific),
                                 str(tpm_f),
                                 str(tpm_r),
                                 str(directionality),
                                 group])

