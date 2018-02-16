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

    def get_t(self, t_list, rl):
        T = 0
        for feature in t_list:
            rg, flg = feature
            T += (float(rg)*float(rl))/float(flg)
        return T

    def main(self, dirseq_output, rl, output_file):

        t_list = []
        
        for idx, line in enumerate(dirseq_output):
            if idx == 0:
                header = line[:10] + ['TPM', 'feature_length'] + [line[10]]
                continue

            gene, contig, type, start, end, strand, rg_f, rg_r, pvalue, rg, directionality, annotation \
                = line
            flg \
                = float(end) - float(start) 
            t_list.append([float(rg), float(flg)])
        
        T = self.get_t(t_list, rl)

        with open(output_file, 'w') as out_io:
            # Print header
            out_io.write('\t'.join(header) + '\n')
            for dirseq_line, t_line in zip(dirseq_output[1:], t_list):
                if dirseq_line[0]=='gene':continue
                flg     = float(dirseq_line[4]) - float(dirseq_line[3]) 
                top     = float(dirseq_line[9]) * rl * 1e6
                bottom  = flg * T
                tpm     = top / bottom
                if(tpm>0):
                    output_line = dirseq_line[:10] + [str(tpm) , str(t_line[1]), dirseq_line[10]]
                    # Print output line
                    out_io.write('\t'.join(output_line) + '\n')
