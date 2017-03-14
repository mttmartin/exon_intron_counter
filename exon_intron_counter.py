#!/usr/bin/env python2.7
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import HTSeq
import collections
import argparse

def main():
    #gtf_filename = "/home/matthew/lab_root/sorghum/test_root/data/GTF/Sorghum_bicolor.Sorbi1.34.gtf"

    #alignment_filename = "/home/matthew/lab_root/sorghum/test_root/data/aligned_reads/20170209_1838_B!_light_20170209_1730_quality_filtered_minq20_minp80_R2Aligned.sortedByCoord.out.bam"

    parser = argparse.ArgumentParser()
    parser._action_groups.pop()
    required = parser.add_argument_group('required argumnets')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-g","--gtf",help="GTF file",required=True)
    optional.add_argument("-i","--input",help="Input BAM file(Default standard input)")
    optional.add_argument("-o","--output",help="Output file(Default standard out)")



    args = parser.parse_args()

    gtf_filename = args.gtf

    if args.output == None:
        output_file = sys.stdout
    else:
        output_filename = args.output
        try:
            output_file = open(output_filename,"w")
        except IOError:
            sys.stderr.write("Error opening output file\n")
            return

    if args.input == None:
        alignment_filename = sys.stdin
    else:
        alignment_filename = args.input

    gtf_file = HTSeq.GFF_Reader(gtf_filename)
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    introns = HTSeq.GenomicArrayOfSets("auto", stranded=True)

    gene_id_not_found_n = 0
    try:
        for feature in gtf_file:
            if "gene_id" not in feature.attr:
                gene_id_not_found_n += 1
                continue
            if feature.type == "exon":
                exons[feature.iv ] += feature.attr["gene_id"]
            elif feature.type == "gene":
                introns[feature.iv] += feature.attr["gene_id"]
    except IOError:
        sys.stderr.write("Error opening GTF file\n")
        return

    if gene_id_not_found_n > 0:
        sys.stderr("Warning: Could not find gene_id in " + str(gene_id_not_found_n) + " entries in GFF")

    exon_counts = collections.Counter()
    intron_counts = collections.Counter()
    almnt_file = HTSeq.BAM_Reader(alignment_filename)

    try:
        for almnt in almnt_file:
            if not almnt.aligned:
                exon_counts["_unmapped"] += 1
                intron_counts["_unmapped"] += 1
                continue
            gene_ids = set()

            for iv, val in exons[almnt.iv].steps():
                gene_ids |= val

            if len(gene_ids) == 1:
                gene_id = list(gene_ids)[0]
                exon_counts[gene_id] += 1
            elif len(gene_ids) == 0:
                intron_gene_ids = set()
                for iv,val in introns[almnt.iv].steps():
                    intron_gene_ids |= val
                if len(intron_gene_ids) == 1:
                    gene_id = list(intron_gene_ids)[0]
                    intron_counts[gene_id] += 1
                elif len(intron_gene_ids) == 0:
                    intron_counts["_no_feature"] += 1
                    exon_counts["_no_feature"] += 1
            else:
                intron_counts["_ambiguous"] += 1
                exon_counts["_ambiguous"] += 1
    except IOError:
        sys.stderr.write("Error opening input BAM file\n")
        return

    total_exons=0
    total_introns=0
    output_file.write("gene ID,exon count,intron count\n")
    for gene_id in exon_counts:
                 output_file.write( str(gene_id) + ',' + str(exon_counts[gene_id]) + ',' + str(intron_counts[gene_id]) + '\n')
                 if gene_id != "_no_feature" and gene_id != "_ambiguous":
                     total_exons += exon_counts[gene_id]
                     total_introns += intron_counts[gene_id]


if __name__ == "__main__":
    main()
