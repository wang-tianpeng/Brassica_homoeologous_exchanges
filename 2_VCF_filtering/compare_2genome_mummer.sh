#!/bin/bash
set -ueo pipefail

##############
#
# run nucmer and delta-filter to compare 2 reference genome
# Author Tianpeng Wang 2022-05
#
#############

ref=$1
query=$2

name_ref=${ref/.fasta/}
name_ref=${name_ref/*\//}
name_query=${query/.fasta/}
name_query=${name_query/*\//}

name=${name_query}_2_${name_ref}
#echo $name
## run nucmer

nucmer -p $name $ref $query -g 1000 -c 90 -l 40 -t 60
delta-filter -r -q -l 5000 $name.delta >$name.delta.filter.rq
show-coords -c $name.delta.filter.rq >$name.delta.filter.rq.coords
#Rscript ~/scripts/mummerCoordsDotPlotly.tp.R -i ragtag.scaffold.zs11pb.contig.review.rename_2_B_rapa_ol_genome.delta.filter.rq.coords -o ragtag.scaffold.zs11pb.contig.review.rename_2_B_rapa_ol_genome.delta.filter.rq.coords.10k -m 10000 -l -x -r A01,A02,A03,A04,A05,A06,A07,A08,A09,A10,C01,C02,C03,C04,C05,C06,C07,C08,C09

delta-filter -m -l 5000 $name.delta >$name.delta.filter.m
show-coords -c $name.delta.filter.m >$name.delta.filter.m.coords
#Rscript ~/scripts/mummerCoordsDotPlotly.tp.R -i ragtag.scaffold.zs11pb.contig.review.rename_2_B_rapa_ol_genome.delta.filter.m.coords -o ragtag.scaffold.zs11pb.contig.review.rename_2_B_rapa_ol_genome.delta.filter.m.coords.10k -m 10000 -l -x -r A01,A02,A03,A04,A05,A06,A07,A08,A09,A10,C01,C02,C03,C04,C05,C06,C07,C08,C09
