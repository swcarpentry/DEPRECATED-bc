#! /bin/bash


for i in {1..5}
do
    bedtools shuffle -chrom -i cpg.bed -g genome.txt > shuffled_cpg/cpg.shuffled${i}.bed
done

