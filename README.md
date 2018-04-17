# metagg

The aim of our bioinformatics pipeline is to provide a reliable tool that assigns taxonomy and gene abundance both at community and taxon level to metagenomics data. The approach relies on assembling all the fastq file reads, DNA short sequences output of the sequencing technology, to create longer sequences (contigs) on which all DNA sequences can be mapped and quantified.  All the scripts are downloadable from this page and run in command line environment. 
We suggest executing the scripts
METAgg_mapping.pl
METAgg_annotation.pl
on the command line and then upload the XX folder, output of the first steps, on the AMAZON server, as the other steps of the pipeline have been optimized to be used on the AMAZON server as they involve an interactive visualisation and screening of taxonomic and functional profiles at community and taxon level. 
Here is a description of the script that are downloadable from this page:
1) METAgg_mapping.pl
uses the bwa algorithm to map the fastq file reads from all the samples to the contigs. In the script all the parameters are already set and optimized for this kind of data, but it is possible for the user to change some of those. The input files for this step are an assembly file and all the fastq  
2) METAgg_annotation.pl
uses a Diamond blastx search to align all the contigs against a non-redundant protein database (UniRef100). The output file reports only proteins that aligned to a contig with bit scores higher and e-value lower than a threshold. All the parameters are pre-set but the user can change them.
3) METAggS_gene.pl assigns protein coding regions to each contig. It parses the output of METAgg_annotation.pl and looks for non-overlapping proteins on each contig. Whenever more proteins matched the same contig region, only the protein with the highest bit score is kept. 
4) METAggS_lca.pl assigns taxonomy to each contig using the Lowest Common Ancestor (LCA) method. The LCA is determined looking at the taxon associated with each protein that has been assigned to each contig. The LCA assignment is also weighted with the bit score of each alignment. 
5) Once taxon (METAggS_lca.pl) and coding regions (METAggS_gene.pl) have been assigned to each contig, METAggS_assignment.pl uses this information and the output files of METAggS_mapping.pl to attribute taxon and gene relative abundances to all the samples. The relative abundances are calculated in reads per kilobase (RPK) and are weighted with the proportion of mapped reads in each sample.
6) METAggS_SNP.pl detects SNPs (Single Nucleotide Polymorphisms) in different samples. The output files report synonymous and non-synonymous SNPs, amino acid and codon substitutions, deletions and insertions. 
