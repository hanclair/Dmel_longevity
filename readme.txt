Site polarization ancestral state inference in three steps.  

Requirements:

1. A tab-delimited input file, which is a table of focal SNPs in D. melanogaster in the following format
CHR	START	END	QVAL LTR	CHR OVERLAP	START GENE	END GENE	GENE ID	GENE NAME	NT OVERLAP
chr2L	499719	499720	0.038673216	chr2L	476220	540560	FBgn0003963	ush	1
chr2L	645185	645186	5.79E-07	chr2L	640013	714983	FBgn0284247	ds	1
chr2L	838002	838003	0.026840972	chr2L	833584	851071	FBgn0020304	drongo	1

Note- only the first 3 columns are required. 

2. BLAST and MUSCLE executables in the following directories:
	- /MUSCLE/muscle3.8.31_i86darwin64
	- /blast/blastn

3. Genome reference sequences that have have been BLAST formatted . You need both an uncompressed (for BLAST) and a gzipped version of a fasta file. Mom OS X requires that gzipped files have a ".Z" at the end. Annoying I know.  Examples:
	- dmel-all-chromosome-r6.27_w60.fasta.gz
	- dmel-all-chromosome-r6.27_w60.fasta.gz.Z
	- Dsim_w501.fasta.gz
	- Dsim_w501.fasta.gz.Z
	- Dyak_Tai18E2.fasta
	- Dyak_Tai18E2.fasta.gz.Z
The Dmel reference used in this paper is r6.27 from Flybase. The Dsim w501 and Dyak Tai18E2 references are unpublished at this time. 

4. Five perl scripts
	- get_Dmel_SNP_regions_in_Dsim_gz.pl
	- get_Dmel_SNP_regions_in_Dyak_gz.pl
	- get_Dsim_SNP_regions_list_gz.pl
	- get_Dyak_SNP_regions_list_gz.pl
	- snp_ancestral_state.pl


Step 1: takes a list of SNPs and uses blast to generate a bedfile coordinates in Dsim and Dyak.
fn="SNPlist.txt"
zcat dmel-all-chromosome-r6.27_w60.fasta.gz.Z | perl get_Dmel_SNP_regions_in_Dsim_gz.pl $fn 50
zcat dmel-all-chromosome-r6.27_w60.fasta.gz.Z | perl get_Dmel_SNP_regions_in_Dyak_gz.pl $fn 50


Step 2: Use the above bedfiles to extract fasta sequences from Dsim and Dyak.
zcat Dyak_Tai18E2.fasta.gz.Z | perl get_Dyak_SNP_regions_list_gz.pl Dyak_snp_coordinates.bed
zcat Dsim_w501_chromosome_arms_mtDNA_Pilon_w60.fasta.gz.Z | perl get_Dsim_SNP_regions_list_gz.pl Dsim_snp_coordinates.bed


Step 3. Use the above fasta to generate MUSCLE alignments and infer ancestral state at focal SNP
perl snp_ancestral_state.pl Dmel_snp_sequences.fa Dsim_snp_sequences.fa Dyak_snp_sequences.fa


Generates the following 9 files:

	check_blast_results_Dsim - logs blast results for each SNP in Dsim
	check_blast_results_Dyak - logs blast results for each SNP in Dyak

	Dsim_snp_coordinates.bed - a .bed file formatted summary of best blast hits coordinates in Dsim
	Dyak_snp_coordinates.bed - a .bed file formatted summary of best blast hits coordinates in Dyak

	Dmel_snp_sequences.fa - a fasta of Dmel regions surrounding focal SNP (each 100 bp with focal SNP at position 50)
	Dsim_snp_sequences.fa - a fasta of best blast hit Dsim sequences
	Dyak_snp_sequences.fa - a fasta of best blast hit Dyak sequences

	snp_ancestral_state_log -  detailed analysis of each SNP ancestry state. This includes the alignment of the three species, removing insertions in the Dmel sequence. A gap in the alignment is included just after position 50 in Dmel, this allows you to quickly blast the first segment and confirm that it ends on the focal SNP.  

	snp_ancestral_state_results - the main results, including states in mel, sim and yak (MSY) and the ancestral state call. Dmel SNP regions with no blast hit in either Dsim, dyak or both are logged as "no_hit". These columns should align with (and thus can be fused to) the input file. 

fn="SNPlist.txt"
perl remcarM.pl $fn
grep . $fn >temp1
paste temp1 snp_ancestral_state_results >"$fn"_results.xls
rm temp1

Caveats:
- currently hard-coded to blast 100bp window surrounding focal SNP, with focal SNP at position 50
- currently hard-coded for Dmel, Dyak, Dsim













