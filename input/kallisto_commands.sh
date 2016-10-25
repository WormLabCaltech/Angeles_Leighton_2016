#!/bin/bash
# make transcript index
kallisto index -i input/transcripts.idx input/c_elegans_WBcel235.rel79.cdna.all.fa;
# kallisto command for ../input/rawseq/JKY1
kallisto quant -i input/transcripts.idx -o input/kallisto_all/JKY1/kallisto --single -s 60 -l 180 -b 200 -t 4 input/rawseq/JKY1/15044_GTTTCG_L001_R1_001.fastq.gz input/rawseq/JKY1/15044_GTTTCG_L001_R1_002.fastq.gz input/rawseq/JKY1/15044_GTTTCG_L001_R1_003.fastq.gz input/rawseq/JKY1/15044_GTTTCG_L001_R1_004.fastq.gz input/rawseq/JKY1/15044_GTTTCG_L001_R1_005.fastq.gz input/rawseq/JKY1/15044_GTTTCG_L002_R1_001.fastq.gz input/rawseq/JKY1/15044_GTTTCG_L002_R1_002.fastq.gz input/rawseq/JKY1/15044_GTTTCG_L002_R1_003.fastq.gz input/rawseq/JKY1/15044_GTTTCG_L002_R1_004.fastq.gz input/rawseq/JKY1/15044_GTTTCG_L002_R1_005.fastq.gz;
# kallisto command for ../input/rawseq/JKY2
kallisto quant -i input/transcripts.idx -o input/kallisto_all/JKY2/kallisto --single -s 60 -l 180 -b 200 -t 4 input/rawseq/JKY2/15045_CGTACG_L001_R1_001.fastq.gz input/rawseq/JKY2/15045_CGTACG_L001_R1_002.fastq.gz input/rawseq/JKY2/15045_CGTACG_L001_R1_003.fastq.gz input/rawseq/JKY2/15045_CGTACG_L001_R1_004.fastq.gz input/rawseq/JKY2/15045_CGTACG_L001_R1_005.fastq.gz input/rawseq/JKY2/15045_CGTACG_L001_R1_006.fastq.gz input/rawseq/JKY2/15045_CGTACG_L002_R1_001.fastq.gz input/rawseq/JKY2/15045_CGTACG_L002_R1_002.fastq.gz input/rawseq/JKY2/15045_CGTACG_L002_R1_003.fastq.gz input/rawseq/JKY2/15045_CGTACG_L002_R1_004.fastq.gz input/rawseq/JKY2/15045_CGTACG_L002_R1_005.fastq.gz input/rawseq/JKY2/15045_CGTACG_L002_R1_006.fastq.gz;
# kallisto command for ../input/rawseq/JKY3
kallisto quant -i input/transcripts.idx -o input/kallisto_all/JKY3/kallisto --single -s 60 -l 180 -b 200 -t 4 input/rawseq/JKY3/15046_GAGTGG_L001_R1_001.fastq.gz input/rawseq/JKY3/15046_GAGTGG_L001_R1_002.fastq.gz input/rawseq/JKY3/15046_GAGTGG_L001_R1_003.fastq.gz input/rawseq/JKY3/15046_GAGTGG_L001_R1_004.fastq.gz input/rawseq/JKY3/15046_GAGTGG_L001_R1_005.fastq.gz input/rawseq/JKY3/15046_GAGTGG_L002_R1_001.fastq.gz input/rawseq/JKY3/15046_GAGTGG_L002_R1_002.fastq.gz input/rawseq/JKY3/15046_GAGTGG_L002_R1_003.fastq.gz input/rawseq/JKY3/15046_GAGTGG_L002_R1_004.fastq.gz input/rawseq/JKY3/15046_GAGTGG_L002_R1_005.fastq.gz;
# kallisto command for ../input/rawseq/N2Y1
kallisto quant -i input/transcripts.idx -o input/kallisto_all/N2Y1/kallisto --single -s 60 -l 180 -b 200 -t 4 input/rawseq/N2Y1/14945_CAGATC_L001_R1_001.fastq.gz input/rawseq/N2Y1/14945_CAGATC_L001_R1_002.fastq.gz input/rawseq/N2Y1/14945_CAGATC_L001_R1_003.fastq.gz input/rawseq/N2Y1/14945_CAGATC_L001_R1_004.fastq.gz input/rawseq/N2Y1/14945_CAGATC_L001_R1_005.fastq.gz input/rawseq/N2Y1/14945_CAGATC_L002_R1_001.fastq.gz input/rawseq/N2Y1/14945_CAGATC_L002_R1_002.fastq.gz input/rawseq/N2Y1/14945_CAGATC_L002_R1_003.fastq.gz input/rawseq/N2Y1/14945_CAGATC_L002_R1_004.fastq.gz input/rawseq/N2Y1/14945_CAGATC_L002_R1_005.fastq.gz;
# kallisto command for ../input/rawseq/N2Y2
kallisto quant -i input/transcripts.idx -o input/kallisto_all/N2Y2/kallisto --single -s 60 -l 180 -b 200 -t 4 input/rawseq/N2Y2/14946_ACTTGA_L001_R1_001.fastq.gz input/rawseq/N2Y2/14946_ACTTGA_L001_R1_002.fastq.gz input/rawseq/N2Y2/14946_ACTTGA_L001_R1_003.fastq.gz input/rawseq/N2Y2/14946_ACTTGA_L001_R1_004.fastq.gz input/rawseq/N2Y2/14946_ACTTGA_L001_R1_005.fastq.gz input/rawseq/N2Y2/14946_ACTTGA_L002_R1_001.fastq.gz input/rawseq/N2Y2/14946_ACTTGA_L002_R1_002.fastq.gz input/rawseq/N2Y2/14946_ACTTGA_L002_R1_003.fastq.gz input/rawseq/N2Y2/14946_ACTTGA_L002_R1_004.fastq.gz input/rawseq/N2Y2/14946_ACTTGA_L002_R1_005.fastq.gz;
# kallisto command for ../input/rawseq/N2Y3
kallisto quant -i input/transcripts.idx -o input/kallisto_all/N2Y3/kallisto --single -s 60 -l 180 -b 200 -t 4 input/rawseq/N2Y3/14947_GATCAG_L001_R1_001.fastq.gz input/rawseq/N2Y3/14947_GATCAG_L001_R1_002.fastq.gz input/rawseq/N2Y3/14947_GATCAG_L001_R1_003.fastq.gz input/rawseq/N2Y3/14947_GATCAG_L002_R1_001.fastq.gz input/rawseq/N2Y3/14947_GATCAG_L002_R1_002.fastq.gz input/rawseq/N2Y3/14947_GATCAG_L002_R1_003.fastq.gz;