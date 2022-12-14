Workflow followed when processing bacterial V4 region sequencing data
<br> - raw sequences preprocess <br> - zOTU construction and taxonomy
assignment <br> - preparing for downstream analyses

<br><br><br>

# **1) Raw sequences preprocess & curation**

| Info              |                                                 |
|-------------------|-------------------------------------------------|
| Working dir (=WD) | `podruh2:/mnt/data/chomic/seqme2021/eroSOM/v4/` |
| Samples           | \# 1-200                                        |
| File name example | 1_R1.fq                                         |

### **1.1.) Paired-end merging and raw combined fastq construction**

-   merge paired-end reads (max 10 mismatches, min % id 90), rename seq
    header acc. to sample name

<!-- -->

    WD$
    usearch11 -fastq_mergepairs raw/*R1.fq -fastqout merged.fq -fastq_maxdiffs 10 -fastq_pctid 90 -relabel @ -report merged_log.txt -tabbedout merged_tab.txt

    Totals:
       5184872  Pairs (5.2M)
       4996540  Merged (5.0M, 96.37%)
       3351086  Alignments with zero diffs (64.63%)
        131475  Too many diffs (> 10) (2.54%)
         56857  No alignment found (1.10%)
             0  Alignment too short (< 16) (0.00%)
           388  Staggered pairs (0.01%) merged & trimmed
        198.16  Mean alignment length
        301.83  Mean merged length
          0.49  Mean fwd expected errors
          0.45  Mean rev expected errors
          0.14  Mean merged expected errors

-   cut adapters & primers

<!-- -->

    WD$
    cutadapt -g GTGYCAGCMGCCGCGGTAA --cores=150 --discard-untrimmed merged.fq | cutadapt -a ATTAGAWACCCBNGTAGTCC --cores=150 --discard-untrimmed - > noadapt.fq
     
    === Summary ===

    Total reads processed:               4,984,492
    Reads with adapters:                 4,971,594 (99.7%)
    Reads written (passing filters):     4,971,594 (99.7%)

    Total basepairs processed: 1,385,667,811 bp
    Total written (filtered):  1,258,671,910 bp (90.8%)

------------------------------------------------------------------------

### **1.2.) Quality curation and length filtering**

-   #### Discard reads \>259 bp, trim to length 252 bp, discard shorter reads, no ambiguous bases

<!-- -->

    WD$
    prinseq -fastq noadapt.fq -out_good l259 -out_bad null -max_len 259 && prinseq -fastq l259.fastq -out_good l252 -out_bad null -min_len 252 -trim_to_len 252 -ns_max_n 0

    Input and filter stats:
            Input sequences: 4,971,594
            Input bases: 1,258,671,910
            Input mean length: 253.17
            Good sequences: 4,942,362 (99.41%)
            Good bases: 1,250,724,283
            Good mean length: 253.06
            Bad sequences: 29,232 (0.59%)
            Bad bases: 7,947,627
            Bad mean length: 271.88
            Sequences filtered by specified parameters:
            max_len: 29232
    Input and filter stats:
            Input sequences: 4,942,362
            Input bases: 1,250,724,283
            Input mean length: 253.06
            Good sequences: 4,933,577 (99.82%)
            Good bases: 1,243,261,404
            Good mean length: 252.00
            Bad sequences: 8,785 (0.18%)
            Bad bases: 2,134,975
            Bad mean length: 243.03
            Sequences filtered by specified parameters:
            min_len: 7474
            ns_max_n: 1311

<br>

-   #### Discard reads with error probability \>1

<!-- -->

    WD$
    usearch11 -fastq_filter l252.fastq -fastq_maxee 1 -fastqout l252_e1.fastq -threads 150

       4933577  Reads (4.9M)
         47447  Discarded reads with expected errs > 1.00
       4886130  Filtered reads (4.9M, 99.0%)

------------------------------------------------------------------------

# **2) zOTUs generation**

-   USEARCH needs sequences in FASTA format

<!-- -->

    WD$
    prinseq -fastq l252_e1.fastq  -out_good l252_e1 -out_bad null -out_format 1 -line_width 0

-   Dereplication

<!-- -->

    WD$
    usearch11 -fastx_uniques l252_e1.fasta -fastaout uniq.fa -sizeout -relabel Uniq

    00:20 1.6Gb   100.0% Reading l252_e1.fasta
    00:20 1.6Gb  CPU has 192 cores, defaulting to 10 threads
    00:24 2.7Gb   100.0% DF
    00:24 2.8Gb  4886130 seqs, 1013541 uniques, 663918 singletons (65.5%)
    00:24 2.8Gb  Min size 1, median 1, max 47899, avg 4.82
    00:29 2.8Gb   100.0% Writing uniq.fa

-   zOTU construction

<!-- -->

    WD$
    usearch11 -unoise3 uniq.fa -zotus zotus.fa -minsize 8

    00:01 367Mb   100.0% Reading uniq.fa
    00:05 504Mb   100.0% 15762 amplicons, 465760 bad (size >= 8)
    05:37 518Mb   100.0% 15554 good, 208 chimeras
    05:37 518Mb   100.0% Writing zotus    

# **3) Taxonomy assignment**

-   BLAST against UNITE database

<!-- -->

    WD$
    parallel_assign_taxonomy_blast.py -i zotus.fa -o . -r /mnt/data/chomic/databases/silva_138_qiime_release/silva-138-ssu-nr99-seqs-515f-806r-majority.fasta -t /mnt/data/chomic/databases/silva_138_qiime_release/silva-138-ssu-nr99-tax-515f-806r-derep-majority.txt -O 190

-   Convert taxonomy syntax from Qiime to utax compatible
    -   in the otus file `zotus.fna`
        -   get rid of extra l within sequences
        -   sort sequences by OTU name
    -   within BLAST output `zotus_tax_assignments.txt`
        -   sort OTUs ascending (by OTU #)
        -   add empty rows to allow merging with `zotus.fna`, get rid of
            OTU \#
        -   reformat taxonomy syntax for utax; delete ???No blast hit???
            fields
    -   merge OTU sequences `zotus.fna` and UTAX comaptible taxonomy

<!-- -->

    WD$
    prinseq -fasta zotus.fa -out_bad null -out_good zotus_line -line_width 0

    cat zotus_line.fasta | paste - -| sed 's/>Zotu//1' | sort -n | sed 's/^/>ZOTU/g; s/\t/\n/g' > zotus_line_sorted.fa

    sed 's/Zotu//1' zotus_tax_assignments.txt | sort -g | sed 's/^/ZOTU/' > zotus_tax_assignments_sorted.txt

    sed -e 's/unidentified//g; s/d__/tax=d:/g; s/__/:/g; s/;/,/g; s/,p:,c:,o:,f:,g:,s://g; s/,c:,o:,f:,g:,s://g; s/,o:,f:,g:,s://g; s/,f:,g:,s://g; s/,g:,s://g; s/,s:$//g; s/tax=d/;tax=d/g; s/No blast hit/;tax=d:No_blast_hit/g; s/$/;/g; s/, /,/g' zotus_tax_assignments_sorted.txt > zotus_tax_assignments_sorted_reformat.txt

    sed -e 'G' zotus_tax_assignments_sorted_reformat.txt | cut -f 2 > blast_taxonomy.txt

    paste -d "" zotus_line_sorted.fa blast_taxonomy.txt > zotus_tax_blast.fna

# **4) OTU table construction**

-   map initial sequences to OTUs

<!-- -->

    WD$
    usearch8 -usearch_global l252_e1.fasta -db zotus_tax_blast.fna -strand plus -id 0.99 -otutabout otu_table.txt

    4647665 / 4886130 mapped to OTUs (95.1%)    

-   convert to phyloseq compatible taxonomy format

<!-- -->

    WD$
    cut -f 202- otu_table.txt > tax.txt
    cut -f 1-201 otu_table.txt > abund.txt

    sed -E 's/d://g; s/p://g; s/c://g; s/o://g; s/f://g; s/g://g; s/s://g' tax.txt| 
    sed 's/,uncultured.*$//g; s/,gamma_proteobacterium.*$//g; s/,_marine.*$//g; s/,metagenome.*$//g;
    s/,groundwater_metagenome.*$//g; s/,marine_metagenome.*$//g; s/,soil_.*$//g; s/,g:Unknown_Family.*$//g; s/,bacterium_enrichment.*$//g' | 
    awk -F"," '{ if ($7 ~ /_sp./) print $1","$2","$3","$4","$5","$6",u_"$6;  else print $0}' | 
    awk -F"," '{ if ($7 ~ /_bacterium/) print $1","$2","$3","$4","$5","$6",u_"$6;  else print $0}' | 
    awk -F"," '{ if ($6 ~ /ceae/) print $1","$2","$3","$4","$5",u_"$5;  else print $0}' | 
    awk -F"," '{ if ($5 ~ /ales/) print $1","$2","$3","$4",u_"$4;  else print $0}' | 
    awk -F"," '{ if ($2 == "") print $1",u_"$1;  else print $0}' | 
    awk -F"," '{ if ($3 == "") print $1","$2",u_"$2;  else print $0}' | 
    awk -F"," '{ if ($4 == "") print $1","$2","$3",u_"$3;  else print $0}' | 
    awk -F"," '{ if ($5 == "") print $1","$2","$3","$4",u_"$4;  else print $0}' | 
    awk -F"," '{ if ($6 == "") print $1","$2","$3","$4","$5",u_"$5;  else print $0}' | 
    awk -F"," '{ if ($7 == "") print $1","$2","$3","$4","$5","$6",u_"$6;  else print $0}' | 
    sed -E '1s/^.*$/taxonomy/; s/,/;/g' | sed 's/\(u_\)\1\{1,\}/u_/g'  > tax_phyloseq.txt

    paste abund.txt tax_phyloseq.txt > zotu_table_v4.txt

-   discard OTUs assigned to mitochondria or chloroplasts and without
    blast hit

<!-- -->

    WD$
    sed -i -E '/Mitochondria/d; /Chloroplast/d; /No_blast_hit/d' zotu_table_v4.txt

-   create and summarize biom file

<!-- -->

    WD$
    biom convert -i zotu_table_v4.txt --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy -o zotu_table_v4.biom

    biom summarize-table -i zotu_table_v4.biom -o zotu_table_v4_summary.txt

    Num samples: 200
    Num observations: 15197
    Total count: 4633752
    Table density (fraction of non-zero values): 0.181

    Counts/sample summary:
     Min: 7259.000
     Max: 40401.000
     Median: 22441.000
     Mean: 23168.770
     Std. dev.: 6126.075
     Sample Metadata Categories: None provided
     Observation Metadata Categories: taxonomy

    Counts/sample detail:
    4: 7259.000
    12: 8283.000
    76: 10727.000
    7: 11421.000
    10: 12767.000
    69: 13289.000
    6: 13355.000
    94: 13767.000
    52: 13810.000
    1: 13897.000
    153: 14093.000
    8: 14214.000
    91: 14277.000
    45: 14323.000
    129: 14427.000
    57: 14588.000
    5: 14609.000
    89: 14708.000
    9: 14882.000
    177: 14987.000
    164: 15238.000
    132: 15502.000
    130: 15519.000
    165: 15665.000
    11: 16004.000
    58: 16099.000
    2: 16192.000
    92: 16582.000
    21: 16692.000
    168: 16869.000
    125: 17044.000
    40: 17347.000
    163: 17371.000
    131: 17485.000
    162: 17530.000
    115: 17541.000
    44: 17561.000
    120: 17826.000
    60: 17858.000
    55: 17879.000
    96: 17970.000
    128: 18080.000
    156: 18141.000
    124: 18242.000
    86: 18258.000
    56: 18396.000
    127: 18660.000
    20: 18752.000
    18: 18994.000
    77: 19129.000
    70: 19146.000
    46: 19165.000
    160: 19169.000
    161: 19282.000
    24: 19309.000
    90: 19321.000
    154: 19328.000
    67: 19369.000
    48: 19391.000
    54: 19457.000
    88: 19525.000
    22: 19569.000
    167: 19611.000
    166: 19674.000
    180: 19745.000
    79: 19941.000
    65: 20072.000
    3: 20168.000
    141: 20188.000
    148: 20234.000
    152: 20338.000
    95: 20374.000
    174: 20439.000
    126: 20463.000
    66: 20552.000
    19: 20560.000
    93: 20679.000
    43: 20791.000
    80: 20834.000
    178: 20860.000
    33: 20934.000
    72: 20962.000
    158: 20978.000
    87: 21016.000
    123: 21230.000
    1182: 21279.000
    103: 21364.000
    117: 21416.000
    42: 21485.000
    53: 21549.000
    172: 21640.000
    17: 21737.000
    151: 21798.000
    139: 21870.000
    136: 22090.000
    28: 22133.000
    175: 22182.000
    81: 22262.000
    102: 22347.000
    16: 22395.000
    142: 22487.000
    122: 22936.000
    82: 23003.000
    101: 23050.000
    78: 23159.000
    41: 23271.000
    68: 23423.000
    155: 23608.000
    119: 23637.000
    50: 23652.000
    1195: 23744.000
    118: 23756.000
    150: 23901.000
    138: 24047.000
    140: 24056.000
    34: 24088.000
    75: 24305.000
    108: 24456.000
    179: 24531.000
    157: 24661.000
    64: 24799.000
    159: 24820.000
    1194: 24857.000
    116: 25007.000
    74: 25176.000
    121: 25280.000
    105: 25282.000
    144: 25291.000
    171: 25305.000
    51: 25585.000
    59: 25659.000
    36: 25687.000
    1196: 25700.000
    137: 25974.000
    170: 26276.000
    98: 26298.000
    1197: 26368.000
    173: 26406.000
    84: 26435.000
    112: 26558.000
    149: 26596.000
    38: 26663.000
    1193: 26768.000
    71: 26923.000
    100: 27044.000
    110: 27139.000
    30: 27353.000
    1199: 27377.000
    134: 27500.000
    104: 27542.000
    63: 27582.000
    135: 27778.000
    1192: 27886.000
    176: 27941.000
    47: 27980.000
    111: 28322.000
    146: 28406.000
    49: 28446.000
    1186: 28468.000
    114: 28604.000
    147: 28721.000
    182: 28730.000
    107: 28831.000
    1188: 29092.000
    1185: 29221.000
    35: 29399.000
    73: 29400.000
    37: 29427.000
    23: 29442.000
    62: 29563.000
    1189: 29587.000
    1198: 29605.000
    32: 29754.000
    39: 29830.000
    143: 30357.000
    25: 30657.000
    85: 30709.000
    106: 30909.000
    99: 30972.000
    97: 31028.000
    14: 31131.000
    109: 31219.000
    169: 31995.000
    13: 32375.000
    133: 32454.000
    113: 32808.000
    83: 32889.000
    29: 33628.000
    1191: 33876.000
    1190: 33883.000
    31: 33914.000
    1183: 34137.000
    61: 34152.000
    15: 34251.000
    26: 34310.000
    1184: 35830.000
    1187: 36210.000
    181: 37694.000
    27: 38181.000
    145: 40401.000

<!-- # **5) OTU functional potential** -->
<!-- ~~~ -->
<!-- WD/picrust -->
<!-- place_seqs.py -s ../zotu/zotus.fa -o out.tre -p 190 --verbose -->
<!-- hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p 190 -->
<!-- ~~~~ -->

------------------------------------------------------------------------

#### SOFTWARE USED

<!-- * gunzip (gzip) 1.6 -->

-   PRINSEQ-lite 0.20.4
-   QIIME 1.9.1 <!-- * FastQC v0.11.8 -->
    <!-- * multiqc, version 1.8 -->
-   sed (GNU sed) 4.4
-   usearch v11.0.667_i86linux64 (and usearch v8.1.1861_i86linux64)
-   GNU Awk 4.1.4, API: 1.1 (GNU MPFR 4.0.1, GNU MP 6.1.2)
-   cut (GNU coreutils) 8.28
-   paste (GNU coreutils) 8.28
-   biom, version 2.1.5 <!-- * PICRUSt 2.3.0-b -->
