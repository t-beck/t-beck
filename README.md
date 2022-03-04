_Primer Design Pipeline Use:_

The primary use of this program has been to compare _Alu_ loci via BLAT of user-specified genomes. There are several things to note before you begin:

1. You will need an input text file (coordinate_file) with one line for each locus of interest that includes an additional 600 bp on either side of the insertions. The flanking sequence is necessary to properly obtain and align orthologous sequence from the comparison genome assemblies. Each locus should be on one line, in the format chromosome:start-stop, e.g.:  
   \
    chr5:163793654-1637948654 </br>
    chr6:678090-679920907 </br>
    …….

2. Paths to several programs have been hard-coded. Therefore, these lines need to be changed before running this pipeline. The lines in each program are as follows:
   1. Line 153 of blatGenomes1.py currently indicates the muscle program is: muscle3.8.31_i86linux64. Not only should the MUSCLE alignment be in your PATH, but should also be the appropriate name for the version of MUSCLE you have installed on your computer.
   2. The current Primer3 parameters can be found on lines 86-99 of the program basetobasecompare.py. Any desired deviations from these can be edited by changing these lines.
   3. Line 102 requires the location of the Primer3 program as well as the Primer3 settings. These are currently set as, &#39;/home/batzer/Downloads/primer3-2.3.6/src/primer3_core&#39; and, &#39;-p3_settings_file=/home/batzer/Downloads/primer3-2.3.6/primer3web_v4_0_0_default_settings.txt&#39;, respectively. Please change these to the appropriate path for your computer.

The following dependencies are required prior to running this pipeline to ensure a successful run:

1. BLAT
2. Primer3
3. Python2
4. Linux/Unix OS
5. MUSCLE
6. 2bit genome assemblies. These can be generated by running the following on the command line:

   `$ faToTwoBit genome.fa genome.2bit`

7. An ooc file for each genome assembly. The ooc file tells BLAT the over-occuring k-mers in the genome assembly, which will increase the speed of the program. The default size for k is 11. This file can be generated by running the following on the command line:

   `$ blat genome /dev/null /dev/null -makeOoc=genome.11.ooc`

Please note that this pipeline has been written in Python2, and utilizes a GUI.

The first step of this pipeline is to obtain fasta sequences from the appropriate genome assembly by running the following on the command line:

`$ python BlatPipe.py coordinate_file 2bit_genome_file > outfile.seqs`

The file named outfile.seqsare all of the FASTA sequences pulled from the genome you originally specifyic with the input coordinates.

The fasta sequences can alternatively be obtained by utilizing bedtools&#39; getfasta utility:

`$ bedtools getfasta -fi genome_assembly -bed bedfile -fo output_fasta_file`

NOTE: This step is not necessary if you already have fasta sequences!

The second step of the pipeline is a genome comparison and alignment via BLAT and MUSCLE, respectively.

`$ python blatGenomes1.py outfile.seqs outfile.blatout`

This program uses the output from the previous step (outfile.seqs) as the input. Again, the first step is not required if you already have a FASTA file, which would take the place of the outfile.seqsfile. The outfile.blatoutis the user-defined name for the output file of this genome comparison step.


Figure 1 indicates the type of output you should expect from running blatGenomes1.py.

| ![image](https://user-images.githubusercontent.com/73801486/156855305-641ff926-19f4-4528-9543-929fd62cd123.png) |
|:--:|
| **Figure 1:** Example of the outfile.blatout file. Note the gaps that have been added as the result of alignment via MUSCLE.|

NOTE: When you run this program you will be faced with several prompts on the terminal window:

1. Input the name of your organism of interest
2. Type &quot;add&quot; if you want to add a 2bit format of a genome for comparison
   1. A popup window shows up of your files; you need to first choose the 2bit genome file and then it&#39;s corresponding ooc file

IMPORTANT: You will need at least 3 genomes for comparison or there will be downstream problems in this primer design pipeline. If you do not have 3 genomes for comparison, you can enter the same genome multiple times.

The third step is to calculate PCR product sizes and primer pairs based on the parameters hard-coded in the basetobasecompare.py program:

`$ python basetobasecompare.py > outfile.step3out`

NOTE: Once you start this program, a popup window will ask you to specify your &quot;comparison&quot; file. This is the outfile.blatoutfile from the previous step. Example output can be found in Figure 2.

| ![image](https://user-images.githubusercontent.com/73801486/156855354-a6f87a60-58c7-476a-a852-a9c5858b7459.png) |
|:--:|
| **Figure 2:** Example of the outfile.step3out. The red bracket indicates the expected product size when a PCR followed by gel electrophoresis is performed. Each entry is separated by an &quot;=&quot; symbol (black arrow). Any loci that fail to produce a primer pair based on the defined parameters will have the word, &quot;FAILED&quot; beneath the seq_ID (red box.). |





You can optionally use a mispriming library; this library states which sequences Primer3 should avoid when generating primers (e.g., repetitive sequences). Information and downloads regarding mispriming libraries can be found here: [https://bioinfo.ut.ee/primer3-0.4.0/primer3/input-help.htm#PRIMER_MISPRIMING_LIBRARY](https://bioinfo.ut.ee/primer3-0.4.0/primer3/input-help.htm#PRIMER_MISPRIMING_LIBRARY)

The fourth step is to reformat the output from the third step in order to facilitate ordering primers for PCR analysis.

`$ python primerOrderGenerator.py outfile.step3out outfile.laststepout primerPrefix`

NOTE: The &quot;primerPrefix&quot; is what is going to start the name of the output primer, e.g. if I type &quot;TestPrimer&quot; as my primerPrefix, the name of the output primer is going to be:

```
TestPrimer_0F\_#

TestPrimer_0R\_#

TestPrimer_1F\_#.....etc.
```

The purpose of this program is to generate a format appropriate for ordering primers. Example output can be found in Figure 3.

| ![image](https://user-images.githubusercontent.com/73801486/156855370-86b8a42f-7ddf-4368-90af-200191d1af59.png) |
|:--:|
| **Figure 3:** Example output of primerOrderGenerator.py. The primer number starts with &quot;0&quot;, and shows each primer pair name and sequence. |

This output file can be opened and saved in Excel.

Troubleshooting:

1. What to do if your outfile.blatout is empty:</br>
   You may have hidden characters in your file that the program doesn&#39;t like. For example, there should be a newline hidden character to denote new lines (&#39;\n&#39;) that is not seen when you open a file, but stored so that a file can have proper formatting. One hidden character in particular that this pipeline does not like is &#39;\r&#39;, or the &quot;raw&quot; character. To remove unwanted characters from your lines, use the following program: reveal_characters_and_replace.py
   </br>

2. Make sure that you have used the correct version of Python when calling the programs! The above assumes that Python2 is your default. However, you may need to specifically refer to Python2 on the command line.





Thomas Beckstrom;
Resident, Dept of Oral and Maxillofacial Surgery;
University of Washington;
1959 NE Pacific Street;
Health Sciences Building B-241;
Seattle, WA 98195-7134;
206.543.7496
