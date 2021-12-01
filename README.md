# CircAST2
CircAST2 is an effective tool for the assembly and quantification of linear, circular transcripts based on total transcriptome sequencing data.

CircAST2 based on python3.8，there are 2 versions，we recommend running CircAST2.py directly.

CircAST2 need a pair-ended sequencing SAM file by Hisat2,circRNA junction TXT file by CIRI2,Reference genome GTF file and FASTA file,also StringTie needs to be built in and provided with a run path,reads length is the Length of pair-end sequencing reads,threshold means the minimum read supportting CircRNA BSj(we recommend 1).


The running parameters are as follows：

python CircAST2.py -L (reads length) -T(threshold) -G (gtf file) -I (input SAM file) -J (circRNA junction list file) -F(genome fasta file path) -S (StringTieDir)


Circast2-sim is a whole-transcriptome simulation sequencing data generation tool,Users are required to provide three input files, which are respectively the annotation file in GTF format, the reference genome sequence file in FASTA format and the circular RNA list information file in TXT format.
