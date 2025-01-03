# c_demultiplex
Program for demultiplexing barcoded illumina data -written in C

## Use: 

`gcc main.c` creates executable `a.out`. 

`./a.out <filepath_R1.fq> <filepath_R2.fq> <filepath_R3.fq> <filepath_R4.fq> <filepath_indices.txt> <(barcode)quality-cutoff>`

User-specified quality cutoff for the barcodes should be an integer value.

Must use command `mkdir -p output` prior to running this program, 
