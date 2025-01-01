/* Program for demultiplexing dual-indexed Illumina reads */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define HEADER_LEN 300
#define SEQUENCE_LEN 500
#define BARCODE_LEN 10 
#define PLUS_LEN 200

int count_lines(char file_name[]);
char *rev_complement(char *seq);

int main(int argc, char *argv[])
{
	/* argv [1] = R1
	   argv [2] = R2
	   argv [3] = R3
	   argv [4] = R4
	   argv [5] = list of barcodes */ 

	/* use count_lines divided by 4 to get entry count. */
	int num_reads = count_lines(argv[1]) / 4;

	/* read in list of barcodes */ 
	FILE *barcodelist;
	int num_barcodes = count_lines(argv[5]) - 1; // there is a header line in tsv of barcodes
	barcodelist = fopen(argv[5], "r");

	char barcodearray[num_barcodes+1][100];
	char rc_barcodearray[num_barcodes+1][100];
	char *temp;
	int i, j, k;

	fgets(barcodearray[0], sizeof(barcodearray[i]), barcodelist); // first line is overwritten  
	for (i = 0; i < num_barcodes; i++) {
		fgets(barcodearray[i], sizeof(barcodearray[i]), barcodelist);
		barcodearray[i][strcspn(barcodearray[i], "\n")] = 0;

		temp = strtok(barcodearray[i], " ");
		for (j = 0; j < 4; j++) {
		temp = strtok(NULL, " ");
		if (j == 3) strcpy(barcodearray[i], temp);
		}
	}
	fclose(barcodelist);

	// rev_complement will overwrite, copy is needed
	memcpy(rc_barcodearray, barcodearray, sizeof(barcodearray));

	/* print out the barcodes taken from the barcode list file */
	printf("Barcodes:\n");
	for (i = 0; i < num_barcodes; i++) { 
		printf("%d_fwd: %s\n", i+1, barcodearray[i]);
		printf("%d_rev: %s\n", i+1, rev_complement(rc_barcodearray[i]));
	}

	// start reading through files for demultiplexing
	FILE *read1, *read2, *read3, *read4;

	read1 = fopen(argv[1], "r");
	read2 = fopen(argv[2], "r");
	read3 = fopen(argv[3], "r");
	read4 = fopen(argv[4], "r");

	/* output files */
	FILE *unkr1, *unkr2;
	unkr1 = fopen("output/unk_R1.fq", "w");
	unkr2 = fopen("output/unk_R2.fq", "w");
	FILE *hoppedr1, *hoppedr2;
	hoppedr1 = fopen("output/hopped_R1.fq", "w");
	hoppedr2 = fopen("output/hopped_R2.fq", "w");

	/* create new files for each barcode */ 
  char filenameR1[BARCODE_LEN+13+1], filenameR2[BARCODE_LEN+13+1];
  FILE *foutR1, *foutR2;
	for (i = 0; i < num_barcodes; i++){
    snprintf(filenameR1, sizeof(filenameR1), "output/%s_R1.fq", barcodearray[i]);
    snprintf(filenameR2, sizeof(filenameR2), "output/%s_R2.fq", barcodearray[i]);
    foutR1 = fopen(filenameR1, "w");
    foutR2 = fopen(filenameR2, "w");
    fclose(foutR1);
    fclose(foutR2);
	}

	char r1name[HEADER_LEN+1], r1seq[SEQUENCE_LEN+1];
	char r1plus[PLUS_LEN+1], r1qual[SEQUENCE_LEN+1]; 
	char r2name[HEADER_LEN+1], r2seq[SEQUENCE_LEN+1];
	char r2plus[PLUS_LEN+1], r2qual[SEQUENCE_LEN+1]; 
	char r3name[HEADER_LEN+1], r3seq[SEQUENCE_LEN+1];
	char r3plus[PLUS_LEN+1], r3qual[SEQUENCE_LEN+1]; 
	char r4name[HEADER_LEN+1], r4seq[SEQUENCE_LEN+1];
	char r4plus[PLUS_LEN+1], r4qual[SEQUENCE_LEN+1]; 

	int hopped, unk, matched;

	for (i=0; i < num_reads; i++) {
		hopped = 0;
		unk = 0;
		matched = 0;

		fgets(r1name, sizeof(r1name), read1);
		fgets(r1seq, sizeof(r1seq), read1);
		fgets(r1plus, sizeof(r1plus), read1);
		fgets(r1qual, sizeof(r1qual), read1);

		fgets(r2name, sizeof(r2name), read2);
		fgets(r2seq, sizeof(r2seq), read2);
		fgets(r2plus, sizeof(r2plus), read2);
		fgets(r2qual, sizeof(r2qual), read2);

		fgets(r3name, sizeof(r3name), read3);
		fgets(r3seq, sizeof(r3seq), read3);
		fgets(r3plus, sizeof(r3plus), read3);
		fgets(r3qual, sizeof(r3qual), read3);
		
		fgets(r4name, sizeof(r4name), read4);
		fgets(r4seq, sizeof(r4seq), read4);
		fgets(r4plus, sizeof(r4plus), read4);
		fgets(r4qual, sizeof(r4qual), read4);
		
		/* WIP: add quality filtering for the barcodes */
		
		for (j = 0; j < num_barcodes; j++) { // check each barcode in list against the r2 barcode
			if(strncmp(barcodearray[j], r2seq, 8) == 0) { // if r2seq is a valid barcode
				if (strncmp(rc_barcodearray[j], r3seq, 8) == 0) { // if r3seq is matching (revcomp) 
					matched++; // matched == 1 
				}
				break;
			}
			else if (j == num_barcodes - 1) unk++;
		}

		if (matched == 0 && unk == 0) { // check the second barcode (r3) if necessary 
			for (k = 0; k < num_barcodes; k++) {
				if (strncmp(rc_barcodearray[k], r3seq, 8) == 0) {
					break;
				}
				else if (k == num_barcodes - 1) unk++;
			}
			if (j == k) matched++;
			else if (unk == 0 && j != k) hopped++;
		}

		if (matched + unk + hopped != 1) {
			printf("Error:\n");
			printf("matched: %d\n", matched);
			printf("unk: %d\n", unk);
			printf("hopped: %d\n", hopped);
			break;
		}

		if (matched != 0) {
			// this is a variable filename 
    	snprintf(filenameR1, sizeof(filenameR1), "output/%s_R1.fq", barcodearray[j]);
    	snprintf(filenameR2, sizeof(filenameR2), "output/%s_R2.fq", barcodearray[j]);
    	foutR1 = fopen(filenameR1, "a");
    	foutR2 = fopen(filenameR2, "a");

			/* write to these files here */
			fprintf(foutR1, r1name);
			fprintf(foutR1, r1seq);
			fprintf(foutR1, r1plus);
			fprintf(foutR1, r1qual);
			fprintf(foutR2, r4name);
			fprintf(foutR2, r4seq);
			fprintf(foutR2, r4plus);
			fprintf(foutR2, r4qual);

    	fclose(foutR1);
    	fclose(foutR2);
		}

		else if (unk != 0) { // write to files unk_R1.fq and unk_R2.fq	
			fprintf(unkr1, r1name);
			fprintf(unkr1, r1seq);
			fprintf(unkr1, r1plus);
			fprintf(unkr1, r1qual);
			fprintf(unkr2, r4name);
			fprintf(unkr2, r4seq);
			fprintf(unkr2, r4plus);
			fprintf(unkr2, r4qual);
		}

		else if (hopped != 0) { // write to files hopped_R1.fq and hopped_R2.fq
			fprintf(hoppedr1, r1name);
			fprintf(hoppedr1, r1seq);
			fprintf(hoppedr1, r1plus);
			fprintf(hoppedr1, r1qual);
			fprintf(hoppedr2, r4name);
			fprintf(hoppedr2, r4seq);
			fprintf(hoppedr2, r4plus);
			fprintf(hoppedr2, r4qual);
		}

	}

	fclose(read1);
	fclose(read2);
	fclose(read3);
	fclose(read4);
	fclose(unkr1);
	fclose(unkr2);
	fclose(hoppedr1);
	fclose(hoppedr2);

	return 0;
}

int count_lines(char file_name[])
{
	FILE *file;
	int ch, line_num = 0;

	file = fopen(file_name, "r");
	if (file == NULL) {
		printf("Can't open %s\n", file_name);
		exit(EXIT_FAILURE);
	}

	for (ch = getc(file); ch != EOF; ch = getc(file))
		if (ch == '\n') line_num++;

	fclose(file);
	//printf("%d\n", line_num);
	return line_num;
}

char *rev_complement(char *sequence)
{
		char complement_base(char base);

    if (!sequence || ! *sequence)
        return sequence;

    int i = strlen(sequence) - 1, j = 0;

    char temp;
    while (i > j)
    {
        temp = complement_base(sequence[i]);
        sequence[i] = complement_base(sequence[j]);
        sequence[j] = temp;
        i--;
        j++;
    }
    return sequence;
}

char complement_base(char base)
{
	switch(base)
	{
		case 'A':
			return 'T';
		case 'a':
			return 't';
		case 'T': case 'U':
			return 'A';
		case 't': case 'u':
			return 'a';
		case 'C':
			return 'G';
		case 'c':
			return 'g';
		case 'G':
			return 'C';
		case 'g':
			return 'c';
		case 'N':
			return 'N';
		case 'n':
			return 'n';
		case ' ':
			return ' ';
		default:
			return '\0';
	}
}
