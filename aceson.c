#include <stdio.h>
#include <io_lib/Read.h>
#include "contig_handler.h"
#include "alignment.h"

void help () {
    puts("aceson is a program to convert ACE and/or AB1 files to a JSON format.");
}

int main ( int argc, char* argv[]) {
    char output_format = 'j'; // j for jason a for ace f for fastq
    char input_format = 'C'; // C for aCe, B for aB1
    _Bool align = 0;
    char** filename = 0;
    unsigned int n_files = 0;
    FILE* f = stdin;
    FILE* out = stdout;
    //struct contig_node* contigs;
    for (unsigned int i=1; i < argc; ++i) {
	if (argv[i][0] == '-') {
	    if (!strcmp(argv[i],"-f") || !strcmp(argv[i],"--file")) {
		unsigned int first = i+1;
		unsigned int last = first;
		while (i+1 < argc && argv[i+1][0] != '-') {
		    ++last;
		    ++i;
		}
		if (last != first) {
		    if (filename == 0)			
			filename = (char**)malloc((last-first)*sizeof(char*));
		    else
			filename = (char**)realloc(filename,(last-first+n_files)*sizeof(char*));
		    for (; first < last; ++first) {
			filename[n_files] = argv[first];
			++n_files;
		    }
		}
		else {
		    fputs("-f / --file require a file name as next argument.\n", stderr);
		    exit(1);
		}
	    }
	    else if (!strcmp(argv[i],"-I") || !strcmp(argv[i],"--input_format")) {
		if (i+1 < argc && argv[++i][0] != '-') {
                    if (!strcmp(argv[i],"ace"))
                        input_format = 'C';
                    else if (!strcmp(argv[i],"ab1")){
                        input_format = 'B';
		    }
                    else {
                        fprintf(stderr, "%s is not a recognized input format\n", argv[i]);
                        exit(1);
                    }
                }
	    }
	    else if (!strcmp(argv[i],"-O") || !strcmp(argv[i],"--output_format")) {
                if (i+1 < argc && argv[++i][0] != '-') {
		    if (!strcmp(argv[i],"aceson"))
			output_format = 'j';
		    else if (!strcmp(argv[i],"ace"))
			output_format = 'a';
		    else if (!strcmp(argv[i],"fastq"))
                        output_format = 'f';
		    else if (!strcmp(argv[i],"fasta"))
                        output_format = 'F';
		    else {
			fprintf(stderr, "%s is not a recognized output format\n", argv[i]);
			exit(1);
		    }
		}
	    }
	    else if (!strcmp(argv[i],"-o") || !strcmp(argv[i],"--output_file")) {
		if (i+1 < argc && argv[++i][0] != '-') {
		    out = fopen(argv[i],"w");
		}
		else {
		    fputs("-o / --output require a file name as next argument.\n", stderr);
		    exit(1);
		}
	    }
	    else if (!strcmp(argv[i],"-A") || !strcmp(argv[i],"--align")) {
		align = 1;
            }
	    else if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) {
		help();
                exit(0);
            }
	    else {
		fprintf(stderr, "Do not recognize argument: %s.\n", argv[i]);
		exit(1);
	    }
	}
       	else if (filename == 0 && n_files==0) {
	    filename=(char**)malloc(sizeof(char*));
	    filename[0] = argv[i];
	    ++n_files;
	}
	else {
	    fprintf(stderr, "Do not recognize argument: %s.\n", argv[i]);
	    exit(1);
	}
    }
    if (align) {
	fputs("Aligning...\n",stderr);
	Read** read_data = (Read**)calloc(n_files,sizeof(Read*));
	struct alignment** sequences = (struct alignment**)calloc(n_files,sizeof(struct alignment*));
	// Read sequences
	for (unsigned int i=0; i < n_files; ++i) {
	    fprintf(stderr, "Reading file: %s. ...\n", filename[i]);
	    read_data[i] = read_reading(filename[i], 0);
	    //fputs(read_data[i]->trace_name,stderr);
	    struct multiple_sequences* temp = multiple_sequences_alocate(1);
	    add_read(temp, read_data[i], 0);
	    sequences[i] = alignment_alocate(temp);
	    //fprint_fasta_alignment(sequences[i], stderr);
	}
	// Align first two sequences
	struct alignment* ali = align_pair(sequences[0], sequences[1]);
	//fprint_fasta_alignment(ali, out);
	//fprintf(stderr, "Made it.\n");
	if (output_format == 'F')
	    fprint_fasta_alignment(ali, out);
	else if (output_format == 'j' || output_format == 'a') {
	    struct contig_node* contigs = (struct contig_node*)calloc(sizeof(struct contig_node),1);
	    contigs->data = alignment_to_contig("contig1",ali,read_data);
	    if (output_format == 'j')
		printJSON(out, contigs);
	    else if (output_format == 'a')
		printACE(out, contigs);
	    else
		fputs("Unrecognised output format\n",stderr);
	    fputs("Going to delete contigs\n", stderr);
	    delete_contigs_x(contigs, 0);
	    fputs("Deleted contig\n", stderr);
	}
	// Clean up
	if (ali) {
	    alignment_ms_dealocate(ali);
	    fputs("Deleted ali\n", stderr);
	}
	for (unsigned int i=0; i < n_files; ++i) {
	    read_deallocate(read_data[i]);
	    fputs("Check C\n", stderr);
	    alignment_ms_dealocate(sequences[i]);
	    fprintf(stderr, "Deleted file and seq %u\n", i);
	}
	if (sequences)
	    free(sequences);
	if (read_data)
	    free(read_data);
	fputs("Freed seqs and reads\n", stderr);	
	
    }
    else {
	fputs("Not aligning\n",stderr);
	if (input_format == 'C') {
	    if (filename != 0)
		f = fopen(filename[0],"r");
	    if (!f) {
		fputs("Could not open input file\n", stderr);
		exit(1);
	    }
	    if (!out) {
		fputs("Could not open output file\n", stderr);
		exit(1);
	    }
	    struct contig_node* contigs = parsACE(f, NULL);
	    add_traces(contigs);
	    if (output_format == 'a')
		printACE(out, contigs);
	    else if (output_format == 'j')
		printJSON(out, contigs);
	    delete_contigs(contigs);
	    fclose(f);
	    fclose(out);
	}
	else if (input_format == 'B') {
	    if (output_format == 'j') {
		fputs("Printing iaceson\n",stderr);
		fputs("{\"reads\":[", out);
	    }
	    for (unsigned int i=0; i < n_files; ++i) {
		Read* read_data = read_reading(filename[i], 0);
		if (output_format == 'f') {
		    fputs("Printing fastq\n",stderr);
		    fprintf(out,"@%s\n", read_data->trace_name);
		    for (unsigned int i = 0; i < read_data->NBases; ++i) {
			fputc(read_data->base[i], out);
		    }
		    fputs("\n+\n", out);
		    for (unsigned int i = 0; i < read_data->NBases; ++i) {
			if (read_data->base[i] == 'A')
			    fputc(read_data->prob_A[i]+0x21, out);
			else if (read_data->base[i] == 'C')
			    fputc(read_data->prob_C[i]+0x21, out);
			else if (read_data->base[i] == 'G')
			    fputc(read_data->prob_G[i]+0x21, out);
			else if (read_data->base[i] == 'T')
			    fputc(read_data->prob_T[i]+0x21, out);
		    }
		    fputc('\n', out);
		}
		else if (output_format == 'j') {
		    if (i)
			fputc(',',out);
		    fputs("{\"chrom\":", out);
		    print_read_dataJSON ( out, read_data);
		    fputc('}', out);
		}
		read_deallocate(read_data);
	    }
	    if (output_format == 'j') {
    		fputs("]}\n", out);
	    }
	}
    }
    if (filename)
	free(filename);
}
