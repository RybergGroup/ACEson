#include <stdio.h>
#include <io_lib/Read.h>
#include "contig_handler.h"

void help () {
    puts("aceson is a program to convert ACE and/or AB1 files to a JSON format.");
}

int main ( int argc, char* argv[]) {
    char output_format = 'j';
    char input_format = 'C'; // C for aCe, B for aB1
    char* filename;
    FILE* f = stdin;
    FILE* out = stdout;
    //struct contig_node* contigs;
    for (unsigned int i=1; i < argc; ++i) {
	if (argv[i][0] == '-') {
	    if (!strcmp(argv[i],"-f") || !strcmp(argv[i],"--file")) {
		if (i+1 < argc && argv[++i][0] != '-') {
		    filename = argv[i];
		    f = fopen(argv[i],"r");
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
	}
       	else if (f == stdin) {
	    f = fopen(argv[i],"r");
	}
	else {
	    fprintf(stderr, "Do not recognize argument: %s.\n", argv[i]);
	}
    }
    if (!f) {
	fputs("Could not open input file\n", stderr);
	exit(1);
    }
    if (!out) {
	fputs("Could not open output file\n", stderr);
	exit(1);
    }
    /*
    fputs(filename,stderr);
    fputc('\n',stderr);
    fputc(input_format,stderr);
    fputc('\n',stderr);
    */
    if (input_format == 'C') {
	struct contig_node* contigs = parsACE(f, NULL);
	add_traces(contigs);
	//fprintf(stderr, "N contigs %u | n reads %u\n", n_contigs(contigs), n_reads(contigs));
	if (output_format == 'a')
	    printACE(out, contigs);
	else if (output_format == 'j')
	    printJSON(out, contigs);
	delete_contigs(contigs);
    }
    else if (input_format == 'B') {
	Read* read_data = fread_reading(f, filename, 0);
	if (output_format == 'f') {
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
	    print_read_dataJSON ( out, read_data);
	    fputc('\n', out);
	}
	read_deallocate(read_data);
    }
    fclose(f);
    fclose(out);
}
