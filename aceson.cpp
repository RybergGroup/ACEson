#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <io_lib/Read.h>
//#include "contig_handler.h"
#include "alignment.h"
#include "sequence.h"

using namespace std;

void help () {
    cerr << "aceson is a program to convert ACE and/or AB1 files to a JSON format." << endl;
}

int main ( int argc, char* argv[]) {
    char output_format = 'j'; // j for jason a for ace f for fastq
    char input_format = 'C'; // C for aCe, B for aB1
    bool align = 0;
    char** filename = 0;
    unsigned int n_files = 0;
    char q_cut = 0x00;
    istream* f = &cin;
    ifstream in_file;
    ostream* out = &cout;
    ofstream out_file;
    //struct contig_node* contigs;
    for (int i=1; i < argc; ++i) {
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
		    cerr << "-f / --file require a file name as next argument." << endl;
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
                        cerr << argv[i] << " is not a recognized input format" << endl;
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
			cerr << argv[i] <<  " is not a recognized output format" << endl;
			exit(1);
		    }
		}
	    }
	    else if (!strcmp(argv[i],"-o") || !strcmp(argv[i],"--output_file")) {
		if (i+1 < argc && argv[++i][0] != '-') {
		    out_file.open(argv[i]);
		    out = &out_file;
		}
		else {
		    cerr << "-o / --output require a file name as next argument." << endl;
		    exit(1);
		}
	    }
	    else if (!strcmp(argv[i],"-q") || !strcmp(argv[i],"--qual_cut_off")) {
                if (i+1 < argc && argv[++i][0] != '-') {
                    q_cut = (char)atoi(argv[i]);
                }
                else {
                    cerr << "-q / --qual_cut_off require a number as next argument." << endl;
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
		cerr <<  "Do not recognize argument: " <<  argv[i] << '.' << endl;
		exit(1);
	    }
	}
       	else if (filename == 0 && n_files==0) {
	    filename=(char**)malloc(sizeof(char*));
	    filename[0] = argv[i];
	    ++n_files;
	}
	else {
	    cerr <<  "Do not recognize argument: " << argv[i] << '.' << endl;
	    exit(1);
	}
    }
    if (align) {
	cerr << "Aligning..." << endl;
	Read** read_data = new Read*[n_files];
	sequence* sequences = new sequence[n_files];
	vector<alignment> alignments;
	// Read sequences
	for (unsigned int i=0; i < n_files; ++i) {
	    cerr << "Reading first file: " << filename[i] << endl;
	    read_data[i] = read_reading(filename[i], 0);
	    if (q_cut) {
		cerr << "Quality clipping based on score " << (int)q_cut << endl;
		for (int j=0; j < read_data[i]->NBases; ++j) {
		    if (!read_data[i]->leftCutoff) {
			if ((read_data[i]->base[j] == 'A' || read_data[i]->base[j] == 'a') && read_data[i]->prob_A[j] > q_cut)
			    read_data[i]->leftCutoff = j;
			else if ((read_data[i]->base[j] == 'C' || read_data[i]->base[j] == 'C') && read_data[i]->prob_C[j] > q_cut)
			    read_data[i]->leftCutoff = j;
			else if ((read_data[i]->base[j] == 'G' || read_data[i]->base[j] == 'g') && read_data[i]->prob_G[j] > q_cut)
			    read_data[i]->leftCutoff = j;
			else if ((read_data[i]->base[j] == 'T' || read_data[i]->base[j] == 't') && read_data[i]->prob_T[j] > q_cut)
			    read_data[i]->leftCutoff = j;
		    }
		    if (!read_data[i]->rightCutoff) {
                        if ((read_data[i]->base[read_data[i]->NBases-1-j] == 'A' || read_data[i]->base[read_data[i]->NBases-1-j] == 'a') && read_data[i]->prob_A[read_data[i]->NBases-1-j] > q_cut)
                            read_data[i]->rightCutoff = read_data[i]->NBases-1-j;
                        else if ((read_data[i]->base[read_data[i]->NBases-1-j] == 'C' || read_data[i]->base[read_data[i]->NBases-1-j] == 'C') && read_data[i]->prob_C[read_data[i]->NBases-1-j] > q_cut)
                            read_data[i]->rightCutoff = read_data[i]->NBases-1-j;
                        else if ((read_data[i]->base[read_data[i]->NBases-1-j] == 'G' || read_data[i]->base[read_data[i]->NBases-1-j] == 'g') && read_data[i]->prob_G[read_data[i]->NBases-1-j] > q_cut)
                            read_data[i]->rightCutoff = read_data[i]->NBases-1-j;
                        else if ((read_data[i]->base[read_data[i]->NBases-1-j] == 'T' || read_data[i]->base[read_data[i]->NBases-1-j] == 't') && read_data[i]->prob_T[read_data[i]->NBases-1-j] > q_cut)
                            read_data[i]->rightCutoff = read_data[i]->NBases-1-j;
                    }
		}
		cerr << filename[i] << " n bases: " << read_data[i]->NBases << ", qual start: " << read_data[i]->leftCutoff << ", qual end: " << read_data[i]->rightCutoff << endl;
	    }
	    sequences[i].add_read(read_data[i]);
	    cerr << "Prepare for alignment" << endl;
	    alignments.emplace_back(&sequences[i]);
	    cerr << "Prepare for alignment" << endl;
	    cerr << "Length: " << sequences[i].length() << " (" << alignments.back().alignment_length() << ")" << " [" << alignments.front().alignment_length() << "]" << endl;
	}
	// Align first two sequences
	cerr << "Aligning first pair..." << endl;
	cerr << "Length " << alignments[0].alignment_length() << " " << alignments[1].alignment_length() << endl;
	alignment ali(2);
       	ali.align_pair(alignments[0], alignments[1]);
	if (output_format == 'F')
	    ali.print_alignment_fasta(*out);
	else if (output_format == 'j') {
	    *out << "{\"contigs\":[";
	    //bool first = true;
	    //fputs("Check 1\n", stderr);
	    //while (contigs != NULL) {
		//fputs("Check 2\n", stderr);
	    //  if (first)
	    //     first = false;
	    //  else
	    //      fputc(',',out);
	    ali.set_name("contig1");
	    ali.print_alignment_json(*out);
	           //if (contigs->next != NULL)
        //    fputc(',',out);
        //contigs = contigs->next;
    //}
	    *out << "]}" << endl;
	}
       	//else if (output_format == 'a')
	// Clean up
	delete[] sequences;
	for (unsigned int i=0; i < n_files; ++i) {
	    read_deallocate(read_data[i]);
	}
	delete read_data;
    }
    else {
/*	fputs("Not aligning\n",stderr);
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
	}*/
    }
    if (filename)
	free(filename);
}
