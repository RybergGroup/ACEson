#include "contig_handler.h"

///////////////////////////
/// assisting functions ///
///////////////////////////
_Bool isWhitespace ( char c ) {
    if (c == ' ' || c == '\t' || c == '\n' || c == '\r') return true;
    return false;
}
_Bool isNewlineChar ( char c ) {
    if (c == '\n' || c == '\r') return true;
    return false;
}

void delete_contigs ( struct contig_node* node ) {
    if (node->next !=0)
        delete_contigs(node->next);
    for (unsigned int i=0; i < node->data.n_reads; ++i) {
        if (node->data.reads[i].seq != NULL)
            free(node->data.reads[i].seq);
        if (node->data.reads[i].name != NULL)
            free(node->data.reads[i].name);
        if (node->data.reads[i].description != NULL)
            free(node->data.reads[i].description);
        if (node->data.reads[i].read_data!= NULL)
            read_deallocate(node->data.reads[i].read_data);
    }
    if (node->data.reads != NULL)
        free(node->data.reads);
    if (node->data.contig_seq != NULL)
        free(node->data.contig_seq);
    if (node->data.contig_qual != NULL)
        free(node->data.contig_qual);
    if (node->data.segment != NULL)
	free(node->data.segment);
    if (node != NULL)
        free(node);
}

struct contig_node* add_contig(struct contig_node* parent) {
    if (parent == NULL) {
        parent = (struct contig_node*)calloc(1,sizeof(struct contig_node));
        return parent;
    }
    else if (parent->next != NULL)
        return add_contig(parent->next);
    else {
        parent->next = (struct contig_node*)calloc(1,sizeof(struct contig_node));
        return parent->next;
    }
}

unsigned int n_contigs (struct contig_node* start) {
    unsigned int n = 0;
    while (start != NULL) {
        ++n;
        start = start->next;
    }
    return n;
}

unsigned int n_reads ( struct contig_node* start ) {
    unsigned int n=0;
    while (start != NULL) {
        n += start->data.n_reads;
        start = start->next;
    }
    return n;
}

////////////////////////

struct contig_node* parsACE ( FILE* in, struct contig_node* contigs ) {
    char word[WORDSIZE] = "\0";
    char mode = '0';
    struct contig_node* first_contig;
    struct contig_node* present_contig;
    _Bool newline = true;
    unsigned int contigNo = 0;
    unsigned int readNo = 0;
    int baseNo = 0;
    char c;
    unsigned int expected_contigs;
    unsigned int expected_reads;
    while ((c = fgetc(in)) != EOF) {
	if (isWhitespace(c)) {
	    //fprintf(stderr, "Word: %s | mode: %c\n", word, mode);
	    if (word[0] != '\0') {
		/// Identify ACE entry ///
		if (newline && !strcmp(word,"AS")) {
		    mode = 'B';
		    fprintf(stderr, "At start\n");
		}
		else if (newline && !strcmp(word,"CO")) {
		    mode = 'C';
		    baseNo= 0;
		    present_contig = add_contig(contigs);
		    if (first_contig == NULL)
			first_contig = present_contig;
		    if (contigs == NULL)
			contigs = first_contig;
		    --expected_contigs;
		    readNo = 0;
		    fprintf(stderr, "Found start of contig\n");
		}
		else if (newline && contigNo < n_contigs(first_contig) &&  present_contig != NULL && !strcmp(word,"BQ")) {
		    mode = 'Q';
		    baseNo= 0;
		    present_contig->data.contig_qual = (unsigned int*)calloc(present_contig->data.contig_length - present_contig->data.n_gaps,sizeof(unsigned int));
		    //fprintf(stderr, "Found q scores\n");
		}
		else if (newline && contigNo < n_contigs(first_contig) &&  present_contig != NULL &&
				    !strcmp(word,"AF")) { mode = 'A'; fprintf(stderr, "Found read %d\n", readNo);}
		else if (newline && !strcmp(word,"RD")) { mode = 'R'; }
		else if (newline && !strcmp(word,"QA")) { mode = 'S'; }
		else if (newline && !strcmp(word,"DS")) { mode = 'D'; fprintf(stderr, "Found read description for read %d\n", readNo);}
		else if (newline && !strcmp(word,"BS")) { mode = 'E'; }
		/// read data ///
		else if (mode == 'B') {
		    expected_contigs = atoi(word);
		    mode = 'b';
		}
		else if (mode == 'b') {
		    expected_reads = atoi(word);
		    mode = '0';
		}
		else if (mode == 'C') {
		    if (present_contig->data.name == NULL ) {
			present_contig->data.name = (char*)malloc((strlen(word)+1)*sizeof(char));
			strcpy(present_contig->data.name,word);
			fprintf(stderr, "Name %s\n", present_contig->data.name );
		    }
		    else if (present_contig->data.contig_seq == NULL) {
			present_contig->data.contig_length = atoi(word);
			present_contig->data.contig_seq = (char*)calloc((present_contig->data.contig_length+1),sizeof(char));
			present_contig->data.contig_seq[0] = '\0';
		    }
		    else if (present_contig->data.reads == NULL) {
			present_contig->data.n_reads = atoi(word);
			present_contig->data.reads = (contig_read*)calloc(present_contig->data.n_reads,sizeof(contig_read));
		    }
		    else if (!strcmp(word,"C"))
			present_contig->data.revcomp = 1;
		    if (isNewlineChar (c))
			mode = 'c';
		}
		else if (mode == 'c') {
		    unsigned int length = strlen(word);
		    unsigned int i = 0;
		    for (; i < length && baseNo < present_contig->data.contig_length; ++i) {
			present_contig->data.contig_seq[baseNo] = word[i];
			if (word[i] == '*')
			    ++present_contig->data.n_gaps;
			++baseNo;
		    }
		    present_contig->data.contig_seq[baseNo] = '\0';
		    if (i < length)
			fprintf(stderr, "Sequence length overflow! %d characters read as part of the contig sequence of %s did not fit in the given sequence length (%d)\n", length - i, present_contig->data.name, present_contig->data.contig_length);
		}
		else if (mode == 'Q') {
		    if (baseNo < present_contig->data.contig_length - present_contig->data.n_gaps) present_contig->data.contig_qual[baseNo] = atoi(word);
		    else
			fprintf(stderr, "Quality length overflow! %d words read as quality scores did not fit in the given length of the sequence\n", baseNo - (present_contig->data.contig_length-present_contig->data.n_gaps));
		    ++baseNo;
		}
		else if (mode == 'A') {
		    if (readNo >= present_contig->data.n_reads) {
			continue;
		    }
		    else if (present_contig->data.reads[readNo].name == NULL) {
			present_contig->data.reads[readNo].name = (char*)malloc((strlen(word)+1)*sizeof(char));
			strcpy(present_contig->data.reads[readNo].name,word);
		    }
		    else if (!strcmp(word,"C")) {
			present_contig->data.reads[readNo].revcomp = true;
		    }
		    else if (!strcmp(word,"U")) {
			present_contig->data.reads[readNo].revcomp = false;
		    }
		    else {
			present_contig->data.reads[readNo].start = atoi(word);
		    }
		    if (isNewlineChar (c)) {mode = '0'; ++readNo;};
		}
		else if (mode == 'R') {
		    for (unsigned int i =0; i < present_contig->data.n_reads; ++i) {
			if (present_contig->data.reads[i].name == NULL) {
			    present_contig->data.reads[readNo].name = (char*)malloc((strlen(word)+1)*sizeof(char));
			    strcpy(present_contig->data.reads[readNo].name,word);
			    readNo=i;
			    mode = 'l';
			    break;
			}
        		else if (!strcmp(word,present_contig->data.reads[i].name)) {
			    if ( present_contig->data.reads[i].seq != NULL && strlen(present_contig->data.reads[i].seq) > 0 ) {
				fprintf(stderr, "Sequence already read for %s, may be duplicate in read names.\n", word);
			    }
			    else {
				readNo=i;
				mode = 'l';
				break;
			    }
			}
		    }
		    if (readNo >= present_contig->data.n_reads) {
			fprintf(stderr,"Could not find read %s, and the given number of reads have been read.\n", word);
		    }
		}
		else if (mode == 'l') {
		    if (readNo < present_contig->data.n_reads && present_contig->data.reads[readNo].seq == NULL) {
			present_contig->data.reads[readNo].seq_length = atoi(word);
			present_contig->data.reads[readNo].qual_begin = 0;
			present_contig->data.reads[readNo].qual_end = present_contig->data.reads[readNo].seq_length +1;
			present_contig->data.reads[readNo].align_begin = 0;
			present_contig->data.reads[readNo].align_end = present_contig->data.reads[readNo].seq_length +1;
			present_contig->data.reads[readNo].seq = (char*)calloc(present_contig->data.reads[readNo].seq_length+1,sizeof(char));
		    }
		    if (isNewlineChar (c)) { mode = 'r'; }
		}
		else if (mode == 'r') {
		    if (readNo < present_contig->data.n_reads) {
			unsigned int length_before = strlen(present_contig->data.reads[readNo].seq);
			strncat(present_contig->data.reads[readNo].seq,word,present_contig->data.reads[readNo].seq_length - length_before);
			if (strlen(word) > present_contig->data.reads[readNo].seq_length - length_before)
			    fprintf(stderr, "Sequence length overflow! %ld characters read as part of the read sequence of %s did not fit in the given sequence length (%d)\n", strlen(word) - (present_contig->data.reads[readNo].seq_length - length_before), present_contig->data.reads[readNo].name, present_contig->data.reads[readNo].seq_length);

		    }
		    else {
			fprintf(stderr,"Read reading error. Number of read sequences out of bound!\n");
		    }
		}
		else if (mode == 'S') {
		    if (readNo < present_contig->data.n_reads) {
			if (present_contig->data.reads[readNo].qual_begin == 0) present_contig->data.reads[readNo].qual_begin = atoi(word);
			else if (present_contig->data.reads[readNo].qual_end > present_contig->data.reads[readNo].seq_length)  present_contig->data.reads[readNo].qual_end = atoi(word);
			else if (present_contig->data.reads[readNo].align_begin == 0) present_contig->data.reads[readNo].align_begin = atoi(word);
			else if (present_contig->data.reads[readNo].align_end < 0) present_contig->data.reads[readNo].align_end = atoi(word);
		    }
		}
		else if (mode == 'D') {
		    if (present_contig->data.reads[readNo].description == NULL) {
			present_contig->data.reads[readNo].description = (char*)malloc((strlen(word)+1)*sizeof(char));
			strcpy(present_contig->data.reads[readNo].description,word);
		    }
		    else {
			char * temp = present_contig->data.reads[readNo].description;
			present_contig->data.reads[readNo].description = (char*)malloc((strlen(word)+strlen(temp)+2)*sizeof(char));
			strcpy(present_contig->data.reads[readNo].description,temp);
			strcat(present_contig->data.reads[readNo].description," ");
			strcat(present_contig->data.reads[readNo].description,word);
			free(temp);
		    }
		    //fprintf(stderr,"Description is: %s\n", present_contig->data.reads[readNo].description);
		    if (isNewlineChar (c))
	    		mode = '0';
		    
		}
		else if (mode == 'E' || mode == 'e') {
		    if (present_contig->data.segment == NULL) {
			present_contig->data.segment = (char*)malloc((strlen(word)+4)*sizeof(char));
			strcpy(present_contig->data.segment,"BS ");
			strcat(present_contig->data.segment,word);
		    }
		    else {
			char * temp = present_contig->data.segment;
			if (mode == 'E') {
			    present_contig->data.segment = (char*)malloc((strlen(word)+strlen(temp)+5)*sizeof(char));
			    strcpy(present_contig->data.segment,temp);
			    strcat(present_contig->data.segment,"\nBS ");
			}
			else {
			    present_contig->data.segment = (char*)malloc((strlen(word)+strlen(temp)+2)*sizeof(char));
			    strcpy(present_contig->data.segment,temp);
			    strcat(present_contig->data.segment," ");
			}
			strcat(present_contig->data.segment,word);
			free(temp);
		    }
		    if (isNewlineChar (c))
			mode = '0';
		    else
			mode = 'e';
		}
		else
		    mode = '0';
		//System.err.println(word);
		word[0] = '\0';
	    }
	    newline = isNewlineChar(c);
	}
	else {
	    for (unsigned int i=0; i < WORDSIZE; ++i) {
		if (i == WORDSIZE-1) {
		    fprintf(stderr,"Word buffer owerflow!: %c.\n", c);
		}
		else if (word[i] == '\0') {
		    word[i] = c;
		    word[i+1] = '\0';
		    break;
		}
	    }
	}
    }
    return first_contig;
}

void print_seq ( FILE* out, unsigned int length, char* seq) {
    for (unsigned int i=0; i < length; ++i) {
	fputc(seq[i], out);
	if (!((i+1)%50) || i == length-1)
	    fputc('\n', out);
    }
}

void add_traces( struct contig_node* contigs ) {
    while (contigs != NULL) {
	for (unsigned int i=0; i < contigs->data.n_reads; ++i) {
	    char filename[WORDSIZE];
	    filename[0] = '\0';
	    if (contigs->data.reads[i].description != NULL) {
		// find CHROMAT_FILE
		unsigned int start = 0;
		char mode = '0';
		fprintf(stderr, "Looking for file name in DS: %s\n", contigs->data.reads[i].description);
		for (unsigned int p=0; contigs->data.reads[i].description[p] != '\0'; ++p) {
		    if (p > 2 &&
			contigs->data.reads[i].description[p] == 'E' &&
	    		contigs->data.reads[i].description[p-1] == 'L' &&
    			contigs->data.reads[i].description[p-2] == 'I' &&
			contigs->data.reads[i].description[p-3] == 'F') {
		    	if (start == 0) start = p;
			if (p > 10 &&
			    contigs->data.reads[i].description[p-4] == '_' &&
			    contigs->data.reads[i].description[p-5] == 'T' &&
			    contigs->data.reads[i].description[p-6] == 'A' &&
			    contigs->data.reads[i].description[p-7] == 'M' &&
			    contigs->data.reads[i].description[p-8] == 'O' &&
			    contigs->data.reads[i].description[p-9] == 'R' &&
			    contigs->data.reads[i].description[p-10] == 'H' &&
			    contigs->data.reads[i].description[p-11] == 'C') {
				start = p;
				mode = 'A';
			}
		    }
		}
		if (start == 0)
		    fprintf(stderr,"Did not find file anotation for: %s\n", contigs->data.reads[i].name);
		unsigned int pos = 0;
		for (unsigned int p=start+1; ; ++p) {
		    //fprintf(stderr, "reading pos: %u | char: %c\n", p, contigs->data.reads[i].description[p]);
		    if ( contigs->data.reads[i].description[p] != ' ' &&
			contigs->data.reads[i].description[p] != ':' &&
			contigs->data.reads[i].description[p] != '\0' &&
			pos < WORDSIZE-1 &&
			contigs->data.reads[i].description[p] != '\n' &&
			contigs->data.reads[i].description[p] != '\r') {
			    filename[pos] = contigs->data.reads[i].description[p];
			    filename[++pos] = '\0';
		    }
		    else if (filename[0] != '\0' || contigs->data.reads[i].description[p] == '\0' )
			break;
		}
		if (mode != 'A' && filename[0] != '\0') {
		    fprintf(stderr,"Did not find chromatogram file name, but other file name: %s.\n",filename); 
		    //unsigned int p = 0;
		    unsigned int last_dot = 0;
		    for (unsigned int p = 0; filename[p] != '\0'; ++p) {
			if (filename[p] == '.')
			    last_dot = p;
		    }
		    if (last_dot+3 < WORDSIZE) {
			if (filename[last_dot+1] == 'a' &&
			    filename[last_dot+2] == 'b' &&
			    ( filename[last_dot+3] == '1' ||
			    filename[last_dot+3] == 'i'))
				mode = 'A';
			else
			    if (filename[last_dot+1] == '1' && last_dot > 3 && filename[last_dot-1] == 'd' && filename[last_dot-2] == 'h' && filename[last_dot-3]  == 'p' && filename[last_dot-4] == '.')
				last_dot -= 4;
			if (mode != 'A') {
			    filename[last_dot+1] = 'a';
			    filename[last_dot+2] = 'b';
			    filename[last_dot+3] = '1';
			    filename[last_dot+4] = '\0';
			    if (access( filename, F_OK ) != 0)
				filename[last_dot+3] = 'i';
			    fprintf(stderr, "Not ab1 or abi file name ending. Trying with %s\n", filename);
			}
		    }
		}
		if (filename[0] == '\0' || access( filename, F_OK ) != 0) {
		    fputs("No name of existing file with ab1 or abi file ending found. Trying with read name and ab1 ending.\n",stderr);
		    if (strlen(contigs->data.reads[i].name) < WORDSIZE-4) {
			strcpy(filename,contigs->data.reads[i].name);
			strcat(filename,".ab1");
			if ( access( filename, F_OK ) != 0 ) {
			    filename[strlen(filename)-1] = 'i';
			    if ( access( filename, F_OK ) != 0 )
				filename[0] = '\0';
			}
		    }
		}
		if (filename[0] == '\0')
		    fprintf(stderr,"Could not find chromatogram filename for: %s.\n", contigs->data.reads[i].name);
	       	else if ( access( filename, F_OK ) != 0 )
		    fprintf(stderr,"No tracefile for: %s.\n", filename);
		else
		    contigs->data.reads[i].read_data = read_reading(filename, 0);
	    }
	}
	contigs = contigs->next;
    }
}

void printACE( FILE* out, struct contig_node* contigs ) {
    fprintf(out,"AS %u %u\n", n_contigs(contigs), n_reads(contigs));
    while (contigs != NULL) {
	char rev = 'U';
	if (contigs->data.revcomp)
	    rev = 'C';
	fprintf(out,"\nCO %s %u %u 0 %c\n", contigs->data.name, contigs->data.contig_length, contigs->data.n_reads, rev);
	print_seq(out,contigs->data.contig_length,contigs->data.contig_seq);
	if (contigs->data.contig_qual != NULL) {
	    fputs("\nBQ\n",out);
	    for (unsigned int i=0; i < contigs->data.contig_length; ++i) {
		fprintf(out, " %u", contigs->data.contig_qual[i]);
		if (!((i+1)%50) || i == contigs->data.contig_length-1)
		    fputc('\n', out);
	    }

	}
	if (contigs->data.reads != NULL) {
	    fputc('\n', out);
	    for (unsigned int i = 0; i < contigs->data.n_reads; ++i) {
		fprintf(out, "AF %s", contigs->data.reads[i].name);
		if (contigs->data.reads[i].revcomp)
		    fputs(" C", out);
		else
		    fputs(" U", out);
		fprintf(out, " %d\n", contigs->data.reads[i].start);
	    }
	    if (contigs->data.segment !=  NULL)
		fprintf(out,"%s\n", contigs->data.segment);
	    for (unsigned int i = 0; i < contigs->data.n_reads; ++i) {
		fprintf(out, "\nRD %s %u 0 0\n", contigs->data.reads[i].name, contigs->data.reads[i].seq_length);
		if (contigs->data.reads[i].seq != NULL)
		    print_seq(out, contigs->data.reads[i].seq_length, contigs->data.reads[i].seq);
		fprintf(out, "\nQA %u %u %u %u\n", contigs->data.reads[i].qual_begin, contigs->data.reads[i].qual_end, contigs->data.reads[i].align_begin, contigs->data.reads[i].align_end);
		if (contigs->data.reads[i].description != NULL)
		    fprintf(out, "DS %s\n", contigs->data.reads[i].description);
	    }
	}
	contigs = contigs->next;
    }
}

void print_read_dataJSON ( FILE* out, Read* read ) {
    if (read != NULL) {
	fprintf(out,"{\"file_name\":\"%s\", \"n_bases\":%d, \"left_cut_off\":%d, \"right_cut_off\":%d, \"traces\":{\"n_points\":%d, \"max_val\":%u, \"base_line\":%d,\"A\":[", read->trace_name, read->NBases, read->leftCutoff, read->rightCutoff, read->NPoints, read->maxTraceVal, read->baseline);
	for (unsigned int i=0; i < read->NPoints; ++i) {
	    if (i!=0) fputc(',',out);
	    fprintf(out,"%u",read->traceA[i]);
	}
	fputs("],\"C\":[",out);
	for (unsigned int i=0; i < read->NPoints; ++i) {
            if (i!=0) fputc(',',out);
            fprintf(out,"%u",read->traceC[i]);
        }
	fputs("],\"G\":[",out);
	for (unsigned int i=0; i < read->NPoints; ++i) {
            if (i!=0) fputc(',',out);
            fprintf(out,"%u",read->traceG[i]);
        }
	fputs("],\"T\":[",out);
	for (unsigned int i=0; i < read->NPoints; ++i) {
            if (i!=0) fputc(',',out);
            fprintf(out,"%u",read->traceT[i]);
        }
	fputs("]}", out);
	fputs(", \"calls\":{ \"seq\":[",out);
	for (unsigned int i=0; i < read->NBases; ++i) {
	    if (i!=0) fputc(',',out);
	    fprintf(out,"\"%c\"", read->base[i]);// , read->basePos[i], read->prob_A[i], read->prob_C[i], read->prob_G[i], read->prob_T[i]);
	}
	fputs("], \"trace_pos\":[",out);
	for (unsigned int i=0; i < read->NBases; ++i) {
	    if (i!=0) fputc(',',out);
	    fprintf(out,"%u", read->basePos[i]);
	}
	fputs("], \"base_prob\":[",out);
	for (unsigned int i=0; i < read->NBases; ++i) {
	    if (i!=0) fputc(',',out);
	    fprintf(out,"{\"A\":%u,\"C\":%u,\"G\":%u,\"T\":%u}", read->prob_A[i], read->prob_C[i], read->prob_G[i], read->prob_T[i]);
	}
	fputs("]}", out);
	fputc('}', out);
    }
}

void print_readJSON ( FILE* out, contig_read* read ) {
    fprintf(out,"{\"name\":\"%s\", \"length\":%u, \"start\":%d, ", read->name, read->seq_length, read->start);
    fputs("\"reverse\":",out);
	if (read->revcomp)
	    fputs("true",out);
        else
            fputs("false",out);
    fprintf(out,", \"qual_begin\":%u, \"qual_end\":%u, \"align_begin\":%u, \"align_end\":%u", read->qual_begin, read->qual_end, read->align_begin, read->align_end);
    if (read->read_data !=0) {
	fputs(", \"chrom\":", out);
	print_read_dataJSON( out, read->read_data);
    }
    else {
	fputs(", \"seq\":[", out);
	for (unsigned int i=0; i < read->seq_length; ++i) {
	    if (read->seq[i] != '*') { 
		if (i != 0)
		    fputc(',',out);
		fprintf(out,"\"%c\"",read->seq[i]);
	    }
	}
	fputc(']', out);
    }
    fprintf(out,", \"description\":\"%s\"", read->description);
    fputc('}', out);
}

void printJSON ( FILE* out, struct contig_node* contigs ) {
    fputs("{\"contigs\":[",out);
    _Bool first = true;
    while (contigs != NULL) {
	if (first)
	    first = false;
	else
	    fputc(',',out);
	fputc('{',out);
	fprintf(out,"\"name\":\"%s\", ", contigs->data.name);
	fprintf(out,"\"length\":%u, ", contigs->data.contig_length);
	//fprintf(out,"\"n_gaps\":%u, ", contigs->data.n_gaps);
	fprintf(out,"\"n_reads\":%u, ", contigs->data.n_reads);
	fputs("\"reverse\":",out);
	if (contigs->data.revcomp)
	    fputs("true, ",out);
	else
	    fputs("false, ",out);
	fputs("\"seq\":[", out);
	for (unsigned int i=0; i < contigs->data.contig_length; ++i) {
	    if (contigs->data.contig_seq[i] != '*') {
		if (i != 0)
		    fputc(',',out);
		fprintf(out,"\"%c\"",contigs->data.contig_seq[i]);
	    }
	}
	fputs("],\"qual\":[ ", out);
	for (unsigned int i=0; i < contigs->data.contig_length; ++i) {
	    if (contigs->data.contig_seq[i] != '*') {
		if (i != 0)
		    fputc(',',out);
		fprintf(out,"%u", contigs->data.contig_qual[i]);
	    }
	}
	fputs("],\"reads\":[ ", out);
	int align_start = 1;
	for (unsigned int i=0; i < contigs->data.n_reads; ++i) {
	    if (i != 0)
                fputs(", ",out);
	    print_readJSON(out, &(contigs->data.reads[i]));
	    if ( contigs->data.reads[i].start < align_start)
		align_start = contigs->data.reads[i].start;;
	}
	fprintf(out, " ], \"alignment\":{\"contig\":{ \"start\":%i, \"gaps\":[", align_start*-1+1);
	unsigned int n_gap = 0;
	unsigned int ali_length = contigs->data.contig_length + align_start;
	for (unsigned int i=0; i < contigs->data.contig_length; ++i) {
	    if (contigs->data.contig_seq[i] == '*') {
		if (n_gap)
		    fputc(',',out);
		++n_gap;
		fprintf(out,"{\"pos\":%u, \"qual\":[%u", i-n_gap, contigs->data.contig_qual[i]);
		unsigned int length = 1;
		while (i+1 < contigs->data.contig_length && contigs->data.contig_seq[i+1] == '*') {
		    ++i;
		    ++length;
		    fprintf(out,",%i",contigs->data.contig_qual[i]);
		}
		n_gap += length-1;
		fprintf(out,"], \"length\":%u}", length);
	    }
	}
	fputs("]},\"reads\":[", out);
	for (unsigned int i=0; i < contigs->data.n_reads; ++i) {
	    if (i != 0)
		fputc(',',out);
	    n_gap = 0;
	    fprintf(out, "{ \"start\":%i, \"gaps\":[", contigs->data.reads[i].start-align_start);
	    if (contigs->data.reads[i].start-align_start + contigs->data.reads[i].seq_length > ali_length)
		ali_length = contigs->data.reads[i].start-align_start + contigs->data.reads[i].seq_length;
	    for (unsigned int j=0; j < contigs->data.reads[i].seq_length; ++j) {
		if (contigs->data.reads[i].seq[j] == '*') {
		    if (n_gap)
			fputc(',',out);
		    ++n_gap;
		    fprintf(out,"{\"pos\":%u", j-n_gap);
		    unsigned int length = 1;
		    while (j+1 < contigs->data.contig_length && contigs->data.contig_seq[j+1] == '*') {
			++j;
			++length;
		    }
		    n_gap += length-1;
		    fprintf(out,", \"length\":%u}", length);
		}
	    }
	    fputs("]}", out);
	}
	fprintf(out, "],\"length\":%u}}", ali_length);
	
	if (contigs->next != NULL)
	    fputc(',',out);
	contigs = contigs->next;
    }
    fputs("]}\n",out);
}
