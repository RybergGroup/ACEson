#include "alignment.h"

struct multiple_sequences* multiple_sequences_alocate( const unsigned int num_seq ) {
    struct multiple_sequences* ms = (struct multiple_sequences*)calloc(1, sizeof (struct multiple_sequences));
    ms->Nseq = num_seq;
    ms->name = (char**)calloc(num_seq,sizeof(char*));
    ms->Nbases = (unsigned int*)calloc(num_seq,sizeof(unsigned int));
    ms->leftCutoff = (unsigned int*)calloc(num_seq,sizeof(unsigned int));
    ms->rightCutoff = (unsigned int*)calloc(num_seq,sizeof(unsigned int));
    ms->base = (char**)calloc(num_seq,sizeof(char*));
    ms->prob = (struct base_prob*)calloc(num_seq,sizeof(struct base_prob));
    return ms;
}

void multiple_sequences_dealocate ( struct multiple_sequences* ms ) {
    if (ms) {
	if (ms->name)
	    free(ms->name);
	if (ms->Nbases)
	    free(ms->Nbases);
	if (ms->leftCutoff)
	    free(ms->leftCutoff);
	if (ms->rightCutoff)
	    free(ms->rightCutoff);
	if (ms->base)
	    free(ms->base);
	if (ms->prob)
	    free(ms->prob);
	free(ms);
    }
}

struct alignment* alignment_alocate( struct multiple_sequences* multiple_seq ) {
    struct alignment* ali = (struct alignment*)calloc(1, sizeof (struct alignment));
    ali->seq = multiple_seq;
    ali->Nseq = multiple_seq->Nseq;
    unsigned int max_seq_length = 0;
    for (unsigned int i = 0; i < multiple_seq->Nseq; ++i)
	if (multiple_seq->Nbases[i] > max_seq_length)
	    max_seq_length = multiple_seq->Nbases[i];
    ali->length = max_seq_length;
    ali->startPos = (unsigned int*)calloc(ali->Nseq,sizeof(unsigned int));
    ali->revcomp = (_Bool*)calloc(ali->Nseq,sizeof(_Bool));
    ali->gap = (struct gap_stack**)calloc(ali->Nseq,sizeof(struct gap_stack*));
    return ali;
}

void gap_stack_dealocate ( struct gap_stack* present ) {
    if (present && present->next)
	gap_stack_dealocate(present->next);
    if (present)
	free(present);
}

struct gap_stack* copy_gap_stack ( struct gap_stack* present ) {
    if (present == 0)
	return 0;
    struct gap_stack* copy = (struct gap_stack*)malloc(sizeof(struct gap_stack));
    copy->pos = present->pos;
    copy->length = present->length;
    copy->next = copy_gap_stack(present->next);
    return copy;
}

void alignment_dealocate ( struct alignment* ali ) {
    if (ali) { 
	if (ali->startPos)
	    free(ali->startPos);
	if (ali->gap) {
	    for (unsigned int i=0; i < ali->Nseq; ++i)
		gap_stack_dealocate(ali->gap[i]);
	    free(ali->gap);
	}
	free(ali);
    }
}

void alignment_ms_dealocate ( struct alignment* ali) {
    fputs("Check A\n", stderr);
    multiple_sequences_dealocate(ali->seq);
    fputs("Check B\n", stderr);
    alignment_dealocate(ali);
}

struct score_matrix* score_matrix_alocate ( unsigned int length_seq_one, unsigned int length_seq_two ) {
    struct score_matrix* matrix = (struct score_matrix*)malloc(sizeof(struct score_matrix));
    matrix->Ncols = length_seq_one;
    matrix->Nrows = length_seq_two;
    matrix->score = (int*)calloc(length_seq_one*length_seq_two,sizeof(int));
    matrix->prevOne = (int*)calloc(length_seq_one*length_seq_two,sizeof(int));
    matrix->prevTwo = (int*)calloc(length_seq_one*length_seq_two,sizeof(int));
    matrix->maxScore = 0;
    matrix->endI = 0;
    matrix->endJ = 0;
    return matrix;
} 

void score_matrix_dealocate( struct score_matrix* matrix ) {
    if (matrix) {
	if (matrix->score)
	    free(matrix->score);
	if (matrix->prevOne)
	    free(matrix->prevOne);
	if (matrix->prevTwo)
	    free(matrix->prevTwo);
	free(matrix);
    }
}

unsigned int get_matrix_pos (struct score_matrix* matrix,  unsigned int x, unsigned int y ) {
    if (y > matrix->Nrows || x > matrix->Ncols)
	return matrix->Ncols * matrix->Nrows;
    return y*matrix->Ncols + x;
}

void add_sequence ( struct multiple_sequences* ms, const unsigned int seq_num, char* name, const unsigned int length, char* seq, char* qual_G, char* qual_C, char* qual_T, char* qual_A, const unsigned int trim_from_start, const unsigned int trim_from_end ) {
    if (seq_num >= ms->Nseq)
	return;
    ms->name[seq_num] = name;
    ms->Nbases[seq_num] = length;
    ms->leftCutoff[seq_num] = trim_from_start;
    ms->rightCutoff[seq_num] = trim_from_end;
    ms->base[seq_num] = seq;
    ms->prob[seq_num].G = qual_G;
    ms->prob[seq_num].C = qual_C;
    ms->prob[seq_num].T = qual_T;
    ms->prob[seq_num].A = qual_A;
}

void add_read ( struct multiple_sequences* ms, Read* read, const unsigned int seq_num) {
    add_sequence( ms, seq_num, read->trace_name, read->NBases, read->base, read->prob_G, read->prob_C, read->prob_T, read->prob_A, read->leftCutoff, read->rightCutoff);
}

void add_gap ( struct alignment* ali, const unsigned int seq_num, const unsigned int pos, const unsigned int length ) {
    unsigned int ali_seq_length = ali->seq->Nbases[seq_num]+ali->startPos[seq_num];
    if (ali->gap[seq_num] == 0 ) {
	ali->gap[seq_num] = (struct gap_stack*)malloc(sizeof(struct gap_stack));
	ali->gap[seq_num]->pos = pos;
	ali->gap[seq_num]->length = length;
	ali->gap[seq_num]->next = 0;
    }
    else {
	struct gap_stack* present = ali->gap[seq_num];
	struct gap_stack* prev = 0;
	// Loop while in a position after the present gap, and not at end
	while ( pos > present->pos && present->next !=0 ) {
	    ali_seq_length += present->length;
	    prev = present;
	    present = present->next;
	}
	ali_seq_length += present->length; // Why?
	if (pos < present->pos) { // if before a gap, new gap is needed
	    struct gap_stack* new_gap_pos = (struct gap_stack*)malloc(sizeof(struct gap_stack));
	    new_gap_pos->pos = pos;
	    new_gap_pos->length = length;
	    new_gap_pos->next = present;
	    if (prev == 0)
		ali->gap[seq_num] = new_gap_pos;
	    else
		prev->next = new_gap_pos;
	}
	else if (pos == present->pos) { // if in a gap extend it
	    present->length += length;
	    present = present->next;
	}
	else if (present->next == 0) { // if at end, add new gap at end
	    present->next = (struct gap_stack*)malloc(sizeof(struct gap_stack));
	    present->next->pos = pos;
	    present->next->length = length;
	    present->next->next = 0;
	    present = 0;
	}
	else
	    fputs("Unforeseen event when inserting gap.\n", stderr);
	while (present != 0) {
	    ali_seq_length += present->length;
	    present = present->next;
	}
    }
    if (ali_seq_length > ali->length)
	ali->length = ali_seq_length;
	
}

void insert_gap(struct alignment* ali, const unsigned int ali_pos, const unsigned int gap_length) {
    for (unsigned int i=0; i < ali->Nseq; ++i) {
	if ((ali->revcomp[i] && ali_pos < ali->startPos[i]+ali->seq->rightCutoff[i]) || ali_pos < ali->startPos[i]+ali->seq->leftCutoff[i]) {
	    ali->startPos[i] += gap_length;
	    //fprintf(stderr,"Insert %u gap at beginning", gap_length);
	}
	else {
	    unsigned int pos = get_seq_pos( ali, i, ali_pos);
	    //fprintf(stderr,"%u %u %u\n", pos, ali->seq->Nbases[i]-1, gap_length);
	    if (pos < ali->seq->Nbases[i]-1 && gap_length > 0) {
		add_gap(ali, i, pos, gap_length);
		//fprintf(stderr,"Insert %u gap at %u (%u)", gap_length, ali_pos, pos);
	    }
	}
    }
}

void move_start_pos (struct alignment* ali, const int n) {
    fputs("Changing start pos.\n", stderr);
    for (unsigned int i=0; i < ali->Nseq; ++i) {
	if ((int)ali->startPos[i] > n*-1)
	    ali->startPos[i] += n;
	else
	    ali->startPos[i] = 0;
    }
    ali->length = alignment_length(ali);
}

unsigned int get_seq_pos ( struct alignment* ali, const unsigned int seq_num, const unsigned int align_pos) {
    unsigned int seq_pos = align_pos - ali->startPos[seq_num];
    //fprintf(stderr,"X: %u %u %u : %u\n", align_pos, ali->startPos[seq_num], seq_pos, ali->seq->Nbases[seq_num]-ali->seq->rightCutoff[seq_num]);
    struct gap_stack* present = ali->gap[seq_num];
    while (present != 0) {
        if (present->pos > seq_pos)
            break;
        else if (present->pos+present->length > seq_pos) {
	    //fprintf(stderr, "---- %u --- %u\n", present->pos+present->length, seq_pos);
            return ali->seq->Nbases[seq_num];
	}
        else
            seq_pos -= present->length;
	present = present->next;
    }
    if (seq_pos > ali->seq->Nbases[seq_num]-ali->seq->rightCutoff[seq_num])
        return ali->seq->Nbases[seq_num];
    if (ali->revcomp[seq_num])
	seq_pos = ali->seq->Nbases[seq_num]-1-seq_pos;
    return seq_pos;
}

char get_call_prob ( struct alignment* ali, const unsigned int seq_num, const unsigned int seq_pos ) {
    if (ali->seq->base[seq_num][seq_pos] == 'G' || ali->seq->base[seq_num][seq_pos] == 'g')
	return ali->seq->prob[seq_num].G[seq_pos];
    else if (ali->seq->base[seq_num][seq_pos] == 'C' || ali->seq->base[seq_num][seq_pos] == 'c')
	return ali->seq->prob[seq_num].C[seq_pos];
    else if (ali->seq->base[seq_num][seq_pos] == 'T' || ali->seq->base[seq_num][seq_pos] == 't')
	return ali->seq->prob[seq_num].T[seq_pos];
    else if (ali->seq->base[seq_num][seq_pos] == 'A' || ali->seq->base[seq_num][seq_pos] == 'a')
	return ali->seq->prob[seq_num].A[seq_pos];
    else
	return 0x00;
}

struct base_call get_base_call ( struct alignment* ali, const unsigned int seq_num, const unsigned int align_pos ) {
    struct base_call base;
    base.base = '-';
    base.prob = 0x00;
    if (seq_num >= ali->Nseq)
        return base;
    unsigned int seq_pos = get_seq_pos ( ali, seq_num, align_pos );
    if (seq_pos >= ali->seq->Nbases[seq_num])
	return base;
    base.base = ali->seq->base[seq_num][seq_pos];
    if (ali->seq->prob != 0)
	base.prob = get_call_prob(ali,seq_num,seq_pos);
    return base;
}

int get_gap_score ( struct alignment* ali, unsigned int seq_num, unsigned int align_pos) {
    align_pos -= ali->startPos[seq_num];
    struct gap_stack* present = ali->gap[seq_num];
    while (present != 0) {
	if (present->pos+present->length > align_pos) {
	    align_pos = present->pos;
	    break;
	}
	else if (present->pos < align_pos)
	    return 0;
	present = present->next;
    }
    if (!present)
	return 0;
    if (align_pos > ali->seq->Nbases[seq_num]-2)
	return 0;
    if (ali->revcomp)
	align_pos = ali->seq->Nbases[seq_num]-2-align_pos;
    return ((int)(get_call_prob(ali,seq_num,align_pos))+(int)(get_call_prob(ali,seq_num,align_pos+1)))/2;
}

struct base_call get_con_base ( struct alignment* ali, const unsigned int align_pos ) {
    struct base_call base;
    int G=0, C=0, T=0, A=0, gap=0;
    for (unsigned int i=0; i<ali->Nseq; ++i) {
	base = get_base_call(ali, i, align_pos);
	if (ali->revcomp[i])
	    base.base = comp_base(base.base);
	if (base.base=='G' || base.base=='g')
	    G+=base.prob;
	else if (base.base=='C' || base.base=='c')
	    C+=base.prob;
	else if (base.base=='T' || base.base=='t')
	    T+=base.prob;
	else if (base.base=='A' || base.base=='a')
	    A+=base.prob;
	else if (base.base=='-')
	    gap += get_gap_score(ali, i, align_pos);
    }
    if (gap >= G && gap >= C && gap >= T && gap >= A) {
	base.base = '-';
     	base.prob = gap-G-T-A-C;
    }
    else if (G >= C && G >=T && G >= A && G>gap) {
	base.base = 'G';
	base.prob = G-T-A-C-gap;
    }
    else if (C >=T && C >= A && C > G && C > gap) {
	base.base = 'C';
	base.prob = C-T-G-A-gap;
    }
    else if (T >=A && T > C && T > A && T > gap) {
	base.base = 'T';
	base.prob = T-A-C-G-gap;
    }
    else {
	base.base = 'A';
	base.prob = A-G-T-C-gap;
    }
    return base;
}

unsigned int alignment_length ( struct alignment* ali ) {
    unsigned int max_length = 0;
    for (unsigned int i=0; i < ali->Nseq; ++i) {
	unsigned int length = ali->seq->Nbases[i]+ali->startPos[i];
	struct gap_stack* present = ali->gap[i];
	while (present != 0) {
	    length += present->length;
	    present = present->next;
	}
	if (length > max_length)
	    max_length = length;
    }
    return max_length;
}

char comp_base ( char base ) {
    if (base == 'G')
	base = 'C';
    else if (base == 'A')
	base = 'T';
    else if (base == 'C')
	base = 'G';
    else if (base == 'T')
	base = 'A';
    return base;
}

void set_max_score ( struct score_matrix* matrix, int score, unsigned int i, unsigned int j) {
    if (score > matrix->maxScore) {
	matrix->maxScore = score;
	matrix->endI = i;
	matrix->endJ = j;
    }
}

struct multiple_sequences* join_sequences ( struct multiple_sequences* seq1, struct multiple_sequences* seq2 ) {
    struct multiple_sequences* sequences = multiple_sequences_alocate(seq1->Nseq+seq2->Nseq);
    for (unsigned int i=0; i <  seq1->Nseq; ++i) {
	sequences->name[i] = seq1->name[i];
	sequences->Nbases[i] = seq1->Nbases[i];
	sequences->leftCutoff[i] = seq1->leftCutoff[i];
	sequences->rightCutoff[i] = seq1->rightCutoff[i];
	sequences->base[i] = seq1->base[i];
	sequences->prob[i] = seq1->prob[i];
    }
    for (unsigned int i=0; i <  seq2->Nseq; ++i) {
	sequences->name[i+seq1->Nseq] = seq2->name[i];
	sequences->Nbases[i+seq1->Nseq] = seq2->Nbases[i];
	sequences->leftCutoff[i+seq1->Nseq] = seq2->leftCutoff[i];
	sequences->rightCutoff[i+seq1->Nseq] = seq2->rightCutoff[i];
	sequences->base[i+seq1->Nseq] = seq2->base[i];
	sequences->prob[i+seq1->Nseq] = seq2->prob[i];
    }
    return sequences;
}

struct alignment* merged_alignment(struct alignment* ali_one, struct alignment* ali_two) {
    struct alignment* merged_alignment = alignment_alocate(join_sequences(ali_one->seq,ali_two->seq));
    fprintf(stderr, "length one: %d, length two: %d\n", ali_one->length, ali_two->length);
    if (ali_one->length > ali_two->length)
	merged_alignment->length = ali_one->length;
    else
	merged_alignment->length = ali_two->length;
    for (unsigned int i=0; i <  ali_one->Nseq; ++i) {
	merged_alignment->startPos[i] = ali_one->startPos[i];
	merged_alignment->revcomp[i] = ali_one->revcomp[i];
	merged_alignment->gap[i] = copy_gap_stack(ali_one->gap[i]);
    }
    for (unsigned int i=0; i <  ali_one->Nseq; ++i) {
	merged_alignment->startPos[i+ali_one->Nseq] = ali_two->startPos[i];
	merged_alignment->revcomp[i+ali_one->Nseq] = ali_two->revcomp[i];
	merged_alignment->gap[i+ali_one->Nseq] = copy_gap_stack(ali_two->gap[i]);
    }
    fprintf(stderr, "length mearged: %d\n", merged_alignment->length);
    return merged_alignment;
}


void revcomp (struct alignment* ali) {
    for (unsigned int i = 0; i < ali->Nseq; ++i) {
	fprintf(stderr,"Start: %d\n", ali->startPos[i]);
	ali->startPos[i] = ali->length - ali->startPos[i] -  ali->seq->Nbases[i];
	fprintf(stderr,"Start rev: %d\n", ali->startPos[i]);
	ali->revcomp[i] = !ali->revcomp[i];
	struct gap_stack* present = ali->gap[i];
	struct gap_stack* prev = 0;
	struct gap_stack* temp = 0;
	while (present != 0) {
	    present->pos = ali->seq->Nbases[i]- 1 - present->pos;
	    temp = present;
	    present = present->next;
	    if (present == 0)
		ali->gap[i] = temp;
	    else {
		present->next = prev;
		prev = temp;
	    }
	    //fputc('m', stderr);
	}
    }
}

struct alignment* align_pair ( struct alignment* ali_one, struct alignment* ali_two ) {
    struct score_matrix* for_score = score_matrix_alocate(ali_one->length, ali_two->length);
    struct score_matrix* rev_score = score_matrix_alocate(ali_one->length, ali_two->length);
    for (unsigned int i=1; i < ali_one->length; ++i) { // loop over alignment pos in "seq1"
	struct base_call one_base = get_con_base(ali_one, i);
	for (unsigned int j=1; j < ali_two->length; ++j) { // loop over alignment pos in "seq2"
	    int match = for_score->score[get_matrix_pos(for_score, i-1, j-1)];
	    int gap_j = for_score->score[get_matrix_pos(for_score, i, j-1)];
	    int gap_i = for_score->score[get_matrix_pos(for_score,i-1,j)];
	    int rev_match = rev_score->score[get_matrix_pos(for_score,i-1,j-1)];
	    int rev_gap_j = rev_score->score[get_matrix_pos(for_score,i,j-1)];
	    int rev_gap_i = rev_score->score[get_matrix_pos(for_score,i-1,j)];
	    struct base_call two_base = get_con_base(ali_two, j);
	    struct base_call two_base_rev = get_con_base(ali_two, ali_two->length-1-j);
	    two_base_rev.base = comp_base(two_base_rev.base);
	    // Calc prob if aligning agains each other
	    if (one_base.base == two_base.base) {
	     	match += one_base.prob + two_base.prob;
	    }
     	    else {
 		match -= (one_base.prob + two_base.prob);
	    }
	    if (one_base.base == two_base_rev.base) {
		rev_match += one_base.prob + two_base_rev.prob;
	    }
	    else {
		rev_match -= (one_base.prob + two_base_rev.prob);
	    }
	    // Calc prob if inserting gaps
	    gap_j -= 2*two_base.prob;
	    gap_i -= 2*one_base.prob;
	    rev_gap_j -= 2*two_base_rev.prob;
	    rev_gap_i -= 2*one_base.prob;
	    // Check which option is best and set score
	    if (match >= gap_j && match >= gap_i) {
		for_score->score[get_matrix_pos(for_score,i,j)] = match;
		for_score->prevOne[get_matrix_pos(for_score,i,j)] = i-1;
		for_score->prevTwo[get_matrix_pos(for_score,i,j)] = j-1;
		set_max_score(for_score, match, i, j);
	    }
	    else if (gap_j > match && gap_j >= gap_i) {
		for_score->score[get_matrix_pos(for_score,i,j)] = gap_j;
		for_score->prevOne[get_matrix_pos(for_score,i,j)] = i;
		for_score->prevTwo[get_matrix_pos(for_score,i,j)] = j-1;
	    }
	    else if (gap_i > match && gap_i > gap_j) {
		for_score->score[get_matrix_pos(for_score,i,j)] = gap_j;
		for_score->prevOne[get_matrix_pos(for_score,i,j)] = i-1;
		for_score->prevTwo[get_matrix_pos(for_score,i,j)] = j;
	    }
	    else {
		for_score->score[get_matrix_pos(for_score,i,j)] = 0;
		for_score->prevOne[get_matrix_pos(for_score,i,j)] = -1;
		for_score->prevTwo[get_matrix_pos(for_score,i,j)] = -1;
	    }
	    if (rev_match >= rev_gap_j && rev_match >= rev_gap_i) {
		rev_score->score[get_matrix_pos(rev_score,i,j)] = rev_match;
		rev_score->prevOne[get_matrix_pos(rev_score,i,j)] = i-1;
		rev_score->prevTwo[get_matrix_pos(rev_score,i,j)] = j-1;
		set_max_score(rev_score, rev_match, i, j);
	    }
	    else if (rev_gap_j > rev_match && rev_gap_j >= rev_gap_i) {
		rev_score->score[get_matrix_pos(rev_score,i,j)] = rev_gap_j;
		rev_score->prevOne[get_matrix_pos(rev_score,i,j)] = i;
		rev_score->prevTwo[get_matrix_pos(rev_score,i,j)] = j-1;
	    }
	    else if (rev_gap_i > rev_match && rev_gap_i > rev_gap_j) {
		rev_score->score[get_matrix_pos(rev_score,i,j)] = rev_gap_i;
		rev_score->prevOne[get_matrix_pos(rev_score,i,j)] = i-1;
		rev_score->prevTwo[get_matrix_pos(rev_score,i,j)] = j;
	    }
	    else {
		rev_score->score[get_matrix_pos(rev_score,i,j)] = 0;
		rev_score->prevOne[get_matrix_pos(rev_score,i,j)] = -1;
		rev_score->prevTwo[get_matrix_pos(rev_score,i,j)] = -1;
	    }
	}
    }
    struct score_matrix* best_score = 0;
    _Bool revcomp_two = 0;
    // get pos in seqs for best score
    if (for_score->maxScore >= rev_score->maxScore)
	best_score = for_score;
    else {
	best_score = rev_score;
	revcomp_two = 1;
    }
    unsigned int i = best_score->endI;
    unsigned int j = best_score->endJ;
    // 
    while ( 1 ) {
	//fprintf(stderr, "p: %d - %d\n", i, j);
	unsigned int next_i = best_score->prevOne[get_matrix_pos(best_score,i,j)];
	unsigned int next_j = best_score->prevTwo[get_matrix_pos(best_score,i,j)];
	/*/// DEBUG /////
	struct base_call first = get_con_base(ali_one,i);
    	struct base_call second;
	if (revcomp_two) {
	    second=get_con_base(ali_two,ali_two->length-j);
	    second.base = comp_base(second.base);
	}
	else
	    second=get_con_base(ali_two,j);
	////////////////*/
	if (next_i == i && next_j == j) {
	    fprintf(stderr, "Error adding gap to both seqs - i:  %d j: %d\n", i, j);
	}
	else if (next_i == i) {
	    insert_gap(ali_one, i, 1);
	    //fprintf(stderr, "Gap: - - %c\n", second.base);
	}
	else if (next_j == j) {
	    if (revcomp_two) {
		insert_gap(ali_two, ali_two->length-1-j, 1);
	    }
	    else {
		insert_gap(ali_two, j, 1);
		//fprintf(stderr, "j:  %d\n", j);
	    }
	    //fprintf(stderr, "Gap: %c - -\n", first.base);
	}
	/*else {
    	    fprintf(stderr, "Match: %c - %c\n", first.base, second.base);
	}*/
	if (next_i == 0 || next_j == 0)
	    break;
	i = next_i;
	j = next_j;
    }
    //fprintf(stderr, "f: %u - %u\n", i, j);
    //fprintf(stderr, "L: %u - %u\n", ali_one->length, ali_two->length);
    //fprintf(stderr, "S: %u - %u\n", ali_one->startPos[0], ali_two->startPos[0]);
    if (j>0)
	move_start_pos(ali_one, j);
    else if (i>0) {
	if (revcomp_two)
	    insert_gap(ali_two, ali_two->length-1, i);
	else
	    move_start_pos(ali_two, i);
    }	    
    //fprintf(stderr, "L: %ui - %u\n", ali_one->length, ali_two->length);
    //fprintf(stderr, "S: %u - %u\n", ali_one->startPos[0], ali_two->startPos[0]);
    if (revcomp_two)
	revcomp(ali_two);
    fputs("Made it here.\n", stderr);
    return merged_alignment(ali_one, ali_two);
}

char get_seq_base ( struct alignment* ali, unsigned int seq_num, unsigned int pos ) {
    if (ali->revcomp[seq_num])
	return comp_base(ali->seq->base[seq_num][ali->seq->Nbases[seq_num]-1-pos]);
    else
	return ali->seq->base[seq_num][pos];
}

void fprint_fasta_alignment (struct alignment* ali, FILE* out ) {
    for (unsigned int i = 0; i < ali->Nseq; ++i) {
	fprintf(out, ">%s\n", ali->seq->name[i]);
	unsigned int j;
	unsigned int length = ali->length;
	//fprintf(stderr, "Alignment length: %u\n", length);
	for (j=0; j < ali->startPos[i]; ++j) {
	    fputc('-', out);
	    if (length >0)
		--length;
	}
	j=0;
	//fprintf(stderr, "Left to print: %u\n", length);
	struct gap_stack* present = ali->gap[i];
	while (present != 0) {
	    for (;j<=present->pos && j < ali->seq->Nbases[i]; ++j) {
		fputc(get_seq_base(ali,i,j),out);
		if (length >0)
		    --length;
	    }
	    for (unsigned int k=0; k < present->length; ++k) {
		fputc('-', out);
		if (length >0)
		    --length;
	    }
	    present = present->next;
	}
	//fprintf(stderr, "Left to print: %u\n", length);
	for (;j < ali->seq->Nbases[i]; ++j) {
	    fputc(get_seq_base(ali,i,j),out);
	    if (length >0)
		--length;
	}
	//fprintf(stderr, "Left to print: %u\n", length);
	while (length > 0) {
	    fputc('-', out);
	    --length;
	}
	fputc('\n', out);
	//fprintf(stderr, "Nothing left to print: %d\n", length);
    }
}

contig alignment_to_contig (char* name, struct alignment* ali, Read** read_data ) {
    contig data;
    data.name = (char*)malloc(sizeof(char)*(strlen(name)+1));
    strcpy(data.name, name); //
    unsigned int contig_length = ali->length;
    unsigned int i=0;
    while (i < ali->length) {
	struct base_call nuc = get_con_base(ali,i);
	if (nuc.base != '-')
	    break;
	++i;
    }
    unsigned int contig_start = i;
    contig_length -= i;
    i = 0;
    while (i < ali->length) {
	struct base_call nuc = get_con_base(ali,ali->length-i-1);
	if (nuc.base != '-')
            break;
        ++i;
    }
    contig_length -= i;
    data.contig_length = contig_length; //
    data.contig_seq = (char*)malloc(sizeof(char)*contig_length);
    data.contig_qual = (unsigned int*)malloc(sizeof(unsigned int)*contig_length);
    unsigned int n_gaps;
    for (i=0; i < contig_length; ++i) {
	struct base_call nuc = get_con_base(ali,i+contig_start);
	if (nuc.base == '-') {
	    ++n_gaps;
	    data.contig_seq[i] = '*'; //-
	}
	else
	    data.contig_seq[i] = nuc.base; //-
	data.contig_qual[i] = nuc.prob; //
	
    }
    data.n_gaps = n_gaps;
    data.revcomp = 0;
    data.n_reads = ali->seq->Nseq;
    data.reads = (contig_read*)calloc(data.n_reads,sizeof(contig_read));
    for (i=0; i < ali->seq->Nseq; ++i) {
	data.reads[i].name = (char*)malloc(sizeof(char)*(strlen(ali->seq->name[i])+1));
	strcpy(data.reads[i].name,ali->seq->name[i]); //
	data.reads[i].seq_length = ali->seq->Nbases[i]; //
	data.reads[i].seq = (char*)malloc(sizeof(char)*data.reads[i].seq_length); // ali->seq.base[i]; //
	for (unsigned int j=0; j < ali->seq->Nbases[i]; ++j) {
	    if (ali->seq->base[i][j] == '-')
		data.reads[i].seq[j] = '*';
	    else
		data.reads[i].seq[j] = ali->seq->base[i][j];
	} //<-
	data.reads[i].start = ali->startPos[i]-(int)contig_start+1; //
	data.reads[i].qual_begin = ali->seq->leftCutoff[i]; //
	data.reads[i].qual_end = ali->seq->rightCutoff[i]; //
	data.reads[i].align_begin = ali->seq->leftCutoff[i]; //
	data.reads[i].align_end = ali->seq->leftCutoff[i]; //
	data.reads[i].description = NULL; //
	data.reads[i].revcomp = ali->revcomp[i]; //
	if (read_data) {
	    if (!strcmp(ali->seq->name[i],read_data[i]->trace_name))
                    data.reads[i].read_data = read_data[i];
    	    else {
		for ( unsigned int k=0; k < ali->seq->Nseq; ++k) {
		    if (k != i && !strcmp(ali->seq->name[i],read_data[k]->trace_name)) {
			data.reads[i].read_data = read_data[k];
			break;
		    }
		}
	    }
	}
    }
    return data;
}
