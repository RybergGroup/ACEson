#include "alignment.h"

/*void add_sequence ( struct multiple_sequences* ms, const unsigned int seq_num, char* name, const unsigned int length, char* seq, char* qual_G, char* qual_C, char* qual_T, char* qual_A, const unsigned int trim_from_start, const unsigned int trim_from_end ) {
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
}*/

void alignment::add_gap ( const unsigned int seq_num, unsigned int ali_pos, const unsigned int length ) {
    ali_pos -= start_pos[seq_num];
    //cerr << "Ali pos: " << ali_pos << endl;
    vector<pair<unsigned int, unsigned int>>::iterator pos;
    for (pos = gaps[seq_num].begin(); pos != gaps[seq_num].end(); ++pos) {
	//cerr << " " << pos->first << " " << pos->second << ";";
	if (ali_pos < pos->first) {
	    cerr << "Add gap: " << ali_pos << " " << length << endl;
	    gaps[seq_num].emplace(pos,ali_pos,length);
	    return;
	}
	else {
	    if (ali_pos < pos->first + pos->second) {
		pos->second += length;
		cerr << endl << "Add to existing gap: " << pos->first << " " << pos->second << endl;
		return;
	    }
	}
    }
    cerr << "Add gap at end: " << ali_pos << " " << length << endl;
    gaps[seq_num].emplace_back(ali_pos,length);
}

/*unsigned int alignment::get_seq_pos ( const unsigned int seq_num, const unsigned int ali_pos) {
    if (ali_pos < start_pos[seq_num])
	return seqs[seq_num]->get_length();
    unsigned int seq_pos = ali_pos - start_pos[seq_num];
    struct gap_stack* present = gap[seq_num];
    while (present != 0) {
        if (present->pos > seq_pos)
            break;
        else if (present->pos+present->length > seq_pos) {
            return seqs[seq_num].length();
	}
        else
            seq_pos -= present->length;
	present = present->next;
    }
    if (revcomp[seq_num]) {
	if (seq_pos > seqs[seq_num].length() || seq_pos > seqs[seq_num].length_left_trimmed() || seq_pos < seqs[seq_num].n_trimmed_right())
	    seq_pos = seq[seq_num].length();
	else    
	    seq_pos = seqs[seq_num].length()-1-seq_pos;
    }
    else if (seq_pos > seqs[seq_num].length() || seq_pos < seqs[seq_num].n_trimed_left() || seq_pos > seqs[seq_num].length_right_trimmed())
        seq_pos = seq[seq_num].length();
    return seq_pos;
}*/
/*
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
*/

pair<char,int> alignment::get_base_call ( const unsigned int seq_num, const unsigned int ali_pos ) {
    pair<char,int> base(0x00,0);
    if (seq_num > n_seq)
        return base;
    if (ali_pos < start_pos[seq_num]) {
	base.first = '-';
	return base;
    }
    unsigned int seq_pos = ali_pos - start_pos[seq_num];
    for (vector<pair<unsigned int, unsigned int>>::const_iterator i=gaps[seq_num].begin(); i != gaps[seq_num].end(); ++ i) {
        if (i->first > seq_pos)
            break;
        else if (i->first+i->second > seq_pos) {
	    base.first = '-';
	    base.second = (seqs[seq_num]->get_qual(seq_pos)+seqs[seq_num]->get_qual(seq_pos+1))/2;
            return base;
	}
        else
            seq_pos -= i->second;
    }
    if (revcomp[seq_num]) {
	if (seq_pos > seqs[seq_num]->length() || seq_pos > seqs[seq_num]->length_left_trimmed() || seq_pos < seqs[seq_num]->n_trimmed_right()) {
	    base.first = '-';
	    return base;
	}
	else    
	    seq_pos = seqs[seq_num]->length()-1-seq_pos;
    }
    else if (seq_pos > seqs[seq_num]->length() || seq_pos < seqs[seq_num]->n_trimmed_left() || seq_pos > seqs[seq_num]->length_right_trimmed()) {
        base.first = '-';
	return base;

    }
    base.first = seqs[seq_num]->get_base(seq_pos);
    base.second = seqs[seq_num]->get_qual(seq_pos);
    return base;
}
/*
int get_gap_score ( unsigned int seq_num, unsigned int align_pos) {
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
*/

pair<char,int> alignment::get_con_base ( const unsigned int ali_pos ) {
    pair<char,int> base;
    if (has_consensus())
	base = get_base_call(n_seq,ali_pos);
    else {
	int G=0, C=0, T=0, A=0, gap=0;
	for (unsigned int i=0; i<n_seq; ++i) {
	    base = get_base_call(i, ali_pos);
	    if (revcomp[i])
		base.first = comp_base(base.first);
	    if (base.first=='G' || base.first=='g')
		G+=base.second;
	    else if (base.first=='C' || base.first=='c')
		C+=base.second;
	    else if (base.first=='T' || base.first=='t')
		T+=base.second;
	    else if (base.first=='A' || base.first=='a')
		A+=base.second;
	    else if (base.first=='-')
	     	gap += base.second;
	}
	if (gap >= G && gap >= C && gap >= T && gap >= A) {
	    base.first = '-';
	    base.second = gap-G-T-A-C;
	}
	else if (G >= C && G >=T && G >= A && G>gap) {
	    base.first = 'G';
	    base.second = G-T-A-C-gap;
	}
	else if (C >=T && C >= A && C > G && C > gap) {
	    base.first = 'C';
	    base.second = C-T-G-A-gap;
	}
	else if (T >=A && T > C && T > A && T > gap) {
	    base.first = 'T';
	    base.second = T-A-C-G-gap;
	}
	else {
	    base.first = 'A';
	    base.second = A-G-T-C-gap;
	}
    }
    return base;
}

unsigned int alignment::alignment_length ( ) {
    unsigned int max_length = 0;
    for (unsigned int i=0; i < n_seq; ++i) {
        if (seqs[i] != 0) {
            unsigned int length = get_seq_length(i)+start_pos[i];
            if (length > max_length)
                max_length = length;
        }
    }
    return max_length;
}

char alignment::comp_base ( char base ) {
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

void score_matrix::set_max_score ( int score, int x, int y ) {
    if (score > max_score) {
	max_score = score;
	end.first = x;
	end.second = y;
    }
}

void alignment::reverse_complement (unsigned int first_seq, unsigned int until_seq) {
    unsigned int ali_length = alignment_length();
    for (unsigned int i = first_seq; i < until_seq; ++i) {
	start_pos[i] = ali_length - start_pos[i] - seqs[i]->length();
	revcomp[i] = !revcomp[i];
	reverse(gaps[i].begin(),gaps[i].end());
	for (vector<pair<unsigned int, unsigned int>>::iterator j = gaps[i].begin(); j != gaps[i].end(); ++j) {
	    j->first = seqs[i]->length()-1-j->first;
	}
    }
}

unsigned int alignment::get_first_call () {
    unsigned int first = alignment_length();
    for (unsigned int i=0; i < n_seq; ++i) {
	unsigned int call(0);
	if (revcomp[i]) {
	    call = start_pos[i] + seqs[i]->n_trimmed_right();
	}
	else {
	    call = start_pos[i]+seqs[i]->n_trimmed_left(); 
	}
	if (call < first) first = call;
    }
    return first;
}

unsigned int alignment::get_last_call ( ) {
    unsigned int last = 0;
    for (unsigned int i=0; i < n_seq; ++i) {
	unsigned int call(0);
	if (revcomp[i]) {
	    call = start_pos[i]+get_seq_length(i)-seqs[i]->n_trimmed_left(); 
	}
	else {
		call = start_pos[i]+get_seq_length(i)-seqs[i]->n_trimmed_right(); 
	}
	if (call > last)
	    last = call;
    }
    return last;
}

void alignment::add_consensus() {
    if (n_seq == 1)
	seqs[1]=seqs[0];
    else if (n_seq > 1) {
	vector<char> seq;
	vector<int> probs;
	unsigned int first=get_first_call(), last=get_last_call();
	start_pos[n_seq] = first;
	for (unsigned int i = first; i < last; ++i) {
	    pair<char,int> base = get_con_base(i);
	    if (base.first == '-' || base.first == 0x00) {
		if (seq.empty())
		    ++start_pos[n_seq];
		else {
		    cerr << "Add gap to concensus " << seq.size() << " " << n_seq << endl;
		    add_gap(n_seq,seq.size(),1);
		}
	    }
	    else {
		seq.push_back(base.first);
		probs.push_back(base.second);
	    }
	}
	if (seqs[n_seq] != 0 && seqs[n_seq] != seqs[0]) {
	    delete seqs[n_seq];
	    seqs[n_seq] = 0;
	}
	seqs[n_seq] = new sequence("con", seq, probs);
    }
}

void alignment::align_pair (alignment& ali_one, alignment& ali_two) {
    n_seq = ali_one.n_sequences();
    seqs.resize(n_seq+1,0);
    revcomp.resize(n_seq+1);
    start_pos.resize(n_seq+1);
    gaps.resize(n_seq+1);
    for (unsigned int i=0; i < n_seq; ++i) {
	seqs[i] = ali_one.get_sequence(i);
	revcomp[i] = ali_one.get_revcomp(i);
	start_pos[i] = ali_one.get_start_pos(i);
	gaps[i] = ali_one.get_gaps(i);
    }
    add_alignment(ali_two);
}

void alignment::add_alignment ( alignment& ali ) {
    //cerr << "In pair_align" << endl;
    //cerr << "Lengths " << ali_one.alignment_length() << " " << ali_two.alignment_length() << endl;
    const unsigned int first_one=get_first_call(), first_two = ali.get_first_call(), last_one=get_last_call(), last_two=ali.get_last_call();
    score_matrix for_score(last_one, last_two, first_one, first_two);
    score_matrix rev_score(last_one, last_two, first_one, first_two);
    if (!has_consensus())
	add_consensus();
    if (!ali.has_consensus())
	ali.add_consensus();
    unsigned int length_two = ali.alignment_length();
    //cerr << "First base: " << first_one << " " << first_two << " Last: " << last_one << " " << last_two << endl;
    for (int i=first_one; i < last_one; ++i) { // loop over alignment pos in "seq1"
	//cerr << "Seq 1 base: " << i << endl;
	pair<char,int> one_base = get_con_base(i);
	// save in consensus?
	for (int j=first_two; j < last_two; ++j) { // loop over alignment pos in "seq2"
	    //cerr << "  Seq 2 base " << j << endl;
	    int match = for_score.get_score(i-1, j-1);
	    int gap_j = for_score.get_score(i,j-1);
	    int gap_i = for_score.get_score(i-1,j);
	    int rev_match = rev_score.get_score(i-1,j-1);
	    int rev_gap_j = rev_score.get_score(i,j-1);
	    int rev_gap_i = rev_score.get_score(i-1,j);
	    //cerr << "Pair align check point 1" << endl;
	    pair<char,int> two_base = ali.get_con_base(j);
	    pair<char,int> two_base_rev = ali.get_con_base(length_two-1-j);
	    two_base_rev.first = comp_base(two_base_rev.first);
	    //cerr << "Pair align check point 2" << endl;
	    // Calc prob if aligning agains each other
	    if (one_base.first == two_base.first || one_base.first == '-' || two_base.first=='-') {
	     	match += one_base.second + two_base.second;
	    }
     	    else {
 		match -= (one_base.second + two_base.second);
	    }
	    if (one_base.first == two_base_rev.first || one_base.first == '-' || two_base_rev.first=='-') {
		rev_match += one_base.second + two_base_rev.second;
	    }
	    else {
		rev_match -= (one_base.second + two_base_rev.second);
	    }
	    //cerr << "Pair align check point 3" << endl;
	    // Calc prob if inserting gaps
	    gap_j -= 2*one_base.second;
	    gap_i -= 2*two_base.second;
	    rev_gap_j -= 2*one_base.second;
	    rev_gap_i -= 2*two_base_rev.second;
	    // Check which option is best and set score
	    //cerr << "Pair align check point 4" << endl;
	    if (match >= gap_j && match >= gap_i) {
		for_score.set_score(i,j, match);
		for_score.set_prev_one(i,j,i-1);
		for_score.set_prev_two(i,j,j-1);
		for_score.set_max_score(match, i, j);
	    }
	    else if (gap_j > match && gap_j >= gap_i) {
		for_score.set_score(i,j,gap_j);
		for_score.set_prev_one(i,j,i);
		for_score.set_prev_two(i,j,j-1);
	    }
	    else if (gap_i > match && gap_i > gap_j) {
		for_score.set_score(i,j,gap_j);
		for_score.set_prev_one(i,j,i-1);
		for_score.set_prev_two(i,j,j);
	    }
	    else {
		for_score.set_score(i,j,0);
		for_score.set_prev_one(i,j,-1);
		for_score.set_prev_two(i,j,-1);
	    }
	    //cerr << "Pair align check point 5" << endl;
	    if (rev_match >= rev_gap_j && rev_match >= rev_gap_i) {
		rev_score.set_score(i,j,rev_match);
		rev_score.set_prev_one(i,j,i-1);
		rev_score.set_prev_two(i,j,j-1);
		rev_score.set_max_score(rev_match, i, j);
	    }
	    else if (rev_gap_j > rev_match && rev_gap_j >= rev_gap_i) {
		rev_score.set_score(i,j,rev_gap_j);
		rev_score.set_prev_one(i,j,i);
		rev_score.set_prev_two(i,j,j-1);
	    }
	    else if (rev_gap_i > rev_match && rev_gap_i > rev_gap_j) {
		rev_score.set_score(i,j,rev_gap_i);
		rev_score.set_prev_one(i,j,i-1);
		rev_score.set_prev_two(i,j,j);
	    }
	    else {
		rev_score.set_score(i,j,0);
		rev_score.set_prev_one(i,j,-1);
		rev_score.set_prev_two(i,j,-1);
	    }
	    //cerr << "Pair align check point 6 prev 1: " << for_score.get_prev_one(i-first_one,j-first_two) << " (" << i << ") prev 2: " << for_score.get_prev_two(i-first_one,j-first_two) << " (" << j <<")" << endl;
	    //cerr << "Pair align check point 6 (rev) prev 1: " << rev_score.get_prev_one(i-first_one,j-first_two) << " (" << i << ") prev 2: " << rev_score.get_prev_two(i-first_one,j-first_two) << " (" << j <<")" << endl;
	}
    }
    //cerr << "Pair align check point 7: " << rev_score.get_prev_one(949-first_one,731-first_two) << " " << rev_score.get_prev_two(949-first_one,730-first_two) << endl;
    struct score_matrix* best_score = 0;
    bool revcomp_two(false);
    // get pos in seqs for best score
    if (for_score.get_max_score() >= rev_score.get_max_score())
	best_score = &for_score;
    else {
	best_score = &rev_score;
	revcomp_two = true;
    }
    // initialize new alignment
    //cerr << "Pair align check point 8" << endl;
    //if (revcomp_two)
	//cerr << "Revcomp two" << endl;
    unsigned int n_seq_o = n_seq;
    if (seqs[n_seq] != 0 && seqs[n_seq] != seqs[0]) // clear any consensus seq
	delete seqs[n_seq];
    n_seq = n_sequences()+ali.n_sequences();
    seqs.resize(n_seq+1);
    for (unsigned int i=0; i < ali.n_sequences(); ++i) {
	seqs[i+n_seq_o] = ali.get_sequence(i);
    }
    revcomp.resize(n_seq+1);
    for (unsigned int i=0; i < ali.n_sequences(); ++i) {
        revcomp[i+n_seq_o] = ali.get_revcomp(i);
    }
    //cerr << "Pair align check point 9" << endl;
    start_pos.resize(n_seq+1);
    gaps.resize(n_seq+1);
    for (unsigned int i=0; i < ali.n_sequences(); ++i) {
        gaps[i+n_seq_o] = ali.get_gaps(i);
    }
    score = best_score->get_max_score();
    //cerr << "Pair align check point 10: " << score << endl;
    // 
    int i = best_score->get_end_one();
    int j = best_score->get_end_two();
    while ( 1 ) {
	const int next_i(best_score->get_prev_one(i,j));
	const int next_j(best_score->get_prev_two(i,j));
	//cerr << "i: " << i <<  " " << ali_one.get_con_base(i).first << " (" << next_i << " [" << rev_score.get_prev_one(i,j) << "]) j: " << j << " " << comp_base(ali_two.get_con_base(length_two-1-j).first) << " (" << next_j << " [" << rev_score.get_prev_two(i,j) << "])" << " - " << best_score->get_score(i, j) << endl;
	if (next_i == i && next_j == j) {
	    //cerr << "Error adding gap to both seqs - i:  " << i << " j: " << j << endl;
	    break;
	}
	else if (next_i == i) {
	    add_gap(i, 1, 0, n_seq_o);
	    cerr << "Add gap to Ali One pos: " << i << endl;
	}
	else if (next_j == j) {
	    cerr << "Add gap to Ali Two pos: " << j << endl;
	    /*if (revcomp_two) {
		add_gap(length_two-1-j, 1, ali_one.n_sequences(), n_seq);
	    }
	    else {*/
		add_gap(j, 1, n_seq_o, n_seq);
	    //}
	}
	if (next_i < 0 || next_j < 0)
	    break;
	i = next_i;
	j = next_j;
    }
    //cerr << "Pair align check point 11. Length one:" << get_seq_length(0) << " (" << length_gaps(0) << ") [" << start_pos[0] << "] length two " << get_seq_length(1) << " (" << length_gaps(1) << ") [" << start_pos[1] << "]" << endl;
    if (revcomp_two) {
	unsigned int length = 0;
	for (unsigned int i=n_seq_o; i < n_seq; ++i) {
	    unsigned int seq_length = get_seq_length(i)+start_pos[i];
	    if (seq_length > length) length=seq_length;
	}
	for (unsigned int i=n_seq_o; i < n_seq; ++i) {
	    revcomp[i] = !revcomp[i];
	    //cerr << "Seq " << i << " is revcomp: " << revcomp[i] << endl;
	    start_pos[i] = length-get_seq_length(i)-start_pos[i];
	}
    }
    //cerr << "Pair align check point 11. Length one:" << get_seq_length(0) << " (" << length_gaps(0) << ") [" << start_pos[0] << "] length two " << get_seq_length(1) << " (" << length_gaps(1) << ") [" << start_pos[1] << "]" << endl;
    if (j>i)
	move_start_pos(j-i,0, n_seq_o);
    else if (i>j) {
	move_start_pos(i-j, n_seq_o, n_seq);
    }	    
    //cerr << "Pair align check point 11. Length one:" << get_seq_length(0) << " (" << length_gaps(0) << ") [" << start_pos[0] << "] length two " << get_seq_length(1) << " (" << length_gaps(1) << ") [" << start_pos[1] << "]" << endl;
    //cerr << "N gaps: " << gaps[0].size() << " " << gaps[1].size() << endl;
}

char alignment::get_seq_base ( unsigned int seq_num, unsigned int seq_pos ) {
    if (revcomp[seq_num])
	return comp_base(seqs[seq_num]->get_base(seqs[seq_num]->length()-1-seq_pos));
    else
	return seqs[seq_num]->get_base(seq_pos);
}

void alignment::print_seq(ostream& out, const unsigned int seq, unsigned int first_ali_pos, unsigned int last_ali_pos, const char gap_char) {
    if (seq > n_seq || last_ali_pos < start_pos[seq] || first_ali_pos > get_seq_length(seq)+start_pos[seq]) {
	//cerr << "Out of bound" << endl;
	last_ali_pos -= first_ali_pos;
	while (last_ali_pos) {
	    out << gap_char;
	    --last_ali_pos;
	}
    }
    else {
	for (unsigned int i=first_ali_pos; i < start_pos[seq]; ++i) {
    	    out << gap_char;
	    //cerr << "Print starting gap" << endl;
	}
	if (first_ali_pos > start_pos[seq])
	    first_ali_pos -= start_pos[seq];
	else
	    first_ali_pos = 0;
	last_ali_pos -= start_pos[seq];
	// Add to trim low qual bases
	for (vector<pair<unsigned int, unsigned int>>::iterator gap = gaps[seq].begin(); gap != gaps[seq].end(); ++gap) {
	    if (gap->first < last_ali_pos) {
		//cerr << "Print until and gap at pos " << gap->first << endl;;
		for (;first_ali_pos <= gap->first; ++first_ali_pos) {
		    out << seqs[seq]->get_base(first_ali_pos, revcomp[seq]);
		}
		unsigned int n_gaps = last_ali_pos - gap->first;
		if (n_gaps > gap->second) {
		    n_gaps = gap->second;
		    last_ali_pos -= gap->second;
		}
		else
		    last_ali_pos = 0;
		for (unsigned int i=0; i < n_gaps; ++i)
		    out << gap_char;
	    }
	    else
		break;
	}
	unsigned int last_base = seqs[seq]->length();
	// Add to trim low qual bases
	if (last_ali_pos < last_base)
	    last_base = last_ali_pos;
	for (; first_ali_pos < last_base; ++first_ali_pos)
	    out << seqs[seq]->get_base(first_ali_pos,revcomp[seq]);
	for (; first_ali_pos < last_ali_pos; ++first_ali_pos)
	    out << gap_char;
    }	    
}

void alignment::print_alignment_fasta ( ostream& out ) {
    unsigned int ali_length = alignment_length();
    for (unsigned int i = 0; i < n_seq; ++i) {
	out << '>' << seqs[i]->get_name() <<endl;
	print_seq(out, i, 0, ali_length, '-');
	out << endl;
    }
}

void alignment::print_alignment_json ( ostream& out ) {
    if (!seqs[n_seq]) {
	add_consensus();
    }
    out << "{\"name\":\"" << name << "\",\"length\":" << get_seq_length(n_seq) << ", ";
    out << "\"n_reads\":" << n_seq << ",\"reverse\":false,";
    out << "\"seq\":[";
    for (unsigned int i=0; i < seqs[n_seq]->length(); ++i) {
	if (i != 0)
	    out << ',';
	out << '"' << seqs[n_seq]->get_base(i) << '"';
    }
    out << "],\"qual\":[ ";
    for (unsigned int i=0; i < seqs[n_seq]->length(); ++i) {
	if (i != 0)
	    out << ',';
	out << seqs[n_seq]->get_qual(i);
    }
    out << "],\"reads\":[ ";
    for (unsigned int i=0; i < n_seq; ++i) {
	if (i != 0)
	    out << ", ";
	seqs[i]->print_json(out, revcomp[i], "");
    }
    out << " ], \"alignment\":{\"contig\":{ \"start\":" << start_pos[n_seq] << ", \"gaps\":[";
    for (vector<pair<unsigned int, unsigned int>>::const_iterator i=gaps[n_seq].begin(); i != gaps[n_seq].end(); ++i) {
	if (i != gaps[n_seq].begin())
	    out << ',';
	out << "{\"pos\":" << i->first << ", \"qual\":[";
	for (unsigned j=0; j < i->second; ++j) {
	    if (j !=0)
		out << ',';
	    out << get_con_base(i->first+1+j).second;
	}
	out << "], \"length\":" << i->second << "}";
    }
    out << "]},\"reads\":[";
    for (unsigned int i=0; i < n_seq; ++i) {
	if (i != 0)
	    out << ',';
	out << "{ \"start\":" << start_pos[i] << ", \"gaps\":[";
	for (vector<pair<unsigned int, unsigned int>>::const_iterator j=gaps[i].begin(); j != gaps[i].end(); ++j) {
	    if (j != gaps[i].begin())
		out << ',';
	    out << "{\"pos\":" << j->first << ",\"length\":" << j->second << "}";
	}
	out << "]}";
    }
    out << "],\"length\":" << alignment_length() << "}}";

}
/*contig alignment_to_contig (char* name, struct alignment* ali, Read** read_data ) {
    fputs("Convert into contig\n", stderr);
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
    fputs("Getting contig\n", stderr);
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
	data.reads[i].seq_length = get_seq_length(ali,i); //
	data.reads[i].seq = (char*)malloc(sizeof(char)*data.reads[i].seq_length); // ali->seq.base[i]; //
	unsigned int pos = 0;
	//struct gap_stack* present = ali->gap[i];
	fprintf(stderr, "Read of length: %u", data.reads[i].seq_length);
	for (unsigned int j=ali->startPos[i]; j < ali->length && pos < data.reads[i].seq_length; ++j) {
	    unsigned int base_pos = get_seq_pos(ali, i, j); 
	    if (base_pos >= ali->seq->Nbases[i])
		data.reads[i].seq[pos] = '*';
	    else
		data.reads[i].seq[pos] = ali->seq->base[i][base_pos];
	    ++pos;
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
}*/
