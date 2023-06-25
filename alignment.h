/********************************************************************
Copyright (C) 2022 Martin Ryberg

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

contact: martin.ryberg@ebc.uu.se
*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <utility>
#include <vector>
#include <algorithm>
#include <io_lib/Read.h>
#include "sequence.h"
#include <string>

#ifndef ALIGNMENT
#define ALIGNMENT

class score_matrix {
public:
    score_matrix (const int x, const int y, const int x_off, const int y_off) {
	n_rows = y;
	n_cols = x;
	if (y_off < y)
	    rows_offset=y_off;
	else
	    rows_offset = n_rows;
	if (x_off < x)
	    cols_offset = x_off;
	else
	    cols_offset = n_cols;
	score.resize(n_rows-rows_offset);
	for (vector<vector<int>>::iterator i=score.begin(); i != score.end(); ++i)
	    i->resize(n_cols-cols_offset);	    
	prev.resize(n_rows-rows_offset);
	for (vector<vector<pair<int,int>>>::iterator i=prev.begin(); i != prev.end(); ++i)
	    i->resize(n_cols-cols_offset);	    
	max_score = 0;
	end.first = 0;
	end.second = 0;
    }
    void set_max_score(int score, int i, int j);
    int get_max_score() { return max_score; };
    void set_score ( int x, int y, int value ) {
	if (x-cols_offset >= 0 && x < n_cols-cols_offset && y-rows_offset >= 0 && y < n_rows-rows_offset)
	    score[y-rows_offset][x-cols_offset] = value;
    }
    int get_score(int x, int y) {
	if (x-cols_offset >= 0 && x < n_cols-cols_offset && y-rows_offset >= 0 && y < n_rows-rows_offset)
	    return score[y-rows_offset][x-cols_offset];
	else
	    return 0;
    }
    int get_end_one() { return end.first; }
    int get_end_two() { return end.second; }
    void set_prev_one ( int x, int y, int value ) {
	if (x-cols_offset >= 0 && x < n_cols-cols_offset && y-rows_offset >= 0 && y < n_rows-rows_offset)
	    prev[y-rows_offset][x-cols_offset].first = value;
    }
    int get_prev_one ( int x, int y ) {
	if (x-cols_offset >= 0 && x < n_cols-cols_offset && y-rows_offset >= 0 && y < n_rows-rows_offset)
	    return prev[y-rows_offset][x-cols_offset].first;
	else
	    return -1;
    }
    void set_prev_two ( int x, int y, int value ) {
	if (x-cols_offset >= 0 && x < n_cols-cols_offset && y-rows_offset >= 0 && y < n_rows-rows_offset)
	    prev[y-rows_offset][x-cols_offset].second = value;
    }
    int get_prev_two ( int x, int y ) {
	if (x-cols_offset >= 0 && x < n_cols-cols_offset && y-rows_offset >= 0 && y < n_rows-rows_offset)
	    return prev[y-rows_offset][x-cols_offset].second;
	else
	    return -1;
    }
private:
    int n_rows;
    int n_cols;
    int rows_offset;
    int cols_offset;
    vector<vector<int>> score;
    vector<vector<pair<int,int>>> prev;
    int max_score;
    pair<int,int> end;
    
   /* unsigned int get_matrix_pos (unsigned int x, unsigned int y ) {
	if (y > n_rows || x > n_cols)
	    return n_cols * n_rows;
	return y*n_cols + x;
    }*/

};

class alignment {
/*
Positions:
There are two main position coordinates, alignment position (ali_pos) which is
the column in the alignment and sequence position (seq_pos) which is the base in
the sequence. Gaps are given by the sequence position, but in reversed sequences
this has to be counted from the back of the sequence. The start position
(start_pos) is the position in the alignment where the sequence start. The
sequence start include any quality trimed bases even if they may not be aligned.

*/
public:
    alignment (unsigned int N) {
	n_seq = N;
	seqs.resize(n_seq+1,0);
	revcomp.resize(n_seq+1);
	start_pos.resize(n_seq+1);
	gaps.resize(n_seq+1);
	score = 0;
    }
    alignment ( sequence* seq ) {
	n_seq = 1;
	seqs.resize(2,0);
	seqs[0] = seq;
	revcomp.resize(2);
	start_pos.resize(2);
	gaps.resize(2);
	score = 0;
    }
    ~alignment ( ) {
	if (seqs[n_seq] && seqs[n_seq] != seqs[0])
	    delete seqs[n_seq];
    }
    alignment ( const alignment& o ) {
	cerr << "Copy alignment" << endl;
	n_seq=o.n_seq;
	seqs = o.seqs;
	cerr << "Length seq 0: " << seqs.at(0)->length() << " (" << o.seqs.at(0)->length() << endl;
	revcomp = o.revcomp;
	start_pos = o.start_pos;
	gaps = o.gaps;
	score = score;
    }
    unsigned int n_sequences() { return n_seq;  }
    bool has_consensus() { return bool(seqs[n_seq]); }
    sequence* get_sequence( unsigned int seq ) { return seqs[seq]; }
    unsigned int get_seq_pos ( const unsigned int seq_num, const unsigned int align_pos);
    char get_seq_base ( unsigned int seq_num, unsigned int pos );
    pair<char,int> get_base_call ( const unsigned int seq_num, const unsigned int ali_pos );
    pair<char,int> get_con_base ( const unsigned int ali_pos );
    unsigned int get_first_call();
    unsigned int get_last_call();
    bool get_revcomp( unsigned int seq ) { return revcomp[seq]; }
    vector<pair<unsigned int, unsigned int>> get_gaps(unsigned int seq ) { return gaps[seq]; }
    unsigned int alignment_length ( );
    void set_name (string ali_name) { name = ali_name; }
    void add_gap ( const unsigned int seq_num, const unsigned int ali_pos, const unsigned int length );
    void add_gap ( const unsigned int ali_pos, const unsigned int length) {
	add_gap(ali_pos,length, 0, n_seq);
    }
    void add_gap ( const unsigned int ali_pos, const unsigned int length, unsigned int first_seq, unsigned int until_seq) {
	for (unsigned int i=first_seq; i < until_seq; ++i) {
	    add_gap(i, ali_pos, length);
	}
    }
    void add_consensus();
    void reverse_complement (unsigned int first_seq, unsigned int until_seq);
    void align_pair ( alignment& ali_one, alignment& ali_two );
    void print_alignment_fasta ( ostream& out );
    void print_alignment_json ( ostream& out );

private:
    /// Variables
    string name;
    vector<sequence*> seqs;
    unsigned int n_seq;
    vector<bool> revcomp;
    vector<unsigned int> start_pos;
    vector<vector<pair<unsigned int, unsigned int>>> gaps;
    int score;
/*
struct base_prob {
    int* G;
    int* C;
    int* T;
    int* A;
    int* gap;
};*/

// Data struct handling //
    unsigned int length_gaps ( int seq ) {
	unsigned int length(0);
	for (vector<pair<unsigned int, unsigned int>>::const_iterator i= gaps[seq].begin(); i != gaps[seq].end(); ++i)
	    length += i->second;
	return length;
    }

    unsigned int get_seq_length (unsigned int seq_num) {
	return seqs[seq_num]->length() + length_gaps(seq_num);
    }

    char comp_base ( char base );

    void print_seq(ostream& out, const unsigned int seq, unsigned int first_ali_pos, unsigned int last_ali_pos, const char gap_char);
    void move_start_pos( const unsigned int n_pos, const unsigned int first_seq, const unsigned int until_seq) {
	for (unsigned int i=first_seq; i < until_seq; ++i)
	    start_pos[i] += n_pos;
    }
    void rev_seqs(const unsigned int first_seq, const unsigned int until_seq) {
	unsigned int length = alignment_length();
        for (unsigned int i=first_seq; i < until_seq; ++i) {
	    revcomp[i] = !revcomp[i];
	    start_pos[i] = length-get_seq_length(i)-start_pos[i];
	}
    }
	    
/*void add_sequence ( struct multiple_sequences* ms, const unsigned int seq_num, char* name, const unsigned int length, char* seq, char* qual_G, char* qual_C, char* qual_T, char* qual_A, const unsigned int trim_from_start, const unsigned int trim_from_end );
void add_read ( struct multiple_sequences* ms, Read* read, const unsigned int seq_num);
*/

/*void insert_gap(struct alignment* ali, const unsigned int ali_pos, const unsigned int gap_length);
void move_start_pos (struct alignment* ali, const int n);
void set_max_score ( struct score_matrix* matrix, int score, unsigned int i, unsigned int j);
*/
//struct multiple_sequences* join_sequences ( struct multiple_sequences* seq1, struct multiple_sequences* seq2 );
//struct alignment* merged_alignment(struct alignment* ali_one, struct alignment* ali_two);
//void revcomp (struct alignment* ali);

// Get infor from data structs //
/*char get_call_prob ( struct alignment* ali, const unsigned int seq_num, const unsigned int seq_pos );
struct base_call get_base_call ( struct alignment* ali, const unsigned int seq_num, const unsigned int align_pos );
int get_gap_score ( struct alignment* ali, unsigned int seq_num, unsigned int align_pos);
unsigned int get_first_call (struct alignment* ali);
unsigned int get_last_call (struct alignment* ali);
unsigned int get_seq_length (struct alignment* ali, unsigned int seq_num);
*/

// Alignment

// Pring and conversion to other struct complex
//void fprint_fasta_alignment (struct alignment* ali, FILE* out );
//contig alignment_to_contig (char* name, struct alignment* ali, Read** read_data );
};
#endif //ALIGNMENT
