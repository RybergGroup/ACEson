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
#include <io_lib/Read.h>
#include "contig_handler.h"

#ifndef ALIGNMENT
#define ALIGNMENT

struct base_call {
    char base;
    int prob;
};

struct score_matrix {
    unsigned int Nrows;
    unsigned int Ncols;
    int* score;
    int* prevOne;
    int* prevTwo;
    int maxScore;
    int endI;
    int endJ;
};

struct gap_stack {
    unsigned int pos;
    unsigned int length;
    struct gap_stack* next;
};

struct base_prob {
    char* G;
    char* C;
    char* T;
    char* A;
};

struct multiple_sequences {
    char** name;
    unsigned int Nseq;
    unsigned int* Nbases;
    unsigned int* leftCutoff;
    unsigned int* rightCutoff;
    char** base;
    struct base_prob* prob;

};

struct alignment {
    struct multiple_sequences* seq;
    unsigned int Nseq;
    unsigned int length;
    unsigned int* startPos;
    _Bool* revcomp;
    struct gap_stack** gap;
};

struct multiple_sequences* multiple_sequences_alocate( const unsigned int num_seq );
void multiple_sequences_dealocate ( struct multiple_sequences* ms );

struct alignment* alignment_alocate( struct multiple_sequences* multiple_seq );
void alignment_dealocate ( struct alignment* ali );
void alignment_ms_dealocate ( struct alignment* ali);

void gap_stack_dealocate ( struct gap_stack* present );
struct gap_stack* copy_gap_stack ( struct gap_stack* present );

struct score_matrix* score_matrix_alocate ( unsigned int length_seq_one, unsigned int length_seq_two );
void score_matrix_dealocate( struct score_matrix* matrix );

void add_sequence ( struct multiple_sequences* ms, const unsigned int seq_num, char* name, const unsigned int length, char* seq, char* qual_G, char* qual_C, char* qual_T, char* qual_A, const unsigned int trim_from_start, const unsigned int trim_from_end );
void add_read ( struct multiple_sequences* ms, Read* read, const unsigned int seq_num);
void add_gap ( struct alignment* ali, const unsigned int seq_num, const unsigned int pos, const unsigned int length );
void insert_gap(struct alignment* ali, const unsigned int ali_pos, const unsigned int gap_length);
void move_start_pos (struct alignment* ali, const int n);
void set_max_score ( struct score_matrix* matrix, int score, unsigned int i, unsigned int j);

unsigned int get_matrix_pos (struct score_matrix* matrix,  unsigned int x, unsigned int y );
unsigned int get_seq_pos ( struct alignment* ali, const unsigned int seq_num, const unsigned int align_pos);
char get_call_prob ( struct alignment* ali, const unsigned int seq_num, const unsigned int seq_pos );
struct base_call get_base_call ( struct alignment* ali, const unsigned int seq_num, const unsigned int align_pos );
char get_seq_base ( struct alignment* ali, unsigned int seq_num, unsigned int pos );
int get_gap_score ( struct alignment* ali, unsigned int seq_num, unsigned int align_pos);
struct base_call get_con_base ( struct alignment* ali, const unsigned int align_pos );
unsigned int alignment_length ( struct alignment* ali );

char comp_base ( char base );
struct multiple_sequences* join_sequences ( struct multiple_sequences* seq1, struct multiple_sequences* seq2 );
struct alignment* merged_alignment(struct alignment* ali_one, struct alignment* ali_two);
void revcomp (struct alignment* ali);
struct alignment* align_pair ( struct alignment* ali_one, struct alignment* ali_two );

void fprint_fasta_alignment (struct alignment* ali, FILE* out );
contig alignment_to_contig (char* name, struct alignment* ali, Read** read_data );
//void fprint_JSON_alignment (struct alignment* ali, FILE* out );
#endif //ALIGNMENT
