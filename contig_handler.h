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
#include <stdbool.h>
#include <string.h>
#include <unistd.h>

#ifndef CONTIGHEADER
#define CONTIGHEADER

#define WORDSIZE 1024

typedef struct {
    char* name;
    char* seq;
    unsigned int seq_length;
    int start;
    unsigned int qual_begin;
    unsigned int qual_end;
    unsigned int align_begin;
    unsigned int align_end;
    char* description;
    Read* read_data;
    _Bool revcomp;
} contig_read;

typedef struct {
    char* name;
    char* contig_seq;
    unsigned int* contig_qual;
    unsigned int contig_length;
    unsigned int n_gaps;
    char* segment;
    contig_read* reads;
    unsigned int n_reads;
    _Bool revcomp;
} contig;

struct contig_node {
    contig data;
    struct contig_node* next;
};

unsigned int n_contigs (struct contig_node* start);

unsigned int n_reads ( struct contig_node* start );

struct contig_node* add_contig(struct contig_node* parent);

void add_traces( struct contig_node* contigs );

void delete_contigs ( struct contig_node* node );

_Bool isWhitespace ( char c );

_Bool isNewlineChar ( char c );

struct contig_node* parsACE ( FILE* in, struct contig_node* contigs );

void printACE ( FILE* out, struct contig_node* contigs );

void printJSON ( FILE* out, struct contig_node* contigs );

void print_readJSON ( FILE* out, contig_read* read );

void print_seq ( FILE* out, unsigned int length, char* seq);

#endif //CONTIGHEADER
