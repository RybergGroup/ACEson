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
#include <cstring>
#include <limits.h>
#include <vector>
#include <utility>
#include <io_lib/Read.h>

using namespace std;

#ifndef SEQUENCE
#define SEQUENCE
class sequence {
public:
    sequence () {
	name = 0;
	read = 0;
	seq = 0;
	prob = 0;
	n_bases = 0;
	qual_start=0;
	qual_end = 0;
    }
    sequence (const char* seq_name) {
	name = new char[strlen(name)+1];
	strcpy(name,seq_name);
	read = 0;
	seq = 0;
	prob = 0;
	n_bases = 0;
	qual_start=0;
	qual_end = 0;
    }
    sequence( const char* seq_name, unsigned int length, bool qual) {
	name = new char[strlen(name)+1];
	strcpy(name,seq_name);
	read = 0;
	seq = new char[length];
	if (qual)
	    prob = new int[length];
	else
	    prob = 0;
	n_bases = length;
	qual_start = 0;
	qual_end = length;
    }
    sequence ( const char* seq_name, vector<char>& sequence, vector<int>& probs) {
	name = new char[strlen(seq_name)+1];
	strcpy(name,seq_name);
	read = 0;
	seq = new char[int(sequence.size())];
	copy(sequence.begin(),sequence.end(),seq);
	if (!probs.empty() && probs.size() == sequence.size()) {
	    prob = new int[int(sequence.size())];
	    copy(probs.begin(),probs.end(),prob);
	}
	else
	    prob = 0;
	n_bases = sequence.size();
	qual_start = 0;
	qual_end = n_bases;
    }
    sequence (Read* seq_read) {
	cerr << "Read: " << seq_read << endl;
	name = new char[strlen(seq_read->trace_name)+1];
	strcpy(name,seq_read->trace_name);
	read = seq_read;
	n_bases = seq_read->NBases;
	seq = seq_read->base;
	prob = new int[n_bases];
	for (int i=0; i < seq_read->NBases; ++i) {
	    if (seq[i] == 'G' || seq[i] == 'g')
		prob[i] = (int)seq_read->prob_G[i];
	    else if (seq[i] == 'C' || seq[i] == 'c')
		prob[i] = (int)seq_read->prob_C[i];
	    else if (seq[i] == 'T' || seq[i] == 't')
		prob[i] = (int)seq_read->prob_T[i];
	    else if (seq[i] == 'A' || seq[i] == 'a')
		prob[i] = (int)seq_read->prob_A[i];
	}
	//revcomp = 0;
	cerr << "Cut offs " << seq_read->leftCutoff << " " << seq_read->rightCutoff << endl;
	qual_start = seq_read->leftCutoff;
	if (seq_read->rightCutoff != 0)
	    qual_end = seq_read->rightCutoff;
	else
	    qual_end = n_bases;
	cerr << "Set cut off " << qual_start << " " << qual_end << endl;
    }
    sequence (const sequence& o) {
	cerr << "Copy seq" << endl;
	name = new char[strlen(o.name)+1];
	strcpy(name, o.name);
	cerr << "copied name " << o.name << " " << name << endl;
	read = o.read;
	cerr << "copied read " << o.read << endl;
	n_bases = o.n_bases;
	cerr << "copied n bases " << o.n_bases << " " << n_bases << endl;
	if (o.seq && o.read && o.seq != o.read->base) {
	    cerr << "Copy sequence" << endl;
	    seq = new char[n_bases];
	    for (unsigned int i=0; i < n_bases; ++i) {
		cerr << o.seq[i] << " to " << seq[i] << endl;
		seq[i] = o.seq[i];
	    }
	}
	else
	    seq = o.seq;
	cerr << "Copied seq but not prob" << endl;
	if (o.prob != 0) {
	    prob = new int[n_bases];
	    for (unsigned int i=0; i < n_bases; ++i)
		prob[i] = o.prob[i];
	}
	else
	    prob = 0;
	qual_start = o.qual_start;
	qual_end = o.qual_end;
	cerr << "Set cut off " << qual_start << " " << qual_end << endl;
	cerr << "Copied seq" << endl;
    }
    ~sequence() {
	delete[] name;
	if (seq != 0 && read != 0 && seq != read->base)
	    delete[] seq;
	if (prob != 0)
	    delete[] prob;
    }
    void add_read (Read* seq_read);
    unsigned int length() { return n_bases; }
    unsigned int length_trimmed () { return qual_end - qual_start; }
    unsigned int length_left_trimmed (){ return n_bases - qual_start; }
    unsigned int length_right_trimmed () { return qual_end; }
    unsigned int n_trimmed_right () { return n_bases-qual_end; }
    unsigned int n_trimmed_left () { return qual_start; }
    char get_base( unsigned int n ) {
	if (n < n_bases)
	    return seq[n];
	else
	    return 0x00;
    }
    char get_revcomp_base (unsigned int n) {
	if (n < n_bases)
	    return comp_base(seq[n_bases-1-n]);
	else
	    return 0x00;
    }
    char get_base(unsigned int n, bool revcomp) {
	if (revcomp)
	    return get_revcomp_base(n);
	else
	    return get_base(n);
    }
    int get_qual( unsigned int n ) {
	if (prob && n < n_bases)
	    return prob[n];
	else
	    return 0;
    }
    Read* get_read() { return read; }
    const char* get_name () { return name; }
    void shorten( ) { --n_bases; }
    bool shorten( unsigned int n ) {
	if (n <= n_bases) {
	    n_bases -= n;;
	    return true;
	}
	else {
	    return false;
	}
    }
    void print_read_json(ostream& out);
    void print_json(ostream& out, bool revcomp, const char* description);
private:
    char* name;
    Read* read;
    unsigned int n_bases;
    char* seq;
    int* prob;
    unsigned int qual_start;
    unsigned int qual_end;
    // bool revcomp;
// data management
    sequence operator= ( const sequence& x);
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

};
#endif //SEQUENCE
