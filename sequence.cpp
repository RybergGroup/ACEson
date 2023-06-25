#include "sequence.h"

void sequence::add_read (Read* seq_read) {
    if (name)
	delete[] name;
    name = new char[strlen(seq_read->trace_name)+1];
    strcpy(name,seq_read->trace_name);
    read = seq_read;
    n_bases = seq_read->NBases;
    if (seq)
	delete[] seq;
    seq = seq_read->base;
    if (prob)
	delete[] prob;
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
    qual_start = seq_read->leftCutoff;
    if (seq_read->rightCutoff != 0)
	qual_end = seq_read->rightCutoff;
    else
	qual_end = n_bases;
}

void sequence::print_read_json(ostream& out) {
    if (read != NULL) {
	out << "{\"file_name\":\"" << read->trace_name << "\", \"n_bases\":" << read->NBases << ", \"left_cut_off\":" << read->leftCutoff <<", \"right_cut_off\":" << read->rightCutoff << ", \"traces\":{\"n_points\":" << read->NPoints << ", \"max_val\":" << read->maxTraceVal << ", \"base_line\":" << read->baseline << ",\"A\":[";
	for (int i=0; i < read->NPoints; ++i) {
	    if (i!=0) out << ',';
	    out << read->traceA[i];
	}
	out << "],\"C\":[";
	for (int i=0; i < read->NPoints; ++i) {
	    if (i!=0) out << ',';
	    out << read->traceC[i];
	}
	out << "],\"G\":[";
	for (int i=0; i < read->NPoints; ++i) {
	    if (i!=0) out << ',';
	    out << read->traceG[i];
	}
	out << "],\"T\":[";
	for (int i=0; i < read->NPoints; ++i) {
	    if (i!=0) out << ',';
	    out << read->traceT[i];
	}
	out << "]}";
	out << ", \"calls\":{ \"seq\":[";
	for (int i=0; i < read->NBases; ++i) {
	    if (i!=0) out << ',';
	    out << '"' << read->base[i] << '"';
	}
	out << "], \"trace_pos\":[";
	for (int i=0; i < read->NBases; ++i) {
	    if (i!=0) out << ',';
	    out << read->basePos[i];
	}
	out << "], \"base_prob\":[";
	for (int i=0; i < read->NBases; ++i) {
	    if (i!=0) out << ',';
	    out << "{\"A\":" << (int)read->prob_A[i] << ",\"C\":" << (int)read->prob_C[i] << ",\"G\":" << (int)read->prob_G[i] << ",\"T\":" << (int)read->prob_T[i] << "}";
	}
	out << "]}";
	out << '}';
    }
}

void sequence::print_json(ostream& out, bool revcomp, const char* description) {
    out <<"{\"name\":\"" << name << "\", \"length\":" << n_bases << ", \"start\":" << qual_start << ",\"reverse\":";
    if (revcomp)
	out << "true";
    else
	out << "false";
    out << ", \"qual_begin\":" << qual_start << ", \"qual_end\":" << qual_end << ", \"align_begin\":" << qual_start << ", \"align_end\":" << qual_end;
    if (read !=0) {
	out << ", \"chrom\":";
	print_read_json( out );
    }
    else {
	out << ", \"seq\":[";
	for (unsigned int i=0; i < n_bases; ++i) {
	    if (i != 0)
		out << ',';
	    out << '"' << seq[i] << '"';
	}
	out << ']';
    }
    out << ", \"description\":\"" << description << "\"}";
}

