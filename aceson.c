#include <stdio.h>
#include <io_lib/Read.h>
#include "contig_handler.h"

int main ( int argc, char* argv[]) {
    FILE* f = fopen(argv[1],"r");
    struct contig_node* contigs = parsACE(f, NULL);
    add_traces(contigs);
    //fprintf(stderr, "N contigs %u | n reads %u\n", n_contigs(contigs), n_reads(contigs));
    //printACE(stdout, contigs);
    printJSON(stdout, contigs);
    delete_contigs(contigs);
}
