CC:=gcc
ACESON:=aceson.c 
CONTIG:=contig_handler.c
IOLIB:=staden-read

aceson: $(ACESON) $(CONTIG)
	$(CC) $(ACESON) $(CONTIG) -l$(IOLIB) -o $@
