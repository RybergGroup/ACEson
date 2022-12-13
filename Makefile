CC:=gcc
ACESON:=aceson.c 
CONTIG:=contig_handler.c
ALIGN:=alignment.c
IOLIB:=staden-read

aceson: $(ACESON) $(CONTIG)
	$(CC) $(ACESON) $(CONTIG) $(ALIGN) -Wall -l$(IOLIB) -o $@
