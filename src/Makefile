CC:=g++
ACESON:=aceson.cpp
ALIGN:=alignment.cpp
SEQ:=sequence.cpp
IOLIB:=staden-read

aceson: $(ACESON) $(CONTIG)
	$(CC) $(ACESON) $(SEQ) $(ALIGN) -Wall -l$(IOLIB) -o $@
