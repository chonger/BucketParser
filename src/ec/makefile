.C.o:
	g++ $(CFLAGS) -c -O5 $<

PARSERTEST_OBJS = \
	ECArgs.o \
	Tree.o \
	Rule.o \
	Unkify.o \
	Trellis.o \
	CellClos.o \
	Parser.o \
	parser.o 
#	parserTest.o



ecall: $(PARSERTEST_OBJS)

ecclean:
	rm -f *.o  *~ parser

parserTest: parserTest.o $(PARSERTEST_OBJS)
	g++ $(CFLAGS) -g $(PARSERTEST_OBJS) parserTest.o -o parserTest -D_REENTRANT -D_XOPEN_SOURCE=600 -lpthread