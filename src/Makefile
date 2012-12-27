PARSERTEST_OBJS = \
	ECArgs.o \
	Tree.o \
	Rule.o \
	Unkify.o \
	Trellis.o \
	CellClos.o \
	Parser.o \
	parser.o 

BINDIR=../bin/
OBJDIR=../obj/
CORE=BucketParser BucketGrammar BucketEM
PROGS=parse em trim maketag
ECDIR=ec/
COREOBJS=$(addsuffix .o,$(addprefix $(OBJDIR),$(CORE))) $(addprefix $(ECDIR),$(PARSERTEST_OBJS))
LIBS=-lgsl -lgslcblas -lm -pthread -lpthread #-lboost_thread
INCLUDE=../lib/sparsehash/include

#FLAGS=-g
FLAGS=-O5
CXX=g++

$(OBJDIR)%.o: %.cpp %.hpp
	$(CXX) $(FLAGS) -I$(INCLUDE)  -w -c $< -o $@

$(OBJDIR)%.o: %.cpp $(COREOBJS)
	$(CXX) $(FLAGS) -I$(INCLUDE)  -w -c $< -o $@

bsall: parse em trim maketag

clean:
	rm -f $(OBJDIR)*
	rm -f $(BINDIR)*
	make -C ec ecclean

parse: $(OBJDIR)BucketParserMain.o $(COREOBJS)
	$(CXX) $(FLAGS)  $(LIBS) -w -o $(BINDIR)$@ $< $(COREOBJS)

em: $(OBJDIR)BucketEMMain.o $(COREOBJS)
	$(CXX) $(FLAGS)  $(LIBS) -w -o $(BINDIR)$@ $< $(COREOBJS)

trim: $(OBJDIR)BucketTrim.o $(COREOBJS)
	$(CXX) $(FLAGS)  $(LIBS) -w -o $(BINDIR)$@ $< $(COREOBJS)

maketag: $(OBJDIR)MakeTAG.o $(COREOBJS)
	$(CXX) $(FLAGS)  $(LIBS) -w -o $(BINDIR)$@ $< $(COREOBJS)

