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
CORE=BucketParser BucketGrammar BucketEM PCFGParser BucketVB VBGrammar
PROGS=parse em trim maketag
ECDIR=ec/
COREOBJS=$(addsuffix .o,$(addprefix $(OBJDIR),$(CORE))) #$(addprefix $(ECDIR),$(PARSERTEST_OBJS))
LIBS=-lgsl -lgslcblas -lm -pthread -lpthread #-lboost_thread
INCLUDE=-I../lib/sparsehash/include -I/home/chonger/Downloads/maxent-master/src 

FLAGS=-g
#FLAGS=-O3
CXX=g++

$(OBJDIR)%.o: %.cpp %.hpp
	$(CXX) $(FLAGS) $(INCLUDE)  -w -c $< -o $@

$(OBJDIR)%.o: %.cpp $(COREOBJS)
	$(CXX) $(FLAGS) $(INCLUDE)  -w -c $< -o $@
 

bsall: vb parse em trim maketag pcfgtest 

clean:
	rm -f $(OBJDIR)*
	rm -f $(BINDIR)*
	make -C ec ecclean

trainme : $(OBJDIR)TrainME.o $(COREOBJS)
	$(CXX) $(FLAGS)  -w -o $(BINDIR)$@ $< $(COREOBJS)  

vb : $(OBJDIR)VBMain.o $(COREOBJS)
	$(CXX) $(FLAGS)  -w -o $(BINDIR)$@ $< $(COREOBJS)

testme : $(OBJDIR)TestME.o $(COREOBJS)
	$(CXX) $(FLAGS)  -w -o $(BINDIR)$@ $< $(COREOBJS)  -L /usr/local/lib -l maxent -lm -lz -lgfortran

parse: $(OBJDIR)BucketParserMain.o $(COREOBJS)
	$(CXX) $(FLAGS)  $(LIBS) -w -o $(BINDIR)$@ $< $(COREOBJS)

em: $(OBJDIR)BucketEMMain.o $(COREOBJS)
	$(CXX) $(FLAGS)  $(LIBS) -w -o $(BINDIR)$@ $< $(COREOBJS)

trim: $(OBJDIR)BucketTrim.o $(COREOBJS)
	$(CXX) $(FLAGS)  $(LIBS) -w -o $(BINDIR)$@ $< $(COREOBJS)

maketag: $(OBJDIR)MakeTAG.o $(COREOBJS)
	$(CXX) $(FLAGS)  $(LIBS) -w -o $(BINDIR)$@ $< $(COREOBJS)

pcfgtest: $(OBJDIR)PCFGParserMain.o $(COREOBJS)
	$(CXX) $(FLAGS)  $(LIBS) -w -o $(BINDIR)$@ $< $(COREOBJS)

