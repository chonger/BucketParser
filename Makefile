OSTAG=/home/chonger/data/osTAG/
TSGGRAM=$(OSTAG)eGrammar.txt
TAGGRAM=$(OSTAG)tagGrammar.txt
TAGGRAME=$(OSTAG)tagGrammar-estimated.txt
TAGGRAMT=$(OSTAG)tagGrammar-trimmed.txt
TOPARSE=$(OSTAG)23.yld
GOLD=$(OSTAG)23.unk
GOLDB=$(OSTAG)23.binarized.unk
TSGOUT=$(OSTAG)tsgout.txt
TAGOUT=$(OSTAG)tagout.txt
TRAIN=$(OSTAG)train.txt.unk
PCFG=$(OSTAG)train.pcfg.txt
MEMODEL=$(OSTAG)memodel

clean:
	make -C src clean

all:
#	make -C src/ec/. ecall
	make -C src/. bsall

pcfg:
	bin/pcfgtest $(PCFG) $(TOPARSE) $(GOLDB)

tsgparse:
	bin/parse $(TSGGRAM) 0 $(TOPARSE) $(TSGOUT) $(PCFG)

vg:
	valgrind --leak-check=full bin/parse $(TSGGRAM) 0 $(TOPARSE) $(TSGOUT) $(PCFG)

tsgeval:
	EVALB/evalb $(GOLD) $(TSGOUT)

tagparse:
	bin/parse $(TAGGRAMT) 1 $(TOPARSE) $(TAGOUT) $(PCFG)

tageval:
	EVALB/evalb $(GOLD) $(TAGOUT)

maketag:
	bin/maketag $(TSGGRAM) $(TAGGRAM)

em:
	bin/em $(TAGGRAM) 1 10 $(TRAIN) $(TAGGRAME)

em2:
	bin/em $(TAGGRAMT) 1 100 $(TRAIN) $(TAGGRAME)2

trim:
	bin/trim $(TAGGRAME)2 $(TAGGRAMT)2

trainme:
	bin/trainme $(PCFG) $(TRAIN) $(MEMODEL)

testme:
	bin/testme $(PCFG) $(MEMODEL) $(GOLDB)
