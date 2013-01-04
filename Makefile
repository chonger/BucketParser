OSTAG=/home/chonger/data/osTAG/
TSGGRAM=$(OSTAG)eGrammar.txt
TAGGRAM=$(OSTAG)tagGrammar.txt
TAGGRAME=$(OSTAG)tagGrammar-estimated.txt
TOPARSE=$(OSTAG)23.yld
GOLD=$(OSTAG)23.unk
GOLDB=$(OSTAG)23.binarized.unk
TSGOUT=$(OSTAG)tsgout.txt
TRAIN=$(OSTAG)train.txt.unk
PCFG=$(OSTAG)train.pcfg.txt

all:
	make -C src/ec/. ecall
	make -C src/. bsall


pcfg:
	bin/pcfgtest $(PCFG) $(TOPARSE) $(GOLDB)

tsgparse:
	bin/parse $(TSGGRAM) 0 $(TOPARSE) $(TSGOUT) $(PCFG)

tsgeval:
	EVALB/evalb $(GOLD) $(TSGOUT)

maketag:
	bin/maketag $(TSGGRAM) $(TAGGRAM)

em:
	bin/em $(TAGGRAM) 1 10 $(TRAIN) $(TAGGRAME)

trim:
	bin/trim data/WSJ/TreebankTAG-EST.txt data/WSJ/TreebankTAG-TRIM.txt

