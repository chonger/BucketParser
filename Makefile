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


all:
#	make -C src/ec/. ecall
	make -C src/. bsall
clean:
	make -C src clean

PTB=/home/chonger/data/PTB/
TESTCTF=$(OSTAG)testCTF
PTBTRAIN=$(PTB)train.txt.unk
PTBPCFG=$(PTB)train.pcfg.txt
PTBYLD=$(PTB)23.yld
PTBGOLD=$(PTB)23.txt.unk
PTBGRAM=$(PTB)elifPTSG.txt

stanparse:
	bin/parse ~/data/PTB/stanfordTrain.pcfgRAW.txt 0 ~/data/PTB/stan23.yld ~/data/PTB/stanout.txt

staneval:
	EVALB/evalb ~/data/PTB/23.txt.unk ~/data/PTB/stanout2.txt

pcfgptb:
	bin/pcfgtest $(PTBPCFG) $(PTBYLD) $(PTBGOLD)

pcfg:
	bin/pcfgtest $(PCFG) $(TOPARSE) $(GOLDB)

tsgparse:
	bin/parse $(PTBGRAM) 0 $(PTBYLD) $(TSGOUT) $(PTBPCFG)

vg:
	valgrind --leak-check=full bin/parse $(TSGGRAM) 0 $(TOPARSE) $(TSGOUT) $(PCFG)

tsgeval:
	EVALB/evalb $(PTBGOLD) $(TSGOUT)

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
	bin/trainme $(PTBPCFG) $(PTB)23.txt.b.unk $(TESTCTF)

testme:
	bin/testme $(TESTCTF) 
