MOD=

PTB=/home/chonger/data/PTB/
TSGGRAM=$(PTB)trainPTSG.txt$(MOD)
TAGGRAM=$(PTB)tagGrammarSYM.txt$(MOD)X
TOPARSE=$(PTB)23.yld$(MOD)
GOLD=$(PTB)23.txt.unk$(MOD)
TSGOUT=$(PTB)tsgout.txt$(MOD)
TAGOUT=$(PTB)tagout.txt$(MOD)
TRAIN=$(PTB)train.txt.unk$(MOD)
PCFG=$(PTB)pcfg$(MOD).txt


all:
#	make -C src/ec/. ecall
	make -C src/. bsall
clean:
	make -C src clean

pcfgptb:
	bin/pcfgtest $(PTBPCFG) $(PTBYLD) $(PTBGOLD)

pcfg:
	bin/pcfgtest $(PCFG) $(TOPARSE) $(GOLDB)

tsgparse:
	bin/parse $(TSGGRAM) 0 $(TOPARSE) $(TSGOUT) 

vg:
	valgrind --leak-check=full bin/parse $(TSGGRAM) 0 $(TOPARSE) $(TSGOUT) $(PCFG)

tsgeval:
	EVALB/evalb $(GOLD) $(TSGOUT)

tagparse:
	bin/parse $(TAGGRAM)E2 1 $(TOPARSE) $(TAGOUT) $(PCFG)

tageval:
	EVALB/evalb $(GOLD) $(TAGOUT)

maketag:
	bin/maketag $(TSGGRAM) $(TAGGRAM)

em:
	bin/em $(TAGGRAM) 2 10 $(TRAIN) $(TAGGRAM)E

em2:
	bin/em $(TAGGRAM)T 2 100 $(TRAIN) $(TAGGRAM)E2

trim2:
	bin/trim $(TAGGRAM)E2 $(TAGGRAM)T2

trim:
	bin/trim $(TAGGRAM)E $(TAGGRAM)T

trainme:
	bin/trainme $(PTBPCFG) $(PTB)23.txt.b.unk $(TESTCTF)

testme:
	bin/testme $(TESTCTF) 

batch: maketag em trim em2 


DATA=/home/chonger/data/

vb:
	bin/vb $(DATA)PTB/pcfg10.txt $(DATA)PTB/train.txt.unk 5 $(DATA)PTB/vbout

valvb:
	valgrind --leak-check=full bin/vb $(DATA)toyG $(DATA)toyT2 5 $(DATA)toyO

#	valgrind --leak-check=full bin/vb $(DATA)PTB/pcfg10.txt $(DATA)PTB/train.txt.unk 5 $(DATA)PTB/vbout



vb2:
	bin/vb $(DATA)toyG $(DATA)toyT2 5 $(DATA)toyO
