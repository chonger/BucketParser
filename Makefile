MOD=20

PTB=/home/chonger/data/PTB/
TSGGRAM=$(PTB)trainPTSG.txt$(MOD)
TAGGRAM=$(PTB)tagGrammar.txt$(MOD)
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
	bin/parse $(TSGGRAM) 0 $(TOPARSE) $(TSGOUT) $(PCFG)

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
	bin/em $(TAGGRAM) 1 10 $(TRAIN) $(TAGGRAM)E

em2:
	bin/em $(TAGGRAM)T 1 100 $(TRAIN) $(TAGGRAM)E2

trim2:
	bin/trim $(TAGGRAM)E2 $(TAGGRAM)T2

trim:
	bin/trim $(TAGGRAM)E $(TAGGRAM)T

trainme:
	bin/trainme $(PTBPCFG) $(PTB)23.txt.b.unk $(TESTCTF)

testme:
	bin/testme $(TESTCTF) 

batch: maketag em trim em2 tagparse
