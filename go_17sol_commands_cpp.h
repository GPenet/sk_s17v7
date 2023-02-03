
//========================================
const char * zh_g_cpt[10] = { "npuz", "guess", "close_d ", "upd1 ", "upd2 ",
"fupd ", "hpair ", "htripl ", " ", " " };

const char * libs_c17_00_cpt2g[100] = {
	"0 bands1+2 processed entries M10",//0
	"1 total bands 3",//1
	"2 steps external loop ",//2
	"3 3 clues ",	
	"4 6 clues ",	
	"5 7p clues last ",
	"6 active 6 clues",
	"7 set b12 ",
	"8 go b3",
	"9 min ok g2 g3",
	"10 uass b3 loaded go",
	"11 after uas b3 size 4 ",
	"12 entry after 11 miss0", 
	"13 entry after 11 others",
	"14 same after ua size>4",
	"15 GoNotMiss0() ",
	"16 still not miss o ",
	"17 miss 1 after [16]",
	"18 miss 1 after [17] ",
	"19 add guam",
	"20 addg2", 
	"21 addg3",
	"22  ",
	"23 ","24  ","25  ",
	"26  ","27 ","28","29 ",	
	"30 check valid2",
	"31 check valid3",
	"32 ","33 ","34 ",	"35 ",	"36 ",	"37 ",
	"38 subset seen", 
	"39 ",
	"40 miss0",	"41 miss1",	"42 miss2",	"43 miss3",
	"44 ",
	"45 miss0 to expand",	
	"46 miss0 go expand",
	"47 ",
	"48 miss1 before guam",
	"49 ",
	"50 temp test more1A",
	"51 full","52 below",
	"53 nfull<20",
	"54 size 7","55 size 8","56 size 9",
	"57 size 10","58 size 11",	"59  ",
	"60 nmiss0","61 nmiss1","62 nmiss2",
	"63 nmiss3","64 ","65 ","66 ","67 ","68 ","69 ",
	"70 entry goendall",
	"71 ",
	"72 exp 0-31",
	"73 exp 32-63",
	"74 exp 64-95",
	"75 exp 96_128",	
	"76 ","77 ","78 ","79 ",
	"80 ntoass<5 ","81 ntoass >=5",
	"82 D only b2","83  D only b1 ",
	"84 D 0-31","85 D 32_63","86 D 64_95 ",
	"87 D >95","88 ","89 ",
	"90 entry 10_12 ",
	"91active 10_12","92 10_12 avec 12","93 ","94 ",
	"95 n10",	"96 n11","97 n12","98 maxn12","99 ",

};
void Go_c17_00( ) {// p2 process
	cout << "Go_c17_00 search batch 17 clues  " << endl;
	op.SetUp(0);
	int it16_start = sgo.vx[0];
	g17b.aigstop=0;
	if (sgo.vx[2] < 0) {
		cerr << "invalid value for skip" << endl;
		return;
	}
	if (sgo.vx[3] < sgo.vx[2]) {
		cerr << "invalid value for last to process" << endl;
		return;
	}
	zh_g.modevalid = 1;
	zh_g2.grid0 = genb12.grid0;
	zh_g2.zsol = zh_g2.stdfirstsol;
	memset(p_cptg, 0, sizeof p_cptg);// used in debugging sequences only
	memset(p_cpt1g, 0, sizeof p_cpt1g);// used in debugging sequences only
	memset(p_cpt2g, 0, sizeof p_cpt2g);// used in debugging sequences only
	genb12.Start(0);
	genb12.NewBand1(op.b1);
	cout << "print final stats" << endl;
	for (int i = 0; i < 100; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
	cout << "exit final stats" << endl;
}

