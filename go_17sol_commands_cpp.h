
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
	"6 9 clues",
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
	"22  build9",
	"23 goa call",
	"24 after goa ",
	"25 gob call ",
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
	"91active 10_12","92 10_12 processed","93 ","94 ",
	"95 n10",	"96 n11","97 n12","98 maxn12","99 ",

};
uint8_t Gethexa(char c) {
	char ret = 0;
	if (c >= '0' && c <= '9')ret = c - '0';
	else if (c >= 'A' && c <= 'F')ret = c - 'A' + 10;
	return ret;
}
char* Char8out(uint8_t w) {
	static char ws[9];
	strncpy(ws, empty_puzzle, 8);
	ws[8] = 0;
	for (int j = 7; j >= 0; j--) if (w & (1 << j))			ws[j] = '1';
	return ws;

}
char* Char4out(uint8_t w) {
	static char ws[5];
	strncpy(ws, empty_puzzle, 4);
	ws[4] = 0;
	for (int j = 3; j >= 0; j--) if (w & (1 << j))			ws[j] = '1';
	return ws;

}





void Go_c17_00( ) {// p2 process
	cout << "Go_c17_00 search batch 17/18 clues  " << endl;
	op.SetUp(0);
	if (op.f18_status)if (genb12.F18_Init()) return;
	cout << "end test f18" << endl;
	//return;
	//int it16_start = sgo.vx[0];
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

void Go_c17_09() {// p2 process locate 

}
struct CLBS {// clues band stack
	uint16_t b[3];
	uint16_t s[3];
	inline void Init(uint64_t bf54, uint16_t n) {
		register uint64_t U = bf54;
		b[0] = (uint16_t)_popcnt64(U & BIT_SET_27);
		b[1] = n - b[0];
		b[2] = 0;
		register uint64_t S = stack1_54;
		s[0] = (uint16_t)_popcnt64(U & S);
		S <<= 3;
		s[1] = (uint16_t)_popcnt64(U & S);
		s[2] = n - s[0] - s[1];
	}
	inline void Add(uint32_t cell) {
		b[cell / 27]++;
		s[C_stack[cell]]++;
	}
	
	void Status() {
		cout << " bx " << b[0] << b[1] << b[2]
			<< " sx " << s[0] << s[1] << s[2];

	}
};


//========================= known s17 file 10/19
void Go_c17_10() {
	//op.SetUp(10,0);
	cout << "Go_10() search 17/18 using a file having known 17 656 " << endl;
	op.SetUp(10, 1);// setup with known
	cout << "back setup " << endl;

	zh_g.modevalid = 1;
	zh_g2.grid0 = genb12.grid0;
	zh_g2.zsol = zh_g2.stdfirstsol;
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	char* ze = finput.ze;
	int* zs0 = genb12.grid0, npuz = 0;
	char* ze2 = &ze[82];
	//if (op.t18) return;
	//if (op.p1) return;
	//cout << "Go_10() search 17/18 using a file having known 17 656 " << endl;

	//op.ton = 1;

	while (finput.GetLigne()) {
		if (strlen(ze) < 160) continue;// skip blank lines
		npuz++;
		g17b.npuz = npuz;
		g17b.a_17_found_here = 0;
		g17b.aigstop = 0;
		if (npuz < op.first) continue;
		if (npuz > op.last) break;
		CLBS cbs;
		g17b.p17diag.SetAll_0();
		int  nclues = 0;
		memset(&cbs, 0, sizeof cbs);
		for (int i = 0; i < 81; i++) if (ze2[i] != '.') {
			g17b.p17diag.Set_c(i);
			cbs.Add(i);
			nclues++;
		}

		if (0) {
			if (cbs.b[2] < 7) continue;
			fout1 << ze << endl;
			cout << &ze[82] << endl;
			continue;
		}

		if (op.ton)cout << ze << " to process  n=" << npuz
			<< " bands " << cbs.b[0] << cbs.b[1] << cbs.b[2]
			<< " stacks " << cbs.s[0] << cbs.s[1] << cbs.s[2] << endl;
		if (op.t18 && nclues == 17) {
			if (op.ton)cout << " see later 17 clues for a 18" << endl;
			continue;
		}
		if (op.t18) {
			op.p2b = 0;
		}
		if ((!op.t18) && op.p2) {
			op.p2b = 0;
			//================================ to avoid the 665 case
			if (cbs.b[2] == 5) {// change band3 <-> band2
				op.p2b = 1;
				for (int i = 0; i < 27; i++) {
					char temp = ze[i + 27];	ze[i + 27] = ze[i + 54];	ze[i + 54] = temp;
					temp = ze[i + 109];	ze[i + 109] = ze[i + 136];	ze[i + 136] = temp;
				}
				uint32_t w = g17b.p17diag.bf.u32[1];
				g17b.p17diag.bf.u32[1] = g17b.p17diag.bf.u32[2];
				g17b.p17diag.bf.u32[2] = w;

			}


		}

		cout << "\n\nto process  n=" << dec << npuz << " debug=" << op.ton << endl;
		if (op.ton)		cout << ze << " to process  n=" << npuz
			<< " bands " << cbs.b[0] << cbs.b[1] << cbs.b[2]
			<< " stacks " << cbs.s[0] << cbs.s[1] << cbs.s[2] << endl;
		// =======================morph entry 
		for (int i = 0; i < 81; i++)zs0[i] = ze[i] - '1';
		BANDMINLEX::PERM perm_ret;
		bandminlex.Getmin(zs0, &perm_ret);
		int ib1 = perm_ret.i416, ib1a = t416_to_n6[ib1];
		myband1.InitBand2_3(ib1, ze, perm_ret, 0);
		bandminlex.Getmin(&zs0[27], &perm_ret);
		int ib2 = perm_ret.i416, ib2a = t416_to_n6[ib2];
		myband2.InitBand2_3(ib2, &ze[27], perm_ret, 1);
		bandminlex.Getmin(&zs0[54], &perm_ret);
		int ib3 = perm_ret.i416, ib3a = t416_to_n6[ib3];
		genb12.bands3[0].InitBand3(ib3, &ze[54], perm_ret);
		genb12.nband3 = 1;
		genb12.i1t16 = ib1a;
		genb12.i2t16 = ib2a;
		ze[81] = 0;
		if (op.ton)
			cout << Char2Xout(g17b.p17diag.bf.u64[0]) << " b12 pattern for the 17" << endl;
		register uint64_t U = g17b.p17diag.bf.u64[0];
		g17b.pk54 = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);
		genb12.ValidInitGang();
		g17b.npuz = npuz;
		if (sgo.bfx[2] & 1)g17b.StartKnown();
		else	g17b.Start();
		//g17b.a_17_found_here = 1;
		if (!g17b.a_17_found_here) {
			cout << "puz=" << npuz << " failed to find the searched 17" << endl;
			cerr << "puz=" << npuz << " failed to find the searched 17" << endl;
			break;
		}
	}
	cout << "print final stats" << endl;
	for (int i = 0; i < 100; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
}


void Go_c17_11() {// extract cout xx5 file1 xx6
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	char * ze = finput.ze;
	int * zs0 = genb12.grid0, npuz = 0;
	cout << "Go_c17_11() split known 17 b3_5 b3_6 " << endl;
	while (finput.GetLigne()) {
		// =======================morph entry to have min n6 count in first
		for (int i = 0; i < 81; i++)zs0[i] = ze[i] - '1';
		BANDMINLEX::PERM perm_ret;
		bandminlex.Getmin(zs0, &perm_ret);
		int ib1 = perm_ret.i416, ib1a = t416_to_n6[ib1];
		bandminlex.Getmin(&zs0[27], &perm_ret);
		int ib2 = perm_ret.i416, ib2a = t416_to_n6[ib2];
		bandminlex.Getmin(&zs0[54], &perm_ret);
		int ib3 = perm_ret.i416, ib3a = t416_to_n6[ib3];
		if (ib1a > ib2a) {// change band2 <-> band1
			for (int i = 0; i < 27; i++) {
				char temp = ze[i];	ze[i] = ze[i + 27];	ze[i + 27] = temp;
				temp = ze[i + 82];	ze[i + 82] = ze[i + 109]; ze[i + 109] = temp;
				int w = ib1a;	ib1a = ib2a;	ib2a = w;
			}
		}
		if (ib1a > ib3a) {// change band3 <-> band1
			for (int i = 0; i < 27; i++) {
				char temp = ze[i];	ze[i] = ze[i + 54];		ze[i + 54] = temp;
				temp = ze[i + 82];	ze[i + 82] = ze[i + 136];	ze[i + 136] = temp;
				int w = ib1a;	ib1a = ib3a;	ib3a = w;
			}
		}
		int ncb3 = 0;
		for (int i = 0; i < 27; i++) {
			if (ze[i + 136] - '.')ncb3++;
		}
		if (ncb3 == 6) {
			fout1 <<&ze[82] << ";" << ib2a 
				<< ";" << ib3a <<endl;
				continue;
		}
		// now 5 clues in band 3 exchange bands
		for (int i = 0; i < 27; i++) {
			char temp = ze[i + 27];	ze[i + 27] = ze[i + 54];	ze[i + 54] = temp;
				temp = ze[i + 109];	ze[i + 109] = ze[i + 136];	ze[i + 136] = temp;
		}
		int w = ib2a;	ib2a = ib3a;	ib3a = w;
		cout << &ze[82]  << ";" << ib2a
			<< ";" << ib3a << endl;

	}

}
void BandReShapewithknown(int* s, int* d, char* sk, char* dk, BANDMINLEX::PERM p);
void Go_c17_12() {// check diagonal status in a 665
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	char* ze = finput.ze;
	int  zs0[81], zs0d[81],
		npuz = 0;
	cout << "Go_c17_12() band analysis option=" << sgo.vx[2] << endl;
	while (finput.GetLigne()) {
		if (strlen(ze) < 81)continue;
		// =======================morph entry to have min n6 count in first
		for (int i = 0; i < 81; i++) {
			zs0[i] = ze[i] - '1';
			zs0d[i] = ze[C_transpose_d[i]] - '1';
		}
		BANDMINLEX::PERM perm_ret;
		bandminlex.Getmin(zs0, &perm_ret);
		int ib1 = perm_ret.i416, ib1a = t416_to_n6[ib1];
		bandminlex.Getmin(&zs0[27], &perm_ret);
		int ib2 = perm_ret.i416, ib2a = t416_to_n6[ib2];
		bandminlex.Getmin(&zs0[54], &perm_ret);
		int ib3 = perm_ret.i416, ib3a = t416_to_n6[ib3];
		bandminlex.Getmin(zs0d, &perm_ret);
		int ib1d = perm_ret.i416, ib1ad = t416_to_n6[ib1d];
		bandminlex.Getmin(&zs0d[27], &perm_ret);
		int ib2d = perm_ret.i416, ib2ad = t416_to_n6[ib2d];
		bandminlex.Getmin(&zs0d[54], &perm_ret);
		int ib3d = perm_ret.i416, ib3ad = t416_to_n6[ib3d];

		if (0) {
			cout << ze << ";\t" << ib1 << ";" << ib2 << ";" << ib3
				<< ";\t" << ib1a << ";" << ib2a << ";" << ib3a
				<< ";\t" << ib1d << ";" << ib2d << ";" << ib3d
				<< ";\t" << ib1ad << ";" << ib2ad << ";" << ib3ad
				<< endl;
		}
		switch (sgo.vx[2]) {
		case 0: {cout << ze << ";\t" << ib1 << ";" << ib2 << ";" << ib3
			<< ";\t" << ib1a << ";" << ib2a << ";" << ib3a
			<< ";\t" << ib1d << ";" << ib2d << ";" << ib3d
			<< ";\t" << ib1ad << ";" << ib2ad << ";" << ib3ad
			<< endl;
			break;
		}
		case 1: {// extract 18 (8 in [87]
			if (strlen(ze) < 88)break;
			if (ze[87] == '8')
				fout1 << ze << endl;
			break;
		}

		case 2: {// extract 18 (8 in [87]) plus info
			if (strlen(ze) < 88)break;
			if (ze[87] == '8')
				fout1 << ze << ";\t" << ib1 << ";" << ib2 << ";" << ib3
				<< ";\t" << ib1a << ";" << ib2a << ";" << ib3a
				<< ";\t" << ib1d << ";" << ib2d << ";" << ib3d
				<< ";\t" << ib1ad << ";" << ib2ad << ";" << ib3ad
				<< endl;
			break;
		}
		case 3: {// sol+18 to analyze get pass2
			if (strlen(ze) < 162)break;
			if (!sgo.vx[3])sgo.vx[3] = 82;
			char* ze2 = &ze[sgo.vx[3]];
			CLBS cbs;
			int  nclues = 0;
			memset(&cbs, 0, sizeof cbs);
			for (int i = 0; i < 81; i++) if (ze2[i] != '.') {
				cbs.Add(i);
			}
			if (cbs.b[0] == 6 && cbs.b[1] == 6 )
				fout1 << ze2 << ";\t" << ib1a << ";" << ib2a << ";" << ib3a
				<< ";\t" << ib1ad << ";" << ib2ad << ";" << ib3ad
				<< "; bands; " << cbs.b[0] << cbs.b[1] << cbs.b[2]
				<< "; stacks; " << cbs.s[0] << cbs.s[1] << cbs.s[2] << endl;
			break;
		}
		case 4: {// sol+18 to analyze get pass1
			if (strlen(ze) < 162)break;
			if (!sgo.vx[3])sgo.vx[3] = 82;
			char* ze2 = &ze[sgo.vx[3]];
			CLBS cbs;
			int  nclues = 0;
			memset(&cbs, 0, sizeof cbs);
			for (int i = 0; i < 81; i++) if (ze2[i] != '.') {
				cbs.Add(i);
			}

			if (cbs.b[0] == 6 && cbs.b[1] == 6 && cbs.s[0] == 6 && cbs.s[1] == 6)
				break;// pass2
			uint16_t x = cbs.s[0];
			for (int i = 1; i < 3; i++) if (x < cbs.s[i])x = cbs.s[i];
			if (cbs.b[2] < cbs.b[0]) break;
			if (cbs.b[2] < cbs.b[1]) break;
			if (cbs.b[2] < x) break;
			fout1 << ze2 << ";\t" << ib1a << ";" << ib2a << ";" << ib3a
				<< ";\t" << ib1ad << ";" << ib2ad << ";" << ib3ad
				<< "; bands; " << cbs.b[0] << cbs.b[1] << cbs.b[2]
				<< "; stacks; " << cbs.s[0] << cbs.s[1] << cbs.s[2] << endl;
			break;
		}
		case 5: {// sol+18 extract  pass2
			if (strlen(ze) < 162)break;
			if (!sgo.vx[3])sgo.vx[3] = 82;
			char* ze2 = &ze[sgo.vx[3]];
			CLBS cbs;
			int  nclues = 0;
			memset(&cbs, 0, sizeof cbs);
			for (int i = 0; i < 81; i++) if (ze2[i] != '.') {
				cbs.Add(i);
			}
			if (cbs.b[0] == 6 && cbs.b[1] == 6)
				fout1 << ze << ";bands 666 "  << endl;
			break;
		}
		case 6: {// morph to CF sol+known
			cout << ze << endl;
			if (strlen(ze) < 162)break;
			bandminlex.Getmin(zs0, &perm_ret);
			int ib1 = perm_ret.i416, ib1a = t416_to_n6[ib1],
				gw[81];
			char zks[82]; zks[81] = 0;
			myband1.Initstd();
			myband1.InitG12(ib1);
			memcpy(gw, myband1.band0, sizeof myband1.band0);
			char* zke = &ze[82];
			BandReShapewithknown(zs0, gw, zke, zks, perm_ret);
			BandReShapewithknown(&zs0[27], &gw[27], &zke[27], &zks[27], perm_ret);
			BandReShapewithknown(&zs0[54], &gw[54], &zke[54], &zks[54], perm_ret);
			for (int i = 0; i < 81; i++)fout1 << gw[i] + 1;
			fout1 << ";" << zks << endl;
	
			break;
		}
		case 7: {// add info at the end dll output
			if (strlen(ze) < 162)break;		char* ze2 = &ze[strlen(ze)];
			CLBS cbs;	int  nclues = 0;	memset(&cbs, 0, sizeof cbs);
			for (int i = 0; i < 81; i++) if (ze[i+82] != '.') 	cbs.Add(i);			
			fout1<<&ze[82] << ";" << ib1a << ";" << ib2a << ";" << ib3a
				<< "; bands; " << cbs.b[0] << ";" << cbs.b[1] << ";" << cbs.b[2]
				 << endl;
			break;
		}
		case 8: {// same want 6 items in output
			if (strlen(ze) < 162)break;		char* ze2 = &ze[strlen(ze)];
			CLBS cbs;	int  nclues = 0;	memset(&cbs, 0, sizeof cbs);
			for (int i = 0; i < 81; i++) if (ze[i + 82] != '.') 	cbs.Add(i);
			fout1 << &ze[82] << ";" << ib1a << ";" << ib2a << ";" << ib3a
				<< "; " << cbs.b[0] << cbs.b[1] << cbs.b[2]		<< endl;
			break;
		}
		}
	}

}
void BandReShapewithknown(int* s, int* d, char* sk, char*dk, BANDMINLEX::PERM p) {
	int* pc = p.cols;
	//uint8_t* pc = p.cols;
	for (int irow = 0; irow < 3; irow++) {
		int drow = 9 * irow;
		for (int icol = 0; icol < 9; icol++) {
			int x = p.map[s[drow + pc[icol]]];
			char c= sk[drow + pc[icol]],
				cx = (char)x + '1';
			d[drow + icol] =x;
			if (c == '.')dk[drow + icol] = '.';
			else dk[drow + icol] = cx;
		}
	}
	int temp[9];// now reorder 
	char tempk[9];
	if (d[0] > d[9]) {
		memcpy(temp, &d[0], sizeof temp);
		memcpy(&d[0], &d[9], sizeof temp);
		memcpy(&d[9], temp, sizeof temp);
		memcpy(tempk, &dk[0], sizeof tempk);
		memcpy(&dk[0], &dk[9], sizeof tempk);
		memcpy(&dk[9], tempk, sizeof tempk);
	}
	if (d[0] > d[18]) {
		memcpy(temp, &d[0], sizeof temp);
		memcpy(&d[0], &d[18], sizeof temp);
		memcpy(&d[18], temp, sizeof temp);
		memcpy(tempk, &dk[0], sizeof tempk);
		memcpy(&dk[0], &dk[18], sizeof tempk);
		memcpy(&dk[18], tempk, sizeof tempk);
	}
	if (d[9] > d[18]) {
		memcpy(temp, &d[9], sizeof temp);
		memcpy(&d[9], &d[18], sizeof temp);
		memcpy(&d[18], temp, sizeof temp);
		memcpy(tempk, &dk[9], sizeof tempk);
		memcpy(&dk[9], &dk[18], sizeof tempk);
		memcpy(&dk[18], tempk, sizeof tempk);
	}
}
void BandReShape(int* s, int* d, BANDMINLEX::PERM p) {
	int * pc = p.cols;
	//uint8_t* pc = p.cols;
	for (int irow = 0; irow < 3; irow++) {
		int drow = 9 * irow;
		for (int icol = 0; icol < 9; icol++)
			d[drow + icol] = p.map[s[drow + pc[icol]]];
	}
	int temp[9];// now reorder 
	if (d[0] > d[9]) {
		memcpy(temp, &d[0], sizeof temp);
		memcpy(&d[0], &d[9], sizeof temp);
		memcpy(&d[9], temp, sizeof temp);
	}
	if (d[0] > d[18]) {
		memcpy(temp, &d[0], sizeof temp);
		memcpy(&d[0], &d[18], sizeof temp);
		memcpy(&d[18], temp, sizeof temp);
	}
	if (d[9] > d[18]) {
		memcpy(temp, &d[9], sizeof temp);
		memcpy(&d[9], &d[18], sizeof temp);
		memcpy(&d[18], temp, sizeof temp);
	}
}
void BandReOrder(int* d) {
	int temp[9];// now reorder 
	if (d[0] > d[9]) {
		memcpy(temp, &d[0], sizeof temp);
		memcpy(&d[0], &d[9], sizeof temp);
		memcpy(&d[9], temp, sizeof temp);
	}
	if (d[0] > d[18]) {
		memcpy(temp, &d[0], sizeof temp);
		memcpy(&d[0], &d[18], sizeof temp);
		memcpy(&d[18], temp, sizeof temp);
	}
	if (d[9] > d[18]) {
		memcpy(temp, &d[9], sizeof temp);
		memcpy(&d[9], &d[18], sizeof temp);
		memcpy(&d[18], temp, sizeof temp);
	}
}
BANDMINLEX::PERM t_auto20[108];

void Go_c17_20() {// Morph to clear redundancy in pass1 entry file
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	char* ze = finput.ze;
	int  zs0[81],  grid0[81],gw[81],
		m1t16, m2t16, m3t16,		// min index 0_415
		i1t16, i2t16, i3t16, npuz = 0,n_auto20=0; // index n6
	STD_B3 bands3[512],mb1;
	i1t16 = sgo.vx[0];	
	cout << "Go_c17_20() band analysis band2 index=" << sgo.vx[2] << endl;
	if (i1t16>415) {
		cout << "invalid band 3 to study" << endl;
	}
	m1t16 = tn6_to_416[i1t16];
	myband1.Initstd();
	myband1.InitG12(m1t16);
	memcpy(grid0, myband1.band0, sizeof myband1.band0);
	n_auto20 = bandminlex.GetAutoMorphs(m1t16, t_auto20);
	if (!n_auto20) {
		cout << "no auto morphism stop" << endl;
		return;
	}
	cout << myband1.band << "index=" << i1t16 << " 0_415=" << m1t16
		<< " n auto morphs=" << n_auto20 << endl;
	while (finput.GetLigne()) {
		int tmini[20], nmini = 1; tmini[0] = -1;
		if (strlen(ze) < 81)continue;
		npuz++;
		//cout << ze << endl;
		// =======================morph entry to have min n6 count in first
		for (int i = 0; i < 81; i++) {
			zs0[i] = ze[i] - '1';
		}
		BANDMINLEX::PERM perm_retb3,perm_retb2,perm_retb1;
		bandminlex.Getmin(&zs0[54], &perm_retb3);

		if (perm_retb3.i416 != m1t16) {
			cout << "illegal band 3 stop" << endl;
			return;
		}
		//continue;
		bandminlex.Getmin(zs0, &perm_retb1);
		m2t16 = perm_retb1.i416;
		i2t16 = t416_to_n6[m2t16];
		bandminlex.Getmin(&zs0[27], &perm_retb2);
		m3t16 = perm_retb2.i416;
		i3t16 = t416_to_n6[m3t16];
		memcpy(gw, myband1.band0, sizeof myband1.band0);
		//for (int i = 0; i < 27; i++) cout << gw[i] + 1;
		//cout <<endl;
		BandReShape(zs0, &gw[27], perm_retb3);
		BandReShape(&zs0[27], &gw[54], perm_retb3);
		for (int i = 0; i < 81; i++) cout << gw[i] + 1;		
		cout << ";"<< i1t16<<";" << i2t16 << ";" << i3t16 << endl;
		// look for the slowest lexical band 2
		int band2min[27];
		memcpy(band2min, &gw[27], sizeof band2min);
		for (int i = 0; i < 27; i++) cout << band2min[i] + 1;
		cout << " band2min" << endl;
		for (int imorph = 0; imorph < n_auto20; imorph++) {
			BANDMINLEX::PERM& p = t_auto20[imorph];
			int band[27], * z0 = &gw[27];// morph the band
			for (int i = 0; i < 9; i++) {
				band[i] = p.map[z0[p.cols[i]]];
				band[i + 9] = p.map[z0[p.cols[i] + 9]];
				band[i + 18] = p.map[z0[p.cols[i] + 18]];
			}
			int tsort[3], bandw[27], aig = 0;;
			G17BuildSort(band, tsort);
			for (int i = 0,k=0; i < 3; i++) {
				int* b = &band[9 * (tsort[i] & 3)];
				for (int j = 0; j < 9; j++,k++) bandw[k] = b[j];
			}
			for (int i = 0; i < 27; i++) 
				if ((aig = bandw[i] - band2min[i])) 	break;	
			if (aig > 0)continue; // not minimal
			for (int i = 0; i < 27; i++) cout << bandw[i] + 1;
			cout << "imorph=" << imorph << " aig" << aig << endl;
			if (!aig) {		tmini[nmini++] = imorph;		continue;	}
			// now a lower 
			nmini = 0;
			tmini[nmini++] = imorph;
			memcpy(band2min, bandw, sizeof band2min);
		}
		if (tmini[0] < 0) {//keep source
			for(int i=27;i<81;i++) fout1 << gw[i] + 1;
			fout1<<";" << i1t16 << ";" << i2t16 << ";" << i3t16 << ";np=" << npuz <<"; source "<<ze << endl;
			continue;
		}
		// morph band 3 to perm
		for (int im = 0; im < nmini; im ++ ) {
			BANDMINLEX::PERM& p = t_auto20[tmini[im]];
			int band3[27], * z0 = &gw[54];// morph the band
			for (int i = 0; i < 9; i++) {
				band3[i] = p.map[z0[p.cols[i]]];
				band3[i + 9] = p.map[z0[p.cols[i] + 9]];
				band3[i + 18] = p.map[z0[p.cols[i] + 18]];
			}

			int tsort[3], bandw[27], aig = 0;;
			G17BuildSort(band3, tsort);

			for (int i = 0, k = 0; i < 3; i++) {
				int* b = &band3[9 * (tsort[i] & 3)];
				for (int j = 0; j < 9; j++, k++) bandw[k] = b[j];
			}
			for (int j = 0; j < 27; j++) fout1 << band2min[j] + 1;
			for (int j = 0; j < 27; j++) fout1 << bandw[j] + 1;
			fout1 <<";" << i1t16 << ";" << i2t16 << ";" << i3t16 << ";np=" << npuz << " morph " << im << " " << ze << endl;
		}

		// is now <= print for analysis purpose 

		/*
int G17ComparedOrderedBand(int * zs0, int * band){
	int tsort[3];
	G17BuildSort(band, tsort);
	// compare now the morph to the source exit if lower
	for (int i = 0; i < 3; i++){
		int *b = &band[9 * (tsort[i] & 3)], *z = &zs0[9 * i];
		for (int j = 0; j < 9; j++){
			if (b[j] < z[j]) return 1;// not a minimal band 1+2
			if (b[j] > z[j]) return 2;; // try next
		}
	}
	return 0;
}		*/
	}
}

//=========================== enumeration test


void Go_c17_80() {// enumeration test
	//return; // revise command line
	cout << "Go_c17_80 phase 2a enumeration test " << endl;
	cout << sgo.vx[0] << " -v0- first id 0_415" << endl;
	cout << sgo.vx[1] << " -v1- second id 0_415" << endl;
	cout << sgo.vx[2] << " -v2- if 1 printout asked" << endl;
	cout << sgo.vx[4] << " -v4- 0 mode p2a 1 mode p2b 2 mode p1" << endl;
	cout << sgo.vx[5] << " -v5- band 2 index if not 9990" << endl;
	cout << sgo.vx[11] << " -vx- band 3 index if not 999" << endl;
	int it16_start = sgo.vx[0], it16_end = sgo.vx[1];
	memset(&op, 0, sizeof op);
	op.first = 0;
	op.last =200000;
	op.t18 = 1;
	op.ton= sgo.vx[2];
	op.b2 = sgo.vx[5];
	op.bx3 = sgo.vx[11];
	op.out_entry = op.ton;
	int x = (int)sgo.vx[4];
	if (x == 2) op.p1 = 1;
	else {
		op.p2 = 1; if (x == 1) op.p2b = 1;
	}
	if (it16_start > 415 || it16_end > 415 || it16_start > it16_end) {
		cerr << "invalid it16_start it16_end" << endl;
		return;
	}
	memset(p_cptg, 0, sizeof p_cptg);// used in debugging sequences only

	genb12.Start(11);// mode p2a or p2b or p1
	for (int i1t16 = it16_start; i1t16 <= it16_end; i1t16++) {
		memset(p_cpt, 0, sizeof p_cpt);// band2 and band 3 count
		genb12.nb12 = 0;
		genb12.NewBand1(i1t16);
		cout << genb12.i1t16 << "\t" << genb12.it16 << "\t" << p_cpt[0]
			<< "\t" << p_cpt[1] << endl;
		p_cptg[0] += p_cpt[0];
		p_cptg[1] += p_cpt[1];
	}
	cout << "total\t\t" << p_cptg[0]
		<< "\t" << p_cptg[1] << "\t" << p_cpt2g[10] << "\t" << p_cpt2g[11] 
		<< "\t" << p_cpt2g[12] << "\t" << p_cpt2g[13] 
		<< "\t" << p_cpt2g[14] << "\t" << p_cpt2g[15]
		<<" max nsgcheck " << p_cpt2g[20] << endl;
}

void Go_c17_81() {// enumeration split of known
	cout << "Go_c17_84 split in slices " << endl;
	cout << sgo.vx[0] << " -v0-  id 0_415" << endl;
	if (sgo.vx[0] > 415) {
		cerr << "invalid it16_start it16_end" << endl;
		return;
	}
	cout << sgo.vx[4] << " -v4- 0 mode p2a 1 mode p2b 2 mode p1" << endl;
	memset(&op, 0, sizeof op);
	op.last = 200000;
	op.t18 = 1;
	op.out_entry = -1;
	op.b2 = sgo.vx[5];
	op.bx3 = sgo.vx[11];
	int x = (int)sgo.vx[4];
	if (x == 2) op.p1 = 1;
	else {	op.p2 = 1; if (x == 1) op.p2b = 1;	}

	// Load in memory 
	//uint8_t bitfield_sgs[200000000 / 8];
	//uint64_t nbitfield_sgs;
	char* ze = finput.ze;
	nbitfield_sgs = 0;
	uint64_t nl=0,idl=0,nbits=0,ntot=0,n;
	while (finput.GetLigne()) {
		nl++;
		idl = nbitfield_sgs;
		n =(int) strlen(ze);
		if (n & 1) {	ze[n] = '0'; n++; ze[n] = 0;}
		for (int i = 0; i < n; i += 2) {
			uint8_t c1 = Gethexa(ze[i]), c2 = Gethexa(ze[i + 1]); // one byte
			c1 |= (c2 << 4);
			bitfield_sgs[nbitfield_sgs++] = c1;// store byte
			ntot += 8;
			nbits += _popcnt32(c1);
		}
	}
	cout << "nl= " << nl << " idl=" << idl <<" ntot="<<ntot<<" nbits="<<nbits
		<< " nbitfield_sgs=" << nbitfield_sgs << endl;
	memset(p_cptg, 0, sizeof p_cptg);// used in debugging sequences only
	memset(p_cpt2g, 0, sizeof p_cpt2g);// used in debugging sequences only
	genb12.Start(11);// mode p2a or p2b or p1
	memset(p_cpt, 0, sizeof p_cpt);// band2 and band 3 count
	genb12.nb12 = 0;
	genb12.NewBand1(sgo.vx[0]);
	cout <<"last " << p_cpt[0] << "\t" << p_cpt[1] << endl;

}

void Go_c17_90() {  // JIM

	// Process a list of SG's

	char pgrd[128];

	op.SetUp(10, 0);    // get -b0- (pass 1 or 2)

	zh_g2.grid0 = genb12.grid0;
	zh_g2.zsol = zh_g2.stdfirstsol;
	zh_g.modevalid = 1;

	char ze[128];
	//char ws[128];

	int* zs0 = genb12.grid0;
	int ib1, ib2, ib3;
	int ib1a, ib2a, ib3a;
	//int rlen, nsg;
	int nb3 = 0;    // nband3

	BANDMINLEX::PERM perm_ret;

	FILE* SGF = fopen(sgo.finput_name, "r");

	if (!SGF) {
		cout << "[c90] not found: " << sgo.finput_name << endl;
		return;
	}

	cout << "[c90] Processing " << sgo.finput_name << endl;

	while (fgets(ze, 120, SGF)) {
		if (strlen(ze) < 81) continue;  // skip blank lines
		if (ze[0] != '1') continue;     // skip grids with 1st byte not '1'

		g17b.aigstop = 0;
		genb12.nb12 = atoi(&ze[94]) << 6;  // reproduce original slice number

		// ======================= morph entry 

		if (memcmp(ze, pgrd, 27)) {    // new band1
			if (nb3) {   // don't change grid0/zsol yet!
				genb12.nband3 = nb3;
				genb12.ValidInitGang();
				g17b.Start();
			}

			nb3 = 0;
			for (int i = 0; i < 81; i++) zs0[i] = ze[i] - '1';
			bandminlex.Getmin(zs0, &perm_ret);
			ib1 = perm_ret.i416; ib1a = t416_to_n6[ib1];
			genb12.i1t16 = ib1a;
			memset(&pgrd[27], 0, 27);  // new b1 always forces new b2 
		}

		if (memcmp(&ze[27], &pgrd[27], 27)) { // new band2
			if (nb3) {
				genb12.nband3 = nb3;
				genb12.ValidInitGang();
				g17b.Start();
			}
			nb3 = 0;
			for (int i = 0; i < 81; i++) zs0[i] = ze[i] - '1';
			myband1.InitBand2_3(ib1, ze, perm_ret, 0);
			bandminlex.Getmin(&zs0[27], &perm_ret);
			ib2 = perm_ret.i416; ib2a = t416_to_n6[ib2];
			myband2.InitBand2_3(ib2, &ze[27], perm_ret, 1);
			genb12.i2t16 = ib2a;
		}

		for (int i = 0; i < 81; i++) zs0[i] = ze[i] - '1';
		bandminlex.Getmin(&zs0[54], &perm_ret);
		ib3 = perm_ret.i416; ib3a = t416_to_n6[ib3];
		genb12.bands3[nb3++].InitBand3(ib3, &ze[54], perm_ret);
		memcpy(pgrd, ze, 81);
	}

	fclose(SGF);

	if (nb3) {   // process the last b12 set
		genb12.nband3 = nb3;
		genb12.ValidInitGang();
		g17b.Start();
	}

	cout << "print final stats" << endl;
	for (int i = 0; i < 100; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
}


void BuildAutoMorph(); // in go_17sol_tables.cpp

//====== 15  extraction des 17 6 6 5
void Go_c17_15() {
	cout << "process 15" << endl;
	BuildAutoMorph();
	return;

	char * ze = finput.ze;

	while (finput.GetLigne()) {
		int count[6];
		memset(count, 0, sizeof count);
		for (int i = 0; i < 81; i++) {
			if (ze[i] - '.') {
				int band = i / 27, stack = C_stack[i] + 3;
				++count[band];
				++count[stack];
			}
		}
		if (count[0] > 6 || count[1] > 6 ||
			count[2] > 6 || count[3] > 6 ||
			count[4] > 6 || count[5] > 6)
			cout << ze << endl;
		else fout1 << ze << endl;
	}
}

void Msp_ReorderBand(char * ze,char * zep=0)
{
	char temp[9];
	if (ze[0] > ze[9]) {
		memmove(temp, ze, sizeof temp);
		memmove(ze, &ze[9], sizeof temp);
		memmove(&ze[9], temp,sizeof temp);
		if (zep) {
			memmove(temp, zep, sizeof temp);
			memmove(zep, &zep[9], sizeof temp);
			memmove(&zep[9], temp, sizeof temp);

		}
	}
	if (ze[0] > ze[18]) {
		memmove(temp, ze, sizeof temp);
		memmove(ze, &ze[18], sizeof temp);
		memmove(&ze[18], temp, sizeof temp);
		if (zep) {
			memmove(temp, zep, sizeof temp);
			memmove(zep, &zep[18], sizeof temp);
			memmove(&zep[18], temp, sizeof temp);
		}
	}
	if (ze[9] > ze[18]) {
		memmove(temp, &ze[9], sizeof temp);
		memmove(&ze[9], &ze[18], sizeof temp);
		memmove(&ze[18], temp, sizeof temp);
		if (zep) {
			memmove(temp, &zep[9], sizeof temp);
			memmove(&zep[9], &zep[18], sizeof temp);
			memmove(&zep[18], temp, sizeof temp);
		}
	}
}

void Go_c17_16() {
	char * ze = finput.ze,zout[82],zdiag[82],zew[82],zes[82], zesf[82];
	zout[81] = zdiag[81] = zes[81] = 0;
	zh_g2.zsol = zout;
	while (finput.GetLigne()) {
		if (zhou[0].CheckValidityQuick(ze) == 1) {
			int count[6];
			memset(count, 0, sizeof count);
			for (int i = 0; i < 81; i++) {
				if (ze[i] - '.') {
					int band = i / 27, stack = C_stack[i] + 3;
					++count[band];
					++count[stack];
				}
			}
			for (int i = 0; i < 81; i++)
				zdiag[C_transpose_d[i]] = zout[i];
			int zw[27], indexr[6], indexmin = 500;
			BANDMINLEX::PERM pr;
			char * zz = zout;
			for (int i = 0; i < 27; i++)
				zw[i] = zz[i] - '1';
			bandminlex.Getmin(zw, &pr);
			indexmin = indexr[0] = t416_to_n6[pr.i416];
			zz = &zout[27];
			for (int i = 0; i < 27; i++)
				zw[i] = zz[i] - '1';
			bandminlex.Getmin(zw, &pr);
			indexr[1] = t416_to_n6[pr.i416];
			if (indexmin > indexr[1])indexmin = indexr[1];
			zz = &zout[54];
			for (int i = 0; i < 27; i++)
				zw[i] = zz[i] - '1';
			bandminlex.Getmin(zw, &pr);
			indexr[2] = t416_to_n6[pr.i416];
			if (indexmin > indexr[2])indexmin = indexr[2];
			zz = zdiag;
			for (int i = 0; i < 27; i++)
				zw[i] = zz[i] - '1';
			bandminlex.Getmin(zw, &pr);
			indexr[3] = t416_to_n6[pr.i416];
			if (indexmin > indexr[3])indexmin = indexr[3];
			zz = &zdiag[27];
			for (int i = 0; i < 27; i++)
				zw[i] = zz[i] - '1';
			bandminlex.Getmin(zw, &pr);
			indexr[4] = t416_to_n6[pr.i416];
			if (indexmin > indexr[4])indexmin = indexr[4];
			zz = &zdiag[54];
			for (int i = 0; i < 27; i++)
				zw[i] = zz[i] - '1';
			bandminlex.Getmin(zw, &pr);
			indexr[5] = t416_to_n6[pr.i416];
			if (indexmin > indexr[5])indexmin = indexr[5];
			// sort and keep track of the count or split the file


			if (sgo.vx[6] && indexmin != sgo.vx[6]) continue;
			uint32_t tx[3], tdx[3], temp;
			tx[0] = (indexr[0] << 8) | count[0];
			tx[1] = (indexr[1] << 8) | count[1]|010;
			tx[2] = (indexr[2] << 8) | count[2]|020;
			tdx[0] = (indexr[3] << 8) | count[3];
			tdx[1] = (indexr[4] << 8) | count[4]|010;
			tdx[2] = (indexr[5] << 8) | count[5]|020;
			if (tx[0] > tx[1]) { temp = tx[0]; tx[0] = tx[1]; tx[1] = temp; }
			if (tx[0] > tx[2]) { temp = tx[0]; tx[0] = tx[2]; tx[2] = temp; }
			if (tx[1] > tx[2]) { temp = tx[1]; tx[1] = tx[2]; tx[2] = temp; }
			if (tdx[0] > tdx[1]) { temp = tdx[0]; tdx[0] = tdx[1]; tdx[1] = temp; }
			if (tdx[0] > tdx[2]) { temp = tdx[0]; tdx[0] = tdx[2]; tdx[2] = temp; }
			if (tdx[1] > tdx[2]) { temp = tdx[1]; tdx[1] = tdx[2]; tdx[2] = temp; }
			BF64 va, vb; va.bf.u64 = vb.bf.u64 = 0;
			va.bf.u16[2] = tx[0] >> 8; va.bf.u16[1] = tx[1] >> 8; va.bf.u16[0] = tx[2] >> 8;
			vb.bf.u16[2] = tdx[0] >> 8; vb.bf.u16[1] = tdx[1] >> 8; vb.bf.u16[0] = tdx[2] >> 8;
			if (sgo.vx[6]) {// check bands 2 3 for test
				cout << "va1=" << va.bf.u16[2] << " va2=" << va.bf.u16[1]
					<< " va3=" << va.bf.u16[0] << "\tvb1=" << vb.bf.u16[2] << " vb2=" << vb.bf.u16[1]
					<< " vb3=" << vb.bf.u16[0];
				if (va.bf.u64 <= vb.bf.u64)			cout << "\tdirect"<<endl;
				else {
					cout << "\tdiag"<<endl;
				}
			}
			strcpy(zew, ze);// to morph the expected 17 puzzle
			if (va.bf.u64 > vb.bf.u64)	{
				va = vb;
				memcpy(tx, tdx, sizeof tx);
				strcpy(zout, zdiag);
				for (int i = 0; i < 81; i++)zew[C_transpose_d[i]] = ze[i];
			}
			if (sgo.vx[6]) {// check bands 2 3 for test
				if (va.bf.u16[2] != sgo.vx[6]) continue;
				if (va.bf.u16[1] != sgo.vx[7]) continue;
			}
			char zs[82], zsa[82]; zsa[81] = zs[81] = 0;
			int i1 = (tx[0] & 070) >> 3, i2 = (tx[1] & 070) >> 3, i3 = (tx[2] & 070) >> 3;
			memcpy(zsa, &zout[27 * i1], 27);
			memcpy(&zsa[27], &zout[27 * i2], 27);
			memcpy(&zsa[54], &zout[27 * i3], 27);
			memcpy(zes, &zew[27 * i1], 27);
			memcpy(&zes[27], &zew[27 * i2], 27);
			memcpy(&zes[54], &zew[27 * i3], 27);
			int zs0[81];
			// redo id to build tables
			int ib1, ib1a, ib2, ib2a, ib3, ib3a;
			for (int i = 0; i < 81; i++)zs0[i] = zsa[i] - '1';
			BANDMINLEX::PERM pb1, prw;
			bandminlex.Getmin(zs0, &pb1);
			ib1 = pb1.i416;
			ib1a = t416_to_n6[ib1];
			//myband1.InitBand2_3(ib1, zsa, perm_ret, 0);
			bandminlex.Getmin(&zs0[27], &prw);
			ib2 = prw.i416;
			ib2a = t416_to_n6[ib2];
			//myband2.InitBand2_3(ib2, &zsa[27], perm_ret, 1);
			bandminlex.Getmin(&zs0[54], &prw);
			ib3 = prw.i416;
			ib3a = t416_to_n6[ib3];
			strcpy(zs, "12345678945");
			strncpy(&zs[11], t416[ib1], 16);
			strcpy(zesf, empty_puzzle);
			// relabel all
			for (int i = 3, ij = 27; i < 9; i++)
				for (int j = 0; j < 9; j++, ij++) {
					int is = 9 * i + pb1.cols[j];
					int c = zsa[is] - '1';
					zs[ij] = pb1.map[c] + '1';
					if (zes[is] != '.')zesf[ij] = zs[ij];
				}
			if (sgo.vx[6]) {
				cout << "\ti1,2,3 " << i1 << i2 << i3 << endl;
				cout << ze << " ze" << endl;
				cout << zes << " zes" << endl;

				cout << zout << endl;
				cout << zsa << endl;
				//genb12.bands3[0].InitBand3(ib3, &zsa[54], perm_ret);
				
			}
			Msp_ReorderBand(&zs[27], &zesf[27]);
			Msp_ReorderBand(&zs[54], &zesf[54]);

			// reorder band1 p17
			for (int i = 0, ij = 0; i < 3; i++)
				for (int j = 0; j < 9; j++, ij++) {
					int is = 9 * pb1.rows[i] + pb1.cols[j],
						c = zes[is];
					if (c != '.')zesf[ij] = zs[ij];
				}
			if (sgo.vx[6]) {
				cout << "ib1a=" << ib1a << " ib2a=" << ib2a << " ib3a=" << ib3a << endl;
				cout << zs << "zsa morphed relabeled" << endl;
				cout << zs << "zsa reordered" << endl;
				cout << zesf << "zes reordered relabelled" << endl;

			}
			fout1 << zs << ";" << zesf << ";" << indexmin << ";";
			fout1  << va.bf.u16[1] << ";"
				<< va.bf.u16[0] << ";" << (tx[2] & 07)  << endl;


		}
		else cout << ze << "non resolu" << endl;
	}
}


/*

void GEN_BANDES_12::Morph_ED_B12_known17(char * z, int ib1){
	i1t16 = ib1;  it16 = tn6_to_416[ib1];
	char zw[164], *zb1p = &z[82], *zb2p = &z[27 + 82], *zb3p = &z[54 + 82];
	strcpy(zw, z);
	unsigned long  g0[81], *z02 = &g0[27], *z03 = &g0[54],zcomp[27];
	for (int i = 0; i < 81; i++)g0[i] = z[i] - '1';
	__movsd(zcomp,z02,27);
	n_auto_b1 = bandminlex.GetAutoMorphs(it16, t_auto_b1);
	cout << "processing band rank=" << i1t16 << " N° 1_416=" << it16 + 1
		<< " n auto morphs=" << n_auto_b1 << endl;
	if (n_auto_b1){
		for (int imorph = 0; imorph < n_auto_b1; imorph++){
			BANDMINLEX::PERM &p = t_auto_b1[imorph];
			unsigned long band[27];// morph the band
			for (int i = 0; i < 9; i++){
				band[i] = p.map[z02[p.cols[i]]];
				band[i + 9] = p.map[z02[p.cols[i] + 9]];
				band[i + 18] = p.map[z02[p.cols[i] + 18]];
			}
			if( G17ComparedOrderedBand((int *)zcomp,(int *) band)==1){//morph to p orderd
				int tsort[3];
				G17BuildSort((int *)band, tsort);
				cout << p.rows[0] << p.rows[1] << p.rows[2] << " rows"<<endl;
				for (int i = 0; i < 3; i++){
					int ir = tsort[i] &= 3,irb1=p.rows[i];
					__movsd(&zcomp[9*i], &band[9 * ir], 9);// new zs0
					for (int j = 0; j < 9; j++){// new entry
						int wc = band[9 * ir + j]+ '1',wc2=wc;
						if (zb2p[ 9 * ir + p.cols[j]] == '.') wc2 = '.';
						zw[27 + 9 * i + j] = wc;
						zw[27+82 + 9 * i + j] = wc2;
						wc = p.map[z03[9 * i + p.cols[j]] ] + '1';
						wc2 = wc;
						if (zb3p[ 9 * i + p.cols[j]] == '.') wc2 = '.';
						zw[54 + 9 * i + j] = wc;
						zw[54 + 82 + 9 * i + j] = wc2;
						wc2 = p.map[g0[9 * irb1 + p.cols[j]]] + '1';
						if (z[82+9 * irb1 + p.cols[j]] == '.') wc2 = '.';
						zw[82 + 9 * i + j] = wc2;

					}
				}
				cout << zw << "après morph une etape" << endl;
			}
		}
	}
	cout << zw << "après morph ED band 12" << endl;
}
*/

void msp_change12(char * ze, int &ib1a, int & ib2a) {
	for (int i = 0; i < 27; i++) {
		char temp = ze[i];	ze[i] = ze[i + 27];	ze[i + 27] = temp;
		temp = ze[i + 82];	ze[i + 82] = ze[i + 109]; ze[i + 109] = temp;
		int w = ib1a;	ib1a = ib2a;	ib2a = w;
	}
}
void msp_change13(char * ze, int &ib1a, int & ib3a) {
	for (int i = 0; i < 27; i++) {
		char temp = ze[i];	ze[i] = ze[i + 54];		ze[i + 54] = temp;
		temp = ze[i + 82];	ze[i + 82] = ze[i + 136];	ze[i + 136] = temp;
		int w = ib1a;	ib1a = ib3a;	ib3a = w;
	}
}
void msp_change23(char * ze, int &ib2a, int & ib3a) {
	msp_change12(&ze[27], ib2a, ib3a);
}
void M1_S17_Morph(char * z, BANDMINLEX::PERM & p1) {
	char zw[164]; strcpy(zw, z); // pose 0 et ;
	// map band 1
	for (int ir = 0, i = 0; ir < 3; ir++)for (int ic = 0; ic < 9; ic++, i++) {
		int iw = 9 * p1.rows[ir] + p1.cols[ic];
		char cw1 = p1.map[z[iw] - '1'] + '1', cw2 = cw1;
		if (z[iw + 82] == '.') cw2 = '.';
		zw[i] = cw1; zw[i + 82] = cw2;
	}
	{// map band2 to band 1 reorder columns
		int iw4 = 27 + p1.cols[0], iw5 = 36 + p1.cols[0], iw6 = 45 + p1.cols[0];
		int cw4 = p1.map[z[iw4] - '1'], cw5 = p1.map[z[iw5] - '1'], cw6 = p1.map[z[iw6] - '1'];
		int tsort[3], w;// sort bands / stacks increasing order of the  id
		tsort[0] = cw4 << 8;
		tsort[1] = 1 | (cw5 << 8);
		tsort[2] = 2 | (cw6 << 8);
		for (int i = 0; i < 2; i++) for (int j = i + 1; j < 3; j++)
			if (tsort[i] > tsort[j]) { w = tsort[i]; tsort[i] = tsort[j]; tsort[j] = w; }
		tsort[0] &= 3;  tsort[1] &= 3; tsort[2] &= 3;
		for (int ir = 0, i = 27; ir < 3; ir++)for (int ic = 0; ic < 9; ic++, i++) {
			int iw = 9 * tsort[ir] + p1.cols[ic] + 27;
			char cw1 = p1.map[z[iw] - '1'] + '1', cw2 = cw1;
			if (z[iw + 82] == '.') cw2 = '.';
			zw[i] = cw1; zw[i + 82] = cw2;
		}

	}

	{// map band3 to band 1 reorder columns
		int iw5 = 54 + p1.cols[0], iw6 = 63 + p1.cols[0], iw7 = 72 + p1.cols[0];
		int cw5 = p1.map[z[iw5] - '1'], cw6 = p1.map[z[iw6] - '1'], cw7 = p1.map[z[iw7] - '1'];
		int tsort[3], w;// sort bands / stacks increasing order of the  id
		tsort[0] = cw5 << 8;
		tsort[1] = 1 | (cw6 << 8);
		tsort[2] = 2 | (cw7 << 8);
		for (int i = 0; i < 2; i++) for (int j = i + 1; j < 3; j++)
			if (tsort[i] > tsort[j]) { w = tsort[i]; tsort[i] = tsort[j]; tsort[j] = w; }
		tsort[0] &= 3;  tsort[1] &= 3; tsort[2] &= 3;
		for (int ir = 0, i = 54; ir < 3; ir++)for (int ic = 0; ic < 9; ic++, i++) {
			int iw = 9 * tsort[ir] + p1.cols[ic] + 54;
			char cw1 = p1.map[z[iw] - '1'] + '1', cw2 = cw1;
			if (z[iw + 82] == '.') cw2 = '.';
			zw[i] = cw1; zw[i + 82] = cw2;
		}

	}
	strcpy(z, zw);
	//zw[54] = 0; zw[54 + 82] = 0;
	//fout1 << zw  << endl;
}
void M1_S17(char * finput_name, char * foutput_name, uint32_t * vx, uint32_t * bx) {
	// temporary code extract solution with a known 17 6 6 5
	finput.open(finput_name);
	if (!finput.is_open()) { cerr << "error open file " << finput_name << endl; return; }
	fout1.open(foutput_name);
	char * ze = finput.ze, *ze2 = &ze[82];
	char zdiag[164], zw[164];
	zdiag[81] = ';'; zdiag[163] = 0;
	int count[6], zwi[81], zei[81], zdiagi[81];
	while (finput.GetLigne()) {
		if (strlen(ze) < 163)continue;
		ze[163] = 0; // forget more info
		memset(count, 0, sizeof count);
		for (int i = 0; i < 81; i++) {
			if (ze2[i] - '.') {
				int band = i / 27, stack = C_stack[i] + 3;
				++count[band];
				++count[stack];
			}
			zdiag[C_transpose_d[i]] = ze[i];
			zdiag[C_transpose_d[i] + 82] = ze[i + 82];
		}
		BANDMINLEX::PERM pr[6];
		for (int ib = 0; ib < 6; ib++)if (count[ib] > 6) {
			if (0) goto next;
			int x;
			// morph the puzzle to the high band in band
			if (ib > 2)memmove(zw, zdiag, sizeof zw);
			else memmove(zw, ze, sizeof zw);
			// put the >6 clues band as band3
			if (ib == 0 || ib == 3) msp_change13(zw, x, x);
			else if (ib == 1 || ib == 4) msp_change23(zw, x, x);
			for (int i = 0; i < 54; i++) zwi[i] = zw[i] - '1';
			bandminlex.Getmin(zwi, &pr[0]);
			int ib1 = pr[0].i416, ib1a = t416_to_n6[ib1];
			bandminlex.Getmin(&zwi[27], &pr[1]);
			int ib2 = pr[1].i416, ib2a = t416_to_n6[ib2],
				i2 = ib2a, i1 = ib1a;
			if (ib1a > ib2a) {
				msp_change12(zw, i1, i2);
				M1_S17_Morph(zw, pr[1]);
			}
			else 	M1_S17_Morph(zw, pr[0]);

			fout1 << zw << ";" << i1 << ";" << i2 << endl;
			cout << zw << endl;
			//genb12.Morph_ED_B12_known17(zw, i1);
			if (0) {// mode look for 566 656 no 7 in stack
				char zout[164];
				for (int i = 0; i < 81; i++) {
					zei[i] = ze[i] - '1';
					zdiagi[i] = zdiag[i] - '1';
				}
				BANDMINLEX::PERM perm_reti[6];
				int ibi[6], ibri[6];
				for (int i = 0; i < 3; i++) {
					bandminlex.Getmin(&zei[27 * i], &perm_reti[i]);
					ibi[i] = perm_reti[i].i416;
					ibri[i] = t416_to_n6[ibi[i]];
					bandminlex.Getmin(&zdiagi[27 * i], &perm_reti[i + 3]);
					ibi[i + 3] = perm_reti[i + 3].i416;
					ibri[i + 3] = t416_to_n6[ibi[i + 3]];
				}
				int tsort[3], w;// sort bands / stacks increasing order of the  id
				tsort[0] = ibri[0] << 8;
				tsort[1] = 1 | (ibri[1] << 8);
				tsort[2] = 2 | (ibri[2] << 8);
				for (int i = 0; i < 2; i++) for (int j = i + 1; j < 3; j++)
					if (tsort[i] > tsort[j]) { w = tsort[i]; tsort[i] = tsort[j]; tsort[j] = w; }
				int ib1 = tsort[0] & 3, ib2 = tsort[1] & 3, ib3 = tsort[2] & 3;
				tsort[0] = 3 | ibri[3] << 8;
				tsort[1] = 4 | (ibri[4] << 8);
				tsort[2] = 5 | (ibri[5] << 8);
				for (int i = 0; i < 2; i++) for (int j = i + 1; j < 3; j++)
					if (tsort[i] > tsort[j]) { w = tsort[i]; tsort[i] = tsort[j]; tsort[j] = w; }
				int ibd1 = tsort[0] & 3, ibd2 = tsort[1] & 3, ibd3 = tsort[2] & 3;
				char *z = ze;
				int * cpt = count, iout = ibri[ib1];
				if (ibri[ib1] < ibri[ibd1]) goto bands;
				if (ibri[ib1] > ibri[ibd1]) goto stacks;
				if (ibri[ib2] < ibri[ibd2]) goto bands;
				if (ibri[ib2] > ibri[ibd2]) goto stacks;
				goto bands;// not a key point, use band
			stacks:
				z = zdiag;
				cpt = &count[3];
				ib1 = ibd1 - 3; ib2 = ibd2 - 3; ib3 = ibd3 - 3;
				iout = ibri[ibd1];
			bands: {
				strcpy(zout, z);
				if (ib1) {
					memcpy(zout, &ze[27 * ib1], 27);
					memcpy(&zout[82], &ze[27 * ib1 + 82], 27);
				}
				if (count[ib3] == 5) {//  must be pass2b
					if (ib3 != 2) {
						memcpy(&zout[27], &z[27 * ib3], 27);
						memcpy(&zout[82 + 27], &z[27 * ib3 + 82], 27);
					}
					if (ib2 != 3) {
						memcpy(&zout[54], &z[27 * ib2], 27);
						memcpy(&zout[82 + 54], &z[27 * ib2 + 82], 27);
					}
					fout1 << z << ";" << iout << endl;
				}
				else {//  pass2a
					if (ib2 != 2) {
						memcpy(&zout[27], &z[27 * ib2], 27);
						memcpy(&zout[82 + 27], &z[27 * ib2 + 82], 27);
					}
					if (ib3 != 3) {
						memcpy(&zout[54], &z[27 * ib3], 27);
						memcpy(&zout[82 + 54], &z[27 * ib3 + 82], 27);
					}
					cout << z << ";" << iout << endl;

				}
				}
			}
		}
	next:;
	}
}


