
//========================================
const char * zh_g_cpt[10] = { "npuz", "guess", "close_d ", "upd1 ", "upd2 ",
"fupd ", "hpair ", "htripl ", " ", " " };

const char * libs_c17_00_cpt2g[60] = {
	"0 bands1+2 processed entries M10",//0
	"1 total bands 3",//1
	"2 steps external loop ",//2
	"3 loop band1 ",//3
	"4 loop band 2  ",//4
	"5 clean entries    ",//5
	"6 pass filters base",//6
	"7 calls brute force ",//7
	"8 valid brute force",//8
	"9 3B nb3 tot miss",//9
	"10 miss 0",//10
	"11 miss 1",//11
	"12 miss 2",//12
	"13 expand",//13
	"14 miss 2 expand",//14
	"15 3A nb ok",//15
	"16 maxuas",//16
	"17 maxuas b1",//17
	"18 switch mv2",//18
	"19 ",//19
	"20 ",//20
	"21 pass mincount",//21
	"22 3a mincount <=6",//22
	"23 check b12 first",//23
	"24 check b12 first ok  ",//24
	"25 ",
	"26 chunk count ",
	"27 clean count ",	
	"28  need to assign", 
	"29  b3 other call check",
	"30  b3_expand call check",	
	"31 ajouts uas base",
	"32 ajouts uas base+1 ",
	"33  ajouts uas base+2",
	"34  ajouts uas base+3",
	"35 max ngua stepb2",	
	"36 ng2_2 + ng3_2",	
	"37 setv nguan",	
	"38 setv nguas", 
	"39 all steps",
	"40 b3 new size 2",
	"41 b3 new size 3",
	"42 b3 new size>3",
	"43 ",
	"44 ",
	"45 step dead  b2",
	"46 ",
	"47",
	"48 n clean b",
	"49 matrix ",
	"50 count clean",
	"51  check g2 ok",
	"52  check g2 ",
	"53  ",
	"54  ",
	"55 zhou[1].CallMultipleB3",
	"56 nmiss1 ua of found",
	"57 nmiss1 go subcritical",
	"58 b3 expand critical 2 pairs ",
	"59 test",
};
void Go_c17_00( ) {// p2 process
	cout << "Go_c17_00 search batch 17 clues 656 566 " << endl;
	cout << sgo.vx[0] << " -v0- band 0_415" << endl;
	cout << sgo.vx[2] << " -v2- skip first nnn restart after batch failure" << endl;
	cout << sgo.vx[3] << " -v3- last entry number for this batch must be > vx[2]" << endl;
	cout << sgo.vx[4] << " -v4- 0 if p2a 1 if p2b" << endl;
	int it16_start = sgo.vx[0];
	g17b.aigstop=0;
	//g17b.diag = sgo.vx[6];
	genb12.skip = sgo.vx[2];
	genb12.last = sgo.vx[3];
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
	memset(p_cpt2g, 0, sizeof p_cpt2g);// used in stats
	genb12.Start(0);
	genb12.NewBand1(sgo.vx[0]);
	cout << "print final stats" << endl;
	for (int i = 0; i < 60; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
	cout << "exit final stats" << endl;
}
//========================= known s17 file 10/19
void Go_c17_10() {
	zh_g.modevalid = 1;
	zh_g2.grid0 = genb12.grid0;
	zh_g2.zsol = zh_g2.stdfirstsol;
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	char * ze = finput.ze;
	int * zs0 = genb12.grid0, npuz = 0;
	cout << "Go_c17_10() search 17 using a file having known 17 656 " << endl;
	while (finput.GetLigne()) {
		if (strlen(ze) < 160) continue;// skip blank lines
		npuz++;
		g17b.npuz = npuz;
		g17b.a_17_found_here = 0;
		g17b.aigstop = 0;
		if (npuz <= (int)sgo.vx[2]) continue;
		if (npuz > (int)sgo.vx[3]) break;
		//if (npuz >5) return;
		cout << "\n\nto process  n=" << dec << npuz  << endl;
		long tdeb = GetTimeMillis();
		//================================ to avoid the 665 case
		int ncb3 = 0;
		for (int i = 0; i < 27; i++)
			if (ze[i + 136] - '.')ncb3++;
		if (ncb3 == 5) {// change band3 <-> band2
			for (int i = 0; i < 27; i++) {
				char temp = ze[i + 27];	ze[i + 27] = ze[i + 54];	ze[i + 54] = temp;
				temp = ze[i + 109];	ze[i + 109] = ze[i + 136];	ze[i + 136] = temp;
			}
		}

		cout << ze << " to process  n=" << npuz << endl;

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
		ze[81] = 0;
		char * ze2 = &ze[82];
		g17b.p17diag.SetAll_0();
		for (int i = 0; i < 81; i++) if (ze2[i] != '.')
			g17b.p17diag.Set_c(i);
		cout << g17b.p17diag.bf.u64[0] << " b12 pattern for the 17" << endl;
		genb12.ValidInitGang();
		g17b.npuz = npuz;
		g17b.ExpandB1();// expand band1
		g17b.GoM10();
		//g17b.a_17_found_here = 1;
		if (!g17b.a_17_found_here) {
			cout << "puz=" << npuz << " failed to find the searched 17" << endl;
			cerr << "puz=" << npuz << " failed to find the searched 17" << endl;
			break;
		}
	}
	if (0) {
		cout << "print final stats" << endl;
		for (int i = 0; i < 60; i++) {
			if (!p_cpt2g[i])continue;
			cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
		}

	}
}
