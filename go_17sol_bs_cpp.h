

G17B::G17B() {
	memset(this, 0, sizeof(*this));
	gang27 = gang[0];
	aigstop = 0;
	// Init SGUAs tables load permanent data
	int tpcol[3][2] = { {1,2},{0,2},{0,1} };
	for (int ist = 0; ist < 3; ist++) {//stack
		for (int ircol = 0; ircol < 3; ircol++) {//first column
			int* p = tpcol[ircol], dcol = 3 * ist;
			for (int irdig = 0; irdig < 3; irdig++) {//digit1
				for (int irdig2 = 0; irdig2 < 3; irdig2++) {//digit 2
					int i = 27 * ist + 9 * ircol + 3 * irdig + irdig2;
					SGUA2& w = tsgua2[i];
					w.i_81 = i; w.i9= dcol + ircol;
					w.col1 = dcol + p[0];		w.col2 = dcol + p[1];
					w.id1 = irdig + 3 * w.col1;
					w.id2 = irdig2 + 3 * w.col2;
					if (0) {
						cout << "sua2=" << w.i_81 //<< " i9=" << w.i9 + 1
							<< " cols12 " << w.col1 + 1 << w.col2 + 1
							<< " id1;id2 " << w.id1 << ";" << w.id2 << endl;
					}
				}
			}
		}
	}
	for (int ist = 0; ist < 3; ist++) {//stack
		for (int irdig1 = 0; irdig1 < 3; irdig1++) {//digit1
			for (int irdig2 = 0; irdig2 < 3; irdig2++) {//digit 2
				for (int irdig3 = 0; irdig3 < 3; irdig3++) {//digit 2
					int i = 27 * ist + 9 * irdig1 + 3 * irdig2 + irdig3;
					SGUA3& w = tsgua3[i];
					w.i_81 = i;
					w.stack = ist;
					w.col1 = 3 * ist;// minirow first column in gangster
					// pointers to gangster digits
					w.id1 = 9 * ist+irdig1;
					w.id2 = 9 * ist + 3+irdig2;
					w.id3 = 9 * ist + 6 + irdig3;
					if (0) {
						cout << "sgua3 " << w.i_81 << " col1=" << w.col1 + 1
							<< " id1;id2,id3 " << w.id1
							<< ";" << w.id2 << ";" << w.id3 << endl;
					}
				}
			}
		}
	}
}
void G17B::StartInit() {
	memcpy(grid0, genb12.grid0, sizeof grid0);// use first b3
	memcpy(&grid0[54], genb12.bands3[0].band0, sizeof genb12.bands3[0].band0);
	zh2b_g.Init_g0(grid0);// init brute force bands 1+2
	if (op.p1)mincluesb3 = 7;
	else mincluesb3 = 6;
	knownt = 0;
	for (int i = 0; i < 9; i++) {// 9 cols to build out of gangcols
		int istack = C_box[i];
		int* d = gang[i], c = genb12.gangcols[i];
		uint32_t bit;
		bitscanforward(bit, c);
		c ^= (1 << bit);
		d[0] = bit;
		gang_digits_cols[bit][istack] = i;
		bitscanforward(bit, c);
		c ^= (1 << bit);
		d[1] = bit;
		gang_digits_cols[bit][istack] = i;
		bitscanforward(bit, c);
		d[2] = bit;
		gang_digits_cols[bit][istack] = i;
	}

	//finish sguas2/3 set up with the gangster 
	gsock2.SetAll_0();  gsock3.SetAll_0();
	//register int* gang27=genb12.gang27; // redefines gang[9][3] as 27 integer
	//for (int i = 0; i < 27; i++) cout << gang27[i] + 1;
	//cout << endl;
	int tpcol[3][2] = { {1,2},{0,2},{0,1} };
	for (register int ist = 0,i=0,stack=07007007; ist < 3; ist++,stack<<=3) {
		//cout << Char27out(stack) << endl;
		for (register int irst = 0; irst < 27; irst++, i++) {//stack
			{//________________________ guas2
				SGUA2& w = tsgua2[i];
				w.dig1 = gang27[w.id1];
				w.dig2 = gang27[w.id2];
				w.digs = (1 << w.dig1) | (1 << w.dig2);
				w.valid = 0;
				//cout << i << " " << w.dig1+1 <<" "<< w.dig2+1 << endl;
				//continue;
				for (register int ib3 = 0; ib3 < genb12.nband3; ib3++) {
					STD_B3& b3 = genb12.bands3[ib3];
					register uint32_t Bf = b3.dpat[w.dig1] | b3.dpat[w.dig2],
						Bfc = (Bf | (Bf >> 9) | (Bf >> 18)) & 0777;// colummn pattern
					Bf &= stack;
					int nr = 0, bfr = 0, irow;
					if (Bf & 0777) { nr++; bfr |= Bf & 0777; irow = 0; }// row1
					if (Bf & 0777000) { nr++; bfr |= Bf & 0777000; irow = 1; }
					if (Bf & 0777000000) { nr++; bfr |= Bf & 0777000000; irow = 2; }
					if (nr == 1) {// we have a gua2
						//cout << ib3 << " " << i << endl;
						gsock2.setBit(i);
						w.valid = 1;
						b3.g.gsocket2.setBit(i);
						b3.g.pat2[i] = Bf;
						int imini = 3 * irow + ist, mask = 7 << (3 * imini);
						b3.g.ua2_imini[i] = imini;
						int bit27= mask ^ Bf,i27;
						bitscanforward(i27, bit27);
						b3.g.ua2bit27[i] = bit27;
						b3.i_27_to_81[i27] = i;
						b3.i_81_to_27[i] = i27;
					}
					else {// gua46 except if 6 columns
						int ncol = _popcnt32(Bfc);
						if (ncol == 5) {// gua2_4
							b3.g.pat2[i] = bfr;
						}
						else b3.g.pat2[i] = Bf;
					}
				}
			}
			
			
			{ //________________________ guas3

				SGUA3& w = tsgua3[i];
				w.dig1 = gang27[w.id1];
				w.dig2 = gang27[w.id2];
				w.dig3 = gang27[w.id3];
				w.digs = 1 << w.dig1 | 1 << w.dig2 | 1 << w.dig3;
				w.valid = 0;
				for (register int ib3 = 0; ib3 < genb12.nband3; ib3++) {
					STD_B3& b3 = genb12.bands3[ib3];
					register uint32_t Bf = b3.dpat[w.dig1]
						| b3.dpat[w.dig2] | b3.dpat[w.dig3];// colummn pattern
					Bf &= stack;
					int nr = 0,  irow=0;
					if (Bf & 0777) { nr++;  irow =0; }// row1
					if (Bf & 0777000) { nr++;  irow = 1; }
					if (Bf & 0777000000) { nr++;  irow = 2; }
					if (nr == 1) {// we have a gua3
						gsock3.setBit(i);
						w.valid = 1;
						b3.g.gsocket3.setBit(i);
						b3.g.pat3[i] = Bf;
						int imini = 3 * irow + ist;
						b3.g.ua3_imini[i] = imini;
						b3.i_9_to_81[imini] = i;

					}
				}
			}
		}
	}
	if (op.ton > 1|| (sgo.bfx[3] & 1)) genb12.bands3[0].DumpIndex();
}
void G17B::StartInitDebug() {
	cout << "check guas2 status" << endl;
	for (register int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& b3 = genb12.bands3[ib3];
		cout << Char64out(b3.g.gsocket2.bf.u64[0]);
		cout << Char27out(b3.g.gsocket2.bf.u32[2])<<" "
			<< b3.g.gsocket2 .Count96() << endl;
	}
	cout << Char64out(gsock2.bf.u64[0]);
	cout << Char27out(gsock2.bf.u32[2]) << " all guas2 "
		<< gsock2.Count96() << endl;

	cout << "check guas3 status" << endl;
	for (register int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& b3 = genb12.bands3[ib3];
		cout << Char64out(b3.g.gsocket3.bf.u64[0]);
		cout << Char27out(b3.g.gsocket3.bf.u32[2]) 
			<< b3.g.gsocket3.Count96() << endl;
	}
	cout << Char64out(gsock3.bf.u64[0]);
	cout << Char27out(gsock3.bf.u32[2]) << " all guas3 " 
		<< gsock3.Count96() << endl;

}
void G17B::StartPrint() {
	cout << myband2.band << " go band2 id=" << myband2.i416 
		<< "index="<< t416_to_n6[myband2.i416] << " nb12=" << genb12.nb12
		<< " nb3=" << genb12.nband3 << " p_cpt2g[0]=" << p_cpt2g[0] << endl;
	cout << Char64out(gsock2.bf.u64[0]);
	cout << Char27out(gsock2.bf.u32[2]) << " all guas2 "
		<< gsock2.Count96() << endl;
	cout << Char64out(gsock3.bf.u64[0]);
	cout << Char27out(gsock3.bf.u32[2]) << " all guas3 "
		<< gsock3.Count96() << endl;
	if (op.ton >1) {
		cout << myband2.band << endl << endl;
		for (int i = 0; i < genb12.nband3; i++) {
			STD_B3& b3 = genb12.bands3[i];
			cout << b3.band << " b3 i=" << i << endl;
			//for (int id = 0; id < 9; id++)
			//	cout << Char27out(b3.dpat[id]) << endl;
		}
	}
	for (int i = 0; i < 81; i++) cout << grid0[i] + 1;
	cout << " grid to use in ua collector" << endl;

	//StartInitDebug();
	//for (int i = 0; i < genb12.nband3; i++) {
		//STD_B3& b3 = genb12.bands3[i];
		//cout << b3.band << " b3 i=" << i << endl;
		//for (int i = 0; i < (int)b3.nua; i++)
		//	cout << Char27out(b3.tua[i]) << endl;
	//}

}

void G17B::Start() {// processing an entry 
	if (aigstop)return;
	p_cpt2g[0] ++;
	p_cpt2g[1] += genb12.nband3;

	if(sgo.vx[5])		cout <<myband2.band<<" [0] "<< p_cpt2g[0] << endl;
	//if (p_cpt2g[0] > 1) {		aigstop = 1; return;	}

	StartInit();// do all preliminary setups
	if(op.ton)	StartPrint();
	UaCollector();//_do uas/guas2_3 initial harvest
	StartAfterUasHarvest();
}

void G17B::StartKnown() {// processing an entry  with a known o follow3
	if (aigstop)return;
	p_cpt2g[0] ++;
	StartInit();// do all preliminary setups
	if(op.ton)StartPrint();
	UaCollector();//_do uas/guas2_3 initial harvest
	StartAfterUasHarvest();
}

struct EXTUAS {// get a sorted subset of uas for next "collect" step
	uint64_t  t2[300];
	uint32_t  nt2;
	void GetSortB(uint64_t* t, uint32_t n, uint64_t filter, uint32_t nold = 0) {
		BF128 vsize[25][4];
		uint64_t  t2a[512];
		uint32_t  nt2a = 0;
		nt2 = nold;
		register uint64_t Fn = (~filter) & BIT_SET_2X, w;
		memset(vsize, 0, sizeof vsize);
		for (uint32_t iua = 0; iua < tuasb12.nua; iua++) {
			w = tuasb12.tua[iua];
			uint64_t cc = w >> 59;
			w &= BIT_SET_2X;
			if (!(w & Fn)) {
				int bloc = nt2a >> 7, ir = nt2a - (bloc << 7);
				vsize[cc][bloc].setBit(ir);
				t2a[nt2a++] = w;// clean the count 
			}
			if (nt2a >= 512)break;
		}
		//cout << nt2a << " nt2a" << endl;
		uint32_t nbl64 = (nt2a + 63) >> 6, x;
		for (int i1 = 4; i1 < 25; i1++) {
			uint64_t* tb64 = vsize[i1]->bf.u64;
			for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
				register uint64_t V = tb64[i2];
				while (bitscanforward64(x, V)) {// switch to 54 mode here
					V ^= (uint64_t)1 << x;
					t2[nt2++] = t2a[x + (i2 << 6)];
					if (nt2 >= 100)return;

				}
			}
		}


	}
	void GetSortA(uint64_t filter, uint32_t nold = 0) {
		nt2 = nold;
		register uint64_t Fn = (~filter) & BIT_SET_2X, w;
		for (uint32_t iua = 0; iua < tuasb12.nua; iua++) {
			w = tuasb12.tua[iua];
			if (!(w & Fn))t2[nt2++] = w;
			if (nt2 >= 100)break;
		}
		if (nt2) {
			if (nt2 > 1)
				for (uint32_t i = 0; i < nt2 - 1; i++) {
					register uint64_t Ri = t2[i];
					for (uint32_t j = i + 1; j < nt2; j++) {
						register uint64_t  Rj = t2[j];
						if ((Ri >> 59) > (Rj >> 59)) { // smaller size first 
							t2[i] = Rj; t2[j] = Ri; Ri = Rj;
						}
					}
				}
			for (uint32_t i = 0; i < nt2; i++)
				t2[i] &= BIT_SET_2X;// clean the count

		}
	}
	void Dump() {
		cout << "dump extuas nt2=" << nt2 << endl;
		for (uint32_t i = 0; i < nt2; i++)
			cout << Char2Xout(t2[i]) << "i=" << i << endl;
	}
}extuas;

void G17B::UaCollector() {
	//for (int ib3 = 0; ib3 < genb12.nband3; ib3++) genb12.bands3[ib3].ntguam = 0;
	zh2gxn.SetupFsol(grid0);
	zh2b[0].InitBands12(grid0);
	zh2gxn.InitKnown(extuas.t2, &extuas.nt2);
	tuasb12.nua = 0;
	FirstUasCollect();// uas 2 digits bands and stacks
	//return;
//	aigstop = 1; return;
	SecondUasCollect();// uas 345 bands 1+2
	UasCollect6_7();
	Guas2Collect();
	UasCollect4box();// 4 box in bands 1+2
	zh2b[0].InitBands12(grid0);// restore 
	//________________________________ insert stack known uas 
	STD_B3 bw;
	BANDMINLEX::PERM pp;
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& mb3 = genb12.bands3[ib3];
		memcpy(&grid0[54], mb3.band0, sizeof genb12.bands3[0].band0);
		int grid_diag[81];
		char ze[81];
		for (int i = 0; i < 81; i++)grid_diag[C_transpose_d[i]] = grid0[i];
		for (int i = 0; i < 81; i++)  ze[i] = grid_diag[i] + '1';
		for (int istack = 0; istack < 3; istack++) {
			int ir = bandminlex.Getmin(&grid_diag[27 * istack], &pp, 0);
			if (ir < 0) {//would be bug} 
				cout << "look for stack bad return" << endl;
				return;
			}
			bw.InitBand2_3(pp.i416, &ze[27 * istack], pp, istack);
			//bw.PrintStatus();
			for (uint32_t i = 0; i < bw.nua; i++) {
				BF128 wd; wd.SetAll_0();// 3x mode to build guam
				uint32_t  zed = 0, ua = bw.tua[i]&BIT_SET_27,ccua=_popcnt32(ua);
				for (uint32_t j = 0; j < 27; j++) {
					if (ua & (1 << j)) {
						int cell = C_transpose_d[j + 27 * istack];
						wd.Set_c(cell);
					}
				}
				int cc3 =  _popcnt32  (wd.bf.u32[2]);
				//cout << Char27out(zed) << endl;
				//cout << Char2Xout(wd.bf.u64[0]) << " ";
				//cout << Char27out(wd.bf.u32[2])<<" "<<wd.Count()<<" "<<cc3 << endl;
				if (cc3 > 3)mb3.Addguam(wd);
				else if(cc3==2  &&  ccua>10) { // must build g2
					int i81=mb3.GetI81_2(wd.bf.u32[2]);
					guah.tg2[i81].AddIf(wd.bf.u64[0]);
				}
				else if (cc3 == 3  &&  ccua>9) { // must build g3
					int i81 = mb3.GetI81_3(wd.bf.u32[2]);
					guah.tg3[i81].AddIf(wd.bf.u64[0]);
				}

			}
		}
		//cout << "band summary ib3 =" << ib3 << "\t";
		//mb3.DumpGuam();


	}
	//guah.DumpOne2(41);
	//guah.Dump2all2();
}

inline void G17B::Adduab12( uint32_t digs, uint32_t nd) {

	uint64_t* t = zh2gxn.tua;
	uint32_t n = zh2gxn.nua,n2=0;
	BF128 tw[100];
	for (uint32_t i = 0; i < n; i++) {
		register uint64_t w = t[i], cc = _popcnt64(w);
		if (cc < 25) {	
			BF128 x; x.bf.u64[0] =  w; x.bf.u64[1] = cc;	
			tw[n2++] = x;
		}
	}
	if (n2 > 1) {// sort to check redundancy
		for (uint32_t i = 0; i < n2 - 1; i++) {
			for (uint32_t j=i+1; j < n2 ; j++) {
				if (tw[i].bf.u64[1] > tw[j].bf.u64[1]) {
					BF128 x = tw[i]; tw[i] = tw[j]; tw[j] = x;
				}
			}
		}
	}
	for (uint32_t i = 0; i < n2; i++) {
		BF128 x = tw[i];
		register uint64_t cc = x.bf.u64[1],
			w= x.bf.u64[0],nw=~w;

		if (i) {
			for (uint32_t j = 0; j < i; j++) {
				BF128 y = tw[j];
				register uint64_t ycc = y.bf.u64[1];
				if (ycc == cc) break;
				if(! (y.bf.u64[0] & nw)) { cc = 0; break; }
			}
			if (!cc)continue;// subset found
		}
		if (cc > 15) cc = 15;
		w |= (cc << 59);
		if (tuasb12.nua < UA12SIZE)tuasb12.AddInit(w, digs, nd);
	}

}
void G17B::FirstUasCollect() {// produce all uas 2/3 digits

	struct TUAS81 {// used to collect uas 2 digits
		BF128 	tall[200];// first collect
		uint32_t ntall;// , ntold;
		int Add(BF128 w, uint32_t floor) {
			uint32_t cc = w.Count96(), nfloor = ~floor;
			for (uint32_t iua = 0; iua < ntall; iua++) {
				BF128 wt = tall[iua];
				uint32_t cct = wt.bf.u32[3] >> 16;
				uint32_t floort = wt.bf.u32[3] & 0777;
				wt.bf.u32[3] = 0;
				if (cct <= cc)continue; 
				// insert here and check super set
				for (uint32_t jua = ntall; jua > iua; jua--)
					tall[jua] = tall[jua - 1];
				tall[iua] = w;// new inserted
				tall[iua].bf.u32[3] = floor | (cc << 16);
				ntall++;
				return 2;
			}
			w.bf.u32[3] = floor | (cc << 16);
			tall[ntall++] = w;// new added
			return 1;
		}
	}tuas81;
#ifdef TEST_ON
	cout << "entry FirstUasCollect()" << endl;
#endif
	for (int i = 0; i < 36; i++) {
		uint32_t myf = floors_2d[i];
		zh2_2[0].GoZ2D(myf);
		//continue;
		if (zh2gxn.nua) {
			//cout << "found UAs for i=" << i << endl;
			for (uint32_t i = 0; i < zh2gxn.nua; i++) {
				register uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);
				//cout << Char2Xout(w) << " " << cc << endl;
			if (cc > 15)cc = 15;
			w |= cc << 59;
			tuasb12.AddInit(w, myf, 2);
			}
		}

	}

	// do zhou2[0].GoZ2(myf) for al bands looing for guas
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& mb3 = genb12.bands3[ib3];
		memcpy(&grid0[54], mb3.band0, sizeof genb12.bands3[0].band0);
		zhgxn.SetupFsol(grid0);
		//cout << "ib3 =" << ib3 << endl;
		for (int i = 0; i < 36; i++) {
			uint32_t myf = floors_2d[i];
			zhou2[0].GoZ2(myf);
			//cout << "back nua=" << zhgxn.nua << endl;
			if (zhgxn.nua == 1) {// this is a 18 cells 12+6
				BF128 w = zhgxn.tua[0];
				register uint32_t C = w.bf.u32[2];
				uint32_t  ncol = _popcnt32((C | (C >> 9) | (C >> 18)) & 0777);
				if (ncol != 4)mb3.Addguam(w);
				else {// this is a 6 clues gua2
					//cout << Char2Xout(w.bf.u64[0]) << "socket 2 12+6" << endl;
					//for (int i = 0; i < 3; i++) {
					//	cout << Char9out(C) << endl;
					//	C >>= 9;
					//}
				}
			}
			else {
				for (uint32_t j = 0; j < zhgxn.nua; j++) {
					BF128 w = zhgxn.tua[j];
					int cc = w.Count();
					int aig = 1;// check for subsets
					if (cc == 18)continue; // has a subset
					for (uint32_t k = 0; k < zhgxn.nua; k++) if (k != j) {
						BF128 w2 = zhgxn.tua[k];
						if (w2.isSubsetOf(w)) { aig = 0; break; }
					}
					if (aig) {
						uint64_t w1 = w.bf.u64[0], cc = _popcnt64(w1);
						register uint32_t C = w.bf.u32[2];
						uint32_t cc2 = _popcnt32(C),
							ncol = _popcnt32((C | (C >> 9) | (C >> 18)) & 0777),
							nplug = 2 * ncol - cc2;
						if ((!cc) || (!cc2)) continue;//ua bands 1+2 or ua band 3
						if (nplug > 2) { mb3.Addguam(w); continue; }
						//cout << Char2Xout(w.bf.u64[0]) << "socket 2 "<< cc<<" + "<<cc2 << endl;
						//for (int i = 0; i < 3; i++) {
						//	cout << Char9out(C) << endl;
						//	C >>= 9;
						//}
					}
				}
			}		
		}
		//mb3.DumpGuam2d();
	}



}
void G17B::SecondUasCollect() {// collect 345 digits in bands 1+2
	//cout << "seconuacollect" << endl;
	//zh2gxn.InitKnown(tuasb12.t2, &tuasb12.nt2);
	//cout << "3 digits" << endl;
	for (int i = 0; i < 84; i++) {// look for 3 digits
		uint32_t myf = floors_3d[i];
		int ir = zh2_3[0].GoZ3(myf);// cells unsolved count
		if (ir < 6) continue;// minimum for a fresh ua 3 digits
		uint64_t F = zh2_3[0].cells_unsolved.bf.u64;
		extuas.GetSortA(F);
		zh2_3[0].DoZ3Go();
		if (zh2gxn.nua) Adduab12(myf,3);
	}
	if (op.ton)tuasb12.DumpShort("3 digits");
	for (int i = 0; i < 126; i++) {// look for 4 digits
		int ir = zh2_4[0].GoZ4(floors_4d[i]);
		if (ir < 8) continue;// minimum for a fresh ua 4 digits
		uint64_t F = zh2gxn.unsolved_field;
		extuas.GetSortA(F);
		zh2_4[0].DoZ4Go();
		if (zh2gxn.nua) Adduab12(floors_4d[i],4);
	}
	if (op.ton)tuasb12.DumpShort("4 digits");
	for (int i = 0; i < 126; i++) {// look for 5 digits
		uint32_t myf = 0777 ^ floors_4d[i];
		int ir = zh2_5[0].GoZ5(myf);
		uint64_t F = zh2gxn.unsolved_field;
		extuas.GetSortB(tuasb12.tua, tuasb12.nua, F);
		zh2_5[0].DoZ5Go();
		if (zh2gxn.nua) Adduab12(myf, 5);
	}
	if (op.ton)tuasb12.DumpShort("5 digits");
}
void G17B::UasCollect6_7() {
	for (int i = 0; i < 84; i++) {
		uint32_t ass= floors_3d[i], myf = 0777^ass;
		uint64_t bf = zh2gxn.Getsol(ass);
		uint64_t F = bf^BIT_SET_2X;
		//cout << Char2Xout(bf) << " ";cout << Char9out(myf) << endl;
		extuas.GetSortB(tuasb12.tua, tuasb12.nua, F);
		//extuas.GetSort(F);
		zh2b[0].InitBf(bf);
		zh2b[0].Dob();
		if (zh2gxn.nua) Adduab12(myf, 6);
	}
	if (op.ton)tuasb12.DumpShort("6 digits");
	return;// 7 seems too expensive 
	for (int i = 0; i < 36; i++) {
		uint32_t ass = floors_2d[i], myf = 0777 ^ ass;
		uint64_t bf = zh2gxn.Getsol(ass);
		uint64_t F = bf ^ BIT_SET_2X;
		extuas.GetSortB(tuasb12.tua, tuasb12.nua,F);
		zh2b[0].InitBf(bf);
		zh2b[0].Dob();
		if (zh2gxn.nua) Adduab12(myf, 7);
	}
}
void G17B::UasCollect4box() {
	//cout << "4 box collect" << endl;
	//___________________ box 1245
	{
		uint64_t stack = units3xBM[5].u64[0], b4 = BIT_SET_2X & (~stack);
		extuas.GetSortB(tuasb12.tua, tuasb12.nua, b4);
		//extuas.Dump();
		zh2b[0].InitBf(stack);

	}
	zh2b[0].Dob(22);
	if (zh2gxn.nua) {
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);	
			//cout << Char2Xout(w) << " to add b4 size " << cc << endl;
			if (tuasb12.nua < 2550)tuasb12.AddInit(w | (cc << 59), 0,0);
		}
	}
	
	//___________________ box 1346
	{
		uint64_t stack = units3xBM[4].u64[0], b4 = BIT_SET_2X & (~stack);
		extuas.GetSortB(tuasb12.tua, tuasb12.nua,b4);
		//extuas.Dump();
		zh2b[0].InitBf(stack);

	}
	zh2b[0].Dob(22);
	if (zh2gxn.nua) {
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);
			//cout << Char2Xout(w) << " to add b4 size " << cc << endl;
			if (tuasb12.nua < UA12SIZE)tuasb12.AddInit(w | (cc << 59), 0,0);
		}
	}
	//___________________ box 2356
		{
		uint64_t stack = units3xBM[5].u64[0], b4 = BIT_SET_2X & (~stack);
		//extuas.GetSort(b4);
		extuas.GetSortB(tuasb12.tua, tuasb12.nua, b4);
		//extuas.Dump();
		zh2b[0].InitBf(stack);

	}
	zh2b[0].Dob(22);
	if (zh2gxn.nua) {
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);	
			//cout << Char2Xout(w) << " to add b4 size " << cc << endl;
			if (tuasb12.nua < UA12SIZE)tuasb12.AddInit(w | (cc << 59), 0,0);
		}
	}
	if(op.ton)tuasb12.DumpShort("4 box ");
}

struct TTT1 {//tuab12;tua <=12 per size, then per floor
#define NB 3
	uint64_t tua0[NB*128], tua[NB * 128];
	uint32_t nua, tdigs0[NB * 128], tdigs[NB * 128],nbl;
	BF128 vs [9][NB] ;// vectors size 4/12
	struct VF {
		BF128 v2[36][NB], v3[84][NB], v4[126][NB], v5[126][NB], v6[84][NB];
	}vf;
	void Init(uint64_t* tu, uint32_t* td, uint32_t n) {
		//cout << "uas <=12 out of n="<<n << endl;
		memset(vs, 0, sizeof vs);
		nua = 0;
		for (uint32_t i = 0; i < n; i++) {
			register uint64_t U = tu[i],cc=U>>59;
			if (cc>=4 && cc <= 12) {
				tua0[nua] = U & BIT_SET_2X;
				tdigs0[nua] = td[i];
				uint32_t bloc = nua / 128, ir = nua - 128 * bloc;
				vs[cc-4][bloc].setBit(ir);
				nua++;
				if (nua >= NB * 128) break;// overflow control
			}
		}
		//cout << "uas <=12 got nua=" << nua << endl;
		// copy per size 
		nbl = (nua + 63) / 64;
		uint32_t x, n2 = 0;
		for (int i = 0; i < 9; i++) {// size 4/12
			uint64_t* v = vs[i][0].bf.u64;
			for (uint32_t j = 0; j < nbl; j++) {// uas for this size
				uint64_t V = v[j];
				while (bitscanforward64(x, V)) {
					V ^= (uint64_t)1 << x;
					register uint32_t ii = x + 64 * j;
					tua[n2] = tua0[ii]  ;
					tdigs[n2++] = tdigs0[ii];
				}
			}
		}
		//DumpUas();
		// build vectors per floor
		memset(&vf, 0, sizeof vf);
		for (int i = 0; i < 36; i++) {// ______________2 digits
			register uint32_t fl = floors_2d[i];
			BF128* pv = vf.v2[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc =j/ 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
		}
		for (int i = 0; i < 84; i++) {// ______________3 digits
			register uint32_t fl = floors_3d[i],fln=~fl;
			BF128* pv = vf.v3[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc = j / 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
			for (int i = 0; i < 36; i++) {// subset___2 digits
				uint32_t fl2 = floors_2d[i];
				if (fl2 & fln) continue;
				BF128* pv2 = vf.v2[i];
				for (uint32_t k = 0; k < nbl; k++) pv[k] |= pv2[k];
			}
		}
		for (int i = 0; i < 126; i++) {// ______________4 digits
			register uint32_t fl = floors_4d[i], fln = ~fl;
			BF128* pv = vf.v4[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc = j / 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
			for (int i = 0; i < 84; i++) {// subset___3 digits
				uint32_t fl2 = floors_3d[i];
				if (fl2 & fln) continue;
				BF128* pv2 = vf.v3[i];
				for (uint32_t k = 0; k < nbl; k++) pv[k] |= pv2[k];
			}
		}
		for (int i = 0; i < 126; i++) {// ______________5 digits
			register uint32_t fl = 0777^floors_4d[i], fln = ~fl;
			BF128* pv = vf.v5[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc = j / 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
			for (int i = 0; i < 126; i++) {// subset___4 digits
				uint32_t fl2 = floors_4d[i];
				if (fl2 & fln) continue;
				BF128* pv2 = vf.v4[i];
				for (uint32_t k = 0; k < nbl; k++) pv[k] |= pv2[k];
			}
		}
		for (int i = 0; i < 84; i++) {// ______________6 digits
			register uint32_t fl = 0777 ^ floors_3d[i], fln = ~fl;
			BF128* pv = vf.v6[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc = j / 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
			for (int i = 0; i < 126; i++) {// subset___5 digits
				uint32_t fl2 = 0777^floors_4d[i];
				if (fl2 & fln) continue;
				BF128* pv2 = vf.v5[i];
				for (uint32_t k = 0; k < nbl; k++) pv[k] |= pv2[k];
			}
		}
		//DumpV();
	}
	int Load(uint64_t* td, int iv, int ix) {
		BF128* v;
		switch (iv) {
		case 2: v = vf.v2[ix]; break;
		case 3: v = vf.v3[ix]; break;
		case 4: v = vf.v4[ix]; break;
		case 5: v = vf.v5[ix]; break;
		case 6: v = vf.v6[ix]; break;
		default: return 0;
		}
		int n = 0,x;
		uint64_t* pv = v->bf.u64;
		for (uint32_t i = 0; i < nbl; i++) {
			uint64_t V = pv[i];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				register uint32_t ii = x + 64 * i;
				td[n++] = tua[ii];
			}
		}
		return n;
	}
	void DumpUas() {
		cout << "uas <=12 sorted" << endl;
		if (nua < 256) for (uint32_t i = 0; i < nua; i++) {
			cout << Char2Xout(tua[i]) << "|";
			cout  << Char9out(tdigs[i]) << " i="<<i << endl;
		}
	}
	void DumpV() {
		cout << " dump v nbl=" << nbl << endl;
		cout << "v2" << endl;
		for (int i = 0; i < 36; i++) {
			uint64_t *pv=vf.v2[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k])<< "|";
			cout << Char9out(floors_2d[i]) <<" i="<<i << endl;

		}
		cout << "v3" << endl;
		for (int i = 0; i < 84; i++) {
			uint64_t* pv = vf.v3[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k])<< "|";
			cout << Char9out(floors_3d[i]) << " i=" << i << endl;

		}
		cout << "v4" << endl;
		for (int i = 0; i < 126; i++) {
			uint64_t* pv = vf.v4[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k]) << "|";
			cout << Char9out(floors_4d[i]) << " i=" << i << endl;

		}
		cout << "v5" << endl;
		for (int i = 0; i < 126; i++) {
			uint64_t* pv = vf.v5[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k]) << "|";
			cout << Char9out(0777^floors_4d[i]) << " i=" << i << endl;

		}
		cout << "v6" << endl;
		for (int i = 0; i < 84; i++) {
			uint64_t* pv = vf.v6[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k]) << "|";
			cout << Char9out(0777^floors_3d[i]) << " i=" << i << endl;

		}

	}
	//36;84;126;126 for 2;3;4;5; 

}ttt1;


void G17B::Guas2Collect() {
	ttt1.Init(tuasb12.tua, tuasb12.tdigs, tuasb12.nua);
	guah.Init();

	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {// initial gua2s 4 cells
		SGUA2 w = tsgua2[i81];
		uint64_t bf = zh2b_g.fd_sols[0][w.dig1].bf.u64;
		bf |= zh2b_g.fd_sols[0][w.dig2].bf.u64;
		uint64_t bf2 = units3xBM[9 + w.col1].u64[0];
		bf2 |= units3xBM[9 + w.col2].u64[0];
		bf &= bf2;
		int xcell1, xcell2;
		bitscanforward64(xcell1, bf);
		bitscanreverse64(xcell2, bf);
		int cell1=From_128_To_81[xcell1], cell2= From_128_To_81[xcell2];
		int r1 = C_row[cell1], r2 = C_row[cell2];
		if (r1 == r2) {
			//cout << Char54out(bf) << " gua2_4 i_81=" << i81 << endl;
			guah.Add2(bf, i81);
		}
	}

	Guas2CollectG2();
	Guas2CollectG3();

}
void G17B::Guas2CollectG2() {
	//___ find guas2 2 digits
	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {
		if (guah.IsUa4(i81)) continue;// nothing to do
		SGUA2 w = tsgua2[i81];
		GUAH::GUA& gt = guah.tg2[i81];
		int i = 0;
		for (i; i < 36; i++)if (floors_2d[i] == w.digs) break;
		extuas.nt2 = ttt1.Load(extuas.t2, 2, i);
		//if (extuas.nt2 && HasUaHitFalse(extuas.t2,0, extuas.nt2, i81))
		//	continue;

		int ir = zh2_2[0].GoZ2G2(w.digs, w.col1, w.dig1, w.col2, w.dig2);
		if (ir < 0)continue; // locked
		if (ir == 1) {//	solved)
			uint64_t U = zh2gxn.tua[0], cc = _popcnt64(U);
			guah.Add2(U, i81);
		}
		uint64_t F = zh2gxn.unsolved_field;
		zh2_2[0].DoZ2Go();
		if (zh2gxn.nua)
			for (uint32_t i = 0; i < zh2gxn.nua; i++) {
				uint64_t U = zh2gxn.tua[i];
				guah.Add2(U, i81);
			}
	}
	guah.SortClean();
	//___ find guas2 3 digits

	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {
		if (guah.IsUa4(i81)) continue;// nothing to do
		SGUA2 w = tsgua2[i81];
		GUAH::GUA& gt = guah.tg2[i81];
		for (int i = 0; i < 84; i++) {// find UAs 3 digits
			int fl = floors_3d[i];
			if (!((fl & w.digs) == w.digs)) continue;
			uint64_t isfl = 0;// cells ok of the floors
			for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
				if (fl & bit)	isfl |= zh2gxn.fsol[i2];
			extuas.nt2=gt.Load(extuas.t2, isfl);
			zh2gxn.nkguas = extuas.nt2;
			uint32_t istart = extuas.nt2;
			extuas.nt2 += ttt1.Load(&extuas.t2[istart], 3, i);
			// check no ua hitting cells "false"
			//if (extuas.nt2>istart && HasUaHitFalse(extuas.t2, istart,extuas.nt2, i81))
			//	continue;
			int ir = zh2_3[0].GoZ3G2(fl, w.col1, w.dig1, w.col2, w.dig2),ir2;
			if (ir < 0)continue; // locked
			if (ir == 1) {//	solved)
				if (zh2gxn.nua) {
					uint64_t U = zh2gxn.tua[0];
					guah.Add2(U, i81);
				}
				continue;
			}
			ir2 = zh2_3[0].DoZ3Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					guah.Add2(U, i81);
				}
		}
	}// end 3 digits
	//___ find guas2 4 digits
	guah.SortClean();

	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {
		if (guah.IsUa4(i81)) continue;// nothing to do
		SGUA2 w = tsgua2[i81];
		GUAH::GUA& gt = guah.tg2[i81];
		for (int i = 0; i < 126; i++) { 
			int fl = floors_4d[i];
			if (!((fl & w.digs) == w.digs)) continue;
			uint64_t isfl = 0;// cells ok of the floors
			for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
				if (fl & bit)	isfl |= zh2gxn.fsol[i2];
			extuas.nt2 = gt.Load(extuas.t2, isfl);
			zh2gxn.nkguas = extuas.nt2;
			uint32_t istart = extuas.nt2;
			extuas.nt2 += ttt1.Load(&extuas.t2[istart], 3, i);
			// check no ua hitting cells "false"
			//if (extuas.nt2 > istart && HasUaHitFalse(extuas.t2, istart, extuas.nt2, i81))
			//	continue;
			int ir = zh2_4[0].GoZ4G2(fl, w.col1, w.dig1, w.col2, w.dig2), ir2;
			if (ir < 0)continue; // locked
			if (ir == 1) {//	solved)
				if (zh2gxn.nua) {
					uint64_t U = zh2gxn.tua[0];
					guah.Add2(U, i81);
				}
				continue;
			}
			ir2 = zh2_4[0].DoZ4Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					guah.Add2(U, i81);
				}
		}
	}// end 4 digits
	//guah.DumpOne2(41);

	guah.SortClean();

	//___ find guas2 5 digits

	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {
		if (guah.IsUa4(i81)) continue;// nothing to do
		SGUA2 w = tsgua2[i81];
		GUAH::GUA& gt = guah.tg2[i81];
		for (int i = 0; i < 126; i++) { 
			int fl =0777^ floors_4d[i];
			if (!((fl & w.digs) == w.digs)) continue;
			uint64_t isfl = 0;// cells ok of the floors
			for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
				if (fl & bit)	isfl |= zh2gxn.fsol[i2];
			extuas.nt2 = gt.Load(extuas.t2, isfl);
			//if (diag2 == 2)extuas.Dump();
			zh2gxn.nkguas = extuas.nt2;
			uint32_t istart = extuas.nt2;
			extuas.nt2 += ttt1.Load(&extuas.t2[istart], 3, i);
			//if (diag2 == 2)extuas.Dump();
			// check no ua hitting cells "false"
			//if (extuas.nt2 > istart && HasUaHitFalse(extuas.t2, istart, extuas.nt2, i81))
			//	continue;
			int ir = zh2_5[0].GoZ5G2(fl, w.col1, w.dig1, w.col2, w.dig2), ir2;
			if (ir < 0)continue; // locked
			//if (diag2 == 2)cout << "back goz5g2 ir=" << ir << endl;
			if (ir == 1) {//	solved)
				if (zh2gxn.nua) {
					uint64_t U = zh2gxn.tua[0];
					guah.Add2(U, i81);
				}
				continue;
			}
			ir2 = zh2_5[0].DoZ5Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					guah.Add2(U, i81);
				}
		}
	}// end 5 digits
	guah.SortClean();

	ng2=guah.CutG2(30);
	if (op.ton)cout << "ng2 after cut30=" << ng2 << endl;
	//guah.DumpOne2(41);

	//guah.Dump2all2();

}
void G17B::Guas2CollectG3() {
	for (int i81 = 0; i81 < 81; i81++) if (gsock3.On(i81)) {

		//if (i81)continue;
		SGUA3 w = tsgua3[i81];
		GUAH::GUA& gt = guah.tg3[i81];

		// Setup the perms for gangsters in minirow
		int bita = 1 << w.dig1, bitb = 1 << w.dig2, bitc = 1 << w.dig3,
			digs = w.digs,triplet_perms[2][3];

		int* p = triplet_perms[0];// for per abc -> bca
		p[0] = bita | bitb; p[1] = bitb | bitc; p[2] = bitc | bita;

		p = triplet_perms[1];// for per abc -> cab
		p[0] = bita | bitc; p[1] = bitb | bita; p[2] = bitc | bitb;
		int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
		for (int ip = 0; ip < 2; ip++) {
			// build revised gangster
			int rgangcols[9];// revised gangster
			memcpy(rgangcols, zh2gxn.gangsters, sizeof rgangcols);
			p = triplet_perms[ip];
			int c1 = w.col1, c2 = c1 + 1, c3 = c1 + 2;
			rgangcols[c1] ^= p[0];
			rgangcols[c2] ^= p[1];
			rgangcols[c3] ^= p[2];
			int i = 0;
			for (i; i < 84; i++)if (floors_3d[i] == w.digs) break;
			extuas.nt2 = ttt1.Load(extuas.t2, 3, i);
			//extuas.Dump();
			// find UAs 3 digits one 3 digits 
			int ir= zh2_3[0].GoZ3G3(w.digs, rgangcols ),ir2;
			if (ir < 0)continue; // locked
			//if (extuas.nt2 && HasUaHitFalse3(extuas.t2, 0, extuas.nt2, i81))
			//	continue;
			//cout << "i81=" << i81 << " digs=" << w.dig1 + 1 << w.dig2 + 1 << w.dig3 + 1
			//	<< " col1=" << c1 + 1<<" ir ="<<ir << endl;
			//zh2_3[0].ImageCandidats();
			if (ir == 1) {//	solved)
				if (zh2gxn.nua) {
					uint64_t U = zh2gxn.tua[0];
					guah.Add3(U, i81);
				}
				continue;
			}
			ir2 = zh2_3[0].DoZ3Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					//cout << Char2Xout(U) << "to add " << _popcnt64(U) << endl;
					guah.Add3(U, i81);
				}
		}
	}
	//guah.DumpOne3(46);
	//guah.Dump2all3();
	Guas2CollectG3_4d();
}
void G17B::Guas2CollectG3_4d() {
	for (int i81 = 0; i81 < 81; i81++) if (gsock3.On(i81)) {
		SGUA3 w = tsgua3[i81];
		GUAH::GUA& gt = guah.tg3[i81];
		if (guah.IsUamin(i81)) continue;// nothing to do

		// Setup the perms for gangsters in minirow
		int bita = 1 << w.dig1, bitb = 1 << w.dig2, bitc = 1 << w.dig3,
			digs = w.digs, triplet_perms[2][3];

		int* p = triplet_perms[0];// for per abc -> bca
		p[0] = bita | bitb; p[1] = bitb | bitc; p[2] = bitc | bita;

		p = triplet_perms[1];// for per abc -> cab
		p[0] = bita | bitc; p[1] = bitb | bita; p[2] = bitc | bitb;
		int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
		for (int ip = 0; ip < 2; ip++) {
			//cout << i81 << " perm " << ip << endl;
			// build revised gangster
			int rgangcols[9];// revised gangster
			memcpy(rgangcols, zh2gxn.gangsters, sizeof rgangcols);
			p = triplet_perms[ip];
			int c1 = w.col1, c2 = c1 + 1, c3 = c1 + 2;
			rgangcols[c1] ^= p[0];
			rgangcols[c2] ^= p[1];
			rgangcols[c3] ^= p[2];
			for (int i = 0; i < 126; i++) {
				int fl = floors_4d[i];
				if (!((fl & w.digs) == w.digs)) continue;
				uint64_t isfl = 0;// cells ok of the floors
				for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
					if (fl & bit)	isfl |= zh2gxn.fsol[i2];
				extuas.nt2 = gt.Load(extuas.t2, isfl);
				zh2gxn.nkguas = extuas.nt2;
				uint32_t istart = extuas.nt2;
				extuas.nt2 += ttt1.Load(&extuas.t2[istart], 4, i);
				//cout << "nolds=" << istart << " ";
				int ir = zh2_4[0].GoZ4G3(fl, rgangcols), ir2;
				if (ir < 0)continue; // locked
				//extuas.Dump();
				//zh2_4[0].ImageCandidats();
				if (ir == 1) {//	solved)
					if (zh2gxn.nua) {
						uint64_t U = zh2gxn.tua[0];
						if (_popcnt64(U) < 18) {
							//cout << Char2Xout(U) << "to add " << _popcnt64(U) << endl;
							guah.Add3(U, i81);
						}
					}
					continue;
				}
				ir2 = zh2_4[0].DoZ4Go();
				if (ir2 < 0) continue;
				if (zh2gxn.nua)
					for (uint32_t i = 0; i < zh2gxn.nua; i++) {
						uint64_t U = zh2gxn.tua[i];
						if (_popcnt64(U) < 18) {
							//cout << Char2Xout(U) << "to add " << _popcnt64(U) << endl;
							guah.Add3(U, i81);
						}
					}
			}
		}
	}
	guah.SortClean3();
	//guah.DumpOne3(46);
	//guah.Dump2all3();
	Guas2CollectG3_5d();
}
void G17B::Guas2CollectG3_5d() {
	for (int i81 = 0; i81 < 81; i81++) if (gsock3.On(i81)) {
		SGUA3 w = tsgua3[i81];
		GUAH::GUA& gt = guah.tg3[i81];
		if (guah.IsUamin(i81)) continue;// nothing to do
		if (gt.nua >= 10) continue;

		// Setup the perms for gangsters in minirow
		int bita = 1 << w.dig1, bitb = 1 << w.dig2, bitc = 1 << w.dig3,
			digs = w.digs, triplet_perms[2][3];

		int* p = triplet_perms[0];// for per abc -> bca
		p[0] = bita | bitb; p[1] = bitb | bitc; p[2] = bitc | bita;

		p = triplet_perms[1];// for per abc -> cab
		p[0] = bita | bitc; p[1] = bitb | bita; p[2] = bitc | bitb;
		int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
		for (int ip = 0; ip < 2; ip++) {
			//if(i81==12) cout << i81 << " perm 5d" << ip << endl;
			// build revised gangster
			int rgangcols[9];// revised gangster
			memcpy(rgangcols, zh2gxn.gangsters, sizeof rgangcols);
			p = triplet_perms[ip];
			int c1 = w.col1, c2 = c1 + 1, c3 = c1 + 2;
			rgangcols[c1] ^= p[0];
			rgangcols[c2] ^= p[1];
			rgangcols[c3] ^= p[2];
			for (int i = 0; i < 126; i++) {
				int fl = 0777^floors_4d[i];
				if (!((fl & w.digs) == w.digs)) continue;
				uint64_t isfl = 0;// cells ok of the floors
				for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
					if (fl & bit)	isfl |= zh2gxn.fsol[i2];
				extuas.nt2 = gt.Load(extuas.t2, isfl);
				zh2gxn.nkguas = extuas.nt2;
				uint32_t istart = extuas.nt2;
				extuas.nt2 += ttt1.Load(&extuas.t2[istart], 5, i);
				int ir = zh2_5[0].GoZ5G3(fl, rgangcols), ir2;
	//if (i81 == 46)cout<<"i="<<i << " nolds=" << istart << " nt2=" << extuas.nt2 << "ir=" << ir << endl;

				if (ir < 0)continue; // locked
				//if (i81 == 12)extuas.Dump();
				//zh2_5[0].ImageCandidats();
				if (ir == 1) {//	solved)
					if (zh2gxn.nua) {
						uint64_t U = zh2gxn.tua[0];
						if (_popcnt64(U) < 18) {
	//if (i81 == 46)cout << Char2Xout(U) << "to add1 " << _popcnt64(U) << endl;
							guah.Add3(U, i81);
						}
					}
					continue;
				}

				ir2 = zh2_5[0].DoZ5Go();
				if (ir2 < 0) continue;
				if (zh2gxn.nua)
					for (uint32_t i = 0; i < zh2gxn.nua; i++) {
						uint64_t U = zh2gxn.tua[i];
						if (_popcnt64(U) < 18) {
		//if (i81 == 46)cout << Char2Xout(U) << "to add 2 " << _popcnt64(U) << endl;
							guah.Add3(U, i81);
						}
					}
			}
		}
	}
	guah.SortClean3();
	ng3 = guah.CutG3(20);
	if (op.ton)cout << "ng3 after cut20=" << ng3 << endl;
	//guah.DumpOne3(46);
	//guah.Dump2all3();
}


//__________ end of initial uas harvest switch to 54 mode

void G17B::StartAfterUasHarvest() {
	guah54.Build();
	t54b12.Build_ta128(tuasb12.tua, tuasb12.nua);
	if (op.ton) {
		if (sgo.bfx[1] & 8)t54b12.DebugA();
		//if (sgo.bfx[3] & 1) { guah54.Dumpall2(); guah54.Dumpall3(); }
	}
	//if (op.known)genb12.bands3[0].DumpGuam(1);
	Expand_03();
}

void T54B12::Build_ta128(uint64_t* t, uint32_t n) {

	// here insert all missing uas one band (1 or 2)
	for (int i = 0; i < 2; i++) {
		STD_B416& wb = (i) ? myband2 : myband1;
		BF64 wt; wt.bf.u64 = 0;
		uint32_t* tu = wb.tua, nu = wb.nua;
		for (uint32_t j = 0; j < nu; j++) {
			register uint32_t u = tu[j] & BIT_SET_27, cc = _popcnt32(u);
			wt.bf.u32[i] = u;
			if (cc > 12)if(n< UA12SIZE)t[n++] = wt.bf.u64;
		}
	}	
	InitA();
	uint32_t nbl64 = (n + 63) >> 6, x;
	memset(vsize, 0, sizeof vsize);
	for (uint32_t i = 0; i < n; i++) {
		int bloc = i >> 7, ir = i - (bloc << 7);
		uint64_t cc = t[i] >> 59;
		if (cc > 24) continue;//safety
		vsize[cc][bloc].setBit(ir);
	}
	uint32_t ntw = 0;
	for (int i1 = 4; i1 < 25; i1++) {
		uint64_t* tb64 = vsize[i1]->bf.u64;
		for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
			register uint64_t V = tb64[i2];
			while (bitscanforward64(x, V)) {// switch to 54 mode here
				V ^= (uint64_t)1 << x;
				register uint64_t R = t[x + (i2 << 6)];
				R = (R & BIT_SET_27) | ((R & BIT_SET_B2) >> 5);
				if (1) {// check redundancy and subsets
					for (uint32_t i = 0; i < ntw; i++) {
						if (!(R & (~tw[i]))) {
							//cout << Char54out(R) << " erased" << endl;;
							//cout << Char54out(tw[i]) <<" i="<<i << endl;
							R = 0; break;						}
					}
					if (R)tw[ntw++] = R;
					else continue;
				}
				AddA(R);
			}
		}
	}
	if (op.ton > 1)cout << "build ta128 ntw=" << ntw << " nold" << n << endl;
}
int T54B12::Build_tb128(SPB03A& s) {
	int tc[3], ntc = 0;
	{// // build table of clues 
		int cell;
		register uint64_t U = s.all_previous_cells;
		while (bitscanforward64(cell, U)) {
			U ^= (uint64_t)1 << cell;
			tc[ntc++] = cell;
		}
	}

	uint32_t lastbloc = t54b12.nablocs;
	tvw[0] = s.v;
	for (uint32_t i = 1; i <= lastbloc; i++) {
		TUVECT& vv = ta128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 0; ic < 3; ic++)
			v &= vc[tc[ic]];
		tvw[i] = v;
	}

	// apply active on still valid uas and flag by size
	memset(vsize, 0, sizeof vsize);
	uint32_t ntw = 0;
	{
		register uint64_t Ac = s.active_cells;
		for (uint32_t i = 0; i <= lastbloc; i++) {
			register uint64_t* t = ta128[i].t;
			BF128 V = tvw[i];
			while (1) {
				register int ir = V.getFirst128();
				if (ir >= 0) {
					V.clearBit(ir);
					register uint64_t R = t[ir] & Ac;
					if (!R)return 1; //dead
					register uint64_t cc = _popcnt64(R);
					if (cc > 18)cc = 18;
					int bloc = ntw >> 7, ir = ntw - (bloc << 7);
					vsize[cc][bloc].setBit(ir);
					tw[ntw++] = R;
					if (ntw >= 128 * UABNBLOCS) break;
				}
				else break;
			}
			if (ntw >= 128 * UABNBLOCS) break;
		}
	}
	//__Build the reduced set of UAs vectors
	InitB();
	uint32_t nbl64 = (ntw + 63) >> 6, x;
	for (int i1 = 1; i1 < 19; i1++) {
		uint64_t* tb64 = vsize[i1]->bf.u64;
		for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
			register uint64_t V = tb64[i2];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				AddB(tw[x + (i2 << 6)]);
			}
		}
	}
	if (nb128 < 128) nb128 = 128; //forcing add outside first bloc
	if (op.ton > 1)cout << "ua3c ntw=" << ntw << " nb128=" << nb128 << endl;

	return 0;
}
int T54B12::Build_tc128(SPB03A& s3, SPB03A& s6) {
	int tc[3], ntc = 0;
	{// // build table of clues 
		int cell;
		register uint64_t U = s6.all_previous_cells
			& (~s3.all_previous_cells);// fresh clues 
		while (bitscanforward64(cell, U)) {
			U ^= (uint64_t)1 << cell;
			tc[ntc++] = cell;
		}
	}

	uint32_t lastbloc = t54b12.nbblocs;
	tvw[0] = s6.v;
	for (uint32_t i = 1; i <= lastbloc; i++) {
		TUVECT& vv = tb128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 0; ic < 3; ic++)
			v &= vc[tc[ic]];
		tvw[i] = v;
	}

	// apply active on still valid uas and flag by size
	memset(vsize, 0, sizeof vsize);
	uint32_t ntw = 0;
	{
		register uint64_t Ac = s6.active_cells;
		for (uint32_t i = 0; i <= lastbloc; i++) {
			register uint64_t* t = tb128[i].t;
			BF128 V = tvw[i];
			while (1) {
				register int ir = V.getFirst128();
				if (ir >= 0) {
					V.clearBit(ir);
					register uint64_t R = t[ir] & Ac;
					if (!R)return 1; //dead
					register uint64_t cc = _popcnt64(R);
					if (cc > 18)cc = 18;
					int bloc = ntw >> 7, ir = ntw - (bloc << 7);
					vsize[cc][bloc].setBit(ir);
					tw[ntw++] = R;
					if (ntw >= 128 * UACNBLOCS) break;
				}
				else break;
			}
			if (ntw >= 128 * UACNBLOCS) break;
		}
	}
	//__Build the reduced set of UAs vectors clean redundancy
	InitC(); InitD();// be sure to have D ready for fresh uas
	uint32_t nbl64 = (ntw + 63) >> 6, x;
	for (int i1 = 1; i1 < 19; i1++) {
		uint64_t* tb64 = vsize[i1]->bf.u64;
		for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
			register uint64_t V = tb64[i2];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				register uint64_t U = tw[x + (i2 << 6)];
				// check redundancy in tc128[0]
				if (t54b12.IsNotRedundant(U))
					AddC(U);
			}
		}
	}
	if (op.ton > 2)cout << "ua6c ntw=" << ntw << " nc128=" << nc128 << endl;
	return 0;
}
int T54B12::Build_td128(SPB03B &s9) {
	uint32_t lastbloc = t54b12.ncblocs;
	// apply active on still valid uas and flag by size
	memset(vsize, 0, sizeof vsize);
	uint32_t ntw = 0;
	{
		register uint64_t Ac = s9.active_cells;
		for (uint32_t i = 0; i <= lastbloc; i++) {
			register uint64_t* t = tc128[i].t;
			BF128 V = s9.v[i];
			while (1) {
				register int ir = V.getFirst128();
				if (ir >= 0) {
					V.clearBit(ir);
					register uint64_t R = t[ir] & Ac;
					if (!R)return 1; //dead
					register uint64_t cc = _popcnt64(R);
					if (cc > 18)cc = 18;
					int bloc = ntw >> 7, ir = ntw - (bloc << 7);
					vsize[cc][bloc].setBit(ir);
					tw[ntw++] = R;
					if (ntw >= 128 * UADNBLOCS) break;
				}
				else break;
			}
			if (ntw >= 128 * UADNBLOCS) break;
		}
	}
	//__Build the reduced set of UAs vectors clean redundancy
	InitD();
	uint32_t nbl64 = (ntw + 63) >> 6, x;
	for (int i1 = 1; i1 < 19; i1++) {
		uint64_t* tb64 = vsize[i1]->bf.u64;
		for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
			register uint64_t V = tb64[i2];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				register uint64_t U = tw[x + (i2 << 6)];
				// check redundancy in tc128[0]
				if (t54b12.IsNotRedundantD(U))
					AddD(U);
			}
		}
	}
	return 0;
}

void GUAH54::Build() {// cut to 30 switch to 54 killer

	uint64_t* pbuf = gbuf;
	for (int i81 = 0; i81 < 81; i81++) {
		tg2[i81].Init(pbuf, 0, i81);
		if (g17b.gsock2.On(i81)) {
			GUAH::GUA& g0 = guah.tg2[i81];
			GUA54& gd = tg2[i81];
			uint32_t n = g0.nua;
			if (n > 30) n = 30;
			gd.nua = n;
			gd.nuamax = gd.nua + 50;
			pbuf += gd.nuamax;// lock the storing place
			register uint64_t K = ~0;
			for (uint32_t i = 0; i < n; i++) {
				register uint64_t U = g0.tua[i];
				U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
				K &= U;
				gd.tua[i] = U;
			}
			gd.killer = K;
		}
	}
	for (int i81 = 0; i81 < 81; i81++) {
		tg3[i81].Init(pbuf, 1, i81);
		if (g17b.gsock3.On(i81)) {
			GUAH::GUA& g0 = guah.tg3[i81];
			GUA54& gd = tg3[i81];
			uint32_t n = g0.nua;
			if (n > 30) n = 30;
			gd.nua = n;
			gd.nuamax = gd.nua + 10;
			pbuf += gd.nuamax;// lock the storing place
			register uint64_t K = ~0;
			for (uint32_t i = 0; i < n; i++) {
				register uint64_t U = g0.tua[i];
				U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
				K &= U;
				gd.tua[i] = U;
			}
			gd.killer = K;
		}
	}
	//DumpOne2(41);
	//Dumpall2();
}
void GUAH54::Build2(uint64_t filter, uint64_t active) {
	register uint64_t F = filter,
		A = active;

	uint64_t* pbuf = gbuf;
	uint64_t tw[60];// valid
	for (int i81 = 0; i81 < 81; i81++) {
		tg2[i81].Init(pbuf, 0, i81);
		if (g17b.gsock2.On(i81)) {
			GUA54& g0 = guah54.tg2[i81];
			GUA54& gd = tg2[i81];
			if (g0.killer & F) {
				pbuf += gd.nuamax;// lock the storing place
				continue; // killed or not valid
			}
			uint32_t n1 = g0.nua;
			if (n1 == 1) {// killer is enough
				register uint64_t U = g0.tua[0] & A;
				gd.Add(U);
				pbuf += gd.nuamax;// lock the storing place
				//pbuf++;
				continue;
			}
			uint64_t vsize[25];// sorting tw2 by size
			memset(vsize, 0, sizeof vsize);
			uint32_t ntw = 0, iw;
			// store still valid for active cells
			for (uint32_t i = 0; i < n1; i++) {
				register uint64_t U = g0.tua[i];
				if (U & F)continue; // hit
				U &= A;// keep active cells
				uint64_t cc = _popcnt64(U);
				if (cc > 20)cc = 20;
				vsize[cc] |= (uint64_t)1 << ntw;
				tw[ntw++] = U;
				if (ntw >= 40) break;
			}
			if (vsize[0]) {// force true (one not hit)
				//cout << "0 g2 found i81=" << i81 << endl;
				gd.Add(0);
				pbuf++;
				continue;
			}			// take back per size and send to main table
			uint32_t istart = ntw; // start for redundancy
			for (int i = 0; i <= 20; i++) {
				register uint64_t V = vsize[i];
				while (bitscanforward64(iw, V)) {
					V ^= (uint64_t)1 << iw;
					gd.AddCheck(tw[iw]);
				}
				if (gd.nua >= 30) break;
			}
			gd.nuamax = gd.nua + 20;
			pbuf += gd.nuamax;// lock the storing place
		}
	}
	for (int i81 = 0; i81 < 81; i81++) {
		tg3[i81].Init(pbuf, 1, i81);
		if (g17b.gsock3.On(i81)) {
			GUA54& g0 = guah54.tg3[i81];
			GUA54& gd = tg3[i81];
			if (g0.killer & F) {
				pbuf += gd.nuamax;// lock the storing place
				continue; // killed or not valid
			}
			uint32_t n1 = g0.nua;
			if (n1 == 1) {// killer is enough
				register uint64_t U = g0.tua[0] & A;
				gd.Add(U);
				pbuf += gd.nuamax;// lock the storing place
				//pbuf++;
				continue;
			}
			uint64_t vsize[25];// sorting tw2 by size
			memset(vsize, 0, sizeof vsize);
			uint32_t ntw = 0, iw;
			// store still valid for active cells
			for (uint32_t i = 0; i < n1; i++) {
				register uint64_t U = g0.tua[i];
				if (U & F)continue; // hit
				U &= A;// keep active cells
				uint64_t cc = _popcnt64(U);
				if (cc > 20) cc = 20;
				vsize[cc] |= (uint64_t)1 << ntw;
				tw[ntw++] = U;
				if (ntw >= 40) break;
			}
			if (vsize[0]) {// force true (one not hit)
				//cout << "0 g3 found i81=" << i81 << endl;
				gd.Add(0);
				pbuf++;
				continue;
			}			// take back per size and send to main table
			uint32_t istart = ntw; // start for redundancy
			for (int i = 0; i <= 20; i++) {
				register uint64_t V = vsize[i];
				while (bitscanforward64(iw, V)) {
					V ^= (uint64_t)1 << iw;
					gd.AddCheck(tw[iw]);
				}
				if (gd.nua >= 30) break;
			}
			gd.nuamax = gd.nua + 10;
			pbuf += gd.nuamax;// lock the storing place
		}
	}
	//DumpOne2(41);

	//Dumpall2();
	//Dumpall3();
}
void GUAH54::Build9(uint64_t filter, uint64_t active) {
	register uint64_t F = filter,		A = active;
	uint64_t* pbuf = gbuf;
	uint64_t tw[60];// valid
	for (int i81 = 0; i81 < 81; i81++) {
		tg2[i81].Init(pbuf, 0, i81);
		if (g17b.gsock2.On(i81)) {
			GUA54& g0 = guah54_2.tg2[i81];
			GUA54& gd = tg2[i81];
			if (g0.killer & F) {
				pbuf += gd.nuamax;// lock the storing place
				continue; // killed or not valid
			}
			uint32_t n1 = g0.nua;
			if (n1 == 1) {// killer is enough
				register uint64_t U = g0.tua[0] & A;
				gd.Add(U);
				pbuf += gd.nuamax;// lock the storing place
				//pbuf++;
				continue;
			}
			uint64_t vsize[25];// sorting tw2 by size
			memset(vsize, 0, sizeof vsize);
			uint32_t ntw = 0, iw;
			// store still valid for active cells
			for (uint32_t i = 0; i < n1; i++) {
				register uint64_t U = g0.tua[i];
				if (U & F)continue; // hit
				U &= A;// keep active cells
				uint64_t cc = _popcnt64(U);
				if (cc > 20)cc = 20;
				vsize[cc] |= (uint64_t)1 << ntw;
				tw[ntw++] = U;
				if (ntw >= 40) break;
			}
			if (vsize[0]) {// force true (one not hit)
				//cout << "0 g2 found i81=" << i81 << endl;
				gd.Add(0);
				pbuf++;
				continue;
			}			// take back per size and send to main table
			uint32_t istart = ntw; // start for redundancy
			for (int i = 0; i <= 20; i++) {
				register uint64_t V = vsize[i];
				while (bitscanforward64(iw, V)) {
					V ^= (uint64_t)1 << iw;
					gd.AddCheck(tw[iw]);
				}
				if (gd.nua >= 30) break;
			}
			gd.nuamax = gd.nua + 20;
			pbuf += gd.nuamax;// lock the storing place
		}
	}
	for (int i81 = 0; i81 < 81; i81++) {
		tg3[i81].Init(pbuf, 1, i81);
		if (g17b.gsock3.On(i81)) {
			GUA54& g0 = guah54_2.tg3[i81];
			GUA54& gd = tg3[i81];
			if (g0.killer & F) {
				pbuf += gd.nuamax;// lock the storing place
				continue; // killed or not valid
			}
			uint32_t n1 = g0.nua;
			if (n1 == 1) {// killer is enough
				register uint64_t U = g0.tua[0] & A;
				gd.Add(U);
				pbuf += gd.nuamax;// lock the storing place
				//pbuf++;
				continue;
			}
			uint64_t vsize[25];// sorting tw2 by size
			memset(vsize, 0, sizeof vsize);
			uint32_t ntw = 0, iw;
			// store still valid for active cells
			for (uint32_t i = 0; i < n1; i++) {
				register uint64_t U = g0.tua[i];
				if (U & F)continue; // hit
				U &= A;// keep active cells
				uint64_t cc = _popcnt64(U);
				if (cc > 20) cc = 20;
				vsize[cc] |= (uint64_t)1 << ntw;
				tw[ntw++] = U;
				if (ntw >= 40) break;
			}
			if (vsize[0]) {// force true (one not hit)
				//cout << "0 g3 found i81=" << i81 << endl;
				gd.Add(0);
				pbuf++;
				continue;
			}			// take back per size and send to main table
			uint32_t istart = ntw; // start for redundancy
			for (int i = 0; i <= 20; i++) {
				register uint64_t V = vsize[i];
				while (bitscanforward64(iw, V)) {
					V ^= (uint64_t)1 << iw;
					gd.AddCheck(tw[iw]);
				}
				if (gd.nua >= 30) break;
			}
			gd.nuamax = gd.nua + 10;
			pbuf += gd.nuamax;// lock the storing place
		}
	}
	//DumpOne2(41);

	//Dumpall2();
	//Dumpall3();
}


BF128 GUAH54::GetG2(uint64_t bf) {
	register uint64_t F = bf;
	BF128 g2;
	g2.SetAll_0();
	for (int i81 = 0; i81 < 81; i81++) {
		GUA54& g0 = tg2[i81];
		if (g0.killer & F)continue; // killed or not valid
		uint32_t n1 = g0.nua;
		if (n1 == 1) {	g2.setBit(i81);	continue;}
		for (uint32_t i = 0; i < n1; i++) {
			register uint64_t U = g0.tua[i];
			if (U & F)continue; // hit
			g2.setBit(i81);	// first not hit is ok
			break;
		}
	}
	return g2;
}
BF128 GUAH54::GetG3(uint64_t bf) {
	register uint64_t F = bf;
	BF128 g3;
	g3.SetAll_0();
	for (int i81 = 0; i81 < 81; i81++) {
		GUA54& g0 = tg3[i81];
		if (g0.killer & F)continue; // killed or not valid
		uint32_t n1 = g0.nua;
		if (n1 == 1) { g3.setBit(i81);	continue; }
		for (uint32_t i = 0; i < n1; i++) {
			register uint64_t U = g0.tua[i];
			if (U & F)continue; // hit
			g3.setBit(i81);	// first not hit is ok
			break;
		}
	}
	return g3;
}

//_______________________ expand common step 1_9 clues

///________ start expand uas bands 12
/// target 10 clues 3+3  +4
/// target 11 clues 3+3  +5
/// target 12 clues 3+3  +3+3
/// end of early steps skrinking and restructuring uas bands 1+2
/// and building a reduced GUAs table 

void G17B::Expand_03() {
	if (aigstop) return;
	SPB03A sp0, sp1,sp2, sp3;
	T54B12::TUVECT& tuv128 = t54b12.ta128[0];
	//tuv128.Dump(30);
	uint64_t* twu = tuv128.t;
	memset(&sp0, 0, sizeof sp0);
	sp0.active_cells = maskLSB[54].u64[0];
	sp0.possible_cells = twu[0];
	sp0.v = tuv128.v0;// initial nothing done
next1:
	{
		register int cell;
		uint64_t p = sp0.possible_cells;
		if (!p) return;
		bitscanforward64(cell, p);
		register uint64_t bit = (uint64_t)1 << cell;
		sp0.possible_cells ^= bit;
		sp0.active_cells ^= bit;
		sp1 = sp0; 
		sp1.all_previous_cells |= bit;
		sp1.v &= tuv128.vc[cell];
		if (!(sp1.possible_cells = twu[sp1.v.getFirst128()] & sp0.active_cells))goto next1;
	}
next2:
	{
		register int cell2;
		uint64_t p = sp1.possible_cells;
		if (!p) goto next1;
		bitscanforward64(cell2, p);
		register uint64_t bit = (uint64_t)1 << cell2;
		sp1.possible_cells ^= bit;
		sp1.active_cells ^= bit;
		sp2 = sp1;
		sp2.all_previous_cells |= bit;
		sp2.v &= tuv128.vc[cell2];
		if (!(sp2.possible_cells = twu[sp2.v.getFirst128()] & sp1.active_cells))goto next2;
	}
next3:
	{
		register int cell3;
		uint64_t p = sp2.possible_cells;
		if (!p) goto next2;
		bitscanforward64(cell3, p);
		register uint64_t bit = (uint64_t)1 << cell3;
		sp2.possible_cells ^= bit;
		sp2.active_cells ^= bit;
		sp3 = sp2;
		sp3.all_previous_cells |= bit;
		sp3.v &= tuv128.vc[cell3];
		//if (!(sp3.possible_cells = twu[sp3.v.getFirst128()] & sp0.active_cells))goto next3;
		p_cpt2g[3]++;
		if (t54b12.Build_tb128(sp3)) goto next3;
		if (op.known > 1) {
			//if (op.known)genb12.bands3[0].DumpGuam(1);
			if (!((~pk54) & sp3.all_previous_cells)) {
				cout << Char54out(sp3.all_previous_cells) << " expected 3" << endl;
				Expand_46(sp3);
				aigstop = 1;
				goto next3;
			}
		}
		Expand_46(sp3);
		if (knownt >= 9)return;
		if (aigstop) return;
		goto next3;

	}
}
void G17B::Expand_46(SPB03A& s3) {
	if (aigstop) return;
	//if (op.known > 1) {
		//cout << Char54out(s3.all_previous_cells) << " entry 3 [3] "
			//<< p_cpt2g[3] << endl;
	//}
	T54B12::TUVECT& tuv128 = t54b12.tb128[0];
	uint64_t* twu = tuv128.t;
	SPB03A sp3, sp4, sp5, sp6;
	sp3 = s3;
	sp3.possible_cells = twu[0];
	sp3.v = tuv128.v0;// initial nothing done
	int locdiag = 0;
	if (op.ton) {
		if (op.ton>1)cout << Char54out(s3.all_previous_cells) << " 3clues [3]" << p_cpt2g[3]  << endl;
		if (op.f3) {
			if (p_cpt2g[3] == op.f3) {
				cout << "call 4_6 good path" << endl;
				if (op.ton > 2) 	tuv128.Dump(30);
				locdiag = 1;
			}
			else {
				if (p_cpt2g[3] > op.f3) { cout << "stop" << endl;	aigstop = 1; return; }
				if (!(op.upto3)) return;
			}
		}
	}
next4:
	{
		register int cell;
		uint64_t p = sp3.possible_cells;
		if (!p) return;
		bitscanforward64(cell, p);
		register uint64_t bit = (uint64_t)1 << cell;
		sp3.possible_cells ^= bit;
		sp3.active_cells ^= bit;
		sp4 = sp3;
		sp4.all_previous_cells |= bit;
		sp4.v &= tuv128.vc[cell];
		if (!(sp4.possible_cells = twu[sp4.v.getFirst128()] & sp3.active_cells))goto next4;
	};
next5:
	{
		register int cell2;
		uint64_t p = sp4.possible_cells;
		if (!p) goto next4;
		bitscanforward64(cell2, p);
		register uint64_t bit = (uint64_t)1 << cell2;
		sp4.possible_cells ^= bit;
		sp4.active_cells ^= bit;
		sp5 = sp4;
		sp5.all_previous_cells |= bit;
		sp5.v &= tuv128.vc[cell2];
		if (!(sp5.possible_cells = twu[sp5.v.getFirst128()] & sp4.active_cells))goto next4;
	}
next6:
	{
		register int cell3;
		uint64_t p = sp5.possible_cells;
		if (!p) goto next5;
		bitscanforward64(cell3, p);
		register uint64_t bit = (uint64_t)1 << cell3;
		sp5.possible_cells ^= bit;
		sp5.active_cells ^= bit;
		sp6 = sp5;
		sp6.all_previous_cells |= bit;
		sp6.v &= tuv128.vc[cell3];
		if (!(sp6.possible_cells = twu[sp6.v.getFirst128()] & sp5.active_cells)) {
			// no possible valid 6 clues should not come
			if (zh2b[1].IsValid(sp6.all_previous_cells)) {
				uint32_t i = zh2gxn.nua - 1;
				register uint64_t ua = zh2gxn.tua[i],ua54 = 
					(ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
				t54b12.AddA(ua54);		t54b12.AddB(ua54);
			}
			else {
				cout << "bug exp 4-7 lower 7" << endl;aigstop = 1; return;			}
		}
		if (t54b12.Build_tc128(s3,sp6)) goto next6;
		p_cpt2g[4]++;
		if (op.known > 1) {
			if (knownt >= 9)return;
			if (!((~pk54) & sp6.all_previous_cells)) {
				cout << Char54out(sp6.all_previous_cells) << " expected 6 [4] " 
					<<p_cpt2g[4] << endl;
				Expand_7_9(sp6);
				knownt=6;
				aigstop = 1;
				goto next6;
			}			
		}
		if (knownt >= 6) { aigstop = 1; return; }
		Expand_7_9(sp6);
		goto next6;

	}

}
void G17B::Expand_7_9(SPB03A& s6) {
	if (aigstop) return;
	if (knownt >= 9)return;
	if (t54b12.nc128<128)t54b12.nc128=128;// be sure to have fresh uas outside
	//if (op.known > 1){
		//cout << Char54out(s6.all_previous_cells) << " entry 6 [4] "
		//	<< p_cpt2g[4] << endl;
		//cout << Char54out(s6.active_cells) << " active "<< endl;
	//}
	// ____ build a reduced table of uas/guas for band 3
	guah54_2.Build2(s6.all_previous_cells, s6.active_cells);
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
		genb12.bands3[ib3].BuildGuam2(s6.all_previous_cells);

	SPB03A sp6,sp7,sp8,sp9,sp10;
	T54B12::TUVECT& tuv128 = t54b12.tc128[0];
	uint64_t* twu = tuv128.t;
	uint32_t tcells[5];
	int locdiag = 0;
	if (op.ton) {
		if (op.f3) {
			if (p_cpt2g[3] == op.f3) {
				cout << Char54out(s6.all_previous_cells) << " 6clues [4]" << p_cpt2g[4]
					<< " nc128=" << t54b12.nc128;
				cout << Char54out(twu[0]) << " [7]"<< p_cpt2g[7] << endl;

			}
		}
		if (op.f4) {
			if (p_cpt2g[4] == op.f4) {
				cout << "call 7_9 good path" << endl;
				if (op.ton > 2) {tuv128.Dump(30); locdiag = 1;	}
			}
			else {
				if (p_cpt2g[4] > op.f4) { cout << "stop" << endl;	aigstop = 1; return; }
				if (!(op.upto4)) return;
			}
		}
	}
	sp6 = s6;
	sp6.possible_cells = twu[0];
	sp6.v= tuv128.v0;
	memset(ntbelow, 0, sizeof ntbelow);//7 8 9 10 11 full
	memset(tandbelow, 255, sizeof tandbelow);//7 8 9 10 11 full


next7:	//_______ add clue 7
	{
		uint64_t p = sp6.possible_cells;
		if (!p) {	
			if (p_cpt2g[4] == op.f4) cout << " 7 8 9 10 " << ntbelow[0] << " " << ntbelow[1] << " "
				<< ntbelow[2] << " " << ntbelow[3] << endl;

			EndExpand_7_9();	return;	}
		bitscanforward64(tcells[0], p);
		register uint64_t bit = (uint64_t)1 << tcells[0];
		sp6.possible_cells ^= bit;
		sp6.active_cells ^= bit;
		sp7 = sp6; 
		sp7.all_previous_cells |= bit;
		if (op.known > 1) {
			if (knownt >= 9)return;
			if (!((~pk54) & sp7.all_previous_cells))
				cout << Char54out(sp7.all_previous_cells) << " expected 7 " << endl;
		}
		if(locdiag)cout << Char54out(sp7.all_previous_cells) << "  7 " << endl;
		// update vectors and get next ua
		register uint64_t U = 0;
		sp7.v &= tuv128.vc[tcells[0]];
		if (sp7.v.isNotEmpty()) U = twu[sp7.v.getFirst128()	];
		else {
			if (t54b12.ncblocs) {
				register T54B12::TUVECT& tv1 = t54b12.tc128[1];
				BF128 w1=tv1.v0 & tv1.vc[tcells[0]];
				if (w1.isNotEmpty()) U = tv1.t[w1.getFirst128()];
			}
		}
		if (!U) {// this is a valid 7 (4 known cases
			if (!IsValid7pbf(sp7.all_previous_cells)) {
				BF128 w;// valid 7
				w.bf.u64[0] = sp7.all_previous_cells;
				w.bf.u64[1] = sp7.active_cells ;
				tbelow7[ntbelow[0]++] = w;
				tandbelow[0] &= w.bf.u64[0];
				p_cpt2g[54]++;
				goto next7;
			}
			else U = ua_ret7p;
		}
		sp7.possible_cells = U& sp7.active_cells;// never empty
	}

next8: //add clue 8 and find valid below size 8
	{
		uint64_t p = sp7.possible_cells;
		if (!p) goto next7;
		bitscanforward64(tcells[1], p);
		register uint64_t bit = (uint64_t)1 << tcells[1];
		sp7.possible_cells ^= bit;
		sp7.active_cells ^= bit;
		sp8 = sp7; 
		sp8.all_previous_cells |= bit;
		if (op.known > 1) {
			if (knownt >= 9)return;
			if (!((~pk54) & sp8.all_previous_cells))
				cout << Char54out(sp8.all_previous_cells) << " expected 8 " << endl;
		}
		if(locdiag)cout << Char54out(sp8.all_previous_cells) << "  8 " << endl;

		register uint64_t U = 0;
		sp8.v &= tuv128.vc[tcells[1]];
		if (sp8.v.isNotEmpty()) U = twu[sp8.v.getFirst128()];
		else {
			if (t54b12.ncblocs) {
				register T54B12::TUVECT& tv1 = t54b12.tc128[1];
				BF128 w1 = tv1.v0 & 
					(tv1.vc[tcells[1]]&tv1.vc[tcells[0]]);
				if (w1.isNotEmpty()) U = tv1.t[w1.getFirst128()];
				else if (t54b12.ncblocs > 1) {// safety code not likely
					register T54B12::TUVECT& tv2 = t54b12.tc128[2];
				register T54B12::TUVECT& tv1 = t54b12.tc128[1];
					BF128 w2 = tv2.v0 &
						( tv2.vc[tcells[1]] & tv2.vc[tcells[0]]);
					if (w2.isNotEmpty()) U = tv2.t[w2.getFirst128()];
					// no chance to need more than 384 uas here
				}
			}
		}
		if (!U) {// this is a valid 8 not often
			if (!IsValid7pbf(sp8.all_previous_cells)) {
				BF128 w;// valid 8
				w.bf.u64[0] = sp8.all_previous_cells;
				w.bf.u64[1] = sp8.active_cells;
				tbelow8[ntbelow[1]++] = w;
				tandbelow[1] &= w.bf.u64[0];
				p_cpt2g[55]++;
				goto next8;
			}
			else U = ua_ret7p;
		}
		sp8.possible_cells = U& sp8.active_cells;// never empty
	}
	 
next9: // add clue 9 build td and call next step (min 10 clues)
	{
		uint64_t p = sp8.possible_cells;
		if (!p) goto next8;
		p_cpt2g[5]++;
		bitscanforward64(tcells[2], p);
		register uint64_t bit = (uint64_t)1 << tcells[2];
		sp8.possible_cells ^= bit;
		sp8.active_cells ^= bit;
		sp9 = sp8; 
		sp9.all_previous_cells |= bit;
		sp9.v &= tuv128.vc[tcells[2]];
		if (op.known > 1) {
			if (knownt >= 9)return;
			if (!((~pk54) & sp9.all_previous_cells)) {
				cout << Char54out(sp9.all_previous_cells) << " expected 9 " << endl;
				knownt = 9;
				if (op.ton > 1) {
					t54b12.DebugC();
					sp9.v.Print(" v ");
				}
			}
		}
		if (locdiag)cout << Char54out(sp9.all_previous_cells) << "  9 " << endl;

		register uint64_t U = 0;
		if (sp9.v.isNotEmpty()) U = twu[sp9.v.getFirst128()];
		else {
			if (t54b12.ncblocs) {
				register T54B12::TUVECT& tv1 = t54b12.tc128[1];
				BF128 w1 = (tv1.v0 & tv1.vc[tcells[2]])
					& (tv1.vc[tcells[1]] & tv1.vc[tcells[0]]);
				if (w1.isNotEmpty()) U = tv1.t[w1.getFirst128()];
				else if (t54b12.ncblocs > 1) {// safety code not likely
					register T54B12::TUVECT& tv2 = t54b12.tc128[2];
					BF128 w2 = (tv2.v0 & tv2.vc[tcells[2]])
						& (tv2.vc[tcells[1]] & tv2.vc[tcells[0]]);
					if (w2.isNotEmpty()) U = tv2.t[w2.getFirst128()];
					else if (t54b12.ncblocs > 2) {// safety code not likely
						register T54B12::TUVECT& tv3 = t54b12.tc128[3];
						BF128 w3 = (tv3.v0 & tv3.vc[tcells[2]])
							& (tv3.vc[tcells[1]] & tv3.vc[tcells[0]]);
						if (w3.isNotEmpty()) U = tv3.t[w3.getFirst128()];
					}
					// no chance to need more than 512 uas here
				}
			}
		}
		if (!U) {// this can be a valid 9
			if (!IsValid7pbf(sp9.all_previous_cells)) {
				BF128 w;// valid 9
				w.bf.u64[0] = sp9.all_previous_cells;
				w.bf.u64[1] = sp9.active_cells;
				tbelow9[ntbelow[2]++] = w;
				tandbelow[2] &= w.bf.u64[0];
				p_cpt2g[56]++;
				goto next9;
			}
			else U = ua_ret7p;
		}
		sp9.possible_cells = U & sp9.active_cells;// never empty
		if (knownt == 9)cout << Char54out(sp9.possible_cells) << " after 9" << endl;
	}
	if (op.t18 && op.p2) {// see later 
		goto next9;
	}
	if ((!op.t18) && op.p1) {// 10 is the last clue
	next10last://this is the last
		{
			uint64_t p = sp9.possible_cells;
			if (!p) goto next9;
			bitscanforward64(tcells[3], p);
			register uint64_t bit = (uint64_t)1 << tcells[3];
			sp10.possible_cells ^= bit;
			sp10 = sp9;
			sp10.all_previous_cells |= bit;
			if (op.known > 8) {
				if (!((~pk54) & sp10.all_previous_cells)) {
					cout << Char54out(sp10.all_previous_cells) << " expected 10 " << endl;
				}
			}
			sp10.v &= tuv128.vc[tcells[3]];
			if (sp10.v.isNotEmpty())goto next10last;
			for (uint32_t i = 1; i <= t54b12.ncblocs; i++) {
				register T54B12::TUVECT& tvi = t54b12.tc128[i];
				BF128 wi = tvi.v0 & (tvi.vc[tcells[3]] & tvi.vc[tcells[2]]);
				wi &= (tvi.vc[tcells[1]] & tvi.vc[tcells[0]]);
				if (wi.isNotEmpty())goto next10last;
			}
			// this is a 10 to process
			goto next10last;
		}
	}
	// now 11 clues 48 p1  or 17 p2
next10:
	{
		uint64_t p = sp9.possible_cells;
		if (!p) goto next9;
		bitscanforward64(tcells[3], p);
		register uint64_t bit = (uint64_t)1 << tcells[3];
		sp9.possible_cells ^= bit;
		sp9.active_cells ^= bit;
		sp10 = sp9;
		sp10.all_previous_cells |= bit;
		if (op.known>1) {
			if (knownt > 9) return;
			if (knownt==9)cout << Char54out(sp10.all_previous_cells) << " 10" << endl;
			if (!((~pk54) & sp10.all_previous_cells)) {
				cout << Char54out(sp10.all_previous_cells) << " expected 10 " << endl;
				knownt = 10;
			}
		}
		// next if last, find and of remaining uas
		register uint64_t Uand = ~0,U=0;
		sp10.v &= tuv128.vc[tcells[3]];
		if (sp10.v.isNotEmpty()) {	U = 1;	tuv128.DoAnd(sp10.v, Uand);		}
		if (knownt == 10)cout << Char54out(Uand) << " after first bloc " << endl;
		if (U &&(!Uand)) goto next10;
		for (uint32_t i = 1; i <= t54b12.ncblocs; i++) {
			register T54B12::TUVECT& tvi = t54b12.tc128[i];
			BF128 wi = tvi.v0 & (tvi.vc[tcells[3]] & tvi.vc[tcells[2]]);
			wi &= (tvi.vc[tcells[1]] & tvi.vc[tcells[0]]);
			if (wi.isNotEmpty()) {	tvi.DoAnd(wi, Uand); U = 1;	}
			if (U && (!Uand)) goto next10;
		}

		if (knownt == 10)cout << Char54out(Uand) << " final U= "<<U << endl;

		if (!U) {//can be a valid 10 uaand empty
			//if (locdiag) cout << "can be a valid 10" << endl;
			if (!IsValid7pbf(sp10.all_previous_cells)) {
				BF128 w;// valid 10
				w.bf.u64[0] = sp10.all_previous_cells;
				w.bf.u64[1] = sp10.active_cells;
				tbelow10[ntbelow[3]++] = w;
				tandbelow[3] &= w.bf.u64[0];
				p_cpt2g[57]++;
				if (knownt == 10) { EndExpand_7_9(); return; }
				goto next10;
			}
			else Uand = anduab12;
		}
		if (knownt == 10)cout << Char54out(Uand) << " final U= " << U << endl;


		// this is a 10 to push to 11
		if (op.t18) {
			register uint64_t P = Uand;
			uint64_t cc = _popcnt64(sp10.all_previous_cells & BIT_SET_27);// band 1 count
			if (cc > 7 || cc < 2) goto next10;
			if (cc == 7)P &= ~(uint64_t)BIT_SET_27;
			if (cc == 3)P &= BIT_SET_27;
			Uand = P;
		}
		if (knownt == 10) {
			cout << Char54out(Uand) << " after bands filter " << endl;
			cout << Char54out(sp10.active_cells) << " active " << endl;
		}

		Uand &= sp10.active_cells;
		if(!Uand)goto next10;
		while (bitscanforward64(tcells[4], Uand)) {
			uint64_t bit2 = (uint64_t)1 << tcells[4];
			Uand ^= bit2; //clear bit
			myb12 = cb3.bf12 = sp10.all_previous_cells | bit2;
			if (knownt >= 10) {
				cout << Char54out(myb12) << " 11" << endl;
				if (!((~pk54) & myb12)) {
					cout << Char54out(myb12) << " expected 11 " << endl;
					knownt = 11;
				}
			}
			clean_valid_done = 0;
			p_cpt2g[7]++;
			cb3.ncl = 11;
			cb3.cbs.Init(myb12, 11);
			cb3.g2t = guah54_2.GetG2(cb3.bf12);
			cb3.g3t = guah54_2.GetG3(cb3.bf12);
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
				STD_B3& b3 = genb12.bands3[ib3];
				b3.Go(cb3);
			}
			if (knownt == 11) return;
		}
		goto next10;
	}
/*

	t54b12.Build_td128(sp9);
	if (t54b12.nd128 > p_cpt2g[59]) p_cpt2g[59] = t54b12.nd128;
	if (op.t18) {
		if (op.p1)ExpandTo11_18(s6, sp9);
		else;

	}
	else {
		if (op.p1)ExpandTo10(s6,sp9);
		else;
	}*/

	goto next9;

/*

if (t54b12.tandd)

			{// adjust active to the band constraint
				register uint64_t nb1 = _popcnt64(bf & BIT_SET_27),
					nb2 = 9 - nb1;
				if (op.p2) {
					if (nb1 > 6 || nb2 > 6) continue;;
					if (nb1 == 6) P &= ~(uint64_t)BIT_SET_27;
					if (nb2 == 6) P &= BIT_SET_27;
					w.bf.u64[1] = P;
				}
				else if ((!nb1) || (!nb2)) continue;
			}

	*/

}
void G17B::EndExpand_7_9() {
	int locdiag = 0;
	if (p_cpt2g[4] == op.f4) {
		cout << "call 7_9 good path end expand" << endl;
		locdiag = 1; 
	}

	if (op.t18) {
		if (op.p1) {
			if (knownt > 10)return;
			if (ntbelow[0]) Go_7_11_18(); // do 7 clues then more
			if (ntbelow[1]) Go_8_11_18(); // do 8 clues then more
			if (ntbelow[2]) Go_9_11_18(); // do 9 clues then more
			if (ntbelow[3]) Go_10_11_18(); // do 10 clues then more
		}
		else {
			if (knownt > 11)return;
			if (ntbelow[0]) Go_7_12(); // push to 12 clues 666
			if (ntbelow[1]) Go_8_12();
			if (ntbelow[2]) Go_9_12();

		}
	}
	else {
		if (op.p1) {
			if (ntbelow[0]) Go_7_11_18(); // do 7 clues then more
			if (ntbelow[1]) Go_8_11_18(); // do 8 clues then more
			if (ntbelow[2]) Go_9_11_18(); // do 9 clues then more
		}
		else {
			if (ntbelow[0]) Go_7_11_17(); // push to 11 clues 656
			if (ntbelow[1]) Go_8_11_17();
			if (ntbelow[2]) Go_9_11_17();
		}
	}
}



//________________expand and go search 17 pass 1

void G17B::Go_9_10() {// 9 clues limit 10 clues 
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[2]; iv++) {
		BF128 ww = tbelow9[iv];
		myb12 = cb3.bf12 = ww.bf.u64[0];
		myac= ww.bf.u64[1];
		guah54_9.Build9(myb12, myac);
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
			genb12.bands3[ib3].BuildGuam9(myb12);
		cb3.cbs.Init(myb12, 9);
		cb3.ncl = 9;
		// try direct
		cb3.g2t = guah54_9.GetG2(cb3.bf12);
		cb3.g3t = guah54_9.GetG3(cb3.bf12);
		p_cpt2g[7]++;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
			STD_B3& b3 = genb12.bands3[ib3];
			b3.Go(cb3);
		}
		// try now one more clue in bands 1+2
		uint64_t Ac = myac;
		{
			register uint64_t nb1 = cb3.cbs.b[0], nb2 = cb3.cbs.b[1];
			if (nb1 > 7 || nb2 > 7) continue; // should have been seen earlier
			if (nb1 == 7) Ac &= ~(uint64_t)BIT_SET_27; // only band 2
			if (nb2 == 7) Ac &= BIT_SET_27; // only band 2

		}
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.ncl = 10;
			cb3n.cbs.Add(cell);
			cb3n.g2t = guah54_9.GetG2(myb12);
			cb3n.g3t = guah54_9.GetG3(myb12);
			p_cpt2g[7]++;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
				genb12.bands3[ib3].Go(cb3n);
		}
	}
}
void G17B::Go_8_10() {// 8 clues limit 10 clues 
	return; // not ready
	clean_valid_done = 1;
	//SPB03* sn = &spb_0_15[7];
	for (uint32_t iv = 0; iv < ntbelow[1]; iv++) {
		BF128 ww = tbelow8[iv];
		cb3.ncl = 8;
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 8);

		// try direct
		cb3.g2t = guah54_2.GetG2(cb3.bf12);
		cb3.g3t = guah54_2.GetG3(cb3.bf12);
		p_cpt2g[7]++;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
			genb12.bands3[ib3].Go(cb3);

		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		{
			register uint64_t nb1 = cb3.cbs.b[0], nb2 = cb3.cbs.b[1];
			if (nb1 > 7 || nb2 > 7) continue; // should have been seen earlier
			if (nb1 == 7) Ac &= ~(uint64_t)BIT_SET_27; // only band 2
			if (nb2 == 7) Ac &= BIT_SET_27; // only band 2

		}

		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.ncl = 9;
			cb3n.cbs.Add(cell);
			cb3n.g2t = guah54_2.GetG2(myb12);
			cb3n.g3t = guah54_2.GetG3(myb12);
			p_cpt2g[7]++;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
				genb12.bands3[ib3].Go(cb3n);

			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			{
				register uint64_t nb1 = cb3n.cbs.b[0], nb2 = cb3n.cbs.b[1];
				if (nb1 > 7 || nb2 > 7) continue; // should have been seen earlier
				if (nb1 == 7) Ac &= ~(uint64_t)BIT_SET_27; // only band 2
				if (nb2 == 7) Ac &= BIT_SET_27; // only band 2

			}
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.ncl = 10;
				cb3n2.cbs.Add(cell2);
				cb3n2.g2t = guah54_2.GetG2(myb12);
				cb3n2.g3t = guah54_2.GetG3(myb12);
				p_cpt2g[7]++;
				for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
					genb12.bands3[ib3].Go(cb3n2);
			}
		}
	}
}
void G17B::Go_7_10() {// 7 clues limit 10 clues 
	return; // not ready copied from 8_10
	cout << " entry 7 clues for 10 clues" << endl;
	clean_valid_done = 1;
	//SPB03* sn = &spb_0_15[7];
	for (uint32_t iv = 0; iv < ntbelow[1]; iv++) {
		BF128 ww = tbelow8[iv];
		cb3.ncl = 8;
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 7);
		// try direct (expected empty 
		cb3.g2t = guah54_2.GetG2(cb3.bf12);
		cb3.g3t = guah54_2.GetG3(cb3.bf12);
		p_cpt2g[7]++;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
			genb12.bands3[ib3].Go(cb3);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		{
			register uint64_t nb1 = cb3.cbs.b[0], nb2 = cb3.cbs.b[1];
			if (nb1 > 7 || nb2 > 7) continue; // should have been seen earlier
			if (nb1 == 7) Ac &= ~(uint64_t)BIT_SET_27; // only band 2
			if (nb2 == 7) Ac &= BIT_SET_27; // only band 2

		}
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.ncl = 9;
			cb3n.cbs.Add(cell);
			cb3n.g2t = guah54_2.GetG2(myb12);
			cb3n.g3t = guah54_2.GetG3(myb12);
			p_cpt2g[7]++;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
				genb12.bands3[ib3].Go(cb3n);

			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			{
				register uint64_t nb1 = cb3n.cbs.b[0], nb2 = cb3n.cbs.b[1];
				if (nb1 > 7 || nb2 > 7) continue; // should have been seen earlier
				if (nb1 == 7) Ac2 &= ~(uint64_t)BIT_SET_27; // only band 2
				if (nb2 == 7) Ac2 &= BIT_SET_27; // only band 2

			}
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.ncl = 10;
				cb3n2.cbs.Add(cell2);
				cb3n2.g2t = guah54_2.GetG2(myb12);
				cb3n2.g3t = guah54_2.GetG3(myb12);
				p_cpt2g[7]++;
				for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
					genb12.bands3[ib3].Go(cb3n2);
				// try now a third clue in bands 1+2
				uint64_t Ac3 = Ac2;// others are not active now
				{
					register uint64_t nb1 = cb3n2.cbs.b[0], nb2 = cb3n2.cbs.b[1];
					if (nb1 > 7 || nb2 > 7) continue; // should have been seen earlier
					if (nb1 == 7) Ac3 &= ~(uint64_t)BIT_SET_27; // only band 2
					if (nb2 == 7) Ac3 &= BIT_SET_27; // only band 2

				}
				int cell3;
				while (bitscanforward64(cell3, Ac3)) {
					CALLBAND3 cb3n3 = cb3n2;
					uint64_t bit3 = (uint64_t)1 << cell3;
					Ac2 ^= bit2; //clear bit
					cb3n3.bf12 |= bit3;
					myb12 = cb3n3.bf12;
					cb3n3.ncl = 10;
					cb3n3.cbs.Add(cell3);
					cb3n2.g2t = guah54_2.GetG2(myb12);
					cb3n2.g3t = guah54_2.GetG3(myb12);
					p_cpt2g[7]++;
					for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
						genb12.bands3[ib3].Go(cb3n3);

				}

			}
		}
	}
}
/*
void G17B::ExpandTo10(SPB03A& s6, SPB03B& s9) {//search 17 last clue and below
	// do all below first
	int nx = ntbelow[0] + ntbelow[1] + ntbelow[2];
	if ((!nx) && (!t54b12.tandd)) return;
	return;//<<<<<<<<<<<<<<<<<<<<<<<<<<
		//nclues = 6;// see why needed
	int locdiag = 0;
	p_cpt2g[6]++;
	myandall = tandbelow[0] & tandbelow[1] & tandbelow[2] & s9.all_previous_cells;

	if (sgo.bfx[3] & 1) return;// debugging phase1
	// process 10 clues  (full)
	if (ntbelow[5] == 1) {// get active g2 g3 from guah54_2 direct
		//guah54_2.Dumpall2();
		p_cpt2g[7]++;
		cb3.ncl = 10;
		myb12 = cb3.bf12 = tfull[0];
		cb3.cbs.Init(myb12, 10);
		cb3.g2t = guah54_2.GetG2(cb3.bf12);

		cb3.g3t = guah54_2.GetG3(cb3.bf12);
		clean_valid_done = 0;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
			STD_B3& b3 = genb12.bands3[ib3];
			b3.Go(cb3);
			if (clean_valid_done == 2)break;
		}

	}
}
*/

//_________________________
uint64_t G17B::GetLastD() {
	register uint64_t   And = ~0;
	for (uint32_t i = 0; i <= t54b12.ndblocs; i++) {
		BF128 w = dvect[i];
		if (w.isEmpty()) continue;
		int ir;
		T54B12::TUVECT& vv = t54b12.td128[i];
		while ((ir = w.getFirst128()) >= 0) {
			w.clearBit(ir);
			And &= vv.t[ir];
		}
	}
	return And;
}


//________________expand and go search 18 pass 1


void G17B::Go_10_11_18() {// 10 clues limit 11 clues 
	int locdiag = 0;
	if (p_cpt2g[4] == op.f4) {
		cout << "go 10 to 11 18 " << ntbelow[3] << endl;
		locdiag = 1;
	}


	//cout << "valid 10 to go to 11 status ntbelow[3]="<< ntbelow[3] << endl;
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[3]; iv++) {
		BF128 ww = tbelow10[iv];
		cb3.ncl = 10;
		myb12 = cb3.bf12 = ww.bf.u64[0];
		if (op.known >1) {
			if (!((~pk54) & myb12)) {
				cout << Char54out(myb12) << " expected 10 to push to 11 " << endl;
			}
		}
		if (locdiag) {
			cout << Char54out(myb12) << " to push to 11 iv=" << iv << endl;

		}

		//cout << Char54out(myb12) << "10 to test" << endl;
		cb3.cbs.Init(myb12, 10);
		// try direct
		cb3.g2t = guah54_2.GetG2(cb3.bf12);
		cb3.g3t = guah54_2.GetG3(cb3.bf12);
		p_cpt2g[7]++;
		// can be here 288 or 828
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
			genb12.bands3[ib3].Go(cb3);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		{
			register uint64_t nb1 = cb3.cbs.b[0], nb2 = cb3.cbs.b[1];
			if (nb1 > 7 || nb2 > 7) continue;
			if (nb1 == 7) Ac &= ~(uint64_t)BIT_SET_27;
			if (nb2 == 7) Ac &= BIT_SET_27;
		}
		if (op.known > 1) {
			if (!((~pk54) & myb12)) {
				cout << Char54out(Ac) << "  active to try " << endl;
			}
		}
		if (locdiag) {
			cout << Char54out(Ac) << "  active to try " << endl;

		}
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.ncl = 11;
			cb3n.cbs.Add(cell);
			cb3n.g2t = guah54_2.GetG2(myb12);
			cb3n.g3t = guah54_2.GetG3(myb12);
			p_cpt2g[7]++;
			if (op.known > 1) {
				if (!((~pk54) & myb12)) {
					cout << Char54out(myb12) << "  expected 11 [7] "
						<< p_cpt2g[7] << endl;
					knownt = 11;
				}
			}
			if (locdiag) 
				cout << Char54out(myb12) << " 11 [7] "	<< p_cpt2g[7] << endl;

			

			for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
				genb12.bands3[ib3].Go(cb3n);
			if (knownt >= 11)return;
		}
	}
}
void G17B::Go_9_11_18() {// 9 clues limit 11 clues 
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[2]; iv++) {
		BF128 ww = tbelow9[iv];
		cb3.ncl = 9;
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 9);
		// try direct
		cb3.g2t = guah54_2.GetG2(cb3.bf12);
		cb3.g3t = guah54_2.GetG3(cb3.bf12);
		p_cpt2g[7]++;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
			genb12.bands3[ib3].Go(cb3);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		int cell;// maxi is 7+2 possible 8+2
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.ncl = 10;
			cb3n.cbs.Add(cell);
			cb3n.g2t = guah54_2.GetG2(myb12);
			cb3n.g3t = guah54_2.GetG3(myb12);
			// could be 288 828
			p_cpt2g[7]++;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
				genb12.bands3[ib3].Go(cb3n);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			{
				register uint64_t nb1 = cb3n.cbs.b[0], nb2 = cb3n.cbs.b[1];
				if (nb1 > 7 || nb2 > 7) continue;
				if (nb1 == 7) Ac2 &= ~(uint64_t)BIT_SET_27;
				if (nb2 == 7) Ac2 &= BIT_SET_27;
			}
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				p_cpt2g[7]++;
				cb3n2.ncl = 11;
				cb3n2.cbs.Add(cell2);
				cb3n2.g2t = guah54_2.GetG2(myb12);
				cb3n2.g3t = guah54_2.GetG3(myb12);
				for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
					genb12.bands3[ib3].Go(cb3n2);
			}
		}
	}

}
void G17B::Go_8_11_18() {// 8 clues limit 11 clues 
	//cout << " entry 8 clues for 11 clues" << endl;
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[1]; iv++) {
		BF128 ww = tbelow8[iv];
		cb3.ncl = 8;
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 8);
		// try direct
		cb3.g2t = guah54_2.GetG2(cb3.bf12);
		cb3.g3t = guah54_2.GetG3(cb3.bf12);
		p_cpt2g[7]++;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
			genb12.bands3[ib3].Go(cb3);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		int cell;// maxi is 6+2 possible 8+2
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.ncl = 9;
			cb3n.cbs.Add(cell);
			cb3n.g2t = guah54_2.GetG2(myb12);
			cb3n.g3t = guah54_2.GetG3(myb12);
			p_cpt2g[7]++;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
				genb12.bands3[ib3].Go(cb3n);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.ncl = 10;
				cb3n2.cbs.Add(cell2);
				cb3n2.g2t = guah54_2.GetG2(myb12);
				cb3n2.g3t = guah54_2.GetG3(myb12);
				// could be 288 828
				for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
					genb12.bands3[ib3].Go(cb3n2);
				// try now a third clue in bands 1+2
				uint64_t Ac3 = Ac2;// others are not active now
				{
					register uint64_t nb1 = cb3n2.cbs.b[0], nb2 = cb3n2.cbs.b[1];
					if (nb1 > 7 || nb2 > 7) continue;
					if (nb1 == 7) Ac3 &= ~(uint64_t)BIT_SET_27;
					if (nb2 == 7) Ac3 &= BIT_SET_27;
				}
				int cell3;
				while (bitscanforward64(cell3, Ac3)) {
					CALLBAND3 cb3n3 = cb3n2;
					uint64_t bit3 = (uint64_t)1 << cell3;
					Ac3 ^= bit3; //clear bit
					cb3n3.bf12 |= bit3;
					myb12 = cb3n3.bf12;
					cb3n3.ncl = 11;
					cb3n3.cbs.Add(cell3);
					cb3n3.g2t = guah54_2.GetG2(myb12);
					cb3n3.g3t = guah54_2.GetG3(myb12);
					for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
						genb12.bands3[ib3].Go(cb3n3);
				}
			}
		}
	}

}
void G17B::Go_7_11_18() {// 7 clues limit 11 clues 
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[0]; iv++) {
		BF128 ww = tbelow7[iv];
		cb3.ncl = 7;
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 7);
		// try direct
		cb3.g2t = guah54_2.GetG2(cb3.bf12);
		cb3.g3t = guah54_2.GetG3(cb3.bf12);
		p_cpt2g[7]++;
		//cout << Char54out(myb12) << " 7 direct" << endl;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
			genb12.bands3[ib3].Go(cb3);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		int cell;// maxi is 5+2 possible 8+2
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.ncl = 8;
			cb3n.cbs.Add(cell);
			cb3n.g2t = guah54_2.GetG2(myb12);
			cb3n.g3t = guah54_2.GetG3(myb12);
			p_cpt2g[7]++;
			//cout << Char54out(myb12) << " 7 to 8" << endl;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
				genb12.bands3[ib3].Go(cb3n);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.ncl = 9;
				cb3n2.cbs.Add(cell2);
				cb3n2.g2t = guah54_2.GetG2(myb12);
				cb3n2.g3t = guah54_2.GetG3(myb12);
				// could be 288 828
				p_cpt2g[7]++;
				//cout << Char54out(myb12) << " 7 to 9" << endl;
				for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
					genb12.bands3[ib3].Go(cb3n2);
				// try now a third clue in bands 1+2
				uint64_t Ac3 = Ac2;
				int cell3;
				while (bitscanforward64(cell3, Ac3)) {
					CALLBAND3 cb3n3 = cb3n2;
					uint64_t bit3 = (uint64_t)1 << cell3;
					Ac3 ^= bit3; //clear bit
					cb3n3.bf12 |= bit3;
					myb12 = cb3n3.bf12;
					cb3n3.ncl = 10;
					cb3n3.cbs.Add(cell3);
					cb3n3.g2t = guah54_2.GetG2(myb12);
					cb3n3.g3t = guah54_2.GetG3(myb12);
					p_cpt2g[7]++;
					//cout << Char54out(myb12) << " 7 to 10" << endl;
					for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
						genb12.bands3[ib3].Go(cb3n3);
					// try now a fourth  clue in bands 1+2
					uint64_t Ac4 = Ac3;
					{
						register uint64_t nb1 = cb3n3.cbs.b[0], nb2 = cb3n3.cbs.b[1];
						if (nb1 > 7 || nb2 > 7) continue;
						if (nb1 == 7) Ac4 &= ~(uint64_t)BIT_SET_27;
						if (nb2 == 7) Ac4 &= BIT_SET_27;
					}
					int cell4;
					while (bitscanforward64(cell4, Ac4)) {
						CALLBAND3 cb3n4 = cb3n3;
						uint64_t bit4 = (uint64_t)1 << cell4;
						Ac4 ^= bit4; //clear bit
						cb3n4.bf12 |= bit4;
						myb12 = cb3n4.bf12;
						cb3n4.ncl = 11;
						cb3n4.cbs.Add(cell4);
						cb3n4.g2t = guah54_2.GetG2(myb12);
						cb3n4.g3t = guah54_2.GetG3(myb12);
						p_cpt2g[7]++;
						//cout << Char54out(myb12) << " 7 to 11" << endl;
						for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
							genb12.bands3[ib3].Go(cb3n4);
					}
				}
			}
		}
	}
}


//________________ expand and go search 17 pass 2 

void G17B::Go_10_11_17() {// 10 clues limit 11 clues 
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[3]; iv++) {
		BF128 ww = tbelow10[iv];
		cb3.ncl = 10;
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 10);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		{
			register uint64_t nb1 = cb3.cbs.b[0], nb2 = cb3.cbs.b[1];
			if (nb1 > 6 || nb2 > 6) continue; // should have been seen earlier
			// 656,566 in p2a 566 in p2b
			if (nb1 == 6) Ac &= ~(uint64_t)BIT_SET_27;
			if (nb2 == 6) Ac &= BIT_SET_27;
			if (op.p2b) {
				if (nb2 > 5) continue;
				if (nb2 == 5)  Ac &= BIT_SET_27;
			}
		}
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.ncl = 11;
			cb3n.cbs.Add(cell);
			cb3n.g2t = guah54_2.GetG2(myb12);
			cb3n.g3t = guah54_2.GetG3(myb12);
			p_cpt2g[7]++;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
				genb12.bands3[ib3].Go(cb3n);

		}
	}

}
void G17B::Go_9_11_17() {// 9 clues limit 11 clues 
	//cout << " entry 9 clues for 11 clues" << endl;
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[2]; iv++) {
		BF128 ww = tbelow9[iv];
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 9);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			//cb3n.ncl = 10;
			cb3n.cbs.Add(cell);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			{
				register uint64_t nb1 = cb3n.cbs.b[0], nb2 = cb3n.cbs.b[1];
				if (nb1 > 6 || nb2 > 6) continue;
				// 656,566 in p2a 566 in p2b
				if (nb1 == 6) Ac2 &= ~(uint64_t)BIT_SET_27;
				if (nb2 == 6) Ac2 &= BIT_SET_27;
				if (op.p2b) {
					if (nb2 > 5) continue;
					if (nb2 == 5)  Ac2 &= BIT_SET_27;
				}
			}
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.ncl = 11;
				cb3n2.cbs.Add(cell2);
				cb3n2.g2t = guah54_2.GetG2(myb12);
				cb3n2.g3t = guah54_2.GetG3(myb12);
				p_cpt2g[7]++;
				for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
					genb12.bands3[ib3].Go(cb3n2);
			}

		}
	}
}
void G17B::Go_8_11_17() {// 8 clues limit 11 clues 
	//cout << " entry 8 clues for 11 clues" << endl;
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[1]; iv++) {
		BF128 ww = tbelow8[iv];
		//cb3.ncl = 8;
		myb12 = cb3.bf12 = ww.bf.u64[0];
		//if (t54b12.NotValid(myb12)) continue;// check more uas b12
		cb3.cbs.Init(myb12, 8);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.cbs.Add(cell);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			{
				register uint64_t nb1 = cb3n.cbs.b[0], nb2 = cb3n.cbs.b[1];
				if (nb1 > 6 || nb2 > 6) continue;
				// 656,566 in p2a 566 in p2b
				if (nb1 == 6) Ac2 &= ~(uint64_t)BIT_SET_27;
				if (nb2 == 6) Ac2 &= BIT_SET_27;
				if (op.p2b) {
					if (nb2 > 5) continue;
					if (nb2 == 5)  Ac2 &= BIT_SET_27;
				}
			}
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.cbs.Add(cell2);
				// try now a third clue in bands 1+2
				uint64_t Ac3 = Ac2;// others are not active now
				{
					register uint64_t nb1 = cb3n2.cbs.b[0], nb2 = cb3n2.cbs.b[1];
					if (nb1 > 6 || nb2 > 6) continue;
					// 656,566 in p2a 566 in p2b
					if (nb1 == 6) Ac3 &= ~(uint64_t)BIT_SET_27;
					if (nb2 == 6) Ac3 &= BIT_SET_27;
					if (op.p2b) {
						if (nb2 > 5) continue;
						if (nb2 == 5)  Ac3 &= BIT_SET_27;
					}
				}
				int cell3;
				while (bitscanforward64(cell3, Ac3)) {
					CALLBAND3 cb3n3 = cb3n2;
					uint64_t bit3 = (uint64_t)1 << cell3;
					Ac3 ^= bit3; //clear bit
					cb3n3.bf12 |= bit3;
					myb12 = cb3n3.bf12;
					cb3n3.ncl = 11;
					cb3n3.cbs.Add(cell2);
					cb3n3.g2t = guah54_2.GetG2(myb12);
					cb3n3.g3t = guah54_2.GetG3(myb12);
					p_cpt2g[7]++;
					for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
						genb12.bands3[ib3].Go(cb3n3);
				}
			}
		}
	}
}
void G17B::Go_7_11_17() {
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[0]; iv++) {
		BF128 ww = tbelow7[iv];
		//cb3.ncl = 8;
		myb12 = cb3.bf12 = ww.bf.u64[0];
		//if (t54b12.NotValid(myb12)) continue;// check more uas b12
		cb3.cbs.Init(myb12, 7);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.cbs.Add(cell);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			{
				register uint64_t nb1 = cb3n.cbs.b[0], nb2 = cb3n.cbs.b[1];
				if (nb1 > 6 || nb2 > 6) continue;
				// 656,566 in p2a 566 in p2b
				if (nb1 == 6) Ac2 &= ~(uint64_t)BIT_SET_27;
				if (nb2 == 6) Ac2 &= BIT_SET_27;
				if (op.p2b) {
					if (nb2 > 5) continue;
					if (nb2 == 5)  Ac2 &= BIT_SET_27;
				}
			}
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.cbs.Add(cell2);
				// try now a third clue in bands 1+2
				uint64_t Ac3 = Ac2;// others are not active now
				{
					register uint64_t nb1 = cb3n2.cbs.b[0], nb2 = cb3n2.cbs.b[1];
					if (nb1 > 6 || nb2 > 6) continue;
					// 656,566 in p2a 566 in p2b
					if (nb1 == 6) Ac3 &= ~(uint64_t)BIT_SET_27;
					if (nb2 == 6) Ac3 &= BIT_SET_27;
					if (op.p2b) {
						if (nb2 > 5) continue;
						if (nb2 == 5)  Ac3 &= BIT_SET_27;
					}
				}
				int cell3;
				while (bitscanforward64(cell3, Ac3)) {
					CALLBAND3 cb3n3 = cb3n2;
					uint64_t bit3 = (uint64_t)1 << cell3;
					Ac3 ^= bit3; //clear bit
					cb3n3.bf12 |= bit3;
					myb12 = cb3n3.bf12;
					cb3n3.cbs.Add(cell2);
					cb3n3.g2t = guah54_2.GetG2(myb12);
					cb3n3.g3t = guah54_2.GetG3(myb12);
					cb3n2.cbs.Add(cell2);
					// try now a third clue in bands 1+2
					uint64_t Ac4 = Ac3;// others are not active now
					{
						register uint64_t nb1 = cb3n3.cbs.b[0], nb2 = cb3n3.cbs.b[1];
						if (nb1 > 6 || nb2 > 6) continue;
						// 656,566 in p2a 566 in p2b
						if (nb1 == 6) Ac4 &= ~(uint64_t)BIT_SET_27;
						if (nb2 == 6) Ac4 &= BIT_SET_27;
						if (op.p2b) {
							if (nb2 > 5) continue;
							if (nb2 == 5)  Ac4 &= BIT_SET_27;
						}
					}
					int cell4;
					while (bitscanforward64(cell4, Ac3)) {
						CALLBAND3 cb3n4 = cb3n3;
						uint64_t bit4 = (uint64_t)1 << cell4;
						Ac4 ^= bit4; //clear bit
						cb3n4.bf12 |= bit4;
						myb12 = cb3n4.bf12;
						cb3n4.ncl = 11;
						cb3n4.cbs.Add(cell2);
						cb3n3.g2t = guah54_2.GetG2(myb12);
						cb3n3.g3t = guah54_2.GetG3(myb12);
						p_cpt2g[7]++;
						for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
							genb12.bands3[ib3].Go(cb3n4);
					}
				}
			}
		}
	}
}

//________________expand and go search 18 pass 2


void G17B::Go_11_12() {
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[4]; iv++) {
		BF128 ww = tbelow11[iv];
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 11);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		Ac &= cb3.cbs.NextActive();
		int locdiag = 0;
		if (op.known) {
			if (!((~pk54) & myb12)) {
				cout << Char54out(myb12) << " expected 12 go 11 12 iv=" << iv << endl;
				cout << Char54out(Ac) << " active 666 ";
				cb3.cbs.Status();
				locdiag = 1;
			}
		}
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.cbs.Add(cell);
			cb3n.ncl = 12;
			cb3n.g2t = guah54_9.GetG2(myb12);
			cb3n.g3t = guah54_9.GetG3(myb12);
			p_cpt2g[7]++;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
				genb12.bands3[ib3].Go(cb3n);
		}
	}
}
void G17B::Go_10_12() {
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[3]; iv++) {
		BF128 ww = tbelow10[iv];
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 10);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.cbs.Add(cell);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			Ac2 &= cb3n.cbs.NextActive();
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.cbs.Add(cell2);
				cb3n2.ncl = 12;
				cb3n2.g2t = guah54_9.GetG2(myb12);
				cb3n2.g3t = guah54_9.GetG3(myb12);
				p_cpt2g[7]++;
				for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
					genb12.bands3[ib3].Go(cb3n2);
			}
		}
	}
}
void G17B::Go_9_12() {
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[2]; iv++) {
		BF128 ww = tbelow9[iv];
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 9);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.cbs.Add(cell);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			{
				register uint64_t nb1 = cb3n.cbs.b[0], nb2 = cb3n.cbs.b[1];
				if (nb1 > 6 || nb2 > 6) continue;
				if (nb1 == 6) Ac2 &= ~(uint64_t)BIT_SET_27;
				if (nb2 == 6) Ac2 &= BIT_SET_27;
			}
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.cbs.Add(cell2);
				// try now a third clue in bands 1+2
				uint64_t Ac3 = Ac2;// others are not active now
				{
					{
						register uint64_t nb1 = cb3n2.cbs.b[0], nb2 = cb3n2.cbs.b[1];
						if (nb1 > 6 || nb2 > 6) continue;
						if (nb1 == 6) Ac3 &= ~(uint64_t)BIT_SET_27;
						if (nb2 == 6) Ac3 &= BIT_SET_27;
					}
				}
				int cell3;
				while (bitscanforward64(cell3, Ac3)) {
					CALLBAND3 cb3n3 = cb3n2;
					uint64_t bit3 = (uint64_t)1 << cell3;
					Ac3 ^= bit3; //clear bit
					cb3n3.bf12 |= bit3;
					myb12 = cb3n3.bf12;
					cb3n3.ncl = 12;
					cb3n3.cbs.Add(cell2);
					cb3n3.g2t = guah54_2.GetG2(myb12);
					cb3n3.g3t = guah54_2.GetG3(myb12);
					p_cpt2g[7]++;
					for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
						genb12.bands3[ib3].Go(cb3n3);
				}
			}
		}
	}
}
void G17B::Go_8_12() {
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[1]; iv++) {
		BF128 ww = tbelow8[iv];
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 8);
		uint64_t Ac = ww.bf.u64[1];
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.cbs.Add(cell);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			{
				register uint64_t nb1 = cb3n.cbs.b[0], nb2 = cb3n.cbs.b[1];
				if (nb1 > 6 || nb2 > 6) continue;
				if (nb1 == 6) Ac2 &= ~(uint64_t)BIT_SET_27;
				if (nb2 == 6) Ac2 &= BIT_SET_27;
			}
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.cbs.Add(cell2);
				// try now a third clue in bands 1+2
				uint64_t Ac3 = Ac2;// others are not active now
				{
					register uint64_t nb1 = cb3n2.cbs.b[0], nb2 = cb3n2.cbs.b[1];
					if (nb1 > 6 || nb2 > 6) continue;
					if (nb1 == 6) Ac3 &= ~(uint64_t)BIT_SET_27;
					if (nb2 == 6) Ac3 &= BIT_SET_27;
				}
				int cell3;
				while (bitscanforward64(cell3, Ac3)) {
					CALLBAND3 cb3n3 = cb3n2;
					uint64_t bit3 = (uint64_t)1 << cell3;
					Ac3 ^= bit3; //clear bit
					cb3n3.bf12 |= bit3;
					myb12 = cb3n3.bf12;
					cb3n3.cbs.Add(cell3);
					// try now a fourth  clue in bands 1+2
					uint64_t Ac4 = Ac3;
					{
						register uint64_t nb1 = cb3n3.cbs.b[0], nb2 = cb3n3.cbs.b[1];
						if (nb1 > 6 || nb2 > 6) continue;
						if (nb1 == 6) Ac4 &= ~(uint64_t)BIT_SET_27;
						if (nb2 == 6) Ac4 &= BIT_SET_27;
					}
					int cell4;
					while (bitscanforward64(cell4, Ac4)) {
						CALLBAND3 cb3n4 = cb3n3;
						uint64_t bit4 = (uint64_t)1 << cell4;
						Ac4 ^= bit4; //clear bit
						cb3n4.bf12 |= bit4;
						myb12 = cb3n4.bf12;
						cb3n4.ncl = 12;
						cb3n4.cbs.Add(cell4);
						cb3n4.g2t = guah54_2.GetG2(myb12);
						cb3n4.g3t = guah54_2.GetG3(myb12);
						p_cpt2g[7]++;
						for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
							genb12.bands3[ib3].Go(cb3n4);
					}
				}
			}
		}
	}
}
void G17B::Go_7_12() {

}


int G17B::IsValid7pbf(uint64_t bf) {
	if (zh2b[1].IsValid(bf)) {
		anduab12 = ~0;
		register uint64_t ua54;
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t ua = zh2gxn.tua[i], cc = _popcnt64(ua);
			ua54 = (ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
			anduab12 &= ua54;
			if (cc < 21) {
				t54b12.AddA(ua54);
				t54b12.AddB(ua54);
				t54b12.AddC(ua54);
			}
		}
		ua_ret7p = ua54;// return last (smaller)
		return 1;
	}
	return 0;
}
inline int G17B::GetNextCell(SPB03* s) {
	SPB03* sn;
	register int cell;
	uint64_t p = s->possible_cells;
	if (!p)return 1;
	bitscanforward64(cell, p);
	register uint64_t bit = (uint64_t)1 << cell;
	s->possible_cells ^= bit;
	tclues[s->ncl] = cell;
	s->active_cells ^= bit;
	sn = s + 1; *sn = *s; sn->ncl++;
	//sn->cbs.Add(cell);
	sn->all_previous_cells |= bit;
	sn->v[0] &= t54b12.tc128[0].vc[cell];
	return 0;
}
inline void G17B::GetNextUa(SPB03* sn) {
	register uint64_t  V;
	if ((V = sn->v[0].bf.u64[0])) {// next ua
		register uint32_t ir;
		bitscanforward64(ir, V);//relative index first active
		ua_ret7p = t54b12.tc128[0].t[ir];
	}
	else {// next ua must be here
		V = sn->v[0].bf.u64[1];
		register uint32_t ir;
		bitscanforward64(ir, V);//relative index first active
		ua_ret7p = t54b12.tc128[0].t[ir + 64];
	}
}
inline void G17B::GetNextUaAdd(SPB03* sn) {
	if (t54b12.nc128 <= 128) return;
	// more uas to check
	for (uint32_t i = 1; i <= t54b12.ncblocs; i++) {
		T54B12::TUVECT& vv = t54b12.tc128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 6; ic < sn->ncl; ic++)
			v &= vc[tclues[ic]];
		if (v.isNotEmpty()) {
			int ir2 = v.getFirst128();
			ua_ret7p = vv.t[ir2];
			return;
		}
	}
}
inline int G17B::GetLastAndUa(SPB03* sn, int d) {
	int aig = 0;
	register uint64_t  V, And = ~0;
	register uint32_t ir;
	if ((V = sn->v[0].bf.u64[0])) {
		aig = 1;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			And &= t54b12.tc128[0].t[ir];
		}
	}
	if (!And) return 1;
	if ((V = sn->v[0].bf.u64[1])) {
		aig = 1;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			And &= t54b12.tc128[0].t[ir + 64];
			if (!And) return 1;
		}
	}
	if (!And) return 1;
	if (t54b12.ncblocs) {
		// more uas to check
		for (uint32_t i = 1; i <= t54b12.ncblocs; i++) {
			T54B12::TUVECT& vv = t54b12.tc128[i];
			BF128 v = vv.v0, * vc = vv.vc;
			for (uint32_t ic = 6; ic < sn->ncl; ic++)
				v &= vc[tclues[ic]];
			if (v.isNotEmpty()) {
				aig = 1;
				int ir2;
				while ((ir2 = v.getFirst128()) >= 0) {
					v.clearBit(ir2);
					And &= vv.t[ir2];
				}
			}
		}
	}
	if (aig) ua_ret7p = And;
	return aig;
}

/// last step in bands 1+2 expansion critical code
/// expand builds tables of potential solutions 
/// one with target number of clues stored 
///		 bit field ;and of all items; number of clues/bands stacks)
/// one per number of clues below the target 
///      storing also the active cells downstream

void DumpPotential(int det=0) {
	cout <<  "DumpPotential() " ;
	for (int i = 0; i < 6; i++) cout << ntbelow[i] << " ";	
	cout << endl;
	if (!det) return;
	for (int i = 0; i < 6; i++)if (ntbelow[i]) {
		cout << Char54out(tandbelow[i]) << " and i" << i 
			<<" "<<_popcnt64(tandbelow[i]) << endl;
		if (det > 1) {
			if (i < 5) {
				BF128* tb[5] = { tbelow7 ,tbelow8,tbelow9,tbelow10, tbelow11 };
				cout << "tbelow nclues=<i+7 ncl=" << i + 7 << endl;
				BF128* tbw = tb[i];
				for (uint32_t j = 0; j < ntbelow[i]; j++)
					cout << Char54out(tbw[j].bf.u64[0]) << " " << j << endl;

			}
			else {
				cout << "t  nclues= full"   << endl;
				for (uint32_t j = 0; j < ntbelow[5]; j++)
					cout << Char54out(tfull[j]) << " " << j << endl;
			}
		}
	}

}

//void Check(uint64_t )
inline int G17B::GetNextCell7_9(SPB03* s) {
	SPB03* sn;
	register int cell;
	uint64_t p = s->possible_cells;
	if (!p)return 1;
	bitscanforward64(cell, p);
	register uint64_t bit = (uint64_t)1 << cell;
	s->possible_cells ^= bit;
	tclues[s->ncl] = cell;
	s->active_cells ^= bit;
	sn = s + 1; *sn = *s; sn->ncl++;
	//sn->cbs.Add(cell);
	sn->all_previous_cells |= bit;
	sn->v[0] &= t54b12.tc128[0].vc[cell];
	return 0;
}
/*

void G17B::GoExpand_7_9B() {// entry p_cpt2g[4 times ]
	SPB03* sn = &spb_0_15[7];
	//Expand_7_9B();
	if (p_cpt2g[72] < nt_7_9)p_cpt2g[72] = nt_7_9;
	int locdiag = 0;
	if (op.f4) {
		if (p_cpt2g[4] == op.f4) {
			cout << Char54out(sn->all_previous_cells) << "[4] good path expand 7_9 " << p_cpt2g[4] << endl;
			cout << "nt7_9" << nt_7_9 << " "; 
			DumpPotential(op.ton > 1);			
			locdiag = 1;
		}
		else {
			if (p_cpt2g[4] > op.f4) { cout << "stop" << endl;	aigstop = 1; return; }
			if (!(op.upto4)) return;
		}
	}
	else if (op.f3 && p_cpt2g[3] == op.f3) {
		cout << Char54out(sn->all_previous_cells) << "[4] " << p_cpt2g[4] << " " << "nt7_9" << nt_7_9 << endl;
	}

	if (sgo.bfx[3] & 2) 		return;// debugging phase1

	myandall = tandbelow[1] & tandbelow[2] & tand_7_9;
	guah54_2.Build2(myandall, sn->active_cells);
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
		genb12.bands3[ib3].BuildGuam2(myandall);

	// can be valid 7 or valid 8 

	if (ntbelow[0] ) { 
		if (op.t18) {
			if (op.p1)Go_7_11_18();
			else;// Go_7_12();
		}
		else {
			if (op.p1);// Go_7_10_17();
			else;// Go_7_11_17();
		}
	}
	if (ntbelow[1]) {
		if (op.t18) {
			if (op.p1) Go_8_11_18();
			else Go_8_12();
		}
		else {
			if (op.p1) Go_8_10();
			else Go_8_11_17();
		}
	}
	for (uint32_t i1 = 0; i1 < nt_7_9; i1++) {
		BF128 w9=t_7_9[i1];
		myb12_9 = w9.bf.u64[0];
		myac_9 = w9.bf.u64[1];
		if (op.known) {
			if (!((~pk54) & myb12_9)) {
				cout << Char54out(myb12_9) << " expected 9 [70]" << p_cpt2g[70] << endl;
				locdiag = 1;
			}
		}
		CBS cbs;
		//if (t54b12.Build_td128()) return;
		cbs.Init(myb12_9, 9);
	}

	//for (uint32_t i1 = 0; i1 < nt_7_9; i1++)
	//	GoExpand_10_12(t_7_9[i1]);
}
void G17B::GoExpand_10p(BF128 ww) {// from sp8 to end 
	// build 
}
void G17B::GoExpand_7_9() {
	SPB03* sn = &spb_0_15[7];
	p_cpt2g[70]++;
	int locdiag = 0;
	if (op.f4) {
		if (p_cpt2g[70] == op.f4) {
			cout << Char54out(sn->all_previous_cells) << "[4] good path expand 7_12 -v8- [70] " << p_cpt2g[70] << endl;
			locdiag = 1;
		}
		else {
			cout << Char54out(sn->all_previous_cells) << "[70] " << p_cpt2g[70] << endl;
			if (p_cpt2g[70] > op.f4) { cout << "stop" << endl;	aigstop = 1; return; }
			if (!(op.upto4)) return;
		}
	}
	Expand_7_9();

	if (ntbelow[0]) {
		cout << "unexpected 7 clues pass2 stop [4]" << p_cpt2g[4] << endl;
		aigstop = 1;		return;
	}
	p_cpt2g[71]++;
	if (p_cpt2g[72] < nt_7_9)p_cpt2g[72] = nt_7_9;

	if (sgo.bfx[3] & 2) 		return;// debugging phase1

	myandall = tandbelow[1] & tandbelow[2] & tand_7_9;
	guah54_2.Build2(myandall, sn->active_cells);
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
		genb12.bands3[ib3].BuildGuam2(myandall);

	if (locdiag && op.ton) {
		cout << "test 7_12 [71]" << p_cpt2g[71];
		cout << " nt_7_9=" << nt_7_9 << "\t\t ";
		DumpPotential(0);
		//cout << Char54out(myandall) << " myandall " << endl;
		//cout << Char54out(sn->active_cells) << " active " << endl;
		//cout << "guah54_2 status" << endl;
		//guah54_2.Dumpall2(); guah54_2.Dumpall3();
	}

	if (ntbelow[1]) Go_8_12();
	if (ntbelow[2]) Go_9_12();
	for (uint32_t i1 = 0; i1 < nt_7_9; i1++)
		GoExpand_10_12(t_7_9[i1]);
}
int G17B::Expand_7_10() {
	memset(ntbelow, 0, sizeof ntbelow);//7 8 9 10 11 full
	memset(tandbelow, 255, sizeof tandbelow);//7 8 9 10 11 full
	uint32_t nfull = 0, nbelow = 0;
	SPB03* sl = &spb_0_15[8], * s = sl, * sn;
	T54B12::TUVECT& tuv128 = t54b12.tc128[0];
	uint64_t* twu = tuv128.t;
	if (t54b12.nc128 < 128)t54b12.nc128 = 128;// force adds outside
	*s = spb_0_15[7];	// duplicate 6 for new vector
	s->possible_cells = twu[0];
	s->v[0] = tuv128.v0;// initial nothing done

	//_______ start search clues 7_10
next:	// catch and apply cell in bitfields
	if (GetNextCell(s))	if (--s >= sl)goto next;	
	else {	return (nfull | nbelow);	}
	sn = s + 1;
	ua_ret7p = 0;
	if (sn->ncl == 10) {// 10 cells
		tfull[ntbelow[5]++] = sn->all_previous_cells;
		tandbelow[5]&= sn->all_previous_cells;
		nfull++;
		goto next;
	}

	if (sn->ncl == 9) {// last step
		p_cpt2g[5]++;
		if (!GetLastAndUa(sn)) {// can be a valid 10 or unknown ua(s)
			if (IsValid7pbf(sn->all_previous_cells)) // got uas to use
				ua_ret7p = anduab12;
			else {// valid 9
				p_cpt2g[57]++;
				BF128 w;
				w.bf.u64[0] = sn->all_previous_cells;
				w.bf.u64[1] = sn->active_cells;
				tbelow9[ntbelow[2]++] = w;
				tandbelow[2] &= sn->all_previous_cells;
				nbelow++;
				goto next;
			}

		}
		{// last not empty direct or after check valid
			if (!ua_ret7p) goto next;
			register uint64_t nb1 = _popcnt64(sn->all_previous_cells & BIT_SET_27),
				nb2 = 9 - nb1, P = ua_ret7p & s->active_cells;
			if (nb1 > 7 || nb2 > 7) goto next;
			if (nb1 == 7) P &= ~(uint64_t)BIT_SET_27;
			if (nb2 == 7) P &= BIT_SET_27; 
			if (P) { sn->possible_cells = P; s++; }
			goto next;
		}

	}
	else {
		if (sn->v[0].isNotEmpty())GetNextUa(sn);// first 128
		else GetNextUaAdd(sn);
	}
	if (!ua_ret7p) {// add ua or go below
		if (IsValid7pbf(sn->all_previous_cells)) {// got uas to use
			if (sn->ncl == 10)ua_ret7p = anduab12;
		}
		else {// valid below
			clean_valid_done = 1;
			int isize = sn->ncl - 7;//table index 0 7 clues
			p_cpt2g[54+isize]++;
			BF128 w;
			w.bf.u64[0] = sn->all_previous_cells;
			w.bf.u64[1] = sn->active_cells;
			switch (isize) {
			case 0:   tbelow7[ntbelow[0]++] = w;
				tandbelow[0] &= sn->all_previous_cells; break;
			case 1:  tbelow8[ntbelow[1]++] = w;
				tandbelow[1] &= sn->all_previous_cells; break;
			case 2:  tbelow9[ntbelow[2]++] = w;
				tandbelow[2] &= sn->all_previous_cells; break;
			}			
			nbelow++;
			goto next;
		}
	}
	sn->possible_cells = ua_ret7p & s->active_cells;
	if (sn->possible_cells) s++;  // switch to next spot
	goto next;
}
void G17B::GoExpand_7_10() {
	if (op.t18) {
		if ( op.p1) { return; }
		else { GoExpand_7_12(); return; }
	}
	else if(op.p2) { return;  } 
	SPB03 * sn= &spb_0_15[7];
	//nclues = 6;// see why needed
	int ire = Expand_7_10();
	int locdiag = 0;
	if (op.ton > 1) {
		cout << Char54out(sn->all_previous_cells) << " 6clues [4]" << p_cpt2g[4] << endl;
		DumpPotential(0);
		if (p_cpt2g[3] == op.f3) {
			if (p_cpt2g[4] == op.f4) {
				cout << "this is the expected path" << endl;
				DumpPotential(2);
				if (op.ton > 2) { t54b12.DebugA(); t54b12.DebugB(); t54b12.DebugC(); }
				locdiag = op.ton;
			}
		}
	}
	if (!ire) return;


	p_cpt2g[6]++;
	int nbelow = ntbelow[4] + ntbelow[3] + ntbelow[2] + ntbelow[1] + ntbelow[0];
	p_cpt2g[51] += ntbelow[5];
	p_cpt2g[52] += nbelow;
	if (p_cpt2g[62] < ntbelow[5])p_cpt2g[62] = ntbelow[5];

	myandall = tandbelow[0] & tandbelow[1] & tandbelow[2] & tandbelow[3] & tandbelow[4] & tandbelow[5];

	if (sgo.bfx[3] & 1) return;// debugging phase1

	guah54_2.Build2(myandall, sn->active_cells);
	if (locdiag) {
		guah54.DumpOne3(76);
		guah54_2.DumpOne3(76);
	}
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
		genb12.bands3[ib3].BuildGuam2(myandall);

	if (ntbelow[0]) Go_7_10(); // do 7 clues then more
	if (ntbelow[1]) Go_8_10(); // do 8 clues then more
	if (ntbelow[2]) Go_9_10(); // do 9 clues then more

		// process 10 clues  (full)
	if (ntbelow[5] ==1) {// get active g2 g3 from guah54_2 direct
		//guah54_2.Dumpall2();
		p_cpt2g[7]++;
		cb3.ncl = 10;
		myb12 = cb3.bf12 = tfull[0];
		cb3.cbs.Init(myb12, 10);
		cb3.g2t = guah54_2.GetG2(cb3.bf12);

		cb3.g3t = guah54_2.GetG3(cb3.bf12);
		clean_valid_done = 0;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
			STD_B3& b3 = genb12.bands3[ib3];
			b3.Go(cb3);
			if(clean_valid_done == 2)break;
		}
	}
	else { 
		for (uint32_t iv = 0; iv < ntbelow[5]; iv ++ ) {
			cb3.ncl = 10;
			myb12 = cb3.bf12 = tfull[iv];
			cb3.cbs.Init(myb12, 10);
			p_cpt2g[7]++;
			if (locdiag) {
				cout << Char54out(myb12) << " [7] " << p_cpt2g[7] << endl;
				guah54_2.DumpOne3(76);
				if (p_cpt2g[7] == sgo.vx[8]) {
					cout << "this is the expected path to go b3" << endl;
					cb3.Dump();
					//guah54_2.Dumpall2();
					//guah54_2.Dumpall3();
				}

			}
			cb3.g2t = guah54_2.GetG2(cb3.bf12);
			cb3.g3t = guah54_2.GetG3(cb3.bf12);
			if (locdiag)guah54_2.DumpOne3(76);
			clean_valid_done = 0;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
				STD_B3& b3 = genb12.bands3[ib3];
				b3.Go(cb3);
				if (clean_valid_done == 2)break;
			}
		}
	}
}
*/
/// 7_11 is for 18 pass1 or 17 pass 2 (A and B) 
/// stack limit bands limit is 6 in pass 2
/// in pass 2B band 2 has 5 clues 

/*


//_________________________ processing 18 clues pass 2 666 666

inline int G17B::GetNextCellD(SPB03* s) {
	SPB03* sn;
	register int cell;
	uint64_t p = s->possible_cells;
	if (!p)return 1;
	bitscanforward64(cell, p);
	register uint64_t bit = (uint64_t)1 << cell;
	s->possible_cells ^= bit;
	tclues[s->ncl] = cell;
	s->active_cells ^= bit;
	sn = s + 1; *sn = *s; sn->ncl++;
	sn->cbs.Add(cell);
	sn->all_previous_cells |= bit;
	sn->v[0] &= t54b12.td128[0].vc[cell];
	return 0;
}
inline void G17B::GetNextUaD(SPB03* sn) {
	register uint64_t  V;
	if ((V = sn->v[0].bf.u64[0])) {// next ua
		register uint32_t ir;
		bitscanforward64(ir, V);//relative index first active
		ua_ret7p = t54b12.td128[0].t[ir];
	}
	else {// next ua must be here
		V = sn->v[0].bf.u64[1];
		register uint32_t ir;
		bitscanforward64(ir, V);//relative index first active
		ua_ret7p = t54b12.td128[0].t[ir + 64];
	}
}
inline void G17B::GetNextUaAddD(SPB03* sn) {
	if (t54b12.nd128 <= 128) return;
	// more uas to check
	for (uint32_t i = 1; i <= t54b12.ndblocs; i++) {
		T54B12::TUVECT& vv = t54b12.td128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 9; ic < sn->ncl; ic++)
			v &= vc[tclues[ic]];
		if (v.isNotEmpty()) {
			int ir2 = v.getFirst128();
			ua_ret7p = vv.t[ir2];
			return;
		}
	}
}
inline int G17B::GetLastAndUaD(SPB03* sn, int d) {
	int aig = 0;
	register uint64_t  V, And = sn->active_cells;
	register uint32_t ir;
	if ((V = sn->v[0].bf.u64[0])) {
		aig = 1;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			And &= t54b12.td128[0].t[ir];
		}
	}
	if (!And) return 1;
	if ((V = sn->v[0].bf.u64[1])) {
		aig = 1;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			And &= t54b12.td128[0].t[ir + 64];
			if (!And) return 1;
		}
	}
	if (!And) return 1;
	if (t54b12.ndblocs) {
		// more uas to check
		for (uint32_t i = 1; i <= t54b12.ndblocs; i++) {
			T54B12::TUVECT& vv = t54b12.td128[i];
			BF128 v = vv.v0, * vc = vv.vc;
			for (uint32_t ic = 9; ic < sn->ncl; ic++)
				v &= vc[tclues[ic]];
			if (v.isNotEmpty()) {
				aig = 1;
				int ir2;
				while ((ir2 = v.getFirst128()) >= 0) {
					v.clearBit(ir2);
					And &= vv.t[ir2];
				}
			}
		}
	}
	if (aig) ua_ret7p = And;
	return aig;
}



void G17B::Expand_7_9() {
	memset(ntbelow, 0, sizeof ntbelow);//7 8 9 10 11 full
	memset(tandbelow, 255, sizeof tandbelow);//7 8 9 10 11 full
	tand_7_9 = ~0;
	nt_7_9 = 0;
	SPB03* sl = &spb_0_15[8], * s = sl, * sn;
	T54B12::TUVECT& tuv128 = t54b12.tc128[0];
	uint64_t* twu = tuv128.t;
	if (t54b12.nc128 < 128)t54b12.nc128 = 128;// force adds outside
	*s = spb_0_15[7];	// duplicate 6 for new vector
	s->possible_cells = twu[0];
	s->v[0] = tuv128.v0;// initial nothing done


	//_______ start search clues 7_9
next:	// catch and apply cell in bitfields
	if (GetNextCell(s))	if (--s >= sl)goto next;	else 	return ;
	sn = s + 1;
	{// adjust active to the band constraint
		register uint64_t nb1 = _popcnt64(sn->all_previous_cells & BIT_SET_27),
			nb2 = sn->ncl - nb1, P =  sn->active_cells;
		if (nb1 > 6 || nb2 > 6) goto next;
		if (nb1 == 6) P &= ~(uint64_t)BIT_SET_27;
		if (nb2 == 6) P &= BIT_SET_27;
		sn->active_cells=P;
	}
	ua_ret7p = 0;
	if (sn->v[0].isNotEmpty())GetNextUa(sn);// first 128
	else GetNextUaAdd(sn);
	if (!ua_ret7p) {// add ua or go below
		if (!IsValid7pbf(sn->all_previous_cells)) {// valid below
			int isize = sn->ncl - 7;//table index 0 7 clues
			p_cpt2g[54 + isize]++;
			BF128 w;
			w.bf.u64[0] = sn->all_previous_cells;
			w.bf.u64[1] = sn->active_cells;
			switch (isize) {
			case 0: tbelow7[ntbelow[0]++] = w;
				tandbelow[0] &= sn->all_previous_cells; break;
			case 1: tbelow8[ntbelow[1]++] = w;
				tandbelow[1] &= sn->all_previous_cells; break;
			case 2: tbelow9[ntbelow[2]++] = w;
				tandbelow[2] &= sn->all_previous_cells; break;
			}
			goto next;
		}
	}
	// now a new ua is granted 
	if (sn->ncl == 9) {// 9 cells store it
		BF128 w;
		w.bf.u64[0] = sn->all_previous_cells;
		w.bf.u64[1] = sn->active_cells;
		t_7_9[nt_7_9++] = w;
		tand_7_9 &= sn->all_previous_cells;
		goto next;
	}
	sn->possible_cells = ua_ret7p & s->active_cells;
	if (sn->possible_cells) s++;  // switch to next spot
	goto next;
}
void G17B::Expand_10_12() {
	SPB03* sl = &spb_0_15[12], * s = sl, * sn;
	T54B12::TUVECT& tuv128 = t54b12.td128[0];
	uint64_t* twu = tuv128.t;
	//*s initial done by the caller
	s->possible_cells = twu[0];
	s->v[0] = tuv128.v0;// initial nothing done
	//cout << Char54out(twu[0]) << "first ua" << endl;
	//_______ start search clues 7_11
next:	// catch and apply cell in bitfields
	if (GetNextCellD(s))	if (--s >= sl)goto next;else  return ;
	sn = s + 1;
	//cout << Char54out(sn->all_previous_cells);
	//sn->cbs.Status(); cout << " ncl=" << sn->ncl<<endl;
	ua_ret7p = 0;
	// adjust active to the band constraint
	if (!(sn->active_cells &= sn->cbs.NextActive())) goto next;
	if (sn->ncl == 11) {// last step
		p_cpt2g[5]++;
		if (!GetLastAndUaD(sn)) {// can be a valid 10 or unknown ua(s)
			if (IsValid7pbf(sn->all_previous_cells)) {// got uas to use add it to D	
				t54b12.AddD(ua_ret7p);// re use smallest
				//cout << Char54out(ua_ret7p) << "added last nd128= "<< t54b12.nd128<< endl;
				ua_ret7p = anduab12;
			}
			else {// valid 11
				p_cpt2g[58]++;
				BF128 w;
				w.bf.u64[0] = sn->all_previous_cells;
				w.bf.u64[1] = sn->active_cells;
				tbelow11[ntbelow[4]++] = w;
				tandbelow[4] &= sn->all_previous_cells;
				goto next;
			}
		}
		{// last not empty direct or after check valid
			register uint64_t  P = ua_ret7p & sn->active_cells;
			//cout << Char54out(P) << "ua for last step " << endl;
			// expand directly P to 12, these are potential 12 666 666
			int cell;
			register uint64_t pc = sn->all_previous_cells;
			while (bitscanforward64(cell, P)) {
				register uint64_t bit  = (uint64_t)1 << cell;
				P ^= bit;
				bit |= pc;
				tfull[ntbelow[5]++] = bit;
				tandbelow[5] &= bit;
			}
			goto next;
		}
	}
	else {
		if (sn->v[0].isNotEmpty())GetNextUaD(sn);// first 128
		else GetNextUaAddD(sn);
	}

	if (!ua_ret7p) {// can only be 10 clues
		if (!IsValid7pbf(sn->all_previous_cells)) {// valid below
			p_cpt2g[56]++;
			BF128 w;
			w.bf.u64[0] = sn->all_previous_cells;
			w.bf.u64[1] = sn->active_cells;
			tbelow10[ntbelow[3]++] = w;
			tandbelow[3] &= sn->all_previous_cells;
			goto next;
		}
		else {
			//cout << Char54out(ua_ret7p) << "added 10 clues nd128= " << t54b12.nd128 << endl;
			t54b12.AddD(ua_ret7p);// re use smallest
		}
	}
	sn->possible_cells = ua_ret7p & sn->active_cells;
	//cout << Char54out(ua_ret7p) << "next ua" << endl;
	//cout << Char54out(sn->active_cells) << "sn->active_cells" << endl;
	//cout << Char54out(sn->possible_cells) << "next possible" << endl;
	if (sn->possible_cells) s++;  // switch to next spot
	goto next;
}



void G17B::GoExpand_7_12() {
	SPB03* sn = &spb_0_15[7];
	p_cpt2g[70]++;
	int locdiag = 0;
	if (op.f4) {
		if (p_cpt2g[70] == op.f4) {
			cout<< Char54out(sn->all_previous_cells) << "[4] good path expand 7_12 -v8- [70] "<< p_cpt2g[70] << endl;
			locdiag = 1;
		}
		else {
			cout << Char54out(sn->all_previous_cells) << "[70] " << p_cpt2g[70] << endl;
			if (p_cpt2g[70] > op.f4) { cout << "stop" << endl;	aigstop = 1; return; }
			if (!(op.upto4)) return;
		}
	}
	Expand_7_9();

	if (ntbelow[0]) {
		cout << "unexpected 7 clues pass2 stop [4]" << p_cpt2g[4] << endl;
		aigstop = 1;		return;
	}
	p_cpt2g[71]++;
	if (p_cpt2g[72] < nt_7_9)p_cpt2g[72] = nt_7_9;

	if (sgo.bfx[3] & 2) 		return;// debugging phase1

	myandall = tandbelow[1] & tandbelow[2]& tand_7_9;
	guah54_2.Build2(myandall, sn->active_cells);
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
		genb12.bands3[ib3].BuildGuam2(myandall);

	if (locdiag && op.ton) {
		cout << "test 7_12 [71]" << p_cpt2g[71];
		cout << " nt_7_9=" << nt_7_9 << "\t\t ";
		DumpPotential(0);
		//cout << Char54out(myandall) << " myandall " << endl;
		//cout << Char54out(sn->active_cells) << " active " << endl;
		//cout << "guah54_2 status" << endl;
		//guah54_2.Dumpall2(); guah54_2.Dumpall3();
	}

	if (ntbelow[1]) Go_8_12();
	if (ntbelow[2]) Go_9_12();
	for (uint32_t i1 = 0; i1 < nt_7_9; i1++)  
		GoExpand_10_12(t_7_9[i1]);
}
void G17B::GoExpand_10_12(BF128 ww){
	//p_cpt2g[79]++;
	//if (p_cpt2g[79] > 10) return;
	myb12_9 = ww.bf.u64[0];
	myac_9 = ww.bf.u64[1];
	int locdiag = 0;
	if (op.known) {
		if (!((~pk54) & myb12_9)) {
			cout << Char54out(myb12_9) << " expected 9 [70]"<< p_cpt2g[70] << endl;
			locdiag = 1;
		}
	}
	CBS cbs;
	//if (t54b12.Build_td128()) return;
	cbs.Init(myb12_9, 9);
	if (locdiag) {
		//t54b12.DebugC();
		//cout << Char54out(myb12_9) << " myb12_9 "  << endl;
		//cout << Char54out(myac_9) << " myac_9 " << endl;
		//t54b12.DebugD();
	}
	spb_0_15[12].Init9(ww, cbs);
	// clear processed 7,8,9 init the count for next expansion step
	memset(ntbelow, 0, sizeof ntbelow);//7 8 9 10 11 full
	memset(tandbelow, 255, sizeof tandbelow);//7 8 9 10 11 full
	if (t54b12.nd128 < 128)t54b12.nd128 = 128;// force adds outside
	Expand_10_12();
	p_cpt2g[73] += ntbelow[5];
	if (locdiag) 		DumpPotential(0);
	if (p_cpt2g[74] < ntbelow[5]) p_cpt2g[74] = ntbelow[5];

	if ((!ntbelow[5]) && (!ntbelow[4]) && (!ntbelow[3])) return;
	//if (p_cpt2g[79] > 1) return;

	myandall_9 = tandbelow[3] & tandbelow[4] & tandbelow[5];
	guah54_9.Build9(myandall_9, myac_9);
	if (locdiag && op.ton) {
		//cout << Char54out(myandall_9) << " myandall_9 " <<endl;
		//cout << Char54out(myac_9) << " myac_9 " << endl;
		for (uint32_t iv = 0; iv < ntbelow[3]; iv++) {
			uint64_t U = tbelow10[iv].bf.u64[0];
			cout << Char54out(U) << " 10 " << iv << endl;
			if (!((~pk54) & U))
				cout << Char54out(U) << " expected 10" << endl;
		}
		for (uint32_t iv = 0; iv < ntbelow[4]; iv++) {
			uint64_t U = tbelow11[iv].bf.u64[0];
			cout << Char54out(U) << " 11 " << iv << endl;
			if (!((~pk54) & U))
				cout << Char54out(U) << " expected 11" << endl;
		}

		for (uint32_t iv = 0; iv < ntbelow[5]; iv++) {
			uint64_t U = tfull[iv];
			cout << Char54out(U) << " " << iv << endl;
			if (!((~pk54) & U))
				cout << Char54out(U) << " expected 12 iv="<<iv << endl;
		}
		//cout << "guah54_2 status" << endl;
		//guah54_2.Dumpall2(); guah54_2.Dumpall3();
		//cout << "guah54_9 status" << endl;
		//guah54_9.Dumpall2(); guah54_9.Dumpall3();
	}
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
		genb12.bands3[ib3].BuildGuam9(myandall);

	if (ntbelow[4])Go_11_12();
	if (ntbelow[3])Go_10_12();


	if (ntbelow[5] == 1) {// get active g2 g3 from guah54_2 direct
		p_cpt2g[7]++;
		cb3.ncl = 12;
		myb12 = cb3.bf12 = tfull[0];
		cb3.cbs.Init(myb12, 12);
		cb3.g2t = guah54_9.GetG2(cb3.bf12);
		cb3.g3t = guah54_9.GetG3(cb3.bf12);
		clean_valid_done = 0;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
			genb12.bands3[ib3].Go(cb3);
			if (clean_valid_done == 2)break;
		}

	}
	else {
		for (uint32_t iv = 0; iv < ntbelow[5]; iv++) {
			cb3.ncl = 12;
			myb12 = cb3.bf12 = tfull[iv];
			cb3.cbs.Init(myb12, 12);
			cb3.g2t = guah54_9.GetG2(cb3.bf12);
			cb3.g3t = guah54_9.GetG3(cb3.bf12);
			p_cpt2g[7]++;
			if (locdiag && op.ton) {
				cout << Char54out(myb12) << " [7] " << p_cpt2g[7] << endl;
				cb3.Dump();
				if (p_cpt2g[7] == op.f7) {
					cout << "this is the expected path to go b3" << endl;
					//guah54_2.Dumpall2();
					//guah54_2.Dumpall3();
				}

			}
			clean_valid_done = 0;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
				genb12.bands3[ib3].Go(cb3);
				if (clean_valid_done == 2)break;
			}
		}
	}
}

*/


//____________processing band 3 to the end


//int tstack27[3] = { 07007007 ,070070070,0700700700 };


int G17B::IsValid_myb12() {
	if (zh2b[1].IsValid(myb12)) {
		if (op.dv12)cout << Char54out(myb12) << " isnot valid [7] " << p_cpt2g[7]
			<< " [3]" << p_cpt2g[3] << " [4]" << p_cpt2g[4]
			<<" n A B C " << t54b12.na128 <<" " << t54b12.nb128 <<" " << t54b12.nc128 << endl;
		anduab12 = ~0;
		register uint64_t ua54;
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t ua = zh2gxn.tua[i], cc = _popcnt64(ua);
			ua54 = (ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
			anduab12 &= ua54;
			if (cc < 12) {
				cout << Char54out(ua54) << " seen add b12 size " << cc << endl;
				aigstop = 1; cout << " stop uab12<12" << endl; return 1;
			}
			if (cc < 22) {
				if (op.dv12) 					
					cout << Char54out(ua54) << " to add b12 size " << cc<< endl;				
				t54b12.AddA(ua54);
				t54b12.AddB(ua54);
				t54b12.AddC(ua54);
			}
		}
		ua_ret7p = ua54;// return last (smaller)
		return 1;
	}
	return 0;
}
uint32_t G17B::IsValidB3(uint32_t bf,int debug) {
	p_cpt2g[31]++;
	if (nt3more) {
		for (uint32_t i = 0; i < nt3more; i++)if (bf == t3more[i]) return t3more[i];
	}
	if (debug) {
		cout <<Char27out(bf) << "IsValidB3 debugging mode " << endl;
	}
	int nret;
	if ((nret=zhou[0].CallCheckB3( bf,debug))) {
		anduab3 = BIT_SET_27;
		for (int iadd = 0; iadd < nret; iadd++) {
			BF128 w = zhgxn.tua[iadd];
			int cc = _popcnt32(w.bf.u32[2]);
			if (!cc) {
				/*
				cout << Char2Xout(w.bf.u64[0]) << " ";
				cout << Char27out(w.bf.u32[2]) << " iadd=" << iadd << endl;

				cout << Char54out(myb12) << " ";
				zhou[0].CallCheckB3(bf, 1);
				cout << " back count =" << zhgxn.nua << endl;
				for (uint32_t i = 0; i < zhgxn.nua; i++) {
					BF128 ww = zhgxn.tua[i];
					cout << Char2Xout(ww.bf.u64[0]) << " ";
					cout <<Char27out(ww.bf.u32[2])<<" i="<<i<<endl;
				}
				*/
				cout << Char27out(bf) << " bug no b3 IsValidB3 [7]" << p_cpt2g[7]
					<< "   [8]" << p_cpt2g[8] << "nr=" << nret << endl;
				aigstop = 1;	return 1;
			}
			register uint32_t ua = w.bf.u32[2];
			t3more[nt3more++] = ua;
			anduab3 &= ua;
			if (ua & ~t3infield)  t3outseen &= ua;

			register uint64_t U = w.bf.u64[0];
			uint64_t cc0 = _popcnt64(U);
			if (cc0 > 16)continue;

			U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
			if (cc <= 6 && cc>1) {
				if (cc > 3) {					
					if (!myband3->Check(w))myband3->Addguam(w, 1);
				}
				else {
					if (cc == 2) {
						int i81 = myband3->GetI81_2(w.bf.u32[2]);
						if (sgo.bfx[3] & 1) {
							if ( guah54.Check2(U, i81)) {
								cout << Char54out(U) << "\t";
								cout << Char27out(ua) << " bug redundant gua2 [7]" << p_cpt2g[7]
									<< "   [8]" << p_cpt2g[8] << endl;
								guah54.DumpOne2(i81);
								guah54_2.DumpOne2(i81);
								aigstop = 1;	return 1;
							}							
						}
						guah54.Add2(U, i81); guah54_2.Add2(U, i81);
						cb3.g2t.setBit(i81);
						p_cpt2g[20]++;
					}
					else {
						int i81 = myband3->GetI81_3(w.bf.u32[2]);
						if (sgo.bfx[3] & 1) {
							if (guah54.Check3(U, i81)) {
								cout << Char54out(U) << "\t";
								cout << Char27out(ua) << " bug redundant gua3 [7]" << p_cpt2g[7]
									<< "   [8]" << p_cpt2g[8] << endl;
								guah54.DumpOne3(i81);
								guah54_2.DumpOne3(i81);
								aigstop = 1;	return 1;
							}
						}
						guah54.Add3(U, i81); guah54_2.Add3(U, i81);
						cb3.g3t.setBit(i81);
						p_cpt2g[21]++;
					}
				}

			}

		}
		return 1;
	}
	else return 0;
}

void STD_B3::Go(CALLBAND3& cb3) {
	p_cpt2g[8]++;
	if (g17b.knownt == 11)
		cout << "entry go expected ok clean=" << g17b.clean_valid_done 
		<<" aigstop= "<< g17b.aigstop << endl;
	//if (op.f7 < p_cpt2g[8]) { g17b.aigstop = 1; return; }
	if (g17b.clean_valid_done == 2) return;
	if (g17b.aigstop) return;
	CRITB3 scritb3;
	BF128 g2 = cb3.g2t & g.gsocket2,	g3 = cb3.g3t & g.gsocket3;
	register uint32_t  vmini, Mg2 = 0, Mg3 = 0;
	int x;
	while ((x = g2.getFirst128()) >= 0) {
		g2.clearBit(x);		Mg2 |= g.ua2bit27[x];	}
	while ((x = g3.getFirst128()) >= 0) {
		g3.clearBit(x);		Mg3 |= 1<<g.ua3_imini[x];	}
	scritb3.Init(cb3.ncl, cb3.cbs);
	for (int imini = 0; imini < 9; imini++, Mg2 >>= 3, Mg3 >>= 1) {
		if (!(vmini = Mg2 & 7))	if (Mg3 & 1)vmini = 8;
		if (vmini)scritb3.MinVmini(imini, vmini);
	}
	if (g17b.knownt == 11) {
		cout << "scritb3.mincount=" << scritb3.mincount << "scritb3.nb3=" << scritb3.nb3 << endl;
		scritb3.Status("aaa");
	}

	if ((int)scritb3.mincount >   scritb3.nb3) return;
	scritb3.nmiss = scritb3.nb3 - scritb3.mincount;
	
	//memcpy(&genb12.grid0[54], band0, sizeof band0);// used in brute force
	memcpy(&g17b.grid0[54], band0, sizeof band0);// used in brute force
	g17b.myband3 = this;
	g17b.pbufuas3 = g17b.bufuas3;
	CRITHANDLER crh; crh.Init(g17b.pbufuas3);
	crh.mycritb3 = scritb3;
	int locdiag = 0;
	if (p_cpt2g[7] <= op.f7) {
		//cout << "go b3 locdiag=1 [8]" << p_cpt2g[8] << endl;
		locdiag = 1;
	}
	//if (locdiag)crh.mycritb3.Status(" before build t3  ");
	g17b.t3infield = crh.mycritb3.critbf;
	g17b.nt3more = 0;
	g17b.ncluesb3 = crh.mycritb3.nb3;
	p_cpt2g[9]++;

	//if(locdiag)cout << "go b3 locdiag=1 [8]" << p_cpt2g[8]<<" nmiss = "<< crh.mycritb3.nmiss<< endl;
	{ // add g2 g2
		register uint32_t V = scritb3.pairs27, i27;
		while (bitscanforward64(i27, V)) {
			V ^= 1 << i27;	crh.Add(g.pat2[i_27_to_81[i27]]);		}
		V = scritb3.minix[0];// triplets
		if (V) {// add guas3 if any (usually 0)
			for (int i = 0; i < 9; i++)if (V & (1 << i))
				crh.Add(7 << (3 * i));		}
	}
	register uint64_t F = g17b.myb12;
	register uint32_t  		Ac = scritb3.critbf, U;// if critical all must be in
	g17b.critical_done = 0;
	if (g17b.knownt == 11) {
		scritb3.Status("bbb");
	}

	g17b.Init_t3o(scritb3);
	if (!scritb3.nmiss) {//critical status uas reduced no out
		for (uint32_t i = 0; i < ntguam2; i++) {
			GUAM&  w = tguam2[i];
			if (F & w.bf12) continue;
			if (!(U= w.bf3 & Ac)) return; // dead branch
			crh.AddIf(U);	
		}
		for (uint32_t i = 0; i < nua; i++) {// now band 3 uas
			if (!(U = tua[i] & Ac)) return; // dead branch
			crh.AddIf(U);
		}
		if (scritb3.minix[2]) {// better try direct
			g17b.pbufuas3 = crh.Lock();
			p_cpt2g[10]++;
			g17b.GoB3Critical(crh);
		}
		else g17b.CheckCritical(cb3, crh);
		return;
	}

	// not miss 0 no "active" split in/out
	crh.mycritb3.active = BIT_SET_27 & ~crh.mycritb3.critbf;
	for (uint32_t i = 0; i < ntguam2; i++) {
		GUAM& w = tguam2[i];
		if (F & w.bf12) continue;
		U = w.bf3;
		if (!(U & Ac))g17b.AddT3o(U);  
		else crh.AddIf(U);
	}
	for (uint32_t i = 0; i < nua; i++) {// now band 3 uas
		U = tua[i]&BIT_SET_27;
		if (!(U & Ac)) g17b.AddT3o(U);
		else crh.AddIf(U);
	}
	g17b.SortT3o(~0);
	//check critical if used later
	g17b.CheckCritical(cb3, crh);
}
inline void G17B::CheckCritical(CALLBAND3& cb3, CRITHANDLER& crh) {
	if (knownt == 11) {
		crh.Dump(); Dumpt3o();
		cout << Char27out(p17diag.bf.u32[2]) << endl;
		cout << myband3->band << endl;
	}
	critical_done = 0;
	CRITB3& mcr = crh.mycritb3;
	if (p_cpt2g[7] == op.f7)
		cout << "\t\t\t\t\t entry check GoCommonB3 mcr.nmiss=" << mcr.nmiss 
		<< " [8] " << p_cpt2g[8] 
		<<" start clean_valid_done="<< clean_valid_done << endl;
	if (mcr.nmiss > 10) {
		cout << "CheckCritical 1 nmiss invalide" << mcr.nmiss << " [3] " << p_cpt2g[3]
			<< " [4] " << p_cpt2g[4]   << " [7] " << p_cpt2g[7]
			<< " [8] " << p_cpt2g[8] << endl;
		aigstop = 1; return;
	}
	if (p_cpt2g[8] == sgo.vx[9]) {
		mcr.Status("entry check critical");
		crh.Dump();
		Dumpt3o();
	}
	if (mcr.nmiss > 3 || mcr.mincount < 5) {	GoCommonB3(crh); return;	}
	if (VerifyValid())return;
	critical_done = 1;
	if (!zhou[0].CallCheckB3(mcr.critbf)) { GoCommonB3(crh); return; }
	if (!mcr.nmiss)return; // can not add out here
	critical_done = 1;// force redo check later
	// critbf not valid, check fresh uas
	for (uint32_t iadd = 0; iadd < zhgxn.nua; iadd++) {
		BF128 w = zhgxn.tua[iadd];
		int cc = _popcnt32(w.bf.u32[2]);
		register uint32_t ua = w.bf.u32[2];
		register uint64_t U = w.bf.u64[0];
		uint64_t cc0 = _popcnt64(U);
		U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
		//cout << Char54out(U) << " \t";
		//cout << Char27out(ua) << " added " << cc0 << " " << cc << endl;


		if (cc == 2) {
			int i81 = myband3->GetI81_2(ua),
				imini = myband3->g.ua2_imini[i81],
				mask = 7 << (3 * imini), bitmini = 1 << imini;
			cb3.g2t.setBit(i81);
			if (mcr.critbf & mask) {// must be a triplet coming before
				if (!(mcr.minix[0] & bitmini)) continue;
				// clear previous
				mcr.minix[0] ^= bitmini;
				mcr.critbf &= ~mask;
				mcr.mincount--; mcr.nmiss++;
			}
			if (!mcr.nmiss) return;
			mcr.mincount++; mcr.nmiss--;
			mcr.minix[1] |= bitmini;
			mcr.pairs27 |= (mask ^ ua);
			mcr.critbf |= ua;
			crh.Add(ua);
			if (cc0 > 16)continue;
			guah54.Add2(U, i81); guah54_2.Add2(U, i81);
			p_cpt2g[20]++;
			continue;
		}
		if (cc == 3) {
			int i81 = myband3->GetI81_3(ua),
				imini = myband3->g.ua3_imini[i81],
				bitmini = 1 << imini;
			cb3.g3t.setBit(i81);
			if (!(mcr.critbf & ua)) {//not  gua2 added
				if (!mcr.nmiss) return;
				mcr.mincount++; mcr.nmiss--;
				mcr.minix[0] |= bitmini;
				mcr.critbf |= ua;
				crh.Add(ua);
			}
			if (cc0 > 16)continue;
			guah54.Add3(U, i81); guah54_2.Add3(U, i81);
			p_cpt2g[21]++;
			continue;
		}
		AddT3o(ua);
		if (cc0 < 17)if (!myband3->Check(w))myband3->Addguam(w, 1);
	}
	if (mcr.nmiss > 10) {
		cout << "CheckCritical ff nmiss invalide" << mcr.nmiss << " [3] " << p_cpt2g[3]
			<< " [4] " << p_cpt2g[4]   << " [7] " << p_cpt2g[7]
			<< " [8] " << p_cpt2g[8] << endl;
		aigstop = 1; return;
	}
	//need to revise out if critbf added
	int n = nt3o;
	mcr.active = BIT_SET_27 & ~mcr.critbf;
	Init_t3o(mcr);
	for (int i = 0; i < n; i++) {
		register uint32_t U = t3o[i];
		if (U & t3infield)crh.Add(U);
		else AddT3o(U);
	}
	//mcr.Status(" after check criticalvalid");
	GoCommonB3(crh);
}
void G17B::GoCommonB3(CRITHANDLER& crh) {
	g17b.pbufuas3 = crh.Lock();
	CRITB3& mcr = crh.mycritb3;
	nt3more = 0;
	if (knownt == 11)
		mcr.Status("entry go commonb3=");

	if (!mcr.nmiss) {
		if (!nt3o) { p_cpt2g[10]++; GoB3Critical(crh); }
		return;
	}
	if (mcr.nmiss == 1) {
		if (nt3o) {	if (!t3ando) return;}
		else 	t3ando = 0;
		p_cpt2g[11]++;
		GoB3Miss1(crh);
		return;
	}
	if (mcr.nmiss == 2) {
		GoB3Miss2(crh);
		return;
	}
	if (nt3o)nt3o = SortT3o();
	if (mcr.nmiss == 3) {	GoB3Miss3(crh); return;	}
	else GoB3MissMore(crh);

}

int CRITHANDLER::CleanOne(CRITHANDLER& hout, int debug) {// assign all singles in uas

	//Build  reduced table in hout
	{
		register uint32_t AC = mycritb3.critbf, F = mycritb3.assigned;
		hout.nuas = 0;
		for (int i = 0; i < (int)nuas; i++) {
			register uint32_t U = tuas[i];
			if (!(U & F)) {
				if (!(U &= AC)) return 1;
				hout.tuas[hout.nuas++] = U;

			}
		}
	}
	if (debug) {
		hout.Dump();
		//return 0;
	}
	return hout.CleanGo(debug);
}
int CRITHANDLER::CleanGo(int debug ) {// assign all singles in uas
	while (1) {
		int is1 = 0;
		for (int i = 0; i < (int)nuas; i++) {
			register uint32_t U = tuas[i];
			if (_popcnt32(U) == 1)  is1 |= U;
		}

		if (is1) {
			if(debug)cout << Char27out(is1) << "is1 p_cpt2g[8]= " << p_cpt2g[8] << endl;
			if (mycritb3.AddAssign(is1, 0)){
				if (debug)cout << " add assign Cleango conflict or dead" << endl;
				return 1;// conflict
			}
			if(debug)mycritb3.Status(" after add assign ");
			register uint32_t  F = mycritb3.assigned;
			register uint32_t AC = mycritb3.active;
			if (debug) {
				cout << Char27out(F) << "F" << endl;
				cout << Char27out(AC) << "AC" << endl;
			}

			int n = nuas; nuas = 0;
			for (int i = 0; i < n; i++) {
				register uint32_t U = tuas[i];
				if (!(U & F)) {
					U &= AC;
					if (!U) return 1; // dead branch
					tuas[nuas++] = U;
				}
			}
			if (debug) cout << "loop nuas=" << nuas << endl;
		}
		else return 0;
	}
}
void CRITHANDLER::SortShrink(int shrink) {// no single here 
	uint32_t t2[100], nt2 = 0, t3[100], nt3 = 0, t4[400], nt4 = 0, 
		tm[400], ntm = 0;
	for (uint32_t i = 0; i < nuas; i++) {
		register uint32_t U = tuas[i], cc = _popcnt32(U);
		switch (cc) {
		case 2:t2[nt2++] = U; break;
		case 3:t3[nt3++] = U; break;
		case 4:t4[nt4++] = U; break;
		default:tm[ntm++] = U; break;
		}
	}
	nuas = 0;
	// add now more clues 
	for (uint32_t i = 0; i < nt2; i++) tuas[nuas++] = t2[i];
	for (uint32_t i = 0; i < nt3; i++) AddIf(t3[i]);
	for (uint32_t i = 0; i < nt4; i++) AddIf(t4[i]);
	for (uint32_t i = 0; i < ntm; i++) AddIf(tm[i]);
}
void G17B::BuildSortT3b() {
	uint32_t tx[7][100], ntx[7];
	memset(ntx, 0, sizeof ntx);
	for (uint32_t i = 0; i < nt3exp; i++) {
		register uint32_t U = t3exp[i] , cc = _popcnt32(U);
		if (cc > 6)cc = 6;
		tx[cc][ntx[cc]++] = U;
	}
	nt3b = 0;
	for (uint32_t i = 0; i < 7; i++) if (ntx[i]) {
		register uint32_t* tu = tx[i];
		for (uint32_t j = 0; j < ntx[i]; j++) {
			register uint32_t U = tu[j], nU = ~U;
			for (uint32_t k = 0; k < nt3b; k++) {
				if (!(t3b[k] & nU)) { U = 0; break; }
			}
			if (U)t3b[nt3b++] = U;
		}
	}
}
int G17B::SortT3o(uint32_t active) {// order t3o on size
	if (!active) {
		uint32_t t4[200], t6[100], tm[100], nt4 = 0, nt6 = 0, ntm = 0;
		for (uint32_t i = 0; i < nt3o; i++) {
			register uint32_t U = t3o[i], cc = _popcnt32(U);
			if (cc < 5)t4[nt4++] = U;
			else if (cc < 7)t6[nt6++] = U;
			else tm[ntm++] = U;
		}
		int n = 0;
		for (uint32_t i = 0; i < nt4; i++)t3o[n++] = t4[i];
		for (uint32_t i = 0; i < nt6; i++)t3o[n++] = t6[i];
		for (uint32_t i = 0; i < ntm; i++)t3o[n++] = tm[i];
		return n;
	}
	// apply active and sort 
	uint32_t tx[7][50], ntx[7];
	memset(ntx, 0, sizeof ntx);
	register uint32_t A = active;
	for (uint32_t i = 0; i < nt3o; i++) {
		register uint32_t U = t3o[i]&A, cc = _popcnt32(U);
		if (cc >6)cc=6;
		tx[cc][ntx[cc]++] = U;
	}
	nt3o = 0;
	if (ntx[0]) return 1;
	for (uint32_t i = 1; i < 7; i++) if (ntx[i]) {
		register uint32_t *tu = tx[i];
		for (uint32_t j = 0; j < ntx[i]; j++) {
			register uint32_t U = tu[j],nU=~U;
			for (uint32_t k = 0; k < nt3o; k++) {
				if (!(t3o[k] & nU)) { U = 0; break; }
			}
			if (U)t3o[nt3o++] = U;
		}
	}
	return 0;
}

void G17B::GoB3Miss1(CRITHANDLER& crh) {
	//int backmiss1 = ~0;
	int locdiag = 0;//  2 * (p_cpt2g[8] == sgo.vx[9]);
	//if (p_cpt2g[7] == 4449)locdiag = 2;
	if (knownt == 11)locdiag = 2;
	if (op.ton > 1)if (p_cpt2g[7] == op.f7) {
		locdiag = op.ton;
	}
	if (locdiag) {
		cout << Char27out(crh.mycritb3.assigned) << "gomiss1" << endl;
		cout << Char27out(t3ando) << "t3ando" << endl;
		crh.mycritb3.Status("ee");
		crh.Dump();
	}
	if (!t3ando) {
		if (VerifyValid()) return;
		if (IsValidB3(crh.mycritb3.critbf| crh.mycritb3.assigned))	t3ando = anduab3;
		if (!t3ando) {
			GoB3SubCritical(crh);	
			t3ando=BIT_SET_27 & ~(t3infield | crh.mycritb3.assigned);
		}
	}
	if (locdiag) 
		cout << Char27out(t3ando) << "t3ando after subcritical" << endl;



	// now t3ando is the outfiekd to use and critical
	if (VerifyValid()) return;
	crh.mycritb3.AssignCritical();
	if (locdiag)crh.mycritb3.Status("ff");
	uint32_t res1;
	uint32_t active = BIT_SET_27 & ~t3infield;
	//nt3o = 0;
	while (bitscanforward(res1, t3ando)) {
		int bit1 = 1 << res1;
		t3ando ^= bit1; active ^= bit1;
		CRITHANDLER crhn = crh;
		crhn.mycritb3.assigned |= bit1 ;
		CRITHANDLER crhn2 = crhn;
		crhn2.Init(pbufuas3);
		if (locdiag)crhn.mycritb3.Status(" avt cleanone");
		if (knownt >= 11) {
			if (IspathB3(crhn.mycritb3.assigned)) {
				cout << "on the path" << endl;
				knownt = 20;
			}
		}

		crhn2.mycritb3.nmiss = 0;// needed in cleanone
		if (crhn.CleanOne(crhn2))continue;;// conflict or dead	
		EndGoB3Critical(crhn2);
		if (knownt >= 20) return;
	}
}
void  G17B::GoB3Miss2(CRITHANDLER& crh) {
	p_cpt2g[12]++;
	nt3more = 0;
	int locdiag = 0;// 1;// 2 * (p_cpt2g[8] == 851);
	//if (p_cpt2g[7] == 749)locdiag = 2;
	//if (p_cpt2g[12] == 24)locdiag = 2;
	if (!crh.mycritb3.minix[2]) {// nothing better than direct
		GoB3MissMore2A(crh); return;}
	if (locdiag)cout << Char9out( crh.mycritb3.minix[2])<<" bf2 miss2 entry [12] "<< p_cpt2g[12] <<endl;
	if (locdiag > 1) { crh.mycritb3.Status(" miss2 entry ");	crh.Dump(); Dumpt3o(); }

	if (!nt3o) {
		if (VerifyValid()) return;
		if (IsValidB3(crh.mycritb3.critbf)) {
			t3ando = anduab3;
			if (locdiag)cout << Char27out(t3ando) << " not valid gob3miss2  " << endl;
			// store uas in t3o
			for (uint32_t i = 0; i < nt3more; i++)		t3o[i] = t3more[i];
			nt3o = nt3more;
			nt3more = 0;
		}
		else t3ando = BIT_SET_27 & (~t3infield) & (~crh.mycritb3.assigned);
	}
	if (locdiag) 	cout << Char27out(t3ando) << " t3ando  " << endl;			

	uint32_t res;
	if (t3ando) {// do in priority hitting all out
		uint32_t v = t3ando,vr=v;
		while (bitscanforward(res, v)) {
			int bit = 1 << res; v ^= bit;
			CRITHANDLER crhn = crh;
			crhn.mycritb3.Addoutbfone(bit);
			uint32_t rmore = nt3more;
			t3ando = 0;
			if (locdiag) 	cout << Char27out(crhn.mycritb3.assigned) << " call miss1  " << endl;

			GoB3Miss1(crhn);
			v &= t3outseen;
			if (nt3more>rmore) {
				for (uint32_t i = rmore; i < nt3more; i++) {
					register uint32_t U = t3more[i];
					if (!(U & t3infield)) {
						if (locdiag)cout << Char27out(U) << " more" << endl;
						v &= U;
						t3o[nt3o++] = U;
						
					}
				}
			}
			nt3more = rmore;
		}
		if((!nt3o) ||t3o[0]==vr)return;
		crh.mycritb3.active &= ~vr; 
	}
	if (SortT3o(crh.mycritb3.active)) {
		if (locdiag)cout << "seen empty ua" << endl;return;	}
	if (locdiag) { cout << "sorted" << endl; Dumpt3o(); }

	if (nt3o<2) { GoB3MissMore2A(crh); return; }
	crh.mycritb3.AssignCritical();
	uint32_t res1, res2;
	uint32_t active = BIT_SET_27 & ~t3infield, u1 = t3o[0];
	if (locdiag)cout << Char27out(active) << " active " << endl;
	while (bitscanforward(res1, u1)) {
		int bit1 = 1 << res1;
		u1 ^= bit1; active ^= bit1;
		uint32_t wand = ~0;// 1 other uas to have 3 clues
		for (uint32_t j =1; j < nt3o; j++) {
			register uint32_t U = t3o[j];
			if (!(U & bit1)) wand &= U & active;
		}
		if (locdiag)cout << Char27out(wand) << " loop2 " << endl;

		while (bitscanforward(res2, wand)) {
			int bit2 = 1 << res2; wand ^= bit2;
			CRITHANDLER crhn = crh;
			crhn.mycritb3.assigned |= (bit1 | bit2);
			CRITHANDLER crhn2 = crhn;
			crhn2.Init(pbufuas3);
			crhn2.mycritb3.nmiss = 0;// needed in cleanone
			if (crhn.CleanOne(crhn2))continue;;// conflict or dead	
			if (locdiag)crhn2.mycritb3.Status("at call endgob3");
			EndGoB3Critical(crhn2);
		}
	}

}
void  G17B::GoB3Miss3(CRITHANDLER& crh) {
	p_cpt2g[13]++;
	int locdiag = 0;// 1;// 2 * (p_cpt2g[8] == 851);
	if (p_cpt2g[8] == 2343933)locdiag = 2;
	if (locdiag) {
		crh.mycritb3.Status("entry miss3");
		Dumpt3o();
	}
	if (crh.mycritb3.minix[2]<2) {// keep that for later
		GoB3MissMore2A(crh); return;	}
	if (t3ando) { GoB3MissMore2A(crh); return; }// too complex see later 
	// keep for later if less than 3 clues possible 
	uint32_t res1, res2, res3;
	{
		uint32_t active = BIT_SET_27 & ~t3infield, u1 = t3o[0],and2;
		while (bitscanforward(res1, u1)) {
			int bit1 = 1 << res1;
			u1 ^= bit1; active ^= bit1;
			and2 = active;
			uint32_t wand = ~0;// 1 check if 2 clues is enough
			for (uint32_t j = 1; j < nt3o; j++) {
				register uint32_t U = t3o[j];
				if (!(U & bit1))and2&= U; 
			}
			if (and2)break;
		}
		if(and2) { GoB3MissMore2A(crh); return; }
	}

	if(locdiag)cout<<Char9out(crh.mycritb3.minix[2])
		<< "entry go miss 3 to try nt3o="<<nt3o << endl;
	if (VerifyValid()) return; 
	// we have outfield potential, what about infield 
	crh.mycritb3.AssignCritical();
	uint32_t active = BIT_SET_27 & ~t3infield, u1 = t3o[0];
	if (locdiag) {
		cout << Char27out(u1)<<"first t3o"<<endl;
		cout << Char27out(active) << "active" << endl;
	}
	while (bitscanforward(res1, u1)) {
		int bit1 = 1 << res1,u2=0; 
		u1 ^= bit1; active ^= bit1;
		uint32_t wand = ~0;// 1 check if 2 clues is enough
		for (uint32_t j = 1; j < nt3o; j++) {
			register uint32_t U = t3o[j];
			if (!(U & bit1)) {u2 = U & active; break;}
		}
		while (bitscanforward(res2, u2)) {
			int bit2 = 1 << res2; u2 ^= bit2;
			uint32_t wand = ~0;// 1 other uas to have 3 clues
			for (uint32_t j = 2; j < nt3o; j++) {
				register uint32_t U = t3o[j];
				if (!(U & bit1)) wand &= U & active;
			}
			while (bitscanforward(res3, wand)) {
				int bit3 = 1 << res3; wand ^= bit3;
				CRITHANDLER crhn = crh;
				crhn.mycritb3.assigned |= (bit1 | bit2 | bit3);
				CRITHANDLER crhn2 = crhn;
				if (locdiag)crhn2.mycritb3.Status("call critical");
				crhn2.Init(pbufuas3);
				crhn2.mycritb3.nmiss = 0;// needed in cleanone
				if (crhn.CleanOne(crhn2))continue;;// conflict or dead	
				EndGoB3Critical(crhn2);
			}
		}
	}
 
}
void  G17B::GoB3MissMore(CRITHANDLER& crh) {
	p_cpt2g[14]++;
	if (p_cpt2g[14] < 10) {
		cout << "GoB3MissMore  " << endl;
		//crh.mycritb3.Status("missmore ");
		//crh.Dump();
		//Dumpt3o();
	}
	for (uint32_t i = 0; i < nt3o; i++)crh.Add(t3o[i]);
	{ GoB3MissMore2A(crh); return; }
	return;
	// go direct if no critical assign
	if (!crh.mycritb3.minix[2]) { GoB3MissMore2A(crh); return; }

}
void  G17B::GoB3MissMore2A(CRITHANDLER& crh) {// add t3o to t3 and go
	//for (uint32_t i = 0; i < nt3o; i++)crh.Add(t3o[i]);
	p_cpt2g[22 + (crh.nuas >> 5)]++;
	ntoassb3 = crh.mycritb3.GetToAss();
	if (knownt > 10) {
		cout << " entry more2a ntoassb3=" << ntoassb3 << " knownt=" << knownt << endl;
		crh.mycritb3.Status(" entry 2A");
		crh.Dump();
	}
	if (ntoassb3 < 3)p_cpt2g[26]++;
	else if(ntoassb3 > 4)p_cpt2g[29]++;
	else p_cpt2g[24+ntoassb3]++;
	t3exp = crh.tuas; nt3exp = crh.nuas;
	if (ntoassb3 >= 4) {
		crh.SortShrink(1);
		if (knownt >= 20) {
			crh.Dump();
		}
		GoB3Expand_1_3(crh.mycritb3);
		return;
	}
	// go direct with t3b
	BuildSortT3b(); 
	if (knownt >= 20) {
		Dumpt3b();
	}
	SP3 sp3;
	sp3.active = BIT_SET_27 ^ crh.mycritb3.assigned;
	sp3.all_previous = crh.mycritb3.assigned;
	GoB3Expand_4_x(sp3);
}


void G17B::GoB3Critical(CRITHANDLER& crh) {
	p_cpt2g[17]++;
	t3outseen = ~0;
	CRITHANDLER crhn = crh;
	crhn.mycritb3.AssignCritical();
	CRITHANDLER crhn2 = crhn;
	crhn2.Init(pbufuas3);
	if (crhn.CleanOne(crhn2))return;// conflict or dead
	EndGoB3Critical(crhn2);
}
void G17B::EndGoB3Critical(CRITHANDLER& crh) {
	t3exp = crh.tuas; nt3exp = crh.nuas;
	int ntoass = crh.mycritb3.GetToAss();
	if (ntoass < 0)return; //safety
	if ((int)crh.nuas < ntoass) {
		cout << " EndGoB3Critical bug manque uas [3] " << p_cpt2g[3]
			<< " [4] " << p_cpt2g[4] << " [7] " << p_cpt2g[7] << endl;
		crh.Dump();		aigstop = 1; return;	}
	p_cpt2g[18]++;

	//<<<<<<<<<<<<<<<<<<<<<<<<<<<
	// here test sort of t3
	//<<<<<<<<<<<<<<<<<<<<<<<<<
	if (!ntoass) {// job is done check if valid
		if (nt3exp)return;// should never be
		if(VerifyValid()) return; 
		if (IsValidB3(crh.mycritb3.assigned)) return;
		Out17(crh.mycritb3.assigned);
		return;
	}
	//if (p_cpt2g[18] > 500) { aigstop = 1; return; }

	//cout << " EndGoB3Critical [3] " << p_cpt2g[3]
		//<< " [4] " << p_cpt2g[4] << " [7] " << p_cpt2g[7]
		//<< " [17] " << p_cpt2g[17] << " [18] " << p_cpt2g[18]
		//<< "  ntoass" << ntoass << endl;		

	if (ntoass > 1) {	GoB3MissMore2A( crh); return;	}
	// ______________________only one to ass must hit all uas 
	register int A = t3exp[0];
	for (int i = 1; i < (int)nt3exp; i++)
		if (!(A &= t3exp[i])) return;
	if (VerifyValid()) return;
	register int cell;
	while (bitscanforward(cell, A)) {
		int ass = crh.mycritb3.assigned,bit= 1 << cell;
		A ^=bit; //clear bit
		ass |= bit;
		if (IsValidB3(ass)) A &= anduab3;
		else 	Out17(ass);
	}
}
void G17B::GoB3SubCritical(CRITHANDLER& crh) {
	p_cpt2g[19]++;
	if(p_cpt2g[19]<10)cout << "entry subcritical [3]"<<p_cpt2g[3]
		<<" [4]" << p_cpt2g[4] << " [7]" << p_cpt2g[7] << endl;
	for (int i = 0, mask = 7, bit = 1; i < 9; i++, mask <<= 3, bit <<= 1) {// 9 mini rows
		CRITHANDLER crhn = crh;
		if (!(crh.mycritb3.critbf & mask)) continue;
		//	uint32_t minix[4],// triplet bf1 bf2 bf3  
		if (crh.mycritb3.minix[3] & bit) {
			crhn.mycritb3.AddAssign(mask);
			GoB3Critical(crhn);
			continue;
		}
		if (crh.mycritb3.minix[1] & bit) {
			crhn.mycritb3.AddAssign(mask& crh.mycritb3.critbf);
			GoB3Critical(crhn);
			continue;
		}
		// now 2 out of 3 (bf2 or triplet)
		register uint32_t A = mask;
		int cell;
		while (bitscanforward(cell, A)) {
			uint32_t b = 1 << cell;
			A ^= b;
			b ^= mask;// 2 cells to assign
			crhn = crh;
			crhn.mycritb3.AddAssign(b);
			GoB3Critical(crhn);
		}
	}
}


void G17B::GoB3Expand_1_3(CRITB3& ecritb3) {
	int locdiag = 0;// = 2 * (p_cpt2g[8] == 933);
	//if (p_cpt2g[18]==43) locdiag = 1;
	if (op.ton)if (p_cpt2g[7] == op.f7) locdiag = op.ton;
	if (knownt >= 20)locdiag = 2;
	if (locdiag) {
		cout << "expand entry nto ass =" << ntoassb3
			<< " clean_valid_done=" << clean_valid_done
			<< " nt3expand=" << nt3exp << endl;
		Dumpt3exp();
	}
	SP3 spb[5];

	spb[0].active =BIT_SET_27^ ecritb3.assigned;
	spb[0].all_previous = ecritb3.assigned;
	spb[0].indtw3 = 0;// initial first ua
	spb[0].possible_cells = t3exp[0];
	if (locdiag) {
		cout << Char27out(spb[0].all_previous) << " initial assigned" << endl;
		//cout << Char27out(spb[0].active) << " initial active" << endl;
	}
	if (locdiag) cout << Char27out(spb[0].possible_cells) << "start" << endl;

next1:
	{// catch and apply cell in bitfields
		register int cell;
		register uint32_t p = spb[0].possible_cells;
		if (!p) {
			if (locdiag) cout << "back 1_3" << endl;
			return;
		}
		bitscanforward(cell, p);
		register uint32_t bit = 1 << cell;
		spb[0].possible_cells ^= bit;
		spb[0].active ^= bit;
		spb[1] = spb[0];
		spb[1].all_previous |= bit;
		register uint32_t F = spb[1].all_previous;
		for (uint32_t i = spb[0].indtw3 + 1; i < nt3exp; i++) {
			register uint32_t U = t3exp[i];
			if (!(U & F)) {
				U &= spb[1].active;
				//if (!U)goto next1;//dead branch
				spb[1].possible_cells = U;
				spb[1].indtw3 = i;
				if (locdiag) cout << Char27out(spb[1].possible_cells) << "first" << endl;
				goto next2;
			}
		}
		if (VerifyValid()) return;
		if (Valid3_1_3(spb[1].all_previous))goto next1;
		spb[1].possible_cells = 0; // see later 
	}
next2:
	{// catch and apply cell in bitfields
		register int cell;
		register uint32_t p = spb[1].possible_cells;
		if (!p) goto next1;
		bitscanforward(cell, p);
		register uint32_t bit = 1 << cell;
		spb[1].possible_cells ^= bit;
		spb[1].active ^= bit;
		spb[2] = spb[1];
		spb[2].all_previous |= bit;
		register uint32_t F = spb[2].all_previous;
		for (uint32_t i = spb[1].indtw3 + 1; i < nt3exp; i++) {
			register uint32_t U = t3exp[i];
			if (!(U & F)) {
				U &= spb[2].active;
				if (!U)goto next2;// dead branch 
				spb[2].possible_cells = U;
				spb[2].indtw3 = i;
				if (locdiag) cout << Char27out(spb[2].possible_cells) << "second" << endl;
				goto next3;
			}
		}
		if (VerifyValid()) return;
		if (Valid3_1_3(spb[2].all_previous))goto next2;
		spb[2].possible_cells = 0; // see later 
	}
next3:
	{// catch and apply cell in bitfields
		register int cell;
		register uint32_t p = spb[2].possible_cells;
		if (!p) goto next2;
		bitscanforward(cell, p);
		register uint32_t bit = 1 << cell;
		spb[2].possible_cells ^= bit;
		spb[2].active ^= bit;
		spb[3] = spb[2];
		spb[3].all_previous |= bit;
		if (locdiag) cout << Char27out(spb[3].all_previous) << " 3 cells " << endl;

	}
	if (ntoassb3 != 4) goto endnext3;
	// this is the last step for 4
	{
		register uint32_t F = spb[3].all_previous,
			wa = spb[2].active, n = 0;
		for (uint32_t i = spb[2].indtw3 + 1; i < nt3exp; i++) {
			register uint32_t U = t3exp[i];
			if (!(U & F)) {
				n++;
				if (!(wa &= U)) goto next3; // dead
			}
		}
		if (locdiag) cout << Char27out(wa) << " wa " << endl;
		if (VerifyValid()) return;// not sure it is done
		if (!n) {// could be a valid
			if (Valid3_1_3(spb[3].all_previous))goto next3;
			wa = anduab3;
		}
		if (locdiag) cout <<  " go for  wa " << endl;

		uint32_t cell;
		while (bitscanforward(cell, wa)) {
			register uint32_t bit = 1 << cell,
				bf= spb[3].all_previous | bit;
			wa ^= bit; //clear bit
			if (Valid3_1_3(bf))
				wa &= anduab3;
			else 	Out17(bf);//can be valid 17 in 18 mode

		}
	}
	goto next3;
endnext3:
	{// build a reduced table for the next clues
		nt3b = 0;
		uint32_t tx[8][30], ntx[7];// gives 60 for 6 and more 
		memset(ntx, 0, sizeof ntx);
		register uint32_t A = spb[2].active, F = spb[3].all_previous;
		for (uint32_t i = spb[2].indtw3 + 1; i < nt3exp; i++) {
			register uint32_t U = t3exp[i];
			if (U & F) continue;
			if (!(U &= A)) goto next3;
			register uint32_t cc = _popcnt32(U);
			if (cc > 6)cc = 6;
			tx[cc][ntx[cc]++] = U;
		}
		for (uint32_t i = 1; i < 7; i++) if (ntx[i]) {
			register uint32_t* tu = tx[i];
			for (uint32_t j = 0; j < ntx[i]; j++) {
				register uint32_t U = tu[j], nU = ~U;
				for (uint32_t k = 0; k < nt3b; k++) {
					if (!(t3b[k] & nU)) { U = 0; break; }
				}
				if (U)t3b[nt3b++] = U;
			}
		}
		if (!nt3b) {
			if (VerifyValid()) return;
			if (!Valid3_1_3(spb[3].all_previous))goto next2;
			// fresh uas see what to do
			cout << " fresh uas last step 1_3 seee what to do" << endl;
			aigstop = 1;
			return;
		}
		//Dumpt3b();
		//SP3 sp3;
		//sp3.active = BIT_SET_27 ^ crh.mycritb3.assigned;
		//sp3.all_previous = crh.mycritb3.assigned;
		GoB3Expand_4_x(spb[3]);
		if (clean_valid_done == 2) return;
		goto next3;
	}

	goto next3;// not expected
}

void G17B::GoB3Expand_4_x(SP3 spe) {
	int ntoass = ncluesb3 - _popcnt32(spe.all_previous);
	int locdiag = 0;// = 2 * (p_cpt2g[8] == 933);
	//if (ntoassb3 > 4) locdiag = 1;
	//if (op.ton)if (p_cpt2g[7] == op.f7) locdiag = op.ton;
	if (knownt >= 20)locdiag = 2;
	uint64_t limspot = ntoass - 1, limm1 = limspot - 1;
	if (locdiag) {
		cout <<Char27out(spe.all_previous) << " expand entry nto ass =" << ntoass
			<< " limspot= " << limspot << " " << limm1
			<< " clean_valid_done=" << clean_valid_done
			<< " ncluesb3=" << ncluesb3
			<< " nt3b=" << nt3b << endl;
		Dumpt3b();
		cout << Char27out(spe.active) << " active entry" << endl;
		if (ntoass < 1) {	cout << " bug ntoass" << endl; return;	}
	}
	SP3 spb[12];
	register SP3* s, * sn;
	register uint64_t ispot;
	s = spb;
	spb[0] = spe;
	s->indtw3 = 0;// initial first ua
	s->possible_cells = t3b[0];
	if (locdiag)
		cout << Char27out(t3b[0]) << " ua start"  << endl;

next:
	ispot = s - spb;
	{// catch and apply cell in bitfields
		register int cell;
		register uint32_t p = s->possible_cells;
		if (!p)if (--s >= spb)goto next; else return;
		bitscanforward(cell, p);
		register uint32_t bit = 1 << cell;
		s->possible_cells ^= bit;
		s->active ^= bit;
		sn = s + 1; *sn = *s;
		sn->all_previous |= bit;
	}
	if (locdiag)
		cout << Char27out(sn->all_previous) << " ispot =" << ispot << endl;

	if (ispot < limm1) { 
		register uint32_t F = sn->all_previous;
		for (uint32_t i = s->indtw3 + 1; i < nt3b; i++) {
			register uint32_t U = t3b[i];
			if (!(U & F)) {
				U &= s->active;
				if (!U)goto next;//dead branch
				sn->possible_cells = U;
				if (locdiag)	cout << Char27out(U) << "next ua"  << endl;
				sn->indtw3 = i;
				s++; // switch to next spot
				goto next;
			}
		}
		if (locdiag)	cout  << "no  ua" << endl;

		// no ua available must check( in 18 mode can  be valid)
		if (VerifyValid())  return; 
		if (IsValidB3(sn->all_previous)) {
			sn->possible_cells = anduab3&sn->active;
			sn->indtw3 = nt3b;// no more ua later
			s++;	goto next;// switch to next spot
		}
		else 	Out17(sn->all_previous);//can be valid 17 in 18 mode
	}

	// this is the last step must hit all pending uas
	{ // find next cells hitting all uas
		int aig = 1;
		register uint32_t andw = sn->active;
		register uint32_t F = sn->all_previous;
		{// and still valid   uas
			for (uint32_t i = s->indtw3 + 1; i < nt3b; i++) {
				register uint32_t U = t3b[i];
				if (!(U & F)) { // not hit ua
					if (!(andw &= U))goto next;//dead branch
					aig = 0;
				}
			}
		}
		if (locdiag)
			cout << Char27out(andw) << " andw aig=" << aig << endl;

		// no more ua or "and all uas" not empty
		if (VerifyValid())  return;
		if (aig) {// no ua could be a 17 valid in 18 mode 
			if (Valid3mm(F)) 	andw &= anduab3;			 
			else { Out17(F);	goto next; }// this is a valid 17
		}
		if (locdiag)
			cout << Char27out(andw) << " andw to go nt3more="<< nt3more << endl;
		register int cell;
		while (bitscanforward(cell, andw)) {
			register uint32_t bit = 1 << cell;
			andw ^= bit; //clear bit
			uint32_t nrt3 = nt3more;// don't touch old ??
			if (Valid3mm(F | bit)) 	andw &= anduab3;			
			else 	Out17(F | bit);			
		}
		goto next;
	}
	goto next;// never here
}




//=========brute force specific to this
int ZHOU::CallCheckB3( uint32_t bf, int nogo) {// 17 search mode
	zh_g.diag = nogo;
	memcpy(this, zhoustart, sizeof zhoustart);
	misc.SetAll_0();
	BF128 dca[9];
	memset(dca, 0, sizeof dca);
	int digitsbf = 0;
	if (nogo) {
		for(int i=0;i<81;i++) cout<< g17b.grid0[i]+1;
		cout << endl;
		cout << Char54out(g17b.myb12); cout << Char27out(bf);
	}

	{// cells in bands 1+2		
		int cell;
		register uint64_t U = g17b.myb12;
		while (bitscanforward64(cell, U)) {
			U ^= (uint64_t)1 << cell;
			int  digit = g17b.grid0[cell];
			int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
			digitsbf |= 1 << digit;
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}

	{
		uint32_t cc;
		register int x = bf;
		while (bitscanforward(cc, x)) {
			x ^= 1 << cc; //clear bit
			int cell = cc + 54, digit = g17b.grid0[cell];
			digitsbf |= 1 << digit;
			int xcell = cc + 64; // the cell value in 3x32 of a 128 bits map
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}
	zh_g2.s17_b3_mini = 1;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled

	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	if (nogo) {
		ImageCandidats();
		//return 0;
	}

	//__________end assign last lot start solver
	zhgxn.nua = 0;
	zh_g.go_back = 0;	zh_g.nsol = 0; zh_g.lim = 1;// modevalid is set to  1
	if(nogo)ImageCandidats();
	int ir = Full17Update();
	if (nogo) {
		ImageCandidats();// return 0;
	}

	if (ir == 2) return 0;// solved can not be multiple
	Guess17(0);
	return zhgxn.nua;
}

int ZHOU::PartialInitSearch17(uint32_t* t, int n) {
	zh_g2.digitsbf = 0;
	memset(zh_g2.Digit_cell_Assigned, 0, sizeof zh_g2.Digit_cell_Assigned);
	memcpy(this, zhoustart, sizeof zhoustart);
	for (int icell = 0; icell < n; icell++) {
		int cell = t[icell], digit = g17b.grid0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
		zh_g2.digitsbf |= 1 << digit;
		if (FD[digit][0].Off(xcell))  return 1;// check not valid entry
		Assign(digit, cell, xcell);
		zh_g2.Digit_cell_Assigned[digit].Set(xcell);
	}
	//cout << "ok init" << endl;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | zh_g2.Digit_cell_Assigned[i];
	return 0;
}

int ZHOU::Apply17SingleOrEmptyCellsB12() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched
	register uint64_t R1, R2, R3, R4;
	{
		register uint64_t* P = FD[0][0].bf.u64, M = *P;
		R1 = M;
		P += 4; M = *P;	                            R2 = R1 & M;  R1 |= M;
		P += 4; M = *P;              R3 = R2 & M;   R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 = R3 & M;  R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P;	R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
	}
	if (cells_unsolved.bf.u64[0] & (~R1)) 	return 1; // empty cells
	R1 &= ~R2; // now true singles
	R1 &= cells_unsolved.bf.u64[0]; // these are new singles
	if (R1) {// singles to apply
		uint32_t xcell;
		while (bitscanforward64(xcell, R1)) {
			R1 ^= (uint64_t)1 << xcell;
			uint32_t cell = From_128_To_81[xcell];
			for (int idig = 0; idig < 9; idig++) {
				if (FD[idig][0].On(xcell)) {
					Assign(idig, cell, xcell);
					goto nextr1;
				}
			}
			return 1; // conflict with previous assign within this lot
		nextr1:;
		}
		zh_g.single_applied = 1;
		return 0;
	}
	// no single store apply  pair in priority ??
	R2 &= ~R3; // now true singles
	if (!R2) {
		R3 &= ~R4;
		if (R3) R2 = R3;
		else R2 = R4;
	}
	bitscanforward64(zh_g2.xcell_to_guess, R2);
	return 0;
}
int ZHOU::Apply17SingleOrEmptyCellsB3() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched
	register uint32_t R1, R2, R3, R4;
	{
		register uint32_t* P = &FD[0][0].bf.u32[2], M = *P;
		R1 = M;
		P += 8; M = *P;	                            R2 = R1 & M;  R1 |= M;
		P += 8; M = *P;              R3 = R2 & M;   R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 = R3 & M;  R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P;	R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
	}
	if (cells_unsolved.bf.u32[2] & (~R1)) 	return 1; // empty cells
	R1 &= ~R2; // now true singles
	R1 &= cells_unsolved.bf.u32[2]; // these are new singles
	if (R1) {// singles to apply in band 3
		uint32_t xcell;
		while (bitscanforward(xcell, R1)) {
			R1 ^= (uint64_t)1 << xcell;
			uint32_t cell = xcell + 54;
			xcell += 64;
			for (int idig = 0; idig < 9; idig++) {
				if (FD[idig][0].On(xcell)) {
					Assign(idig, cell, xcell);
					goto nextr1;
				}
			}
			return 1; // conflict with previous assign within this lot
		nextr1:;
		}
		zh_g.single_applied = 1;
		return 0;
	}
	// no single store apply  pair in priority ??
	R2 &= ~R3; // now true pairs
	if (!R2) {
		R3 &= ~R4; // now true tripletss
		if (R3) R2 = R3;
		else R2 = R4;
	}
	bitscanforward(zh_g2.xcell_to_guess, R2);
	zh_g2.xcell_to_guess += 64;
	return 0;
}
int ZHOU::Full17Update() {
	if (zh_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		//if(diag &&cells_unsolved.bf.u32[2]==0)ImageCandidats();
		if (!Unsolved_Count()) return 2;
		if (cells_unsolved.bf.u32[2]) {// fill B3 first
			if (Apply17SingleOrEmptyCellsB3())	return 0; //  empty cell or conflict singles in cells
		}
		else {
			if ((!ISFALSEON))return 0;
			if (Apply17SingleOrEmptyCellsB12())	return 0;
		}
		if (!zh_g.single_applied)	break;
	}
	return 1;
}
void ZHOU::Compute17Next(int index) {
	int ir = Full17Update();
	if (zh_g.diag) {
		cout << "index=" << index << endl;
		ImageCandidats();
	}
	if (!ir) return;// locked 
	if (ir == 2) {//solved
		if (index) {// store false as ua
			BF128  wua;
			int* sol = g17b.grid0;
			wua.SetAll_0();
			for (int i = 0; i < 81; i++) {
				int d = sol[i];
				if (FD[d][0].Off_c(i))	wua.Set_c(i);
			}
			if (wua.isNotEmpty()) {
				//if (zh_g.diag) {
					//cout << "zhgxn.nua=" << zhgxn.nua << endl;
					//wua.Print3(" ");
				//}
				int cc = _popcnt32(wua.bf.u32[2]);
				if ((!zhgxn.nua) || cc < 3)
					zhgxn.tua[zhgxn.nua++] = wua;
				if (cc < 3 || zhgxn.nua>5)	zh_g.go_back = 1;
			}
		}
		return;
	}
	Guess17(index);// continue the process
}
void ZHOU::Guess17(int index) {
	if (zh_g.go_back) return;
	int xcell = zh_g2.xcell_to_guess,
		cell = From_128_To_81[xcell],
		digit = g17b.grid0[cell];
	// true first if possible
	if (FD[digit][0].On(xcell)) {
		ZHOU* mynext = (this + 1);
		*mynext = *this;
		mynext->SetaCom(digit, cell, xcell);
		mynext->Compute17Next(index + 1);
		if (zh_g.go_back) return;
	}
	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue;
		if (ISFALSEON && zhgxn.nua) continue;
		if (FD[idig][0].On(xcell)) {
			ZHOU* mynext = (this + 1);
			*mynext = *this;
			if (cell >= 54)mynext->ISFALSEON++;
			mynext->SetaCom(idig, cell, xcell);
			mynext->Compute17Next(index + 1);
			if (zh_g.go_back) return;
		}
	}
}



void G17B::Out17(uint32_t bfb3) {
	if ((_popcnt32(bfb3) + _popcnt64(myb12)) > 18) {
		cout << Char54out(myb12) << "\t";
		cout <<Char27out(bfb3) << " more than 18 [3] " << p_cpt2g[3]
			<< " [4] " << p_cpt2g[4] << " [7] " << p_cpt2g[7] << endl;

		DumpPotential();
		aigstop = 1;
		return;
	}
	if (op.out_one) if (myband3->poutdone++) return;	
	char ws[128];
	strcpy(ws, empty_puzzle);
	{// cells in bands 1+2		
		int cell;
		register uint64_t U = myb12;
		while (bitscanforward64(cell, U)) {
			U ^= (uint64_t)1 << cell;
			ws[cell] = grid0[cell] + '1';
		}
	}

	for (int i = 0, bit = 1; i < 27; i++, bit <<= 1)if (bfb3 & bit)
		ws[i+54]= grid0[i+54] + '1';
	if (op.ton ) {
		cout << ws << "\t\t one sol   "  << "[3] " << p_cpt2g[3] 
			<< " [4] " << p_cpt2g[4] << " [7] " << p_cpt2g[7] << endl;
		if (op.ton > 1)if (p_cpt2g[7] == op.f7) {
			cout << "nt3_2=" << nt3_2 << endl;
			for (uint32_t i = 0; i < nt3_2; i++)
				cout << Char27out(t3_2[i]) << endl;
		}
	}



	sprintf(&ws[81], ";%5d;%3d;%3d;%3d",(int)( genb12.nb12 / 64), genb12.i1t16, genb12.i2t16, t416_to_n6[myband3->i416]);     
	fout1 << ws << endl;
	/*
	fout1 << ws << ";";
	fout1.width(10);
	fout1 << genb12.nb12 / 64 << ";";
	fout1.width(5);
	fout1 << genb12.i1t16 << ";";
	fout1.width(5);
	fout1 << genb12.i2t16 << ";";
	fout1.width(5);
	fout1 << t416_to_n6[myband3->i416] << endl;
	*/
	a_17_found_here++;

}