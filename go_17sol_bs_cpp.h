

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
						Bfs = Bf & stack,
						Bfc = (Bf | (Bf >> 9) | (Bf >> 18)) & 0777;// colummn pattern

					uint32_t c1, c2, r1, r2;
					bitscanforward(c1, Bfs);	bitscanreverse(c2, Bfs);
					r1 = c1/9; r2 = c2 / 9;
					if (r1 == r2) {// socket 2
						gsock2.setBit(i);
						w.valid = 1;
						b3.g.gsocket2.setBit(i);
						b3.g.pat2[i] = Bfs;
						int imini = 3 * r1 + ist, mask = 7 << (3 * imini);
						b3.g.ua2_imini[i] = imini;
						int bit27 = mask ^ Bfs, i27;
						bitscanforward(i27, bit27);
						b3.g.ua2bit27[i] = bit27;
						b3.i_27_to_81[i27] = i;
						b3.i_81_to_27[i] = i27;
						b3.g.pat2_27[i27] = Bfs;
						continue;
					}
					// can be a 4/6
					int ncol = _popcnt32(Bfc);
					if (ncol > 5) continue;
					if (ncol == 4) {// this is a ua6
						b3.g.gsocket6.setBit(i);
						gsock2.setBit(i);
						b3.g.pat2[i] = Bf;
						continue;
					}
					// 5 column find the ua4
					for (int is2 = 0,mask=7, st2 = 07007007; is2 < 3; is2++,mask<<=3,st2<<=3)
						if (is2 != ist) 		if (_popcnt32(mask & Bfc) != 1)Bf &= ~st2;
					// must also be 2 rows 
					int maskr = 0777 << (9 * r1) | 0777 << (9 * r2);
					if (Bf & ~maskr) continue;
					b3.g.gsocket4.setBit(i);
					gsock2.setBit(i);
					b3.g.pat2[i] = Bf;

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
			if (op.ton > 2) {
				for (int j = 0; j < 81; j++)if (b3.g.gsocket4.On(j))
					cout << j << "\t" <<Char27out( b3.g.pat2[j]) << " ua4" << endl;
			}

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
//#define PUZTEST 1483
void G17B::Start() {// processing an entry 
#ifdef PUZTEST
	if (npuz > PUZTEST) return;
	if (npuz == PUZTEST) {
		memset(p_cpt2g, 0, sizeof p_cpt2g);
		op.ton = 3;
		op.f3 = 1;
		op.f4 = 3;
		op.f7 = 0;
		op.f10 = sgo.vx[9] = 0;
		op.f10 = sgo.vx[9] = 15;
		op.dv12 = op.dv3 = 1;
	}
	else {
		op.ton = 0;
		op.f3 = 0;
		op.f4 = 0;
		op.f7 = 0;
		op.f10=sgo.vx[9] = 0;
		op.dv12 = op.dv3 = 0;

	}
#endif
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
				if (cc3 > 3) {
					if(cc3==4)	mb3.Addm4(wd);
					else mb3.Addmm(wd);
				}
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
	}

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
				if (ncol != 4)  mb3.Addmm(w); 
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
						if (nplug > 2) { 
							if(cc2==4)mb3.Addm4(w); 
							else mb3.Addmm(w);
							continue; 
						}
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

	//ng2=guah.CutG2(30);
	//if (op.ton)cout << "ng2 after cut30=" << ng2 << endl;
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
	//ng3 = guah.CutG3(20);
	//if (op.ton)cout << "ng3 after cut20=" << ng3 << endl;
	//guah.DumpOne3(46);
	//guah.Dump2all3();
}


//__________ end of initial uas harvest switch to 54 mode

void G17B::StartAfterUasHarvest() {

	guah54.Build();
	cout << "guas first ng2=" << guah54.nag2 << " ng3=" << guah54.nag3 << endl;
	//guah54.DumpA2();
	//return;
	t54b12.Build_ta128(tuasb12.tua, tuasb12.nua);
	if (op.ton) {
		if (sgo.bfx[1] & 8)t54b12.DebugA();
		/*
		STD_B3& b = genb12.bands3[0];
		cout <<b.band<< "status  pour pat 2" << endl;
		((b.g.gsocket2 | b.g.gsocket4) | b.g.gsocket6).Print("g2 4 6");
		for (int i = 0; i < 81; i++) {
			int p = b.g.pat2[i];
			if (p) cout << Char27out(p) << " i=" << i << endl;
		}
		b.DumpTgmm();
		*/
	}
	//if (op.known > 1) genb12.bands3[0].DumpTgm();
	//return;
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
	/*
	int tc[3], ntc = 0;
	{// // build table of clues
		int cell;
		register uint64_t U = s.all_previous_cells;
		while (bitscanforward64(cell, U)) {
			U ^= (uint64_t)1 << cell;
			tc[ntc++] = cell;
		}
	}
	*/
	uint32_t* tc = g17b.tc_1_3;


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
int T54B12::Build_tc128( SPB03A& s6) {
	/*
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
	*/
	uint32_t *tc = g17b.tc_4_6;

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
	return 0;
}
int T54B12::Build_td128(SPB03A& s9) {
	uint32_t* tc = g17b.tc_7_9;
	uint32_t lastbloc = t54b12.ncblocs;
	tvw[0] = s9.v;
	for (uint32_t i = 1; i <= lastbloc; i++) {
		TUVECT& vv = tc128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 0; ic < 3; ic++)
			v &= vc[tc[ic]];
		tvw[i] = v;
	}
	// apply active on still valid uas and flag by size
	memset(vsize, 0, sizeof vsize);
	uint32_t ntw = 0;
	{
		register uint64_t Ac = s9.active_cells,
			nb1=_popcnt64(s9.all_previous_cells &BIT_SET_27);
		if (nb1 > 6 || nb1 < 3) return 1; // no 66
		if (nb1 == 6) Ac &= ~(uint64_t)BIT_SET_27;// only b2
		if (nb1 == 3) Ac &= BIT_SET_27;// only b1
		for (uint32_t i = 0; i <= lastbloc; i++) {
			register uint64_t* t = tc128[i].t;
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
					if (ntw >= 128 * UADNBLOCS) break;
				}
				else break;
			}
			if (ntw >= 128 * UADNBLOCS) break;
		}
	}
	//__Build the reduced set of UAs vectors clean redundancy
	InitD();// be sure to have D ready for fresh uas
	uint32_t nbl64 = (ntw + 63) >> 6, x;
	for (int i1 = 1; i1 < 19; i1++) {
		uint64_t* tb64 = vsize[i1]->bf.u64;
		for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
			register uint64_t V = tb64[i2];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				register uint64_t U = tw[x + (i2 << 6)];
				// check redundancy in td128[0]
				if (t54b12.IsNotRedundantD(U))
					AddD(U);
			}
		}
	}
	int ii = nd128 >> 5; if (ii > 3)ii = 3;
	p_cpt2g[84+ii]++;
	//if (ii>2) { cout << p_cpt2g[90] << " ";	DebugD(); }

	uint64_t nb1 = _popcnt64(s9.all_previous_cells & BIT_SET_27);
	if (nb1 == 6) p_cpt2g[82]++;;// only b2
	if (nb1 == 3) p_cpt2g[83]++;;// only b1
	if (!nd128) {
		//cout << "empty D p_cpt2g[91]=" << p_cpt2g[91] << endl;
		return 1;
	}

	if (nd128 < 128)nd128 = 128;// fresh uas in add
	return 0;
}

void GUAH54::Build() {// cut to 30 switch to 54 killer
	InitA();
	for (int i81 = 0; i81 < 81; i81++) {
		if (g17b.gsock2.On(i81)) {
			GUAH::GUA& g0 = guah.tg2[i81];
			uint32_t n = g0.nua;
			//if (n > 30) n = 30;
			for (uint32_t i = 0; i < n; i++) {
				register uint64_t U = g0.tua[i];
				if (n > 30 && _popcnt64(U) > 16) break;
				U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
				AddA2(U, i81, 0);// take all of them
			}
		}
	}
	for (int i81 = 0; i81 < 81; i81++) {
		if (g17b.gsock3.On(i81)) {
			GUAH::GUA& g0 = guah.tg3[i81];
			uint32_t n = g0.nua;
			//if (n > 30) n = 30;
			for (uint32_t i = 0; i < n; i++) {
				register uint64_t U = g0.tua[i];
				if (n > 30 && _popcnt64(U) > 16) break;
				U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
				AddA3(U, i81, 0);// take all of them
			}
		}
	}
	//cout << " gua2 table" << endl; DumpA2();
	//cout << " gua3 table" << endl; DumpA3();

}
void GUAH54::Build6(int debug) {
	wc = 0;// no count filter
	uint32_t *tc=g17b.tc_1_6;// build table of cells 
	register uint64_t 	A = g17b.ac_cl6;
	if (debug) {
		cout << "debug Build6 for gua2 21" << endl;
		cout << Char54out(g17b.bf_cl6) << " cl" << endl;
		cout << Char54out(g17b.ac_cl6) << " ac" << endl;
	}
	InitB();
	for (uint32_t i = 0; i <= nbag2; i++) {
		G64& g = ag2[i];
		register uint64_t v = g.v0;
		for (uint32_t j = 0; j < 6; j++)
			v &= g.vc[tc[j]];
		register int k;
		while (bitscanforward64(k, v)) {
			v ^= (uint64_t)1 << k;
			wbf = g.t12[k] & A;
			wi = g.ti81[k];
			Add2B();
		}
	}
	for (uint32_t i = 0; i <= nbag3; i++) {
		G64& g = ag3[i];
		register uint64_t v = g.v0;
		for (uint32_t j = 0; j < 6; j++)
			v &= g.vc[tc[j]];
		register int k;
		while (bitscanforward64(k, v)) {
			v ^= (uint64_t)1 << k;
			wbf = g.t12[k] & A;
			wi = g.ti81[k];
			Add3B();
		}
	}
	//if (!(p_cpt2g[4] & 127)) { DumpB2(); DumpB3(); }
}
void GUAH54::GetG2G3(BF128& g2,BF128& g3) {
	uint32_t *tc=g17b.tclues6p,
		ntc = g17b.nclues6p ;// build table of cells 
	g2.SetAll_0();
	for (uint32_t i = 0; i <= nbbg2; i++) {
		G64B& g = bg2[i];
		register uint64_t v = g.v0;
		for (uint32_t j = 0; j < ntc; j++)
			v &= g.vc[tc[j]];
		register int k;
		while (bitscanforward64(k, v)) {
			v ^= (uint64_t)1 << k;
			g2.setBit(g.ti81[k]);
		}
	}		
	g3.SetAll_0();
	for (uint32_t i = 0; i <= nbbg3; i++) {
		G64B& g = bg3[i];
		register uint64_t v = g.v0;
		for (uint32_t j = 0; j < ntc; j++)
			v &= g.vc[tc[j]];
		register int k;
		while (bitscanforward64(k, v)) {
			v ^= (uint64_t)1 << k;
			g3.setBit(g.ti81[k]);
		}
	}
}
void GUAH54::GetG2G3_12(BF128& g2, BF128& g3) {
	uint32_t* tc = g17b.tc_10_12;//  table of cells 
	g2.SetAll_0();
	for (uint32_t i = 0; i <= nbbg2; i++) {
		G64B& g = bg2[i];
		register uint64_t v = g.v9;
		for (uint32_t j = 0; j <3; j++)
			v &= g.vc[tc[j]];
		register int k;
		while (bitscanforward64(k, v)) {
			v ^= (uint64_t)1 << k;
			g2.setBit(g.ti81[k]);
		}
	}
	g3.SetAll_0();
	for (uint32_t i = 0; i <= nbbg3; i++) {
		G64B& g = bg3[i];
		register uint64_t v = g.v9;
		for (uint32_t j = 0; j < 3; j++)
			v &= g.vc[tc[j]];
		register int k;
		while (bitscanforward64(k, v)) {
			v ^= (uint64_t)1 << k;
			g3.setBit(g.ti81[k]);
		}
	}
}

void GUAH54::SetupV9() {
	uint32_t* tc = g17b.tc_7_9;// table of cells 
	for (uint32_t i = 0; i <= nbbg2; i++) {
		G64B& g = bg2[i];
		register uint64_t v = g.v0;
		for (uint32_t j = 0; j <3; j++)
			v &= g.vc[tc[j]];
		g.v9 = v;
	}
	for (uint32_t i = 0; i <= nbbg3; i++) {
		G64B& g = bg3[i];
		register uint64_t v = g.v0;
		for (uint32_t j = 0; j <3; j++)
			v &= g.vc[tc[j]];
		g.v9 = v;
	}
}

//_______________________ expand common step 1_9 clues

///________ start expand uas bands 12
/// target 10 clues 3+3  +4
/// target 11 clues 3+3  +5
/// target 12 clues 3+3  +3+3
/// end of early steps skrinking and restructuring uas bands 1+2
/// and building a reduced GUAs table 

void DumpPotential(int det = 0) {
	cout << "DumpPotential() ";
	for (int i = 0; i < 6; i++) cout << ntbelow[i] << " ";
	cout << endl;

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
		tc_1_3[0] = cell;
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
		tc_1_3[1] = cell2;
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
		tc_1_3[2] = cell3;
		sp2.possible_cells ^= bit;
		sp2.active_cells ^= bit;
		sp3 = sp2;
		sp3.all_previous_cells |= bit;
		sp3.v &= tuv128.vc[cell3];
		//if (!(sp3.possible_cells = twu[sp3.v.getFirst128()] & sp0.active_cells))goto next3;
		p_cpt2g[3]++;
		if (t54b12.Build_tb128(sp3)) goto next3;
		Set3(sp3);
		if (op.known > 1) {
			//if (op.known)genb12.bands3[0].DumpGuam(1);
			if (!((~pk54) & sp3.all_previous_cells)) {
				cout << Char54out(sp3.all_previous_cells) << " expected 3 cells "
					<<tc_1_3[0] <<" " << tc_1_3[1] << " " <<tc_1_3[2] << endl;
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
				cout << "call 4_6 good path cells "<< tc_1_6[0] << " " << tc_1_6[1] << " " << tc_1_6[2] << endl;
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
		tc_4_6[0] = cell;
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
		tc_4_6[1] = cell2;
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
		tc_4_6[2] = cell3;
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
		//if (t54b12.Build_tc128(s3,sp6)) goto next6;
		if (t54b12.Build_tc128( sp6)) goto next6;
		p_cpt2g[4]++;
		Set6(sp6);
		if (op.known > 1) {
			if (knownt >= 9)return;
			if (op.ton > 2)cout << "ua6c  nc128=" << t54b12.nc128
				<< "[4] " << p_cpt2g[4] << endl;
			if (!((~pk54) & sp6.all_previous_cells)) {
				cout << Char54out(bf_cl6) << " expected 6 [4] " 
					<<p_cpt2g[4] << endl;
				cout << Char54out(ac_cl6) << " expected 6 [4] "
					<< p_cpt2g[4] << endl;
				cout << "cells " << tc_1_6[0] << " " << tc_1_6[1] << " "
					<< tc_1_6[2] << " " << tc_1_6[3] << " "
					<< tc_1_6[4] << " " << tc_1_6[5] << endl;
				knownt = 6;
				Expand_7_9(sp6);
				aigstop = 1;
				return;;
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
	// ____ build a reduced table of uas/guas for band 3
	guah54.Build6();
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		genb12.bands3[ib3].Set6clues(s6.all_previous_cells);
	}
	//if (knownt == 6)guah54.DumpB2(1);
	SPB03A sp6,sp7,sp8,sp9,sp10;
	T54B12::TUVECT& tuv128 = t54b12.tc128[0];
	uint64_t* twu = tuv128.t;
	uint32_t * tcells=g17b.tc_7_9;
	int locdiag = 0;
	if (op.ton) {
		if (op.f3) {
			if (p_cpt2g[3] == op.f3&& (!op.f4)) {
				cout << Char54out(s6.all_previous_cells) << " 6clues [4]" << p_cpt2g[4]
					<< " nc128=" << t54b12.nc128;
				cout << Char54out(twu[0]) << " [7]"<< p_cpt2g[7] << endl;

			}
		}
		if (op.f4) {
			if (p_cpt2g[4] == op.f4) {
				cout << Char54out(s6.all_previous_cells) << "call 7_9 good path" << endl;
				cout << Char54out(s6.active_cells) << " active " << endl;
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
	//memset(tandbelow, 255, sizeof tandbelow);//7 8 9 10 11 full


next7:	//_______ add clue 7
	{
		uint64_t p = sp6.possible_cells;
		if (!p) {	
			if (p_cpt2g[4] == op.f4) cout << " 7;8;9;10 " << ntbelow[0] << " " << ntbelow[1] << " "
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
				//tandbelow[0] &= w.bf.u64[0];
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
				//tandbelow[1] &= w.bf.u64[0];
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
				//tandbelow[2] &= w.bf.u64[0];
				p_cpt2g[56]++;
				goto next9;
			}
			else U = ua_ret7p;
		}
		sp9.possible_cells = U & sp9.active_cells;// never empty
		if (knownt == 9)cout << Char54out(sp9.possible_cells) << "possible after 9" << endl;
	}
	if (op.t18 && op.p2) {// must expand 10_12
		if(!sp9.possible_cells)goto next9;// dead branch
		Set9(sp9);
		if (locdiag)cout << Char54out(sp9.all_previous_cells) << "  9 call 10 12 before build " << endl;
		if (! t54b12.Build_td128(sp9))Expand_10_12( sp9);
		if (aigstop) return;
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
				//tandbelow[3] &= w.bf.u64[0];
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
			cb3.ncl = 11;
			cb3.cbs.Init(myb12, 11);
			GoCallB3( cb3);
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
			return;
			if (knownt > 11)return;
			if (ntbelow[0]) Go_7_12(); // push to 12 clues 666
			if (ntbelow[1]) Go_8_12();
			if (ntbelow[2]) Go_9_12();

		}
	}
	else {
		if (op.p1) {
			if (ntbelow[0]) Go_7_10(); // do 7 clues then more
			if (ntbelow[1]) Go_8_10(); // do 8 clues then more
			if (ntbelow[2]) Go_9_10(); // do 9 clues then more
		}
		else {
			if (ntbelow[0]) Go_7_11_17(); // push to 11 clues 656
			if (ntbelow[1]) Go_8_11_17();
			if (ntbelow[2]) Go_9_11_17();
		}
	}
}
//#define T90 162366
void G17B::Expand_10_12(SPB03A& s9) {
	if (aigstop) return;
	p_cpt2g[90]++;
	int locdiag = 0;
	if (!((~pk54) & bf_cl9)) 		locdiag = 1;
#ifdef  T90
	if (op.f4) {
		if (p_cpt2g[4] == op.f4) {
			if (op.ton > 2) {
				cout << "  [90] " << p_cpt2g[90] << endl;
				if (p_cpt2g[90] > T90) { aigstop = 1; return; }
				if (p_cpt2g[90] == T90) {
					t54b12.DebugD();
					locdiag = 2;
				}
			}
		}
	}
#endif


	//if (op.f4 == p_cpt2g[4])cout << Char54out(s9.all_previous_cells) << "  [90] "<< p_cpt2g[90] << endl;
	SPB03A   sp9, sp10,sp11,sp12;
	T54B12::TUVECT& tuv128 = t54b12.td128[0];
	uint64_t* twu = tuv128.t;
	sp9 = s9;
	sp9.possible_cells = twu[0];
	sp9.v = tuv128.v0;
	memset(&ntbelow[3], 0, 3* sizeof ntbelow[0]);// 10 11  12
next10:	//_______ add clue 10
	{
		uint64_t p = sp9.possible_cells;
		if (!p) {	
#ifdef  T90
		//	if (locdiag) return;
#endif
			EndExpand_10_12();	return;		}
		bitscanforward64(tc_10_12[0], p);
		register uint64_t bit = (uint64_t)1 << tc_10_12[0];
		sp9.possible_cells ^= bit;
		sp9.active_cells ^= bit;
		sp10 = sp9;
		sp10.all_previous_cells |= bit;
		if (op.known > 1) {
			if (!((~pk54) & sp10.all_previous_cells)) {
				cout << Char54out(sp10.all_previous_cells) << " expected 10 " << endl;
				knownt = 10;
			}
		}
		if (locdiag > 1) {
			cout << Char54out(sp10.all_previous_cells) << " new 10 " << endl;
		}
		// update vectors and get next ua
		register uint64_t U = 0;
		sp10.v &= tuv128.vc[tc_10_12[0]];
		if (sp10.v.isNotEmpty()) U = twu[sp10.v.getFirst128()];
		else if(t54b12.ndblocs) {// usually "more" not if stored  
			register T54B12::TUVECT& tv1 = t54b12.td128[1];
			BF128 w1 = tv1.v0 & tv1.vc[tc_10_12[0]];
			if (w1.isNotEmpty()) U = tv1.t[w1.getFirst128()];
		}
		if (!U) {// this is a valid 10  
			if (locdiag>1)cout << " potential 10 n="<< ntbelow[3] << endl;

			if (!IsValid7pbf(sp10.all_previous_cells)) {
				BF128 w;// valid 10
				w.bf.u64[0] = sp10.all_previous_cells;
				w.bf.u64[1] = sp10.active_cells;
				tbelow10[ntbelow[3]++] = w;
				p_cpt2g[57]++;
				goto next10;
			}
			U = ua_ret7p;
			if (locdiag>1)cout<<Char54out(ua_ret7p) << " not  10 " << endl;
		}
		sp10.possible_cells = U & sp10.active_cells;// never empty
		uint64_t nb1 = _popcnt64(sp10.all_previous_cells & BIT_SET_27);
		if (nb1 == 6) sp10.possible_cells&= ~(uint64_t)BIT_SET_27;// only b2
		if (nb1 == 4) sp10.possible_cells &= BIT_SET_27;// only b1

	}

next11: //add clue 11  
	{
		uint64_t p = sp10.possible_cells;
		if (!p) goto next10;
		bitscanforward64(tc_10_12[1], p);
		register uint64_t bit = (uint64_t)1 << tc_10_12[1];
		sp10.possible_cells ^= bit;
		sp10.active_cells ^= bit;
		sp11 = sp10;
		sp11.all_previous_cells |= bit;
		if (op.known > 1) {
			if (!((~pk54) & sp11.all_previous_cells)) {
				cout << Char54out(sp11.all_previous_cells) << " expected 11 " << endl;
				knownt = 11;
			}
		}


		if (locdiag > 1) {
			cout << Char54out(sp11.all_previous_cells) << " new 11 " << endl;
		}
		register uint64_t Uand = ~0, U = 0;
		sp11.v &= tuv128.vc[tc_10_12[1]];
		if (sp11.v.isNotEmpty()) { U = 1;	tuv128.DoAnd(sp11.v, Uand); }
		if (U && (!Uand)) goto next11;
		for (uint32_t i = 1; i <= t54b12.ndblocs; i++) {
			register T54B12::TUVECT& tvi = t54b12.td128[i];
			BF128 wi = tvi.v0 & (tvi.vc[tc_10_12[1]] & tvi.vc[tc_10_12[0]]);
			if (wi.isNotEmpty()) { tvi.DoAnd(wi, Uand); U = 1; }
			if (U && (!Uand)) goto next11;
		}

		if (!U) {////can be a valid 11 uaand empty
			if (locdiag>1|| knownt == 11)cout << " potential 11 n=" << ntbelow[4] << endl;
			if (!IsValid7pbf(sp11.all_previous_cells)) {
				BF128 w;// valid 11
				w.bf.u64[0] = sp11.all_previous_cells;
				w.bf.u64[1] = sp11.active_cells;
				tbelow11[ntbelow[4]++] = w;
				p_cpt2g[58]++;
				goto next11;
			}
			Uand = anduab12;
			if (locdiag)cout << Char54out(Uand) << " not  11 " << endl;
		}
		//this is a 11 to push to 12
		if (locdiag>1 || knownt == 11)cout << Char54out(Uand) << "  <-11 to push to 12" << endl;
		uint64_t nb1 = _popcnt64(sp11.all_previous_cells & BIT_SET_27);
		if (nb1 == 6)Uand &= ~(uint64_t)BIT_SET_27;// only b2
		else Uand &= BIT_SET_27;// only b1

		Uand &= sp11.active_cells;// never empty
		if (!Uand)goto next11;


		while (bitscanforward64(tc_10_12[2], Uand)) {
			uint64_t bit2 = (uint64_t)1 << tc_10_12[2];
			Uand ^= bit2; //clear bit
			register uint64_t u12= sp11.all_previous_cells | bit2;
			if (op.known > 1) {
				if (!((~pk54) & u12)) {
					cout << Char54out(u12) << " expected 12 " << endl;
					knownt = 12;
				}
				else continue;
			}
			if (locdiag)
				cout << " add n=" << ntbelow[5] << endl;
			tfull[ ntbelow[5]++]= sp11.all_previous_cells | bit2;
		}
		goto next11;
	}
}
void G17B::EndExpand_10_12() {
	if (!(ntbelow[3] | ntbelow[4] | ntbelow[5])) return;
	p_cpt2g[91]++;	p_cpt2g[95] += ntbelow[3];
	p_cpt2g[96] += ntbelow[4];	p_cpt2g[97] += ntbelow[5];
	if(p_cpt2g[98] < ntbelow[5])p_cpt2g[98] = ntbelow[5];
	int locdiag = 0;	//if (p_cpt2g[91] == 225393) locdiag = 1;
	//if (p_cpt2g[91] > 4) return;
	if (op.known &&!((~pk54) & bf_cl9)) 		locdiag = 1;
#ifdef  T90
	if (p_cpt2g[90] == T90) {
		DumpPotential();
		//return;
		locdiag = 2;
	}
#endif

	if (locdiag) {
		cout << Char54out(bf_cl9) << " end after 9 [91] " << p_cpt2g[91]
			<< " [3] " << p_cpt2g[3] << " [4] " << p_cpt2g[4] << endl;
		DumpPotential();// ntbelow[5]=4;
	}

	guah54.SetupV9();
	if (locdiag) {
		cout << " dump of still valid 9 in gua2s cells "
			<< tc_7_9[0] << " " << tc_7_9[1] << " " << tc_7_9[2] << endl;
		cout << Char54out(bf_cl9) << "cl9" << endl;
		//guah54.DumpB2_9();
	}

	
	if (ntbelow[5]) {
		p_cpt2g[92]++;
		for (uint32_t i = 0; i < ntbelow[5]; i++) {
			myb12 = cb3.bf12 = tfull[i];
			if (locdiag)cout << Char54out(myb12) << " 12" << endl;
			if (knownt >= 11) {
				cout << Char54out(myb12) << " 12" << endl;
				if (!((~pk54) & myb12)) {
					cout << Char54out(myb12) << " expected 12 " << endl;
					knownt = 12;
				}
			}
			{
				register uint64_t F = myb12 & ~bf_cl9;
				int n = 0;// cells 9_12
				register int cell;// build table of cells 
				while (bitscanforward64(cell, F)) {
					F ^= (uint64_t)1 << cell;
					tc_10_12[n++] = cell;
				}
			}
			clean_valid_done = 0;
			cb3.ncl = 12;
			cb3.cbs.Init(myb12, 12);
			//if (locdiag && ntbelow[5]> 4) continue;		
			GoCallB3_12( cb3);
			if (knownt == 12) return;
		}
	}
	if (locdiag > 1) {
		//cout << "stop test " << endl; return;
	}
	if (ntbelow[4]) Go_11_12();	
	
	return;
	// Set up 9 clues vector in guas
	// 
	//cout <<Char54out(bf_cl9)<<" 9 common " << p_cpt2g[90] << " ";

	if (ntbelow[3]) Go_10_12();
}



//________________expand and go search 17 pass 1

void G17B::Go_9_10() {// 9 clues limit 10 clues 
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[2]; iv++) {
		BF128 ww = tbelow9[iv];
		myb12 = cb3.bf12 = ww.bf.u64[0];
		myac= ww.bf.u64[1];
		//guah54_9.Build9(myb12, myac);
		//for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
		//	genb12.bands3[ib3].BuildGuam9(myb12);
		cb3.cbs.Init(myb12, 9);
		cb3.ncl = 9;
		// try direct
		GoCallB3(cb3);
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
			GoCallB3(cb3n);
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
		GoCallB3(cb3);

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
			GoCallB3(cb3n);

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
				GoCallB3(cb3n2);
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
		GoCallB3(cb3);
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
			GoCallB3(cb3n);

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
				GoCallB3(cb3n2);
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
					GoCallB3(cb3n3);
				}

			}
		}
	}
}

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
			cout << Char54out(myb12) << "10 to test to push to 11 iv=" << iv << endl;

		}
		cb3.cbs.Init(myb12, 10);
		// try direct
		GoCallB3(cb3);
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
			if (locdiag) 
				cout << Char54out(myb12) << " 11 [7] "	<< p_cpt2g[7] << endl;
			GoCallB3(cb3n);
			if (aigstop==1)return;
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
		GoCallB3(cb3);
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
			GoCallB3(cb3n);
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
				cb3n2.ncl = 11;
				GoCallB3(cb3n2);
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
		GoCallB3(cb3);
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
			GoCallB3(cb3n);
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
				GoCallB3(cb3n2);
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
					GoCallB3(cb3n3);
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
		GoCallB3(cb3);
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
			GoCallB3(cb3n);
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
				GoCallB3(cb3n2);
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
					GoCallB3(cb3n3);
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
						GoCallB3(cb3n4);
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
			GoCallB3(cb3n);
			if (aigstop == 1) return;
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
				GoCallB3(cb3n2);
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
					GoCallB3(cb3n3);
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
						GoCallB3(cb3n4);
					}
				}
			}
		}
	}
}

//________________expand and go search 18 pass 2


void G17B::Go_11_12() {
	int locdiag = 0;
#ifdef  T90
	if (p_cpt2g[90] == T90) {
		locdiag = 2;
	}
#endif
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[4]; iv++) {
		if (locdiag) cout << "Go_11_12() in test iv=" << iv << endl;
		BF128 ww = tbelow11[iv];
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 11);
		{
			register uint64_t F = myb12 & ~bf_cl9;
			int n = 0;// cells 9_12
			register int cell;// build table of cells 
			while (bitscanforward64(cell, F)) {
				F ^= (uint64_t)1 << cell;
				tc_10_12[n++] = cell;
			}
		}



		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		Ac &= cb3.cbs.NextActive();
		if (op.known) {
			if (!((~pk54) & myb12)) {
				cout << Char54out(myb12) << " expected 11 go 11 12 iv=" << iv << endl;
				cout << Char54out(Ac) << " active  ";
				cb3.cbs.Status();
				locdiag = 1;
			}
		}
		if (locdiag>1) {
			cout << Char54out(myb12) << "  11 go 11 12 iv=" << iv << endl;
			cout << Char54out(Ac) << " active  ";
		}

		
		while (bitscanforward64(tc_10_12[2], Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << tc_10_12[2];
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			if (op.known > 1) {
				if (!((~pk54) & cb3n.bf12)) {
					cout << Char54out(cb3n.bf12) << " expected 12 in 11_12 added clue  "
						<< tc_10_12[2] << endl;
					knownt = 12;
				}
			}
			myb12 = cb3n.bf12;
			cb3n.cbs.Add(tc_10_12[2]);
			cb3n.ncl = 12;
			if (locdiag > 1) {
				cout << Char54out(myb12) << "  12 previous [7]" << p_cpt2g[7] << endl;
				cout << Char54out(Ac) << " active  ";
			}

			GoCallB3_12(cb3n);
		}
	}
}
void G17B::Go_10_12() {
	int locdiag = 0;
#ifdef  T90
	if (p_cpt2g[90] == T90) {
		locdiag = 2;
	}
#endif
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[3]; iv++) {
		if (locdiag) cout << "Go_10_12() in test iv=" << iv << endl;
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
				GoCallB3_12(cb3n2);
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
					GoCallB3(cb3n3);
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
						GoCallB3(cb3n4);
					}
				}
			}
		}
	}
}
void G17B::Go_7_12() {

}



//____________processing band 3 to the end
void GUAH54::AddA2(uint64_t bf, int i81, int cc) {
	/*
	if (i81 == 37) {
		cout << Char54out(bf) << ";";
		cout << Char54out(genb12.bands3[0].g.pat2[i81])
			<< "add g2=37 [3]" << p_cpt2g[3] << " [4] " << p_cpt2g[4]
			<< " [7] " << p_cpt2g[7] << " [10] " << p_cpt2g[10]
			<< "<<<<<<<<<<<<<<<<<<<<<<<" << endl;
	}
	*/
	if (op.dv3&& p_cpt2g[10]) {
		cout << Char54out(bf) << "\ti81="<<i81  << " added g2  [10] " << p_cpt2g[10]  << endl;
/*
		if (i81 == 13 || i81 == 71) {
			cout << Char54out(bf) << ";";
			cout << Char54out(genb12.bands3[0].g.pat2[i81])
				<< "add g2=13/71 [3]" << p_cpt2g[3] << " [4] " << p_cpt2g[4]
				<< " [7] " << p_cpt2g[7] << " [10] " << p_cpt2g[10]
				<< "<<<<<<<<<<<<<<<<<<<<<<<" << endl;
			xq.Status();
			if (p_cpt2g[7] > 103) g17b.aigstop = 1;
		}
*/
	}



	wbf = bf; wi = i81, wc = cc;
	Add2A();
	if (cc)Add2B();
}

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
				cout << Char27out(bf) << " bug no b3 IsValidB3 [7]"					
					<< p_cpt2g[7]	<< "   [8]" << p_cpt2g[8] 
					<< "nr=" << nret << endl;
				aigstop = 1;	return 1;
			}
			register uint32_t ua = w.bf.u32[2];
			t3more[nt3more++] = ua;
			anduab3 &= ua;
			if (ua & ~t3infield)  t3outseen &= ua;

			register uint64_t U = w.bf.u64[0];
			uint64_t cc0 = _popcnt64(U);
			if (cc0 > 16)continue;
			if (!cc0) {
				cout << Char27out(bf) << " bug ua nob12 [10]" << p_cpt2g[10] << endl;
				cout << Char54out(myb12) << " " << endl;
				cout << "list of uas found" << endl;
				for (uint32_t i = 0; i < zhgxn.nua; i++) {
					BF128 ww = zhgxn.tua[i];
					cout << Char2Xout(ww.bf.u64[0]) << " ";
					cout << Char27out(ww.bf.u32[2]) << " i=" << i << endl;
				}
				zhou[0].CallCheckB3(bf, 1);		aigstop = 1;	return 1;
			}

			U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
			if (cc <= 6 && cc > 1) {
				if (cc > 3) {// is it a 2 digits??
					int ndigs = 0;
					for (int i = 0; i < 9; i++)
						if (ua & myband3->dpat[i]) ndigs++;
					if (ndigs == 2) {
						int i81 = myband3->GetI81_x(ua);
						if (i81 >= 0) {
							guah54.AddA2(U, i81, (int)cc0);
							cb3.g2t.setBit(i81);
							p_cpt2g[20]++;
							continue;
						}
					}
					p_cpt2g[19]++;
					if (cc == 4) {
						myband3->Addm4(w);
						if (op.dv3) {
							cout << Char54out(U) <<"\t";
							cout << Char27out(w.bf.u32[2]) << " added 4 " << myband3->nbgm 
								<<" [10] "<< p_cpt2g[10] << endl;
						}
					}
					else {
						myband3->Addmm(w);
						if (op.dv3) {
							cout << Char54out(U) << "\t";
							cout << Char27out(w.bf.u32[2]) << " added >4 " 
								<< " [10] " << p_cpt2g[10] << myband3->nbgmm << endl;
						}
					}
					continue;
				}

				else {
					if (cc == 2) {
						int i81 = myband3->GetI81_2(w.bf.u32[2]);
						guah54.AddA2(U, i81, (int)cc0);
						cb3.g2t.setBit(i81);
						p_cpt2g[20]++;
					}
					else {
						int i81 = myband3->GetI81_3(w.bf.u32[2]);
						guah54.AddA3(U, i81, (int)cc0);
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
int G17B::Valid3_1_3(uint32_t bf) {
	if (!IsValidB3(bf))return 0;
	for (uint32_t i = 0; i < nt3more; i++) {
		xq.tout[xq.nout++] = t3more[i];
	}
	nt3more = 0;
	return 1;
}
int  G17B::Valid3mm(uint32_t bf) {
	if (!IsValidB3(bf))return 0;
	for (uint32_t i = 0; i < nt3more; i++) {
		t3b[nt3b++] = t3more[i];
		xq.tout[xq.nout++] = t3more[i];
		//cout << Char27out(t3more[i]) << " add fresh tb" << endl;
	}
	nt3more = 0;
	return 1;
}


//__________ new potential valid bands 1+2

//#define MYTEST  1512000

#define VTEST sgo.vx[9]
void G17B::GoCallB3(CALLBAND3& cb3w) {
	p_cpt2g[7]++;
	if (op.t18 && op.p2) return;// should not be called 
	if(op.f4== p_cpt2g[4] )cout << Char54out(myb12) << "  [7] "
		<< p_cpt2g[7] << endl;
	int locdiag = 0;
	if ( p_cpt2g[7]== op.f7){
		cout << Char54out(myb12) << "GoCallB3 in diag [7]  " << p_cpt2g[7] << endl;
		locdiag = 1;
		if (op.known) {
			cout << "genb12.bands3[0].nbgm "<< genb12.bands3[0].nbgm 
				<< " genb12.bands3[0].nbgmm " << genb12.bands3[0].nbgmm << endl;
			//genb12.bands3[0].DumpTgmm();
		}
	}
	if (op.f7 && p_cpt2g[7] >op.f7) { aigstop = 1; return; }
	//if (op.f7 > p_cpt2g[7]) return;
	if (op.known > 1) {
		if (!((~pk54) & myb12)) {
			cout << Char54out(myb12) << "  expected 11 [7] "<< p_cpt2g[7] << endl;
			knownt = 11;
		}
	}
	{
		register uint64_t F = g17b.myb12 & ~g17b.bf_cl6; // cells 6 to x
		nclues6p = 0;
		register int cell;// build table of cells 
		while (bitscanforward64(cell, F)) {
			F ^= (uint64_t)1 << cell;
			tclues6p[nclues6p++] = cell;
		}
		if (g17b.knownt == 12|| locdiag ) {
			cout <<Char54out(g17b.myb12 & ~g17b.bf_cl6)
				<< " entry go expected ok clean=" << g17b.clean_valid_done
				<< " aigstop= " << g17b.aigstop << endl;

		}
	}
	guah54.GetG2G3(cb3w.g2t, cb3w.g3t);
	if (locdiag) {
		cout << "g2t g3t done" << endl; //return;
		if (op.ton > 1) {
			cb3w.g2t.Print("g2");
			//guah54.DumpA2(); guah54.DumpB2(1);
		}
	}

	//p_cpt2g[16] += cb3.g2t.Count(); p_cpt2g[18] += cb3.g3t.Count();
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& b3 = genb12.bands3[ib3];
		b3.Go(cb3w);
	}

}
void G17B::GoCallB3_12(CALLBAND3& cb3w) {
	p_cpt2g[7]++;

	//if (op.f4 == p_cpt2g[4])cout << Char54out(myb12) << "  [7] "
	//	<< p_cpt2g[7] << endl;
	int locdiag = 0;

	if (op.known) {
		if (knownt > 12) return;
		if (knownt == 12) locdiag = 1;
	}	



	if (locdiag) {
		cout << "b3_12 [7]" << p_cpt2g[7] << endl;
		cout << Char54out(bf_cl9) << " clues " << tc_10_12[0]
			<< "  " << tc_10_12[1] << "  " << tc_10_12[2] << endl;
		locdiag = 1; 
	}
	if (p_cpt2g[7] == op.f7) {
		cout << Char54out(myb12) << "GoCallB3_12 in diag [7]  " << p_cpt2g[7] << endl;
		locdiag = 1;
		if (op.known) {
			cout << "genb12.bands3[0].nbgm " << genb12.bands3[0].nbgm
				<< " genb12.bands3[0].nbgmm " << genb12.bands3[0].nbgmm << endl;
			//genb12.bands3[0].DumpTgmm();
		}
	}
	if (op.f7 && p_cpt2g[7] > op.f7) { aigstop = 1; return; }

	guah54.GetG2G3_12(cb3w.g2t, cb3w.g3t);
	if (locdiag) {
		cb3w.g2t.Print("g2");
	}
	//return;
	if (locdiag) {
		cout << "g2t g3t done" << endl; //return;
		if (op.ton > 1) {
			cb3w.g2t.Print("g2");
			//guah54.DumpA2(); guah54.DumpB2(1);
		}
	}
	{
		register uint64_t F = g17b.myb12 & ~g17b.bf_cl6; // cells 6 to x
		nclues6p = 0;
		register int cell;// build table of cells 
		while (bitscanforward64(cell, F)) {
			F ^= (uint64_t)1 << cell;
			tclues6p[nclues6p++] = cell;
		}
	}
	//p_cpt2g[16] += cb3.g2t.Count(); p_cpt2g[18] += cb3.g3t.Count();
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& b3 = genb12.bands3[ib3];
		b3.Go(cb3w);
	}
	if (op.known && knownt == 12) knownt = 13;
}

void STD_B3::Go(CALLBAND3& cb3e) {
	if (VTEST && p_cpt2g[10] > VTEST) return;
	int locdiag = 0;
	if (op.known&& g17b.knownt == 12) locdiag = 1;
	
	if (p_cpt2g[7] == op.f7) {
		cout<<band << "Go(CALLBAND3& cb3e) [7]  " << p_cpt2g[7] << endl;
		locdiag = 1;
		cout << "xq Dump2 entry go (end of previous)" << endl;
		xq.Dump2();
	}


	p_cpt2g[8]++;
	if (g17b.clean_valid_done == 2 || g17b.aigstop) return;
	memcpy(&g17b.grid0[54], band0, sizeof band0);// used in brute force
	g17b.myband3 = this;
	// note : this is a critical code 
	register uint32_t  Mg2 = 0, Mg3 = 0;
	{
		register uint32_t r;
		register uint64_t V = cb3e.g2t.bf.u64[0] & g.gsocket2.bf.u64[0];
		while (bitscanforward64(r, V)) {
			V ^= (uint64_t)1 << r;//clear bit
			Mg2 |= g.ua2bit27[r];
		}
		V = cb3e.g2t.bf.u64[1] & g.gsocket2.bf.u64[1];
		while (bitscanforward64(r, V)) {
			V ^= (uint64_t)1 << r;//clear bit
			Mg2 |= g.ua2bit27[r + 64];
		}
		V = cb3e.g3t.bf.u64[0] & g.gsocket3.bf.u64[0];
		while (bitscanforward64(r, V)) {
			V ^= (uint64_t)1 << r;//clear bit
			Mg3 |= 1 << g.ua3_imini[r];
		}
		V = cb3e.g3t.bf.u64[1] & g.gsocket3.bf.u64[1];
		while (bitscanforward64(r, V)) {
			V ^= (uint64_t)1 << r;//clear bit
			Mg3 |= 1 << g.ua3_imini[r + 64];
		}
	}
	scritb3.Init(cb3e.ncl, cb3e.cbs);
	g17b.ncluesb3 = scritb3.nb3;
	{
		register uint32_t mask = 7, imini = 0, vmini, bitmini = 1;
		for (; imini < 9; imini++, bitmini <<= 1, mask <<= 3) {
			if (!(vmini = Mg2 & mask)) {
				if (Mg3 & bitmini) {
					scritb3.minix[0] |= bitmini;// mini triplet
					scritb3.critbf |= mask;
					scritb3.mincount++;
				}
			}
			else {
				scritb3.pairs27 |= vmini;
				uint32_t cc = _popcnt32(vmini);
				scritb3.minix[cc] |= bitmini;
				if (cc > 1)scritb3.critbf |= mask;
				else scritb3.critbf |= (mask ^ vmini);
				scritb3.mincount++;
				if (cc == 3) scritb3.mincount++;
				if (cc == 2) {// setup data for later
					register int i27, bit;
					scritb3.t1a |= (mask ^ vmini);// common cell
					bitscanforward(i27, vmini);
					bit = 1 << i27;
					scritb3.t2a[scritb3.nt2a++] = mask ^ bit;
					bit ^= vmini;// other i27
					scritb3.t2a[scritb3.nt2a++] = mask ^ bit;

				}
			}
		}
	}
	if (locdiag) {
		cout << " before  [9] mincount= " << scritb3.mincount
			<< " [10]" << p_cpt2g[10] << endl;

	}
	if ((int)scritb3.mincount > scritb3.nb3) {
		if (locdiag) {
			scritb3.Status("exit a");
			(cb3e.g2t & g.gsocket2).Print("g2");
			(cb3e.g3t & g.gsocket3).Print("g3");
			Dumpg2();
		}
		return;
	}
	p_cpt2g[9]++;
	//if (locdiag) cout << "   [9] " << p_cpt2g[9] << endl;

	Go2();
	if (Go3(cb3e)) return;
	//if (p_cpt2g[10] == VTEST) cout << "  VTEST après go3 "  << endl;


	//cout << " after g03 [9]" << p_cpt2g[9] << " [10]" << p_cpt2g[10] << endl;
	if (xq.nmiss == 1 && xq.nout)	if (xq.Miss1ToMiss0()) return;
	p_cpt2g[10]++;
	if (locdiag) cout << "   [10] " << p_cpt2g[10] << endl;
	if (p_cpt2g[10] == VTEST	||g17b.knownt == 12){
		cout << Char54out(g17b.myb12) << "b12 in test [10]"<< p_cpt2g[10] << endl;
		cout << band << "band3 in test" << endl;
	}
	GoAfter10(cb3e);
}
void STD_B3::Checkkown4() {

	cout << "check tgm ngm=" << nbgm << endl;
	for (uint32_t i = 0; i <= nbbgm; i++)
		tgm64[i].Dump(i << 6);

}

int XQ::Min1_4Disjoint() {
	if (nout < 4)return 0;
	register uint32_t r1 = tout[0]|tout[1], r2=tout[0] & tout[1],r3,r4,x;
	x = tout[2]; r3 = r2 & x; r2 |= r1 & x; r1 |= x;
	x = tout[3]; r4 = r3 & x; r3 |= r2 & x; r2 |= r1 & x; r1 |= x;
	for (uint32_t i = 4; i < nout; i++) {
		x = tout[i];
		r4 |= r3 & x; r3 |= r2 & x; r2 |= r1 & x; r1 |= x;
	}
	if (r4)return 3; else if (r3) return 2; else if (r2) return 1; else return 0;
}
int XQ::Miss1ToMiss0(int debug) {
	if (debug) cout << "Miss1ToMiss0 debugging" << endl;
	if (!nout) {
		cout << "bug miss1 to miss0" << endl;
		g17b.aigstop = 1;		return 1;	}
	{
		register uint32_t A; if (!(A = GetAndout())) return 1;
		t2b[n2b++] = A;
		critbf |= A;
		nout = nmiss = 0;
		SetFilters();
	}
	if (!nin) return 0;
	while (1) {// open the door tomore than one assign
		int as = 0;
		{
			register uint32_t F = fa, A = fb;
			if (debug) {
				cout << Char32out(F) << " F cycle" << endl;
				cout << Char32out(A) << "A cycle" << endl;
			}
			uint32_t n = nin;
			nin = 0;
			for (uint32_t i = 0; i < n; i++) {
				register uint32_t U = tin[i];
				if (F & U) continue;
				if (!(U &= A)) return 1;
				if (debug) {
					cout << Char54out(U) << " U nin="<<nin << endl;
				}
				tin[nin++] = U;
				if (_popcnt32(U) == 1) as = U;
			}			
		}
		if (!as) return 0;
		AddAssigned(as);
	}
}
int XQ::MissxToMiss0(uint32_t ubf) {
	if (!ubf) return 0;
	for (uint32_t i = 0, bit = 1; i < nout; i++, bit <<= 1) {
		register uint32_t u = tout[i];
		if (bit & ubf) {// new crit
			t2b[n2b++] = u;
			critbf |= u;
		}
		else tin[nin++] = u; // put it in to check later
	}

	nout = nmiss = 0;
	SetFilters();
	if (!nin) return 0;
	while (1) {// open the door tomore than one assign
		int as = 0;
		{
			register uint32_t F = fa, A = fb;
			uint32_t n = nin;
			nin = 0;
			for (uint32_t i = 0; i < n; i++) {
				register uint32_t U = tin[i];
				if (F & U) continue;
				if (!(U &= A)) return 1;
				tin[nin++] = U;
				if (_popcnt32(U) == 1) as = U;
			}
		}
		if (!as) return 0;
		AddAssigned(as);
	}
}
int XQ::Miss0CheckTin() {// see if new invalid or new assigned 
	if (nmiss || (!nin)) return 0;
	// check here if fresh potential in "in table"
	while (1) {// open the door tomore than one assign
		int as = 0;
		{
			register uint32_t F = fa, A = fb;
			int n = 0;
			for (uint32_t i = 0; i < nin; i++) {
				register uint32_t U = tin[i];
				if (F & U) continue;
				if (!(U &= A)) return 1;
				tin[n++] = U;
				if (_popcnt32(U) == 1) as = U;
			}
			nin = n;
		}
		if (!as) return 0;
		AddAssigned(as);
	}
}
int XQ::Isoutsize2() {
	if (nout < 2) return 0;
	int bf = 0;
	for (uint32_t i1 = 0; i1 < nout - 1; i1++) {
		register uint32_t U1 = tout[i1];
		for (uint32_t i2 = i1+1; i2 < nout; i2++) {
			register uint32_t U2 = tout[i2];
			if (U1 & U2) continue;
			return (1 << i1) | (1 << i2);
		}
	}
	return 0;
}
int XQ::Isoutsize3() {
	if (nout < 3) return 0;
	int bf = 0;
	for (uint32_t i1 = 0; i1 < nout - 2; i1++) {
		register uint32_t U1 = tout[i1];
		for (uint32_t i2 = i1 + 1; i2 < nout-1; i2++) {
			register uint32_t U2 = tout[i2];
			if (U1 & U2) continue;
			for (uint32_t i3 = i2 + 1; i3 < nout; i3++) {
				register uint32_t U3 = tout[i3];
				if ((U1 | U2) & U3) continue;
				return (1 << i1) | (1 << i2) | (1 << i3);
			}
		}
	}
	return 0;
}
int XQ::Isoutsize4() {
	if (nout < 4) return 0;
	int bf = 0;
	uint32_t limi1 = nout - 3; if (limi1 > 3)limi1 = 3;
	uint32_t limi2 = nout - 2; if (limi2 > 4)limi2 = 4;
	uint32_t limi3 = nout - 1; if (limi3 > 6)limi3 = 6;
	for (uint32_t i1 = 0; i1 < limi1; i1++) {
		register uint32_t U1 = tout[i1];
		for (uint32_t i2 = i1 + 1; i2 < limi2; i2++) {
			register uint32_t U2 = tout[i2];
			if (U1 & U2) continue;
			for (uint32_t i3 = i2 + 1; i3 < limi3; i3++) {
				register uint32_t U3 = tout[i3];
				if ((U1 | U2) & U3) continue;
				register uint32_t uo3 = U1 | U2 |U3;
				for (uint32_t i4 = i3 + 1; i4 < nout ; i4++) {
					register uint32_t U4 = tout[i4];
					if (uo3 & U4) continue;
					return (1 << i1) | (1 << i2) | (1 << i3) | (1 << i4);
				}
			}
		}
	}
	return 0;
}
int XQ::Isoutsize5() {
	if (nout < 5) return 0;
	int bf = 0;
	uint32_t limi1 = nout - 4; if (limi1 > 3)limi1 = 3;
	uint32_t limi2 = nout - 3; if (limi2 > 4)limi2 = 4;
	uint32_t limi3 = nout - 2; if (limi3 > 6)limi3 = 6;
	for (uint32_t i1 = 0; i1 < limi1; i1++) {
		register uint32_t U1 = tout[i1];
		for (uint32_t i2 = i1 + 1; i2 < limi2; i2++) {
			register uint32_t U2 = tout[i2];
			if (U1 & U2) continue;
			for (uint32_t i3 = i2 + 1; i3 < limi3; i3++) {
				register uint32_t U3 = tout[i3];
				if ((U1 | U2) & U3) continue;
				register uint32_t uo3 = U1 | U2 | U3;
				for (uint32_t i4 = i3 + 1; i4 < nout-1; i4++) {
					register uint32_t U4 = tout[i4];
					if (uo3 & U4) continue;
					register uint32_t uo4 = uo3 | U4;
					for (uint32_t i5= i4 + 1; i5 < nout; i5++) {
						register uint32_t U5 = tout[i5];
						if (uo4 & U5) continue;
						//cout << "found size5 " << endl;
						return (1 << i1) | (1 << i2) | (1 << i3) | 
							(1 << i4) | (1 << i5);
					}
				}
			}
		}
	}
	return 0;
}
void XQ::CleanIn() {// usually miss1 last step
	uint32_t nn = nin; nin = 0;
	uint32_t tx[7][100],ntx[7];
	memset(ntx, 0, sizeof ntx);
	for (uint32_t i = 0; i < nn; i++) {
		register uint32_t u = tin[i],nu=~u;
		for (uint32_t j = 0; j < n2a; j++) {
			if (!(t2a[j] & nu)) { u = 0; break; }
		}
		if (!u) continue;
		for (uint32_t j = 0; j < n2b; j++) {
			if (!(t2b[j] & nu)) { u = 0; break; }
			// see also reverse ?? reducing t2b
		}
		if (u) {
			register uint32_t cc = _popcnt32(u);
			if (cc > 6)cc = 6;
			tx[cc][ntx[cc]++] = u;
		}
	}
	nin = 0;
	for (uint32_t i = 0; i < 7; i++) if (ntx[i]) {
		register uint32_t* tu = tx[i];
		for (uint32_t j = 0; j < ntx[i]; j++) {
			register uint32_t U = tu[j], nU = ~U;
			for (uint32_t k = 0; k < nin; k++) {
				if (!(tin[k] & nU)) { U = 0; break; }
			}
			if (U)tin[nin++] = U;
		}
	}
}
void XQ::CleanOut() {// usually miss1 last step
	uint32_t nn = nout; nout = 0;
	uint32_t tx[7][100], ntx[7];
	memset(ntx, 0, sizeof ntx);
	for (uint32_t i = 0; i < nn; i++) {
		register uint32_t u = tout[i],  cc = _popcnt32(u);
		if (cc > 6)cc = 6;
		tx[cc][ntx[cc]++] = u;		
	}
	for (uint32_t i = 0; i < 7; i++) if (ntx[i]) {
		register uint32_t* tu = tx[i];
		for (uint32_t j = 0; j < ntx[i]; j++) {
			register uint32_t U = tu[j], nU = ~U;
			for (uint32_t k = 0; k < nout; k++) {
				if (!(tout[k] & nU)) { U = 0; break; }
			}
			if (U)tout[nout++] = U;
		}
	}
}
void XQ::BuildAllOut() {
	uint32_t a;
	memcpy(tout, t2a, n2a * sizeof a);
	memcpy(&tout[n2a], t2b, n2b * sizeof a);
	memcpy(&tout[n2a+n2b], tin, nin * sizeof a);
	nout = n2a+n2b+nin;
}
void XQ::BuildAllOutMiss0() {
	uint32_t a;
	memcpy(tout,  t2b, n2b * sizeof a);
	memcpy(&tout[ n2b], tin, nin * sizeof a);
	nout =  n2b + nin;
}


inline void STD_B3::Go2(int debug ) {
	if (p_cpt2g[10] == VTEST-1) {
		cout << "test [10] status entry go2" << endl;
		debug = 1;
	}
	xq.Init(scritb3.critbf );
	xq.nb3 = scritb3.nb3;
	xq.nmiss = scritb3.nb3 - scritb3.mincount;
	// Load g2 g3
	{// get back stored status for bf2
		xq.t1a = scritb3.t1a;
		xq.n2a= scritb3.nt2a;
		memcpy(xq.t2a, scritb3.t2a, scritb3.nt2a* sizeof xq.t2a[0]);
	}
	{// get other of mincount
		register uint32_t a;
		if ((a = scritb3.minix[3])) {// now bf3 2 clues
			xq.bf3p=a;// count of minirows 3 pairs in xq
			register uint32_t  i;
			while (bitscanforward(i, a)) {
				a ^= 1 << i;
				i *= 3;// now first of 3 cells
				xq.t2b[xq.n2b++] = g.pat2_27[i++];
				xq.t2b[xq.n2b++] = g.pat2_27[i++];
				xq.t2b[xq.n2b++] = g.pat2_27[i];
			}	
		}
		register uint32_t p27 = scritb3.pairs27;
		if ((a = scritb3.minix[1])) {// now one pair
			register uint32_t  i = 0, bit = 1, mask = 7;
			for (; i < 9; i++, bit <<= 1, mask <<= 3)if (a & bit) {
				register uint32_t i27, mask2 = mask & p27;
				bitscanforward(i27, mask2);
				xq.t2b[xq.n2b++] = g.pat2_27[i27];
			}
		}
		if ((a = scritb3.minix[0])) {// now triplet 
			register uint32_t a = scritb3.minix[0], i = 0, bit = 1, mask = 7;
			for (; i < 9; i++, bit <<= 1, mask <<= 3)if (a & bit) {
				register uint32_t i27, mask2 = mask & p27;
				bitscanforward(i27, mask2);
				xq.t2b[xq.n2b++] = mask;
			}
		}
	}
	if (debug) {
		cout << "\n\ngo2 [9]=" << p_cpt2g[9] << endl;
		scritb3.Status("go2 ");
		xq.Status();
	}

}
int STD_B3::Go3(CALLBAND3& cb3) {
	int debug = 0;
	xq.SetFilters();
	if ( p_cpt2g[10] == VTEST-1) {
		cout << "test [10] status entry go3" << endl;
		xq.Status();
		debug = 2;
	}

	{// add band 3 UA size 4
		register uint32_t F = xq.fa, C = xq.critbf, A = xq.fb; 
		//if (debug)xq.Dump2();
		for (uint32_t i = 0; i < nua; i++) {
			register uint32_t U = tua[i], cc = U >> 27;
			U &= BIT_SET_27;
			if (cc > 4) { xq.iuas4 = i; break; }
			if(debug)cout << Char27out(U) << "u" << endl;
			if (F & U) continue;
			if (!(U &=A)) return 1; // dead
			cc = _popcnt32(U); // remaining cells if critical
			if (cc == 1) {// new assign 
				F |= U;
				xq.t1a |= U;
				A = C = xq.AddAssigned(U);
				if (debug) {
					cout << Char27out(U) << "u new assign" << endl;
					xq.Dump2();
				}
				continue;
			}
			if (U & C) {
				if (debug)cout << "this is infield" << endl;
				xq.Addin(U);
			}
			else {
				if (debug)cout << "this is outfield" << endl;
				xq.Addout(U);
			}
		}
		{// here to insert guas 4/6 

			uint32_t tw[80], ntw = 0;
			register BF128 v4 = cb3.g2t & g.gsocket4;
			if (debug>1 ) {
				v4.Print(" potential 4 clues in g4 ");
				cout << band << endl;
				//guah54.DumpB2(1);
			}
			{
				register uint32_t r;
				register uint64_t V = cb3.g2t.bf.u64[0] & g.gsocket4.bf.u64[0];
				while (bitscanforward64(r, V)) {
					V ^= (uint64_t)1 << r;//clear bit
					if (debug)cout << Char27out(g.pat2[r]) << endl;
					tw[ntw++] = g.pat2[r];
				}
				V = cb3.g2t.bf.u64[1] & g.gsocket4.bf.u64[1];
				while (bitscanforward64(r, V)) {
					V ^= (uint64_t)1 << r;//clear bit
					if (debug)cout << Char27out(g.pat2[r + 64]) << endl;;
					tw[ntw++] = g.pat2[r + 64];
				}
			}
			if (debug > 1)cout << "check ntw=" << ntw << endl;
			for (uint32_t i = 0; i < ntw; i++) {
				register uint32_t U = tw[i];
				if (F & U) continue;
				if (!(U &= A)) return 1; // dead
				register uint32_t cc = _popcnt32(U);
				if (cc == 1) {// new assign
					F |= U;
					xq.t1a |= U;
					A = C = xq.AddAssigned(U);
					if (debug) {
						cout << Char27out(U) << "u new assign" << endl;
						xq.Dump2();
					}
					continue;
				}
				if (U & C) {
					if (debug)cout << "this is infield" << endl;
					xq.Addin(U);
				}
				else {
					if (debug)cout << "this is outfield" << endl;
					xq.Addout(U);
				}

			}
		}
	}
	if (p_cpt2g[10] == VTEST-1) {
		cout << "test [10] status end go3 before Miss0CheckTin" << endl;
		xq.Status();
	}

	return xq.Miss0CheckTin();

}
void STD_B3::GoAfter10(CALLBAND3& cb3) {
	//if (xq.nmiss > 1) return;
	int locdiag = 0;
	if (p_cpt2g[10] == VTEST || g17b.knownt == 12) {
		cout << "test [10] status entry after 10" << endl;
		if(p_cpt2g[10] == VTEST)xq.Status();
		locdiag = 1;
	}

	if (locdiag) {
		cout << " continue after [10] nbgm=" << nbgm << endl;
		//xq.Status();
	}
	xq.BuildCheckRedundant();
	uint32_t tw[1000], ntw = 0;
	{// Add expensive 4 
		register uint32_t  		A = xq.fb, F = xq.fa, U;// if critical all must be in
		for (uint32_t i = 0; i <= nbbgm; i++) {
			GUM64& gw = tgm64[i];
			register uint64_t V = gw.Getv(g17b.tclues6p, g17b.nclues6p);
			register uint32_t r;
			while (bitscanforward64(r, V)) {
				V ^= (uint64_t)1 << r;
				U = gw.tb3[r];
				if (U & F)continue;
				if (!U) {
					cout << "bug tguam 0 [10]" << p_cpt2g[10] << endl;
					g17b.aigstop = 1;
					//DumpTgm();
					return;
				}
				if (!(U &= A)) {
					//if (locdiag)cout << "dead guam4" << endl;
					return;// dead
				}
				else	tw[ntw++] = U;
			}
		}
	}
	//if (ntw > 500) cout << "lim500 ntw=" << ntw << " [10] " << p_cpt2g[10] << endl;
	{ // check and clear redundancy
		if (locdiag)cout << Char27out(xq.critbf) << " end ntw=" << ntw << endl;
		uint32_t nn = 0;
		for (uint32_t i = 0; i < ntw; i++) {
			register uint32_t U = tw[i], nU = ~U;
			//if (locdiag)cout << Char27out(U) << "u to check" << endl;
			// internal redundancy or equal
			for (uint32_t j = 0; j < nn; j++)
				if (!(tw[j] & nU)) { nU = 0; break; }
			if (!nU) {
				//if (locdiag)cout << "internally redundant" << endl; 
				continue;
			}
			for (uint32_t j = 0; j < xq.nred; j++)
				if (!(xq.t2b[j] & nU)) { nU = 0; break; }
			if (!nU) {
				//if (locdiag)cout << Char27out(U) << "externally redundant" << endl; 
				continue;
			}
			tw[nn++] = U;
		}
		if (!(ntw = nn)) goto end10; // nothing more to do
	}
	if (locdiag)cout << " ntw residual count " << ntw << endl;
	if(!xq.nmiss) {// nmiss0 look for one clue
		int aig = 0;
		while (1) {
			uint32_t nn = ntw; ntw = 0;
			register uint32_t A = xq.fb, F = xq.fa;
			for (uint32_t i = 0; i < nn; i++) {
				register uint32_t U = tw[i];
				if (U & F)continue;
				if (!(U &= A)) return;
				register uint32_t cc = _popcnt32(U);
				if (cc > 1)tw[ntw++] = U;
				else {
					aig = 1;
					xq.AddAssigned(U);
					A = xq.fb; F = xq.fa;// reload
				}
			}
			if (aig) {
				if (xq.Miss0CheckTin()) return;
				//cout << " continue after [10] miss0 with one" << p_cpt2g[10] << endl;
				//xq.Status();
				aig = 0;
			}
			else break;
		}
		memcpy(&xq.tin[xq.nin], tw, ntw * sizeof xq.tin[0]);
		xq.nin += ntw;
		goto end10;
	}
	//___________ look for out disjoints leading to miss0
	{ //split in/out 
		register uint32_t C = xq.critbf;
		for (uint32_t i = 0; i < ntw; i++) {
			register uint32_t U = tw[i];
			if (U & C)xq.Addin(U);
			else  xq.Addout(U);
		}
	}
	//if (xq.nin > 200) cout << "lim200 nin=" << xq.nin << " [10] " << p_cpt2g[10] << endl;
	//if (xq.nout > 200) cout << "lim200 nout=" << xq.nout << " [10] " << p_cpt2g[10] << endl;
	if (xq.nout < xq.nmiss) goto end10;
	if (xq.nout == xq.nmiss && xq.NoDisjoint()) goto end10;
	{ // look for disjoints of required size 
		switch (xq.nmiss) {
		case 1:if (xq.Miss1ToMiss0()) 	return;
			goto end10;
		case 2:
			if(xq.MissxToMiss0(xq.Isoutsize2())) return;
			goto end10;
		case 3:
			if (xq.MissxToMiss0(xq.Isoutsize3()))return;
			goto end10;
		case 4:
			if (xq.MissxToMiss0(xq.Isoutsize4()))return;
			goto end10;
		}
	}
	end10:
	p_cpt2g[11]++;
	if (locdiag) {
		cout << " continue after [11] xq.nmiss= "<< xq.nmiss << endl;
		if (p_cpt2g[10] == VTEST)xq.Status();
	}
	if(xq.nmiss)	GoAfter11(cb3);
	else GoAfter11Miss0(cb3);
}
void STD_B3::GoAfter11Miss0(CALLBAND3& cb3) {// add now size 6 and more 
	int locdiag = 0;
	if (p_cpt2g[10] == VTEST) {
		cout << "test  status entry after 11 xq.iuas4="<< xq.iuas4 << endl;
		locdiag = 1;
		xq.Status();
	}
	p_cpt2g[12]++;
	// add band 3 UA size > 4
	{
		register uint32_t F = xq.fa, A = xq.fb;
		for (uint32_t i = xq.iuas4; i < nua; i++) {
			register uint32_t U = tua[i], cc = U >> 27;
			U &= BIT_SET_27;
			if (F & U) continue;
			if (!(U &= A)) return; // dead
			cc = _popcnt32(U); // remaining cells if critical
			if (cc == 1) {// new assign 
				F |= U;
				xq.t1a |= U;
				A  = xq.AddAssigned(U);
				continue;
			}
			xq.Addin(U);
		}
		if (locdiag) xq.Status();
		{//  insert guas 6 
			uint32_t tw[80], ntw = 0;
			register BF128 v4 = cb3.g2t & g.gsocket4;
			{
				register uint32_t r;
				register uint64_t V = cb3.g2t.bf.u64[0] & g.gsocket6.bf.u64[0];
				while (bitscanforward64(r, V)) {
					V ^= (uint64_t)1 << r;//clear bit
					tw[ntw++] = g.pat2[r];
				}
				V = cb3.g2t.bf.u64[1] & g.gsocket6.bf.u64[1];
				while (bitscanforward64(r, V)) {
					V ^= (uint64_t)1 << r;//clear bit
					tw[ntw++] = g.pat2[r + 64];
				}
			}
			for (uint32_t i = 0; i < ntw; i++) {
				register uint32_t U = tw[i];
				if (F & U) continue;
				if (!(U &= A)) return; // dead
				register uint32_t cc = _popcnt32(U);
				if (cc == 1) {// new assign
					F |= U;
					xq.t1a |= U;
					A =  xq.AddAssigned(U);
					continue;
				}
				xq.Addin(U);
			}
		}
	}
	if (locdiag) {
		cout << "miss0 bbb nbbgmm="<< nbbgmm << endl; xq.Status(); //return;
	}
	// add now guamm from the band (6 clues and more) and add t2b tin
	uint32_t tw[500], ntw = 0;
	{
		register uint32_t  	A = xq.fb, F = xq.fa, U;// miss0 all must be in
		if (locdiag) {
			//cout << Char27out(F) << "F27" << endl;
			//cout << Char32out(F) << "F" << endl;
			//cout << Char32out(A) << "A" << endl;
		}	
		for (uint32_t i = 0; i <= nbbgmm; i++) {
			GUM64& gw = tgm64m[i];
			register uint64_t V = gw.Getv(g17b.tclues6p, g17b.nclues6p);
			if (locdiag) {
				cout << Char64out(V) << " v  imm=" << i << endl;
				//gw.Dump(0);
			}

			register uint32_t r;
			while (bitscanforward64(r, V)) {
				V ^= (uint64_t)1 << r;
				U = gw.tb3[r];
				if (U & F)continue;
				if (!U) {
					cout << "bug tguamm 0 [10]" << p_cpt2g[10] << endl;
					g17b.aigstop = 1;
					return;
				}
				if (locdiag)cout << Char32out(U) << " uamm" << endl;
				if (!(U &= A)) return;// dead				
				{
					register uint32_t cc = _popcnt32(U);
					if (cc == 1) {// new assign
						F |= U;
						xq.t1a |= U;
						A = xq.AddAssigned(U);
					}
					else tw[ntw++] = U;
				}

			}
		}
		if (locdiag) {
			cout << "ntwa=" << ntw << " xq.n2b=" << xq.n2b << endl;
			//return;
			//cout << Char27out(F) << "F27" << endl;
			//cout << Char32out(F) << "F" << endl;
			//cout << Char32out(F) << "A" << endl;
		}
		//if (locdiag) return;

		for (uint32_t i = 0; i < xq.n2b; i++) {
			//if (locdiag)continue;
			U = xq.t2b[i];
			if (U & F)continue;
			if (!(U &= A))return;
			if (_popcnt32(U) == 1) {// new assign
				F |= U;
				xq.t1a |= U;
				A = xq.AddAssigned(U);
				if (locdiag) cout<< Char27out(U) << " new assign"  << endl;
			}
			else tw[ntw++] = U;
		}
		if (locdiag) cout << "ntwb=" << ntw << endl;
		//if (locdiag) return;

		// add in to clean all single
		if(1)for (uint32_t i = 0; i < xq.nin; i++) {
			U = xq.tin[i];
			if (U & F)continue;
			if (!(U &= A))return;
			if (_popcnt32(U) == 1) {// new assign
				F |= U;
				xq.t1a |= U;
				A = xq.AddAssigned(U);
			}
			else tw[ntw++] = U;
		}
		xq.nin = 0;
	}
	// clean size 1 in tw
	while (1) {
		register uint32_t  	A = xq.fb, F = xq.fa, U;// miss0 all must be in
		uint32_t bf = 0, nn = ntw; ntw = 0;
		for (uint32_t i = 0; i < nn; i++) {
			U = tw[i];
			if (U & F)continue;
			if (!(U &= A))return;// new dead
			if(_popcnt32(U)==1)bf=U ;
			tw[ntw++] = U;
		}
		if (bf)  xq.AddAssigned(bf);
		else break;
	}

	// last check and expand builder  this is miss0, t2a is covered

	// sort by size
	{ 
		register uint32_t  	A = xq.fb, F = xq.fa, U,cc;// last status
		uint32_t tws[9][100], ntws[9];// size 1 to 6 plus margin
		memset(ntws, 0, sizeof ntws);
		for (uint32_t i = 0; i < ntw; i++) {
			U = tw[i];
			if (U & F)continue;
			if (!(U &= A))return;
			cc = _popcnt32(U); if (cc > 6)cc = 6;
			tws[cc][ntws[cc]++] = U;
		}
		if (locdiag) {
			cout << "sort n 1-6 ";
			for (int i = 1; i < 7; i++) cout << ntws[i] << " ";
			cout << endl;
		}
		//if (locdiag) return;


		xq.nout = 0;
		for (int i = 1; i < 7; i++) if (ntws[i]) {
			uint32_t* t = tws[i];
			for (uint32_t j = 0; j < ntws[i]; j++) {
				U = t[j];
				register uint32_t  nU = ~U;
				for (uint32_t k = 0; k < xq.nout; k++) {
					if (!(xq.tout[k] & nU)) { U = 0; break; }
				}
				if (U)xq.tout[xq.nout++] = U;
			}
		}
	}
	g17b.GoEndMiss0();
}


void STD_B3::GoAfter11(CALLBAND3& cb3) {// add now size 6 and more 
	p_cpt2g[13]++;
	int locdiag = 0;
	if (p_cpt2g[10] == VTEST || g17b.knownt == 12) {
		cout << "test [10]  entry GoAfter11" << endl;
		if (p_cpt2g[10] == VTEST)xq.Status();
		locdiag = 1;
	}
	// Add first some 4 clues disjoints
	{ // look for disjoints of required size 
		register uint32_t bf = 0;
		switch (xq.nmiss) {
		default: if (xq.nout >= 4)
			if ((bf = xq.Isoutsize4())) goto endaok;
		case 4:if (xq.nout >= 3)
			if ((bf = xq.Isoutsize3())) goto endaok;
		case 3:if (xq.nout >= 2)
			if ((bf = xq.Isoutsize2())) goto endaok;
		case 2:if (xq.nout) { bf = 1; goto endaok; }
		case 1:goto endb;// do nothing if not more than one 
		}
	endaok: {// insert the disjoints in t2b 
		uint32_t nn = xq.nout, nn2; xq.nout = 0;
		bitscanreverse(nn2, bf);// last bit in bf

		for (uint32_t i = 0, bit = 1; i <= nn2; i++, bit <<= 1) {
			register uint32_t u = xq.tout[i];
			if (bit & bf) {// new crit
				xq.t2b[xq.n2b++] = u;
				xq.critbf |= u;
			}
		}
		register uint32_t A = xq.critbf;
		for (uint32_t i = 0, bit = 1; i < nn; i++, bit <<= 1) {
			register uint32_t u = xq.tout[i];
			if (!(bit & bf)) {// not a disjoint
				if (u & A)xq.tin[xq.nin++] = u;
				else xq.tout[xq.nout++] = u; // put it in to check later
			}
		}
		xq.nmiss -= _popcnt32(bf);
		}
endb:	;
	}
	// -----------------------  add uas size>4
	{// add band 3 UA size 4
		register uint32_t F = xq.fa, C = xq.critbf;
		for (uint32_t i = xq.iuas4; i < nua; i++) {
			register uint32_t U = tua[i], cc = U >> 27;
			U &= BIT_SET_27;
			if (F & U) continue;
			cc = _popcnt32(U); // remaining cells if critical
			if (U & C) 	xq.Addin(U);
			else	xq.Addout(U);
		}
		uint32_t tw[300], ntw = 0;
		{// here to insert guas 4/6 
			register BF128 v4 = cb3.g2t & g.gsocket4;
			register uint32_t r;
			register uint64_t V = cb3.g2t.bf.u64[0] & g.gsocket4.bf.u64[0];
			while (bitscanforward64(r, V)) {
				V ^= (uint64_t)1 << r;//clear bit
				tw[ntw++] = g.pat2[r];
			}
			V = cb3.g2t.bf.u64[1] & g.gsocket4.bf.u64[1];
			while (bitscanforward64(r, V)) {
				V ^= (uint64_t)1 << r;//clear bit
				tw[ntw++] = g.pat2[r + 64];
			}
		}
		for (uint32_t i = 0; i < ntw; i++) {
			register uint32_t U = tw[i];
			//if (F & U) continue;// can not be here
			if (U & C)	xq.Addin(U);
			else xq.Addout(U);
		}
		if (locdiag && op.ton) {
			cout << "  GoAfter11 bbbb  " << endl;
			xq.Status();
			//cout << Char27out(F) << "F" << endl;
			//cout << Char27out(C) << "C" << endl;
			//DumpTgmm();
		}

		{
			ntw = 0;
			for (uint32_t i = 0; i <= nbbgmm; i++) {
				GUM64& gw = tgm64m[i];
				register uint64_t V = gw.Getv(g17b.tclues6p, g17b.nclues6p);
				//if (locdiag) {
					//cout << Char64out(V) << " ibgmm=" << i << endl;
				//}
				register uint32_t r;
				while (bitscanforward64(r, V)) {
					V ^= (uint64_t)1 << r;
					register uint32_t U = gw.tb3[r],nU=~U;
					//if (U & F)continue;
					for (uint32_t j = 0; j < ntw; j++) {// clean redundancy
						if (!(tw[j] & nU)) {U = 0; break;	}						
					}
					if (U) tw[ntw++] = U;
					//if (U && g17b.knownt == 11) {
						//cout << Char27out(U) << " index=" << r + 64 * i << endl;
					//}
				}
			}
		}
		for (uint32_t i = 0; i < ntw; i++) {
			register uint32_t U = tw[i];
			if (F & U) continue; 
			if (U & C)	xq.Addin(U);
			else xq.Addout(U);
		}
	}
	if (xq.nin > 200) cout << "lim200 nin=" << xq.nin << " [10] " << p_cpt2g[10] << endl;
	if (xq.nout > 200) cout << "lim200 nout=" << xq.nout << " [10] " << p_cpt2g[10] << endl;


	// try to push to lower "miss" using out disjoints
	if (xq.nout > xq.nmiss || (xq.nout == xq.nmiss && (!xq.NoDisjoint()))) {
		switch (xq.nmiss) {
		case 0: if (xq.nout) return; break;;
		case 1:if (xq.Miss1ToMiss0(locdiag)) 	return;
			break;
		case 2:
			if (xq.MissxToMiss0(xq.Isoutsize2())) return;
			break;
		case 3:
			if (xq.MissxToMiss0(xq.Isoutsize3()))return;
			break;
		}
	}
	// do the best to clean xq.out
	for (uint32_t i = 0; i < xq.nout; i++) {// good chance to get disjoints 
		register uint32_t U = xq.tout[i];
		if (!(xq.critbf & U)) {
			if (!xq.nmiss) return;
			xq.t2b[xq.n2b++] = U;
			xq.nmiss--;
			xq.critbf |= U;
		}
		else xq.Addin(U);
	}
	xq.nout = 0;
	// final status, 
	p_cpt2g[14]++;
	if (locdiag) {
		cout << "test [10] status  after 14" << endl;
		if (p_cpt2g[10] == VTEST)xq.Status();
	}
	p_cpt2g[60 + xq.nmiss]++;

	if ((!xq.nmiss)) {
		//if (!(p_cpt2g[60] & 1023))) xq.Status();
		g17b.BridgeEndMiss0();
	} 
	else  g17b.GoNotMiss0();
}

void G17B::GoNotMiss0() {
	p_cpt2g[15]++;
	if (VerifyValid()) return;
	if (xq.nout) {
		cout << "\n\nGoNotMiss0() entry with xqnout, see why" << endl;
		xq.Status();
		aigstop = 1;
		return;
	}
	int locdiag = 0;
	if (p_cpt2g[10] == VTEST || g17b.knownt == 12)  locdiag = 1;
	//___ test global till valid critbf
	{
		while (1) {
			nt3more = 0;
			if (locdiag)cout << Char27out(xq.critbf) << "chek valid b3" << endl;
			if (IsValidB3(xq.critbf)) {// not valid, new outfield
				 if (locdiag)for (uint32_t i = 0; i < nt3more; i++)
					cout << Char27out(t3more[i]) << " fresh ua i=" << i << endl;
				if (xq.nmiss == 1) {// now dummy ua and miss 0
					xq.Addout(anduab3);
					if (xq.Miss1ToMiss0()) return;
					if (locdiag) {
						cout << "after xq.Miss1ToMiss0() " << endl;
						xq.Status();
					}
					BridgeEndMiss0();
					return;
				}
				// missing here possibility to have 2 gu2 same mini row
				if (nt3more == 1) {// most often
					xq.t2b[xq.n2b++] = anduab3;
					xq.nmiss--;
					xq.critbf|= anduab3;
					continue;
				}
				for (uint32_t i = 0; i < nt3more; i++) {// good chance to get disjoints 
					register uint32_t U= t3more[i];
					if (!(xq.critbf & U)) {
						if (!xq.nmiss) return;
						xq.t2b[xq.n2b++] = U;
						xq.nmiss--;
						xq.critbf |= U;
					}
					else xq.Addin(U);
				}
			}
			else break;
		}
	}
	//if (locdiag) { cout << "exit xhile" << endl; xq.Status(); }
	if (!xq.nmiss) { BridgeEndMiss0(); return; }
	// now in field always possible any outlied can be added 
	p_cpt2g[16]++;
	//if (p_cpt2g[10] == 3510150) locdiag = 1;
	xq.CleanIn();// hit by critical or subset/equal
	if (locdiag) { cout << "after clean in" << endl; 
		if (p_cpt2g[10] == VTEST)xq.Status();
	}

	if (xq.nmiss > 1) {// direct to expand 
		xq.BuildAllOut();
		//if (locdiag) { cout << "after Buildout" << endl; xq.Status(); }
		GoEndAll(0, BIT_SET_27);
		return;
	}
	p_cpt2g[17]++;
	// now miss1 no out try all in  
	XQ xqr = xq;
	TryMiss1Subcritical();
	if (locdiag)  cout << "back subcritical" << endl;
	// then add out field as dummy ua leading to miss 2
	xq = xqr;
	xq.t2b[xq.n2b++] = BIT_SET_27 ^ xq.critbf;
	xq.critbf= BIT_SET_27;
	xq.nmiss = 0;
	xq.SetFilters();
	if (xq.Miss0CheckTin()) return;
	xq.BuildAllOutMiss0();
	xq.CleanOut();
	if (locdiag) {
		cout << " [17]after cleanout" << endl;
		if (p_cpt2g[10] == VTEST)xq.Status();
	}
	int ntoass = xq.nb3 - _popcnt32(xq.t1a);
	xq.nin = 0;
	switch (ntoass) {
	case 0:if (xq.nout) return;
		if (!IsValidB3(xq.t1a)) Out17(xq.t1a);
		return;
	case 1:if (xq.Isoutsize2()) return;		break;
	case 2:if (xq.Isoutsize3()) return;		break;
	case 3:if (xq.Isoutsize4()) return;		break;
	case 4:if (xq.Isoutsize5()) return;		break;
	}
	if (locdiag) { cout << " [18]do miss1" << endl;
		if (p_cpt2g[10] == VTEST)xq.Status();
	}
	p_cpt2g[18]++;
	GoEndAll(xq.fa, xq.fb);

}
void G17B::TryMiss1Subcritical() {
	int locdiag = 0;
	if (p_cpt2g[10] == VTEST || g17b.knownt == 12)  locdiag = 1;
	xq.BuildAllOut(); xq.nin = 0;
	if (locdiag) { 
		cout << "Subcritical() after buildout" << endl; 
		if (p_cpt2g[10] == VTEST)xq.Status();
	}
	register uint32_t F = 0, A = xq.critbf;

	{ // could see assigned
		while (1) {
			uint32_t nn = xq.nout, ass = 0; xq.nout = 0;
			for (uint32_t i = 0; i < nn; i++) {
				register uint32_t u = xq.tout[i];
				if (F & u) continue;
				if (!(u &= A)) return;// dead
				xq.tout[xq.nout++] = u;
				if (_popcnt32(u) == 1)ass |= u;
			}
			if (!ass) break;
			F |= ass;xq.fa |= ass;
			A &= ~ass;xq.fb &= ~ass;
		}
		p_cpt2g[89]++;
		if (locdiag) { 
			cout << "after [89]" << endl;
			if (p_cpt2g[10] == VTEST) {
				cout << Char27out(F) << "F" << endl;
				cout << Char27out(A) << "A" << endl;
				xq.Status();
			}
		}
		// look for clues "more"
		int nmiss = xq.nmiss, ff = F;
		if (F) {// can be out of the min count
			for (uint32_t i = 0; i < xq.n2a; i += 2) {
				register uint32_t a = xq.t2a[i], b = xq.t2a[i + 1],
					mask = a | b,f=mask&F;
				if (f) {
					ff &= ~(mask & F);
					if (f & ~(a & b)) {
						nmiss--;
						if (_popcnt32(f) > 2)return;
					}
				}
			}
			if (ff) { // hit in t2b where all are disjoints
				// mini rows with 3 pairs
				int n2b0 = 0;
				if (xq.bf3p) {
					for (int i = 0, bit = 1,mask=7; i < 9; i++, bit <<= 1,mask<<=3)
						if (bit & xq.bf3p) {
							n2b0 += 3; // skip 3 pairs
							register uint32_t a = mask & ff;
							if (_popcnt32(a) > 2)	nmiss--;
						}
				}
				// others are disjoints
				for (uint32_t i = n2b0; i < xq.n2b; i++) {
					register uint32_t a = xq.t2b[i] & ff;
					if (a) 	nmiss -= (_popcnt32(a) - 1);
				}
			}
		}
		if (nmiss < 0) return;
		xq.nmiss = nmiss;
		if (locdiag) { 
			cout << "after [90] aaa nmiss="<<nmiss << endl;
			if (p_cpt2g[10] == VTEST) {
				cout << Char27out(F) << "F" << endl;
				cout << Char27out(A) << "A" << endl;
				xq.Status();
			}
		}

		if (!nmiss) {
			// rebuild F and A in critical mode 
			for (uint32_t i = 0; i < xq.n2a; i++) 
				 A &= ~xq.t2a[i];// all are hit
			int n2b0 = 0;
			if (xq.bf3p) {
				for (int i = 0, bit = 1, mask = 7; i < 9; i++, bit <<= 1, mask <<= 3)
					if (bit & xq.bf3p) {
						n2b0 += 3; // skip 3 pairs
						register uint32_t a = mask & ff;
						if (_popcnt32(a) > 1)	A &= ~mask;
					}
			}
			// other are disjoints
			for (uint32_t i = n2b0; i < xq.n2b; i++) {
				register uint32_t a = xq.t2b[i];
				if (a & F) A &= ~a;
			}
			if (locdiag) { 
				cout << "after [90] bbb" << endl; 
				if (p_cpt2g[10] == VTEST) {
					cout << Char27out(F) << "F" << endl;
					cout << Char27out(A) << "A" << endl;
					xq.Status();
				}
			}
			// new attempt to assign with the updated status
			while (1) {
				uint32_t nn = xq.nout, ass = 0; xq.nout = 0;
				for (uint32_t i = 0; i < nn; i++) {
					register uint32_t u = xq.tout[i];
					if (F & u) continue;
					if (!(u &= A)) return;// dead
					xq.tout[xq.nout++] = u;
					if (_popcnt32(u) == 1)ass = u;
				}
				if (!ass) break;// ass must be one ua in t2b
				F |= ass;
				A &= ~ass;
			}
			if (locdiag) {
				cout << "after [90] cccc" << endl;
				if (p_cpt2g[10] == VTEST) {
					cout << Char27out(F) << "F" << endl;
					cout << Char27out(A) << "A" << endl;
					xq.Status();
				}
			}



			uint32_t nass = _popcnt32(F);
			if (nass > xq.nb3) return; // should not be
			if (nass == xq.nb3) {// no expand
				if (xq.nout) return;
				//if (VerifyValid()) return;
				if(!IsValidB3(F)) Out17(F);
				return;
			}

			if (0) {
				cout << Char27out(F) << "expected miss1 in  end nmiss= 0" << endl;
				cout << Char27out(A) << "active" << endl;
				cout << Char27out(F) << "in miss0 mode0 "  << endl;
				cout << Char27out(A) << "active" << endl;
				xq.Status();
				cout << Char27out(F) << "final assign  " << endl;
			}
			//xq.DumpOut();
			GoEndAll(F, A,locdiag);
			return;
		}			
	}
	// here nmiss=1 if no out take A
	if (!xq.nout)xq.tout[xq.nout++] = A;
	GoEndAll(F, A,locdiag);// no "more clue" found
}

void G17B::BridgeEndMiss0() {// build tout mis0 exit  GoAfter11() 
	int locdiag = 0;
	if (p_cpt2g[10] == VTEST) {
		cout << "entry BridgeEndMiss0() in test" << endl;
		locdiag = 1;
		xq.Status();
	}
	// sort by size
	{
		register uint32_t  	A = xq.fb, F = xq.fa, U, cc;// last status
		uint32_t tws[8][100], ntws[8];// size 1 to 6
		memset(ntws, 0, sizeof ntws);
		for (uint32_t i = 0; i < xq.n2b; i++) {
			U = xq.t2b[i];
			if (U & F)continue;
			if (!(U &= A))return;
			cc = _popcnt32(U); if (cc > 6)cc = 7;
			tws[cc][ntws[cc]++] = U;
		}
		for (uint32_t i = 0; i < xq.nin; i++) {
			U = xq.tin[i];
			if (U & F)continue;
			if (!(U &= A))return;
			cc = _popcnt32(U); if (cc > 6)cc = 7;
			tws[cc][ntws[cc]++] = U;
		}

		for (uint32_t i = 0; i < xq.nout; i++) {
			U = xq.tout[i];
			if (U & F)continue;
			if (!(U &= A))return;
			cc = _popcnt32(U); if (cc > 6)cc = 7;
			tws[cc][ntws[cc]++] = U;
		}
		xq.nout = xq.nin = 0;
		for (int i = 1; i < 8; i++) if (ntws[i]) {
			uint32_t* t = tws[i];
			for (uint32_t j = 0; j < ntws[i]; j++) {
				U = t[j];
				register uint32_t  nU = ~U;
				for (uint32_t k = 0; k < xq.nout; k++) {
					if (!(xq.tout[k] & nU)) { U = 0; break; }
				}
				if (U)xq.tout[xq.nout++] = U;
			}
		}
	}

	if (locdiag) {
		cout << "call GoEndMiss0() in test" << endl;
		xq.Status();
	}
	GoEndMiss0();
}
void G17B::GoEndMiss0() {
	int locdiag = 0;
	if (p_cpt2g[10] == VTEST) {
		cout << "GoEndMiss0() in test" << endl;
		locdiag = 1;
		//return;
	}	
	// if one in "tout" is subset of tb, replace it
	{	
		int aig = 0;
		for (uint32_t i = 0; i < xq.nout; i++) {
			register uint32_t U = xq.tout[i];
			for (uint32_t j = 0; j < xq.n2b; j++) {
				register uint32_t U2 = xq.t2b[j];
				if (U == U2)  break;
				if (!(U & ~U2)) {// subset replace U2
					aig = 1;
					xq.t2b[j] = U;
				}
			}
		}
		if (aig) {// could have more to clear
			p_cpt2g[38]++;
			xq.SetFreshCrit();
			while (1) {
				register uint32_t A = xq.fb, F = xq.t1a;
				uint32_t bf = 0, nn = xq.nout; xq.nout = 0;
				for (uint32_t i = 0; i < nn; i++) {
					register uint32_t U = xq.tout[i];
					if (U & F) continue;
					if (!(U &= A)) return;
					if (_popcnt32(U) == 1)bf = U;
					xq.tout[xq.nout++] = U;
				}
				if (bf) xq.AddAssigned(bf);
				else break;
			}
			p_cpt2g[39]++;
		}
	}
	p_cpt2g[40]++;
	nt3more = 0;
	int ntoass = xq.nb3 - _popcnt32(xq.t1a);
	if (locdiag) {
		cout << "GoEndMiss0() in test [40] ntoass="<< ntoass << endl;
		xq.Status();
	}
	if (ntoass < 0) {
		cout << "bug miss0 ntoass<0  [8]" << p_cpt2g[8]
			<< " na,b " << scritb3.nb3 << " " << _popcnt32(xq.t1a) << endl;
		xq.Status();
		aigstop = 1; return;
	}
	if (VerifyValid())return;
	if (!ntoass) {
		if(!IsValidB3(xq.t1a))Out17(xq.t1a);
		return;
	}
	p_cpt2g[45]++;
	if (IsValidB3(xq.fb | xq.t1a)) return;// possible 18??
	p_cpt2g[46]++; // around 6000 in sample test
	GoEndAll(xq.t1a, xq.fb,locdiag);
	return;
	if (p_cpt2g[7] <= op.f7) {	}
	if (knownt == 11) 		scritb3.Status("bbb");
}
void  G17B::GoEndAll(uint32_t bf, uint32_t ac, int debug) {//  call the relevant expand b3
	//if (debug) return;
	p_cpt2g[70]++;
	p_cpt2g[72 + (xq.nout >> 5)]++;
	ac &= BIT_SET_27;// be sure to have no extra digit
	/*
	if (xq.nout > 80) {
		cout << " high count uas to expand [10]" << p_cpt2g[10]
			<< endl;
		cout << Char27out(bf) << "assigned" << endl;
		cout << Char27out(ac) << "active " << endl;
		xq.Status();
	}
	*/
	if (debug) {
		cout << " GoEndAll debug [70]" << p_cpt2g[70] << endl;
		cout << Char27out(bf) << "assigned" << endl;
		cout << Char27out(ac) << "active " << endl;
		xq.DumpOut();
		//return;
	}

	int ass = _popcnt32(bf);
	ntoassb3 = xq.nb3-ass;
	int locdiag = 0;
	if (p_cpt2g[70] == VTEST || debug) {
		cout << " entry GoEndAll ntoassb3=" << ntoassb3 << " knownt=" << knownt << endl;
		xq.Status();		locdiag = 1;
	}
	if (ntoassb3 == 1) {
		if (locdiag) cout << "process for last clue in b3" << endl;
		int x,u;
		if (!(u= xq.GetAndout())) return; // no single clue
		while (bitscanforward(x, u)) {
			register int bit= 1 << x;
			u ^= bit;
			if (IsValidB3(bf | bit)) u &= anduab3;
			else Out17(bf | bit);
		}
		return;
	}

	if (ntoassb3 < 5)p_cpt2g[80]++;else p_cpt2g[81]++;
	if (ntoassb3 >= 4) {
		//if (locdiag)xq.Status();
		GoB3Expand_1_3(bf,ac,locdiag);
		//if (locdiag)cout << " back expand 1_3" << endl;
		return;
	}
	// go direct with t3b
	BuildSortT3b();
	if (knownt >= 20) 	Dumpt3b();
	
	SP3 sp3;
	sp3.active = ac;
	sp3.all_previous = bf;
	GoB3Expand_4_x(sp3);
}

//________________________________________________________

void G17B::BuildSortT3b() {
	uint32_t tx[7][100], ntx[7];
	memset(ntx, 0, sizeof ntx);
	for (uint32_t i = 0; i < xq.nout; i++) {
		register uint32_t U = xq.tout[i] , cc = _popcnt32(U);
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
	uint32_t tx[7][100], ntx[7];
	memset(ntx, 0, sizeof ntx);
	register uint32_t A = active;
	for (uint32_t i = 0; i < nt3o; i++) {
		register uint32_t U = t3o[i]&A, cc = _popcnt32(U);
		if (cc >6)cc=6;
		uint32_t* t = tx[cc], n = ntx[cc];
		for (uint32_t j = 0; j < n; j++) {// clear =
			if (t[j] ==U) { U = 0; break; }
		}
		if(U)tx[cc][ntx[cc]++] = U;
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

void G17B::GoB3Expand_1_3(uint32_t bf, uint32_t ac, int debug ) {
	int locdiag = debug;// = 2 * (p_cpt2g[8] == 933);

	if (locdiag) {
		cout << "expand entry nto ass =" << ntoassb3
			<< " clean_valid_done=" << clean_valid_done << endl;
		xq.DumpOut();
	}
	SP3 spb[5];
	spb[0].all_previous = bf;	
	spb[0].active = ac&BIT_SET_27;	
	spb[0].indtw3 = 0;// initial first ua
	spb[0].possible_cells = xq.tout[0];
	if (locdiag) cout << Char27out(xq.tout[0]) << " initial assigned" << endl;
//	if (debug) return;
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
		for (uint32_t i = spb[0].indtw3 + 1; i < xq.nout; i++) {
			register uint32_t U = xq.tout[i];
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
		for (uint32_t i = spb[1].indtw3 + 1; i <xq.nout; i++) {
			register uint32_t U = xq.tout[i];
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
		for (uint32_t i = spb[2].indtw3 + 1; i < xq.nout; i++) {
			register uint32_t U = xq.tout[i];
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
	//if (debug) goto next3;
	{// build a reduced table for the next clues
		nt3b = 0;
		uint32_t tx[8][50], ntx[7];// gives 60 for 6 and more 
		memset(ntx, 0, sizeof ntx);
		register uint32_t A = spb[2].active, F = spb[3].all_previous;
		for (uint32_t i = spb[2].indtw3 + 1; i < xq.nout; i++) {
			register uint32_t U =xq.tout[i];
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
		GoB3Expand_4_x(spb[3],locdiag);
		if (clean_valid_done == 2) return;
		goto next3;
	}

	goto next3;// not expected
}

void G17B::GoB3Expand_4_x(SP3 spe, int debug ) {
	//if (debug) return;
	int ntoass = ncluesb3 - _popcnt32(spe.all_previous);
	int locdiag = debug;// = 2 * (p_cpt2g[8] == 933);
	//if (ntoassb3 > 4) locdiag = 1;
	//if (op.ton)if (p_cpt2g[7] == op.f7) locdiag = op.ton;
	//if (knownt >= 20)locdiag = 2;
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
		register uint32_t andw = sn->active ;
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
			cout << Char27out(andw) << " andw to go nt3more="<< nt3more<< " andw="<<andw << endl;
		if (!andw)goto next;  
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
			<< " [4] " << p_cpt2g[4] << " [10] " << p_cpt2g[10] << endl;

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
			<< " [4] " << p_cpt2g[4] << " [10] " << p_cpt2g[10] << endl;

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