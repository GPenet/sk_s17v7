

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
	for (register int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& b3 = genb12.bands3[ib3];
		b3.g.gsocket2all = b3.g.gsocket2 | (b3.g.gsocket4 | b3.g.gsocket6);
		memset(b3.isg2, 0, sizeof b3.isg2);
		memset(b3.isg3, 0, sizeof b3.isg3);
		for (int i = 0; i < 81; i++) {
			if (b3.g.gsocket2all.On(i)) {
				if (b3.g.gsocket2.On(i))b3.isg2[i] = 2;
				else if (b3.g.gsocket4.On(i))b3.isg2[i] = 4;
				else b3.isg2[i] = 6;
			}
			if (b3.g.gsocket3.On(i))b3.isg3[i] = 1;
		}
	}

	if (op.ton > 1) genb12.bands3[0].DumpIndex();
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
			cout << b3.band << " b3 i=" << i 
				<<" bn="<<b3.i416<<" bx="<< t416_to_n6[ b3.i416] << endl;
			if (op.ton > 3) {
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
	if(op.f0)if (p_cpt2g[0] < op.f0) return;
	p_cpt2g[1] += genb12.nband3;
	nb3_not_found = genb12.nband3;

	if(sgo.vx[5])		cout <<myband2.band<<" [0] "<< p_cpt2g[0] << endl;
	//if (p_cpt2g[0] > 1) {		aigstop = 1; return;	}

	StartInit();// do all preliminary setups
	if(op.ton)	StartPrint();
	UaCollector();//_do uas/guas2_3 initial harvest
	StartAfterUasHarvest();
	if (sgo.bfx[4] & 8) {
		aigstop = 1;
		guah54n.Status(2);
	}
	if (op.out_one) aigstop = 0;// reinit stop at first

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

	guah.SortClean();//60?

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
	guah.SortClean(); //80?

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
	guah54n.Build();
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
		// p1+p2 must now accept 7
		if (nb1 >7 || nb1 < 2) return 1; // no 66
		if (nb1 == 7) Ac &= ~(uint64_t)BIT_SET_27;// only b2
		if (nb1 == 2) Ac &= BIT_SET_27;// only b1
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

void GUAH54N::Build() {// cut to 30 switch to 54 killer
	Init();
	for (int i81 = 0; i81 < 81; i81++) {
		if (g17b.gsock2.On(i81)) {
			GUAH::GUA& g0 = guah.tg2[i81];
			uint32_t n = g0.nua;
			if (0 && n == 1) {
				register uint64_t U = g0.tua[0];
				U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
				BLOC0& b0 = bloc0[nbl0++];
				b0.ua = U; b0.id = i81;
				continue;
			}
			indg2[i81] = nzzused;
			Z128& myz = zz[nzzused++];
			myz.Init(0, i81);
			for (uint32_t i = 0; i < n; i++) {
				register uint64_t U = g0.tua[i];
				//if (n > 30 && nn > 16) break;
				U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
				myz.Enter(U);
			}
		}
	}
	nzzg2= nzzused;
	for (int i81 = 0; i81 < 81; i81++) {
		if (g17b.gsock3.On(i81)) {
			GUAH::GUA& g0 = guah.tg3[i81];
			uint32_t n = g0.nua;
			if (0 && n == 1) {
				register uint64_t U = g0.tua[0];
				U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
				//cout << "using bloc0 for " << nbl0 << endl;
				BLOC0& b0 = bloc0[nbl0++];
				b0.ua = U; b0.id = i81+100;
				continue;
			}
			//cout << "using nzzused for " << nzzused << endl;
			indg3[i81] = nzzused;
			Z128& myz = zz[nzzused++];
			myz.Init(1, i81);
			for (uint32_t i = 0; i < n; i++) {
				register uint64_t U = g0.tua[i];
				U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
				myz.Enter(U);
			}
		}
	}
	if(op.ton)Status();
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

//#define T90 162366
// 


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
			}
			t54b12.AddC(ua54);
		}
		ua_ret7p = ua54;// return last (smaller)
		return 1;
	}
	return 0;
}


#define SKTEXA(I,J,C) \
bitscanforward64(C, p);\
register uint64_t bit = (uint64_t)1 << C;\
I.possible_cells ^= bit;\
I.active_cells ^= bit;\
J = I;\
J.all_previous_cells |= bit

#define CALLB3(I,J) \
ncluesb12 = I; ncluesb3 = J;\
GoCallB3Com();if (aigstop)return

#define KNOWNX(I) \
if (op.known) {	if (knownt == I)return;\
if (!((~pk54) & myb12)) { \
cout << Char54out(myb12) << " expected "<<I  << endl; \
knownt = I;	}}


void G17B::Expand_03() {
	if (aigstop) return;
	SPB03A sp0, sp1,sp2, sp3;
	T54B12::TUVECT& tuv128 = t54b12.ta128[0];
	uint64_t* twu = tuv128.t;
	memset(&sp0, 0, sizeof sp0);
	sp0.active_cells = maskLSB[54].u64[0];
	sp0.possible_cells = twu[0];
	sp0.v = tuv128.v0;// initial nothing done
next1:
	{
		uint64_t p = sp0.possible_cells;
		if (!p) return;
		SKTEXA(sp0, sp1, tc_1_3[0]);
		sp1.v &= tuv128.vc[tc_1_3[0]];
		if (!(sp1.possible_cells = twu[sp1.v.getFirst128()] & sp0.active_cells))goto next1;
	}
next2:
	{
		uint64_t p = sp1.possible_cells;
		if (!p) goto next1;
		SKTEXA(sp1, sp2, tc_1_3[1]);
		sp2.v &= tuv128.vc[tc_1_3[1]];
		if (!(sp2.possible_cells = twu[sp2.v.getFirst128()] & sp1.active_cells))goto next2;
	}
next3:
	{
		uint64_t p = sp2.possible_cells;
		if (!p) goto next2;
		SKTEXA(sp2, sp3, tc_1_3[2]);
		sp3.v &= tuv128.vc[tc_1_3[2]];
		p_cpt2g[3]++;
		if (t54b12.Build_tb128(sp3)) goto next3;
		Set3(sp3);
		if (op.known > 1) {
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
	int locdiag = 0;
	if (op.known ) {
		cout << Char54out(s3.all_previous_cells) << " entry 3 [3] "	<< p_cpt2g[3] << endl;
		if (!((~pk54) & s3.all_previous_cells)) {
			cout << Char54out(s3.all_previous_cells) << " expected 3 [3] "
				<< p_cpt2g[3] << endl;
			knownt = 3;
			locdiag = 1;
		}
	}


	T54B12::TUVECT& tuv128 = t54b12.tb128[0];
	uint64_t* twu = tuv128.t;
	SPB03A sp3, sp4, sp5, sp6;
	sp3 = s3;
	sp3.possible_cells = twu[0];
	sp3.v = tuv128.v0;// initial nothing done
	if (op.ton) {
		if (op.ton>1)cout << Char54out(s3.all_previous_cells) << " 3clues [3]" << p_cpt2g[3]  << endl;
		if (op.f3) {
			if (p_cpt2g[3] == op.f3) {
				cout << "call 4_6 good path cells "<< tc_1_6[0] << " " << tc_1_6[1] << " " << tc_1_6[2]
					<< " [3]"<< p_cpt2g[3] << " start [4]" << p_cpt2g[4] << endl;

				if (op.ton > 2) {
					cout << Char54out(s3.active_cells) << " active"  << endl;
					tuv128.Dump(30);
				}
				locdiag = 1;
			}
			else {
				if (p_cpt2g[3] > op.f3) { cout << "stop  [3]" << p_cpt2g[3] << "  [4]" << p_cpt2g[4] << endl;	aigstop = 1; return; }
				if (!(op.upto3)) return;
			}
		}
	}
next4:
	{
		uint64_t p = sp3.possible_cells;
		if (!p) return;
		SKTEXA(sp3,sp4,tc_4_6[0]);
		sp4.v &= tuv128.vc[tc_4_6[0]];
		if (!(sp4.possible_cells = twu[sp4.v.getFirst128()] & sp3.active_cells))goto next4;
	};
next5:
	{
		uint64_t p = sp4.possible_cells;
		if (!p) goto next4;
		SKTEXA(sp4, sp5, tc_4_6[1]);
		sp5.v &= tuv128.vc[tc_4_6[1]];
		if (!(sp5.possible_cells = twu[sp5.v.getFirst128()] & sp4.active_cells))goto next4;
	}
next6:
	{
		uint64_t p = sp5.possible_cells;
		if (!p) goto next5;
		SKTEXA(sp5, sp6, tc_4_6[2]);
		sp6.v &= tuv128.vc[tc_4_6[2]];
		if (!(sp6.possible_cells = twu[sp6.v.getFirst128()] & sp5.active_cells)) {
			// no possible valid 6 clues should not come
			if (zh2b[1].IsValid(sp6.all_previous_cells)) {
				uint32_t i = zh2gxn.nua - 1;
				register uint64_t ua = zh2gxn.tua[i],ua54 = 
					(ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
				t54b12.AddA(ua54);		t54b12.AddB(ua54);
			}
			else {cout << "bug exp 4-7 lower 7" << endl;aigstop = 1; return;}
		}
		if (t54b12.Build_tc128( sp6)) goto next6;
		p_cpt2g[4]++;
		Set6(sp6);
		if (op.known ) {
			if (knownt >= 9)return;
			if (op.ton > 2)cout << "ua6c  nc128=" << t54b12.nc128
				<< "[4] " << p_cpt2g[4] << endl;
			if (!((~pk54) & sp6.all_previous_cells)) {
				cout << Char54out(bf_cl6) << " expected 6 [4] " <<p_cpt2g[4]  << endl;
				knownt = 6;
				Expand_7_9(sp6);	aigstop = 1;	return;;
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
	guah54n.Build6(tc_1_6);
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		genb12.bands3[ib3].Set6clues(tc_1_6);
	}
	//if (knownt == 6)guah54.DumpB2(1);
	SPB03A sp6,sp7,sp8,sp9,sp10;
	T54B12::TUVECT& tuv128 = t54b12.tc128[0];
	uint64_t* twu = tuv128.t;
	uint32_t * tcells=tc_7_9;
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
				cout << Char54out(s6.all_previous_cells) << "call 7_9 good path[4]" << p_cpt2g[4] << endl;
				cout << Char54out(s6.active_cells) << " active " << endl;
				if (op.ton > 2) {tuv128.Dump(30); locdiag = 1;	}
				locdiag = 1;
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
		SKTEXA(sp6, sp7, tcells[0]);
		if (op.known > 1) {
			if (knownt >= 9)return;
			if (!((~pk54) & sp7.all_previous_cells))
				cout << Char54out(sp7.all_previous_cells) << " expected 7 " << endl;
		}
		if(locdiag)cout << Char54out(sp7.all_previous_cells) << "  7 " << endl;
		// update vectors and get next ua
		register uint64_t U = 0;
		sp7.v &= tuv128.vc[tcells[0]];
		if (sp7.v.isNotEmpty()) U = twu[sp7.v.getFirst128()];
		else {
			if (t54b12.ncblocs) {
				register T54B12::TUVECT& tv1 = t54b12.tc128[1];
				BF128 w1 = tv1.v0 & tv1.vc[tcells[0]];
				if (w1.isNotEmpty()) U = tv1.t[w1.getFirst128()];
			}
		}
		if (!U) {// unlikely can be a valid 7 (4 known cases)
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
		SKTEXA(sp7, sp8, tcells[1]);
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
		SKTEXA(sp8, sp9, tcells[2]);
		sp9.v &= tuv128.vc[tcells[2]];
		Set9(sp9);
		build9done = 0;
		myb12 = sp9.all_previous_cells;
		KNOWNX(9)
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
		clean_valid_done = 0;
		if (!U) {// this can be a valid 9
			if (!IsValid7pbf(sp9.all_previous_cells)) {
				Valid9(sp9);
				if (aigstop) return;	goto next9;
			}
			else U = ua_ret7p;
		}
		sp9.possible_cells = U & sp9.active_cells;// never empty
		if (knownt == 9)cout << Char54out(sp9.possible_cells) << "possible after 9" << endl;
	}
	if (!sp9.possible_cells)goto next9;// dead branch
	if (op.t18) {
		if (op.p1) {
			Expand_10_11_18(sp9);
			if (aigstop) return;	goto next9;
		}
		else {// must expand 10_12
			if (!t54b12.Build_td128(sp9)) {
				Expand_10_12(sp9);
				if (aigstop) return;	goto next9;
			}
			else {
				if (knownt == 9) {
					cout << " whith known empty td after 9 see what to do" << endl;
					t54b12.DebugC();
				}
				if (aigstop) return;	goto next9;
			}
		}
	}
	//now 17 p1 max 10 in b1+2 or p2 11 in b1+2
	if (op.p2) {
		Expand_10_11_17(sp9);
		if (aigstop) return;	goto next9;
	}
next10last://this is the last 17 p1 (10 + 7)
	{
		uint64_t p = sp9.possible_cells;
		if (!p) goto next9;
		SKTEXA(sp9, sp10, tcells[3]);
		sp10.v &= tuv128.vc[tcells[3]];
		if (sp10.v.isNotEmpty())goto next10last;
		for (uint32_t i = 1; i <= t54b12.ncblocs; i++) {
			register T54B12::TUVECT& tvi = t54b12.tc128[i];
			BF128 wi = tvi.v0 & (tvi.vc[tcells[3]] & tvi.vc[tcells[2]]);
			wi &= (tvi.vc[tcells[1]] & tvi.vc[tcells[0]]);
			if (wi.isNotEmpty())goto next10last;
		}
		// this is a 10 + 7 to process
		myb12 = sp10.all_previous_cells ;
		guah54n.GetG2G3_10(myb12, tcells[3]);
		CALLB3(10, 8);
		goto next10last;
	}

}
void G17B::Valid9(SPB03A& s9) {// can be any process
	if (op.f4 && p_cpt2g[4] == op.f4)cout << "valid 9" << endl;
	p_cpt2g[56]++;
	clean_valid_done = 1;
	DoBuild9();
	if (op.t18) {
		if (op.p1) {// 9;10;11 clues  b12
			myb12 = s9.all_previous_cells;
			guah54n.GetG2G3_9(myb12);// try direct
			CALLB3(9, 9);
			uint64_t Ac = s9.active_cells;
			while (bitscanforward64(tc_10_12[0], Ac)) {
				uint64_t bit1 = (uint64_t)1 << tc_10_12[0];
				Ac ^= bit1; //clear bit
				myb12 = s9.all_previous_cells | bit1;
				guah54n.GetG2G3_10(myb12, tc_10_12[0]);
				CALLB3(10, 8);
				uint64_t Ac2 = Ac;// others are not active now
				{
					register uint64_t  A = s9.all_previous_cells | bit1,
						nb1 = _popcnt64(A & BIT_SET_27),
						nb2 = 10 - nb1;// here 10 clues 
					if (nb1 > 7 || nb2 > 7) return;
					if (nb1 == 7) Ac2 &= ~(uint64_t)BIT_SET_27;
					if (nb2 == 7) Ac2 &= BIT_SET_27;
				}
				while (bitscanforward64(tc_10_12[1], Ac2)) {
					{
						register uint64_t bit2 = (uint64_t)1 << tc_10_12[1];
						Ac2 ^= bit2; //clear bit
						myb12 = s9.all_previous_cells | bit1 | bit2;
						guah54n.GetG2G3_11(myb12, tc_10_12[0], tc_10_12[1]);
						CALLB3(11, 7);
					}
				}
			}
		}
		else {// 12 clues 666 + p1
			myb12 = s9.all_previous_cells;
			guah54n.GetG2G3_9(myb12);// try direct
			CALLB3(9, 9);
			uint64_t Ac = s9.active_cells;
			while (bitscanforward64(tc_10_12[0], Ac)) {
				uint64_t bit1 = (uint64_t)1 << tc_10_12[0];
				Ac ^= bit1; //clear bit
				myb12 = s9.all_previous_cells | bit1;
				guah54n.GetG2G3_10(myb12, tc_10_12[0]);
				CALLB3(10, 8);
				uint64_t Ac2 = Ac;// others are not active now
				{
					register uint64_t A = s9.all_previous_cells | bit1,
						nb1 = _popcnt64(A & BIT_SET_27),
						nb2 = 10 - nb1;// here 10 clues 
					if (nb1 > 7 || nb2 > 7) return;
					if (nb1 == 7) Ac2 &= ~(uint64_t)BIT_SET_27;
					if (nb2 == 7) Ac2 &= BIT_SET_27;
				}
				while (bitscanforward64(tc_10_12[1], Ac2)) {
					{
						register uint64_t bit2 = (uint64_t)1 << tc_10_12[1];
						Ac2 ^= bit2; //clear bit
						bit2 |= bit1;
						myb12 = s9.all_previous_cells  | bit2;
						guah54n.GetG2G3_11(myb12, tc_10_12[0], tc_10_12[1]);
						CALLB3(11, 7);
						uint64_t Ac3 = Ac2;// others are not active now
						{
							register uint64_t A = s9.all_previous_cells | bit2,
								nb1 = _popcnt64(A & BIT_SET_27),
								nb2 = 10 - nb1;// here 10 clues 
							if (nb1 > 6 || nb2 > 6) continue;
							if (nb1 == 6) Ac3 &= ~(uint64_t)BIT_SET_27;
							if (nb2 == 6) Ac3 &= BIT_SET_27;
						}
						while (bitscanforward64(tc_10_12[2], Ac3)) {
							register uint64_t bit3 = (uint64_t)1 << tc_10_12[2];
							Ac3 ^= bit3; //clear bit
							bit3 |= bit2;
							myb12 = s9.all_previous_cells | bit3;
							guah54n.GetG2G3_12(myb12, tc_10_12);
							CALLB3(12, 6);
						}
					}
				}
			}
		}
	}
	else if(op.p1) {// t17 p1 9 10
		myb12 = s9.all_previous_cells;
		guah54n.GetG2G3_9(myb12);// try direct
		CALLB3(9, 8);
		uint64_t Ac = s9.active_cells;
		{
			register uint64_t nb1 = _popcnt64(myb12 & BIT_SET_27),
				nb2 = 10 - nb1;// here 10 clues 
			if (nb1 > 7 || nb2 > 7) return;
			if (nb1 == 7) Ac &= ~(uint64_t)BIT_SET_27;
			if (nb2 == 7) Ac &= BIT_SET_27;	
		}
		while (bitscanforward64(tc_10_12[0], Ac)) {
			uint64_t bit1 = (uint64_t)1 << tc_10_12[0];
			Ac ^= bit1; //clear bit
			myb12 = s9.all_previous_cells | bit1 ;
			guah54n.GetG2G3_10(myb12, tc_10_12[0]);
			CALLB3(10, 7);
		}
	}
	else {// t17  p2 11+6
		uint64_t Ac =s9.active_cells;
		myb12 = s9.all_previous_cells;
		guah54n.GetG2G3_9(myb12);// try direct
		if (!op.p2b)CALLB3(9, 8);
		{
			register uint64_t nb1 = _popcnt64(myb12 & BIT_SET_27),
				nb2 = 10 - nb1;// here 10 clues 
			if (nb1 > 7 || nb2 > 7) return;
			if (nb1 == 7) Ac &= ~(uint64_t)BIT_SET_27;
			if (nb2 == 7) Ac &= BIT_SET_27;
		}
		while (bitscanforward64(tc_10_12[0], Ac)) {
			uint64_t bit1 = (uint64_t)1 << tc_10_12[0];
			Ac ^= bit1; //clear bit
			myb12 = s9.all_previous_cells | bit1 ;
			guah54n.GetG2G3_10(myb12, tc_10_12[0]);
			if (!op.p2b)CALLB3(10, 7);
			uint64_t Ac2 = Ac;// others are not active now
			{
				register uint64_t A = s9.all_previous_cells | bit1,
					nb1 = _popcnt64(A & BIT_SET_27),
					nb2 = 10 - nb1;// here 10 clues 
				if (nb1 > 6 || nb2 > 6) continue;
				// 656,566 in p2a 566 in p2b
				if (nb1 == 6) Ac2 &= ~(uint64_t)BIT_SET_27;
				if (nb2 == 6) Ac2 &= BIT_SET_27;
				if (op.p2b) {
					if (nb2 > 5) continue;
					if (nb2 == 5)  Ac2 &= BIT_SET_27;
				}
			}
			while (bitscanforward64(tc_10_12[1], Ac2)) {
				{
					register uint64_t bit2 = (uint64_t)1 << tc_10_12[1];
					Ac2 ^= bit2; //clear bit
					myb12 = s9.all_previous_cells | bit1 | bit2;
					guah54n.GetG2G3_11(myb12, tc_10_12[0], tc_10_12[1]);
					CALLB3(11, 6);
				}
			}
		}
	}
}
void G17B::EndExpand_7_9() {
	if (!(ntbelow[0] | ntbelow[1])) return;// nothing to do
	clean_valid_done = 1;// always done here

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
		}
		else {
			if (knownt > 11)return;
			if (ntbelow[0]) Go_7_12(); // push to 12 clues 666
			if (ntbelow[1]) Go_8_12();
		}
	}
	else {
		if (op.p1) {
			if (ntbelow[0]) Go_7_10(); // do 7 clues then more
			if (ntbelow[1]) Go_8_10(); // do 8 clues then more
		}
		else {
			if (ntbelow[0]) Go_7_11_17(); // push to 11 clues 656
			if (ntbelow[1]) Go_8_11_17();
		}
	}
}


void G17B::Expand_10_11_17(SPB03A& s9) {// 17p2 11+6 (65 56) + 6
	if (aigstop) return;
	if (knownt >= 9)return;
	SPB03A  sp9, sp10;
	T54B12::TUVECT& tuv128 = t54b12.tc128[0];
	uint64_t* twu = tuv128.t;
	uint32_t* tcells = tc_7_9;
	sp9 = s9;
	// now 11 clues 48 p1  or 17 p2
next10:
	{
		uint64_t p = sp9.possible_cells;
		if (!p) return;
		SKTEXA(sp9, sp10, tcells[3]);

		if (op.known > 1) {
			if (knownt > 9) return;
			if (knownt == 9)cout << Char54out(sp10.all_previous_cells) << " 10" << endl;
			if (!((~pk54) & sp10.all_previous_cells)) {
				cout << Char54out(sp10.all_previous_cells) << " expected 10 " << endl;
				knownt = 10;
			}
		}
		// next if last, find and of remaining uas
		register uint64_t Uand = ~0, U = 0;
		sp10.v &= tuv128.vc[tcells[3]];
		if (sp10.v.isNotEmpty()) { U = 1;	tuv128.DoAnd(sp10.v, Uand); }
		if (knownt == 10)cout << Char54out(Uand) << " after first bloc " << endl;
		if (U && (!Uand)) goto next10;
		for (uint32_t i = 1; i <= t54b12.ncblocs; i++) {
			register T54B12::TUVECT& tvi = t54b12.tc128[i];
			BF128 wi = tvi.v0 & (tvi.vc[tcells[3]] & tvi.vc[tcells[2]]);
			wi &= (tvi.vc[tcells[1]] & tvi.vc[tcells[0]]);
			if (wi.isNotEmpty()) { tvi.DoAnd(wi, Uand); U = 1; }
			if (U && (!Uand)) goto next10;
		}
		if (!U) {//can be a valid 10 uaand empty
			if (!IsValid7pbf(sp10.all_previous_cells)) {
				BF128 w;// valid 10
				w.bf.u64[0] = sp10.all_previous_cells;
				w.bf.u64[1] = sp10.active_cells;
				tbelow10[ntbelow[3]++] = w;

				p_cpt2g[57]++;
				if (knownt == 10)  return;
				goto next10;
			}
			else Uand = anduab12;
		}
		// this is a 10 to push to 11 656 566
		{
			register uint64_t P = Uand;
			uint64_t cc = _popcnt64(sp10.all_previous_cells & BIT_SET_27);// band 1 count
			if (cc > 6 || cc < 3) goto next10;
			if (cc == 6)P &= ~(uint64_t)BIT_SET_27;
			if (cc == 3)P &= BIT_SET_27;
			Uand = P;
		}

		Uand &= sp10.active_cells;
		if (!Uand)goto next10;
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
			guah54n.GetG2G3_11(myb12, tcells[3], tcells[4]);
			clean_valid_done = 0;
			CALLB3(11, 6);	
			if (knownt == 11) return;
		}
		goto next10;
	}

}


void G17B::Expand_10_12(SPB03A& s9) {
	if (aigstop) return;
	p_cpt2g[90]++;
	int locdiag = 0;
	if (!((~pk54) & bf_cl9)) 		locdiag = 1;


	SPB03A   sp9, sp10,sp11,sp12;
	T54B12::TUVECT& tuv128 = t54b12.td128[0];
	uint64_t* twu = tuv128.t;
	sp9 = s9;
	sp9.possible_cells = twu[0];
	sp9.v = tuv128.v0;
	memset(&ntbelow[3], 0, 3* sizeof ntbelow[0]);// 10 11  12
	if (op.f4) {
		if (p_cpt2g[4] == op.f4) {
			cout  << "call 10_12 good path[90]" << p_cpt2g[90] << endl;
			if (op.ton > 2)  tuv128.Dump(30); 
			locdiag = 1;
		}

	}
next10:	//_______ add clue 10
	{
		uint64_t p = sp9.possible_cells;
		if (!p) {	EndExpand_10_12();	return;		}
		SKTEXA(sp9, sp10, tc_10_12[0]);
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
		if ((nb1 >7)||(nb1 <3)) goto next10;
		if (nb1 == 7) sp10.possible_cells &= ~(uint64_t)BIT_SET_27;// only b2
		if (nb1 == 3) sp10.possible_cells &= BIT_SET_27;// only b1

		//if (nb1 == 6) sp10.possible_cells&= ~(uint64_t)BIT_SET_27;// only b2
		//if (nb1 == 4) sp10.possible_cells &= BIT_SET_27;// only b1

	}

next11: //add clue 11  
	{
		uint64_t p = sp10.possible_cells;
		if (!p) goto next10;
		SKTEXA(sp10, sp11, tc_10_12[1]);
		if (op.known > 1) {
			if (!((~pk54) & sp11.all_previous_cells)) {
				cout << Char54out(sp11.all_previous_cells) << " expected 11 " << endl;
				knownt = 11;
			}
		}
		if (locdiag > 1) 
			cout << Char54out(sp11.all_previous_cells) << " new 11 " << endl;
		
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
			//if (locdiag)cout << Char54out(Uand) << " not  11 " << endl;
		}
		//this is a 11 to push to 12
		if (locdiag>1 || knownt == 11)cout << Char54out(Uand) << "  <-11 to push to 12" << endl;
		uint64_t nb1 = _popcnt64(sp11.all_previous_cells & BIT_SET_27);

		if ((nb1 >6)||(nb1 <5)) goto next11;
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
			tfull[ ntbelow[5]++]= sp11.all_previous_cells | bit2;
		}
		goto next11;
	}
}
void G17B::EndExpand_10_12() {// this is 18 pass 2 + pass1 first NED
	if (op.f4) {
		if (p_cpt2g[4] == op.f4) {
			DumpPotential();
			cout << "end 10_12 good path[90]" << p_cpt2g[90] << endl;
		}

	}
	if (!(ntbelow[3] | ntbelow[4] | ntbelow[5])) return;
	p_cpt2g[91]++;	p_cpt2g[95] += ntbelow[3];
	p_cpt2g[96] += ntbelow[4];	p_cpt2g[97] += ntbelow[5];
	if(p_cpt2g[98] < ntbelow[5])p_cpt2g[98] = ntbelow[5];
	int locdiag = 0;	//if (p_cpt2g[91] == 225393) locdiag = 1;
	//if (p_cpt2g[91] > 4) return;
	if (op.known &&!((~pk54) & bf_cl9)) 		locdiag = 1;
	guah54n.Build9(tc_7_9);


	if (ntbelow[5]) {
		for (uint32_t i = 0; i < ntbelow[5]; i++) {
			p_cpt2g[92]++;
			myb12 = tfull[i];
			KNOWNX(12)	if (locdiag)cout << Char54out(myb12) << " 12" << endl;
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
			guah54n.GetG2G3_12(myb12, tc_10_12);
			CALLB3(12, 6);
			if (knownt == 12) return;
		}
	}

	clean_valid_done = 1;// done for valid 10 11

	if (ntbelow[4]) 	for (uint32_t iv = 0; iv < ntbelow[4]; iv++) {
		if (locdiag) cout << "Go_11_12() in test iv=" << iv << endl;
		BF128 ww = tbelow11[iv];
		uint64_t bf11 =  ww.bf.u64[0],y= bf11 & ~bf_cl9;
		bitscanforward64(tc_10_12[0],y);
		bitscanreverse64(tc_10_12[1], y);
		myb12 = bf11 ;
		KNOWNX(11)
		guah54n.GetG2G3_11(myb12, tc_10_12[0],tc_10_12[1]);
		CALLB3(11,7);

		uint64_t Ac = ww.bf.u64[1];//  one more clue in bands 1+2
		{
			register uint64_t nb1 = _popcnt64(bf11 & BIT_SET_27),
				nb2 = 11 - nb1;
			if (nb1 > 6 || nb2 > 6) continue;
			if (nb1 == 6) Ac &= ~(uint64_t)BIT_SET_27;
			if (nb2 == 6) Ac &= BIT_SET_27;
		}

		while (bitscanforward64(tc_10_12[2], Ac)) {
			uint64_t bit = (uint64_t)1 << tc_10_12[2];
			Ac ^= bit; //clear bit
			myb12 = bf11 | bit;
			KNOWNX(12)
			guah54n.GetG2G3_12(myb12, tc_10_12);
			CALLB3(12, 6);
		}

	}
	if (ntbelow[3]) 	for (uint32_t iv = 0; iv < ntbelow[3]; iv++) {
		if (locdiag) cout << "Go_10_12() in test iv=" << iv << endl;
		BF128 ww = tbelow10[iv];
		uint64_t bf10 = ww.bf.u64[0], y = bf10 & ~bf_cl9;
		bitscanforward64(tc_10_12[0], y);
		myb12 = bf10 ;
		KNOWNX(10)
		guah54n.GetG2G3_10(myb12, tc_10_12[0]);
		CALLB3(10, 8);
		uint64_t Ac = ww.bf.u64[1];//  one more clue in bands 1+2
		{
			register uint64_t nb1 = _popcnt64(bf10 & BIT_SET_27),
				nb2 = 10 - nb1;
			if (nb1 > 7 || nb2 > 7) continue;
			if (nb1 == 7) Ac &= ~(uint64_t)BIT_SET_27;
			if (nb2 == 7) Ac &= BIT_SET_27;
		}

		while (bitscanforward64(tc_10_12[1], Ac)) {
			uint64_t bit = (uint64_t)1 << tc_10_12[1];
			Ac ^= bit; //clear bit
			uint64_t bf11 = bf10 | bit;
			myb12 = bf11 ;
			KNOWNX(11)
			guah54n.GetG2G3_11(myb12, tc_10_12[0],tc_10_12[1]);
			CALLB3(11, 7);
			uint64_t Ac2 = Ac;//  one more clue in bands 1+2
			{
				register uint64_t nb1 = _popcnt64(bf11 & BIT_SET_27),
					nb2 = 11 - nb1;
				if (nb1 > 6 || nb2 > 6) continue;
				if (nb1 == 6) Ac2 &= ~(uint64_t)BIT_SET_27;
				if (nb2 == 6) Ac2 &= BIT_SET_27;
			}
			while (bitscanforward64(tc_10_12[2], Ac2)) {
				uint64_t bit2 = (uint64_t)1 << tc_10_12[2];
				Ac2 ^= bit2; //clear bit
				myb12 = bf11 | bit2;
				KNOWNX(12)
				guah54n.GetG2G3_12(myb12, tc_10_12);
				CALLB3(12, 6);
			}
		}
	}
}

//__________________________________________expand and go search 18 pass 1
void G17B::Expand_10_11_18(SPB03A& s9) {
	if (aigstop) return;
	if (knownt > 9)return;
	//build9done = 0;
	int locdiag = 0;
	SPB03A  sp9, sp10;
	T54B12::TUVECT& tuv128 = t54b12.tc128[0];
	uint64_t* twu = tuv128.t;
	uint32_t* tcells = tc_7_9;
	sp9 = s9;
	if (knownt == 9) {
		cout << Char54out(sp9.possible_cells )<< " possible start known Expand_10_11_18" << endl;
		locdiag = 1;
	}
	// now 11 clues 48 p1  or 17 p2
next10:
	{
		uint64_t p = sp9.possible_cells;
		if (!p) return;
		//bitscanforward64(tcells[3], p);
		//register uint64_t bit = (uint64_t)1 << tcells[3];
		//sp9.possible_cells ^= bit;
		//sp9.active_cells ^= bit;
		//sp10 = sp9;
		//sp10.all_previous_cells |= bit;
		SKTEXA(sp9, sp10, tcells[3]);
		myb12 = sp10.all_previous_cells;
		KNOWNX(10)
		// next if last, find and of remaining uas
		register uint64_t Uand = ~0, U = 0;
		sp10.v &= tuv128.vc[tcells[3]];
		if (sp10.v.isNotEmpty()) { U = 1;	tuv128.DoAnd(sp10.v, Uand); }
		if (U && (!Uand)) goto next10;
		for (uint32_t i = 1; i <= t54b12.ncblocs; i++) {
			register T54B12::TUVECT& tvi = t54b12.tc128[i];
			BF128 wi = tvi.v0 & (tvi.vc[tcells[3]] & tvi.vc[tcells[2]]);
			wi &= (tvi.vc[tcells[1]] & tvi.vc[tcells[0]]);
			if (wi.isNotEmpty()) { tvi.DoAnd(wi, Uand); U = 1; }
			if (U && (!Uand)) goto next10;
		}
		if (!U) {//can be a valid 10 uaand empty
			if (!IsValid7pbf(sp10.all_previous_cells)) {
				p_cpt2g[57]++;  // valid 10
				clean_valid_done = 1;
				myb12 = sp10.all_previous_cells;
				DoBuild9();
				guah54n.GetG2G3_10(myb12, tcells[3]);
				CALLB3(10, 8);
				uint64_t Ac = sp10.active_cells;
				{
					register uint64_t nb1 = _popcnt64(myb12 & BIT_SET_27),
						nb2 = 10-nb1;
					if (nb1 > 7 || nb2 > 7) goto next10;
					if (nb1 == 7) Ac &= ~(uint64_t)BIT_SET_27;
					if (nb2 == 7) Ac &= BIT_SET_27;
				}
				while (bitscanforward64(tcells[4], Ac)) {
					uint64_t bit = (uint64_t)1 << tcells[4];
					Ac ^= bit; //clear bit
					myb12 = sp10.all_previous_cells|bit;
					KNOWNX(11)
					if (knownt == 11) {
						cout << "dump check 50" << endl;
						int iz = guah54n.indg2[50];
						guah54n.zz[iz].Dump(4);
					}
					guah54n.GetG2G3_11(myb12, tcells[3], tcells[4]);
					CALLB3(11, 7);
				}
				clean_valid_done = 0;
				if (knownt == 10)  return;
				goto next10;
			}
			else Uand = anduab12;
		}
		// this is a 10 to push to 11
		{
			register uint64_t P = Uand;
			uint64_t cc = _popcnt64(sp10.all_previous_cells & BIT_SET_27);// band 1 count
			if (cc > 7 || cc < 2) goto next10;
			if (cc == 7)P &= ~(uint64_t)BIT_SET_27;
			if (cc == 3)P &= BIT_SET_27;
			Uand = P;
		}
		Uand &= sp10.active_cells;

		if (knownt == 10) cout<<Char54out(Uand) << " push 10 to 11 in known" << endl;
		if (!Uand)goto next10;
		DoBuild9();
		while (bitscanforward64(tcells[4], Uand)) {
			uint64_t bit2 = (uint64_t)1 << tcells[4];
			Uand ^= bit2; //clear bit
			myb12 =  sp10.all_previous_cells | bit2;
			KNOWNX(11)
			guah54n.GetG2G3_11(myb12, tcells[3], tcells[4]);
			clean_valid_done = 0;
			CALLB3(11, 7);
			if (knownt == 11) return;
		}
		clean_valid_done = 0;
		goto next10;
	}

}
void G17B::Go_8_11_18() {// 8 clues limit 11 clues 
	if (p_cpt2g[4] == op.f4) 	cout << "go 8 to 11 18 " << ntbelow[1] << endl;

	for (uint32_t iv = 0; iv < ntbelow[1]; iv++) {
		BF128 ww = tbelow8[iv];		
		bf_clc =  ww.bf.u64[0];
		{
			register uint64_t U = bf_clc & ~bf_cl6;
			bitscanforward64(tclx[0], U);
			bitscanreverse64(tclx[1], U);
		}
		ac_clc = ww.bf.u64[1];
		Go_x_11_18();
	}
}
void G17B::Go_7_11_18() {// 7 clues limit 11 clues 
	if (p_cpt2g[4] == op.f4) 	cout << "go 7 to 11 18 " << ntbelow[0] << endl;
	
	for (uint32_t iv = 0; iv < ntbelow[0]; iv++) {
		BF128 ww = tbelow7[iv];
		uint64_t bf0 = ww.bf.u64[0];
		myb12 = bf0;
		uint64_t U = myb12 & ~bf_cl6;
		bitscanforward64(tclx[0], U);
		guah54n.GetG2G3_7(myb12, tclx[0]);// try direct
		CALLB3(7, 11);
		uint64_t Ac = ww.bf.u64[1];
		// maxi is 5+2 possible 8+2
		while (bitscanforward64(tclx[1], Ac)) {
			uint64_t bit = (uint64_t)1 << tclx[1];
			Ac ^= bit; //clear bit
			bf_clc = bf0 | bit;
			ac_clc = Ac;
			Go_x_11_18();
		}
	}
}
void G17B::Go_x_11_18() {// 8 clues limit 11 clues 
	int locdiag = 0;
	if (p_cpt2g[4] == op.f4) {
		cout << Char54out(bf_clc) << "entry Go_x_11_18 cells " << tclx[0] << " " << tclx[1] << endl;
		locdiag = 1;
	}
	myb12 = bf_clc; 
	guah54n.GetG2G3_8(myb12, tclx);// try 8
	CALLB3(8, 10);
	// try now a second clue in bands 1+2

	uint64_t Ac2 = ac_clc;
	while (bitscanforward64(tclx[2], Ac2)) {
		uint64_t bit2 = (uint64_t)1 << tclx[2];
		Ac2 ^= bit2; //clear bit
		uint64_t bf2 = bf_clc | bit2;

		myb12 = bf2;  
		guah54n.Build9(tclx);
		guah54n.GetG2G3_9(myb12);
		CALLB3(9, 9);
		// try now a third clue in bands 1+2
		uint64_t Ac3 = Ac2;
		while (bitscanforward64(tclx[3], Ac3)) {
			uint64_t bit3 = (uint64_t)1 << tclx[3];
			Ac3 ^= bit3; //clear bit
			uint64_t bf3 = bf2 | bit3;
			myb12 = bf3;
			guah54n.GetG2G3_10(myb12, tclx[3]);
			CALLB3(10, 8);

			// try now a fourth  clue in bands 1+2
			uint64_t Ac4 = Ac3;
			{
				register uint64_t nb1 = _popcnt64(myb12 & BIT_SET_27),
					nb2 = 10 - nb1;
				if (nb1 > 7 || nb2 > 7) continue;
				if (nb1 == 7) Ac4 &= ~(uint64_t)BIT_SET_27;
				if (nb2 == 7) Ac4 &= BIT_SET_27;
			}
			while (bitscanforward64(tclx[4], Ac4)) {
				uint64_t bit4 = (uint64_t)1 << tclx[4];
				Ac4 ^= bit4; //clear bit
				uint64_t bf4 = bf3 | bit4;
				myb12 = bf4; 
				guah54n.GetG2G3_11(myb12, tclx[3], tclx[4]);
				CALLB3(11, 7);
			}
		}
	}

}


//________________expand and go search 17 pass 1
void G17B::Go_8_10() {// 8 clues limit 10 clues 
	if (p_cpt2g[4] == op.f4) 	cout << "go 8 to 10 " << ntbelow[1] << endl;

	for (uint32_t iv = 0; iv < ntbelow[1]; iv++) {
		BF128 ww = tbelow8[iv];
		bf_clc = ww.bf.u64[0];
		{
			register uint64_t U = bf_clc & ~bf_cl6;
			bitscanforward64(tclx[0], U);
			bitscanreverse64(tclx[1], U);
		}
		ac_clc = ww.bf.u64[1];
		Go_x_10();
	}
}
void G17B::Go_7_10() {// 7 clues limit 10 clues 
	if (p_cpt2g[4] == op.f4) 	cout << "go 7 to 10 " << ntbelow[0] << endl;

	for (uint32_t iv = 0; iv < ntbelow[0]; iv++) {
		BF128 ww = tbelow7[iv];
		uint64_t bf0 = ww.bf.u64[0];
		myb12 = bf0;
		uint64_t U = myb12 & ~bf_cl6;
		bitscanforward64(tclx[0], U);
		guah54n.GetG2G3_7(myb12, tclx[0]);// try direct
		CALLB3(7, 10);
		uint64_t Ac = ww.bf.u64[1];
		// maxi is 5+2 possible 8+2
		while (bitscanforward64(tclx[1], Ac)) {
			uint64_t bit = (uint64_t)1 << tclx[1];
			Ac ^= bit; //clear bit
			bf_clc = bf0 | bit;
			ac_clc = Ac;
			Go_x_10();
		}
	}
}
void G17B::Go_x_10() {// 8 clues limit 10 clues 
	int locdiag = 0;
	if (p_cpt2g[4] == op.f4) {
		cout << Char54out(bf_clc) << "entry Go_x_10 cells " << tclx[0] << " " << tclx[1] << endl;
		locdiag = 1;
	}
	myb12 = bf_clc;  
	guah54n.GetG2G3_8(myb12, tclx);// try 8
	CALLB3(8, 9);
	// try now a second clue in bands 1+2
	uint64_t Ac2 = ac_clc;
	while (bitscanforward64(tclx[2], Ac2)) {
		uint64_t bit2 = (uint64_t)1 << tclx[2];
		Ac2 ^= bit2; //clear bit
		uint64_t bf2 = bf_clc | bit2;

		myb12 = bf2; 
		guah54n.Build9(tclx);
		guah54n.GetG2G3_9(myb12);
		CALLB3(9, 8);
		// try now a third clue in bands 1+2
		uint64_t Ac3 = Ac2;
		{
			register uint64_t nb1 = _popcnt64(myb12 & BIT_SET_27),
				nb2 = 10 - nb1;
			if (nb1 > 7 || nb2 > 7) continue;
			if (nb1 == 7) Ac3 &= ~(uint64_t)BIT_SET_27;
			if (nb2 == 7) Ac3 &= BIT_SET_27;
		}

		while (bitscanforward64(tclx[3], Ac3)) {
			uint64_t bit3 = (uint64_t)1 << tclx[3];
			Ac3 ^= bit3; //clear bit
			uint64_t bf3 = bf2 | bit3;
			myb12 = bf3;
			guah54n.GetG2G3_10(myb12, tclx[3]);
			CALLB3(9, 8);
		}
	}
}


//________________expand and go search 18 pass 2
void G17B::Go_8_12() {
	if (p_cpt2g[4] == op.f4) 	cout << "go 8 to 12 " << ntbelow[1] << endl;
	for (uint32_t iv = 0; iv < ntbelow[1]; iv++) {
		BF128 ww = tbelow8[iv];
		bf_clc = ww.bf.u64[0];
		{
			register uint64_t U = bf_clc & ~bf_cl6;
			bitscanforward64(tclx[0], U);
			bitscanreverse64(tclx[1], U);
		}
		ac_clc = ww.bf.u64[1];
		Go_x_12();
	}
}
void G17B::Go_7_12() {
	if (p_cpt2g[4] == op.f4) 	cout << "go 7 to 12 " << ntbelow[0] << endl;
	for (uint32_t iv = 0; iv < ntbelow[0]; iv++) {
		BF128 ww = tbelow7[iv];
		uint64_t bf0 = ww.bf.u64[0];
		uint64_t U = bf0 & ~bf_cl6;
		bitscanforward64(tclx[0], U);
		myb12 = bf0;//  7/ 11;
		guah54n.GetG2G3_7(myb12, tclx[0]);// try direct
		CALLB3(7, 11);
		uint64_t Ac = ww.bf.u64[1];
		while (bitscanforward64(tclx[1], Ac)) {
			uint64_t bit = (uint64_t)1 << tclx[1];
			Ac ^= bit; //clear bit
			bf_clc = bf0 | bit;
			ac_clc = Ac;
			Go_x_12();
		}
	}

}
void G17B::Go_x_12() {// 8 clues limit 12 clues 666
	int locdiag = 0;
	if (p_cpt2g[4] == op.f4) {
		cout << Char54out(bf_clc) << "entry Go_x_12 cells " << tclx[0] << " " << tclx[1] << endl;
		locdiag = 1;
	}
	myb12 = bf_clc;//  8/ 10;
	KNOWNX(8)
	guah54n.GetG2G3_8(myb12, tclx);// try 8
	if (locdiag) cout << "CALLB3(8, 10); [7]" << p_cpt2g[4] << endl;
	CALLB3(8, 10);
	uint64_t Ac2 = ac_clc;
	while (bitscanforward64(tclx[2], Ac2)) { 
		uint64_t bit2 = (uint64_t)1 << tclx[2];
		Ac2 ^= bit2; //clear bit
		uint64_t bf2 = bf_clc | bit2;	 // 9  9;
		guah54n.Build9(tclx);
		myb12 = bf2;
		KNOWNX(9)
		guah54n.GetG2G3_9(myb12);
		if (locdiag) cout << "CALLB3(9, 9); [7]" << p_cpt2g[4] << endl;
		CALLB3(9, 9);
		uint64_t Ac3 = Ac2;
		while (bitscanforward64(tclx[3], Ac3)) {
			uint64_t bit3 = (uint64_t)1 << tclx[3];
			Ac3 ^= bit3; //clear bit
			uint64_t bf3 = bf2 | bit3;	 //10 8 
			myb12 = bf3;  
			KNOWNX(10)
			guah54n.GetG2G3_10(myb12, tclx[3]);
			if (locdiag) cout << "CALLB3(10, 8); [7]" << p_cpt2g[4] << endl;
			CALLB3(10, 8);
			uint64_t Ac4 = Ac3;
			{
				register uint64_t nb1 = _popcnt64(bf3 & BIT_SET_27),
					nb2 =10-nb1;
				if (nb1 > 7 || nb2 > 7) continue;
				if (nb1 == 7) Ac2 &= ~(uint64_t)BIT_SET_27;
				if (nb2 == 7) Ac2 &= BIT_SET_27;
			}
			while (bitscanforward64(tclx[4], Ac4)) {
				uint64_t bit4 = (uint64_t)1 << tclx[4];
				Ac4 ^= bit4; //clear bit
				uint64_t bf4 = bf3 | bit4;	 // 11 7
				myb12 = bf4; 
				KNOWNX(11)
				guah54n.GetG2G3_11(myb12, tclx[3], tclx[4]);
				if (locdiag) cout << "CALLB3(11, 7); [7]" << p_cpt2g[4] << endl;
				CALLB3(11, 7);
				uint64_t Ac5 = Ac4;
				{
					register uint64_t nb1 = _popcnt64(bf4 & BIT_SET_27),
						nb2 = 11 - nb1;
					if (nb1 > 6 || nb2 > 6) continue;
					if (nb1 == 6) Ac5 &= ~(uint64_t)BIT_SET_27;
					if (nb2 == 6) Ac5 &= BIT_SET_27;
				}
				while (bitscanforward64(tclx[5], Ac5)) {
					uint64_t bit5 = (uint64_t)1 << tclx[5];
					Ac5 ^= bit5; //clear bit
					myb12 = bf4 | bit5;		// 12 6
					guah54n.GetG2G3_12(myb12, tclx);
					if (locdiag) cout << "CALLB3(12, 6); [7]" << p_cpt2g[4] << endl;
					CALLB3(12, 6);
				}			
			}
		}
	}
}

//________________ expand and go search 17 pass 2 
void G17B::Go_8_11_17() {// 8 clues limit 11 clues 
	//cout << " entry 8 clues for 11 clues" << endl;
	if (p_cpt2g[4] == op.f4) 	cout << "go 8 to 11_17 " << ntbelow[1] << endl;
	for (uint32_t iv = 0; iv < ntbelow[1]; iv++) {
		BF128 ww = tbelow8[iv];
		bf_clc = ww.bf.u64[0];
		{
			register uint64_t U = bf_clc & ~bf_cl6;
			bitscanforward64(tclx[0], U);
			bitscanreverse64(tclx[1], U);
		}
		ac_clc = ww.bf.u64[1];
		Go_x_11_17();
	}
}
void G17B::Go_7_11_17() {
	if (p_cpt2g[4] == op.f4) 	cout << "go 7 to 11 17 " << ntbelow[0] << endl;
	for (uint32_t iv = 0; iv < ntbelow[0]; iv++) {
		BF128 ww = tbelow7[iv];
		uint64_t bf0 = ww.bf.u64[0];
		uint64_t U = bf0 & ~bf_cl6;
		bitscanforward64(tclx[0], U);
		myb12 = bf0;//  7/ 10;
		guah54n.GetG2G3_7(myb12, tclx[0]);// try direct
		CALLB3(7, 10);
		uint64_t Ac = ww.bf.u64[1];
		while (bitscanforward64(tclx[1], Ac)) {
			uint64_t bit = (uint64_t)1 << tclx[1];
			Ac ^= bit; //clear bit
			bf_clc = bf0 | bit;
			ac_clc = Ac;
			Go_x_11_17();
		}
	}
}
void G17B::Go_x_11_17() {// 8 clues limit 11 clues 656 566
	int locdiag = 0;
	if (p_cpt2g[4] == op.f4) {
		cout << Char54out(bf_clc) << "entry Go_x_11_17 cells " << tclx[0] << " " << tclx[1] << endl;
		locdiag = 1;
	}
	myb12 = bf_clc;//  8/ 9;
	guah54n.GetG2G3_8(myb12, tclx);// try 8
	CALLB3(8, 9);
	uint64_t Ac2 = ac_clc;
	while (bitscanforward64(tclx[2], Ac2)) {
		uint64_t bit2 = (uint64_t)1 << tclx[2];
		Ac2 ^= bit2; //clear bit
		uint64_t bf2 = bf_clc | bit2;	 // 9  9;
		guah54n.Build9(tclx);
		myb12 = bf2;
		guah54n.GetG2G3_9(myb12);
		CALLB3(9, 8);
		uint64_t Ac3 = Ac2;
		{
			register uint64_t nb1 = _popcnt64(bf2 & BIT_SET_27),
				nb2 = 10 - nb1;
			if (nb1 > 7 || nb2 > 7) continue;
			if (nb1 == 7) Ac3 &= ~(uint64_t)BIT_SET_27;
			if (nb2 == 7) Ac3 &= BIT_SET_27;
		}
		while (bitscanforward64(tclx[3], Ac3)) {
			uint64_t bit3 = (uint64_t)1 << tclx[3];
			Ac3 ^= bit3; //clear bit
			uint64_t bf3 = bf2 | bit3;	 //10 8 
			myb12 = bf3;
			guah54n.GetG2G3_10(myb12, tclx[3]);
			CALLB3(10, 7);
			uint64_t Ac4 = Ac3;
			{
				register uint64_t nb1 = _popcnt64(bf3 & BIT_SET_27),
					nb2 = 10 - nb1;
				if (nb1 == 6) Ac4 &= ~(uint64_t)BIT_SET_27;
				if (nb2 == 6) Ac4 &= BIT_SET_27;
			}
			while (bitscanforward64(tclx[4], Ac4)) {
				uint64_t bit4 = (uint64_t)1 << tclx[4];
				Ac4 ^= bit4; //clear bit
				myb12 = bf3 | bit4;		// 11 7
				guah54n.GetG2G3_11(myb12, tclx[0], tclx[1]);
				CALLB3(11, 6);
			}
		}
	}
}



#define VTEST sgo.vx[9]



void G17B::GoCallB3Com() {
	p_cpt2g[7]++;
	int locdiag = 0;
	if (knownt >= 11) locdiag = 1;	
	if (op.f7 ){
		if(op.f4 && p_cpt2g[4] == op.f4)
		if (p_cpt2g[7] == op.f7)			locdiag = 1;		
	}
	if (locdiag) {
		cout << Char54out(myb12) << " :GoCallB3Com() in test clean_valid_done" << clean_valid_done << endl;
		//guah54n.StatusFiltered(myb12);
		//guah54n.g2.Print("g2");
		//guah54n.g3.Print("g3");
		//guah54n.Statusv6();
		//guah54n.Statusv9();
	}
	tcluesxpdone = 0;
	guah54n.InitCom();//open the door for g2 g3 seen later 
	int tb3[256], ntb3 = 0;
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& b3 = genb12.bands3[ib3];
		if (b3.aigskip)continue;// op first and one 18 found
		b3.ntg2_4 = b3.ntg2_6 = 0;
		p_cpt2g[23]++;

		register uint32_t  Mg2 = 0, Mg3 = 0, minc = 0,
			R1 = 0, R2 = 0, R3 = 0;
		{
			register uint32_t mini;
			for (uint32_t i = 0; i < guah54n.ntg2; i++) {
				register uint32_t i81 = guah54n.tg2[i], a;
				if ((a = b3.isg2[i81])) { 
					if (locdiag) cout << Char27out(b3.g.ua2bit27[i81]) << " ib3" << ib3 << " g2 i81 " << i81 
						<< " mini  " << b3.g.ua2_imini[i81] << endl;
					if (a == 2) {
						Mg2 |= b3.g.ua2bit27[i81];
						mini =1<< b3.g.ua2_imini[i81];
						R3 |= R2 & mini; R2 |= R1 & mini; R1 |= mini;
					}
					else if (a == 4)b3.tg2_4[b3.ntg2_4++] = i81;
					else b3.tg2_6[b3.ntg2_6++] = i81;
				}
			}
			minc = _popcnt32(R3) + _popcnt32(R1);
			if (locdiag) {
				cout << "after guas2 ib3=" << ib3 << " minc=" << minc << endl;
				cout << Char9out(R1) << " R1" << endl;
				cout << Char9out(R2) << " R2" << endl;
				cout << Char9out(R3) << " R3" << endl;
				cout << Char27out(Mg2) << " Mg2" << endl;
			}


			if (minc > (uint32_t)ncluesb3) continue;
			for (uint32_t i = 0; i < guah54n.ntg3; i++) {
				register uint32_t i81 = guah54n.tg3[i];
				if (b3.isg3[i81]) 	Mg3 |= 1 << b3.g.ua3_imini[i81];
			}
			Mg3 &= (R1 ^ 0x1ff);// free mini rows
			minc += _popcnt32(Mg3);
		}
		if (locdiag) cout << " ib3=" << ib3 << " minc=" << minc << endl;
		if (minc <= (uint32_t)ncluesb3) {
			b3.ming2_3 = R3; 			
			b3.ming2_2 = R2&~R3;
			b3.ming2_1 = R1 & ~R2;
			b3.mg3 = Mg3;
			b3.mg2 = Mg2;  
			b3.minc = minc;
			tb3[ntb3++] = ib3;
		}
	}
	if (!ntb3) return;
	p_cpt2g[24]++;
	p_cpt2g[25] += ntb3;
	if (locdiag) 
		cout  << " :GoCallB3Com() in test loop ntb3=" << ntb3 << endl;
	
	for (int i = 0; i < ntb3; i++) {
		genb12.bands3[tb3[i]].GoA();
		//if (locdiag && i > 9) break;
	}

	if (locdiag) 
		cout << " :GoCallB3Com() in test endloop ntb3"   << endl;
	

}

int maskminirow[9] = { 07,070,0700,07000,070000,0700000,07000000,070000000,0700000000 };

void STD_B3::GoAg23(int debug) {
	for (uint32_t i = 0; i < guah54n.ntmg2; i++) {
		register uint32_t i81 = guah54n.tmg2[i], a;
		if (debug)cout << "GoAg23 add g2 i81=" << i81 << endl;
			if ((a = isg2[i81])) {
			if (a == 2) {
				register uint32_t  bitmini = 1 << g.ua2_imini[i81];
				mg2 |= g.ua2bit27[i81];
				if (bitmini & ming2_2) {
					ming2_2 ^= bitmini;
					ming2_3 |= bitmini;
				}
				else if (bitmini & ming2_1) {
					ming2_1 ^= bitmini;
					ming2_2 |= bitmini;
				}
				else ming2_1 |= bitmini;
			}
			else if (a == 4)tg2_4[ntg2_4++] = i81;
			else tg2_6[ntg2_6++] = i81;
		}
	}
	for (uint32_t i = 0; i < guah54n.ntmg3; i++) {
		register uint32_t i81 = guah54n.tmg3[i];
		if (debug)cout << "GoAg23 add g3 i81=" << i81 << endl;
		if (isg3[i81]) mg3 |= 1 << g.ua3_imini[i81];
	}
	register uint32_t all = ming2_1 | ming2_2 | ming2_3;
	mg3 &= ~all;
	minc = _popcnt32(all | mg3) + _popcnt32(ming2_3);
}
void STD_B3::GoA() {
	g17b.nt3more = 0;
	CALLBAND3& cb3e = g17b.cb3;
	p_cpt2g[8]++;
	if (VTEST && p_cpt2g[8] > VTEST) return;
	int locdiag = 0;
	if (op.known && g17b.knownt >= 10) {
		cout << band << "GoA() in diag known [7]  " << p_cpt2g[7] << " [8]  " << p_cpt2g[8] << endl;
		locdiag = 1;
	}

	if (p_cpt2g[7] == op.f7) {
		cout << band << "GoA() [7]  " << p_cpt2g[7] << " [8]  " << p_cpt2g[8] << endl;
		locdiag = 1;
		cout << "xq Dump2 entry goA clean_valid_done  "<< g17b.clean_valid_done
			<< " nmore "<< guah54n.ntmg2 + guah54n.ntmg3 << endl;
		xq.Dump2();
	}

	if (g17b.clean_valid_done == 2 || g17b.aigstop) return;
	memcpy(&g17b.grid0[54], band0, sizeof band0);// used in brute force
	g17b.myband3 = this;
	if (guah54n.ntmg2 | guah54n.ntmg3) {// possible change in mincout
		if (locdiag) {

		}
		GoAg23(locdiag);
		if ((int)minc > g17b.ncluesb3) return;
	}


	// note : this is a critical code 
	xq.Init();
	xq.nb3 = g17b.ncluesb3;
	xq.nmiss = xq.nb3 - minc;

	{
		register uint32_t  Mg2 = mg2;
		// load 2 pairs in minirow

		register uint32_t M = ming2_2, x;
		while (bitscanforward(x, M)) {
			M ^= 1 << x;
			register int mask = maskminirow[x],
				v = mask & Mg2,
				i27, bit;
			xq.t1a |= (mask ^ v);// common cell
			xq.critbf |= mask;
			bitscanforward(i27, v);
			bit = 1 << i27;
			xq.t2a[xq.n2a++] = mask ^ bit;
			bit ^= v;// other i27
			xq.t2a[xq.n2a++] = mask ^ bit;

		}
		xq.t2b3 = xq.t2a + xq.n2a;
		M = ming2_3;
		while (bitscanforward(x, M)) {
			M ^= 1 << x;
			register int i0 = 3 * x;
			xq.t2b3[xq.n2b3++] = g.pat2_27[i0++];
			xq.t2b3[xq.n2b3++] = g.pat2_27[i0++];
			xq.t2b3[xq.n2b3++] = g.pat2_27[i0];
			xq.critbf |= maskminirow[x];
		}
		xq.t2b = xq.t2b3 + xq.n2b3;// now one pair
		M = ming2_1;
		while (bitscanforward(x, M)) {
			M ^= 1 << x;
			register int mask = maskminirow[x],
				v = mask & Mg2,
				a = mask ^ v;
			xq.critbf |= a;
			xq.t2b[xq.n2b++] = a;
		}
		M = mg3;// now triplets
		while (bitscanforward(x, M)) {
			M ^= 1 << x;
			register int mask = maskminirow[x];
			xq.critbf |= mask;
			xq.t2b[xq.n2b++] = mask;
		}
	}
	if (locdiag) {
		cout << "goa in diag status before call b0 b1 bmore clean_valid_done  " << g17b.clean_valid_done << endl;
		xq.Status();
	}
	if (!xq.nmiss) { GoB0(); return; }
	p_cpt2g[9]++;
	if (xq.nmiss == 1) { GoB1(); return; }
	GoBMore1();

}

//____________________ band3 miss0 at start

uint32_t  XQ::AssignMiss0(uint32_t bf) {
	nadded++;
	register uint32_t ret = 0;
	// is it a bf3 assigned
	if (n2b3) {
		register uint32_t  n = n2b3; n2b3 = 0;
		for (uint32_t i = 0; i < n; i += 3) {
			register uint32_t wa = t2b3[i] | t2b3[i + 1] | t2b3[i + 2];
			if (!(wa & bf)) {// no hit keep it
				ret |= wa;	continue;
			}
			// push the last in n2b
			wa = t2b3[i];
			if (wa & bf) {
				wa = t2b3[i+1];
				if (wa & bf) {	wa = t2b3[i + 2];}
			}
			t2b[n2b++] =wa;	ret |= wa;	n2b3 = i;
			// reload others   n2b and exit
			for (uint32_t j = i + 3; j < n; j++) {
				ret |= t2b3[j];
				t2b3[n2b3++] = t2b3[j];
			}
			for (uint32_t j = 0; j < n2b; j++) 	ret |= t2b[j];			
			return ret;
		}
		n2b3 = n; // no hit
	}
	//bf3 hits a single clue ua assigned
	register  uint32_t n = n2b;
	n2b = 0;
	for (uint32_t i = 0; i < n; i++) {
		register uint32_t U = t2b[i];
		if (U & bf) {
			for (uint32_t j = i + 1; j < n; j++) {
				register uint32_t U = t2b[j];
				t2b[n2b++] = U;		ret |=U;
			}
			break;
		}
		else {	ret |= U; n2b++;}
	}
	return ret;// should never be
}

#define MISS0ADD \
if (F & U) continue;\
if (!(U &= A)) return ; \
if (_popcnt32(U) == 1) {F |= U;	\
A =  xq.AssignMiss0(U);}\
else xq.Addin(U)

void STD_B3::GoB0() {
	p_cpt2g[40]++;
	int locdiag = 0;
	if (VTEST && p_cpt2g[8] == VTEST) locdiag = 1;
	xq.SetFilters();
	{// add band 3 UA size 4
		register uint32_t F = xq.fa, A = xq.fb;
		if (locdiag) {
			cout<<Char54out(F) << " F GoB0" << endl;
			cout << Char54out(A) << " A" << endl;
		}
		for (uint32_t i = 0; i < nua; i++) {
			register uint32_t U = tua[i];
			if( (U>>27) > 4) { xq.iuas4 = i; break; }//
			U &= BIT_SET_27;
			MISS0ADD;
		}
		for (uint32_t i = 0; i < ntg2_4; i++) {//guas 4
			int i81 = tg2_4[i];
			register uint32_t U= g.pat2[i81];
			MISS0ADD;
		}
		for (uint32_t i = 0; i < ntg2_6; i++) {//guas 6
			int i81 = tg2_6[i];
			register uint32_t U = g.pat2[i81];
			MISS0ADD;
		}
		for (uint32_t i = xq.iuas4; i < nua; i++) {// uas b3 >4
			register uint32_t U = tua[i] & BIT_SET_27;
			MISS0ADD;
		}
		if (locdiag) xq.Statusmiss0(F,A);
		if (xq.nin) {// check fresh potential 
			while (1) {// open the door tomore than one assign
				register uint32_t Fr = F, nn = xq.nin;
				xq.nin = 0;
				for (uint32_t i = 0; i < nn; i++) {
					register uint32_t U = xq.tin[i];
					MISS0ADD;
				}
				if (Fr == F) break;;
			}
		}
		if (locdiag) {
			cout << "still valid end of goB0 clean_valid_done  " << g17b.clean_valid_done << endl;
			xq.Statusmiss0(F, A);
		}

		if (_popcnt32 (F) == xq.nb3) {	GoC0F(F); return;}
		GoC0(F,A);
		//cout << "back GoC0(F,A) end of goB0" << endl;
	}
}
void STD_B3::GoC0F(uint32_t bf) {// misso all known
	register uint32_t F = bf;
	g17b.Do7x();//  expensive 4 guam

	for (uint32_t i = 0; i <= nbbgm; i++) {
		GUM64& gw = tgm64[i];
		register uint64_t V = gw.Getv(g17b.tcluesxp, g17b.ncluesxp);
		register uint32_t r;
		while (bitscanforward64(r, V)) {
			V ^= (uint64_t)1 << r;
			/*
			if (0) {
				uint64_t bit= (uint64_t)1 << r;
				cout << Char27out(gw.tb3[r]) << "\t";
				for (int j = 0; j < 54; j++)
					if (gw.vc[j] & bit) cout << ".";
					else cout << "1";
				cout  << endl;
			}
			*/
			if (!( gw.tb3[r] & F)) return;
		}
	}
	for (uint32_t i = 0; i <= nbbgmm; i++) {
		GUM64& gw = tgm64m[i];
		register uint64_t V = gw.Getv(g17b.tcluesxp, g17b.ncluesxp);
		register uint32_t r;
		while (bitscanforward64(r, V)) {
			V ^= (uint64_t)1 << r;
			/*
			if (0) {
				uint64_t bit = (uint64_t)1 << r;
				cout << Char27out(gw.tb3[r]) << "\t";
				for (int j = 0; j < 54; j++)
					if (gw.vc[j] & bit) cout << ".";
					else cout << "1";
				cout << endl;
			}
			*/
			if (!(gw.tb3[r] & F)) return;
		}
	}
	// solution to test
	//cout << Char27out(F) << " this is a solution to test" << endl;
	g17b.GoD0F(F);

}
void STD_B3::GoC0(uint32_t bf, uint32_t a) {
	int locdiag = 0;
	if (VTEST && p_cpt2g[8] == VTEST) locdiag = 1;
	register uint32_t F = bf,A=a;
	g17b.Do7x();
	xq.BuildMiss0Redundant();
	if (locdiag) {
		cout << "goC0 no redundant" << endl;
		xq.Status();
	}
	for (uint32_t i = 0; i <= nbbgm; i++) {
		GUM64& gw = tgm64[i];
		register uint64_t V = gw.Getv(g17b.tcluesxp, g17b.ncluesxp);
		register uint32_t r;
		while (bitscanforward64(r, V)) {
			V ^= (uint64_t)1 << r;
			register uint32_t U = gw.tb3[r];
			if (U & F) continue;;
			//cout << Char27out(U) << " u gm" << endl;;
			if (!(U &= A)) return;
			if (_popcnt32(U) == 1) {
				F |= U;
				A = xq.AssignMiss0(U);
			}
			else if (xq.AddRedundant(U))
				xq.tin[xq.nin++] = U;

		}
	}
	for (uint32_t i = 0; i <= nbbgmm; i++) {
		GUM64& gw = tgm64m[i];
		register uint64_t V = gw.Getv(g17b.tcluesxp, g17b.ncluesxp);
		register uint32_t r;
		while (bitscanforward64(r, V)) {
			V ^= (uint64_t)1 << r;
			register uint32_t U = gw.tb3[r];
			if (U & F) continue;;
			//cout << Char27out(U) << " u gmm" << endl;;
			if (!(U &= A)) return;
			if (_popcnt32(U) == 1) {
				F |= U;
				A = xq.AssignMiss0(U);
			}
			else if (xq.AddRedundant(U))
				xq.tin[xq.nin++] = U;
		}
	}
	if (locdiag) {
		cout << Char27out(F) << "all add done" << endl;
		xq.Status();
	}
	if (xq.nin) {// check again fresh potential 
		while (1) {// open the door tomore than one assign
			register uint32_t Fr = F, nn = xq.nin;
			xq.nin = 0;
			for (uint32_t i = 0; i < nn; i++) {
				register uint32_t U = xq.tin[i];
				//cout << Char27out(U&A) << " uin to see" << endl;
				MISS0ADD;
				//cout << Char27out(F) << " F new nin" << xq.nin << endl;
				//cout << Char27out(A) << " A" << endl;

			}
			if (Fr == F) break;;
		}
	}
	if (locdiag) {
		cout <<Char27out(F)<< " F final C0" << endl;
		xq.Status();
	}
	g17b.ntoassb3 = xq.nb3- _popcnt32(F) ;
	//cout << "[8] " << p_cpt2g[8] << " minc " << minc<<" ntoass "<< ntoass << endl;
	if (g17b.ntoassb3 < 0) return; // should never be
	if (!g17b.ntoassb3) {
		if (xq.nin)return;
		g17b.GoD0F(F); 
		return; 
	}
	// must now go to expand with ordered and cleaned file
	xq.BuildMiss0Out();
	if (g17b.ntoassb3 == 1) {
		int x, u;
		if (!(u = xq.GetAndout())) return; // no single clue
		if (g17b.VerifyValidb3())return;
		if (locdiag) 	cout << Char27out(u) << " 1 to assign try this nmore ="<<g17b.nt3more << endl;
		while (bitscanforward(x, u)) {
			register int bit = 1 << x;
			u ^= bit;
			if (g17b.IsValidB3(F | bit,locdiag)) u &= g17b.anduab3;
			else g17b.Out17(F | bit);
		}
		return;
	}
	if (g17b.ntoassb3 >= 4) {		g17b.GoB3Expand_1_3(F, A);		return;	}
	memcpy(g17b.t3b, xq.tout, xq.nout * sizeof g17b.t3b[0]);
	g17b.nt3b = xq.nout;
	SP3 sp3;
	sp3.active = A;
	sp3.all_previous = F;
	g17b.GoB3Expand_4_x(sp3);
}
void G17B::GoD0F(uint32_t bf) {
	if (VerifyValidb3())return;
	if (!IsValidB3(bf))Out17(bf);
}

/*
	cout << "[8] " << p_cpt2g[8] << " minc " << minc << " ntoass " << g17b.ntoassb3 << endl;
	cout << "final status to expand nout=" << xq.nout << endl;
	cout << Char27out(F) << " F " << endl;
	cout << Char27out(A) << " A" << endl;
	xq.DumpOut();
*/
//____________________ band3 miss1 at start
#define MISS1ADD \
if (!(U & A)){ \
wa&=U;xq.Addout(U);}\
else xq.Addin(U)

#define MISS1ADDIF if(xq.AddRedundant(U)) \
if (!(U & A)){ \
wa&=U;xq.Addout(U);}\
else xq.Addin(U)

void STD_B3::GoB1() {
	p_cpt2g[41]++;
	int locdiag = 0;
	if (VTEST && p_cpt2g[8] == VTEST) locdiag = 1;
	//xq.Status();
	{// add band 3 UA size 4
		register uint32_t  A = xq.critbf, wa = BIT_SET_27;
		for (uint32_t i = 0; i < nua; i++) {
			register uint32_t U = tua[i];
			if ((U >> 27) > 4) { xq.iuas4 = i; break; }//
			U &= BIT_SET_27;
			MISS1ADD;
		}
		for (uint32_t i = 0; i < ntg2_4; i++) {//guas 4
			int i81 = tg2_4[i];
			register uint32_t U = g.pat2[i81];
			MISS1ADD;
		}
		for (uint32_t i = 0; i < ntg2_6; i++) {//guas 6
			int i81 = tg2_6[i];
			register uint32_t U = g.pat2[i81];
			MISS1ADD;
		}
		for (uint32_t i = xq.iuas4; i < nua; i++) {// uas b3 >4
			register uint32_t U = tua[i] & BIT_SET_27;
			MISS1ADD;
		}
		if (locdiag) {
			cout << Char27out(wa) << "Gob1 aa wa  end after first adds nin="<<xq.nin << endl;
			//xq.Status();
		}
		if (!wa) return;
		if (xq.nout) {	GoB1toMiss0(wa); return; }
		g17b.Do7x();
		xq.BuildMissxRedundant();
		for (uint32_t i = 0; i <= nbbgm; i++) {
			GUM64& gw = tgm64[i];
			register uint64_t V = gw.Getv(g17b.tcluesxp, g17b.ncluesxp);
			register uint32_t r;
			while (bitscanforward64(r, V)) {
				V ^= (uint64_t)1 << r;
				register uint32_t U = gw.tb3[r];
				MISS1ADDIF;
			}
		}
		for (uint32_t i = 0; i <= nbbgmm; i++) {
			GUM64& gw = tgm64m[i];
			register uint64_t V = gw.Getv(g17b.tcluesxp, g17b.ncluesxp);
			register uint32_t r;
			while (bitscanforward64(r, V)) {
				V ^= (uint64_t)1 << r;
				register uint32_t U = gw.tb3[r];
				MISS1ADDIF;
			}
		}
		if (locdiag) {
			cout << Char27out(wa) << "Gob1 bb wa  end after adds nin=" << xq.nin << endl;
			xq.Status();
		}
		if (!wa) return;
		if (xq.nout) { GoB1toMiss0(wa); return; }
		//___ test global  
		if (g17b.VerifyValidb3())return;

		g17b.nt3more = 0;
		if (g17b.IsValidB3(xq.critbf)) {// not valid, new outfield
			register uint32_t U = g17b.anduab3;
			if (!U) return;
			if (locdiag) cout << Char27out(U) << "Not valid critbf go to miss0" << endl;
			GoB1toMiss0(U); return;
		}
		if (locdiag) 		cout  << "Gob1 try subcritical"  << endl;
		// now miss1 no out try all in
		XQ xqr = xq;
		g17b.TryMiss1Subcritical();
		if (locdiag) 		cout << "Gob1 go miss0 with all out" << endl;
		xq = xqr;
		uint32_t U = BIT_SET_27 ^ xq.critbf;// dummy UA for any outfield
		GoB1toMiss0(U);
	}
}
void STD_B3::GoB1toMiss0(uint32_t wa) {
	int locdiag = 0;
	//cout << " now GoB1toMiss0 miss 0" << endl;
	xq.t2b[xq.n2b++] = wa;
	xq.critbf |= wa;
	xq.nmiss = 0;
	xq.nout = 0;
	//xq.Status();
	xq.SetFilters();
	register uint32_t F = xq.fa, A = xq.fb;
	if (xq.nin) {// check xq.in potential 
		while (1) {// open the door tomore than one assign
			register uint32_t Fr = F, nn = xq.nin;
			xq.nin = 0;
			for (uint32_t i = 0; i < nn; i++) {
				register uint32_t U = xq.tin[i];
				MISS0ADD;
			}
			if (Fr == F) break;;
		}
	}
	//cout << "still valid end of GoB1toMiss0" << endl;
	//xq.Statusmiss0(F, A);
	if (_popcnt32(F) == xq.nb3) { GoC0F(F); return; }
	GoC0(F, A);
}
void G17B::TryMiss1Subcritical() {
	int locdiag = 0;
	if (p_cpt2g[8] == VTEST || g17b.knownt == 12)  locdiag = 1;
	if(locdiag)cout << "Subcritical() after buildout" << endl;
	memcpy(xq.tout, xq.t2a, xq.n2a * sizeof xq.tout[0]);
	memcpy(&xq.tout[xq.n2a], xq.t2b3, xq.n2b3 * sizeof xq.tout[0]);
	memcpy(&xq.tout[xq.n2a + xq.n2b3], xq.t2b, xq.n2b * sizeof xq.tout[0]);
	xq.nout = xq.n2a + xq.n2b3 + xq.n2b;// +nin;
	register uint32_t F = 0, A = xq.critbf;

	{
		uint32_t ass = 0;// load xq.in in the limit of A
		for (uint32_t i = 0; i < xq.nin; i++) {
			register uint32_t u = xq.tin[i];
			if (F & u) continue;
			if (!(u &= A)) return;// dead
			if (_popcnt32(u) == 1)ass |= u;
			else xq.tout[xq.nout++] = u;
		}
		F |= ass; xq.fa |= ass;
		A &= ~ass; xq.fb &= ~ass;
	}
	xq.nin = 0;
	if (locdiag) xq.Statusmiss0(F, A);
	if (!F) { 
		xq.CleanOut(F, A);
		if (locdiag) xq.Statusmiss0(F, A);
		GoEndAll(F, A,locdiag);
		return; 
	}
	int more = xq.SubCritMore(F);
	if (locdiag) cout << "nmore="<<more<<endl;
	if (more > 1) return;
	if (_popcnt32(F) == ncluesb3) {
		if (xq.NonHitOut(F))return;
	}


	if (more) {//now miss0
		xq.SubCritMoreDo(F);
		A = xq.SubCritActive(F);
		GoSubcritToMiss0(F, A);
		return;
	}
	xq.CleanOut(F, A);
	if (locdiag) {
		cout << "after cleanout" << endl;
		xq.Statusmiss0(F, A);
	}
	GoEndAll(F, A,locdiag);
}
void G17B::GoSubcritToMiss0(uint32_t bf, uint32_t ac) {
	int locdiag = 0;
	if (p_cpt2g[8] == VTEST || g17b.knownt == 12)  locdiag = 1;
	xq.nmiss = 0;
	//xq.CleanOut(bf, ac);
	if (locdiag) {
		cout << " GoSubcritToMiss0" << endl;
		xq.Status();
	}
	register uint32_t F = bf, A = 0;
	F |= xq.t1a;
	for (uint32_t i = 0; i < xq.n2b3; i++)
		if (!(F & xq.t2b3[i]))A |= xq.t2b3[i];

	for (uint32_t i = 0; i < xq.n2b; i++)
		if (!(F & xq.t2b[i]))A |= xq.t2b[i];
	if (locdiag) {
		cout << " GoSubcritToMiss0 while process tout" << endl;
		xq.Statusmiss0(F, A);
	}
	while (1) {// open the door tomore than one assign
		register uint32_t Fr = F, nn = xq.nout;
		xq.nout = 0;
		for (uint32_t i = 0; i < nn; i++) {
			register uint32_t U = xq.tout[i];
			if (F & U) continue;
			if (!(U &= A)) return;
			if (_popcnt32(U) == 1) {
				F |= U;
				A = xq.AssignMiss0(U);
			}
			else xq.Addout(U);
		}
		if (Fr == F) break;;
	}
	if (_popcnt32(F) == xq.nb3)  GoD0F(F);
	else GoEndAll(F, A);

}


//______________________ band3 missing >1 at start

void STD_B3::GoBMore1() {
	p_cpt2g[42]++;
	int locdiag = 0;
	if (op.known && g17b.knownt >= 10) locdiag = 1;
	if (VTEST && p_cpt2g[8] == VTEST) locdiag = 1;
	//xq.Status();
	register uint32_t  A = xq.critbf, wa = BIT_SET_27;
	xq.BuildMissxRedundant();
	{// add band 3 UA size 4
		for (uint32_t i = 0; i < nua; i++) {
			register uint32_t U = tua[i];
			if ((U >> 27) > 4) { xq.iuas4 = i; break; }//
			U &= BIT_SET_27;
			MISS1ADDIF;
		}
		for (uint32_t i = 0; i < ntg2_4; i++) {//guas 4
			int i81 = tg2_4[i];
			register uint32_t U = g.pat2[i81];
			MISS1ADDIF;
		}
		for (uint32_t i = 0; i < ntg2_6; i++) {//guas 6
			int i81 = tg2_6[i];
			register uint32_t U = g.pat2[i81];
			MISS1ADDIF;
		}
		for (uint32_t i = xq.iuas4; i < nua; i++) {// uas b3 >4
			register uint32_t U = tua[i] & BIT_SET_27;
			MISS1ADDIF;
		}
	}

	if(locdiag)	cout << Char27out(wa) << "GoBMore1 end after first adds nin="<<xq.nin << endl;
	g17b.Do7x();// be sure to hae the list of clues
	for (uint32_t i = 0; i <= nbbgm; i++) {
		GUM64& gw = tgm64[i];
		register uint64_t V = gw.Getv(g17b.tcluesxp, g17b.ncluesxp);
		register uint32_t r;
		while (bitscanforward64(r, V)) {
			V ^= (uint64_t)1 << r;
			register uint32_t U = gw.tb3[r];
			MISS1ADDIF;
		}
	}

	for (uint32_t i = 0; i <= nbbgmm; i++) {
		GUM64& gw = tgm64m[i];
		register uint64_t V = gw.Getv(g17b.tcluesxp, g17b.ncluesxp);
		register uint32_t r;
		while (bitscanforward64(r, V)) {
			V ^= (uint64_t)1 << r;
			register uint32_t U = gw.tb3[r];
			MISS1ADDIF;
		}
	}
	if (locdiag) {
		cout  << "GoBMore1 all loaded"  << endl;
		xq.Status();
	}
	if (!xq.nout) {//go direct
		xq.BuildMissxOut();
		xq.CleanOut();
		if (locdiag) {
			cout << "GoBMore1 buildout done" << endl;
			xq.Status();
		}
		g17b.GoEndAll(0, BIT_SET_27,locdiag);
		return;
	}
	if (xq.nmiss == 2) {// exit or miss1
		if (!wa) {// try miss0 or exit 
			//cout << "miss2 try miss0 or exit   [8] " << p_cpt2g[8] << endl;
			//xq.DumpOut();
//			if(xq.MissxToMiss0(xq.Isoutsize2())) return;
			uint64_t x = xq.Isoutsize2();
			if (x) {
				GoBMoretoMiss0(x);
				return;
			}
		}
		GoBMoretoMiss1(1);// first as disjoint
		return;
	}
	xq.CleanOut();
	if (locdiag) {
		cout << "GoBMore1 cleaned out" << endl;
		xq.Status();
	}
	if (xq.nout<=xq.nmiss-2 || wa) {//go direct
		xq.BuildMissxOut();
		xq.CleanOut();
		if (locdiag) {
			cout << "GoBMore1 all out" << endl;
			xq.Status();
		}
		g17b.GoEndAll(0, BIT_SET_27,locdiag);
		return;
	}
	if (xq.nmiss == 3) {
		uint64_t x = xq.Isoutsize2();
		if (!x) {//go direct
			xq.BuildMissxOut();
			xq.CleanOut();
			g17b.GoEndAll(0, BIT_SET_27);
			return;
		}
		uint64_t y = xq.Isoutsize3();
		if (!y) {// push to miss 1
			GoBMoretoMiss1(x);
			return;
		}
		else {
			GoBMoretoMiss0(y);
			return;
		}
	}
	if (xq.nmiss == 4) {
		uint64_t x = xq.Isoutsize3();
		if (locdiag) {
			cout << Char64out(x) << "x size 3=" << endl;
		}
		if (!x) {//go direct
			xq.BuildMissxOut();
			g17b.GoEndAll(0, BIT_SET_27);
			return;
		}
		uint64_t y = xq.Isoutsize4();
		if (locdiag)	cout << Char64out(y) << "y size 4=" << endl;
		if (!y) {// push to miss 1
			if (locdiag) {
				cout << "<miss4 + size3 push to miss1" << endl;
				for (uint32_t i = 0, bit = 1; i < xq.nout; i++, bit <<= 1) {
					register uint32_t u = xq.tout[i];
					if (bit & x) {// new crit
						cout << Char2Xout(xq.tout[i]) << " bit " << i << endl;
					}
				}
			}
			GoBMoretoMiss1(x);
			return;
		}
		else {
			GoBMoretoMiss0(y);
			return;
		}
	}
	//cout << "miss>2  [8] " << p_cpt2g[8] << " miss= "<< xq.nmiss << endl;
	//xq.DumpOut();
	xq.BuildMissxOut();
	g17b.GoEndAll(0, BIT_SET_27);

}
void STD_B3::GoBMoretoMiss0(uint64_t ubf) {
	int locdiag = 0;
	if (!ubf) return ;
	for (uint64_t i = 0, bit = 1; i < xq.nout; i++, bit <<= 1) {
		register uint32_t u = xq.tout[i];
		if (bit & ubf) {// new crit
			xq.t2b[xq.n2b++] = u;
			xq.critbf |= u;
		}
		else xq.tin[xq.nin++] = u; // put it in to check later
	}
	xq.nmiss = xq.nout = 0;
	//xq.Status();
	xq.SetFilters();
	register uint32_t F = xq.fa, A = xq.fb;
	if (xq.nin) {// check xq.in potential 
		while (1) {// open the door tomore than one assign
			register uint32_t Fr = F, nn = xq.nin;
			xq.nin = 0;
			for (uint32_t i = 0; i < nn; i++) {
				register uint32_t U = xq.tin[i];
				MISS0ADD;
			}
			if (Fr == F) break;;
		}
	}
	//cout << "still valid end of GoB1toMiss0" << endl;
	//xq.Statusmiss0(F, A);
	if (_popcnt32(F) == xq.nb3) { GoC0F(F); return; }
	GoC0(F, A);
}
void STD_B3::GoBMoretoMiss1(uint64_t ubf) {
	if (g17b.VerifyValidb3())return;
	int locdiag = 0;
	if (op.known && g17b.knownt >= 10) locdiag = 1;
	if (VTEST && p_cpt2g[8] == VTEST) locdiag = 1;
	if (!ubf) return;
	for (uint64_t i = 0, bit = 1; i < xq.nout; i++, bit <<= 1) {
		register uint32_t u = xq.tout[i];
		if (bit & ubf) {// new crit
			xq.t2b[xq.n2b++] = u;
			xq.critbf |= u;
		}
		else xq.tin[xq.nin++] = u; // put it in to check later
	}
	xq.nmiss = 1;
	xq.nout = 0;
	if (locdiag) xq.Status();
	g17b.nt3more = 0;
	if (g17b.IsValidB3(xq.critbf)) {// not valid, new outfield
		register uint32_t U = g17b.anduab3;
		if (!U) return;
		if (locdiag)cout << Char27out(U) << "Not valid critbf go to miss0" << endl;
		GoB1toMiss0(U);
		return;
	}
	// now miss1 no out try all in
	XQ xqr = xq;
	if (locdiag)cout  << "valid critbf continue miss1 mode " << endl;
	g17b.TryMiss1Subcritical();
	xq = xqr;
	if (locdiag)cout << "miss1 mode  go to miss0" << endl;
	uint32_t U = BIT_SET_27 ^ xq.critbf;// dummy UA for any outfield
	GoB1toMiss0(U);
	return;

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


//____________processing band 3 to the end
void OutFreshUa(uint64_t u,uint32_t u3, const char * lib) {
	if (op.known)if (g17b.pk54 &u) return;
	cout << Char54out(u) << "\t";
	cout << Char27out(u3) << lib
		<< "   [3]" << p_cpt2g[3] << "   [4]" << p_cpt2g[4]
		<< "[7]" << p_cpt2g[7] << "   [8]" << p_cpt2g[8]
		<< endl;

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
	if (debug)  
		cout <<Char27out(bf) << "IsValidB3 debugging mode [8]"<< p_cpt2g[8] << endl;
	 
	int nret;
	if ((nret=zhou[0].CallCheckB3( bf,(debug>1)))) {
		anduab3 = BIT_SET_27;
		for (int iadd = 0; iadd < nret; iadd++) {
			BF128 w = zhgxn.tua[iadd];
			if (debug) {
				cout << Char2Xout(w.bf.u64[0]) << " ";
				cout << Char27out(w.bf.u32[2]) << "fresh ua" << endl;
			}
			int cc = _popcnt32(w.bf.u32[2]);
			if (!cc) {
				cout << Char27out(bf) << " bug no b3 IsValidB3 [7]"					
					<< p_cpt2g[7]	<< "   [8]" << p_cpt2g[8] 
					<< "   [3]" << p_cpt2g[3] << "   [4]" << p_cpt2g[4]
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
				cout << Char27out(bf) << " bug ua nob12 [8]" << p_cpt2g[8] 
					<< "   [3]" << p_cpt2g[3] << "   [4]" << p_cpt2g[4] 
					<<" [7] "	<< p_cpt2g[7] << endl;
				cout << Char54out(myb12) << " " << endl;
				xq.Status();
				cout << "list of uas found" << endl;
				for (uint32_t i = 0; i < zhgxn.nua; i++) {
					BF128 ww = zhgxn.tua[i];
					//cout << Char2Xout(ww.bf.u64[0]) << " ";
					cout << Char27out(ww.bf.u32[2]) << " i=" << i << endl;
				}
				//zhou[0].CallCheckB3(bf, 1);		
				aigstop = 1;	return 1;
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
							if (op.dv3)OutFreshUa(U, w.bf.u32[2], " added 0A");
							guah54n.AddA2(U, i81,(int)cc0);
							guah54n.g2.setBit(i81);
							p_cpt2g[20]++;
							continue;
						}
					}
					p_cpt2g[19]++;
					if (cc == 4) {
						myband3->Addm4(w);
						if (op.dv3)OutFreshUa(U, w.bf.u32[2], " added 4");
					}
					else {
						myband3->Addmm(w);
						if (op.dv3)OutFreshUa(U, w.bf.u32[2], " added >4");
					}
					continue;
				}

				else {
					if (cc == 2) {
						int i81 = myband3->GetI81_2(w.bf.u32[2]);
						if (op.dv3) OutFreshUa(U, w.bf.u32[2], " added ob");
						guah54n.AddA2(U, i81,(int)cc0);
						guah54n.g2.setBit(i81);
						p_cpt2g[20]++;
					}
					else  {
						int i81 = myband3->GetI81_3(w.bf.u32[2]);						
						guah54n.AddA3(U, i81);
						guah54n.g3.setBit(i81);
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
	nt3more = 0;
	if (!IsValidB3(bf))return 0;
	for (uint32_t i = 0; i < nt3more; i++) {
		xq.tout[xq.nout++] = t3more[i];
	}
	nt3more = 0;
	return 1;
}
int  G17B::Valid3mm(uint32_t bf) {
	nt3more = 0;
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

void G17B::GoCallB3_12(CALLBAND3& cb3w) {
	/*
	return; // not yet adjusted to gua54hn
	p_cpt2g[7]++;
	if (sgo.bfx[3] & 1) return; // only ohase 1

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

	//guah54n.GetG2G3_12(cb3w.g2t, cb3w.g3t);
	if (sgo.bfx[3] & 2) return; // stop here

	if (locdiag) {
		guah54n.g2.Print("g2");
	}
	//return;
	if (locdiag) {
		cout << "g2t g3t done" << endl; //return;
		if (op.ton > 1) {
			guah54n.g2.Print("g2");
			//guah54.DumpA2(); guah54.DumpB2(1);
		}
	}
	{
		register uint64_t F = g17b.myb12 & ~g17b.bf_cl9; // cells 6 to x
		ncluesxp = 0;
		register int cell;// build table of cells
		while (bitscanforward64(cell, F)) {
			F ^= (uint64_t)1 << cell;
			tcluesxp[ncluesxp++] = cell;
		}
	}
	//p_cpt2g[16] += cb3.g2t.Count(); p_cpt2g[18] += cb3.g3t.Count();
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& b3 = genb12.bands3[ib3];
		b3.Go(cb3w);
	}
	if (op.known && knownt == 12) knownt = 13;	*/

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
uint64_t XQ::Isoutsize2() {
	if (nout < 2) return 0;
	register uint32_t n = nout;
	if (n > 64) n = 64;
	for (uint32_t i1 = 0; i1 < n - 1; i1++) {
		register uint32_t U1 = tout[i1];
		for (uint32_t i2 = i1+1; i2 < n; i2++) {
			register uint32_t U2 = tout[i2];
			if (U1 & U2) continue;
			return ((uint64_t)1 << i1) | ((uint64_t)1 << i2);
		}
	}
	return 0;
}
uint64_t XQ::Isoutsize3() {
	if (nout < 3) return 0;
	register uint32_t n = nout;
	if (n > 64) n = 64;
	for (uint32_t i1 = 0; i1 < n - 2; i1++) {
		register uint32_t U1 = tout[i1];
		for (uint32_t i2 = i1 + 1; i2 < n-1; i2++) {
			register uint32_t U2 = tout[i2];
			if (U1 & U2) continue;
			for (uint32_t i3 = i2 + 1; i3 < n; i3++) {
				register uint32_t U3 = tout[i3];
				if ((U1 | U2) & U3) continue;
				return ((uint64_t)1 << i1) | ((uint64_t)1 << i2) | ((uint64_t)1 << i3);
			}
		}
	}
	return 0;
}
uint64_t XQ::Isoutsize4() {
	if (nout < 4) return 0;
	register uint32_t n = nout;  
	if (n > 64) n = 64;
	for (uint32_t i1 = 0; i1 < n-3; i1++) {
		register uint32_t U1 = tout[i1];
		for (uint32_t i2 = i1 + 1; i2 < n-2; i2++) {
			register uint32_t U2 = tout[i2];
			if (U1 & U2) continue;
			for (uint32_t i3 = i2 + 1; i3 < n-1; i3++) {
				register uint32_t U3 = tout[i3];
				if ((U1 | U2) & U3) continue;
				register uint32_t uo3 = U1 | U2 |U3;
				for (uint32_t i4 = i3 + 1; i4 < n ; i4++) {
					register uint32_t U4 = tout[i4];
					if (uo3 & U4) continue;
					return ((uint64_t)1 << i1) | ((uint64_t)1 << i2) |
						((uint64_t)1 << i3) | ((uint64_t)1 << i4);
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
	uint32_t tx[8][100], ntx[7];
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
void XQ::CleanOut(uint32_t F, uint32_t A) {// usually miss1 last step
	uint32_t nn = nout; nout = 0;
	uint32_t tx[7][100], ntx[7];
	memset(ntx, 0, sizeof ntx);
	for (uint32_t i = 0; i < nn; i++) {
		register uint32_t u = tout[i];
		if (u & F) continue;
		u &= A;
		register uint32_t cc = _popcnt32(u);
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


void  G17B::GoEndAll(uint32_t bf, uint32_t ac, int debug) {//  call the relevant expand b3
	ac &= BIT_SET_27;// be sure to have no extra digit
	int ass = _popcnt32(bf);
	ntoassb3 = xq.nb3-ass;
	if (VerifyValidb3())return;
	if (ntoassb3 == 1) {
		int x,u;
		if (!xq.nout) {
			if (!IsValidB3(bf))u = ac;
			else u= anduab3;
		}
		else // must try a possible "17"belox"
			if (!(u = xq.GetAndout())) return; // no single clue
		while (bitscanforward(x, u)) {
			register int bit= 1 << x;
			u ^= bit;
			if (IsValidB3(bf | bit)) u &= anduab3;
			else Out17(bf | bit);
		}
		return;
	}
	if (ntoassb3 >= 4) {
		GoB3Expand_1_3(bf,ac,debug);
		return;
	}
	BuildSortT3b();
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
	//if (p_cpt2g[8] == 3475969) locdiag = 1;
	if (locdiag) {
		cout << "expand entry 1-3 nto ass =" << ntoassb3
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
		if (VerifyValidb3()) return;
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
		if (VerifyValidb3()) return;
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
		if (locdiag) cout << Char27out(spb[3].all_previous) << " next3 all cells " << endl;

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
		if (VerifyValidb3()) return;// not sure it is done
		if (!n) {// could be a valid
			if (Valid3_1_3(spb[3].all_previous))
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
	if (locdiag) {
		cout << " endnext3 nout=" << xq.nout << endl;
		xq.DumpOut();
	}

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
		if (!nt3b) {// still un assigned >=2 can not be valid
			if (VerifyValidb3()) return;
			if (!Valid3_1_3(spb[3].all_previous))goto next3;
			for (uint32_t i = 0; i < zhgxn.nua; i++) {
				BF128 ww = zhgxn.tua[i];
				t3b[nt3b++] = ww.bf.u32[2];
			}
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
		cout <<Char27out(spe.all_previous) << " expand entry 4_x nto ass =" << ntoass
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
		if (VerifyValidb3())  return; 
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
		if (VerifyValidb3())  return;
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
			//uint32_t nrt3 = nt3more;// don't touch old ??
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
			<< " [4] " << p_cpt2g[4] << " [7] " << p_cpt2g[7] << " [8] " << p_cpt2g[8] << endl;

		DumpPotential();
		aigstop = 1;
		return;
	}
	if (op.out_one) {
		if (myband3->aigskip) return;
		myband3->aigskip = 1;
		nb3_not_found--;
		if (nb3_not_found <= 0) aigstop = 1;
	}
	myband3->poutdone++;
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
			<< " [4] " << p_cpt2g[4] << " [7] " << p_cpt2g[7]
			 << " [8] " << p_cpt2g[8] << endl;

	}
	sprintf(&ws[81], ";%3d;%3d;%3d;%5d", genb12.i1t16, genb12.i2t16, t416_to_n6[myband3->i416],(int)( genb12.nb12 >> 6));
	fout1 << ws << endl;
	
	a_17_found_here++;

}