struct CBS {// clues band stack
	uint16_t b[3];
	uint16_t s[3];
	inline void Add(uint32_t cell) {
		b[cell / 27]++;
		s[C_stack[cell]]++;
	}
	inline int IsFilt10_17() {// 737 377 is possible
		if (b[0] > 7 || b[1] > 7)return 1;
		return 0;
	}
	inline uint64_t LimBand6(uint64_t bf) {
		if (b[0] >6 || b[1] > 67)return 0;
		if (b[0] == 6) return bf & (~BIT_SET_27);
		if (b[1] == 6) return bf & BIT_SET_27;
		return bf;
	}

	inline int IsFilt11_17() {// 566 656 no stack>6
		if (b[0] > 6 || b[1] > 6)return 1;
		if (s[0] > 6 || s[1] > 6 || s[2] > 6)return 1;
		return 0;
	}

	inline int IsFilt11_18() {
		if (b[0] > 7 || b[1] > 6)return 1;
		if (s[0] > 7 || s[1] > 7 || s[2] > 7)return 1;
		return 0;
	}
	inline int IsFilt12_18() {
		if (b[0] != 6)return 1;
		if (s[0] > 6 || s[1] > 6 || s[2] > 6)return 1;
		return 0;
	}
};
struct SPB03 {// spots to first 7 clues
	BF128 v;
	uint64_t  possible_cells, all_previous_cells, active_cells;
	CBS cbs;
	uint32_t ncl;
	void Dump(uint64_t x) {
		cout << "spb03 status ncl=" << ncl << " " << x << endl;
		cout << Char54out(all_previous_cells) << " assigned" << endl;
		cout << Char54out(active_cells) << " active" << endl;
		cout << Char54out(possible_cells) << " possible" << endl;
		cout << Char64out(v.bf.u64[0]) << " 64 v" << endl;
	}
}spb_0_15[16]; 
struct CALLBAND3 {
	BF128 g2t, g3t;
	uint64_t bf12,bfcom;
	uint32_t ncl;
	CBS cbs; 
	void Dump() {
		cout << Char54out(bf12) << " b12 to process ncl=" <<ncl
			<<" "<<cbs.b[0]<<cbs.b[1]<<" "
			<<cbs.s[0] << cbs.s[1] << cbs.s[2] << endl;
		cout << Char64out(g2t.bf.u64[0]);
		cout << Char27out(g2t.bf.u32[2]) << " g2t" << endl;;
		cout << Char64out(g3t.bf.u64[0]);
		cout << Char27out(g3t.bf.u32[2]) << " g3t" << endl;;

	}
}cb3;

struct SGUA2 {// 81 possible UA2 sockets
	// permanent data
	//uint64_t* tua;
	int col1, col2;// columns of the socket
	int i_81; // index 0_80 for this 
	int i9;// complementary column in minirow
	int id1, id2; // index of digits in gang 27 
	// Current band1+2 data
	int digs, dig1, dig2;// depending on gang27 status
	int valid, // valid if guas 
		validuas,// gua2s found
		used;// if needed in bands3
	int gangcols[9];// revised gangster
	//uint32_t nua;// nua_start, nua_end;

}tsgua2[81];
struct SGUA3 {// 81 possible UA3 sockets
	// permanent data
	//uint64_t* tua;// , killer;
	int col1;// first columns 0-9 
	int i_81, stack;// , iguan; // index 0_80 for this 
	int id1, id2, id3; // index of digits in gang 27 
	// Current band1+2 data
	int  dig1, dig2, dig3, digs;// depending on gang27 status
	int valid, // valid if guas 
		validuas,// gua2s found
		used;// if needed in bands3
	//uint32_t nua;// nua_start, nua_end, nua;
}tsgua3[81];
#define UA12SIZE 3840
#define UA12BLOCS 30

struct TUASB12 {//  initial set of uas bands 1+2 
	uint64_t tua[UA12SIZE]; // 30x128
	uint32_t nua,  tdigs[UA12SIZE],ndigs[UA12SIZE];
	inline void AddInit(
		int64_t ua, uint32_t digs,  uint32_t nd) {
		tdigs[nua] = digs;
		ndigs[nua] = nd;
		tua[nua++] = ua;
	}
	void DumpInit() {
		cout << "dumpinit  tua nua=" << nua << endl;
		for (uint32_t iua = 0; iua < nua; iua++) {
			cout << Char2Xout(tua[iua]) << " i="
				<< iua << " " << (tua[iua] >> 59)
				<< " " <<_popcnt64(tua[iua] & BIT_SET_2X);
			cout << " digs " << Char9out(tdigs[iua]) << " " << ndigs[iua] << endl;;
		}
	}
	void DumpShort(const char * lib) {
		cout << lib  << "tuasb12 tua nua=" << nua << endl;
	}

}tuasb12;
struct T54B12 {//   uas bands 1+2 in 54 mode
	struct TUVECT {//  128 uas and vectors
		BF128 v0, vc[54];
		uint64_t t[128];
		void Init() {
			v0.SetAll_0();
			memset(vc, 255, sizeof vc);
		}
		void Dumpv(int lim = 32) {
			cout << Char32out(v0.bf.u32[0]) << " v0" << endl;
			for (int i = 0; i < 54; i++)
				cout << Char32out(vc[i].bf.u32[0]) << " cell=" << i << endl;
		}
		void Dump(int lim = 128) {
			for (int i = 0; i < lim; i++)
				if (v0.On(i))
					cout << Char54out(t[i]) << " " << i <<" "<<_popcnt64(t[i]) << endl;
				else return;
		}
	};

	// initial status after harvest plus fresh uas B12
	TUVECT ta128[UA12BLOCS];// max start 20*128=2560 uas 
	uint32_t na128, nablocs, nta128[UA12BLOCS];
	void InitA() {
		memset(nta128, 0, sizeof nta128);
		na128 = nablocs = 0;
		for (int i = 0; i < UA12BLOCS; i++)ta128[i].Init();
	}

	inline void AddA(uint64_t u) {
		if (na128 >= UA12BLOCS * 128) return;
		register uint32_t bloc = na128 >> 7, ir = na128 - (bloc << 7);
		na128++; nablocs = bloc; nta128[bloc]++;
		ta128[bloc].v0.setBit(ir);
		BF128* myvc = ta128[bloc].vc;
		register uint64_t R = u;
		ta128[bloc].t[ir] = R;
		R &= BIT_SET_54;// clear extra bits
		uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvc[cell].clearBit(ir);
		}
	}
	void Build_ta128(uint64_t* t, uint32_t n);

#define UABNBLOCS 20
	// status after 3 clues (B)
	TUVECT tb128[UABNBLOCS];// max start 10*128=1280 uas 
	uint32_t nb128, nbblocs, ntb128[UABNBLOCS];
	void InitB() {
		memset(ntb128, 0, sizeof ntb128);
		nb128 = nbblocs = 0;
		for (int i = 0; i < UABNBLOCS; i++)tb128[i].Init();
	}
	inline void AddB(uint64_t u) {
		if (nb128 >= UABNBLOCS*128) return;
		register uint32_t bloc = nb128 >> 7,
			ir = nb128 - 128 * bloc;
		nb128++; nbblocs = bloc; ntb128[bloc]++;
		tb128[bloc].v0.setBit(ir);
		BF128* myvc = tb128[bloc].vc;
		register uint64_t R = u;
		tb128[bloc].t[ir] = R;
		register uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvc[cell].clearBit(ir);
		}
	}
	int Build_tb128();

#define UACNBLOCS 15
	// status after 6 clues (C)
	TUVECT tc128[UACNBLOCS];// max start 10*128=1280 uas 
	uint32_t nc128, ncblocs, ntc128[UACNBLOCS];
	void InitC() {
		memset(ntc128, 0, sizeof ntc128);
		nc128 = ncblocs = 0;
		for (int i = 0; i < UACNBLOCS; i++)tc128[i].Init();
	}
	inline void AddC(uint64_t u) {
		if (nc128 >= UACNBLOCS * 128) return;
		register uint32_t bloc = nc128 >> 7,
			ir = nc128 - 128 * bloc;
		nc128++; ncblocs = bloc; ntc128[bloc]++;
		tc128[bloc].v0.setBit(ir);
		BF128* myvc = tc128[bloc].vc;
		register uint64_t R = u;
		tc128[bloc].t[ir] = R;
		register uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvc[cell].clearBit(ir);
		}
	}
	int Build_tc128();// after 6 clues
	int Build_tc128_7();//after 7 clues 
	inline int IsNotRedundant(uint64_t u) {
		register uint64_t nu = ~u;
		for (uint32_t i = 0; i < ntc128[0]; i++)
			if (!(tc128[0].t[i] & nu)) return 0;
		return 1;
	}

	// more in a chunk expand
	uint64_t tm[200], ntm;
	void AddM(uint64_t u) {
		if (ntm < 200)tm[ntm++] = u;
	}
	inline int NotValid(uint64_t u) {
		for (uint64_t i = 0; i < ntm; i++)
			if (!(u & tm[i])) return 1;
		return 0;
	}

	void DebugA() {
		cout << " debugA na128=" << na128<<" nablocs=" <<nablocs<< endl;
		for (uint32_t i = 0; i <= nablocs; i++) {
			cout << "+ " << 128 * i<< " count " <<nta128[i]  << endl;
			ta128[i].Dump();
		}
	}

	void DebugB() {
		cout << " debugB nb128=" << nb128 << " nbblocs=" << nbblocs << endl;
		for (uint32_t i = 0; i <= nbblocs; i++) {
			cout << "+ " << 128 * i << endl;
			tb128[i].Dump();
		}
	}
	void DebugC() {
		cout << " debugC nc128=" << nc128 << endl;
		tc128[0].Dump();
	}

}t54b12;

struct T54G2 {//
	struct G2VECT {//  128 guas2 and vectors  
		BF128 v0, vc[54], v81[81];
		//vectors 128 {base, cells, i81s}
		uint64_t t[128];// ua54 + i81
		void Init() {
			v0.SetAll_0();
			memset(vc, 255, sizeof vc);
			memset(v81, 255, sizeof v81);
		}	

	};
#define NG2BLOCS6 30
	G2VECT t128[NG2BLOCS6];// max start 20*128=2560 uas 
	BF128 vcl[7][NG2BLOCS6];//6 clues  

	uint32_t n128, nblocs, nt128[NG2BLOCS6];
	BF128 g2;
	void Init() {
		memset(nt128, 0, sizeof nt128);
		n128 = nblocs = 0;
		for (int i = 0; i < NG2BLOCS6; i++)t128[i].Init();
	}
	inline void Add(uint64_t u) {
		if (n128 >= NG2BLOCS6 * 128) return;
		register uint32_t bloc = n128 >> 7, ir = n128 - (bloc << 7);
		n128++; nblocs = bloc; nt128[bloc]++;
		register G2VECT& myb = t128[bloc];
		myb.v0.setBit(ir);
		BF128* myvc = myb.vc,*myvi=myb.v81;
		register uint64_t R = u;
		myb.t[ir] = R;
		myvi[R>>56].clearBit(ir);
		R &= BIT_SET_54;// clear i81
		uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvc[cell].clearBit(ir);
		}
	}
	void BuildG2();
	inline void NewVcl(int ncl, int clue) {
		register int i1 = ncl - 7;
		register BF128* V1 = vcl[i1], * V2 = vcl[i1 + 1];
		for (uint32_t i = 0; i <= nblocs; i++)
			V2[i] = V1[i] & t128[i].vc[clue];
	}
	void GetActive(int nclues) {
		BF128* V2 = vcl[nclues - 6];
		g2.SetAll_0();
		for (uint32_t ib = 0; ib < nblocs; ib++) {
			BF128 w128 = V2[ib];
			if (w128.isNotEmpty()) {
				G2VECT& gv = t128[ib];
				int x;
				while ((x = w128.getFirst128()) >= 0) {
					register uint64_t U = gv.t[x],
						i81=U>>56;
					// Apply i81 here and downstream
					w128 &= gv.v81[i81];
					for (uint32_t ib2 = ib+1; ib2 < nblocs; ib2++)
						V2[ib2]&= t128[ib2].v81[i81];
					g2.setBit((int)i81);
				}
			}
		}
	}

	void Add54b3(uint64_t u) {

	}
}t54g2;
struct T54G3 {//
	struct G3VECT {//  128 guas2 and vectors 
		BF128 v0, vc[54], v81[81];
		uint64_t t[128];// ua54 + i81
		void Init() {
			v0.SetAll_0();
			memset(vc, 255, sizeof vc);
			memset(v81, 255, sizeof v81);
		}
	};

#define NG3BLOCS6 15
	// initial status after harvest plus fresh guas
	G3VECT t128[NG3BLOCS6];// max start 15*128=1920 uas 
	BF128 vcl[7][NG3BLOCS6];//6 clues  
	uint32_t n128, nblocs, nt128[NG3BLOCS6];
	BF128 g3;
	void Init() {
		memset(nt128, 0, sizeof nt128);
		n128 = nblocs = 0;
		for (int i = 0; i < NG3BLOCS6; i++)t128[i].Init();
	}
	inline void Add(uint64_t u) {
		if (n128 >= NG3BLOCS6 * 128) return;
		register uint32_t bloc = n128 >> 7, ir = n128 - (bloc << 7);
		n128++; nblocs = bloc; nt128[bloc]++;
		register G3VECT& myb = t128[bloc];
		myb.v0.setBit(ir);
		BF128* myvc = myb.vc, * myvi = myb.v81;
		register uint64_t R = u;
		myb.t[ir] = R;
		myvi[R >> 56].clearBit(ir);
		R &= BIT_SET_54;// clear i81
		uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvc[cell].clearBit(ir);
		}
	}
	void BuildG3();
	inline void NewVcl(int ncl, int clue) {
		register int i1 = ncl - 7;
		register BF128* V1 = vcl[i1], * V2 = vcl[i1 + 1];
		for (uint32_t i = 0; i <= nblocs; i++)
			V2[i] = V1[i] & t128[i].vc[clue];
	}
	void GetActive(int nclues) {
		BF128* V2 = vcl[nclues - 6];
		g3.SetAll_0();
		for (uint32_t ib = 0; ib < nblocs; ib++) {
			BF128 w128 = V2[ib];
			if (w128.isNotEmpty()) {
				G3VECT& gv = t128[ib];
				int x;
				while ((x = w128.getFirst128()) >= 0) {
					register uint64_t U = gv.t[x],
						i81 = U >> 56;
					// Apply i81 here and downstream
					w128 &= gv.v81[i81];
					for (uint32_t ib2 = ib + 1; ib2 < nblocs; ib2++)
						V2[ib2] &= t128[ib2].v81[i81];
					g3.setBit((int)i81);
				}
			}
		}
	}
	void Add54b3(uint64_t u) {

	}
}t54g3;

struct GUA54 {
	uint64_t* tua, killer;
	uint32_t nua, nuamax, type, i81;
	inline void Init(uint64_t* p, uint32_t t, uint32_t i) {
		tua = p; type = t; i81 = i; killer = ~0;
		nua = nuamax = 0;
	}
	inline void Add(uint64_t u) {
		if (nua >= nuamax) return;
		killer &= u;	tua[nua++] = u;
	}
	inline void AddCheck(uint64_t u) {// no redundancy
		register uint64_t nU = ~u;
		for (uint32_t j = 0; j < nua; j++)
			if (!(tua[j] & nU)) return;
		killer &= u;	tua[nua++] = u;
	}

	void Debug(int nodet = 1) {
		if (!nua)return;
		cout << "gua54 type=" << type << " i81=" << i81
			<< "\tnua=" << nua << "\tnuamax=" << nuamax << endl;
		cout << Char54out(killer) << " K" << endl;
		if (nodet) return;
		for (uint32_t i = 0; i < nua; i++) {
			cout << Char54out(tua[i]) << " " << _popcnt64(tua[i])
				<< " " << i << endl;
		}
	}
};
struct GUAH54 {// handler guas 2 3 in 54 mode
	uint64_t gbuf[162 * 60]; // room for cut 30+10 in average
	GUA54 tg2[81], tg3[81];

	void Build();
	void Build2(uint64_t filter, uint64_t active);
	BF128 GetG2(uint64_t bf);
	BF128 GetG3(uint64_t bf);
	void Add2(uint64_t bf, int i81) { tg2[i81].Add(bf); }
	void Add3(uint64_t bf, int i81) { tg3[i81].Add(bf); }
	void Dumpall2() {
		for (int i = 0; i < 81; i++) {
			tg2[i].Debug(0);
		}
	}
	void DumpOne2(int i81) {
		tg2[i81].Debug(0);
	}
	void Dumpall3() {
		for (int i = 0; i < 81; i++) {
			tg3[i].Debug(0);
		}
	}
}guah54, guah54_2;

struct CHUNK1B {//storing 64 uas and vectors band 3s
	uint64_t v0, vc[27], nt;
	uint32_t tua[64];
	inline void Init() {
		nt = 0, v0 = 0;
		memset(vc, 255, sizeof vc);
	}
	inline void Add(uint32_t v) {//add a new ua
		if (nt >= 64) return;// safety should never be

		{//__ check redundancy
			register uint32_t vn = ~v;
			for (uint64_t i = 0; i < nt; i++)
				if (!(tua[i] & vn))return; // == or subset
		}

		uint64_t bit = (uint64_t)1 << nt;
		v &= BIT_SET_27;//no extra bit
		tua[nt++] = v;
		v0 |= bit;
		uint32_t cc27;// build cells vectors
		register  uint32_t Rw = v;
		while (bitscanforward(cc27, Rw)) {
			Rw ^= 1 << cc27;// clear bit
			vc[cc27] ^= bit;
		}
	}


	inline uint64_t ApplyXY(uint32_t* tcells, uint32_t ntcells) {
		if (!nt) return 0;
		uint64_t w = v0;
		for (uint32_t i = 0; i < ntcells; i++)
			w &= vc[tcells[i]];
		return w;
	}
	void Debug(int all = 1) {
		for (uint32_t i = 0; i < nt; i++) {
			cout << Char27out(tua[i]) << " i=" << i << endl;
		}
		if (!all) return;
		cout << Char64out(v0) << " v0" << endl;
		for (int i = 0; i < 27; i++)
			if ((vc[i] & v0) != v0)
				cout << Char64out(vc[i] & v0) << " cell=" << i << endl;
	}
}b3direct;


// standard first band (or base any band band)
struct STD_B416 {
	char band[28];// band in char mode
	int i416,// id 0-415 in minlex mode of the band
		map[27],// mapping cells from minlex to solution grid
		band0[27], // band in 0-8 integer mode
		gangster[9], // bit field digits per column
		dpat[9],// solution pat per digit
		dpmi[9],// initial map gangster per digit
		dband;
	uint32_t tua[82], nua;//  maximum 81 morphed uas 
	uint32_t fd_sols[2][9];//start puzzle/ solution
	void Initstd();// first band initial (once)
	void GetBandTable(int i416e);//pick up band 1 from table
	void InitG12(int i416e);// end init band 1 
	void InitBand2_3(int i16, char* ze, BANDMINLEX::PERM& p
		, int iband = 1);

	void SetGangster();
	inline void GetUAs() {
		nua = t16_nua[i416];
		memcpy(tua, &t16_UAs[t16_indua[i416]], 4 * nua);
	}
	void MorphUas()	;
	void InitC10(int i);// known mode


	void PrintStatus();
};
struct STD_B3 :STD_B416 {// data specific to bands 3
	// permanent gangster information
	struct G {
		BF128 gsocket2, gsocket3;// active i81 mode 81 bits
		int pat2[81], pat3[81]; // storing ua bitfields
		int ua2_imini[81], ua3_imini[81],	ua2bit27[81];
	}g;
	int minirows_bf[9];
	int triplet_perms[9][2][3];
	uint32_t i_27_to_81[27], i_9_to_81[9]; //band i81 for guas guas3
	uint32_t i_81_to_27[81]; //band i81 for guas guas3

	struct GUAM {
		uint64_t bf12;// mode 54
		uint32_t bf3;
	}tguam[300],tguam2[200];
	uint32_t ntguam,ntguam2,guam2done;
	//_______________________

	void InitBand3(int i16, char * ze, BANDMINLEX::PERM & p);
	void Go(CALLBAND3& cb3);
	inline void BuildGuam2(uint64_t known) {
		register uint64_t F = known, n = 0;
		for (uint32_t i = 0; i < ntguam; i++)
			if (!(F & tguam[i].bf12))
				if(n<200)tguam2[n++] = tguam[i];
		ntguam2 = (int)n;
		guam2done = 1;
	}
	void Addguam(BF128 w) {// entry mode 3x32
		if (ntguam >= 300) return;
		tguam[ntguam].bf3 = w.bf.u32[2];
		register uint64_t U = w.bf.u64[0];
		U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
		tguam[ntguam++].bf12 = U;
	}

	uint32_t Get2d(int d1, int d2) {
		return fd_sols[0][d1] | fd_sols[0][d2];
	}
	inline int GetI81_2(int bf) {
		for (uint32_t i = 0; i < 27; i++) {
			register uint32_t i81 = i_27_to_81[i];
			if (g.pat2[i81] == bf) return i81;
		}
		return -1;
	}
	inline int GetI81_3(int bf) {
		for (uint32_t i = 0; i < 9; i++) {
			register uint32_t i81 = i_9_to_81[i];
			if (g.pat3[i81] == bf) return i81;
		}
		return -1;
	}
	void Pat_to_digitsx(int bf,int* tcells, int* tdigits, int& nt) {
		register int cell;
		nt = 0;
		while (bitscanforward(cell, bf)) {
			bf ^= 1 << cell;
			tcells[nt] = cell;
			tdigits[nt++] = band0[cell];
		}
	}
	int Is_Pat_For_Mex(int* tcells, int* tdigits, int nt) {
		for (int i = 0; i < nt; i++) {
			if (band0[tcells[nt]] != tdigits[nt]) return 0;
		}
		return 1;
	}
	void DumpGuam2d() {// only if 2 digits (check socket
		for (uint32_t i = 0; i < ntguam; i++) {
			GUAM w = tguam[i];
			register uint32_t C = w.bf3;
			uint32_t cc = _popcnt32(C),
				ncol = _popcnt32((C | (C >> 9) | (C >> 18)) & 0777),
				nsock = 2 * ncol - cc;
			cout << Char54out(w.bf12) << " ";
			cout << Char27out(w.bf3) << "  i=" << i << " nsock=" <<nsock 
				<< " ncells="<< cc  << endl;

		}
	}
	void DumpGuam(int det=0) {// only if 2 digits (check socket
		cout << "tguam n =" << ntguam << endl;
		if (!det) return;
		for (uint32_t i = 0; i < ntguam; i++) {
			GUAM w = tguam[i];
			cout << Char54out(w.bf12) << " ";
			cout << Char27out(w.bf3) << "  i=" << i
				<< " " << _popcnt64(w.bf12)
				<<" "<< _popcnt32(w.bf3) << endl;

		}
	}

	void DumpIndex() {
		cout << "index 27 to 81 "  << endl;
		for (int i = 0; i < 27; i++) {
			register uint32_t i81 = i_27_to_81[i];
			cout << i << " i81=" << i81 << " " << Char27out(g.pat2[i81]) << endl;
		}
		cout << "index 9 to 81 " << endl;
		for (int i = 0; i < 9; i++) {
			register uint32_t i81 = i_9_to_81[i];
			cout << i << " i81=" << i81 << " " << Char27out(g.pat3[i81]) << endl;
		}

	}
};

struct GEN_BANDES_12 {// encapsulating global data 
	STD_B3 bands3[512];
	int modeb12, go_back, diagmore, diagbug, ip20,
		it16, it16_2, imin16_1, imin16_2, imin16_3;
	int i1t16, i2t16, i3t16, maxnb3; // index 416 ordered in increasing size of valid clues 6
	char zsol[82], rband2[28];
	int grid0[81], tc[6], ntc;
	int gcheck[82], ib2check, ib3check;
	int skip, last;// restart point; last entry in the batch
	uint64_t   nb12;
	BANDMINLEX::PERM t_auto_b1[108], // maxi is 107excluding start position
		t_auto_b1b2[108], t_auto_b2b1[108],
		pband2, pband3, pcheck2, pcheck3;
	int n_auto_b1, n_auto_b1b2, n_auto_b2b1;
	int cold[9], coldf[9], rowd[6], boxd[6], rowdb3[3], boxdb3[3]; //free digits 
	//_________________ gangster 
	int gangcols[9];// 9 bits fields band3 for each colum (see valid band2)
	//int gang[9][3]; // gangcols expanded (buildgang ) 3 digits
	int gangb12[9]; // digit bf bands 12 per column
	//int   *gang27; // redefines gang[9][3] as 27 integer
	//int   gang_digits_cols[9][3];// active cols for a given digit
	//____________structs hosting the 81 GUA entries


	int GET_I81_G2(int digs, uint32_t pat) {
		register  uint32_t A = pat,
			B = (A | (A >> 9) | (A >> 18)) & 0777; // all columns
		uint32_t dstack;
		bitscanforward(dstack, B);
		dstack = (dstack / 3) * 27;// now stack0 0, 27 57 chunks of 27 i81
		for (uint32_t i = dstack; i < dstack + 27; i++)
			if (digs == tsgua2[i].digs) return i;
		return 0; // would be a bug
	}
	int GET_I81_G2_4(int digs, uint32_t pat) {
		register  uint32_t A = pat,
			B = (A | (A >> 9) | (A >> 18)) & 0777; // all columns
		// keep the stack with 2 columns
		for (int i = 0, mask = 7; i < 3; i++, mask <<= 3) {
			if (_popcnt32(B&mask) == 2) {
				B = mask;
				break;
			}
		}
		uint32_t dstack; // now same as pat 2 cells
		bitscanforward(dstack, B);
		dstack = (dstack / 3) * 27;// now stack0 0, 27 57 chunks of 27 i81
		for (uint32_t i = dstack; i < dstack + 27; i++)
			if (digs == tsgua2[i].digs) return i;
		return 0; // would be a bug
	}

	int GET_I81_G3(int digs, uint32_t pat) {
		register  uint32_t A = pat,
			B = (A | (A >> 9) | (A >> 18)) & 0777; // all columns
		uint32_t dstack;
		bitscanforward(dstack, B);
		dstack = (dstack / 3) * 27;// now stack0 0, 27 57 chunks of 27 i81
		for (uint32_t i = dstack; i < dstack + 27; i++)
			if (digs == tsgua3[i].digs) return i;
		return 0; // would be a bug
	}

	// __________________________  primary UAs tables and creation of such tables
	uint64_t  *ptua2;// pointer to current table cycle search 2/3
	uint32_t   nua2; // nua2 for cycle search 2/3  
	//================== bands 3 and gangster band 3 analysis
	int nband3;
	int   tcolok[2], ncolok;

	int ngua6_7, c1, c2, band, floors, digp, i81;
	uint64_t wua0, ua;// partial gua ua to check
	uint64_t tuacheck[100], tua_for_check[500];
	uint32_t uadigs_for_check[500], nua_for_check, nua_check;
	//================ A creating a catalogue for the 17 search 
	//sorted increasing number of valid bands 6 clues

	GEN_BANDES_12() {
		//gang27 = gang[0];
		//InitialSockets2Setup();
		//InitialSockets3Setup();
	}
	//void InitialSockets2Setup();// batch level
	//void InitialSockets3Setup();// batch level
	//================================= functions
	void GetStartB2(int i); // one of the 20 starts 
	void Start(int mode = 0);
	void NewBand1(int iw);
	int F17Novalid1_2();
	int Band2Check();
	int Band3Check();
	void Find_band2B();
	int ValidBand2();
	void ValidInitGang();
	void Find_band3B(int m10 = 1);
	void Find_band3B_pass1(int m10=1);
	//============= loops control for UAs 5;6;7 digits collection (collect more=
	int iband, ibox, iminirow, ibox2, iminirow2, pat1, pat2, ncells;
	int tcells[6], tcols[6];
	int bcols[2][9], mycols[9], myfloors;
	uint64_t mybf;
};

struct G17B {// hosting the search in 6 6 5 mode combining bands solutions
	G17B();// initial tasks all commands

	BF128 p17diag;// known 17 pattern for tests
	uint64_t pk54;
	int b3lim,	 aigstop, aigstopxy,
		npuz, a_17_found_here,nsearched ;
	int  debug17, debug17_check, diag, diagbug, debugb3,
		is_test_on,ng2,ng3;
	int grid0[81];

	//____gangsters, brute force,sockets setup
	//_________________ gangster 
	int gang[9][3]; // gangcols expanded (buildgang ) 3 digits
	int* gang27; // redefines gang[9][3] as 27 integer
	int   gang_digits_cols[9][3];// active cols for a given digit

	//______sockets common to  all bands 3  
	BF128 gsock2, gsock3;
	
	//============================ b12 no more uas to test
	uint64_t ua_ret7p, myb12, myb12f,
		myac, myacf,
		myb12add, myacadd, 
		anduab12, clean_valid_done;

	uint32_t tclues[40],tb3[512],*tcluesxy; 
	int nclues_step, nclues,ntb3;
	//BF128 bands_active_pairs, bands_active_triplets,
	//	valid_vect;
	//============  go band3
	int nclgo, nmiss;
	int  ncluesb3, mincluesb3;
	uint64_t tuaddb12[50];
	uint32_t tadd[50], ntadd, ntuaddb12;// add b1 b2
	uint32_t anduab3, stopexpandb3;// b3 expand


	STD_B3* myband3;
	uint32_t fstk, andmiss1, noutmiss1, wactive0;
	uint32_t free1, free2, free3;
	uint32_t tua128_b3[1000], ntua128_b3;
	BF128 wg46;
	//uint32_t cur_ib;
	uint32_t tcluesb12[20], ncluesb3x;
	uint32_t t3[2000], nt3,
		t3_2[1000], nt3_2,
		uasb3_1[2000], uasb3_2[2000], uas_in[2000],
		nuasb3_1, nuasb3_2, nuas_in, b3_andout;
	
	inline void AddB3_2_to_B3_1() {
		memcpy(&uasb3_1[nuasb3_1], uasb3_2, nuasb3_2 * 4);
		nuasb3_1 += nuasb3_2;
	}
	
	
	
	
	//=====================process for a new band 2 / set of bands 3
	void Start();// standard entry
	void StartPrint();// standard entry
	void StartKnown();// entry for a known 17
	void StartInit();// initial task gangster set up 
	void StartInitDebug();// initial task gangster set up 
	void UaCollector();
	inline void Adduab12(uint32_t digs, uint32_t nd);
	void FirstUasCollect();
	void SecondUasCollect();
	void UasCollect4box();
	void UasCollect6_7();

	void StartAfterUasHarvest();
	//inline int BuildGua(BF128& w);
	//inline void BuildGua(BF128& w, int cc);
	void Guas2Collect();
	void Guas2CollectG2();
	void Guas2CollectG3();
	void Guas2CollectG3_4d();
	void Guas2CollectG3_5d();
	void Expand_03();
	void Expand_46();
	void Expand_47();

	//int SetupExpand_7p();
	//void Init7p_guas();

	int IsValid7p(SPB03* sn);
	int IsValid_myb12();
	uint32_t IsValidB3(uint32_t bf);
	inline int GetNextCell(SPB03* s);
	inline void GetNextUa(SPB03* sn);
	inline void GetNextUaAdd(SPB03* sn);
	inline int GetLastAndUa(SPB03* sn, int diag = 0);
	int Expand_7_10();
	void GoExpand_7_10();
	void Go_9_10();
	void Go_8_10();

	int Expand_8_11_18();
	int Expand_8_11_17();
	void GoExpand_8_11();
	void Go_10_11();
	void Go_9_11();
	void Go_8_11();




	void Expand_7_11x();
	void Expand_7_12x();

	void GoBelow(SPB03* sn);
	void ExpandAddB1B2(SPB03* sn);
	void GoAfterExpand(SPB03* sn, uint32_t nadd=0);
	void GoB3CleanOne();
	void GoB3Miss0();
	void GoB3Miss1();
	void GoB3MissMore();

	void GoB3small(STD_B3* mb3);
	void GoB3Big(STD_B3* mb3);
	void ExpandB3Direct(int ntoass);
	void Out17(uint32_t bfb3);

	int nt4ok, okcheck;// for known
	// bands 1+2 valid epansion


};