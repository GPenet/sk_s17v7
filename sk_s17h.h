#define stack1_54 07007007007007007


struct OPCOMMAND {// decoding command line option for this rpocess
	// processing options 
	int opcode;
	int t18,p1,p2,p2b,//17 of 18 clues, pass or 2 (2a or 2b)
		p2c,//asking for list of attached ED grids (coded as t18 p2b)
		b2slice; // runing a slice of bands 2 in 18 mode 
	int b1;//band 1 in process 
	// debugging options
	int b2,b2_is ;//bands usually b2 not given
	char* b2start;
	int skip, last;
	int ton;//test on and test level
	int known; // 1 if known process 2 if known filter active 
	uint64_t f3, f4, f7; // filters p_cpt2g [3] [4) [7]
	int upto3, upto4; // active below f3 below f4
	int dv12, dv3;// print fresh uas bands 1 2 band 3
	int dumpa;
	void SetUp(int opcod,int k = 0,int p=1) {// init known or not
		memset(this, 0, sizeof * this);
		opcode = opcod;
		known = k;
		if (sgo.bfx[0] & 1)t18 = 1;
		if (sgo.bfx[0] & 6) {// pass 1 2a 2b
			p2 = 1;
			if (sgo.bfx[0] & 4) p2b = 1;
			if (p2b && t18) { p2b = 0; p2c = 1; }
		}
		else p1 = 1;
		if(t18 && (sgo.bfx[0] & 8)) {// slice of bands 
			b2slice = 1; b2_is= sgo.vx[4];
		}

		if (known)if (sgo.bfx[0] & 8) known = 2;
		b1 = sgo.vx[0];
		skip = sgo.vx[2];
		last = sgo.vx[3];
		b2 = sgo.vx[5];
		if (sgo.s_strings[0])	if(strlen(sgo.s_strings[0]))
			b2start = sgo.s_strings[0];
		ton= sgo.vx[1];
		f3 = sgo.vx[6];
		if (sgo.bfx[1] & 1)upto3 = 1;
		f4 = sgo.vx[7];
		if (sgo.bfx[1] & 2)upto4 = 1;
		f7 = sgo.vx[8];
		if (sgo.bfx[1] & 4)dv12 = 1;
		if (sgo.bfx[1] & 8)dumpa = 1;

/*
-b1- 1 mode up to 3	2 mode upto 4		4 cout adds 12 		8 debug A

	int it16_start = sgo.vx[0];
	g17b.debug17_check =   g17b.aigstop=0;
	g17b.diag = g17b.debug17 = sgo.vx[1];
	genb12.skip = sgo.vx[2];
	genb12.last = sgo.vx[3];
*/


		if (p) {
			cout << "standard processing commands_______________" << endl;
			if(t18) cout <<"\t\tsearch 18 clues via -b0-x."<<endl;
			else cout << "\t\tsearch 17 clues via -b0-x." << endl;
			if(p1)cout << "\t\tpass1 via -b0-.x." << endl;
			if (p2)cout << "\t\tpass2 via -b0-.x." << endl;
			if (p2b)cout << "\t\tpass2b via -b0-..x." << endl;
			if (p2c) cout << " file1 contains attached solution grids" << endl;
			cout << sgo.vx[0] << " b1  -v0- band 0_415" << endl;
			cout << sgo.vx[2] << " skip  -v2- skip first nnn restart after batch failure" << endl;
			cout << sgo.vx[3] << " last  -v3- last entry number for this batch must be > vx[2]" << endl;
			if (b2slice) {
				cout << "running a slice of bands 2 index from="
					<< b2_is << " to=" << b2 << endl;
			}
			
			cout << "debugging commands___________________" << endl;
			if (known) {
				cout << "processing solution grids with known" << endl;
				if (known > 1)cout << "\tfilter on path active" << endl; 
			}
			cout << sgo.vx[5] << " b2 -v5- filter band 2 index" << endl;
			if (b2start)	cout << b2start << " filter band 2 start" << endl;

			if (ton)cout << ton << "  test on  -v1- verbose mode " << endl;
			if (f3)cout << f3 << "  f3  -v6- diag filter 3 clues [3]" << endl;
			if (f4)cout << f4 << "  f4  -v7- diag filter 6 clues [6]" << endl;
			if (f7)cout << f7 << "  f7  -v8- diag filter go full [7]" << endl;
			if (dv12)cout << "  -b1-..x  dump add in valid b12" << endl;
			if (dumpa)cout << "  -b1-...x  dump uas b12 at the start" << endl;
			if (sgo.bfx[1])cout << "upto debugging  through -b1-"
				<< Char9out(sgo.bfx[1]) << endl;
			if (sgo.bfx[2])cout << "partial process through -b2-"
				<< Char9out(sgo.bfx[2]) << endl;
			if (sgo.bfx[3])cout << "add b3 diag through -b3-"
				<< Char9out(sgo.bfx[3]) << endl;

		}
	}
}op;
struct CBS {// clues band stack
	uint16_t b[3];
	uint16_t s[3];
	inline void Init(uint64_t bf54, uint16_t n) {
		register uint64_t U = bf54;
		b[0] =(uint16_t) _popcnt64(U & BIT_SET_27);
		b[1] = n - b[0];
		b[2]=0;
		register uint64_t S = stack1_54;
		s[0]=(uint16_t) _popcnt64(U & S);
		S <<= 3;
		s[1] =(uint16_t) _popcnt64(U & S);
		s[2] = n - s[0] - s[1];
	}
	inline void Add(uint32_t cell) {
		b[cell / 27]++;
		s[C_stack[cell]]++;
	}
	inline int IsFilt10_17() {// 737 377 is possible
		if (b[0] > 7 || b[1] > 7)return 1;
		return 0;
	}
	inline uint64_t LimBand6(uint64_t bf) {
		if (b[0] >6 || b[1] > 6)return 0;
		if (b[0] == 6) return bf & (~BIT_SET_27);
		if (b[1] == 6) return bf & BIT_SET_27;
		return bf;
	}
	inline int StackMore6() {
		if (s[0] > 6 || s[1] > 6 || s[2] > 6) return 1;
		return 0;
	}
	inline uint64_t NextActive() {// called in expand 10_12
		if (b[0] > 6 || b[1] > 6) return 0;
		if(b[0]==6) return ~(uint64_t)BIT_SET_27;
		if(b[1]==6) return  BIT_SET_27;
		return ~0;
	}
	inline uint64_t NextActiveStack() {// called in expand 10_12
		if (s[0] > 6 || s[1] > 6 || s[2] > 6) return 0;
		if (s[0] == 6) return ~stack1_54;
		if (s[1] == 6) return ~(stack1_54 <<3);
		if (s[2] == 6) return ~(stack1_54 << 6);
		return ~0;
	}
	inline int IsFilt11_17(int p2b) {// 566 656 no stack>6
		if (b[0] > 6 || b[1] > 6)return 1;
		if (s[0] > 6 || s[1] > 6 || s[2] > 6)return 1;
		if(p2b && (b[1] > 5))return 1;
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
	void Status() {
		cout << " bx " << b[0] << b[1] << b[2]
			<< " sx " << s[0] << s[1] << s[2];

	}
};
struct SPB03 {// spots to first 7 clues
	BF128 v;
	uint64_t  possible_cells, all_previous_cells, active_cells;
	CBS cbs;
	uint32_t ncl;
	void Init9(BF128 w, CBS & c){
		all_previous_cells= w.bf.u64[0];
		active_cells = w.bf.u64[1];
		ncl = 9;
		cbs = c;
	}
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
	uint64_t bf12;
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
	inline int IsNotRedundant(uint64_t u) {
		register uint64_t nu = ~u;
		for (uint32_t i = 0; i < ntc128[0]; i++)
			if (!(tc128[0].t[i] & nu)) return 0;
		return 1;
	}

#define UADNBLOCS 10
	// status after 9 clues (D)
	TUVECT td128[UADNBLOCS];// max start 10*128=1280 uas 
	uint32_t nd128, ndblocs, ntd128[UADNBLOCS];
	void InitD() {
		memset(ntd128, 0, sizeof ntd128);
		nd128 = ndblocs = 0;
		for (int i = 0; i < UADNBLOCS; i++)td128[i].Init();
	}
	inline void AddD(uint64_t u) {
		if (nd128 >= UADNBLOCS * 128) return;
		register uint32_t bloc = nd128 >> 7,
			ir = nd128 - ( bloc<<7);
		nd128++; ndblocs = bloc; ntd128[bloc]++;
		td128[bloc].v0.setBit(ir);
		BF128* myvc = td128[bloc].vc;
		register uint64_t R = u;
		td128[bloc].t[ir] = R;
		register uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvc[cell].clearBit(ir);
		}
	}
	int Build_td128(); 
	inline int IsNotRedundantD(uint64_t u) {
		register uint64_t nu = ~u;
		for (uint32_t i = 0; i < ntd128[0]; i++)
			if (!(td128[0].t[i] & nu)) return 0;
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
	void DebugD() {
		cout << " debugD nd128=" << nd128 << endl;
		td128[0].Dump();
	}

}t54b12;
/*
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
*/

struct GUA54 {
	uint64_t* tua, killer;
	uint32_t nua, nuamax, type, i81;
	inline void Init(uint64_t* p, uint32_t t, uint32_t i) {
		tua = p; type = t; i81 = i; killer = ~0;
		nua = 0;
		nuamax = 10;
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

	inline int Check(uint64_t u) {// no redundancy
		register uint64_t nU = ~u;
		for (uint32_t j = 0; j < nua; j++)
			if (!(tua[j] & nU)) return 1;
		return 0;
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
	void Build9(uint64_t filter, uint64_t active);
	BF128 GetG2(uint64_t bf);
	BF128 GetG3(uint64_t bf);
	int  Check2(uint64_t bf, int i81)	{
		if (tg2[i81].Check(bf)) return 1;
		return 0;
	}
	int  Check3(uint64_t bf, int i81) {
		if (tg3[i81].Check(bf)) return 1;
		return 0;
	}


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
	void DumpOne3(int i81) {
		tg3[i81].Debug(0);
	}


	void Dumpall3() {
		for (int i = 0; i < 81; i++) {
			tg3[i].Debug(0);
		}
	}
}guah54, guah54_2, guah54_9;

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
		inline int Count() {
			return ((int)_popcnt64(bf12) + _popcnt32(bf3));
		}
		void Print(){
			cout <<Char54out(bf12)<<"\t";
			cout << Char27out(bf3) 
				<< " " << Count() << endl;
		}
	}tguam[384],tguam2[300], tguam9[300];
	uint32_t ntguam,ntguam2, ntguam9, guam2done, guam9done;
	//_______________________

	void InitBand3(int i16, char * ze, BANDMINLEX::PERM & p);
	void Go(CALLBAND3& cb3);
	inline void BuildGuam2(uint64_t known) {
		register uint64_t F = known, n = 0;
		for (uint32_t i = 0; i < ntguam; i++)
			if (!(F & tguam[i].bf12))
				if(n<300)tguam2[n++] = tguam[i];
		ntguam2 = (int)n;
		guam2done = 1;
		guam9done = 0;
	}
	inline void BuildGuam9(uint64_t known) {
		register uint64_t F = known, n = 0;
		for (uint32_t i = 0; i < ntguam2; i++)
			if (!(F & tguam[i].bf12))
				if (n < 300)tguam9[n++] = tguam2[i];
		ntguam9 = (int)n;
		guam9done = 1;
	}

	void Pack() {// keep only the best  
		GUAM  tt[384];
		BF128 vsize[14][3];
		uint32_t ntt = 0;
		memset(vsize, 0, sizeof vsize);
		for (uint32_t i = 0; i < ntguam; i++) {
			int cc = tguam[i].Count(),bloc=ntt>>7,ir=ntt-(bloc<<7);
			if (cc >= 20 || cc<6) continue;
			tt[ntt++] = tguam[i];
			vsize[cc-6][bloc].setBit(ir);
		}
		ntguam=0;
		for (int  i = 0; i <14; i++) {
			BF128 * w  = vsize[i];
			for (int j = 0; j < 3; j++) {
				BF128 V = w[j];
				while (1) {
					register int ir = V.getFirst128();
					if (ir >= 0) {
						V.clearBit(ir);
						tguam[ntguam++] = tt[ir+128*j];
					}
					else break;
				}
			}
			if (ntguam >= 256) break;
		}
		if (ntguam > 256) ntguam = 256;

		//cout << "check pack guam" << endl;
		//for (uint32_t i = 0; i < 10; i++)
			//tguam[i].Print();
		//for (uint32_t i =250; i < 256; i++)
			//tguam[i].Print();

	}
	void Addguam(BF128 w,int go2=0) {// entry mode 3x32
		if (ntguam >= 384) Pack();
		tguam[ntguam].bf3 = w.bf.u32[2];
		register uint64_t U = w.bf.u64[0];
		U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
		tguam[ntguam++].bf12 = U;
		if (go2&& ntguam2<300) {
			tguam2[ntguam2].bf3 = w.bf.u32[2];
			tguam2[ntguam2++].bf12 = U;

		}
	}
	int Check(BF128 w) {
		register uint64_t U = w.bf.u64[0];
		U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
		register uint32_t u3 = w.bf.u32[2];
		for (uint32_t i = 0; i < ntguam; i++) {
			if (U == tguam[i].bf12 && u3 == tguam[i].bf3) return 1;
		}
		return 0;
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
	//int skip, last;// restart point; last entry in the batch
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
		npuz, a_17_found_here ;
	int ng2,ng3;
	int grid0[81];

	//____gangsters, brute force,sockets setup
	//_________________ gangster 
	int gang[9][3]; // gangcols expanded (buildgang ) 3 digits
	int* gang27; // redefines gang[9][3] as 27 integer
	int   gang_digits_cols[9][3];// active cols for a given digit

	//______sockets common to  all bands 3  
	BF128 gsock2, gsock3;
	
	//============================ b12 no more uas to test
	uint64_t ua_ret7p, myb12, myandall,	myac,
		myb12_9,myac_9, myandall_9,
		anduab12, clean_valid_done;

	uint32_t tclues[20]; 
	int ncluesx;
	//============  go band3
	int nclgo, nmiss;
	int  ncluesb3, mincluesb3;
	uint32_t anduab3, stopexpandb3;// b3 expand


	STD_B3* myband3;
	uint32_t t3[1000], nt3,		t3_2[1000], nt3_2,
		uasb3_1[1000], uasb3_2[1000], uas_in[1000],
		nuasb3_1, nuasb3_2, nuas_in, b3_andout;
	inline void AddT3If(uint32_t u) {
		register uint32_t nu = ~u;
		for (uint32_t i = 0; i < nt3; i++)
			if (!(nu & t3[i])) return; // subset or eqal
		t3[nt3++] = u;
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

	uint32_t IsValidB3(uint32_t bf);
	inline int GetNextCell(SPB03* s);
	inline void GetNextUa(SPB03* sn);
	inline void GetNextUaAdd(SPB03* sn);
	inline int GetLastAndUa(SPB03* sn, int diag = 0);

	// option no table of clues, no live count per band stack

	int IsValid7pbf(SPB03* sn);

	int IsValid_myb12();

	int Expand_7_10();
	void GoExpand_7_10();
	void Go_9_10();
	void Go_8_10();
	void Go_7_10();

	void Go_10_11_18();
	void Go_9_11_18();
	void Go_8_11_18();
	void Go_7_11_18();
	void Go_10_11_17();
	void Go_9_11_17();
	void Go_8_11_17();
	int Expand_7_11();
	void GoExpand_7_11();


	inline int GetNextCellD(SPB03* s);
	inline void GetNextUaD(SPB03* sn);
	inline void GetNextUaAddD(SPB03* sn);
	inline int GetLastAndUaD(SPB03* sn, int diag = 0);


	void Go_11_12();
	void Go_10_12();
	void Go_9_12();
	void Go_8_12();


	void Expand_7_9();// 18 clues pass2
	void Expand_10_12();// 18 clues pass2
	void GoExpand_7_12();
	void GoExpand_10_12(BF128 ww);


	void GoB3CleanOne();
	void GoB3Miss0();
	void GoB3Miss1();
	void GoB3MissMore();

	void GoB3End(int ntoass);
	void GoB3Expand (int ntoass);
	void Out17(uint32_t bfb3);

	int nt4ok, okcheck;// for known
	// bands 1+2 valid epansion


};