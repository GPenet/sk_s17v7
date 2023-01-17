#define stack1_54 07007007007007007
struct OPCOMMAND {// decoding command line option for this rpocess
	// processing options 
	int opcode;
	int t18, p1, p2, p2b,//17 of 18 clues, pass or 2 (2a or 2b)
		p2c,//asking for list of attached ED grids (coded as t18 p2b)
		b2slice, // runing a slice of bands 2 in 18 mode bfx[0] & 8
		b3low, // running band 1 pass1 for slices in pass2 bfx[0] & 16
		out_one,// limit output to one per band 3 .bfx[2] & 1
		out_entry, //output of the entry file for test DLL .bfx[2] & 2
	    known; // 1 if known process 2 if known filter active .bfx[2] & 4
	// bfx[2] & 8 special use b2_is as limit b3
	int b1;//band 1 in process 
	int b2,b2_is ;//bands b2  forced
	char* b2start;
	int skip, last;
	int ton;//test on and test level
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
			b2slice = 1; 
		}
		if (p1 && (sgo.bfx[0] & 16))b3low = 1;

		b1 = sgo.vx[0];
		skip = sgo.vx[2];
		last = sgo.vx[3];
		b2_is = sgo.vx[4];
		b2 = sgo.vx[5];
		if (sgo.s_strings[0])	if(strlen(sgo.s_strings[0]))
			b2start = sgo.s_strings[0];
		ton= sgo.vx[1];

		f3 = sgo.vx[6];		f4 = sgo.vx[7];		f7 = sgo.vx[8];

		if (sgo.bfx[1] & 1)upto3 = 1;		if (sgo.bfx[1] & 2)upto4 = 1;
		if (sgo.bfx[1] & 4)dv12 = 1;
		if (sgo.bfx[1] & 8)dumpa = 1;

		if (sgo.bfx[2] & 1) out_one = 1;
		if (sgo.bfx[2] & 2) out_entry = 1;
		if (known)if (sgo.bfx[2] & 4) known = 2;

		// sgo.bfx[3] is for partial process 

		if (p) {
			cout << Char9out(sgo.bfx[0]) << " sgo.bfx[0 " << endl;
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
			if (b3low)
				cout << " pass1 with limit in band 3 index <= band1 index " << endl;
			if (out_one) cout << " max one out per band 3 sgo.bfx[2] & 1 " << endl;
			if (out_entry)  cout << " file1 contains attached solution grids" << endl;
			cout << "debugging commands___________________" << endl;
			if (known) {
				cout << "processing solution grids with known" << endl;
				if (known > 1)cout << "\tfilter on path active  sgo.bfx[2] & 2" << endl; 
			}
			cout << sgo.vx[5] << " b2 -v5- filter band 2 index" << endl;
			if (b2start)	cout << b2start << " filter band 2 start" << endl;

			if (ton)cout << ton << "  test on  -v1- verbose mode " << endl;
			if (f3)cout << f3 << "  f3  -v6- diag filter 3 clues [3]" << endl;
			if (f4)cout << f4 << "  f4  -v7- diag filter 6 clues [6]" << endl;
			if (f7)cout << f7 << "  f7  -v8- diag filter go full [7]" << endl;
			if (dv12)cout << "  -b1-..x  dump add in valid b12" << endl;
			if (dumpa)cout << "  -b1-...x  dump uas b12 at the start" << endl;
			if (upto3)cout << "upto debugging [3]  sgo.bfx[1] & 1 " << endl;
			if (upto4)cout << "upto debugging [4]  sgo.bfx[1] & 2 " << endl;
			if (dv12)cout << " print fresh adds sgo.bfx[1] & 4 " << endl;
			if (dumpa)cout << " print initial uas 12 sgo.bfx[1] & 8 " << endl;

		}
	}
}op;
struct CBS {// clues band stack
	uint16_t b[3];
	inline void Init(uint64_t bf54, uint16_t n) {
		register uint64_t U = bf54;
		b[0] =(uint16_t) _popcnt64(U & BIT_SET_27);
		b[1] = n - b[0];
		b[2]=0;
	}
	inline void Add(uint32_t cell) {b[cell / 27]++;	}
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

	inline uint64_t NextActive() {// called in expand 10_12
		if (b[0] > 6 || b[1] > 6) return 0;
		if(b[0]==6) return ~(uint64_t)BIT_SET_27;
		if(b[1]==6) return  BIT_SET_27;
		return ~0;
	}

	inline int IsFilt11_17(int p2b) {// 566 656 no stack>6
		if (b[0] > 6 || b[1] > 6)return 1;
		if(p2b && (b[1] > 5))return 1;
		return 0;
	}

	inline int IsFilt18p1() {
		if (b[0] > 8 || b[1] > 8)return 1;
		return 0;
	}
	inline int IsFilt12_18() {
		if (b[0] != 6)return 1;
		return 0;
	}
	void Status() {
		cout << " bx " << b[0] << b[1] << b[2];

	}
};
#define UACNBLOCS 15
//___ expand uas in bands 1+2
struct SPB03A {// spots 6 first clues
	BF128 v;	
	uint64_t  possible_cells, all_previous_cells, active_cells;
};
struct SPB03B {// spots 6 first clues
	BF128 v[UACNBLOCS];
	uint64_t  possible_cells, all_previous_cells, active_cells;
};
struct SPB03 {// spots to first 7 clues
	BF128 v[UACNBLOCS];// v2 v3 in expand 7_9 max 384 uas
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
		cout << Char64out(v[0].bf.u64[0]) << " 64 v" << endl;
	}
};
// expand uas in band 3
struct SP3 {
	uint32_t  possible_cells, all_previous, active,
		indtw3;
};

//___ permanent data for guas2 guas3
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

//_____ uas bands 1+2
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
		inline void DoAnd(BF128 v, uint64_t & wa) {
			int ir;
			while ((ir = v.getFirst128()) >= 0) {
				v.clearBit(ir);
				wa &= t[ir];
			}
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
	// working area for "build"
	uint64_t tw[UA12SIZE]; // to check redundancy
	BF128 vsize[25][UA12BLOCS];
	BF128 tvw[UA12BLOCS];

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
	int Build_tb128(SPB03A & s);

//   #define UACNBLOCS 15 see befor SPB 
	// status after 6 clues (C)
	TUVECT tc128[UACNBLOCS];// max start 10*128=1280 uas 
	uint64_t tandc;
	uint32_t nc128, ncblocs, ntc128[UACNBLOCS];
	void InitC() {
		memset(ntc128, 0, sizeof ntc128);
		nc128 = ncblocs = 0;
		for (int i = 0; i < UACNBLOCS; i++)tc128[i].Init();
		tandc = ~0;
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
		tandc &= R;
		register uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvc[cell].clearBit(ir);
		}
	}
	int Build_tc128(SPB03A& s3,SPB03A& s6);// after 6 clues
	inline int IsNotRedundant(uint64_t u) {
		register uint64_t nu = ~u;
		for (uint32_t i = 0; i < ntc128[0]; i++)
			if (!(tc128[0].t[i] & nu)) return 0;
		return 1;
	}

#define UADNBLOCS 10
	// status after 9 clues (D)
	TUVECT td128[UADNBLOCS];// max start 10*128=1280 uas 
	uint64_t tandd;
	uint32_t nd128, ndblocs, ntd128[UADNBLOCS];
	void InitD() {
		memset(ntd128, 0, sizeof ntd128);
		nd128 = ndblocs = 0;
		for (int i = 0; i < UADNBLOCS; i++)td128[i].Init();
		tandd = ~0;
	}
	inline void AddD(uint64_t u) {
		if (nd128 >= UADNBLOCS * 128) return;
		register uint32_t bloc = nd128 >> 7,
			ir = nd128 - ( bloc<<7);
		nd128++; ndblocs = bloc; ntd128[bloc]++;
		td128[bloc].v0.setBit(ir);
		BF128* myvd = td128[bloc].vc;
		register uint64_t R = u;
		td128[bloc].t[ir] = R;
		tandd &= R;
		register uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvd[cell].clearBit(ir);
		}
	}
	int Build_td128(SPB03B& s9);
	inline int IsNotRedundantD(uint64_t u) {
		register uint64_t nu = ~u;
		for (uint32_t i = 0; i < ntd128[0]; i++)
			if (!(td128[0].t[i] & nu)) return 0;
		return 1;
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
		for (uint32_t i = 0; i <= ncblocs; i++) {
			cout << "+ " << 128 * i << endl;
			tc128[i].Dump();
		}
	}
	void DebugD() {
		cout << " debugD nd128=" << nd128 << endl;
		td128[0].Dump();
	}

}t54b12;

//____ guas (uas with band3)
struct GUAH {// handler for initial collection of guas 2 3

	struct GUA {
		uint64_t tua[128];
		uint32_t nua, type, i81;
		inline void Add(uint64_t ua) {
			if (nua < 128)tua[nua++] = ua;
		}
		inline void AddIf(uint64_t ua) {
			for (uint32_t i = 0; i < nua; i++) {
				if (ua == tua[i]) return;
			}
			tua[nua++] = ua;
		}

		inline void Init(uint32_t n, uint32_t t, uint32_t i) {
			nua = n; type = t; i81 = i;
		}
		int Load(uint64_t* tu, uint64_t bf) {
			int n = 0;
			register uint64_t nbf = ~bf;
			for (uint32_t i = 0; i < nua; i++) {
				register uint64_t U = tua[i], cc = _popcnt64(U);
				if (n > 10) break;// 
				if (n > 5 && cc > 12) break;// 
				if (!(U & nbf))tu[n++] = U;
			}
			return n;
		}
		void SortClean() {
			if (nua < 2)return;
			GUA w = *this;
			BF128 vsize[30];
			memset(vsize, 0, sizeof vsize);
			uint64_t tcopy[128];
			memcpy(tcopy, w.tua, sizeof tcopy);
			for (uint32_t i = 0; i < nua; i++) {
				register uint64_t  cc = _popcnt64(tua[i] & BIT_SET_2X);
				vsize[cc].setBit(i);
			}
			nua = 0;
			for (int i = 0; i < 20; i++) if (vsize[i].isNotEmpty()) {
				int j;
				BF128 v = vsize[i];
				while ((j = v.getFirst128()) >= 0) {
					v.clearBit(j);
					register uint64_t U = tcopy[j];
					for (uint32_t k = 0; k < nua; k++) {
						register uint64_t U2 = tua[k];
						if (!(U & (~U2))) { U = 0; break; }
					}
					if (U) tua[nua++] = U;
					if (nua > 40) break;
				}
				if (nua > 40) break;
			}

		}
		void Debug(int nodet = 1) {
			cout << "gua type=" << type << " i81=" << i81
				<< "\tnua=" << nua << endl;
			if (nodet) return;
			for (uint32_t i = 0; i < nua; i++) {
				cout << Char2Xout(tua[i]) << " " << _popcnt64(tua[i] & BIT_SET_2X)
					<< " " << i << endl;
			}
		}
	}tg2[81], tg3[81], guaw;
	void Init() {
		for (int i = 0; i < 81; i++) {
			tg2[i].Init(0, 0, i); tg3[i].Init(0, 1, i);
		}
	}
	inline void Add2(uint64_t u, uint32_t i) { if (_popcnt64(u) < 18)		tg2[i].Add(u); }
	inline void Add3(uint64_t u, uint32_t i) { if (_popcnt64(u) < 18)		tg3[i].Add(u); }

	int IsUa4(int i81) {
		GUA& g = tg2[i81];
		if (g.nua == 1 && _popcnt64(g.tua[0]) == 2) return 1;
		return 0;
	}
	int IsUamin(int i81) {
		GUA& g = tg3[i81];
		if (g.nua == 1 && _popcnt64(g.tua[0]) < 5) return 1;
		return 0;
	}

	void SortClean3() {
		for (int i = 0; i < 81; i++)
			if (tg3[i].nua > 1) tg3[i].SortClean();
	}
	void SortClean() {
		for (int i = 0; i < 81; i++)
			if (tg2[i].nua > 1) tg2[i].SortClean();
	}

	int CutG2(int lim) {
		int n = 0;
		for (int i = 0; i < 81; i++) {
			if ((int)tg2[i].nua > lim)tg2[i].nua = lim;
			n += tg2[i].nua;
		}
		return n;
	}
	int CutG3(int lim) {
		int n = 0;
		for (int i = 0; i < 81; i++) {
			if ((int)tg3[i].nua > lim)tg3[i].nua = lim;
			n += tg3[i].nua;
		}
		return n;
	}
	void Dump2all2() {
		for (int i = 0; i < 81; i++) if (tg2[i].nua) {
			tg2[i].Debug(0);
		}
	}
	void Dump2all3() {
		for (int i = 0; i < 81; i++) if (tg3[i].nua) {
			tg3[i].Debug(0);
		}
	}
	void DumpOne2(int i) {
		cout << "debug2 i_1=" << i << endl;
		if (tg2[i].nua) 			tg2[i].Debug(0);
	}

	void DumpOne3(int i) {
		cout << "debug3 i_1=" << i << endl;
		if (tg3[i].nua) 			tg3[i].Debug(0);
	}
}guah;
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
#define G2ANBLOCS 100
#define G3ANBLOCS 100
	uint32_t nag2, nag3, nbag2, nbag3;
	uint64_t wbf;	uint32_t  wi, wc;
	struct G64 {
		uint64_t t12[64],v0,vc[54];
		uint32_t ti81[64],n;
		void Init() {
			n = 0;		v0 = 0;
			memset(vc, 255, sizeof vc);
		}
		int Add(uint64_t bf, uint32_t i81) {
			if (n >= 64) return 0;//safety
			t12[n] = bf;
			ti81[n] = i81;
			uint64_t bit= (uint64_t)1 << n++;
			v0 |= bit;
			register uint64_t B = bf, nbit = ~bit;
			int cell;
			while (bitscanforward64(cell, B)) {
				B^= (uint64_t)1 << cell;
				vc[cell] &= nbit;
			}
			return 1;
		}
		void Dump(int i0) {
			for (uint32_t i = 0; i < n; i++)
				cout << ti81[i] << "\t" << Char54out(t12[i])
				<< "  " << i + i0 << endl;
		}
	}ag2[G2ANBLOCS], ag3[G3ANBLOCS];
	// room for 6400 g2 or g3
	//___________________ start g54 mode 
	void InitA() {
		nag2 = nag3 = nbag2 = nbag3 = 0;
		for (int i = 0; i < G2ANBLOCS; i++) ag2[i].Init();
		for (int i = 0; i < G3ANBLOCS; i++) ag3[i].Init();
	}
	inline void Add2A() {
		if (nag2 >= 64 * G2ANBLOCS) return;
		if ((G2ANBLOCS - nbag2) < 3 && wc > 12)return;
		nbag2 = nag2++ >> 6; // new last bloc 
		ag2[nbag2].Add(wbf, wi);
	}
	inline void Add3A() {
		if (nag3 >= 64 * G3ANBLOCS) return;
		if ((G3ANBLOCS - nbag3) < 3 && wc > 12)return;
		nbag3 = nag3++ >> 6; // new last bloc 
		ag3[nbag3].Add(wbf, wi);
	}

	//_______________ after 6 clues

#define G2BNBLOCS 50
#define G3BNBLOCS 50
	uint32_t nbg2, nbg3, nbbg2, nbbg3;
	struct G64B {
		uint64_t v0, vc[54];
		uint32_t ti81[64], n;
		void Init() {
			n = 0;		v0 = 0;
			memset(vc, 255, sizeof vc);
		}
		int Add(uint64_t bf, uint32_t i81) {
			if (n >= 64) return 0;//safety
			ti81[n] = i81;
			uint64_t bit = (uint64_t)1 << n++;
			v0 |= bit;
			register uint64_t B = bf, nbit = ~bit;
			int cell;
			while (bitscanforward64(cell, B)) {
				B ^= (uint64_t)1 << cell;
				vc[cell] &= nbit;
			}
			return 1;
		}
		void Dump(int i0) {
			for (uint32_t i = 0; i < n; i++) {
				uint64_t bit = (uint64_t)1 << i;
				cout << ti81[i] << "\t";
				for (int j = 0; j < 54; j++)
					if (vc[j] & bit) cout << ".";
					else cout << "1";
				cout << "  " << i + i0 << endl;
			}
		}
	}bg2[G2BNBLOCS], bg3[G3BNBLOCS];
	void InitB() {
		nbg2 = nbg3 = nbbg2 = nbbg3 = 0;
		for (int i = 0; i < G2BNBLOCS; i++) bg2[i].Init();
		for (int i = 0; i < G3BNBLOCS; i++) bg3[i].Init();
	}
	inline void Add2B() {
		if (nbg2 >= 64 * G2BNBLOCS) return;
		if ((G2BNBLOCS - nbbg2) < 3 && wc > 12)return;
		nbbg2 = nbg2++ >> 6; // new last bloc 
		bg2[nbbg2].Add(wbf, wi);
	}
	inline void Add3B() {
		if (nbg3 >= 64 * G3BNBLOCS) return;
		if ((G3BNBLOCS - nbbg3) < 3 && wc > 12)return;
		nbbg3 = nbg3++ >> 6; // new last bloc 
		bg3[nbbg3].Add(wbf, wi);
	}

	//_______________ add fresh uas in check b3
	void AddA2(uint64_t bf, int i81, int cc);
	/*
	{
		wbf = bf; wi = i81, wc = cc;

		Add2A();
		if(cc)Add2B();
	}	*/

	void AddA3(uint64_t bf, int i81, int cc) {
		wbf = bf; wi = i81, wc = cc;
		Add3A();
		if(cc)Add3B();
	}

	void Build();
	void Build2(uint64_t filter, uint64_t active);
	void Build9(uint64_t filter, uint64_t active);
	void GetG2G3(BF128& g2, BF128& g3);

	void DumpA2() {
		cout << "dumpA2 ng2=" << nag2 << endl;
		for (uint32_t i = 0; i <= nbag2; i++)  
			ag2[i].Dump(i << 6);		 
	}
	void DumpA3() {
		cout << "dumpA3 ng3=" << nag3 << endl;
		for (uint32_t i = 0; i <= nbag3; i++)
			ag3[i].Dump(i << 6);
	}
	void DumpB2(int det=0) {
		cout << "dumpB2 ng2=" << nbg2 << endl;
		if(det)
		for (uint32_t i = 0; i <= nbbg2; i++)
			bg2[i].Dump(i << 6);
	}
	void DumpB3(int det=0) {
		cout << "dumpB3 ng3=" << nbg3 << endl;
		if(det)
		for (uint32_t i = 0; i <= nbbg3; i++)
			bg3[i].Dump(i << 6);
	}


}guah54;


//_____ processing band3
struct XQ {//to build the UAs b3 to expand
	uint32_t t1a, t1b; //27 bits field assigned 
	uint32_t  critbf,fa,fb;
	uint32_t t2a[12], t2b[30];// pairs bf2 other pairs and triplet
	uint32_t n2a, n2b, nb3,nmiss,nadded;
	uint32_t tin[400], tout[400];
	uint32_t nin, nout,nred;
	uint32_t iuas4;
	void Init(uint32_t cbf) { 
		t1a = t1b = n2a = n2b //= n4 = nm 
			=nin=nout=nadded= 0; 
		critbf = cbf;
	}
	inline void SetFilters () {
		if (!nmiss) {
			if ((fa = t1a)) {// initial safety
				critbf = 0;
				for (uint32_t i = 0; i < n2b; i++)
					critbf |= t2b[i];
			}
			fb = critbf;
		}
		else {	fa = 0; fb = BIT_SET_27;}
	}
	int Miss1ToMiss0();
	int MissxToMiss0(uint32_t ubf);
	int Miss0CheckTin();
	int NoRoomToAssign() {
		return(_popcnt32(t1a) >= nb3);
	}
	int NToAssign() {
		return( nb3- _popcnt32(t1a));
	}
	inline void SetFreshCrit() {
		fb = 0;
		for (uint32_t i = 0; i < n2b; i++) {
			fb |= t2b[i];
		}
		critbf = fb;
	}
	uint32_t  AddAssigned(uint32_t bf) {
		nadded++;
		t1a |= bf;
		fa = t1a;
		register uint32_t F = t1a, C = 0,n=n2b;
		n2b = 0;
		for (uint32_t i = 0; i < n; i++) {
			register uint32_t U = t2b[i];
			if (!(U & F)) {
				t2b[n2b++]=U;
				C |= U;
			}
		}
		fb=critbf = C;
		return C;
	}
	inline void Addin(uint32_t bf) { tin[nin++] = bf; }
	inline void Addout(uint32_t bf) { tout[nout++] = bf; }
	inline uint32_t GetAndout() {
		register uint32_t A = tout[0];
		for (uint32_t i = 1; i < nout; i++)
			A&=tout[i];
		return A;
	}
	void BuildCheckRedundant() {
		nred = n2b;
		if (nmiss) {
			nred += n2a;
			memcpy(&t2b[n2b], t2a, n2a * sizeof t2b[0]);
		}

	}
	int Isoutsize2();
	int Isoutsize3();
	int Isoutsize4();
	int Isoutsize5();
	int NoDisjoint() {
		register uint32_t r1 = tout[0] , r2 ;
		for (uint32_t i = 1; i < nout; i++) {
			register uint32_t u = tout[i];
			if (r2 |= r1 & u)return 1; 
				r1 |= u;
		}
		return 0;
	}
	int Min1_4Disjoint();
	void CleanIn();
	void CleanOut();
	void BuildAllOut();
	void BuildAllOutMiss0();
	void Dump1() {
		cout << "xq nmiss= " << nmiss<< " nb3="<<nb3 << endl;
		cout << Char27out(t1a) << "bf2 assigned" << endl;
		cout << Char27out(critbf) << " critbf" << endl;
		for(uint32_t i=0;i<n2a;i++)
			cout << Char27out(t2a[i])  << endl;
		cout << "t2b uas" << endl;
		for (uint32_t i = 0; i < n2b; i++)
			cout << Char27out(t2b[i]) << endl;

	}
	void Dump2() {
		cout << Char27out(fa) << " F assigned" << endl;
		cout << Char27out(fb) << " fb active" << endl;
		cout << Char27out(critbf) << " critbf" << endl;
	}
	void DumpOut() {
		cout << "out status nout=" <<nout<< endl;
		for (uint32_t j = 0; j < nout; j++)
			cout << Char27out(tout[j]) << endl;

	}

	void Status() {
		Dump1(); Dump2();
		cout << " nadded= " << nadded 
			<<" nout="<<nout << endl;
		for (uint32_t i = 0; i < nout; i++)
			cout << Char27out(tout[i]) << endl;
		cout << " nin=" << nin << endl;
		for (uint32_t i = 0; i < nin; i++)
			cout << Char27out(tin[i]) << endl;
		cout << "end xq status \n" << endl;
	}
}xq;
struct CALLBAND3 {
	BF128 g2t, g3t;
	uint64_t bf12;
	uint32_t ncl;
	CBS cbs;
	void Dump() {
		cout << Char54out(bf12) << " b12 to process ncl=" << ncl
			<< " " << cbs.b[0] << cbs.b[1] << endl;
		cout << Char64out(g2t.bf.u64[0]);
		cout << Char27out(g2t.bf.u32[2]) << " g2t" << endl;;
		cout << Char64out(g3t.bf.u64[0]);
		cout << Char27out(g3t.bf.u32[2]) << " g3t" << endl;;

	}
};
struct CRITB3 {
	uint32_t minix[4],// triplet bf1 bf2 bf3  
		critbf, pairs27, mincount,
		t1a,t2a[27],nt2a,
		assigned, active,
		ncl, nb3, nmiss;
	inline void Init(int ncb12, CBS& cbsx) {
		memset(this, 0, sizeof(*this));
		ncl = ncb12;
		if (op.t18) nb3 = 18 - ncl;
		else nb3 = 17 - ncl;
	}

	inline void Addoutbfone(uint32_t bf) { assigned |= bf; nmiss--; }
	inline int Addone(uint32_t i27) {// back 1 if not possible
		int bit27 = 1 << i27;
		int imini = i27 / 3, bitmini = 1 << imini;
		assigned |= bit27;
		if (!(bit27 & critbf)) {// clue added outfield
			if (!nmiss)return 1;// not possible
			nmiss--;
			if (!nmiss) active &= critbf;// no more outfield
			return 0;;
		}
		// now add in field within mincount
		if (minix[3] & bitmini) {// 2 clues expected
			critbf ^= bit27; minix[3] ^= bitmini;
			minix[1] ^= bitmini;//send the bit in one pair to assign
			if (!nmiss)active &= critbf; 
			return 0;
		}
		if (minix[2] & bitmini) {
			if (pairs27 & bit27) {// 2 clues not common clue one more clue
				if (!nmiss)return 1;// not possible
				critbf ^= bit27;	minix[2] ^= bitmini;
				nmiss--;
				if (!nmiss) active &= critbf;// no more outfield
				return 0;
			}
		}
		// now the last clue in the mini row
		register int mask = ~(7 << (3 * imini));// clear minirow
		critbf &= mask;
		if( !nmiss )active &= critbf; 
		if (minix[1] & bitmini) minix[1] ^= bitmini; // one clue expected 		
		else if (minix[2] & bitmini)minix[2] ^= bitmini;//common cell in 2 pairs 
		else if (minix[0] & bitmini)minix[0] ^= bitmini;// triplet
		return 0;
	}
	inline int AddAssign(uint32_t bfa, int debug = 0) {// back 1 if not possible
		if (assigned & bfa)return 1; //should never be
		if (debug)Status(" entry add assign ");
		active &= ~bfa; // minimum is to kill new assign
		register int i27, X = bfa;
		while (bitscanforward(i27, X)) {
			if (debug) {
				cout << "i27 =" << i27 << " ";
				Status(" assign i27 ");
			}
			X ^= 1 << i27;// clear bit
			if (Addone(i27))return 1;
			if (debug) cout << " back addone" << endl;
		}
		return 0;
	}

	void AssignBf2(int bf) {// all or part of the bf2
		minix[2] &= ~bf;// clean the bf2 
		int imini;
		while (bitscanforward(imini, bf)) {
			int  mask = 7 << (3 * imini), m27 = pairs27 & mask;
			bf ^= 1 << imini;// clean the bf
			critbf &= ~mask;// clean the mini row as crit field
			active &= ~mask;// no more clue in this minirow
			pairs27 ^= m27;// clean the 27 pairs to assign
			assigned |= mask ^ m27;// assign the third cell of the mini row
		}
	}
	inline void AssignCritical() {//always called with nmiss=0
		active = critbf;
		if (minix[2])  AssignBf2(minix[2]);
	}
	inline int GetToAss() { return (nb3 - _popcnt32(assigned)); }
	void Status(const char* lib) {
		cout << lib << "critical Status mincount =" << mincount << " nmiss=" << nmiss
			<< " nb3 = " << nb3 << endl;
		cout << Char27out(critbf) << " critical bf" << endl;
		cout << Char27out(pairs27) << " pairs 27" << endl;
		cout << Char27out(assigned) << " assigned" << endl;
		cout << Char27out(active) << " active" << endl;
		if (minix[1])cout << Char9out(minix[1]) << "     minis bf1" << endl;
		if (minix[2])cout << Char9out(minix[2]) << "     minis bf2" << endl;
		if (minix[3])cout << Char9out(minix[3]) << "     minis bf3" << endl;
		if (minix[0])cout << Char9out(minix[0]) << " mini triplets" << endl << endl;
	}

}scritb3;
struct CRITHANDLER {
	CRITB3 mycritb3;
	uint32_t* tuas, nuas;
	void Init(uint32_t* t) { tuas = t; nuas = 0; }
	inline void Add(uint32_t u) { tuas[nuas++] = u; }
	inline void AddIf(uint32_t u) {
		register uint32_t nu = ~u;
		for (uint32_t i = 0; i < nuas; i++)
			if (!(nu & tuas[i])) return; // subset or eqal
		tuas[nuas++] = u;
	}
	uint32_t* Lock() { return &tuas[nuas]; }
	inline uint32_t GetAnd() {
		if (!nuas) return 0;
		uint32_t wa = BIT_SET_27;
		for (uint32_t i = 0; i < nuas; i++)wa &= tuas[i];
		return wa;
	}
	void Dump() {
		cout << "uasb3 main table " << nuas << endl;
		for (uint32_t i = 0; i < nuas; i++) {
			if (tuas[i] & ~BIT_SET_27)
				cout << Char32out(tuas[i]) << "  t3infield with out bits" << endl;
			cout << Char27out(tuas[i]) << " " << i 
				<< " " << _popcnt32(tuas[i]) << endl;

		}
	}
};


// standard  bands  
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
	void MorphUas();
	void InitC10(int i);// known mode


	void PrintStatus();
};
struct STD_B3 :STD_B416 {// data specific to bands 3
	// permanent gangster information
	struct G {
		BF128 gsocket2, gsocket3;// active i81 mode 81 bits
		BF128 gsocket4, gsocket6;//g2 with 4 cells in band 3
		int pat2[81], pat3[81]; // storing ua bitfields
		int pat2_27[27]; // storing ua bitfields
		int ua2_imini[81], ua3_imini[81], ua2bit27[81];
	}g;
	int minirows_bf[9];
	int triplet_perms[9][2][3];
	uint32_t i_27_to_81[27], i_9_to_81[9]; //band i81 for guas guas3
	uint32_t i_81_to_27[81]; //band i81 for guas guas3
	uint32_t  poutdone;
	//_______________________
	void InitBand3(int i16, char* ze, BANDMINLEX::PERM& p);
	void Go(CALLBAND3& cb3);
	inline void Go2(int debug=0);
	int  Go3(CALLBAND3& cb3);
	//int  GoMiss1Out(int debug = 0);
	void GoAfter10(CALLBAND3& cb3);
	void GoAfter11(CALLBAND3& cb3);
	void GoAfter11Miss0(CALLBAND3& cb3);

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
	inline int GetI81_x(int bf) {
		for (uint32_t i = 0; i < 81; i++) {
			if (g.pat2[i] == bf) return i;
		}
		return -1;
	}


	void Pat_to_digitsx(int bf, int* tcells, int* tdigits, int& nt) {
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
	void DumpIndex() {
		cout << "index 27 to 81 " << endl;
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
	void Dumpg2() {
		cout << "g2 status " << endl;
		for (int i = 0; i < 81; i++) 
			if(g.gsocket2.On(i))	{
			cout << i << "\t"  << " " << Char27out(g.pat2[i]) << endl;
		}

	}

#define GM_NB4 50
#define GM_NBM 50

	struct GUM64 {
		uint64_t v0,  vc[55];//vc[54] is for 6 clues
		uint32_t tb3[64], n;
		void Init() {
			n = 0;		v0 = 0;
			memset(vc, 255, sizeof vc);
			vc[54] = 0;
		}
		int Add(uint64_t bf, uint32_t patb3) {
			if (n >= 64) return 0;//safety
			tb3[n] = patb3;
			uint64_t bit = (uint64_t)1 << n++;
			v0 |= bit; vc[54] |= bit;
			register uint64_t B = bf, nbit = ~bit;
			int cell;
			while (bitscanforward64(cell, B)) {
				B ^= (uint64_t)1 << cell;
				vc[cell] &= nbit;
			}
			return 1;
		}
		inline void Set6(uint32_t* t) {
			register uint64_t v = v0;
			for (int i = 0; i < 6; i++) v &= vc[t[i]];
			vc[54] = v;
		}
		inline uint64_t Getv(uint32_t* t, int nt) {
			register uint64_t v = vc[54];
			for (int i = 0; i < nt; i++) v &= vc[t[i]];
			return v;
		}

		void Dump(int i0) {
			for (uint32_t i = 0; i < n; i++) {
				uint64_t bit = (uint64_t)1 << i;
				cout <<Char27out(tb3[i]) << "\t";
				for (int j = 0; j < 54; j++)
					if (vc[j] & bit) cout << ".";
					else cout << "1";
				cout << "  " << i + i0 << endl;
			}
		}

		uint64_t bf12;// mode 54
		uint32_t bf3;
		inline int Count() {
			return ((int)_popcnt64(bf12) + _popcnt32(bf3));
		}
		void Print() {
			cout << Char54out(bf12) << "\t";
			cout << Char27out(bf3)
				<< " " << Count() << endl;
		}
	}tgm64[GM_NB4], tgm64m[GM_NBM];
	uint32_t nbgm,  nbbgm, nbgmm, nbbgmm;
	void InitTg() {
		nbgm = nbbgm =  0;
		for (int i = 0; i < GM_NB4; i++) tgm64[i].Init();
		for (int i = 0; i < GM_NBM; i++) tgm64m[i].Init();
	}
	inline void Addm4(BF128 &w) {// entry mode 3x
		if (nbgm >= 64 * GM_NB4) return;
		nbbgm = nbgm++ >> 6; // new last bloc 
		register uint64_t U = w.bf.u64[0];
		U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
		tgm64[nbbgm].Add(U, w.bf.u32[2]);
	}
	inline void Addmm(BF128& w) {// entry mode 3x
		if (nbgmm >= 64 * GM_NBM) return;
		nbbgmm = nbgmm++ >> 6; // new last bloc 
		register uint64_t U = w.bf.u64[0];
		U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
		tgm64m[nbbgmm].Add(U, w.bf.u32[2]);
	}

	inline void Set6clues(uint64_t known) {
		uint32_t tc[10], ntc = 0, cell;// build table of cells 
		register uint64_t F = known;
		while (bitscanforward64(cell, F)) {
			F ^= (uint64_t)1 << cell;
			tc[ntc++] = cell;
		}
		for (uint32_t i = 0; i <= nbbgm; i++) 
			tgm64[i].Set6(tc);
		for (uint32_t i = 0; i <= nbbgmm; i++) 
				tgm64m[i].Set6(tc);
	}
	void DumpTgm() {
		cout << "dumptgm ngm=" << nbgm << endl;
		for (uint32_t i = 0; i <= nbbgm; i++)
			tgm64[i].Dump(i << 6);
	}
	void DumpTgmm() {
		cout << "dumptgm ngmm=" << nbgmm << endl;
		for (uint32_t i = 0; i <= nbbgmm; i++)
			tgm64m[i].Dump(i << 6);
	}
	void Checkkown4();

};


//____ entry builder

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
	BANDMINLEX::PERM t_auto_b1b[108], // maxi is 107excluding start position
		t_auto_b1b2b[108], t_auto_b2b1b[108];
	int n_auto_b1b, n_auto_b1b2b, n_auto_b2b1b;



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
			if (_popcnt32(B & mask) == 2) {
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
	uint64_t* ptua2;// pointer to current table cycle search 2/3
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
	void Find_band3B_pass1(int m10 = 1);
	void Find_band3B_pass1B(int m10 = 1);
	//============= loops control for UAs 5;6;7 digits collection (collect more=
	int iband, ibox, iminirow, ibox2, iminirow2, pat1, pat2, ncells;
	int tcells[6], tcols[6];
	int bcols[2][9], mycols[9], myfloors;
	uint64_t mybf;
};


// main process 

struct G17B {// hosting the search in 6 6 5 mode combining bands solutions
	G17B();// initial tasks all commands

	BF128 p17diag;// known 17 pattern for tests

	inline int IspathB3(uint32_t bf) {
		if (!((~p17diag.bf.u32[2]) & bf)) return 1;
		return 0;
	}

	uint64_t pk54;
	int b3lim,	 aigstop, aigstopxy,knownt,
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
	BF128 dvect[2];
	uint64_t ua_ret7p, myb12, myandall,	myac,
		myb12_6,myb12_9,myac_9, myandall_9,
		anduab12, clean_valid_done,critical_done;

	uint32_t tclues6p[20]; 
	int nclues6p;
	//============  go band3
	int nclgo, nmiss;
	int  ncluesb3, mincluesb3;
	uint32_t anduab3;// b3 expand
	CALLBAND3 cb3;

	STD_B3* myband3;

	uint32_t t3o[1000], nt3o, t3ando, t3infield, t3outseen;
	uint32_t t3b[200], nt3b, t3c[100], nt3c;// after 3/6 clues
	uint32_t t3more[200], nt3more;
	uint32_t ntoassb3;
	void Dumpt3o() {
		cout << "t3o n=" << nt3o << endl;
		for (uint32_t i = 0; i < nt3o; i++) {
			if (t3o[i] & ~BIT_SET_27)
				cout << Char32out(t3o[i]) << "  t3outfield with out bits" << endl;
			cout << Char27out(t3o[i]) << " " << i
				<<" "<<_popcnt32(t3o[i]) << endl;

		}
	}
	void Dumpt3b() {
		cout << "t3b n=" << nt3b << endl;
		for (uint32_t i = 0; i < nt3b; i++)
			cout << Char27out(t3b[i]) << endl;
	}

	uint32_t t3[1000], nt3,		t3_2[1000], nt3_2,
		uasb3_1[1000], uasb3_2[1000], uas_in[1000],
		nuasb3_1, nuasb3_2, nuas_in, b3_andout;
	CRITHANDLER mycrh;
	inline void Init_t3o(CRITB3 & cr) {
		nt3o = 0; t3ando = ~0; 
		t3infield = cr.critbf;
	}
	inline void AddT3o(uint32_t u) {
		if (nt3o >= 1000)return;
		t3o[nt3o++] = u; 
		t3ando &= u;
	}
	int SortT3o(uint32_t active=0);
	void BuildSortT3b();
	inline uint32_t GetAndT3o() {
		uint32_t wa = BIT_SET_27;
		for (uint32_t i = 0; i < nt3o; i++)
			wa &= t3o[i];
		return wa;
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

	int IsValid7pbf(uint64_t bf);

	void Expand_03();
	void Expand_46(SPB03A& s3);
	void Expand_7_9(SPB03A& s6);
	void EndExpand_7_9();
	void GoCallB3(CALLBAND3& cb3);
	uint64_t GetLastD();


	void CheckCritical(CALLBAND3& cb3, CRITHANDLER & crh);
	uint32_t IsValidB3(uint32_t bf,int debug=0);

	// option no table of clues, no live count per band stack


	int IsValid_myb12();

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
	void Go_7_11_17();


	inline int GetNextCellD(SPB03* s);
	inline void GetNextUaD(SPB03* sn);
	inline void GetNextUaAddD(SPB03* sn);
	inline int GetLastAndUaD(SPB03* sn, int diag = 0);


	void Go_11_12();
	void Go_10_12();
	void Go_9_12();
	void Go_8_12();
	void Go_7_12();

	inline int VerifyValid() {
		if (clean_valid_done) return 0;
		if (clean_valid_done==2) return 1;
		clean_valid_done = 1;
		if (!IsValid_myb12()) return 0;
		clean_valid_done = 2; return 1;
	}
	int Valid3_1_3(uint32_t bf);
	int Valid3mm(uint32_t bf);
	void BridgeEndMiss0();
	void GoEndMiss0();
	void GoEndAll(uint32_t bf, uint32_t ac,int debug=0);
	void GoNotMiss0();
	void TryMiss1Subcritical();


	void GoB3Expand_1_3(uint32_t bf, uint32_t ac, int debug = 0);
	void GoB3Expand_4_x(SP3 spe, int debug = 0);

	void Out17(uint32_t bfb3);

	int nt4ok, okcheck;// for known
	// bands 1+2 valid epansion


};