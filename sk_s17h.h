struct MINCOUNT {
	uint32_t mini_bf1, mini_bf2, mini_bf3, mini_triplet,
		critbf, pairs27;
	uint32_t all_used_minis, mincount, minplus;
	inline void SetMincount() {// after direct setting minis
		all_used_minis = mini_bf1 | mini_triplet;
		mini_triplet &= ~mini_bf1;// count only triplets with no pair
		mincount = _popcnt32(all_used_minis) + _popcnt32(mini_bf3);
		mini_bf1 &= ~mini_bf2;// now pure one pair
		mini_bf2 &= ~mini_bf3;// now pure 2 pairs
		if (mini_triplet) {// must add triplet minirow
			for (int i = 0, bit = 1, field = 7; i < 9; i++, bit <<= 1, field <<= 3)
				if (mini_triplet&bit)
					critbf |= field;
		}
		minplus = mincount;
	}
	GINT64_t  Count_per_stack() {
		GINT64 cc; cc.u64 = 0;
		for (int i = 0, st = 0111; i < 3; i++, st <<= 1) {
			cc.u16[i] = _popcnt32(all_used_minis&st) +
				_popcnt32(mini_bf3&st);
		}
		return cc;
	}
};
struct BI2 {//Index 2 valid in band
	uint64_t bf, // 2 cells in bit fiekd
		active; // remaining possible cells
	uint32_t tval[2], // 2 cells in int mode
		istart, iend;// in the VALIB table
};
struct ZS128 {// main loop 128 bits 
	BF128 v;
	uint64_t bf, filler;
};

struct VALIDB {//  valid band 1+2
	uint64_t bf; // 2 cells in bit fiekd
	uint32_t tval[4]; // 3-4 cells in int mode
	uint64_t nval;// n cells over the 2
	inline void Enter(uint64_t ebf, uint32_t * tc2) {
		bf = ebf;
		nval = _popcnt64(bf) - 2;
		memcpy(tval, tc2, sizeof tval);
	}

};
struct VALIDB64 {// minimal valid band no index
	uint64_t bf; // 2 cells in bit fiekd
	uint32_t tval[6]; // 0-6 cells in int mode
	uint64_t nval;// n cells over the 2
};


// standard first band (or unique band)
struct STD_B416 {
	BI2 * my_bi2;
	VALIDB * my_validb;
	uint32_t nbi2, nvalidb; 
	char band[28];
	int band0[27], i416, gangster[9], map[27], dband;
	uint32_t tua[100], nua;//   maximum 81  
	uint32_t fd_sols[2][9];//start puzzle/ solution
	void Initstd();
	void GetBandTable(int i);
	void SetGangster();
	inline void GetUAs() {
		nua = t16_nua[i416];
		memcpy(tua, &t16_UAs[t16_indua[i416]], 4 * nua);
	}
	void MorphUas()	;
	void InitC10(int i);
	void InitG12(int i);
	void InitBand2_3(int i16, char * ze, BANDMINLEX::PERM & p
		, int iband = 1);
};
struct STD_B1_2 :STD_B416 {
	// row solution pattern in digit
	int mini_digs[9], mini_pairs[27],
		revised_g[9];// after false forced in 1/2 minirows
	int  tv_pairs[27], nvpairs; //9 or 27 bits 
	void FillMiniDigsMiniPairs(STD_B1_2 & bb);
	inline void InitRevisedg() {
		memcpy(revised_g, gangster, sizeof gangster);
	}
	int ReviseG_triplet(int imini, int ip, STD_B1_2 * bb);
	uint32_t GetMiniData(int index, uint32_t & bcells, STD_B1_2 *bb);
};
struct STD_B3 :STD_B416 {// data specific to bands 3
	BF128 tua128[1000], tua128_b1[1000], tua128_b2[1000];
	struct GUAs {
		BF128 isguasocket2, isguasocket3, isguasocket2_46;// active i81
		BF128 isguasocketc2, isguasocketc3, isguasocketc2_46;// active i81
		int triplet[9];//same gua3s
		int triplet_imini[81];
		int ua_pair[81], ua_triplet[81]; // storing ua bitfields
		int ua2_imini[81], ua3_imini[81],
			ua2pair27[81], ua2bit[81], ua3bit[81];
	}guas;
	BF128 isguasocketc246;//all linked to a socket 2
	uint32_t ntua128, ntua128_b1, ntua128_b2;
	int minirows_bf[9];
	int triplet_perms[9][2][3];
	uint32_t i_27_to_81[27], i_9_to_81[9]; //band i81 for guas guas3
	//_______________ handling mincount
	MINCOUNT smin, sminr;
	GINT64  stack_count;
	//_______________________
	void InitBand3(int i16, char * ze, BANDMINLEX::PERM & p);
	void InitStack(int i16, int  * z0, BANDMINLEX::PERM & p, int iband);
	int IsGua(int i81);
	int IsGua3(int i81);
	inline int GetI81_2(int bf) {
		for (uint32_t i = 0; i < 27; i++) {
			register uint32_t i81 = i_27_to_81[i];
			if (guas.ua_pair[i81] == bf) return i81;
		}
		return -1;
	}
	inline int GetI81_3(int bf) {
		for (uint32_t i = 0; i < 9; i++) {
			register uint32_t i81 = i_9_to_81[i];
			if (guas.ua_triplet[i81] == bf) return i81;
		}
		return -1;
	}
	int CleanG2(); //check guas2 
	int CleanG3(); //check guas3 
	inline void Insert2(uint32_t i81) {
		register uint32_t	bit = guas.ua2bit[i81];
		smin.mini_bf3 |= smin.mini_bf2&bit;
		smin.mini_bf2 |= smin.mini_bf1&bit;
		smin.mini_bf1 |= bit;
		smin.critbf |= guas.ua_pair[i81];
		smin.pairs27 |= guas.ua2pair27[i81];
	}
	inline void Insert3(uint32_t i81) {
		smin.mini_triplet |= guas.ua3bit[i81];
	}
};

struct BINDEXN {
	BI2 * t2;
	VALIDB * tvb;
	uint32_t nt2, ntvb;
	inline void Attach(BI2 * t2e, VALIDB * tvbe) { t2 = t2e; tvb = tvbe; }
	void Copy(STD_B1_2 & b);
	void Copy(BINDEXN & b);
}bin_b1, bin_b2, bin_b1yes, bin_b2yes,
bin2_b1, bin2_b2, bin2_b1yes, bin2_b2yes;

struct MORE64VECT {// FIFO table of more for bands 1+2
	BF128 vect,vc[54];
	uint64_t  t[128];
	uint32_t nt;
	inline void Init() {  
		nt = 0; 
		memset(&vect, 0, sizeof vect);
		memset(vc, 255, sizeof vc);

	}

	inline void Add(uint64_t v) {//add a new more if <128 
		if (nt < 128) {
			uint32_t cc64;// build cells vectors 
			vect.Set(nt);
			register uint64_t Rw = v;
			while (bitscanforward64(cc64, Rw)) {
				Rw ^= (uint64_t)1 << cc64;// clear bit
				if (cc64 > 26)cc64 -= 5;
				vc[cc64].clearBit(nt);
			}
			t[nt++] = v;
		}
	}
	inline void Add2(uint64_t v,uint64_t ac) {
		uint32_t cc64;// build cells vectors 
		vect.Set(nt);
		register uint64_t Rw = v&ac;
		while (bitscanforward64(cc64, Rw)) {
			Rw ^= (uint64_t)1 << cc64;// clear bit
			if (cc64 > 26)cc64 -= 5;
			vc[cc64].clearBit(nt);
		}
		t[nt++] = v;
	}
	inline int ApplyXY(uint32_t *tcells, uint32_t ntcells) {
		if(!nt) return 0;
		BF128 w = vect;
		for (uint32_t i = 0; i < ntcells; i++)
			w &= vc[tcells[i]];
		return (w.isNotEmpty()) ;
	}


};

struct MOREV2 {// 2 more64vect paired
	MORE64VECT mv1, mv2;
	uint32_t sw12;
	inline void Init() { mv1.Init(); mv2.Init(); sw12 = 0; }
	inline void Add(uint64_t v) {//add a new more if <128 
		if (sw12 &&(mv2.nt == 128)) {
			mv1.Init();
			sw12 = 0;// mv1 active
		}
		if ((!sw12) &&(mv1.nt == 128)) {
			mv2.Init();
			sw12 = 1;// mv2 active
		}
		if (sw12)mv2.Add(v);
		else mv1.Add(v);
	}
	inline int ApplyXY(uint32_t *tcells, uint32_t ntcells) {
		if (mv1.ApplyXY(tcells, ntcells)) return 1;
		return mv2.ApplyXY(tcells, ntcells);
	}

}morev2a,morev2b, morev2c;

struct MORE32 {// FIFO table of more for band b
	uint32_t  t[32];
	int nt, maxt, curt;
	inline void Init() { maxt = 32; nt = 0; }
	inline void Add(uint32_t v) {//add a new more in FIFO 
		if (nt < maxt) {// use next location
			curt = nt;
			t[nt++] = v;
		}
		else {// replace the oldest
			curt++;
			if (curt == maxt)curt = 0;
			t[curt] = v;
		}

	}

	inline int Check(uint32_t v) {// check oldest first
		if (!nt) return 0;
		register uint32_t V = v, *Rt = &t[curt], *Rtl;
	loop1:// form curt to 0
		if (!((*Rt) & V))return 1;
		if (--Rt >= t)goto loop1;
		if (nt < maxt) return 0;
		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt > Rtl)if (!((*Rt) & V))return 1;
		return 0;
	}
	inline uint32_t CheckNew(uint32_t bf, uint32_t & ua) {
		if (nt) {
			register uint32_t V = bf, *Rt = &t[nt];
			while (--Rt >= t) {
				if ((*Rt) & V) continue;
				ua &= *Rt;
			}
		}
		return ua;
	}
};

//================== UA collector 2 bands 
struct GENUAS_B12 {// for uas collection in bands 1+2 using brute force 
	int dig_cells[9][9], gangbf[9],// columns for band 3 in bit field
		revised_gangbf[9],// same revised UA2s UA3s ***
		mini_digs[9], mini_pairs[27], // UA2s UA3  ***
		nfloors, limstep, map[9], cptdebug, modemore;
	BF128 valid_sockets;
	//=============== uas collector 
	int limsize, floors;
	uint64_t  tuaold[TUA64_12SIZE],tua[TUA64_12SIZE], 
		 tuab1b2[200];// collecting bands uas in 2x mode
	uint32_t nuaold, nua, nuab1b2,		tuamore[500];
	//_______________uamore control
	STD_B1_2 *ba, *bb;
	uint32_t patb, ib, digp, colb, cola;
	uint64_t w0, ua, p12;
	//_____________________ functions collect UAs bands 1+2
	int Initgen();
	void BuildFloorsAndCollectOlds(int fl);
	inline void AddUA(uint64_t v) {	ua = v; AddUA64(tua, nua, ua);	}
	inline void AddUACheck(uint64_t v) {
		if (nua >= TUA64_12SIZE) nua = TUA64_12SIZE - 1;
		ua = v; AddUA64(tua, nua, ua);
	}
	int BuilOldUAs(uint32_t r0);
	int CheckOld();
	int CheckMain(uint64_t wua);
	void CollectMore2digits();
	void Collect2digits2_4_cols();
	void CollectMoreTry6_7();
	void EndCollectMoreStep();
	void CollectTriplets();
	void CollectMore2minirows();
	void ProcessSocket2(int i81);
};

struct GEN_BANDES_12 {// encapsulating global data 
	STD_B3 bands3[512];
	int modeb12, go_back,  ip20,
		it16, it16_2, imin16_1, imin16_2, imin16_3;
	int i1t16, i2t16, i3t16, maxnb3; // index 416 ordered in increasing size of valid clues 6
	char zsol[82], rband2[28];
	int grid0[82], tc[6], ntc;
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
	int gang[9][3]; // gangcols expanded (buildgang ) 3 digits
	int gangb12[9]; // digit bf bands 12 per column
	int   *gang27; // redefines gang[9][3] as 27 integer
	int   gang_digits_cols[9][3];// active cols for a given digit
	//____________structs hosting the 81 GUA entries
	struct SGUA2 {// 81 possible UA2 sockets
		// permanent data
		uint64_t * tua;
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
		uint32_t nua;// nua_start, nua_end;
	}tsgua2[81];

	//______ guas2 gua3s killed per cell assigned in band 1 2
	BF128 kguas2[54], kguas3[54];
	void SetUpkg();
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
	struct SGUA3 {// 81 possible UA3 sockets
		// permanent data
		uint64_t * tua;// , killer;
		int col1;// first columns 0-9 
		int i_81, imini;// , iguan; // index 0_80 for this 
		int id1, id2, id3; // index of digits in gang 27 
		// Current band1+2 data
		int  dig1, dig2, dig3,digs;// depending on gang27 status
		int valid, // valid if guas 
			validuas,// gua2s found
			used;// if needed in bands3
		uint32_t nua;// nua_start, nua_end, nua;
	}tsgua3[81];

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
		gang27 = gang[0];
		InitialSockets2Setup();
		InitialSockets3Setup();
	}
	void InitialSockets2Setup();// batch level
	void InitialSockets3Setup();// batch level
	void BuildGang9x3();
	void Build_CheckUAs_Subsets_Table();
	void Build_tuacheck(int fl);
	int Have_tuacheck_subset();
	void SecondSockets2Setup();// band1+2 level
	void SecondSockets2MoreUAs();// band1+2 level
	void GuaCollectMore();
	void SecondSockets3Setup();// band1+2 level
	void GuaCollect(int fl);
	//================================= functions
	void GetStartB2(int i); // one of the 20 starts 
	void Start(int mode = 0);
	void NewBand1(int iw);
	int Band2Check();
	int Band3Check();
	void Find_band2B();
	int ValidBand2();
	void ValidInitGang();
	void Find_band3B(int m10 = 1);

	//============= loops control for UAs 5;6;7 digits collection (collect more=
	int iband, ibox, iminirow, ibox2, iminirow2, pat1, pat2, ncells;
	int tcells[6], tcols[6];
	int bcols[2][9], mycols[9], myfloors;
	uint64_t mybf;
};


// guas table size 64 but SIZETGUA as initial limit
#define SIZETGUA 35

/*
a GUA can be 
type 
0 one of the 81 gangster pairs
1 one of the 81 gangster triplets

for a given band we can have in BF128 pattern
a  UA produced in band 3 process not gua2 gua3
a stack ua for the band 

in bit fields "active" 
all 81 types are encoded in 3x27 mode in a 128 bits field



At the start, a working GUA struct is open

for each  step 
all GUAs tables are first reduced to still active GUAs
then vectors are built to apply later xy clues
in a first pass, only the guas2 will be checked
this is usually enough to pass the limit count in band 3
*/
struct GUA {
	uint64_t tua[64];
	uint32_t nua,i162 ;
	inline void Add(uint64_t ua) {
		if (nua < 64)tua[nua++] = ua;
	}
	inline void Adduacheck(uint64_t ua) {
		if (nua > 63) nua = 63;
		AddUA64(tua, nua, ua);
	}
	inline void Init(uint32_t i) {
		nua = 0; i162 = i; ;
	}
};
struct TVG128 {// gua vector for 128 bits
	BF128 v, cells[54];
	uint32_t ti162[128];
	inline void SetVect54(uint64_t ua, uint32_t i128, uint32_t i162) {
		register uint64_t 	 W = ua;;
		if (!i128) {//new bloc to create
			v = maskLSB[1];
			memset(cells, 255, sizeof cells);
		}
		else v .Set(i128);
		ti162[i128] = i162;
		register uint32_t cc54;
		while (bitscanforward64(cc54, W)) {// look for  possible cells
			W ^= (uint64_t)1 << cc54;// clear bit
			cells[cc54].clearBit(i128);
		}
	}
	inline void ApplyXYcells(uint32_t *tc, uint32_t ntc) {
		for (uint32_t i = 0; i < ntc; i++)
			v &= cells[tc[i]];
	}

} tvg128g2[10], tvg128g3[5];// designed for 1280/640  uas
struct TGUAS {
	GUA tgua_start[162], tgua_b1[162];// max is 81+81
	uint32_t 	nvg2, nvg3, nguasb2,
		nb64_1, nb64_2, nb128_1, nb128_2;
	BF128 g2ok, g3ok, g2moreok;
	TGUAS() {
		for (uint32_t i = 0; i < 162; i++) {
			tgua_start[i].i162 = i;
			tgua_b1[i].i162 = i;
		}
	}
	inline void AddVG2_128(uint64_t ua, uint32_t i162) {
		if (nvg2 >= 1280) return;
		uint32_t ibloc = nvg2 >> 7, ir = nvg2 - 128 * ibloc;
		tvg128g2[ibloc].SetVect54(ua, ir, i162);
		nvg2++;
	}
	inline void AddVG3_128(uint64_t ua, uint32_t i162) {
		if (nvg3 >= 640) return;
		uint32_t ibloc = nvg3 >> 7, ir = nvg3 - 128 * ibloc;
		tvg128g3[ibloc].SetVect54(ua, ir, i162);
		nvg3++;
	}
	void ApplyLoopB1();	void ApplyLoopB2();
	int ApplyG2_128();	int ApplyG3_128();
}tguas;

struct G17B3HANDLER {
	MINCOUNT smin;
	int known_b3, rknown_b3, active_b3,	 ndead, 
		wactive0, nmiss, irloop, stack;
	uint32_t *uasb3if, nuasb3if, *uasb3of, nuasb3of, andoutf,
		active_sub, wua;
	GINT64 stack_count;
	int diagh;
	// ================== entry in the proces
	void GoMiss0();	
	void GoMiss1(uint32_t andout);	void Do_miss1();
	void GoMiss2Init();	void GoMiss2( uint32_t uamin);
	void AddCell_Miss2(uint32_t cell, uint32_t wand);
	uint32_t IsMultiple(int bf);
	//=============== process critical
	void CriticalAssignCell(int Ru);
	void Critical2pairs(int modesub=0);
	void CriticalLoop();
	//==================== process subcritical no cell added outside the GUAs field
	void SubMini(int M, int mask);
	void Go_Subcritical();
	void Go_SubcriticalMiniRow();

};

struct G17B {// hosting the search in 6 6 5 mode combining bands solutions
	G17B();// initial tasks all commands
	int b3lim,  aigstop, aigstopxy,
		iretb1,doloopnotok,		npuz, a_17_found_here;
	uint32_t	iband1,iband2, step1count;
	G17B3HANDLER hh0;
	//====== data for band expansion
	uint32_t nexp, bnua;
	uint64_t btua[300],start_active, b1cpt[8], b2cpt[8],b1cptdiag;
	BI2 * mybi2t,wbi_1,wbi_2;
	VALIDB * myokt, validb_known,validb_known_b1, validb_known_b2;
	uint32_t nmybi2t, nmyokt,nbi2_1,nbi2_2,
		nvb1,nvb1steps[10],
		nzs1_6,nzs2_6, nzs1_5, nzs2_5;
	//======= status after band 1 and step
	uint64_t tusb1[TUA64_12SIZE], tusb1_128[128];
	uint32_t ntusb1, ntusb1_128,ntusb2;
	uint64_t fb1,acb1, fb2, fb12, acb2, acb12;

	//_____ data and functions of the external loop

	int loopb1;
	struct EXTL {
		uint64_t noxyes;
		uint32_t bfx, tbfy[20], ntbfy, mode, ratio;
		void Init(uint32_t bf, int mod) { bfx = bf; ntbfy = 0; mode = mod; }
	}extl1[10], extl2[10], extlw, extlr;
	uint32_t nextl1, nextl2;
	uint64_t n_yesb1, n_yesb2, n_nob1, n_nob2,minratio;
	uint64_t FindSockets(uint64_t active, uint64_t lim);
	void ExtractMin(uint64_t active, BINDEXN & bin1, BINDEXN & bin2);
	void ExtSplitY(BINDEXN & binw, uint32_t *tbf, uint32_t ntbf,
		uint32_t & activer,int bande=1);
	void ExtSplitX(BINDEXN & bin1no, BINDEXN & bin1yes,
		uint32_t bf, uint32_t & activer,int bande=1);
	void Go2_Ext_Loop();
	void Go2b_Ext_Loop(uint64_t activeloop, uint32_t mode2);

	//___ process sub lots 
	uint32_t nbi5_1, nbi5_2, nbi6_1, nbi6_2;
	void Go3_Build_Band1(uint32_t ib1, BINDEXN & binw);
	void Go3_Build_Band2(uint32_t ib2, BINDEXN & binw);
	void Go3_Apply_B1_V();
	int Go3_Apply_B2_V();
	void Go3(BINDEXN & bin1, BINDEXN & bin2);

	//============ vectors 64 128 bits 
	struct UB2 {// for 2560 uas over 128 in step b1b2
		BF128 vx[20], vcx[20][54];
		void Init() {
			memset(vx, 0, sizeof vx);
			memset( vcx, 255, sizeof vcx);
		}
		inline void Add(uint64_t & ua, uint32_t i) {
			uint32_t ibloc = i >> 7, ir = i & 127;
			vx[ibloc].Set(ir);
			uint32_t cc64;// build cells vectors 
			register uint64_t Rw = ua;
			BF128 * vc = vcx[ibloc];
			while (bitscanforward64(cc64, Rw)) {
				Rw ^= (uint64_t)1 << cc64;// clear bit
				if (cc64 > 26)cc64 -= 5;
				vc[cc64].clearBit(ir);
			}
		}
		void ApplyStep(uint64_t & bf, uint32_t nuas) {
			// put bf in table
			uint32_t tcells[54], ntcells = 0,cc64;
			register uint64_t Rw = bf;
			while (bitscanforward64(cc64, Rw)) {
				Rw ^= (uint64_t)1 << cc64;// clear bit
				if (cc64 > 26)cc64 -= 5;
				tcells[ntcells++]=cc64;
			}
			uint32_t ibloc = 0;
			while (1) {
				BF128 w = vx[ibloc], *vc = vcx[ibloc];
				for (uint32_t i = 0; i < ntcells; i++)
					w &= vc[tcells[i]];
				vx[ibloc] = w;
				if (nuas > 128) { nuas -= 128; ibloc++; }
				else break;
			}
		}
		inline int ApplyXY(uint32_t *tcells, uint32_t ntcells, uint32_t nuas) {
			uint32_t ibloc = 0;
			while (1) {
				BF128 w = vx[ibloc], *vc = vcx[ibloc];
				for (uint32_t i = 0; i < ntcells; i++)
					w &= vc[tcells[i]];
				if (w.isNotEmpty()) return 1;
				if (nuas > 128) { nuas -= 128; ibloc++; }
				else break;
			}
			return 0;
		}
	}ub2;

	//============ main loop process
	BF128 v128uas, vc128[54];// 0 to 128
	// parameter for chunk 128
	ZS128 *za, *zb;
	uint32_t nza, nzb;
	uint64_t  n_to_clean, n_to_clean2,nwc;

	//___________ extract potential valid bands 1+2 (no more uas)

	void Do128uas();
	void Do128Chunk();
	void Do128Go(ZS128 * a, ZS128 * b, uint32_t na, uint32_t nb);

	//============================ b12 no more uas to test
	uint64_t wb12bf, wb12active,myua;
	GINT64  stack_count_step, stack_count, stack_countf;
	uint32_t tclues[40],*tcluesxy; 
	int nclues_step, nclues,  clean_valid_done;
	uint32_t	ua_of_seen;
	//============  go band3
	STD_B3* myband3;
	uint32_t fstk, andmiss1, noutmiss1, wactive0;
	BF128 wg46;
	int ib3_current;
	uint32_t   nmiss, ua_out_seen;
	uint32_t uasb3_1[2000], uasb3_2[2000],
		nuasb3_1, nuasb3_2, b3_andout;
	MINCOUNT smin;
	MORE32 moreuas_b3, moreuas_b3_small;

		
	//=====================process for a new band 2 / set of bands 3
	void GoM10();// standard entry
	void GoM10Uas();// collect uas guas
	void StackUas();
	void ExpandB1();
	void ExpandB2();
	void ExpandB3();
	void ExpandOneBand(int ib);// find 2-6 valid bands no redundant clue
	//_______ processing potential valid bands 1+2
	int Clean_Valid();// test first valid

	int Clean_2a();// test first valid
	int Clean_2b();// test first guas2
	int Clean_2c();// test only guas2
	void CleanAll();

	void NewUaB12();
	void DebugAdd12();
	void GoB3(STD_B3 & b);

	void FinalCheckB3(uint32_t bfb3);
	void Out17(uint32_t bfb3);
	void NewUaB3();
};