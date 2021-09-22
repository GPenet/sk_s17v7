
//#define DODEBUG
#ifdef DODEBUG

//___ start process expand bands collect uas guas ...
const char * diagband = "274965318316847925598231467"; // la bonne
const char * diagband3 = "635714892782593641941628573";
const char * diagpuz = ".2.45......7....36..................3.6..7......2..4........8.......3..194....5..";
#endif
// initial taks all commands
G17B::G17B() {
	bin_b1.Attach(bi2_b1w, vab1w);
	bin_b2.Attach(bi2_b2w, vab2w);
	bin_b1yes.Attach(bi2_b1yes, vab1yes);
	bin_b2yes.Attach(bi2_b2yes, vab2yes);
	bin2_b1.Attach(bi2_b1w2, vab1w2);
	bin2_b2.Attach(bi2_b2w2, vab2w2);
	bin2_b1yes.Attach(bi2_b1yes2, vab1yes2);
	bin2_b2yes.Attach(bi2_b2yes2, vab2yes2);
	debug17 = debug17_check=aigstop = 0;
}


void G17B::GoM10() {// processing an entry 656 566 with the relevant table of ba,ds3
	//if(1)return;
	if (aigstop)return;

#ifdef DODEBUG
	diag = diagbug = 0;
	if (!strcmp(diagband, myband2.band)) {
		cout << myband2.band << "go band2 id=" << myband2.i416 << " nb12=" << genb12.nb12
			<< " nb3=" << genb12.nband3 << endl;
		DebugGetPuz(diagpuz);
		for (int i = 0; i < genb12.nband3; i++)
			if (!strcmp(diagband3, genb12.bands3[i].band))
			cout << genb12.bands3[i].band << " b3 i=" << i << " myband i416=  " << genb12.bands3[i].i416 << endl;
		//diagbug = 2;
		cout << Char27out(p17diag.bf.u32[0]) << endl;
		cout << Char27out(p17diag.bf.u32[1]) << endl;
		cout << Char27out(p17diag.bf.u32[2]) << endl;
	}
	else return;
#endif
	//if (p_cpt2g[0]) return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	p_cpt2g[0] ++;
	p_cpt2g[1] += genb12.nband3;
	cout << myband2.band << " go band2 id=" << myband2.i416 << " nb12=" << genb12.nb12
		<< " nb3=" << genb12.nband3 << " p_cpt2g[0]=" << p_cpt2g[0] << endl;

	//______ true start band 2 expand
	ExpandB2();// expand band2 272 ms
	if (!myband2.nvalidb)return ; // mode 656 only no 5
	b1cptdiag = b1cpt[1];
	GoM10Uas();//expand bands 3  collect UAs 350 ms
	StackUas();
	bin_b1.Copy(myband1);
	bin_b2.Copy(myband2);
	Go2_Ext_Loop();//next is external (outer) loop 
}
#ifdef DEBUGKNOWN

void G17B::GoM10Known() {// processing an entry 656 566 with the relevant table of ba,ds3
	if (aigstop)return;
	p_cpt2g[0] ++;
	if (g17b.debug17 > 1) diag = diagbug = g17b.debug17;
	else diag = diagbug = 0;
	ExpandB1();// expand band1
	ExpandB2();// expand band2
	GoM10Uas();//collect UAs 
	StackUas();
	//GodebugInit(1);
	bin_b1.Copy(myband1);
	bin_b2.Copy(myband2);
	Go2_Ext_Loop(); 
}
#else
void G17B::GoM10Known() {// processing an entry 656 566 with the relevant table of ba,ds3
}
#endif
 
//__________________ bands expansion collect valid 5/6


void G17B::ExpandB2() {
	STD_B1_2 &b = myband2;
	bnua = b.nua;
	for (uint32_t i = 0; i < bnua; i++) btua[i] = (uint64_t)b.tua[i] << 32;
	start_active = (uint64_t)BIT_SET_27 << 32;
	nexp = 6;
	mybi2t = bi2_2;
	memset(bi2_2, 0, sizeof bi2_2);
	myokt = vab_2;
	ExpandOneBand(2);
	nbi2_2 = nmybi2t;
	memcpy(b2cpt, p_cpt, sizeof b2cpt);
	//cout << p_cpt2g[0] << " end expand b2 count=n=" << nbi2_2 << " " << b2cpt[0] << " " << b2cpt[1] << endl;;
	b2count = totb2 = b2cpt[0]; // debugging for max expected
	b.my_bi2 = bi2_2;
	b.nbi2 = nbi2_2;
	b.my_validb = vab_2;
	b.nvalidb = (uint32_t)(p_cpt[0] + p_cpt[1]);
}
void G17B::ExpandB1() {
	STD_B1_2 &b = myband1;
	bnua = b.nua;
	for (uint32_t i = 0; i < bnua; i++) btua[i] = b.tua[i];// just expand to 64
	start_active = BIT_SET_27;
	nexp = 6;
	mybi2t = bi2_1;
	memset(bi2_1, 0, sizeof bi2_1);
	myokt = vab_1;
	ExpandOneBand(1);
	nbi2_1 = nmybi2t;
	memcpy(b1cpt, p_cpt, sizeof b1cpt);
	b.my_bi2 = bi2_1;
	//cout  << " end expand b1 count=n=" << nbi2_1 << " " << b1cpt[0] << " " << b1cpt[1] << endl;;
	b.nbi2 = nbi2_1;
	b.my_validb = vab_1;
	b.nvalidb = (uint32_t)(p_cpt[0] + p_cpt[1]);
}
struct SPOT_E64 {// spots to find band12 valid solutions n clues
	SPOT_E64 * sp;
	uint64_t  all_previous_cells, active_cells;
	uint32_t * start_possibles, n_possibles, ipos, ispot;
	uint64_t * tua;
	uint32_t stack[3], bands[2], missing_clues, nua;
	inline void Copy(SPOT_E64 * old) {
		*this = *old;
		start_possibles += n_possibles;
		ispot++;
		missing_clues--;
		ipos = 0;
		tua += nua;
	}
	inline void AddCellBandStack(int cell, uint32_t * ncb) {
		// if the stack is limit update sn active
		int st = C_stack[cell];
		stack[st]++;
		if (stack[st] > 5) {
			//cout << "stack pleine" << st << endl;
			active_cells &= ~band3xBM[st + 3].u64[0];

		}
		// if the band is limit update sn active
		int b = C_div27[cell];
		bands[b]++;
		if (bands[b] > 5) {// more in mode b 656 only
			//cout << "bande pleine" << b<< endl;
			active_cells &= ~band3xBM[b].u64[0];
		}
	}
	inline void GetOne(uint64_t v) {
		n_possibles = 1;
		bitscanforward64(start_possibles[0], v);
		start_possibles[0] = From_128_To_81[start_possibles[0]];
	}
	inline void AddPossibles(uint64_t v) {
		uint32_t cc;
		while (bitscanforward64(cc, v)) {// look for  possible cells
			register uint64_t bit2 = (uint64_t)1 << cc;
			v ^= (uint64_t)1 << cc;// clear bit
			start_possibles[n_possibles++] = From_128_To_81[cc];
		}
	}
	inline int GetUa(uint64_t v) {
		n_possibles = 0;
		AddPossibles(v);
		return n_possibles;
	}
	uint32_t GetPossibles() {
		// UAs are limited to active cells no empty or single ua
		if (missing_clues < 2) return 0; // minimum to call this process
		uint32_t cells_count[64],
			min_elims = (nua + missing_clues - 1) / missing_clues;
		memset(cells_count, 0, sizeof cells_count);
		uint32_t cc;
		for (uint32_t iua = 0; iua < nua; iua++) {
			register uint64_t Rw = tua[iua] & BIT_SET_2X;
			while (bitscanforward64(cc, Rw)) {// look for  possible cells
				register uint64_t bit2 = (uint64_t)1 << cc;
				Rw ^= (uint64_t)1 << cc;// clear bit
				cells_count[cc]++;
				//if (!cc) cout << "cell 0 pour i=" << iua << " cpt="<< cells_count[cc] << endl;
			}
		}
		//cout << "compte brut critical "<< min_elims << endl;
		//for (int i = 0; i < 64; i++)if (cells_count[i])
		//	cout << From_128_To_81[i] << "\t" << cells_count[i] << endl;

		// collect cells over critical count
		GINT64 tcells[64], temp;
		uint32_t ntcells = 0;
		for (int i = 0; i < 54; i++) {
			register uint32_t my_cell_count = cells_count[C_To128[i]];
			if (my_cell_count >= min_elims) {
				GINT64 & w = tcells[ntcells++];
				w.u32[0] = i;
				w.u32[1] = my_cell_count;
			}
		}
		if (!ntcells) return 0;
		if (ntcells > 1) {// sort in decreasing order
			for (uint32_t i = 0; i < ntcells - 1; i++) {
				for (uint32_t j = i + 1; j < ntcells; j++) {
					if (tcells[i].u64 < tcells[j].u64) {
						temp.u64 = tcells[i].u64;
						tcells[i].u64 = tcells[j].u64;
						tcells[j].u64 = temp.u64;
					}
				}
			}
			//if (ntcells > 64 - missing_clues) ntcells = 64 - missing_clues;
		}
		// load the final table of cells to consider
		//cout << "final count" << endl;
		for (uint32_t i = 0; i < ntcells; i++) {
			start_possibles[i] = tcells[i].u32[0];
			//cout <<i<<"\t"<< tcells[i].u32[0]<<"\t"<< tcells[i].u32[1] <<endl;
		}
		n_possibles = ntcells;
		return ntcells;
	}
	inline int GetLast() {
		n_possibles = 0;
		uint64_t andx = (uint64_t)BIT_SET_2X;
		for (uint32_t iua = 0; iua < nua; iua++) {
			andx &= tua[iua];
			if (!andx)return 0;
		}
		uint32_t cc;
		while (bitscanforward64(cc, andx)) {// look for  possible cells
			register uint64_t bit2 = (uint64_t)1 << cc;
			andx ^= (uint64_t)1 << cc;// clear bit
			start_possibles[n_possibles++] = From_128_To_81[cc];
		}
		return n_possibles;
	}
	inline uint32_t GetPossiblesP() {// using parallel
		uint32_t lim = (nua + missing_clues - 1) / missing_clues,
			lim2 = lim + 5;
		//uint64_t tval[17];
		uint64_t tval[14];
		memset(tval, 0, sizeof tval);
		tval[1] = tua[0] & BIT_SET_2X;
		uint32_t tend = 1; // last used so far
		for (uint32_t iua = 1; iua < nua; iua++) {
			register uint64_t Rw = tua[iua] & BIT_SET_2X,
				Rmore = tval[tend] & Rw;
			if (tend < lim2 && Rmore) tend++;
			switch (tend) {
				//case 15:tval[15]|=tval[14] & Rw;
				//case 14:tval[14]|=tval[13] & Rw;
			case 13:tval[13] |= tval[12] & Rw;
			case 12:tval[12] |= tval[11] & Rw;
			case 11:tval[11] |= tval[10] & Rw;
			case 10:tval[10] |= tval[9] & Rw;
			case 9:tval[9] |= tval[8] & Rw;
			case 8:tval[8] |= tval[7] & Rw;
			case 7:tval[7] |= tval[6] & Rw;
			case 6:tval[6] |= tval[5] & Rw;
			case 5:tval[5] |= tval[4] & Rw;
			case 4:tval[4] |= tval[3] & Rw;
			case 3:tval[3] |= tval[2] & Rw;
			case 2:tval[2] |= tval[1] & Rw;
			}
			tval[1] |= Rw;
		}
		if (tend < lim)return 0;// nothing to do
		n_possibles = 0;
		for (uint32_t i = lim; i < tend; i++) tval[i] &= ~tval[i + 1];
		for (uint32_t i = tend; i >= lim; i--)AddPossibles(tval[i]);
		return n_possibles;
	}
	inline uint32_t GetPossibles2() {// using parallel
		register uint64_t R5 = 0, R4 = 0, R3 = 0, R2 = 0, R1 = 0, R;
		for (uint32_t iua = 1; iua < nua; iua++) {
			R = tua[iua] & BIT_SET_2X;
			R5 |= R4 & R; R4 |= R3 & R; R3 |= R2 & R;	R2 |= R1 & R;
			R1 |= R;
		}
		//if (!R2) return 0;
		R2 &= ~R3; R3 &= ~R4; R4 &= ~R5;
		n_possibles = 0;
		if (R5)					AddPossibles(R5);
		if (R4)					AddPossibles(R4);
		if (R3)					AddPossibles(R3);
		if (R2)					AddPossibles(R2);
		return n_possibles;
	}
	inline uint32_t GetPossibles3() {// using parallel
		register uint64_t R6 = 0, R5 = 0, R4 = 0, R3 = 0, R2 = 0, R1 = 0, R;
		for (uint32_t iua = 1; iua < nua; iua++) {
			R = tua[iua] & BIT_SET_2X;
			R6 |= R5 & R; R5 |= R4 & R; R4 |= R3 & R; R3 |= R2 & R;	R2 |= R1 & R;
			R1 |= R;
		}
		//if (!R3) return 0;
		R3 &= ~R4; R4 &= ~R5; R5 &= ~R6;
		n_possibles = 0;
		if (R6)	AddPossibles(R6);
		if (R5)	AddPossibles(R5);
		if (R4)	AddPossibles(R4);
		if (R3)	AddPossibles(R3);
		return n_possibles;
	}
	inline uint32_t GetPossibles4() {// using parallel
		register uint64_t R6 = 0, R5 = 0, R4 = 0, R3 = 0, R2 = 0, R1 = 0, R;
		for (uint32_t iua = 1; iua < nua; iua++) {
			R = tua[iua] & BIT_SET_2X;
			R6 |= R5 & R; R5 |= R4 & R; R4 |= R3 & R; R3 |= R2 & R;	R2 |= R1 & R;
			R1 |= R;
		}
		//if (!R4) return 0;
		R4 &= ~R5; R5 &= ~R6;
		n_possibles = 0;
		if (R6)	AddPossibles(R6);
		if (R5)	AddPossibles(R5);
		if (R4)	AddPossibles(R4);
		return n_possibles;
	}
	inline uint32_t GetPossibles5() {// using parallel
		register uint64_t R7 = 0, R6 = 0, R5 = 0, R4 = 0, R3 = 0, R2 = 0, R1 = 0, R;
		for (uint32_t iua = 1; iua < nua; iua++) {
			R = tua[iua] & BIT_SET_2X;
			R7 |= R6 & R; R6 |= R5 & R; R5 |= R4 & R; R4 |= R3 & R; R3 |= R2 & R;	R2 |= R1 & R;
			R1 |= R;
		}
		//if (!R5) return 0;
		R4 &= ~R5; R5 &= ~R6; R6 &= ~R7;
		n_possibles = 0;
		if (R7)	AddPossibles(R7);
		if (R6)	AddPossibles(R6);
		if (R5)	AddPossibles(R5);
		return n_possibles;
	}
	void D1() {
		cout << "D1\t" << ispot << "\t" << ipos << endl;
		//cout << Char2Xout(active_cells) << " spot=" << ispot << " pos=" << ipos
		//	<< " npos=" << n_possibles << endl;
	}
	void D2(int all = 0) {
		//cout << Char2Xout(active_cells) << " Shrink active"  << endl;
		cout << Char2Xout(all_previous_cells) << " known nuas=" << nua << endl;
		if (!all) return;
		for (uint32_t iua = 0; iua < nua; iua++)
			cout << Char2Xout(tua[iua]) << endl;
	}

	void D3() {
		cout << Char2Xout(all_previous_cells) << " known" << endl;
		cout << "get possibles  nua=" << nua << "  missing_clues " << missing_clues
			<< "n_possibles" << n_possibles << endl;
	}
	inline void Ddead(uint32_t iua) {
		//			cout <<"\t\t"<< ispot <<" "<<ipos<<" dead branch iua=" << iua << endl;
	}
	inline void Dass(uint32_t iua, uint64_t Ru) {
		//			cout<<"\t\t" << ispot << " " << ipos << " assign iua=" << iua
		//				<< " cell=" << start_possibles[0] << endl;
		//			cout << ispot << " " << ipos << " assign iua=" << iua <<endl
		//				<< Char2Xout(Ru) <<" cell="<< start_possibles[0] << endl;
	}
	inline void DNoMore() {
		cout << Char2Xout(all_previous_cells) << "no more uas "
			<< bands[0] << bands[1] << " "
			<< stack[0] << stack[1] << stack[2] << endl;

	}
};

void G17B::ExpandOneBand(int ib) {
	int notfirst = 0,
		c56 = 011;// collect 5_6
	if (sgo.vx[4]) {// this is pass b only 656
		if (ib == 1)c56 = 010; // must be 6
		else c56 = 1;// must be 5
	}
	memset(p_cpt, 0, sizeof p_cpt);
	uint32_t tclues[12], bufp[1000], lastspot = nexp - 1;
	uint64_t bufua12[30000];
	SPOT_E64 spt[8], *s, *sn;
	nmybi2t = nmyokt = 0;
	s = spt;
	memset(s, 0, sizeof spt[0]);// init the stack status ...
	s->missing_clues = nexp;
	s->active_cells = start_active;// all cells active
	s->start_possibles = bufp;
	s->tua = bufua12;// copy init table to the buffer
	s->nua = bnua;
	memcpy(s->tua, btua, 8 * bnua);
	s->GetUa(s->tua[0] & BIT_SET_2X);
next:
	{// catch and apply cell in bitfields
		uint32_t iw = s->ipos++;
		if (iw >= s->n_possibles)goto back;
		register uint32_t cell = s->start_possibles[iw];
		tclues[s->ispot] = cell;
		register uint64_t bit = (uint64_t)1 << C_To128[cell];
		register  uint64_t filter = s->all_previous_cells | bit,
			ac = s->active_cells ^ bit;
		sn = s + 1; sn->Copy(s); // prepare next spot
		sn->all_previous_cells = filter;
		sn->active_cells = s->active_cells = ac;
		{//level>0 shrink the ua table in new
			register uint32_t nua1 = s->nua, iua;
			sn->nua = 0;
			register uint64_t Ra = sn->active_cells;
			for (iua = 0; iua < nua1; iua++) {
				register uint64_t Ru = s->tua[iua];
				if (Ru&filter)continue;
				Ru &= Ra;
				if (!Ru) goto next;// dead branch  
				register uint64_t cc = _popcnt64(Ru);
				Ru |= (cc << 59);
				AddUA64(sn->tua, sn->nua, Ru);
			}
		}
		if (s->ispot == 1) {// open a new index2
			if (notfirst) {// save previous if active
				BI2 & pr = mybi2t[nmybi2t];
				if (pr.istart != pr.iend) {
					nmybi2t++;
					BI2 & pn = mybi2t[nmybi2t];
					pn.istart = pn.iend = pr.iend;
				}
			}
			notfirst = 1;
			BI2 & pn = mybi2t[nmybi2t];// init the ne status
			pn.bf = sn->all_previous_cells;
			pn.active = sn->active_cells;
			memcpy(pn.tval, tclues, sizeof pn.tval);
		}
		if (!sn->nua)goto no_more_uas;
		else if (sn->missing_clues == 1) { if (!sn->GetLast())goto next; }
		else sn->GetUa(sn->tua[0] & BIT_SET_2X);
		if (s->ispot < lastspot)s++;// {	s++; s->D3(); }

		goto next;
	}
no_more_uas:
	{	BI2 & pi = mybi2t[nmybi2t];
	register uint64_t R0 = sn->all_previous_cells;// ^pi.bf;
	if (s->ispot == 5) {
		if (c56 & 010) {
			myokt[pi.iend++].Enter(R0, &tclues[2]);
			p_cpt[1]++;
		}
	}
	else {//if below 6 loop  for redundant clues
		int tc[64], nt = 0;
		uint64_t tbit[64];
		{	uint32_t register xcell;
		register uint64_t ac = s->active_cells&BIT_SET_2X;
		while (bitscanforward64(xcell, ac)) {// put active cells in table
			uint64_t bit = (uint64_t)1 << xcell;
			ac ^= bit;
			tc[nt] = From_128_To_81[xcell];
			tbit[nt++] = bit;
		}
		if (s->ispot == 4) { //5 clues
			if (c56 & 1) {
				myokt[pi.iend++].Enter(R0, &tclues[2]);
				p_cpt[0]++;
			}
			for (int i6 = 0; i6 < nt; i6++) {
				tclues[5] = tc[i6];
				if (c56 & 010) {
					myokt[pi.iend++].Enter(tbit[i6] | R0, &tclues[2]);
					p_cpt[1]++;
				}
			}
		}
		else if (s->ispot == 3) { // valid 4 clues
			for (int i5 = 0; i5 < nt; i5++) {
				register uint64_t R5 = tbit[i5] | R0;
				tclues[4] = tc[i5];
				if (c56 & 1) {
					myokt[pi.iend++].Enter(R5, &tclues[2]);
					p_cpt[0]++;
				}
				if (c56 & 010)for (int i6 = i5 + 1; i6 < nt; i6++) {
					tclues[5] = tc[i6];
					myokt[pi.iend++].Enter(tbit[i6] | R5, &tclues[2]);
					p_cpt[1]++;
				}
			}
		}
		else if (s->ispot == 2) { // valid 3 clues
			for (int i4 = 0; i4 < nt; i4++) {
				register uint64_t R4 = tbit[i4] | R0;
				tclues[3] = tc[i4];
				for (int i5 = i4 + 1; i5 < nt; i5++) {
					register uint64_t R5 = tbit[i5] | R4;
					tclues[4] = tc[i5];
					if (c56 & 1) {
						myokt[pi.iend++].Enter(R5, &tclues[2]);
						p_cpt[0]++;
					}
					if (c56 & 010) for (int i6 = i5 + 1; i6 < nt; i6++) {
						tclues[5] = tc[i6];
						myokt[pi.iend++].Enter(tbit[i6] | R5, &tclues[2]);
						p_cpt[1]++;
					}
				}
			}
		}
		else if (s->ispot == 1) { // valid 2 clues
			for (int i3 = 0; i3 < nt; i3++) {
				register uint64_t R3 = tbit[i3] | R0;
				tclues[2] = tc[i3];
				for (int i4 = i3 + 1; i4 < nt; i4++) {
					register uint64_t R4 = tbit[i4] | R3;
					tclues[3] = tc[i4];
					for (int i5 = i4 + 1; i5 < nt; i5++) {
						register uint64_t R5 = tbit[i5] | R4;
						tclues[4] = tc[i5];
						if (c56 & 1) {
							myokt[pi.iend++].Enter(R5, &tclues[2]);
							p_cpt[0]++;
						}
						if (c56 & 010)	for (int i6 = i5 + 1; i6 < nt; i6++) {
							tclues[5] = tc[i6];
							myokt[pi.iend++].Enter(tbit[i6] | R5, &tclues[2]);
							p_cpt[1]++;
						}
					}
				}

			}
		}
		}
	}
	}
	goto next;

back:
	if (--s >= spt)goto next;
	// save the last index if 
	BI2 & pr = mybi2t[nmybi2t];
	if (pr.istart != pr.iend) 	nmybi2t++;
}


//____ ___________________________________initial tasks collect uas 12;  guas
void G17B::GoM10Uas() {
	//=========================== collect UAs  old process 
	zh1b_g.modegua = 0;//must be to activate filter in UAs b12 more
	if (genuasb12.Initgen()) return;
	genb12.BuildGang9x3();
	// _____ GUAs 
	zh1b_g.modegua = 1;//must be to kill  filter in GUAs 6_7 more
	genb12.SecondSockets2Setup();// collect GUA2s 
	genb12.SecondSockets3Setup();// collect GUA3s 

	// setupsockets common to all band3
	isguasocket2all.SetAll_0();
	isguasocket3all.SetAll_0();
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3 & b = genb12.bands3[ib3];
		b.guas.isguasocketc2.Convert3X27to81(b.guas.isguasocket2);
		b.guas.isguasocketc3.Convert3X27to81(b.guas.isguasocket3);
		b.guas.isguasocketc2_46.Convert3X27to81(b.guas.isguasocket2_46);
		b.isguasocketc246 = b.guas.isguasocketc2 | b.guas.isguasocketc2_46;
		isguasocket2all |= b.guas.isguasocket2;
		isguasocket3all |= b.guas.isguasocket3;
	}

}


void G17B::StackUas() {// band 3 uas used as gangsters via  C_transpose_d[81]
	STD_B3 wbs;

	// transpose bands 1+2 transpose bands 3
	for (int ib = 0; ib < genb12.nband3; ib++) {// bands 3
		STD_B3 &myb = genb12.bands3[ib];
		wbs.ntua128 = 0;
		memcpy(&genb12.grid0[54], myb.band0, 4 * 27);
		int zt[81];
		for (int i = 0; i < 81; i++) {
			zt[i] = genb12.grid0[C_transpose_d[i]];
			//cout << genb12.grid0[i] +1;
		}
		//cout << endl;
		BANDMINLEX::PERM perm_ret;
		//___ stack 1
		bandminlex.Getmin(zt, &perm_ret);
		int ib1 = perm_ret.i416;
		wbs.InitStack(ib1, zt, perm_ret,0);

		//___ stack 2
		bandminlex.Getmin(&zt[27], &perm_ret);
		int ib2 = perm_ret.i416;
		wbs.InitStack(ib2, &zt[27], perm_ret,1);	
		//___ stack 3
		bandminlex.Getmin(&zt[54], &perm_ret);
		int ib3 = perm_ret.i416;
		wbs.InitStack(ib3, &zt[54], perm_ret,2);

		// check the final status and copy to band3
		myb.ntua128 = 0;
		for (uint32_t i = 0; i < wbs.ntua128; i++) {
			BF128 w = wbs.tua128[i];
			int cc3 = _popcnt32(w.bf.u32[2]);
			if(cc3>3)	myb.tua128[myb.ntua128++] = wbs.tua128[i];
			else if (cc3 ==2) {// be sure to have gua 2 clues in table
				uint64_t cc12 = _popcnt64(w.bf.u64[0]);
				if (cc12 < 12) continue; // must be there
				// find the digits pattern from the current band 3
				int * cur_b3 = myb.band0, wdigs = 0, c27;
				{
					register uint32_t wua = w.bf.u32[2];
					while (bitscanforward(c27, wua)) {// look for  possible cells
						wua ^= 1 << c27;// clear bit
						wdigs |= 1 << cur_b3[c27];
					}
				}
				uint32_t my_i81 = genb12.GET_I81_G2(wdigs, w.bf.u32[2]);
				GUA & g = tguas.tgua_start[my_i81];
				AddUA64(g.tua, g.nua, w.bf.u64[0]);
			}
		}
	}
}
void STD_B3::InitStack(int i16, int * z0, BANDMINLEX::PERM & p, int iband) {
	i416 = i16;
	int dstack = 3 * iband;
	GetUAs();
	// create the cell map in output
	for (int i = 0; i < 3; i++) {
		int vr = 9 * p.rows[i], vr0 = 9 * i;
		for (int j = 0; j < 9; j++)
			map[vr0 + j] = vr + p.cols[j];
	}
	// morph all uas
	for (uint32_t i = 0; i < nua; i++) {
		register int uao = tua[i] & BIT_SET_27;
		BF128  ua; ua.SetAll_0();
		for (int i=0, bit = 1; i < 27; i++, bit <<= 1)
			if (bit & uao) {// set the bit in the band form
				int cell = C_transpose_d[map[i]];// +dstack;
				ua.Set_c(cell+dstack);
			}
		tua128[ntua128++] = ua;
	}
}


//_________________ external loop on bands 1+2 small uas

/*external loop control
 the external loop cut the process using small UAs bands 1+2
 if 'yes' is the part hiiting a ua in a band
 yes1 hit in band 1 ; yes 2 hit in band2
 the process can then be cut in 2 chunks
 yes1 * all2		all1-yes1 * yes2
 this is of interest if the count is significantly lower than
		all1 * all 2
 and if all1*all2 is big enough

 the process can be split several times with different uas

 here, the first chunk yes * all is split again
 the main loop continues to split    all1-yes1 * yes2
 Note the first band can be band 1 or band 2

 Control of the main loop

 EXLRATIO is the minimal wanted reduction of the count
 EXLNLOOP1 the maximum number of steps in the main loop
 EXLNLOOP2 the maximum number of steps in the first chunk split
 EXLBREAK the count minimum to loop
 */
#define EXLRATIO 80
#define EXLNLOOP1 4
#define EXLNLOOP2 3
#define EXLBREAK 4000000


void BINDEXN::Copy(STD_B1_2 & b) {
	ntvb = b.nvalidb;
	for (uint32_t i = 0; i < ntvb; i++) tvb[i] = b.my_validb[i];
	nt2 = b.nbi2;
	memcpy(t2, b.my_bi2, nt2 * sizeof t2[0]);
}
void BINDEXN::Copy(BINDEXN & b) {
	ntvb = b.ntvb;
	memcpy(tvb, b.tvb, ntvb * sizeof tvb[0]);
	nt2 = b.nt2;
	memcpy(t2, b.t2, (nt2 + 1) * sizeof t2[0]);
}

uint64_t G17B::FindSockets(uint64_t active, uint64_t lim) {
	nextl1 = 0, nextl2 = 0;
	uint64_t *t = genuasb12.tua;
	uint32_t n = genuasb12.nua;
	for (uint32_t i = 0; i < n; i++) {
		register uint64_t U = t[i] & active,
			U1 = U & BIT_SET_27, U2 = U & BIT_SET_B2;
		if ((!U1) || (!U2))continue;
		uint64_t n = _popcnt64(U& BIT_SET_27);
		if ((lim == 2 && n < lim) || n == lim) {
			uint32_t i1, aig = 1;

			for (i1 = 0; i1 < nextl1; i1++)
				if (U1 == extl1[i1].bfx) { aig = 0; break; }
			if (aig) {
				if (nextl1 < 6) { // open a new window (limit 10
					i1 = nextl1++; aig = 0;
					extl1[i1].Init((uint32_t)U1, 1);
				}
			}
			if (!aig)if (extl1[i1].ntbfy < 20)
				extl1[i1].tbfy[extl1[i1].ntbfy++] = (uint32_t)(U2 >> 32);

		}
		n = _popcnt64(U& BIT_SET_B2);
		if ((lim == 2 && n < lim) || n == lim) {
			register uint32_t U2a = (uint32_t)(U2 >> 32);
			//if (nt2sk < 5)t2sk[nt2sk++] = (uint32_t)(U2 >> 32);
			uint32_t i2, aig = 1;
			for (i2 = 0; i2 < nextl2; i2++)
				if (U2a == extl2[i2].bfx) { aig = 0; break; }
			if (aig) {
				if (nextl2 < 6) { // open a new window (limit 10
					i2 = nextl2++; aig = 0;
					extl2[i2].Init(U2a, 2);
				}
			}
			if (!aig)if (extl2[i2].ntbfy < 20)
				extl2[i2].tbfy[extl2[i2].ntbfy++] = (uint32_t)U1;
		}
	}
#ifdef DEBUGEXL
	//cout << "FindSockets lim=" << lim << "\tnt1=" << nextl1 << "\tnt2=" << nextl2 << endl;
	//for (uint32_t i = 0; i < nextl1; i++) extl1[i].Debug();
	//for (uint32_t i = 0; i < nextl2; i++) extl2[i].Debug();
#endif
	if (nextl1 | nextl2) {
		return lim - 1;
	}
	else	return 0;
}

uint64_t GetYes_b1( VALIDB *vb, uint32_t nvb, uint32_t bf) {
	uint64_t n = 0;
	register uint64_t F = bf;
	for (uint32_t i = 0; i < nvb; i++) 
		if (F& vb[i].bf)n++;	
	return n;
}
uint64_t GetYes_b1(VALIDB *vb, uint32_t nvb, uint32_t *tbf, uint32_t ntbf) {
	uint64_t n = 0;
	for (uint32_t i = 0; i < nvb; i++) {
		VALIDB &w = vb[i];
		uint32_t aig = 0;
		register uint32_t U =(uint32_t) w.bf;
		for (uint32_t j = 0; j < ntbf; j++) {
			if ((tbf[j] & U)) { aig = 1; break; }
		}
		n += aig;
	}
	return n;
}
uint64_t GetYes_b2(VALIDB *vb, uint32_t nvb, uint32_t bf) {
	uint64_t n = 0;
	register uint64_t F = (uint64_t)bf << 32;
	for (uint32_t i = 0; i < nvb; i++) 
		if (F& vb[i].bf)n++;	
	return n;
}
uint64_t GetYes_b2(VALIDB *vb, uint32_t nvb, uint32_t *tbf, uint32_t ntbf) {
	uint64_t n = 0;
	for (uint32_t i = 0; i < nvb; i++) {
		VALIDB &w = vb[i];
		uint32_t aig = 0;
		register uint32_t U = (uint32_t)(w.bf >> 32);
		for (uint32_t j = 0; j < ntbf; j++) {
			if ((tbf[j] & U)) { aig = 1; break; }
		}
		n += aig;
	}
	return n;
}

void G17B::ExtractMin(uint64_t active, BINDEXN & bin1, BINDEXN & bin2) {
	uint64_t limb1 = _popcnt64(active & BIT_SET_27) - 3,
		limb2 = _popcnt64(active & BIT_SET_B2) - 3;
	for (uint32_t i = 0; i < nextl1; i++) {
		uint32_t orw = 0, *t = extl1[i].tbfy;
		for (uint32_t j = 0; j < extl1[i].ntbfy; j++) {
			orw |= t[j];
		}
		if (_popcnt32(orw) > (uint32_t)limb1) continue;//min 3 clues will end as 100% ratio
		n_yesb1= GetYes_b1(bin1.tvb, bin1.ntvb, extl1[i].bfx);
		if (extl1[i].ntbfy == 1)
			n_yesb2 = GetYes_b2(bin2.tvb, bin2.ntvb, extl1[i].tbfy[0]);
		else n_yesb2 = GetYes_b2(bin2.tvb, bin2.ntvb, extl1[i].tbfy, extl1[i].ntbfy);
		n_nob1 = bin1.ntvb - n_yesb1;
		n_nob2 = bin2.ntvb - n_yesb2;

		uint64_t ratio1 = n_yesb1 * bin2.ntvb + n_yesb2 * n_nob1,
			ratio = (uint64_t)100 * ratio1 / bin1.ntvb / bin2.ntvb;

		if (ratio < minratio) {
			minratio = ratio;
			extl1[i].ratio = (uint32_t)ratio;
			extlw = extl1[i];
			extlw.noxyes = n_yesb2 * n_nob1;
#ifdef DEBUGEXL
			cout << Char27out(extl1[i].bfx) << " t1";
			cout << "\t" << n_nob1 << "\t" << bin1.ntvb << "\t" << n_nob2 << "\t" << bin2.ntvb
				 << "\t" << ratio << endl;
#endif
		}

		}
	for (uint32_t i = 0; i < nextl2; i++) {
		uint32_t orw = 0, *t = extl2[i].tbfy;
		for (uint32_t j = 0; j < extl2[i].ntbfy; j++) {
			orw |= t[j];
		}
		if (_popcnt32(orw) > (uint32_t)limb2) continue;// will end as 100% ratio
		n_yesb2 = GetYes_b2(bin2.tvb, bin2.ntvb, extl2[i].bfx);

		if (extl2[i].ntbfy == 1)
			n_yesb1 = GetYes_b1(bin1.tvb, bin1.ntvb, extl2[i].tbfy[0]);
		else n_yesb1 = GetYes_b1(bin1.tvb, bin1.ntvb, extl2[i].tbfy, extl2[i].ntbfy);
		n_nob1 = bin1.ntvb - n_yesb1;
		n_nob2 = bin2.ntvb - n_yesb2;
		uint64_t ratio1 = n_yesb2 * bin1.ntvb + n_yesb1 * n_nob2,
			ratio = (uint64_t)100 * ratio1 / bin1.ntvb / bin2.ntvb;
		if (ratio < minratio) {
			minratio = ratio;
			extl2[i].ratio = (uint32_t)ratio;
			extlw = extl2[i];
			extlw.noxyes = n_yesb1 * n_nob2;
#ifdef DEBUGEXL
			cout << Char27out(extl2[i].bfx) << " t2";
			cout << "\t" << n_nob1 << "\t" << bin1.ntvb << "\t" << n_nob2 << "\t" << bin2.ntvb
				<< "\t" << ratio << endl;
#endif
		}
	}
	}



void G17B::ExtSplitY(BINDEXN & binw, uint32_t *tbf, uint32_t ntbf,
	uint32_t & activer, int bande) {// source binw exit yes in binw
	uint32_t lim_source = binw.nt2;// table no is also source
	binw.nt2 = binw.ntvb = 0;
	VALIDB *vb = binw.tvb;
	activer = 0; // or status in new vb1 (no) start null
	for (uint32_t i = 0; i < lim_source; i++) { // all index 2
		BI2 wi = binw.t2[i], win = wi;
		win.istart = binw.ntvb;
		for (uint32_t j = wi.istart; j < wi.iend; j++) {
			VALIDB wj = binw.tvb[j];
			register uint64_t Bf = wj.bf;

			for (uint32_t k = 0; k < ntbf; k++) {
				register uint64_t F = tbf[k]; // UA to hit to be yes
				if (bande == 2)F <<= 32;
				if (Bf & F) {// a new "yes" 
					vb[binw.ntvb++] = wj;
					activer |= wj.bf;
					break;// one hit is enough 
				}
			}
		}
		if (binw.ntvb != win.istart) {// new group in "yes"
			win.iend = binw.ntvb;
			binw.t2[binw.nt2++] = win;
		}
	}
}
void  G17B::ExtSplitX(BINDEXN & bin1no, BINDEXN & bin1yes,
	uint32_t bf, uint32_t & activer, int bande) {
	uint32_t lim_source = bin1no.nt2;// table no is also source

	// initial status empty for table yes and no
	bin1no.nt2 = bin1no.ntvb = 0;
	bin1yes.nt2 = bin1yes.ntvb = 0;
	activer = 0; // or status in new vb1 (no) start null
	VALIDB *vb1 = bin1no.tvb, *vb2 = bin1yes.tvb;

	register uint64_t F = bf; // UA to hit to be yes
	if (bande == 2)F <<= 32;
	for (uint32_t i = 0; i < lim_source; i++) { // all index 2
		BI2 wi = bin1no.t2[i], win = wi, win1;
		win.istart = bin1yes.ntvb;
		if (F&wi.bf) {// all the group is "yes"
			uint32_t n = wi.iend - wi.istart;
			win.iend = bin1yes.ntvb + n;
			memcpy(&vb2[bin1yes.ntvb], &vb1[wi.istart], n * sizeof  vb2[0]);
			bin1yes.t2[bin1yes.nt2++] = win;
			bin1yes.ntvb += n;
			continue;
		}
		// must check each validb of the group
		win1 = wi;
		win1.istart = bin1no.ntvb;
		for (uint32_t j = wi.istart; j < wi.iend; j++) {
			VALIDB wj = vb1[j];
			if (F&wj.bf) 	vb2[bin1yes.ntvb++] = wj;
			else { vb1[bin1no.ntvb++] = wj; activer |= wj.bf; }
		}
		if (bin1no.ntvb != win1.istart) {// new group in "no"
			win1.iend = bin1no.ntvb;
			bin1no.t2[bin1no.nt2++] = win1;
		}
		if (bin1yes.ntvb != win.istart) {// new group in "yes"
			win.iend = bin1yes.ntvb;
			bin1yes.t2[bin1yes.nt2++] = win;
		}
	}

}

int CheckBf(BINDEXN & binw, uint64_t bfw) {
	register VALIDB * tvb = binw.tvb;
	for (uint32_t i = 0; i < binw.ntvb; i++) {
		if (tvb[i].bf == bfw)return i;
	}
	return -1;
}

void Filter_initial(BINDEXN & bin1, BINDEXN & bin2) { //must have the valid b1b2 
	cout << "initial check locate puz" << endl;
	uint64_t f = g17b.p17diag.bf.u32[0];
	int kn_ir1 = CheckBf(bin1, f);
	cout << "kn_ir1=" << kn_ir1 << endl;
	if (kn_ir1 < 0)return;
	f = g17b.p17diag.bf.u32[1]; f <<= 32;// switch to mode 64
	int kn_ir2 = CheckBf(bin2, f);
	cout << "kn_ir2=" << kn_ir2 << endl;
	if (kn_ir2 < 0)return;

}
void Filter_Go3(BINDEXN & bin1, BINDEXN & bin2) { //must have the valid b1b2 
	if (g17b.aigstop) return;
	uint64_t f = g17b.p17diag.bf.u32[0];
	int kn_ir1 = CheckBf(bin1, f);
	//cout << "kn_ir1=" << kn_ir1 << endl;
	if (kn_ir1 < 0)return;
	f = g17b.p17diag.bf.u32[1]; f <<= 32;// switch to mode 64
	int kn_ir2 = CheckBf(bin2, f);
	//cout << "kn_ir2=" << kn_ir2 << endl;
	if (kn_ir2 < 0)return;
	cout << "got the right Go3 set" << endl;
	cout << Char2Xout(bin1.tvb[kn_ir1].bf) << " bfb1 ir1=" << kn_ir1 << endl;
	cout << Char2Xout(bin2.tvb[kn_ir2].bf) << " bfb2 ir2=" << kn_ir2 << endl;
	g17b.kn_ir1 = kn_ir1;
	g17b.kn_ir2 = kn_ir2;
	g17b.Go3(bin1, bin2);
	g17b.aigstop = 1;// stop after this
}
void Filter_Go2b_Ext_Loop(uint64_t activeloop, uint32_t mode2) {
	//cout << "entry Filter_Go2b_Ext_Loop mode=" <<mode2 << endl;
	if (g17b.aigstop) return;
	uint64_t f = g17b.p17diag.bf.u32[0];
	int kn_ir1;
	if (mode2 == 1)  kn_ir1=CheckBf(bin_b1yes, f);
	else kn_ir1 = CheckBf(bin_b1, f);
	//cout << "kn_ir1=" << kn_ir1 << endl;
	if (kn_ir1 < 0) return;

	int kn_ir2;
	f = g17b.p17diag.bf.u32[1]; f <<= 32;// switch to mode 64
	if (mode2 == 1)  kn_ir2 = CheckBf(bin_b2, f);
	else kn_ir2 = CheckBf(bin_b2yes, f);
	//cout << "kn_ir2=" << kn_ir2 << endl;
	if (kn_ir2 < 0)return;

	cout << "the 17 searched is in the first lot" << endl;
	g17b.Go2b_Ext_Loop(activeloop, mode2);
	g17b.aigstop = 1;// stop after this

}


void G17B::Go2_Ext_Loop() {	//_____________ outer loop
	loopb1 = 0;
	uint32_t activerb1, activerb2;
	uint64_t activeloop = BIT_SET_2X;
	activerb1= activerb2 = BIT_SET_27;
	//Filter_initial(bin_b1, bin_b2);
	//_________________ external loop
	while (++loopb1 << EXLNLOOP1) {
#ifdef DEBUGEXL
		if (aigstop)cout << "seen aigstop=1" << endl;
		cout << Char2Xout(activeloop) << "========= loop " << loopb1 << endl;
#endif
		if (aigstop)break;
		minratio = extlr.ratio=1000;
		uint64_t ir = FindSockets(activeloop,2);
		if (ir) 	ExtractMin(activeloop, bin_b1, bin_b2);
		if (!ir || minratio > EXLRATIO) {
			ir = FindSockets(activeloop,3);
			if (ir)ExtractMin(activeloop, bin_b1, bin_b2);
			if (!ir || minratio > EXLRATIO) {
				ir = FindSockets(activeloop,4);
				if (ir)ExtractMin(activeloop, bin_b1, bin_b2);
			}
		}
#ifdef DEBUGEXL
		//else 	cout << "finf2 minratio= " << minratio << endl;
		cout << "final selection minratio=" << minratio << endl;
#endif
		if (minratio > EXLRATIO) break;
		else {
			extlr = extlw;
#ifdef DEBUGEXL
			cout << "entry count=" << extlw.noxyes << endl;
			extlw.Debug();
			cout << "shrink mode=" << extlw.mode << endl;
#endif
			if (extlw.mode == 1) {// this is a band1 X band2 Y
				ExtSplitX(bin_b1, bin_b1yes, extlw.bfx, activerb1);
#ifdef DEBUGKNOWN
				if (loopb1 == 1)
					Filter_Go2b_Ext_Loop(BIT_SET_2X | activerb1, 1);
				else  Filter_Go3(bin_b1yes, bin_b2);
#else
				if (loopb1 == 1)
					Go2b_Ext_Loop(BIT_SET_2X | activerb1, 1);
				else  Go3(bin_b1yes, bin_b2);
#endif

				ExtSplitY(bin_b2, extlr.tbfy, extlr.ntbfy, activerb2, 2);
			}
			else {// this is a band2 X band1 Y
				ExtSplitX(bin_b2, bin_b2yes, extlw.bfx, activerb2, 2);
#ifdef DEBUGKNOWN
				if (loopb1 == 1)
					Filter_Go2b_Ext_Loop(BIT_SET_2X, 2);
				else Filter_Go3(bin_b1, bin_b2yes);
#else
				if (loopb1 == 1)
					Go2b_Ext_Loop(BIT_SET_2X, 2);
				else  Go3(bin_b1, bin_b2yes);
#endif
				ExtSplitY(bin_b1,	extlr.tbfy, extlr.ntbfy,  activerb1);
			}
		}
		activeloop = activerb2; activeloop <<= 32; activeloop |= activerb1;
		if (extlr.noxyes < EXLBREAK) break;
	}
#ifdef DEBUGKNOWN
	Filter_Go3(bin_b1, bin_b2);// last call
#else
	Go3(bin_b1, bin_b2);// last call
#endif

}
void G17B::Go2b_Ext_Loop(uint64_t activeloop, uint32_t mode2) {
		//___________init the working areas bin2_b1, bin2_b2
#ifdef DEBUGEXL
		cout << Char2Xout(activeloop) << "entry Go2b_Ext_Loop mode " << mode2 << endl;
#endif
		if (mode2 == 1) {// b1 yes b2 all
			bin2_b1.Copy(bin_b1yes);
			bin2_b2.Copy(bin_b2);
		}
		else {// b1 all b2 yes
			bin2_b1.Copy(bin_b1);
			bin2_b2.Copy(bin_b2yes);
		}
		uint32_t loopb2 = 0;
		uint32_t activerb1, activerb2;
		while (++loopb2 <= EXLNLOOP2) {
#ifdef DEBUGEXL
			if (aigstop)cout << "seen aigstop=1" << endl;
			cout << Char2Xout(activeloop) << "==== loop first " << loopb2 << endl;
#endif
			if (aigstop) return;
			minratio = 1000;
			uint64_t ir = FindSockets(activeloop, 2);
			if (ir)  ExtractMin(activeloop, bin2_b1, bin2_b2);
			if (!ir || minratio > EXLRATIO) {
				ir = FindSockets(activeloop, 3);
				if (ir)ExtractMin(activeloop, bin2_b1, bin2_b2);
				if (!ir || minratio > EXLRATIO) {
					ir = FindSockets(activeloop, 4);
					if (ir)ExtractMin(activeloop, bin2_b1, bin2_b2);
				}
			}
#ifdef DEBUGEXL
			cout << "final selection minratio=" << minratio << endl;
#endif
			if (minratio >= EXLRATIO) break;
			else {
#ifdef DEBUGEXL
				extlw.Debug();
#endif
				if (extlw.mode == 1) {// this is a band1 X band2 Y
					ExtSplitX(bin2_b1, bin2_b1yes, extlw.bfx, activerb1);

#ifdef DEBUGKNOWN
					Filter_Go3(bin2_b1yes, bin2_b2);
#else
					Go3(bin2_b1yes, bin2_b2);
#endif
					ExtSplitY(bin2_b2, extlw.tbfy, extlw.ntbfy, activerb2,2);
				}
				else {// this is a band2 X band1 Y
					ExtSplitX(bin2_b2, bin2_b2yes, extlw.bfx, activerb2, 2);

#ifdef DEBUGKNOWN
					Filter_Go3(bin2_b1, bin2_b2yes);
#else
					Go3(bin2_b1, bin2_b2yes);
#endif
					ExtSplitY(bin2_b1, extlw.tbfy, extlw.ntbfy, activerb1);
				}
			}
			activeloop = activerb2; activeloop <<= 32; activeloop |= activerb1;
			if (extlw.noxyes < EXLBREAK) break;
		}
#ifdef DEBUGKNOWN
		Filter_Go3(bin2_b1, bin2_b2);// last call
#else
		Go3(bin2_b1, bin2_b2);// last call
#endif
	}


//_________"step" _________ process a subset of valid {band1;band2}

	/* the sub lot is made of 2 pieces of the expansion
		usually the smaller piece is in band 1

		here 2 loops, outer band1 inner band 2
		for each loop, the process in cut in sub lots
		  depending on the expansion index (2 common cells)
		  valid bands are split by size
		  a pair of 2 cells index gives a step

		 UAs and GUAs tables are reduced to reach a "step size"
		 common cells and possible cells of the step are identified to optimize the process


	*/


void G17B::Go3_Build_Band1(uint32_t ib1, BINDEXN & binw) {
	nbi5_1 = nbi6_1 = 0;
	BI2 biw = binw.t2[ib1];
	uint32_t id = biw.istart, iend = biw.iend;
	uint64_t	andw = BIT_SET_2X, orw = 0;
	for (uint32_t iv = id; iv < iend; iv++) {
		VALIDB & wv = binw.tvb[iv];
		register uint64_t Bf = wv.bf;
		andw &= Bf;	orw |= Bf;
		if (wv.nval == 3)Z128_5_1[nbi5_1++].bf = Bf; // 5 clues			
		else Z128_6_1[nbi6_1++].bf = Bf;
	}
	fb1 = andw;
	acb1 = orw;
}
void G17B::Go3_Apply_B1_V() {
	ntusb1 = ntusb1_128 = 0;
	register uint64_t * tua = genuasb12.tua;
	register uint32_t nua = genuasb12.nua;
	register uint64_t filter = fb1, Ra = acb1 | BIT_SET_B2;
	//g4t_b1.Shrink(g4t_start, filter);
	for (uint32_t iua = 0; iua < nua; iua++) {
		register uint64_t Ru = tua[iua];
		if (Ru&filter) continue;
		Ru &= Ra;
		if (ntusb1_128 < 128) {
			tusb1_128[ntusb1_128++] = Ru;
		}
		else if(ntusb1<2000)tusb1[ntusb1++] = Ru;
	}
	// Vector for the first 128
	v128uas = maskLSB[ntusb1_128];// Uas vector
	memset(vc128, 255, sizeof vc128);// all bits to 1
	uint32_t cc64;// build cells vectors A
	for (uint32_t i = 0; i < ntusb1_128; i++) {
		register uint64_t Rw = tusb1_128[i];
		while (bitscanforward64(cc64, Rw)) {// look for  possible cells
			Rw ^= (uint64_t)1 << cc64;// clear bit
			vc128[From_128_To_81[cc64]].clearBit(i);
		}
	}
	// Apply vect to b1 5 clues 
	for (uint32_t i5b1 = 0; i5b1 < nbi5_1; i5b1++) {
		ZS128 & w128 = Z128_5_1[i5b1];
		BF128 wvect = v128uas;// only uas not hit by fb1
		register uint64_t W = w128.bf^fb1;
		while (bitscanforward64(cc64, W)) {// look for  possible cells
			W ^= (uint64_t)1 << cc64;// clear bit
			wvect &= vc128[cc64];// b1 cc64<27
		}
		w128.v = wvect;
	}
	// Apply vect to b1 6 clues 
	for (uint32_t i6b1 = 0; i6b1 < nbi6_1; i6b1++) {
		ZS128 & w128 = Z128_6_1[i6b1];
		BF128 wvect = v128uas;// only uas not hit by fb1
		register uint64_t W = w128.bf^fb1;
		while (bitscanforward64(cc64, W)) {// look for  possible cells
			W ^= (uint64_t)1 << cc64;// clear bit
			wvect &= vc128[cc64];
		}
		w128.v = wvect;
	}


}

void G17B::Go3_Build_Band2(uint32_t ib2, BINDEXN & binw) {
	nbi5_2 = nbi6_2 = 0;
	BI2 biw = binw.t2[ib2];
	uint32_t id = biw.istart, iend = biw.iend;
	uint64_t	andw = BIT_SET_2X, orw = 0;
	for (uint32_t iv = id; iv < iend; iv++) {
		VALIDB & wv = binw.tvb[iv];
		register uint64_t Bf = wv.bf;
		andw &= Bf;	orw |= Bf;
		if (wv.nval == 3)Z128_5_2[nbi5_2++].bf = Bf; //  5 clues
		else 	Z128_6_2[nbi6_2++].bf = Bf;
		
	}
	fb2 = andw;
	acb2 = orw;
}
int G17B::Go3_Apply_B2_V() {
	register uint64_t Ra = acb12, filter = fb2;
	register uint64_t * tua = tusb1;
	//g4t_b2.Shrink(g4t_b1, filter);
	ntusb2 = 0;
	for (uint32_t iua = 0; iua < ntusb1; iua++) {
		register uint64_t Ru = tua[iua];
		if (!(Ru&filter)) {
			Ru &= Ra;
			if (!Ru) { return 1; }
			if(ntusb2<1000)tusb2[ntusb2++] = Ru;
		}
	}
	uint32_t cc64;// build cells vectors A
	BF128 wvectc = v128uas;// setup common cells
	register uint64_t F = fb2;
	while (bitscanforward64(cc64, F)) {// look for  possible cells
		F ^= (uint64_t)1 << cc64;// clear bit
		wvectc &= vc128[From_128_To_81[cc64]];
	}

	// Apply vect to b2 5 clues 
	for (uint32_t i5b2 = 0; i5b2 < nbi5_2; i5b2++) {
		ZS128 & w128 = Z128_5_2[i5b2];
		BF128 wvect = wvectc;//uas not common
		register uint64_t W = w128.bf^fb2;
		while (bitscanforward64(cc64, W)) {// look for  possible cells
			W ^= (uint64_t)1 << cc64;// clear bit
			wvect &= vc128[From_128_To_81[cc64]];
		}
		w128.v = wvect;
	}
	// Apply vect to b2 6 clues 
	for (uint32_t i6b2 = 0; i6b2 < nbi6_2; i6b2++) {
		ZS128 & w128 = Z128_6_2[i6b2];
		BF128 wvect = wvectc;// only uas not hit by fb1
		register uint64_t W = w128.bf^fb2;
		while (bitscanforward64(cc64, W)) {// look for  possible cells
			W ^= (uint64_t)1 << cc64;// clear bit
			wvect &= vc128[From_128_To_81[cc64]];
		}
		w128.v = wvect;
	}

	//_________ setup the brute force start
	nclues_step = 0;
	uint64_t w = fb12;
	uint32_t xcell;
	while (bitscanforward64(xcell, w)) {
		w ^= (uint64_t)1 << xcell;
		tclues[nclues_step++] = From_128_To_81[xcell];
	}
	tcluesxy = &tclues[nclues_step];
	zh2b_i1.ValidXY_Step(tclues, nclues_step);

	return 0;
}

void G17B::Go3(BINDEXN & bin1, BINDEXN & bin2) {
	p_cpt2g[2]++;
#ifdef DEBUGL1L2 
	cout << p_cpt2g[2] << " go3\t" << bin1.ntvb << "\t " << bin2.ntvb << endl;
#endif
	if (aigstop) return;
	if ((!bin1.ntvb) || (!bin2.ntvb)) return;


	//______________________________________ loop b1 main loop
	for (uint32_t ib1 = 0; ib1 < bin1.nt2; ib1++) {
		if (aigstop) return;
		Go3_Build_Band1(ib1, bin1);
		Go3_Apply_B1_V();
#ifdef DEBUGL1 
		if (ib1 != DEBUGL1) continue;
		cout << Char27out((uint32_t)fb1) << " ntusb1=" << ntusb1
			<< " ntusb1_128=" << ntusb1_128
			<< "\tn5=" << nbi5_1 << " n6=" << nbi6_1 << "\t" << _popcnt64(acb1)
			<< "\tdebufl1 ib1=" << ib1 << endl;
		//for (uint32_t iu = 0; iu < 32; iu++)
		//	cout << Char2Xout(tusb1_128[iu]) << " " << iu << endl;
		//for (uint32_t iv = 0; iv < 54; iv++)
		//	cout << Char32out(vc128[iv].bf.u32[0]) << " v=" << iv << endl;

#endif
		//for (uint32_t i61 = 0; i61 < nbi6_1; i61++) Z128_6_1[i61].dump();
		tguas.ApplyLoopB1();
		for (uint32_t ib2 = 0; ib2 < bin2.nt2; ib2++) {
			if (aigstop) return;
			//if (ib2 != 32) continue;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			Go3_Build_Band2(ib2, bin2);
			fb12 = fb1 | fb2;
			acb12 = acb1 | acb2;
			if (Go3_Apply_B2_V()) continue;
#ifdef DEBUGL1 
			if (ib2 != DEBUGL2) continue;
			cout << Char2Xout(fb2) << "\tn5=" << nbi5_2 << " n6=" << nbi6_2
				<< "\t" << _popcnt64(acb2) << "\tib2=" << ib2 << endl;
#endif
			//for (uint32_t i52 = 0; i52 < nbi5_2; i52++) Z128_5_2[i52].dump();
			Do128uas();
		}

	}
}

//___________________ gangster specific

void TGUAS::ApplyLoopB1() {
	register uint64_t Bf = g17b.fb1;
	register uint64_t Ac = g17b.acb1 | BIT_SET_B2;

	for (uint32_t i = 0; i < 162; i++) {
		GUA & wg = tgua_start[i], &wgd = tgua_b1[i];
		wgd.nua=0;
		for (uint32_t j = 0; j < wg.nua; j++) {// apply new subsets
			register uint64_t Ua = wg.tua[j];
			if (Ua & Bf) continue;
			Ua &= Ac;// could be empty
			uint64_t cc = _popcnt64(Ua);
			wgd.Adduacheck(Ua | (cc << 59)); // no redundancy
		}
	}
	// reduce the tua128 table in bands 3
	for (int ib = 0; ib < genb12.nband3; ib++) {// bands 3
		STD_B3 &myb = genb12.bands3[ib];
		myb.ntua128_b1 = 0;
		for (uint32_t i = 0; i < myb.ntua128; i++) {
			BF128 w = myb.tua128[i];
			if (!(w.bf.u64[0] & Bf)) {
				myb.tua128_b1[myb.ntua128_b1++] =w;
			}
		}
		//cout << "ib=" << ib << " n128=" << myb.ntua128_b1 << endl;
	}
	//for (uint32_t i = 0; i < 162; i++)tgua_b1[i].Debug(1);
	

}
void TGUAS::ApplyLoopB2() {
	// relay table per size
	uint32_t ntt[38]; // count for tt
	BF128 tt[38][2000];// guas 54 ua12 plus code 0_161 for the gangster
	memset(ntt, 0, sizeof ntt);
	register uint64_t Bf = g17b.fb12c;
	register uint64_t Ac = g17b.acb12c;

	//cout << Char2Xout(Bf) << " Bf" << endl;
	//cout << Char2Xout(Ac) << " Ac" << endl;

	for (uint32_t igu = 0; igu < 162; igu++) {
		GUA & wg = tgua_b1[igu];
		uint64_t tua[64];
		uint32_t nua=0 ;
		for (uint32_t j = 0; j < wg.nua; j++) {// apply new subsets
			register uint64_t Ua = wg.tua[j];
			if (Ua & Bf) continue;
			Ua &= Ac;// could be empty
			if (!Ua) {
				tua[0] = 0;
				nua = 1;
				break;
			}
			tua[nua++]=Ua;
		}
		if (nua) {	// split per size in 54 + index 0-162
			BF128 w;
			w.bf.u64[1] =igu;
			for (uint32_t j = 0; j < nua; j++) {
				register uint64_t Ua = tua[j],
					cc = _popcnt64(Ua);
				if (cc > 18) continue; // should not be
				if (igu > 80) {// this is a g3
					cc += 19;
				}
				w.bf.u64[0] = (Ua & BIT_SET_27) | ((Ua & BIT_SET_B2) >> 5);
				tt[cc][ntt[cc]++] = w;
			}
		}
	}
	//_____________ create vectors
	nvg2 = 0;
	for (uint32_t i = 0; i < 19; i++) {// all guas2
		uint32_t n = ntt[i];
		BF128 * tw = tt[i];
		for (uint32_t j = 0; j < n; j++) {
			BF128 w = tw[j];
			AddVG2(w.bf.u64[0], w.bf.u32[2]);
			//if (i < 7)cout<<i<<" " << Char54out(w.bf.u64[0]) << " bb " << w.bf.u32[2]<<" "<< nvg2 << endl;
			if (nvg2 >= 1024) break;
		}
		if (nvg2 >= 1024) break;
	}
	nvg3 = 0;
	for (uint32_t i = 19; i < 38; i++) {// all guas2
		uint32_t n = ntt[i];
		BF128 * tw = tt[i];
		for (uint32_t j = 0; j < n; j++) {
			BF128 w = tw[j];
			AddVG3(w.bf.u64[0], w.bf.u32[2]);
			if (nvg3 >= 512) break;
		}
		if (nvg3 >= 512) break;
	}
	// reduce the tua128 table in bands 3
	for (int ib = 0; ib < genb12.nband3; ib++) {// bands 3
		STD_B3 &myb = genb12.bands3[ib];
		myb.ntua128_b2 = 0;
		for (uint32_t i = 0; i < myb.ntua128_b1; i++) {
			BF128 w = myb.tua128_b1[i];
			if (!(w.bf.u64[0] & Bf)) {
				myb.tua128_b2[myb.ntua128_b2++] = w;
			}
		}
		//cout << "ib=" << ib << " n128_b2=" << myb.ntua128_b2 << endl;
	}

}
void TGUAS::ApplyG2() {
	//cout << "ApplyG2() nclues=" << g17b.nclues << endl;
	g2ok.SetAll_0();
	nb64_1 = (nvg2 + 63) >> 6;
	for (uint32_t iv = 0; iv < nb64_1; iv++) {
		TVG64 &vv = tvg64g2[iv];
		register uint64_t V = vv.v;
		for (int j = 0; j < g17b.nclues; j++)
			V &= vv.cells[g17b.tcluesxy[j]];
		register uint32_t cc64;
		while (bitscanforward64(cc64, V)) {
			V ^= (uint64_t)1 << cc64;// clear bit
			g2ok.Set(vv.ti162[cc64]);
		}
	}
}
int TGUAS::ApplyG3() {
	g3ok.SetAll_0();
	nb64_2 = (nvg3 + 63) >> 6;
	for (uint32_t iv = 0; iv < nb64_2; iv++) {
		TVG64 &vv = tvg64g3[iv];
		register uint64_t V = vv.v;
		for (int j = 0; j < g17b.nclues; j++)
			V &= vv.cells[g17b.tcluesxy[j]];
		register uint32_t cc64;
		while (bitscanforward64(cc64, V)) {
			V ^= (uint64_t)1 << cc64;// clear bit
			g3ok.Set(vv.ti162[cc64]-81);
		}
	}
	return g3ok.isNotEmpty();
}



//______________main loop 128 first uas
void G17B::Do128uas() {//apply first 128 filter
	// reset more uas tables
	moreuas_12_13.Init();	moreuas_14.Init();
	moreuas_15.Init();	moreuas_AB_small.Init();
	moreuas_AB.Init();	moreuas_AB_big.Init();
	nza = nbi5_1;
	if (nza) { // 5 b1  6 b2  
		za = Z128_5_1;
		zb = Z128_6_2;
		nzb = nbi6_2;
		Do128Chunk();
	}
	nza = nbi6_1;
	if (nza) { // 6 b1  5 b2  
		za = Z128_6_1;
		zb = Z128_5_2;
		nzb = nbi5_2;
		Do128Chunk();
	}
}

void G17B::Do128Chunk() {

	if (!nzb) return;
	if (aigstop) return;
	p_cpt2g[27]++;
	p_cpt2g[54] += (nza * nzb);
#ifdef DEBUGCHUNK
	cout << "Do128Chunk() " << nza << " " << nzb << endl;
	cout << " a list" << endl;
	if (DEBUGCHUNK == p_cpt2g[27]) {
		for (uint32_t i = 0; i < nza; i++) {
			cout << Char2Xout(za[i].bf) << " ia=" << i << endl;
			if (i == 4)				za[i].v.Print("v");
		}
		cout << " b list" << endl;
		for (uint32_t i = 0; i < nzb; i++) {
			cout << Char2Xout(zb[i].bf) << " ib=" << i << endl;
			if (i ==78)				zb[i].v.Print("v");
		}
		v128uas.Print(" base vector");
		for (int i = 0; i < 54; i++) {
			cout << i << "\t";
			vc128[i].Print("v");
		}
	}

#endif
	n_to_clean = 0;
	if ((nza * nzb) < 5000) {
		if (nza < nzb)Do128Go(za, zb, nza, nzb);
		else Do128Go(zb, za, nzb, nza);
		if (n_to_clean) CleanAll();
		return;
	}
	// cut in chunks max Xchunk Ychunk
	ZS128 * z1 = za, *z2 = zb;
	uint32_t n1 = nza, n2 = nzb;
	uint32_t  ideb2 = 0, iend2 = YCHUNK64;
	if (iend2 > n2)iend2 = n2;
	while (ideb2 < n2) { //Y chunk
		uint32_t ny = iend2 - ideb2;
		uint32_t ideb1 = 0, iend1 = XCHUNK64;
		if (iend1 > n1)iend1 = n1;

		while (ideb1 < n1) {// X chunk
			uint32_t nx = iend1 - ideb1;
			if (nx < ny)Do128Go(&z1[ideb1], &z2[ideb2], nx, ny);
			else Do128Go(&z2[ideb2], &z1[ideb1], ny, nx);
			ideb1 = iend1; iend1 += XCHUNK64;
			if (iend1 > n1)iend1 = n1;
		}
		ideb2 = iend2; iend2 += YCHUNK64;
		if (iend2 > n2)iend2 = n2;
	}
	if (n_to_clean) CleanAll();
}
void G17B::Do128Go(ZS128 * a, ZS128 * b, uint32_t na, uint32_t nb) {
	if (aigstop)return;
	//check a matrix band 1 band2 for potential 2 bands valid 11 clues
	register ZS128 * Ra = &a[na - 1];
	register uint64_t * Rs = &to_clean[n_to_clean];
	for (; Ra >= a; Ra--) {
		register ZS128 * Rb = &b[nb - 1];
		BF128 va = (*Ra).v;
		register uint64_t bfa = (*Ra).bf;
		for (; Rb >= b; Rb--)
			if ((Rb->v&va).isEmpty()) 		*Rs++ = bfa | Rb->bf;
	}
	n_to_clean = Rs - to_clean;
	if (n_to_clean > 10000)CleanAll();
}


//___________ process potential valid bands 1+2 called once per step

void G17B::CleanAll() {
	//cout << " entry clean " << n_to_clean << " p_cpt2g[4]=" << p_cpt2g[4] << endl;
	p_cpt2g[5]++;	p_cpt2g[7] += n_to_clean;
	if (CleanAll1()) return;// filter all uas
	n_to_clean = n_to_clean2; // remaining valid
	if (CleanAllFifo()) return;// filter Fifo tables 
	n_to_clean = n_to_clean2; // remaining valid
	//cout << " clean2 " << n_to_clean  << endl;
	p_cpt2g[6]++;	p_cpt2g[8] += n_to_clean;
	// no more uas 12 filters, stack control and sockets settings
#ifdef DEFPHASE
	if (DEFPHASE == -2) {
		n_to_clean = 0; return;
	}
#endif

	Clean_2();
	n_to_clean = 0;
}

int G17B::CleanAll1_a() {// <= 10 remaining uas 
	if (!ntusb_clean)	{
		n_to_clean2= n_to_clean;
		return 0;
	}
	for (uint64_t i = 0; i < n_to_clean; i++) {// direct check
		register uint64_t Bf = to_clean[i];
		int aig = 1;
		for (uint32_t iua = 0; iua < ntusb_clean; iua++)
			if (!(Bf & tusb_clean[iua])) {	aig = 0; break;	}
		if (aig)  to_clean[n_to_clean2++] = Bf;
	}
	return (n_to_clean2==0);
}
int G17B::CleanAll1_64(uint64_t bfc) {// <= 64 uas
	uint64_t ve= maskLSB[ntusb_clean].u64[0],vc[54];
	memset(vc, 255, sizeof vc);// all bits to 1
	uint32_t cc64;// build cells vectors A
	uint64_t biti = 1;
	for (uint32_t i = 0; i < ntusb_clean; i++, biti <<= 1) {
		register uint64_t Rw = tusb_clean[i] & BIT_SET_2X;
		while (bitscanforward64(cc64, Rw)) {// look for  possible cells
			Rw ^= (uint64_t)1 << cc64;// clear bit
			vc[From_128_To_81[cc64]] ^= biti;
		}
	}
	{// apply common cells to ve
		register uint64_t Bf = bfc;
		while (bitscanforward64(cc64, Bf)) {
			Bf ^= (uint64_t)1 << cc64;
			ve &= vc[From_128_To_81[cc64]];
		}
	}
	for (uint64_t i = 0; i < n_to_clean; i++) { 
		register uint64_t Bf = to_clean[i]^bfc,
			v=ve;// initial vector
		while (bitscanforward64(cc64, Bf)) {
			Bf ^= (uint64_t)1 << cc64;
			v &= vc[From_128_To_81[cc64]];
		}
		if (!v) to_clean[n_to_clean2++] = to_clean[i];
	}
	return (n_to_clean2 == 0);
}
int G17B::CleanAll1_128(uint64_t bfc, uint32_t nu) {
	BF128 ve , vc[54];
	ve = maskLSB[nu];// Uas vector
	memset(vc, 255, sizeof vc);// all bits to 1
	uint32_t cc64;
	for (uint32_t i = 0; i < nu; i++) {
		register uint64_t Rw = tusb_clean[i] & BIT_SET_2X;
		while (bitscanforward64(cc64, Rw)) {// look for  possible cells
			Rw ^= (uint64_t)1 << cc64;// clear bit
			vc[From_128_To_81[cc64]].clearBit(i);
		}
	}
	   
	{// apply common cells to ve
		register uint64_t Bf = bfc;
		while (bitscanforward64(cc64, Bf)) {
			Bf ^= (uint64_t)1 << cc64;
			ve &= vc[From_128_To_81[cc64]];
		}
	}
	for (uint64_t i = 0; i < n_to_clean; i++) {
		register uint64_t Bf = to_clean[i]^bfc;
		BF128 	v=ve;// initial vector
		while (bitscanforward64(cc64, Bf)) {
			Bf ^= (uint64_t)1 << cc64;
			v &= vc[From_128_To_81[cc64]];
		}
		if (v.isEmpty()) to_clean[n_to_clean2++] = to_clean[i];
	}
	return (n_to_clean2 == 0);
}
int G17B::CleanAll1() {
	diagbugclean = 0;
#ifdef DEBUGSTEP
	if (p_cpt2g[4] == DEBUGSTEP) {
		cout << " entry clean " << n_to_clean << " p_cpt2g[4]=" << p_cpt2g[4] << endl;
	}
#endif
#ifdef DEBUGKNOWN
	//cout << " entry clean " << nw << endl;
	int aig = 1;
	for (uint64_t i = 0; i < n_to_clean; i++) {
		uint64_t bf = to_clean[i];
		if (bf == p17diag.bf.u64[0]) {
			cout << Char2Xout(bf)
				<< "\t\tclean all seen bf i=" << i << " forced to one clean " << endl;
			cout << Char2Xout(fb1 | fb2) << " step applied" << endl;
			to_clean[0] = bf;
			n_to_clean = 1;
#ifdef DEBUGONECLEAN
			diagbugclean = 1;
#endif
			aig = 0;
			break;
		}
	}
	if (aig) return 1;
#endif	

	register uint64_t And = BIT_SET_2X, Or = 0;
	for (uint64_t i = 0; i < n_to_clean; i++) {// setup and / or for this set
		register uint64_t bf = to_clean[i];
		And &= bf; Or |= bf;
	}
	{ //  collect still valid uas and check dead branch  
		ntusb_clean = 0;
		for (uint32_t iua = 0; iua < ntusb2; iua++) {
			register uint64_t Ru = tusb2[iua];
			if (!(Ru&And)) {
				Ru &= Or;
				if (!Ru) 	return 1;
				else tusb_clean[ntusb_clean++] = Ru;
			}
		}
	}
	while (1) {// loop if > 128 uas in tusb2
		n_to_clean2 = 0;
		if (ntusb_clean <= 10 || n_to_clean < 10) return CleanAll1_a();
		if (ntusb_clean <= 64) return CleanAll1_64(And);
		if (ntusb_clean <= 128) return CleanAll1_128(And, ntusb_clean);
		// more than 128, do the first 128
		if (CleanAll1_128(And, 128)) return 1; // at least one ua not hit
		n_to_clean = n_to_clean2; // new status for clean (>0)
		And = BIT_SET_2X; Or = 0;// redo and/or 
		for (uint64_t i = 0; i < n_to_clean; i++) {
			register uint64_t bf = to_clean[i];
			And &= bf; Or |= bf;
		}
		for (uint32_t i = 128; i < ntusb_clean; i++)
			tusb_clean[i - 128]= tusb_clean[i] & Or;
		ntusb_clean -= 128;
	}
}

int G17B::CleanAllFifo() {
	n_to_clean2=0;
	for (uint64_t i = 0; i < n_to_clean; i++) {
		register uint64_t bf = to_clean[i];
		if (moreuas_12_13.Check(bf))continue;
		if (moreuas_14.Check(bf))continue;
		if (moreuas_15.Check(bf))continue;
		if (moreuas_AB_small.Check(bf))continue;
		if (moreuas_AB.Check(bf)) continue;
		if (moreuas_AB_big.Check(bf)) continue;
		to_clean[n_to_clean2++] = bf;
	}
	return (n_to_clean2 == 0);
}

void G17B::Clean_2() {
	if (diagbugclean) cout << "entry clean 2" << endl;
	register uint64_t And = BIT_SET_2X, Or = 0; // re do and or
	for (uint64_t i = 0; i < n_to_clean2; i++) {// setup and / or for this set
		register uint64_t bf = to_clean[i];
		And &= bf; Or |= bf;
	}
	fb12c = And;	acb12c = Or;
	tguas.ApplyLoopB2();// create vectors sockets 2 3
	//g4t_clean.Shrink(g4t_b2, fb12c);

	//cout << "dump first 64 gua2 nvg2="<<tguas.nvg2 << endl;
	//tvg64g2[0].Dump();
	//cout << "dump first 64 gua3 nvg2=" << tguas.nvg3 << endl;
	//tvg64g3[0].Dump();

	{	//_________ setup the brute force start
		nclues_step = 0;
		register uint32_t xcell;
		stack_count_step.u64 = 0;
		while (bitscanforward64(xcell, And)) {
			And^= (uint64_t)1 << xcell;
			uint32_t cell = From_128_To_81[xcell];
			tclues[nclues_step++] = cell;
			stack_count_step.u16[C_stack[cell]]++;
		}
		tcluesxy = &tclues[nclues_step];

		// apply the common cells G2
		tguas.nb64_1 = (tguas.nvg2 + 63) >> 6;
		for (uint32_t iv = 0; iv < tguas.nb64_1; iv++) {
			tvg64g2[iv].ApplyXYcells(tclues, nclues_step);
		}
		tguas.nb64_2 = (tguas.nvg3 + 63) >> 6;
		for (uint32_t iv = 0; iv < tguas.nb64_2; iv++) {
			tvg64g3[iv].ApplyXYcells(tclues, nclues_step);
		}
	}
	//cout << Char2Xout(fb12c) << " fb12c nclues_step="<< nclues_step
		//<< " nb64_1=" << tguas.nb64_1 << " nb64_2=" << tguas.nb64_2
		//<< endl;
	// chek the minimum clues in guas2 guas3
	//n_to_clean2 = 0;
	moreuasxy.Init();
	for (uint64_t iclean = 0; iclean < n_to_clean; iclean++) {// now
		if (aigstop) {
			cout << "aigstop=1 wbf=" << wb12bf << endl;
			return;
		}
		wb12bf = to_clean[iclean];
		if (moreuasxy.Check(wb12bf))continue;
		Clean_3();
	}

}

//_________________________________ apply guas in band3


uint32_t STD_B3::CleanG2() {// critical code guas2
	memset(&smin, 0, sizeof smin);
	BF128 g2 = tguas.g2ok&guas.isguasocketc2;
	int cc81;
	while ((cc81 = g2.getFirst128()) >= 0) {
		g2.Clear(cc81);
		Insert2(cc81);
	}
	sminr = smin;
	smin.SetMincount();
	//cout << "min=" << smin.mincount << endl;
	if (smin.mincount > 6)return 0;
	stack_count.u64 = g17b.stack_count.u64 + smin.Count_per_stack().u64;
	if (stack_count.u16[0] > 6 || stack_count.u16[1] > 6 ||
		stack_count.u16[2] > 6) return 0; // not ok
	//smin.Status("clean0");
	return smin.mincount;
}
uint32_t STD_B3::CleanG3() {// guas 3
	smin = sminr;
	BF128 g3 = tguas.g3ok & guas.isguasocketc3;
	int cc81;
	while ((cc81 = g3.getFirst128()) >= 0) {
		g3.Clear(cc81);
		Insert3(cc81);
	}
	smin.SetMincount();

	if (smin.mincount > 6)return 0;
	stack_count.u64 = g17b.stack_count.u64 + smin.Count_per_stack().u64;
	if (stack_count.u16[0] > 6 || stack_count.u16[1] > 6 ||
		stack_count.u16[2] > 6) return 0; // not ok
	return smin.mincount;
}
void STD_B3::CleanBuildIfOf() {
	register uint32_t Fstk = BIT_SET_27;
	if (stack_count.u16[0] == 6) Fstk ^= 07007007;
	if (stack_count.u16[1] == 6) Fstk ^= 070070070;
	if (stack_count.u16[2] == 6) Fstk ^= 0700700700;
	wactive0 = fstk = Fstk;
	g17b.nuasb3_1 = g17b.nuasb3_2 = 0;
	//nu_if = nu_of = 0;
	// load uas guas 2 guas 3 in field
	int cc81;
	{
		BF128 g2 = tguas.g2ok & guas.isguasocketc2;
		while ((cc81 = g2.getFirst128()) >= 0) {
			g2.Clear(cc81);
			g17b.uasb3_1[g17b.nuasb3_1++] = guas.ua_pair[cc81];
		}
	}
	{
		BF128 g3 = tguas.g3ok&guas.isguasocketc3;
		while ((cc81 = g3.getFirst128()) >= 0) {
			g3.Clear(cc81);
			g17b.uasb3_1[g17b.nuasb3_1++] = guas.ua_triplet[cc81];
		}
	}
	// build out field and And out
	register uint32_t If = smin.critbf; // in field 
	register uint32_t andout = Fstk & (~If);
	for (uint32_t i = 0; i < ntua128_b2; i++) {
		BF128 w = tua128_b2[i];
		if (w.bf.u64[0] & g17b.wb12bf) continue;
		register uint32_t u3 = w.bf.u32[2];
		if (u3&If)g17b.uasb3_1[g17b.nuasb3_1++] = u3;
		else {
			g17b.uasb3_2[g17b.nuasb3_2++] = u3 & Fstk;
			andout &= u3;
		}
	}
	BF128 w = tguas.g2ok &guas.isguasocketc2_46;
	int i81;
	while ((i81 = w.getFirst128()) >= 0) {
		w.Clear(i81);
		register uint32_t u3 = guas.ua_pair[i81];
		if (u3&If)g17b.uasb3_1[g17b.nuasb3_1++] = u3;
		else {
			g17b.uasb3_2[g17b.nuasb3_2++] = u3 & Fstk;
			andout &= u3;
		}
	}
	for (uint32_t i = 0; i < nua; i++) {
		register uint32_t u3 = tua[i];
		if (u3&If)g17b.uasb3_1[g17b.nuasb3_1++] = u3;
		else {
			g17b.uasb3_2[g17b.nuasb3_2++] = u3 & Fstk;
			andout &= u3;
		}
	}

	and_out = andout;
	if (g17b.diagbugclean) {
		cout << "end build if og nif= "<< g17b.nuasb3_1 << " n_of=" << g17b.nuasb3_2 << endl;
		//cout << Char27out(smin.critbf) << "critbf" << endl;
		for(uint32_t i=0;i< g17b.nuasb3_1;i++)
			cout << Char27out(g17b.uasb3_1[i]) << endl;
	}
}


void G17B::Clean_3() {// this is for a given band 1+2
	if (diagbugclean) cout << "entry clean 3" << endl;

	p_cpt2g[42]++;
	//if (p_cpt2g[42]>300) return;//<<<<<<<<<<<<<<<<<<<<<<

	if (aigstop) return;
	clean_valid_done = 0;
	aigstopxy = 0;
	nclues = 0;
	{
		register uint64_t w = wb12bf ^ fb12c;
		register uint32_t xcell, cell;
		stack_count = stack_count_step;
		while (bitscanforward64(xcell, w)) {
			w ^= (uint64_t)1 << xcell;
			cell = From_128_To_81[xcell];
			tcluesxy[nclues++] = cell;
			stack_count.u16[C_stack[cell]]++;
		}
	}
	if (stack_count.u16[0] > 6 || stack_count.u16[1] > 6 ||
		stack_count.u16[2] > 6) return;
	p_cpt2g[43]++;

	if (genb12.nband3 > 10) {
		p_cpt2g[28]++;		
		clean_valid_done = 1;
		myua = zh2b[0].ValidXY(tclues, nclues + nclues_step);
		if (myua) { NewUaB12();		return; }
	}
#ifdef DEFPHASE
	if (DEFPHASE == -3)  return;	
#endif

	tguas.ApplyG2();// all guas2 
	//tguas.g2ok.Print("g2ok");
	uint32_t tvb3[256], nvb3 = 0, minf=10,minr; //for Bands 3 still valid
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		minr = genb12.bands3[ib3].CleanG2();
		if (minr) { 
			tvb3[nvb3++] = ib3; 
			if (minr<minf) minf = minr;
		}

	}
	if (!nvb3) return;
	//cout << Char2Xout(wb12bf) << " wb12bf p_cpt2g[50]=" << p_cpt2g[50]
	//	<< " nb3=" << genb12.nband3 << endl;
	//cout << "after clean G2 nb3=" << nvb3 << " minf=" << minf << endl;
	/*
	if (minf < 5) {// do check bands 1 2 here
		if (!clean_valid_done) {
			p_cpt2g[44]++;
			clean_valid_done = 1;
			myua = zh2b[0].ValidXY(tclues, nclues + nclues_step);
			if (myua) {
				NewUaB12();		aigstopxy = 1; return;
			}
		}
	}
	*/
	p_cpt2g[45]++;
	if (tguas.ApplyG3()) {
		uint32_t n = nvb3;
		nvb3 = 0;
		for (uint32_t iw = 0; iw < n; iw++) {
			uint32_t ib3 = tvb3[iw];
			minr = genb12.bands3[ib3].CleanG3();
			if (minr) {
				tvb3[nvb3++] = ib3;
				if (minr < minf) minf = minr;
			}
		}
	}
	if (!nvb3) return;
	//cout << "after clean G3 nb3=" << nvb3 << " minf=" << minf << endl;
	
	if (minf < 5 && (!clean_valid_done)) {// do check bands 1 2 here
		p_cpt2g[46]++;
		clean_valid_done = 1;
		myua = zh2b[0].ValidXY(tclues, nclues + nclues_step);
		if (myua) {
			NewUaB12();		aigstopxy = 1; return;
		}
	}
	
	p_cpt2g[47]++;		p_cpt2g[53] += nvb3;
#ifdef DEFPHASE
	if (DEFPHASE == -4)  return;
#endif
	if (diagbugclean) cout << Char2Xout(wb12bf) << " wb12bf p_cpt2g[47]=" << p_cpt2g[47]
		<< " nb3=" << nvb3 << endl;
	//for (uint32_t iw = 0; iw < nvb3; iw++)// pick up uas b3 to use
	//	genb12.bands3[tvb3[iw]].CleanBuildIfOf();
	zhou[0].PartialInitSearch17(tclues, nclues + nclues_step);
	//__________no more guas2 guas3 process bands 3	

	for (uint32_t iw = 0; iw < nvb3; iw++) {
		uint32_t ib3 = tvb3[iw];
		p_cpt2g[48]++;// debugging one band
		GoB3(genb12.bands3[ib3]);
		if (aigstopxy)break;// added 
	}

}



//______________________ start process final b3
void G17B::GoB3(STD_B3 & b) {
	if (diagbugclean)cout << "gob3 p_cpt2g[48]=" << p_cpt2g[48] << endl;
#ifdef DEBUGB3
	int locdiag = 0;
	if (p_cpt2g[48] == DEBUGB3) locdiag = 1;
	if (locdiag) {
		cout << b.band << "band3 in diag" << endl;
		b.smin.Status(" ");
	}
#endif
	b.CleanBuildIfOf();
	myband3 = &b;
	moreuas_b3.Init();
	ua_out_seen = 0;
	memcpy(&genb12.grid0[54], b.band0, 4 * 27);
	memset(&hh0, 0, sizeof hh0);
	//hh0.diagh = diagbug;
	stack_countf = b.stack_count;
	smin = b.smin;
	uint32_t nmiss = 6 - smin.mincount;
	wactive0 = fstk =b.fstk;
#ifdef DEBUGB3
	if (locdiag) {
		cout << Char27out(wactive0) << " active ";
		cout << stack_countf.u16[0] << stack_countf.u16[1] << stack_countf.u16[2] << endl;
	}
	//Debug_If_Of_b3();
#endif
#ifdef DEBUGKNOWN
	cout <<b.band<< " gob3 smincount =" << b.smin.mincount << endl;
#endif
	if (!nmiss) {
		if (nuasb3_2)return;
		else {
			p_cpt2g[10]++;
			hh0.GoMiss0((*myband3)); 
			return; 
		}
	}
	if (nmiss == 1) {
		if (!b.and_out)return;
		else {	
			p_cpt2g[11]++;
			hh0.GoMiss1((*myband3));
			return;	
		}
	}
	if (nmiss == 2 && nuasb3_2>1 && (!myband3->and_out)) {// want always 2 clues 
		// start with the smallest ua next will be "and" of remaining uas
		p_cpt2g[12]++;
		hh0.GoMiss2Init((*myband3));
		uint32_t  uamin = uasb3_2[0];
		{
			register uint32_t min = _popcnt32(uamin);
			for (uint32_t i = 1; i < nuasb3_2; i++) {
				register uint32_t Ru = uasb3_2[i], cc = _popcnt32(Ru);
				if (cc < min) {
					min = cc;
					uamin = Ru;
				}
			}
		}
		hh0.GoMiss2((*myband3), uamin);
		return;
	}
	if (nuasb3_2) {// all uas b3 same table
		memcpy(&uasb3_1[nuasb3_1], uasb3_2, nuasb3_2 * sizeof uasb3_1[0]);
		nuasb3_1 += nuasb3_2;
	}
	ExpandB3();
}



//__________ phase 2___ find band 3 clues for one band 3

void G17B3HANDLER::GoMiss0(STD_B3 & b3) {
	smin = b3.smin;
	uasb3if = g17b.uasb3_1;
	nuasb3if = g17b.nuasb3_1;
	active_b3 = b3.smin.critbf;
	known_b3 = rknown_b3 = 0;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	CriticalLoop();
}
void G17B3HANDLER::GoMiss1(STD_B3 & b3) {
	nmiss = 1;
	smin = b3.smin;
	stack_count = b3.stack_count;
	stack_count_b12 = g17b.stack_count;
	uasb3if = g17b.uasb3_1;
	nuasb3if = g17b.nuasb3_1;
	nuasb3of = g17b.nuasb3_2;
	active_b3 = b3.smin.critbf;
	known_b3 = rknown_b3 = 0;
	wactive0 = b3.wactive0;
	//wua = g17b.andmiss1&wactive0;
	wua = b3.and_out;
	Do_miss1();
}
void G17B3HANDLER::Do_miss1() {
#ifdef DEBUGB3
	cout << Char27out(wua) << " wua entry miss1   nout=" << nuasb3of << endl;
	cout << stack_count_b12.u16[0] << stack_count_b12.u16[1] << stack_count_b12.u16[2]
		<< " count per stack bands 1 2" << endl;
	cout << stack_count.u16[0] << stack_count.u16[1] << stack_count.u16[2] << endl;
	smin.Status("entry miss1");
#endif	
	/*
	if (diagh) {
		cout << Char27out(wua) << " wua entry miss1   nout=" << nuasb3of << endl;
		cout << stack_count_b12.u16[0] << stack_count_b12.u16[1] << stack_count_b12.u16[2]
			<< " count per stack bands 1 2" << endl;
		cout << stack_count.u16[0] << stack_count.u16[1] << stack_count.u16[2] << endl;
		cout << Char27out(wactive0) << " wactive0" << endl;
		cout << Char27out(known_b3) << " known b3" << endl;
		cout << Char27out(active_b3) << " active_b3" << endl;
	}
	*/
	if (!nuasb3of) {// subcritical in hn if solved
		G17B3HANDLER hn = *this;
		hn.Critical2pairs(1);// only in critical stacks
#ifdef DEBUGB3
		hn.smin.Status("after crit 1 2pairs");
		cout << Char27out(hn.known_b3) << " known b3" << endl;
		cout << Char27out(hn.active_b3) << " active_b3" << endl;

#endif


		hn.SubCriticalLoop();// try direct in field
		if (g17b.aigstop || g17b.aigstopxy)return;
		if (g17b.ua_out_seen) {
			wua &= g17b.ua_out_seen;
			g17b.ua_out_seen = 0;
		}
		else if (g17b.moreuas_b3.nt) {// Uas to insert as uasb1
			g17b.moreuas_b3.InsertIn(uasb3if, nuasb3if);
			// here room for a revision of the critical status
			// if a new Gua2 appeared, the game is over
		}
	}
	uint32_t res;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	//cout << Char27out(wua) << " wua go out" << endl;
	while (bitscanforward(res, wua)) {
		int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
		G17B3HANDLER hn = *this;
		hn.known_b3 |= bit;
		hn.stack_count.u16[C_stack[res]]++;
		hn.CriticalLoop();
		if (g17b.aigstop || g17b.aigstopxy)return;
		if (g17b.ua_out_seen) {
			cout << Char27out(g17b.ua_out_seen) << " g17b.ua_out_seen" << endl;
			wua &= g17b.ua_out_seen;
			g17b.ua_out_seen = 0;
		}
	}
}



void G17B3HANDLER::GoMiss2Init(STD_B3 & b3) {
	smin = b3.smin;
	active_b3 = smin.critbf;
	known_b3 = rknown_b3 = 0;
	wactive0 = g17b.wactive0 & (BIT_SET_27 ^ active_b3);//  active cells out field
	stack_count = b3.stack_count;
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	for (int istack = 0, stp = 0111; istack < 3; istack++, stp <<= 1)
		if (stack_count.u16[istack] > 5) {// critical stack
			register int m2stack = stp & smin.mini_bf2, shrink = TblShrinkMask[m2stack];
			if (m2stack) {// common cell(s) to assign
				register int Mask = tbitsrows[shrink] << (3 * istack);
				//adjust count and known
				known_b3 |= Mask & (~smin.pairs27);// and set the common cell as assigned
				smin.mini_bf2 &= ~stp; // clear the 2pairs bit(s) in stack
				active_b3 &= (~Mask);// clear the  field bf
				smin.critbf &= (~Mask);
				smin.pairs27 &= (~Mask);
				smin.mincount -= _popcnt32(shrink);
			}
		}
}

void G17B3HANDLER::GoMiss2(STD_B3 & b3, uint32_t uamin) {
	nmiss = 2;
	uasb3if = g17b.uasb3_1;	nuasb3if = g17b.nuasb3_1;
	uasb3of = g17b.uasb3_2;	nuasb3of = g17b.nuasb3_2;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	wua = uamin;
	// cells added must produce cells hitting all remaining uas
	uint32_t res, tcellsok[27][3], ntcellsok = 0 ;
	while (bitscanforward(res, wua)) {
		uint32_t   nout=0;
		register uint32_t  bit = 1 << res;
		wua ^= bit; 
		wactive0 ^= bit;
		//if (g17b.debug17 && (!(bit & g17b.p17diag.bf.u32[2])))continue;
		//if (diagh) cout << Char27out(bit) << "seen" << endl;
		register uint32_t andx = wactive0, s = C_stack[res];
		if (stack_count.u16[s] == 5) {
			s = ~(07007007 << (3 * s));// mask
			andx &= s;
		}
		for (uint32_t i = 0; i < nuasb3of; i++) {
			register uint32_t ua = uasb3of[i];
			if (!(ua&bit)) {
				nout = 1;	andx &= ua;	if (!andx) break;
			}
		}
		if (andx|| (!nuasb3of)) {
			tcellsok[ntcellsok][0] = res;
			tcellsok[ntcellsok][1] = andx;
			tcellsok[ntcellsok++][2] = nout; 

		}
	}
	//if (diagh) 	cout <<  "  miss2 uamin ntcellsok=" << ntcellsok << endl;
	for (uint32_t i = 0; i < ntcellsok; i++) {// now call 
		G17B3HANDLER hn = *this;
		hn.AddCell_Miss2(tcellsok[i]);
		if (g17b.aigstop || g17b.aigstopxy)return;
	}
}
inline void G17B3HANDLER::AddCell_Miss2(uint32_t * t) {//uint32_t cell, int bit) {
	if (g17b.diagbug == 2) {
		cout << "  miss2 add cell " << t[0] << " nuaof=" << t[2] << endl;
		cout << Char27out(t[1]) << " wua" << endl;
	}
	nuasb3of = t[2];
	wua = t[1];
	{
		register uint32_t cell = t[0];
		register int s = C_stack[cell];
		stack_count.u16[s]++;
		if (stack_count.u16[s] > 5) {
			s = ~(07007007 << (3 * s));// mask
			wua &= s;
			wactive0 &= s;
		}
	}
	nmiss--;
	known_b3 |= 1 << t[0];
	uint32_t res;
	while (bitscanforward(res, wua)) {
		int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
		G17B3HANDLER hn = *this;
		hn.known_b3 |= bit;
		hn.stack_count.u16[C_stack[res]]++;
		hn.CriticalLoop();
		if (g17b.aigstop || g17b.aigstopxy)return;
		if (g17b.ua_out_seen) {
			wua &= g17b.ua_out_seen;
			g17b.ua_out_seen = 0;
		}
	}
}
void G17B3HANDLER::Critical2pairs(int modesub) {// assign 2 pairs in minirow to common cell
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	if (smin.mini_bf2) {// and 2 pairs in minirow forced to common cell
		register int Rst = 07007007;// stack 0 pattern
		for (int ist = 0; ist < 3; ist++) {
			if (modesub && stack_count.u16[ist] != 6) continue;
			int shrink = TblShrinkMask[smin.mini_bf2 & (0111 << ist)];
			if (shrink) {// minirows 2 pairs in that stack
				register int Mask = tbitsrows[shrink] << (3 * ist);
				active_b3 &= (~Mask); // clear the minirow
				known_b3 |= Mask & (~smin.pairs27);// and set the common cell as assigned
			}
		}
		if(!modesub)smin.mini_bf2 = 0;
	}
}
void G17B3HANDLER::CriticalAssignCell(int Ru) {// assign a cell within the critical cells
	// Ru is usually a regidster containing a 27 bits field with one bit on
	// 2 pairs in a miniriow have already been applied
	known_b3 |= Ru;
	uint32_t cell;
	bitscanforward(cell, Ru); // catch the cell
	register int mini = C_minirow[cell],// minirow to clear
		bit = 1 << mini,
		Mask = 7 << (3 * mini);
	if (bit & smin.mini_bf3) {// the cell is in a minirow with 3 pairs active
		active_b3 &= ~Ru; //clear the cell
		smin.mini_bf3 ^= bit; // now only a pairto hit
		smin.mini_bf1 |= bit;
	}
	else {// either one pair or a triplet in the minirow
		active_b3 &= (~Mask); // kill the minirow as active
		smin.mini_bf1 &= ~bit;
		smin.mini_triplet &= ~bit;
	}
}

void G17B3HANDLER::CriticalLoop() {// after optional assignment
	while (1) {// first shrink uas in field
		//if (g17b.diagbug>1)cout << "CriticalLoop() cycle" << endl;
		//if (1)DebugCycle();
		irloop = 0;
		uint32_t * tn = &uasb3if[nuasb3if], n = 0;
		register uint32_t Ra = active_b3,
			Rfilt = known_b3;
		for (uint32_t iua = 0; iua < nuasb3if; iua++) {
			register int Ru = uasb3if[iua];
			if (Ru & Rfilt) continue;// already hit, forget it
			Ru &= Ra;
			if (!Ru) return;// dead branch
			if (_popcnt32(Ru) == 1) {// assign it and reduce the active cells
				CriticalAssignCell(Ru);
				Rfilt = known_b3;
				Ra = active_b3; //can be  modified
				irloop = 1;// should loop for new singles
			}
			else tn[n++] = Ru;			
		}
		uasb3if = tn;
		nuasb3if = n;
		if (!n) irloop = 0;// no need to loop again
		if (!irloop) break;
	}
	if (_popcnt32(known_b3) > 6) return;
	if (!active_b3) {// must be here expected number of clues
		if (nuasb3if) return; //can not be valid
		//cout << Char27out(known_b3) << "final check crit" << endl;
		//zh_g2.cells_assigned.bf.u32[2] = 0;
		g17b.FinalCheckB3(known_b3);
		//if (zh_g2.cells_assigned.bf.u32[2])
			//cout << Char27out(zh_g2.cells_assigned.bf.u32[2]) << " new ua" << endl;
		return; // branch closed
	}
	int wua = uasb3if[0] & active_b3, cell;
	while (bitscanforward(cell, wua)) {
		register int bit = 1 << cell;
		wua ^= bit;// clear bit

		// clean the bit in active_b3, this is now a dead cell downstream
		active_b3 ^= bit;
		//if (g17b.debug17 && (!(bit & g17b.p17diag.bf.u32[2])))continue;
		G17B3HANDLER hn = *this;
		hn.CriticalAssignCell(bit);
		//if (g17b.debug17) hn.smin.Status("after assign");
		hn.CriticalLoop();
	}
}


void G17B3HANDLER::SubCriticalLoop() {// after optional assignment
	// loop1 reduce uas to in field


	while (1) {// first shrink uas in field
		//if (g17b.diagbug>1)cout << "CriticalLoop() cycle" << endl;
		//if (1)DebugCycle();
		irloop = 0;
		uint32_t * tn = &uasb3if[nuasb3if], n = 0;
		register uint32_t Ra = active_b3,
			Rfilt = known_b3;
		for (uint32_t iua = 0; iua < nuasb3if; iua++) {
			register int Ru = uasb3if[iua];
			if (Ru & known_b3) continue;// already hit, forget it
			Ru &= Ra;
			if (!Ru) return;// dead branch
			if (_popcnt32(Ru) == 1) {// assign it and reduce the active cells
				known_b3 |= Ru;
				active_b3 ^= Ru;// dead for next steps
				uint32_t cell;
				bitscanforward(cell, Ru); // catch the cell
				register int s = C_stack[cell];
				stack_count_b12.u16[s]++;
				if (stack_count_b12.u16[s] > 5) {
					s = ~(07007007 << (3 * s));// mask
					active_b3 &= s;
					Rfilt = known_b3;
					Ra = active_b3; //can be  modified
				}
				irloop = 1;// should loop for new singles
			}
			else tn[n++] = Ru;
		}
		uasb3if = tn;
		nuasb3if = n;
		if (!n) irloop = 0;// no need to loop again
		if (!irloop) break;
	}
	int cc = _popcnt32(known_b3);
	if (cc > 6) return;
	if (cc == 6) {
		if (nuasb3if) return; //can not be valid
		p_cpt2g[56]++;
#ifdef DEBUGB3
		cout << "final check sub crit" << endl;
		zh_g2.cells_assigned.bf.u32[2] = 0;
#endif
		g17b.FinalCheckB3(known_b3);
#ifdef DEBUGB3
		if (zh_g2.cells_assigned.bf.u32[2])
			cout << Char27out(zh_g2.cells_assigned.bf.u32[2]) << " new ua" << endl;
#endif
		return; // branch closed
	}
	uint32_t wua, cell;
	
	if (cc == 5) {// must hit all residual uas
		wua =  active_b3;
		for (uint32_t i = 0; i < nuasb3if; i++)
			wua &= uasb3if[i];
		if (!wua) return;
	}
	else wua = uasb3if[0];
	while (bitscanforward(cell, wua)) {
		register int bit = 1 << cell;
		wua ^= bit;// clear bit
		active_b3 ^= bit;// dead for next steps
		//if (g17b.debug17 && (!(bit & g17b.p17diag.bf.u32[2])))continue;
		G17B3HANDLER hn = *this;
		//hn.SubCriticalAssignCell(bit);
		hn.known_b3 |= bit;
		{
			register int s = C_stack[cell];
			hn.stack_count_b12.u16[s]++;
			if (hn.stack_count_b12.u16[s] > 5) {
				s = ~(07007007 << (3 * s));// mask
				hn.active_b3 &= s;
			}
		}
		//if (g17b.debug17) hn.smin.Status("after assign");
		hn.SubCriticalLoop();
		if (g17b.moreuas_b3.nt) {
			//cout << "should check the new ua constraint" << endl;
			if(_popcnt32(known_b3)>4) 
			 g17b.moreuas_b3.CheckNew(known_b3, wua);
		}
	}
}

void G17B::ExpandB3(){// uint32_t *tua, uint32_t nua) {// find all 5 and 6 clues solutions
	// Build tua
	p_cpt2g[13]++;
	uint32_t *tuaw=uasb3_1, nuaw = nuasb3_1;// use  in field pre loaded in the right way
#ifdef DEBUGKNOWN
	if (p_cpt2g[13]<2 && DEBUGKNOWN>1) {
		cout << "ExpandB3() nuaw=" << nuaw
			<< "\tstacks b12 " << stack_count.u16[0] << stack_count.u16[1] << stack_count.u16[2]
			<< "\tstacks " << stack_countf.u16[0] << stack_countf.u16[1] << stack_countf.u16[2] << endl;
	}
#endif

	struct SPB3 {// spots to find band 3 minimum valid solutions
		// ====================== constant after initialization
		uint32_t  possible_cells, all_previous_cells, active_cells, iuab3,
			stack[3];
	}spb3[7], *s3, *sn3;
	s3 = spb3;
	s3->all_previous_cells = 0;
	s3->active_cells = BIT_SET_27;// all cells active

	// init the stack status
	for (int i = 0; i < 3; i++) {
		s3->stack[i] = stack_count.u16[i];// count before band 3 min count
		if (s3->stack[i] == 6)s3->active_cells &= ~(07007007 << (3 * i));
	}
	s3->iuab3 = 0; // copy the start table
	s3->possible_cells = tuaw[0] & s3->active_cells;
	int tcells[10];

	//____________ here start the search 6 clues
next:
	uint64_t ispot = s3 - spb3;
	{// catch and apply cell in bitfields
		register uint32_t cell, p = s3->possible_cells;
		if (!p)goto back;
		bitscanforward(cell, p);
		register int bit = 1 << cell;
		tcells[ispot] = cell;
		s3->possible_cells ^= bit;// clear bit
		register int filter = s3->all_previous_cells | bit,
			ac = s3->active_cells ^ bit;
		sn3 = s3 + 1; *sn3 = *s3; // (copy the stack count)
		sn3->all_previous_cells = filter;
		sn3->active_cells = s3->active_cells = ac;
		{// if the stack is limit update sn3 active
			int st = C_stack[cell];
			sn3->stack[st]++;
			if (sn3->stack[st] > 5)
				sn3->active_cells &= ~(07007007 << (3 * st));
		}
		// nextspot:take the next available ua to loop
		for (uint32_t i = s3->iuab3 + 1; i < nuaw; i++) {
			if (tuaw[i] & filter)continue;
			if (ispot >= 5) 	goto next;//passing the limit
			sn3->iuab3 = i;
			register uint32_t Ru = tuaw[i] & sn3->active_cells;
			if (!Ru)goto next;
			if (ispot == 4) {// last must hit all remaining uas
				for (uint32_t i2 = i + 1; i2 < nuaw; i2++) {
					if (tuaw[i2] & filter)continue;
					Ru &= tuaw[i2];
					if (!Ru)goto next;
				}
			}
			sn3->possible_cells = Ru;
			s3 = sn3; // switch to next spot
			goto next;
		}
	}
	if (ispot < 5) {// no more uas use active as possible
		sn3->possible_cells = sn3->active_cells;
		sn3->iuab3 = nuaw;
		s3 = sn3; // switch to next spot
		goto next;
	}
	p_cpt2g[30]++;	// this is a possible 17 do final check
	if (!clean_valid_done) {
		clean_valid_done = 1;
		myua = zh2b[0].ValidXY(tclues, nclues + nclues_step);
		if (myua) {
			NewUaB12();		aigstopxy = 1; return;
		}
	}

	if (zhou[1].CallMultipleB3(zhou[0], sn3->all_previous_cells, 0)) {
		register uint32_t ua = zh_g2.cells_assigned.bf.u32[2];
		if (nuaw < 300)tuaw[nuaw++] = ua;
		NewUaB3();		
		s3->possible_cells &= ua;
	}
	else Out17(sn3->all_previous_cells);
	goto next;

	// going back, for a non empty index, count it back
back:
	if (--s3 >= spb3)goto next;
}



//=============== part 2  band 3 processing using guas2/3

int ZHOU::CallMultipleB3(ZHOU & o, uint32_t bf, int diag) {
	*this = o;
	//if (diag) cout << Char27out(bf) << " call multipleb3" << endl;
	BF128 dca[9];
	int digitsbf = zh_g2.digitsbf;
	memcpy(dca, zh_g2.Digit_cell_Assigned, sizeof dca);
	{	
		uint32_t cc;
		register int x = bf;
		while (bitscanforward(cc, x)) {
			x ^= 1 << cc; //clear bit
			int cell = cc + 54, digit = genb12.grid0[cell];
			digitsbf |= 1 << digit;
			int xcell = cc + 64; // the cell value in 3x32 of a 128 bits map
			if (FD[digit][0].Off(xcell))  return 0;// check not valid entry
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}
	if (_popcnt32(digitsbf < 8)) 	return 1;// not  one solution
	zh_g2.s17_b3_mini = 1;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	//if (diag) ImageCandidats();
	//__________end assign last lot start solver
	zh_g.go_back = 0;	zh_g.nsol = 0; // modevalid is set to  1
	int ir = Full17Update();
	//if (diag) {	cout << "after update" << endl;	ImageCandidats();	}
	if (ir == 2) return 0;// solved can not be multiple
	Guess17(0,diag);
	return zh_g.nsol;  
}
int ZHOU::Apply17SingleOrEmptyCells() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched 
	BF128 Map;
	BF128 R1 = FD[0][0], R2 = R1 & FD[1][0]; 	R1 |= FD[1][0];
	BF128 R3 = R2 & FD[2][0]; R2 |= R1 & FD[2][0]; R1 |= FD[2][0];
	Map = FD[3][0];	BF128 R4 = R3 & Map; R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = FD[4][0]; BF128 R5 = R4 & Map; R4 |= R3 & Map; R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

	Map = FD[5][0];  BF128 R6 = R5 & Map;
	R5 |= R4 & Map; R4 |= R3 & Map; R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

	Map = FD[6][0];  BF128 R7 = R6 & Map; R6 |= R5 & Map;
	R5 |= R4 & Map; R4 |= R3 & Map; R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

	Map = FD[7][0];  R7 |= R6 & Map; R6 |= R5 & Map;
	R5 |= R4 & Map; R4 |= R3 & Map; R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

	Map = FD[8][0]; R7 |= R6 & Map; R6 |= R5 & Map;
	R5 |= R4 & Map; R4 |= R3 & Map; R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	if ((cells_unsolved - R1).isNotEmpty()) 	return 1; // empty cells
	R1 -= R2; // now true singles
	R1 &= cells_unsolved; // these are new singles
	if (R1.isEmpty()) {// no single store pairs and more
		if (cells_unsolved.bf.u32[2]) {// use only b3
			if (R7.bf.u32[2])zh_g2.cells_for_guess = R7;
			else if (R6.bf.u32[2])zh_g2.cells_for_guess = R6;
			else if (R5.bf.u32[2])zh_g2.cells_for_guess = R5;
			else if (R4.bf.u32[2])zh_g2.cells_for_guess = R4;
			else if (R3.bf.u32[2])zh_g2.cells_for_guess = R3;
			else zh_g2.cells_for_guess = R2;
		}
		else {
			if (R7.bf.u64[0])zh_g2.cells_for_guess = R7;
			else if (R6.bf.u64[0])zh_g2.cells_for_guess = R6;
			else if (R5.bf.u64[0])zh_g2.cells_for_guess = R5;
			else if (R4.bf.u64[0])zh_g2.cells_for_guess = R4;
			else if (R3.bf.u64[0])zh_g2.cells_for_guess = R3;
			else zh_g2.cells_for_guess = R2;
		}
		return 0;
	}
	int tcells[80], ntcells = R1.Table3X27(tcells);
	for (int i = 0; i < ntcells; i++) {
		int cell = tcells[i];
		for (int idig = 0; idig < 9; idig++) {
			if (FD[idig][0].On_c(cell)) {
				Assign(idig, cell, C_To128[cell]);
				goto nextr1;
			}
		}
		return 1; // conflict with previous assign within this lot
	nextr1:;
	}
	zh_g.single_applied = 1;
	return 0;
}
int ZHOU::Full17Update() {
	if (zh_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (Apply17SingleOrEmptyCells())	return 0; // locked empty cell or conflict singles in cells
		if (!zh_g.single_applied)	break;
	}
	return 1;
}
void ZHOU::Guess17(int index, int diag) {
	if (zh_g.go_back) return;
	if (index > 20) return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< debugging temporary code
	if (diag) ImageCandidats();
	uint32_t xcell, cell, digit;
	if (cells_unsolved.bf.u32[2]) {// fill in priority band 3
		uint32_t w3 = zh_g2.cells_for_guess.bf.u32[2];
		bitscanforward(cell, w3);
		cell += 54;
		xcell = cell + 10;
		if (diag) {
			cout << Char27out(w3) << " w3 go for " << zh_g2.grid0[cell] + 1 << cellsFixedData[cell].pt << endl;
		}
		//if (!zh_g2.isfalse_on) {// skip if band 3 full
			//if (_popcnt32(cells_unsolved.bf.u32[2]) < 2) {
				//digit = zh_g2.grid0[cell];
				//FD[digit][0].Clear(xcell);// force false
			//}
		//}
	}
	// stop if all true 
	else {// fill band 12 highest first 
		if (!zh_g2.isfalse_on)	return; // already checked
		uint64_t w12 = zh_g2.cells_for_guess.bf.u64[0];
		bitscanforward64(xcell, w12);
		cell = From_128_To_81[xcell];
		if (diag) cout << zh_g2.grid0[cell] + 1 << cellsFixedData[cell].pt << " in band 12 index=" << index << endl;
	}
	digit = zh_g2.grid0[cell];
	// true first if possible
	if (FD[digit][0].On(xcell)) {
		if (diag) cout << digit + 1 << cellsFixedData[cell].pt << " true index=" << index << endl;
		ZHOU * mynext = (this + 1);
		*mynext = *this;
		mynext->SetaCom(digit, cell, xcell);
		mynext->Compute17Next(index + 1, diag);
		if (zh_g.go_back) return;

	}
	// if first step try first false
	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue;
		if (FD[idig][0].On(xcell)) {
			if (diag) cout << idig + 1 << cellsFixedData[cell].pt << " false " << endl;
			if (xcell >= 64)zh_g2.isfalse_on = 1;
			ZHOU * mynext = (this + 1);
			*mynext = *this;
			mynext->SetaCom(idig, cell, xcell);
			mynext->Compute17Next(index + 1, diag);
			if (zh_g.go_back) return;
		}
	}
}
void ZHOU::Compute17Next(int index, int diag) {
	int ir = Full17Update();
	if (!ir) return;// locked 
	if (ir == 2) {//solved
		if (index) {// store false as ua
			BF128 & wua = zh_g2.cells_assigned;
			int * sol = genb12.grid0;
			wua.SetAll_0();;
			for (int i = 0; i < 81; i++) {
				int d = sol[i];
				if (FD[d][0].Off_c(i))	wua.Set_c(i);
			}
			if (wua.isNotEmpty()) {// ignore true solution
				zh_g.nsol++;
				zh_g.go_back = 1;// closed anyway
			}
		}
		return;
	}
	Guess17(index, diag);// continue the process
}
/*
void G17B3HANDLER::Init( ) {
	G17B & bab = g17b;
	smin = bab.smin;
	uasb3of = bab.uasb3_2;	nuasb3of = bab.nuasb3_2;
	uasb3if = bab.uasb3_1;	nuasb3if = bab.nuasb3_1;
	andoutf = bab.b3_andout;
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	known_b3 = rknown_b3 = 0;
	ndead = BIT_SET_27;
	active_b3 =smin.critbf;// active cells in field
	wactive0 = BIT_SET_27 ^ active_b3;//  active cells out field
	nmiss = 6 - smin.mincount;
	nb3 = 6;
	stack_count = bab.stack_countf;
	for (int istack = 0, stp = 0111; istack < 3; istack++, stp <<= 1)
		if (stack_count.u16[istack] > 5) {// critical stack
			wactive0 &= ~(07007007 << (3 * istack));// clear outfield
			register int m2stack = stp & smin.mini_bf2, shrink = TblShrinkMask[m2stack];
			if (m2stack) {// common cell(s) to assign
				register int Mask = tbitsrows[shrink] << (3 * istack);
				//adjust count and known
				known_b3 |= Mask & (~smin.pairs27);// and set the common cell as assigned
				smin.mini_bf2 &= ~stp; // clear the 2pairs bit(s) in stack
				active_b3 &= (~Mask);// clear the  field bf
				smin.critbf &= (~Mask);
				smin.pairs27 &= (~Mask);
				smin.mincount -= _popcnt32(shrink);
			}
		}

}
*/

uint32_t G17B3HANDLER::IsMultiple(int bf) {
	if (bf == rknown_b3) return 0;
	uint32_t ua = 0;
	rknown_b3 = bf;
	G17B & bab = g17b;
	// check first if all tuab3 is hit
	p_cpt2g[55] ++;
	int ir = zhou[1].CallMultipleB3(zhou[0], bf, 0);
	//if (g17b.diag >= 2)	cout << Char27out(ir) << "ir retour multiple ir=" << ir << endl;
	if (ir) {
		ua = zh_g2.cells_assigned.bf.u32[2];
		g17b.moreuas_b3.Add(ua);
		g17b.NewUaB3();
		//if (g17b.diag >= 2)cout << Char27out(ua) << " b3 ua to add from is multiple" << endl;
	}
	return ua;
}

//================= critical process


void G17B3HANDLER::Go_Critical(){// critical situation all clues in pairs tripl:ets
	//if (g17b.debug17 > 1 && known_b3)cout << Char27out(known_b3) << " entry critical" << endl;
	//if (g17b.diag >= 2 || diagh) 	cout << Char27out(known_b3)
	//	<< "entry critical nb3if= " <<	nuasb3if << endl;

	active_b3 = smin.critbf;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	CriticalLoop();
}


//=============== sub critical process   missing(s)  in the critical area
void G17B3HANDLER::Go_SubcriticalMiniRow() {
	//if (diagh)cout << "entry Go_SubcriticalMiniRow() ndead=" << ndead << endl;
	int c2[3] = { 3, 5, 6 };// 2 cells in a mini row
	int bit = 1 << ndead, mask = 7 << (3 * ndead);
	for (int i = ndead; i < 9; i++,  bit <<= 1, mask <<= 3) {
		stack = i % 3;
		register int M = active_sub & mask;
		if (!M)continue;
		ndead = i;
		//if (diagh)cout << Char27out(M) << " mini row to process i=" << i << endl;
		if (bit & smin.mini_bf1) {// it was a gua2 pair assign both
			G17B3HANDLER hn = *this;
			hn.smin.mini_bf1 ^= bit;
			hn.SubMini( M, mask);
		}
		else if (bit & smin.mini_bf2)// it was 2 gua2 pair assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				G17B3HANDLER hn = *this;
				hn.smin.mini_bf2 ^= bit;
				hn.SubMini(M, mask);
			}
		else if (bit & smin.mini_bf3) {// it was 3 gua2 pair assign 3 out of 3
			G17B3HANDLER hn = *this;
			hn.smin.mini_bf3 ^= bit;
			hn.SubMini(M, mask);
		}
		else if (bit & smin.mini_triplet)// it was a gua3 triplet assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				G17B3HANDLER hn = *this;
				hn.smin.mini_triplet ^= bit;
				hn.SubMini(M, mask);
			}
		else {// second add in the mini row one residual cell take it
			G17B3HANDLER hn = *this;
			hn.SubMini(M, mask);
		}
	}
}
void G17B3HANDLER::SubMini( int M, int mask) {
	known_b3 |= M;// assign 1 or 2
	nmiss--;// one added
	active_b3 &= ~mask;
	active_sub ^= M;
	// now adjust the stack count
	stack_count.u16[stack]++;
	if (stack_count.u16[stack] > 5)active_sub &= ~(07007007 << (3 * stack));
	/*
	if (diagh) {
		cout << "SubMini go nmiss="<<nmiss << endl;
		cout << Char27out(known_b3) << "known b3" << endl;
		cout << Char27out(active_b3) << "active_b3" << endl;
	}
	*/
	if (nmiss) Go_SubcriticalMiniRow();// continue till a"no missing clue condition"
	else {	// leave sub critical mode and enter the critical mode
		Critical2pairs();// assign 2 pairs in minirow to common cell
		/*
		if (diagh) {
			cout << "SubMini go after assign 2 pairs" << endl;
			cout << Char27out(known_b3) << "known b3" << endl;
			cout << Char27out(active_b3) << "active_b3" << endl;
			smin.Status("SubMini go after assign 2 pairs");
		}
		*/
		CriticalLoop();
	}
}
void G17B3HANDLER::Go_Subcritical() {// nmiss to select in the critical field
	active_b3 = active_sub = smin.critbf;
	// check first if a global solution  is still possible
	for (int ist = 0; ist < 3; ist++) {// check stacks
		if (stack_count.u16[ist] > 5)active_sub &= ~(07007007 << (3 * ist));// kill the stack for more clues
	}
	/*
	if (diagh) {
		cout << Char27out(known_b3) << " entry Go_Subcritical() nmiss= "<<nmiss << endl;
		cout << Char27out(smin.critbf) << " smin.critbf " << endl;
		cout << Char27out(active_sub) << " active_sub " << endl;
	}
	*/
	ndead = 0;
	Go_SubcriticalMiniRow();// find the first miss
}


//________ final called by all branches
void G17B::FinalCheckB3(uint32_t bfb3) {
	p_cpt2g[29]++;
	if (aigstopxy) return;
	if (_popcnt32(bfb3) > 6) 		return;
	if (moreuas_b3.Check(bfb3))return;
	if (!clean_valid_done) {
		clean_valid_done = 1;
		myua = zh2b[0].ValidXY(tclues, nclues + nclues_step);
		if (myua) {
			NewUaB12();		aigstopxy = 1; return;
		}
	}	
	//cout << " call zhou p_cpt2g[29]=" << p_cpt2g[29] << endl;
	register uint32_t ir = zhou[1].CallMultipleB3(zhou[0], bfb3, 0);
	if (ir) {
		register uint32_t ua = zh_g2.cells_assigned.bf.u32[2];
		if (ua && (!(ua & hh0.smin.critbf)))	ua_out_seen = ua;
		NewUaB3();
		moreuas_b3.Add(ua);// if empty lock the call 
		return;
	}
	Out17(bfb3);
}
void G17B::Out17(uint32_t bfb3) {
	cout << Char27out(bfb3) << "\t\tone sol to print final check "<< p_cpt2g[48] << endl;
	char ws[82];
	strcpy(ws, empty_puzzle);
	for (int i = 0; i < (nclues_step+nclues); i++) {
		int cell = tclues[i];
		ws[cell] = genb12.grid0[cell] + '1';
	}
	for (int i = 0, bit = 1; i < 27; i++, bit <<= 1)if (bfb3 & bit)
		ws[54 + i] = genb12.grid0[54 + i] + '1';
	fout1 << ws << ";" << genb12.nb12 / 64 << ";" << genb12.i1t16 << ";" << genb12.i2t16 << endl;
	a_17_found_here++;
}

void G17B::NewUaB12() {
	uint64_t cc64 = _popcnt64(myua&BIT_SET_2X);
	if (cc64 < 12) {// this should never be check for a bug
		DebugAdd12();		return;
	}
	//if (cc64 < 18)cout << Char2Xout(myua) << " addua cc=" << cc64 << endl;
	moreuasxy.Add(myua);
	if (cc64 < 18) {
		if (cc64 < 14)moreuas_12_13.Add(myua);
		else if (cc64 == 14)moreuas_14.Add(myua);
		else if (cc64 == 15)moreuas_15.Add(myua);
		else moreuas_AB_small.Add(myua);
		if (genuasb12.nua < 1500 || cc64 < 16) {
			register uint64_t ua_add = myua | (cc64 << 59);
			genuasb12.AddUACheck(ua_add);
			if (ntusb1 < 1000)tusb1[ntusb1++] = myua;
			p_cpt2g[31]++;
		}
	}
	else if (cc64 < 21)			moreuas_AB.Add(myua);
	else moreuas_AB_big.Add(myua);
}


void G17B::DebugAdd12() {
	aigstop = 1;
	cerr << "ua < 12 to add clean" << endl;
	cout << endl << endl << Char2Xout(myua) << " ua < 12 to add   clean" << endl;
	cout << "bug location band 2 id=" << genb12.nb12 << endl;
	cout << myband1.band << endl;
	cout << myband2.band << endl;
	cout << Char2Xout(wb12bf) << " b12 at call" << endl;
	cout << "ntusb1=" << ntusb1 << " n11=" << n_to_clean << endl;
	for (int i = 0; i < nclues_step; i++) cout << tclues[i] << " ";
	cout << "\t";
	for (int i = 0; i < nclues; i++) cout << tcluesxy[i] << " ";
	cout << endl;
	zh2b_i.ImageCandidats();
	zh2b_i1.ImageCandidats();
	zh2b[0].ImageCandidats();
	cout << "table uas" << endl;
	uint64_t *t = genuasb12.tua;
	uint32_t n = genuasb12.nua;
	for (uint32_t i = 0; i < n; i++) {
		uint64_t cc = _popcnt64(t[i] & BIT_SET_2X);
		if (cc > 12)break;
		cout << Char2Xout(t[i]) << " " << i << " " << cc << endl;

	}
	for (uint32_t i = 0; i < ntusb1; i++) {
		uint64_t cc = _popcnt64(tusb1[i] & BIT_SET_2X);
		if (cc > 12)break;
		cout << Char2Xout(tusb1[i]) << " b1 i=" << i << " " << cc << endl;

	}

}

void G17B::NewUaB3() {// new ua from final check zh_g2.cells_assigned
	BF128 ua128 = zh_g2.cells_assigned;
	register uint64_t ua12 = ua128.bf.u64[0];
	register uint32_t ua = ua128.bf.u32[2],
		cc = _popcnt32(ua),
		cc0 = (uint32_t)_popcnt64(ua12);
	if (!cc) {// bands 1+2 not valid should never be here
		myua = ua12;		NewUaB12();		aigstopxy = 1;
		return;
	}

	if (cc0 > 18) return;//18 max  TGUAS::ApplyLoopB2()
	if (cc > 3) {
		if ((cc0 + cc) > 15) return;
		if (cc == 4 && cc0 > 15) return;
	}

	p_cpt2g[32]++;// 2;3 or 4 cells in band 3
	uint64_t ua54 = (ua12 & BIT_SET_27) | ((ua12 & BIT_SET_B2) >> 5);
	// find the digits pattern from the current band 3
	int * cur_b3 = &genb12.grid0[54], wdigs = 0,c27;
	{
		register uint32_t wua = ua;
		while (bitscanforward(c27, wua)) {// look for  possible cells
			wua ^= 1 << c27;// clear bit
			wdigs |=1<<cur_b3[c27];
		}	
	}
	uint32_t my_i81;
	if (cc == 2)my_i81 = genb12.GET_I81_G2(wdigs, ua);
	if (cc == 4) {// can be gua2
		register  uint32_t A = ua,
			B = (A | (A >> 9) | (A >> 18)) & 0777; // all columns
		uint32_t ncols = _popcnt32(B), ndigs = _popcnt32(wdigs);
		if (ncols == 3 && ndigs == 2) {
			my_i81 = genb12.GET_I81_G2_4(wdigs, ua);
			cc = 2;
		}
	}
	if (cc >3) { 
		p_cpt2g[33]++;
		// could also be added later to other bands 3 where it is valid
		if (myband3->ntua128 < 1000) {
			myband3->tua128[myband3->ntua128++] = ua128;
			myband3->tua128_b1[myband3->ntua128_b1++] = ua128;
			myband3->tua128_b2[myband3->ntua128_b2++] = ua128;
		}
		return;
	}

	//if(cc0<16)cout << Char27out(ua) << " " << Char2Xout(ua12) << " new uab3 to store " << endl;

	if (cc == 2) {// one of the 27 GUA2s add to the table
		p_cpt2g[40]++;
		///NewUaB3_g2(my_i81, ua12);
		tguas.tgua_start[my_i81].Adduacheck(ua12);// for new steps
		tguas.tgua_b1[my_i81].Adduacheck(ua12);// for new steps
		tguas.g2ok.Set(my_i81);
		if (tguas.nvg2 <1024) {
			uint32_t ibloc = tguas.nvg2 >> 6, ir = tguas.nvg2 - 64 * ibloc;
			tvg64g2[ibloc].SetVect54(ua54, ir, my_i81);
			tguas.nvg2++;
		}
		return;
	}
	if (cc == 3) {// one of the 3 GUA3s add to the table
		p_cpt2g[41]++;
		my_i81 = genb12.GET_I81_G3(wdigs, ua);
		tguas.g3ok.Set(my_i81);
		my_i81 += 81;// now mode 162
		tguas.tgua_start[my_i81].Adduacheck(ua12);// for new steps
		tguas.tgua_b1[my_i81].Adduacheck(ua12);// for new steps
		if (tguas.nvg3 <512) {
			uint32_t ibloc = tguas.nvg3 >> 6, ir = tguas.nvg3 - 64 * ibloc;
			tvg64g3[ibloc].SetVect54(ua54, ir, my_i81);
			tguas.nvg3++;
		}
		return;
	}
}


//=============================== debugging sequences
int VALIDB::DebugFindKnown17() {
	register uint64_t R = g17b.p17diag.bf.u64[0];
	if (~R & bf) return 0;// must be both in the 17
	// check also the rest
	uint64_t target = bf;// | active | g17b.fb12;
	if(bf& BIT_SET_27)target |= BIT_SET_B2;
	else target |= BIT_SET_27;
	if (R & ~target)return 0;
	//cout << Char2Xout(bf) << " bf" << endl;
	//cout << Char2Xout(active) << " active" << endl;
	g17b.validb_known = *this;
	return 1;
}
uint64_t G17B::DebugFindKnown17_valids_2bands() {
	register uint64_t R = g17b.p17diag.bf.u64[0];
	if (n_to_clean)
		for (uint64_t i = 0; i < n_to_clean; i++)
			if(! (R & ~to_clean[i]))return to_clean[i];

	return 0;
}
void G17B::Debug_b1b2cpt() {
	cout << "number of entries downstream" << endl;
	uint64_t w = b1cpt[6] * b2cpt[5] + b1cpt[5] * b2cpt[6];
	if (w) cout << "11 clues\t" << w << endl;
	w = b1cpt[6] * b2cpt[4] + b1cpt[4] * b2cpt[6]+	b1cpt[5] * b2cpt[5];
	if (w) cout << "10 clues\t" << w << endl;
	w = b1cpt[5] * b2cpt[4] + b1cpt[4] * b2cpt[5] + 
		b1cpt[6] * b2cpt[3] + b1cpt[3] * b2cpt[6];
	if (w) cout << "9 clues\t" << w << endl;
	w = b1cpt[4] * b2cpt[4] + b1cpt[5] * b2cpt[3] + b1cpt[3] * b2cpt[5]
		+ 	b1cpt[6] * b2cpt[2] + b1cpt[2] * b2cpt[6];
	if (w) cout << "8 clues\t" << w << endl;

}
void G17B::GodebugInit(int mode) {// called from DebugK17M10
	cout << "n bands3      \t" << genb12.nband3  << "\tua bands1+2   \t" << genuasb12.nua;

	for (int i = 0; i < genb12.nband3; i++) {
		cout  <<genb12.bands3[i].band<<" " << i << "  band3 nstacks=" << genb12.bands3[i].ntua128 << endl;
	}


	if (mode & 1) {// status of guas2 3 first band
		tguas.DebugStart(0);
	}


/*
	if (debug17 > 1) {
		STD_B3 &b3 = genb12.bands3[0];
		cout << "band3 permanent data" << endl;
		b3.guas.isguasocket2.Print3("guasocket2");
		b3.guas.isguasocket3.Print3("guasocket3");
		b3.guas.isguasocket2_46.Print3("guasocket2_46");
	}



	if (mode & 2) {
		cout << "table uas" << endl;
		uint64_t *t = genuasb12.tua;
		uint32_t n = genuasb12.nua;
		for (uint32_t i = 0; i < n; i++) cout << Char2Xout(t[i]) << " " <<i<< endl;

	}*/


	/*
	if (mode & 4) {// status of guas2 3 first band
		cout << "debug status of guas2 3 first band" << endl;
		STD_B3::GUAs & bguas = genb12.bands3[0].guas;
		BF128 ws = bguas.isguasocket2 | bguas.isguasocket3;
		for (uint32_t i = 0; i < 128; i++) {
			GUAN wg = tuguan.tguan[i];
			if (wg.ncol > 3) continue;
			int i81 = wg.i81;
			if (wg.ncol == 2 && ws.On_c(i81)) {// here 3X 27
				if (bguas.isguasocket2.On_c(i81))
					cout << Char27out(bguas.ua_pair[i81]) << " sock2";
				else cout << Char27out(bguas.ua_triplet[i81]) << " sock3";

				cout <<"i = " << i	<< " i81=" << i81 << " "
					<< Char2Xout(wg.killer) << "killer" << endl;
				uint64_t *tua = wg.tua;
				uint32_t nua = wg.nua;
				for (uint32_t iu = 0; iu < nua; iu++)
					cout << Char2Xout(tua[iu]) << endl;
			}

		}

	}
	*/
	

}
int G17B::GodebugCheckUas(const char * lib) {
	uint32_t nua = genuasb12.nua;
	uint64_t * tua= genuasb12.tua;
	for (uint32_t i = 0; i < nua; i++) {
		if (tua[i] & g17b.p17diag.bf.u64[0]) continue;
		cout << lib << "check ua failed" << endl;
		cout << Char2Xout(tua[i]) << " not hit by the known 17" << endl;
		return 1;
	}
	return 0;
}
void G17B::Debug_If_Of_b3() {
	cout << "nif=" << nuasb3_1 << endl;
	for (uint32_t i = 0; i < nuasb3_1; i++)
		cout << Char27out(uasb3_1[i] )<< endl;
	cout << "nof=" << nuasb3_2 << endl;
	for (uint32_t i = 0; i < nuasb3_2; i++)
		cout << Char27out(uasb3_2[i]) << endl;
}

void G17B3HANDLER::DebugCycle() {
	cout << "CriticalLoop() cycle nif=  nuasb3if=" << nuasb3if
		<< " n_new_uas="<< g17b.moreuas_b3.nt << endl;
	cout << Char27out(known_b3) << " known" << endl;
	cout << Char27out(active_b3) << " active_b3" << endl;
	//smin.Status("");
	for (uint32_t i = 0; i < nuasb3if; i++)
		if(!(uasb3if[i]& known_b3))
		cout << Char27out(uasb3if[i]& active_b3) << endl;

}
