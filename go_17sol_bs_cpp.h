

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
	aigstop = 0;
}
void G17B::GoM10() {// processing an entry 656 566 with the relevant table of ba,ds3
	if (aigstop)return;
	p_cpt2g[0] ++;
	p_cpt2g[1] += genb12.nband3;
	ExpandB2();// expand band2 272 ms
	if (!myband2.nvalidb)return ; // mode 656 only no 5
	b1cptdiag = b1cpt[1];
	GoM10Uas();//expand bands 3  collect UAs 350 ms
	StackUas();
	bin_b1.Copy(myband1);
	bin_b2.Copy(myband2);
	Go2_Ext_Loop();//next is external (outer) loop 
}
 
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
	b.nbi2 = nbi2_1;
	b.my_validb = vab_1;
	b.nvalidb = (uint32_t)(p_cpt[0] + p_cpt[1]);
}
struct SPOT_E64 {// spots to find band12 valid solutions n clues
	uint64_t  all_previous_cells, active_cells;
	uint32_t * start_possibles, n_possibles, ipos, ispot;
	uint64_t * tua;
	uint32_t  missing_clues, nua;
	inline void Copy(SPOT_E64 * old) {
		*this = *old;
		start_possibles += n_possibles;
		ispot++;
		missing_clues--;
		ipos = 0;
		tua += nua;
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
	morev2a.Init();	morev2b.Init();	morev2c.Init();
	zh1b_g.modegua = 0;//must be to activate filter in UAs b12 more
	zh2b_g.diag = 0;
	if (genuasb12.Initgen()) return;
	genb12.BuildGang9x3();
	// _____ GUAs 
	zh1b_g.modegua = 1;//must be to kill  filter in GUAs 6_7 more
	genb12.SecondSockets2Setup();// collect GUA2s 
	genb12.SecondSockets3Setup();// collect GUA3s 
	genb12.SetUpkg();// collect killed per cell
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3 & b = genb12.bands3[ib3];
		b.guas.isguasocketc2.Convert3X27to81(b.guas.isguasocket2);
		b.guas.isguasocketc3.Convert3X27to81(b.guas.isguasocket3);
		b.guas.isguasocketc2_46.Convert3X27to81(b.guas.isguasocket2_46);
		b.isguasocketc246 = b.guas.isguasocketc2 | b.guas.isguasocketc2_46;
	}

}
void G17B::StackUas() {// band 3 uas used as gangsters via  C_transpose_d[81]
	STD_B3 wbs;
	// transpose bands 1+2  transpose bands 3
	for (int ib = 0; ib < genb12.nband3; ib++) {// bands 3
		STD_B3 &myb = genb12.bands3[ib];
		wbs.ntua128 = 0;
		memcpy(&genb12.grid0[54], myb.band0, 4 * 27);
		int zt[81];
		for (int i = 0; i < 81; i++) 
			zt[i] = genb12.grid0[C_transpose_d[i]];
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
 if 'yes' is the part hitting a ua in a band
 yes1 hit in band 1 ; yes2 hit in band2
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
	if (nextl1 | nextl2) {	return lim - 1;	}
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

void G17B::Go2_Ext_Loop() {	//_____________ outer loop
	loopb1 = 0;
	uint32_t activerb1, activerb2;
	uint64_t activeloop = BIT_SET_2X;
	activerb1= activerb2 = BIT_SET_27;
	//_________________ external loop
	while (++loopb1 << EXLNLOOP1) {
		if (aigstop)return;;
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
		if (minratio > EXLRATIO) break;
		else {
			extlr = extlw;
			if (extlw.mode == 1) {// this is a band1 X band2 Y
				ExtSplitX(bin_b1, bin_b1yes, extlw.bfx, activerb1);
				if (loopb1 == 1)
					Go2b_Ext_Loop(BIT_SET_2X | activerb1, 1);
				else  Go3(bin_b1yes, bin_b2);
				ExtSplitY(bin_b2, extlr.tbfy, extlr.ntbfy, activerb2, 2);
			}
			else {// this is a band2 X band1 Y
				ExtSplitX(bin_b2, bin_b2yes, extlw.bfx, activerb2, 2);
				if (loopb1 == 1)
					Go2b_Ext_Loop(BIT_SET_2X, 2);
				else  Go3(bin_b1, bin_b2yes);
				ExtSplitY(bin_b1,	extlr.tbfy, extlr.ntbfy,  activerb1);
			}
		}
		activeloop = activerb2; activeloop <<= 32; activeloop |= activerb1;
		if (extlr.noxyes < EXLBREAK) break;
	}
	Go3(bin_b1, bin_b2);// last call
}

void G17B::Go2b_Ext_Loop(uint64_t activeloop, uint32_t mode2) {
		//___________init the working areas bin2_b1, bin2_b2
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
			if (minratio >= EXLRATIO) break;
			else {
				if (extlw.mode == 1) {// this is a band1 X band2 Y
					ExtSplitX(bin2_b1, bin2_b1yes, extlw.bfx, activerb1);
					Go3(bin2_b1yes, bin2_b2);
					ExtSplitY(bin2_b2, extlw.tbfy, extlw.ntbfy, activerb2,2);
				}
				else {// this is a band2 X band1 Y
					ExtSplitX(bin2_b2, bin2_b2yes, extlw.bfx, activerb2,2);
					Go3(bin2_b1, bin2_b2yes);
					ExtSplitY(bin2_b1, extlw.tbfy, extlw.ntbfy, activerb1);
				}
			}
			activeloop = activerb2; activeloop <<= 32; activeloop |= activerb1;
			if (extlw.noxyes < EXLBREAK) break;
		}
		Go3(bin2_b1, bin2_b2);// last call
	}


//_________"step" _________ process a subset of valid {band1;band2}

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
	for (uint32_t iua = 0; iua < nua; iua++) {
		register uint64_t Ru = tua[iua];
		if (Ru&filter) continue;
		Ru &= Ra;
		if (ntusb1_128 < 128)		tusb1_128[ntusb1_128++] = Ru;
		else tusb1[ntusb1++] = Ru;
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
	stack_count_step.u64 = 0;
	uint64_t w = fb12;
	uint32_t xcell,cell;
	while (bitscanforward64(xcell, w)) {
		w ^= (uint64_t)1 << xcell;
		cell = From_128_To_81[xcell];
		tclues[nclues_step++] = cell;
		stack_count_step.u16[C_stack[cell]]++;
	}
	tcluesxy = &tclues[nclues_step];
	ub2.Init();// all to "no ua over 128"
	ntusb2 = 0;
	for (uint32_t iua = 0; iua < ntusb1; iua++) {
		register uint64_t Ru = tua[iua];
		if (!(Ru&filter)) {
			Ru &= Ra;
			if (!Ru) { return 1; }
			ub2.Add(Ru, ntusb2++);
			if (ntusb2 >= 2560)break;
		}
	}
	return 0;
}

void G17B::Go3(BINDEXN & bin1, BINDEXN & bin2) {
	if (aigstop) return;
	p_cpt2g[2]++;
	if (aigstop) return;
	if ((!bin1.ntvb) || (!bin2.ntvb)) return;
	//______________________________________ loops b1b2 main loops
	for (uint32_t ib1 = 0; ib1 < bin1.nt2; ib1++) {
		if (aigstop) return;
		p_cpt2g[3]++;
		morev2a.Init();	morev2b.Init();	morev2c.Init();
		Go3_Build_Band1(ib1, bin1);
		Go3_Apply_B1_V();
		tguas.ApplyLoopB1();
		for (uint32_t ib2 = 0; ib2 < bin2.nt2; ib2++) {
			if (aigstop) return;
			p_cpt2g[4]++;
			if (aigstop) return;
			Go3_Build_Band2(ib2, bin2);
			fb12 = fb1 | fb2;
			acb12 = acb1 | acb2;
			tguas.ApplyLoopB2();
			if (Go3_Apply_B2_V()) continue;
			Do128uas();
		}

	}
}

//______________main loop 128 first uas
void G17B::Do128uas() {//apply first 128 filter
	//morev2a.Skrink(fb12,acb12);	morev2b.Skrink(fb12,acb12); morev2c.Skrink(fb12, acb12);
	nza = nbi5_1;
	if (nza) { // 5 b1  6 b2  
		za = Z128_5_1;		zb = Z128_6_2;
		nzb = nbi6_2;		Do128Chunk();
	}
	nza = nbi6_1;
	if (nza) { // 6 b1  5 b2  
		za = Z128_6_1;		zb = Z128_5_2;
		nzb = nbi5_2;		Do128Chunk();
	}
}
void G17B::Do128Chunk() {
	if (!nzb) return;
	if (aigstop) return;
	p_cpt2g[26]++;
	p_cpt2g[49] += (nza * nzb);
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
	uint32_t  ideb2 = 0, iend2 = YCHUNK128;
	if (iend2 > n2)iend2 = n2;
	while (ideb2 < n2) { //Y chunk
		uint32_t ny = iend2 - ideb2;
		uint32_t ideb1 = 0, iend1 = XCHUNK128;
		if (iend1 > n1)iend1 = n1;

		while (ideb1 < n1) {// X chunk
			uint32_t nx = iend1 - ideb1;
			if (nx < ny)Do128Go(&z1[ideb1], &z2[ideb2], nx, ny);
			else Do128Go(&z2[ideb2], &z1[ideb1], ny, nx);
			ideb1 = iend1; iend1 += XCHUNK128;
			if (iend1 > n1)iend1 = n1;
		}
		ideb2 = iend2; iend2 += YCHUNK128;
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


//___________________ gangster specific

void TGUAS::ApplyLoopB1() {
	register uint64_t Bf = g17b.fb1;
	register uint64_t Ac = g17b.acb1 | BIT_SET_B2;
	// killed by known
	BF128 kk2, kk3; kk2.SetAll_0(); kk3.SetAll_0();
	{
		register uint64_t w = Bf;// known cells 
		uint32_t x;
		while (bitscanforward64(x, w)) {// here bf is <27
			w ^= (uint64_t)1 << x;
			kk2 |= genb12.kguas2[x];
			kk3 |= genb12.kguas3[x];
		}
	}
	for (uint32_t i = 0; i < 162; i++) {
		GUA & wg = tgua_start[i], &wgd = tgua_b1[i];
		wgd.nua = 0;
		if (i < 81) {
			if (kk2.On(i)) continue;// killed guas2
		}
		else {
			if (kk3.On(i - 81)) continue;// killed guas3
		}
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
				myb.tua128_b1[myb.ntua128_b1++] = w;
			}
		}
	}
}
void TGUAS::ApplyLoopB2() {
	register uint64_t Bf = g17b.fb12;
	register uint64_t Ac = g17b.acb12;
	// killed by known
	BF128 kk2, kk3; kk2.SetAll_0(); kk3.SetAll_0();
	{
		register uint64_t w = Bf >> 5;// known cells band 2 in 81 mode
		uint32_t x;
		while (bitscanforward64(x, w)) {// here bf is >=27
			w ^= (uint64_t)1 << x;
			kk2 |= genb12.kguas2[x];
			kk3 |= genb12.kguas3[x];
		}
	}
	nvg2 = nvg3 = 0;

	for (uint32_t igu = 0; igu < 162; igu++) {
		GUA & wg = tgua_b1[igu];
		uint64_t tua[64];
		uint32_t nua = 0;
		if (igu < 81) { if (kk2.On(igu)) continue; }// killed guas2
		else 	if (kk3.On(igu - 81)) continue;// killed guas3

		for (uint32_t j = 0; j < wg.nua; j++) {// apply new subsets
			register uint64_t Ua = wg.tua[j];
			if (Ua & Bf) continue;
			Ua &= Ac;// could be empty
			if (!Ua) { tua[0] = 0;		nua = 1;	break; }
			tua[nua++] = Ua;
		}
		if (nua) {	// split per size in 54 + index 0-162
			for (uint32_t j = 0; j < nua; j++) {
				register uint64_t Ua = tua[j],
					cc = _popcnt64(Ua),
					Ua54 = (Ua & BIT_SET_27) | ((Ua & BIT_SET_B2) >> 5);
				if (igu < 81)AddVG2_128(Ua54, igu);// this is a g2					
				else AddVG3_128(Ua54, igu);
			}
		}
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
	}
}
int TGUAS::ApplyG2_128() {
	g2ok.SetAll_0();
	nb64_1 = (nvg2 + 127) >> 7;
	for (uint32_t iv = 0; iv < nb64_1; iv++) {
		TVG128 &vv = tvg128g2[iv];
		BF128 V = vv.v;
		for (int j = 0; j < g17b.nclues; j++)
			V &= vv.cells[g17b.tcluesxy[j]];
		int bit;
		while ((bit = V.getFirst128()) >= 0) {
			V.Clear(bit);// clear bit
			g2ok.Set(vv.ti162[bit]);
		}
	}
	int cc = g2ok.Count();
	return cc;
}
int TGUAS::ApplyG3_128() {
	g3ok.SetAll_0();
	nb64_1 = (nvg3 + 127) >> 7;
	for (uint32_t iv = 0; iv < nb64_1; iv++) {
		TVG128 &vv = tvg128g3[iv];
		BF128 V = vv.v;
		for (int j = 0; j < g17b.nclues; j++)
			V &= vv.cells[g17b.tcluesxy[j]];
		int bit;
		while ((bit = V.getFirst128()) >= 0) {
			V.Clear(bit);// clear bit
			g3ok.Set(vv.ti162[bit] - 81);
		}
	}
	return g3ok.isNotEmpty();
}


//___________ potential valid bands 1+2  after 128 uas

void G17B::CleanAll() {
	nwc = n_to_clean;	n_to_clean = 0;
	if (aigstop)  return;
	p_cpt2g[5]++;	
	if (!nwc) return;

	for (uint64_t i = 0; i < nwc; i++) {
		if (aigstop) return;
		aigstopxy = 0;
		wb12bf = to_clean[i];
		nclues = 0;
		{
			register uint64_t w = wb12bf ^ fb12;
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
			stack_count.u16[2] > 6) continue;
		if (ub2.ApplyXY(tcluesxy, nclues, ntusb2)) continue;
		p_cpt2g[6]++;
		if (morev2a.ApplyXY(tclues, nclues + nclues_step)) continue;
		if (morev2b.ApplyXY(tclues, nclues + nclues_step)) continue;
		if (morev2c.ApplyXY(tclues, nclues + nclues_step)) continue;
		clean_valid_done = 0;
		if (genb12.nband3 > 20) { if (Clean_2a()) continue; }// valid first
		else if (genb12.nband3 > 5) {	if (Clean_2b()) continue;}// guas2 then valid
		else if (Clean_2c()) continue; // only guas2
		int isg3 = tguas.ApplyG3_128(), minr;
		for (int ib3= ib3_current; ib3 < genb12.nband3; ib3++) {
			STD_B3 & b = genb12.bands3[ib3];
			if(ib3!= ib3_current)	if (b.CleanG2() < 0) continue;
			if(isg3)if ((minr = b.CleanG3()) < 0) continue;
			zhou[0].PartialInitSearch17(tclues, nclues + nclues_step);
			GoB3(b);
			if (aigstopxy)break;// ua 12 added
		}
	}
}
int G17B::Clean_Valid() {
	clean_valid_done = 1;
	myua = zh2b[0].ValidXY(tclues, nclues + nclues_step);
	if (myua) { NewUaB12();	return 1; }
	return 0;
}
int G17B::Clean_2a() {// this is for a given band 1+2
	p_cpt2g[23]++;
	if (g17b.Clean_Valid()) return 1;  
	p_cpt2g[24]++;
	tguas.ApplyG2_128();
	for (ib3_current = 0; ib3_current < genb12.nband3; ib3_current++) 
		if (genb12.bands3[ib3_current].CleanG2() >= 0) return 0;
	return 1;
}
int G17B::Clean_2b() {// this is for a given band 1+2
	tguas.ApplyG2_128();
	for (ib3_current = 0; ib3_current < genb12.nband3; ib3_current++)
		if (genb12.bands3[ib3_current].CleanG2() >= 0) {
			p_cpt2g[23]++;
			if (g17b.Clean_Valid()) return 1;
			p_cpt2g[24]++;
			return 0;
		}
	return 1;
}
int G17B::Clean_2c() {// this is for a given band 1+2
	tguas.ApplyG2_128();
	for (ib3_current = 0; ib3_current < genb12.nband3; ib3_current++)
		if (genb12.bands3[ib3_current].CleanG2() >= 0) 
			return 0;
	return 1;
}

int STD_B3::CleanG2() {// critical code guas2
	memset(&smin, 0, sizeof smin);
	BF128 g2 = tguas.g2ok&guas.isguasocketc2;
	if (g2.Count() >12)return -1;// 
	int cc81;
	while ((cc81 = g2.getFirst128()) >= 0) {
		g2.Clear(cc81);		Insert2(cc81);
	}
	sminr = smin;
	smin.SetMincount();
	if (smin.mincount > 6)return -1;
	stack_count.u64 = g17b.stack_count.u64 + smin.Count_per_stack().u64;
	if (stack_count.u16[0] > 6 || stack_count.u16[1] > 6 ||
		stack_count.u16[2] > 6) return -1; // not ok
	return smin.mincount;
}
int STD_B3::CleanG3() {// guas 3
	BF128 g3 = tguas.g3ok & guas.isguasocketc3;
	if(g3.isEmpty()) return smin.mincount;
	smin = sminr;
	int cc81;
	while ((cc81 = g3.getFirst128()) >= 0) {
		g3.Clear(cc81);
		Insert3(cc81);
	}
	smin.SetMincount();

	if (smin.mincount > 6)return -1;
	stack_count.u64 = g17b.stack_count.u64 + smin.Count_per_stack().u64;
	if (stack_count.u16[0] > 6 || stack_count.u16[1] > 6 ||
		stack_count.u16[2] > 6) return -1; // not ok
	return smin.mincount;
}

//______________________ start process final b3
void G17B::GoB3(STD_B3 & b) {
	myband3 = &b;
	moreuas_b3.Init();
	memcpy(&genb12.grid0[54], b.band0, 4 * 27);
	memset(&hh0, 0, sizeof hh0);
	stack_countf = b.stack_count;
	smin = b.smin;
	uint32_t nmiss = 6 - smin.mincount;
	register uint32_t Fstk = BIT_SET_27;
	if (stack_countf.u16[0] == 6) Fstk ^= 07007007;
	if (stack_countf.u16[1] == 6) Fstk ^= 070070070;
	if (stack_countf.u16[2] == 6) Fstk ^= 0700700700;
	wactive0 = fstk = Fstk;
	nuasb3_1 = nuasb3_2 = 0;
	// load uas guas 2 guas 3 in field
	int cc81;
	{
		BF128 g2 = tguas.g2ok & b.guas.isguasocketc2;
		while ((cc81 = g2.getFirst128()) >= 0) {
			g2.Clear(cc81);
			uasb3_1[nuasb3_1++] =b.guas.ua_pair[cc81];
		}
	}
	{
		BF128 g3 = tguas.g3ok&b.guas.isguasocketc3;
		while ((cc81 = g3.getFirst128()) >= 0) {
			g3.Clear(cc81);
			uasb3_1[nuasb3_1++] = b.guas.ua_triplet[cc81];
		}
	}
	// build out field and And out
	register uint32_t If = smin.critbf; // in field
	register uint32_t andout = Fstk & (~If);
	for (uint32_t i = 0; i < b.ntua128_b2; i++) {
		BF128 w = b.tua128_b2[i];
		if (w.bf.u64[0] & wb12bf) continue;
		register uint32_t u3 = w.bf.u32[2];
		if (u3&If)uasb3_1[nuasb3_1++] = u3;
		else {
			uasb3_2[nuasb3_2++] = u3 & Fstk;
			andout &= u3;
		}
	}
	BF128 w = tguas.g2ok &b.guas.isguasocketc2_46;
	int i81;
	while ((i81 = w.getFirst128()) >= 0) {
		w.Clear(i81);
		register uint32_t u3 = b.guas.ua_pair[i81];
		if (u3&If)uasb3_1[nuasb3_1++] = u3;
		else {
			uasb3_2[nuasb3_2++] = u3 & Fstk;
			andout &= u3;
		}
	}
	for (uint32_t i = 0; i < b.nua; i++) {
		register uint32_t u3 = b.tua[i];
		if (u3&If)uasb3_1[nuasb3_1++] = u3;
		else {
			uasb3_2[nuasb3_2++] = u3 & Fstk;
			andout &= u3;
		}
	}
	if (!nmiss) {
		if (nuasb3_2)return;
		else {p_cpt2g[10]++;hh0.GoMiss0(); return;}
	}
	if (nmiss == 1) {
		if (!andout)return;
		else {
			p_cpt2g[11]++;
			hh0.GoMiss1(andout);
			return;
		}
	}

	if (nmiss == 2 && nuasb3_2 > 1 && (!andout)) {// want always 2 clues 
		// start with the smallest ua next will be "and" of remaining uas
		p_cpt2g[12]++;
		hh0.GoMiss2Init();
		uint32_t  uamin = uasb3_2[0];
		{
			register uint32_t min = _popcnt32(uamin);
			for (uint32_t i = 1; i < nuasb3_2; i++) {
				register uint32_t Ru = uasb3_2[i], cc = _popcnt32(Ru);
				if (cc < min) {	min = cc;	uamin = Ru;		}
			}
		}
		hh0.GoMiss2( uamin);
		return;
	}
	if (nuasb3_2) {// all uas b3 same table
		memcpy(&uasb3_1[nuasb3_1], uasb3_2, nuasb3_2 * sizeof uasb3_1[0]);
		nuasb3_1 += nuasb3_2;
	}
	ExpandB3();
}

//__________ phase 2___ find band 3 clues for one band 3
void G17B3HANDLER::GoMiss0() {
	smin = g17b.smin;
	uasb3if = g17b.uasb3_1;
	nuasb3if = g17b.nuasb3_1;
	active_b3 = smin.critbf;
	known_b3 = rknown_b3 = 0;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	CriticalLoop();
}

void G17B3HANDLER::GoMiss1(uint32_t andout) {
	nmiss = 1;
	smin = g17b.smin;
	stack_count = g17b.stack_countf;
	uasb3if = g17b.uasb3_1;
	nuasb3if = g17b.nuasb3_1;
	nuasb3of = g17b.nuasb3_2;
	active_b3 = smin.critbf;
	known_b3 = rknown_b3 = 0;
	wactive0 = g17b.wactive0;
	wua = andout;
	Do_miss1();
}
void G17B3HANDLER::Do_miss1(){
	if (!nuasb3of) {// subcritical in hn if solved
		if (!g17b.clean_valid_done) 
			if (g17b.Clean_Valid()) {
				g17b.aigstopxy = 1;
				return; // b 12 not valid
			}
		int uabr = IsMultiple(active_b3);

		if (uabr)wua &= uabr ; // one ua outfield seen
		else {// confirmed subcritical possible
			G17B3HANDLER hn = *this;
			hn.Critical2pairs(1);// only in critical stacks
			hn.Go_Subcritical();
			wua &= ~active_b3; // don't re use this as first cell
		}
	}
	Critical2pairs();// assign 2 pairs in minirow to common cell
	uint32_t res;
	if (!g17b.moreuas_b3.CheckNew(known_b3, wua))return;
	while (bitscanforward(res, wua)) {
		int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
		G17B3HANDLER hn = *this;
		hn.known_b3 |= bit;
		g17b.ua_of_seen = 0;
		hn.CriticalLoop();
		if (g17b.aigstopxy ) return ;
		if (g17b.ua_of_seen)wua &= g17b.ua_of_seen;
	}
}

void G17B3HANDLER::GoMiss2Init() {
	smin = g17b.smin;
	stack_count = g17b.stack_countf;
	active_b3 = smin.critbf;
	known_b3 = rknown_b3 = 0;
	wactive0 = g17b.wactive0 & (BIT_SET_27 ^ active_b3);//  active cells out field
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
void G17B3HANDLER::AddCell_Miss2(uint32_t cell, uint32_t wand) {
	nuasb3of = 1;
	wua = wand;
	{
		register int s = C_stack[cell];
		stack_count.u16[s]++;
		if (stack_count.u16[s] > 5) {
			s = ~(07007007 << (3 * s));// mask
			wua &= s;
			wactive0 &= s;
		}
	}
	nmiss--;
	known_b3 |= 1 << cell;
	Do_miss1();
}
void G17B3HANDLER::GoMiss2( uint32_t uamin) {
	nmiss = 2;
	uasb3if = g17b.uasb3_1;
	nuasb3if = g17b.nuasb3_1;
	uasb3of = g17b.uasb3_2;
	wua = uamin;
	// cells added must produce cells hitting all remaining uas
	uint32_t res ;
	while (bitscanforward(res, wua)) {
		register uint32_t  bit = 1 << res;
		wua ^= bit; 
		wactive0 ^= bit;
		register uint32_t andx = wactive0, s = C_stack[res];
		if (stack_count.u16[s] == 5) {
			s = ~(07007007 << (3 * s));// mask
			andx &= s;
		}
		for (uint32_t i = 0; i < nuasb3of; i++) {
			register uint32_t ua = uasb3of[i];
			if (!(ua&bit)) {andx &= ua;	if (!andx) break;	}
		}
		if (andx) {
			G17B3HANDLER hn = *this;
			AddCell_Miss2(res, andx);
			if (g17b.aigstopxy) return;
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
		if (!modesub)smin.mini_bf2 = 0;
	}
}
void G17B3HANDLER::CriticalLoop() {// after optional assignment
	while (1) {// first shrink uas in field
		irloop = 0;
		uint32_t * tn = &uasb3if[nuasb3if], n = 0;
		register uint32_t Ra = active_b3,
			Rfilt = known_b3;
		for (uint32_t iua = 0; iua < nuasb3if; iua++) {
			register int Ru = uasb3if[iua];
			if (Ru & known_b3) continue;// already hit, forget it
			Ru &= active_b3;
			if (!Ru) return;// dead branch
			if (_popcnt32(Ru) == 1) {// assign it and reduce the active cells
				CriticalAssignCell(Ru);
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
		g17b.FinalCheckB3(known_b3);
		return; // branch closed
	}
	uint32_t wua = uasb3if[0] & active_b3, cell;
	while (bitscanforward(cell, wua)) {
		register int bit = 1 << cell;
		wua ^= bit;// clear bit

		// clean the bit in active_b3, this is now a dead cell downstream
		active_b3 ^= bit;
		G17B3HANDLER hn = *this;
		hn.CriticalAssignCell(bit);
		hn.CriticalLoop();
		if (g17b.aigstopxy || g17b.ua_of_seen) return;
	}
}

void G17B::ExpandB3(){
	p_cpt2g[13]++;
	uint32_t *tuaw=uasb3_1, nuaw = nuasb3_1;// use  in field pre loaded in the right way
	struct SPB3 {// spots to find band 3 minimum valid solutions
		uint32_t  possible_cells, all_previous_cells, active_cells, iuab3,
			stack[3];
	}spb3[7], *s3, *sn3;
	s3 = spb3;
	memset(s3, 0, sizeof spb3[0]); // previous, possible,iuab3
	s3->active_cells = BIT_SET_27;// all cells active
	// init the stack status
	int tcells[10];
	{
		register uint32_t bf2 = smin.mini_bf2, kb3 = 0;
		for (int i = 0; i < 3; i++) {
			s3->stack[i] = stack_count.u16[i];// count before band 3 min count
			if (s3->stack[i] == 6)s3->active_cells &= ~(07007007 << (3 * i));			
			else if (bf2) {
				uint32_t st = bf2 & (0111 << i);// stack bf2
				bf2 ^= st; // clear stack in bf2
				if (stack_countf.u16[i] != 6) continue;
				if (st) {// stack to assign
					int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000,
						07000007, 07007000, 07007007 };
					int shrink = TblShrinkMask[smin.mini_bf2 & (0111 << i)];
					if (shrink) {// minirows 2 pairs in that stack
						register int Mask = tbitsrows[shrink] << (3 * i);
						//s3->active_cells &= ~(07007007 << (3 * i));
						kb3 |= Mask & (~smin.pairs27);// and set the common cell as assigned
					}
				}
			}
			
		}
		if (kb3) {//cells to assign (critical 2 pairs)
			// shrink the uas table and assign
			uint32_t *tuawn =&tuaw[nuaw], nuawn = 0,cell,ispot=0;
			for (uint32_t i = 0; i < nuaw; i++) {
				register uint32_t w = tuaw[i];
				if (!(w & kb3)) {
					w &= s3->active_cells;
					tuawn[nuawn++] = w ;
				}
			}
			tuaw = tuawn;
			nuaw = nuawn;
			while (bitscanforward(cell, kb3)) {
				register int bit = 1 << cell;
				kb3 ^= bit;
				tcells[ispot++ ] = cell;
				sn3 = s3 + 1; *sn3 = *s3; // (copy the stack count)
				int st = C_stack[cell];
				sn3->stack[st]++;
				if (sn3->stack[st] > 5)// stack is limit update sn3 active
					sn3->active_cells &= ~(07007007 << (3 * st));
				sn3->all_previous_cells |= bit;
				s3 = sn3;
			}
		}
	}
	s3->possible_cells = tuaw[0] & s3->active_cells;

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
	if (!clean_valid_done) 
		if (Clean_Valid()) {
			moreuas_b3.Add(0);//lock the call 
			aigstopxy = 1;		return;
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

int ZHOU::CallMultipleB3(ZHOU & o, uint32_t bf, int diagx) {
	*this = o;
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
	if (_popcnt32(digitsbf < 8)) 	return 1;// can not be one solution

	zh_g2.s17_b3_mini = 1;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	//__________end assign last lot start solver
	zh_g.go_back = 0;	zh_g.nsol = 0; zh_g.lim = 1;// modevalid is set to  1
	zh_g2.isfalse_on = -1;
	int ir = Full17Update();
	if (ir == 2) return 0;// solved can not be multiple
	Guess17(0);

	return zh_g.nsol;
}
int ZHOU::Apply17SingleOrEmptyCellsB12() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched
	register uint64_t R1, R2, R3, R4;
	{
		register uint64_t * P = FD[0][0].bf.u64, M = *P;
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
	bitscanforward64(zh_g2.xcell_to_guess,R2);
	return 0;
}
int ZHOU::Apply17SingleOrEmptyCellsB3() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched
	register uint32_t R1, R2, R3, R4;
	{
		register uint32_t * P =& FD[0][0].bf.u32[2], M = *P;
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
			uint32_t cell = xcell+54;
			xcell+= 64;
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
		R3 &= ~R4; // now true singles
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
			if ((!zh_g2.isfalse_on))return 0;
			if (Apply17SingleOrEmptyCellsB12())	return 0;
		}
		if (!zh_g.single_applied)	break;
	}
	return 1;
}
void ZHOU::Compute17Next(int index) {
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
			if (wua.isNotEmpty()) {
				zh_g.nsol++;
				zh_g.go_back = 1;// closed anyway
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
		digit = zh_g2.grid0[cell];
	// true first if possible
	if (FD[digit][0].On(xcell)) {
		ZHOU * mynext = (this + 1);
		*mynext = *this;
		mynext->SetaCom(digit, cell, xcell);
		mynext->Compute17Next(index + 1);
		if (zh_g.go_back) return;
	}
	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue;
		if (FD[idig][0].On(xcell)) {
			if (cell >= 54)zh_g2.isfalse_on = 1;
			ZHOU * mynext = (this + 1);
			*mynext = *this;
			mynext->SetaCom(idig, cell, xcell);
			mynext->Compute17Next(index + 1);
			if (zh_g.go_back) return;
		}
	}
}


//================= critical process
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
uint32_t G17B3HANDLER::IsMultiple(int bf) {
	if (bf == rknown_b3) return 0;
	uint32_t ua = 0;
	rknown_b3 = bf;
	G17B & bab = g17b;
	// check first if all tuab3 is hit
	int ir = zhou[1].CallMultipleB3(zhou[0], bf, 0);
	if (ir) {
		ua = zh_g2.cells_assigned.bf.u32[2];
		g17b.moreuas_b3.Add(ua);
		g17b.NewUaB3();
	}
	return ua;
}



//=============== sub critical process   missing(s)  in the critical area
void G17B3HANDLER::Go_SubcriticalMiniRow() {
	int c2[3] = { 3, 5, 6 };// 2 cells in a mini row
	int bit = 1 << ndead, mask = 7 << (3 * ndead);
	for (int i = ndead; i < 9; i++,  bit <<= 1, mask <<= 3) {
		stack = i % 3;
		register int M = active_sub & mask;
		if (!M)continue;
		ndead = i;
		if (bit & smin.mini_bf1) {// it was a gua2 pair assign both
			G17B3HANDLER hn = *this;
			hn.smin.mini_bf1 ^= bit;
			hn.SubMini( M, mask);
			if(!g17b.moreuas_b3.CheckNew(known_b3, active_sub))return;

		}
		else if (bit & smin.mini_bf2)// it was 2 gua2 pair assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				G17B3HANDLER hn = *this;
				hn.smin.mini_bf2 ^= bit;
				hn.SubMini(M, mask);
				if (!g17b.moreuas_b3.CheckNew(known_b3, active_sub))return;
			}
		else if (bit & smin.mini_bf3) {// it was 3 gua2 pair assign 3 out of 3
			G17B3HANDLER hn = *this;
			hn.smin.mini_bf3 ^= bit;
			hn.SubMini(M, mask);
			if (!g17b.moreuas_b3.CheckNew(known_b3, active_sub))return;
		}
		else if (bit & smin.mini_triplet)// it was a gua3 triplet assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				G17B3HANDLER hn = *this;
				hn.smin.mini_triplet ^= bit;
			if(!g17b.moreuas_b3.CheckNew(known_b3, active_sub))return;
				hn.SubMini(M, mask);
				if (!g17b.moreuas_b3.CheckNew(known_b3, active_sub))return;
			}
		else {// second add in the mini row one residual cell take it
			G17B3HANDLER hn = *this;
			hn.SubMini(M, mask);
			if (!g17b.moreuas_b3.CheckNew(known_b3, active_sub))return;
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
	if (nmiss) Go_SubcriticalMiniRow();// continue till a"no missing clue condition"
	else {	// leave sub critical mode and enter the critical mode
		Critical2pairs();// assign 2 pairs in minirow to common cell
		CriticalLoop();
	}
}
void G17B3HANDLER::Go_Subcritical() {// nmiss to select in the critical field
	active_b3 = active_sub = smin.critbf;
	// check first if a global solution  is still possible
	for (int ist = 0; ist < 3; ist++) {// check stacks
		if (stack_count.u16[ist] > 5)active_sub &= ~(07007007 << (3 * ist));// kill the stack for more clues
	}
	ndead = 0;
	Go_SubcriticalMiniRow();// find the first miss
}


//________ final called by all branches
void G17B::FinalCheckB3(uint32_t bfb3) {
	p_cpt2g[29]++;
	if (!clean_valid_done) {
		if (Clean_Valid()) {
			moreuas_b3.Add(0);//lock the call 
			aigstopxy = 1;		return;
		}
	}
	if (moreuas_b3.Check(bfb3))return;
	register uint32_t ir = zhou[1].CallMultipleB3(zhou[0], bfb3, 0);
	if (ir) { NewUaB3();	return; }
	Out17(bfb3);
}
void G17B::Out17(uint32_t bfb3) {
	cout << Char27out(bfb3) << "\t\tone sol to print final check " << endl;
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
	if (cc64 < 12) {// this should never be this is a check for a bug
		DebugAdd12();		return;
	}
	if (cc64 < UALIMSIZE) {
		if (ntusb2 < 2560)ub2.Add(myua, ntusb2++);
		if (genuasb12.nua < (TUA64_12SIZE - 500) || cc64 < 16) {
			register uint64_t ua_add = myua | (cc64 << 59);
			genuasb12.AddUACheck(ua_add);
			if (ntusb1 < TUA64_12SIZE)tusb1[ntusb1++] = myua;
			if (ntusb1 > p_cpt2g[17])p_cpt2g[17] = ntusb1;
			if (genuasb12.nua > p_cpt2g[16])p_cpt2g[16] = genuasb12.nua;
			p_cpt2g[31]++;
		}
		else morev2a.Add(myua);
	}
	else if (cc64 < (UALIMSIZE + 1)) { morev2a.Add(myua); p_cpt2g[32]++; }
	else if (cc64 < (UALIMSIZE + 2)) { morev2b.Add(myua); p_cpt2g[33]++; }
	else {
		p_cpt2g[34]++;
		morev2c.Add(myua);
	}
}
void G17B::NewUaB3() {// new ua from final check zh_g2.cells_assigned
	BF128 ua128 = zh_g2.cells_assigned;
	register uint64_t ua12 = ua128.bf.u64[0];
	register uint32_t ua = ua128.bf.u32[2],
		cc = _popcnt32(ua),		cc0 = (uint32_t)_popcnt64(ua12);
	if (!cc) {// bands 1+2 not valid
		cerr << " uab12 in newuab3" << endl;
		cout << " uab12 in newuab3 nmiss="<<hh0.nmiss << endl;
		DebugAdd12();		aigstop = 1;
		return;
	}
	if (!(ua&hh0.smin.critbf))ua_out_seen = ua;
	moreuas_b3.Add(ua);
	if (cc0 > GUALIMSIZE) return; 
	if (cc > 3) {
		if ((cc0 + cc) > 15) return;
		if (cc == 4 && cc0 > 15) return;
	}
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
	if (cc > 3) {
		p_cpt2g[42]++;
		// could also be added later to other bands 3 where it is valid
		if (myband3->ntua128 < 1000) {
			myband3->tua128[myband3->ntua128++] = ua128;
			myband3->tua128_b1[myband3->ntua128_b1++] = ua128;
			myband3->tua128_b2[myband3->ntua128_b2++] = ua128;
		}
		return;
	}

	if (cc == 2) {// one of the 27 GUA2s add to the table
		p_cpt2g[40]++;
		///NewUaB3_g2(my_i81, ua12);
		tguas.tgua_start[my_i81].Adduacheck(ua12);// for new steps
		tguas.tgua_b1[my_i81].Adduacheck(ua12);// for new steps
		tguas.g2ok.Set(my_i81);
		if (tguas.nvg2 < 1280) {
			uint32_t ibloc = tguas.nvg2 >> 7, ir = tguas.nvg2 - 128 * ibloc;
			tvg128g2[ibloc].SetVect54(ua54, ir, my_i81);
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
		if (tguas.nvg3 < 640) {
			uint32_t ibloc = tguas.nvg3 >> 7, ir = tguas.nvg3 - 128 * ibloc;
			tvg128g3[ibloc].SetVect54(ua54, ir, my_i81);
			tguas.nvg3++;
		}
		return;
	}

}


void G17B::DebugAdd12() {
	aigstop = 1;
	cerr << "ua < 12 to add clean" << endl;
	cout << endl << endl << Char2Xout(myua) << " ua < 12 to add   clean" << endl;
	cout << "bug location band 2 id=" << genb12.nb12 << endl;
	cout << "p_cpt2g[0]=" << p_cpt2g[0] << "\tp_cpt2g[2]=" << p_cpt2g[2];
	cout << "\tp_cpt2g[3]=" << p_cpt2g[3] << "\tp_cpt2g[6]=" << p_cpt2g[6] << endl;
	cout << Char2Xout(wb12bf) << " b12 at call" << endl;
	cout << "ntusb1=" << ntusb1 << " ntoclean=" << n_to_clean << endl;
	for (int i = 0; i < nclues_step; i++) cout << tclues[i] << " ";
	cout << "\t";
	for (int i = 0; i < nclues; i++) cout << tcluesxy[i] << " ";
	cout << endl;
	zh2b[0].ImageCandidats();
}
