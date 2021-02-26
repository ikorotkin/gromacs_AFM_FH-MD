/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Sat Jul 30 16:41:20 EDT 2016 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_notw.native -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -n 16 -name n2sv_16 -with-ostride 1 -include n2s.h -store-multiple 4 */

/*
 * This function contains 144 FP additions, 40 FP multiplications,
 * (or, 104 additions, 0 multiplications, 40 fused multiply/add),
 * 110 stack variables, 3 constants, and 72 memory accesses
 */
#include "n2s.h"

static void n2sv_16(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DVK(KP414213562, +0.414213562373095048801688724209698078569671875);
     {
	  INT i;
	  for (i = v; i > 0; i = i - (2 * VL), ri = ri + ((2 * VL) * ivs), ii = ii + ((2 * VL) * ivs), ro = ro + ((2 * VL) * ovs), io = io + ((2 * VL) * ovs), MAKE_VOLATILE_STRIDE(64, is), MAKE_VOLATILE_STRIDE(64, os)) {
	       V T2p, T2q, T2r, T2s, T2x, T2y, T2z, T2A, T1M, T1N, T1L, T1P, T2F, T2G, T2H;
	       V T2I, T1O, T1Q;
	       {
		    V T1l, T1H, T1R, T7, T1x, TN, TC, T25, T1E, T1b, T1Z, Tt, T2h, T22, T1D;
		    V T1g, T1n, TQ, T11, Ti, Te, T26, T1m, TT, T1S, TJ, TZ, T1V, TW, Tl;
		    V T12, T13;
		    {
			 V Tq, T1c, Tp, T20, T1a, Tr, T1d, T1e;
			 {
			      V T1, T2, Tw, Tx, T4, T5, Tz, TA;
			      T1 = LD(&(ri[0]), ivs, &(ri[0]));
			      T2 = LD(&(ri[WS(is, 8)]), ivs, &(ri[0]));
			      Tw = LD(&(ii[0]), ivs, &(ii[0]));
			      Tx = LD(&(ii[WS(is, 8)]), ivs, &(ii[0]));
			      T4 = LD(&(ri[WS(is, 4)]), ivs, &(ri[0]));
			      T5 = LD(&(ri[WS(is, 12)]), ivs, &(ri[0]));
			      Tz = LD(&(ii[WS(is, 4)]), ivs, &(ii[0]));
			      TA = LD(&(ii[WS(is, 12)]), ivs, &(ii[0]));
			      {
				   V Tn, TL, T3, T1k, Ty, T1j, T6, TM, TB, To, T18, T19;
				   Tn = LD(&(ri[WS(is, 15)]), ivs, &(ri[WS(is, 1)]));
				   TL = VSUB(T1, T2);
				   T3 = VADD(T1, T2);
				   T1k = VSUB(Tw, Tx);
				   Ty = VADD(Tw, Tx);
				   T1j = VSUB(T4, T5);
				   T6 = VADD(T4, T5);
				   TM = VSUB(Tz, TA);
				   TB = VADD(Tz, TA);
				   To = LD(&(ri[WS(is, 7)]), ivs, &(ri[WS(is, 1)]));
				   T18 = LD(&(ii[WS(is, 15)]), ivs, &(ii[WS(is, 1)]));
				   T19 = LD(&(ii[WS(is, 7)]), ivs, &(ii[WS(is, 1)]));
				   Tq = LD(&(ri[WS(is, 3)]), ivs, &(ri[WS(is, 1)]));
				   T1l = VADD(T1j, T1k);
				   T1H = VSUB(T1k, T1j);
				   T1R = VSUB(T3, T6);
				   T7 = VADD(T3, T6);
				   T1x = VADD(TL, TM);
				   TN = VSUB(TL, TM);
				   TC = VADD(Ty, TB);
				   T25 = VSUB(Ty, TB);
				   T1c = VSUB(Tn, To);
				   Tp = VADD(Tn, To);
				   T20 = VADD(T18, T19);
				   T1a = VSUB(T18, T19);
				   Tr = LD(&(ri[WS(is, 11)]), ivs, &(ri[WS(is, 1)]));
				   T1d = LD(&(ii[WS(is, 3)]), ivs, &(ii[WS(is, 1)]));
				   T1e = LD(&(ii[WS(is, 11)]), ivs, &(ii[WS(is, 1)]));
			      }
			 }
			 {
			      V Tb, Ta, TF, Tc, TG, TH, TP, TO;
			      {
				   V T8, T9, TD, TE;
				   T8 = LD(&(ri[WS(is, 2)]), ivs, &(ri[0]));
				   T9 = LD(&(ri[WS(is, 10)]), ivs, &(ri[0]));
				   TD = LD(&(ii[WS(is, 2)]), ivs, &(ii[0]));
				   TE = LD(&(ii[WS(is, 10)]), ivs, &(ii[0]));
				   Tb = LD(&(ri[WS(is, 14)]), ivs, &(ri[0]));
				   {
					V T17, Ts, T21, T1f;
					T17 = VSUB(Tq, Tr);
					Ts = VADD(Tq, Tr);
					T21 = VADD(T1d, T1e);
					T1f = VSUB(T1d, T1e);
					TP = VSUB(T8, T9);
					Ta = VADD(T8, T9);
					TO = VSUB(TD, TE);
					TF = VADD(TD, TE);
					T1E = VSUB(T1a, T17);
					T1b = VADD(T17, T1a);
					T1Z = VSUB(Tp, Ts);
					Tt = VADD(Tp, Ts);
					T2h = VADD(T20, T21);
					T22 = VSUB(T20, T21);
					T1D = VADD(T1c, T1f);
					T1g = VSUB(T1c, T1f);
					Tc = LD(&(ri[WS(is, 6)]), ivs, &(ri[0]));
				   }
				   TG = LD(&(ii[WS(is, 14)]), ivs, &(ii[0]));
				   TH = LD(&(ii[WS(is, 6)]), ivs, &(ii[0]));
			      }
			      T1n = VADD(TP, TO);
			      TQ = VSUB(TO, TP);
			      {
				   V Tg, Th, TX, TR, Td, TS, TI, TY, Tj, Tk;
				   Tg = LD(&(ri[WS(is, 1)]), ivs, &(ri[WS(is, 1)]));
				   Th = LD(&(ri[WS(is, 9)]), ivs, &(ri[WS(is, 1)]));
				   TX = LD(&(ii[WS(is, 1)]), ivs, &(ii[WS(is, 1)]));
				   TR = VSUB(Tb, Tc);
				   Td = VADD(Tb, Tc);
				   TS = VSUB(TG, TH);
				   TI = VADD(TG, TH);
				   TY = LD(&(ii[WS(is, 9)]), ivs, &(ii[WS(is, 1)]));
				   Tj = LD(&(ri[WS(is, 5)]), ivs, &(ri[WS(is, 1)]));
				   T11 = VSUB(Tg, Th);
				   Ti = VADD(Tg, Th);
				   Tk = LD(&(ri[WS(is, 13)]), ivs, &(ri[WS(is, 1)]));
				   Te = VADD(Ta, Td);
				   T26 = VSUB(Td, Ta);
				   T1m = VSUB(TR, TS);
				   TT = VADD(TR, TS);
				   T1S = VSUB(TF, TI);
				   TJ = VADD(TF, TI);
				   TZ = VSUB(TX, TY);
				   T1V = VADD(TX, TY);
				   TW = VSUB(Tj, Tk);
				   Tl = VADD(Tj, Tk);
				   T12 = LD(&(ii[WS(is, 5)]), ivs, &(ii[WS(is, 1)]));
				   T13 = LD(&(ii[WS(is, 13)]), ivs, &(ii[WS(is, 1)]));
			      }
			 }
		    }
		    {
			 V T2f, Tf, T2j, TK, Tm, T1U, T10, T1B, T14, T1W;
			 T2f = VSUB(T7, Te);
			 Tf = VADD(T7, Te);
			 T2j = VADD(TC, TJ);
			 TK = VSUB(TC, TJ);
			 Tm = VADD(Ti, Tl);
			 T1U = VSUB(Ti, Tl);
			 T10 = VADD(TW, TZ);
			 T1B = VSUB(TZ, TW);
			 T14 = VSUB(T12, T13);
			 T1W = VADD(T12, T13);
			 {
			      V T29, T1T, T27, T2d, T2b, T23, T15, T1A, T2l, T2m, T2n, T2o, T2i, T2k, T1Y;
			      V T2a;
			      {
				   V Tv, Tu, T1X, T2g;
				   T29 = VSUB(T1R, T1S);
				   T1T = VADD(T1R, T1S);
				   T27 = VSUB(T25, T26);
				   T2d = VADD(T26, T25);
				   T2b = VADD(T1Z, T22);
				   T23 = VSUB(T1Z, T22);
				   Tv = VSUB(Tt, Tm);
				   Tu = VADD(Tm, Tt);
				   T1X = VSUB(T1V, T1W);
				   T2g = VADD(T1V, T1W);
				   T15 = VSUB(T11, T14);
				   T1A = VADD(T11, T14);
				   T2l = VSUB(TK, Tv);
				   STM4(&(io[12]), T2l, ovs, &(io[0]));
				   T2m = VADD(Tv, TK);
				   STM4(&(io[4]), T2m, ovs, &(io[0]));
				   T2n = VADD(Tf, Tu);
				   STM4(&(ro[0]), T2n, ovs, &(ro[0]));
				   T2o = VSUB(Tf, Tu);
				   STM4(&(ro[8]), T2o, ovs, &(ro[0]));
				   T2i = VSUB(T2g, T2h);
				   T2k = VADD(T2g, T2h);
				   T1Y = VADD(T1U, T1X);
				   T2a = VSUB(T1X, T1U);
			      }
			      {
				   V T1I, T1y, T1t, T16, T1v, TV, T1r, T1p, T2t, T2u, T2v, T2w, T1h, T1s, TU;
				   V T1o;
				   T1I = VADD(TQ, TT);
				   TU = VSUB(TQ, TT);
				   T1o = VSUB(T1m, T1n);
				   T1y = VADD(T1n, T1m);
				   T1t = VFNMS(LDK(KP414213562), T10, T15);
				   T16 = VFMA(LDK(KP414213562), T15, T10);
				   T2p = VADD(T2f, T2i);
				   STM4(&(ro[4]), T2p, ovs, &(ro[0]));
				   T2q = VSUB(T2f, T2i);
				   STM4(&(ro[12]), T2q, ovs, &(ro[0]));
				   T2r = VADD(T2j, T2k);
				   STM4(&(io[0]), T2r, ovs, &(io[0]));
				   T2s = VSUB(T2j, T2k);
				   STM4(&(io[8]), T2s, ovs, &(io[0]));
				   {
					V T28, T24, T2e, T2c;
					T28 = VSUB(T23, T1Y);
					T24 = VADD(T1Y, T23);
					T2e = VADD(T2a, T2b);
					T2c = VSUB(T2a, T2b);
					T1v = VFNMS(LDK(KP707106781), TU, TN);
					TV = VFMA(LDK(KP707106781), TU, TN);
					T1r = VFMA(LDK(KP707106781), T1o, T1l);
					T1p = VFNMS(LDK(KP707106781), T1o, T1l);
					T2t = VFNMS(LDK(KP707106781), T28, T27);
					STM4(&(io[14]), T2t, ovs, &(io[0]));
					T2u = VFMA(LDK(KP707106781), T28, T27);
					STM4(&(io[6]), T2u, ovs, &(io[0]));
					T2v = VFMA(LDK(KP707106781), T24, T1T);
					STM4(&(ro[2]), T2v, ovs, &(ro[0]));
					T2w = VFNMS(LDK(KP707106781), T24, T1T);
					STM4(&(ro[10]), T2w, ovs, &(ro[0]));
					T2x = VFNMS(LDK(KP707106781), T2e, T2d);
					STM4(&(io[10]), T2x, ovs, &(io[0]));
					T2y = VFMA(LDK(KP707106781), T2e, T2d);
					STM4(&(io[2]), T2y, ovs, &(io[0]));
					T2z = VFMA(LDK(KP707106781), T2c, T29);
					STM4(&(ro[6]), T2z, ovs, &(ro[0]));
					T2A = VFNMS(LDK(KP707106781), T2c, T29);
					STM4(&(ro[14]), T2A, ovs, &(ro[0]));
					T1h = VFNMS(LDK(KP414213562), T1g, T1b);
					T1s = VFMA(LDK(KP414213562), T1b, T1g);
				   }
				   {
					V T1z, T1J, T1K, T1G, T2B, T2C, T2D, T2E, T1C, T1F;
					T1M = VFNMS(LDK(KP414213562), T1A, T1B);
					T1C = VFMA(LDK(KP414213562), T1B, T1A);
					T1F = VFNMS(LDK(KP414213562), T1E, T1D);
					T1N = VFMA(LDK(KP414213562), T1D, T1E);
					{
					     V T1q, T1i, T1w, T1u;
					     T1q = VADD(T16, T1h);
					     T1i = VSUB(T16, T1h);
					     T1w = VADD(T1t, T1s);
					     T1u = VSUB(T1s, T1t);
					     T1L = VFNMS(LDK(KP707106781), T1y, T1x);
					     T1z = VFMA(LDK(KP707106781), T1y, T1x);
					     T1P = VFMA(LDK(KP707106781), T1I, T1H);
					     T1J = VFNMS(LDK(KP707106781), T1I, T1H);
					     T1K = VSUB(T1F, T1C);
					     T1G = VADD(T1C, T1F);
					     T2B = VFMA(LDK(KP923879532), T1q, T1p);
					     STM4(&(io[15]), T2B, ovs, &(io[1]));
					     T2C = VFNMS(LDK(KP923879532), T1q, T1p);
					     STM4(&(io[7]), T2C, ovs, &(io[1]));
					     T2D = VFMA(LDK(KP923879532), T1i, TV);
					     STM4(&(ro[3]), T2D, ovs, &(ro[1]));
					     T2E = VFNMS(LDK(KP923879532), T1i, TV);
					     STM4(&(ro[11]), T2E, ovs, &(ro[1]));
					     T2F = VFMA(LDK(KP923879532), T1w, T1v);
					     STM4(&(ro[15]), T2F, ovs, &(ro[1]));
					     T2G = VFNMS(LDK(KP923879532), T1w, T1v);
					     STM4(&(ro[7]), T2G, ovs, &(ro[1]));
					     T2H = VFMA(LDK(KP923879532), T1u, T1r);
					     STM4(&(io[3]), T2H, ovs, &(io[1]));
					     T2I = VFNMS(LDK(KP923879532), T1u, T1r);
					     STM4(&(io[11]), T2I, ovs, &(io[1]));
					}
					{
					     V T2J, T2K, T2L, T2M;
					     T2J = VFNMS(LDK(KP923879532), T1G, T1z);
					     STM4(&(ro[9]), T2J, ovs, &(ro[1]));
					     STN4(&(ro[8]), T2o, T2J, T2w, T2E, ovs);
					     T2K = VFMA(LDK(KP923879532), T1G, T1z);
					     STM4(&(ro[1]), T2K, ovs, &(ro[1]));
					     STN4(&(ro[0]), T2n, T2K, T2v, T2D, ovs);
					     T2L = VFNMS(LDK(KP923879532), T1K, T1J);
					     STM4(&(io[13]), T2L, ovs, &(io[1]));
					     STN4(&(io[12]), T2l, T2L, T2t, T2B, ovs);
					     T2M = VFMA(LDK(KP923879532), T1K, T1J);
					     STM4(&(io[5]), T2M, ovs, &(io[1]));
					     STN4(&(io[4]), T2m, T2M, T2u, T2C, ovs);
					}
				   }
			      }
			 }
		    }
	       }
	       T1O = VSUB(T1M, T1N);
	       T1Q = VADD(T1M, T1N);
	       {
		    V T2N, T2O, T2P, T2Q;
		    T2N = VFMA(LDK(KP923879532), T1Q, T1P);
		    STM4(&(io[1]), T2N, ovs, &(io[1]));
		    STN4(&(io[0]), T2r, T2N, T2y, T2H, ovs);
		    T2O = VFNMS(LDK(KP923879532), T1Q, T1P);
		    STM4(&(io[9]), T2O, ovs, &(io[1]));
		    STN4(&(io[8]), T2s, T2O, T2x, T2I, ovs);
		    T2P = VFMA(LDK(KP923879532), T1O, T1L);
		    STM4(&(ro[5]), T2P, ovs, &(ro[1]));
		    STN4(&(ro[4]), T2p, T2P, T2z, T2G, ovs);
		    T2Q = VFNMS(LDK(KP923879532), T1O, T1L);
		    STM4(&(ro[13]), T2Q, ovs, &(ro[1]));
		    STN4(&(ro[12]), T2q, T2Q, T2A, T2F, ovs);
	       }
	  }
     }
     VLEAVE();
}

static const kdft_desc desc = { 16, XSIMD_STRING("n2sv_16"), {104, 0, 40, 0}, &GENUS, 0, 1, 0, 0 };

void XSIMD(codelet_n2sv_16) (planner *p) {
     X(kdft_register) (p, n2sv_16, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_notw.native -simd -compact -variables 4 -pipeline-latency 8 -n 16 -name n2sv_16 -with-ostride 1 -include n2s.h -store-multiple 4 */

/*
 * This function contains 144 FP additions, 24 FP multiplications,
 * (or, 136 additions, 16 multiplications, 8 fused multiply/add),
 * 74 stack variables, 3 constants, and 72 memory accesses
 */
#include "n2s.h"

static void n2sv_16(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT i;
	  for (i = v; i > 0; i = i - (2 * VL), ri = ri + ((2 * VL) * ivs), ii = ii + ((2 * VL) * ivs), ro = ro + ((2 * VL) * ovs), io = io + ((2 * VL) * ovs), MAKE_VOLATILE_STRIDE(64, is), MAKE_VOLATILE_STRIDE(64, os)) {
	       V T7, T1R, T25, TC, TN, T1x, T1H, T1l, Tt, T22, T2h, T1b, T1g, T1E, T1Z;
	       V T1D, Te, T1S, T26, TJ, TQ, T1m, T1n, TT, Tm, T1X, T2g, T10, T15, T1B;
	       V T1U, T1A;
	       {
		    V T3, TL, Ty, T1k, T6, T1j, TB, TM;
		    {
			 V T1, T2, Tw, Tx;
			 T1 = LD(&(ri[0]), ivs, &(ri[0]));
			 T2 = LD(&(ri[WS(is, 8)]), ivs, &(ri[0]));
			 T3 = VADD(T1, T2);
			 TL = VSUB(T1, T2);
			 Tw = LD(&(ii[0]), ivs, &(ii[0]));
			 Tx = LD(&(ii[WS(is, 8)]), ivs, &(ii[0]));
			 Ty = VADD(Tw, Tx);
			 T1k = VSUB(Tw, Tx);
		    }
		    {
			 V T4, T5, Tz, TA;
			 T4 = LD(&(ri[WS(is, 4)]), ivs, &(ri[0]));
			 T5 = LD(&(ri[WS(is, 12)]), ivs, &(ri[0]));
			 T6 = VADD(T4, T5);
			 T1j = VSUB(T4, T5);
			 Tz = LD(&(ii[WS(is, 4)]), ivs, &(ii[0]));
			 TA = LD(&(ii[WS(is, 12)]), ivs, &(ii[0]));
			 TB = VADD(Tz, TA);
			 TM = VSUB(Tz, TA);
		    }
		    T7 = VADD(T3, T6);
		    T1R = VSUB(T3, T6);
		    T25 = VSUB(Ty, TB);
		    TC = VADD(Ty, TB);
		    TN = VSUB(TL, TM);
		    T1x = VADD(TL, TM);
		    T1H = VSUB(T1k, T1j);
		    T1l = VADD(T1j, T1k);
	       }
	       {
		    V Tp, T17, T1f, T20, Ts, T1c, T1a, T21;
		    {
			 V Tn, To, T1d, T1e;
			 Tn = LD(&(ri[WS(is, 15)]), ivs, &(ri[WS(is, 1)]));
			 To = LD(&(ri[WS(is, 7)]), ivs, &(ri[WS(is, 1)]));
			 Tp = VADD(Tn, To);
			 T17 = VSUB(Tn, To);
			 T1d = LD(&(ii[WS(is, 15)]), ivs, &(ii[WS(is, 1)]));
			 T1e = LD(&(ii[WS(is, 7)]), ivs, &(ii[WS(is, 1)]));
			 T1f = VSUB(T1d, T1e);
			 T20 = VADD(T1d, T1e);
		    }
		    {
			 V Tq, Tr, T18, T19;
			 Tq = LD(&(ri[WS(is, 3)]), ivs, &(ri[WS(is, 1)]));
			 Tr = LD(&(ri[WS(is, 11)]), ivs, &(ri[WS(is, 1)]));
			 Ts = VADD(Tq, Tr);
			 T1c = VSUB(Tq, Tr);
			 T18 = LD(&(ii[WS(is, 3)]), ivs, &(ii[WS(is, 1)]));
			 T19 = LD(&(ii[WS(is, 11)]), ivs, &(ii[WS(is, 1)]));
			 T1a = VSUB(T18, T19);
			 T21 = VADD(T18, T19);
		    }
		    Tt = VADD(Tp, Ts);
		    T22 = VSUB(T20, T21);
		    T2h = VADD(T20, T21);
		    T1b = VSUB(T17, T1a);
		    T1g = VADD(T1c, T1f);
		    T1E = VSUB(T1f, T1c);
		    T1Z = VSUB(Tp, Ts);
		    T1D = VADD(T17, T1a);
	       }
	       {
		    V Ta, TP, TF, TO, Td, TR, TI, TS;
		    {
			 V T8, T9, TD, TE;
			 T8 = LD(&(ri[WS(is, 2)]), ivs, &(ri[0]));
			 T9 = LD(&(ri[WS(is, 10)]), ivs, &(ri[0]));
			 Ta = VADD(T8, T9);
			 TP = VSUB(T8, T9);
			 TD = LD(&(ii[WS(is, 2)]), ivs, &(ii[0]));
			 TE = LD(&(ii[WS(is, 10)]), ivs, &(ii[0]));
			 TF = VADD(TD, TE);
			 TO = VSUB(TD, TE);
		    }
		    {
			 V Tb, Tc, TG, TH;
			 Tb = LD(&(ri[WS(is, 14)]), ivs, &(ri[0]));
			 Tc = LD(&(ri[WS(is, 6)]), ivs, &(ri[0]));
			 Td = VADD(Tb, Tc);
			 TR = VSUB(Tb, Tc);
			 TG = LD(&(ii[WS(is, 14)]), ivs, &(ii[0]));
			 TH = LD(&(ii[WS(is, 6)]), ivs, &(ii[0]));
			 TI = VADD(TG, TH);
			 TS = VSUB(TG, TH);
		    }
		    Te = VADD(Ta, Td);
		    T1S = VSUB(TF, TI);
		    T26 = VSUB(Td, Ta);
		    TJ = VADD(TF, TI);
		    TQ = VSUB(TO, TP);
		    T1m = VSUB(TR, TS);
		    T1n = VADD(TP, TO);
		    TT = VADD(TR, TS);
	       }
	       {
		    V Ti, T11, TZ, T1V, Tl, TW, T14, T1W;
		    {
			 V Tg, Th, TX, TY;
			 Tg = LD(&(ri[WS(is, 1)]), ivs, &(ri[WS(is, 1)]));
			 Th = LD(&(ri[WS(is, 9)]), ivs, &(ri[WS(is, 1)]));
			 Ti = VADD(Tg, Th);
			 T11 = VSUB(Tg, Th);
			 TX = LD(&(ii[WS(is, 1)]), ivs, &(ii[WS(is, 1)]));
			 TY = LD(&(ii[WS(is, 9)]), ivs, &(ii[WS(is, 1)]));
			 TZ = VSUB(TX, TY);
			 T1V = VADD(TX, TY);
		    }
		    {
			 V Tj, Tk, T12, T13;
			 Tj = LD(&(ri[WS(is, 5)]), ivs, &(ri[WS(is, 1)]));
			 Tk = LD(&(ri[WS(is, 13)]), ivs, &(ri[WS(is, 1)]));
			 Tl = VADD(Tj, Tk);
			 TW = VSUB(Tj, Tk);
			 T12 = LD(&(ii[WS(is, 5)]), ivs, &(ii[WS(is, 1)]));
			 T13 = LD(&(ii[WS(is, 13)]), ivs, &(ii[WS(is, 1)]));
			 T14 = VSUB(T12, T13);
			 T1W = VADD(T12, T13);
		    }
		    Tm = VADD(Ti, Tl);
		    T1X = VSUB(T1V, T1W);
		    T2g = VADD(T1V, T1W);
		    T10 = VADD(TW, TZ);
		    T15 = VSUB(T11, T14);
		    T1B = VADD(T11, T14);
		    T1U = VSUB(Ti, Tl);
		    T1A = VSUB(TZ, TW);
	       }
	       {
		    V T2l, T2m, T2n, T2o, T2p, T2q, T2r, T2s;
		    {
			 V Tf, Tu, T2j, T2k;
			 Tf = VADD(T7, Te);
			 Tu = VADD(Tm, Tt);
			 T2l = VSUB(Tf, Tu);
			 STM4(&(ro[8]), T2l, ovs, &(ro[0]));
			 T2m = VADD(Tf, Tu);
			 STM4(&(ro[0]), T2m, ovs, &(ro[0]));
			 T2j = VADD(TC, TJ);
			 T2k = VADD(T2g, T2h);
			 T2n = VSUB(T2j, T2k);
			 STM4(&(io[8]), T2n, ovs, &(io[0]));
			 T2o = VADD(T2j, T2k);
			 STM4(&(io[0]), T2o, ovs, &(io[0]));
		    }
		    {
			 V Tv, TK, T2f, T2i;
			 Tv = VSUB(Tt, Tm);
			 TK = VSUB(TC, TJ);
			 T2p = VADD(Tv, TK);
			 STM4(&(io[4]), T2p, ovs, &(io[0]));
			 T2q = VSUB(TK, Tv);
			 STM4(&(io[12]), T2q, ovs, &(io[0]));
			 T2f = VSUB(T7, Te);
			 T2i = VSUB(T2g, T2h);
			 T2r = VSUB(T2f, T2i);
			 STM4(&(ro[12]), T2r, ovs, &(ro[0]));
			 T2s = VADD(T2f, T2i);
			 STM4(&(ro[4]), T2s, ovs, &(ro[0]));
		    }
		    {
			 V T2t, T2u, T2v, T2w, T2x, T2y, T2z, T2A;
			 {
			      V T1T, T27, T24, T28, T1Y, T23;
			      T1T = VADD(T1R, T1S);
			      T27 = VSUB(T25, T26);
			      T1Y = VADD(T1U, T1X);
			      T23 = VSUB(T1Z, T22);
			      T24 = VMUL(LDK(KP707106781), VADD(T1Y, T23));
			      T28 = VMUL(LDK(KP707106781), VSUB(T23, T1Y));
			      T2t = VSUB(T1T, T24);
			      STM4(&(ro[10]), T2t, ovs, &(ro[0]));
			      T2u = VADD(T27, T28);
			      STM4(&(io[6]), T2u, ovs, &(io[0]));
			      T2v = VADD(T1T, T24);
			      STM4(&(ro[2]), T2v, ovs, &(ro[0]));
			      T2w = VSUB(T27, T28);
			      STM4(&(io[14]), T2w, ovs, &(io[0]));
			 }
			 {
			      V T29, T2d, T2c, T2e, T2a, T2b;
			      T29 = VSUB(T1R, T1S);
			      T2d = VADD(T26, T25);
			      T2a = VSUB(T1X, T1U);
			      T2b = VADD(T1Z, T22);
			      T2c = VMUL(LDK(KP707106781), VSUB(T2a, T2b));
			      T2e = VMUL(LDK(KP707106781), VADD(T2a, T2b));
			      T2x = VSUB(T29, T2c);
			      STM4(&(ro[14]), T2x, ovs, &(ro[0]));
			      T2y = VADD(T2d, T2e);
			      STM4(&(io[2]), T2y, ovs, &(io[0]));
			      T2z = VADD(T29, T2c);
			      STM4(&(ro[6]), T2z, ovs, &(ro[0]));
			      T2A = VSUB(T2d, T2e);
			      STM4(&(io[10]), T2A, ovs, &(io[0]));
			 }
			 {
			      V T2B, T2C, T2D, T2E, T2F, T2G, T2H, T2I;
			      {
				   V TV, T1r, T1p, T1v, T1i, T1q, T1u, T1w, TU, T1o;
				   TU = VMUL(LDK(KP707106781), VSUB(TQ, TT));
				   TV = VADD(TN, TU);
				   T1r = VSUB(TN, TU);
				   T1o = VMUL(LDK(KP707106781), VSUB(T1m, T1n));
				   T1p = VSUB(T1l, T1o);
				   T1v = VADD(T1l, T1o);
				   {
					V T16, T1h, T1s, T1t;
					T16 = VFMA(LDK(KP923879532), T10, VMUL(LDK(KP382683432), T15));
					T1h = VFNMS(LDK(KP923879532), T1g, VMUL(LDK(KP382683432), T1b));
					T1i = VADD(T16, T1h);
					T1q = VSUB(T1h, T16);
					T1s = VFNMS(LDK(KP923879532), T15, VMUL(LDK(KP382683432), T10));
					T1t = VFMA(LDK(KP382683432), T1g, VMUL(LDK(KP923879532), T1b));
					T1u = VSUB(T1s, T1t);
					T1w = VADD(T1s, T1t);
				   }
				   T2B = VSUB(TV, T1i);
				   STM4(&(ro[11]), T2B, ovs, &(ro[1]));
				   T2C = VSUB(T1v, T1w);
				   STM4(&(io[11]), T2C, ovs, &(io[1]));
				   T2D = VADD(TV, T1i);
				   STM4(&(ro[3]), T2D, ovs, &(ro[1]));
				   T2E = VADD(T1v, T1w);
				   STM4(&(io[3]), T2E, ovs, &(io[1]));
				   T2F = VSUB(T1p, T1q);
				   STM4(&(io[15]), T2F, ovs, &(io[1]));
				   T2G = VSUB(T1r, T1u);
				   STM4(&(ro[15]), T2G, ovs, &(ro[1]));
				   T2H = VADD(T1p, T1q);
				   STM4(&(io[7]), T2H, ovs, &(io[1]));
				   T2I = VADD(T1r, T1u);
				   STM4(&(ro[7]), T2I, ovs, &(ro[1]));
			      }
			      {
				   V T1z, T1L, T1J, T1P, T1G, T1K, T1O, T1Q, T1y, T1I;
				   T1y = VMUL(LDK(KP707106781), VADD(T1n, T1m));
				   T1z = VADD(T1x, T1y);
				   T1L = VSUB(T1x, T1y);
				   T1I = VMUL(LDK(KP707106781), VADD(TQ, TT));
				   T1J = VSUB(T1H, T1I);
				   T1P = VADD(T1H, T1I);
				   {
					V T1C, T1F, T1M, T1N;
					T1C = VFMA(LDK(KP382683432), T1A, VMUL(LDK(KP923879532), T1B));
					T1F = VFNMS(LDK(KP382683432), T1E, VMUL(LDK(KP923879532), T1D));
					T1G = VADD(T1C, T1F);
					T1K = VSUB(T1F, T1C);
					T1M = VFNMS(LDK(KP382683432), T1B, VMUL(LDK(KP923879532), T1A));
					T1N = VFMA(LDK(KP923879532), T1E, VMUL(LDK(KP382683432), T1D));
					T1O = VSUB(T1M, T1N);
					T1Q = VADD(T1M, T1N);
				   }
				   {
					V T2J, T2K, T2L, T2M;
					T2J = VSUB(T1z, T1G);
					STM4(&(ro[9]), T2J, ovs, &(ro[1]));
					STN4(&(ro[8]), T2l, T2J, T2t, T2B, ovs);
					T2K = VSUB(T1P, T1Q);
					STM4(&(io[9]), T2K, ovs, &(io[1]));
					STN4(&(io[8]), T2n, T2K, T2A, T2C, ovs);
					T2L = VADD(T1z, T1G);
					STM4(&(ro[1]), T2L, ovs, &(ro[1]));
					STN4(&(ro[0]), T2m, T2L, T2v, T2D, ovs);
					T2M = VADD(T1P, T1Q);
					STM4(&(io[1]), T2M, ovs, &(io[1]));
					STN4(&(io[0]), T2o, T2M, T2y, T2E, ovs);
				   }
				   {
					V T2N, T2O, T2P, T2Q;
					T2N = VSUB(T1J, T1K);
					STM4(&(io[13]), T2N, ovs, &(io[1]));
					STN4(&(io[12]), T2q, T2N, T2w, T2F, ovs);
					T2O = VSUB(T1L, T1O);
					STM4(&(ro[13]), T2O, ovs, &(ro[1]));
					STN4(&(ro[12]), T2r, T2O, T2x, T2G, ovs);
					T2P = VADD(T1J, T1K);
					STM4(&(io[5]), T2P, ovs, &(io[1]));
					STN4(&(io[4]), T2p, T2P, T2u, T2H, ovs);
					T2Q = VADD(T1L, T1O);
					STM4(&(ro[5]), T2Q, ovs, &(ro[1]));
					STN4(&(ro[4]), T2s, T2Q, T2z, T2I, ovs);
				   }
			      }
			 }
		    }
	       }
	  }
     }
     VLEAVE();
}

static const kdft_desc desc = { 16, XSIMD_STRING("n2sv_16"), {136, 16, 8, 0}, &GENUS, 0, 1, 0, 0 };

void XSIMD(codelet_n2sv_16) (planner *p) {
     X(kdft_register) (p, n2sv_16, &desc);
}

#endif				/* HAVE_FMA */
