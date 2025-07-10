#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Derivation for flux density levels in airgap due to only PM *)

Clear["Global`*"]
Circuit = Fpm/km == sigmap*2*Rgr + sigmap*Rpm; 
Fpm     = Hc * lpm;             (*PM MMF*)
Rgr     = (lgr-x)/(mu_o*Agr);   (*Airgap reluctance*)
mupm    = Br/Hc;                (*From PM demagnetization line*)
Rpm     = lpm/(mupm*Am);        (*PM reluctance*)
sigmap  = Bgr * Agr;             (*Total flux in the circuit*)

FullSimplify[Solve[Circuit,Bgr]]

(*lawcons = Bgr*Agr == Bm*Am/kl;                     (*As used in Betschon thesis Eq 5.29*)*)
(*Hgr = Bgr/mu_o; (*x = 0;*)                        (*Kimman pg34*)

demag = Bm -> (Br*Hm + Br*Hc)/Hc;                  (*Eq. of demagnetization curve from Fig3.7 Kimman thesis- Derive simply from y=mx+c*)
lawcons = lawcons /. demag;                        (*We are substituting values of Bm and Hm in load line. Here we do Bm. Below, we do Hm.*)

hmm = Flatten[Solve[FourthMax,Hm]];
lawcons = lawcons /. hmm;

eqBgr = FullSimplify[Solve[lawcons, Bgr]];
eqBgr = Bgr/.%       *)          (*This is actually Bgrm, same as in Kimman's paper*)


(*Clear["Global`*"]*)
(*Isolating the factor to be used in MATLAB*)       
eqBFEM = BFEM == eqBgr
Solve[eqBFEM,kl]
Solve[eqBFEM,km]


(* RADIAL BEARING FORCE LINEARIZATION AND DEPENDENCIES *)
Bg4Manual  = (ix mu_o nr (2 lgr-x))/(kc (2 lgr^2-x^2));                       (*Control flux density in opposite airgap to what derived above*)
BgrmManual = (Am Br Hc lm mu_o)/(2 Am Br km (lgr+x) + Agr Hc kl lm mu_o);
Fx = (eqBgr + eqBg2)^2 Agr/mu_o - (BgrmManual - Bg4Manual)^2 Agr/mu_o
Series[Fx,{ix,0,1},{x,0,1}]    (*2 variable Taylor series expansion*)

x =0;
Fx = (eqBgr + eqBg2)^2 Agr/mu_o - (BgrmManual - Bg4Manual)^2 Agr/mu_o
FullSimplify[Fx]


(* Axial bearing linearization *)
(*F = mu n^2 Ap /4 ((ib+ip)^2/(lo - z)^2 - (ib-ip)^2/(lo + z)^2 )
Series[F,{ip,0,1},{z,0,1}]*)


(* ::InheritFromParent:: *)
(**)


(*Han's calculations*)
(*Gxp = mu Agr / lgr;
Gxm = Gxp; Gyp = Gxp; Gym = Gxp;
Gzp = Gxp; Gzm = Gxp;
Gr = Gxp + Gxm + Gyp + Gym;
Ga = Gzp + Gzm;
Gm = mu Apm / lm;
G  = Gr Ga + Gm Ga + Gm Gr;
Fm = Hc lm / sigma;

Bgrm = Gm Ga Fm Gxp / G / Agr*)


(*Determining the kc using loop equations and getting Bgr values from FEM. We end up with the same kc as we got from comparison b/w Bgrc from FEM and from analytical.*)
(*Clear["Global`*"]
loop1 = Hg1*(lgr-y) - Hg2*(lgr-x) == nr*iy/kc - nr*ix/kc;
loop2 = Hg2*(lgr-x) + Hg3*(lgr+y) == nr*ix/kc + nr*iy/kc;
Hg1=Bg1/mu; Hg2=Bg2/mu; Hg3=Bg3/mu; Hg4=Bg4/mu;
x = 0; y=0;

nr = 30; ix = 3; lgr = 4*10^(-4); mu= 1.25663753*10^(-6); Bg1=0.2475; Bg3=0.2475; iy=3;
(*Manually:*) kc = 2 nr iy /  ( (Bg1 lgr + Bg3 lgr) /mu )*)


(* ::InheritFromParent:: *)
(**)


(* ::Input:: *)
(**)
