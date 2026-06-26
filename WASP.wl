(* ::Package:: *)

(* ::Section:: *)
(*Begin package*)


BeginPackage["WASP`",
  {
  "HSCSolverHomogeneousRadialTeukolsky`",
  "SpinWeightedSpheroidalHarmonics`",
  "AnalyticSpinOrbitsHamiltonJacobiSpectral`"
  }
];


(* ::Subsubsection:: *)
(*Functions for eccentric equatorial orbits*)


WASPEccEq0PA::usage = "WASPEccEq0PA[l,m,n,orbit,tol] calculates the waveform amplitudes and fluxes for a non-spinnning test-particle on eccentric, equatorial orbits.";


WASPEccEqFF::usage = "WASPEccEqFF[l,m,n,orbit,tol] calculates the waveform amplitudes and fluxes for a spinning test-particle on eccentric, equatorial orbits. The calculations are performed in the fixed frequency spin-gauge.";


WASPEccEqGenericGauge::usage = "WASPEccEqGenericGauge[l,m,n,orbit,wr0,tol] calculates the waveform amplitudes and fluxes for a spinning test-particle on eccentric, equatorial orbits. Valid for any spin-gauge.";


(* ::Subsubsection:: *)
(*Functions for circular orbits*)


WASPCircEqFF::usage = "WASPCircEqFF[l,m,a,rp] calculates the waveform amplitudes and fluxes for a spinning test-particle on circular, equatorial orbits. The calculations are performed in the fixed frequency spin-gauge.";


WASPCircEqFE::usage = "WASPCircEqFE[l,m,a,rp] calculates the waveform amplitudes and fluxes for a spinning test-particle on circular, equatorial orbits. The calculations are performed in the fixed eccentricity spin-gauge.";


WASPCircEqFT::usage = "WASPCircEqFT[l,m,a,rp] calculates the waveform amplitudes and fluxes for a spinning test-particle on circular, equatorial orbits. The calculations are performed in the fixed radius spin-gauge.";


(* ::Subsubsection:: *)
(*Begin package*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Auxiliary functions source term f_ab^(i)*)


fabi[\[Zeta]_,\[Zeta]bar_,sin\[Theta]_,\[CapitalDelta]_,d\[CapitalDelta]_,K_,dK_,S_,L2S_,L1L2S_,a_]:={-2*\[Zeta]^2/\[CapitalDelta]^2*(L1L2S-2I*a/\[Zeta]*sin\[Theta]*L2S),
                                                       4/Sqrt[2]*\[Zeta]^2/\[Zeta]bar/\[CapitalDelta]*((I*K/\[CapitalDelta]+1/\[Zeta]+1/\[Zeta]bar)*L2S-a*sin\[Theta]*K/\[CapitalDelta]*(1/\[Zeta]bar-1/\[Zeta])*S),
                                                       4/Sqrt[2]*\[Zeta]^2/\[Zeta]bar/\[CapitalDelta]*(L2S+I*a*sin\[Theta]*(1/\[Zeta]bar-1/\[Zeta])*S),
                                                       \[Zeta]^2/\[Zeta]bar^2*(I*(dK*\[CapitalDelta]-K*d\[CapitalDelta])-2/\[Zeta]*I*K*\[CapitalDelta]+K^2)/\[CapitalDelta]^2*S,
                                                       -2*\[Zeta]^2/\[Zeta]bar^2*(1/\[Zeta]+I*K/\[CapitalDelta])*S,
                                                       -\[Zeta]^2/\[Zeta]bar^2*S};
                                                       
dfabidS[\[Zeta]_,\[Zeta]bar_,sin\[Theta]_,\[CapitalDelta]_,d\[CapitalDelta]_,K_,dK_,K1_,dK1_,S_,L2S_,L1L2S_,S1_,L2S1_,L1L2S1_,a_]:={
      -2*\[Zeta]^2/\[CapitalDelta]^2*(L1L2S1-2I*a/\[Zeta]*sin\[Theta]*L2S1),
      4/Sqrt[2]*\[Zeta]^2/\[Zeta]bar/\[CapitalDelta]*((I*K/\[CapitalDelta]+1/\[Zeta]+1/\[Zeta]bar)*L2S1+(I*K1/\[CapitalDelta])*L2S-a*sin\[Theta]*K/\[CapitalDelta]*(1/\[Zeta]bar-1/\[Zeta])*S1-a*sin\[Theta]*K1/\[CapitalDelta]*(1/\[Zeta]bar-1/\[Zeta])*S),
      4/Sqrt[2]*\[Zeta]^2/\[Zeta]bar/\[CapitalDelta]*(L2S1+I*a*sin\[Theta]*(1/\[Zeta]bar-1/\[Zeta])*S1),
      \[Zeta]^2/\[Zeta]bar^2*((I*(dK*\[CapitalDelta]-K*d\[CapitalDelta])-2/\[Zeta]*I*K*\[CapitalDelta]+K^2)/\[CapitalDelta]^2*S1+(I*(dK1*\[CapitalDelta]-K1*d\[CapitalDelta])-2/\[Zeta]*I*K1*\[CapitalDelta]+2*K*K1)/\[CapitalDelta]^2*S),
      -2*\[Zeta]^2/\[Zeta]bar^2*((1/\[Zeta]+I*K/\[CapitalDelta])*S1+(I*K1/\[CapitalDelta])*S),
      -\[Zeta]^2/\[Zeta]bar^2*S1};

dfabidr[\[Zeta]_,\[Zeta]bar_,sin\[Theta]_,\[CapitalDelta]_,d\[CapitalDelta]_,K_,dK_,S_,L2S_,L1L2S_,a_,\[Omega]_]:={
      1/\[CapitalDelta]^3 (\[CapitalDelta] (4 I a L2S sin\[Theta]-4 L1L2S \[Zeta])+4 \[Zeta] (-2 I a L2S sin\[Theta]+L1L2S \[Zeta]) d\[CapitalDelta]),
      1/(\[Zeta]bar^3 \[CapitalDelta]^3) 2 Sqrt[2] (a S (\[Zeta]-\[Zeta]bar) sin\[Theta] (2 d\[CapitalDelta] \[Zeta] \[Zeta]bar K-(dK \[Zeta] \[Zeta]bar+(-2 \[Zeta]+\[Zeta]bar) K) \[CapitalDelta])+L2S (-I \[Zeta] \[Zeta]bar K (2 d\[CapitalDelta] \[Zeta] \[Zeta]bar+(\[Zeta]-2 \[Zeta]bar) \[CapitalDelta])+\[CapitalDelta] (\[Zeta] \[Zeta]bar (I dK \[Zeta] \[Zeta]bar-d\[CapitalDelta] (\[Zeta]+\[Zeta]bar))-(\[Zeta]-\[Zeta]bar) (2 \[Zeta]+\[Zeta]bar) \[CapitalDelta]))),
      1/(\[Zeta]bar^3 \[CapitalDelta]^2) 2 Sqrt[2] (-L2S \[Zeta] \[Zeta]bar (d\[CapitalDelta] \[Zeta] \[Zeta]bar+(\[Zeta]-2 \[Zeta]bar) \[CapitalDelta])-I a S (\[Zeta]-\[Zeta]bar) sin\[Theta] (d\[CapitalDelta] \[Zeta] \[Zeta]bar+(2 \[Zeta]-\[Zeta]bar) \[CapitalDelta])),
      -(1/(\[Zeta]bar^3 \[CapitalDelta]^3))2 I S (-I \[Zeta] K^2 (d\[CapitalDelta] \[Zeta] \[Zeta]bar+(\[Zeta]-\[Zeta]bar) \[CapitalDelta])+\[Zeta]^2 \[CapitalDelta] (dK d\[CapitalDelta] \[Zeta]bar+(dK-\[Zeta]bar \[Omega]) \[CapitalDelta])+K (-d\[CapitalDelta]^2 \[Zeta]^2 \[Zeta]bar+\[Zeta]^2 (-d\[CapitalDelta]+\[Zeta]bar+I dK \[Zeta]bar) \[CapitalDelta]+(-2 \[Zeta]+\[Zeta]bar) \[CapitalDelta]^2)),
      1/(\[Zeta]bar^3 \[CapitalDelta]^2) 2 S (I \[Zeta] K (d\[CapitalDelta] \[Zeta] \[Zeta]bar+2 (\[Zeta]-\[Zeta]bar) \[CapitalDelta])+\[CapitalDelta] (-I dK \[Zeta]^2 \[Zeta]bar+(2 \[Zeta]-\[Zeta]bar) \[CapitalDelta])),
      (2 S \[Zeta] (\[Zeta]-\[Zeta]bar))/\[Zeta]bar^3};

dfabid\[Theta][\[Zeta]_,\[Zeta]bar_,sin\[Theta]_,\[CapitalDelta]_,d\[CapitalDelta]_,K_,dK_,S_,L2S_,L1L2S_,dSd\[Theta]_,dL2Sd\[Theta]_,dL1L2Sd\[Theta]_,a_]:={
      (-4 a^2 L2S sin\[Theta]^2+4 I a (dL2Sd\[Theta]-L1L2S) sin\[Theta] \[Zeta]+2 \[Zeta] (-((dL1L2Sd\[Theta]+L2S) \[Zeta])+L2S \[Zeta]bar))/\[CapitalDelta]^2,
      Sqrt[2] (-2 I a^2 K S sin\[Theta]^2 (2 \[Zeta]-\[Zeta]bar) (\[Zeta]+\[Zeta]bar)+\[Zeta] \[Zeta]bar (-I K S (\[Zeta]-\[Zeta]bar)^2+2 I dL2Sd\[Theta] K \[Zeta] \[Zeta]bar+2 dL2Sd\[Theta] \[CapitalDelta] (\[Zeta]+\[Zeta]bar))+2 a sin\[Theta] (dSd\[Theta] K \[Zeta] \[Zeta]bar (-\[Zeta]+\[Zeta]bar)+I L2S \[CapitalDelta] (\[Zeta]+\[Zeta]bar) (2 \[Zeta]+\[Zeta]bar)-K L2S \[Zeta] \[Zeta]bar (\[Zeta]+2 \[Zeta]bar)))/(\[CapitalDelta]^2 \[Zeta]bar^3),
      Sqrt[2] (-2 a^2 S sin\[Theta]^2 (2 \[Zeta]-\[Zeta]bar) (\[Zeta]+\[Zeta]bar)-\[Zeta] \[Zeta]bar (S (\[Zeta]-\[Zeta]bar)^2-2 dL2Sd\[Theta] \[Zeta] \[Zeta]bar)+2 I a sin\[Theta] \[Zeta] \[Zeta]bar (dSd\[Theta] (\[Zeta]-\[Zeta]bar)+L2S (\[Zeta]+2 \[Zeta]bar)))/(\[CapitalDelta] \[Zeta]bar^3),
      1/(\[Zeta]bar^3 \[CapitalDelta]^2) (\[Zeta] (dSd\[Theta] \[Zeta] \[Zeta]bar+2 I a S sin\[Theta] (\[Zeta]+\[Zeta]bar)) K^2+dK \[Zeta] (I dSd\[Theta] \[Zeta] \[Zeta]bar-2 a S sin\[Theta] (\[Zeta]+\[Zeta]bar)) \[CapitalDelta]+K (d\[CapitalDelta] \[Zeta] (-I dSd\[Theta] \[Zeta] \[Zeta]bar+2 a S sin\[Theta] (\[Zeta]+\[Zeta]bar))+2 (-I dSd\[Theta] \[Zeta] \[Zeta]bar+a S sin\[Theta] (2 \[Zeta]+\[Zeta]bar)) \[CapitalDelta])),
      -(1/(\[Zeta]bar^3 \[CapitalDelta]))2 (\[Zeta] (I dSd\[Theta] \[Zeta] \[Zeta]bar-2 a S sin\[Theta] (\[Zeta]+\[Zeta]bar)) K+(dSd\[Theta] \[Zeta] \[Zeta]bar+I a S sin\[Theta] (2 \[Zeta]+\[Zeta]bar)) \[CapitalDelta]),
      -((\[Zeta] (dSd\[Theta] \[Zeta] \[Zeta]bar+2 I a S sin\[Theta] (\[Zeta]+\[Zeta]bar)))/\[Zeta]bar^3)};



(* ::Subsection::Closed:: *)
(*Computation phase in the source term*)


(* ::Subsubsection::Closed:: *)
(*Reference*)


\[CapitalDelta]\[Psi]funref[{\[Omega]_,m_},{EEg_,Lzg_},a_,p_,e_,x_]:=Module[{rp,rm,r1g,r2g,r3g,krg,\[Gamma]r,ellK,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]r=1-r1g/r2g;
	ellK=2EllipticK[krg];
	ellE=2EllipticE[krg];
	ellPi=2EllipticPi[\[Gamma]r,krg];

	\[ScriptCapitalI]=(2ellK)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]rover=2/Sqrt[r2g(r1g-r3g)]r1g*ellPi;
	\[ScriptCapitalI]r2over=-Sqrt[r2g(r1g-r3g)](-ellE+r1g*ellK/(r1g-r3g)-r1g(r1g+r2g+r3g)ellPi/(r2g(r1g-r3g)));
	\[ScriptCapitalI]p=2/Sqrt[(r1g-r3g)r2g] 1/rp (r1g/(r1g-rp) 2EllipticPi[(rp(r1g-r2g))/((r1g-rp)r2g),krg]- ellK);
	\[ScriptCapitalI]mreg=2/Sqrt[(r1g-r3g)r2g] (r1g/(r1g-rm) 2EllipticPi[(rm(r1g-r2g))/((r1g-rm)r2g),krg]- ellK);

	1/(2\[Pi])1/Sqrt[1-EEg^2](-EEg*(\[ScriptCapitalI]r2over+2\[ScriptCapitalI]rover)\[Omega]+\[ScriptCapitalI](Lzg*m-4*EEg*\[Omega])+(rp*\[ScriptCapitalI]p(Lzg*m*rm-4EEg(-2+rm)\[Omega]-2a(EEg*m+Lzg*\[Omega])))/(rm-rp)+(\[ScriptCapitalI]mreg(-Lzg*m*rp+4EEg(-2+rp)\[Omega]+2a(EEg*m+Lzg*\[Omega])))/(rm-rp))
]


(* ::Subsubsection::Closed:: *)
(*Periastron as initial condition*)


\[CapitalDelta]\[Psi]ICr2gfun[wr_,{\[Omega]_,m_},{EEg_,Lzg_},\[Psi]ref_,a_,p_,e_,x_]:=Module[{rp,rm,r1g,r2g,r3g,krg,krghold,\[Gamma]r,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]p,\[ScriptCapitalI]m},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]r=(r1g-r2g)/(r1g-r3g);
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result. Having \[Phi]=\[Pi]/2 with MachinePrecision cause some problems with the integrals \[ScriptCapitalI]p and \[ScriptCapitalI]m*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]rover=2/Sqrt[r2g(r1g-r3g)]((r2g-r3g)ellPi+r3g ellF);
	\[ScriptCapitalI]r2over=2r2g^2/Sqrt[r2g(r1g-r3g)]((r3g/r2g)^2*ellF+2 r3g/r2g(1-r3g/r2g)ellPi+(1-r3g/r2g)^2 1/(2(\[Gamma]r-1)(krg-\[Gamma]r)) (\[Gamma]r*ellE+(krg-\[Gamma]r)ellF+(2\[Gamma]r*krg+2\[Gamma]r-\[Gamma]r^2-3krg)ellPi-(\[Gamma]r^2Sin[\[Phi]]Cos[\[Phi]]Sqrt[1-krg Sin[\[Phi]]^2])/(1-\[Gamma]r*Sin[\[Phi]]^2)));

	\[ScriptCapitalI]p=2/(Sqrt[r2g(r1g-r3g)](r3g-rp))(ellF-(r2g-r3g)/(r2g-rp)EllipticPi[(r3g-rp)(r1g-r2g)/((r1g-r3g)(r2g-rp)),\[Phi],krg]);
	\[ScriptCapitalI]m=2/(Sqrt[r2g(r1g-r3g)](r3g-rm))(ellF-(r2g-r3g)/(r2g-rm)EllipticPi[(r3g-rm)(r1g-r2g)/((r1g-r3g)(r2g-rm)),\[Phi],krg]);

	Remove[Evaluate[ToExpression[Pick[Names["WASP`Private`*"],StringMatchQ[#,"WASP`Private`krghold$"~~__]&/@Names["WASP`Private`*"],True]]]];

	1/Sqrt[1-EEg^2](EEg(\[ScriptCapitalI]r2over+2*\[ScriptCapitalI]rover)\[Omega]+\[ScriptCapitalI](4EEg*\[Omega]-Lzg*m)+1/(rm-rp)\[ScriptCapitalI]m(8EEg*rm*\[Omega]+a^2(Lzg*m-4EEg*\[Omega])-2a*rm(EEg*m+Lzg*\[Omega]))+1/(rm-rp)\[ScriptCapitalI]p(-8EEg*rp*\[Omega]+a^2(-Lzg*m+4EEg*\[Omega])+2a*rp(EEg*m+Lzg*\[Omega])))+\[Psi]ref*wr
]


(* ::Subsubsection::Closed:: *)
(*Apoastron as initial condition*)


\[CapitalDelta]\[Psi]ICr1gfun[wr_,{\[Omega]_,m_},{EEg_,Lzg_},\[Psi]ref_,a_,p_,e_,x_]:=Module[{rp,rm,r1g,r2g,r3g,krg,krghold,\[Gamma]r,\[Phi],ellF,ellE,ellPi,\[ScriptCapitalI],\[ScriptCapitalI]rover,\[ScriptCapitalI]r2over,\[ScriptCapitalI]p,\[ScriptCapitalI]mreg},
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=2/(1-EEg^2)-(r1g+r2g);

	krg=((r1g-r2g)r3g)/((r1g-r3g)r2g);
	\[Gamma]r=1-r1g/r2g;
	(*The following "trick" ensures that, for wr=\[Pi], \[Phi] is exaclty \[Pi]/2 and not an approximate result. Having \[Phi]=\[Pi]/2 with MachinePrecision cause some problems with the integrals \[ScriptCapitalI]p and \[ScriptCapitalI]m*)
	\[Phi]=JacobiAmplitude[1/\[Pi] EllipticK[krghold] wr,krghold]/.krghold->krg;
	ellF=EllipticF[\[Phi],krg];
	ellE=EllipticE[\[Phi],krg];
	ellPi=EllipticPi[\[Gamma]r,\[Phi],krg];

	\[ScriptCapitalI]=(2ellF)/Sqrt[r2g(r1g-r3g)];
	\[ScriptCapitalI]rover=2/Sqrt[r2g(r1g-r3g)]r1g*ellPi;
	\[ScriptCapitalI]r2over=-Sqrt[r2g(r1g-r3g)](-ellE+r1g*ellF/(r1g-r3g)-r1g(r1g+r2g+r3g)ellPi/(r2g(r1g-r3g))+((r1g-r2g)Cos[\[Phi]]*Sin[\[Phi]]Sqrt[1-krg*Sin[\[Phi]]^2])/(-r2g+(-r1g+r2g)Sin[\[Phi]]^2));
	\[ScriptCapitalI]p=2/Sqrt[(r1g-r3g)r2g] 1/rp (r1g/(r1g-rp) EllipticPi[(rp(r1g-r2g))/((r1g-rp)r2g),\[Phi],krg]- ellF);
	\[ScriptCapitalI]mreg=2/Sqrt[(r1g-r3g)r2g](r1g/(r1g-rm) EllipticPi[(rm(r1g-r2g))/((r1g-rm)r2g),\[Phi],krg]- ellF);
	
	Remove[Evaluate[ToExpression[Pick[Names["WASP`Private`*"],StringMatchQ[#,"WASP`Private`krghold$"~~__]&/@Names["WASP`Private`*"],True]]]];

	1/Sqrt[1-EEg^2](EEg(\[ScriptCapitalI]r2over+2*\[ScriptCapitalI]rover)\[Omega]+\[ScriptCapitalI](4EEg*\[Omega]-Lzg*m)+1/(rm-rp)\[ScriptCapitalI]mreg(8EEg*\[Omega]+rp(Lzg*m-4EEg*\[Omega])-2a(EEg*m+Lzg*\[Omega]))+1/(rm-rp)\[ScriptCapitalI]p(-8EEg*rp*\[Omega]+a^2(-Lzg*m+4EEg*\[Omega])+2a*rp(EEg*m+Lzg*\[Omega])))+\[Psi]ref*wr
]


(* ::Section::Closed:: *)
(*Functions to compute linearized angular Teukolsky equations*)


(* ::Subsection::Closed:: *)
(*SWSH eigenvalues and eigenfunctions*)


(* ::Text:: *)
(*No linear in spin corrections*)


angpar0[s_,l_,m_,a_,\[Omega]0_]:=Module[{\[Lambda]0,\[Gamma]0,swshY0,swsh0all,swshcoeff,swshS},
	\[Gamma]0 = a*\[Omega]0;
	If[a==0,
	 \[Lambda]0=l(1+l)-s(1+s);
	 swshY0 = Evaluate[SpinWeightedSphericalHarmonicY[s,l,m,#,0]]&;
	 swsh0all = Function[{\[Theta]},{swshY0[\[Theta]],swshY0'[\[Theta]]}];
	 ,
	 \[Lambda]0 = SpinWeightedSpheroidalEigenvalue[s,l,m,\[Gamma]0]//Quiet;
	 swshcoeff = SpinWeightedSpheroidalHarmonicS[s,l,m,\[Gamma]0,Method->{"SphericalExpansion"}]//Quiet;
	 (*swshcoeff = SpinWeightedSpheroidalHarmonicS[s,l,m,\[Gamma]0]//Quiet;*)
	 swshS = swshcoeff[#,0]&;
	 swsh0all = Function[{\[Theta]},{swshS[\[Theta]]//Quiet,swshS'[\[Theta]]//Quiet}];
	];
	 {\[Lambda]0,swsh0all}
]


(* ::Text:: *)
(*With linear in spin corrections*)


angpar[s_,l_,m_,a_,\[Omega]0_,\[Omega]1_]:=Module[{\[Lambda]0,\[Lambda]1,\[Gamma]0,\[Gamma]1,swshY0,swsh0all,swsh1all},
	\[Gamma]0 = a*\[Omega]0;
	\[Gamma]1 = a*\[Omega]1;
	If[a==0,
		\[Lambda]0=l(1+l)-s(1+s);
		swshY0 = Evaluate[SpinWeightedSphericalHarmonicY[s,l,m,#,0]]&;
		swsh0all = Function[{\[Theta]},{swshY0[\[Theta]],swshY0'[\[Theta]]}];
		{\[Lambda]1,swsh1all} = {0,{0,0,0}&};
		,
		{\[Lambda]0,\[Lambda]1,swsh0all,swsh1all}= angparKerr[s,l,m,\[Gamma]0,\[Gamma]1]//Quiet;
	];
	{\[Lambda]0,\[Lambda]1,swsh0all,swsh1all}
]


angparKerr[s_,l_,m_,\[Gamma]0_,\[Gamma]1_]:=Module[{swshcoeff,\[Lambda]0,\[Lambda]1,swsh0,swsh1},
	\[Lambda]0 = SpinWeightedSpheroidalEigenvalue[s,l,m,\[Gamma]0]//Quiet;
	swshcoeff = SpinWeightedSpheroidalHarmonicS[s,l,m,\[Gamma]0,Method->{"SphericalExpansion"}]//Quiet;
	swsh0 = swshcoeff[#,0]&;
	\[Lambda]1=angeigenvalue1corrint[s,l,m,\[Gamma]0,\[Gamma]1];
	swsh1=angeigenfunction1corrint[s,l,m,\[Gamma]0,\[Gamma]1];

	(*Remove[swshcoeff];*)
	{\[Lambda]0,\[Lambda]1,Function[{\[Theta]},Evaluate[{swsh0[\[Theta]]//Quiet,swsh0'[\[Theta]]//Quiet}]], Function[{\[Theta]},Evaluate[{swsh1[\[Theta]],swsh1'[\[Theta]]}]]}
]


(* ::Subsubsection::Closed:: *)
(*Kernel linear corrections SWSH and their eigenvalues - DEPRECATED*)


angpar1corr[s_,l_,m_,\[Gamma]0_,\[Gamma]1_,\[Lambda]0_,coeff_]:=Module[{km, kp, ks,an1, nmin,nmax,\[Alpha]n,an0Tab,an1Tab,an0square,an1square,an1an0prod,\[Beta]n0,\[Beta]n1,\[Gamma]n0,\[Gamma]n1,\[Lambda]1,swsh1,normterm0,norm0,\[ScriptCapitalA],\[ScriptCapitalB],coeff1order,sign,sum0,sum1,nprec1,\[Theta]},
	km = Abs[m-s]/2;
	kp = Abs[m+s]/2;
	ks = 2(kp+km);
	
	an0Tab = coeff[[1]];
	nmin = coeff[[2]];
	nmax = coeff[[3]];
	
	an0square = Table[Sum[an0Tab[[j]]an0Tab[[i+2-nmin-j]], {j, 1, i-nmin+1}],{i,nmin,nmax}];
	normterm0[i_]:= 2^i Pochhammer[2+ks+i,-2kp-1] Hypergeometric1F1[1+2 km+i,2+ks+i,4\[Gamma]0]an0square[[i+1]];
	sum0 = Total[Table[normterm0[i],{i,nmin,nmax}],Method->"CompensatedSummation"];
	
	(* Normalisation *)
	norm0=Sqrt[2\[Pi]](2^(1+ks)E^(-2 \[Gamma]0) Gamma[1+2 kp]sum0)^(1/2);
	
	(* Overall sign such that the \[Gamma]\[Rule]0 limit is continuous *)
	sign = If[(OddQ[l]&&EvenQ[m+s])||(OddQ[l]&&OddQ[m+s]&&m>=s)||(EvenQ[l]&&OddQ[m-s]&&m<=s), -1, 1];
	
	coeff1order[n_]:= 2\[Pi](Gamma[1+2 km+n]/ Gamma[3+ks+n])2^(2+ks+n)E^(-2\[Gamma]0)Gamma[1+2kp]((1+2kp)(2+ks+n+2s)Hypergeometric1F1[1+2 km+n,3+ks+n,4 \[Gamma]0]-(2+ks+n)(1+2kp-m+s)Hypergeometric1F1[1+2 km+n,2+ks+n,4\[Gamma]0]);
	
	(* First order correction to the eigenvalues: \[Lambda]1 = \[Gamma]1(Alm\[Gamma]1 + 2(\[Gamma]0 - m)) *)
	\[Lambda]1 =-(1/norm0)^2*Total[Table[coeff1order[i]an0square[[i+1]], {i,nmin, nmax}],Method->"CompensatedSummation"];
	
	(* The interpolation method introduces an exta \[Gamma]1 that is not needed in the computation of the eigenfunctions. This factor is reintroduced at the end *)
	If[Precision[\[Lambda]1]<Min[Precision[\[Gamma]1],15],
		\[Lambda]1 = (1/\[Gamma]1)angeigenvalue1corrint[s,l,m,\[Gamma]0,\[Gamma]1];
	];
	
	(* Coefficients of first order expansion in \[Sigma] \[Gamma]1*);
	\[Alpha]n[n_]:= -2(n+1)(n+2 km+1);
	\[Beta]n0[n_]:= n(n-1)+2n(km+kp+1-2\[Gamma]0)-(2\[Gamma]0(2km+s+1)-(km+kp)(km+kp+1))-(s(s+1)+\[Lambda]0+2m*\[Gamma]0);
	\[Beta]n1[n_]:= -2(1+2 km+m+2 n+s)-\[Lambda]1;
	\[Gamma]n0[n_]:= 2\[Gamma]0(n+kp+km+s);
	\[Gamma]n1[n_]:= (ks+2 (n+s));
	
	an1[nmin] = 0;
	an1[nmin+1] = -(\[Beta]n1[nmin]/\[Alpha]n[nmin])an0Tab[[nmin+1]];
	an1[n_]:= an1[n] =(-an1[n-1] \[Beta]n0[n-1]-an0Tab[[n]] \[Beta]n1[n-1]-an1[n-2] \[Gamma]n0[n-1]-an0Tab[[n-1]] \[Gamma]n1[n-1])/\[Alpha]n[n-1];     
	an1Tab = Table[an1[i],{i,nmin,nmax}];
	
	(*This algorithm is taken from "Convergence acceleration for the numerical solution of second-order linear recurrence relations"
	   by P. Levrie and A. Bultheel *)
	Block[{rs1,rs2,fs},
	 (*It is possible to reuse the rs[n] of the zeroth order correction *)
	  rs1[n_]:= rs1[n]= -(\[Gamma]n0[n+1]/\[Alpha]n[n+1])/((\[Beta]n0[n+1]/\[Alpha]n[n+1])+rs1[n+1]);
	  rs1[nmax] = 0;
	  Table[rs1[j],{j,nmax,1,-1}];
	  rs2[n_]:= rs2[n]= -((an0Tab[[n+2]] \[Beta]n1[n+1]+an0Tab[[n+1]] \[Gamma]n1[n+1])/\[Alpha]n[n+1]+rs2[n+1])/((\[Beta]n0[n+1]/\[Alpha]n[n+1])+rs1[n+1]);
	  rs2[nmax] = 0;
	  Table[rs2[j],{j,nmax,1,-1}];
	  (*Number of high precision elements obtained by solving the three-terms recurrence relation directly. These terms are the seeds for the backward algorithm*)
	  fs[n_] := fs[n] = rs1[n-1]fs[n-1]+rs2[n-1];
	  nprec1 = Length[Position[an1Tab,_?(Precision[#]>(Precision[\[Gamma]0]-5)&)]];
	  Table[fs[i] = an1[i],{i,nmin,nprec1}];
	  an1Tab = Table[fs[j],{j,nmin,nmax}];
	];
	an1an0prod = Table[Sum[an1Tab[[j]]an0Tab[[i+2-nmin-j]], {j,1,i-nmin+1}],{i,nmin,nmax}];
	
	\[ScriptCapitalA][i_]:= 2^i Pochhammer[2+i+ks,-1-2kp](4(1+i+2km)/(2+i+ks))Hypergeometric1F1[2+i+2 km,3+i+ks,4\[Gamma]0];
	\[ScriptCapitalB][i_]:= 2^i Pochhammer[2+i+ks,-1-2kp]Hypergeometric1F1[1+i+2km,2+i+ks,4\[Gamma]0];
	
	sum1 = Total[Table[\[ScriptCapitalB][i]2an1an0prod[[i+1]]+\[ScriptCapitalA][i] an0square[[i+1]],{i,nmin,nmax}],Method->"CompensatedSummation"];
	
	swsh1=Function[{\[Theta]},Evaluate[\[Gamma]1(sign/norm0)E^(\[Gamma]0  Cos[\[Theta]])If[km==0,1,(2Cos[\[Theta]/2]^2)^km] If[kp==0,1,(2Sin[\[Theta]/2]^2)^kp] ((-sum1/(2sum0)+2Cos[\[Theta]/2]^2)an0Tab+an1Tab) . (2Cos[\[Theta]/2]^2)^Range[0,nmax]]];

	Remove[\[Alpha]n,\[Beta]n0,\[Beta]n1,\[Gamma]n0,\[Gamma]n1,sum0,sum1,norm0,sign,an0Tab,an1Tab,kp,km];
	{\[Gamma]1*\[Lambda]1,swsh1}
]


(* ::Subsubsection::Closed:: *)
(*Corrections to eigenvalues and eigenfunctions through numerical derivative*)


angeigenvalue1corrint[s_,l_,m_,\[Gamma]0_,\[Gamma]1_]:=Module[{qratio,\[Lambda]p,\[Lambda]m},
	qratio = 1 10^(-5);
	\[Lambda]p = SpinWeightedSpheroidalEigenvalue[s,l,m,\[Gamma]0+qratio*\[Gamma]1];
	\[Lambda]m = SpinWeightedSpheroidalEigenvalue[s,l,m,\[Gamma]0-qratio*\[Gamma]1];
	(\[Lambda]p-\[Lambda]m)/(2qratio)
]


angeigenfunction1corrint[s_,l_,m_,\[Gamma]0_,\[Gamma]1_]:=Module[{qratio,swshcoeffp,swshcoeffm,swshSp,swshSm},
	qratio = 1 10^(-5);
	swshcoeffp=SpinWeightedSpheroidalHarmonicS[s,l,m,\[Gamma]0+qratio*\[Gamma]1,Method->{"SphericalExpansion"}]//Quiet;
	swshcoeffm=SpinWeightedSpheroidalHarmonicS[s,l,m,\[Gamma]0-qratio*\[Gamma]1,Method->{"SphericalExpansion"}]//Quiet;
	swshSp=swshcoeffp[#,0]&;
	swshSm=swshcoeffm[#,0]&;

	Function[{\[Theta]},Evaluate[(swshSp[\[Theta]]-swshSm[\[Theta]])/(2qratio)//Quiet]]
]


(* ::Section::Closed:: *)
(*Source term for a spinning particle - eccentric, equatorial orbits*)


(* ::Subsection::Closed:: *)
(*Geodesic*)


(* ::Text:: *)
(*The sign of the amplitudes now follow the same convention of the BHPToolkit (the sign convention was changed after the first version of this code).*)


SourceTermEccEq0PA[wr_,l_,m_,n_,orbit_,{\[Lambda]_,SWSH_},{RInsol0_,dRInsol0_},{RUpsol0_,dRUpsol0_}]:=
Module[{wr0=\[Pi],sign,a,p,e,x,En,Lz,Kc,\[CapitalOmega]r,\[CapitalOmega]\[Phi],rg,\[CapitalGamma],\[Omega],W,sumPlus0,sumMin0,rp,Ur,\[Psi]ref,phase,expr,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,V,dV,RIn0,dRIn0,ddRIn0,RUp0,dRUp0,ddRUp0,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,\[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,vl,vn,vm,vmb,Cmnn,Cmnmb,Cmmbmb,CPlus0,CMin0,A2g,A1g,A0g},
    
     {a,p,e,x}=orbit["OrbitalElements"];
    (* Geodesic constants of motion *)
    En = orbit["Eg"]; 
    Lz = orbit["Lzg"];
    Kc = (Lz-a*En)^2;
    (* Geodesic BL frequencies *)
    \[CapitalOmega]r = orbit["BLFrequenciesGeo"][[1]];
    \[CapitalOmega]\[Phi] = orbit["BLFrequenciesGeo"][[2]];
    (* Linear corrections to the frequencies *)
    \[CapitalGamma] = orbit["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
    \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r;
       
    rp = orbit["rg"][wr+wr0];
    Ur = sign*Sqrt[((rp^2+a^2)*En-a*Lz)^2-(rp^2-2rp+a^2)*(rp^2+Kc)];
   (* phase = \[Omega]*orbit["\[CapitalDelta]trg"][wr+wr0]-m*orbit["\[CapitalDelta]\[Phi]rg"][wr+wr0];*)
    \[Psi]ref=\[CapitalDelta]\[Psi]funref[{\[Omega],m},{En,Lz},a,p,e,x];
    phase=\[CapitalDelta]\[Psi]ICr1gfun[wr+wr0,{\[Omega],m},{En,Lz},\[Psi]ref,a,p,e,x];
    
    expr = Exp[sign*I*(phase+n*wr)];(* Geodesic exponential term with \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    \[CapitalDelta] = rp^2-2rp+a^2;
    KTeu = (rp^2+a^2)*\[Omega]-a*m;
    d\[CapitalDelta] = 2*(rp-1);
    dKTeu = 2*rp*\[Omega];
    V = -(KTeu^2+4I*(rp-1)*KTeu)/\[CapitalDelta]+8*I*\[Omega]*rp+\[Lambda]; (* Potential in radial Teukolsky equation *)
    dV = -((2KTeu*dKTeu+4I*KTeu+4I*(rp-1)*dKTeu)*\[CapitalDelta]-(KTeu^2+4I*(rp-1)*KTeu)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
    RIn0 = RInsol0[rp]; (* radial function *)
    dRIn0 = dRInsol0[rp];
    ddRIn0 = (V*RIn0+d\[CapitalDelta]*dRIn0)/\[CapitalDelta]; (* second derivative of radial function from Teukolsky equation *)
    RUp0 = RUpsol0[rp];
    dRUp0 = dRUpsol0[rp];
    ddRUp0 = (V*RUp0+d\[CapitalDelta]*dRUp0)/\[CapitalDelta];  
    (* Spin-weighted spheroidal harmonics S(\[Theta](z)) *)
    S = SWSH[\[Pi]/2][[1]]; 
    dSd\[Theta] = SWSH[\[Pi]/2][[2]]; (* First derivative of S wrt \[Theta] *)
    d2Sd\[Theta]2 = S((a*\[Omega])^2-2(a*\[Omega])*m+m^2+2-\[Lambda]); (* second derivative of S from Teukolsky equation *)
    L2S = dSd\[Theta]-S(m-a*\[Omega]); (* Operators acting on S(\[Theta]) *)
    L1L2S = d2Sd\[Theta]2+(dSd\[Theta](-2m+2a*\[Omega]))+S(-2+m^2-2a*m*\[Omega]+a^2*\[Omega]^2);   
    
    \[Zeta] = rp;
    \[Zeta]bar = rp;
    \[CapitalSigma] = rp^2;
    {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[rp,rp,1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,S,L2S,L1L2S,a];
    (* Four-velocity in Kinnersley tetrad *)
    vn = -((rp^2+a^2)*En - a*Lz + Ur)/(2*\[CapitalSigma]); 
    vl = -((rp^2+a^2)*En - a*Lz - Ur)/(\[CapitalDelta]);
    vm = (I*(a*En - Lz))/(-Sqrt[2]*\[Zeta]bar);
    vmb = (I*(a*En - Lz))/(Sqrt[2]*\[Zeta]bar);
    Cmnn = vn^2;
    Cmnmb = vn*vmb;
    Cmmbmb = vmb^2;
    
    A2g = Cmmbmb*fmbmb2;
    A1g = (Cmnmb*fnmb1+Cmmbmb*fmbmb1);
    A0g = (Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0);
   
    sumPlus0 = (\[CapitalSigma]*(A0g*RIn0 - A1g*dRIn0 + A2g*ddRIn0)*expr/.sign->1) + (\[CapitalSigma]*(A0g*RIn0 - A1g*dRIn0 + A2g*ddRIn0)*expr/.sign->-1);
    sumMin0 = (\[CapitalSigma]*(A0g*RUp0 - A1g*dRUp0 + A2g*ddRUp0)*expr/.sign->1) + (\[CapitalSigma]*(A0g*RUp0 - A1g*dRUp0 + A2g*ddRUp0)*expr/.sign->-1);  
    
    W = (RInsol0[p/(1-e)]*dRUpsol0[p/(1-e)] - dRInsol0[p/(1-e)]*RUpsol0[p/(1-e)])/((p/(1-e))^2-2(p/(1-e))+a^2); (* Invariant Wronskian *)
    
    CPlus0 = 2*Pi*sumPlus0/(\[CapitalGamma]*W); (* Geodesic amplitudes *)
    CMin0 = 2*Pi*sumMin0/(\[CapitalGamma]*W);
    
    -{CPlus0,CMin0}
]


(* ::Subsection::Closed:: *)
(*Fixed frequency*)


(* ::Text:: *)
(*The sign of the amplitudes now follow the same convention of the BHPToolkit (the sign convention was changed after the first version of this code).*)


SourceTermEccEqFF[wr_,l_,m_,n_,orbit_,{\[Lambda]_,SWSH_},{RInsol0_,dRInsol0_},{RUpsol0_,dRUpsol0_}]:=
Module[{wr0=\[Pi],sign,a,p,e,x,En,Lz,Kc,En1,Lz1,\[CapitalOmega]r,\[CapitalOmega]\[Phi],rg,\[CapitalGamma],\[CapitalGamma]1,\[Omega],W,sumPlus0,sumMin0,sumPlus1,sumMin1,
    rp,Ur,\[Psi]ref,phase,expr,expr1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,V,dV,RIn0,dRIn0,ddRIn0,dddRIn0,RUp0,dRUp0,ddRUp0,dddRUp0,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,\[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    vl,vn,vm,vmb,Sln,Snm,Snmb,Slmb,Smmb,rho,beta,pi,tau,mu,gamma,alpha,Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cmnn,Cmnmb,Cmmbmb,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,rp1,Urp1,\[CapitalSigma]1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,CPlus0,CMin0,CPlus1,CMin1,A2g,A1g,A0g,A2par,A1par,A0par,B3par,B2par,B1par},
    
     {a,p,e,x}=orbit["OrbitalElements"];
    (* Geodesic constants of motion *)
    En = orbit["Eg"]; 
    Lz = orbit["Lzg"];
    Kc = (Lz-a*En)^2;
    (* Linear corrections to the constants of motion. These terms only appears in the projections to the Kinnersley tetrad of the velocity vn1 and vmb1*)
    En1 = orbit["Es"]; 
    Lz1 = orbit["Jzs"]; 
    (* Geodesic BL frequencies *)
    \[CapitalOmega]r = orbit["BLFrequenciesGeo"][[1]];
    \[CapitalOmega]\[Phi] = orbit["BLFrequenciesGeo"][[2]];
    (* Linear corrections to the frequencies *)
    \[CapitalGamma] = orbit["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
    \[CapitalGamma]1 = orbit["MinoFrequenciesCorrection"][[1]];
    \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r;
       
    rp = orbit["rg"][wr+wr0];
    Ur = sign*Sqrt[((rp^2+a^2)*En-a*Lz)^2-(rp^2-2rp+a^2)*(rp^2+Kc)];
    (* phase = \[Omega]*orbit["\[CapitalDelta]trg"][wr+wr0]-m*orbit["\[CapitalDelta]\[Phi]rg"][wr+wr0];*)
    \[Psi]ref=\[CapitalDelta]\[Psi]funref[{\[Omega],m},{En,Lz},a,p,e,x];
    phase=\[CapitalDelta]\[Psi]ICr1gfun[wr+wr0,{\[Omega],m},{En,Lz},\[Psi]ref,a,p,e,x];
    expr = Exp[sign*I*(phase+n*wr)];(* Geodesic exponential term with \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    \[CapitalDelta] = rp^2-2rp+a^2;
    KTeu = (rp^2+a^2)*\[Omega]-a*m;
    d\[CapitalDelta] = 2*(rp-1);
    dKTeu = 2*rp*\[Omega];
    V = -(KTeu^2+4I*(rp-1)*KTeu)/\[CapitalDelta]+8*I*\[Omega]*rp+\[Lambda]; (* Potential in radial Teukolsky equation *)
    dV = -((2KTeu*dKTeu+4I*KTeu+4I*(rp-1)*dKTeu)*\[CapitalDelta]-(KTeu^2+4I*(rp-1)*KTeu)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
    RIn0 = RInsol0[rp]; (* radial function *)
    dRIn0 = dRInsol0[rp];
    ddRIn0 = (V*RIn0+d\[CapitalDelta]*dRIn0)/\[CapitalDelta]; (* second derivative of radial function from Teukolsky equation *)
    dddRIn0 = 1/\[CapitalDelta] (dV*RIn0+(V+2)*dRIn0); (* third derivative of radial function from Teukolsky equation *)
    RUp0 = RUpsol0[rp];
    dRUp0 = dRUpsol0[rp];
    ddRUp0 = (V*RUp0+d\[CapitalDelta]*dRUp0)/\[CapitalDelta];  
    dddRUp0 = 1/\[CapitalDelta](dV*RUp0+(V+2)*dRUp0);
    (* Spin-weighted spheroidal harmonics S(\[Theta](z)) *)
    S = SWSH[\[Pi]/2][[1]]; 
    dSd\[Theta] = SWSH[\[Pi]/2][[2]]; (* First derivative of S wrt \[Theta] *)
    d2Sd\[Theta]2 = S((a*\[Omega])^2-2(a*\[Omega])*m+m^2+2-\[Lambda]); (* second derivative of S from Teukolsky equation *)
    L2S = dSd\[Theta]-S(m-a*\[Omega]); (* Operators acting on S(\[Theta]) *)
    L1L2S = d2Sd\[Theta]2+(dSd\[Theta](-2m+2a*\[Omega]))+S(-2+m^2-2a*m*\[Omega]+a^2*\[Omega]^2);   
    
    \[Zeta] = rp;
    \[Zeta]bar = rp;
    \[CapitalSigma] = rp^2;
    {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[rp,rp,1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,S,L2S,L1L2S,a];
    {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,S,L2S,L1L2S,a,\[Omega]];
    (* Four-velocity in Kinnersley tetrad *)
    vn = -((rp^2+a^2)*En - a*Lz + Ur)/(2*\[CapitalSigma]); 
    vl = -((rp^2+a^2)*En - a*Lz - Ur)/(\[CapitalDelta]);
    vm = (I*(a*En - Lz))/(-Sqrt[2]*\[Zeta]bar);
    vmb = (I*(a*En - Lz))/(Sqrt[2]*\[Zeta]bar);
    Cmnn = vn^2;
    Cmnmb = vn*vmb;
    Cmmbmb = vmb^2;
    (* Spin tensor in Kinnersley tetrad *)
    Sln  = -Sqrt[Kc]/rp;
    Snm  = (\[Zeta]/Sqrt[Kc])*vm*vn;
    Snmb = Conjugate[Snm];
    Slmb = (-(\[Zeta]/Sqrt[Kc]))*vl*vmb;
    Smmb = 0;
	(* Spin coefficients *)
    rho = 1/\[Zeta];
    beta = 0;
    pi = -((I*a)/(\[Zeta]^2*Sqrt[2]));
    tau = (I*a)/(Sqrt[2]*\[CapitalSigma]);
    mu = \[CapitalDelta]/(2\[Zeta]^2*\[Zeta]bar);
    gamma = (a^2-rp)/(2\[Zeta]^2*\[Zeta]bar);
    alpha = -I*a/(Sqrt[2]\[Zeta]^2);
    Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
    Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
    Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
    Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
    Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
        -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
    St\[Phi]n  = -I*KTeu/(2\[CapitalSigma])*Sln+(a*\[Omega]-m)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
    St\[Phi]mb = -I*KTeu*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]-m)/(Sqrt[2]*\[Zeta])*Smmb;
    Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
    Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
    
    rp1 = orbit["\[Delta]rpar"][wr+wr0];
    (*The complex exponential only serves to change the sign of the spin correction to the radial velocity with different initial conditions (encoded by wr0). 
    It is a makeshift solution, by it is good enough for my purpose. Exp[I*wr0] should never take complex values. *)
    Urp1 = sign*Exp[I*wr0]*orbit["\[Delta]vrpar"][wr+wr0];
    \[CapitalSigma]1 = 2*rp*rp1; (* Linear correction to \[CapitalSigma] *)
    expr1 = I*sign*(\[Omega]*orbit["\[CapitalDelta]\[Delta]tpar"][wr+wr0]-m*orbit["\[CapitalDelta]\[Delta]\[Phi]par"][wr+wr0]);(* Linear parts of the exponential terms with \[CapitalDelta]t and \[CapitalDelta]\[Phi] *)
  
    (* Linear parts of the four-velocity in Kinnersley tetrad *)
    vn1 =  1/(2Sqrt[Kc]*rp^3)(Lz^2+a^2(-En1*Sqrt[Kc]*rp+En^2(1+rp)+2En*Sqrt[Kc]*rp1)-a(En*Lz(2+rp)+Sqrt[Kc](-Lz1*rp+2Lz*rp1))-Sqrt[Kc](En1*rp^3-2rp1*Ur+rp*Urp1));
    vmb1 = -I(En*Sqrt[Kc]*rp-(a*En-Lz)(a*En1-Lz1)rp+Kc*rp1)/(Sqrt[2](a*En-Lz)rp^2);
    
    Ann0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn + 2*vn1)*vn + Cdnn;
    Annt\[Phi]S = (St\[Phi]n + expr1*vn)*vn;
    AnnrS  = (Srn + rp1*vn)*vn;
    Anmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb);
    Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + expr1*vn*vmb);
    AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
    Ambmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vmb + 2*vmb1)*vmb + Cdmbmb;
    Ambmbt\[Phi]S = (St\[Phi]mb + expr1*vmb)*vmb;
    AmbmbrS  = (Srmb + rp1*vmb)*vmb;
    
    A2g = Cmmbmb*fmbmb2;
    A1g = (Cmnmb*fnmb1+Cmmbmb*fmbmb1);
    A0g = (Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0);
    
    (*These coefficients are different respect to the ones define in the companion paper. Here I put together everthing, including the differentials \[Delta]r, \[Delta]z*)
    B3par = -AmbmbrS*fmbmb2;
    B2par = -AnmbrS*fnmb1-AmbmbrS*fmbmb1;
    B1par = -AmbmbrS*fmbmb0-AnmbrS*fnmb0-AnnrS*fnn0;
    
	A2par = AmbmbrS*dfmbmb2dr+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb2;
    A1par = AmbmbrS*dfmbmb1dr+AnmbrS*dfnmb1dr+(Ambmb0S+Ambmbt\[Phi]S)fmbmb1+(Anmb0S+Anmbt\[Phi]S)fnmb1;
    A0par = AmbmbrS*dfmbmb0dr+AnmbrS*dfnmb0dr+AnnrS*dfnn0dr+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb0+(Anmb0S+Anmbt\[Phi]S)*fnmb0+(Ann0S+Annt\[Phi]S)*fnn0;
    
    sumPlus0 = (\[CapitalSigma]*(A0g*RIn0 - A1g*dRIn0 + A2g*ddRIn0)*expr/.sign->1) + (\[CapitalSigma]*(A0g*RIn0 - A1g*dRIn0 + A2g*ddRIn0)*expr/.sign->-1);
    sumMin0 = (\[CapitalSigma]*(A0g*RUp0 - A1g*dRUp0 + A2g*ddRUp0)*expr/.sign->1) + (\[CapitalSigma]*(A0g*RUp0 - A1g*dRUp0 + A2g*ddRUp0)*expr/.sign->-1);  
    sumPlus1 = (\[CapitalSigma]*(A0par*RIn0 - (A1par+B1par)*dRIn0 + (A2par+B2par)ddRIn0 - B3par*dddRIn0)*expr/.sign->1) + (\[CapitalSigma]*(A0par*RIn0 - (A1par+B1par)*dRIn0 + (A2par+B2par)ddRIn0 - B3par*dddRIn0)*expr/.sign->-1);
    sumMin1 = (\[CapitalSigma]*(A0par*RUp0 - (A1par+B1par)*dRUp0 + (A2par+B2par)ddRUp0 - B3par*dddRUp0)*expr/.sign->1) + (\[CapitalSigma]*(A0par*RUp0 - (A1par+B1par)*dRUp0 + (A2par+B2par)ddRUp0 - B3par*dddRUp0)*expr/.sign->-1);      
    
    W = (RInsol0[p/(1-e)]*dRUpsol0[p/(1-e)] - dRInsol0[p/(1-e)]*RUpsol0[p/(1-e)])/((p/(1-e))^2-2(p/(1-e))+a^2); (* Invariant Wronskian *)
    
    CPlus0 = 2*Pi*sumPlus0/(\[CapitalGamma]*W); (* Geodesic amplitudes *)
    CMin0 = 2*Pi*sumMin0/(\[CapitalGamma]*W);
    CPlus1 = 2*Pi*(sumPlus1 - \[CapitalGamma]1/\[CapitalGamma]*sumPlus0)/(\[CapitalGamma]*W); (* Linear parts of the amplitudes *)
    CMin1 = 2*Pi*(sumMin1 - \[CapitalGamma]1/\[CapitalGamma]*sumMin0)/(\[CapitalGamma]*W);
    
    -{CPlus0,CMin0,CPlus1,CMin1}
]


(* ::Subsection::Closed:: *)
(*Generic spin-gauge (FF, FT, FTE)*)


(* ::Text:: *)
(*The sign of the amplitudes now follow the same convention of the BHPToolkit (the sign convention was changed after the first version of this code).*)
(*This version is valid for any i.c. *)


SourceTermEccEqGenericGauge[wr_,l_,m_,n_,orbit_,wr0_,{\[Lambda]_,\[Lambda]1_,SWSH_,SWSHS_},{{RInsol0_,dRInsol0_},{RInsol1_,dRInsol1_}},{{RUpsol0_,dRUpsol0_},{RUpsol1_,dRUpsol1_}}]:=
Module[{sign,a,p,e,x,En,Lz,Kc,En1,Lz1,\[CapitalOmega]r,\[CapitalOmega]\[Phi],\[CapitalOmega]r1,\[CapitalOmega]\[Phi]1,\[CapitalGamma],\[CapitalGamma]1,\[Omega],\[Omega]1,W,W1,sumPlus0,sumMin0,sumPlus1,sumMin1,
    rp,Ur,\[Psi]ref,phase,expr,expr1,\[CapitalDelta],d\[CapitalDelta],KTeu,KTeu1,dKTeu,dKTeu1,V,V1,dV,RIn0,dRIn0,ddRIn0,RIn1,dRIn1,ddRIn1,dddRIn0,
    RUp0,dRUp0,ddRUp0,RUp1,dRUp1,ddRUp1,dddRUp0,S,S1,L2S,L2S1,L1L2S,L1L2S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,\[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    vl,vn,vm,vmb,Sln,Snm,Snmb,Slmb,Smmb,rho,beta,pi,tau,mu,gamma,alpha,Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cmnn,Cmnmb,Cmmbmb,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,rp1,Urp1,\[CapitalSigma]1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,CPlus0,CMin0,CPlus1,CMin1,A2g,A1g,A0g,A2par,A1par,A0par,B3par,B2par,B1par},
    
    {a,p,e,x}=orbit["OrbitalElements"];
    (* Geodesic constants of motion *)
    En = orbit["Eg"]; 
    Lz = orbit["Lzg"];
    Kc = (Lz-a*En)^2;
    (* Linear corrections to the constants of motion. These terms only appears in the projections to the Kinnersley tetrad of the velocity vn1 and vmb1*)
    En1 = orbit["Es"]; 
    Lz1 = orbit["Jzs"]; 
    (* Geodesic BL frequencies *)
    \[CapitalOmega]r = orbit["BLFrequenciesGeo"][[1]];
    \[CapitalOmega]\[Phi] = orbit["BLFrequenciesGeo"][[2]];
    (* Linear corrections to the frequencies *)
    {\[CapitalOmega]r1,\[CapitalOmega]\[Phi]1} = orbit["BLFrequenciesCorrection"]; 
    \[CapitalGamma] = orbit["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
    \[CapitalGamma]1 = orbit["MinoFrequenciesCorrection"][[1]];
    \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r;
    \[Omega]1 = m*\[CapitalOmega]\[Phi]1 + n*\[CapitalOmega]r1;
       
    rp = orbit["rg"][wr+wr0];
    Ur = sign*Sqrt[((rp^2+a^2)*En-a*Lz)^2-(rp^2-2rp+a^2)*(rp^2+Kc)];
    (* phase = \[Omega]*orbit["\[CapitalDelta]trg"][wr+wr0]-m*orbit["\[CapitalDelta]\[Phi]rg"][wr+wr0];*)
    \[Psi]ref=\[CapitalDelta]\[Psi]funref[{\[Omega],m},{En,Lz},a,p,e,x];
    phase=\[CapitalDelta]\[Psi]ICr1gfun[wr+\[Pi],{\[Omega],m},{En,Lz},\[Psi]ref,a,p,e,x];
    expr = Exp[sign*I*(phase+n*wr)];(* Geodesic exponential term with \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    \[CapitalDelta] = rp^2-2rp+a^2;
    KTeu = (rp^2+a^2)*\[Omega]-a*m;
    KTeu1 = (rp^2+a^2)*\[Omega]1;
    d\[CapitalDelta] = 2*(rp-1);
    dKTeu = 2*rp*\[Omega];
    dKTeu1 = 2*rp*\[Omega]1;
    V = -(KTeu^2+4I*(rp-1)*KTeu)/\[CapitalDelta]+8*I*\[Omega]*rp+\[Lambda]; (* Potential in radial Teukolsky equation *)
    V1 = -(2*KTeu*KTeu1+4I*(rp-1)*KTeu1)/\[CapitalDelta]+8*I*\[Omega]1*rp+\[Lambda]1; (* Linear part of the potential in radial Teukolsky equation *)
    dV = -((2KTeu*dKTeu+4I*KTeu+4I*(rp-1)*dKTeu)*\[CapitalDelta]-(KTeu^2+4I*(rp-1)*KTeu)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
    RIn0 = RInsol0[rp]; (* radial function *)
    dRIn0 = dRInsol0[rp];
    ddRIn0 = (V*RIn0+d\[CapitalDelta]*dRIn0)/\[CapitalDelta]; (* second derivative of radial function from Teukolsky equation *)
    RIn1 = RInsol1[rp]; (* Linear part of the radial function *)
    dRIn1 = dRInsol1[rp];
    ddRIn1 = (V*RIn1+V1*RIn0+d\[CapitalDelta]*dRIn1)/\[CapitalDelta];
    dddRIn0 = 1/\[CapitalDelta] (dV*RIn0+(V+2)*dRIn0); (* third derivative of radial function from Teukolsky equation *)
    RUp0 = RUpsol0[rp];
    dRUp0 = dRUpsol0[rp];
    ddRUp0 = (V*RUp0+d\[CapitalDelta]*dRUp0)/\[CapitalDelta];  
    RUp1 = RUpsol1[rp];
    dRUp1 = dRUpsol1[rp];
    ddRUp1 = (V*RUp1+V1*RUp0+d\[CapitalDelta]*dRUp1)/\[CapitalDelta];
    dddRUp0 = 1/\[CapitalDelta](dV*RUp0+(V+2)*dRUp0);
    (* Spin-weighted spheroidal harmonics S(\[Theta](z)) *)
    S = SWSH[\[Pi]/2][[1]]; 
    S1 = SWSHS[\[Pi]/2][[1]]; (* Linear part of S(\[Theta](z)) *)
    dSd\[Theta] = SWSH[\[Pi]/2][[2]]; (* First derivative of S wrt \[Theta] *)
    dS1d\[Theta] = SWSHS[\[Pi]/2][[2]]; (* Linear part of first derivative of S wrt \[Theta] *)
    d2Sd\[Theta]2 = S((a*\[Omega])^2-2(a*\[Omega])*m+m^2+2-\[Lambda]); (* second derivative of S from Teukolsky equation *)
    d2S1d\[Theta]2 = -(-a^2*\[Omega]^2-m^2-2+2m*a*\[Omega]+\[Lambda])S1-(2a*\[Omega]1(-a*\[Omega]+m)+\[Lambda]1)S; (* Linear part of second derivative of S from Teukolsky equation *)
    L2S = dSd\[Theta]-S(m-a*\[Omega]); (* Operators acting on S(\[Theta]) *)
    L2S1 = dS1d\[Theta]-S1(m-a*\[Omega])+S(a*\[Omega]1);(* Operators acting on S(\[Theta]) *)
    L1L2S = d2Sd\[Theta]2+(dSd\[Theta](-2m+2a*\[Omega]))+S(-2+m^2-2a*m*\[Omega]+a^2*\[Omega]^2);  
    L1L2S1 = d2S1d\[Theta]2+(dS1d\[Theta](-2m+2a*\[Omega])+dSd\[Theta](2a*\[Omega]1)) + S1(-2+m^2-2a*m*\[Omega]+a^2*\[Omega]^2) + S(-2a*m*\[Omega]1 + a^2*2*\[Omega]*\[Omega]1);  
    
    \[Zeta] = rp;
    \[Zeta]bar = rp;
    \[CapitalSigma] = rp^2;
    {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[rp,rp,1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,S,L2S,L1L2S,a];
    {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,S,L2S,L1L2S,a,\[Omega]];
    {fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21} = dfabidS[\[Zeta],\[Zeta]bar,1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,KTeu1,dKTeu1,S,L2S,L1L2S,S1,L2S1,L1L2S1,a];
    (* Four-velocity in Kinnersley tetrad *)
    vn = -((rp^2+a^2)*En - a*Lz + Ur)/(2*\[CapitalSigma]); 
    vl = -((rp^2+a^2)*En - a*Lz - Ur)/(\[CapitalDelta]);
    vm = (I*(a*En - Lz))/(-Sqrt[2]*\[Zeta]bar);
    vmb = (I*(a*En - Lz))/(Sqrt[2]*\[Zeta]bar);
    Cmnn = vn^2;
    Cmnmb = vn*vmb;
    Cmmbmb = vmb^2;
    (* Spin tensor in Kinnersley tetrad *)
    Sln  = -Sqrt[Kc]/rp;
    Snm  = (\[Zeta]/Sqrt[Kc])*vm*vn;
    Snmb = Conjugate[Snm];
    Slmb = (-(\[Zeta]/Sqrt[Kc]))*vl*vmb;
    Smmb = 0;
	(* Spin coefficients *)
    rho = 1/\[Zeta];
    beta = 0;
    pi = -((I*a)/(\[Zeta]^2*Sqrt[2]));
    tau = (I*a)/(Sqrt[2]*\[CapitalSigma]);
    mu = \[CapitalDelta]/(2\[Zeta]^2*\[Zeta]bar);
    gamma = (a^2-rp)/(2\[Zeta]^2*\[Zeta]bar);
    alpha = -I*a/(Sqrt[2]\[Zeta]^2);
    Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
    Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
    Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
    Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
    Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
        -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
    St\[Phi]n  = -I*KTeu/(2\[CapitalSigma])*Sln+(a*\[Omega]-m)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
    St\[Phi]mb = -I*KTeu*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]-m)/(Sqrt[2]*\[Zeta])*Smmb;
    Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
    Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
    
    rp1 = orbit["\[Delta]rpar"][wr+wr0];
    (*The complex exponential only serves to change the sign of the spin correction to the radial velocity with different initial conditions (encoded by wr0). 
    It is a makeshift solution, by it is good enough for my purpose. Exp[I*wr0] should never take complex values. *)
    Urp1 = sign*Exp[I*wr0]*orbit["\[Delta]vrpar"][wr+wr0];
    \[CapitalSigma]1 = 2*rp*rp1; (* Linear correction to \[CapitalSigma] *)
    expr1 = I*sign*(\[Omega]*orbit["\[CapitalDelta]\[Delta]tpar"][wr+wr0]-m*orbit["\[CapitalDelta]\[Delta]\[Phi]par"][wr+wr0])+I*sign*\[Omega]1*orbit["\[CapitalDelta]trg"][wr+wr0];(* Linear parts of the exponential terms with \[CapitalDelta]t and \[CapitalDelta]\[Phi] *)
    (* Linear parts of the four-velocity in Kinnersley tetrad *)
    vn1 = 1/(2Sqrt[Kc]*rp^3)(Lz^2+a^2(-En1*Sqrt[Kc]*rp+En^2(1+rp)+2En*Sqrt[Kc]*rp1)-a(En*Lz(2+rp)+Sqrt[Kc](-Lz1*rp+2Lz*rp1))-Sqrt[Kc](En1*rp^3-2rp1*Ur+rp*Urp1));
    vmb1 = -I(En*Sqrt[Kc]*rp-(a*En-Lz)(a*En1-Lz1)rp+Kc*rp1)/(Sqrt[2](a*En-Lz)rp^2);
    
    Ann0S = ((\[CapitalSigma]1/\[CapitalSigma])*vn + 2*vn1)*vn + Cdnn;
    Annt\[Phi]S = (St\[Phi]n + expr1*vn)*vn;
    AnnrS = (Srn + rp1*vn)*vn;
    Anmb0S = ((\[CapitalSigma]1/\[CapitalSigma])*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb);
    Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + expr1*vn*vmb);
    AnmbrS = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
    Ambmb0S = ((\[CapitalSigma]1/\[CapitalSigma])*vmb + 2*vmb1)*vmb + Cdmbmb;
    Ambmbt\[Phi]S = (St\[Phi]mb + expr1*vmb)*vmb;
    AmbmbrS = (Srmb + rp1*vmb)*vmb;
    
    A2g = Cmmbmb*fmbmb2;
    A1g = (Cmnmb*fnmb1+Cmmbmb*fmbmb1);
    A0g = (Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0);
    
    (*These coefficients are different respect to the ones define in the companion paper. Here I put together everthing, including the differentials \[Delta]r, \[Delta]z*)
    B3par = -AmbmbrS*fmbmb2;
    B2par = -AnmbrS*fnmb1-AmbmbrS*fmbmb1;
    B1par = -AmbmbrS*fmbmb0-AnmbrS*fnmb0-AnnrS*fnn0;
    
    A2par = AmbmbrS*dfmbmb2dr+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb2 +Cmmbmb*fmbmb21;
    A1par = AmbmbrS*dfmbmb1dr+AnmbrS*dfnmb1dr+(Ambmb0S+Ambmbt\[Phi]S)fmbmb1+(Anmb0S+Anmbt\[Phi]S)fnmb1+(Cmmbmb*fmbmb11+Cmnmb*fnmb11);
    A0par = AmbmbrS*dfmbmb0dr+AnmbrS*dfnmb0dr+AnnrS*dfnn0dr+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb0+(Anmb0S+Anmbt\[Phi]S)*fnmb0+(Ann0S+Annt\[Phi]S)*fnn0+(Cmmbmb*fmbmb01+Cmnmb*fnmb01+Cmnn*fnn01);
    
    sumPlus0 = (\[CapitalSigma]*(A0g*RIn0 - A1g*dRIn0 + A2g*ddRIn0)*expr/.sign->1) + (\[CapitalSigma]*(A0g*RIn0 - A1g*dRIn0 + A2g*ddRIn0)*expr/.sign->-1);
    sumMin0 = (\[CapitalSigma]*(A0g*RUp0 - A1g*dRUp0 + A2g*ddRUp0)*expr/.sign->1) + (\[CapitalSigma]*(A0g*RUp0 - A1g*dRUp0 + A2g*ddRUp0)*expr/.sign->-1);
    sumPlus1 = (\[CapitalSigma]*(A0par*RIn0 - (A1par+B1par)*dRIn0 + (A2par+B2par)ddRIn0 - B3par*dddRIn0)*expr + \[CapitalSigma]*(A0g*RIn1 - A1g*dRIn1 + A2g*ddRIn1)*expr/.sign->1) + (\[CapitalSigma]*(A0par*RIn0 - (A1par+B1par)*dRIn0 + (A2par+B2par)ddRIn0 - B3par*dddRIn0)*expr + \[CapitalSigma]*(A0g*RIn1 - A1g*dRIn1 + A2g*ddRIn1)*expr/.sign->-1);
    sumMin1 = (\[CapitalSigma]*(A0par*RUp0 - (A1par+B1par)*dRUp0 + (A2par+B2par)ddRUp0 - B3par*dddRUp0)*expr + \[CapitalSigma]*(A0g*RUp1 - A1g*dRUp1 + A2g*ddRUp1)*expr/.sign->1) + (\[CapitalSigma]*(A0par*RUp0 - (A1par+B1par)*dRUp0 + (A2par+B2par)ddRUp0 - B3par*dddRUp0)*expr + \[CapitalSigma]*(A0g*RUp1 - A1g*dRUp1 + A2g*ddRUp1)*expr/.sign->-1);      
    
    W = (RInsol0[p/(1-e)]*dRUpsol0[p/(1-e)] - dRInsol0[p/(1-e)]*RUpsol0[p/(1-e)])/((p/(1-e))^2-2(p/(1-e))+a^2); (* Invariant Wronskian *)
    W1 = (RInsol1[p/(1-e)]*dRUpsol0[p/(1-e)] + RInsol0[p/(1-e)]*dRUpsol1[p/(1-e)] - dRInsol1[p/(1-e)]*RUpsol0[p/(1-e)] -dRInsol0[p/(1-e)]*RUpsol1[p/(1-e)])/((p/(1-e))^2-2(p/(1-e))+a^2); (* Linear part of the invariant Wronskian *)
    
    CPlus0 = 2*Pi*sumPlus0/(\[CapitalGamma]*W); (* Geodesic amplitudes *)
    CMin0 = 2*Pi*sumMin0/(\[CapitalGamma]*W);
    CPlus1 = 2*Pi*(sumPlus1 - \[CapitalGamma]1/\[CapitalGamma]*sumPlus0 - W1/W*sumPlus0)/(\[CapitalGamma]*W); (* Linear parts of the amplitudes *)
    CMin1 = 2*Pi*(sumMin1 - \[CapitalGamma]1/\[CapitalGamma]*sumMin0 - W1/W*sumMin0)/(\[CapitalGamma]*W);
    
    -{CPlus0,CMin0,CPlus1,CMin1}
]


(* ::Section::Closed:: *)
(*Fluxes for a spinning particle - eccentric orbits*)


(* ::Subsection::Closed:: *)
(*Geodesic*)


(* ::Text:: *)
(*The sign of the amplitudes now follow the same convention of the BHPToolkit (the sign convention was changed after the first version of this code).*)


WASPEccEq0PA[l_,m_,n_,orbit_,tol_]:=
Module[{a,p,e,x,\[CapitalOmega]r,\[CapitalOmega]\[Phi],\[CapitalGamma],\[Omega],\[Lambda],SWSH,RInsol0,dRInsol0,RUpsol0,dRUpsol0,\[ScriptCapitalC]2,rplus,\[Kappa],\[Epsilon],\[Alpha],nInt,wrlist,sampledFunc,sumFun,i,CPlus0,CMin0,relerr},
	{a,p,e,x}=orbit["OrbitalElements"];
    (* Geodesic BL frequencies *)
    \[CapitalOmega]r = orbit["BLFrequenciesGeo"][[1]];
    \[CapitalOmega]\[Phi] = orbit["BLFrequenciesGeo"][[2]];
    \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r;

	{\[Lambda],SWSH} = angpar0[-2,l,m,a,\[Omega]]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
	{RInsol0,dRInsol0} = TsolverIn[-2,l,m,a,\[Omega],\[Lambda],p/(1-e),All][[1;;2]];
    {RUpsol0,dRUpsol0} = TsolverUpFromTSidentity[-2,l,m,a,\[Omega],\[Lambda],p/(1+e),All];
    
    \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (* TS constant *)
    rplus = 1+Sqrt[1-a^2]; (* horizon r_+  *)
    \[Kappa] = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
    \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
    \[Alpha] = 256*(2*rplus)^5*\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
	
	nInt[nsamp_]:=2^nsamp*2^2;
	wrlist[nsamp_]:=N[Table[wr,{wr,Pi/(nInt[nsamp]),Pi,2Pi/nInt[nsamp]}]];
	
	sampledFunc[0]=SourceTermEccEq0PA[wrlist[0],l,m,n,orbit,{\[Lambda],SWSH},{RInsol0,dRInsol0},{RUpsol0,dRUpsol0}];
	sampledFunc[nsamp_]:=sampledFunc[nsamp]=SourceTermEccEq0PA[wrlist[nsamp],l,m,n,orbit,{\[Lambda],SWSH},{RInsol0,dRInsol0},{RUpsol0,dRUpsol0}];
	
    sumFun[0]=1/nInt[0]*{Total[sampledFunc[0][[1]],Method->"CompensatedSummation"],Total[sampledFunc[0][[2]],Method->"CompensatedSummation"]};
	sumFun[1]=1/nInt[1]*{Total[sampledFunc[1][[1]],Method->"CompensatedSummation"],Total[sampledFunc[1][[2]],Method->"CompensatedSummation"]};
	sumFun[nsamp_]:=sumFun[nsamp]=1/nInt[nsamp]*{Total[sampledFunc[nsamp][[1]],Method->"CompensatedSummation"],Total[sampledFunc[nsamp][[2]],Method->"CompensatedSummation"]};
	
	relerr[0]=Abs[1-sumFun[1][[1]]/sumFun[0][[1]]];
	relerr[nsamp_]:=relerr[nsamp]=Abs[1-sumFun[nsamp][[1]]/sumFun[nsamp-1][[1]]];
	
	i=1;
	While[relerr[i]>10^(-tol)&&i<=10,i++];
	
	CPlus0=sumFun[i][[1]];
	CMin0=sumFun[i][[2]];
	
	Remove[RInsol0,dRInsol0,RUpsol0,dRUpsol0,\[Lambda],SWSH,sampledFunc,wrlist,sumFun];

  Association[
    "gg"->{nInt[i],relerr[i]},
    "Amplitudes"->
    <|
      "I"->CPlus0,
      "H"->CMin0
    |>,
    "Fluxes"->
    <|
      "En0PA"-><|
        "I"->Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "H"->\[Alpha]*Abs[CMin0]^2/(4Pi*\[Omega]^2)
      |>,
      "Jz0PA"->
      <|
        "I"->Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "H"->\[Alpha]*Abs[CMin0]^2*m/(4Pi*\[Omega]^3)
      |>
    |>
  ]
]


(* ::Subsection::Closed:: *)
(*Fixed frequency*)


(* ::Text:: *)
(*The sign of the amplitudes now follow the same convention of the BHPToolkit (the sign convention was changed after the first version of this code).*)


WASPEccEqFF[l_,m_,n_,orbit_,tol_]:=
Module[{a,p,e,x,\[CapitalOmega]r,\[CapitalOmega]\[Phi],\[CapitalGamma],\[Omega],\[Lambda],SWSH,RInsol0,dRInsol0,RUpsol0,dRUpsol0,\[ScriptCapitalC]2,rplus,\[Kappa],\[Epsilon],\[Alpha],nInt,wrlist,sampledFunc,sumFun,i,CPlus0,CMin0,CPlus1,CMin1,relerr},
	{a,p,e,x}=orbit["OrbitalElements"];
    (* Geodesic BL frequencies *)
    \[CapitalOmega]r = orbit["BLFrequenciesGeo"][[1]];
    \[CapitalOmega]\[Phi] = orbit["BLFrequenciesGeo"][[2]];
    \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r;

	{\[Lambda],SWSH} = angpar0[-2,l,m,a,\[Omega]]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
	{RInsol0,dRInsol0} = TsolverIn[-2,l,m,a,\[Omega],\[Lambda],p/(1-e),All][[1;;2]];
    {RUpsol0,dRUpsol0} = TsolverUpFromTSidentity[-2,l,m,a,\[Omega],\[Lambda],p/(1+e),All][[1;;2]];
    
    \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (* TS constant *)
    rplus = 1+Sqrt[1-a^2]; (* horizon r_+  *)
    \[Kappa] = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
    \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
    \[Alpha] = 256*(2*rplus)^5*\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
	
	nInt[nsamp_]:=2^nsamp*2^2;
	wrlist[nsamp_]:=N[Table[wr,{wr,Pi/(nInt[nsamp]),Pi,2Pi/nInt[nsamp]}]];
	
	sampledFunc[0]=SourceTermEccEqFF[wrlist[0],l,m,n,orbit,{\[Lambda],SWSH},{RInsol0,dRInsol0},{RUpsol0,dRUpsol0}];
	sampledFunc[nsamp_]:=sampledFunc[nsamp]=SourceTermEccEqFF[wrlist[nsamp],l,m,n,orbit,{\[Lambda],SWSH},{RInsol0,dRInsol0},{RUpsol0,dRUpsol0}];
	
    sumFun[0]=1/nInt[0]*{Total[sampledFunc[0][[1]],Method->"CompensatedSummation"],Total[sampledFunc[0][[2]],Method->"CompensatedSummation"],Total[sampledFunc[0][[3]],Method->"CompensatedSummation"],Total[sampledFunc[0][[4]],Method->"CompensatedSummation"]};
	sumFun[1]=1/nInt[1]*{Total[sampledFunc[1][[1]],Method->"CompensatedSummation"],Total[sampledFunc[1][[2]],Method->"CompensatedSummation"],Total[sampledFunc[1][[3]],Method->"CompensatedSummation"],Total[sampledFunc[1][[4]],Method->"CompensatedSummation"]};
	sumFun[nsamp_]:=sumFun[nsamp]=1/nInt[nsamp]*{Total[sampledFunc[nsamp][[1]],Method->"CompensatedSummation"],Total[sampledFunc[nsamp][[2]],Method->"CompensatedSummation"],Total[sampledFunc[nsamp][[3]],Method->"CompensatedSummation"],Total[sampledFunc[nsamp][[4]],Method->"CompensatedSummation"]};
	
	relerr[0]=Abs[1-sumFun[1][[1]]/sumFun[0][[1]]];
	relerr[nsamp_]:=relerr[nsamp]=Abs[1-sumFun[nsamp][[1]]/sumFun[nsamp-1][[1]]];
	
	i=1;
	While[relerr[i]>10^(-tol)&&i<=10,i++];
	
	CPlus0=sumFun[i][[1]];
	CMin0=sumFun[i][[2]];
	CPlus1=sumFun[i][[3]];
	CMin1=sumFun[i][[4]];
	
	Remove[RInsol0,dRInsol0,RUpsol0,dRUpsol0,\[Lambda],SWSH,sampledFunc,wrlist,sumFun];

  Association[
    "gg"->{nInt[i],relerr[i]},
    "Amplitudes"->
    <|
      "I"->CPlus0,
      "H"->CMin0
    |>,
    "AmplitudesCorrection"->
    <|
      "I"->CPlus1,
      "H"->CMin1
    |>,
    "Fluxes"->
    <|
      "En0PA"-><|
        "I"->Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "H"->\[Alpha]*Abs[CMin0]^2/(4Pi*\[Omega]^2)
      |>,
      "Jz0PA"->
      <|
        "I"->Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "H"->\[Alpha]*Abs[CMin0]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection"->
    <|
      "En1PA"-><|
        "I"->(2*Re[CPlus1*Conjugate[CPlus0]])/(4Pi*\[Omega]^2),
        "H"->\[Alpha]*(2*Re[CMin1*Conjugate[CMin0]])/(4Pi*\[Omega]^2)
      |>,
      "Jz1PA"->
      <|
        "I"->(2*Re[CPlus1*Conjugate[CPlus0]])*m/(4Pi*\[Omega]^3),
        "H"->\[Alpha]*(2*Re[CMin1*Conjugate[CMin0]])*m/(4Pi*\[Omega]^3)
      |>
    |>
  ]
]


(* ::Subsection::Closed:: *)
(*Generic spin-gauge (FT, FE)*)


(* ::Text:: *)
(*The sign of the amplitudes now follow the same convention of the BHPToolkit (the sign convention was changed after the first version of this code).*)


WASPEccEqGenericGauge[l_,m_,n_,orbit_,wr0_,tol_]:=
Module[{a,p,e,x,\[CapitalOmega]r0,\[CapitalOmega]\[Phi]0,\[CapitalOmega]r1,\[CapitalOmega]\[Phi]1,\[Omega],\[Omega]1,\[Lambda],\[Lambda]1,SWSH,SWSHS,RInsol0,dRInsol0,RInsol1,dRInsol1,RUpsol0,dRUpsol0,RUpsol1,dRUpsol1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,\[Kappa],\[Epsilon],\[Alpha],\[Alpha]1,nInt,wrlist,sampledFunc,sumFun,i,CPlus0,CMin0,CPlus1,CMin1,relerr},
	{a,p,e,x}=orbit["OrbitalElements"];
	(* Geodesic BL frequencies *)
    \[CapitalOmega]r0 = orbit["BLFrequenciesGeo"][[1]];
    \[CapitalOmega]\[Phi]0 = orbit["BLFrequenciesGeo"][[2]];
    (* Linear corrections to the frequencies *)
    {\[CapitalOmega]r1,\[CapitalOmega]\[Phi]1} = orbit["BLFrequenciesCorrection"];
    \[Omega] = m*\[CapitalOmega]\[Phi]0 + n*\[CapitalOmega]r0;
    \[Omega]1 = m*\[CapitalOmega]\[Phi]1 + n*\[CapitalOmega]r1;
	
	{\[Lambda],\[Lambda]1,SWSH,SWSHS} = angpar[-2,l,m,a,\[Omega],\[Omega]1]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
    {RInsol0,dRInsol0,RInsol1,dRInsol1} = TsolverIn1spin[-2,l,m,a,\[Omega],\[Omega]1,\[Lambda],\[Lambda]1,p/(1-e),All][[1;;4]];
    {RUpsol0,dRUpsol0,RUpsol1,dRUpsol1} = TsolverUpFromTSidentity1spin[-2,l,m,a,\[Omega],\[Omega]1,\[Lambda],\[Lambda]1,p/(1+e),All][[1;;4]];
    
    \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (* TS constant *)
    \[ScriptCapitalC]21 = 4\[Lambda]^3*\[Lambda]1+4\[Lambda]^2(3\[Lambda]1+10a(m-2a*\[Omega])\[Omega]1)+8\[Lambda](\[Lambda]1(1+10a*m*\[Omega]-10a^2*\[Omega]^2)+6a(m+2a*\[Omega])\[Omega]1)+48\[Omega](a*m*\[Lambda]1+6\[Omega]1-18a^3*m*\[Omega]*\[Omega]1+12a^4*\[Omega]^2*\[Omega]1+a^2(\[Lambda]1*\[Omega]+6m^2*\[Omega]1));  (* linear part of the TS constant *)
    rplus = 1+Sqrt[1-a^2]; (* horizon r_+  *)
    \[Kappa] = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
    \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
    \[Alpha] = 256*(2*rplus)^5*\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
    \[Alpha]1 = -256*(2*rplus)^5*\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*\[ScriptCapitalC]21+256*(2*rplus)^5*(\[Omega]1*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3+\[Kappa]*(2*\[Kappa]*\[Omega]1)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3+\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(2*\[Kappa]*\[Omega]1)*\[Omega]^3+\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*3*\[Omega]^2*\[Omega]1)/\[ScriptCapitalC]2; (* linear part of the constant for horizon fluxes *)

	nInt[nsamp_]:=2^nsamp*2^2;
	wrlist[nsamp_]:=N[Table[wr,{wr,Pi/(nInt[nsamp]),Pi,2Pi/nInt[nsamp]}]];
	
	sampledFunc[0]=SourceTermEccEqGenericGauge[wrlist[0],l,m,n,orbit,wr0,{\[Lambda],\[Lambda]1,SWSH,SWSHS},{{RInsol0,dRInsol0},{RInsol1,dRInsol1}},{{RUpsol0,dRUpsol0},{RUpsol1,dRUpsol1}}];
	sampledFunc[nsamp_]:=sampledFunc[nsamp]=SourceTermEccEqGenericGauge[wrlist[nsamp],l,m,n,orbit,wr0,{\[Lambda],\[Lambda]1,SWSH,SWSHS},{{RInsol0,dRInsol0},{RInsol1,dRInsol1}},{{RUpsol0,dRUpsol0},{RUpsol1,dRUpsol1}}];
	
    sumFun[0]=1/nInt[0]*{Total[sampledFunc[0][[1]],Method->"CompensatedSummation"],Total[sampledFunc[0][[2]],Method->"CompensatedSummation"],Total[sampledFunc[0][[3]],Method->"CompensatedSummation"],Total[sampledFunc[0][[4]],Method->"CompensatedSummation"]};
	sumFun[1]=1/nInt[1]*{Total[sampledFunc[1][[1]],Method->"CompensatedSummation"],Total[sampledFunc[1][[2]],Method->"CompensatedSummation"],Total[sampledFunc[1][[3]],Method->"CompensatedSummation"],Total[sampledFunc[1][[4]],Method->"CompensatedSummation"]};
	sumFun[nsamp_]:=sumFun[nsamp]=1/nInt[nsamp]*{Total[sampledFunc[nsamp][[1]],Method->"CompensatedSummation"],Total[sampledFunc[nsamp][[2]],Method->"CompensatedSummation"],Total[sampledFunc[nsamp][[3]],Method->"CompensatedSummation"],Total[sampledFunc[nsamp][[4]],Method->"CompensatedSummation"]};
	
	relerr[0]=Abs[1-sumFun[1][[1]]/sumFun[0][[1]]];
	relerr[nsamp_]:=relerr[nsamp]=Abs[1-sumFun[nsamp][[1]]/sumFun[nsamp-1][[1]]];
	
	i=1;
	While[relerr[i]>10^(-tol)&&i<=10,i++];
	
	CPlus0=sumFun[i][[1]];
	CMin0=sumFun[i][[2]];
	CPlus1=sumFun[i][[3]];
	CMin1=sumFun[i][[4]];
	
	Remove[RInsol0,dRInsol0,RUpsol0,dRUpsol0,RInsol1,dRInsol1,RUpsol1,dRUpsol1,\[Lambda],\[Lambda]1,SWSH,SWSHS,sampledFunc,wrlist,sumFun];
	
  Association[
    "gg"->{nInt[i],relerr[i]},
    "Amplitudes"->
    <|
      "I"->CPlus0,
      "H"->CMin0
    |>,
    "AmplitudesCorrection"->
    <|
      "I"->CPlus1,
      "H"->CMin1
    |>,
    "Fluxes"->
    <|
      "En0PA"-><|
        "I"->Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "H"->\[Alpha]*Abs[CMin0]^2/(4Pi*\[Omega]^2)
      |>,
      "Jz0PA"->
      <|
        "I"->Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "H"->\[Alpha]*Abs[CMin0]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection"->
    <|
      "En1PA"-><|
        "I"->(2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2),
        "H"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMin0]^2 + 2*Re[CMin1*Conjugate[CMin0]] - 2*Abs[CMin0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2)
      |>,
      "Jz1PA"->
      <|
        "I"->(2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3),
        "H"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMin0]^2 + 2*Re[CMin1*Conjugate[CMin0]] - 3*Abs[CMin0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3)
      |>
    |>
  ]
]


(* ::Section::Closed:: *)
(*Fluxes for a spinning particle - circular orbits*)


(* ::Subsection::Closed:: *)
(*Fixed frequency*)


WASPCircEqFF[l_,m_,a_,rp_]:=
Module[{h1,h2,En,Lz,Kc,En1,En1tot,Enrp1,Lz1,Lzrp1,Lz1tot,\[CapitalOmega]\[Phi],\[CapitalGamma],\[Omega],
    SWSH,\[Lambda],\[ScriptCapitalC]2,rplus,\[Kappa],\[Epsilon],\[Alpha],W,sumPlus0,sumMin0,sumPlus1,sumMin1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,dKTeu1,V,dV,RInsol0,dRInsol0,ddRInsol0,dddRInsol0,
    RUpsol0,dRUpsol0,ddRUpsol0,dddRUpsol0,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Cmnn,Cmnmb,Cmmbmb,rho,pi,alpha,mu,gamma,tau,
    Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,rp1,\[CapitalSigma]1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,CPlus0,CMin0,CPlus1,CMin1,A2g,A1g,A0g,A2par,A1par,A0par,B3par,B2par,B1par,Rp,\[CapitalGamma]1,\[CapitalGamma]1rp1,\[CapitalGamma]1tot},
    
   \[CapitalDelta]  = rp^2-2rp+a^2;
   Rp = (rp^2+a^2)En - a*Lz;
   (* Geodesic constants of motion *)
   En = (a*rp+rp^(3/2)(-2+rp))/Sqrt[2a*rp^(7/2)+rp^4(-3+rp)]; 
   Lz = (rp(a^2-2a*Sqrt[rp]+rp^2))/Sqrt[2a*rp^(7/2)+rp^4(-3+rp)];
   Kc = (Lz-a*En)^2;
   (* Linear corrections to the constants of motion. These terms only appears in the projections to the Kinnersley tetrad of the velocity vn1 and vmb1*)
   En1 = ((a-Sqrt[rp])rp^(5/2)(3a^2-4a*Sqrt[rp]+rp^2))/(2(2a*rp^(7/2)+rp^4(-3+rp))^(3/2)); 
   Lz1 = ((3a^4+(9-2rp)(2-rp)rp^3+2a^2rp(2+rp)+a^3Sqrt[rp](-7+3rp)+3a*rp^(5/2)(-7+3rp)))/(2(2a*rp+rp^(3/2)(-3+rp))Sqrt[2a*rp^(7/2)+rp^4(-3+rp)]);
   Enrp1 = (-3a^2+8a*Sqrt[rp]+(-6+rp)rp)/(2 (2a+(-3+rp)Sqrt[rp])^(3/2)rp^(7/4)); 
   Lzrp1 = -(3a^2-8a*Sqrt[rp]-(-6+rp)rp)(a+rp^(3/2))/(2(2a+(-3+rp)Sqrt[rp])^(3/2)rp^(7/4));
   En1tot = En1 + Enrp1*rp1;
   Lz1tot = Lz1 + Lzrp1*rp1;
   
   rp1 = (a-Sqrt[rp])/rp;
   \[CapitalSigma]1 = 2*rp*rp1; 
   
   \[CapitalOmega]\[Phi] =  1/(rp^(3/2)+a); (* Geodesic BL frequencies *)
   \[CapitalGamma]  =  -((2a*Lz*rp)/\[CapitalDelta])+En(-a^2+(a^2+rp^2)^2/\[CapitalDelta]); (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
   \[CapitalGamma]1  = -((a(a^2+3rp^2)Rp(a^2(2En^2-rp)+a*En*Lz(-4+rp)-(-2+rp)(Lz^2+rp^2)))/(rp^2((Lz-a*En)^2+rp^2)\[CapitalDelta]^2))+((a^2+rp^2)Rp(-En*Lz*rp^3-a^2En*Lz(4+rp)+a^3(-rp+En^2(2+rp))+a (2Lz^2+rp^2(2+(-1+En^2)rp))))/(rp^2((Lz-a*En)^2+rp^2)\[CapitalDelta]^2)+(rp^4+a^2*rp*(2+rp))/\[CapitalDelta](En1)-(2a*rp)/\[CapitalDelta](Lz1);
   \[CapitalGamma]1rp1 = (rp^4+a^2*rp*(2+rp))/\[CapitalDelta](Enrp1)-(2a*rp)/\[CapitalDelta](Lzrp1) + ((2a*Lz*rp(-2+2rp))/\[CapitalDelta]^2-(2a*Lz)/\[CapitalDelta]+En(-(((-2+2rp)(a^2+rp^2)^2)/\[CapitalDelta]^2)+(4rp(a^2+rp^2))/\[CapitalDelta]));
   \[CapitalGamma]1tot = \[CapitalGamma]1 + \[CapitalGamma]1rp1*rp1;
   \[Omega]  = m*\[CapitalOmega]\[Phi]; (* Geodesic frequency and the linear correction *)  
   (*These functions are needed for the corrections to the coordinate time and azimuthal velocities. These terms only appears in the projections to the Kinnersley tetrad of the velocity vn1 and vmb1*)    
   h1 = Sqrt[Kc]/rp^3;
   h2 = -Sign[Lz-a*En]1/rp^3(-a^2*En+a*Lz+En*rp^3);
   (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
   {\[Lambda],SWSH} = angpar0[-2,l,m,a,\[Omega]]; 
   {RInsol0,dRInsol0} = Through[TsolverIn[-2,l,m,a,\[Omega],\[Lambda],rp,Automatic][[1;;2]][rp]];
   {RUpsol0,dRUpsol0} = Through[TsolverUpFromTSidentity[-2,l,m,a,\[Omega],\[Lambda],rp,Automatic][[1;;2]][rp]];
   \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
   rplus = 1+Sqrt[1-a^2]; 
   \[Kappa] = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
   \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
   \[Alpha] = 256*(2*rplus)^5*\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
      
   KTeu  = (rp^2+a^2)*\[Omega]-a*m;
   d\[CapitalDelta] = 2*(rp-1);
   dKTeu = 2*rp*\[Omega];
   V  = -(KTeu^2 + 4I*(rp-1)*KTeu)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
   dV = -((2KTeu*dKTeu+4I*KTeu+4I*(rp-1)*dKTeu)*\[CapitalDelta]-(KTeu^2+4I*(rp-1)*KTeu)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
   ddRInsol0  = (V*RInsol0+d\[CapitalDelta]*dRInsol0)/\[CapitalDelta]; (* second derivative of radial function from Teukolsky equation *)
   dddRInsol0 = 1/\[CapitalDelta](dV*RInsol0+(V+2)*dRInsol0);  
   ddRUpsol0  = (V*RUpsol0+d\[CapitalDelta]*dRUpsol0)/\[CapitalDelta];  
   dddRUpsol0 = 1/\[CapitalDelta](dV*RUpsol0+(V+2)*dRUpsol0);
       
   S = SWSH[\[Pi]/2][[1]];  (*  Spin-weighted spheroidal harmonics S(\[Theta](z))  *)
   dSd\[Theta] = SWSH[\[Pi]/2][[2]]; (* First derivative of S wrt \[Theta] *)
   d2Sd\[Theta]2 = S((a*\[Omega])^2-2(a*\[Omega])*m+m^2+2-\[Lambda]); (* second derivative of S from Teukolsky equation *)
   L2S = dSd\[Theta]-S(m-a*\[Omega]);(* Operators acting on S(\[Theta]) *)
   L1L2S = d2Sd\[Theta]2-2dSd\[Theta]*(m-a \[Omega])+S(-2+(m-a*\[Omega])^2);  

   \[Zeta] = rp;
   \[Zeta]bar = rp;
   \[CapitalSigma] = rp^2;
   {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,S,L2S,L1L2S,a];
   {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,S,L2S,L1L2S,a,\[Omega]];
   (* Four-velocity in Kinnersley tetrad *)
   vn = -((rp^2+a^2)*En - a*Lz)/(2*\[CapitalSigma]); 
   vl = -((rp^2+a^2)*En - a*Lz)/(\[CapitalDelta]);
   vm = (I*(a*En - Lz))/(-Sqrt[2]*\[Zeta]bar);
   vmb = Conjugate[vm];
   Cmnn = vn^2;
   Cmnmb = vn*vmb;
   Cmmbmb = vmb^2;
   (* Spin tensor in Kinnersley tetrad *)
   Sln  = -(rp*Sqrt[Kc])/(\[CapitalSigma]); 
   Snm  = (\[Zeta]/Sqrt[Kc])*vm*vn;
   Snmb = Conjugate[Snm];
   Slmb = -(\[Zeta]/Sqrt[Kc])*vl*vmb;
   (* Spin coefficients *)
   rho = 1/\[Zeta];
   pi = -(I*a)/(\[Zeta]^2Sqrt[2]);
   tau = (I*a)/(Sqrt[2]\[CapitalSigma]);
   mu = \[CapitalDelta]/(2\[Zeta]^2*\[Zeta]bar);
   gamma = (a^2-rp)/(2\[Zeta]^2*\[Zeta]bar);
   alpha = -(I*a)/(\[Zeta]^2Sqrt[2]);
   Scd\[Gamma]ndc = -Sln*2*Re[gamma];
   Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]+Slmb*mu;
   Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha])*vn-mu*vm)]);
   Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn));
   Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha])*vmb)
        -Snm*(alpha*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn))/2;
   St\[Phi]n  = -I*KTeu/(2\[CapitalSigma])*Sln+(a*\[Omega]-m)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
   St\[Phi]mb = -I*KTeu*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb);
   Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
   Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
   (* Linear parts of the four-velocity in Kinnersley tetrad *)          
   vn1  = -((rp^2+a^2)*(En1tot-h1) - a*(Lz1tot+h2))/(2*\[CapitalSigma]) - (2*rp*rp1*En-((rp^2+a^2)*En-a*Lz)/\[CapitalSigma]*\[CapitalSigma]1)/(2*\[CapitalSigma]);
   vmb1 = I (-(a*En - Lz)*(1/\[Zeta])*rp1+(a*(En1tot-h1) - (Lz1tot+h2)))/(Sqrt[2]*\[Zeta]);
                            
   Ann0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn + 2*vn1)*vn + Cdnn;
   Annt\[Phi]S = St\[Phi]n*vn;
   AnnrS  = (Srn + rp1*vn)*vn;
   Anmb0S  = (\[CapitalSigma]1/\[CapitalSigma])*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb;
   Anmbt\[Phi]S = (St\[Phi]n*vmb + St\[Phi]mb*vn)/2;
   AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
   Ambmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vmb + 2*vmb1)*vmb + Cdmbmb;
   Ambmbt\[Phi]S = St\[Phi]mb*vmb;
   AmbmbrS  = (Srmb + rp1*vmb)*vmb;
      
   A2g = Cmmbmb*fmbmb2;
   A1g = (Cmnmb*fnmb1+Cmmbmb*fmbmb1);
   A0g = (Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0);         
   (*These coefficients are different respect to the ones define in the companion paper. Here I put together everthing, including the differentials \[Delta]r, \[Delta]z*)
   B3par = -AmbmbrS*fmbmb2;
   B2par = -AnmbrS*fnmb1-AmbmbrS*fmbmb1;
   B1par = -AmbmbrS*fmbmb0-AnmbrS*fnmb0-AnnrS*fnn0;
      
   A2par = AmbmbrS*dfmbmb2dr+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb2;
   A1par = AmbmbrS*dfmbmb1dr+AnmbrS*dfnmb1dr+(Ambmb0S+Ambmbt\[Phi]S)fmbmb1+(Anmb0S+Anmbt\[Phi]S)fnmb1;
   A0par = AmbmbrS*dfmbmb0dr+AnmbrS*dfnmb0dr+AnnrS*dfnn0dr+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb0+(Anmb0S+Anmbt\[Phi]S)*fnmb0+(Ann0S+Annt\[Phi]S)*fnn0;
      
   sumPlus0  = \[CapitalSigma]*(A0g*RInsol0 - A1g*dRInsol0 + A2g*ddRInsol0);
   sumMin0 = \[CapitalSigma]*(A0g*RUpsol0 - A1g*dRUpsol0 + A2g*ddRUpsol0);   
   sumPlus1 = \[CapitalSigma]*(A0par*RInsol0 - (A1par+B1par)*dRInsol0 + (A2par+B2par)ddRInsol0 - B3par*dddRInsol0); 
   sumMin1 = \[CapitalSigma]*(A0par*RUpsol0 - (A1par+B1par)*dRUpsol0 + (A2par+B2par)ddRUpsol0 - B3par*dddRUpsol0);

   W = (RInsol0*dRUpsol0 - dRInsol0*RUpsol0)/\[CapitalDelta]; (* Invariant Wronskian *)
   CPlus0  = 2*Pi*sumPlus0/(\[CapitalGamma]*W); (* Geodesic amplitudes *)
   CMin0 = 2*Pi*sumMin0/(\[CapitalGamma]*W);
   CPlus1  = 2*Pi*(sumPlus1 - \[CapitalGamma]1tot/\[CapitalGamma]*sumPlus0)/(\[CapitalGamma]*W); (* Linear parts of the amplitudes *)
   CMin1 = 2*Pi*(sumMin1 - \[CapitalGamma]1tot/\[CapitalGamma]*sumMin0)/(\[CapitalGamma]*W);
    
    Association[
    "Amplitudes"->
    <|
      "I"->-CPlus0,
      "H"->-CMin0
    |>,
    "AmplitudesCorrection"->
    <|
      "I"->-CPlus1,
      "H"->-CMin1
    |>,
    "Fluxes"->
    <|
      "En0PA"-><|
        "I"->Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "H"->\[Alpha]*Abs[CMin0]^2/(4Pi*\[Omega]^2)
      |>,
      "Jz0PA"->
      <|
        "I"->Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "H"->\[Alpha]*Abs[CMin0]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection"->
    <|
      "En1PA"-><|
        "I"->(2*Re[CPlus1*Conjugate[CPlus0]])/(4Pi*\[Omega]^2),
        "H"->\[Alpha]*(2*Re[CMin1*Conjugate[CMin0]])/(4Pi*\[Omega]^2)
      |>,
      "Jz1PA"->
      <|
        "I"->(2*Re[CPlus1*Conjugate[CPlus0]])*m/(4Pi*\[Omega]^3),
        "H"->\[Alpha]*(2*Re[CMin1*Conjugate[CMin0]])*m/(4Pi*\[Omega]^3)
      |>
    |>
  ]
]


(* ::Subsection::Closed:: *)
(*Fixed eccentricity*)


WASPCircEqFE[l_,m_,a_,rp_]:=
Module[{s=-2,h1,h2,En,Lz,Kc,En1,Enrp1,En1tot,Lz1,Lzrp1,Lz1tot,\[CapitalOmega]\[Phi],\[CapitalOmega]\[Phi]1,\[CapitalOmega]\[Phi]1rp1,\[CapitalOmega]\[Phi]1tot,\[CapitalGamma],\[Omega],\[Omega]1,
    SWSH,SWSHS,\[Lambda],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,\[Kappa],\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumMin0,sumPlus1,sumMin1,\[CapitalDelta],d\[CapitalDelta],KTeu,KTeu1,dKTeu,dKTeu1,V,V1,dV,RInsol0,dRInsol0,ddRInsol0,dddRInsol0,
    RUpsol0,dRUpsol0,ddRUpsol0,RInsol1,dRInsol1,RUpsol1,dRUpsol1,dddRUpsol0,RIn0,RUp0,dRIn0,dRUp0,ddRIn0,dddRIn0,ddRUp0,dddRUp0,RIn1,RUp1,dRIn1,dRUp1,ddRIn1,ddRUp1,
    S,S1,L2S,L2S1,L1L2S,L1L2S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Cmnn,Cmnmb,Cmmbmb,rho,pi,alpha,mu,gamma,tau,
    Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,rp1,\[CapitalSigma]1,vn1,vn1rp1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,CPlus0,CMin0,CPlus1,CMin1,A2g,A1g,A0g,A2par,A1par,A0par,B3par,B2par,B1par,Rp,\[CapitalGamma]1,\[CapitalGamma]1rp1,\[CapitalGamma]1tot},
    
    \[CapitalDelta]  = rp^2-2rp+a^2;
    Rp = (rp^2+a^2)En - a*Lz;
    (* Geodesic constants of motion *)
    En = (a*rp+rp^(3/2)(-2+rp))/Sqrt[2a*rp^(7/2)+rp^4(-3+rp)]; 
    Lz = (rp(a^2-2a*Sqrt[rp]+rp^2))/Sqrt[2a*rp^(7/2)+rp^4(-3+rp)];
    Kc = (Lz-a*En)^2;
    
    rp1 = -((a (1-En^2)^2(a^2(1-En^2)+Lz^2)+8En*Kc*En1)/(2(4*Enrp1*En-3(1-En^2)^2)Kc));
    \[CapitalSigma]1 = 2*rp*rp1;     
    (* Linear corrections to the constants of motion. These terms only appears in the projections to the Kinnersley tetrad of the velocity vn1 and vmb1*)
    En1 = ((a-Sqrt[rp])rp^(5/2)(3a^2-4a*Sqrt[rp]+rp^2))/(2 (2a*rp^(7/2)+rp^4(-3+rp))^(3/2));
    Lz1 = ((3a^4+(9-2rp)(2-rp)rp^3+2a^2*rp(2+rp)+a^3*Sqrt[rp](-7+3rp)+3a*rp^(5/2)(-7+3rp)))/(2(2a*rp+rp^(3/2)(-3+rp))Sqrt[2a*rp^(7/2)+rp^4(-3+rp)]);
    Enrp1 = (-3a^2+8a*Sqrt[rp]+(-6+rp)rp)/(2 (2a+(-3+rp)Sqrt[rp])^(3/2)rp^(7/4)); 
    Lzrp1 = -(3a^2-8a*Sqrt[rp]-(-6+rp)rp)(a+rp^(3/2))/(2(2a+(-3+rp)Sqrt[rp])^(3/2)rp^(7/4));
    En1tot = En1 + Enrp1*rp1;
    Lz1tot = Lz1 + Lzrp1*rp1;
    
    \[CapitalOmega]\[Phi] = 1/(rp^(3/2)+a); (* Geodesic BL frequency *)
    \[CapitalOmega]\[Phi]1 = -3(-a+Sqrt[rp])/(2Sqrt[rp](a+rp^(3/2))^2); (* Linear correction to the BL frequency *)
    \[CapitalOmega]\[Phi]1rp1 = -3Sqrt[rp]/(2(a+rp^(3/2))^2);
    \[CapitalOmega]\[Phi]1tot = \[CapitalOmega]\[Phi]1 + \[CapitalOmega]\[Phi]1rp1*rp1;
    \[CapitalGamma]  =  -((2a*Lz*rp)/\[CapitalDelta])+En(-a^2+(a^2+rp^2)^2/\[CapitalDelta]); (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
    \[CapitalGamma]1 = (-Lz(a^2+rp^2)+a*En(a^2+3rp^2))/(rp*\[CapitalDelta])+(rp(2a^2+a^2*rp+rp^3)En1)/\[CapitalDelta]-(2a*rp*Lz1)/\[CapitalDelta];
    \[CapitalGamma]1rp1  = (rp^4+a^2*rp*(2+rp))/\[CapitalDelta](Enrp1)-(2a*rp)/\[CapitalDelta](Lzrp1) + ((2a*Lz*rp(-2+2rp))/\[CapitalDelta]^2-(2a*Lz)/\[CapitalDelta]+En(-(((-2+2rp)(a^2+rp^2)^2)/\[CapitalDelta]^2)+(4rp(a^2+rp^2))/\[CapitalDelta]));
    \[CapitalGamma]1tot = \[CapitalGamma]1 + \[CapitalGamma]1rp1*rp1;
    \[Omega] = m*\[CapitalOmega]\[Phi];
    \[Omega]1 = m*\[CapitalOmega]\[Phi]1tot;   
    (*These functions are needed for the corrections to the coordinate time and azimuthal velocities. 
     These terms only appears in the projections to the Kinnersley tetrad of the velocity vn1 and vmb1*)    
    h1 = Sqrt[Kc]/rp^3;
    h2 = -Sign[Lz-a*En]1/rp^3(-a^2*En+a*Lz+En*rp^3);
    (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
    {\[Lambda],\[Lambda]1,SWSH,SWSHS} = angpar[-2,l,m,a,\[Omega],\[Omega]1]; 
    {RInsol0,dRInsol0,RInsol1,dRInsol1} = TsolverIn1spin[-2,l,m,a,\[Omega],\[Omega]1,\[Lambda],\[Lambda]1,rp,Automatic][[1;;4]];
    {RUpsol0,dRUpsol0,RUpsol1,dRUpsol1} = TsolverUpFromTSidentity1spin[-2,l,m,a,\[Omega],\[Omega]1,\[Lambda],\[Lambda]1,rp,Automatic][[1;;4]];
    
    \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2);
    \[ScriptCapitalC]21 = 4\[Lambda]^3*\[Lambda]1+4\[Lambda]^2(3\[Lambda]1+10a(m-2a*\[Omega])\[Omega]1)+8\[Lambda](\[Lambda]1(1+10a*m*\[Omega]-10a^2 \[Omega]^2)+6a(m+2a*\[Omega])\[Omega]1)+48\[Omega](a*m*\[Lambda]1+6\[Omega]1-18a^3*m*\[Omega]*\[Omega]1+12a^4*\[Omega]^2*\[Omega]1+a^2(\[Lambda]1*\[Omega]+6m^2*\[Omega]1));
    rplus = 1+Sqrt[1-a^2]; 
    \[Kappa] = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
    \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
    (* constant for horizon fluxes and linear part*)
    \[Alpha] = 256*(2*rplus)^5*\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; 
    \[Alpha]1 = -256*(2*rplus)^5*\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*\[ScriptCapitalC]21+256*(2*rplus)^5*(\[Omega]1*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3+\[Kappa]*(2*\[Kappa]*\[Omega]1)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3+\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(2*\[Kappa]*\[Omega]1)*\[Omega]^3+\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*3*\[Omega]^2*\[Omega]1)/\[ScriptCapitalC]2;
    
    KTeu = (rp^2+a^2)*\[Omega]-a*m;
    KTeu1 = (rp^2+a^2)*\[Omega]1;
    d\[CapitalDelta] = 2*(rp-1);
    dKTeu = 2*rp*\[Omega];
    dKTeu1 = 2*rp*\[Omega]1;
    V = -(KTeu^2+4I*(rp-1)*KTeu)/\[CapitalDelta]+8*I*\[Omega]*rp+\[Lambda]; (* Potential in radial Teukolsky equation *)
    V1 = -(2*KTeu*KTeu1+4I*(rp-1)*KTeu1)/\[CapitalDelta]+8*I*\[Omega]1*rp+\[Lambda]1; (* Linear part of the potential in radial Teukolsky equation *)
    dV = -((2KTeu*dKTeu+4I*KTeu+4I*(rp-1)*dKTeu)*\[CapitalDelta]-(KTeu^2+4I*(rp-1)*KTeu)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
    RIn0 = RInsol0[rp]; (* radial function *)
    dRIn0 = dRInsol0[rp];
    ddRIn0 = (V*RIn0+d\[CapitalDelta]*dRIn0)/\[CapitalDelta]; (* second derivative of radial function from Teukolsky equation *)
    RIn1 = RInsol1[rp]; (* Linear part of the radial function *)
    dRIn1 = dRInsol1[rp];
    ddRIn1 = (V*RIn1+V1*RIn0+d\[CapitalDelta]*dRIn1)/\[CapitalDelta];
    dddRIn0 = 1/\[CapitalDelta](dV*RIn0+(V+2)*dRIn0); (* third derivative of radial function from Teukolsky equation *)
    RUp0 = RUpsol0[rp];
    dRUp0 = dRUpsol0[rp];
    ddRUp0 = (V*RUp0+d\[CapitalDelta]*dRUp0)/\[CapitalDelta];  
    RUp1 = RUpsol1[rp];
    dRUp1 = dRUpsol1[rp];
    ddRUp1 = (V*RUp1+V1*RUp0+d\[CapitalDelta]*dRUp1)/\[CapitalDelta];
    dddRUp0 = 1/\[CapitalDelta](dV*RUp0+(V+2)*dRUp0);
    (* Spin-weighted spheroidal harmonics S(\[Theta](z))*)
    S = SWSH[\[Pi]/2][[1]]; 
    S1 = SWSHS[\[Pi]/2][[1]]; (* Linear part of S(\[Theta](z)) *)
    dSd\[Theta] = SWSH[\[Pi]/2][[2]]; (* First derivative of S wrt \[Theta] *)
    dS1d\[Theta] = SWSHS[\[Pi]/2][[2]]; (* Linear part of first derivative of S wrt \[Theta] *)
    d2Sd\[Theta]2 = S((a*\[Omega])^2-2(a*\[Omega])*m+m^2+2-\[Lambda]); (* second derivative of S from Teukolsky equation *)
    d2S1d\[Theta]2 = -(-a^2*\[Omega]^2-m^2-2+2m*a*\[Omega]+\[Lambda])S1-(2a*\[Omega]1(-a*\[Omega]+m)+\[Lambda]1)S; (* Linear part of second derivative of S from Teukolsky equation *)
    L2S = dSd\[Theta]-S(m-a*\[Omega]); (* Operators acting on S(\[Theta]) *)
    L2S1 = dS1d\[Theta]-S1(m-a*\[Omega])+S(a*\[Omega]1);(* Operators acting on S(\[Theta]) *)
    L1L2S = d2Sd\[Theta]2+(dSd\[Theta](-2m+2a*\[Omega]))+S(-2+m^2-2a*m*\[Omega]+a^2*\[Omega]^2);  
    L1L2S1 = d2S1d\[Theta]2+(dS1d\[Theta](-2m+2a*\[Omega])+dSd\[Theta](2a*\[Omega]1)) + S1(-2+m^2-2a*m*\[Omega]+a^2*\[Omega]^2) + S(-2a*m*\[Omega]1 + a^2*2*\[Omega]*\[Omega]1);
    
    \[Zeta] = rp;
    \[Zeta]bar = rp;
    \[CapitalSigma] = rp^2;
    {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[rp,rp,1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,S,L2S,L1L2S,a];
    {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,S,L2S,L1L2S,a,\[Omega]];
    {fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21} = dfabidS[\[Zeta],\[Zeta]bar,1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,KTeu1,dKTeu1,S,L2S,L1L2S,S1,L2S1,L1L2S1,a];
    vn = -((rp^2+a^2)*En - a*Lz)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
    vl = -((rp^2+a^2)*En - a*Lz)/(\[CapitalDelta]);
    vm = (I*(a*En - Lz))/(-Sqrt[2]*\[Zeta]bar);
    vmb = Conjugate[vm];
    Cmnn   = vn^2;
    Cmnmb  = vn*vmb;
    Cmmbmb = vmb^2;
    (* Spin tensor in Kinnersley tetrad *)
    Sln  = -(rp*Sqrt[Kc])/(\[CapitalSigma]);
    Snm  = (\[Zeta]/Sqrt[Kc])*vm*vn;
    Snmb = Conjugate[Snm];
    Slmb = -(\[Zeta]/Sqrt[Kc])*vl*vmb;
    (* Spin coefficients *)
    rho = 1/\[Zeta]; 
    pi = -(I*a)/(\[Zeta]^2Sqrt[2]);
    tau = (I*a)/(Sqrt[2]\[CapitalSigma]);
    mu = \[CapitalDelta]/(2\[Zeta]^2*\[Zeta]bar);
    gamma = (a^2-rp)/(2\[Zeta]^2*\[Zeta]bar);
    alpha = -(I*a)/(\[Zeta]^2Sqrt[2]);
    Scd\[Gamma]ndc = -Sln*2*Re[gamma];
    Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]+Slmb*mu;
    Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha])*vn-mu*vm)]);
    Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn));
    Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha])*vmb)-Snm*(alpha*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn))/2;
    St\[Phi]n  = -I*KTeu/(2\[CapitalSigma])*Sln+(a*\[Omega]-m)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
    St\[Phi]mb = -I*KTeu*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb);
    Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
    Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
    (* Linear parts of the four-velocity in Kinnersley tetrad *)          
    vn1  = -((rp^2+a^2)*(En1tot-h1) - a*(Lz1tot+h2))/(2*\[CapitalSigma]) - (2*rp*rp1*En-((rp^2+a^2)*En-a*Lz)/\[CapitalSigma]*\[CapitalSigma]1)/(2*\[CapitalSigma]);
    vmb1 = I (-(a*En - Lz)*(1/\[Zeta])*rp1+(a*(En1tot-h1) - (Lz1tot+h2)))/(Sqrt[2]*\[Zeta]);
    
    Ann0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn + 2*vn1)*vn + Cdnn;
    Annt\[Phi]S = St\[Phi]n*vn;
    AnnrS  = (Srn + rp1*vn)*vn;
    Anmb0S  = (\[CapitalSigma]1/\[CapitalSigma])*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb;
    Anmbt\[Phi]S = (St\[Phi]n*vmb + St\[Phi]mb*vn)/2;
    AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
    Ambmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vmb + 2*vmb1)*vmb + Cdmbmb;
    Ambmbt\[Phi]S = St\[Phi]mb*vmb;
    AmbmbrS  = (Srmb + rp1*vmb)*vmb;

    A2g = Cmmbmb*fmbmb2;
    A1g = (Cmnmb*fnmb1+Cmmbmb*fmbmb1);
    A0g = (Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0);
    (*These coefficients are different respect to the ones define in the companion paper. Here I put together everthing, including the differentials \[Delta]r, \[Delta]z*)
    B3par = -AmbmbrS*fmbmb2;
    B2par = -AnmbrS*fnmb1-AmbmbrS*fmbmb1;
    B1par = -AmbmbrS*fmbmb0-AnmbrS*fnmb0-AnnrS*fnn0;
    
    A2par = AmbmbrS*dfmbmb2dr+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb2+Cmmbmb*fmbmb21;
    A1par = AmbmbrS*dfmbmb1dr+AnmbrS*dfnmb1dr+(Ambmb0S+Ambmbt\[Phi]S)fmbmb1+(Anmb0S+Anmbt\[Phi]S)fnmb1+(Cmmbmb*fmbmb11+Cmnmb*fnmb11); 
    A0par = AmbmbrS*dfmbmb0dr+AnmbrS*dfnmb0dr+AnnrS*dfnn0dr+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb0+(Anmb0S+Anmbt\[Phi]S)*fnmb0+(Ann0S+Annt\[Phi]S)*fnn0+(Cmmbmb*fmbmb01+Cmnmb*fnmb01+Cmnn*fnn01);
    
    sumPlus0 = \[CapitalSigma]*(A0g*RIn0 - A1g*dRIn0 + A2g*ddRIn0);
    sumMin0 = \[CapitalSigma]*(A0g*RUp0 - A1g*dRUp0 + A2g*ddRUp0);
    
    sumPlus1 = \[CapitalSigma]*(A0par*RIn0 - (A1par+B1par)*dRIn0 + (A2par+B2par)ddRIn0 - B3par*dddRIn0) + \[CapitalSigma]*(A0g*RIn1 - A1g*dRIn1 + A2g*ddRIn1);
    sumMin1 = \[CapitalSigma]*(A0par*RUp0 - (A1par+B1par)*dRUp0 + (A2par+B2par)ddRUp0 - B3par*dddRUp0) + \[CapitalSigma]*(A0g*RUp1 - A1g*dRUp1 + A2g*ddRUp1);      
    
    W = (RIn0*dRUp0 - dRIn0*RUp0)/\[CapitalDelta]; (* Invariant Wronskian *)
    W1 = (RIn1*dRUp0 + RIn0*dRUp1 - dRIn1*RUp0 -dRIn0*RUp1)/\[CapitalDelta]; (* Linear part of the invariant Wronskian *)
    
    CPlus0 = 2*Pi*sumPlus0/(\[CapitalGamma]*W); (* Geodesic amplitudes *)
    CMin0 = 2*Pi*sumMin0/(\[CapitalGamma]*W);
    CPlus1 = 2*Pi*(sumPlus1 - \[CapitalGamma]1tot/\[CapitalGamma]*sumPlus0 - W1/W*sumPlus0)/(\[CapitalGamma]*W); (* Linear parts of the amplitudes *)
    CMin1 = 2*Pi*(sumMin1 - \[CapitalGamma]1tot/\[CapitalGamma]*sumMin0 - W1/W*sumMin0)/(\[CapitalGamma]*W);
                         
	Association[
	"\[Delta]r"->rp1,
	"Amplitudes"->
	<|
     "I"->-CPlus0,
     "H"->-CMin0
    |>,
    "AmplitudesCorrection"->
    <|
    "I"->-CPlus1,
    "H"->-CMin1
    |>,
    "Fluxes"->
    <|
     "En0PA"-><|
       "I"->Abs[CPlus0]^2/(4Pi*\[Omega]^2),
       "H"->\[Alpha]*Abs[CMin0]^2/(4Pi*\[Omega]^2)
     |>,
     "Jz0PA"->
     <|
       "I"->Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
       "H"->\[Alpha]*Abs[CMin0]^2*m/(4Pi*\[Omega]^3)
     |>
    |>,
    "FluxesCorrection"->
    <|
    "En1PA"-><|
      "I"->(2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2),
      "H"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMin0]^2 + 2*Re[CMin1*Conjugate[CMin0]] - 2*Abs[CMin0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2)
    |>,
    "Jz1PA"->
    <|
     "I"->(2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3),
     "H"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMin0]^2 + 2*Re[CMin1*Conjugate[CMin0]] - 3*Abs[CMin0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3)
    |>
   |>
  ]
]


(* ::Subsection::Closed:: *)
(*Fixed turning points*)


WASPCircEqFT[l_,m_,a_,rp_]:=
Module[{s=-2,En,Lz,Kc,En1,Lz1,\[CapitalOmega]\[Phi],\[CapitalOmega]\[Phi]1,\[CapitalGamma],\[CapitalGamma]1,\[Omega],\[Omega]1,\[Lambda],\[Lambda]1,SWSH,SWSHS,RInsol0,dRInsol0,RUpsol0,dRUpsol0,RInsol1,dRInsol1,RUpsol1,dRUpsol1,
    \[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,\[Kappa],\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumMin0,sumPlus1,sumMin1,
    \[CapitalDelta],d\[CapitalDelta],KTeu,KTeu1,dKTeu,dKTeu1,V,V1,dV,RIn0,dRIn0,ddRIn0,RIn1,dRIn1,ddRIn1,dddRIn0,
    RUp0,dRUp0,ddRUp0,RUp1,dRUp1,ddRUp1,dddRUp0,S,S1,L2S,L2S1,L1L2S,L1L2S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,\[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    vl,vn,vm,vmb,Sln,Cmnn,Cmnmb,Cmmbmb,
    Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,Srn,Srmb,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Ambmb0S,AmbmbrS,CPlus0,CMin0,CPlus1,CMin1,\[CapitalDelta]tspar,\[Delta]rpar,\[Delta]Vrpar,\[CapitalDelta]\[Phi]spar,A2g,A1g,A0g,A2par,A1par,A0par,B3par,B2par,B1par},
    
    \[CapitalDelta] = rp^2-2rp+a^2;
    (* Geodesic constants of motion *)   
    En = (a*rp+rp^(3/2)(-2+rp))/Sqrt[2a*rp^(7/2)+rp^4(-3+rp)]; 
    Lz = (rp(a^2-2a*Sqrt[rp]+rp^2))/Sqrt[2a*rp^(7/2)+rp^4(-3+rp)];
    Kc = (Lz-a*En)^2;
    (* Linear corrections to the constants of motion. These terms only appears in the projections to the Kinnersley tetrad of the velocity vn1 and vmb1*)
    En1 = ((a-Sqrt[rp])rp^(5/2)(3a^2-4a*Sqrt[rp]+rp^2))/(2 (2a*rp^(7/2)+rp^4(-3+rp))^(3/2));
    Lz1 = ((3a^4+(9-2rp)(2-rp)rp^3+2a^2*rp(2+rp)+a^3*Sqrt[rp](-7+3rp)+3a*rp^(5/2)(-7+3rp)))/(2(2a*rp+rp^(3/2)(-3+rp))Sqrt[2a*rp^(7/2)+rp^4(-3+rp)]);
    
    \[CapitalOmega]\[Phi] = 1/(rp^(3/2)+a); (* Geodesic BL frequency *)
    \[CapitalOmega]\[Phi]1 = -3(-a+Sqrt[rp])/(2Sqrt[rp](a+rp^(3/2))^2); (* Linear correction to the BL frequency *)
    \[CapitalGamma] = -((2a*Lz*rp)/\[CapitalDelta])+En(-a^2+(a^2+rp^2)^2/\[CapitalDelta]); (* Geodesic average rate of change of BL time in Mino time *)
    \[CapitalGamma]1 = (-Lz(a^2+rp^2)+a*En(a^2+3rp^2))/(rp*\[CapitalDelta])+(rp(2a^2+a^2*rp+rp^3)En1)/\[CapitalDelta]-(2a*rp*Lz1)/\[CapitalDelta];
    \[Omega] = m*\[CapitalOmega]\[Phi];
    \[Omega]1 = m*\[CapitalOmega]\[Phi]1;
    (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
    {\[Lambda],\[Lambda]1,SWSH,SWSHS} = angpar[-2,l,m,a,\[Omega],\[Omega]1]; 
    {RInsol0,dRInsol0,RInsol1,dRInsol1} = TsolverIn1spin[-2,l,m,a,\[Omega],\[Omega]1,\[Lambda],\[Lambda]1,rp,Automatic][[1;;4]];
    {RUpsol0,dRUpsol0,RUpsol1,dRUpsol1} = TsolverUpFromTSidentity1spin[-2,l,m,a,\[Omega],\[Omega]1,\[Lambda],\[Lambda]1,rp,Automatic][[1;;4]];
    
    \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2);
    \[ScriptCapitalC]21 = 4\[Lambda]^3*\[Lambda]1+4\[Lambda]^2(3\[Lambda]1+10a(m-2a*\[Omega])\[Omega]1)+8\[Lambda](\[Lambda]1 (1+10a*m*\[Omega]-10a^2 \[Omega]^2)+6a(m+2a*\[Omega])\[Omega]1)+48\[Omega](a*m*\[Lambda]1+6\[Omega]1-18a^3*m*\[Omega]*\[Omega]1+12a^4*\[Omega]^2*\[Omega]1+a^2(\[Lambda]1*\[Omega]+6m^2*\[Omega]1));
    rplus = 1+Sqrt[1-a^2]; 
    \[Kappa] = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
    \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
    (* constant for horizon fluxes and linear part*)
    \[Alpha] = 256*(2*rplus)^5*\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; 
    \[Alpha]1 = -256*(2*rplus)^5*\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*\[ScriptCapitalC]21+256*(2*rplus)^5*(\[Omega]1*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3+\[Kappa]*(2*\[Kappa]*\[Omega]1)*(\[Kappa]^2+16*\[Epsilon]^2)*\[Omega]^3+\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(2*\[Kappa]*\[Omega]1)*\[Omega]^3+\[Kappa]*(\[Kappa]^2+4*\[Epsilon]^2)*(\[Kappa]^2+16*\[Epsilon]^2)*3*\[Omega]^2*\[Omega]1)/\[ScriptCapitalC]2;
    
    KTeu = (rp^2+a^2)*\[Omega]-a*m;
    KTeu1 = (rp^2+a^2)*\[Omega]1;
    d\[CapitalDelta] = 2*(rp-1);
    dKTeu = 2*rp*\[Omega];
    dKTeu1 = 2*rp*\[Omega]1;
    V = -(KTeu^2+4I*(rp-1)*KTeu)/\[CapitalDelta]+8*I*\[Omega]*rp+\[Lambda]; (* Potential in radial Teukolsky equation *)
    V1 = -(2*KTeu*KTeu1+4I*(rp-1)*KTeu1)/\[CapitalDelta]+8*I*\[Omega]1*rp+\[Lambda]1; (* Linear part of the potential in radial Teukolsky equation *)
    dV = -((2KTeu*dKTeu+4I*KTeu+4I*(rp-1)*dKTeu)*\[CapitalDelta]-(KTeu^2+4I*(rp-1)*KTeu)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
    RIn0 = RInsol0[rp]; (* radial function *)
    dRIn0 = dRInsol0[rp];
    ddRIn0 = (V*RIn0+d\[CapitalDelta]*dRIn0)/\[CapitalDelta]; (* second derivative of radial function from Teukolsky equation *)
    RIn1 = RInsol1[rp]; (* Linear part of the radial function *)
    dRIn1 = dRInsol1[rp];
    ddRIn1 = (V*RIn1+V1*RIn0+d\[CapitalDelta]*dRIn1)/\[CapitalDelta];
    dddRIn0 = 1/\[CapitalDelta](dV*RIn0+(V+2)*dRIn0); (* third derivative of radial function from Teukolsky equation *)
    RUp0 = RUpsol0[rp];
    dRUp0 = dRUpsol0[rp];
    ddRUp0 = (V*RUp0+d\[CapitalDelta]*dRUp0)/\[CapitalDelta];  
    RUp1 = RUpsol1[rp];
    dRUp1 = dRUpsol1[rp];
    ddRUp1 = (V*RUp1+V1*RUp0+d\[CapitalDelta]*dRUp1)/\[CapitalDelta];
    dddRUp0 = 1/\[CapitalDelta](dV*RUp0+(V+2)*dRUp0);
    (* Spin-weighted spheroidal harmonics S(\[Theta](z))*)
    S = SWSH[\[Pi]/2][[1]]; 
    S1 = SWSHS[\[Pi]/2][[1]]; (* Linear part of S(\[Theta](z)) *)
    dSd\[Theta] = SWSH[\[Pi]/2][[2]]; (* First derivative of S wrt \[Theta] *)
    dS1d\[Theta] = SWSHS[\[Pi]/2][[2]]; (* Linear part of first derivative of S wrt \[Theta] *)
    d2Sd\[Theta]2 = S((a*\[Omega])^2-2(a*\[Omega])*m+m^2+2-\[Lambda]); (* second derivative of S from Teukolsky equation *)
    d2S1d\[Theta]2 = -(-a^2*\[Omega]^2-m^2-2+2m*a*\[Omega]+\[Lambda])S1-(2a*\[Omega]1(-a*\[Omega]+m)+\[Lambda]1)S; (* Linear part of second derivative of S from Teukolsky equation *)
    L2S = dSd\[Theta]-S(m-a*\[Omega]); (* Operators acting on S(\[Theta]) *)
    L2S1 = dS1d\[Theta]-S1(m-a*\[Omega])+S(a*\[Omega]1);(* Operators acting on S(\[Theta]) *)
    L1L2S = d2Sd\[Theta]2+(dSd\[Theta] (-2m+2a*\[Omega]))+S(-2+m^2-2a*m*\[Omega]+a^2*\[Omega]^2);  
    L1L2S1 = d2S1d\[Theta]2+(dS1d\[Theta](-2m+2a*\[Omega])+dSd\[Theta](2a*\[Omega]1)) + S1(-2+m^2-2a*m*\[Omega]+a^2*\[Omega]^2) + S(-2 a*m*\[Omega]1 + a^2*2*\[Omega]*\[Omega]1);  
    
    \[Zeta] = rp;
    \[Zeta]bar = rp;
    \[CapitalSigma] = rp^2;  
    {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[rp,rp,1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,S,L2S,L1L2S,a];
    {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,S,L2S,L1L2S,a,\[Omega]];
    {fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21} = dfabidS[\[Zeta],\[Zeta]bar,1,\[CapitalDelta],d\[CapitalDelta],KTeu,dKTeu,KTeu1,dKTeu1,S,L2S,L1L2S,S1,L2S1,L1L2S1,a];
    (* Four-velocity in Kinnersley tetrad *)
    vn = -((rp^2+a^2)*En - a*Lz)/(2*\[CapitalSigma]); 
    vl = -((rp^2+a^2)*En - a*Lz)/(\[CapitalDelta]);
    vm = (I*(a*En - Lz))/(-Sqrt[2]*\[Zeta]bar);
    vmb = (I*(a*En - Lz))/(Sqrt[2]*\[Zeta]bar);
    Cmnn = vn^2;
    Cmnmb = vn*vmb;
    Cmmbmb = vmb^2;
    (* Spin tensor in Kinnersley tetrad *)
    Sln = -Sqrt[Kc]/rp;
    Cdnn = (Lz-a*En)(a^2*En-a*Lz+En*rp^2)(3a^3En-3a^2*Lz+Lz(4rp-\[CapitalDelta])-a(En*rp(4+rp)-En*\[CapitalDelta]))/(4Sqrt[Kc]*rp^6);
    Cdmbmb = (a En-Lz)(a((a^2*En-a*Lz+En*rp^2)^2)+(a*En-Lz)(a^2*En-a*Lz-En*rp^2)\[CapitalDelta])/(2Sqrt[Kc]*rp^4*\[CapitalDelta]);
    Cdnmb = I (4Kc(a*En-Lz)(a^2-rp)+3a*Kc(a^2*En-a*Lz+En*rp^2)+2(a*En-Lz)(a^2*En-a*Lz-En*rp^2)(a^2*En-a*Lz+En*rp^2)+(a*En-Lz)(a^2*En-a*Lz)(a^2*En-a*Lz+En*rp^2)+2(a*En-Lz)(a^2-rp)(a^2*En-a*Lz+En*rp^2)(a^2*En-a*Lz+En*rp^2)/\[CapitalDelta])/(8Sqrt[2]*Sqrt[Kc]rp^5);
    St\[Phi]n = I (Kc*KTeu+(Lz-a*En)(a^2*En-a*Lz+En*rp^2)(a*\[Omega]-m))/(2Sqrt[Kc]rp^3);
    Srn = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
    Srmb = I(a*En-Lz)(a^2*En-a*Lz+En*rp^2)/(Sqrt[2]*Sqrt[Kc]*rp^2);
    (* Linear parts of the four-velocity in Kinnersley tetrad *)
    vn1 = 1/(2Sqrt[Kc]*rp^3) (Lz^2+a^2(-En1*Sqrt[Kc]*rp+En^2(1+rp))-a(En*Lz(2+rp)+Sqrt[Kc](-Lz1*rp))-Sqrt[Kc](En1*rp^3));
    vmb1 = -I(En*Sqrt[Kc]*rp-(a*En-Lz)(a*En1-Lz1)rp)/(Sqrt[2](a*En-Lz)rp^2);
    
    Ann0S = (2*vn1)*vn + Cdnn;
    Annt\[Phi]S = (St\[Phi]n)*vn;
    AnnrS = (Srn)*vn;
    Anmb0S = (vn1*vmb + vn*vmb1 + Cdnmb);
    Anmbt\[Phi]S = ((St\[Phi]n*vmb)/2);
    AnmbrS = ((Srn*vmb + Srmb*vn)/2);
    Ambmb0S = (2*vmb1)*vmb + Cdmbmb;
    AmbmbrS = (Srmb)*vmb;
    
    A2g = Cmmbmb*fmbmb2;
    A1g = (Cmnmb*fnmb1+Cmmbmb*fmbmb1);
    A0g = (Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0);   
    (*These coefficients are different respect to the ones define in the companion paper. Here I put together everthing, including the differentials \[Delta]r, \[Delta]z*)
    B3par = -AmbmbrS*fmbmb2;
    B2par = -AnmbrS*fnmb1-AmbmbrS*fmbmb1;
    B1par = -AmbmbrS*fmbmb0-AnmbrS*fnmb0-AnnrS*fnn0;
    
    A2par = AmbmbrS*dfmbmb2dr+(Ambmb0S)*fmbmb2+Cmmbmb*fmbmb21;
    A1par = AmbmbrS*dfmbmb1dr+AnmbrS*dfnmb1dr+(Ambmb0S)fmbmb1+(Anmb0S+Anmbt\[Phi]S)fnmb1+(Cmmbmb*fmbmb11+Cmnmb*fnmb11);
    A0par = AmbmbrS*dfmbmb0dr+AnmbrS*dfnmb0dr+AnnrS*dfnn0dr+(Ambmb0S)*fmbmb0+(Anmb0S+Anmbt\[Phi]S)*fnmb0+(Ann0S+Annt\[Phi]S)*fnn0+(Cmmbmb*fmbmb01+Cmnmb*fnmb01+Cmnn*fnn01);
    
    sumPlus0 = \[CapitalSigma]*(A0g*RIn0 - A1g*dRIn0 + A2g*ddRIn0); (* Total of all quadrants *)
    sumMin0 = \[CapitalSigma]*(A0g*RUp0 - A1g*dRUp0 + A2g*ddRUp0);
    
    sumPlus1 = \[CapitalSigma]*(A0par*RIn0 - (A1par+B1par)*dRIn0 + (A2par+B2par)ddRIn0 - B3par*dddRIn0) + \[CapitalSigma]*(A0g*RIn1 - A1g*dRIn1 + A2g*ddRIn1);
    sumMin1 = \[CapitalSigma]*(A0par*RUp0 - (A1par+B1par)*dRUp0 + (A2par+B2par)ddRUp0 - B3par*dddRUp0) + \[CapitalSigma]*(A0g*RUp1 - A1g*dRUp1 + A2g*ddRUp1);      
    
    W = (RIn0*dRUp0 - dRIn0*RUp0)/\[CapitalDelta]; (* Invariant Wronskian *)
    W1 = (RIn1*dRUp0 + RIn0*dRUp1 - dRIn1*RUp0 -dRIn0*RUp1)/\[CapitalDelta]; (* Linear part of the invariant Wronskian *)
    
    CPlus0 = 2*Pi*sumPlus0/(\[CapitalGamma]*W); (* Geodesic amplitudes *)
    CMin0 = 2*Pi*sumMin0/(\[CapitalGamma]*W);
    CPlus1 = 2*Pi*(sumPlus1 - \[CapitalGamma]1/\[CapitalGamma]*sumPlus0 - W1/W*sumPlus0)/(\[CapitalGamma]*W); (* Linear parts of the amplitudes *)
    CMin1 = 2*Pi*(sumMin1 - \[CapitalGamma]1/\[CapitalGamma]*sumMin0 - W1/W*sumMin0)/(\[CapitalGamma]*W);
    
    Association[
    "Amplitudes"->
    <|
      "I"->-CPlus0,
      "H"->-CMin0
    |>,
    "AmplitudesCorrection"->
    <|
      "I"->-CPlus1,
      "H"->-CMin1
    |>,
    "Fluxes"->
    <|
      "En0PA"-><|
        "I"->Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "H"->\[Alpha]*Abs[CMin0]^2/(4Pi*\[Omega]^2)
      |>,
      "Jz0PA"->
      <|
        "I"->Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "H"->\[Alpha]*Abs[CMin0]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection"->
    <|
      "En1PA"-><|
        "I"->(2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2),
        "H"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMin0]^2 + 2*Re[CMin1*Conjugate[CMin0]] - 2*Abs[CMin0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2)
      |>,
      "Jz1PA"->
      <|
        "I"->(2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3),
        "H"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMin0]^2 + 2*Re[CMin1*Conjugate[CMin0]] - 3*Abs[CMin0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3)
      |>
    |>
  ]
]


(* ::Section::Closed:: *)
(*End package*)


End[];


EndPackage[];
