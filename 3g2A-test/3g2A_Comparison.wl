(* ::Package:: *)

Exit[]


<< "InitTwoLoopToolsFF.m"
SetDirectory["/home/francesco/git/myfiniteflowexamples/amplitudes/gggaa"];
<< "setupfiles/Setup_gggaa.m"


(* ::Item:: *)
(*Upload the phase space points for comparison with NJet*)


Get["/home/francesco/git/myfiniteflowexamples/InitTwoLoopToolsFF.m"];
Get[ "/home/francesco/git/myfiniteflowexamples/setupfiles/Setup_5g.m"];
Import["/home/francesco/Downloads/mom.txt"] ;
PStable=ToExpression[StringReplace[%, {"e+" :> "*^", "e-" :> "*^-"}]];

mom[i_]:=MOM4[PStable[[i]]]/. MOM4[{x__}]:>MOM4[x];
mom2[i_]:=PStable[[i]];
PSnumeric1=Table[mom[i], {i,1,5}];
PSnumeric2=Table[mom[i], {i,6,10}];
PSnumeric3=Table[mom[i], {i,11,15}];
PSnumeric4=Table[mom[i], {i,16,20}];
PSnumeric5=Table[mom[i], {i,21,25}];
PSnumeric2Perm=Table[mom[i], {i,26,30}];
EvalMTanalyticRules[PSnumeric1]={};
EvalMTanalyticRules[PSnumeric2]={};
EvalMTanalyticRules[PSnumeric3]={};
EvalMTanalyticRules[PSnumeric4]={};
EvalMTanalyticRules[PSnumeric5]={};
EvalMTanalyticRules[PSnumeric2Perm]={};


subs1= fromMT /.lpS->s /. spA[i_, j_]:>spA[p[i], p[j]] /. spB[i_, j_]:>spB[p[i], p[j]] /. Rule[ex[x_],a_]:>Rule[ex[x],SetPrecision[SimplifyKinematics[a, PSnumeric1],16]];
subs2= fromMT /.lpS->s /. spA[i_, j_]:>spA[p[i], p[j]] /. spB[i_, j_]:>spB[p[i], p[j]] /. Rule[ex[x_],a_]:>Rule[ex[x],SetPrecision[SimplifyKinematics[a, PSnumeric2],16]];
subs3= fromMT /.lpS->s /. spA[i_, j_]:>spA[p[i], p[j]] /. spB[i_, j_]:>spB[p[i], p[j]] /. Rule[ex[x_],a_]:>Rule[ex[x],SetPrecision[SimplifyKinematics[a, PSnumeric3],16]];
subs4= fromMT /.lpS->s /. spA[i_, j_]:>spA[p[i], p[j]] /. spB[i_, j_]:>spB[p[i], p[j]] /. Rule[ex[x_],a_]:>Rule[ex[x],SetPrecision[SimplifyKinematics[a, PSnumeric4],16]];
subs5= fromMT /.lpS->s /. spA[i_, j_]:>spA[p[i], p[j]] /. spB[i_, j_]:>spB[p[i], p[j]] /. Rule[ex[x_],a_]:>Rule[ex[x],SetPrecision[SimplifyKinematics[a, PSnumeric5],16]];
subs2Perm=fromMT /.lpS->s /. spA[i_, j_]:>spA[p[i], p[j]] /. spB[i_, j_]:>spB[p[i], p[j]] /. Rule[ex[x_],a_]:>Rule[ex[x],SetPrecision[SimplifyKinematics[a, PSnumeric2Perm],16]];


(* ::Item:: *)
(*Functions to handle the notation and output*)


ang[i_, j_]:= spA[p[i], p[j]]
sq[i_, j_]:= spB[p[i], p[j]]
exToX[expr_]:= expr /. ex[i_]:> ToExpression[StringJoin["x", ToString[i]]];
Unprotect[Power];
Format[Power[E, a_], CForm] := exp[a]
Format[Power[a_, 1/2], CForm] := sqrt[a]
Format[Power[a_, b_], CForm] := pow[a, b]
Protect[Power];


(* ::Item:: *)
(*Renorm scale= mZ*)


mU2=91.118^2;


(* ::Item:: *)
(*Defining the polylogs in the complex plane (fboxId2 reconstructs the whole box finite part, fboxId the "truncated" one)*)


MyLog[x_,y_]:=Log[Abs[x/y]] - I*Pi*(HeavisideTheta[-Re[x]]-HeavisideTheta[-Re[y]]);
MyDiLog[x_,y_]:=Re[PolyLog[2,1-x/y]] - I*HeavisideTheta[-Re[x/y]]*MyLog[y-x,y]*Im[MyLog[x,y]];
I4[s_,t_,m2_,mur2_]:=Module[{eps2,eps1,eps0},
eps2 = 2/s/t;
eps1 = 2/s/t*(MyLog[mur2,-s]+MyLog[mur2,-t]-MyLog[mur2,-m2]);
eps0 = 1/s/t*(MyLog[mur2,-s]^2+MyLog[mur2,-t]^2-MyLog[mur2,-m2]^2
- 2*MyDiLog[-m2,-s] - 2*MyDiLog[-m2,-t] - MyLog[-s,-t]^2 -Pi^2/2 +Pi^2/6 );
Return[{eps2,eps1,eps0}];
];
analyticSubs= {L[1,s_/t_]:>MyLog[-s,-t]/(1-s/t), Lhat[2, s_/t_]:>MyLog[-s,-t]/(1-s/t)^2+1/(1-s/t)};
I4hat[s_,t_,m2_,mur2_]:=Module[{eps2,eps1,eps0},
eps2 = 2/s/t;
eps1 = 2/s/t*(MyLog[mur2,-s]+MyLog[mur2,-t]-MyLog[mur2,-m2]);
eps0 = 1/s/t*(- 2*MyDiLog[-m2, -s] - 2*MyDiLog[-m2,-t] - MyLog[-s,-t]^2 -Pi^2/3);
Return[{eps2,eps1,eps0}];
];
fboxId=FBOX1M[topo[{{i_},{j_},{k_},{l_,m_}}]]->I4hat[s[i,j], s[j,k], s[l,m], mU2][[3]];
fboxId2=FBOX1M[topo[{{i_},{j_},{k_},{l_,m_}}]]->I4[s[i,j], s[j,k], s[l,m], mU2][[3]];


(* ::Subsubsection::Closed:: *)
(*All-plus*)


Plus@@Integrand["gggaa","1L","trT[1,2,3]","Nc^0","fl[ClosedLoop1[qk,GA,GA]]","dsm2^0","+++++"] /. Rule[cc[n_,T_topo],coeff_]:>coeff*INT[n,Table[1,{ii,1,Length[T[[1]]]}],T];
ampMT= %/. INT[MU[1,1],{1,1,1,1,1},T_topo]:>0 /. INT[MU[1,1]^2,{1,1,1,1},T_topo]:>-1/6 //Factor;
phaseExpr=s[2,4]*s[2,3]/(spA[p[1],p[2]]*spA[p[2],p[3]]*spA[p[3],p[4]]*spA[p[4],p[5]]*spA[p[5],p[1]]);

phaseMT= SimplifyKinematics[phaseExpr, PSanalytic];
phaseSP1= SimplifyKinematics[phaseExpr, PSnumeric1];
phaseSP2= SimplifyKinematics[phaseExpr, PSnumeric2];
phaseSP3= SimplifyKinematics[phaseExpr, PSnumeric3];
phaseSP4= SimplifyKinematics[phaseExpr, PSnumeric4];
phaseSP2Perm= SimplifyKinematics[phaseExpr, PSnumeric2Perm];

checkPS1=ampMT/phaseMT *phaseSP1/.subs1;
checkPS2=ampMT/phaseMT *phaseSP2/.subs2;
checkPS3=ampMT/phaseMT *phaseSP3/.subs3;
checkPS4=ampMT/phaseMT *phaseSP4/.subs4;

checkPS1*Conjugate[checkPS1]/27.295166546462983
checkPS2*Conjugate[checkPS2]/81836.910945227093
checkPS3*Conjugate[checkPS3]/249374.66273155983
checkPS4*Conjugate[checkPS4]/292.25884053205414


ampMT/phaseMT*phaseSP2 /.subs2


ampMT/phaseMT*phaseSP2Perm /.subs2Perm


(* ::Subsubsection:: *)
(*-++++*)


Plus@@Integrand["gggaa","1L","trT[1,2,3]","Nc^0","fl[ClosedLoop1[qk,GA,GA]]","dsm2^0","-++++"] /. Rule[cc[n_,T_topo],coeff_]:>coeff*INT[n,Table[1,{ii,1,Length[T[[1]]]}],T];
ampMT=% /. INT[MU[1,1],{1,1,1,1,1},T_topo]:>0 /. INT[MU[1,1]^2,{1,1,1,1},T_topo]:>-1/6 /. INT[MU[1,1],{1,1,1},T_topo]:>1/2;
phaseExpr=spB[p[2],p[5]]^2*s[1,2]/(spB[p[1],p[2]]*spA[p[2],p[3]]*spA[p[3],p[4]]*spA[p[4],p[5]]*spB[p[5],p[1]]);
phaseMT= SimplifyKinematics[phaseExpr, PSanalytic];
phaseSP1= SimplifyKinematics[phaseExpr, PSnumeric1]
phaseSP2= SimplifyKinematics[phaseExpr, PSnumeric2];
phaseSP3= SimplifyKinematics[phaseExpr, PSnumeric3];
phaseSP2Perm=SimplifyKinematics[phaseExpr, PSnumeric2Perm];

checkPS1=ampMT/phaseMT*phaseSP1 /.subs1;
checkPS2=ampMT/phaseMT *phaseSP2/.subs2;
checkPS3=ampMT/phaseMT *phaseSP3/.subs3;

checkPS1* Conjugate[checkPS1]/ 4060.6938075559292
checkPS2*Conjugate[checkPS2]/115956.22515121225
checkPS3*Conjugate[checkPS3]/456556.80551549117


(* ::Subsubsection:: *)
(*+++-+*)


ampMT=(64 ex[1]^2 ex[3] (ex[2]-ex[4]+ex[5]+ex[3] ex[5]+ex[2] ex[3] ex[5]))/((1+ex[3]) (1+(1+ex[2]) ex[3]));


phaseExpr=spB[p[5],p[3]]^2/(spB[p[4],p[5]]*spA[p[5],p[1]]*spA[p[1],p[2]]*spA[p[2],p[3]]*spB[p[3],p[4]]);
phaseMT=SimplifyKinematics[phaseExpr, PSanalytic]//Simplify;
phaseSP1= SimplifyKinematics[phaseExpr, PSnumeric1];
phaseSP2= SimplifyKinematics[phaseExpr, PSnumeric2];
phaseSP3= SimplifyKinematics[phaseExpr, PSnumeric3];
checkPS1=ampMT/phaseMT* phaseSP1 /.subs1;
checkPS2=ampMT/phaseMT* phaseSP2 /.subs2;
checkPS3=ampMT/phaseMT* phaseSP3 /.subs3;
checkPS1* Conjugate[checkPS1]/60.117660741213406
checkPS2* Conjugate[checkPS2]/108701.59807154817
checkPS3* Conjugate[checkPS3]/249375.34214479593


(* ::Subsubsection:: *)
(*--+++*)


(* ::Item:: *)
(*Defining the analytically continued polylogs*)


Import["--+++.txt"] //ToExpression;
ampMT=% /.analyticSubs;


phaseExpr=(spA[p[1],p[2]]^4)/(spA[p[1],p[2]]*spA[p[2],p[3]]*spA[p[3],p[4]]*spA[p[4],p[5]]*spA[p[5],p[1]]);

ampMT=ampMT /.fboxId  /. s[i_, j_]:>SimplifyKinematics[s[i, j], PSanalytic]; 


phaseMT=SimplifyKinematics[phaseExpr, PSanalytic];
phaseSP1=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric1],16];
phaseSP2=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric2],16];
phaseSP3=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric3],16];
phaseSP4=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric4],16];
phaseSP5=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric5],16];
phaseSP2Perm=SimplifyKinematics[phaseExpr, PSnumeric2Perm];


checkPS1= SetPrecision[ampMT/phaseMT*phaseSP1 /.subs1,16];
checkPS2= ampMT/phaseMT*phaseSP2 /.subs2;
checkPS3= ampMT/phaseMT*phaseSP3 /.subs3;
checkPS4= ampMT/phaseMT*phaseSP4 /.subs4;
checkPS5= ampMT/phaseMT*phaseSP5 /.subs5;

SetPrecision[checkPS1*Conjugate[checkPS1]/33974.020263668674,16]
SetPrecision[checkPS2*Conjugate[checkPS2]/1943228.7278599341, 16]
SetPrecision[checkPS3*Conjugate[checkPS3]/8460000.1664732806, 16]
SetPrecision[checkPS4*Conjugate[checkPS4]/25320.262861295505, 16]
SetPrecision[checkPS5*Conjugate[checkPS5]/1381219.4688844401,16]




(* ::Subsection:: *)
(*+--+-*)
(**)


Import["+--+-.txt"] //ToExpression;
ampMT=% /.analyticSubs;
ampMT2= Import["+--+-Rough.txt"] //ToExpression;
ampMT2= ampMT2 /.Log[s_]:>MyLog[-s,mU2];


phaseExpr=(spB[p[1],p[4]]^4)/(spB[p[1],p[2]]*spB[p[2],p[3]]*spB[p[3],p[4]]*spB[p[4],p[5]]*spB[p[5],p[1]]);

ampMT=ampMT /.fboxId  /. s[i_, j_]:>SimplifyKinematics[s[i, j], PSanalytic]; 


phaseMT=SimplifyKinematics[phaseExpr, PSanalytic];
phaseSP1=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric1],16];
phaseSP2=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric2],16];
phaseSP3=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric3],16];
phaseSP4=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric4],16];
phaseSP5=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric5],16];
phaseSP2Perm=SimplifyKinematics[phaseExpr, PSnumeric2Perm];


checkPS1=ampMT/phaseMT*phaseSP1 /.subs1;
checkPS2= ampMT/phaseMT*phaseSP2 /.subs2;
checkPS3= ampMT/phaseMT*phaseSP3 /.subs3;
checkPS4= ampMT/phaseMT*phaseSP4 /.subs4;
checkPS5= ampMT/phaseMT*phaseSP5 /.subs5;

SetPrecision[checkPS1*Conjugate[checkPS1]/8662.2880538179343,16]
SetPrecision[checkPS2*Conjugate[checkPS2]/1718996.4964392409, 16]
SetPrecision[checkPS3*Conjugate[checkPS3]/6040039.7726440653, 16]


(* ::Subsection:: *)
(*+++--*)


Import["+++--.txt"] //ToExpression;
ampMT=% /.analyticSubs;
ampMT2= am /.Log[s_]:>MyLog[-s,mU2];

phaseExpr=(spA[p[4],p[5]]^4)/(spA[p[1],p[2]]*spA[p[2],p[3]]*spA[p[3],p[4]]*spA[p[4],p[5]]*spA[p[5],p[1]]);

ampMT=ampMT /.fboxId  /. s[i_, j_]:>SimplifyKinematics[s[i, j], PSanalytic]; 
ampMT2=ampMT2 /.fboxId  /. s[i_, j_]:>SimplifyKinematics[s[i, j], PSanalytic]; 

phaseMT=SimplifyKinematics[phaseExpr, PSanalytic];
phaseSP1=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric1],16]; 
phaseSP2=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric2],16];
phaseSP3=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric3],16];
phaseSP4=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric4],16];
phaseSP5=SetPrecision[SimplifyKinematics[phaseExpr, PSnumeric5],16];
phaseSP2Perm=SimplifyKinematics[phaseExpr, PSnumeric2Perm];


checkPS1= SetPrecision[ampMT/phaseMT*phaseSP1 /.subs1,16];
checkPS2= ampMT/phaseMT*phaseSP2 /.subs2;
checkPS3= ampMT/phaseMT*phaseSP3 /.subs3;
checkPS4= ampMT/phaseMT*phaseSP4 /.subs4;
checkPS5= ampMT/phaseMT*phaseSP5 /.subs5;

SetPrecision[checkPS1*Conjugate[checkPS1]/1498.0167865122278,16]
SetPrecision[checkPS2*Conjugate[checkPS2]/1648910.9983746347, 16]
SetPrecision[checkPS3*Conjugate[checkPS3]/4364303.8581917509, 16]
SetPrecision[checkPS4*Conjugate[checkPS4]/261.19087827049003, 16]
