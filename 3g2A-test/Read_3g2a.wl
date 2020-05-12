(* ::Package:: *)

(* ::Item:: *)
(*Setting the directory and importing the amplitude expressions*)


<< "InitTwoLoopToolsFF.m"


SetDirectory["/home/francesco/git/myfiniteflowexamples/amplitudes/gggaa"];


<< "setupfiles/Setup_gggaa.m"


hels = {"+++++", "-++++", "+++-+", "--+++", "---++", "-++-+", "+++--", "+--+-"}


(* ::Item:: *)
(*Function to turn the expression into a format readable by FORM*)


readable[expr_]:= expr/. {ex[i_]:>ToExpression[StringJoin["x" , ToString[i]]], s[i_, j_]:>ToExpression[StringJoin["s", ToString[i], ToString[j]]]};
Unprotect[Power];
Format[Power[E, a_], CForm] := exp[a]
Format[Power[a_, 1/2], CForm] := sqrt[a]
Format[Power[a_, b_], CForm] := pow[a, b]
Protect[Power];


Do[
Get["integrand_1L_gggaa_trT123_qkGAGA_Ncp0_dsm2p0_"<>StringReplace[h,"-"->"m"]<>".m"];,
{h,hels}];


(* ::Subsection::Closed:: *)
(*+++++*)


Plus@@Integrand["gggaa","1L","trT[1,2,3]","Nc^0","fl[ClosedLoop1[qk,GA,GA]]","dsm2^0","+++++"] /. Rule[cc[n_,T_topo],coeff_]:>coeff*INT[n,Table[1,{ii,1,Length[T[[1]]]}],T];
ampMT= %/. INT[MU[1,1],{1,1,1,1,1},T_topo]:>0 /. INT[MU[1,1]^2,{1,1,1,1},T_topo]:>-1/6 //Factor;
% //readable


(* ::Subsection:: *)
(*+--+-*)


Plus@@Integrand["gggaa","1L","trT[1,2,3]","Nc^0","fl[ClosedLoop1[qk,GA,GA]]","dsm2^0","+--+-"] /. Rule[cc[n_,T_topo],coeff_]:>coeff*INT[n,Table[1,{ii,1,Length[T[[1]]]}],T];
% /. INT[MU[1,1],{1,1,1,1,1},T_topo]:>0 /. INT[MU[1,1]^2,{1,1,1,1},T_topo]:>-1/6 /. INT[MU[1,1],{1,1,1},T_topo]:>1/2;
% /. INT[1,{1,1},topo[{x_,y_}]]:>(1/eps-Log[s@@x]+2);
% /. INT[1,{1,1,1},topo[{___,{x1_,x2_,x3_},___}]]:>1/SimplifyKinematics[s[x1,x2,x3],PSanalytic]*(1/eps^2*SS[x1,x2,x3]);

% /. INT[1,{1,1,1},topo[{y1___,{x1_},y2___}]]:>INT[1,{1,1,1},topo[{{x1},y2,y1}]];

% /. INT[1,{1,1,1},topo[{{x1_},y1_,y2_}]]:>1/SimplifyKinematics[s@@y1-s@@y2,PSanalytic]*(1/eps^2*(SS@@y1-SS@@y2));

% /. INT[1,{1,1,1,1},topo[{y1___,{x1_,x2_},y2___}]]:>INT[1,{1,1,1,1},topo[{y2,y1,{x1,x2}}]];

% /. INT[1,{1,1,1,1},topo[{{x1_},{x2_},{x3_},{x4_,x5_}}]]:>2/SimplifyKinematics[s[x1,x2]*s[x2,x3],PSanalytic]*(SS[x1,x2]+SS[x2,x3]-SS[x4,x5])/eps^2+FBOX1M[topo[{{x1},{x2},{x3},{x4,x5}}]];
% /. SS[x1_,x2_,x3_]:>SS@@Complement[{1,2,3,4,5},{x1,x2,x3}];
% /. s[x1_,x2_,x3_]:>s@@Complement[{1,2,3,4,5},{x1,x2,x3}];
% /. (f:(s|SS))[x__]:>f@@Sort[{x}];
tmp= Collect[%,{eps,SS[__],Log[_],INT[__],FBOX1M[__],FBOX2ME[__],FBOX2MH[__]},Factor];
firstRational= tmp /. {FBOX1M[__]->0, Log[__]->0};


coupleList= {{{1,2}, {3,4}},{{1,2}, {4,5}},{{1,3}, {2,4}},{{1,3}, {4,5}},{{1,5}, {3,4}},{{1,5}, {2,4}}};


(*function to write part of the expression in terms of logarithms of ratios*)
setUpSystem[v_List]:=Module[ {tmpMatr,len, intCoeffs1, intCoeffs2, res},
tmpMatr=Position[coupleList, v];
len=Length@tmpMatr;
intCoeffs1= Table[tmpMatr[[i,1]], {i,1,len}];
intCoeffs2=Table[tmpMatr[[i,2]], {i,1,len}];
res= Total[Table[-(-1)^intCoeffs2[[i]]*nameListc1[[intCoeffs1[[i]]]], {i,1,len}]];
Return[res];
];


tb1= Table[coupleList[[i,1]], {i,1, Length@coupleList}];
tb2= Table[coupleList[[i,2]], {i,1, Length@coupleList}];
systemCouples=DeleteDuplicates[Join[tb1,tb2]]
systemCoefficient=Table[ToExpression[StringReplace[ToString[systemCouples[[i]]], { "{"->"c[", "}"->"]"}]] , {i, 1, Length@systemCouples}] 
nameListc2=Table[ToExpression[StringJoin["cB", StringJoin[ToString/@(coupleList[[i]] //Flatten)]]] , {i,1,Length[coupleList]}];
nameListc1=Table[ToExpression[StringJoin["cA", StringJoin[ToString/@(coupleList[[i]]  //Flatten) ]]] , {i,1,Length[coupleList]}];
resVarList= Table[mulSol[SimplifyKinematics[s@@coupleList[[i,1]]-s@@coupleList[[i,2]], PSanalytic]], {i,1,Length@coupleList} ];
varsList=Table[resVarList[[i,1]], {i,1,Length@coupleList}]
resList=Table[resVarList[[i,2]], {i,1,Length@coupleList}]


modRule[Rule[a_, b_]]:={Rule[a,b], Rule[-a,-b]};


shortSys=Drop[systemCoefficient, -1];
vars=Drop[nameListc1,-1];
systemList= Table[setUpSystem[systemCouples[[i]]], {i, 1, Length@systemCoefficient-1}]
Solve[shortSys==systemList,  vars]/. nameListc1[[-1]]->0
l1SubsList= % /. c[i_,  j_]:>Coefficient[tmp, Log[s[i,j]]] //Flatten //Simplify;
l1Coeffs=Table[vars[[i]]/.l1SubsList, {i, 1, Length@vars}] //Simplify;

dSlist={(-ex[2]-ex[2] ex[3]+ex[4]+ex[3] ex[4]+ex[2] ex[3] ex[5])->ex[2]/ex[1]*(s[3,4]-s[1,2]), (-1+ex[5])->(s[4,5]-s[1,2])/ex[1],(ex[2]+ex[2] ex[3]+ex[2]^2 ex[3]-ex[4]-ex[3] ex[4]-ex[2] ex[3] ex[4]-ex[2] ex[5])->ex[2]/ex[1]*(s[2,4]-s[1,3]), (1+ex[4])->(s[4,5]-s[1,3])/ex[1],(-ex[2] ex[3]+ex[4]+ex[2] ex[4]+ex[3] ex[4]+ex[2] ex[3] ex[5])->-ex[2]/ex[1]*(s[2,4]-s[1,5]),(ex[2] ex[3]+ex[2]^2 ex[3]-ex[4]-ex[3] ex[4]-ex[2] ex[3] ex[4])->-ex[2]/ex[1]*(s[3,4]-s[1,5])};
deltaSlist= Flatten[Table[modRule[dSlist[[i]]], {i,1,Length[dSlist]}]];


(* ::InheritFromParent:: *)
(**)


(* ::Subsubsection:: *)
(*Cancellation of Spurious Poles *)


c1234=Apart[l1Coeffs[[1]], ex[4]] /.modRule[(-ex[2]-ex[2] ex[3]+ex[4]+ex[3] ex[4]+ex[2] ex[3] ex[5])->ex[2]/ex[1]*(s[3,4]-s[1,2])] //Simplify;
cB1234=Coefficient[%, 1/(s[1,2]-s[3,4])^2];
l2Hat1234= SimplifyKinematics[1/s[3,4]^2*cB1234, PSanalytic]
l11234=SimplifyKinematics[(c1234-cB1234/(s[1,2]-s[3,4])^2)*(s[1,2]-s[3,4])/s[3,4] , PSanalytic] //Simplify;


c1245=Apart[l1Coeffs[[2]], ex[4]] /. modRule[(-1+ex[5])->(s[4,5]-s[1,2])/ex[1]] //Simplify;
cB1245=Coefficient[%, 1/(s[1,2]-s[4,5])^2];
l2Hat1245= SimplifyKinematics[1/s[4,5]^2*cB1245, PSanalytic]
l11245=SimplifyKinematics[(c1245-cB1245/(s[1,2]-s[4,5])^2)*(s[1,2]-s[4,5])/s[4,5] , PSanalytic] //Simplify;


c1324=Apart[l1Coeffs[[3]], ex[5]] /.modRule[(ex[2]+ex[2] ex[3]+ex[2]^2 ex[3]-ex[4]-ex[3] ex[4]-ex[2] ex[3] ex[4]-ex[2] ex[5])->ex[2]/ex[1]*(s[2,4]-s[1,3])] //Simplify;
cB1324= Coefficient[%, 1/(s[1,3]-s[2,4])^2];
l2Hat1324= SimplifyKinematics[1/s[2,4]^2*cB1324, PSanalytic];
l11324=SimplifyKinematics[(c1324-cB1324/(s[1,3]-s[2,4])^2)*(s[1,3]-s[2,4])/s[2,4] , PSanalytic] //Simplify;


c1345=Apart[l1Coeffs[[4]], ex[4]] /.modRule[(1+ex[4])->(s[4,5]-s[1,3])/ex[1]]  //Simplify;
cB1345= Coefficient[%, 1/(s[1,3]-s[4,5])^2];
l2Hat1345= SimplifyKinematics[1/s[4,5]^2*cB1345, PSanalytic];
l11345=SimplifyKinematics[(c1345-cB1345/(s[1,3]-s[4,5])^2)*(s[1,3]-s[4,5])/s[4,5] , PSanalytic] //Simplify;


c1534=Apart[l1Coeffs[[5]], ex[5]] /. modRule[(ex[2] ex[3]+ex[2]^2 ex[3]-ex[4]-ex[3] ex[4]-ex[2] ex[3] ex[4])->-ex[2]/ex[1]*(s[3,4]-s[1,5])]  //Simplify;
l11534= SimplifyKinematics[c1534*(s[1,3]-s[4,5])/s[4,5] , PSanalytic] //Simplify


c2List={l2Hat1234, l2Hat1245, l2Hat1324, l2Hat1345};
c1List= {l11234* L[1, s[1,2]/s[3,4]], l11245* L[1, s[1,2]/s[4,5]], l11324* L[1, s[1,3]/s[2,4]], l11345* L[1, s[1,3]/s[4,5]], l11534* L[1, s[1,5]/s[3,4]]};
l2HatVec={Lhat[2,s[1,2]/s[3,4]]+SimplifyKinematics[s[3,4]/(s[1,2]-s[3,4]), PSanalytic], Lhat[2,s[1,2]/s[4,5]]+SimplifyKinematics[s[4,5]/(s[1,2]-s[4,5]), PSanalytic],Lhat[2,s[1,3]/s[2,4]]+SimplifyKinematics[s[2,4]/(s[1,3]-s[2,4]), PSanalytic], Lhat[2,s[1,3]/s[4,5]]+SimplifyKinematics[s[4,5]/(s[1,3]-s[4,5]), PSanalytic]};
tmp /.{Log[__]->0};
final= %+c2List.l2HatVec+Total[c1List];


Export["+--+-.txt",final]
