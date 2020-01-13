(* ::Package:: *)

(* ::Subsection:: *)
(*Initialising the 6-particle phase space point*)


minkDot[a_,b_]:={{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}}.a.b 
yRotation[a_]:={{Cos[a],0,Sin[a]},{0,1,0},{-Sin[a],0,Cos[a]}}
zRotation[a_]:={{Cos[a],-Sin[a],0},{Sin[a],Cos[a],0},{0,0,1}}
xRotation[a_]:= {{1,0,0}, {0, Cos[a], -Sin[a]},{0, Sin[a], Cos[a]}}


p1=(sqS/2)*{1,0,0,-1}
p2=(sqS/2)*{1,0,0,1}
nz={0,0,1}
n3=yRotation[\[Theta]1].nz
n4=zRotation[\[Theta]2].n3
n5=xRotation[\[Theta]3].n4
p3=-E3*Flatten[{1,n3}]
p4=-E4*Flatten[{1,n4}]
p5=-E5*lambda*Flatten[{1,n5}]
p6=-p1-p2-p3-p4-p5


(* ::Subsection:: *)
(*Solving for E4 to ensure masslessness of p6*)


sol=Solve[p6[[1]]^2-p6[[2]]^2-p6[[3]]^2-p6[[4]]^2==0, E4] [[1]]



(* ::Subsection:: *)
(*Conmputing numeric values of kinematic invariants*)


subs= {\[Theta]1->\[Pi]/7 , \[Theta]2->\[Pi]/3, \[Theta]3 ->\[Pi]/5, lambda->1, sqS->500, E3-> 100, E5->90,E4->\!\(\*
TagBox[
InterpretationBox["\"\<64.21552491862796\>\"",
64.21552491862796,
AutoDelete->True],
NumberForm[#, 16]& ]\)} 

minkDot[p1,p2]/.subs //N //Simplify
minkDot[p2,p3]/.subs //N //Simplify
minkDot[p3,p4]/.subs //N //Simplify
minkDot[p4,p5]/.subs //N //Simplify
minkDot[p5,p6]/.subs //N //Simplify
minkDot[p6,p6]/. subs //N //Simplify
minkDot[p5,p5]/. subs //N //Simplify
minkDot[p3,p6]+minkDot[p3,p4]+minkDot[p2,p3]+minkDot[p1,p3]+minkDot[p3,p5]/.subs


Unprotect[Power]
Format[Power[a_,b_], CForm]:=pow[a,b]
Protect[Power]

cExpr1=CForm[(E3^2+2 E3 E5 lambda+E5^2 lambda^2-2 E3 sqS-2 E5 lambda sqS+sqS^2-E3^2 Cos[\[Theta]1]^2-2 E3 E5 lambda Cos[\[Theta]1]^2 Cos[\[Theta]3]-E5^2 lambda^2 Cos[\[Theta]1]^2 Cos[\[Theta]3]^2-E3^2 Sin[\[Theta]1]^2-2 E3 E5 lambda Cos[\[Theta]2] Sin[\[Theta]1]^2-E5^2 lambda^2 Cos[\[Theta]2]^2 Sin[\[Theta]1]^2-E5^2 lambda^2 Cos[\[Theta]3]^2 Sin[\[Theta]1]^2 Sin[\[Theta]2]^2-2 E3 E5 lambda Cos[\[Theta]1] Sin[\[Theta]1] Sin[\[Theta]2] Sin[\[Theta]3]-E5^2 lambda^2 Cos[\[Theta]1]^2 Sin[\[Theta]3]^2-E5^2 lambda^2 Sin[\[Theta]1]^2 Sin[\[Theta]2]^2 Sin[\[Theta]3]^2)/(2 (-E3-E5 lambda+sqS+E3 Cos[\[Theta]1]^2+E5 lambda Cos[\[Theta]1]^2 Cos[\[Theta]3]+E3 Cos[\[Theta]2] Sin[\[Theta]1]^2+E5 lambda Cos[\[Theta]2]^2 Sin[\[Theta]1]^2+E5 lambda Cos[\[Theta]3] Sin[\[Theta]1]^2 Sin[\[Theta]2]^2))/.{Cos[\[Theta]1]->ct1, Cos[\[Theta]2]->ct2, Cos[\[Theta]3]->ct3, Sin[\[Theta]1]->st1, Sin[\[Theta]2]->st2, Sin[\[Theta]3]->st3, ct4->Cos[\[Theta]4], st4->Sin[\[Theta]4]}] 


CForm[-p3/E3]/.{Cos[\[Theta]1]->ct1, Cos[\[Theta]2]->ct2, Cos[\[Theta]3]->ct3, Sin[\[Theta]1]->st1, Sin[\[Theta]2]->st2, Sin[\[Theta]3]->st3, ct4->Cos[\[Theta]4], st4->Sin[\[Theta]4]}




(* ::InheritFromParent:: *)
(*List(1,st1,0,ct1)*)


CForm[-p4/E4]/.{Cos[\[Theta]1]->ct1, Cos[\[Theta]2]->ct2, Cos[\[Theta]3]->ct3, Sin[\[Theta]1]->st1, Sin[\[Theta]2]->st2, Sin[\[Theta]3]->st3, ct4->Cos[\[Theta]4], st4->Sin[\[Theta]4]}


CForm[-p5/(E5*lambda)]/.{Cos[\[Theta]1]->ct1, Cos[\[Theta]2]->ct2, Cos[\[Theta]3]->ct3, Sin[\[Theta]1]->st1, Sin[\[Theta]2]->st2, Sin[\[Theta]3]->st3, ct4->Cos[\[Theta]4], st4->Sin[\[Theta]4]}


(* ::Subsubsection:: *)
(*2 to 3*)


(*minkDot[a_,b_]:={{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}}.a.b
q=sqS {1,0,0,0}
p1=1/2 sqS {1,0,0,1}
p2=q-p1
yRotation[a_]:={{Cos[a],0,Sin[a]},{0,1,0},{-Sin[a],0,Cos[a]}}
zRotation[a_]:={{Cos[a],-Sin[a],0},{Sin[a],Cos[a],0},{0,0,1}}
nz={0,0,1}
n3=zRotation[\[Theta]1].yRotation[\[Theta]2].nz
n4=zRotation[\[Theta]1].yRotation[\[Theta]2].zRotation[\[Theta]3].yRotation[\[Theta]4].nz
p3= -E3*Flatten[{1,n3}]
p4=-E4*lambda*Flatten[{1,n4}]
p5=-p1-p2-p3-p4


(*Solve[p5[[1]]^2-p5[[2]]^2-p5[[3]]^2-p5[[4]]^2==0, E3] //Simplify


(*Unprotect[Power]
Format[Power[a_,b_], CForm]:=pow[a,b]
Protect[Power]
cExpr=CForm[(sqS (-2 E4 lambda+sqS))/(2 (-E4 lambda+sqS+E4 lambda Cos[\[Theta]4]))/.{Cos[\[Theta]1]->ct1,Cos[\[Theta]2]->ct2, Cos[\[Theta]3]->ct3,Cos[\[Theta]4]->ct4, Sin[\[Theta]1]->st1, Sin[\[Theta]2]->st2, Sin[\[Theta]3]->st3,  Sin[\[Theta]4]->st4}]








(* ::InheritFromParent:: *)
(**)


(*-((E3 lambda-lambda sqS+E3 lambda Cos[\[Theta]1]^2 Cos[\[Theta]4] Sin[\[Theta]2]^2+E3 lambda Cos[\[Theta]4] Sin[\[Theta]1]^2 Sin[\[Theta]2]^2+E3 lambda Cos[\[Theta]1]^2 Cos[\[Theta]2] Cos[\[Theta]3] Sin[\[Theta]2] Sin[\[Theta]4]+E3 lambda Cos[\[Theta]2] Cos[\[Theta]3] Sin[\[Theta]1]^2 Sin[\[Theta]2] Sin[\[Theta]4]+\[Sqrt](lambda^2 ((E3-sqS+E3 Cos[\[Theta]4] Sin[\[Theta]2]^2+E3 Cos[\[Theta]2] Cos[\[Theta]3] Sin[\[Theta]2] Sin[\[Theta]4])^2+1/2 (-3 E3^2+2 List^2+4 E3 sqS-2 sqS^2+E3^2 Cos[2 \[Theta]2]) (1+Cos[\[Theta]4]^2 Sin[\[Theta]1]^2 Sin[\[Theta]2]^2+2 Cos[\[Theta]2] Cos[\[Theta]3] Cos[\[Theta]4] Sin[\[Theta]1]^2 Sin[\[Theta]2] Sin[\[Theta]4]+Cos[\[Theta]2]^2 Cos[\[Theta]3]^2 Sin[\[Theta]1]^2 Sin[\[Theta]4]^2+Sin[\[Theta]1]^2 Sin[\[Theta]3]^2 Sin[\[Theta]4]^2+Cos[\[Theta]1]^2 (Cos[\[Theta]4]^2 Sin[\[Theta]2]^2+2 Cos[\[Theta]2] Cos[\[Theta]3] Cos[\[Theta]4] Sin[\[Theta]2] Sin[\[Theta]4]+(Cos[\[Theta]2]^2 Cos[\[Theta]3]^2+Sin[\[Theta]3]^2) Sin[\[Theta]4]^2)))))/(lambda^2 (1+Cos[\[Theta]4]^2 Sin[\[Theta]1]^2 Sin[\[Theta]2]^2+2 Cos[\[Theta]2] Cos[\[Theta]3] Cos[\[Theta]4] Sin[\[Theta]1]^2 Sin[\[Theta]2] Sin[\[Theta]4]+Cos[\[Theta]2]^2 Cos[\[Theta]3]^2 Sin[\[Theta]1]^2 Sin[\[Theta]4]^2+Sin[\[Theta]1]^2 Sin[\[Theta]3]^2 Sin[\[Theta]4]^2+Cos[\[Theta]1]^2 (Cos[\[Theta]4]^2 Sin[\[Theta]2]^2+2 Cos[\[Theta]2] Cos[\[Theta]3] Cos[\[Theta]4] Sin[\[Theta]2] Sin[\[Theta]4]+(Cos[\[Theta]2]^2 Cos[\[Theta]3]^2+Sin[\[Theta]3]^2) Sin[\[Theta]4]^2))))




(* ::Subsection:: *)
(*For 4 particle-PS*)


p2New= E2*p2

n3New= xRotation[\[Theta]1].yRotation[\[Theta]2].nz
p3New= -E3*lambda*Flatten[{1,n3New}]
p4New=-p1-p2New-p3New


Solve[minkDot[p4New,p4New]==0,E2] //Simplify


(E3 lambda (1+Cos[\[Theta]1] Cos[\[Theta]2]))/(-E3 lambda+sqS+E3 lambda Cos[\[Theta]1] Cos[\[Theta]2])/.{\[Theta]1->\[Pi]/3, \[Theta]2->\[Pi]/5, lambda->1, sqS->1000, E3->400} //N


CForm[(E3 lambda (1+Cos[\[Theta]1] Cos[\[Theta]2]))/(-E3 lambda+sqS+E3 lambda Cos[\[Theta]1] Cos[\[Theta]2])]/.{Cos[\[Theta]1]->ct1, Cos[\[Theta]2]->ct2, Sin[\[Theta]1]->st1,Sin[\[Theta]2]->st2}
CForm[p3New]/.{Cos[\[Theta]1]->ct1, Cos[\[Theta]2]->ct2, Sin[\[Theta]1]->st1,Sin[\[Theta]2]->st2}
