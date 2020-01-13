#include <iostream>
#include <getopt.h>
#include <cstdio>
#include <math.h>

#include "tools/PhaseSpace.h"
#include "ngluon2/Initialize.h"
#include "ngluon2/Model.h"
#include "ngluon2/NParton2.h"
#include "photonPerms.h"
#include "softPS6.h"

using namespace std;

int* trial(int vec[]) {
int* a= &vec[3];
return a;

    };

  

//The code verifies the factorisation of the colour-dressed AAgggg amplitude in the soft limit for E5 ->0, the soft phase-space point is computed using softPS6.h
int main(int argc, char **argv){

//constructing the vector in colour space, we only need 2 ordered amplitudes because of the overcompleteness

std::vector<int> order1(1);
std::vector<int> order2(1);
std::vector<int> order3(1);

order1.push_back(1);
order1.push_back(2);
order1.push_back(3);
order1.push_back(4);
order1.push_back(5);

order2.push_back(1);
order2.push_back(2);
order2.push_back(4);
order2.push_back(5);
order2.push_back(3);

order3.push_back(1);
order3.push_back(2);
order3.push_back(5);
order3.push_back(3);
order3.push_back(4);

//template <typename T>

//initialising the colour-ordered amplitudes
const Flavour<dd_real> g= StandardModel::G();
const Flavour<dd_real> A= StandardModel::Ax(StandardModel::IL(), StandardModel::IL().C());


const Flavour<dd_real> flavs[]={A,A,g,g,g,g};
const int hels[]={-1,-1,+1,+1,+1,+1};

NParton2<dd_real> amp1;
NParton2<dd_real> amp2;
NParton2<dd_real> amp3;

const dd_real MuR= 91.188;

amp1.setProcess(6,flavs);
amp2.setProcess(6,flavs);
amp3.setProcess(6,flavs);

amp1.setMuR2(MuR*MuR);
amp2.setMuR2(MuR*MuR);
amp3.setMuR2(MuR*MuR);


amp1.setHelicity(hels);
amp2.setHelicity(hels);
amp3.setHelicity(hels);


amp1.setLoopInduced();
amp2.setLoopInduced();
amp3.setLoopInduced();

//constructing the soft phase-space point 

for(int k=0; k<10; ++k){
std::vector<MOM<dd_real> > mom;
mom= getSoftPoint6<dd_real>(pow(10,-k));


amp1.setMomenta(mom);
amp2.setMomenta(mom);
amp3.setMomenta(mom);



complex<dd_real> totAmp1;
complex<dd_real> totAmp2;
complex<dd_real> totAmp3;

// computing the sum of sub-amplitudes over all the photon permutations 
std::vector<std::vector<int> > matr1;
std::vector<std::vector<int> > matr2;
std::vector<std::vector<int> > matr3;

   matr1= photonPerms(order1);
   for(int j=0; j<20; ++j){
      amp1.setOrder(matr1[j]);
      totAmp1+= amp1.eval(1).get0();    
    }     

   matr2= photonPerms(order2);
   for(int j=0; j<20; ++j){
      amp2.setOrder(matr2[j]);
      totAmp2+= amp2.eval(1).get0();
    }

 matr3= photonPerms(order3);
   for(int j=0; j<20; ++j){
      amp3.setOrder(matr3[j]);
      totAmp3+= amp3.eval(1).get0();
    }


    complex<dd_real> ampVec[]={totAmp1, totAmp2, totAmp3};

    complex<dd_real> dotProd1;
    complex<dd_real> dotProd2;
    complex<dd_real> dotProd3;

//constructing the colour matrix
     const dd_real colourMatr[3][3]={
    { 92./3., -16./3., -16./3.},
    {-16./3., 92./3., -16./3.},
    {-16./3., -16./3., 92./3.}};

//computing the squared amplitude
    for(int i=0; i<3; ++i){
       for(int j=0; j<3; ++j){
           dotProd1+=conj(ampVec[j])*colourMatr[j][i]*ampVec[i];

           }
       } 


// perfomring the soft factorisation

const Flavour<dd_real> fact_flavs[]={A,A,g,g,g};

const int fact_hels[]={-1,-1,+1,+1,+1};

std::vector<int> orderSoft(1);

orderSoft.push_back(1);
orderSoft.push_back(2);
orderSoft.push_back(3);
orderSoft.push_back(4);

std::vector<std::vector<int> > permSoft;
permSoft= photonPerms(orderSoft);

std::vector<MOM<dd_real> > fact_mom(5);

fact_mom[0]=mom[0];
fact_mom[1]=mom[1];
fact_mom[2]=mom[2];
fact_mom[3]=mom[3];
fact_mom[4]=mom[5];

//construct the 3 eikonal functions

MOM<dd_real> p1=mom[0];
MOM<dd_real> p2=mom[1];
MOM<dd_real> p3=mom[2];
MOM<dd_real> p4=mom[3];
MOM<dd_real> p5=mom[4];
MOM<dd_real> p6=mom[5];

NParton2<dd_real> fact_amp;

dd_real eikonal1= 2*dot(p4,p6)/(dot(p4,p5)*dot(p5,p6));
dd_real eikonal2= 2*dot(p3,p6)/(dot(p3,p5)*dot(p5,p6));
dd_real eikonal3= 2*dot(p3,p4)/(dot(p3,p5)*dot(p5,p4));

const dd_real col_cc3=9;
const dd_real col_cc1=9;
const dd_real col_cc2=9;

fact_amp.setProcess(5,fact_flavs);
fact_amp.setMomenta(fact_mom);
fact_amp.setHelicity(fact_hels);
fact_amp.setLoopInduced();
fact_amp.setMuR2(MuR*MuR);

complex<dd_real> totAmpSoft;

for(int i=0; i<12; ++i){
  fact_amp.setOrder(permSoft[i]);
  totAmpSoft+= fact_amp.eval(1).get0();
  }

const dd_real prefact= col_cc1*eikonal1+ col_cc3*eikonal3+ col_cc2*eikonal2;

complex<dd_real> sqFact=prefact*totAmpSoft*conj(totAmpSoft);

//printing the ratio of the amplitude and the factorised expression, the ratio tends to 1 for E5->0
cout << sqFact/dotProd1<<endl;
}
return 1;
 }

