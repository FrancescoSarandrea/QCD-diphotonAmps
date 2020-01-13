#ifndef NJET_EXAMPLES_SOFTPS_H
#define NJET_EXAMPLES_SOFTPS_H


#include <iostream>
#include <getopt.h>
#include <cstdio>
#include <math.h>

#include "tools/PhaseSpace.h"
#include "ngluon2/Initialize.h"
#include "ngluon2/Model.h"
#include "ngluon2/NParton2.h"

using namespace std;

//Define function getSoftPoint6 which gives a phase-space point which depends on the energy scale lambda, so to hard-code the p5 single soft limit for AA ->gggg. The code uses the results from 6-particle_Phase_Space.wl

template <typename T>

std::vector<MOM<T> > getSoftPoint6(const T lambda)
{
  const T sqS= 5e2;
  std::vector<MOM<T> > Mom(6);
  
  const T E3= T(100.);
  const T E5= T(90.);
  const T theta1= T(M_PI/7.);
  const T theta2= T(M_PI/3.);
  const T theta3= T(M_PI/5.);
  const T ct1= cos(theta1);
  const T ct2= cos(theta2);
  const T ct3= cos(theta3);
  const T st1= sin(theta1);
  const T st2= sin(theta2);
  const T st3= sin(theta3);


  const T E4=(pow(E3,2) - pow(ct1,2)*pow(E3,2) + 2*E3*E5*lambda - 2*pow(ct1,2)*ct3*E3*E5*lambda + pow(E5,2)*pow(lambda,2) - pow(ct1,2)*pow(ct3,2)*pow(E5,2)*pow(lambda,2) - 2*E3*sqS - 2*E5*lambda*sqS + pow(sqS,2) - pow(E3,2)*pow(st1,2) - 2*ct2*E3*E5*lambda*pow(st1,2) - pow(ct2,2)*pow(E5,2)*pow(lambda,2)*pow(st1,2) - pow(ct3,2)*pow(E5,2)*pow(lambda,2)*pow(st1,2)*pow(st2,2) - 2*ct1*E3*E5*lambda*st1*st2*st3 - pow(ct1,2)*pow(E5,2)*pow(lambda,2)*pow(st3,2) - pow(E5,2)*pow(lambda,2)*pow(st1,2)*pow(st2,2)*pow(st3,2))/(2.*(-E3 + pow(ct1,2)*E3 - E5*lambda + pow(ct1,2)*ct3*E5*lambda + sqS + ct2*E3*pow(st1,2) + pow(ct2,2)*E5*lambda*pow(st1,2) + ct3*E5*lambda*pow(st1,2)*pow(st2,2)));

  Mom[0]=MOM<T>(sqS/2.,0,0,-sqS/2.);
  Mom[1]= MOM<T>(sqS/2.,0,0,sqS/2.);
  Mom[2]= -E3*MOM<T>(1,st1,0,ct1);
  Mom[3]= -E4*MOM<T>(1,ct2*st1,st1*st2,ct1);
  Mom[4]=-E5*lambda*MOM<T>(1,ct2*st1,ct3*st1*st2 - ct1*st3,ct1*ct3 + st1*st2*st3);
  Mom[5]= -Mom[0]-Mom[1]-Mom[2]-Mom[3]-Mom[4];

  

 
  return Mom;
}

#endif

