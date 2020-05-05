#include <cstdio>
#include <getopt.h>
#include <iostream>
#include <math.h>

#include "analytic/0q3gA-analytic.h"

// #include "tools/PhaseSpace.h"
#include "ngluon2/Model.h"

using namespace std;

int main(int argc, char** argv)
{

    //const Flavour<double> A= StandardModel::Ax(StandardModel::IL(), StandardModel::IL().C());
    //const Flavour<double> G= StandardModel::G();
    //Flavour<double> flavs[]={A,A,G,G,G};
    //
    //
    //
    //
    //const double scalefact= 1;
    //
    //const double s= 1;
    //const double t= 1;
    //const double M1= 0;
    //const double M2= 0;
    //const double M3= 0;
    //const double M4= 1.3;
    //const double m1= 0;
    //const double m2= 0;
    //const double m3= 0;
    //const double m4= 0;
    //const double mur= 1;

    Amp0q3gAA_a<double> analytic(StandardModel::Ax(StandardModel::IL(), StandardModel::IL().C()), 1., 1.);

    //NJetAnalytic<double> njetan(scalefact, 5,1);

    return 1;
}
