#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"


using namespace std;

//----------------------------------------------------
//WE ASSUMED THAT IN REACTION A + B = 1 + 2 + 3
//A - VIRTUAL PHOTON
//B - INITIAL PROTON
//1 - PARTICLE WITH MASS m1
//2 - PARTICLE WITH MASS m2
//3 - PARTICLE WITH MASS m3

//We assume as global variables: MP - initial proton mass (is taken from pdg-table), E_beam - beam energy (is taken from input file).

//Here all derivations are performed via special root functions. If you need further level of explanation please see anti_rot_explanation, where all calculations are shown in more details and rotations are performed via matrices.

 void fermi_anti_rot(Float_t &theta_rot2,Float_t &phi_rot2, TLorentzVector &P4_1, TLorentzVector &P4_2, TLorentzVector &P4_3, TLorentzVector &P4_E_prime,TLorentzVector &P4_Eini_fermi) {
 
// cout << "qqqqq2 "<< theta_rot2<< "\n";
 
 
 Float_t theta_fermi, phi_fermi;
 
 
 theta_fermi = acos(pz_fermi/sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
 
  phi_fermi = acos(fabs(px_fermi)/sqrt(px_fermi*px_fermi+py_fermi*py_fermi));
 
 if ((px_fermi < 0.)&&(py_fermi > 0.)) phi_fermi = M_PI-phi_fermi;
 if ((px_fermi < 0.)&&(py_fermi < 0.)) phi_fermi = phi_fermi + M_PI;
 if ((px_fermi > 0.)&&(py_fermi < 0.)) phi_fermi = 2.*M_PI - phi_fermi;
 
 
 //------------------------------------------
 
  Float_t beta;

 beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi);
//  beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/MP;
 //--------------------------------------------
  P4_E_prime.RotateZ(-phi_rot2);
  P4_E_prime.RotateY(-theta_rot2);
  
  P4_1.RotateZ(-phi_rot2);
  P4_1.RotateY(-theta_rot2);
  
  P4_2.RotateZ(-phi_rot2);
  P4_2.RotateY(-theta_rot2);
  
  P4_3.RotateZ(-phi_rot2);
  P4_3.RotateY(-theta_rot2);
  
//  P4_Eini_fermi.RotateZ(-phi_rot2);
//  P4_Eini_fermi.RotateY(-theta_rot2);
  
  
 
  
  P4_E_prime.Boost(0,0,beta);
  P4_1.Boost(0,0,beta);
  P4_2.Boost(0,0,beta);
  P4_3.Boost(0,0,beta);
    
  //  P4_Eini_fermi.Boost(0,0,beta);
    
    P4_E_prime.RotateY(theta_fermi);
   P4_E_prime.RotateZ(phi_fermi);
    
  P4_1.RotateY(theta_fermi);
  P4_1.RotateZ(phi_fermi);
  

  P4_2.RotateY(theta_fermi);
  P4_2.RotateZ(phi_fermi);
  
  P4_3.RotateY(theta_fermi);
  P4_3.RotateZ(phi_fermi);
 
 //cout << "ttt \n";
 
 
// P4_Eini_fermi.RotateY(theta_fermi);
// P4_Eini_fermi.RotateZ(phi_fermi);
 
 
 
 
 };
 
