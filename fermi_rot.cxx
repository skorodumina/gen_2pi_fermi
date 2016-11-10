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

 void fermi_rot(Float_t &E_beam_fermi, Float_t &theta_rot2,Float_t &phi_rot2,Float_t E_beam,TLorentzVector  P4_E_prime, TLorentzVector  &P4_E_prime_boosted) {
 
 Float_t theta_fermi, phi_fermi;
 TLorentzVector P4_EL, P4_in_Prot, P4_gamma;
 
 P4_EL.SetXYZT(0.,0.,E_beam,E_beam);
 P4_in_Prot.SetXYZT(px_fermi,py_fermi,pz_fermi,sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
 
 P4_gamma = P4_EL -  P4_E_prime;
 
  Float_t Q2 = -P4_gamma.Mag2();


 theta_fermi = acos(pz_fermi/sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
 phi_fermi = acos(fabs(px_fermi)/sqrt(px_fermi*px_fermi+py_fermi*py_fermi));
 
 if ((px_fermi < 0.)&&(py_fermi > 0.)) phi_fermi = M_PI-phi_fermi;
 if ((px_fermi < 0.)&&(py_fermi < 0.)) phi_fermi = phi_fermi + M_PI;
 if ((px_fermi > 0.)&&(py_fermi < 0.)) phi_fermi = 2.*M_PI - phi_fermi;
 
Float_t eps_t1 = 1./(1.+Q2*((P4_gamma.Vect()).Mag2())/2./(((P4_EL.Vect()).Cross(P4_E_prime.Vect())).Mag2()));
Float_t eps_l1 = Q2*eps_t1/P4_gamma[3]/P4_gamma[3];

  P4_EL.RotateZ(-phi_fermi);
  P4_EL.RotateY(-theta_fermi);
  
  P4_gamma.RotateZ(-phi_fermi);
  P4_gamma.RotateY(-theta_fermi);
   
  P4_E_prime.RotateZ(-phi_fermi);
  P4_E_prime.RotateY(-theta_fermi);
   
  P4_in_Prot.RotateZ(-phi_fermi);
  P4_in_Prot.RotateY(-theta_fermi);
  

//cout << P4_EL[0] << "   "<<  P4_EL[1] << "   "<<  P4_EL[2] << "   "<< P4_EL[3] << " ioi\n";
//cout << R[0] << "   "<<  R[1] << "   "<< R[2] <<  "\n";

  
// cout << "rrr   "<< P4_in_Prot[0] << "   "<<  P4_in_Prot[1] << "   "<<  P4_in_Prot[2] << "   "<< P4_in_Prot[3]<<"   " <<sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi) << "\n";
 

//------------------------------------
 
 Float_t beta;
 
 beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi);
// beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/MP;
 P4_in_Prot.Boost(0,0,-beta);
 P4_EL.Boost(0,0,-beta);
 P4_E_prime.Boost(0,0,-beta);
 P4_gamma.Boost(0,0,-beta);
 
//cout << "rrr   "<< P4_in_Prot[0] << "   "<<  P4_in_Prot[1] << "   "<<  P4_in_Prot[2] << "   "<< P4_in_Prot[3]<<"   " <<sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi) << "\n";

//-----------------------

//ROT2


  theta_rot2 = acos(P4_EL[2]/sqrt(P4_EL[0]*P4_EL[0]+P4_EL[1]*P4_EL[1]+P4_EL[2]*P4_EL[2]));


  phi_rot2 = acos(fabs(P4_EL[0])/sqrt(P4_EL[0]*P4_EL[0]+P4_EL[1]*P4_EL[1]+P4_EL[2]*P4_EL[2]));
 
 if ((P4_EL[0] < 0.)&&(P4_EL[1] > 0.)) phi_rot2 = M_PI-phi_rot2;
 if ((P4_EL[0] < 0.)&&(P4_EL[1] < 0.)) phi_rot2 = phi_rot2 + M_PI;
 if ((P4_EL[0] > 0.)&&(P4_EL[1] < 0.)) phi_rot2 = 2.*M_PI - phi_rot2;


 P4_EL.RotateY(theta_rot2);
 P4_EL.RotateZ(phi_rot2);
 

 P4_gamma.RotateY(theta_rot2);
 P4_gamma.RotateZ(phi_rot2);
 
 E_beam_fermi = P4_EL[3];

//cout <<"rrrtrrrr                  "<< P4_EL[0] << "   "<<  P4_EL[1] << "   "<<  P4_EL[2] << "   "<< P4_EL[3] << "\n"; 
 
 P4_E_prime.RotateY(theta_rot2);
 P4_E_prime.RotateZ(phi_rot2);
 
  P4_in_Prot.RotateY(theta_rot2);
  P4_in_Prot.RotateZ(phi_rot2);
 //----------------------------
 P4_E_prime_boosted=P4_E_prime;
 
Float_t eps_t2 = 1./(1.+Q2*((P4_gamma.Vect()).Mag2())/2./(((P4_EL.Vect()).Cross(P4_E_prime.Vect())).Mag2()));
Float_t eps_l2 = Q2*eps_t2/P4_gamma[3]/P4_gamma[3];

//cout << eps_l1<< " "<< eps_l2<<" "<< eps_l1/eps_l2<<"  \n";
//cout << eps_l2<<"  f\n";

 /*
  Float_t W_f = (P4_in_Prot + P4_EL - P4_E_prime).Mag();
  Float_t Q2 = -P4_gamma.Mag2();
 
// cout << (W_f*W_f+Q2-MP*MP)/2./MP <<" "<< P4_gamma[3]<<"  \n";
TRotation rot;
TVector3 uz = P4_gamma.Vect().Unit();
 TVector3 ux = (P4_EL.Vect().Cross(P4_E_prime.Vect())).Unit();
 ux.Rotate(3.*M_PI/2,uz);
 rot.SetZAxis(uz,ux).Invert();

 P4_gamma.Transform(rot);


 
P4_gamma.Boost(0,0,-sqrt(P4_gamma[3]*P4_gamma[3]+Q2)/(P4_gamma[3]+MP));
*/
//cout << (W_f*W_f-Q2-MP*MP)/2./W_f << " "<< P4_gamma[3]<< "  cms\n";
 };
 
 void get_EpsL_Ebeam_ferm( Float_t E_beam,Float_t Ep, Float_t phi_e,Float_t theta_e, Float_t &eps_l,Float_t &eps_t, Float_t &E_beam_fermi,Float_t &E_p_fermi) { 
 
 Float_t theta_fermi, phi_fermi,theta_rot2,phi_rot2;
 Float_t nu, E_E_prime, Theta_e_prime;
 TLorentzVector P4_EL, P4_in_Prot, P4_gamma,P4_E_prime;
 
// cout << phi_e<< " "<< theta_e<< " "<< Ep<<" rsoft\n";
 
 P4_EL.SetXYZT(0.,0.,E_beam,E_beam);
 P4_in_Prot.SetXYZT(px_fermi,py_fermi,pz_fermi,sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
  
   
  
 P4_E_prime.SetXYZT(Ep*cos(phi_e)*sin(theta_e),Ep*sin(phi_e)*sin(theta_e),Ep*cos(theta_e),Ep);
// cout << P4_E_prime[0]<<" "<< P4_E_prime[1]<<" " << P4_E_prime[2]<<" "<< P4_E_prime[3]<<" "<< " rc\n";   
//    cout << P4_in_Prot[0]<<" "<< P4_in_Prot[1]<<" " << P4_in_Prot[2]<<" "<< P4_in_Prot[3]<<" "<< " rc\n"; 
      
  //  cout << (P4_in_Prot+ P4_EL-P4_E_prime).Mag()<<" tt\n";  
      
 P4_gamma = P4_EL -  P4_E_prime;
//if (isnan(P4_gamma[3])) cout << P4_gamma[3]<<" "<<P4_EL[3]<<" "<<P4_E_prime[3] <<" "<<Ep<<" \n";
//---------------------------
//cout << P4_E_prime[3]/P4_EL[3]<<" before\n";

 theta_fermi = acos(pz_fermi/sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
 phi_fermi = acos(fabs(px_fermi)/sqrt(px_fermi*px_fermi+py_fermi*py_fermi));
 
 if ((px_fermi < 0.)&&(py_fermi > 0.)) phi_fermi = M_PI-phi_fermi;
 if ((px_fermi < 0.)&&(py_fermi < 0.)) phi_fermi = phi_fermi + M_PI;
 if ((px_fermi > 0.)&&(py_fermi < 0.)) phi_fermi = 2.*M_PI - phi_fermi;
 


  P4_EL.RotateZ(-phi_fermi);
  P4_EL.RotateY(-theta_fermi);
  
  P4_gamma.RotateZ(-phi_fermi);
  P4_gamma.RotateY(-theta_fermi);
   
  P4_E_prime.RotateZ(-phi_fermi);
  P4_E_prime.RotateY(-theta_fermi);
   
  P4_in_Prot.RotateZ(-phi_fermi);
  P4_in_Prot.RotateY(-theta_fermi);
  


//------------------------------------
 
 Float_t beta;
 
 beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi);
// beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/MP;
 P4_in_Prot.Boost(0,0,-beta);
 P4_EL.Boost(0,0,-beta);
 P4_E_prime.Boost(0,0,-beta);
 P4_gamma.Boost(0,0,-beta);
 

//-----------------------

//ROT2


  theta_rot2 = acos(P4_EL[2]/sqrt(P4_EL[0]*P4_EL[0]+P4_EL[1]*P4_EL[1]+P4_EL[2]*P4_EL[2]));

  phi_rot2 = acos(fabs(P4_EL[0])/sqrt(P4_EL[0]*P4_EL[0]+P4_EL[1]*P4_EL[1]+P4_EL[2]*P4_EL[2]));
 
 if ((P4_EL[0] < 0.)&&(P4_EL[1] > 0.)) phi_rot2 = M_PI-phi_rot2;
 if ((P4_EL[0] < 0.)&&(P4_EL[1] < 0.)) phi_rot2 = phi_rot2 + M_PI;
 if ((P4_EL[0] > 0.)&&(P4_EL[1] < 0.)) phi_rot2 = 2.*M_PI - phi_rot2;


 P4_EL.RotateY(theta_rot2);
 P4_EL.RotateZ(phi_rot2);
 
 P4_gamma.RotateY(theta_rot2);
 P4_gamma.RotateZ(phi_rot2);
  
 P4_E_prime.RotateY(theta_rot2);
 P4_E_prime.RotateZ(phi_rot2);
 
 P4_in_Prot.RotateY(theta_rot2);
 P4_in_Prot.RotateZ(phi_rot2);
 
// cout <<  P4_in_Prot[0]<<" "<<P4_in_Prot[1]<<" "<<P4_in_Prot[2]<<" "<<P4_in_Prot[3]<<" "<<" pr\n";
// cout <<  P4_EL[0]<<" "<<P4_EL[1]<<" "<<P4_EL[2]<<" "<<P4_EL[3]<<" "<<" el\n"; 
// cout << P4_E_prime[3]/P4_EL[3]<<" after\n";
 //----------------------------

 E_beam_fermi = P4_EL[3];
 E_p_fermi = P4_E_prime[3];
 
 Float_t Q2 = -P4_gamma.Mag2();
 
 eps_t = 1./(1.+Q2*((P4_gamma.Vect()).Mag2())/2./(((P4_EL.Vect()).Cross(P4_E_prime.Vect())).Mag2()));
 eps_l = Q2*eps_t/P4_gamma[3]/P4_gamma[3];

// cout << theta_rot2<<" "<<phi_rot2<<" "<< P4_E_prime.Theta()<<" "<< P4_E_prime.Phi()<<" rrrrr\n";
// cout << -P4_gamma.Mag2()<<" "<<2.*E_beam_fermi*E_p_fermi*(1. - cos(P4_E_prime.Theta())) <<" q2\n";
// cout << P4_E_prime.Theta()<< " "<<P4_EL.Angle(P4_E_prime.Vect())<<" ang\n";
//cout << eps_l<<" "<<P4_gamma[3] <<" \n";
//cout << (P4_in_Prot+ P4_EL-P4_E_prime).Mag()<<" tt\n";
//cout <<P4_E_prime[3] <<"  ,,\n";
//cout << P4_gamma[3]<< " "<<(P4_in_Prot+ P4_EL-P4_E_prime).Mag() << " \n";
// cout << E_beam_fermi<<"\n";
 };
 
 
 
 //-------------------------------------
 
 
 
 void get_rot2( Float_t E_beam,Float_t Ep, Float_t phi_e,Float_t theta_e,Float_t &theta_rot2, Float_t &phi_rot2, Float_t &ph_e_ferm) { 
 
 Float_t theta_fermi, phi_fermi;
 Float_t nu, E_E_prime, Theta_e_prime;
 TLorentzVector P4_EL, P4_in_Prot, P4_gamma,P4_E_prime;
 
 P4_EL.SetXYZT(0.,0.,E_beam,E_beam);
 P4_in_Prot.SetXYZT(px_fermi,py_fermi,pz_fermi,sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
  if (isnan(theta_e)) cout << theta_e<<"  oo\n";  
   
  
 P4_E_prime.SetXYZT(Ep*cos(phi_e)*sin(theta_e),Ep*sin(phi_e)*sin(theta_e),Ep*cos(theta_e),Ep);
// cout << P4_E_prime[0]<<" "<< P4_E_prime[1]<<" " << P4_E_prime[2]<<" "<< P4_E_prime[3]<<" "<< " rc\n";   
//    cout << P4_in_Prot[0]<<" "<< P4_in_Prot[1]<<" " << P4_in_Prot[2]<<" "<< P4_in_Prot[3]<<" "<< " rc\n"; 
      
     
      
 P4_gamma = P4_EL -  P4_E_prime;
//if (isnan(P4_gamma[3])) cout << P4_gamma[3]<<" "<<P4_EL[3]<<" "<<P4_E_prime[3] <<" "<<Ep<<" \n";
//---------------------------

 theta_fermi = acos(pz_fermi/sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
 phi_fermi = acos(fabs(px_fermi)/sqrt(px_fermi*px_fermi+py_fermi*py_fermi));
 
 if ((px_fermi < 0.)&&(py_fermi > 0.)) phi_fermi = M_PI-phi_fermi;
 if ((px_fermi < 0.)&&(py_fermi < 0.)) phi_fermi = phi_fermi + M_PI;
 if ((px_fermi > 0.)&&(py_fermi < 0.)) phi_fermi = 2.*M_PI - phi_fermi;
 


  P4_EL.RotateZ(-phi_fermi);
  P4_EL.RotateY(-theta_fermi);
  
  P4_gamma.RotateZ(-phi_fermi);
  P4_gamma.RotateY(-theta_fermi);
   
  P4_E_prime.RotateZ(-phi_fermi);
  P4_E_prime.RotateY(-theta_fermi);
   
  P4_in_Prot.RotateZ(-phi_fermi);
  P4_in_Prot.RotateY(-theta_fermi);
  


//------------------------------------
 
 Float_t beta;
 
 beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi);
// beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/MP;
 P4_in_Prot.Boost(0,0,-beta);
 P4_EL.Boost(0,0,-beta);
 P4_E_prime.Boost(0,0,-beta);
 P4_gamma.Boost(0,0,-beta);
 

//-----------------------

//ROT2


  theta_rot2 = acos(P4_EL[2]/sqrt(P4_EL[0]*P4_EL[0]+P4_EL[1]*P4_EL[1]+P4_EL[2]*P4_EL[2]));

  phi_rot2 = acos(fabs(P4_EL[0])/sqrt(P4_EL[0]*P4_EL[0]+P4_EL[1]*P4_EL[1]+P4_EL[2]*P4_EL[2]));
 
 if ((P4_EL[0] < 0.)&&(P4_EL[1] > 0.)) phi_rot2 = M_PI-phi_rot2;
 if ((P4_EL[0] < 0.)&&(P4_EL[1] < 0.)) phi_rot2 = phi_rot2 + M_PI;
 if ((P4_EL[0] > 0.)&&(P4_EL[1] < 0.)) phi_rot2 = 2.*M_PI - phi_rot2;


 P4_EL.RotateY(theta_rot2);
 P4_EL.RotateZ(phi_rot2);
 
 P4_gamma.RotateY(theta_rot2);
 P4_gamma.RotateZ(phi_rot2);
  
 P4_E_prime.RotateY(theta_rot2);
 P4_E_prime.RotateZ(phi_rot2);
 
 P4_in_Prot.RotateY(theta_rot2);
 P4_in_Prot.RotateZ(phi_rot2);
 //----------------------------

 ph_e_ferm = P4_E_prime.Phi();
// cout << theta_rot2<<" "<<phi_rot2<<" "<<P4_E_prime.Theta()<<" "<<ph_e_ferm<<" rtyrty0\n";
 
 };
 
 
 
void from_qulab_to_lab(Float_t theta_rot2,Float_t phi_rot2,Float_t theta_e_ferm, Float_t ph_e_ferm,Float_t Es,Float_t Ep,Float_t &Es_lab,Float_t &Ep_lab) {

 Float_t theta_fermi, phi_fermi;
  TLorentzVector P4_Eini_qualab,P4_E_prime_qualab;
 
 theta_fermi = acos(pz_fermi/sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
 
  phi_fermi = acos(fabs(px_fermi)/sqrt(px_fermi*px_fermi+py_fermi*py_fermi));
 
 if ((px_fermi < 0.)&&(py_fermi > 0.)) phi_fermi = M_PI-phi_fermi;
 if ((px_fermi < 0.)&&(py_fermi < 0.)) phi_fermi = phi_fermi + M_PI;
 if ((px_fermi > 0.)&&(py_fermi < 0.)) phi_fermi = 2.*M_PI - phi_fermi;
 
 
 //------------------------------------------
 
  Float_t beta;

 beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi);
 
 //------------------------------------------
 P4_Eini_qualab.SetXYZT(0.,0.,Es,Es); P4_E_prime_qualab.SetXYZT(Ep*cos(ph_e_ferm)*sin(theta_e_ferm),Ep*sin(ph_e_ferm)*sin(theta_e_ferm),Ep*cos(theta_e_ferm),Ep); 
 
// cout<< P4_Eini_qualab[0]<<" "<< P4_Eini_qualab[1]<<" "<< P4_Eini_qualab[2]<<" "<< P4_Eini_qualab[3]<<"   bef\n ";
  P4_Eini_qualab.RotateZ(-phi_rot2);
  P4_Eini_qualab.RotateY(-theta_rot2);
  
  P4_E_prime_qualab.RotateZ(-phi_rot2);
  P4_E_prime_qualab.RotateY(-theta_rot2);
  
  P4_Eini_qualab.Boost(0,0,beta);
  
  P4_E_prime_qualab.Boost(0,0,beta);
  
  P4_Eini_qualab.RotateY(theta_fermi);
  P4_Eini_qualab.RotateZ(phi_fermi);
  
  P4_E_prime_qualab.RotateY(theta_fermi);
  P4_E_prime_qualab.RotateZ(phi_fermi);
  
 // cout << theta_rot2<<" "<<phi_rot2<<" "<<theta_e_ferm<<" "<<ph_e_ferm<<" rtyrty2\n";

   Ep_lab = P4_E_prime_qualab[3];

   Es_lab = P4_Eini_qualab[3];
   
 //  cout<< P4_Eini_qualab[0]<<" "<< P4_Eini_qualab[1]<<" "<< P4_Eini_qualab[2]<<" "<< P4_Eini_qualab[3]<<" aft\n ";
}; 
 
 
 
