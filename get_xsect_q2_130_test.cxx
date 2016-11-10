#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"
#include "interpol_rip2.h"
#include "interpol_rip3.h"
#include "get_xsect_ripani.h"
#include "get_xsect_gol2.h"

#include "get_xsect_14_18_lowq2_fit.h"
#include "get_xsect_golovach.h"
using namespace std;


Short_t getWbin_rip2_test (Float_t W) {
//return int(W*10000. - 1.4125*10000.)/250;
//return int((W-1.4125)/0.025);

if ((W>=1.837499)&&(W<=1.8625)) return 0;
if ((W>=1.8625)&&(W<=1.8875)) return 1;
if ((W>=1.8875)&&(W<=1.9125)) return 2;
if ((W>=1.9125)&&(W<=1.9375)) return 3;

if ((W>=1.9375)&&(W<=1.9625)) return 4;
if ((W>=1.9625)&&(W<=1.9875)) return 5;
if ((W>=1.9875)&&(W<=2.0125)) return 6;
if ((W>=2.0125)&&(W<=2.0375)) return 7;

if ((W>=2.0375)&&(W<=2.0625)) return 8;
if ((W>=2.0625)&&(W<=2.0875)) return 9;
if ((W>=2.0875)&&(W<=2.1125)) return 10;
if ((W>=2.1125)&&(W<=2.1375)) return 11;

if ((W>=2.1375)&&(W<=2.1875)) return 12;
if ((W>=2.1875)&&(W<=2.2375)) return 13;
if ((W>=2.2375)&&(W<=2.2875)) return 14;
if ((W>=2.2875)&&(W<=2.3375)) return 15;
if ((W>=2.3375)&&(W<=2.3875)) return 16;
if ((W>=2.3875)&&(W<=2.4375)) return 17;
if ((W>=2.4375)&&(W<=2.4875)) return 18;
if ((W>=2.4875)&&(W<=2.5375)) return 19;




if ((W < 1.8125)||(W > 2.538)) {
cout << "Error, wrong W range " << W<< " e" << "\n";
return -100;
}

};


Short_t getWbin_rip3_test (Float_t W) {

if ((W>=2.5875)&&(W<=2.6375)) return 0;
if ((W>=2.6375)&&(W<=2.6875)) return 1;
if ((W>=2.6875)&&(W<=2.7375)) return 2;
if ((W>=2.7375)&&(W<=2.7875)) return 3;
if ((W>=2.7875)&&(W<=2.8375)) return 4;
if ((W>=2.8375)&&(W<=2.8875)) return 5;
if ((W>=2.8875)&&(W<=2.9375)) return 6;
if ((W>=2.9375)&&(W<=2.9875)) return 7;
if ((W>=2.9875)&&(W<=3.0375)) return 8;


if ((W <2.58749 )||(W > 3.038)) {
cout << "Error, wrong W range " << W<< " e" << "\n";
return -100;
}

};



Short_t getsbin_rip3 (Short_t Wbin, Float_t sgen, Float_t Smax, Float_t Smin) {
if ((sgen>=Smin)&&(sgen<=Smax)) return int((sgen-Smin)/((Smax - Smin)/15.));
if (sgen<Smin) return 0;
if (sgen>Smax) return 14;



};



void get_xsect_q2_130_test(Float_t Q2gen, Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Float_t phigen, Float_t &sigma_t_final, Float_t &sigma_l_final,Float_t  &sigma_c2f_final,Float_t  &sigma_s2f_final,Float_t &sigma_cf_final,Float_t  &sigma_sf_final){

Short_t Wleft_bin;
Short_t Wright_bin;
Float_t A_tmp[11];
//cout << getWbin_rip2(2.0875) << "\n";

//cout << Wgen <<" "<< Wleft_bin << " " <<Wright_bin <<"\n";
Float_t sigma_t_wright_gol,sigma_t_wleft_gol, sigma_t_rip2,sigma_l_rip2;
Float_t sigma_t_gol;
Short_t w_left_bin_gol;

Short_t s12left_wleft_bin, s12right_wleft_bin, s12left_wright_bin, s12right_wright_bin, s23left_wleft_bin, s23right_wleft_bin, s23left_wright_bin, s23right_wright_bin;
Float_t sigma_final[6];
Float_t sigma_wright[6];
Float_t sigma_wleft[6];

if ((Q2gen>1.27)&&(Q2gen<1.3)){


Short_t thetaleft_bin = getanglebin(thetagen,THETA_ARR[5]);
Short_t thetaright_bin = thetaleft_bin+1;

Short_t alphaleft_bin = getanglebin(alphagen,ALPHA_ARR[5]);
Short_t alpharight_bin = alphaleft_bin+1;
//cout << Wgen << " "<< Wleft_bin << " t \n"; 
if ((Wgen >1.8375)&&(Wgen< 2.5375)){

Wleft_bin = getWbin_rip2_test(Wgen);
 Wright_bin = Wleft_bin+1;
 s12left_wleft_bin = getsbin(Wleft_bin, s12gen, S12_ARR_RIP2[11][Wleft_bin], S12_ARR_RIP2[0][Wleft_bin]);
 s12right_wleft_bin = s12left_wleft_bin +1;

 s12left_wright_bin = getsbin(Wright_bin, s12gen, S12_ARR_RIP2[11][Wright_bin], S12_ARR_RIP2[0][Wright_bin]);
 s12right_wright_bin = s12left_wright_bin +1;

//cout << Wgen << " "<< Wleft_bin << " t \n"; 

 s23left_wleft_bin = getsbin(Wleft_bin, s23gen, S23_ARR_RIP2[11][Wleft_bin], S23_ARR_RIP2[0][Wleft_bin]);
 s23right_wleft_bin = s23left_wleft_bin +1;

 s23left_wright_bin = getsbin(Wright_bin, s23gen, S23_ARR_RIP2[11][Wright_bin], S23_ARR_RIP2[0][Wright_bin]);
 s23right_wright_bin = s23left_wright_bin +1;
 
 for (Short_t i=0;i<6;i++){
interpol_rip2(4,Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright[i],i);

interpol_rip2(4,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft[i],i);
};

for (Short_t i=0;i<6;i++){
sigma_final[i] = 1./fabs(W_ARR_RIP2[Wright_bin]-W_ARR_RIP2[Wleft_bin]);
sigma_final[i] = sigma_final[i]*(sigma_wright[i]*fabs(W_ARR_RIP2[Wleft_bin]-Wgen)+sigma_wleft[i]*fabs(W_ARR_RIP2[Wright_bin]-Wgen));
//cout <<i<< " "<< sigma_final[i] <<"\n";
};

};



if ((Wgen >2.5875)&&(Wgen< 3.0375)){

Wleft_bin = getWbin_rip3_test(Wgen);
 Wright_bin = Wleft_bin+1;
 
 s12left_wleft_bin = getsbin_rip3(Wleft_bin, s12gen, S12_ARR_RIP3[15][Wleft_bin], S12_ARR_RIP3[0][Wleft_bin]);
 s12right_wleft_bin = s12left_wleft_bin +1;

 s12left_wright_bin = getsbin_rip3(Wright_bin, s12gen, S12_ARR_RIP3[15][Wright_bin], S12_ARR_RIP3[0][Wright_bin]);
 s12right_wright_bin = s12left_wright_bin +1;

//cout << Wgen << " "<< Wleft_bin << " tt \n"; 

 s23left_wleft_bin = getsbin_rip3(Wleft_bin, s23gen, S23_ARR_RIP3[15][Wleft_bin], S23_ARR_RIP3[0][Wleft_bin]);
 s23right_wleft_bin = s23left_wleft_bin +1;

 s23left_wright_bin = getsbin_rip3(Wright_bin, s23gen, S23_ARR_RIP3[15][Wright_bin], S23_ARR_RIP3[0][Wright_bin]);
 s23right_wright_bin = s23left_wright_bin +1;
 
 
 for (Short_t i=0;i<6;i++){
interpol_rip3(4,Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright[i],i);

interpol_rip3(4,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft[i],i);
};
for (Short_t i=0;i<6;i++){
sigma_final[i] = 1./fabs(W_ARR_RIP3[Wright_bin]-W_ARR_RIP3[Wleft_bin]);
sigma_final[i] = sigma_final[i]*(sigma_wright[i]*fabs(W_ARR_RIP3[Wleft_bin]-Wgen)+sigma_wleft[i]*fabs(W_ARR_RIP3[Wright_bin]-Wgen));
//cout <<i<< " "<< sigma_final[i] <<"\n";
};

};


if ((Wgen >2.5375)&&(Wgen< 2.5875)){

for (Short_t i=0;i<6;i++){
sigma_final[i] = 0.;
//cout <<i<< " "<< sigma_final[i] <<"\n";
};
};




};

if ((Wgen >2.5875)&&(Wgen< 3.0375)&&(Q2gen>0.0002)&&(Q2gen<0.02)){
get_xsect_gol2(Wgen, s12gen,s23gen, thetagen, alphagen, w_left_bin_gol, sigma_t_wright_gol,sigma_t_wleft_gol );
Wleft_bin = getWbin_rip3_test(Wgen);
 Wright_bin = Wleft_bin+1;
//cout << sigma_t_wright_gol<<" "<< sigma_t_wleft_gol <<" ";

sigma_t_gol = 1./fabs(W_ARR_RIP3[Wright_bin]-W_ARR_RIP3[Wleft_bin]);
sigma_t_gol = sigma_t_gol*(sigma_t_wright_gol*fabs(W_ARR_RIP3[Wleft_bin]-Wgen)+sigma_t_wleft_gol*fabs(W_ARR_RIP3[Wright_bin]-Wgen));
sigma_final[0] = sigma_t_gol;
//cout << sigma_t_gol <<"\n";
};


if ((Wgen >2.5375)&&(Wgen< 2.5875)&&(Q2gen>0.0002)&&(Q2gen<0.02)){


sigma_final[0] = 0.;
//cout <<i<< " "<< sigma_final[i] <<"\n";

};




/*
if ((Wgen >2.5375)&&(Wgen< 2.5875)){

 s12left_wleft_bin = getsbin(Wleft_bin, s12gen, S12_ARR_RIP2[11][20], S12_ARR_RIP2[0][20]);
 s12right_wleft_bin = s12left_wleft_bin +1;
 
 
 s12left_wright_bin = getsbin_rip3(Wright_bin, s12gen, S12_ARR_RIP3[15][0], S12_ARR_RIP3[0][0]);
 s12right_wright_bin = s12left_wright_bin +1;
 
  s23left_wleft_bin = getsbin(Wleft_bin, s23gen, S23_ARR_RIP2[11][20], S23_ARR_RIP2[0][20]);
 s23right_wleft_bin = s23left_wleft_bin +1;
 
 
 s23left_wright_bin = getsbin_rip3(Wright_bin, s23gen, S23_ARR_RIP3[15][0], S23_ARR_RIP3[0][0]);
 s23right_wright_bin = s23left_wright_bin +1;
 
 

};*/


//cout << s23left_wleft_bin <<" "<< Wgen << " "<<Wleft_bin<< " " << s23gen << " "<< S23_ARR[s23left_wleft_bin][Wleft_bin] << " "<< S23_ARR[s23right_wleft_bin][Wleft_bin] << "\n";




//cout << alphaleft_bin << " "<< alphagen<< " "<< ALPHA_ARR[alphaleft_bin] << " "<< ALPHA_ARR[ alpharight_bin] << "\n";





//then we are doing 4d-interpolation for each (Wleft_bin, Q2_left_bin),  (Wright_bin, Q2_left_bin), (Wright_bin, Q2_right_bin) and (Wleft_bin, Q2_right_bin) and obtain cross-secton in that points (4 GLOBAL 6dim arrays)
//0 - sigma_t, 1 - sigma_l, 2 - sigma_c2f, 3 - sigma_s2f, 4 - sigma_cf, 5 - sigma_sf





//cout << Wgen <<"  "<< Wleft_bin<< " "<<sigma_wright[0]<< " "<< sigma_wleft[0]<< "\n";





//cout << W_ARR_RIP2[Wright_bin] << "\n";
//1dim W-interpolation

//We get explicitly different sigmas from the array
 sigma_t_final = sigma_final[0];
  sigma_l_final = 0.;
// sigma_l_final = sigma_final[1];
 sigma_c2f_final = 0.;
 sigma_s2f_final = 0.;
 sigma_cf_final = 0.;
 sigma_sf_final = 0.;

};
