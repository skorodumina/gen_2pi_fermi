#include "TROOT.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TParticlePDG.h"

#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TMacro.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TText.h"
#include "TStyle.h"
#include "TObject.h"
#include <TRint.h>
#include <stdio.h>
#include <dlfcn.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include <iostream>

#include "inp_file_read.h"
#include "read_xsect_files.h"
#include "read_fit_param_files.h"
#include "out_file_write.h"
#include "out_file_open.h"
#include "out_file_close.h"
#include "anti_rot.h"
#include "get_xsect_ripani.h"
#include "get_xsect_near_threshold.h"

#include "get_xsect_golovach.h"
#include "get_xsect_gol2.h"
#include "get_xsect_fedotov.h"
#include "get_xsect_rip_fed_join.h"
#include "get_xsect_14_18_lowq2_fit.h"
#include "get_xsect_scale_gol_18_25.h"
#include "get_xsect_q2_130_w_gt_18_lt_21.h"
#include "get_xsect_25_30.h"
#include "get_xsect_q2_13_wgt_3.h"
#include "get_xsect_q2_130_test.h"

#include "interpol_int.h"


#include "radcorr.h"
#include "fermi_bonn.h"
#include "fermi_rot.h"
#include "fermi_anti_rot.h"

#include "rot.h"
#include "global.h"
#include <stdlib.h>
#include <time.h>
#include <TLorentzVector.h>
#include <sstream>
#include <TRandom3.h>
#include <fstream>
//for BOS-files creation
/*
extern "C" 
{

#include <signal.h>
#include <errno.h>
#include <ntypes.h>
#include <bostypes.h>
#include <ec.h>
#include <clas_cern.h>
#include <ctype.h>
#include <kinematics.h>
#include <map_manager.h>
#include <trk.h>
#include <clasmdl.h>
#include <utility.h>
#include <pid.h>
#include <makebanks.h>
#include <call.h>
#include <bosddl.h>
#include <tagtnorm.h>
#include <vertex.h>

  void initbos();
  int  getBOS(BOSbank *bcs, int lun, char *list);
  void cleanBanks(BOSbank *bcs);
  void dropAllBanks(BOSbank *bcs, char *list);
  void *getBank(BOSbank *,const char *);
  void open_fpack_unit(char *filename,char *dataname,int unitnum);
  void close_fpack_unit(char *dataname);
  BOSbank bcs_ ;
  BOSbank wcs_ ;

  float ranf_();
  int ranset_( float* );
  void ranlux_( float* , const int* ) ;
 


//  time_t   time( time_t );
//  int   clock();
}

*/
 using namespace std;      
     
     
     


//Byckling function declaration
     Float_t G_BYCKLING(Float_t x, Float_t y, Float_t z, Float_t u, Float_t v, Float_t w) {
     return x*x*y+x*y*y+z*z*u+z*u*u+v*v*w+v*w*w+x*z*w+x*u*v+y*z*v+y*u*w-x*y*(z+u+v+w)-z*u*(x+y+v+w)-v*w*(x+y+z+u);
     };



int main(int argc, char** argv) {




    vector<string> args;
   
    
    Float_t sigma_t_final = 0.;
    Float_t sigma_l_final = 0.; 
    Float_t sigma_c2f_final = 0.;
    Float_t sigma_s2f_final = 0.;
    Float_t sigma_cf_final = 0.;
   Float_t  sigma_sf_final = 0.;
   
    Float_t tmp = 0.;
   
   Float_t sigma_t_final_1, sigma_l_final_1 , sigma_c2f_final_1, sigma_s2f_final_1,sigma_cf_final_1 ,sigma_sf_final_1 ;
   
   Float_t sigma_t_final_2, sigma_l_final_2 , sigma_c2f_final_2, sigma_s2f_final_2,sigma_cf_final_2 ,sigma_sf_final_2; 
   
   Float_t  E_beam_fermi, theta_rot2,phi_rot2,eps_l,eps_t;
    Float_t  eps_l_for_flux,eps_t_for_flux,Ebeam_ferm_for_flux,Ep_ferm_for_flux;
    Float_t  E_beam,sigma_total, sigma_total_1, sigma_total_2,xsect_int_test_t,xsect_int_test_l,p_el_test;
    Float_t Wnew, Q2new, E_beam_new,Ep_new,ph_e_qualab;
    
    Float_t V_flux;
    Float_t alph_const = 1./137.035;
  //  sigma_total = 0.;
    
     Float_t W,Q2,phi_e,W_old,Q2nodata,Q2_old;
     const Float_t phi_e_min = 0;
     const Float_t phi_e_max = 2*M_PI;
     string* file = NULL;
     
    px_fermi =0.;
    py_fermi =0.;
    pz_fermi =0.;
//This needed for taking masses of the particles from pdg_table located in ROOT_DIR


    const char *HOME_ROOT;
    const char *HOME_ROOT1;
//HOME_ROOT = getenv("ROOT");
system("root_home=`root-config --etcdir`");
HOME_ROOT = getenv("root_home");
ostringstream ROOT_DIR;
ROOT_DIR << HOME_ROOT << "/pdg_table.txt";



TDatabasePDG *pdg = new TDatabasePDG();
pdg->ReadPDGTable(ROOT_DIR.str().c_str());
TParticlePDG *part1 = new TParticlePDG();
part1 = pdg->GetParticle("proton");
MP= part1->Mass();
part1 = pdg->GetParticle("pi+");
MPIP= part1->Mass();  
part1 = pdg->GetParticle("pi-");
MPIM= part1->Mass();  
part1 = pdg->GetParticle("e-");
Me= part1->Mass();  


   //Reading input parameters
   inp_file_read(E_beam);

   
   //Reading diff cross section from the tables in .dat files (filling out GLOBAL arrays)
   read_xsect_files();
   
 
   read_fit_param_files();
   

     
   //Reasonably changing kinematical variables if needed
    if (Q2_max > 4.*E_beam*E_beam*sin(Theta_max*M_PI/180./2.)*sin(Theta_max*M_PI/180./2.)) {
    Q2_max = 4.*E_beam*E_beam*sin(Theta_max*M_PI/180./2.)*sin(Theta_max*M_PI/180./2.);
    cout << "maximum Q2 has been changed to " << Q2_max << "\n";
    };
    
     if (Q2_min < 4.*E_beam*E_eprime_min*sin(Theta_min*M_PI/180./2.)*sin(Theta_min*M_PI/180./2.)) {
    Q2_min = 4.*E_beam*E_eprime_min*sin(Theta_min*M_PI/180./2.)*sin(Theta_min*M_PI/180./2.);
    cout << "minimum Q2 has been changed to " << Q2_min << "\n";
    };   
    
    if (W_max*W_max > (MP*MP+2.*MP*(E_beam-E_eprime_min) - Q2_min)) {
    W_max = sqrt(MP*MP+2.*MP*(E_beam-E_eprime_min)- Q2_min);
    cout << "maximum W has been changed to " << W_max << "\n";
    };
    
 if (W_min < (MP + MPIP + MPIM + 0.01)) {
    W_min = MP + MPIP + MPIM + 0.01;
   
    cout << "minimum W has been changed to " << W_min << "\n";
   
    };
    
  //  cout << "ttt1 \n";
    TH1F *h_0_miss = new TH1F("h_0_miss","h_0_miss",1000,-0.05,0.05);
    TH1F *h_pim_miss = new TH1F("h_pim_miss","h_pim_miss",200,MPIM*MPIM-0.05,MPIM*MPIM+0.05);
   
    TH1F *h_0_miss_2 = new TH1F("h_0_miss_2","h_0_miss_2",1000,-0.05,0.05);
    TH1F *h_pim_miss_2 = new TH1F("h_pim_miss_2","h_pim_miss_2",200,MPIM*MPIM-0.05,MPIM*MPIM+0.05);
    
    
    
    TH1F *h_0_miss_fermi = new TH1F("h_0_miss_fermi","h_0_miss_fermi",1000,-0.05,0.05);
    TH1F *h_pim_miss_fermi = new TH1F("h_pim_miss_fermi","h_pim_miss_fermi",200,MPIM*MPIM-0.1,MPIM*MPIM+0.1);
   
    TH1F *h_0_miss_fermi_2 = new TH1F("h_0_miss_fermi_2","h_0_miss_fermi_2",1000,-0.05,0.05);
    TH1F *h_pim_miss_fermi_2 = new TH1F("h_pim_miss_fermi_2","h_pim_miss_fermi_2",200,MPIM*MPIM-0.1,MPIM*MPIM+0.1);
    
     TH1F *h_0_miss_en_fermi = new TH1F("h_0_miss_en_fermi","h_0_miss_en_fermi",200,-0.05,0.05);
    TH1F *h_0_miss_en_fermi_2 = new TH1F("h_0_miss_en_fermi_2","h_0_miss_en_fermi_2",200,-0.05,0.05);
    
     TH1F *h_fermi_bonn = new TH1F("h_fermi_bonn","h_fermi_bonn",100,-0.1,0.9);
    
     TH1F *h_0_miss_en = new TH1F("h_0_miss_en","h_0_miss_en",200,-0.05,0.05);
    TH1F *h_0_miss_en_2 = new TH1F("h_0_miss_en_2","h_0_miss_en_2",200,-0.05,0.05);
    
    TH1F *h_eradgam = new TH1F("h_eradgam","h_eradgam",500,0.,0.6);
   
    
    TH1F *h_W = new TH1F("W","W",100,W_min,W_max);
    TH1F *h_W_2 = new TH1F("W_2","W_2",20,1.3,1.8);
    
    TH1F *h_Q2 = new TH1F("Q2","Q2",250,Q2_min,Q2_max);
    TH1F *h_phi_e = new TH1F("phi_e","phi_e",100,phi_e_min,phi_e_max);
    TH2F *h_Q2vsW = new TH2F("Q2vsW","Q2vsWW",69,W_min,W_max,4,Q2_min,Q2_max);
    TH2F *h_Q2vsW2 = new TH2F("Q2vsW2","Q2vsWW2",69,W_min,W_max,4,Q2_min,Q2_max);
    
     TH2F *h_Q2vsW_t = new TH2F("Q2vsW_t","Q2vsWW_t",71,W_min,W_max,25,Q2_min,Q2_max);
     TH2F *h_Q2vsW_l = new TH2F("Q2vsW_l","Q2vsWW_l",71,W_min,W_max,25,Q2_min,Q2_max);
    
    
     TH2F *h_eps_l = new TH2F("eps_l","eps_l",100,W_min,W_max,100,Q2_min,Q2_max);
    
    TH1F *h_nu = new TH1F("nu","nu",100,-1.*E_beam,E_beam);
    TH1F *h_zel = new TH1F("Z_EL","Z_EL",100,Targ_off-Targ_len/2.-1.,Targ_off+Targ_len/2.+1.);
    
    TH1F *h_inv_m12 = new TH1F("h_inv_m12","h_inv_m12",100,(MPIP+MPIM)*(MPIP+MPIM)-0.02,(1.625-MP)*(1.625-MP)+0.02);
    TH1F *h_inv_m23 = new TH1F("h_inv_m23","h_inv_m23",100,(MPIP+MP)*(MPIP+MP)-0.02,(1.625-MPIM)*(1.625-MPIM)+0.02);
     TH1F *h_th_hadr = new TH1F("h_th_hadr","h_th_hadr",100,0,M_PI);
      TH1F *h_th_hadr_2 = new TH1F("h_th_hadr_2","h_th_hadr_2",100,0,M_PI);
     TH1F *h_ph_hadr = new TH1F("h_ph_hadr","h_ph_hadr",100,0,2.*M_PI);
     TH1F *h_alph_hadr = new TH1F("h_alph_hadr","h_alph_hadr",100,0,2.*M_PI);
    
    TH2F *h_dalitz = new TH2F("dalitz","dalitz",100,MPIM+MPIP,W_max-MP,100,MPIP+MP,W_max-MPIM);
    
    TH1F *h_odn_inv_m12[37];
     TH1F *h_odn_inv_m23[37];
      TH1F *h_odn_alpha[37];
      TH1F *h_odn_theta[37]; 
     TH1F *h_odn_theta_2[37];
     
     
     
      TH1F *h_odn_wwide_inv_m12[18];
     TH1F *h_odn_wwide_inv_m23[18];
      TH1F *h_odn_wwide_alpha[18];
      TH1F *h_odn_wwide_theta[18]; 
     TH1F *h_odn_wwide_theta_2[18];
     
     TH1F *h_odn_wgt3_inv_m12[15];
     TH1F *h_odn_wgt3_inv_m23[15];
     TH1F *h_odn_wgt3_alpha[15];
     TH1F *h_odn_wgt3_theta[15]; 
     TH1F *h_odn_wgt3_theta_2[15];
     
     TH1F *h_odn_q2_dep_t[33];
     TH1F *h_odn_q2_dep_l[33];
     TH1F *h_odn_q2_dep_l2[33];
     TH1F *h_odn_q2_dep_tot[33];
     
     TH1F *h_odn_w_dep_t[12];
      TH1F *h_odn_w_dep_l[12];
      TH1F *h_odn_w_dep_l2[12];
       TH1F *h_odn_w_dep_tot[12];
       
         TH1F *h_int_crsect_t[27];
         TH1F *h_int_crsect_l[27];
	 
	 
	 TTree*t21 = new TTree("h10","Tree h10"); 
  t21->SetDirectory(0);
  t21->Branch("sigma",&sigma_total); 
  t21->Branch("p_el_test",&p_el_test); 
  t21->Branch("px_fermi",&px_fermi); 
  t21->Branch("py_fermi",&py_fermi);
  t21->Branch("pz_fermi",&pz_fermi);
	 
       
       ostringstream qqq;
       
       
    h_int_crsect_t[0] = new TH1F("h_int_crsect_t_0","h_int_crsect_t_0",71, 1.25, 3.025);   
      h_int_crsect_l[0] = new TH1F("h_int_crsect_l_0","h_int_crsect_l_0",71, 1.25, 3.025);  
      
      h_int_crsect_t[26] = new TH1F("h_int_crsect_t_26","h_int_crsect_t_26",71, 1.25, 3.025);   
      h_int_crsect_l[26] = new TH1F("h_int_crsect_l_26","h_int_crsect_l_26",71, 1.25, 3.025);   
      
        
       
       for (Short_t ii=1; ii<=25; ii++) {
       
        qqq << "h_int_crsect_t_" <<ii;
     h_int_crsect_t[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),71, 1.25, 3.025);
       qqq.str("");
       
       qqq << "h_int_crsect_l_" <<ii;
     h_int_crsect_l[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),71, 1.25, 3.025);
       qqq.str("");
       
       
       };
       
       
      for (Short_t ii=0; ii<=11; ii++) {
      
      qqq << "h_odn_w_dep_t_" << 100*(0.15+0.1*ii);
h_odn_w_dep_t[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100, W_min, W_max);
qqq.str("");

  qqq << "h_odn_w_dep_l_" << 100*(0.15+0.1*ii);
h_odn_w_dep_l[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100, W_min, W_max);
qqq.str("");


  qqq << "h_odn_w_dep_l2_" << 100*(0.15+0.1*ii);
h_odn_w_dep_l2[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100, W_min, W_max);
qqq.str("");



 qqq << "h_odn_w_dep_tot_" << 100*(0.15+0.1*ii);
h_odn_w_dep_tot[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100, W_min, W_max);
qqq.str("");


      
      };
      
     for (Short_t ii=0; ii<=32; ii++) { 
     
     qqq << "h_odn_q2_dep_t_" << 10000*(1.2625+0.025*ii);
h_odn_q2_dep_t[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
qqq.str("");

qqq << "h_odn_q2_dep_l_" << 10000*(1.2625+0.025*ii);
h_odn_q2_dep_l[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
qqq.str("");

qqq << "h_odn_q2_dep_l2_" << 10000*(1.2625+0.025*ii);
h_odn_q2_dep_l2[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
qqq.str("");


qqq << "h_odn_q2_dep_tot_" << 10000*(1.2625+0.025*ii);
h_odn_q2_dep_tot[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
qqq.str("");

     
     };
      
      
      
       
    
    for (Short_t ii=0; ii<=36; ii++) {
qqq << "h_odn_inv_m12_" << 10000*(1.2375+0.025*ii);
h_odn_inv_m12[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MPIM-0.02 ,1.2375+0.025*ii-MP+0.02 );
qqq.str("");

qqq << "h_odn_inv_m23_" << 10000*(1.2375+0.025*ii);
h_odn_inv_m23[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MP-0.02 ,1.2375+0.025*ii-MPIM+0.02 );
qqq.str("");

qqq << "h_odn_alpha_" << 10000*(1.2375+0.025*ii);
h_odn_alpha[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,2*M_PI );
qqq.str("");

qqq << "h_odn_theta_" << 10000*(1.2375+0.025*ii);
h_odn_theta[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");

qqq << "h_odn_theta_2_" << 10000*(1.2375+0.025*ii);
h_odn_theta_2[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");

//qqq << "h_odn_q2_dep_" << 10000*(1.4375+0.025*ii);
///h_odn_q2_dep[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
//qqq.str("");


};


 for (Short_t ii=0; ii<18; ii++) {
qqq << "h_odn_wwide_inv_m12_" << 10000*(2.1875+0.05*ii);
h_odn_wwide_inv_m12[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MPIM-0.02 ,2.1875+0.05*ii-MP+0.02 );
qqq.str("");

qqq << "h_odn_wwide_inv_m23_" << 10000*(2.1875+0.05*ii);
h_odn_wwide_inv_m23[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MP-0.02 ,2.1875+0.05*ii-MPIM+0.02 );
qqq.str("");

qqq << "h_odn_wwide_alpha_" << 10000*(2.1875+0.05*ii);
h_odn_wwide_alpha[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,2*M_PI );
qqq.str("");

qqq << "h_odn_wwide_theta_" << 10000*(2.1875+0.05*ii);
h_odn_wwide_theta[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");

qqq << "h_odn_wwide_theta_2_" << 10000*(2.1875+0.05*ii);
h_odn_wwide_theta_2[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");

//qqq << "h_odn_q2_dep_" << 10000*(1.4375+0.025*ii);
///h_odn_q2_dep[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
//qqq.str("");


};


for (Short_t ii=0; ii<15; ii++) {
qqq << "h_odn_wgt3_inv_m12_" << 10000*(3.1375+0.1*ii);
h_odn_wgt3_inv_m12[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MPIM-0.02 ,3.1375+0.1*ii-MP+0.02 );
qqq.str("");

qqq << "h_odn_wgt3_inv_m23_" << 10000*(3.1375+0.1*ii);
h_odn_wgt3_inv_m23[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MP-0.02 ,3.1375+0.1*ii-MPIM+0.02 );
qqq.str("");

qqq << "h_odn_wgt3_alpha_" << 10000*(3.1375+0.1*ii);
h_odn_wgt3_alpha[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,2*M_PI );
qqq.str("");

qqq << "h_odn_wgt3_theta_" << 10000*(3.1375+0.1*ii);
h_odn_wgt3_theta[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");

qqq << "h_odn_wgt3_theta_2_" << 10000*(3.1375+0.1*ii);
h_odn_wgt3_theta_2[ii] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");




};



    
    
    Float_t nu,E_E_prime,Theta_e_prime,E_E_prime_new,Theta_e_prime2,E_E_prime_new2;
    Float_t M1,M2,M3;
    TLorentzVector P4_E_prime, P4_PIP,P4_Pfin,P4_PIM;
    
    TLorentzVector P4_Eini, P4_Eini_new, P4_E_prime_new, P4_Pini, P4_E_prime_qulab_rad;
    TLorentzVector P4_Eini_lab_rad,P4_Eini_qulab_rad;
    TLorentzVector P4_0_miss, P4_pim_miss; 
    TLorentzVector P4_0_miss_2, P4_pim_miss_2; 
    
     TLorentzVector P4_0_miss_fermi, P4_pim_miss_fermi; 
    TLorentzVector P4_0_miss_fermi_2, P4_pim_miss_fermi_2; 
    
   TVector3 P3_P_in_D;
     TLorentzVector P4_E_prime_boosted,P4_gamma_test,P4_gamma_test2, P4_Eini_fermi;
     
     TLorentzVector P4_Pini_fermi;
    
    P4_Pini.SetXYZT(0.,0.,0.,MP); 
    P4_Eini.SetXYZT(0.,0.,E_beam,E_beam);
   
    
    srand (time(NULL));
    Int_t k=0;
    Float_t z_EL, x_EL, y_EL;;
    Float_t r_vert, phi_vert;
     Float_t e_rad_phot, cr_rad_fact;
     
   Float_t E_E_prime_ferm,Theta_e_prime_ferm,W_ferm;
   
   Float_t inv_m12,inv_m23,inv_m13,th_hadr,alph_hadr,ph_hadr,s12,s23,s13,en1,en2,en3,mag1,mag2,mag3; 
   Float_t W_tmp,dummy;
   Int_t dummy_int,ii;
   
     //FOR 1-PIM, 2-PIP, 3-P
     M1 = MPIM;
     M2 = MPIP;
     M3 = MP; 
     //FOR 1-P, 2-PIP, 3-PIM
//     M1 = MP;
//    M2 = MPIP;
//     M3 = MPIM;
     //FOR 1-PIP, 2-PIM, 3-P
 //    M1 = MPIP;
 //    M2 = MPIM;
 //    M3 = MP;
    
  //open input file
     out_file_open();
     
     TRandom3 ph_e_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
     TRandom3 th_hadr_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
     TRandom3 W_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
     TRandom3 Q2_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
     TRandom3 z_EL_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
     TRandom3 alph_hadr_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.)); 
     TRandom3 s12_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
     TRandom3 s23_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));     
     TRandom3 ph_hadr_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));          

     TRandom3 r_vert_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
     TRandom3 phi_vert_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
  // Start to generate electrons    
    for (Int_t i=1; i<=Nevents; i++) {
    
phi_e =ph_e_rndm.Uniform(phi_e_min,phi_e_max);
z_EL = z_EL_rndm.Uniform(Targ_off  - Targ_len/2.,Targ_off  + Targ_len/2.);

r_vert = r_vert_rndm.Uniform(0.,Targ_rad);
phi_vert = phi_vert_rndm.Uniform(0.,6.28318);
x_EL = r_vert*cos(phi_vert);
y_EL = r_vert*sin(phi_vert);

alph_hadr = alph_hadr_rndm.Uniform(0.,6.28318);
ph_hadr = ph_hadr_rndm.Uniform(0.,6.28318);    
    do {
    
    k++;
    if ((k % 1000) == 0) cout  << "N generated = " << k << "; N accepted = " << i << "\n";
    W = W_rndm.Uniform(W_min,W_max);
    Q2 = Q2_rndm.Uniform(Q2_min,Q2_max);
    
    
   
    nu=(W*W+Q2-MP*MP)/2./MP;
     E_E_prime=E_beam-nu;

  
    } while ((nu > E_beam-E_eprime_min)||(isnan(acos(1.-Q2/E_beam/E_E_prime/2.)))||(isnan(acos((Q2+2.*E_beam*nu)/2./E_beam/(sqrt(Q2+nu*nu))))));

    E_E_prime=E_beam-nu;
    Theta_e_prime = acos(1.-Q2/E_beam/E_E_prime/2.);
  // cout<< phi_e<<" "<<Theta_e_prime<<" "<<E_E_prime<<" mp\n";
  // cout << phi_e<<" "<<E_E_prime<<" "<< Theta_e_prime<<" "<<nu<<" mp\n";
  P4_E_prime.SetXYZT(E_E_prime*cos(phi_e)*sin(Theta_e_prime),E_E_prime*sin(phi_e)*sin(Theta_e_prime),E_E_prime*cos(Theta_e_prime),E_E_prime);
   
  // cout << P4_E_prime[0]<<" "<<P4_E_prime[1]<<" "<<P4_E_prime[2]<<" "<<P4_E_prime[3]<<" 1\n";
    
 // cout << E_beam<<" "<<E_E_prime<<" "<< Theta_e_prime<<" "<<phi_e<<" mp\n";
    E_beam_new = E_beam;

    W_old = W;
    Q2_old = Q2;
    
    Wnew = W;
    Q2new = Q2;
 
      

    E_beam_fermi = E_beam_new;
    
// if (flag_fermi == 1) {
  do {
 
    fermi_bonn();
 

                          P4_Pini_fermi.SetXYZT(px_fermi,py_fermi,pz_fermi,sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
	
W_ferm = (P4_Pini_fermi + P4_Eini - P4_E_prime).Mag();
 		
   

// } while (W_ferm < 1.22741 );
 } while (W_ferm < 1.2375 );
 // cout << P4_E_prime[0]<<" "<< P4_E_prime[1]<<" " << P4_E_prime[2]<<" "<< P4_E_prime[3]<<" "<< " mp\n";
 // cout << P4_Pini_fermi[0]<<" "<< P4_Pini_fermi[1]<<" " << P4_Pini_fermi[2]<<" "<< P4_Pini_fermi[3]<<" "<< " mp\n";
    radcorr(E_beam,Q2,W_old,W_ferm,Wnew,Q2new,E_beam_new,Ep_new,E_beam_fermi,eps_l,eps_t,e_rad_phot,cr_rad_fact,phi_e,Theta_e_prime);
    h_eradgam->Fill(e_rad_phot,1.);
    
get_rot2(E_beam_new, Ep_new, phi_e,Theta_e_prime, theta_rot2,  phi_rot2,  ph_e_qualab);
    
    W = Wnew;
    Q2 = Q2new;
   
    P4_Eini_lab_rad.SetXYZT(0.,0.,E_beam_new, E_beam_new);
    
    
    P4_Eini_qulab_rad.SetXYZT(0.,0.,E_beam_fermi, E_beam_fermi);  
    
    E_E_prime_new2 = E_beam_fermi - (W*W+Q2-MP*MP)/2./MP; 
    Theta_e_prime2 = acos(1.-Q2/E_beam_fermi/E_E_prime_new2/2.); 
     
P4_E_prime_qulab_rad.SetXYZT(E_E_prime_new2*cos(ph_e_qualab)*sin(Theta_e_prime2),E_E_prime_new2*sin(ph_e_qualab)*sin(Theta_e_prime2),E_E_prime_new2*cos(Theta_e_prime2),E_E_prime_new2);


 
  // };
   
   
     
  //  W_tmp = W;
   // W = 1.6125;
    dummy_int = 0;
 
    do {

   s12 = s12_rndm.Uniform((M1+M2)*(M1+M2),(W-M3)*(W-M3));
   s23 = s23_rndm.Uniform((M2+M3)*(M2+M3),(W-M1)*(W-M1));
   inv_m12 = sqrt(s12);
   inv_m23 = sqrt(s23); 
   
   
   //this variables are calculated here for the check of correct near-boundaries generation - conditions in while (...)
   s13 = W*W+M1*M1+M2*M2+M3*M3-s12-s23;
   en1 = (W*W+M1*M1-s23)/2./W;
   en2 = (W*W+M2*M2-s13)/2./W;
   en3 = (W*W+M3*M3-s12)/2./W;
   mag1 = sqrt(en1*en1 - M1*M1);
   mag2 = sqrt(en2*en2 - M2*M2);
   mag3 = sqrt(en3*en3 - M3*M3);
   
 
   
   dummy_int++;



    } while ((G_BYCKLING(inv_m12*inv_m12,inv_m23*inv_m23,W*W,M2*M2,M1*M1,M3*M3) > 0.)||(isnan(acos((M1*M1+M2*M2+2*en1*en2-s12)/2./mag1/mag2)))||(isnan(acos((M1*M1+M3*M3+2*en1*en3-s13)/2./mag1/mag3)))||(en1 < M1)||(en2 < M2)||(en3 < M3)||(sqrt(s13) < M1+M3)||(sqrt(s13)>W-M2));
   
  
   h_inv_m12 ->Fill(s12,1.);
   

th_hadr = acos(th_hadr_rndm.Uniform(-1.,1.));
	
Q2nodata = Q2;


if (Q2 > 1.299)Q2 = 1.299;
if (Q2 < 0.0005)Q2 = 0.0005;

if (W<1.2375){

sigma_t_final = 0.;
sigma_l_final = 0.;
sigma_c2f_final = 0.;
sigma_s2f_final = 0.;
sigma_cf_final = 0.;
sigma_sf_final = 0.;
};

//Getting cross section in given generated (W, Q2, s12, s23, theta, alpha)-point
sigma_total = 0.;

 if ((W>=1.4125)&&(W<=1.8125)&&(Q2>=0.65)&&(Q2<=1.3)) {
 get_xsect_ripani(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);

};



 if ((W>=1.6125)&&(W<=1.8125)&&(Q2>0.000001)&&(Q2<0.65)) {

 get_xsect_14_18_lowq2_fit(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);
};



if (((W>=1.3125)&&(W<=1.4375)&&(Q2>0.275)&&(Q2<0.575))||((W>=1.4125)&&(W<=1.4875)&&(Q2>0.275)&&(Q2<0.525))||((W>=1.4125)&&(W<=1.5125)&&(Q2>0.275)&&(Q2<0.425))||((W>=1.4125)&&(W<=1.5375)&&(Q2>0.225)&&(Q2<0.375))||((W>=1.4125)&&(W<=1.5625)&&(Q2>0.225)&&(Q2<0.275))){

 get_xsect_fedotov(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);


};

 if (((W>=1.4125)&&(W<=1.4375)&&(Q2>=0.575)&&(Q2<=0.65))||((W>=1.4375)&&(W<=1.4875)&&(Q2>=0.525)&&(Q2<=0.65))||((W>=1.4875)&&(W<=1.5125)&&(Q2>=0.425)&&(Q2<=0.65))||((W>=1.5125)&&(W<=1.5375)&&(Q2>=0.325)&&(Q2<=0.65))||((W>=1.5375)&&(W<=1.5625)&&(Q2>=0.275)&&(Q2<=0.65))||((W>=1.5625)&&(W<=1.5875)&&(Q2>=0.225)&&(Q2<=0.65))||((W>=1.3125)&&(W<1.4125)&&(Q2>=0.575)&&(Q2<=1.3))||((W>=1.3125)&&(W<1.5125)&&(Q2>=0.0002)&&(Q2<=0.275))||((W>=1.5125)&&(W<=1.5875)&&(Q2>=0.0002)&&(Q2<=0.225))){
 
 get_xsect_rip_fed_join(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final); 

 };
 
  

 if ((W>=1.2375)&&(W<1.3125)&&(Q2>0.00002)&&(Q2<1.3)) {


get_xsect_near_threshold(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final );


};


 if ((W>=1.8375)&&(W<=2.5375)&&(Q2>0.00002)&&(Q2<1.3)) {
get_xsect_q2_130_w_gt_18_lt_21(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final );


};

 if ((W>=2.5875)&&(W<=3.0375)&&(Q2>0.00002)&&(Q2<1.3)) {
get_xsect_25_30(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final );

};






if ((W>=1.5875)&&(W<=1.6125)&&(Q2>0.00001)&&(Q2<0.65)) {
get_xsect_14_18_lowq2_fit(Q2, 1.6125, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_2, sigma_l_final_2,sigma_c2f_final_2,sigma_s2f_final_2,sigma_cf_final_2,sigma_sf_final_2);

 get_xsect_rip_fed_join(Q2,1.5875, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_1, sigma_l_final_1, sigma_c2f_final_1,sigma_s2f_final_1,sigma_cf_final_1,sigma_sf_final_1); 
 
sigma_t_final = 1./0.025;
sigma_t_final = sigma_t_final*(sigma_t_final_2*fabs(1.5875-W)+sigma_t_final_1*fabs(1.6125-W));

sigma_l_final = 1./0.025;
sigma_l_final = sigma_l_final*(sigma_l_final_2*fabs(1.5875-W)+sigma_l_final_1*fabs(1.6125-W));

sigma_c2f_final = 1./0.025;
sigma_c2f_final = sigma_c2f_final*(sigma_c2f_final_2*fabs(1.5875-W)+sigma_c2f_final_1*fabs(1.6125-W));

sigma_s2f_final = 1./0.025;
sigma_s2f_final = sigma_s2f_final*(sigma_s2f_final_2*fabs(1.5875-W)+sigma_s2f_final_1*fabs(1.6125-W));

sigma_cf_final = 1./0.025;
sigma_cf_final = sigma_cf_final*(sigma_cf_final_2*fabs(1.5875-W)+sigma_cf_final_1*fabs(1.6125-W));

sigma_sf_final = 1./0.025;
sigma_sf_final = sigma_sf_final*(sigma_sf_final_2*fabs(1.5875-W)+sigma_sf_final_1*fabs(1.6125-W));


};

if ((W>=1.8125)&&(W<=1.8375)&&(Q2>0.00002)&&(Q2<=1.3)) {


if ((Q2>0.00002)&&(Q2<=0.65)) get_xsect_14_18_lowq2_fit(Q2, 1.8125, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_2, sigma_l_final_2,sigma_c2f_final_2,sigma_s2f_final_2,sigma_cf_final_2,sigma_sf_final_2);


if ((Q2>0.65)&&(Q2<=1.3)) get_xsect_ripani(Q2, 1.8125, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_2, sigma_l_final_2,sigma_c2f_final_2,sigma_s2f_final_2,sigma_cf_final_2,sigma_sf_final_2);



get_xsect_q2_130_w_gt_18_lt_21(Q2, 1.8375, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_1, sigma_l_final_1, sigma_c2f_final_1,sigma_s2f_final_1,sigma_cf_final_1,sigma_sf_final_1); 

sigma_t_final = 1./0.025;
sigma_t_final = sigma_t_final*(sigma_t_final_2*fabs(1.8375-W)+sigma_t_final_1*fabs(1.8125-W));

sigma_l_final = 1./0.025;
sigma_l_final = sigma_l_final*(sigma_l_final_2*fabs(1.8375-W)+sigma_l_final_1*fabs(1.8125-W));

sigma_c2f_final = 1./0.025;
sigma_c2f_final = sigma_c2f_final*(sigma_c2f_final_2*fabs(1.8375-W)+sigma_c2f_final_1*fabs(1.8125-W));

sigma_s2f_final = 1./0.025;
sigma_s2f_final = sigma_s2f_final*(sigma_s2f_final_2*fabs(1.8375-W)+sigma_s2f_final_1*fabs(1.8125-W));

sigma_cf_final = 1./0.025;
sigma_cf_final = sigma_cf_final*(sigma_cf_final_2*fabs(1.8375-W)+sigma_cf_final_1*fabs(1.8125-W));

sigma_sf_final = 1./0.025;
sigma_sf_final = sigma_sf_final*(sigma_sf_final_2*fabs(1.8375-W)+sigma_sf_final_1*fabs(1.8125-W));
};






if ((W>=2.5375)&&(W<=2.5875)&&(Q2>0.00001)&&(Q2<1.3)) {
get_xsect_q2_130_w_gt_18_lt_21(Q2, 2.5375, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_2, sigma_l_final_2,sigma_c2f_final_2,sigma_s2f_final_2,sigma_cf_final_2,sigma_sf_final_2);

 get_xsect_25_30(Q2,2.5875, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_1, sigma_l_final_1, sigma_c2f_final_1,sigma_s2f_final_1,sigma_cf_final_1,sigma_sf_final_1); 
 
sigma_t_final = 1./0.05;
sigma_t_final = sigma_t_final*(sigma_t_final_2*fabs(2.5875-W)+sigma_t_final_1*fabs(2.5375-W));

sigma_l_final = 1./0.05;
sigma_l_final = sigma_l_final*(sigma_l_final_2*fabs(2.5875-W)+sigma_l_final_1*fabs(2.5375-W));

sigma_s2f_final = 1./0.05;
sigma_s2f_final = sigma_s2f_final*(sigma_s2f_final_2*fabs(2.5875-W)+sigma_s2f_final_1*fabs(2.5375-W));

sigma_c2f_final = 1./0.05;
sigma_c2f_final = sigma_c2f_final*(sigma_c2f_final_2*fabs(2.5875-W)+sigma_c2f_final_1*fabs(2.5375-W));

sigma_cf_final = 1./0.05;
sigma_cf_final = sigma_cf_final*(sigma_cf_final_2*fabs(2.5875-W)+sigma_cf_final_1*fabs(2.5375-W));

sigma_sf_final = 1./0.05;
sigma_sf_final = sigma_sf_final*(sigma_sf_final_2*fabs(2.5875-W)+sigma_sf_final_1*fabs(2.5375-W));
};


 if ((W>=3.1375)&&(W<=4.5375)&&(Q2>0.00001)&&(Q2<1.3)) {
get_xsect_q2_13_wgt_3(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final );
//cout << sigma_l_final<< "\n";
};

if ((W>=3.0375)&&(W<=3.1375)&&(Q2>0.00001)&&(Q2<1.3)) {

get_xsect_25_30(Q2,3.0375, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_2, sigma_l_final_2, sigma_c2f_final_2,sigma_s2f_final_2,sigma_cf_final_2,sigma_sf_final_2);

get_xsect_q2_13_wgt_3(Q2, 3.1375, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_1, sigma_l_final_1, sigma_c2f_final_1,sigma_s2f_final_1,sigma_cf_final_1,sigma_sf_final_1 );


sigma_t_final = 1./0.1;
sigma_t_final = sigma_t_final*(sigma_t_final_2*fabs(3.1375-W)+sigma_t_final_1*fabs(3.0375-W));

sigma_l_final = 1./0.1;
sigma_l_final = sigma_l_final*(sigma_l_final_2*fabs(3.1375-W)+sigma_l_final_1*fabs(3.0375-W));

sigma_c2f_final = 1./0.1;
sigma_c2f_final = sigma_c2f_final*(sigma_c2f_final_2*fabs(3.1375-W)+sigma_c2f_final_1*fabs(3.0375-W));

sigma_s2f_final = 1./0.1;
sigma_s2f_final = sigma_s2f_final*(sigma_s2f_final_2*fabs(3.1375-W)+sigma_s2f_final_1*fabs(3.0375-W));

sigma_cf_final = 1./0.1;
sigma_cf_final = sigma_cf_final*(sigma_cf_final_2*fabs(3.1375-W)+sigma_cf_final_1*fabs(3.0375-W));

sigma_sf_final = 1./0.1;
sigma_sf_final = sigma_sf_final*(sigma_sf_final_2*fabs(3.1375-W)+sigma_sf_final_1*fabs(3.0375-W));
};



if (Q2nodata > 1.299){

sigma_t_final = sigma_t_final*Func_q2_dep(Q2nodata)/Func_q2_dep(1.299);
sigma_l_final = sigma_l_final*Func_q2_dep(Q2nodata)/Func_q2_dep(1.299);
sigma_c2f_final = sigma_c2f_final*Func_q2_dep(Q2nodata)/Func_q2_dep(1.299);
sigma_s2f_final = sigma_s2f_final*Func_q2_dep(Q2nodata)/Func_q2_dep(1.299);
sigma_cf_final = sigma_cf_final*Func_q2_dep(Q2nodata)/Func_q2_dep(1.299);
sigma_sf_final = sigma_sf_final*Func_q2_dep(Q2nodata)/Func_q2_dep(1.299);



};

Q2 = Q2nodata;

if (W < 1.2375){
sigma_t_final = 0.;
sigma_l_final = 0.;
sigma_c2f_final = 0.;
sigma_s2f_final = 0.;
sigma_cf_final = 0.;
sigma_sf_final = 0.;
};



if (W < 3.0125) interpol_int(Q2,W,xsect_int_test_t, xsect_int_test_l);

//calculating sigma_total from different sigmas, eps_l and eps_t

sigma_total =0.;

if ((isnan(eps_l))||(isnan(eps_t))) cout << eps_l<< " "<< eps_t<<" eps\n";

if  (!(eps_l>0.)&&!(eps_l<0)) eps_l = 0.;
if  (!(eps_t>0.)&&!(eps_t<0)) eps_t = 0.;

sigma_total = sigma_t_final;
sigma_total = sigma_total + eps_l*sigma_l_final;
sigma_total = sigma_total + eps_t*(sigma_c2f_final*cos(2.*ph_hadr) + sigma_s2f_final*sin(2.*ph_hadr));
sigma_total = sigma_total + sqrt(2.*eps_l*(eps_t+1))*(sigma_cf_final*cos(ph_hadr) + sigma_sf_final*sin(ph_hadr));

//if ((flag_radmod == 1)||(flag_radmod == 2)) {
//cout << eps_l<< " "<<eps_t<< " mp\n";

sigma_total = sigma_total*cr_rad_fact;

//};


//multiply sigma_total by virtual photon flux


if (flag_flux==1){

get_EpsL_Ebeam_ferm(E_beam,E_E_prime, phi_e,Theta_e_prime,eps_l_for_flux,eps_t_for_flux,Ebeam_ferm_for_flux,Ep_ferm_for_flux);
//cout<< E_beam<<" "<<E_E_prime<<" "<< Theta_e_prime<<" p[\n";

V_flux = alph_const/4./M_PI;
V_flux = V_flux/Ebeam_ferm_for_flux/Ebeam_ferm_for_flux/MP/MP;
V_flux = V_flux/(1.-eps_t_for_flux)/Q2_old;
V_flux = V_flux*W_ferm*(W_ferm*W_ferm-MP*MP);

sigma_total = sigma_total*V_flux;    
};

if ((isnan(abs(sigma_total)))||(isnan(abs(cr_rad_fact)))) cout<< sigma_total<<" "<<cr_rad_fact<<" oo\n";
//cout << W_old<< " "<<W<<" "<<px_fermi<<" "<<py_fermi<<" zzzzzzzzzzzzzzzz<<pz_fermi<<"" \n";
if (!(sigma_total>0.)&&!(sigma_total<0.)) cout << sigma_total <<" "<<cr_rad_fact<<" "<<V_flux<<" "<< W_old<< " "<<W_ferm<<" "<<W<<" "<<Q2_old<<" "<<Q2<< " uuuuuuuu\n";
if (!(cr_rad_fact>0.)&&!(cr_rad_fact<0.)) cout<< sigma_total<<" "<<cr_rad_fact<<" oo2\n";


    //FOR 1-PIM, 2-PIP, 3-P 
    
  //  cout <<"qqq1  "<<th_hadr <<"  "<< ph_hadr <<"  "<<alph_hadr<<  "\n";
  anti_rot(W, Q2, ph_e_qualab, E_beam_fermi, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr,  MPIM, MPIP, MP, P4_PIM, P4_PIP,  P4_Pfin);
  //get_rot2( E_beam,Ep, phi_e,theta_e,th_rot2, ph_rot2, ph_e_ferm);
  
 // if (flag_fermi ==1) 
 fermi_anti_rot(theta_rot2,phi_rot2,P4_PIM,P4_PIP, P4_Pfin,P4_E_prime_qulab_rad,P4_Eini_fermi); 
  
//cout << P4_E_prime_qulab_rad[0]<<" "<< P4_E_prime_qulab_rad[1]<<" "<< P4_E_prime_qulab_rad[2]<<" "<< P4_E_prime_qulab_rad[3]<<" 2\n";
   //FOR 1-P, 2-PIP, 3-PIM
     // anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr, MP,MPIP, MPIM, P4_Pfin,P4_PIP, P4_PIM);
   
     //FOR 1-PIP, 2-PIM, 3-P 
     //   anti_rot(W, Q2, phi_e, E_beam, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr, MPIP,MPIM, MP, P4_PIP,P4_PIM, P4_Pfin);
	
//-------------------------------------------------------------------------------------  
    

    //FOR 1-PIM, 2-PIP, 3-P
  //  rot(Q2, E_beam, P4_E_prime,P4_PIM, P4_PIP,  P4_Pfin, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
   
     //FOR 1-P, 2-PIP, 3-PIM
  //   rot(Q2, E_beam, P4_E_prime,P4_Pfin, P4_PIP,  P4_PIM, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
  
  //FOR 1-PIP, 2-PIM, 3-P
     //  rot(Q2, E_beam, P4_E_prime,P4_PIP, P4_PIM,  P4_Pfin, inv_m12, inv_m23, th_hadr, alph_hadr, ph_hadr);
  
   //  cout <<"qqq2  "<<th_hadr <<"  "<< ph_hadr <<"  "<<alph_hadr<<  "\n";

    
    P4_0_miss = P4_Pini + P4_Eini - P4_E_prime - P4_PIP - P4_Pfin - P4_PIM;
    P4_0_miss_2 = P4_Pini + P4_Eini_new - P4_E_prime_new - P4_PIP - P4_Pfin - P4_PIM;
    
    
        P4_pim_miss = P4_Pini + P4_Eini - P4_E_prime - P4_PIP - P4_Pfin;
        P4_pim_miss_2 = P4_Pini + P4_Eini_new - P4_E_prime_new - P4_PIP - P4_Pfin;
	
	//----fermi------
	
	P4_0_miss_fermi = P4_Pini + P4_Eini - P4_E_prime - P4_PIP - P4_Pfin - P4_PIM;
	P4_0_miss_fermi_2 = P4_Pini_fermi + P4_Eini_lab_rad - P4_E_prime_qulab_rad - P4_PIP - P4_Pfin - P4_PIM;
	
	 P4_pim_miss_fermi = P4_Pini + P4_Eini - P4_E_prime - P4_PIP - P4_Pfin;
	 P4_pim_miss_fermi_2 = P4_Pini_fermi + P4_Eini_lab_rad - P4_E_prime_qulab_rad - P4_PIP - P4_Pfin;
	
	
    
        //writing generated events into desired input file
	p_el_test = (P4_E_prime.Vect()).Mag();
	
	
	if (!(P4_E_prime.Mag()>0.)&&!(P4_E_prime.Mag()<0.)) cout << P4_E_prime.Mag() <<" "<< W<< " " << Q2<<" "<< sigma_total<< " tt\n";
	if (!(-P4_E_prime.Mag()>0.)&&!(-P4_E_prime.Mag()<0.)) cout << P4_E_prime.Mag() <<" "<< W<< " " << Q2<<" "<< sigma_total<< " ttr\n";
	
	
	//do not remove. this is the nan-check before writing all momenta to output file
	if ((isnan(P4_E_prime[0]))||(isnan(P4_E_prime[1]))||(isnan(P4_E_prime[2]))||(isnan(P4_E_prime[3]))) cout << P4_E_prime[0]<< " "<< P4_E_prime[1]<< " "<<P4_E_prime[2]<< " "<<P4_E_prime[3]<<" "<<W<<" "<<Q2<< " final electron is nan \n";
	
	if ((isnan(P4_Pfin[0]))||(isnan(P4_Pfin[1]))||(isnan(P4_Pfin[2]))||(isnan(P4_Pfin[3]))) cout << P4_Pfin[0]<< " "<< P4_Pfin[1]<< " "<<P4_Pfin[2]<< " "<<P4_Pfin[3]<<" "<<W<<" "<<Q2<< " final proton is nan \n";
	
	if ((isnan(P4_PIP[0]))||(isnan(P4_PIP[1]))||(isnan(P4_PIP[2]))||(isnan(P4_PIP[3]))) cout << P4_PIP[0]<< " "<< P4_PIP[1]<< " "<<P4_PIP[2]<< " "<<P4_PIP[3]<<" "<<W<<" "<<Q2<< " pip is nan \n";
	
	
	if ((isnan(P4_PIM[0]))||(isnan(P4_PIM[1]))||(isnan(P4_PIM[2]))||(isnan(P4_PIM[3]))) cout << P4_PIM[0]<< " "<< P4_PIM[1]<< " "<<P4_PIM[2]<< " "<<P4_PIM[3]<<" "<<W<<" "<<Q2<< " pim is nan \n";
	
	
	
   out_file_write(i,sigma_total, W, Q2,  P4_E_prime,P4_Pfin, P4_PIP,P4_PIM,z_EL,x_EL,y_EL);
   

    h_zel->Fill(z_EL,1.);


if ((W>1.3)&&(W<1.85)){
 h_W->Fill(W,sigma_total);
   };
//    if ((Q2_old>0.3)&&(Q2_old<0.4)) 
      h_W_2->Fill(W_old,sigma_total);
   
   
   h_eps_l->Fill(W,Q2,eps_l);
   
   h_0_miss->Fill(P4_0_miss.Mag2(),sigma_total);
   h_0_miss_2->Fill(P4_0_miss_2.Mag2(),sigma_total);
   
   
    h_pim_miss->Fill(P4_pim_miss.Mag2(),sigma_total);
   h_pim_miss_2->Fill(P4_pim_miss_2.Mag2(),sigma_total);

 h_0_miss_en->Fill(P4_0_miss[3],sigma_total);
  h_0_miss_en_2->Fill(P4_0_miss_2[3],sigma_total);
  
//-----for fermi  
   h_0_miss_fermi->Fill(P4_0_miss_fermi.Mag2(),sigma_total);
   h_0_miss_fermi_2->Fill(P4_0_miss_fermi_2.Mag2(),sigma_total);
   
   
    h_pim_miss_fermi->Fill(P4_pim_miss_fermi.Mag2(),sigma_total);
   h_pim_miss_fermi_2->Fill(P4_pim_miss_fermi_2.Mag2(),sigma_total);

 h_0_miss_en_fermi->Fill(P4_0_miss_fermi[3],sigma_total);
  h_0_miss_en_fermi_2->Fill(P4_0_miss_fermi_2[3],sigma_total);
  
 h_fermi_bonn->Fill(P4_0_miss_fermi.Vect().Mag(),sigma_total);
//    h_fermi_bonn->Fill(sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
	
  
   
   h_Q2->Fill(Q2,sigma_total);
   
    h_phi_e->Fill(phi_e,1.);
     h_Q2vsW->Fill(W_old,Q2_old,sigma_total);
     h_Q2vsW2->Fill(W_old,Q2_old,1.);
     h_Q2vsW_t->Fill(W,Q2,sigma_t_final);
     h_Q2vsW_l->Fill(W,Q2,sigma_l_final);
     
     
     
     
     if ((Q2>=0.00049)&&(Q2<=0.00051))  {
     h_int_crsect_t[0]->Fill(W,xsect_int_test_t);
     h_int_crsect_l[0]->Fill(W,xsect_int_test_l);
     
     };
     
     
     if ((Q2>=0.0005)&&(Q2<=1.2505))  {
     
     h_int_crsect_t[int((Q2-0.0005)/0.05)+1]->Fill(W,xsect_int_test_t);
     h_int_crsect_l[int((Q2-0.0005)/0.05)+1]->Fill(W,xsect_int_test_l);
     
     
     };
     
     
     if ((Q2>=1.29)&&(Q2<=1.3))  {
     h_int_crsect_t[26]->Fill(W,xsect_int_test_t);
     h_int_crsect_l[26]->Fill(W,xsect_int_test_l);
     
     };
     
     
     
     
    h_nu->Fill(nu,1.);
    h_dalitz->Fill(inv_m12,inv_m23,1.);
    
// if ((W>1.6)&&(W<1.625)&&(Q2>0.7)&&(Q2<0.8)) 
// h_inv_m12 ->Fill(s12,1.);
   //if ((W>1.6)&&(W<1.625)&&(Q2>0.7)&&(Q2<0.8))
    h_inv_m23 ->Fill(s23,1.); 
    h_th_hadr->Fill(th_hadr,1.); 
     for (Int_t ij=1;ij<=100;ij++){
     h_th_hadr_2 ->SetBinContent (ij, h_th_hadr->GetBinContent(ij)/(cos((h_th_hadr->GetBinLowEdge(ij)))-cos((h_th_hadr->GetBinLowEdge(ij)+h_th_hadr->GetBinWidth(ij)))));
     };
     
     h_ph_hadr ->Fill(ph_hadr,sigma_total); 
 
     h_alph_hadr  ->Fill(alph_hadr,1.); 
    
    
    
    if ((Q2>=0.1)&&(Q2<=1.3))  {
 
     h_odn_w_dep_t[int((Q2-0.1)/0.1)]->Fill(W,sigma_t_final);
     h_odn_w_dep_l[int((Q2-0.1)/0.1)]->Fill(W,eps_l*sigma_l_final);
     h_odn_w_dep_l2[int((Q2-0.1)/0.1)]->Fill(W,sigma_l_final);
     
     h_odn_w_dep_tot[int((Q2-0.1)/0.1)]->Fill(W,sigma_total); 
     

    
    };
    
 //   if ((Q2>=0.25)&&(Q2<=0.6))  {
 //    h_odn_w_dep_tot[int((Q2-0.25)/0.05)]->Fill(W,sigma_total); 
  //  };
    
    
    
    
     if ((W>1.25)&&(W<=2.075))  { 
     
     h_odn_q2_dep_t[int((W-1.25)/0.025)]->Fill(Q2,sigma_t_final);
     
     h_odn_q2_dep_l[int((W-1.25)/0.025)]->Fill(Q2,eps_l*sigma_l_final);
      h_odn_q2_dep_l2[int((W-1.25)/0.025)]->Fill(Q2,sigma_l_final);
     
     
     h_odn_q2_dep_tot[int((W-1.25)/0.025)]->Fill(Q2,sigma_total);
     
     
     
     };
    
    
    
    
  if ((W_old>=1.225)&&(W_old<=2.15)&&(Q2>=0.0001)&&(Q2<=1.4))  {
  
 // h_odn_q2_dep[int((W-1.425)/0.025)]->Fill(Q2,diff_xsect_ripani);
  
  h_odn_inv_m12[int((W_old-1.225)/0.025)]->Fill(inv_m12,sigma_total);
  h_odn_inv_m23[int((W_old-1.225)/0.025)]->Fill(inv_m23,sigma_total);
  h_odn_alpha[int((W_old-1.225)/0.025)]->Fill(alph_hadr,sigma_total);
 h_odn_theta[int((W_old-1.225)/0.025)]->Fill(th_hadr,sigma_total);

  };
  
   if ((W>=2.1625)&&(W<=3.0375)&&(Q2>=0.0001)&&(Q2<=1.3))  {
  

  
  h_odn_wwide_inv_m12[int((W-2.1625)/0.05)]->Fill(inv_m12,sigma_total);
  h_odn_wwide_inv_m23[int((W-2.1625)/0.05)]->Fill(inv_m23,sigma_total);
  h_odn_wwide_alpha[int((W-2.1625)/0.05)]->Fill(alph_hadr,sigma_total);
 h_odn_wwide_theta[int((W-2.1625)/0.05)]->Fill(th_hadr,sigma_total);

  };
  
  
  
  if ((W>=3.0875)&&(W<=4.5375)&&(Q2>=0.0001))  {
  

  
  h_odn_wgt3_inv_m12[int((W-3.0875)/0.1)]->Fill(inv_m12,sigma_total);
  h_odn_wgt3_inv_m23[int((W-3.0875)/0.1)]->Fill(inv_m23,sigma_total);
  h_odn_wgt3_alpha[int((W-3.0875)/0.1)]->Fill(alph_hadr,sigma_total);
  h_odn_wgt3_theta[int((W-3.0875)/0.1)]->Fill(th_hadr,sigma_total);

  };
  


       t21->Fill();  
    
    };//end event-loop
    
  
  
   
  for (Int_t ik=0;ik<=36;ik++){
  for (Int_t ij=1;ij<=100;ij++){
     h_odn_theta_2[ik] ->SetBinContent (ij, h_odn_theta[ik]->GetBinContent(ij)/(cos((h_odn_theta[ik]->GetBinLowEdge(ij)))-cos((h_odn_theta[ik]->GetBinLowEdge(ij)+h_odn_theta[ik]->GetBinWidth(ij)))));
     };
  };
  
  for (Int_t ik=0;ik<18;ik++){
  for (Int_t ij=1;ij<=100;ij++){
     h_odn_wwide_theta_2[ik] ->SetBinContent (ij, h_odn_wwide_theta[ik]->GetBinContent(ij)/(cos((h_odn_theta[ik]->GetBinLowEdge(ij)))-cos((h_odn_theta[ik]->GetBinLowEdge(ij)+h_odn_theta[ik]->GetBinWidth(ij)))));
     };
  };
  
  
  for (Int_t ik=0;ik<15;ik++){
  for (Int_t ij=1;ij<=100;ij++){
     h_odn_wgt3_theta_2[ik] ->SetBinContent (ij, h_odn_wgt3_theta[ik]->GetBinContent(ij)/(cos((h_odn_theta[ik]->GetBinLowEdge(ij)))-cos((h_odn_theta[ik]->GetBinLowEdge(ij)+h_odn_theta[ik]->GetBinWidth(ij)))));
     };
  };
    
    //closing output file
    out_file_close();
    
      TFile *outFile = new TFile("tree_sigma.root","recreate");
     outFile->cd();
     t21->Write("", TObject::kOverwrite);
     outFile->Write();
     outFile->Close();
     t21->Delete(); 
     
     
     TFile *file1 = TFile::Open("out_11Apr_test.root","RECREATE");
file1->cd();

h_zel->SetMinimum(0.);
h_W->SetMinimum(0.);
h_W_2->SetMinimum(0.);
h_Q2->SetMinimum(0.);
h_phi_e->SetMinimum(0.);
h_Q2vsW->SetMinimum(0.); 
h_nu->SetMinimum(0.);
h_dalitz->SetMinimum(0.);
 h_inv_m12 ->SetMinimum(0.);
 h_inv_m23 ->SetMinimum(0.); 
 h_th_hadr->SetMinimum(0.);
 h_th_hadr_2->SetMinimum(0.);
 h_ph_hadr ->SetMinimum(0.);
 h_alph_hadr ->SetMinimum(0.); 


h_zel->Write("", TObject::kOverwrite);
h_W->Write("", TObject::kOverwrite);
h_W_2->Write("", TObject::kOverwrite);
h_Q2->Write("", TObject::kOverwrite);
h_phi_e->Write("", TObject::kOverwrite);
h_Q2vsW->Write("", TObject::kOverwrite);
h_Q2vsW2->Write("", TObject::kOverwrite);

h_0_miss->Write("", TObject::kOverwrite);
h_0_miss_2->Write("", TObject::kOverwrite);
h_pim_miss->Write("", TObject::kOverwrite);
h_pim_miss_2->Write("", TObject::kOverwrite);
 h_0_miss_en->Write("", TObject::kOverwrite);
 h_0_miss_en_2->Write("", TObject::kOverwrite);
 
h_eradgam->Write("", TObject::kOverwrite);


h_0_miss_fermi->Write("", TObject::kOverwrite);
h_0_miss_fermi_2->Write("", TObject::kOverwrite);
h_pim_miss_fermi->Write("", TObject::kOverwrite);
h_pim_miss_fermi_2->Write("", TObject::kOverwrite);
 h_0_miss_en_fermi->Write("", TObject::kOverwrite);
 h_0_miss_en_fermi_2->Write("", TObject::kOverwrite);


h_fermi_bonn->Write("", TObject::kOverwrite);



h_eps_l->Write("", TObject::kOverwrite);
h_nu->Write("", TObject::kOverwrite);
h_dalitz->Write("", TObject::kOverwrite);

  h_inv_m12->Write("", TObject::kOverwrite); 
    h_inv_m23->Write("", TObject::kOverwrite);  
     h_th_hadr->Write("", TObject::kOverwrite);
      h_th_hadr_2->Write("", TObject::kOverwrite);
     h_ph_hadr ->Write("", TObject::kOverwrite);
     h_alph_hadr->Write("", TObject::kOverwrite);  
     
     
for (Short_t ii=0; ii<=36; ii++) {
h_odn_inv_m12[ii]->Write("", TObject::kOverwrite); 
};

for (Short_t ii=0; ii<=36; ii++) {
h_odn_inv_m23[ii]->Write("", TObject::kOverwrite); 
};

for (Short_t ii=0; ii<=36; ii++) {
h_odn_alpha[ii]->Write("", TObject::kOverwrite); 
};

for (Short_t ii=0; ii<=36; ii++) {
h_odn_theta[ii]->Write("", TObject::kOverwrite); 
};

for (Short_t ii=0; ii<=36; ii++) {
h_odn_theta_2[ii]->Write("", TObject::kOverwrite); 
};


for (Short_t ii=0; ii<18; ii++) {
h_odn_wwide_inv_m12[ii]->Write("", TObject::kOverwrite); 
};

for (Short_t ii=0; ii<18; ii++) {
h_odn_wwide_inv_m23[ii]->Write("", TObject::kOverwrite); 
};

for (Short_t ii=0; ii<18; ii++) {
h_odn_wwide_alpha[ii]->Write("", TObject::kOverwrite); 
};

for (Short_t ii=0; ii<18; ii++) {
h_odn_wwide_theta[ii]->Write("", TObject::kOverwrite); 
};

for (Short_t ii=0; ii<18; ii++) {
h_odn_wwide_theta_2[ii]->Write("", TObject::kOverwrite); 
};


for (Short_t ii=0; ii<15; ii++) {
h_odn_wgt3_inv_m12[ii]->Write("", TObject::kOverwrite); 
};

for (Short_t ii=0; ii<15; ii++) {
h_odn_wgt3_inv_m23[ii]->Write("", TObject::kOverwrite); 
};

for (Short_t ii=0; ii<15; ii++) {
h_odn_wgt3_alpha[ii]->Write("", TObject::kOverwrite); 
};

for (Short_t ii=0; ii<15; ii++) {
h_odn_wgt3_theta[ii]->Write("", TObject::kOverwrite); 
};

for (Short_t ii=0; ii<15; ii++) {
h_odn_wgt3_theta_2[ii]->Write("", TObject::kOverwrite); 
};



for (Short_t ii=0; ii<12; ii++) {
h_odn_w_dep_tot[ii]->SetMinimum(0.);
h_odn_w_dep_t[ii]->SetMinimum(0.);
h_odn_w_dep_l[ii]->SetMinimum(0.);
h_odn_w_dep_l2[ii]->SetMinimum(0.);

h_odn_w_dep_t[ii]->Write("", TObject::kOverwrite); 
h_odn_w_dep_l[ii]->Write("", TObject::kOverwrite); 
h_odn_w_dep_l2[ii]->Write("", TObject::kOverwrite); 
h_odn_w_dep_tot[ii]->Write("", TObject::kOverwrite);


};


for (Short_t ii=0; ii<33; ii++) {

h_odn_q2_dep_t[ii]->SetMinimum(0.);
h_odn_q2_dep_l[ii]->SetMinimum(0.);
h_odn_q2_dep_l2[ii]->SetMinimum(0.);
h_odn_q2_dep_tot[ii]->SetMinimum(0.);


h_odn_q2_dep_t[ii]->Write("", TObject::kOverwrite); 
h_odn_q2_dep_l[ii]->Write("", TObject::kOverwrite); 
h_odn_q2_dep_l2[ii]->Write("", TObject::kOverwrite); 
h_odn_q2_dep_tot[ii]->Write("", TObject::kOverwrite); 

};



for (Short_t ii=0; ii<27; ii++) {
h_int_crsect_t[ii]->Scale(1./h_int_crsect_t[ii]->GetEntries()*71.);
h_int_crsect_l[ii]->Scale(1./h_int_crsect_l[ii]->GetEntries()*71.);
   
 };  
   
   
for (Short_t ii=0; ii<27; ii++) {
h_int_crsect_t[ii]->Write("", TObject::kOverwrite); 
h_int_crsect_l[ii]->Write("", TObject::kOverwrite); 

};


//for (Short_t ii=0; ii<15; ii++) {
//h_odn_q2_dep[ii]->Write("", TObject::kOverwrite); 
//};
file1->Write(); 
    
// h_Q2vsW->Scale(1./Nevents*71.*25.);
// h_Q2vsW_t->Scale(1./Nevents*71.*25.);
// h_Q2vsW_l->Scale(1./Nevents*71.*25.);
 
 
 
/*  
 
  for (Short_t j=0;j<=26;j++){
  
 qqq.str("");
   //  qqq << "Q2_" << 0.0005+0.05*j <<  ".dat";
     qqq << "intsec_q2_" <<j <<  ".dat";
     
     std::ofstream ofs (qqq.str().c_str(), std::ofstream::out);
     qqq.str(""); 
  
 for (Short_t i=1;i<=71;i++){
 




   ofs  << 1.2375 +0.025*i<< "\n";    
  ofs <<h_int_crsect_t[j]->GetBinContent(i) << "\n";  
  ofs <<h_int_crsect_l[j]->GetBinContent(i) << "\n"; 

     
 };
 
 };

 */
 
 
  
    TApplication *theApp = new TApplication("App", &argc, argv);

  


  delete [] file;
  file = NULL;
  
  return 0;
}
