#include "TFile.h"
#include "TMath.h"
#include <stdio.h>
#include <dlfcn.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <TLorentzVector.h>
#include "global.h"
#include <fstream>
#include <iomanip>
 using namespace std; 
 


void out_file_write(int i,Float_t sigma, Float_t W, Float_t Q2, TLorentzVector &P4_E, TLorentzVector &P4_1, TLorentzVector &P4_2, TLorentzVector &P4_3, Float_t z_EL, Float_t x_EL, Float_t y_EL) {
 
 
// cout <<sigma<< " qqq2 \n";
//cout << i<<"\n";
 //test
// int i=10;
// Float_t z_EL=2.2;
 
 

 

 
 
 //------CREATE BOS OUTPUT WITH DIFFERENT BANKS IF NEEDED
 
 
 if ((flag_bos == 1)||(flag_bos == 2)){
 
 //Create HEAD-bank anyway
      clasHEAD_t* HEAD ;
  if ( ( HEAD = (clasHEAD_t*)makeBank( &bcs_, "HEAD", 0, sizeof(head_t)/4, 4 ) ) )
    {
      HEAD->head[0].version = 1; 
      HEAD->head[0].nevent = i ;
      HEAD->head[0].nrun = 1;
      HEAD->head[0].evtclass = 7 ;
      HEAD->head[0].type = -2 ;
      HEAD->head[0].time = 0;
      HEAD->head[0].trigbits = 1 ;
      
    
    }
    
   clasMCEV_t* MCEV ;
  if ( ( MCEV = (clasMCEV_t*)makeBank( &bcs_, "MCEV", 0, sizeof(mcev_t)/4, 1 ) ) )
    {
      MCEV->mcev[0].i1 = 100000*((float) rand() / (RAND_MAX)); 
      MCEV->mcev[0].i2 =100000*((float) rand() / (RAND_MAX)) ;
      
      
    
    } 
    
//-----------------------------------------------------------------------    
  //Create MCTK & MCVX banks if needed  
    if (flag_bos==1) {
    
       clasMCTK_t* MCTK ;
  if ( ( MCTK = (clasMCTK_t*)makeBank( &bcs_, "MCTK", 0, sizeof(mctk_t)/4, 4 ) ) )
    {
    
    //electron
      MCTK->mctk[0].cx = P4_E[0]/((P4_E.Vect()).Mag()); 
      MCTK->mctk[0].cy = P4_E[1]/((P4_E.Vect()).Mag()); 
      MCTK->mctk[0].cz = P4_E[2]/((P4_E.Vect()).Mag()); 
      MCTK->mctk[0].pmom = (P4_E.Vect()).Mag();
      MCTK->mctk[0].mass = 0.001;
      MCTK->mctk[0].charge = -1.;
      MCTK->mctk[0].id = 11;
      MCTK->mctk[0].flag = 0;
      MCTK->mctk[0].beg_vtx = 1;
      MCTK->mctk[0].end_vtx = 0;
      MCTK->mctk[0].parent = 0;
      
      //pi-
     
      MCTK->mctk[1].cx = P4_3[0]/((P4_3.Vect()).Mag()); 
      MCTK->mctk[1].cy = P4_3[1]/((P4_3.Vect()).Mag()); 
      MCTK->mctk[1].cz = P4_3[2]/((P4_3.Vect()).Mag()); 
      MCTK->mctk[1].pmom = (P4_3.Vect()).Mag();
      MCTK->mctk[1].mass = MPIM;
      MCTK->mctk[1].charge = -1.;
      MCTK->mctk[1].id = -211;
      MCTK->mctk[1].flag = 0;
      MCTK->mctk[1].beg_vtx = 1;
      MCTK->mctk[1].end_vtx = 0;
      MCTK->mctk[1].parent = 0;
      
      //pi+
     
      MCTK->mctk[2].cx = P4_2[0]/((P4_2.Vect()).Mag()); 
      MCTK->mctk[2].cy = P4_2[1]/((P4_2.Vect()).Mag()); 
      MCTK->mctk[2].cz = P4_2[2]/((P4_2.Vect()).Mag()); 
      MCTK->mctk[2].pmom = (P4_2.Vect()).Mag();
      MCTK->mctk[2].mass = MPIP;
      MCTK->mctk[2].charge = 1.;
      MCTK->mctk[2].id = 211;
      MCTK->mctk[2].flag = 0;
      MCTK->mctk[2].beg_vtx = 1;
      MCTK->mctk[2].end_vtx = 0;
      MCTK->mctk[2].parent = 0;
      
     //proton
     
      MCTK->mctk[3].cx = P4_1[0]/((P4_1.Vect()).Mag()); 
      MCTK->mctk[3].cy = P4_1[1]/((P4_1.Vect()).Mag()); 
      MCTK->mctk[3].cz = P4_1[2]/((P4_1.Vect()).Mag()); 
      MCTK->mctk[3].pmom = (P4_1.Vect()).Mag();
      MCTK->mctk[3].mass = MP;
      MCTK->mctk[3].charge = 1.;
      MCTK->mctk[3].id = 2212;
      MCTK->mctk[3].flag = 0;
      MCTK->mctk[3].beg_vtx = 1;
      MCTK->mctk[3].end_vtx = 0;
      MCTK->mctk[3].parent = 0;
      
      
       
      
      
      
      
      
      
      
    } //end MCTK bank creation 

  clasMCVX_t* MCVX ;
  if ( ( MCVX = (clasMCVX_t*)makeBank( &bcs_, "MCVX", 0, sizeof(mcvx_t)/4, 1 ) ) )
    {
      MCVX->mcvx[0].x = x_EL; 
      MCVX->mcvx[0].y = y_EL;
      MCVX->mcvx[0].z = z_EL;
      MCVX->mcvx[0].tof = 0.;
      MCVX->mcvx[0].flag = 0;
    }    //end MCVX bank creation   
};//end if flag_bos==1

//-------------------------------------------------

//Create PART bank if needed 
    if (flag_bos==2) {

 clasPART_t* PART ;
  if ( ( PART = (clasPART_t*)makeBank( &bcs_, "PART", 0, sizeof(part_t)/4, 4) ) )
    {
    
    //electron
    PART->part[0].pid = 3.; 
    PART->part[0].vert.x = 0.;
    PART->part[0].vert.y = 0.;
    PART->part[0].vert.z = z_EL;
    PART->part[0].p.t = P4_E[3];
    PART->part[0].p.space.x = P4_E[0];
    PART->part[0].p.space.y = P4_E[1];
    PART->part[0].p.space.z =P4_E[2] ;
    PART->part[0].q = -1.;
    PART->part[0].trkid = 0.;
    PART->part[0].qpid = 0.;
    PART->part[0].qtrk = 0.;
    PART->part[0].flags = 0.;
       
       //proton
       
    PART->part[1].pid = 14.; 
    PART->part[1].vert.x = 0.;
    PART->part[1].vert.y = 0.;
    PART->part[1].vert.z = 0.;
    PART->part[1].p.t = P4_1[3];
    PART->part[1].p.space.x = P4_1[0];
    PART->part[1].p.space.y = P4_1[1];
    PART->part[1].p.space.z =P4_1[2] ;
    PART->part[1].q = 1.;
    PART->part[1].trkid = 0.;
    PART->part[1].qpid = 0.;
    PART->part[1].qtrk = 0.;
    PART->part[1].flags = 0.;   
       
           //pi+
       
    PART->part[2].pid = 8.; 
    PART->part[2].vert.x = 0.;
    PART->part[2].vert.y = 0.;
    PART->part[2].vert.z = 0.;
    PART->part[2].p.t = P4_2[3];
    PART->part[2].p.space.x = P4_2[0];
    PART->part[2].p.space.y = P4_2[1];
    PART->part[2].p.space.z =P4_2[2] ;
    PART->part[2].q = 1.;
    PART->part[2].trkid = 0.;
    PART->part[2].qpid = 0.;
    PART->part[2].qtrk = 0.;
    PART->part[2].flags = 0.;  
    
    
      //pi-
       
    PART->part[3].pid = 9.; 
    PART->part[3].vert.x = 0.;
    PART->part[3].vert.y = 0.;
    PART->part[3].vert.z = 0.;
    PART->part[3].p.t = P4_3[3];
    PART->part[3].p.space.x = P4_3[0];
    PART->part[3].p.space.y = P4_3[1];
    PART->part[3].p.space.z =P4_3[2] ;
    PART->part[3].q = -1.;
    PART->part[3].trkid = 0.;
    PART->part[3].qpid = 0.;
    PART->part[3].qtrk = 0.;
    PART->part[3].flags = 0.; 
       
       
            
    }  //end PART bank creation     
};//end if flag_bos==2

 putBOS( &bcs_, 1, "E" );
      dropAllBanks( &bcs_, "E");
      cleanBanks( &bcs_ );
};//end if ((flag_bos == 1)||(flag_bos == 2))

      
   if (flag_lund == 1){   
   
      out_lund_stream <<"4  1.  1.  0  0  "<<i<< "  0  "<<std::fixed<<std::setprecision(6)<< W<< "  "<< Q2<< "  "<<sigma << "\n"; 
      
      //electron
       out_lund_stream <<" 1  -1.  1  11  0  0  "<<P4_E[0]<<"  ";
       out_lund_stream <<P4_E[1]<<"  "<<P4_E[2]<<"  "<<P4_E[3]<<"  ";
       out_lund_stream <<"0.0005  0.0000  0.0000  "<<z_EL<<"\n";
       
       
       //proton
       out_lund_stream <<" 2  1.  1  2212  0  0  "<<P4_1[0]<<"  ";
       out_lund_stream <<P4_1[1]<<"  "<<P4_1[2]<<"  "<<P4_1[3]<<"  ";
       out_lund_stream <<"0.9383  0.0000  0.0000  "<<z_EL<<"\n";
       
       
        //pi+
       out_lund_stream <<" 3  1.  1  211  0  0  "<<P4_2[0]<<"  ";
       out_lund_stream <<P4_2[1]<<"  "<<P4_2[2]<<"  "<<P4_2[3]<<"  ";
       out_lund_stream <<"0.1396  0.0000  0.0000  "<<z_EL<<"\n";
       
       
        //pi-
       out_lund_stream <<" 4  -1.  1  -211  0  0  "<<P4_3[0]<<"  ";
       out_lund_stream <<P4_3[1]<<"  "<<P4_3[2]<<"  "<<P4_3[3]<<"  ";
       out_lund_stream <<"0.1396  0.0000  0.0000  "<<z_EL<<"\n";
       
       
       
    };       

 };
 
