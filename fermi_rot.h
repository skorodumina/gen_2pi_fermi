#ifndef FERMI_ROT_H
#define FERMI_ROT_H


void fermi_rot(Float_t &E_beam_fermi, Float_t &theta_rot2,Float_t &phi_rot2, Float_t E_beam,TLorentzVector P4_E_prime,TLorentzVector  &P4_E_prime_boosted);

  void get_EpsL_Ebeam_ferm( Float_t E_beam, Float_t Ep, Float_t phi_e, Float_t theta_e, Float_t &eps_l, Float_t &eps_t,Float_t &E_beam_fermi,Float_t &E_p_fermi); 
  
  
   void get_rot2( Float_t E_beam,Float_t Ep, Float_t phi_e,Float_t theta_e,Float_t &th_rot2, Float_t &ph_rot2, Float_t &ph_e_ferm);
   
 void from_qulab_to_lab(Float_t theta_rot2,Float_t phi_rot2,Float_t theta_e_ferm, Float_t ph_e_ferm,Float_t Es,Float_t Ep,Float_t &Es_lab,Float_t &Ep_lab);  
#endif
