#ifndef GET_XSECT_Q2_TEST_H
#define GET_XSECT_Q2_TEST_H



Short_t getWbin_rip2_test (Float_t W); 
Short_t getWbin_rip3_test (Float_t W); 

Short_t getsbin_rip3 (Short_t Wbin, Float_t sgen, Float_t Smax, Float_t Smin);
void get_xsect_q2_130_test(Float_t Q2gen, Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Float_t phigen,  Float_t &sigma_t_final, Float_t &sigma_l_final,Float_t  &sigma_c2f_final,Float_t  &sigma_s2f_final,Float_t &sigma_cf_final,Float_t  &sigma_sf_final);

#endif
