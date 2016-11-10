#include "TMath.h"
#include <stdio.h>
#include <dlfcn.h>
#include <TGClient.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"

#include <fstream>
 using namespace std;

int out_file_close() {
//Close BOS output if exist

 if ((flag_bos == 1)||(flag_bos == 2)){
close_fpack_unit("BOSOUTPUT");
     
  };  
  
 
   if (flag_lund == 1){
  
   out_lund_stream.close();
   };
     };
