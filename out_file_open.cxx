#include "TMath.h"
#include <stdio.h>
#include <dlfcn.h>
#include <TGClient.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"

#include <sys/stat.h>

#include <sstream>
#include <fstream>
 using namespace std;

int out_file_open() {

//Create BOS output if needed

 if ((flag_bos == 1)||(flag_bos == 2)){
 
  remove(out_bos_file.c_str());
 
char mess[256];
      initbos();
      bankList( &bcs_, "E=", "HEADPARTMCTKMCVXMCEV" );
     sprintf( mess, "OPEN BOSOUTPUT UNIT=1 FILE=\"%s\" WRITE STATUS=NEW RECL=3600", out_bos_file.c_str() );
     fparm_c( mess );
};//end of BOS output flag check
     
//Create BOS output if needed 
  
  if (flag_lund == 1){ 
  
//  if (std::filesystem::exists(out_lund_file) == 0){
    remove(out_lund_file.c_str());
//    cout << "EXISTING LUND-FILE " << out_lund_file.c_str()<<" WAS RECREATED \n";
//  };
 out_lund_stream.open(out_lund_file.c_str(), std::ofstream::out);
  
  
  
  
  };
  
};
