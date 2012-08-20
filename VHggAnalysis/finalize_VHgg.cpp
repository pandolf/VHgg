#include "RedNtpFinalizer_VHgg.h"
#include "TMath.h"
#include <iostream>








int main( int argc, char* argv[] ) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./finalize_VHgg [dataset] [selectionType] [bTaggerType=\"TCHE\"]" <<std::endl;
    return 13;
  }


  std::string dataset(argv[1]);
  std::string selectionType(argv[2]);

  std::string bTaggerType="TCHE";
  if( argc>3 ) {
    std::string bTaggerType_str(argv[3]);
    bTaggerType = bTaggerType_str;
  }




  RedNtpFinalizer_VHgg* nf = new RedNtpFinalizer_VHgg( dataset, selectionType, bTaggerType );
  
  } else {
  
    nf->addFile( dataset );

  }

  nf->finalize();


  return 0;

}



