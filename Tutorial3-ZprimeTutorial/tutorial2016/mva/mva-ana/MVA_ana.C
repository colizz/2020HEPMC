#include "MVA_ana.h"
#include "MVA_ana.cxx"

void MVA_ana()
{

MVA_ana train;
MVA_ana applica;

train.Classification("","BDTG","zprimetraining.root");

applica.ClassApplication("","BDTG","zprimetraining.root","selected.root","result.root");


exit(0);


}

