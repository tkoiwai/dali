#include <stdio.h>

#include "/home/apollo/anarootSAMURAI/sources/Core/include/TArtStoreManager.hh"
#include "/home/apollo/anarootSAMURAI/sources/Core/include/TArtEventStore.hh"
#include "TArtDALIParameters.hh"
#include "TArtCalibDALI.hh"
#include "TArtDALINaI.hh"
#include "TSystem.h"
#include "TH2F.h"
#include "TFile.h"
#include "TClonesArray.h"

#include "iostream"
#include <TRandom.h>
#include "signal.h"

using namespace std;
// function to exit loop at keyboard interrupt. 
bool stoploop = false;
void stop_interrupt()
{
  printf("keyboard interrupt\n");
  stoploop = true;
}

//void MakeDALIHist(char *infile, char *outfile="dalia.root"){
int main(int argc, char *argv[]){
  //gSystem->Load("libXMLParser.so");
  //gSystem->Load("libanaroot.so");

  char* infile = argv[1];
  char* outfile = argv[2];
  char* runnum = argv[3];
  
  TArtStoreManager * sman = TArtStoreManager::Instance();
  TArtEventStore *estore = new TArtEventStore();
  estore->SetInterrupt(&stoploop);
  estore->Open(infile);
  if(estore->Open(infile))
    cout << "infile: " << infile << endl;

  TArtDALIParameters *dpara = TArtDALIParameters::Instance();
  dpara->LoadParameter("../db/DALI.xml");
  TArtCalibDALI *dalicalib= new TArtCalibDALI();

  TFile *fout = new TFile(outfile,"RECREATE");
  TH2F *h101 = new TH2F("h101",runnum,226,0,226,500,0,1000);
  TH2F *h102 = new TH2F("h102",runnum,226,0,226,400,0,4000);
  TH1F *hmulti = new TH1F("hmulti","",10,0,10);

  // define data nodes which are supposed to be dumped to tree 
  TClonesArray * info_array = (TClonesArray *)sman->FindDataContainer("EventInfo");
  std::cout<<info_array->GetName()<<std::endl;

  TClonesArray * dali_array=
     (TClonesArray *)sman->FindDataContainer("DALINaI");

  Int_t dalimultwithoutt = 0;
  Int_t dalimult = 0;
  Int_t dalitimetruemult = 0;
  Int_t dalimultthres = 0;
  Int_t dalitimetruemultthres = 0;

  Float_t fADC,fEnergy;
  TRandom *rand = new TRandom();
  Int_t neve = 0;
  Int_t counter = 0;
  while(estore->GetNextEvent()){
  //while(neve<1001){
    if(neve%10000==0)
      std::cout << "event: " << neve << std::endl;

    dalicalib->ClearData();
    dalicalib->ReconstructData();
    //    dalimultwithoutt = dalicalib->GetMultWithoutT();
    dalimult = dalicalib->GetMult();
    //    dalitimetruemult = dalicalib->GetTimeTrueMult();
    //    dalimultthres = dalicalib->GetMultThres();
    //   dalitimetruemultthres = dalicalib->GetTimeTrueMultThres();
    neve ++;
    hmulti->Fill(dalimult);
    if(dalimult<2) {;} else continue;

    Int_t NumDALI = dalicalib->GetNumNaI();
    for(Int_t j=0;j<NumDALI;j++){
      TArtDALINaI *nai = (TArtDALINaI *)dalicalib->GetNaI(j);
      if(nai->GetID()!=187 && nai->GetID()!=188){
	fADC=nai->GetRawADC()+(Float_t)rand->Rndm();
	fEnergy=nai->GetEnergy();
	h101->Fill((int)nai->GetID(),fADC);
	h102->Fill((int)nai->GetID(),fEnergy);
	counter++;
      }
    }
  }
  cout << "exited loop" << endl;
  fout->Write();
  cout << "file wrote" << endl;
  fout->Close("R");
  cout << "file closed" << endl;
  estore->ClearData();
  cout << "estore cleared" << endl;
  delete dalicalib;
  cout << "dalicalib deleted" << endl;
  //delete dpara;
  //cout << "dpara deleted" << endl;
  delete sman;
  cout << "sman deleted" << endl;
  cout << "Conversion done." << endl;
  cout << "counter = " << counter << endl;
  return 0;
}

