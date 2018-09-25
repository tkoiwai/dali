#include <TSystem.h>
#include <string>
#include <iostream>
#include <sstream> //string stream
#include <fstream>
#include <sys/stat.h> //get the status of files. "st_"
#include <sys/types.h>
#include <unistd.h> //UNIx STanDard Header file
#include <climits> //char limits

#include "TTree.h"
#include "TFile.h"
#include "THashList.h"
#include "TEnv.h" //related to read file like geometry
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h" //moji toka mojiretsu
#include "TVector3.h" //3 vector
#include "TVectorD.h" //?
#include <vector>
#include "TMath.h"
#include "TCutG.h"

//ANAROOT headers
#include "TArtStoreManager.hh" //I/O of data and parameter
#include "TArtEventStore.hh" //convert ridf to raw data
#include "TArtEventInfo.hh"
#include "TArtRawEventObject.hh"
#include "TArtCore.hh" //debug, info, error, warning etc.

#include "TArtBigRIPSParameters.hh"
#include "TArtSAMURAIParameters.hh"
#include "TArtCalibCoin.hh"
#include "TArtRecoRIPS.hh" //added for delta. NOV14
#include "TArtRIPS.hh"
#include "TArtRIPSPara.hh"


//DALI
#include "TArtCalibDALI.hh"
#include "TArtDALINaIPara.hh"
#include "TArtDALINaI.hh"
#include "TArtDALIParameters.hh"

#define WriteOneEnvFile(file) file->Write(TString(TString(file->GetRcName()).ReplaceAll("/","_")).ReplaceAll(".","_")); // / toka . no youna moji wo _ ni replace surudake.
#define WriteAllEnvFiles WriteOneEnvFile(env_par); WriteOneEnvFile(env_geo); //WriteOneEnvFile(env_nebt0); WriteOneEnvFile(env_neut0);

//=====External Function defined at last=======================================
void SortDaliHit(Int_t, Int_t,vector <Int_t> *,vector <Double_t> *, vector <Double_t> *, vector <Double_t> *, vector <Double_t> *);
Double_t DopplerCorrection(Double_t, Double_t, Double_t);
inline bool exists_test(const std::string&);
inline bool exists_test(const TString&);
//=====main Function==========================================================
int main(int argc, char *argv[]){
  Long64_t MaxEventNumber = LLONG_MAX; //signed 8bit int. LLOMG_MAX ~ 2^63-1(const.)
  
  Int_t FileNumber = TString(argv[1]).Atoi();
  TString RidfFileName;

  if(FileNumber==0){
    std::cerr <<  " You should provide either a runnumber" << endl;
  }

  if(FileNumber>1&&FileNumber<36)
    RidfFileName = Form("/home/koiwai/analysis/ridf/smggdaq04/DALI2CalibRun%04d.ridf.gz",FileNumber);  
  else if(FileNumber>=36&&FileNumber<231)
    RidfFileName = Form("/home/koiwai/analysis/ridf/sdaq02/run%04d.ridf.gz",FileNumber);
  else{
    std::cerr << "run number is not correct." << endl;
  }
  
  //=====Load setting parameters=========================================
  TEnv *env_par = new TEnv("./conversion_settings.prm");
  TEnv *env_geo = new TEnv(env_par->GetValue("geometrydata","")); //unit of length:mm
  Double_t Dist_F5F7 = env_geo->GetValue("Dist_F5F7",0.0);
  Double_t Dist_F7F13 = env_geo->GetValue("Dist_F7F13",0.0);
  Double_t BDC1_Width = env_geo->GetValue("BDC1_Width",0.0); //mm
  Double_t FDC1_Width = env_geo->GetValue("FDC1_Width",0.0); //mm
  Double_t Dist_BDC1BDC2 = env_geo->GetValue("Dist_BDC1BDC2",0.0); //mm
  Double_t Dist_BDC1Target = env_geo->GetValue("Dist_BDC1Target",0.0); //mm
  Double_t Dist_BDC1FDC1 = env_geo->GetValue("Dist_BDC1FDC1",0.0); //Distance between the middle of BDC1 and the middle of FDC1 mm
  Double_t Dist_SBTTarget = env_geo->GetValue("Dist_SBTTarget",0.0); //mm
  //=====Load ANAROOT parameters===========================================
  TArtStoreManager *sman = TArtStoreManager::Instance();
  TArtBigRIPSParameters *BigRIPSPara = TArtBigRIPSParameters::Instance();
  TArtSAMURAIParameters *SAMURAIPara = TArtSAMURAIParameters::Instance();
  TArtDALIParameters *DALIPara = TArtDALIParameters::Instance();
  
  TObjString *xmlfile;

  TIter next(TString(env_par->GetValue("BigRIPSPara","")).Tokenize(" ")); //Make array of BigRIPSPara xmlfiles in conversion_settings.prm
  while((xmlfile = (TObjString *)next())){
    BigRIPSPara->LoadParameter((char*)xmlfile->GetName()); //load the xmlfile(?)
  }

  next = TString(env_par->GetValue("SamuraiPara","")).Tokenize(" ");
  while((xmlfile = (TObjString*)next())){
    SAMURAIPara->LoadParameter((char*)xmlfile->GetName());
  }

  next = TString(env_par->GetValue("DALIPara","")).Tokenize(" ");
  while(( xmlfile = (TObjString *)next())) {
      DALIPara->LoadParameter((char*)xmlfile->GetName());
  }
  

  //Load detectors.
  TArtCalibDALI *CalibDALI = new TArtCalibDALI();
  
  TArtRecoRIPS *RecoRIPS = new TArtRecoRIPS(); //added NOV14
  
  //=====Conversion ridf -> raw data======================================
  TArtEventStore *EventStore = new TArtEventStore;

  if(!EventStore->Open(RidfFileName)){
    std::cerr << "cannot open " << RidfFileName << std::endl << "aborting..." << std::endl;
    return 1;
  }

  //===== Load Beam file for PID cut =====
  TString infnameB = Form("/home/koiwai/analysis/rootfiles/ana/beam/ana_beam%04d.root",FileNumber);
  TFile   *infileB = TFile::Open(infnameB);
  TTree   *anatrB;
  infileB->GetObject("anatrB",anatrB);

  Double_t zetBR, aoqBR;
  Double_t betaF7F13, betaF3F13, bammaF7F13, gammaF3F13;

  anatrB->SetBranchAddress("zetBR",&zetBR);
  anatrB->SetBranchAddress("aoqBR",&aoqBR);
  anatrB->SetBranchAddress("betaF7F13",&betaF7F13);
  anatrB->SetBranchAddress("betaF3F13",&betaF3F13);
  anatrB->SetBranchAddress("gammaF7F13",&gammaF7F13);
  anatrB->SetBranchAddress("gammaF3F13",&gammaF3F13);

  //===== Load smri file for PID cut =====
  TString infnameS = Form("/home/koiwai/analysis/rootfiles/ana/smri/ana_smri%04d.root",FileNumber);
  TFile   *infileS = TFile::Open(infnameS);
  TTree   *anatrS;
  infileS->GetObject("anatrS",anatrS);

  Double_t zetSA, aoqSA;
  Double_t beta_minoshodo, gamma_minoshodo;

  anatrS->SetBranchAddress("zetSA",&zetSA);
  anatrS->SetBranchAddress("aoqSA",&aoqSA);
  anatrS->SetBranchAddress("beta_minoshodo",&beta_minoshodo);
  anatrS->SetBranchAddress("gamma_minoshodo",&gamma_minoshodo);

  //===== Load MINOS file =====
  TString infnameM = Form("/home/koiwai/analysis/rootfiles/minos/vertex/vertex%04d.root",FileNumber);
  TFile   *infileM = TFile::Open(infnameM);
  TTree   *trM;
  infileM->GetObject("tr",trM);
  
  Double_t vertexX, vertexY, vertexZ, vertexTheta, vertexPhi, theta2p;

  trM->SetBranchAddress("vertexX",&vertexX);
  trM->SetBranchAddress("vertexY",&vertexY);
  trM->SetBranchAddress("vertexZ",&vertexZ);
  trM->SetBranchAddress("vertexTheta",&vertexTheta);
  trM->SetBranchAddress("vertexPhi",&vertexPhi);
  trM->SetBranchAddress("theta2p",&theta2p);
  
  anatrB->AddFriend(anatrS);
  anatrB->AddFriend(trM);
  
  //===== Load cut files -----
  TFile *brcuts = new TFile("/home/koiwai/analysis/cutfiles/BRpid.root","");

  TCutG *cbr56sc = (TCutG*)brcuts->Get("br56sc");
  TCutG *cbr56ca = (TCutG*)brcuts->Get("br56ca");
  TCutG *cbr55ca = (TCutG*)brcuts->Get("br55ca");
  TCutG *cbr54ca = (TCutG*)brcuts->Get("br54ca");
  TCutG *cbr53ca = (TCutG*)brcuts->Get("br53ca");
  TCutG *cbr52ca = (TCutG*)brcuts->Get("br52ca");

  //brcuts->Close();

  TFile *sacuts = new TFile("/home/koiwai/analysis/cutfiles/SApid.root","");

  TCutG *csa56ca = (TCutG*)sacuts->Get("sa56ca");
  TCutG *csa55ca = (TCutG*)sacuts->Get("sa55ca");
  TCutG *csa54ca = (TCutG*)sacuts->Get("sa54ca");
  TCutG *csa53ca = (TCutG*)sacuts->Get("sa53ca");
  TCutG *csa52ca = (TCutG*)sacuts->Get("sa52ca");

  TCutG *csa55k  = (TCutG*)sacuts->Get("sa55k");

  //sacuts->Close();
  
  //=====ROOT file setting==========================================================
  TString ofname;
  if(FileNumber>1&&FileNumber<36){
  //ofname = Form("/home/koiwai/analysis/rootfiles/dali_calibrun/dali_calibrun%04d.root",FileNumber);
    if(FileNumber==10||FileNumber==19||FileNumber==20||FileNumber==32)
      ofname = Form("/home/koiwai/analysis/rootfiles/dali_calibrun/137cs%04d.root",FileNumber);
    else if(FileNumber==12||FileNumber==18||FileNumber==31)
      ofname = Form("/home/koiwai/analysis/rootfiles/dali_calibrun/88y%04d.root",FileNumber);
    else if(FileNumber==13||FileNumber==33)
      ofname = Form("/home/koiwai/analysis/rootfiles/dali_calibrun/133ba%04d.root",FileNumber);
    else if(FileNumber==14||FileNumber==27)
      ofname = Form("/home/koiwai/analysis/rootfiles/dali_calibrun/BG%04d.root",FileNumber);
    else
      ofname = Form("/home/koiwai/analysis/rootfiles/dali_calibrun/60co%04d.root",FileNumber);
      }
  else if(FileNumber>=36&&FileNumber<231)
    ofname = Form("/home/koiwai/analysis/rootfiles/dali/cal_dali%04d.root",FileNumber);

  TFile *outfile = new TFile(ofname,"RECREATE");
  TTree *tr = new TTree("tr","tr");
  tr->SetAutoSave(1e5); //Autosave every 10,000 events filled.

  //=====Define variables===================================================
  Long64_t EventNumber = 0;
  Int_t RunNumber = -1;

  vector <Double_t> *DALI_Energy = new vector <Double_t>();
  vector <Double_t> *DALI_CosTheta = new vector<Double_t>();
  vector <Double_t> *DALI_Time = new vector<Double_t>();
  vector <Int_t> *DALI_ID = new vector <Int_t>();
  vector <Double_t> *DALI_Layer = new vector<Double_t>();
  vector <Double_t> *DALI_X = new vector<Double_t>();
  vector <Double_t> *DALI_Y = new vector<Double_t>();
  vector <Double_t> *DALI_Z = new vector<Double_t>();
  //vector <TVector3> *DALI_Pos = new vector<TVector3>();

  Int_t DALI_Multi;
  Int_t source;

  Int_t br56sc, br56ca, br55ca, br54ca, br53ca, br52ca;
  Int_t sa56ca, sa55ca, sa54ca, sa53ca, sa52ca, sa55k;
  
  tr->Branch("EventNumber",&EventNumber);
  tr->Branch("RunNumber",&RunNumber);
  tr->Branch("DALI_Energy",&DALI_Energy);
  tr->Branch("DALI_CosTheta",&DALI_CosTheta);
  tr->Branch("DALI_Time",&DALI_Time);
  tr->Branch("DALI_ID",&DALI_ID);
  tr->Branch("DALI_Layer",&DALI_Layer);
  tr->Branch("DALI_X",&DALI_X);
  tr->Branch("DALI_Y",&DALI_Y);
  tr->Branch("DALI_Z",&DALI_Z);
  //tr->Branch("DALI_Pos",&DALI_Pos);
  tr->Branch("DALI_Multi",&DALI_Multi);
  tr->Branch("source",&source);

  tr->Branch("br56sc",&br56sc);
  tr->Branch("br56ca",&br56ca);
  tr->Branch("br55ca",&br55ca);
  tr->Branch("br54ca",&br54ca);
  tr->Branch("br53ca",&br53ca);
  tr->Branch("br52ca",&br52ca);

  tr->Branch("sa56ca",&sa56ca);
  tr->Branch("sa55ca",&sa55ca);
  tr->Branch("sa54ca",&sa54ca);
  tr->Branch("sa53ca",&sa53ca);
  tr->Branch("sa52ca",&sa52ca);
  tr->Branch("sa55k",&sa55k);
  
  
  while(EventStore->GetNextEvent()&&EventNumber<MaxEventNumber){
    //while(EventStore->GetNextEvent()&&EventNumber<5){
    EventNumber++;
    if(EventNumber%100 == 0){
      std::clog << EventNumber/1000 << "k events treated..." << "\r";
    }
    
    

    CalibDALI->ClearData();
    
    TArtRawEventObject *fEvent = (TArtRawEventObject *)sman->FindDataContainer("RawEvent");

    for(Int_t i=0;i<fEvent->GetNumSeg();i++){
      TArtRawSegmentObject *seg = fEvent->GetSegment(i);
      Int_t device = seg->GetDevice();
      Int_t detector = seg->GetDetector();

      if(device== BIGRIPS){
	switch(detector) {
	  //case   PPACT: if(kBEAM) CalibPPAC->LoadData(seg); break;
	  //case   PPACQ: if(kBEAM) CalibPPAC->LoadData(seg); break;
	default:  break;
	}
      } else {
	switch(detector) {
	case    DALIA: CalibDALI->LoadData(seg); break;
	case    DALIT: CalibDALI->LoadData(seg); break;
	default:  break;
	}
      } 
    }

    RunNumber = FileNumber;

    anatrB->GetEntry(EventNumber);
    
    DALI_ID->clear();
    DALI_Time->clear();
    DALI_Energy->clear();
    DALI_CosTheta->clear();
    DALI_Layer->clear();
    DALI_X->clear();
    DALI_Y->clear();
    DALI_Z->clear();
    //DALI_Pos->clear();
    
    DALI_Multi = 0;
    source = 0;

    br56sc = 0;
    br56ca = 0;
    br55ca = 0;
    br54ca = 0;
    br53ca = 0;
    br52ca = 0;

    sa56ca = 0;
    sa55ca = 0;
    sa54ca = 0;
    sa53ca = 0;
    sa52ca = 0;
    sa55k  = 0;
    
    CalibDALI->ReconstructData();
    
    //===== cut events not interested =====
    bool BRcut_bool = false;
    bool SAcut_bool = false;
    if(cbr56ca->IsInside(aoqBR,zetBR)) BRcut_bool = true;
    else if(cbr56sc->IsInside(aoqBR,zetBR)) BRcut_bool = true;
    else if(cbr54ca->IsInside(aoqBR,zetBR)) BRcut_bool = true;
    if(csa55ca->IsInside(aoqSA,zetSA)) SAcut_bool = true;
    else if(csa55k->IsInside(aoqSA,zetSA)) SAcut_bool = true;
    else if(csa53ca->IsInside(aoqSA,zetSA)) SAcut_bool = true;

    if((SAcut_bool==false)||(BRcut_bool==false)){
      tr->Fill();
      continue;
      }
    
    TClonesArray *DALINaIHits = (TClonesArray *)sman->FindDataContainer("DALINaI");
    if(DALINaIHits) {
      Int_t DALI_Mult = 0;
      Int_t NumberOfNaIHit = DALINaIHits->GetEntries();
      Double_t ADC;
      Double_t Time;
      for(Int_t i=0; i<NumberOfNaIHit; i++) {
	TArtDALINaI *DALINaI = (TArtDALINaI *)DALINaIHits->At(i);
	ADC = DALINaI->GetRawADC();
	Time = DALINaI->GetTimeOffseted();
	if(ADC > 0 &&  ADC < 4095 /*&& Time > 0*/) {
	  DALI_ID->push_back(DALINaI->GetID());
	  DALI_Energy->push_back(DALINaI->GetEnergy());
	  DALI_CosTheta->push_back(DALINaI->GetCosTheta());
	  DALI_Time->push_back(Time);
	  DALI_Layer->push_back(DALINaI->GetLayer());
	  DALI_X->push_back(DALINaI->GetXPos());
	  DALI_Y->push_back(DALINaI->GetYPos());
	  DALI_Z->push_back(DALINaI->GetZPos());
	  //DALI_Pos->push_back(TVector3(DALINaI->GetXPos(),DALINaI->GetYPos(),DALINaI->GetZPos()));
	  DALI_Mult++;
	  
	  
	}
      }
      DALI_Multi = DALI_Mult;
      // if (DALI_Mult > 1)
      //   SortDaliHit(0,DALI_Mult-1, DALI_ID, DALI_Energy, DALI_EnergyDopplerCorrected, DALI_Time, DALI_CosTheta);
    }

    //===== PID cut =====
    if(cbr56ca->IsInside(aoqBR,zetBR)) br56ca = 1;
    if(cbr55ca->IsInside(aoqBR,zetBR)) br55ca = 1;
    if(cbr54ca->IsInside(aoqBR,zetBR)) br54ca = 1;
    if(cbr53ca->IsInside(aoqBR,zetBR)) br53ca = 1;
    if(cbr52ca->IsInside(aoqBR,zetBR)) br52ca = 1;
    if(csa56ca->IsInside(aoqSA,zetSA)) sa56ca = 1;
    if(csa55ca->IsInside(aoqSA,zetSA)) sa55ca = 1;
    if(csa54ca->IsInside(aoqSA,zetSA)) sa54ca = 1;
    if(csa53ca->IsInside(aoqSA,zetSA)) sa53ca = 1;
    if(csa52ca->IsInside(aoqSA,zetSA)) sa52ca = 1;
    if(csa55k->IsInside(aoqSA,zetSA)) sa55k = 1;

    
    tr ->Fill();
    
  }//while loop
  std::clog << std::endl;

  tr->BuildIndex("RunNumber","EventNumber");
  outfile->cd();
  WriteAllEnvFiles;
  outfile->Write();
  outfile->Close("R");

  EventStore->ClearData();
  delete CalibDALI;
  
  delete DALI_ID;
  delete DALI_Time;
  delete DALI_Energy;
  delete DALI_CosTheta;
  delete TArtStoreManager::Instance();
  
  return 0;
}//main()

void SortDaliHit(Int_t left, Int_t right,vector <Int_t> *DALI_ID,vector <Double_t> *DALI_Energy, vector <Double_t> *DALI_EnergyDopplerCorrected, vector <Double_t> *DALI_Time, vector <Double_t> *DALI_CosTheta)
{
  Int_t TempID;
  Double_t TempEnergy;
  Double_t TempEnergyDopplerCorrected;
  Double_t TempTime;
  Double_t TempCosTheta;

  int i = left, j = right;
  double pivot = DALI_EnergyDopplerCorrected->at((left + right) / 2);

  /* partition */
  while (i <= j) {
    while (DALI_EnergyDopplerCorrected->at(i) > pivot)
      i++;
    while (DALI_EnergyDopplerCorrected->at(j) < pivot)
      j--;
    if (i <= j) {
      TempID = DALI_ID->at(j);
      TempEnergy = DALI_Energy->at(j);
      TempEnergyDopplerCorrected = DALI_EnergyDopplerCorrected->at(j);
      TempTime = DALI_Time->at(j);
      TempCosTheta = DALI_CosTheta->at(j);

      DALI_ID->at(j) = DALI_ID->at(i);
      DALI_Energy->at(j) = DALI_Energy->at(i);
      DALI_EnergyDopplerCorrected->at(j) = DALI_EnergyDopplerCorrected->at(i);
      DALI_Time->at(j) = DALI_Time->at(i);
      DALI_CosTheta->at(j) = DALI_CosTheta->at(i);

      DALI_ID->at(i) = TempID;
      DALI_Energy->at(i) = TempEnergy;
      DALI_EnergyDopplerCorrected->at(i) = TempEnergyDopplerCorrected;
      DALI_Time->at(i) = TempTime;
      DALI_CosTheta->at(i) = TempCosTheta;

      i++;
      j--;
    }
  };
  /* recursion */
  if (left < j)
    SortDaliHit(left, j, DALI_ID, DALI_Energy, DALI_EnergyDopplerCorrected, DALI_Time, DALI_CosTheta);
  if (i < right)
    SortDaliHit(i, right, DALI_ID, DALI_Energy, DALI_EnergyDopplerCorrected, DALI_Time, DALI_CosTheta);
}

Double_t DopplerCorrection(Double_t GammaDopplerEnergy, Double_t Beta, Double_t CosTheta) {
  Double_t Gamma = 1 / TMath::Sqrt(1 - Beta * Beta);
  Double_t DopplerCorrected = GammaDopplerEnergy * Gamma * (1 - Beta * CosTheta);
  return DopplerCorrected;
}

inline bool exists_test (const std::string& name) {
  return ( access( name.c_str(), F_OK ) != -1 );
}
inline bool exists_test (const TString& name) {
  return ( access( name.Data(), F_OK ) != -1 );
}
