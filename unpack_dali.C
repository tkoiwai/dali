//!===============================================================================
//!  2020 June 12 TK wrote.
//!  unpacking+calibration for DALI2+ date of SEASTAR3-2017
//!  Doppler cor will be done in next step
//!
//!
//!===============================================================================

#define unpack_dali_cxx
#include "/home/koiwai/analysis/include/unpackdalidef.h"

//*=====main Function==========================================================
int main(int argc, char *argv[]) {
  initiate_timer_tk();

  gInterpreter->GenerateDictionary("vector<TVector3>", "TVector3.h");

  Long64_t MaxEventNumber = LLONG_MAX;  //signed 8bit int. LLOMG_MAX ~ 2^63-1(const.)

  bool TestMode = false;
  bool ENum_flag = false;
  int FileNumber = -1;

  int opt;

  while((opt = getopt(argc, argv, "tr:e:")) != -1) {
    switch(opt) {
      case 'r':
        FileNumber = atoi(optarg);
        break;
      case 'e':
        ENum_flag = true;
        MaxEventNumber = atoi(optarg);
        break;
      case 't':
        TestMode = true;
        break;

      default:
        break;
    }
  }

  printf("=======================================\n");

  TString RidfFileName = Form("/home/koiwai/analysis/ridf/sdaq02/run%04d.ridf.gz", FileNumber);

  if(!exists_test(RidfFileName)) {
    cerr << " ERROR - '" << RidfFileName << "' does not exist '" << endl;
    exit(EXIT_FAILURE);
  }

  printf("\n%s %d %s \n\n", "=== Execute unpack_dali for RUN", FileNumber, "===");
  if(ENum_flag)
    printf("=== You will process %lld events\n \n", MaxEventNumber);

  //*=====Load setting parameters=========================================
  TEnv *env_par = new TEnv("/home/koiwai/analysis/conversion_settings.prm");
  TEnv *env_geo = new TEnv(env_par->GetValue("geometrydata", ""));  //unit of length:mm
  Double_t Dist_F5F7 = env_geo->GetValue("Dist_F5F7", 0.0);
  Double_t Dist_F7F13 = env_geo->GetValue("Dist_F7F13", 0.0);
  Double_t BDC1_Width = env_geo->GetValue("BDC1_Width", 0.0);                    //mm
  Double_t FDC1_Width = env_geo->GetValue("FDC1_Width", 0.0);                    //mm
  Double_t Dist_BDC1BDC2 = env_geo->GetValue("Dist_BDC1BDC2", 0.0);              //mm
  Double_t Dist_BDC1Target = env_geo->GetValue("Dist_BDC1Target", 0.0);          //mm
  Double_t Dist_BDC1FDC1 = env_geo->GetValue("Dist_BDC1FDC1", 0.0);              //Distance between the middle of BDC1 and the middle of FDC1 mm
  Double_t Dist_BDCFDC1 = env_geo->GetValue("Dist_BDCFDC1", 0.0);                //Distance between the middle of BDC and FDC1 mm
  Double_t Dist_SBTTarget = env_geo->GetValue("Dist_SBTTarget", 0.0);            //mm
  Double_t Dist_MINOSfrontFDC1 = env_geo->GetValue("Dist_MINOSfrontFDC1", 0.0);  //mm
  Double_t Dist_MINOSfrontBDC = env_geo->GetValue("Dist_MINOSfrontBDC", 0.0);    //mm, negative value

  Double_t MINOSoffsetZ = env_geo->GetValue("MINOSoffsetZ", 0.0);  //offset of vertexZ. [mm]
  Double_t DALIoffset = env_geo->GetValue("DALIoffset", 0.0);

  //*=====Load ANAROOT parameters===========================================
  TArtStoreManager *sman = TArtStoreManager::Instance();
  TArtDALIParameters *DALIPara = TArtDALIParameters::Instance();

  if(FileNumber >= 36 && FileNumber <= 61)
    DALIPara->LoadParameter((char *)"/home/koiwai/analysis/db/DALI/DALI_20170505.xml");
  else if(FileNumber >= 62 && FileNumber <= 124)
    DALIPara->LoadParameter((char *)"/home/koiwai/analysis/db/DALI/DALI_20170508.xml");
  if(FileNumber >= 125 && FileNumber <= 196)
    DALIPara->LoadParameter((char *)"/home/koiwai/analysis/db/DALI/DALI_20170511.xml");
  if(FileNumber >= 197 && FileNumber <= 230)
    DALIPara->LoadParameter((char *)"/home/koiwai/analysis/db/DALI/DALI_20170514.xml");

  //*===Load detectors.
  TArtCalibDALI *CalibDALI = new TArtCalibDALI();

  //=====Conversion ridf -> raw data======================================
  TArtEventStore *EventStore = new TArtEventStore;

  if(!EventStore->Open(RidfFileName)) {
    std::cerr << "cannot open " << RidfFileName << std::endl
              << "aborting..." << std::endl;
    return 1;
  } else {
    cout << "input file: " << RidfFileName << endl;
  }

  //=====ROOT file setting ==========================================================
  TString ofname;
  if(!TestMode)
    ofname = Form("/home/koiwai/analysis/rootfiles/ana/dali/unpack_dali%04d.root", FileNumber);
  else
    ofname = Form("/home/koiwai/analysis/dali/test_unpack_dali%04d.root", FileNumber);

  TFile *outfile = new TFile(ofname, "RECREATE");
  TTree *tr = new TTree("tr", "tr");
  //  tr->SetAutoSave(1e5); //Autosave every 10,000 events filled.

  printf("%-20s %s \n\n", "Output dali file:", ofname.Data());

  Set_Branch_dali(tr);

  //===== LOOP start ================================================================
  int iEntry = 0;
  int AllEntry;
  //printf("\n Number of events to treat: %d\n\n", AllEntry);

  prepare_timer_tk();

  while(EventStore->GetNextEvent() && EventNumber < MaxEventNumber) {
    iEntry = EventNumber;
    // anatrB->GetEntry(EventNumber);
    EventNumber++;

    start_timer_tk(iEntry, 1000);

    CalibDALI->ClearData();

    TArtRawEventObject *fEvent = (TArtRawEventObject *)sman->FindDataContainer("RawEvent");

    for(Int_t i = 0; i < fEvent->GetNumSeg(); i++) {
      TArtRawSegmentObject *seg = fEvent->GetSegment(i);
      Int_t device = seg->GetDevice();
      Int_t detector = seg->GetDetector();

      if(device == BIGRIPS) {
        switch(detector) {
            //case   PPACT: if(kBEAM) CalibPPAC->LoadData(seg); break;
            //case   PPACQ: if(kBEAM) CalibPPAC->LoadData(seg); break;
          default:
            break;
        }
      } else {
        switch(detector) {
          case DALIA:
            CalibDALI->LoadData(seg);
            break;
          case DALIT:
            CalibDALI->LoadData(seg);
            break;
          default:
            break;
        }
      }
    }

    TClonesArray *info_array = (TClonesArray *)sman->FindDataContainer("EventInfo");
    //cout << "RunNumber: " << ((TArtEventInfo *)info_array->At(0))->GetRunNumber() << endl;

    RunNumber = FileNumber;

    DALI_ID->clear();
    DALI_Time->clear();
    DALI_Energy->clear();
    DALI_CosTheta->clear();
    DALI_X->clear();
    DALI_Y->clear();
    DALI_Z->clear();

    DALI_ID_orig->clear();
    DALI_Time_orig->clear();
    DALI_Energy_orig->clear();
    DALI_CosTheta_orig->clear();
    DALI_X_orig->clear();
    DALI_Y_orig->clear();
    DALI_Z_orig->clear();

    DALI_Multi = 0;

    CalibDALI->ReconstructData();

    TClonesArray *DALINaIHits = (TClonesArray *)sman->FindDataContainer("DALINaI");
    if(DALINaIHits) {
      Int_t DALI_Mult = 0;
      Int_t NumberOfNaIHit = DALINaIHits->GetEntries();
      Double_t ADC;
      Double_t Time;

      //=====DEBUG==========================

      //TArtDALINaI *DALINaI = (TArtDALINaI *)DALINaIHits->At(0);
      //DALI_ID->push_back(DALINaI->GetID());
      //cout << DALINaI->GetID();

      //=====DEBUG end========================

      for(Int_t i = 0; i < NumberOfNaIHit; i++) {
        TArtDALINaI *DALINaI = (TArtDALINaI *)DALINaIHits->At(i);
        ADC = DALINaI->GetRawADC();
        Time = DALINaI->GetTimeOffseted();

        if(ADC > 0 && ADC < 4095 && DALINaI->GetEnergy() > 0 /*&& Time > 0*/) {
          DALI_ID->push_back(DALINaI->GetID());

          DALI_Energy->push_back(DALINaI->GetEnergy());
          DALI_CosTheta->push_back(DALINaI->GetCosTheta());
          DALI_Time->push_back(Time);
          DALI_X->push_back(DALINaI->GetXPos());
          DALI_Y->push_back(DALINaI->GetYPos());
          DALI_Z->push_back(DALINaI->GetZPos());

          DALI_ID_orig->push_back(DALINaI->GetID());
          DALI_Energy_orig->push_back(DALINaI->GetEnergy());
          DALI_CosTheta_orig->push_back(DALINaI->GetCosTheta());
          DALI_Time_orig->push_back(Time);
          DALI_X_orig->push_back(DALINaI->GetXPos());
          DALI_Y_orig->push_back(DALINaI->GetYPos());
          DALI_Z_orig->push_back(DALINaI->GetZPos());

          DALI_Mult++;
        }
      }
      DALI_Multi = DALI_ID->size();
    }

    if(DALI_Multi > 1) {
      SortDaliHit(0, DALI_Multi - 1,
                  DALI_ID,
                  DALI_Energy,
                  DALI_Time,
                  DALI_X,
                  DALI_Y,
                  DALI_Z,
                  DALI_CosTheta);
    } else if(DALI_Multi == 0) {
      tr->Fill();
      continue;
    }

    //===== ADDBACK END =======================

    //===== TIMING GATE =======================
    //===== TIMING GATE END ===================

    tr->Fill();

  }  //while loop
  AllEntry = iEntry;
  std::clog << std::endl;

  tr->BuildIndex("RunNumber", "EventNumber");
  outfile->cd();
  outfile->Write();
  outfile->Close("R");

  EventStore->ClearData();
  delete CalibDALI;

  delete DALI_ID;
  delete DALI_Time;
  delete DALI_Energy;
  delete DALI_CosTheta;
  delete DALI_X;
  delete DALI_Y;
  delete DALI_Z;

  delete DALI_ID_orig;
  delete DALI_Time_orig;
  delete DALI_Energy_orig;
  delete DALI_CosTheta_orig;
  delete DALI_X_orig;
  delete DALI_Y_orig;
  delete DALI_Z_orig;

  delete TArtStoreManager::Instance();

  stop_timer_tk(FileNumber, AllEntry);

  return 0;
}  //main()

void SortDaliHit(Int_t left, Int_t right, vector<Int_t> *DALI_ID, vector<Double_t> *DALI_Energy, vector<Double_t> *DALI_Time, vector<Double_t> *DALI_X, vector<Double_t> *DALI_Y, vector<Double_t> *DALI_Z, vector<Double_t> *DALI_CosTheta) {
  Int_t TempID;
  Double_t TempEnergy;
  Double_t TempTime;
  Double_t TempX;
  Double_t TempY;
  Double_t TempZ;
  Double_t TempCosTheta;

  int i = left, j = right;
  double pivot = DALI_Energy->at((left + right) / 2);

  //--partition---/
  while(i <= j) {
    while(DALI_Energy->at(i) > pivot)
      i++;
    while(DALI_Energy->at(j) < pivot)
      j--;
    if(i <= j) {
      TempID = DALI_ID->at(j);
      TempEnergy = DALI_Energy->at(j);
      TempTime = DALI_Time->at(j);
      TempX = DALI_X->at(j);
      TempY = DALI_Y->at(j);
      TempZ = DALI_Z->at(j);
      TempCosTheta = DALI_CosTheta->at(j);

      DALI_ID->at(j) = DALI_ID->at(i);
      DALI_Energy->at(j) = DALI_Energy->at(i);
      DALI_Time->at(j) = DALI_Time->at(i);
      DALI_X->at(j) = DALI_X->at(i);
      DALI_Y->at(j) = DALI_Y->at(i);
      DALI_Z->at(j) = DALI_Z->at(i);
      DALI_CosTheta->at(j) = DALI_CosTheta->at(i);

      DALI_ID->at(i) = TempID;
      DALI_Energy->at(i) = TempEnergy;
      DALI_Time->at(i) = TempTime;
      DALI_X->at(i) = TempX;
      DALI_Y->at(i) = TempY;
      DALI_Z->at(i) = TempZ;
      DALI_CosTheta->at(i) = TempCosTheta;

      i++;
      j--;
    }
  };

  // recursion //
  if(left < j)
    SortDaliHit(left, j, DALI_ID, DALI_Energy, DALI_Time, DALI_X, DALI_Y, DALI_Z, DALI_CosTheta);
  if(i < right)
    SortDaliHit(i, right, DALI_ID, DALI_Energy, DALI_Time, DALI_X, DALI_Y, DALI_Z, DALI_CosTheta);
}

Double_t DopplerCorrection(Double_t GammaDopplerEnergy, Double_t Beta, Double_t CosTheta) {
  Double_t Gamma = 1 / TMath::Sqrt(1 - Beta * Beta);
  Double_t DopplerCorrected = GammaDopplerEnergy * Gamma * (1 - Beta * CosTheta);
  return DopplerCorrected;
}

inline bool exists_test(const std::string &name) {
  return (access(name.c_str(), F_OK) != -1);
}
inline bool exists_test(const TString &name) {
  return (access(name.Data(), F_OK) != -1);
}
