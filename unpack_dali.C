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

  if(argc < 2 || argc > 4) {
    printf("Usage: ./unpack_dali RUNNUMBER\nOR     ./unpack_dali RUNNUMBER MAXEVENTS\nOR     ./unpack_dali RUNNUMBER MAXEVENTS TEST\n");
    exit(EXIT_FAILURE);
  }
  printf("=======================================\n");
  if(argc > 2) {
    MaxEventNumber = TString(argv[2]).Atoi();
    printf(" You will process %lld events\n", MaxEventNumber);
  }

  Int_t FileNumber = TString(argv[1]).Atoi();
  TString RidfFileName = Form("/home/koiwai/analysis/ridf/sdaq02/run%04d.ridf.gz", FileNumber);

  if(!exists_test(RidfFileName)) {
    cerr << " ERROR - '" << RidfFileName << "' does not exist '" << endl;
    exit(EXIT_FAILURE);
  }

  printf("\n%s %d %s \n\n", "=== Execute cal_dali for RUN", FileNumber, "===");

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

  //===== Load Beam file for PID cut =====
  TString infnameB = Form("/home/koiwai/analysis/rootfiles/ana/beam/ana_beam%04d.root", FileNumber);
  TFile *infileB = TFile::Open(infnameB);
  TTree *anatrB;
  infileB->GetObject("anatrB", anatrB);

  printf("%-20s %s \n", "Input beam file:", infnameB.Data());

  anaD_Get_Branch_beam(anatrB);

  //===== Load smri file for PID cut =====
  TString infnameS = Form("/home/koiwai/analysis/rootfiles/ana/smri/ana_smri%04d.root", FileNumber);
  TFile *infileS = TFile::Open(infnameS);
  TTree *anatrS;
  infileS->GetObject("anatrS", anatrS);

  printf("%-20s %s \n", "Input smri file:", infnameS.Data());

  anaD_Get_Branch_smri(anatrS);

  //===== Load MINOS file =====
  TString infnameM = Form("/home/koiwai/analysis/rootfiles/minos/vertex/vertex%04d.root", FileNumber);
  TFile *infileM = TFile::Open(infnameM);
  TTree *trM;
  infileM->GetObject("tr", trM);

  printf("%-20s %s \n\n", "Input minos file:", infnameM.Data());

  anaD_Get_Branch_minos(trM);

  //===== AddFriend =====

  anatrB->AddFriend(anatrS);
  anatrB->AddFriend(trM);

  //===== Load cuts =======================================================================
  TFile *fcutSA_Ca = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_Ca.root");

  TCutG *csa57ca = (TCutG *)fcutSA_Ca->Get("sa57ca");
  TCutG *csa55ca = (TCutG *)fcutSA_Ca->Get("sa55ca");
  TCutG *csa54ca = (TCutG *)fcutSA_Ca->Get("sa54ca");
  TCutG *csa53ca = (TCutG *)fcutSA_Ca->Get("sa53ca");

  //=====ROOT file setting ==========================================================
  TString ofname;
  if(argc < 4)
    ofname = Form("/home/koiwai/analysis/rootfiles/ana/dali/cal_dali%04d.root", FileNumber);
  else if(argc == 4)
    ofname = Form("/home/koiwai/analysis/dali/testcal_dali%04d.root", FileNumber);

  TFile *outfile = new TFile(ofname, "RECREATE");
  TTree *tr = new TTree("tr", "tr");
  //  tr->SetAutoSave(1e5); //Autosave every 10,000 events filled.

  printf("%-20s %s \n\n", "Output dali file:", ofname.Data());

  Set_Branch_dali(tr);

  //===== LOOP start ================================================================
  int iEntry = 0;
  int AllEntry;
  if(argc > 2)
    AllEntry = MaxEventNumber;
  else
    AllEntry = anatrB->GetEntries();

  printf("\n Number of events to treat: %d\n\n", AllEntry);

  prepare_timer_tk();

  while(EventStore->GetNextEvent() && EventNumber < MaxEventNumber) {
    iEntry = EventNumber;
    anatrB->GetEntry(EventNumber);
    EventNumber++;

    start_timer_tk(iEntry, AllEntry, 1000);

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

    vertexZ_cor = Sqrt(-1);
    beta_vertex = Sqrt(-1);
    gamma_vertex = Sqrt(-1);

    dali_edop->clear();

    vertex.SetXYZ(Sqrt(-1), Sqrt(-1), Sqrt(-1));
    dali_pos.clear();

    gamma_pos.clear();
    gamma_cos.clear();

    dali_edop_simple->clear();

    dali_e_ab->clear();
    dali_id_ab->clear();
    dali_cos_ab->clear();
    dali_t_ab->clear();
    dali_x_ab->clear();
    dali_y_ab->clear();
    dali_z_ab->clear();
    dali_multi_ab = 0;

    dali_edop_ab->clear();

    sa57ca = kFALSE;
    sa55ca = kFALSE;
    sa54ca = kFALSE;
    sa53ca = kFALSE;

    CalibDALI->ReconstructData();

    //===== cut events not interested =====
    bool BRcut_bool = false;

    if((br59sc == 1) || (br58sc == 1) || (br57sc == 1) || (br56sc == 1) || (br55sc == 1) || (br54sc == 1) ||
       (br56ca == 1) || (br55ca == 1) || (br54ca == 1) || (br53ca == 1) || (br52ca == 1) ||
       (br54k == 1) || (br53k == 1) || (br52k == 1) || (br51k == 1) || (br50k == 1) || (br49k == 1))
      BRcut_bool = true;

    if(BRcut_bool == false) {
      tr->Fill();
      continue;
    }

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

    vertexZ_cor = vertexZ + MINOSoffsetZ;

    beta_vertex = betaF7F13 - (betaF7F13 - betaTH) * vertexZ_cor / 150.0;
    gamma_vertex = 1 / Sqrt(1 - beta_vertex * beta_vertex);

    vertex.SetXYZ(vertexX, vertexY, vertexZ_cor - DALIoffset);
    //To match the centre of MINOS cell and DALI Z = 0.(DALIOffset)

    for(Int_t i = 0; i < DALI_Multi; i++) {
      dali_pos.push_back(TVector3(10 * DALI_X->at(i), 10 * DALI_Y->at(i), 10 * DALI_Z->at(i)));
      gamma_pos.push_back(dali_pos.at(i) - vertex);
      gamma_cos.push_back((gamma_pos.at(i)).CosTheta());
    }

    //===== DOPPLER CORRECTION =====

    if(DALI_Multi >= 1) {
      for(Int_t i = 0; i < DALI_Multi; i++) {
        Double_t dali_edop_tmp = Sqrt(-1);
        dali_edop_tmp = DALI_Energy->at(i) * gamma_vertex * (1 - beta_vertex * gamma_cos.at(i));
        dali_edop->push_back(dali_edop_tmp);
      }
    }

    //===== DOPPLER CORRECTION END =====

    //===== SINPLE DOPPLER CORRECTIOMN (WITHOUT MINOS )=====

    const Double_t beta_mid = 0.57;
    const Double_t gamma_mid = 1 / Sqrt(1 - beta_mid * beta_mid);
    if(DALI_Multi >= 1) {
      for(Int_t i = 0; i < DALI_Multi; i++) {
        dali_edop_simple->push_back(DALI_Energy->at(i) * gamma_mid * (1 - beta_mid * DALI_CosTheta->at(i)));
      }
    }

    //===== SIMPLE DOPPLER CORRECTION END =====

    //===== ADDBACK ===========================
    if(DALI_Multi > 1) {
      dali_e_ab->push_back(DALI_Energy->at(0));
      dali_id_ab->push_back(DALI_ID->at(0));
      dali_cos_ab->push_back(DALI_CosTheta->at(0));
      dali_t_ab->push_back(DALI_Time->at(0));
      dali_x_ab->push_back(DALI_X->at(0));
      dali_y_ab->push_back(DALI_Y->at(0));
      dali_z_ab->push_back(DALI_Z->at(0));

      //===== ADD BACK =====

      Bool_t AddBack_flag = false;

      for(Int_t i = 1; i < DALI_Multi; i++) {
        AddBack_flag = false;
        dali_multi_ab = dali_id_ab->size();
        for(Int_t k = 0; k < dali_multi_ab; k++) {
          for(Int_t j = 0; j < 7; j++) {
            if(DALI_ID->at(i) == AddBackTable[dali_id_ab->at(k)][j] && DALI_Energy->at(i) > 300) {
              dali_e_ab->at(k) = dali_e_ab->at(k) + DALI_Energy->at(i);
              AddBack_flag = true;
              break;
            }
          }
          if(AddBack_flag == true) break;
        }
        if(AddBack_flag == false) {
          dali_e_ab->push_back(DALI_Energy->at(i));
          dali_id_ab->push_back(DALI_ID->at(i));
          dali_cos_ab->push_back(DALI_CosTheta->at(i));
          dali_t_ab->push_back(DALI_Time->at(i));
          dali_x_ab->push_back(DALI_X->at(i));
          dali_y_ab->push_back(DALI_Y->at(i));
          dali_z_ab->push_back(DALI_Z->at(i));
        }
      }  //for DALI_Multi

      dali_multi_ab = dali_id_ab->size();

    }  //DALI_Multi>1
    else if(DALI_Multi == 1) {
      dali_e_ab->push_back(DALI_Energy->at(0));
      dali_id_ab->push_back(DALI_ID->at(0));
      dali_cos_ab->push_back(DALI_CosTheta->at(0));
      dali_t_ab->push_back(DALI_Time->at(0));
      dali_x_ab->push_back(DALI_X->at(0));
      dali_y_ab->push_back(DALI_Y->at(0));
      dali_z_ab->push_back(DALI_Z->at(0));

      dali_multi_ab = 1;
    } else if(DALI_Multi == 0) {
      tr->Fill();
      continue;
    }

    if(dali_multi_ab >= 1) {
      for(Int_t i = 0; i < dali_multi_ab; i++) {
        Double_t dali_edop_ab_tmp = Sqrt(-1);
        dali_edop_ab_tmp = dali_e_ab->at(i) * gamma_vertex * (1 - beta_vertex * gamma_cos.at(i));
        dali_edop_ab->push_back(dali_edop_ab_tmp);
      }
    }

    //===== ADDBACK END =======================

    //===== TIMING GATE =======================
    //===== TIMING GATE END ===================

    //===== isotope gate =====
    if(csa57ca->IsInside(aoqSA, zetSA)) sa57ca = kTRUE;
    if(csa55ca->IsInside(aoqSA, zetSA)) sa55ca = kTRUE;
    if(csa54ca->IsInside(aoqSA, zetSA)) sa54ca = kTRUE;
    if(csa53ca->IsInside(aoqSA, zetSA)) sa53ca = kTRUE;

    tr->Fill();

  }  //while loop
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

  delete dali_edop;
  delete dali_edop_simple;

  delete dali_e_ab;
  delete dali_id_ab;
  delete dali_cos_ab;
  delete dali_t_ab;
  delete dali_x_ab;
  delete dali_y_ab;
  delete dali_z_ab;

  delete dali_edop_ab;

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

/*
void SortDaliHit(Int_t left, Int_t right,vector <Int_t> *DALI_ID,vector <Double_t> *DALI_Energy, vector <Double_t> *DALI_EnergyDopplerCorrected, vector <Double_t> *DALI_Time, vector <Double_t> *DALI_CosThetav)
{
  Int_t TempID;
  Double_t TempEnergy;
  Double_t TempEnergyDopplerCorrected;
  Double_t TempTime;
  Double_t TempCosTheta;

  int i = left, j = right;
  double pivot = DALI_EnergyDopplerCorrected->at((left + right) / 2);
*/
/* partition */
/*
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
*/
/* recursion */
/*
  if (left < j)
    SortDaliHit(left, j, DALI_ID, DALI_Energy, DALI_EnergyDopplerCorrected, DALI_Time, DALI_CosTheta);
  if (i < right)
    SortDaliHit(i, right, DALI_ID, DALI_Energy, DALI_EnergyDopplerCorrected, DALI_Time, DALI_CosTheta);
}
*/

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
