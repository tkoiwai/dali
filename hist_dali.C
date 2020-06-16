#include "../include/anadalidef.h"

//*=====main Function==========================================================
int main(int argc, char *argv[]) {
  initiate_timer_tk();

  gInterpreter->GenerateDictionary("vector<TVector3>", "TVector3.h");

  Int_t FileNumber = TString(argv[1]).Atoi();

  if(FileNumber == 0) {
    std::cerr << " You should provide either a runnumber" << endl;
  }
  if(argc < 2 || argc > 5) {
    printf("Usage: ./cal_dali RUNNUMBER\nOR     ./cal_dali RUNNUMBER MAXEVENTS\nOR     ./cal_dali RUNNUMBER MAXEVENTS TEST\n");
    exit(EXIT_FAILURE);
  }

  int MaxEventNumber = 0;

  if(argc > 2) {
    MaxEventNumber = TString(argv[2]).Atoi();
    printf(" You will process %d events\n", MaxEventNumber);
  }

  const double MINOSoffsetZ = 12.37;
  const double DALIoffset = 96.5;

  //===== Load input files =====
  TString infname = Form("/home/koiwai/rootfiles/ana/dali/unpack_dali%04d.root", FileNumber);
  TFile *infile = TFile::Open(infname);
  TTree *intr = (TTree *)infile->Get("tr");

  Get_Branch_calD(intr);

  int AddBackTable[226][7] = {
      {10, 19, -1, -1, -1, -1, -1},  //0
      {10, 11, -1, -1, -1, -1, -1},
      {11, 12, -1, -1, -1, -1, -1},
      {12, 13, -1, -1, -1, -1, -1},
      {13, 14, -1, -1, -1, -1, -1},
      {14, 15, -1, -1, -1, -1, -1},
      {15, 16, -1, -1, -1, -1, -1},
      {16, 17, -1, -1, -1, -1, -1},
      {17, 18, -1, -1, -1, -1, -1},
      {18, 19, -1, -1, -1, -1, -1},
      {0, 1, 20, -1, -1, -1, -1},  //10
      {1, 2, 22, -1, -1, -1, -1},
      {2, 3, 23, -1, -1, -1, -1},
      {3, 4, 24, -1, -1, -1, -1},
      {4, 5, 25, -1, -1, -1, -1},
      {5, 6, 26, -1, -1, -1, -1},
      {6, 7, 28, -1, -1, -1, -1},
      {7, 8, 29, -1, -1, -1, -1},
      {8, 9, 30, -1, -1, -1, -1},
      {0, 9, 31, -1, -1, -1, -1},
      {10, 32, -1, -1, -1, -1, -1},  //20
      {33, -1, -1, -1, -1, -1, -1},
      {11, 34, -1, -1, -1, -1, -1},
      {12, 35, 36, -1, -1, -1, -1},
      {13, 37, -1, -1, -1, -1, -1},
      {14, 38, -1, -1, -1, -1, -1},
      {15, 39, -1, -1, -1, -1, -1},
      {40, -1, -1, -1, -1, -1, -1},
      {16, 41, -1, -1, -1, -1, -1},
      {17, 42, 43, -1, -1, -1, -1},
      {18, 44, -1, -1, -1, -1, -1},  //30
      {19, 45, -1, -1, -1, -1, -1},
      {20, -1, -1, -1, -1, -1, -1},
      {21, -1, -1, -1, -1, -1, -1},
      {22, 48, -1, -1, -1, -1, -1},
      {23, 49, 50, -1, -1, -1, -1},
      {23, 51, 52, -1, -1, -1, -1},
      {24, 53, -1, -1, -1, -1, -1},
      {25, 54, -1, -1, -1, -1, -1},
      {26, -1, -1, -1, -1, -1, -1},
      {27, -1, -1, -1, -1, -1, -1},  //40
      {28, 58, -1, -1, -1, -1, -1},
      {29, 59, 60, -1, -1, -1, -1},
      {29, 61, 62, -1, -1, -1, -1},
      {30, 63, -1, -1, -1, -1, -1},
      {31, 64, -1, -1, -1, -1, -1},
      {66, -1, -1, -1, -1, -1, -1},
      {67, -1, -1, -1, -1, -1, -1},
      {34, 49, 68, -1, -1, -1, -1},
      {35, 48, 50, 69, -1, -1, -1},
      {35, 49, 51, 70, -1, -1, -1},  //50
      {36, 50, 52, 71, -1, -1, -1},
      {36, 51, 53, 72, -1, -1, -1},
      {37, 52, 73, -1, -1, -1, -1},
      {38, 74, -1, -1, -1, -1, -1},
      {75, -1, -1, -1, -1, -1, -1},
      {76, -1, -1, -1, -1, -1, -1},
      {77, -1, -1, -1, -1, -1, -1},
      {41, 59, 78, -1, -1, -1, -1},
      {42, 58, 60, 79, -1, -1, -1},
      {42, 59, 61, 80, -1, -1, -1},  //60
      {43, 60, 62, 81, -1, -1, -1},
      {43, 61, 63, 82, -1, -1, -1},
      {44, 62, 83, -1, -1, -1, -1},
      {45, 84, -1, -1, -1, -1, -1},
      {85, -1, -1, -1, -1, -1, -1},
      {46, 86, -1, -1, -1, -1, -1},
      {47, 87, -1, -1, -1, -1, -1},
      {48, 69, 88, -1, -1, -1, -1},
      {49, 68, 70, 89, -1, -1, -1},
      {50, 69, 71, 90, -1, -1, -1},  //70
      {51, 70, 72, 91, -1, -1, -1},
      {52, 71, 73, 92, -1, -1, -1},
      {53, 72, 93, -1, -1, -1, -1},
      {54, 94, -1, -1, -1, -1, -1},
      {55, 95, -1, -1, -1, -1, -1},
      {56, 96, -1, -1, -1, -1, -1},
      {57, 97, -1, -1, -1, -1, -1},
      {58, 79, 98, -1, -1, -1, -1},
      {59, 78, 80, 99, -1, -1, -1},
      {60, 79, 81, 100, -1, -1, -1},  //80
      {61, 80, 82, 101, -1, -1, -1},
      {62, 81, 83, 102, -1, -1, -1},
      {63, 82, 103, -1, -1, -1, -1},
      {64, 104, -1, -1, -1, -1, -1},
      {65, 105, -1, -1, -1, -1, -1},
      {66, 106, 107, -1, -1, -1, -1},
      {67, 108, 109, -1, -1, -1, -1},
      {68, 89, 110, -1, -1, -1, -1},
      {69, 88, 90, 110, 111, -1, -1},
      {70, 89, 91, 111, 112, -1, -1},  //90
      {71, 90, 92, 113, 114, -1, -1},
      {72, 91, 93, 114, 115, -1, -1},
      {73, 92, 115, -1, -1, -1, -1},
      {74, 116, 117, -1, -1, -1, -1},
      {75, 118, 119, -1, -1, -1, -1},
      {76, 120, 121, -1, -1, -1, -1},
      {77, 122, 123, -1, -1, -1, -1},
      {78, 99, 124, -1, -1, -1, -1},
      {79, 98, 100, 124, 125, -1, -1},
      {80, 99, 101, 125, 126, -1, -1},  //100
      {81, 100, 102, 127, 128, -1, -1},
      {82, 101, 103, 128, 129, -1, -1},
      {83, 102, 129, -1, -1, -1, -1},
      {84, 130, 131, -1, -1, -1, -1},
      {85, 132, 133, -1, -1, -1, -1},
      {86, 107, 133, 134, -1, -1, -1},
      {86, 106, 108, 135, -1, -1, -1},
      {87, 107, 109, 136, -1, -1, -1},
      {87, 108, 137, -1, -1, -1, -1},
      {88, 89, 111, 138, -1, -1, -1},  //110
      {89, 90, 110, 112, 139, -1, -1},
      {90, 111, 113, 140, -1, -1, -1},
      {91, 112, 114, 141, -1, -1, -1},
      {91, 92, 113, 115, 142, -1, -1},
      {92, 93, 114, 143, -1, -1, -1},
      {94, 117, 144, -1, -1, -1, -1},
      {94, 116, 118, 145, -1, -1, -1},
      {95, 117, 119, 146, -1, -1, -1},
      {95, 118, 120, 147, -1, -1, -1},
      {96, 119, 121, 148, -1, -1, -1},  //120
      {96, 120, 122, 149, -1, -1, -1},
      {97, 121, 123, 150, -1, -1, -1},
      {97, 122, 151, -1, -1, -1, -1},
      {98, 99, 125, 152, -1, -1, -1},
      {99, 100, 124, 126, 153, -1, -1},
      {100, 125, 127, 154, -1, -1, -1},
      {101, 126, 128, 155, -1, -1, -1},
      {101, 102, 127, 129, 156, -1, -1},
      {102, 103, 128, 157, -1, -1, -1},
      {104, 131, 158, -1, -1, -1, -1},  //130
      {104, 130, 132, 159, -1, -1, -1},
      {105, 131, 133, 160, -1, -1, -1},
      {105, 106, 132, 161, -1, -1, -1},
      {106, 135, 161, -1, -1, -1, -1},
      {107, 134, 136, -1, -1, -1, -1},
      {108, 135, 137, -1, -1, -1, -1},
      {109, 136, 138, -1, -1, -1, -1},
      {110, 137, 139, -1, -1, -1, -1},
      {111, 138, 140, -1, -1, -1, -1},
      {112, 139, 141, -1, -1, -1, -1},  //140
      {113, 140, 142, -1, -1, -1, -1},
      {114, 141, 143, -1, -1, -1, -1},
      {115, 142, 144, -1, -1, -1, -1},
      {116, 143, 145, -1, -1, -1, -1},
      {117, 144, 146, -1, -1, -1, -1},
      {118, 145, 147, -1, -1, -1, -1},
      {119, 146, 148, -1, -1, -1, -1},
      {120, 147, 149, -1, -1, -1, -1},
      {121, 148, 150, -1, -1, -1, -1},
      {122, 149, 151, -1, -1, -1, -1},  //150
      {123, 150, 152, -1, -1, -1, -1},
      {124, 151, 153, -1, -1, -1, -1},
      {125, 152, 154, -1, -1, -1, -1},
      {126, 153, 155, -1, -1, -1, -1},
      {127, 154, 156, -1, -1, -1, -1},
      {128, 155, 157, -1, -1, -1, -1},
      {129, 156, 158, -1, -1, -1, -1},
      {130, 157, 159, -1, -1, -1, -1},
      {131, 158, 160, -1, -1, -1, -1},
      {132, 159, 161, -1, -1, -1, -1},  //160
      {133, 134, 160, -1, -1, -1, -1},
      {163, 171, 172, 191, -1, -1, -1},
      {162, 172, 173, 192, -1, -1, -1},
      {165, 176, 177, 195, -1, -1, -1},
      {164, 177, 178, 196, -1, -1, -1},
      {167, 181, 182, 199, -1, -1, -1},
      {166, 182, 183, 200, -1, -1, -1},
      {169, 186, 187, 203, -1, -1, -1},
      {168, 187, 188, 204, -1, -1, -1},
      {171, 189, 190, -1, -1, -1, -1},  //170
      {162, 170, 172, 190, 191, -1, -1},
      {162, 163, 171, 173, 191, 192, 207},
      {163, 172, 174, 192, 193, -1, -1},
      {173, 175, 193, -1, -1, -1, -1},
      {174, 176, 194, -1, -1, -1, -1},
      {164, 175, 177, 194, 195, -1, -1},
      {164, 165, 176, 178, 195, 196, 210},
      {165, 177, 179, 196, 197, -1, -1},
      {178, 180, 197, -1, -1, -1, -1},
      {179, 181, 198, -1, -1, -1, -1},  //180
      {166, 180, 182, 198, 199, -1, -1},
      {166, 167, 181, 183, 199, 200, 213},
      {167, 182, 184, 200, 201, -1, -1},
      {183, 185, 201, -1, -1, -1, -1},
      {184, 186, 202, -1, -1, -1, -1},
      {168, 185, 187, 202, 203, -1, -1},
      {168, 169, 186, 188, 203, 204, 216},
      {169, 187, 189, 204, 205, -1, -1},
      {170, 188, 205, -1, -1, -1, -1},
      {170, 171, 191, 205, 206, -1, -1},  //190
      {162, 171, 172, 190, 192, 206, 207},
      {163, 172, 173, 191, 193, 207, 208},
      {173, 174, 192, 194, 208, -1, -1},
      {175, 176, 193, 195, 209, -1, -1},
      {164, 176, 177, 194, 196, 209, 210},
      {165, 177, 178, 195, 197, 210, 211},
      {178, 179, 196, 198, 211, -1, -1},
      {180, 181, 197, 199, 212, -1, -1},
      {166, 181, 182, 198, 200, 212, 213},
      {167, 182, 183, 199, 201, 213, 214},  //200
      {183, 184, 200, 202, 214, -1, -1},
      {185, 186, 201, 203, 215, -1, -1},
      {168, 186, 187, 202, 204, 215, 216},
      {169, 187, 188, 203, 205, 216, 217},
      {188, 189, 190, 204, 217, -1, -1},
      {190, 191, 207, 217, 218, -1, -1},
      {172, 191, 192, 206, 208, 218, 219},
      {192, 193, 207, 209, 219, -1, -1},
      {194, 195, 208, 210, 220, -1, -1},
      {177, 195, 196, 209, 211, 220, 221},  //210
      {196, 197, 210, 212, 221, -1, -1},
      {198, 199, 211, 213, 222, -1, -1},
      {182, 199, 200, 212, 214, 222, 223},
      {200, 201, 213, 215, 223, -1, -1},
      {202, 203, 214, 216, 224, -1, -1},
      {187, 203, 204, 215, 217, 224, 225},
      {204, 205, 206, 216, 225, -1, -1},
      {206, 207, 219, 225, -1, -1, -1},
      {207, 208, 218, 220, -1, -1, -1},
      {209, 210, 219, 221, -1, -1, -1},  //220
      {210, 211, 220, 222, -1, -1, -1},
      {212, 213, 221, 223, -1, -1, -1},
      {213, 214, 222, 224, -1, -1, -1},
      {215, 216, 223, 225, -1, -1, -1},
      {216, 217, 218, 224, -1, -1, -1}  //225
  };

  //=====ROOT file setting==========================================================
  //TString ofname =  Form("rootfiles/hist_dali%04d.root",FileNumber);
  TString ofname;

  if(argc < 4)
    ofname = Form("/home/koiwai/analysis/rootfiles/dali_hist/hist_dali%04d.root", FileNumber);
  else if(argc == 4)
    ofname = Form("/home/koiwai/analysis/dali/testhist_dali%04d.root", FileNumber);

  TFile *outfile = new TFile(ofname, "RECREATE");
  //TTree *tr = new TTree("anatrD","anatrD");
  //tr->SetAutoSave(1e5); //Autosave every 10,000 events filled.

  //=====Define variables===================================================
  Int_t EventNumber = 0;
  Int_t RunNumber = -1;

  Int_t dali_multi_ab;
  vector<Double_t> *dali_e_ab = new vector<Double_t>();
  vector<Int_t> *dali_id_ab = new vector<Int_t>();
  vector<Double_t> *dali_cos_ab = new vector<Double_t>();
  vector<Double_t> *dali_t_ab = new vector<Double_t>();
  vector<Double_t> *dali_x_ab = new vector<Double_t>();
  vector<Double_t> *dali_y_ab = new vector<Double_t>();
  vector<Double_t> *dali_z_ab = new vector<Double_t>();

  vector<Double_t> *dali_edop_ab = new vector<Double_t>();

  vector<Double_t> *dali_edop_simple_ab = new vector<Double_t>();
  vector<Double_t> *dali_edop_beta_ab = new vector<Double_t>();
  vector<Double_t> *dali_edop_theta_ab = new vector<Double_t>();

  //===== LOAD CUTS =====================================================================
  TFile *fcutSA_K = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_K.root");
  TFile *fcutSA_Ca_MINOS = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_Ca_wMINOS.root");
  TFile *fcutBR_Sc = TFile::Open("/home/koiwai/analysis/cutfiles/cutBR_Sc.root");
  TFile *fcutSA_50Ar = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_50Ar.root");
  TFile *fcutBR_51K = TFile::Open("/home/koiwai/analysis/cutfiles/cutBR_51K.root");

  TCutG *csa55k = (TCutG *)fcutSA_K->Get("sa55k");

  TCutG *csa53ca_minos = (TCutG *)fcutSA_Ca_MINOS->Get("csa53ca_wminos");

  TCutG *cbr54sc = (TCutG *)fcutBR_Sc->Get("br54sc");
  TCutG *cbr55sc = (TCutG *)fcutBR_Sc->Get("br55sc");
  TCutG *cbr56sc = (TCutG *)fcutBR_Sc->Get("br56sc");
  TCutG *cbr57sc = (TCutG *)fcutBR_Sc->Get("br57sc");
  TCutG *cbr58sc = (TCutG *)fcutBR_Sc->Get("br58sc");
  TCutG *cbr59sc = (TCutG *)fcutBR_Sc->Get("br59sc");

  TCutG *cbr51k = (TCutG *)fcutBR_51K->Get("br51k");
  TCutG *csa50ar = (TCutG *)fcutSA_50Ar->Get("sa50ar");

  //===== DEFINE HIST ====================================================================
  char *cnamebr[10] = {(char *)"br54ca", (char *)"br56ca", (char *)"br56ca", (char *)"br56sc", (char *)"br58sc", (char *)"br59sc", (char *)"br51k", (char *)"", (char *)"", (char *)""};
  char *cnamesa[10] = {(char *)"sa53ca", (char *)"sa55ca", (char *)"sa55k", (char *)"sa55ca", (char *)"sa57ca", (char *)"sa57ca", (char *)"sa50ar", (char *)"", (char *)"", (char *)""};
  char *hnames[10] = {(char *)"all", (char *)"m1", (char *)"m2", (char *)"m3", (char *)"mle3", (char *)"all_ab", (char *)"m1_ab", (char *)"m2_ab", (char *)"m3_ab", (char *)"mle3_ab"};
  TH1F *hdop[100];
  TH1F *hdopsimple[100];

  for(int i = 0; i < 7; i++) {
    for(int j = 0; j < 10; j++) {
      hdop[i * 10 + j] = new TH1F(Form("h_edop_%s_%s_%s", cnamebr[i], cnamesa[i], hnames[j]), Form("h_edop_%s_%s_%s", cnamebr[i], cnamesa[i], hnames[j]), 4000, 0, 4000);
    }
    for(int jj = 0; jj < 7; jj++) {
      hdopsimple[i * 5 + jj] = new TH1F(Form("h_edop_simple_%s_%s_%s", cnamebr[i], cnamesa[i], hnames[jj]), Form("h_edop_simple_%s_%s_%s", cnamebr[i], cnamesa[i], hnames[jj]), 4000, 0, 4000);
    }
  }

  //===== LOOP =========================================================================

  Int_t nEntry = intr->GetEntries();
  int iEntry = 0;
  int AllEntry;
  if(argc > 2 && MaxEventNumber < nEntry)
    AllEntry = MaxEventNumber;
  else
    AllEntry = nEntry;

  prepare_timer_tk();

  for(Int_t iEntry = 0; iEntry < AllEntry; iEntry++) {
    intr->GetEntry(iEntry);

    start_timer_tk(iEntry, AllEntry, 1000);

    RunNumber = FileNumber;
    EventNumber = EventNumber_calD;

    //*===== GATES =================================================================

    //*===== DALI timing gate ======================================================
    for(int i = 0; i < DALI_Multi; i++)
      if(abs(DALI_Time->at(i)) > 15.) continue;

    //*===== PID gate ==============================================================
    // done in cal_dali

    //*===== INIT ==================================================================
    dali_e_ab->clear();
    dali_id_ab->clear();
    dali_cos_ab->clear();
    dali_t_ab->clear();
    dali_x_ab->clear();
    dali_y_ab->clear();
    dali_z_ab->clear();
    dali_multi_ab = 0;

    dali_edop_ab->clear();
    dali_edop_simple_ab->clear();
    dali_edop_beta_ab->clear();
    dali_edop_theta_ab->clear();

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
      //tr->Fill();
      continue;
    }

    //===== ADD BACK END =====

    //===== Timing gate =====
    for(int i = 0; i < dali_multi_ab; i++) {
      if(-13.851 > dali_t_ab->at(i) || dali_t_ab->at(i) > 13.621) continue;  //5 sigma
    }
    if(!goodEvt) continue;

    //for(Int_t i=0;i<dali_multi_ab;i++){
    //  dali_pos->push_back(TVector3(10*dali_x_ab->at(i),10*dali_y_ab->at(i),10*dali_z_ab->at(i)));
    //  gamma_pos.push_back(dali_pos.at(i)-vertex);
    //  gamma_cos.push_back((gamma_pos.at(i)).CosTheta());
    //}

    //===== DOPPLER CORRECTION =====

    if(dali_multi_ab >= 1) {
      for(Int_t i = 0; i < dali_multi_ab; i++) {
        Double_t dali_edop_tmp = Sqrt(-1);
        dali_edop_tmp = dali_e_ab->at(i) * gamma_vertex * (1 - beta_vertex * gamma_cos->at(i));
        dali_edop_ab->push_back(dali_edop_tmp);
      }
    }

    //===== DOPPLER CORRECTION END =====

    //===== SINPLE DOPPLER CORRECTIOMN (WITHOUT MINOS )=====

    //beta_vertex_simple  = 0.5*(betaF7F13 + beta_minoshodo);
    //gamma_vertex_simple = 1/Sqrt(1 - beta_vertex_simple*beta_vertex_simple);
    const Double_t beta_mid = 0.57;
    const Double_t gamma_mid = 1 / Sqrt(1 - beta_mid * beta_mid);
    if(dali_multi_ab >= 1) {
      for(Int_t i = 0; i < dali_multi_ab; i++) {
        dali_edop_simple_ab->push_back(dali_e_ab->at(i) * gamma_mid * (1 - beta_mid * dali_cos_ab->at(i)));
        dali_edop_beta_ab->push_back(dali_e_ab->at(i) * gamma_vertex * (1 - beta_vertex * dali_cos_ab->at(i)));
        dali_edop_theta_ab->push_back(dali_e_ab->at(i) * gamma_mid * (1 - beta_mid * gamma_cos->at(i)));
      }
    }

    //===== SIMPLE DOPPLER CORRECTION END =====

    //===== FILL HIST ==============================================================================
    if(br54ca && csa53ca_minos->IsInside(aoqSA, zetSA)) {
      hdop[0]->Fill(dali_edop->at(0));
      if(DALI_Multi == 1)
        hdop[1]->Fill(dali_edop->at(0));
      if(DALI_Multi == 2)
        hdop[2]->Fill(dali_edop->at(0));
      if(DALI_Multi == 3)
        hdop[3]->Fill(dali_edop->at(0));
      if(DALI_Multi < 4)
        hdop[4]->Fill(dali_edop->at(0));

      hdop[5]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 1)
        hdop[6]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 2)
        hdop[7]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 3)
        hdop[8]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab < 4)
        hdop[9]->Fill(dali_edop_ab->at(0));
    }

    if(br56ca && sa55ca) {
      hdop[10]->Fill(dali_edop->at(0));
      if(DALI_Multi == 1)
        hdop[11]->Fill(dali_edop->at(0));
      if(DALI_Multi == 2)
        hdop[12]->Fill(dali_edop->at(0));
      if(DALI_Multi == 3)
        hdop[13]->Fill(dali_edop->at(0));
      if(DALI_Multi < 4)
        hdop[14]->Fill(dali_edop->at(0));

      hdop[15]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 1)
        hdop[16]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 2)
        hdop[17]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 3)
        hdop[18]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab < 4)
        hdop[19]->Fill(dali_edop_ab->at(0));
    }

    if(br56ca && csa55k->IsInside(aoqSA, zetSA)) {
      hdop[20]->Fill(dali_edop->at(0));
      if(DALI_Multi == 1)
        hdop[21]->Fill(dali_edop->at(0));
      if(DALI_Multi == 2)
        hdop[22]->Fill(dali_edop->at(0));
      if(DALI_Multi == 3)
        hdop[23]->Fill(dali_edop->at(0));
      if(DALI_Multi < 4)
        hdop[24]->Fill(dali_edop->at(0));

      hdop[25]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 1)
        hdop[26]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 2)
        hdop[27]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 3)
        hdop[28]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab < 4)
        hdop[29]->Fill(dali_edop_ab->at(0));
    }

    if(cbr56sc->IsInside(aoqBR, zetBR) && sa55ca) {
      hdop[30]->Fill(dali_edop->at(0));
      if(DALI_Multi == 1)
        hdop[31]->Fill(dali_edop->at(0));
      if(DALI_Multi == 2)
        hdop[32]->Fill(dali_edop->at(0));
      if(DALI_Multi == 3)
        hdop[33]->Fill(dali_edop->at(0));
      if(DALI_Multi < 4)
        hdop[34]->Fill(dali_edop->at(0));

      hdop[35]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 1)
        hdop[36]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 2)
        hdop[37]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 3)
        hdop[38]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab < 4)
        hdop[39]->Fill(dali_edop_ab->at(0));
    }

    if(cbr58sc->IsInside(aoqBR, zetBR) && sa57ca) {
      hdop[40]->Fill(dali_edop->at(0));
      if(DALI_Multi == 1)
        hdop[41]->Fill(dali_edop->at(0));
      if(DALI_Multi == 2)
        hdop[42]->Fill(dali_edop->at(0));
      if(DALI_Multi == 3)
        hdop[43]->Fill(dali_edop->at(0));
      if(DALI_Multi < 4)
        hdop[44]->Fill(dali_edop->at(0));

      hdop[45]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 1)
        hdop[46]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 2)
        hdop[47]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 3)
        hdop[48]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab < 4)
        hdop[49]->Fill(dali_edop_ab->at(0));
    }

    if(cbr59sc->IsInside(aoqBR, zetBR) && sa57ca) {
      hdop[50]->Fill(dali_edop->at(0));
      if(DALI_Multi == 1)
        hdop[51]->Fill(dali_edop->at(0));
      if(DALI_Multi == 2)
        hdop[52]->Fill(dali_edop->at(0));
      if(DALI_Multi == 3)
        hdop[53]->Fill(dali_edop->at(0));
      if(DALI_Multi < 4)
        hdop[54]->Fill(dali_edop->at(0));

      hdop[55]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 1)
        hdop[56]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 2)
        hdop[57]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 3)
        hdop[58]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab < 4)
        hdop[59]->Fill(dali_edop_ab->at(0));
    }

    if(cbr51k->IsInside(aoqBR, zetBR) && csa50ar->IsInside(aoqSA, zetSA)) {
      hdop[60]->Fill(dali_edop->at(0));
      if(DALI_Multi == 1)
        hdop[61]->Fill(dali_edop->at(0));
      if(DALI_Multi == 2)
        hdop[62]->Fill(dali_edop->at(0));
      if(DALI_Multi == 3)
        hdop[63]->Fill(dali_edop->at(0));
      if(DALI_Multi < 4)
        hdop[64]->Fill(dali_edop->at(0));

      hdop[65]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 1)
        hdop[66]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 2)
        hdop[67]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab == 3)
        hdop[68]->Fill(dali_edop_ab->at(0));
      if(dali_multi_ab < 4)
        hdop[69]->Fill(dali_edop_ab->at(0));
    }

    if(br54ca && csa53ca_minos->IsInside(aoqSA, zetSA)) {
      hdopsimple[0]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 1)
        hdopsimple[1]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 2)
        hdopsimple[2]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 3)
        hdopsimple[3]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi < 4)
        hdopsimple[4]->Fill(dali_edop_simple->at(0));
    }

    if(br56ca && sa55ca) {
      hdopsimple[5]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 1)
        hdopsimple[6]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 2)
        hdopsimple[7]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 3)
        hdopsimple[8]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi < 4)
        hdopsimple[9]->Fill(dali_edop_simple->at(0));
    }

    if(br56ca && csa55k->IsInside(aoqSA, zetSA)) {
      hdopsimple[10]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 1)
        hdopsimple[11]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 2)
        hdopsimple[12]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 3)
        hdopsimple[13]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi < 4)
        hdopsimple[14]->Fill(dali_edop_simple->at(0));
    }

    if(cbr56sc->IsInside(aoqBR, zetBR) && sa55ca) {
      hdopsimple[15]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 1)
        hdopsimple[16]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 2)
        hdopsimple[17]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 3)
        hdopsimple[18]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi < 4)
        hdopsimple[19]->Fill(dali_edop_simple->at(0));
    }

    if(cbr58sc->IsInside(aoqBR, zetBR) && sa57ca) {
      hdopsimple[20]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 1)
        hdopsimple[21]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 2)
        hdopsimple[22]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 3)
        hdopsimple[23]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi < 4)
        hdopsimple[24]->Fill(dali_edop_simple->at(0));
    }

    if(cbr59sc->IsInside(aoqBR, zetBR) && sa57ca) {
      hdopsimple[25]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 1)
        hdopsimple[26]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 2)
        hdopsimple[27]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 3)
        hdopsimple[28]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi < 4)
        hdopsimple[29]->Fill(dali_edop_simple->at(0));
    }

    if(cbr51k->IsInside(aoqBR, zetBR) && csa50ar->IsInside(aoqSA, zetSA)) {
      hdopsimple[30]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 1)
        hdopsimple[31]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 2)
        hdopsimple[32]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi == 3)
        hdopsimple[33]->Fill(dali_edop_simple->at(0));
      if(DALI_Multi < 4)
        hdopsimple[34]->Fill(dali_edop_simple->at(0));
    }

    //tr->Fill();

  }  //while loop
  std::clog << std::endl;

  outfile->cd();
  for(int i = 0; i < 70; i++)
    hdop[i]->Write();
  for(int i = 0; i < 35; i++)
    hdopsimple[i]->Write();
  outfile->Write();
  outfile->Close("R");

  delete DALI_ID;
  delete DALI_Time;
  delete DALI_Energy;
  delete DALI_CosTheta;
  delete DALI_X;
  delete DALI_Y;
  delete DALI_Z;

  delete dali_e_ab;
  delete dali_id_ab;
  delete dali_cos_ab;
  delete dali_t_ab;
  delete dali_x_ab;
  delete dali_y_ab;
  delete dali_z_ab;

  delete dali_edop_ab;
  delete dali_edop_simple_ab;
  delete dali_edop_beta_ab;
  delete dali_edop_theta_ab;

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
void SortDaliHit(Int_t left, Int_t right,vector <Int_t> *DALI_ID,vector <Double_t> *DALI_Energy, vector <Double_t> *DALI_EnergyDopplerCorrected, vector <Double_t> *DALI_Time, vector <Double_t> *DALI_CosTheta)
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
