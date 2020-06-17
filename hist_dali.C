#include "../include/anadalidef.h"

//&=====main Function==========================================================
int main(int argc, char *argv[]) {
  initiate_timer_tk();

  gInterpreter->GenerateDictionary("vector<TVector3>", "TVector3.h");

  bool   TestMode       = false;
  bool   ENum_flag      = false;
  int    FileNumber     = -1;
  int    MaxEventNumber = 0;
  double addbackRadius  = -1.;

  struct option longopts[] = {
      {"ab", required_argument, NULL, 'a'},
      {"eventnumber", required_argument, NULL, 'e'},
      {"runnumber", required_argument, NULL, 'r'},
      {"testmode", no_argument, NULL, 't'},
      {0, 0, 0, 0},
  };

  int opt;
  int longindex;

  while((opt = getopt_long(argc, argv, "a:tr:e:", longopts, &longindex)) != -1) {
    switch(opt) {
      case 'a':
        addbackRadius = atoi(optarg);
        break;
      case 'r':
        FileNumber = atoi(optarg);
        break;
      case 'e':
        ENum_flag      = true;
        MaxEventNumber = atoi(optarg);
        break;
      case 't':
        TestMode = true;
        break;

      default:
        break;
    }
  }

  if(argc == 1) {
    printf("Usage: ./hist_dali -r <run number> -e <max event number to treat> -t (to activate test mode)");
    exit(EXIT_FAILURE);
  }

  if(FileNumber == -1) {
    std::cerr << " You should provide a runnumber" << endl;
    exit(EXIT_FAILURE);
  }
  //if(argc < 2 || argc > 5) {
  //  printf("Usage: ./cal_dali RUNNUMBER\nOR     ./cal_dali RUNNUMBER MAXEVENTS\nOR     ./cal_dali RUNNUMBER MAXEVENTS TEST\n");
  //  exit(EXIT_FAILURE);
  //}

  //if(argc > 2) {
  //  MaxEventNumber = TString(argv[2]).Atoi();
  //  printf(" You will process %d events\n", MaxEventNumber);
  //}

  const double MINOSoffsetZ = 12.37;
  const double DALIoffset   = 96.5;

  //===== Load input files =====
  TString infnameB = Form("/home/koiwai/rootfiles/ana/beam/ana_beam%04d.root", FileNumber);
  TFile * infileB  = TFile::Open(infnameB);
  TTree * intrB    = (TTree *)infileB->Get("anatrB");

  Get_Branch_beam(intrB);

  TString infnameS = Form("/home/koiwai/rootfiles/ana/smri/ana_smri%04d.root", FileNumber);
  TFile * infileS  = TFile::Open(infnameS);
  TTree * intrS    = (TTree *)infileS->Get("anatrS");

  Get_Branch_smri(intrS);

  TString infnameV = Form("/home/koiwai/rootfiles/minos/vertex/vertex_frank%04d.root", FileNumber);
  TFile * infileV  = TFile::Open(infnameV);
  TTree * intrV    = (TTree *)infileV->Get("tr");

  Get_Branch_vertex(intrV);

  TString infname = Form("/home/koiwai/rootfiles/ana/dali/unpack_dali%04d.root", FileNumber);
  TFile * infile  = TFile::Open(infname);
  TTree * intr    = (TTree *)infile->Get("tr");

  Get_Branch_calD(intr);

  intr->AddFriend(intrB);
  intr->AddFriend(intrS);
  intr->AddFriend(intrV);

  //+=====ROOT file setting==========================================================
  //TString ofname =  Form("rootfiles/hist_dali%04d.root",FileNumber);
  TString ofname;

  if(!TestMode)
    ofname = Form("/home/koiwai/analysis/rootfiles/dali_hist/hist_dali%04d.root", FileNumber);
  else
    ofname = Form("/home/koiwai/analysis/dali/testhist_dali%04d.root", FileNumber);

  TFile *outfile = new TFile(ofname, "RECREATE");

  //=====Define variables===================================================
  //done in header

  //+===== LOAD CUTS =====================================================================
  TFile *fcutSA_K        = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_K.root");
  TFile *fcutSA_Ca_MINOS = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_Ca_wMINOS.root");
  TFile *fcutBR_Sc       = TFile::Open("/home/koiwai/analysis/cutfiles/cutBR_Sc.root");
  TFile *fcutSA_50Ar     = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_50Ar.root");
  TFile *fcutBR_51K      = TFile::Open("/home/koiwai/analysis/cutfiles/cutBR_51K.root");

  TCutG *csa55k = (TCutG *)fcutSA_K->Get("sa55k");

  TCutG *csa53ca_minos = (TCutG *)fcutSA_Ca_MINOS->Get("csa53ca_wminos");

  TCutG *cbr54sc = (TCutG *)fcutBR_Sc->Get("br54sc");
  TCutG *cbr55sc = (TCutG *)fcutBR_Sc->Get("br55sc");
  TCutG *cbr56sc = (TCutG *)fcutBR_Sc->Get("br56sc");
  TCutG *cbr57sc = (TCutG *)fcutBR_Sc->Get("br57sc");
  TCutG *cbr58sc = (TCutG *)fcutBR_Sc->Get("br58sc");
  TCutG *cbr59sc = (TCutG *)fcutBR_Sc->Get("br59sc");

  TCutG *cbr51k  = (TCutG *)fcutBR_51K->Get("br51k");
  TCutG *csa50ar = (TCutG *)fcutSA_50Ar->Get("sa50ar");

  //===== DEFINE HIST ====================================================================
  char *cnamebr[10] = {(char *)"br54ca", (char *)"br56ca", (char *)"br56ca", (char *)"br56sc", (char *)"br58sc",
                       (char *)"br59sc", (char *)"br51k", (char *)"", (char *)"", (char *)""};
  char *cnamesa[10] = {(char *)"sa53ca", (char *)"sa55ca", (char *)"sa55k", (char *)"sa55ca", (char *)"sa57ca",
                       (char *)"sa57ca", (char *)"sa50ar", (char *)"", (char *)"", (char *)""};
  char *hnames[10]  = {(char *)"all", (char *)"m1", (char *)"m2", (char *)"m3", (char *)"mle3", (char *)"all_ab",
                      (char *)"m1_ab", (char *)"m2_ab", (char *)"m3_ab", (char *)"mle3_ab"};
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
  int   iEntry = 0;
  int   AllEntry;
  if(argc > 2 && MaxEventNumber < nEntry)
    AllEntry = MaxEventNumber;
  else
    AllEntry = nEntry;

  prepare_timer_tk();

  //&===== LOOP =================================================================================
  for(Int_t iEntry = 0; iEntry < AllEntry; iEntry++) {
    intr->GetEntry(iEntry);

    start_timer_tk(iEntry, AllEntry, 1000);

    RunNumber   = FileNumber;
    EventNumber = EventNumber_calD;

    //+===== GATES =================================================================

    //+===== DALI timing gate ======================================================
    for(int i = 0; i < dali_multi; i++)
      if(abs(dali_t->at(i)) > 15.) continue;

    //+===== PID gate ==============================================================
    //- done in cal_dali

    //+===== INIT ==================================================================
    dali_e->clear();
    dali_t->clear();
    dali_cos->clear();
    dali_x->clear();
    dali_y->clear();
    dali_z->clear();
    dali_id->clear();
    dali_multi = 0;

    dali_edop->clear();
    dali_edop_simple->clear();
    dali_edop_beta->clear();
    dali_edop_theta->clear();

    dali_e_ab->clear();
    dali_t_ab->clear();
    dali_cos_ab->clear();
    dali_x_ab->clear();
    dali_y_ab->clear();
    dali_z_ab->clear();
    dali_id_ab->clear();
    dali_multi_ab = 0;

    dali_edop_ab->clear();
    dali_edop_simple_ab->clear();
    dali_edop_beta_ab->clear();
    dali_edop_theta_ab->clear();

    if(dali_multi > 1) {
      dali_e_ab->push_back(dali_e->at(0));
      dali_t_ab->push_back(dali_t->at(0));
      dali_cos_ab->push_back(dali_cos->at(0));
      dali_x_ab->push_back(dali_x->at(0));
      dali_y_ab->push_back(dali_y->at(0));
      dali_z_ab->push_back(dali_z->at(0));
      dali_id_ab->push_back(dali_id->at(0));

      //+===== ADD BACK =====

      Bool_t AddBack_flag = false;

      for(Int_t i = 1; i < dali_multi; i++) {
        AddBack_flag  = false;
        dali_multi_ab = dali_id_ab->size();
        for(Int_t k = 0; k < dali_multi_ab; k++) {
          for(Int_t j = 0; j < 7; j++) {  //TODOfor dali_e.size()
            //  if(dali_id->at(i) == AddBackTable[dali_id_ab->at(k)][j] && dali_Energy->at(i) > 300) {
            //    dali_e_ab->at(k) = dali_e_ab->at(k) + dali_Energy->at(i);
            //    AddBack_flag     = true;
            //    break;
            //  }
            //}
            if(AddBack_flag == true) break;
          }

          if(AddBack_flag == false) {
            dali_e_ab->push_back(dali_e->at(i));
            dali_t_ab->push_back(dali_t->at(i));
            dali_cos_ab->push_back(dali_cos->at(i));
            dali_x_ab->push_back(dali_x->at(i));
            dali_y_ab->push_back(dali_y->at(i));
            dali_z_ab->push_back(dali_z->at(i));
            dali_id_ab->push_back(dali_id->at(i));
          }
        }  //for dali_Multi

        dali_multi_ab = dali_id_ab->size();

      }  //dali_Multi>1

    } else if(dali_multi == 1) {
      dali_e_ab->push_back(dali_e->at(0));
      dali_t_ab->push_back(dali_t->at(0));
      dali_cos_ab->push_back(dali_cos->at(0));
      dali_x_ab->push_back(dali_x->at(0));
      dali_y_ab->push_back(dali_y->at(0));
      dali_z_ab->push_back(dali_z->at(0));
      dali_id_ab->push_back(dali_id->at(0));

      dali_multi_ab = 1;
    } else if(dali_multi == 0) {
      //tr->Fill();
      continue;
    }

    //-===== ADD BACK END =====

    //+===== Timing gate =====
    for(int i = 0; i < dali_multi_ab; i++) {
      if(-13.851 > dali_t_ab->at(i) || dali_t_ab->at(i) > 13.621) continue;  //5 sigma
    }
    //? if(!goodEvt) continue;

    //for(Int_t i=0;i<dali_multi_ab;i++){
    //  dali_pos->push_back(TVector3(10*dali_x_ab->at(i),10*dali_y_ab->at(i),10*dali_z_ab->at(i)));
    //  gamma_pos.push_back(dali_pos.at(i)-vertex);
    //  gamma_cos.push_back((gamma_pos.at(i)).CosTheta());
    //}

    //+===== DOPPLER CORRECTION =====

    //if(dali_multi_ab >= 1) {
    //  for(Int_t i = 0; i < dali_multi_ab; i++) {
    //    Double_t dali_edop_tmp = Sqrt(-1);
    //    dali_edop_tmp          = dali_e_ab->at(i) * gamma_vertex * (1 - beta_vertex * gamma_cos->at(i));
    //    dali_edop_ab->push_back(dali_edop_tmp);
    //  }
    //}

    //-===== DOPPLER CORRECTION END =====

    //+===== SINPLE DOPPLER CORRECTIOMN (WITHOUT MINOS )=====

    //beta_vertex_simple  = 0.5*(betaF7F13 + beta_minoshodo);
    //gamma_vertex_simple = 1/Sqrt(1 - beta_vertex_simple*beta_vertex_simple);
    //const Double_t beta_mid  = 0.57;
    //const Double_t gamma_mid = 1 / Sqrt(1 - beta_mid * beta_mid);
    //if(dali_multi_ab >= 1) {
    //  for(Int_t i = 0; i < dali_multi_ab; i++) {
    //    dali_edop_simple_ab->push_back(dali_e_ab->at(i) * gamma_mid * (1 - beta_mid * dali_cos_ab->at(i)));
    //    dali_edop_beta_ab->push_back(dali_e_ab->at(i) * gamma_vertex * (1 - beta_vertex * dali_cos_ab->at(i)));
    //    dali_edop_theta_ab->push_back(dali_e_ab->at(i) * gamma_mid * (1 - beta_mid * gamma_cos->at(i)));
    //  }
    //}

    //-===== SIMPLE DOPPLER CORRECTION END =====

    //+===== FILL HIST ==============================================================================
    if(br54ca && csa53ca_minos->IsInside(aoqSA, zetSA)) {
      hdop[0]->Fill(dali_edop->at(0));
      if(dali_multi == 1)
        hdop[1]->Fill(dali_edop->at(0));
      if(dali_multi == 2)
        hdop[2]->Fill(dali_edop->at(0));
      if(dali_multi == 3)
        hdop[3]->Fill(dali_edop->at(0));
      if(dali_multi < 4)
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
      if(dali_multi == 1)
        hdop[11]->Fill(dali_edop->at(0));
      if(dali_multi == 2)
        hdop[12]->Fill(dali_edop->at(0));
      if(dali_multi == 3)
        hdop[13]->Fill(dali_edop->at(0));
      if(dali_multi < 4)
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
      if(dali_multi == 1)
        hdop[21]->Fill(dali_edop->at(0));
      if(dali_multi == 2)
        hdop[22]->Fill(dali_edop->at(0));
      if(dali_multi == 3)
        hdop[23]->Fill(dali_edop->at(0));
      if(dali_multi < 4)
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
      if(dali_multi == 1)
        hdop[31]->Fill(dali_edop->at(0));
      if(dali_multi == 2)
        hdop[32]->Fill(dali_edop->at(0));
      if(dali_multi == 3)
        hdop[33]->Fill(dali_edop->at(0));
      if(dali_multi < 4)
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
      if(dali_multi == 1)
        hdop[41]->Fill(dali_edop->at(0));
      if(dali_multi == 2)
        hdop[42]->Fill(dali_edop->at(0));
      if(dali_multi == 3)
        hdop[43]->Fill(dali_edop->at(0));
      if(dali_multi < 4)
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
      if(dali_multi == 1)
        hdop[51]->Fill(dali_edop->at(0));
      if(dali_multi == 2)
        hdop[52]->Fill(dali_edop->at(0));
      if(dali_multi == 3)
        hdop[53]->Fill(dali_edop->at(0));
      if(dali_multi < 4)
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
      if(dali_multi == 1)
        hdop[61]->Fill(dali_edop->at(0));
      if(dali_multi == 2)
        hdop[62]->Fill(dali_edop->at(0));
      if(dali_multi == 3)
        hdop[63]->Fill(dali_edop->at(0));
      if(dali_multi < 4)
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
      if(dali_multi == 1)
        hdopsimple[1]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 2)
        hdopsimple[2]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 3)
        hdopsimple[3]->Fill(dali_edop_simple->at(0));
      if(dali_multi < 4)
        hdopsimple[4]->Fill(dali_edop_simple->at(0));
    }

    if(br56ca && sa55ca) {
      hdopsimple[5]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 1)
        hdopsimple[6]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 2)
        hdopsimple[7]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 3)
        hdopsimple[8]->Fill(dali_edop_simple->at(0));
      if(dali_multi < 4)
        hdopsimple[9]->Fill(dali_edop_simple->at(0));
    }

    if(br56ca && csa55k->IsInside(aoqSA, zetSA)) {
      hdopsimple[10]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 1)
        hdopsimple[11]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 2)
        hdopsimple[12]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 3)
        hdopsimple[13]->Fill(dali_edop_simple->at(0));
      if(dali_multi < 4)
        hdopsimple[14]->Fill(dali_edop_simple->at(0));
    }

    if(cbr56sc->IsInside(aoqBR, zetBR) && sa55ca) {
      hdopsimple[15]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 1)
        hdopsimple[16]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 2)
        hdopsimple[17]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 3)
        hdopsimple[18]->Fill(dali_edop_simple->at(0));
      if(dali_multi < 4)
        hdopsimple[19]->Fill(dali_edop_simple->at(0));
    }

    if(cbr58sc->IsInside(aoqBR, zetBR) && sa57ca) {
      hdopsimple[20]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 1)
        hdopsimple[21]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 2)
        hdopsimple[22]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 3)
        hdopsimple[23]->Fill(dali_edop_simple->at(0));
      if(dali_multi < 4)
        hdopsimple[24]->Fill(dali_edop_simple->at(0));
    }

    if(cbr59sc->IsInside(aoqBR, zetBR) && sa57ca) {
      hdopsimple[25]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 1)
        hdopsimple[26]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 2)
        hdopsimple[27]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 3)
        hdopsimple[28]->Fill(dali_edop_simple->at(0));
      if(dali_multi < 4)
        hdopsimple[29]->Fill(dali_edop_simple->at(0));
    }

    if(cbr51k->IsInside(aoqBR, zetBR) && csa50ar->IsInside(aoqSA, zetSA)) {
      hdopsimple[30]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 1)
        hdopsimple[31]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 2)
        hdopsimple[32]->Fill(dali_edop_simple->at(0));
      if(dali_multi == 3)
        hdopsimple[33]->Fill(dali_edop_simple->at(0));
      if(dali_multi < 4)
        hdopsimple[34]->Fill(dali_edop_simple->at(0));
    }

  }  //-while loop

  std::clog << std::endl;

  outfile->cd();
  for(int i = 0; i < 70; i++)
    hdop[i]->Write();
  for(int i = 0; i < 35; i++)
    hdopsimple[i]->Write();

  outfile->Write();
  outfile->Close("R");

  delete dali_e;
  delete dali_t;
  delete dali_cos;
  delete dali_x;
  delete dali_y;
  delete dali_z;
  delete dali_id;

  delete dali_edop;
  delete dali_edop_simple;
  delete dali_edop_beta;
  delete dali_edop_theta;

  delete dali_e_ab;
  delete dali_t_ab;
  delete dali_cos_ab;
  delete dali_x_ab;
  delete dali_y_ab;
  delete dali_z_ab;
  delete dali_id_ab;

  delete dali_edop_ab;
  delete dali_edop_simple_ab;
  delete dali_edop_beta_ab;
  delete dali_edop_theta_ab;

  stop_timer_tk(FileNumber, AllEntry);

  return 0;
}  //main()

void SortDaliHit(Int_t left, Int_t right, vector<Int_t> *DALI_ID, vector<Double_t> *DALI_Energy, vector<Double_t> *DALI_Time, vector<Double_t> *DALI_X, vector<Double_t> *DALI_Y, vector<Double_t> *DALI_Z, vector<Double_t> *DALI_CosTheta) {
  Int_t    TempID;
  Double_t TempEnergy;
  Double_t TempTime;
  Double_t TempX;
  Double_t TempY;
  Double_t TempZ;
  Double_t TempCosTheta;

  int    i = left, j = right;
  double pivot = DALI_Energy->at((left + right) / 2);

  //--partition---/
  while(i <= j) {
    while(DALI_Energy->at(i) > pivot)
      i++;
    while(DALI_Energy->at(j) < pivot)
      j--;
    if(i <= j) {
      TempID       = DALI_ID->at(j);
      TempEnergy   = DALI_Energy->at(j);
      TempTime     = DALI_Time->at(j);
      TempX        = DALI_X->at(j);
      TempY        = DALI_Y->at(j);
      TempZ        = DALI_Z->at(j);
      TempCosTheta = DALI_CosTheta->at(j);

      DALI_ID->at(j)       = DALI_ID->at(i);
      DALI_Energy->at(j)   = DALI_Energy->at(i);
      DALI_Time->at(j)     = DALI_Time->at(i);
      DALI_X->at(j)        = DALI_X->at(i);
      DALI_Y->at(j)        = DALI_Y->at(i);
      DALI_Z->at(j)        = DALI_Z->at(i);
      DALI_CosTheta->at(j) = DALI_CosTheta->at(i);

      DALI_ID->at(i)       = TempID;
      DALI_Energy->at(i)   = TempEnergy;
      DALI_Time->at(i)     = TempTime;
      DALI_X->at(i)        = TempX;
      DALI_Y->at(i)        = TempY;
      DALI_Z->at(i)        = TempZ;
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
  Double_t Gamma            = 1 / TMath::Sqrt(1 - Beta * Beta);
  Double_t DopplerCorrected = GammaDopplerEnergy * Gamma * (1 - Beta * CosTheta);
  return DopplerCorrected;
}

inline bool exists_test(const std::string &name) {
  return (access(name.c_str(), F_OK) != -1);
}
inline bool exists_test(const TString &name) {
  return (access(name.Data(), F_OK) != -1);
}
