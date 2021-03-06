#include "/home/koiwai/analysis/include/anadalidef.h"

//&=====main Function==========================================================
int main(int argc, char *argv[]) {
  initiate_timer_tk();

  gInterpreter->GenerateDictionary("vector<TVector3>", "TVector3.h");

  bool   TestMode       = false;
  bool   ENum_flag      = false;
  int    FileNumber     = -1;
  int    MaxEventNumber = 0;
  double addbackRadius  = -1.;
  bool   addback_flag   = false;
  double DALITimegateUp = 100.;

  struct option longopts[] = {
      {"ab", required_argument, NULL, 'a'},
      {"eventnumber", required_argument, NULL, 'e'},
      {"runnumber", required_argument, NULL, 'r'},
      {"testmode", no_argument, NULL, 't'},
      {"timegateup", required_argument, NULL, 'u'},
      {0, 0, 0, 0},
  };

  int opt;
  int longindex;

  while((opt = getopt_long(argc, argv, "a:tr:e:u:", longopts, &longindex)) != -1) {
    switch(opt) {
      case 'a':
        addback_flag  = true;
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
      case 'u':
        DALITimegateUp = atoi(optarg);
        break;

      default:
        break;
    }
  }

  if(argc == 1) {
    printf("\nUsage: ./hist_dali \n -r <run number> \n -e <max event number to treat> \n -t (to activate test mode)\n -a <addback radius>\n -u <time gate upper limit> \n");
    exit(EXIT_FAILURE);
  }

  if(FileNumber == -1) {
    std::cerr << " You should provide a runnumber" << endl;
    exit(EXIT_FAILURE);
  }

  //+===== Define constants =====

  const double MINOSoffsetZ = 12.37;
  const double DALIoffset   = 96.5;
  const double DALITimeGate = 15.;
  const double TPClength    = 151.;

  //+===== Load input files =====
  TString infnameB = Form("/home/koiwai/analysis/rootfiles/ana/beam/ana_beam%04d.root", FileNumber);
  TFile * infileB  = TFile::Open(infnameB);
  TTree * intrB    = (TTree *)infileB->Get("anatrB");

  Get_Branch_beam(intrB);

  TString infnameS = Form("/home/koiwai/analysis/rootfiles/ana/smri/ana_smri%04d.root", FileNumber);
  TFile * infileS  = TFile::Open(infnameS);
  TTree * intrS    = (TTree *)infileS->Get("anatrS");

  Get_Branch_smri(intrS);

  TString infnameV = Form("/home/koiwai/analysis/rootfiles/minos/vertex/vertex_frank%04d.root", FileNumber);
  TFile * infileV  = TFile::Open(infnameV);
  TTree * intrV    = (TTree *)infileV->Get("tr");

  Get_Branch_vertex(intrV);

  TString infname = Form("/home/koiwai/analysis/rootfiles/ana/dali/unpack_dali%04d.root", FileNumber);
  TFile * infile  = TFile::Open(infname);
  TTree * intr    = (TTree *)infile->Get("tr");

  Get_Branch_calD(intr);

  intr->AddFriend(intrB);
  intr->AddFriend(intrS);
  intr->AddFriend(intrV);

  //+=====ROOT file setting==========================================================
  TString ofname;

  if(!TestMode)
    ofname = Form("/home/koiwai/analysis/rootfiles/dali_hist/hist_dali%04d.root", FileNumber);
  else
    ofname = Form("/home/koiwai/analysis/dali/testhist_dali%04d.root", FileNumber);

  TFile *outfile = new TFile(ofname, "RECREATE");

  //+=====Define variables===================================================
  //- done in header

  //+===== LOAD CUTS =====================================================================

  TFile *fcutSA_Ca_MINOS = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_Ca_wMINOS.root");
  TFile *fcutSA_Ca       = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_Ca.root");
  TFile *fcutSA_K        = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_K.root");
  TFile *fcutSA_50Ar     = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_50Ar.root");

  TCutG *csa53ca_minos = (TCutG *)fcutSA_Ca_MINOS->Get("csa53ca_wminos");
  TCutG *csa57ca       = (TCutG *)fcutSA_Ca->Get("sa57ca");
  TCutG *csa55ca       = (TCutG *)fcutSA_Ca->Get("sa55ca");
  TCutG *csa55k        = (TCutG *)fcutSA_K->Get("sa55k");
  TCutG *csa53k        = (TCutG *)fcutSA_K->Get("sa53k");
  TCutG *csa51k        = (TCutG *)fcutSA_K->Get("sa51k");
  TCutG *csa50ar       = (TCutG *)fcutSA_50Ar->Get("sa50ar");

  //+===== DEFINE HIST ====================================================================

  //char *cnamech[11] = {(char *)"br54ca_sa53ca",
  //                     (char *)"br51k_sa50ar",
  //                     (char *)"br52ca_sa51k",
  //                     (char *)"br54ca_sa53k",
  //                     (char *)"br56ca_sa55k",
  //                     (char *)"br56ca_sa55ca",
  //                     (char *)"br56sc_sa55ca",
  //                     (char *)"br58sc_sa57ca",
  //                     (char *)"br59sc_sa57ca",
  //                     (char *)"br59ti_sa57ca",
  //                     (char *)"br60ti_sa57ca"};

  char *cnamech[5] = {(char *)"br56ca_sa55k",
                      (char *)"br56ca_sa55ca",
                      (char *)"br56sc_sa55ca",
                      (char *)"br58sc_sa57ca",
                      (char *)"br59sc_sa57ca"};

  int gatenum = 5;

  char *hnames[10] = {(char *)"all",
                      (char *)"m1",
                      (char *)"m2",
                      (char *)"m3",
                      (char *)"mle4",
                      (char *)"all_ab",
                      (char *)"m1_ab",
                      (char *)"m2_ab",
                      (char *)"m3_ab",
                      (char *)"mle4_ab"};

  TH1F *hdop[50];
  TH1F *hdopafter[50];
  TH1F *hdopsimple[50];

  TH2F *hdopggcoin[11];
  TH2F *hdopggsym[11];

  TH2F *hdopmult[11];
  TH2F *hdoptime[11];
  TH2F *hdopID[11];

  TH1F *hbetaF7F13[11];
  TH1F *hbetaTH[11];
  TH1F *hzvertex[11];
  TH1F *hbetavertex[11];

  TH1F *htransbr[11];
  TH1F *htransbrsa[11];

  for(int i = 0; i < gatenum; i++) {
    hdopmult[i] = new TH2F(
        Form("h_edopmult_%s", cnamech[i]),
        Form("Multiplicity vs Doppler corrected E (%s)", cnamech[i]),
        10, 0, 10, 4000, 0, 4000);

    hdoptime[i] = new TH2F(
        Form("h_edoptime_%s", cnamech[i]),
        Form("DALI_Time vs Doppler corrected E (%s)", cnamech[i]),
        300, -50, 50, 4000, 0, 4000);

    hdopID[i] = new TH2F(
        Form("h_edopID_%s", cnamech[i]),
        Form("DALI ID vs Doppler corrected E (%s)", cnamech[i]),
        230, 0, 230, 4000, 0, 4000);

    for(int j = 0; j < 10; j++) {
      hdop[i * 10 + j]      = new TH1F(  //TODO histo number to be re-considered
          Form("h_edop_%s_%s", cnamech[i], hnames[j]),
          Form("h_edop_%s_%s", cnamech[i], hnames[j]),
          4000, 0, 4000);
      hdopafter[i * 10 + j] = new TH1F(  //TODO histo number to be re-considered
          Form("h_edopafter_%s_%s", cnamech[i], hnames[j]),
          Form("h_edopafter_%s_%s", cnamech[i], hnames[j]),
          4000, 0, 4000);

      hdopsimple[i * 10 + j] = new TH1F(
          Form("h_edop_simple_%s_%s", cnamech[i], hnames[j]),
          Form("h_edop_simple_%s_%s", cnamech[i], hnames[j]),
          4000, 0, 4000);
    }

    hdopggcoin[i] = new TH2F(
        Form("h_edop_ggcoin_%s", cnamech[i]),
        Form("GG coin: First hit vs else (%s)", cnamech[i]),
        4000, 0, 4000, 4000, 0, 4000);

    hdopggsym[i] = new TH2F(
        Form("h_edop_ggsym_%s", cnamech[i]),
        Form("GG coin: G1&G2 vs G2&G1 (%s)", cnamech[i]),
        4000, 0, 4000, 4000, 0, 4000);

    hbetaF7F13[i] = new TH1F(
        Form("h_betaF7F13_%s", cnamech[i]),
        Form("h_betaF7F13_%s", cnamech[i]),
        500, 0.3, 0.8);
    hbetaTH[i] = new TH1F(
        Form("h_betaTH_%s", cnamech[i]),
        Form("h_betaTH_%s", cnamech[i]),
        500, 0.3, 0.8);
    hzvertex[i] = new TH1F(
        Form("h_zvertex_%s", cnamech[i]),
        Form("h_zvertex_%s", cnamech[i]),
        150, 0, 150);
    hbetavertex[i] = new TH1F(
        Form("h_betavertex_%s", cnamech[i]),
        Form("h_betavertex_%s", cnamech[i]),
        30, 0.4, 0.7);

    htransbr[i] = new TH1F(
        Form("h_trans_BR_%s", cnamech[i]),
        Form("F5X (gated on BR for %s)", cnamech[i]),
        100, -150, 150);
    htransbrsa[i] = new TH1F(
        Form("h_trans_BRSA_%s", cnamech[i]),
        Form("F5X (gated on BR and SA for %s)", cnamech[i]),
        100, -150, 150);
  }

  TH1F *hdopFW[4];

  hdopFW[0] = new TH1F("h_edopFW31_br56ca_sa55ca", " BR56Ca -> SA55Ca (DALI ID > 31)", 4000, 0, 4000);
  hdopFW[1] = new TH1F("h_edopFW45_br56ca_sa55ca", " BR56Ca -> SA55Ca (DALI ID > 45)", 4000, 0, 4000);
  hdopFW[2] = new TH1F("h_edopFW67_br56ca_sa55ca", " BR56Ca -> SA55Ca (DALI ID > 67)", 4000, 0, 4000);
  hdopFW[3] = new TH1F("h_edopFW85_br56ca_sa55ca", " BR56Ca -> SA55Ca (DALI ID > 85)", 4000, 0, 4000);

  //TH1F *h_minoseff_50ar[3];
  //TH1F *h_minoseff_53ca[3];
  //
  //h_minoseff_50ar[0] = new TH1F(
  //    "h_edop_simple_br51k_sa50ar_all_wovertex",
  //    "50Ar: MINOS effciency (simple edop)",
  //    4000, 0, 4000);
  //h_minoseff_50ar[1] = new TH1F(
  //    "h_edop_simple_br51k_sa50ar_all_wallvertex",
  //    "50Ar: MINOS effciency (simple edop plus all MINOS verteces reco.ed)",
  //    4000, 0, 4000);
  //h_minoseff_50ar[2] = new TH1F(
  //    "h_edop_simple_br51k_sa50ar_all_w1or2vertex",
  //    "50Ar: MINOS effciency (simple edop plus 1 or 2 MINOS vertex reco.ed)",
  //    4000, 0, 4000);
  //h_minoseff_53ca[0] = new TH1F(
  //    "h_edop_simple_br54ca_sa53ca_all_wovertex",
  //    "53Ca: MINOS effciency (simple edop)",
  //    4000, 0, 4000);
  //h_minoseff_53ca[1] = new TH1F(
  //    "h_edop_simple_br54ca_sa53ca_all_wallvertex",
  //    "53Ca: MINOS effciency (simple edop plus all MINOS verteces reco.ed)",
  //    4000, 0, 4000);
  //h_minoseff_53ca[2] = new TH1F(
  //    "h_edop_simple_br54ca_sa53ca_all_w1or2vertex",
  //    "53Ca: MINOS effciency (simple edop plus 1 or 2 MINOS vertex reco.ed)",
  //    4000, 0, 4000);

  TH1F *h_dalit     = new TH1F("h_dalit", "DALI time of first hit", 300, -50, 50);
  TH1F *h_dalit_all = new TH1F("h_dalit_all", "DALI time of all hits", 300, -50, 50);

  TH1F *hMINOSZ      = new TH1F("MINOS_Z_cor", "MINOS_Z_cor", 250, -50, 200);
  TH1F *hbeta_vertex = new TH1F("beta_vertex", "beta_vertex", 100, 0, 1);
  TH1F *hgamma_cos   = new TH1F("gamma_cos", "gamma_cos", 100, -1.1, 1.1);
  TH2F *hIDtheta     = new TH2F("h_IDtheta", "DALI ID vs DALI theta", 230, 0, 230, 100, -1.1, 1.1);
  TH1F *hNumTracks   = new TH1F("MINOS_NumberTracks", "MINOS NumberTracks", 10, 0, 10);

  //&===== LOOP =========================================================================

  Int_t nEntry = intr->GetEntries();
  int   iEntry = 0;
  int   AllEntry;

  if(ENum_flag && MaxEventNumber < nEntry)
    AllEntry = MaxEventNumber;
  else
    AllEntry = nEntry;

  prepare_timer_tk();

  for(Int_t iEntry = 0; iEntry < AllEntry; iEntry++) {
    intr->GetEntry(iEntry);

    start_timer_tk(iEntry, AllEntry, 1000);

    RunNumber   = FileNumber;
    EventNumber = EventNumber_calD;

    //+===== GATES =================================================================
    if(!br59ti && !br60ti && !br59sc && !br58sc && !br56sc && !br56ca && !br54ca && !br52ca && !br51k)
      continue;

    if(dali_multi == 0)
      continue;

    //+===== PID gate ==============================================================
    //bool PIDgates[11] = {
    //    br54ca && csa53ca_minos->IsInside(aoqSA, zetSA),
    //    br51k && csa50ar->IsInside(aoqSA, zetSA),
    //    br52ca && csa51k->IsInside(aoqSA, zetSA),
    //    br54ca && csa53k->IsInside(aoqSA, zetSA),
    //    br56ca && csa55k->IsInside(aoqSA, zetSA),
    //    br56ca && csa55ca->IsInside(aoqSA, zetSA),
    //    br56sc && csa55ca->IsInside(aoqSA, zetSA),
    //    br58sc && csa57ca->IsInside(aoqSA, zetSA),
    //    br59sc && csa57ca->IsInside(aoqSA, zetSA),
    //    br59ti && csa57ca->IsInside(aoqSA, zetSA),
    //    br60ti && csa57ca->IsInside(aoqSA, zetSA),
    //};

    bool PIDgates[5] = {
        br56ca && csa55k->IsInside(aoqSA, zetSA),
        br56ca && csa55ca->IsInside(aoqSA, zetSA),
        br56sc && csa55ca->IsInside(aoqSA, zetSA),
        br58sc && csa57ca->IsInside(aoqSA, zetSA),
        br59sc && csa57ca->IsInside(aoqSA, zetSA),
    };

    //+===== INIT ==================================================================
    dali_edop->clear();
    dali_edopafter->clear();
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
    dali_edopafter_ab->clear();
    dali_edop_simple_ab->clear();
    dali_edop_beta_ab->clear();
    dali_edop_theta_ab->clear();

    gamma_cos.clear();
    gamma_cosafter.clear();
    gamma_cos_ab.clear();
    gamma_cosafter_ab.clear();

    vertex.SetXYZ(Sqrt(-1), Sqrt(-1), Sqrt(-1));

    //+===== calculate velocities =========================================
    MINOS_Z_cor = MINOS_Z + MINOSoffsetZ;
    beta_vertex = betaF7F13 - (betaF7F13 - betaTH) * MINOS_Z_cor / TPClength;
    vertex.SetXYZ(MINOS_X, MINOS_Y, MINOS_Z_cor - DALIoffset);

    //+===== Create DALI crystal vector ==========================================
    for(int i = 0; i < dali_multi; i++) {
      TVector3 gamma_pos_tmp(10 * dali_x->at(i) - MINOS_X, 10 * dali_y->at(i) - MINOS_Y, 10 * dali_z->at(i) - vertex.Z());
      gamma_cos.push_back(gamma_pos_tmp.CosTheta());
      TVector3 gamma_posafter_tmp(10 * dali_x->at(i) - MINOS_X, 10 * dali_y->at(i) - MINOS_Y, 10 * dali_z->at(i) - 150. - DALIoffset);
      gamma_cos.push_back(gamma_pos_tmp.CosTheta());
      gamma_cosafter.push_back(gamma_posafter_tmp.CosTheta());
    }

    //+===== ADD BACK ============================================================
    int       DALI_NClust                                 = 0;
    Double_t  addbackThreshold                            = 40.;  //! keV
    const int NUMBEROFDALICRYSTALS                        = 226;
    bool      crystalUsedForAddback[NUMBEROFDALICRYSTALS] = {false};
    double    DUMM_Energy[NUMBEROFDALICRYSTALS]           = {sqrt(-1.)};
    double    DUMM_Time[NUMBEROFDALICRYSTALS]             = {sqrt(-1.)};
    double    DUMM_Cos[NUMBEROFDALICRYSTALS]              = {sqrt(-1.)};
    double    DUMM_X[NUMBEROFDALICRYSTALS]                = {sqrt(-1.)};
    double    DUMM_Y[NUMBEROFDALICRYSTALS]                = {sqrt(-1.)};
    double    DUMM_Z[NUMBEROFDALICRYSTALS]                = {sqrt(-1.)};
    int       DUMM_ID[NUMBEROFDALICRYSTALS]               = {-1};

    if(addback_flag) {
      for(int j = 0; j < dali_multi; j++) {
        if(crystalUsedForAddback[j] == true)
          continue;
        DUMM_Energy[DALI_NClust] = dali_e->at(j);
        DUMM_Time[DALI_NClust]   = dali_t->at(j);
        DUMM_Cos[DALI_NClust]    = dali_cos->at(j);
        DUMM_X[DALI_NClust]      = dali_x->at(j);
        DUMM_Y[DALI_NClust]      = dali_y->at(j);
        DUMM_Z[DALI_NClust]      = dali_z->at(j);
        DUMM_ID[DALI_NClust]     = dali_id->at(j);
        crystalUsedForAddback[j] = true;
        for(unsigned int k = j + 1; k < dali_id->size(); k++) {
          if(crystalUsedForAddback[k] == true)
            continue;
          TVector3 dali_pos_tmp(dali_x->at(j) - dali_x->at(k), dali_y->at(j) - dali_y->at(k), dali_z->at(j) - dali_z->at(k));
          if(dali_pos_tmp.Mag() < addbackRadius && dali_e->at(k) > addbackThreshold && (-5 < dali_t->at(k) && dali_t->at(k) < DALITimegateUp)) {
            DUMM_Energy[DALI_NClust] += dali_e->at(k);
            crystalUsedForAddback[k] = true;
          }
        }
        DALI_NClust++;
      }

      dali_multi_ab = DALI_NClust;
      for(int j = 0; j < DALI_NClust; j++) {
        dali_e_ab->push_back(DUMM_Energy[j]);
        dali_t_ab->push_back(DUMM_Time[j]);
        dali_x_ab->push_back(DUMM_X[j]);
        dali_y_ab->push_back(DUMM_Y[j]);
        dali_z_ab->push_back(DUMM_Z[j]);
        dali_cos_ab->push_back(DUMM_Cos[j]);
        dali_id_ab->push_back(DUMM_ID[j]);
      }
    }

    if(dali_multi_ab > 1) {
      SortDaliHit(0, dali_multi_ab - 1,
                  dali_id_ab,
                  dali_e_ab,
                  dali_t_ab,
                  dali_x_ab,
                  dali_y_ab,
                  dali_z_ab,
                  dali_cos_ab);
    }

    //-===== ADD BACK END =====

    //? if(!goodEvt) continue;

    //+===== Reconstruct gamma-ray vector (cos theta) =====
    for(int i = 0; i < dali_multi_ab; i++) {
      TVector3 gamma_pos_ab_tmp(10 * dali_x->at(i) - MINOS_X, 10 * dali_y->at(i) - MINOS_Y, 10 * dali_z->at(i) - vertex.Z());
      gamma_cos_ab.push_back(gamma_pos_ab_tmp.CosTheta());
      TVector3 gamma_posafter_ab_tmp(10 * dali_x->at(i) - MINOS_X, 10 * dali_y->at(i) - MINOS_Y, 10 * dali_z->at(i) - 150. - DALIoffset);
      gamma_cosafter_ab.push_back(gamma_posafter_ab_tmp.CosTheta());
    }

    //+===== DOPPLER CORRECTION =====
    if(dali_multi_ab >= 1) {
      for(Int_t i = 0; i < dali_multi; i++) {
        dali_edop->push_back(DopplerCorrection(dali_e->at(i), beta_vertex, gamma_cos.at(i)));
        dali_edopafter->push_back(DopplerCorrection(dali_e->at(i), betaTH, gamma_cosafter.at(i)));
      }
      for(Int_t i = 0; i < dali_multi_ab; i++) {
        dali_edop_ab->push_back(DopplerCorrection(dali_e_ab->at(i), beta_vertex, gamma_cos_ab.at(i)));
        dali_edopafter_ab->push_back(DopplerCorrection(dali_e_ab->at(i), betaTH, gamma_cosafter_ab.at(i)));
      }
    }

    //-===== DOPPLER CORRECTION END =====

    //+===== SINPLE DOPPLER CORRECTIOMN (WITHOUT MINOS )=====
    beta_vertex_simple      = 0.5 * (betaF7F13 + betaTH);
    const Double_t beta_mid = 0.57;

    if(dali_multi_ab >= 1) {
      for(Int_t i = 0; i < dali_multi; i++) {
        dali_edop_simple->push_back(DopplerCorrection(dali_e->at(i), beta_vertex_simple, dali_cos->at(i)));
        dali_edop_beta->push_back(DopplerCorrection(dali_e->at(i), beta_vertex, dali_cos->at(i)));
        dali_edop_theta->push_back(DopplerCorrection(dali_e->at(i), beta_vertex_simple, gamma_cos.at(i)));
      }
      for(Int_t i = 0; i < dali_multi_ab; i++) {
        dali_edop_simple_ab->push_back(DopplerCorrection(dali_e_ab->at(i), beta_vertex_simple, dali_cos_ab->at(i)));
        dali_edop_beta_ab->push_back(DopplerCorrection(dali_e_ab->at(i), beta_vertex, dali_cos_ab->at(i)));
        dali_edop_theta_ab->push_back(DopplerCorrection(dali_e_ab->at(i), beta_vertex_simple, gamma_cos_ab.at(i)));
      }
    }

    //-===== SIMPLE DOPPLER CORRECTION END =====

    //+===== DEBUG =====
    //if(dali_edop->size() != dali_edop_simple->size()) {
    //  cout << "dali_edop->size()!= dali_edop_simple->size()" << endl;
    //  exit(EXIT_FAILURE);
    //}
    //if(dali_edop_ab->size() != dali_edop_simple_ab->size()) {
    //  cout << "dali_edop_ab->size()!= dali_edop_simple_ab->size()" << endl;
    //  exit(EXIT_FAILURE);
    //}

    //+===== FILL HIST ==============================================================================

    h_dalit->Fill(dali_t->at(0));
    for(int i = 0; i < dali_multi; i++) {
      h_dalit_all->Fill(dali_t->at(i));
      hIDtheta->Fill(dali_id->at(i), dali_cos->at(i));
    }

    hMINOSZ->Fill(MINOS_Z_cor);
    hbeta_vertex->Fill(beta_vertex);
    if(gamma_cos.size() > 0)
      hgamma_cos->Fill(gamma_cos.at(0));
    hNumTracks->Fill(MINOS_NumberTracks);

    //for(unsigned int j = 0; j < dali_edop_simple->size(); j++) {
    //  if(br51k && csa50ar->IsInside(aoqSA, zetSA)) {
    //    h_minoseff_50ar[0]->Fill(dali_edop_simple->at(j));
    //    if(MINOS_NumberTracks > 0)
    //      h_minoseff_50ar[1]->Fill(dali_edop_simple->at(j));
    //    if(MINOS_NumberTracks == 1 || MINOS_NumberTracks == 2)
    //      h_minoseff_50ar[2]->Fill(dali_edop_simple->at(j));
    //  }
    //  if(br54ca && csa53ca_minos->IsInside(aoqSA, zetSA)) {
    //    h_minoseff_53ca[0]->Fill(dali_edop_simple->at(j));
    //    if(MINOS_NumberTracks > 0)
    //      h_minoseff_53ca[1]->Fill(dali_edop_simple->at(j));
    //    if(MINOS_NumberTracks == 1 || MINOS_NumberTracks == 2)
    //      h_minoseff_53ca[2]->Fill(dali_edop_simple->at(j));
    //  }
    //}

    if(br56ca->IsInside(aoqBR, zetBR)) {
      htransbr[0]->Fill(F5X);
      htransbr[1]->Fill(F5X);
      if(PIDgate[0])
        htransbrsa[0]->Fill(F5X);
      if(PIDgate[1])
        htransbrsa[1]->Fill(F5X);
    }
    if(br56sc->IsInside(aoqBR, zetBR)) {
      htransbr[2]->Fill(F5X);
      if(PIDgate[2])
        htransbrsa[2]->Fill(F5X);
    }
    if(br58sc->IsInside(aoqBR, zetBR)) {
      htransbr[3]->Fill(F5X);
      if(PIDgate[3])
        htransbrsa[3]->Fill(F5X);
    }
    if(br59sc->IsInside(aoqBR, zetBR)) {
      htransbr[4]->Fill(F5X);
      if(PIDgate[4])
        htransbrsa[4]->Fill(F5X);
    }

    for(int i = 0; i < gatenum; i++) {
      if(PIDgates[i]) {
        hbetaF7F13[i]->Fill(betaF7F13);
        hbetaTH[i]->Fill(betaTH);
        hzvertex[i]->Fill(MINOS_Z_cor);
        hbetavertex[i]->Fill(beta_vertex);

        for(unsigned int j = 0; j < dali_edop->size(); j++) {  //+ w/o Addback
          hdoptime[i]->Fill(dali_t->at(j), dali_edop->at(j));
          if(-5 < dali_t->at(j) && dali_t->at(j) < DALITimegateUp) {
            hdopmult[i]->Fill(dali_multi, dali_edop->at(j));
            hdop[i * 10 + 0]->Fill(dali_edop->at(j));
            hdopafter[i * 10 + 0]->Fill(dali_edopafter->at(j));
            hdopsimple[i * 10 + 0]->Fill(dali_edop_simple->at(j));
            if(dali_multi == 1) {
              hdop[i * 10 + 1]->Fill(dali_edop->at(j));
              hdopafter[i * 10 + 1]->Fill(dali_edopafter->at(j));
              hdopsimple[i * 10 + 1]->Fill(dali_edop_simple->at(j));
            }
            if(dali_multi == 2) {
              hdop[i * 10 + 2]->Fill(dali_edop->at(j));
              hdopafter[i * 10 + 2]->Fill(dali_edopafter->at(j));
              hdopsimple[i * 10 + 2]->Fill(dali_edop_simple->at(j));
            }
            if(dali_multi == 3) {
              hdop[i * 10 + 3]->Fill(dali_edop->at(j));
              hdopafter[i * 10 + 3]->Fill(dali_edopafter->at(j));
              hdopsimple[i * 10 + 3]->Fill(dali_edop_simple->at(j));
            }
            if(dali_multi < 4) {
              hdop[i * 10 + 4]->Fill(dali_edop->at(j));
              hdopafter[i * 10 + 4]->Fill(dali_edopafter->at(j));
              hdopsimple[i * 10 + 4]->Fill(dali_edop_simple->at(j));
            }
          }
        }

        for(unsigned int j = 0; j < dali_edop_ab->size(); j++) {  //+ w/ Addback
          if(-5 < dali_t_ab->at(j) && dali_t_ab->at(j) < DALITimegateUp) {
            hdopID[i]->Fill(dali_id_ab->at(j), dali_edop_ab->at(j));
            hdop[i * 10 + 5]->Fill(dali_edop_ab->at(j));
            hdopafter[i * 10 + 5]->Fill(dali_edopafter_ab->at(j));
            hdopsimple[i * 10 + 5]->Fill(dali_edop_simple_ab->at(j));
            if(j != 0)
              hdopggcoin[i]->Fill(dali_edop_ab->at(0), dali_edop_ab->at(j));
            if(dali_multi_ab > 2 && j == 0) {
              hdopggsym[i]->Fill(dali_edop_ab->at(0), dali_edop_ab->at(1));
              hdopggsym[i]->Fill(dali_edop_ab->at(1), dali_edop_ab->at(0));
            }
            if(dali_multi_ab == 1) {
              hdop[i * 10 + 6]->Fill(dali_edop_ab->at(j));
              hdopafter[i * 10 + 6]->Fill(dali_edopafter_ab->at(j));
              hdopsimple[i * 10 + 6]->Fill(dali_edop_simple_ab->at(j));
            }
            if(dali_multi_ab == 2) {
              hdop[i * 10 + 7]->Fill(dali_edop_ab->at(j));
              hdopafter[i * 10 + 7]->Fill(dali_edopafter_ab->at(j));
              hdopsimple[i * 10 + 7]->Fill(dali_edop_simple_ab->at(j));
            }
            if(dali_multi_ab == 3) {
              hdop[i * 10 + 8]->Fill(dali_edop_ab->at(j));
              hdopafter[i * 10 + 8]->Fill(dali_edopafter_ab->at(j));
              hdopsimple[i * 10 + 8]->Fill(dali_edop_simple_ab->at(j));
            }
            if(dali_multi_ab < 4) {
              hdop[i * 10 + 9]->Fill(dali_edop_ab->at(j));
              hdopafter[i * 10 + 9]->Fill(dali_edopafter_ab->at(j));
              hdopsimple[i * 10 + 9]->Fill(dali_edop_simple_ab->at(j));
            }
          }
        }
      }
    }

    if(PIDgates[1]) {
      for(unsigned int j = 0; j < dali_edop_ab->size(); j++) {  //+ w/ Addback
        if(-5 < dali_t_ab->at(j) && dali_t_ab->at(j) < DALITimegateUp) {
          if(dali_id_ab->at(j) > 31)
            hdopFW[0]->Fill(dali_edop_ab->at(j));
          if(dali_id_ab->at(j) > 45)
            hdopFW[1]->Fill(dali_edop_ab->at(j));
          if(dali_id_ab->at(j) > 67)
            hdopFW[2]->Fill(dali_edop_ab->at(j));
          if(dali_id_ab->at(j) > 85)
            hdopFW[3]->Fill(dali_edop_ab->at(j));
        }
      }
    }

  }  //-while loop

  std::clog << std::endl;

  //+===== Write to the output file =====
  outfile->cd();

  for(int i = 0; i < 50; i++) {
    hdop[i]->Write();
    hdopafter[i]->Write();
  }

  //for(int i = 0; i < 110; i++)
  //  hdopsimple[i]->Write();
  hMINOSZ->Write();
  hbeta_vertex->Write();
  hgamma_cos->Write();
  hIDtheta->Write();
  hNumTracks->Write();
  //for(int i = 0; i < 3; i++) {
  //  h_minoseff_50ar[i]->Write();
  //  h_minoseff_53ca[i]->Write();
  //}

  h_dalit->Write();
  h_dalit_all->Write();
  for(int i = 0; i < gatenum; i++) {
    hdoptime[i]->Write();
    hdopmult[i]->Write();
    hdopID[i]->Write();
    hdopggcoin[i]->Write();
    hdopggsym[i]->Write();
  }

  for(int i = 0; i < gatenum; i++) {
    hbetaF7F13[i]->Write();
    hbetaTH[i]->Write();
    hzvertex[i]->Write();
    hbetavertex[i]->Write();

    htransbr[i]->Write();
    htransbrsa[i]->Write();
  }
  for(int i = 0; i < 4; i++) {
    hdopFW[i]->Write();
  }

  outfile->Write();
  outfile->Close("R");

  //+===== Delete vectors =====

  delete dali_e;
  delete dali_t;
  delete dali_cos;
  delete dali_x;
  delete dali_y;
  delete dali_z;
  delete dali_id;

  delete dali_edop;
  delete dali_edopafter;
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
  delete dali_edopafter_ab;
  delete dali_edop_simple_ab;
  delete dali_edop_beta_ab;
  delete dali_edop_theta_ab;

  //delete gamma_cos;
  //delete gamma_cosafter;
  //delete gamma_cos_ab;
  //delete gamma_cosafter_ab;

  stop_timer_tk(FileNumber, AllEntry);

  return 0;
}  //-main()

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
