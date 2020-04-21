#define cal_dali_cxx
#include "/home/koiwai/analysis/include/dalidef.h"

bool signal_recieved = false;
void signalhandler(int sig){
  if(sig==SIGINT){
    signal_recieved = true;
  }
}

double get_time(){
  struct timeval t;
  gettimeofday(&t,NULL);
  double d = t.tv_sec + (double)t.tv_usec/1000000;
  return d;
}

//=====main Function==========================================================
int main(int argc, char *argv[]){

  double time_start = get_time();
  signal(SIGINT,signalhandler);

  time_t start, stop;
  time(&start);

  gInterpreter->GenerateDictionary("vector<TVector3>","TVector3.h");
  
  Long64_t MaxEventNumber = LLONG_MAX; //signed 8bit int. LLOMG_MAX ~ 2^63-1(const.)

  if (argc < 2 || argc > 4){
    printf("Usage: ./cal_dali RUNNUMBER\nOR     ./cal_dali RUNNUMBER MAXEVENTS\nOR     ./cal_dali RUNNUMBER MAXEVENTS TEST\n");
    exit(EXIT_FAILURE); 
  }
  printf("=======================================\n");
  if (argc > 2) {
    MaxEventNumber = TString(argv[2]).Atoi();
    printf(" You will process %lld events\n",MaxEventNumber);
  }
  
  Int_t FileNumber = TString(argv[1]).Atoi();
  TString RidfFileName = Form("/home/koiwai/analysis/ridf/sdaq02/run%04d.ridf.gz",FileNumber);

  if(!exists_test(RidfFileName)){
    cerr << " ERROR - '" << RidfFileName << "' does not exist '" << endl;
    exit(EXIT_FAILURE); 
  }

  printf("\n%s %d %s \n\n","=== Execute cal_dali for RUN",FileNumber,"===");

   
  //=====Load setting parameters=========================================
  TEnv *env_par = new TEnv("/home/koiwai/analysis/conversion_settings.prm");
  TEnv *env_geo = new TEnv(env_par->GetValue("geometrydata","")); //unit of length:mm
  Double_t Dist_F5F7           = env_geo->GetValue("Dist_F5F7",0.0);
  Double_t Dist_F7F13          = env_geo->GetValue("Dist_F7F13",0.0);
  Double_t BDC1_Width          = env_geo->GetValue("BDC1_Width",0.0); //mm
  Double_t FDC1_Width          = env_geo->GetValue("FDC1_Width",0.0); //mm
  Double_t Dist_BDC1BDC2       = env_geo->GetValue("Dist_BDC1BDC2",0.0); //mm
  Double_t Dist_BDC1Target     = env_geo->GetValue("Dist_BDC1Target",0.0); //mm
  Double_t Dist_BDC1FDC1       = env_geo->GetValue("Dist_BDC1FDC1",0.0); //Distance between the middle of BDC1 and the middle of FDC1 mm
  Double_t Dist_BDCFDC1        = env_geo->GetValue("Dist_BDCFDC1",0.0); //Distance between the middle of BDC and FDC1 mm
  Double_t Dist_SBTTarget      = env_geo->GetValue("Dist_SBTTarget",0.0); //mm
  Double_t Dist_MINOSfrontFDC1 = env_geo->GetValue("Dist_MINOSfrontFDC1",0.0); //mm
  Double_t Dist_MINOSfrontBDC  = env_geo->GetValue("Dist_MINOSfrontBDC",0.0); //mm, negative value
  
  Double_t MINOSoffsetZ = env_geo->GetValue("MINOSoffsetZ",0.0); //offset of vertexZ. [mm]
  Double_t DALIoffset   = env_geo->GetValue("DALIoffset",0.0);

  //===== ADDBACK TABLE ==================================================

  int AddBackTable[226][7] = {
    { 10,  19,  -1,  -1,  -1,  -1,  -1},//0
    { 10,  11,  -1,  -1,  -1,  -1,  -1},
    { 11,  12,  -1,  -1,  -1,  -1,  -1},
    { 12,  13,  -1,  -1,  -1,  -1,  -1},
    { 13,  14,  -1,  -1,  -1,  -1,  -1},
    { 14,  15,  -1,  -1,  -1,  -1,  -1},
    { 15,  16,  -1,  -1,  -1,  -1,  -1},
    { 16,  17,  -1,  -1,  -1,  -1,  -1},
    { 17,  18,  -1,  -1,  -1,  -1,  -1},
    { 18,  19,  -1,  -1,  -1,  -1,  -1},
    {  0,   1,  20,  -1,  -1,  -1,  -1},//10
    {  1,   2,  22,  -1,  -1,  -1,  -1},
    {  2,   3,  23,  -1,  -1,  -1,  -1},
    {  3,   4,  24,  -1,  -1,  -1,  -1},
    {  4,   5,  25,  -1,  -1,  -1,  -1},
    {  5,   6,  26,  -1,  -1,  -1,  -1},
    {  6,   7,  28,  -1,  -1,  -1,  -1},
    {  7,   8,  29,  -1,  -1,  -1,  -1},
    {  8,   9,  30,  -1,  -1,  -1,  -1},
    {  0,   9,  31,  -1,  -1,  -1,  -1}, 
    { 10,  32,  -1,  -1,  -1,  -1,  -1}, //20 
    { 33,  -1,  -1,  -1,  -1,  -1,  -1},
    { 11,  34,  -1,  -1,  -1,  -1,  -1},
    { 12,  35,  36,  -1,  -1,  -1,  -1}, 
    { 13,  37,  -1,  -1,  -1,  -1,  -1},
    { 14,  38,  -1,  -1,  -1,  -1,  -1},
    { 15,  39,  -1,  -1,  -1,  -1,  -1},
    { 40,  -1,  -1,  -1,  -1,  -1,  -1},
    { 16,  41,  -1,  -1,  -1,  -1,  -1},
    { 17,  42,  43,  -1,  -1,  -1,  -1}, 
    { 18,  44,  -1,  -1,  -1,  -1,  -1},//30 
    { 19,  45,  -1,  -1,  -1,  -1,  -1},
    { 20,  -1,  -1,  -1,  -1,  -1,  -1},
    { 21,  -1,  -1,  -1,  -1,  -1,  -1},
    { 22,  48,  -1,  -1,  -1,  -1,  -1},
    { 23,  49,  50,  -1,  -1,  -1,  -1}, 
    { 23,  51,  52,  -1,  -1,  -1,  -1}, 
    { 24,  53,  -1,  -1,  -1,  -1,  -1}, 
    { 25,  54,  -1,  -1,  -1,  -1,  -1},
    { 26,  -1,  -1,  -1,  -1,  -1,  -1},
    { 27,  -1,  -1,  -1,  -1,  -1,  -1},//40
    { 28,  58,  -1,  -1,  -1,  -1,  -1},
    { 29,  59,  60,  -1,  -1,  -1,  -1},
    { 29,  61,  62,  -1,  -1,  -1,  -1},
    { 30,  63,  -1,  -1,  -1,  -1,  -1},
    { 31,  64,  -1,  -1,  -1,  -1,  -1},
    { 66,  -1,  -1,  -1,  -1,  -1,  -1},
    { 67,  -1,  -1,  -1,  -1,  -1,  -1},
    { 34,  49,  68,  -1,  -1,  -1,  -1},
    { 35,  48,  50,  69,  -1,  -1,  -1}, 
    { 35,  49,  51,  70,  -1,  -1,  -1}, //50
    { 36,  50,  52,  71,  -1,  -1,  -1},
    { 36,  51,  53,  72,  -1,  -1,  -1},
    { 37,  52,  73,  -1,  -1,  -1,  -1},
    { 38,  74,  -1,  -1,  -1,  -1,  -1},
    { 75,  -1,  -1,  -1,  -1,  -1,  -1},
    { 76,  -1,  -1,  -1,  -1,  -1,  -1},
    { 77,  -1,  -1,  -1,  -1,  -1,  -1},
    { 41,  59,  78,  -1,  -1,  -1,  -1},
    { 42,  58,  60,  79,  -1,  -1,  -1},
    { 42,  59,  61,  80,  -1,  -1,  -1}, //60
    { 43,  60,  62,  81,  -1,  -1,  -1},
    { 43,  61,  63,  82,  -1,  -1,  -1},
    { 44,  62,  83,  -1,  -1,  -1,  -1},
    { 45,  84,  -1,  -1,  -1,  -1,  -1},
    { 85,  -1,  -1,  -1,  -1,  -1,  -1},
    { 46,  86,  -1,  -1,  -1,  -1,  -1},
    { 47,  87,  -1,  -1,  -1,  -1,  -1},
    { 48,  69,  88,  -1,  -1,  -1,  -1},
    { 49,  68,  70,  89,  -1,  -1,  -1}, 
    { 50,  69,  71,  90,  -1,  -1,  -1},//70
    { 51,  70,  72,  91,  -1,  -1,  -1},
    { 52,  71,  73,  92,  -1,  -1,  -1},
    { 53,  72,  93,  -1,  -1,  -1,  -1},
    { 54,  94,  -1,  -1,  -1,  -1,  -1},
    { 55,  95,  -1 , -1,  -1,  -1,  -1},
    { 56,  96,  -1,  -1,  -1,  -1,  -1},
    { 57,  97,  -1,  -1,  -1,  -1,  -1},
    { 58,  79,  98,  -1,  -1,  -1,  -1},
    { 59,  78,  80,  99,  -1,  -1,  -1},
    { 60,  79,  81, 100,  -1,  -1,  -1},//80
    { 61,  80,  82, 101,  -1,  -1,  -1},
    { 62,  81,  83, 102,  -1,  -1,  -1},
    { 63,  82, 103,  -1,  -1,  -1,  -1},
    { 64, 104,  -1,  -1,  -1,  -1,  -1},
    { 65, 105,  -1,  -1,  -1,  -1,  -1},
    { 66, 106, 107,  -1,  -1,  -1,  -1},
    { 67, 108, 109,  -1,  -1,  -1,  -1},
    { 68,  89, 110,  -1,  -1,  -1,  -1},
    { 69,  88,  90, 110, 111,  -1,  -1},
    { 70,  89,  91, 111, 112,  -1,  -1},//90
    { 71,  90,  92, 113, 114,  -1,  -1},
    { 72,  91,  93, 114, 115,  -1,  -1},
    { 73,  92, 115,  -1,  -1,  -1,  -1},
    { 74, 116, 117,  -1,  -1,  -1,  -1},
    { 75, 118, 119,  -1,  -1,  -1,  -1},
    { 76, 120, 121,  -1,  -1,  -1,  -1},
    { 77, 122, 123,  -1,  -1,  -1,  -1},
    { 78,  99, 124,  -1,  -1,  -1,  -1},
    { 79,  98, 100, 124, 125,  -1,  -1},
    { 80,  99, 101, 125, 126,  -1,  -1},//100
    { 81, 100, 102, 127, 128,  -1,  -1},
    { 82, 101, 103, 128, 129,  -1,  -1}, 
    { 83, 102, 129,  -1,  -1,  -1,  -1},
    { 84, 130, 131,  -1,  -1,  -1,  -1},
    { 85, 132, 133,  -1,  -1,  -1,  -1},
    { 86, 107, 133, 134,  -1,  -1,  -1}, 
    { 86, 106, 108, 135,  -1,  -1,  -1},
    { 87, 107, 109, 136,  -1,  -1,  -1},
    { 87, 108, 137,  -1,  -1,  -1,  -1},
    { 88,  89, 111, 138,  -1,  -1,  -1},//110
    { 89,  90, 110, 112, 139,  -1,  -1},
    { 90, 111, 113, 140,  -1,  -1,  -1},
    { 91, 112, 114, 141,  -1,  -1,  -1},
    { 91,  92, 113, 115, 142,  -1,  -1},
    { 92,  93, 114, 143,  -1,  -1,  -1},
    { 94, 117, 144,  -1,  -1,  -1,  -1},
    { 94, 116, 118, 145,  -1,  -1,  -1}, 
    { 95, 117, 119, 146,  -1,  -1,  -1},
    { 95, 118, 120, 147,  -1,  -1,  -1},
    { 96, 119, 121, 148,  -1,  -1,  -1},//120
    { 96, 120, 122, 149,  -1,  -1,  -1},
    { 97, 121, 123, 150,  -1,  -1,  -1},
    { 97, 122, 151,  -1,  -1,  -1,  -1},
    { 98,  99, 125, 152,  -1,  -1,  -1},
    { 99, 100, 124, 126, 153,  -1,  -1},
    {100, 125, 127, 154,  -1,  -1,  -1},
    {101, 126, 128, 155,  -1,  -1,  -1},
    {101, 102, 127, 129, 156,  -1,  -1},
    {102, 103, 128, 157,  -1,  -1,  -1},
    {104, 131, 158,  -1,  -1,  -1,  -1},//130
    {104, 130, 132, 159,  -1,  -1,  -1},
    {105, 131, 133, 160,  -1,  -1,  -1},
    {105, 106, 132, 161,  -1,  -1,  -1},
    {106, 135, 161,  -1,  -1,  -1,  -1},
    {107, 134, 136,  -1,  -1,  -1,  -1},
    {108, 135, 137,  -1,  -1,  -1,  -1},
    {109, 136, 138,  -1,  -1,  -1,  -1},
    {110, 137, 139,  -1,  -1,  -1,  -1},
    {111, 138, 140,  -1,  -1,  -1,  -1},
    {112, 139, 141,  -1,  -1,  -1,  -1},//140
    {113, 140, 142,  -1,  -1,  -1,  -1},
    {114, 141, 143,  -1,  -1,  -1,  -1},
    {115, 142, 144,  -1,  -1,  -1,  -1},
    {116, 143, 145,  -1,  -1,  -1,  -1},
    {117, 144, 146,  -1,  -1,  -1,  -1},
    {118, 145, 147,  -1,  -1,  -1,  -1},
    {119, 146, 148,  -1,  -1,  -1,  -1},
    {120, 147, 149,  -1,  -1,  -1,  -1},
    {121, 148, 150,  -1,  -1,  -1,  -1},
    {122, 149, 151,  -1,  -1,  -1,  -1},//150
    {123, 150, 152,  -1,  -1,  -1,  -1},
    {124, 151, 153,  -1,  -1,  -1,  -1},
    {125, 152, 154,  -1,  -1,  -1,  -1},
    {126, 153, 155,  -1,  -1,  -1,  -1},
    {127, 154, 156,  -1,  -1,  -1,  -1},
    {128, 155, 157,  -1,  -1,  -1,  -1},
    {129, 156, 158,  -1,  -1,  -1,  -1},
    {130, 157, 159,  -1,  -1,  -1,  -1},
    {131, 158, 160,  -1,  -1,  -1,  -1},
    {132, 159, 161,  -1,  -1,  -1,  -1},//160
    {133, 134, 160,  -1,  -1,  -1,  -1},
    {163, 171, 172, 191,  -1,  -1,  -1},
    {162, 172, 173, 192,  -1,  -1,  -1},
    {165, 176, 177, 195,  -1,  -1,  -1},
    {164, 177, 178, 196,  -1,  -1,  -1},
    {167, 181, 182, 199,  -1,  -1,  -1},
    {166, 182, 183, 200,  -1,  -1,  -1},
    {169, 186, 187, 203,  -1,  -1,  -1},
    {168, 187, 188, 204,  -1,  -1,  -1},
    {171, 189, 190,  -1,  -1,  -1,  -1},//170
    {162, 170, 172, 190, 191,  -1,  -1},
    {162, 163, 171, 173, 191, 192, 207}, 
    {163, 172, 174, 192, 193,  -1,  -1},
    {173, 175, 193,  -1,  -1,  -1,  -1},
    {174, 176, 194,  -1,  -1,  -1,  -1},
    {164, 175, 177, 194, 195,  -1,  -1},
    {164, 165, 176, 178, 195, 196, 210},
    {165, 177, 179, 196, 197,  -1,  -1},
    {178, 180, 197,  -1,  -1,  -1,  -1},
    {179, 181, 198,  -1,  -1,  -1,  -1},//180
    {166, 180, 182, 198, 199,  -1,  -1},
    {166, 167, 181, 183, 199, 200, 213},
    {167, 182, 184, 200, 201,  -1,  -1},
    {183, 185, 201,  -1,  -1,  -1,  -1},
    {184, 186, 202,  -1,  -1,  -1,  -1},
    {168, 185, 187, 202, 203,  -1,  -1},
    {168, 169, 186, 188, 203, 204, 216},
    {169, 187, 189, 204, 205,  -1,  -1},
    {170, 188, 205,  -1,  -1,  -1,  -1},
    {170, 171, 191, 205, 206,  -1,  -1},//190
    {162, 171, 172, 190, 192, 206, 207},
    {163, 172, 173, 191, 193, 207, 208}, 
    {173, 174, 192, 194, 208,  -1,  -1},
    {175, 176, 193, 195, 209,  -1,  -1},
    {164, 176, 177, 194, 196, 209, 210}, 
    {165, 177, 178, 195, 197, 210, 211}, 
    {178, 179, 196, 198, 211,  -1,  -1},
    {180, 181, 197, 199, 212,  -1,  -1},
    {166, 181, 182, 198, 200, 212, 213}, 
    {167, 182, 183, 199, 201, 213, 214},//200 
    {183, 184, 200, 202, 214,  -1,  -1},
    {185, 186, 201, 203, 215,  -1,  -1},
    {168, 186, 187, 202, 204, 215, 216}, 
    {169, 187, 188, 203, 205, 216, 217}, 
    {188, 189, 190, 204, 217,  -1,  -1},
    {190, 191, 207, 217, 218,  -1,  -1},
    {172, 191, 192, 206, 208, 218, 219}, 
    {192, 193, 207, 209, 219,  -1,  -1},
    {194, 195, 208, 210, 220,  -1,  -1},
    {177, 195, 196, 209, 211, 220, 221},//210
    {196, 197, 210, 212, 221,  -1,  -1},
    {198, 199, 211, 213, 222,  -1,  -1},
    {182, 199, 200, 212, 214, 222, 223},
    {200, 201, 213, 215, 223,  -1,  -1},
    {202, 203, 214, 216, 224,  -1,  -1},
    {187, 203, 204, 215, 217, 224, 225},
    {204, 205, 206, 216, 225,  -1,  -1},
    {206, 207, 219, 225,  -1,  -1,  -1},
    {207, 208, 218, 220,  -1,  -1,  -1},
    {209, 210, 219, 221,  -1,  -1,  -1},//220
    {210, 211, 220, 222,  -1,  -1,  -1},
    {212, 213, 221, 223,  -1,  -1,  -1},
    {213, 214, 222, 224,  -1,  -1,  -1},
    {215, 216, 223, 225,  -1,  -1,  -1},
    {216, 217, 218, 224,  -1,  -1,  -1}//225
  };

  
  //=====Load ANAROOT parameters===========================================
  TArtStoreManager *sman = TArtStoreManager::Instance();
  //TArtBigRIPSParameters *BigRIPSPara = TArtBigRIPSParameters::Instance();
  //TArtSAMURAIParameters *SAMURAIPara = TArtSAMURAIParameters::Instance();
  TArtDALIParameters *DALIPara = TArtDALIParameters::Instance();
  
  //TObjString *xmlfile;

  //TIter next(TString(env_par->GetValue("BigRIPSPara","")).Tokenize(" ")); //Make array of BigRIPSPara xmlfiles in conversion_settings.prm
  //while((xmlfile = (TObjString *)next())){
  //  BigRIPSPara->LoadParameter((char*)xmlfile->GetName()); //load the xmlfile(?)
  //}
  //
  //next = TString(env_par->GetValue("SamuraiPara","")).Tokenize(" ");
  //while((xmlfile = (TObjString*)next())){
  //  SAMURAIPara->LoadParameter((char*)xmlfile->GetName());
  //}
  //
  //next = TString(env_par->GetValue("DALIPara","")).Tokenize(" ");
  //while(( xmlfile = (TObjString *)next())) {
  //    DALIPara->LoadParameter((char*)xmlfile->GetName());
  //}

  if(FileNumber>=36&&FileNumber<=61)
    DALIPara->LoadParameter((char*)"/home/koiwai/analysis/db/DALI/DALI_20170505.xml");
  else if(FileNumber>=62&&FileNumber<=124)
    DALIPara->LoadParameter((char*)"/home/koiwai/analysis/db/DALI/DALI_20170508.xml");
  if(FileNumber>=125&&FileNumber<=196)
    DALIPara->LoadParameter((char*)"/home/koiwai/analysis/db/DALI/DALI_20170511.xml");
  if(FileNumber>=197&&FileNumber<=230)
    DALIPara->LoadParameter((char*)"/home/koiwai/analysis/db/DALI/DALI_20170514.xml");
  
  

  //Load detectors.
  TArtCalibDALI *CalibDALI = new TArtCalibDALI();
  
  //TArtRecoRIPS *RecoRIPS = new TArtRecoRIPS(); //added NOV14

  
  
  //=====Conversion ridf -> raw data======================================
  TArtEventStore *EventStore = new TArtEventStore;

  if(!EventStore->Open(RidfFileName)){
    std::cerr << "cannot open " << RidfFileName << std::endl << "aborting..." << std::endl;
    return 1;
  }else{
    cout << "input file: " << RidfFileName << endl;
  }

  //===== Load Beam file for PID cut =====
  TString infnameB = Form("/home/koiwai/analysis/rootfiles/ana/beam/ana_beam%04d.root",FileNumber);
  TFile   *infileB = TFile::Open(infnameB);
  TTree   *anatrB;
  infileB->GetObject("anatrB",anatrB);

  printf("%-20s %s \n","Input beam file:",infnameB.Data());
  
  anaD_Get_Branch_beam(anatrB);
  
  //===== Load smri file for PID cut =====
  TString infnameS = Form("/home/koiwai/analysis/rootfiles/ana/smri/ana_smri%04d.root",FileNumber);
  TFile   *infileS = TFile::Open(infnameS);
  TTree   *anatrS;
  infileS->GetObject("anatrS",anatrS);

  printf("%-20s %s \n","Input smri file:",infnameS.Data());

  anaD_Get_Branch_smri(anatrS);

  //===== Load MINOS file =====
  TString infnameM = Form("/home/koiwai/analysis/rootfiles/minos/vertex/vertex%04d.root",FileNumber);
  TFile   *infileM = TFile::Open(infnameM);
  TTree   *trM;
  infileM->GetObject("tr",trM);

  printf("%-20s %s \n\n","Input minos file:",infnameM.Data());
  
  anaD_Get_Branch_minos(trM);

  //===== AddFriend =====
  
  anatrB->AddFriend(anatrS);
  anatrB->AddFriend(trM);

  //===== Load cuts =======================================================================
  TFile *fcutSA_Ca = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_Ca.root");


  TCutG *csa57ca = (TCutG*)fcutSA_Ca->Get("sa57ca");
  TCutG *csa55ca = (TCutG*)fcutSA_Ca->Get("sa55ca");
  TCutG *csa54ca = (TCutG*)fcutSA_Ca->Get("sa54ca");
  TCutG *csa53ca = (TCutG*)fcutSA_Ca->Get("sa53ca");

  
 
  //=====ROOT file setting ==========================================================
  TString ofname;
  if(argc < 4)
    ofname = Form("/home/koiwai/analysis/rootfiles/ana/dali/cal_dali%04d.root",FileNumber);
  else if(argc == 4)
    ofname = Form("/home/koiwai/analysis/dali/testcal_dali%04d.root",FileNumber);

  TFile *outfile = new TFile(ofname,"RECREATE");
  TTree *tr = new TTree("tr","tr");
  //  tr->SetAutoSave(1e5); //Autosave every 10,000 events filled.

  printf("%-20s %s \n\n","Output dali file:",ofname.Data());

  Set_Branch_dali(tr);

  //===== LOOP start ================================================================
  int iEntry = 0;
  int AllEntry;
  if(argc > 2)
    AllEntry = MaxEventNumber;
  else
    AllEntry= anatrB->GetEntries();

  printf("\n Number of events to treat: %d\n\n",AllEntry);

  double time_prev = 0;
  double time_startloop = get_time();
  
  while(EventStore->GetNextEvent()&&EventNumber<MaxEventNumber){

    iEntry = EventNumber;
    anatrB->GetEntry(EventNumber);
    EventNumber++;

    if(iEntry%1000 == 0){	
      double time_end = get_time();	
      double t_diff_a = time_end - time_start;
      int t_hour_a    = t_diff_a/3600;
      int t_min_a     = (t_diff_a - 3600*t_hour_a)/60;
      double t_sec_a  = t_diff_a -3600*t_hour_a - 60*t_min_a;
      
      printf("%dh %dm %.2fs elapsed: %dk events done: %.2f events/s: about %.2fs to go: current speed: %.2f events/s \n",
	     t_hour_a,
	     t_min_a,
	     t_sec_a,
	     (int)(iEntry/1000),
	     iEntry/(time_end - time_startloop),
	     (AllEntry - iEntry)*(time_end - time_startloop)/(double)iEntry,
	     1000./(time_end - time_prev));	
      time_prev = get_time();
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
    
    vertexZ_cor  = Sqrt(-1);
    beta_vertex  = Sqrt(-1);
    gamma_vertex = Sqrt(-1);

    dali_edop->clear();

    vertex.SetXYZ(Sqrt(-1),Sqrt(-1),Sqrt(-1));
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
    
    if((br59sc==1)||(br58sc==1)||(br57sc==1)||(br56sc==1)||(br55sc==1)||(br54sc==1)||
       (br56ca==1)||(br55ca==1)||(br54ca==1)||(br53ca==1)||(br52ca==1)||
       (br54k ==1)||(br53k ==1)||(br52k ==1)||(br51k ==1)||(br50k ==1)||(br49k ==1))
      BRcut_bool = true;
    
    if(BRcut_bool==false){
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
      
      for(Int_t i=0; i<NumberOfNaIHit; i++) {

	TArtDALINaI *DALINaI = (TArtDALINaI *)DALINaIHits->At(i);
	ADC = DALINaI->GetRawADC();
	Time = DALINaI->GetTimeOffseted();

	if(ADC > 0 &&  ADC < 4095 && DALINaI->GetEnergy() > 0 /*&& Time > 0*/) {
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

    if(DALI_Multi>1){
      SortDaliHit(0,DALI_Multi-1,
		  DALI_ID,
		  DALI_Energy,
		  DALI_Time,
		  DALI_X,
		  DALI_Y,
		  DALI_Z,
		  DALI_CosTheta);
    }else if(DALI_Multi==0){
      tr->Fill();
      continue;
    }

    vertexZ_cor = vertexZ + MINOSoffsetZ;

    beta_vertex  = betaF7F13 - (betaF7F13 - betaTH)*vertexZ_cor/150.0;
    gamma_vertex = 1/Sqrt(1-beta_vertex*beta_vertex);
    
    vertex.SetXYZ(vertexX,vertexY,vertexZ_cor - DALIoffset);
    //To match the centre of MINOS cell and DALI Z = 0.(DALIOffset)

    for(Int_t i=0;i<DALI_Multi;i++){
      dali_pos.push_back(TVector3(10*DALI_X->at(i),10*DALI_Y->at(i),10*DALI_Z->at(i)));
      gamma_pos.push_back(dali_pos.at(i)-vertex);      
      gamma_cos.push_back((gamma_pos.at(i)).CosTheta());
    }
    
    //===== DOPPLER CORRECTION =====
    
    if(DALI_Multi>=1){
      for(Int_t i=0;i<DALI_Multi;i++){
	Double_t dali_edop_tmp = Sqrt(-1);
	dali_edop_tmp = DALI_Energy->at(i)*gamma_vertex*(1-beta_vertex*gamma_cos.at(i));
	dali_edop->push_back(dali_edop_tmp);	
      }
    }

    //===== DOPPLER CORRECTION END =====

    //===== SINPLE DOPPLER CORRECTIOMN (WITHOUT MINOS )=====

    const Double_t beta_mid = 0.57;
    const Double_t gamma_mid = 1/Sqrt(1 - beta_mid*beta_mid);
    if(DALI_Multi>=1){
      for(Int_t i=0;i<DALI_Multi;i++){
	dali_edop_simple->push_back(DALI_Energy->at(i)*gamma_mid*(1-beta_mid*DALI_CosTheta->at(i)));
      }
    }

    //===== SIMPLE DOPPLER CORRECTION END =====

    //===== ADDBACK ===========================
    if(DALI_Multi>1){
      dali_e_ab->push_back(DALI_Energy->at(0));
      dali_id_ab->push_back(DALI_ID->at(0));
      dali_cos_ab->push_back(DALI_CosTheta->at(0));
      dali_t_ab->push_back(DALI_Time->at(0));      
      dali_x_ab->push_back(DALI_X->at(0));
      dali_y_ab->push_back(DALI_Y->at(0));
      dali_z_ab->push_back(DALI_Z->at(0));


      //===== ADD BACK =====
      
      Bool_t AddBack_flag = false;
      
      for(Int_t i=1;i<DALI_Multi;i++){
	AddBack_flag = false;
	dali_multi_ab = dali_id_ab->size();
	for(Int_t k=0;k<dali_multi_ab;k++){
	  for(Int_t j=0;j<7;j++){	    
	    if(DALI_ID->at(i)==AddBackTable[dali_id_ab->at(k)][j]&&DALI_Energy->at(i)>300){
	      dali_e_ab->at(k) = dali_e_ab->at(k) + DALI_Energy->at(i);
	      AddBack_flag = true;
	      break;
	    }
	  }
	  if(AddBack_flag==true) break;
	}
	if(AddBack_flag==false){
	  dali_e_ab->push_back(DALI_Energy->at(i));
	  dali_id_ab->push_back(DALI_ID->at(i));
	  dali_cos_ab->push_back(DALI_CosTheta->at(i));
	  dali_t_ab->push_back(DALI_Time->at(i));
	  dali_x_ab->push_back(DALI_X->at(i));
	  dali_y_ab->push_back(DALI_Y->at(i));
	  dali_z_ab->push_back(DALI_Z->at(i));

	}
      }//for DALI_Multi
      
      
      dali_multi_ab = dali_id_ab->size();
      
    }//DALI_Multi>1
    else if(DALI_Multi==1){
      dali_e_ab->push_back(DALI_Energy->at(0));
      dali_id_ab->push_back(DALI_ID->at(0));
      dali_cos_ab->push_back(DALI_CosTheta->at(0));
      dali_t_ab->push_back(DALI_Time->at(0));
      dali_x_ab->push_back(DALI_X->at(0));
      dali_y_ab->push_back(DALI_Y->at(0));
      dali_z_ab->push_back(DALI_Z->at(0));
      
      dali_multi_ab = 1;
    }else if(DALI_Multi==0){
      tr->Fill();
      continue;
    }

    if(dali_multi_ab>=1){
      for(Int_t i=0;i<dali_multi_ab;i++){
	Double_t dali_edop_ab_tmp = Sqrt(-1);
	dali_edop_ab_tmp = dali_e_ab->at(i)*gamma_vertex*(1-beta_vertex*gamma_cos.at(i));
	dali_edop_ab->push_back(dali_edop_ab_tmp);	
      }
    } 

    //===== ADDBACK END =======================

    //===== TIMING GATE =======================
    //===== TIMING GATE END ===================

    //===== isotope gate =====
    if(csa57ca->IsInside(aoqSA,zetSA))  sa57ca = kTRUE;
    if(csa55ca->IsInside(aoqSA,zetSA))  sa55ca = kTRUE;
    if(csa54ca->IsInside(aoqSA,zetSA))  sa54ca = kTRUE;
    if(csa53ca->IsInside(aoqSA,zetSA))  sa53ca = kTRUE;

   tr->Fill();
    
  }//while loop
  std::clog << std::endl;

  tr->BuildIndex("RunNumber","EventNumber");
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

  time(&stop);

  int t_hour = (int)difftime(stop,start)/3600;
  int t_min  = (int)(difftime(stop,start) - t_hour*3600)/60;
  double t_sec  = difftime(stop,start) - t_hour*3600 - t_min*60;
  printf("Elapsed time: %dh %dm %.1f seconds\n",t_hour,t_min,t_sec);

  double time_end = get_time();
  printf("Average process speed: %f events/s\n",AllEntry/(time_end - time_start));
  printf("RUN%d: Conversion finished!: caldaliok%d\n",FileNumber,FileNumber);
  return 0;
}//main()










void SortDaliHit(Int_t left, Int_t right,vector <Int_t> *DALI_ID,vector <Double_t> *DALI_Energy, vector <Double_t> *DALI_Time, vector <Double_t> *DALI_X, vector <Double_t> *DALI_Y, vector <Double_t> *DALI_Z, vector <Double_t> *DALI_CosTheta)
{
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
  while (i <= j) {
    while (DALI_Energy->at(i) > pivot)
      i++;
    while (DALI_Energy->at(j) < pivot)
      j--;
    if (i <= j) {
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
  if (left < j)
    SortDaliHit(left, j, DALI_ID, DALI_Energy, DALI_Time, DALI_X, DALI_Y, DALI_Z, DALI_CosTheta);
  if (i < right)
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

inline bool exists_test (const std::string& name) {
  return ( access( name.c_str(), F_OK ) != -1 );
}
inline bool exists_test (const TString& name) {
  return ( access( name.Data(), F_OK ) != -1 );
}
