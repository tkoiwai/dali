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
#include "TInterpreter.h"

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

//tk defined header files
#include "/home/koiwai/analysis/include/beamdef.h"
#include "/home/koiwai/analysis/include/smridef.h"
#include "/home/koiwai/analysis/include/mwdcdef.h"
#include "/home/koiwai/analysis/include/minosdef.h"
#include "/home/koiwai/analysis/include/dalidef.h"
#include "/home/koiwai/analysis/include/mychannelsdef.h"

using namespace std;
using namespace TMath;

#define WriteOneEnvFile(file) file->Write(TString(TString(file->GetRcName()).ReplaceAll("/","_")).ReplaceAll(".","_")); // / toka . no youna moji wo _ ni replace surudake.
#define WriteAllEnvFiles WriteOneEnvFile(env_par); WriteOneEnvFile(env_geo); //WriteOneEnvFile(env_nebt0); WriteOneEnvFile(env_neut0);

//=====External Function defined at last=======================================
//void SortDaliHit(Int_t, Int_t,vector <Int_t> *,vector <Double_t> *, vector <Double_t> *, vector <Double_t> *, vector <Double_t> *);
void SortDaliHit(Int_t, Int_t,vector <Int_t> *,vector <Double_t> *, vector <Double_t> *, vector <Double_t> *,vector <Double_t> *, vector <Double_t> *, vector <Double_t> *);
Double_t DopplerCorrection(Double_t, Double_t, Double_t);
inline bool exists_test(const std::string&);
inline bool exists_test(const TString&);
//=====main Function==========================================================
int main(int argc, char *argv[]){

  time_t start, stop;
  time(&start);
  
  gInterpreter->GenerateDictionary("vector<TVector3>","TVector3.h");
  
  Int_t FileNumber = TString(argv[1]).Atoi();
  
  if(FileNumber==0){
    std::cerr <<  " You should provide either a runnumber" << endl;
  }
  
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
  
  //===== Load input files =====
  TString infname = Form("/home/koiwai/analysis/rootfiles/ana/mychannels/ana%04d.root",FileNumber);
  TFile   *infile = TFile::Open(infname);
  TTree   *intr   = (TTree*)infile->Get("tr");

  Get_Branch_beam(intr);
  Get_Branch_smri(intr);
  Get_Branch_mwdc(intr);
  Get_Branch_minos(intr);
  Get_Branch_dali(intr);
  Get_Branch_mych(intr);
  
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
 
  
  //=====ROOT file setting==========================================================
  TString ofname =  Form("/home/koiwai/analysis/rootfiles/ana/dali/ana_dali%04d.root",FileNumber);
  TFile *outfile = new TFile(ofname,"RECREATE");
  TTree *tr = new TTree("anatrD","anatrD");
  tr->SetAutoSave(1e5); //Autosave every 10,000 events filled.

  //=====Define variables===================================================
  Int_t EventNumber = 0;
  Int_t RunNumber = -1;
  Int_t IncrementNumber = 0;

  Int_t br56sc_C;
  Int_t br56ca_C;
  Int_t br54ca_C;
  Int_t br51k_C;
  Int_t sa55ca_C;
  Int_t sa55k_C;
  Int_t sa53ca_C;
  Int_t sa50ar_C;
  
  Double_t vertexZ_cor;
  Double_t beta_vertex, gamma_vertex;

  Double_t beta_vertex_simple, gamma_vertex_simple;

  Int_t    dali_multi_ab;
  vector<Double_t> *dali_e_ab     = new vector<Double_t>();
  vector<Int_t>    *dali_id_ab    = new vector<Int_t>();
  vector<Double_t> *dali_cos_ab   = new vector<Double_t>();
  vector<Double_t> *dali_t_ab     = new vector<Double_t>();
  vector<Double_t> *dali_x_ab     = new vector<Double_t>();
  vector<Double_t> *dali_y_ab     = new vector<Double_t>();
  vector<Double_t> *dali_z_ab     = new vector<Double_t>();

  vector<Double_t> *dali_edop_ab  = new vector<Double_t>();

  TVector3 vertex;
  TVector3 fdc1;
  TVector3 bdc;
  vector<TVector3> dali_pos;

  TVector3 beam;
  vector<TVector3> gamma;

  vector<Double_t> gamma_cos;

  vector<Double_t> *dali_edop_simple_ab  = new vector<Double_t>();
  TVector3 vertex_simple;
  TVector3 beam_simple;
  vector<TVector3> gamma_simple;

  tr->Branch("EventNumber",&EventNumber);
  tr->Branch("RunNumber",&RunNumber);
  tr->Branch("IncrementNumber",&IncrementNumber);

  tr->Branch("br56sc",&br56sc_C);
  tr->Branch("br56ca",&br56ca_C);
  tr->Branch("br54ca",&br54ca_C);
  tr->Branch("br51k",&br51k_C);
  tr->Branch("sa55ca",&sa55ca_C);
  tr->Branch("sa55k",&sa55k_C);
  tr->Branch("sa53ca",&sa53ca_C);
  tr->Branch("sa50ar",&sa50ar_C);
  
  tr->Branch("vertexZ_cor",&vertexZ_cor);
  tr->Branch("beta_vertex",&beta_vertex);
  tr->Branch("gamma_vertex",&gamma_vertex);

  tr->Branch("beta_vertex_simple",&beta_vertex_simple);
  tr->Branch("gamma_vertex_simple",&gamma_vertex_simple);

  tr->Branch("dali_multi_ab",&dali_multi_ab);
  tr->Branch("dali_e_ab",&dali_e_ab);
  tr->Branch("dali_id_ab",&dali_id_ab);
  tr->Branch("dali_cos_ab",&dali_cos_ab);
  tr->Branch("dali_t_ab",&dali_t_ab);
  tr->Branch("dali_x_ab",&dali_x_ab);
  tr->Branch("dali_y_ab",&dali_y_ab);
  tr->Branch("dali_z_ab",&dali_z_ab);

  tr->Branch("dali_edop_ab",&dali_edop_ab);
   
  tr->Branch("vertexZ",&vertexZ);

  tr->Branch("vertex",&vertex);
  tr->Branch("fdc1",&fdc1);
  tr->Branch("bdc",&bdc);
  tr->Branch("dali_pos",&dali_pos);

  tr->Branch("gamma",&gamma);
  tr->Branch("beam",&beam);

  tr->Branch("gamma_cos",&gamma_cos);

  tr->Branch("dali_edop_simple_ab",&dali_edop_simple_ab);
  tr->Branch("vertex_simple",&vertex_simple);
  tr->Branch("beam_simple",&beam_simple);
  tr->Branch("gamma_simple",&gamma_simple);

  Int_t nEntry = intr->GetEntries();

  for(Int_t iEntry=0;iEntry<nEntry;iEntry++){
    //for(Int_t iEntry=0;iEntry<5;iEntry++){

    intr->GetEntry(iEntry);
    IncrementNumber++;

    if(IncrementNumber%100 == 0){
      std::clog << IncrementNumber/100 << " * 100 events treated..." << "\r";
    }
 
    RunNumber = FileNumber;
    EventNumber = EventNumber_mych;
    //cout << "e num mych " << EventNumber_mych << endl;; 

    br56sc_C = br56sc;
    br56ca_C = br56ca;
    br54ca_C = br54ca;
    br51k_C  = br51k;
    sa55ca_C = sa55ca;
    sa55k_C  = sa55k;
    sa53ca_C = sa53ca;
    sa50ar_C = sa50ar;
    
    vertexZ_cor  = Sqrt(-1);
    beta_vertex  = Sqrt(-1);
    gamma_vertex = Sqrt(-1);

    beta_vertex_simple  = Sqrt(-1);
    gamma_vertex_simple = Sqrt(-1);

    dali_e_ab->clear();
    dali_id_ab->clear();
    dali_cos_ab->clear();
    dali_t_ab->clear();
    dali_x_ab->clear();
    dali_y_ab->clear();
    dali_z_ab->clear();
    dali_multi_ab = 0;

    dali_edop_ab->clear();

    vertex.SetXYZ(Sqrt(-1),Sqrt(-1),Sqrt(-1));
    fdc1.SetXYZ(Sqrt(-1),Sqrt(-1),Sqrt(-1));
    bdc.SetXYZ(Sqrt(-1),Sqrt(-1),Sqrt(-1));
    dali_pos.clear();

    beam.SetXYZ(Sqrt(-1),Sqrt(-1),Sqrt(-1));
    gamma.clear();

    gamma_cos.clear();

    dali_edop_simple_ab->clear();
    vertex_simple.SetXYZ(Sqrt(-1),Sqrt(-1),Sqrt(-1));
    beam_simple.SetXYZ(Sqrt(-1),Sqrt(-1),Sqrt(-1));
    gamma_simple.clear();

    
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

    //===== ADD BACK END =====
    
    vertexZ_cor = vertexZ + MINOSoffsetZ;

    beta_vertex  = betaF7F13 - (betaF7F13 - beta_minoshodo)*vertexZ_cor/150.0;
    gamma_vertex = 1/Sqrt(1-beta_vertex*beta_vertex);
    
    vertex.SetXYZ(vertexX,vertexY,vertexZ_cor - DALIoffset); //To match the centre of MINOS cell and DALI Z = 0.(DALIOffset)

    //vertex.SetZ(vertex.Z() + 75.); 
    
    fdc1.SetXYZ(FDC1_X,FDC1_Y,Dist_MINOSfrontFDC1);
    //beam = fdc1 - vertex;
    
    for(Int_t i=0;i<dali_multi_ab;i++){
      dali_pos.push_back(TVector3(10*dali_x_ab->at(i),10*dali_y_ab->at(i),10*dali_z_ab->at(i)));        
      gamma.push_back(dali_pos.at(i)-vertex);      
      gamma_cos.push_back((gamma.at(i)).CosTheta());
    }

    

    //===== DOPPLER CORRECTION =====
    
    if(dali_multi_ab>=1){
      for(Int_t i=0;i<dali_multi_ab;i++){
	Double_t dali_edop_tmp = Sqrt(-1);
	dali_edop_tmp = dali_e_ab->at(i)*gamma_vertex*(1-beta_vertex*gamma_cos.at(i));
	dali_edop_ab->push_back(dali_edop_tmp);	
      }
    }

    //===== DOPPLER CORRECTION END =====

    //===== SINPLE DOPPLER CORRECTIOMN (WITHOUT MINOS )=====

    //beta_vertex_simple  = 0.5*(betaF7F13 + beta_minoshodo);
    //gamma_vertex_simple = 1/Sqrt(1 - beta_vertex_simple*beta_vertex_simple);
    const Double_t beta_mid = 0.57;
    const Double_t gamma_mid = 1/Sqrt(1 - beta_mid*beta_mid);
    //bdc.SetXYZ(BDC_X,BDC_Y,Dist_MINOSfrontBDC);
    //beam_simple = fdc1 - bdc;
    //vertex_simple.SetXYZ((FDC1_X-BDC_X)*-1.*Dist_MINOSfrontBDC/Dist_BDCFDC1,(FDC1_Y-BDC_Y)*-1.*Dist_MINOSfrontBDC/Dist_BDCFDC1,0.);
    //vertex_simple.SetXYZ(0,0,0);
    /*
    for(Int_t i=0;i<dali_multi_ab;i++){
       TVector3 gamma_simple_tmp;
       gamma_simple_tmp = dali_pos.at(i) - vertex_simple;
       gamma_simple.push_back(gamma_simple_tmp);
     }
    */
    if(dali_multi_ab>=1){
      for(Int_t i=0;i<dali_multi_ab;i++){
	//Double_t dali_edop_simple_tmp = Sqrt(-1);
	//dali_edop_simple_tmp = dali_e_ab->at(i)*gamma_vertex_simple*(1-beta_vertex_simple*dali_cos_ab->at(i));
	//dali_edop_simple_ab->push_back(dali_e_ab->at(i)*gamma_vertex_simple*(1-beta_vertex_simple*dali_cos_ab->at(i)));
	dali_edop_simple_ab->push_back(dali_e_ab->at(i)*gamma_mid*(1-beta_mid*dali_cos_ab->at(i)));
      }
    }

    //===== SIMPLE DOPPLER CORRECTION END =====    

   tr->Fill();
    
  }//while loop
  std::clog << std::endl;

  //tr->BuildIndex("RunNumber","EventNumber");
  outfile->cd();
  WriteAllEnvFiles;
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

  cout << "Conversion done." << endl;
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

inline bool exists_test (const std::string& name) {
  return ( access( name.c_str(), F_OK ) != -1 );
}
inline bool exists_test (const TString& name) {
  return ( access( name.Data(), F_OK ) != -1 );
}
