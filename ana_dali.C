#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<bitset>

#include"TROOT.h"
#include"TDirectory.h"
#include"TFile.h"
#include"TH1I.h"
#include"TTree.h"
#include"TCut.h"
#include"TCutG.h"
#include"TMath.h"
#include"TString.h"
#include"TEnv.h"
#include"TVector3.h"

using namespace std;
using namespace TMath;

//===== external func. =====
//void SortDaliHit(Int_t, Int_t,vector <Int_t> *,vector <Double_t> *, vector <Double_t> *,vector <Double_t> *);
void SortDaliHit(Int_t, Int_t,vector <Int_t> *,vector <Double_t> *, vector <Double_t> *, vector <Double_t> *,vector <Double_t> *, vector <Double_t> *, vector <Double_t> *);

//==== main func. =====
int main(int argc, char *argv[]){

  Int_t FileNum = TString(argv[1]).Atoi();
  
  //===== Load input Beam file =====
  TString infnameB = Form("/home/koiwai/analysis/rootfiles/ana/beam/ana_beam%04d.root",FileNum);
  TFile *infileB = TFile::Open(infnameB);
  
  TTree *anatrB;
  infileB->GetObject("anatrB",anatrB);
  
  //===== Declare Beam Variable =====
  Int_t EventNumber_beam, RunNumber_beam;

  Double_t betaF7F13;
  
  Double_t zetBR, aoqBR;

  Int_t BG_flag_beam;
  
  //===== Beam SetBranchAddress =====
  anatrB->SetBranchAddress("EventNumber",&EventNumber_beam);
  anatrB->SetBranchAddress("RunNumber",&RunNumber_beam);

  anatrB->SetBranchAddress("betaF7F13",&betaF7F13); // 
  
  anatrB->SetBranchAddress("zetBR",&zetBR);
  anatrB->SetBranchAddress("aoqBR",&aoqBR);

  anatrB->SetBranchAddress("BG_flag",&BG_flag_beam);

  //@@@@@ Beam end @@@@@

  //===== Load MINOS vertex =====
  TString infnameM = Form("/home/koiwai/analysis/rootfiles/minos/vertex/vertex%04d.root",FileNum);
  TFile  *infileM  = TFile::Open(infnameM);

  TTree *tr_vertex;
  infileM->GetObject("tr",tr_vertex);

  //=====
  Int_t EventNumber_MINOS, RunNumber_MINOS;

  Double_t vertexX, vertexY, vertexZ;
  Double_t vertexR, vertexTheta, vertexPhi;
  Double_t theta2p;

  
  //===== MINOS SetBranchAddress =====
  tr_vertex->SetBranchAddress("eventnum",&EventNumber_MINOS);
  tr_vertex->SetBranchAddress("runnum",&RunNumber_MINOS);
  tr_vertex->SetBranchAddress("vertexX",&vertexX);
  tr_vertex->SetBranchAddress("vertexY",&vertexY);
  tr_vertex->SetBranchAddress("vertexZ",&vertexZ);
  tr_vertex->SetBranchAddress("vertexR",&vertexR);
  tr_vertex->SetBranchAddress("vertexTheta",&vertexTheta);
  tr_vertex->SetBranchAddress("vertexPhi",&vertexPhi);
  tr_vertex->SetBranchAddress("theta2p",&theta2p);

  //@@@@@ MINOS vertex end @@@@@
  
  //===== Load DALI =====
  TString infnameD = Form("/home/koiwai/analysis/rootfiles/dali/cal_dali%04d.root",FileNum);
  TFile  *infileD  = TFile::Open(infnameD);

  TTree *tr_dali;
  infileD->GetObject("tr",tr_dali);

  //=====
  Long64_t EventNumber_DALI;
  Int_t    RunNumber_DALI;

  vector <Double_t> *DALI_Energy = 0;
  vector <Double_t> *DALI_Time = 0;
  vector <Double_t> *DALI_CosTheta = 0;
  vector <Double_t> *DALI_X = 0;
  vector <Double_t> *DALI_Y = 0;
  vector <Double_t> *DALI_Z = 0;
  vector <Double_t> *DALI_Layer = 0;
  vector <Int_t>    *DALI_ID = 0;
  

  Int_t DALI_Multi;

  
  //=====
  tr_dali->SetBranchAddress("EventNumber",&EventNumber_DALI);
  tr_dali->SetBranchAddress("RunNumber",&RunNumber_DALI);
  tr_dali->SetBranchAddress("DALI_Energy",&DALI_Energy);
  tr_dali->SetBranchAddress("DALI_CosTheta",&DALI_CosTheta);
  tr_dali->SetBranchAddress("DALI_Time",&DALI_Time);
  tr_dali->SetBranchAddress("DALI_X",&DALI_X);
  tr_dali->SetBranchAddress("DALI_Y",&DALI_Y);
  tr_dali->SetBranchAddress("DALI_Z",&DALI_Z);
  tr_dali->SetBranchAddress("DALI_Layer",&DALI_Layer);
  tr_dali->SetBranchAddress("DALI_ID",&DALI_ID);
  tr_dali->SetBranchAddress("DALI_Multi",&DALI_Multi);

  //===== AddFriend =====
  anatrB->AddFriend(tr_vertex);
  anatrB->AddFriend(tr_dali);

  //===== Load dat
  /*int crystal_id = -1;
  int ab_crystal_id[226][6] = {-1};
  
  ifstream ab_table;
  ab_table.open("/home/koiwai/analysis/AddbackTabl.out");
  if(!ab_table) cout << "Cannot open AddbackTable.out" << endl;

  while(ab_table.good()){
    ab_table >> crystal_id >> ab_crystal_id[]
    
    }*/ 

  int ab_crystal_id[226][7] = {
    { 10,  19,  -1,  -1,  -1,  -1,  -1},
    { 10,  11,  -1,  -1,  -1,  -1,  -1},
    { 11,  12,  -1,  -1,  -1,  -1,  -1},
    { 12,  13,  -1,  -1,  -1,  -1,  -1},
    { 13,  14,  -1,  -1,  -1,  -1,  -1},
    { 14,  15,  -1,  -1,  -1,  -1,  -1},
    { 15,  16,  -1,  -1,  -1,  -1,  -1},
    { 16,  17,  -1,  -1,  -1,  -1,  -1},
    { 17,  18,  -1,  -1,  -1,  -1,  -1},
    { 18,  19,  -1,  -1,  -1,  -1,  -1},
    {  0,   1,  20,  -1,  -1,  -1,  -1}, //10
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
    { 50,  69,  71,  90,  -1,  -1,  -1},
    { 51,  70,  72,  91,  -1,  -1,  -1},
    { 52,  71,  73,  92,  -1,  -1,  -1},
    { 53,  72,  93,  -1,  -1,  -1,  -1},
    { 54,  94,  -1,  -1,  -1,  -1,  -1},
    { 55,  95,  -1 , -1,  -1,  -1,  -1},
    { 56,  96,  -1,  -1,  -1,  -1,  -1},
    { 57,  97,  -1,  -1,  -1,  -1,  -1},
    { 58,  79,  98,  -1,  -1,  -1,  -1},
    { 59,  78,  80,  99,  -1,  -1,  -1},
    { 60,  79,  81, 100,  -1,  -1,  -1},
    { 61,  80,  82, 101,  -1,  -1,  -1},
    { 62,  81,  83, 102,  -1,  -1,  -1},
    { 63,  82, 103,  -1,  -1,  -1,  -1},
    { 64, 104,  -1,  -1,  -1,  -1,  -1},
    { 65, 105,  -1,  -1,  -1,  -1,  -1},
    { 66, 106, 107,  -1,  -1,  -1,  -1},
    { 67, 108, 109,  -1,  -1,  -1,  -1},
    { 68,  89, 110,  -1,  -1,  -1,  -1},
    { 69,  88,  90, 110, 111,  -1,  -1},
    { 70,  89,  91, 111, 112,  -1,  -1},
    { 71,  90,  92, 113, 114,  -1,  -1},
    { 72,  91,  93, 114, 115,  -1,  -1},
    { 73,  92, 115,  -1,  -1,  -1,  -1},
    { 74, 116, 117,  -1,  -1,  -1,  -1},
    { 75, 118, 119,  -1,  -1,  -1,  -1},
    { 76, 120, 121,  -1,  -1,  -1,  -1},
    { 77, 122, 123,  -1,  -1,  -1,  -1},
    { 78,  99, 124,  -1,  -1,  -1,  -1},
    { 79,  98, 100, 124, 125,  -1,  -1},
    { 80,  99, 101, 125, 126,  -1,  -1},
    { 81, 100, 102, 127, 128,  -1,  -1}, 
    { 82, 101, 103, 128, 129,  -1,  -1}, 
    { 83, 102, 129,  -1,  -1,  -1,  -1},
    { 84, 130, 131,  -1,  -1,  -1,  -1},
    { 85, 132, 133,  -1,  -1,  -1,  -1},
    { 86, 107, 133, 134,  -1,  -1,  -1}, 
    { 86, 106, 108, 135,  -1,  -1,  -1},
    { 87, 107, 109, 136,  -1,  -1,  -1},
    { 87, 108, 137,  -1,  -1,  -1,  -1},
    { 88,  89, 111, 138,  -1,  -1,  -1},
    { 89,  90, 110, 112, 139,  -1,  -1},
    { 90, 111, 113, 140,  -1,  -1,  -1},
    { 91, 112, 114, 141,  -1,  -1,  -1},
    { 91,  92, 113, 115, 142,  -1,  -1},
    { 92,  93, 114, 143,  -1,  -1,  -1},
    { 94, 117, 144,  -1,  -1,  -1,  -1},
    { 94, 116, 118, 145,  -1,  -1,  -1}, 
    { 95, 117, 119, 146,  -1,  -1,  -1},
    { 95, 118, 120, 147,  -1,  -1,  -1},
    { 96, 119, 121, 148,  -1,  -1,  -1},
    { 96, 120, 122, 149,  -1,  -1,  -1},
    { 97, 121, 123, 150,  -1,  -1,  -1},
    { 97, 122, 151,  -1,  -1,  -1,  -1},
    { 98,  99, 125, 152,  -1,  -1,  -1},
    { 99, 100, 124, 126, 153,  -1,  -1},
    {100, 125, 127, 154,  -1,  -1,  -1},
    {101, 126, 128, 155,  -1,  -1,  -1},
    {101, 102, 127, 129, 156,  -1,  -1},
    {102, 103, 128, 157,  -1,  -1,  -1},
    {104, 131, 158,  -1,  -1,  -1,  -1},
    {104, 130, 132, 159,  -1,  -1,  -1}, 
    {105, 131, 133, 160,  -1,  -1,  -1},
    {105, 106, 132, 161,  -1,  -1,  -1},
    {106, 135, 161,  -1,  -1,  -1,  -1},
    {107, 134, 136,  -1,  -1,  -1,  -1},
    {108, 135, 137,  -1,  -1,  -1,  -1},
    {109, 136, 138,  -1,  -1,  -1,  -1},
    {110, 137, 139,  -1,  -1,  -1,  -1},
    {111, 138, 140,  -1,  -1,  -1,  -1},
    {112, 139, 141,  -1,  -1,  -1,  -1},
    {113, 140, 142,  -1,  -1,  -1,  -1},
    {114, 141, 143,  -1,  -1,  -1,  -1},
    {115, 142, 144,  -1,  -1,  -1,  -1},
    {116, 143, 145,  -1,  -1,  -1,  -1},
    {117, 144, 146,  -1,  -1,  -1,  -1},
    {118, 145, 147,  -1,  -1,  -1,  -1},
    {119, 146, 148,  -1,  -1,  -1,  -1},
    {120, 147, 149,  -1,  -1,  -1,  -1},
    {121, 148, 150,  -1,  -1,  -1,  -1},
    {122, 149, 151,  -1,  -1,  -1,  -1},
    {123, 150, 152,  -1,  -1,  -1,  -1},
    {124, 151, 153,  -1,  -1,  -1,  -1},
    {125, 152, 154,  -1,  -1,  -1,  -1},
    {126, 153, 155,  -1,  -1,  -1,  -1},
    {127, 154, 156,  -1,  -1,  -1,  -1},
    {128, 155, 157,  -1,  -1,  -1,  -1},
    {129, 156, 158,  -1,  -1,  -1,  -1},
    {130, 157, 159,  -1,  -1,  -1,  -1},
    {131, 158, 160,  -1,  -1,  -1,  -1},
    {132, 159, 161,  -1,  -1,  -1,  -1},
    {133, 134, 160,  -1,  -1,  -1,  -1},
    {163, 171, 172, 191,  -1,  -1,  -1},
    {162, 172, 173, 192,  -1,  -1,  -1},
    {165, 176, 177, 195,  -1,  -1,  -1},
    {164, 177, 178, 196,  -1,  -1,  -1},
    {167, 181, 182, 199,  -1,  -1,  -1},
    {166, 182, 183, 200,  -1,  -1,  -1},
    {169, 186, 187, 203,  -1,  -1,  -1},
    {168, 187, 188, 204,  -1,  -1,  -1},
    {171, 189, 190,  -1,  -1,  -1,  -1},
    {162, 170, 172, 190, 191,  -1,  -1},
    {162, 163, 171, 173, 191, 192, 207}, 
    {163, 172, 174, 192, 193,  -1,  -1},
    {173, 175, 193,  -1,  -1,  -1,  -1},
    {174, 176, 194,  -1,  -1,  -1,  -1},
    {164, 175, 177, 194, 195,  -1,  -1},
    {164, 165, 176, 178, 195, 196, 210},
    {165, 177, 179, 196, 197,  -1,  -1},
    {178, 180, 197,  -1,  -1,  -1,  -1},
    {179, 181, 198,  -1,  -1,  -1,  -1},
    {166, 180, 182, 198, 199,  -1,  -1},
    {166, 167, 181, 183, 199, 200, 213},
    {167, 182, 184, 200, 201,  -1,  -1},
    {183, 185, 201,  -1,  -1,  -1,  -1},
    {184, 186, 202,  -1,  -1,  -1,  -1},
    {168, 185, 187, 202, 203,  -1,  -1},
    {168, 169, 186, 188, 203, 204, 216},
    {169, 187, 189, 204, 205,  -1,  -1},
    {170, 188, 205,  -1,  -1,  -1,  -1},
    {170, 171, 191, 205, 206,  -1,  -1},
    {162, 171, 172, 190, 192, 206, 207},  
    {163, 172, 173, 191, 193, 207, 208}, 
    {173, 174, 192, 194, 208,  -1,  -1},
    {175, 176, 193, 195, 209,  -1,  -1},
    {164, 176, 177, 194, 196, 209, 210}, 
    {165, 177, 178, 195, 197, 210, 211}, 
    {178, 179, 196, 198, 211,  -1,  -1},
    {180, 181, 197, 199, 212,  -1,  -1},
    {166, 181, 182, 198, 200, 212, 213}, 
    {167, 182, 183, 199, 201, 213, 214}, 
    {183, 184, 200, 202, 214,  -1,  -1},
    {185, 186, 201, 203, 215,  -1,  -1},
    {168, 186, 187, 202, 204, 215, 216}, 
    {169, 187, 188, 203, 205, 216, 217}, 
    {188, 189, 190, 204, 217,  -1,  -1},
    {190, 191, 207, 217, 218,  -1,  -1},
    {172, 191, 192, 206, 208, 218, 219}, 
    {192, 193, 207, 209, 219,  -1,  -1},
    {194, 195, 208, 210, 220,  -1,  -1},
    {177, 195, 196, 209, 211, 220, 221},
    {196, 197, 210, 212, 221,  -1,  -1},
    {198, 199, 211, 213, 222,  -1,  -1},
    {182, 199, 200, 212, 214, 222, 223},
    {200, 201, 213, 215, 223,  -1,  -1},
    {202, 203, 214, 216, 224,  -1,  -1},
    {187, 203, 204, 215, 217, 224, 225},
    {204, 205, 206, 216, 225,  -1,  -1},
    {206, 207, 219, 225,  -1,  -1,  -1},
    {207, 208, 218, 220,  -1,  -1,  -1},
    {209, 210, 219, 221,  -1,  -1,  -1},
    {210, 211, 220, 222,  -1,  -1,  -1},
    {212, 213, 221, 223,  -1,  -1,  -1},
    {213, 214, 222, 224,  -1,  -1,  -1},
    {215, 216, 223, 225,  -1,  -1,  -1},
    {216, 217, 218, 224,  -1,  -1,  -1}
  };
    
  
  //===== output
  TString ofname = Form("/home/koiwai/analysis/rootfiles/ana/dali/ana_dali%04d.root",FileNum);
  TFile  *ofile  = new TFile(ofname,"RECREATE");
  TTree  *anatrD = new TTree("anatrD","anatrD");

  //===== const, variables
  int runnum, eventnum;

  vector <Double_t> *dali_e = new vector <Double_t>();
  vector <Double_t> *dali_t = new vector <Double_t>();
  vector <Double_t> *dali_cos = new vector <Double_t>();
  vector <Double_t> *dali_x = new vector <Double_t>();
  vector <Double_t> *dali_y = new vector <Double_t>();
  vector <Double_t> *dali_z = new vector <Double_t>();
  vector <Double_t> *dali_layer = new vector <Double_t>();
  vector <Int_t> *dali_id = new vector <Int_t>();

  vector <Double_t> *dali_ab = new vector <Double_t>();
  vector <Double_t> *dali_ab_ecor = new vector <Double_t>();

  //vector <TVector3> *dali_pos = new vector <TVector3>();

  Int_t dali_multi = 0;
  
  //===== branch
  anatrD->Branch("runnum",&runnum);
  anatrD->Branch("eventnum",&eventnum);
  anatrD->Branch("dali_e",&dali_e);
  anatrD->Branch("dali_t",&dali_t);
  anatrD->Branch("dali_cos",&dali_cos);
  anatrD->Branch("dali_x",&dali_x);
  anatrD->Branch("dali_y",&dali_y);
  anatrD->Branch("dali_z",&dali_z);
  anatrD->Branch("dali_id",&dali_id);
  anatrD->Branch("dali_multi",&dali_multi);

  anatrD->Branch("dali_ab",&dali_ab);
  anatrD->Branch("dali_ab_ecor",&dali_ab_ecor);
  
  //anatrD->Branch("dali_pos",&dali_pos);

  //===== LOOP
  int nEntry = anatrB->GetEntries();
  for(int iEntry=0;iEntry<nEntry;iEntry++){
    //for(int iEntry=0;iEntry<5;iEntry++){
    
    if(iEntry%100==0&&iEntry!=0) clog << iEntry/1000 << "k events treated..." << "\r";
    
    anatrB->GetEntry(iEntry);
    
    runnum   = RunNumber_beam;
    eventnum = EventNumber_beam;
    
    //=== init
    int dali_multi = 0;
    for(unsigned int j=0;j<DALI_ID->size();j++){
      //koiwai + wimmer, removed negative energy hits
      if(DALI_Energy->at(j)<0)
	continue;
      if(DALI_ID->at(j)>225)
	continue;
      //end
      dali_id->push_back(DALI_ID->at(j));
      dali_e->push_back(DALI_Energy->at(j));
      dali_t->push_back(DALI_Time->at(j));
      dali_cos->push_back(DALI_CosTheta->at(j));
      dali_x->push_back(DALI_X->at(j));
      dali_y->push_back(DALI_Y->at(j));
      dali_z->push_back(DALI_Z->at(j));
      dali_multi++;
    }
    
    dali_ab->clear();
    dali_ab_ecor->clear();
    
    //=== calc
    
    if(dali_multi>1){

      SortDaliHit(0,dali_multi-1, dali_id, dali_e, dali_t, dali_x, dali_y, dali_z, dali_cos);
     

      double tmp_abenergy = 0;
    
      for(unsigned int i=1;i<dali_id->size();i++){
	for(int j=0;j<7;j++){
	  if(dali_id->at(i)==ab_crystal_id[dali_id->at(0)][j]){
	    //dali_ab->at(0) += dali_e->at(i);
	    tmp_abenergy += dali_e->at(i);
	  }
	  else continue;
	}
	dali_ab->push_back(tmp_abenergy);
      }
      
    }    
    
    
    

    



    anatrD->Fill();
  }//for
  ofile->cd();
  anatrD->Write();
  ofile->Close();
}//main()
  /*
void SortDaliHit(Int_t left, Int_t right,vector <Int_t> *dali_id,vector <Double_t> *dali_e,vector <Double_t> *dali_t, vector <Double_t> *dali_cos)
{
  Int_t TempID;
  Double_t TempEnergy;
  Double_t TempTime;
  Double_t TempCosTheta;

  int i = left, j = right;
  double pivot = dali_e->at((left + right) / 2);

  // partition 
  while (i <= j) {
    while (dali_e->at(i) > pivot)
      i++;
    while (dali_e->at(j) < pivot)
      j--;
    if (i <= j) {
      TempID = dali_id->at(j);
      TempEnergy = dali_e->at(j);
      TempTime = dali_t->at(j);
      TempCosTheta = dali_cos->at(j);

      dali_id->at(j) = dali_id->at(i);
      dali_e->at(j) = dali_e->at(i);
      dali_t->at(j) = dali_t->at(i);
      dali_cos->at(j) = dali_cos->at(i);

      dali_id->at(i) = TempID;
      dali_e->at(i) = TempEnergy;
      dali_t->at(i) = TempTime;
      dali_cos->at(i) = TempCosTheta;

      i++;
      j--;
    }
  };
  // recursion 
  if (left < j)
    SortDaliHit(left, j, dali_id, dali_e, dali_t, dali_cos);
  if (i < right)
    SortDaliHit(i, right, dali_id, dali_e, dali_t, dali_cos);
}
  */
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
