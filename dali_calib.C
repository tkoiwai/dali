#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TH2.h"
#include "TH1.h"
#include "TH1F.h"
#include "TLine.h"
#include "TFile.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TVirtualPad.h"
#include "TText.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom.h"
#include "TTree.h"

#define numberOfDetectors 226
 
FILE *fCalibrationResult;
double PeakPositions[2][numberOfDetectors] = {{0.0}};
double PeakErrors[2][numberOfDetectors] = {{0.0}};
TFile* inFile;
TFile* outFile;
TH2F* h2;
TTree* tr;
TCanvas *fCanvas = new TCanvas("Canvas", "Canvas", 900, 700);
TH1F* h[numberOfDetectors];
TH1F* h_new[numberOfDetectors];
int if_fit[numberOfDetectors];

int FindPeaks(int det_id, int num_peaks, float chmin, float chmax, float resolution = 1, double sigma = 2, double threshold = 0.3, int rebin = 1){
   int nFound;
   if(num_peaks!=1 && num_peaks!=2){
      cout<<"Wrong number of gamma peaks!"<<endl;
      return 0;
   }
   TSpectrum* s;
   TF1 *f;
   s = new TSpectrum(num_peaks+1);
   s->SetResolution(resolution);
   h_new[det_id] = (TH1F*)h[det_id]->Clone();
   h_new[det_id]->Rebin(rebin);
   h_new[det_id]->GetXaxis()->SetRangeUser(chmin,chmax);

   nFound = s->Search(h_new[det_id],sigma,"new",threshold);
   printf("In spectra [%i] I found %i candidate peaks to fit\n",det_id,nFound);
   if(nFound!=num_peaks){
      for(int i=0; i<num_peaks; i++){
         PeakPositions[i][det_id]=0;
         PeakErrors[i][det_id]=0;
      }
      cout<<"Strange spectrum"<<endl;
      if_fit[det_id]=1;
      return 0;
   }
   Float_t *xPeaks = s->GetPositionX(); 
   vector<float> Data, height;
   for(int i=0; i<nFound; i++){
      Data.push_back( *(xPeaks+i) );
   }
   sort( Data.begin(), Data.end() );
   Double_t minValue, maxValue;
   Double_t rangeMin, rangeMax;
   for(int i=0; i<nFound; i++){
      int counts = h_new[det_id]->GetBinContent(Data[0]);
      height.push_back(counts);
   }
   Double_t minValue, maxValue;
   Double_t rangeMin, rangeMax;
   minValue = 0.82*Data[0];
   maxValue = 1.12*Data[num_peaks-1];
   rangeMin = 0.80*Data[0];
   rangeMax = 1.22*Data[num_peaks-1];
   fCanvas->cd();
   h_new[det_id]->Draw();
   h_new[det_id]->GetXaxis()->SetRangeUser(rangeMin,rangeMax);
   if(num_peaks==1){
      f = new TF1("f", "gaus(0)+pol2(3)", rangeMin, rangeMax);
      f->SetParameter(0,height[0]);  f->SetParameter(1,Data[0]);  f->SetParameter(2,10);
      f->SetParLimits(1,0.97*Data[0],1.03*Data[0]);
      f->SetParLimits(2,0,15);
   }else{
      f = new TF1("f", "gaus(0)+gaus(3)+pol2(6)", rangeMin, rangeMax);
      f->SetParameter(0,height[0]);  f->SetParameter(1,Data[0]);  f->SetParameter(2,10);
      f->SetParLimits(1,0.97*Data[0],1.03*Data[0]);
      f->SetParLimits(2,0,15);
      f->SetParameter(3,height[1]);  f->SetParameter(4,Data[1]);  f->SetParameter(5,10);
      f->SetParLimits(4,0.97*Data[1],1.03*Data[1]);
      f->SetParLimits(5,0,15);
   }
   h_new[det_id]->Fit("f","QR");
   if(num_peaks==1){
      PeakPositions[0][det_id] = f->GetParameter(1);
      PeakErrors[0][det_id] = f->GetParameter(2);
   }else{
      PeakPositions[0][det_id] = f->GetParameter(1);
      PeakPositions[1][det_id] = f->GetParameter(4);
      PeakErrors[0][det_id] = f->GetParameter(2);
      PeakErrors[1][det_id] = f->GetParameter(5);
   }
   fCanvas->Update();
   if_fit[det_id]=1;
}

int OpenRoot(int calibrationDate, int source){
   char dummyText[1000];
   char rootFile[200];
   ifstream out_calib;
   int id;
   float peak1, peak2;
   float perr1, perr2;
   cout<<"Opening root files."<<endl;
   if(source==1){
      sprintf(rootFile,"./rootfiles/dali_calibrun/137cs00%i.root",calibrationDate);
      sprintf(dummyText,"./PeakPositions_Cs.txt.%i",calibrationDate);
      //for a temp workaround, use first peak of 88y instead of 137cs
      //sprintf(rootFile,"./rootfiles/88y%i.root",calibrationDate);
   }else if(source==2){
      sprintf(rootFile,"./rootfiles/dali_calibrun/60co00%i.root",calibrationDate);
      sprintf(dummyText,"./PeakPositions_Co.txt.%i",calibrationDate);
   }else if(source==3){
      sprintf(rootFile,"./rootfiles/dali_calibrun/88y00%i.root",calibrationDate);
      sprintf(dummyText,"./PeakPositions_Y.txt.%i",calibrationDate);
   }else{
      cout<<"source 1: 137Cs, 2: 60Co, 3: 88Y"<<endl;
      return 0;
   }

   inFile = new TFile(rootFile,"READ");
   //inFile->GetObject("h101",h2);
   inFile->GetObject("tr",tr);
   tr->Draw("DALI_Energy>>h2(300,0,3000)");
   int yBins= h2->GetYaxis()->GetNbins();
   int yMin = h2->GetYaxis()->GetXmin();
   int yMax = h2->GetYaxis()->GetXmax();
   for(int j=0;j<numberOfDetectors;j++)  {
      if_fit[j]=0;
      sprintf(dummyText,"h[%i]",j);
      for(int i=0; i<2; i++){
         PeakPositions[i][j]=0;
         PeakErrors[i][j]=0;
      }
      h[j] = new TH1F(dummyText,"",yBins,yMin,yMax);
      for(int k=0;k<yBins;k++) {
         int bin = h2->GetBinContent(j+1,k);
         h[j]->SetBinContent(k,bin);
      }
   }
   if(source==1){
      sprintf(dummyText,"./PeakPositions_Cs.txt.%i",calibrationDate);
   }else if(source==2){
      sprintf(dummyText,"./PeakPositions_Co.txt.%i",calibrationDate);
   }else if(source==3){
      sprintf(dummyText,"./PeakPositions_Y.txt.%i",calibrationDate);
   }else{
      cout<<"source 1: 137Cs, 2: 60Co, 3: 88Y"<<endl;
      return 0;
   }
   out_calib.open(dummyText, ios::in);
   while(out_calib.good()){
      if(source==1 || source==3){
         out_calib >> id>>peak1>>perr1;
         PeakPositions[0][id] = peak1;
         PeakErrors[0][id] = perr1;
      }else{
         out_calib >> id>>peak1>>perr1>>peak2>>perr2;
         PeakPositions[0][id] = peak1;
         PeakErrors[0][id] = perr1;
         PeakPositions[1][id] = peak2;
         PeakErrors[1][id] = perr2;
      }
   }
   return 1;
}

int SavePeaks(int calibrationDate, int source){
   char rootFile[200];
   char dummyText[1000];
   if(source==1){
      sprintf(dummyText,"./PeakPositions_Cs.txt.%i",calibrationDate);
      sprintf(rootFile,"./rootfiles/FitResult_137cs%i.root",calibrationDate);
   }else if(source==2){
      sprintf(dummyText,"./PeakPositions_Co.txt.%i",calibrationDate);
      sprintf(rootFile,"./rootfiles/FitResult_60co%i.root",calibrationDate);
   }else if(source==3){
      sprintf(dummyText,"./PeakPositions_Y.txt.%i",calibrationDate);
      sprintf(rootFile,"./rootfiles/FitResult_88y%i.root",calibrationDate);
   }else{
      cout<<"source 1: 137Cs, 2: 60Co, 3: 88Y"<<endl;
      return 0;
   }
   double peak1, peak2, error1, error2;

   fCalibrationResult= fopen(dummyText,"w");
   for(int j=0;j<numberOfDetectors;j++){
      if(source==1||source==3){
         peak1 = PeakPositions[0][j];
         error1 = PeakErrors[0][j];
         fprintf(fCalibrationResult,"%i %10.3f %10.3f\n",j,peak1,error1);
      }else{
         peak1 = PeakPositions[0][j];
         error1 = PeakErrors[0][j];
         peak2 = PeakPositions[1][j];
         error2 = PeakErrors[1][j];
         fprintf(fCalibrationResult,"%i %10.3f %10.3f %10.3f %10.3f\n",j,peak1,error1,peak2,error2);
      }
   }
   fclose(fCalibrationResult);

   outFile = new TFile(rootFile,"UPDATE");
   for(int j=0;j<numberOfDetectors;j++)  {
      if(if_fit[j]==1){
         h_new[j]->Write();
      }
   }
   outFile->Write();
   outFile->Close();
   return 1;
}

void FitPeaks(int calibrationDate){
   double y[4]={661.66, 1173.24, 1332.50, 1836.06}, x[numberOfDetectors][4];
   double yerr[4] = {0.0}, xerr[numberOfDetectors][4];
   ifstream cs_calib, co_calib, y_calib;
   char dummyText[1000];
   int id;
   float peak1, peak2, peak3, peak4;
   float perr1, perr2, perr3, perr4;
   sprintf(dummyText,"./PeakPositions_Cs.txt.%i",calibrationDate);
   cs_calib.open(dummyText, ios::in);
   sprintf(dummyText,"./PeakPositions_Co.txt.%i",calibrationDate);
   co_calib.open(dummyText, ios::in);
   sprintf(dummyText,"./PeakPositions_Y.txt.%i",calibrationDate);
   y_calib.open(dummyText, ios::in);
   double offset[numberOfDetectors];
   double gain[numberOfDetectors];
   double chi2[numberOfDetectors];
   while(cs_calib.good()){
      cs_calib >> id>>peak1>>perr1;
      x[id][0] = peak1;
      xerr[id][0] = 3*perr1;
   }
   while(co_calib.good()){
      co_calib >> id>>peak2>>perr2>>peak3>>perr3;
      x[id][1] = peak2;
      x[id][2] = peak3;
      xerr[id][1] = 3*perr2;
      xerr[id][2] = 3*perr3;
   }
   while(y_calib.good()){
      y_calib >> id>>peak4>>perr4;
      x[id][3] = peak4;
      xerr[id][3] = 3*perr4;
   }
   TGraphErrors* gr[numberOfDetectors];
   TCanvas *c1=new TCanvas("c1","c1",500,800);

   sprintf(dummyText,"./CalibrationResults.txt.%i",calibrationDate);
   fCalibrationResult= fopen(dummyText,"w");
   sprintf(dummyText,"./CalibrationResults.%i.root",calibrationDate);
   TFile fRoot(dummyText,"RECREATE");

   for(int m=0; m<numberOfDetectors; m++){
      if(x[m][0]>0 && x[m][1] && x[m][2]>0 && x[m][3]>0){
         //c1->cd();
         gr[m] = new TGraphErrors(4,x[m],y,xerr[m],yerr);
         f = new TF1("f","pol1", 0,2000);
         f->SetParameters(0,3);
         gr[m]->Fit("f","QR");
         f = gr[m]->GetFunction("f");
         f->SetLineColor(kBlue);
         f->SetLineWidth(1);
         offset[m] = f->GetParameter(0);
         gain[m] = f->GetParameter(1);
         chi2[m] = f->GetChisquare();
         cout<<"Good, id = "<<m<<", k = "<<gain[m]<<", b = "<<offset[m]<<", chi2 = "<<chi2[m]/2<<endl;
         gr[m]->Draw("AP");
         c1->Update();
         //c1->WaitPrimitive();
         fprintf(fCalibrationResult,"%i %10.3f %10.3f %10.3f\n",m,gain[m],offset[m],chi2[m]/2);
         gr[m]->Write();
      }else if(x[m][0]>0 && x[m][1] && x[m][2]>0 && x[m][3]==0){
         //c1->cd();
         gr[m] = new TGraphErrors(3,x[m],y,xerr[m],yerr);
         f = new TF1("f","pol1", 0,2000);
         f->SetParameters(0,3);
         gr[m]->Fit("f","QR");
         f = gr[m]->GetFunction("f");
         f->SetLineColor(kBlue);
         f->SetLineWidth(1);
         offset[m] = f->GetParameter(0);
         gain[m] = f->GetParameter(1);
         chi2[m] = f->GetChisquare();
         cout<<"Good, id = "<<m<<", k = "<<gain[m]<<", b = "<<offset[m]<<", chi2 = "<<chi2[m]/2<<endl;
         gr[m]->Draw("AP");
         c1->Update();
         //c1->WaitPrimitive();
         fprintf(fCalibrationResult,"%i %10.3f %10.3f %10.3f\n",m,gain[m],offset[m],chi2[m]/2);
         gr[m]->Write();
      }else if(x[m][0]>0 && x[m][1] && x[m][2]>0 && x[m][4]==0 && x[m][5]==0){
         //c1->cd();
         gr[m] = new TGraphErrors(3,x[m],y,xerr[m],yerr);
         f = new TF1("f","pol1", 0,2000);
         f->SetParameters(0,3);
         gr[m]->Fit("f","QR");
         f = gr[m]->GetFunction("f");
         f->SetLineColor(kBlue);
         f->SetLineWidth(1);
         offset[m] = f->GetParameter(0);
         gain[m] = f->GetParameter(1);
         chi2[m] = f->GetChisquare();
         cout<<"Good, id = "<<m<<", k = "<<gain[m]<<", b = "<<offset[m]<<", chi2 = "<<chi2[m]/2<<endl;
         gr[m]->Draw("AP");
         c1->Update();
         //c1->WaitPrimitive();
         fprintf(fCalibrationResult,"%i %10.3f %10.3f %10.3f\n",m,gain[m],offset[m],chi2[m]/2);
         gr[m]->Write();
      }else{
         gr[m] = new TGraphErrors(4,x[m],y,xerr[m],yerr);
         gr[m]->Draw("A*");
         gr[m]->Write();
         fprintf(fCalibrationResult,"%i %10.3f %10.3f %10.3f\n",m,0,0,10000);
      }
   }
   fRoot.Write();
   fRoot.Close();
}

int FindPeaks_All_crystals(int calibrationDate, int source, float chmin, float chmax, float resolution = 1, double sigma = 2, double threshold = 0.3, int rebin = 1){
   int if_good=OpenRoot(calibrationDate,source);
   if(if_good!=1){
      cout<<"Error"<<endl;
      return 0;
   }

   //start calibration work
   fCanvas->SetFillColor(0);
   fCanvas->SetFrameFillColor(0);
   fCanvas->SetBorderSize(0.);
   fCanvas->SetBorderMode(0);
   fCanvas->SetFrameBorderMode(0);
   fCanvas->SetTopMargin(0.05);
   fCanvas->SetBottomMargin(0.05);
   fCanvas->SetLeftMargin(0.05);
   fCanvas->SetRightMargin(0.05);
   int num_peaks;
   if(source==1||source==3){
      num_peaks=1;
   }else{
      num_peaks=2;
   }

   for(int m=0;m<numberOfDetectors;m++){
      if(PeakPositions[0][m]==0){
         FindPeaks(m,num_peaks,chmin,chmax,resolution,sigma,threshold,rebin);
      }
   }
   SavePeaks(calibrationDate,source);
}
