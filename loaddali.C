#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>

#include "TProof.h"
#include "TFile.h"
#include "TString.h"




void loaddali(Int_t runnum){
  
  if(runnum > 30){
    TString FileNameD = Form("/home/koiwai/analysis/rootfiles/dali/cal_dali%04d.root",runnum);
    TString FileNameM = Form("/home/koiwai/analysis/rootfiles/minos/cal/cal_minos%04d.root",runnum);
    TString FileNameV = Form("/home/koiwai/analysis/rootfiles/minos/vertex/vertex%04d.root",runnum);
    TString FileNameB = Form("/home/koiwai/analysis/rootfiles/ana/beam/ana_beam%04d.root",runnum);
    TString FileNameS = Form("/home/koiwai/analysis/rootfiles/ana/smri/ana_smri%04d.root",runnum);
    TString FileNameDC= Form("/home/koiwai/analysis/rootfiles/ana/mwdc/ana_mwdc%04d.root",runnum);
    TFile *f = TFile::Open(FileNameD);
    
  }else{
    TFile *f = TFile::Open("/home/koiwai/analysis/chain/calDch.root");
    //TString FileNameM = "/home/koiwai/analysis/chain/calMch.root";
    TString FileNameV = "/home/koiwai/analysis/chain/vertexch.root";
    TString FileNameB = "/home/koiwai/analysis/chain/anaBch.root";
    TString FileNameS = "/home/koiwai/analysis/chain/anaSch.root";
    //TString FileNameDC= "/home/koiwai/analysis/chain/anaDCch.root";
  }

  //tr->AddFriend("caltrM",FileNameM);
  
  tr->AddFriend("anatrB",FileNameB);
  tr->AddFriend("tr",FileNameV);
  tr->AddFriend("anatrS",FileNameS);
  //tr->AddFriend("anatrDC",FileNameDC);
  
  TFile *brcuts = new TFile("/home/koiwai/analysis/cutfiles/BRpid.root","");
  TCutG* brcut[50];
  int j=0;
  brcut[j] = (TCutG*)brcuts->Get("br64v"); brcut[j++]->SetName("br64v");
  brcut[j] = (TCutG*)brcuts->Get("br63v"); brcut[j++]->SetName("br63v");
  brcut[j] = (TCutG*)brcuts->Get("br62v"); brcut[j++]->SetName("br62v");
  brcut[j] = (TCutG*)brcuts->Get("br61v"); brcut[j++]->SetName("br61v");
  brcut[j] = (TCutG*)brcuts->Get("br60v"); brcut[j++]->SetName("br60v");

  //brcut[j] = (TCutG*)brcuts->Get("in62Ti");brcut[j++]->SetName("in62Ti");
  //brcut[j] = (TCutG*)brcuts->Get("br61ti");brcut[j++]->SetName("br61ti");
  brcut[j] = (TCutG*)brcuts->Get("br60ti");brcut[j++]->SetName("br60ti");
  brcut[j] = (TCutG*)brcuts->Get("br59ti");brcut[j++]->SetName("br59ti");
  brcut[j] = (TCutG*)brcuts->Get("br58ti");brcut[j++]->SetName("br58ti");
  brcut[j] = (TCutG*)brcuts->Get("br57ti");brcut[j++]->SetName("br57ti");

  brcut[j] = (TCutG*)brcuts->Get("br59sc");brcut[j++]->SetName("br59sc");
  brcut[j] = (TCutG*)brcuts->Get("br58sc");brcut[j++]->SetName("br58sc");
  brcut[j] = (TCutG*)brcuts->Get("br57sc");brcut[j++]->SetName("br57sc");
  brcut[j] = (TCutG*)brcuts->Get("br56sc");brcut[j++]->SetName("br56sc");
  brcut[j] = (TCutG*)brcuts->Get("br55sc");brcut[j++]->SetName("br55sc");
  brcut[j] = (TCutG*)brcuts->Get("br54sc");brcut[j++]->SetName("br54sc");

  //brcut[j] = (TCutG*)brcuts->Get("in57Ca");brcut[j++]->SetName("in57Ca");
  brcut[j] = (TCutG*)brcuts->Get("br56ca");brcut[j++]->SetName("br56ca");
  brcut[j] = (TCutG*)brcuts->Get("br55ca");brcut[j++]->SetName("br55ca");
  brcut[j] = (TCutG*)brcuts->Get("br54ca");brcut[j++]->SetName("br54ca");
  brcut[j] = (TCutG*)brcuts->Get("br53ca");brcut[j++]->SetName("br53ca");
  brcut[j] = (TCutG*)brcuts->Get("br52ca");brcut[j++]->SetName("br52ca");

  brcut[j] = (TCutG*)brcuts->Get("br53k"); brcut[j++]->SetName("br53k");
  brcut[j] = (TCutG*)brcuts->Get("br52k"); brcut[j++]->SetName("br52k");
  brcut[j] = (TCutG*)brcuts->Get("br51k"); brcut[j++]->SetName("br51k");
  brcut[j] = (TCutG*)brcuts->Get("br50k"); brcut[j++]->SetName("br50k");

  //brcut[j] = (TCutG*)brcuts->Get("in51Ar");brcut[j++]->SetName("in51Ar");
  brcut[j] = (TCutG*)brcuts->Get("br50ar");brcut[j++]->SetName("br50ar");
  brcut[j] = (TCutG*)brcuts->Get("br49ar");brcut[j++]->SetName("br49ar");

  //brcut[j] = (TCutG*)brcuts->Get("in49Cl");brcut[j++]->SetName("in49Cl");
  //brcut[j] = (TCutG*)brcuts->Get("in48Cl");brcut[j++]->SetName("in48Cl");
  
  cout << j << " incoming  BR cuts found " << endl;
  int nbr = j;
  
  TFile *sacuts = new TFile("/home/koiwai/analysis/cutfiles/SApid.root","");
  TCutG* sacut[50];
  j=0;
  sacut[j] = (TCutG*)sacuts->Get("sa63ti");sacut[j++]->SetName("sa63ti");
  sacut[j] = (TCutG*)sacuts->Get("sa62ti");sacut[j++]->SetName("sa62ti");
  sacut[j] = (TCutG*)sacuts->Get("sa61ti");sacut[j++]->SetName("sa61ti");
  sacut[j] = (TCutG*)sacuts->Get("sa60ti");sacut[j++]->SetName("sa60ti");
  sacut[j] = (TCutG*)sacuts->Get("sa59ti");sacut[j++]->SetName("sa59ti");
  sacut[j] = (TCutG*)sacuts->Get("sa58ti");sacut[j++]->SetName("sa58ti");
  sacut[j] = (TCutG*)sacuts->Get("sa57ti");sacut[j++]->SetName("sa57ti");
  sacut[j] = (TCutG*)sacuts->Get("sa56ti");sacut[j++]->SetName("sa56ti");
  sacut[j] = (TCutG*)sacuts->Get("sa55ti");sacut[j++]->SetName("sa55ti");
  sacut[j] = (TCutG*)sacuts->Get("sa54ti");sacut[j++]->SetName("sa54ti");
  sacut[j] = (TCutG*)sacuts->Get("sa53ti");sacut[j++]->SetName("sa53ti");

  //sacut[j] = (TCutG*)sacuts->Get("sa61sc");sacut[j++]->SetName("sa61sc");
  sacut[j] = (TCutG*)sacuts->Get("sa60sc");sacut[j++]->SetName("sa60sc");
  sacut[j] = (TCutG*)sacuts->Get("sa59sc");sacut[j++]->SetName("sa59sc");
  sacut[j] = (TCutG*)sacuts->Get("sa58sc");sacut[j++]->SetName("sa58sc");
  sacut[j] = (TCutG*)sacuts->Get("sa57sc");sacut[j++]->SetName("sa57sc");
  sacut[j] = (TCutG*)sacuts->Get("sa56sc");sacut[j++]->SetName("sa56sc");
  sacut[j] = (TCutG*)sacuts->Get("sa55sc");sacut[j++]->SetName("sa55sc");
  sacut[j] = (TCutG*)sacuts->Get("sa54sc");sacut[j++]->SetName("sa54sc");
  sacut[j] = (TCutG*)sacuts->Get("sa53sc");sacut[j++]->SetName("sa53sc");
  sacut[j] = (TCutG*)sacuts->Get("sa52sc");sacut[j++]->SetName("sa52sc");
  sacut[j] = (TCutG*)sacuts->Get("sa51sc");sacut[j++]->SetName("sa51sc");
  sacut[j] = (TCutG*)sacuts->Get("sa50sc");sacut[j++]->SetName("sa50sc");
  sacut[j] = (TCutG*)sacuts->Get("sa49sc");sacut[j++]->SetName("sa49sc");

  //sacut[j] = (TCutG*)sacuts->Get("sa58ca");sacut[j++]->SetName("sa58ca");
  sacut[j] = (TCutG*)sacuts->Get("sa57ca");sacut[j++]->SetName("sa57ca");
  sacut[j] = (TCutG*)sacuts->Get("sa56ca");sacut[j++]->SetName("sa56ca");
  sacut[j] = (TCutG*)sacuts->Get("sa55ca");sacut[j++]->SetName("sa55ca");
  sacut[j] = (TCutG*)sacuts->Get("sa54ca");sacut[j++]->SetName("sa54ca");
  sacut[j] = (TCutG*)sacuts->Get("sa53ca");sacut[j++]->SetName("sa53ca");
  sacut[j] = (TCutG*)sacuts->Get("sa52ca");sacut[j++]->SetName("sa52ca");
  sacut[j] = (TCutG*)sacuts->Get("sa51ca");sacut[j++]->SetName("sa51ca");
  sacut[j] = (TCutG*)sacuts->Get("sa50ca");sacut[j++]->SetName("sa50ca");
  sacut[j] = (TCutG*)sacuts->Get("sa49ca");sacut[j++]->SetName("sa49ca");
  sacut[j] = (TCutG*)sacuts->Get("sa48ca");sacut[j++]->SetName("sa48ca");
  sacut[j] = (TCutG*)sacuts->Get("sa47ca");sacut[j++]->SetName("sa47ca");
  sacut[j] = (TCutG*)sacuts->Get("sa46ca");sacut[j++]->SetName("sa46ca");

  sacut[j] = (TCutG*)sacuts->Get("sa55k"); sacut[j++]->SetName("sa55k"); 
  sacut[j] = (TCutG*)sacuts->Get("sa54k"); sacut[j++]->SetName("sa54k"); 
  sacut[j] = (TCutG*)sacuts->Get("sa53k"); sacut[j++]->SetName("sa53k"); 
  sacut[j] = (TCutG*)sacuts->Get("sa52k"); sacut[j++]->SetName("sa52k"); 
  sacut[j] = (TCutG*)sacuts->Get("sa51k"); sacut[j++]->SetName("sa51k"); 
  sacut[j] = (TCutG*)sacuts->Get("sa50k"); sacut[j++]->SetName("sa50k");
  sacut[j] = (TCutG*)sacuts->Get("sa49k"); sacut[j++]->SetName("sa49k"); 
  sacut[j] = (TCutG*)sacuts->Get("sa48k"); sacut[j++]->SetName("sa48k"); 
  sacut[j] = (TCutG*)sacuts->Get("sa47k"); sacut[j++]->SetName("sa47k"); 
  sacut[j] = (TCutG*)sacuts->Get("sa46k"); sacut[j++]->SetName("sa46k");
  sacut[j] = (TCutG*)sacuts->Get("sa45k"); sacut[j++]->SetName("sa45k"); 
  sacut[j] = (TCutG*)sacuts->Get("sa44k"); sacut[j++]->SetName("sa44k"); 
  sacut[j] = (TCutG*)sacuts->Get("sa43k"); sacut[j++]->SetName("sa43k"); 
 
  //sacut[j] = (TCutG*)sacuts->Get("sa53Ar");sacut[j++]->SetName("sa53Ar");
  //sacut[j] = (TCutG*)sacuts->Get("sa52Ar");sacut[j++]->SetName("sa52Ar");
  //sacut[j] = (TCutG*)sacuts->Get("sa51Ar");sacut[j++]->SetName("sa51Ar");
  //sacut[j] = (TCutG*)sacuts->Get("sa50Ar");sacut[j++]->SetName("sa50Ar");
  //sacut[j] = (TCutG*)sacuts->Get("sa49Ar");sacut[j++]->SetName("sa49Ar");
  //
  //sacut[j] = (TCutG*)sacuts->Get("sa50Cl");sacut[j++]->SetName("sa50Cl");
  //sacut[j] = (TCutG*)sacuts->Get("sa49Cl");sacut[j++]->SetName("sa49Cl");
  //sacut[j] = (TCutG*)sacuts->Get("sa48Cl");sacut[j++]->SetName("sa48Cl");
  //
  //sacut[j] = (TCutG*)sacuts->Get("sa48S"); sacut[j++]->SetName("sa48S"); 
  //sacut[j] = (TCutG*)sacuts->Get("sa47S"); sacut[j++]->SetName("sa47S"); 
  //sacut[j] = (TCutG*)sacuts->Get("sa46S"); sacut[j++]->SetName("sa46S"); 
 
  cout << j << " outgoming SA cuts found " << endl;
 
  int nsa = j;
  
  /*  
  gStyle->SetPalette(1);
  anatrS->Draw("zetSA:aoqSA>>h(1000,2.2,3,1000,18,24)","BG_flag==0&&BG_flag_beam==0&&hodo_id>6&&hodo_id<22","colz");
  for(int i=0;i<nsa;i++) sacut[i]->Draw("same");
  */
  f->cd();

  //TProof::Open("");
  //tr->SetProof();
  
}
