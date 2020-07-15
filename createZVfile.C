#include "TFile.h"
#include "TH1F.h"

void createZVfile() {
  TFile *f = TFile::Open("hist_daliall_AB15_ABth40_wTimeGate_minus5plus3.root");

  TH1F *h56ca55ca = ((TH1F *)f->Get("h_zvertex_br56ca_sa55ca"))->Clone("h");
  TH1F *h56ca55k  = ((TH1F *)f->Get("h_zvertex_br56ca_sa55k"))->Clone("h");
  TH1F *h58sc57ca = ((TH1F *)f->Get("h_zvertex_br58sc_sa57ca"))->Clone("h");
  TH1F *h59sc57ca = ((TH1F *)f->Get("h_zvertex_br59sc_sa57ca"))->Clone("h");
  TH1F *h59ti57ca = ((TH1F *)f->Get("h_zvertex_br59ti_sa57ca"))->Clone("h");
  TH1F *h60ti57ca = ((TH1F *)f->Get("h_zvertex_br60ti_sa57ca"))->Clone("h");

  TFile *f_56ca55ca = new TFile("hz_56Ca_55Ca.root", "RECREATE");
  TFile *f_56ca55k  = new TFile("hz_56Ca_55K.root", "RECREATE");
  TFile *f_58sc57ca = new TFile("hz_58Sc_57Ca.root", "RECREATE");
  TFile *f_59sc57ca = new TFile("hz_59Sc_57Ca.root", "RECREATE");
  TFile *f_59ti57ca = new TFile("hz_59Ti_57Ca.root", "RECREATE");
  TFile *f_60ti57ca = new TFile("hz_60Ti_57Ca.root", "RECREATE");

  f_56ca55ca->cd();
  h56ca55ca->Write();
  f_56ca55ca->Write();
  f_56ca55ca->Close();

  f_56ca55k->cd();
  h56ca55k->Write();
  f_56ca55k->Write();
  f_56ca55k->Close();

  f_58sc57ca->cd();
  h58sc57ca->Write();
  f_58sc57ca->Write();
  f_58sc57ca->Close();

  f_59sc57ca->cd();
  h59sc57ca->Write();
  f_59sc57ca->Write();
  f_59sc57ca->Close();

  f_59ti57ca->cd();
  h59ti57ca->Write();
  f_59ti57ca->Write();
  f_59ti57ca->Close();

  f_60ti57ca->cd();
  h60ti57ca->Write();
  f_60ti57ca->Write();
  f_60ti57ca->Close();

  f->Close();
}