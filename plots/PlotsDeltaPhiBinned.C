#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TPad.h>
#include <TMath.h>
#include <TLegend.h>
#include <TText.h>
#include <TPaveText.h>

using namespace std;

void CorrectStat(TH1D *h)
{
   for(int i = 1; i <= h->GetNbinsX(); ++i) {
        double val = h->GetBinContent(i);
        double err = val*1./sqrt(val * 3e6);
        h->SetBinError(i,err);
   }

}

int Plots(char *Save, TCanvas *c1, TString tag){

    TH1::SetDefaultSumw2();

  TPaveText* pavetext1 = new TPaveText(0.14,0.84,0.2,0.82,"NDC");
  pavetext1->SetFillColor(0); // text is black on white
  pavetext1->SetFillStyle(4000);
  pavetext1->SetTextSize(0.045); 
  pavetext1->SetTextAlign(12);
  pavetext1->SetTextFont(42);
  TText *text = new TText();
  text = pavetext1->AddText("#bf{CMS}");

  TPaveText* pavetext3 = new TPaveText(0.120,0.79,0.28,0.77,"NDC");
  pavetext3->SetFillColor(0); // text is black on white
  pavetext3->SetFillStyle(4000);
  pavetext3->SetTextSize(0.045); 
  pavetext3->SetTextAlign(12);
  pavetext3->SetTextFont(42);
  TText *text3 = new TText();
  text3 = pavetext3->AddText(" #it{Simulation}");

  TPaveText* pavetext4 = new TPaveText(0.4,0.94,0.92,0.92,"NDC");
  pavetext4->SetFillColor(0); // text is black on white
  pavetext4->SetFillStyle(4000);
  pavetext4->SetTextSize(0.045);
  pavetext4->SetTextAlign(12);
  pavetext4->SetTextFont(42);
  TText *text4 = new TText();
  text4 = pavetext4->AddText("HL-LHC 3000 fb^{-1} (14 TeV)");


  TFile* f1=new TFile("../farm/merged/bb3.root");

  TDirectoryFile *s1=(TDirectoryFile*)f1->Get("bjetsNew");

  TH1D *h1,*h2,*h3, *h4, *h5;

  if (tag == "bb") {
      h1=(TH1D*)s1->Get("AK8bbDeltaPhi");
      h2=(TH1D*)s1->Get("AK8bbDeltaPhi1");
      h3=(TH1D*)s1->Get("AK8bbDeltaPhi2");
      h4=(TH1D*)s1->Get("AK8bbDeltaPhi3");
      h5=(TH1D*)s1->Get("AK8bbDeltaPhi4");
      h4->Add(h5);
  }
  else if(tag == "bj") {
      h1=(TH1D*)s1->Get("AK8bDeltaPhi");
      h2=(TH1D*)s1->Get("AK8bDeltaPhi1");
      h3=(TH1D*)s1->Get("AK8bDeltaPhi2");
      h4=(TH1D*)s1->Get("AK8bDeltaPhi3");
      h5=(TH1D*)s1->Get("AK8bDeltaPhi4");
      h4->Add(h5);
  }
  else {
      h1=(TH1D*)s1->Get("AK8DeltaPhi");
      h2=(TH1D*)s1->Get("AK8DeltaPhi1");
      h3=(TH1D*)s1->Get("AK8DeltaPhi2");
      h4=(TH1D*)s1->Get("AK8DeltaPhi3");
      h5=(TH1D*)s1->Get("AK8DeltaPhi4");
      h4->Add(h5);
  }

  cout << h1->Integral() << " "<< h2->Integral() << endl;

  CorrectStat(h1);
  CorrectStat(h2);
  CorrectStat(h3);
  CorrectStat(h4);

  h1->Scale(1.0, "width");
  h2->Scale(1.0, "width");
  h3->Scale(1.0, "width");
  h4->Scale(1.0, "width");

  h1->Scale(1./h1->Integral("width"));
  h2->Scale(1./h2->Integral("width"));
  h3->Scale(1./h3->Integral("width"));
  h4->Scale(1./h4->Integral("width"));
 
  h2->Scale(10);
  h3->Scale(100);
  h4->Scale(1000);



  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0.01, 0.01, 1., 1.);
  gStyle->SetOptTitle(0);

  //   pad1->SetBottomMargin(0.005); // Upper and lower plot are joined
  pad1->SetTickx();
  pad1->SetTicky();
  pad1->SetLogy(); //log
  //pad1->SetLogx();
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  h1->SetStats(0);          // No statistics on upper plot
  h1->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#it{#sigma}}{d#Delta#phi}  (rad^{-1})");
  h1->GetXaxis()->SetTitle("#Delta#phi (rad)");
  h1->Draw();               // Draw h1
  h2->Draw("same");
  h3->Draw("same");
  h4->Draw("same");
  pavetext1->Draw("same");
  pavetext3->Draw("same");
  pavetext4->Draw("same");

  TLegend* leg = new TLegend(0.35, 0.63, .86, .88);
  if(tag =="bb")
      leg->SetHeader("Di-jet production (b+b)","C");
  else if(tag=="bj")
      leg->SetHeader("Di-jet production (b+j)","C");
  else 
      leg->SetHeader("Di-jet production (j+j)","C");

  leg->AddEntry(h1, "400 GeV <  #it{p}_{T} (x10^{0})", "p");
  leg->AddEntry(h2, "400 <  #it{p}_{T} < 800 GeV (x10^{1})", "p");
  leg->AddEntry(h3, "800 <  #it{p}_{T} < 1200 GeV (x10^{2})", "p");
  leg->AddEntry(h4, "1200 GeV <  #it{p}_{T} (x10^{3})", "p");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("same");


  h1->SetLineColor(kBlue-4);
  h1->SetLineWidth(2);
  h1->SetMarkerStyle(8);
  h1->SetMarkerColor(kBlue-4);

  h2->SetLineColor(kGreen-4);
  h2->SetLineWidth(2);
  h2->SetMarkerStyle(22);
  h2->SetMarkerColor(kGreen-4);

  h3->SetLineColor(kOrange+7);
  h3->SetLineWidth(2);
  h3->SetMarkerStyle(23);
  h3->SetMarkerColor(kOrange+7);

  h4->SetLineColor(kRed-4);
  h4->SetLineWidth(2);
  h4->SetMarkerStyle(34);
  h4->SetMarkerColor(kRed-4);

  h1->GetYaxis()->SetTitleSize(20);
  h1->SetTitleFont(42);
  h1->GetYaxis()->SetTitleSize(0.03);
  h1->GetYaxis()->SetTitleOffset(1.6);
  h1->GetYaxis()->SetLabelSize(0.035);
  h1->GetXaxis()->SetLabelSize(0.035);
  h1->GetXaxis()->SetTitleSize(0.04);
  h1->GetXaxis()->SetTitleOffset(.9);
  h1->GetYaxis()->SetNdivisions(6,5,0);
  h1->GetXaxis()->SetTitleOffset(1.2);
   
  h1->GetYaxis()->SetRangeUser(0.00001, 1e9);

  c1->Print((char*)Save);
  return 0;

}

void PlotsDeltaPhiBinned(){

  TCanvas *c1=new TCanvas("c1","c1",800,800);
  Plots((char*)"DeltaPhiBinnedBB.pdf", c1, "bb");
  c1->Clear();
  Plots((char*)"DeltaPhiBinnedBJ.pdf", c1, "bj");
  c1->Clear();
  Plots((char*)"DeltaPhiBinnedJJ.pdf", c1, "jj");
  c1->Clear();


}
