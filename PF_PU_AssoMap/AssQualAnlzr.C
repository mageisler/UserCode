#include "/home/home2/institut_3b/geisler/Phd-Study/CMSSW/Helpers/includes.hpp"
#include "/home/home2/institut_3b/geisler/Phd-Study/CMSSW/Helpers/littleHelpers.hpp"

int main(){

  using namespace std;

  /***************************************************************/
  /* read the Files & get the vector <TH1F>'s                    */
  /***************************************************************/

  TH1F* histo_RT_all;
  TH1F* histo_RT_1st;
  TH1F* histo_RT_all_corr;
  TH1F* histo_RT_1st_corr;

  TH2F* histo_GEN_Pos;
  TH2F* histo_GEN_Dist;
  TH2F* histo_FRT_Pos;
  TH2F* histo_FRT_Dist;

  TFile* File_ref = TFile::Open("/home/home2/institut_3b/geisler/Phd-Study/CMSSW/CMSSW_4_4_2_patch6/src/MGeisler/PF_PU_AssoMap/AssociationQualityAnalysis3.root");
  File_ref->cd();

  gDirectory->GetObject("h_RT_all",histo_RT_all);
  gDirectory->GetObject("h_RT_1st",histo_RT_1st);
  gDirectory->GetObject("h_RT_corr",histo_RT_all_corr);
  gDirectory->GetObject("h_RT_1st_corr",histo_RT_1st_corr);

  gDirectory->GetObject("h_GEN_Pos",histo_GEN_Pos);
  gDirectory->GetObject("h_GEN_Dist",histo_GEN_Dist);
  gDirectory->GetObject("h_FRT_Pos",histo_FRT_Pos);
  gDirectory->GetObject("h_FRT_Dist",histo_FRT_Dist);


  /***************************************************************/
  /* Draw the TH1F's*                                            */
  /***************************************************************/

  int my_colors[4] = {kBlack,kRed+2,kViolet+7,kGreen};

  TCanvas *canv_AQA = new TCanvas("canv_AQA","Guck mal",0,0,1440,900);
  gStyle->SetNumberContours(255);
  gStyle->SetOptStat("ne");
  canv_AQA->Divide(1,1);
   
  gPad->SetBorderMode(0);
  gPad->SetFillColor(0);
   
  canv_AQA->cd(1);
  gPad->SetGrid();  
  gPad->SetLogy();  
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.08);
  gPad->SetFillColor(0);

  histo_RT_all->SetMaximum(histo_RT_all->GetMaximum() + 1000.*sqrt(histo_RT_all->GetMaximum()));

  histo_RT_all->GetXaxis()->CenterTitle();
  histo_RT_all->GetXaxis()->SetTitleOffset(1.1);

  histo_RT_all->GetYaxis()->CenterTitle();
  histo_RT_all->GetYaxis()->SetTitleOffset(1.3);
  histo_RT_all->GetYaxis()->SetTitleSize(0.055);
  histo_RT_all->GetYaxis()->SetLabelSize(0.055);

  histo_RT_all->SetLineColor(my_colors[0]);
  histo_RT_all->SetLineWidth(2);
  histo_RT_all->Draw();
  moveStatBox2(canv_AQA,histo_RT_all,my_colors[0],0);

  histo_RT_1st->SetLineColor(my_colors[1]);
  histo_RT_1st->SetLineWidth(2);
  histo_RT_1st->Draw("sames");
  moveStatBox2(canv_AQA,histo_RT_1st,my_colors[1],1);

  histo_RT_all_corr->SetLineColor(my_colors[2]);
  histo_RT_all_corr->SetLineWidth(2);
  histo_RT_all_corr->Draw("sames");
  moveStatBox2(canv_AQA,histo_RT_all_corr,my_colors[2],2);

  histo_RT_1st_corr->SetLineColor(my_colors[3]);
  histo_RT_1st_corr->SetLineWidth(2);
  histo_RT_1st_corr->Draw("sames");
  moveStatBox2(canv_AQA,histo_RT_1st_corr,my_colors[3],3);
  
  canv_AQA->SaveAs("AssociationQualityAnalysis3.png");

  /*********************************************/

  TCanvas *canv_2d = new TCanvas("canv_2d","Guck mal",0,0,1440,900);
  gStyle->SetNumberContours(255);
  gStyle->SetOptStat("ne");
  canv_2d->Divide(2,2);
   
  gPad->SetBorderMode(0);
  gPad->SetFillColor(0);
   
  canv_2d->cd(1);
  gPad->SetGrid();  
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.08);
  gPad->SetFillColor(0);

  histo_GEN_Pos->GetXaxis()->CenterTitle();
  histo_GEN_Pos->GetYaxis()->CenterTitle();
  histo_GEN_Pos->GetZaxis()->CenterTitle();
  
  histo_GEN_Pos->GetXaxis()->SetTitleOffset(1.1);
  histo_GEN_Pos->GetYaxis()->SetTitleOffset(1.1);
  histo_GEN_Pos->SetMarkerColor(my_colors[0]);
  histo_GEN_Pos->SetStats(0);
  histo_GEN_Pos->Draw();
   
  canv_2d->cd(2);
  gPad->SetGrid();  
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.08);
  gPad->SetFillColor(0);

  histo_GEN_Dist->GetXaxis()->CenterTitle();
  histo_GEN_Dist->GetYaxis()->CenterTitle();
  histo_GEN_Dist->GetZaxis()->CenterTitle();
  
  histo_GEN_Dist->GetXaxis()->SetTitleOffset(1.1);
  histo_GEN_Dist->GetYaxis()->SetTitleOffset(1.1);
  histo_GEN_Dist->SetMarkerColor(my_colors[0]);
  histo_GEN_Dist->SetStats(0);
  histo_GEN_Dist->Draw();
   
  canv_2d->cd(3);
  gPad->SetGrid();  
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.08);
  gPad->SetFillColor(0);

  histo_FRT_Pos->GetXaxis()->CenterTitle();
  histo_FRT_Pos->GetYaxis()->CenterTitle();
  histo_FRT_Pos->GetZaxis()->CenterTitle();
  
  histo_FRT_Pos->GetXaxis()->SetTitleOffset(1.1);
  histo_FRT_Pos->GetYaxis()->SetTitleOffset(1.1);
  histo_FRT_Pos->SetMarkerColor(my_colors[0]);
  histo_FRT_Pos->SetStats(0);
  histo_FRT_Pos->Draw();
   
  canv_2d->cd(4);
  gPad->SetGrid();  
  gPad->SetRightMargin(0.12);
  gPad->SetLeftMargin(0.08);
  gPad->SetFillColor(0);

  histo_FRT_Dist->GetXaxis()->CenterTitle();
  histo_FRT_Dist->GetYaxis()->CenterTitle();
  histo_FRT_Dist->GetZaxis()->CenterTitle();
  
  histo_FRT_Dist->GetXaxis()->SetTitleOffset(1.1);
  histo_FRT_Dist->GetYaxis()->SetTitleOffset(1.1);
  histo_FRT_Dist->SetMarkerColor(my_colors[0]);
  histo_FRT_Dist->SetStats(0);
  histo_FRT_Dist->Draw();
  
  canv_2d->SaveAs("AssociationQualityAnalysis3_2D.png");
} 