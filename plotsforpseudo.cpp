#include "plotstyle.hpp"
#include "cmsstyle.hpp"
#include "helpers.hpp"
#include "binning.hpp"
#include "TPaveStats.h"

#include "TNtuple.h"

int plotsforpseudo()
{
  TFile* inputfile = new TFile(TString("output/")+fileSuffix+"/pseudo.root");

  CMSStyle* cmsstyle = new CMSStyle();
  CMSStyle* cstyle = new CMSStyle();
  //CStyle* cstyle = new CStyle(2);
  //cstyle->cd();  
  

  gStyle->SetOptStat("RM");
  //gStyle->SetStatH(0.5);
  //gStyle->SetStatW(0.2);
  
  TCanvas* canv1 = new TCanvas("canv1","canv1",300,0,800,600);
  //canv1->SetBatch(true);
  //canv1->SetLeftMargin(0.2);
  canv1->SetRightMargin(0.04);
  

  TH1F* h1;
  TH1F* h2;
  
  /// asydifference
  h1 = (TH1F*)inputfile->Get("asydifference");
  cstyle->setup_style(h1, "A_{C,2d} - A_{C,1d}", "N", 1, 0, 1);
  
  h1->SetMaximum(h1->GetMaximum()*1.2);
  h1->Draw("HIST");
  //cmsstyle->CMSSimulation();
    
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_asydifference"));
  
  /// asys 1d overlayed over 2d
  TPaveStats* st;
  
  h1 = (TH1F*)inputfile->Get("asys");
  h2 = (TH1F*)inputfile->Get("asys1d");
  cstyle->setup_style(h1, "A_{C}", "N", 4, 0, 1);
  cstyle->setup_style(h2, "A_{C}", "N", 2, 0, 1);
  
  h1->SetMaximum(h2->GetMaximum()*1.2);
  h1->Draw("HIST");
  h2->Draw("HIST, SAMES");
  
  gPad->Update(); // to make sure the statistics boxes are created
  
  st = (TPaveStats*)h2->FindObject("stats");
  double height = st->GetY2NDC() - st->GetY1NDC();
  double offset = 0.2;
  st->SetY1NDC(st->GetY1NDC()-offset);
  st->SetY2NDC(st->GetY2NDC()-offset);

  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_asys"));
  

  /// to make z-axis stuff show completely in correlation plots
  canv1->SetRightMargin(0.14);
  
  
  /// 1d vs 2d
  TH2F* h2_1 = (TH2F*)inputfile->Get("asy1dvs2d");
  cstyle->setup_style_2D(h2_1, "A_{C} 1d unfolded", "A_{C} 2d unfolded");
  //h2_1->SetMaximum(0.07);
  h2_1->GetXaxis()->SetNdivisions(505, kTRUE);
  h2_1->GetYaxis()->SetNdivisions(505, kTRUE);
  h2_1->Draw("COLZ");
  //cmsstyle->CMSSimulation_2D();
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_asy1dvs2d"));
  
  
  /// asys  
  gStyle->SetOptStat("rm");
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.3);
  double statFontSize = 0.045;
  TNtuple* asytup = (TNtuple*)inputfile->Get("trueasys");
  float trueasy = 0.1234;
  asytup->SetBranchAddress("trueasy", &trueasy);
  
  h1 = (TH1F*)inputfile->Get("asys1d");
  h1->Rebin(2);
  cstyle->setup_style(h1, "A_{C}", "N", 4, 0, 1);
  h1->SetMaximum(h1->GetMaximum()*1.5);
  h1->GetXaxis()->SetRangeUser(-0.1, 0.1);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  asytup->GetEntry(0);
  cmsstyle->TextforCMS("inclusive",0.5);
  cmsstyle->TextforCMS(TString::Format("Simulated A_{C}=%.4f", trueasy).Data(),1.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_asys1d"));
  
  h1 = (TH1F*)inputfile->Get("asys");
  h1->Rebin(2);
  cstyle->setup_style(h1, "A_{C}", "N", 4, 0, 1);
  h1->SetMaximum();
  h1->SetMaximum(h1->GetMaximum()*1.5);
  h1->GetXaxis()->SetRangeUser(-0.1, 0.1);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  cmsstyle->TextforCMS("inclusive",0.5);
  cmsstyle->TextforCMS(TString::Format("Simulated A_{C}=%.4f", trueasy).Data(),1.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_asys2d"));
  
  h1 = (TH1F*)inputfile->Get("diffasybin0");
  h1->Rebin(2);
  cstyle->setup_style(h1, "A_{C}", "N", 4, 0, 1);
  h1->SetMaximum(h1->GetMaximum()*1.5);
  h1->GetXaxis()->SetRangeUser(-0.2, 0.2);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  asytup->GetEntry(1);
  cmsstyle->TextforCMS(TString::Format("%s bin 1", labelOfXAxisVar_noUnit).Data(),0.5);
  cmsstyle->TextforCMS(TString::Format("Simulated A_{C}=%.4f", trueasy).Data(),1.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_asysbin0"));
  
  h1 = (TH1F*)inputfile->Get("diffasybin1");
  h1->Rebin(2);
  cstyle->setup_style(h1, "A_{C}", "N", 4, 0, 1);
  h1->SetMaximum(h1->GetMaximum()*1.5);
  h1->GetXaxis()->SetRangeUser(-0.2, 0.2);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  asytup->GetEntry(2);
  cmsstyle->TextforCMS(TString::Format("%s bin 2", labelOfXAxisVar_noUnit).Data(),0.5);
  cmsstyle->TextforCMS(TString::Format("Simulated A_{C}=%.4f", trueasy).Data(),1.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_asysbin1"));
  
  h1 = (TH1F*)inputfile->Get("diffasybin2");
  h1->Rebin(2);
  cstyle->setup_style(h1, "A_{C}", "N", 4, 0, 1);
  h1->SetMaximum(h1->GetMaximum()*1.5);
  h1->GetXaxis()->SetRangeUser(-0.2, 0.2);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  asytup->GetEntry(3);
  cmsstyle->TextforCMS(TString::Format("%s bin 3", labelOfXAxisVar_noUnit).Data(),0.5);
  cmsstyle->TextforCMS(TString::Format("Simulated A_{C}=%.4f", trueasy).Data(),1.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_asysbin2"));
  
  
  /// pulls  
  h1 = (TH1F*)inputfile->Get("asypull1d");
  h1->Rebin(2);
  cstyle->setup_style(h1, "pull", "N", 4, 0, 1);
  h1->SetMaximum(h1->GetMaximum()*1.5);
  //h1->GetXaxis()->SetRangeUser(-0.1, 0.1);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  cmsstyle->TextforCMS("inclusive",0.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_pull1d"));
  
  h1 = (TH1F*)inputfile->Get("asypull");
  h1->Rebin(2);
  cstyle->setup_style(h1, "pull", "N", 4, 0, 1);
  h1->SetMaximum(h1->GetMaximum()*1.5);
  //h1->GetXaxis()->SetRangeUser(-0.1, 0.1);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  cmsstyle->TextforCMS("inclusive",0.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_pull2d"));
  
  h1 = (TH1F*)inputfile->Get("diffasypullbin0");
  h1->Rebin(2);
  cstyle->setup_style(h1, "pull", "N", 4, 0, 1);
  h1->SetMaximum(h1->GetMaximum()*1.5);
  //h1->GetXaxis()->SetRangeUser(-0.2, 0.2);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  cmsstyle->TextforCMS(TString::Format("%s bin 1", labelOfXAxisVar_noUnit).Data(),0.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_pullbin0"));
  
  h1 = (TH1F*)inputfile->Get("diffasypullbin1");
  h1->Rebin(2);
  cstyle->setup_style(h1, "pull", "N", 4, 0, 1);
  h1->SetMaximum(h1->GetMaximum()*1.5);
  //h1->GetXaxis()->SetRangeUser(-0.2, 0.2);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  cmsstyle->TextforCMS(TString::Format("%s bin 2", labelOfXAxisVar_noUnit).Data(),0.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_pullbin1"));
  
  h1 = (TH1F*)inputfile->Get("diffasypullbin2");
  h1->Rebin(2);
  cstyle->setup_style(h1, "pull", "N", 4, 0, 1);
  h1->SetMaximum(h1->GetMaximum()*1.5);
  //h1->GetXaxis()->SetRangeUser(-0.2, 0.2);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  cmsstyle->TextforCMS(TString::Format("%s bin 3", labelOfXAxisVar_noUnit).Data(),0.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_pullbin2"));
  
  
  
  /// stat errors  
  h1 = (TH1F*)inputfile->Get("asyerrors1d");
  cstyle->setup_style(h1, "#sigma_{A_{C}}", "N", 1, 0, 1);
  h1->SetMaximum(h1->GetMaximum()*1.5);
  //h1->GetXaxis()->SetRangeUser(-0.1, 0.1);
  h1->GetXaxis()->SetRangeUser(0., 0.1);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  cmsstyle->TextforCMS("inclusive",0.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_staterrors1d"));
  
  h1 = (TH1F*)inputfile->Get("asyerrors");
  cstyle->setup_style(h1, "#sigma_{A_{C}}", "N", 1, 0, 1);
  h1->SetMaximum(h1->GetMaximum()*1.5);
  //h1->GetXaxis()->SetRangeUser(-0.1, 0.1);
  h1->GetXaxis()->SetRangeUser(0., 0.5);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  cmsstyle->TextforCMS("inclusive",0.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_staterrors2d"));
  
  h1 = (TH1F*)inputfile->Get("diffasyerrorbin0");
  cstyle->setup_style(h1, "#sigma_{A_{C}}", "N", 1, 0, 1);
  h1->SetMaximum(h1->GetMaximum()*1.5);
  //h1->GetXaxis()->SetRangeUser(-0.2, 0.2);
  h1->GetXaxis()->SetRangeUser(0., 0.5);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  cmsstyle->TextforCMS(TString::Format("%s bin 1", labelOfXAxisVar_noUnit).Data(),0.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_staterrorsbin0"));
  
  h1 = (TH1F*)inputfile->Get("diffasyerrorbin1");
  cstyle->setup_style(h1, "#sigma_{A_{C}}", "N", 1, 0, 1);
  h1->SetMaximum(h1->GetMaximum()*1.5);
  //h1->GetXaxis()->SetRangeUser(-0.2, 0.2);
  h1->GetXaxis()->SetRangeUser(0., 0.5);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  cmsstyle->TextforCMS(TString::Format("%s bin 2", labelOfXAxisVar_noUnit).Data(),0.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_staterrorsbin1"));
  
  h1 = (TH1F*)inputfile->Get("diffasyerrorbin2");
  cstyle->setup_style(h1, "#sigma_{A_{C}}", "N", 1, 0, 1);
  h1->SetMaximum(h1->GetMaximum()*1.5);
  //h1->GetXaxis()->SetRangeUser(-0.2, 0.2);
  h1->GetXaxis()->SetRangeUser(0., 0.5);
  h1->Draw("HIST");
  gPad->Update(); // to make sure the statistics boxes are created
  st = (TPaveStats*)h1->FindObject("stats");
  st->SetY1NDC(1-canv1->GetTopMargin()-0.2);
  st->SetY2NDC(1-canv1->GetTopMargin());
  st->SetX1NDC(1-canv1->GetRightMargin()-0.3);
  st->SetX2NDC(1-canv1->GetRightMargin());
  st->SetTextSize(statFontSize);
  cmsstyle->TextforCMS(TString::Format("%s bin 3", labelOfXAxisVar_noUnit).Data(),0.5);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_staterrorsbin2"));
  
  /// ////////////////////////////////
  inputfile->Close();
  delete inputfile;
   
  inputfile = new TFile(TString("output/")+fileSuffix+"/lincheck.root");
  
  double max, min, border;
  TF1* bisector;
  
  /// linchecks
  h1 = (TH1F*)inputfile->Get("hlincheck");
  max = h1->GetMaximum();
  min = h1->GetMinimum();
  border = (max-min)*0.2;
  cstyle->setup_style(h1, "A_{C,gen}", "A_{C,unfolded}", 1, 0, 0);
  //gStyle->SetOptStat("RM");
  //h1->SetMaximum(max+border);
  //h1->SetMinimum(min-border);
  h1->SetMaximum(max+border);
  h1->SetMinimum(min-border);
  h1->SetAxisRange(min-border, max+border);
  h1->Draw("E");
  bisector = new TF1("bisector","x", h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
  bisector->SetLineColor(4);
  bisector->Draw("SAME");
  cmsstyle->CMSSimulation();
  cmsstyle->TextforCMS("inclusive",1.0);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_lincheck2d"));  
  delete bisector;
  
  h1 = (TH1F*)inputfile->Get("hlincheck1d");
  max = h1->GetMaximum();
  min = h1->GetMinimum();
  border = (max-min)*0.2;
  cstyle->setup_style(h1, "A_{C,gen}", "A_{C,unfolded}", 1, 0, 0);
  //gStyle->SetOptStat("RM");
  //h1->SetMaximum(max+border);
  //h1->SetMinimum(min-border);
  h1->SetMaximum(max+border);
  h1->SetMinimum(min-border);
  h1->SetAxisRange(min-border, max+border);
  h1->Draw("E");
  bisector = new TF1("bisector","x", h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
  bisector->SetLineColor(4);
  bisector->Draw("SAME");
  cmsstyle->CMSSimulation();
  cmsstyle->TextforCMS("inclusive",1.0);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_lincheck1d"));  
  delete bisector;
  
  
  h1 = (TH1F*)inputfile->Get("lincheck_diffasybin0"); 
  max = h1->GetMaximum();
  min = h1->GetMinimum();
  border = (max-min)*1.;
  cstyle->setup_style(h1, "A_{C,gen}", "A_{C,unfolded}", 1, 0, 0);
  //gStyle->SetOptStat("RM");
  //h1->SetMaximum(max+border);
  //h1->SetMinimum(min-border);
  h1->SetMaximum(max+border);
  h1->SetMinimum(min-border);
  h1->SetAxisRange(min-border, max+border);
  h1->Draw("E");
  bisector = new TF1("bisector","x", h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
  bisector->SetLineColor(4);
  bisector->Draw("SAME");
  cmsstyle->CMSSimulation();
  cmsstyle->TextforCMS(TString::Format("%s bin 1", labelOfXAxisVar_noUnit).Data(),1.0);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_lincheckbin0"));
  delete bisector;
  
  
  h1 = (TH1F*)inputfile->Get("lincheck_diffasybin1"); 
  max = h1->GetMaximum();
  min = h1->GetMinimum();
  border = (max-min)*0.7;
  cstyle->setup_style(h1, "A_{C,gen}", "A_{C,unfolded}", 1, 0, 0);
  //gStyle->SetOptStat("RM");
  //h1->SetMaximum(max+border);
  //h1->SetMinimum(min-border);
  h1->SetMaximum(max+border);
  h1->SetMinimum(min-border);
  h1->SetAxisRange(min-border, max+border);
  h1->Draw("E");
  bisector = new TF1("bisector","x", h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
  bisector->SetLineColor(4);
  bisector->Draw("SAME");
  cmsstyle->CMSSimulation();
  cmsstyle->TextforCMS(TString::Format("%s bin 2", labelOfXAxisVar_noUnit).Data(),1.0);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_lincheckbin1"));
  delete bisector;
  
  
  h1 = (TH1F*)inputfile->Get("lincheck_diffasybin2"); 
  max = h1->GetMaximum();
  min = h1->GetMinimum();
  border = (max-min)*0.2;
  cstyle->setup_style(h1, "A_{C,gen}", "A_{C,unfolded}", 1, 0, 0);
  //gStyle->SetOptStat("RM");
  //h1->SetMaximum(max+border);
  //h1->SetMinimum(min-border);
  h1->SetMaximum(max+border);
  h1->SetMinimum(min-border);
  h1->SetAxisRange(min-border, max+border);
  h1->Draw("E");
  bisector = new TF1("bisector","x", h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
  bisector->SetLineColor(4);
  bisector->Draw("SAME");
  cmsstyle->CMSSimulation();
  cmsstyle->TextforCMS(TString::Format("%s bin 3", labelOfXAxisVar_noUnit).Data(),1.0);
  saveAs(canv1, (TString("output/")+fileSuffix+"plotsforpseudo_lincheckbin2"));
  delete bisector;  
  
  return 0;
}


int main()
{
  plotsforpseudo();
  return 0;
}
