#ifndef CSTYLE_HPP
#define CSTYLE_HPP

#include <iostream>
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TColor.h"

// various includes for compiling
#include "TGraph.h"
#include "TFile.h"
#include "TPaveStats.h"
#include <string>

using namespace std;

class CStyle{

public:
  
  CStyle(int width){
      
    LINE_WIDTH = width;

    gROOT->SetStyle("Plain");

    cStyle.SetPalette(1);

    cStyle.SetLineWidth(LINE_WIDTH);
    cStyle.SetFrameLineWidth(LINE_WIDTH);

    // For the canvas:
    cStyle.SetCanvasBorderMode(0);
    cStyle.SetCanvasColor(kWhite);
    cStyle.SetCanvasDefH(600); //Height of canvas  
    cStyle.SetCanvasDefW(600); //Width of canvas
    cStyle.SetCanvasDefX(0);   //POsition on screen
    cStyle.SetCanvasDefY(0);

    // For the Pad:
    cStyle.SetPadBorderMode(0);
    cStyle.SetPadColor(kWhite);
    cStyle.SetPadGridX(false);
    cStyle.SetPadGridY(false);
    cStyle.SetGridColor(0);
    cStyle.SetGridStyle(3);
    cStyle.SetGridWidth(1);

    // For the frame:
    cStyle.SetFrameBorderMode(0);
    cStyle.SetFrameBorderSize(1);
    cStyle.SetFrameFillColor(0);
    cStyle.SetFrameFillStyle(0);
    cStyle.SetFrameLineColor(1);
    cStyle.SetFrameLineStyle(1);
    cStyle.SetFrameLineWidth(2);

    // For the histo:
    cStyle.SetHistLineColor(1);
    cStyle.SetHistLineStyle(0);
    cStyle.SetHistLineWidth(2);

    cStyle.SetEndErrorSize(2);
    cStyle.SetErrorX(0.);
    cStyle.SetMarkerStyle(20);

    //For the date:
    cStyle.SetOptDate(0);

    // For the statistics box:
    cStyle.SetOptFile(0);
    cStyle.SetOptStat("mr"); // To display the mean and RMS:   SetOptStat("mr");
    cStyle.SetOptStat(1100);
    cStyle.SetStatColor(kWhite);
    cStyle.SetStatFont(62);
    cStyle.SetStatFontSize(0.035);
    cStyle.SetStatTextColor(1);
    cStyle.SetStatFormat("6.4g");
    cStyle.SetStatBorderSize(0);
    cStyle.SetStatH(0.25);
    cStyle.SetStatW(0.25); // 0.35

    // Margins:
    cStyle.SetPadTopMargin(0.07);
    cStyle.SetPadBottomMargin(0.2);
    cStyle.SetPadLeftMargin(0.26);
    cStyle.SetPadRightMargin(0.16);

    // For the Global title:
    cStyle.SetOptTitle(0);
    cStyle.SetTitleFont(62);
    cStyle.SetTitleColor(1);
    cStyle.SetTitleTextColor(1);
    cStyle.SetTitleFillColor(kWhite);
    cStyle.SetTitleFontSize(0.05);

    // For the axis titles:
    cStyle.SetTitleColor(1, "XYZ");
    cStyle.SetTitleFont(62, "XYZ");
    cStyle.SetTitleSize(0.07, "XYZ");
    cStyle.SetTitleXOffset(1.);
    cStyle.SetTitleYOffset(1.2);

    // For the axis labels:
    cStyle.SetLabelColor(1, "XYZ");
    cStyle.SetLabelFont(62, "XYZ");
    cStyle.SetLabelOffset(0.015, "XYZ");
    cStyle.SetLabelSize(0.07, "XYZ");

    // For the axis:
    cStyle.SetAxisColor(1, "XYZ");
    cStyle.SetStripDecimals(kTRUE);
    cStyle.SetTickLength(0.03, "XYZ");
    cStyle.SetNdivisions(605, "XYZ");
    //  cStyle.SetNdivisions(405, "X"); // for Njets plots

    // Change for log plots:
    cStyle.SetOptLogx(0);
    cStyle.SetOptLogy(0);
    cStyle.SetOptLogz(0);
    cStyle.SetPalette(1);

    //For the legend
    cStyle.SetLegendBorderSize(1);

    TColor KIT(10000, 255, 255, 255, "KIT", 1);
    KIT_gruen=KIT.GetColor(128, 187, 61);
    KIT_gelb=KIT.GetColor( 244, 235, 22);
    KIT_orange=KIT.GetColor(221, 164, 33);
    KIT_braun=KIT.GetColor(160, 134, 49);
    KIT_rot=KIT.GetColor(155, 29, 45);
    KIT_lila=KIT.GetColor(159, 0, 121);
    KIT_blau=KIT.GetColor(92, 112, 176);
    KIT_cyan=KIT.GetColor(77, 166, 227);
    KIT_logo=KIT.GetColor(0, 154, 131);
        
    cStyle.cd();

  };

 // ~CStyle();

  void cd(){
    cStyle.cd();
  };

  void setup_global_style_ex1()
  {
    gROOT->SetStyle("Plain");
    //cStyle.SetOptStat(0);
    
    cStyle.SetPadTopMargin(0.15);
    cStyle.SetPadBottomMargin(0.2);
    cStyle.SetPadLeftMargin(0.23);
    cStyle.SetPadRightMargin(0.12);

    cStyle.SetPalette(1);

    cStyle.SetLineWidth(LINE_WIDTH);
    cStyle.SetFrameLineWidth(LINE_WIDTH);
  }

  void setup_global_style_fit()
  {
    gROOT->SetStyle("Plain");
    //cStyle.SetOptStat(0);
    
    cStyle.SetPadTopMargin(0.15);
    cStyle.SetPadBottomMargin(0.2);
    cStyle.SetPadLeftMargin(0.23);
    cStyle.SetPadRightMargin(0.05);

    cStyle.SetPalette(1);

    cStyle.SetLineWidth(LINE_WIDTH);
    cStyle.SetFrameLineWidth(LINE_WIDTH);
  }

  void setup_style(TH1* histo, std::string xtitle, std::string ytitle, int lcolor, bool fill, bool StatBox)
  {
    histo->SetTitle("");      
    histo->SetXTitle(xtitle.c_str());
    histo->SetYTitle(ytitle.c_str());
    if(fill){histo->SetFillColor(lcolor);}
    else{
      histo->SetFillColor(0);
      histo->SetLineColor(lcolor);}

    histo->SetNdivisions(505,"x");
    histo->SetNdivisions(505,"y");

    histo->SetLabelOffset(0.01,"x");
    histo->SetLabelSize(0.03,"x"); // 0.07
    histo->SetTitleSize(0.03,"x");
    histo->SetTitleOffset(1.15,"x");
    
    histo->SetLabelOffset(0.01,"y");
    histo->SetLabelSize(0.03,"y");
    histo->SetTitleSize(0.03,"y");
    histo->SetTitleOffset(1.5,"y");
    histo->SetLineWidth(LINE_WIDTH);
    histo->SetMinimum(0);
    
    histo->SetStats(StatBox);

  }
  void setup_style_data(TH1* histo, std::string xtitle, std::string ytitle, int lcolor, bool fill, bool StatBox)
  {
    histo->SetMarkerStyle(20);
    histo->SetMarkerSize(1);
    histo->SetMarkerColor(lcolor);
    histo->SetTitle("");      
    histo->SetXTitle(xtitle.c_str());
    histo->SetYTitle(ytitle.c_str());
    if(fill){histo->SetFillColor(lcolor);}
    else{
      histo->SetFillColor(0);
      histo->SetLineColor(lcolor);}

    histo->SetNdivisions(505,"x");
    histo->SetNdivisions(505,"y");

    histo->SetLabelOffset(0.01,"x");
    histo->SetLabelSize(0.07,"x");
    histo->SetTitleSize(0.07,"x");
    histo->SetTitleOffset(1.15,"x");
    
    histo->SetLabelOffset(0.01,"y");
    histo->SetLabelSize(0.07,"y");
    histo->SetTitleSize(0.07,"y");
    histo->SetTitleOffset(1.5,"y");
    histo->SetLineWidth(LINE_WIDTH);
    histo->SetMinimum(0);
    
    histo->SetStats(StatBox);

  }



  void setup_style_array(TH1F* histo[],int groesse, std::string xtitle, std::string ytitle, int lcolor, bool fill, bool StatBox)
  {
    for(int i=0; i<groesse; i++){

      histo[i]->SetMarkerStyle(1);
      histo[i]->SetTitle("");      
      histo[i]->SetXTitle(xtitle.c_str());
      histo[i]->SetYTitle(ytitle.c_str());
      if(fill){histo[i]->SetFillColor(lcolor);}
      else{
	histo[i]->SetFillColor(0);
	histo[i]->SetLineColor(lcolor);}
      
      histo[i]->SetNdivisions(505,"x");
      histo[i]->SetNdivisions(505,"y");
      
      histo[i]->SetLabelOffset(0.01,"x");
      histo[i]->SetLabelSize(0.07,"x");
      histo[i]->SetTitleSize(0.07,"x");
      histo[i]->SetTitleOffset(1.15,"x");
      
      histo[i]->SetLabelOffset(0.01,"y");
      histo[i]->SetLabelSize(0.07,"y");
      histo[i]->SetTitleSize(0.07,"y");
      histo[i]->SetTitleOffset(1.5,"y");
      histo[i]->SetLineWidth(LINE_WIDTH);
      histo[i]->SetMinimum(0);
      
      histo[i]->SetStats(StatBox);
    }
    
  }

  void setup_TGraph(TGraph* graph, std::string xtitle, std::string ytitle, int lcolor){
    graph->SetLineColor(lcolor);
    graph->SetLineWidth(LINE_WIDTH);
    
    // graph->Draw();
    TAxis *xaxis = graph->GetXaxis();
    TAxis *yaxis = graph->GetYaxis();
    xaxis->SetTitle(xtitle.c_str());
    yaxis->SetTitle(ytitle.c_str());
    xaxis->SetTitleOffset(1.1);
    yaxis->SetTitleOffset(1.1);

    
  }

  void setup_style_2D(TH2* histo, std::string xtitle, std::string ytitle)
  {
    histo->SetTitle("");      
    histo->SetXTitle(xtitle.c_str());
    histo->SetYTitle(ytitle.c_str());
    histo->SetLineWidth(LINE_WIDTH);

    histo->SetNdivisions(505,"x");
    histo->SetNdivisions(505,"y");

    histo->SetLabelOffset(0.01,"x");
    histo->SetLabelSize(0.07,"x");
    histo->SetTitleSize(0.07,"x");
    histo->SetTitleOffset(1.15,"x");
   
    histo->SetLabelOffset(0.01,"y");
    histo->SetLabelSize(0.07,"y");
    histo->SetTitleSize(0.07,"y");
    histo->SetTitleOffset(1.2,"y");

    histo->SetMinimum(0);

    histo->SetStats(0);
  }

  Color_t KIT_gruen;
  Color_t KIT_gelb;
  Color_t KIT_orange;
  Color_t KIT_braun;
  Color_t KIT_rot;
  Color_t KIT_lila;
  Color_t KIT_blau;
  Color_t KIT_cyan;
  Color_t KIT_logo;
  
private:
  TStyle cStyle;
  int LINE_WIDTH;
  
};

#endif
