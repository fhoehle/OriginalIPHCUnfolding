#include <iostream>
#include <assert.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TSpline.h>
#include <TUnfold.h>
#include <TUnfoldSys.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <sstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <TPave.h>
#include <TPaveText.h>
#include <TVirtualPad.h>
#include <TClass.h>
#include "TApplication.h"
#include "TRandom3.h"
#include <fstream>



#include "cmsstyle.hpp"
#include "helpers.hpp"
#include "specialhelpers.hpp"
#include "samplelocations.hpp"


#include "binning.hpp"

// WARNING make sure variables in binning.hpp are set correctly


int migrationmatrix()
{
  TH1::SetDefaultSumw2(true);  
//  TFile* selefffile = new TFile(TString("output/")+fileSuffix+"/"+dimstring+"seleff.root");
  TFile* selefffile = new TFile(TString("output/")+dimstring+"seleff.root");

  CMSStyle* cmsstyle = new CMSStyle();
  gStyle->SetOptStat("RM");

  TCanvas* canv1 = new TCanvas("canv1","canv1",300,0,800,600);
  canv1->SetRightMargin(0.04);
  
  // to make z-axis stuff show completely in 2d plots
  canv1->SetRightMargin(0.14);
  
  /// generate 1d histos for the 2 variables with correct binning, content doesnt matter
  H2gen* genbinning = new H2gen(xedgesgen, yedgesgen); 
  H2rec* recbinning = new H2rec(xedgesrec, yedgesrec); 
  
  H2rec recoplotele(xedgesrec, yedgesrec); 
  //  H2rec recoplotmu(xedgesrec, yedgesrec); 

  TH2F* migmatrixcombined = new TH2F("migmatrixcombined","migmatrixcombined", nbinsafter, 0.5, nbinsafter+0.5, nbins, 0.5, nbins+0.5);
  TH2F* migmatrixele = new TH2F("migmatrixele","migmatrixele", nbinsafter, 0.5, nbinsafter+0.5, nbins, 0.5, nbins+0.5);
  //  TH2F* migmatrixmu = new TH2F("migmatrixmu","migmatrixmu", nbinsafter, 0.5, nbinsafter+0.5, nbins, 0.5, nbins+0.5);
  
  if(nbinsmafter != 1)
  {
    // FIXME uses fixed number of m bins
    cmsstyle->setup_style_2D(migmatrixcombined, (TString("bins of #Delta|y|_{gen} in 3 bins of ")+labelOfXAxisVar_noUnit_gen).Data(),
                                                (TString("bins of #Delta|y|_{rec} in 6 bins of ")+labelOfXAxisVar_noUnit_rec).Data());
    cmsstyle->setup_style_2D(migmatrixele, (TString("bins of #Delta|y|_{gen} in 3 bins of ")+labelOfXAxisVar_noUnit_gen).Data(),
                                           (TString("bins of #Delta|y|_{rec} in 6 bins of ")+labelOfXAxisVar_noUnit_rec).Data());
    //    cmsstyle->setup_style_2D(migmatrixmu, (TString("bins of #Delta|y|_{gen} in 3 bins of ")+labelOfXAxisVar_noUnit_gen).Data(),
    //                                    (TString("bins of #Delta|y|_{rec} in 6 bins of ")+labelOfXAxisVar_noUnit_rec).Data());
 
    TAxis* axis = migmatrixcombined->GetXaxis();
    for(int i=1; i<nbinsetaafter+1; i++)
    {
      axis->SetBinLabel(i,TString::Format("%d",i).Data());
      axis->SetBinLabel(i+nbinsetaafter,TString::Format("%d",i).Data());
      axis->SetBinLabel(i+nbinsetaafter*2,TString::Format("%d",i).Data());
    }
    axis->SetTickLength(0);
    
    axis = migmatrixcombined->GetYaxis();
    for(int i=1; i<nbinseta+1; i++)
    {
      if(i%4) continue;
      axis->SetBinLabel(i,TString::Format("%d",i).Data());
      axis->SetBinLabel(i+nbinseta,TString::Format("%d",i).Data());
      axis->SetBinLabel(i+nbinseta*2,TString::Format("%d",i).Data());
      axis->SetBinLabel(i+nbinseta*3,TString::Format("%d",i).Data());
      axis->SetBinLabel(i+nbinseta*4,TString::Format("%d",i).Data());
      axis->SetBinLabel(i+nbinseta*5,TString::Format("%d",i).Data());      
    }
    axis->SetTickLength(0);
    
  }
  else
  {
    cmsstyle->setup_style_2D(migmatrixcombined, "bin of generated |y_{t}|-|y_{#bar{t}}|", "bin of reconstructed |y_{t}|-|y_{#bar{t}}|");
    cmsstyle->setup_style_2D(migmatrixele, "bin of generated |y_{t}|-|y_{#bar{t}}|", "bin of reconstructed |y_{t}|-|y_{#bar{t}}|");
    //    cmsstyle->setup_style_2D(migmatrixmu, "bin of generated |y_{t}|-|y_{#bar{t}}|", "bin of reconstructed |y_{t}|-|y_{#bar{t}}|");
    // FIXME uses fixed numbers of bins
    TAxis* axis = migmatrixcombined->GetXaxis();
    axis->SetBinLabel(1,"1");
    axis->SetBinLabel(2,"2");
    axis->SetBinLabel(3,"3");
    axis->SetBinLabel(4,"4");
    axis->SetBinLabel(5,"5");
    axis->SetBinLabel(6,"6");
    axis->SetBinLabel(7,"7");
    axis->SetBinLabel(8,"8");
    
    axis = migmatrixcombined->GetYaxis();
    axis->SetBinLabel(1,"1");
    axis->SetBinLabel(2,"2");
    axis->SetBinLabel(3,"3");
    axis->SetBinLabel(4,"4");
    axis->SetBinLabel(5,"5");
    axis->SetBinLabel(6,"6");
    axis->SetBinLabel(7,"7");
    axis->SetBinLabel(8,"8");
    axis->SetBinLabel(9,"9");
    axis->SetBinLabel(10,"10");
    axis->SetBinLabel(11,"11");
    axis->SetBinLabel(12,"12");
    axis->SetBinLabel(13,"13");
    axis->SetBinLabel(14,"14");
    axis->SetBinLabel(15,"15");
    axis->SetBinLabel(16,"16");
  }
  
  TH2F* visuplots[nbinsmafter];
  for(int i=0; i<nbinsmafter; i++)
  {
    visuplots[i] = new TH2F("","", 100, yedgesgen[0][0], yedgesgen[0][nbinsetaafter], 100, yedgesgen[0][0], yedgesgen[0][nbinsetaafter]);
  }
  
  TH1F* sensrecoplotinc = new TH1F("","", 100, yedgesrec[0][0], yedgesrec[0][nbinseta]);
  
  
  float mttbargen, mttbarrec, diffabsetagen, diffabsetarec, weight;
  

  #ifdef POWHEG
    int nPostSamples = 1;
  #endif
  
  for(int j=0; j<nPostSamples; j++) // do this for both ele and mu
  {
    TH2F* migmatrix;
    TNtuple* tuple;
    TH2F* nonselected;
    H2rec* recoplot;
    //double sampleweight = 1;
    
    TFile* inputfile = 0;

    
    #ifdef POWHEG
      if(j==0) // ele
      {
        migmatrix = migmatrixele;
        if(!random_split)
          inputfile = new TFile(sl_ttbarele);
        else
          inputfile = new TFile("FIXME/home/froscher/outputroot/powhegsample/ttbar_ele_combined_subsample1_histos.root");
        nonselected = (TH2F*)selefffile->Get("nonselectedele");
        recoplot = &recoplotele;
      }
      /*
      else // mu
      {
        migmatrix = migmatrixmu;
        if(!random_split)
          inputfile = new TFile(sl_ttbarmu);
        else
          inputfile = new TFile("FIXME/home/froscher/outputroot/powhegsample/ttbar_mu_combined_subsample1_histos.root");
        nonselected = (TH2F*)selefffile->Get("nonselectedmu");
        recoplot = &recoplotmu;
	}*/
    #endif
    
//    tuple = (TNtuple*)inputfile->Get("TTbarSelHyppsiprime/BBB_UnfoldingTuple");
    tuple = (TNtuple*)inputfile->Get("Tuple");
                                        
    tuple->SetBranchAddress(nameOfXAxisVarGen, &mttbargen);
    tuple->SetBranchAddress(nameOfXAxisVarRec, &mttbarrec);
    tuple->SetBranchAddress(nameOfSensVarGen, &diffabsetagen);
    tuple->SetBranchAddress(nameOfSensVarRec, &diffabsetarec);
    tuple->SetBranchAddress("weight", &weight);
    weight *= lumi_factor;
    
    migmatrix->Sumw2();
    
    int nentries = tuple->GetEntries();
    
    // only use part of events so unfolding can be done on statistically independent sample
    int nentries_toconsider = nentries; 
    int nentries_toskip = 0;
    
    // fill main content of migmatrix
    for(int i=nentries_toskip; i<nentries_toconsider; i++)
    {

      
      tuple->GetEntry(i);
      //weight *= sampleweight;
      
      recoplot->fill(fabs(mttbarrec), diffabsetarec, weight);
      
      
      const int genbinm = genbinning->findXBin(fabs(mttbargen));
      const int genbineta = genbinning->findYBin(genbinm, diffabsetagen);
      const int genbin = genbinm * nbinsetaafter + genbineta + 1;
      
      const int recbinm = recbinning->findXBin(fabs(mttbarrec));
      const int recbineta = recbinning->findYBin(recbinm, diffabsetarec);
      const int recbin = recbinm * nbinseta + recbineta + 1; 
      
      migmatrix->Fill(genbin, recbin, weight);
      
      visuplots[genbinm]->Fill(diffabsetagen, diffabsetarec, weight);
      sensrecoplotinc->Fill(diffabsetarec, weight);
      
  
    }
    
    migmatrix->Sumw2();


    // fill overflow bins of migmatrix for selection efficiency
    TH1F nonselected1d("nonselected1d","nonselected1d", nbinsafter, 0.5, nbinsafter+0.5); 

    unwrap2dhisto(nonselected, &nonselected1d);

    for(int i=0; i<nbinsafter; i++)
    {
      double content = nonselected1d.GetBinContent(i+1);
      double error = nonselected1d.GetBinError(i+1);  
      migmatrix->SetBinContent(i+1, 0, content);
      migmatrix->SetBinError(i+1, 0, error);
    }
  
  }
  
  // combine matrices in right proportion (same as in fit to data).
  // overflow bins will be rewritten in next step.
  TH2F* postNges = (TH2F*)selefffile->Get("postNges");
  TH2F* postNele = (TH2F*)selefffile->Get("postNele");
  //  TH2F* postNmu = (TH2F*)selefffile->Get("postNmu");
  
  migmatrixcombined->Add(migmatrixele, betattbarele);
  //  migmatrixcombined->Add(migmatrixmu, betattbarmu);
  
  H2rec recoplot(xedgesrec, yedgesrec);
  recoplot.add(&recoplotele, betattbarele);
  //  recoplot.add(&recoplotmu, betattbarmu);
  TH2F* recovis = recoplot.toVisTH2F();
  cmsstyle->setup_style_2D(recovis, labelOfXAxisVar, labelOfSensVar);
  //recovis->SetMinimum(0);
  
  
  // .. set overflow bins correctly as well
  TH2F* nonselected = (TH2F*)selefffile->Get("nonselected_reweighted");
  TH1F nonselected1d("nonselected1d","nonselected1d", nbinsafter, 0.5, nbinsafter+0.5); 

  unwrap2dhisto(nonselected, &nonselected1d);


  for(int i=0; i<nbinsafter; i++)
  {
    double content = nonselected1d.GetBinContent(i+1);
    double error = nonselected1d.GetBinError(i+1);  
    migmatrixcombined->SetBinContent(i+1, 0, content);
    migmatrixcombined->SetBinError(i+1, 0, error);
  }
  
  double temp_margin = canv1->GetRightMargin();
  
  /// some plots
  
  TH2F* migmatrixclone = (TH2F*) normalizeMigMat(migmatrixcombined);
  migmatrixclone->SetMinimum(0);
  migmatrixclone->Draw("COLZ");
  drawMigMatBinEdges();
  saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring+"migmatrix_"+ nameOfXAxisVarShort));
  canv1->SetRightMargin(temp_margin);
  
  migmatrixclone->GetXaxis()->SetNdivisions(0,0);
  migmatrixclone->GetYaxis()->SetNdivisions(0,0);
  migmatrixclone->GetXaxis()->SetLabelSize(0);
  migmatrixclone->GetYaxis()->SetLabelSize(0);
  migmatrixclone->SetXTitle("");
  migmatrixclone->SetYTitle("");
  migmatrixclone->Draw("COLZ");
  drawMigMatBinEdges();
  saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring+"migmatrix_nolabels"));
  
  for(int i=0; i<nbinsmafter; i++)
  {
    cmsstyle->setup_style_2D(visuplots[i], "gen", "rec");
    visuplots[i]->Draw("COLZ");
    saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring + TString::Format("recoplotbin%d",i) ));
  }
  
  recovis->Draw("COLZ");
  //recovis->GetZaxis()->SetRangeUser(80.,100.);
  drawRecBinEdges();
  saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring+"reconstructed"));
  
  
  cmsstyle->setup_style(sensrecoplotinc, labelOfSensVar, "N", 1, 0, 0);
  sensrecoplotinc->Draw("HIST");
  saveAs(gPad, (TString("output")+fileSuffix+"/"+ dimstring + "sensrecoplotinc" ));
  sensrecoplotinc->SetMaximum(sensrecoplotinc->GetMaximum()*1.2);
  sensrecoplotinc->SetMinimum(0);
  draw1dRecBinEdges(sensrecoplotinc, 2);
  saveAs(gPad, (TString("output")+fileSuffix+"/"+ dimstring + "sensrecoplotinc_withlines" ));
  
  
  /// write output
  
  
//  TFile* file = new TFile(TString("output/")+fileSuffix+"/"+dimstring+"migmatrix.root","recreate");
  TFile* file = new TFile(TString("output/")+dimstring+"migmatrix.root","recreate");
  file->cd();
  
  migmatrixele->Write();
  //  migmatrixmu->Write();
  migmatrixcombined->Write();
  
//  smudges->Write();
//  msmudge->Write();
//  etasmudge->Write();

  
  file->Close();
  
  return 0;
}


int main()
{
  return migrationmatrix();
}
