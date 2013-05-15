#include <iostream>
#include <assert.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
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

#include <sys/stat.h> // for mkdir
#include <sys/types.h>


#include "cmsstyle.hpp"
#include "fitresults.hpp"
#include "helpers.hpp"

#include "binning.hpp"
#include "specialhelpers.hpp"

// WARNING make sure variables in binning.hpp are set correctly


#include "myunfold_class.hpp"
#include "myunfold_class1d.hpp"

#ifndef FINDTAU
  const int numexperiments = 5000; // should be about 50,000: to run pseudo experiments for linearity plots
#else
  const int numexperiments = 100; // to find tau
#endif

int pseudoexp()
{
  TH1::SetDefaultSumw2(true);

//  if(correct_for_lincheck)
//    cout << "WARNING: correcting for lincheck results!" << endl << endl;
//  else
//    cout << "WARNING: NOT correcting for lincheck results!" << endl << endl;
  

  TRandom3* random = new TRandom3();
  random->SetSeed();
  
  //style stuff
  CMSStyle* cmsstyle = new CMSStyle();
  //gStyle->SetOptStat("");
  TCanvas* canv1 = new TCanvas("canv1","canv1",300,0,800,600);
  canv1->SetBatch();
  canv1->SetRightMargin(0.04); // this would be cms recommendation, but we need...
  canv1->SetRightMargin(0.14); //  this - to make z-axis stuff show completely in 2d plots 
  
  // get true generated values from before selection
//  TFile* preselfile = new TFile(TString("output/")+fileSuffix+"/seleff.root");
  TFile* preselfile = new TFile(TString("output/")+"seleff.root");
  TH2F* hpresel = (TH2F*)preselfile->Get("preN");
  TH1F* hpreselunwrapped = new TH1F("hpreselunwrapped","hpreselunwrapped", nbinsafter, 0.5, nbinsafter+0.5); ;
  unwrap2dhisto(hpresel, hpreselunwrapped);
  cmsstyle->setup_style_2D(hpresel, labelOfXAxisVar, labelOfSensVar);
  
  // same for "1d" case where there's only 1 m-bin
//  TFile* preselfile1d = new TFile(TString("output/")+fileSuffix+"/1dseleff.root");
  TFile* preselfile1d = new TFile(TString("output/")+"1dseleff.root");
  TH2F* hpresel1d = (TH2F*)preselfile1d->Get("preN");
  TH1F* hpreselunwrapped1d = new TH1F("hpreselunwrapped1d","hpreselunwrapped1d", nbinsafter1d, 0.5, nbinsafter1d+0.5); ;
  unwrap2dhisto(hpresel1d, hpreselunwrapped1d);
  cmsstyle->setup_style_2D(hpresel1d, labelOfXAxisVar, labelOfSensVar);
  
  MyUnfold myunfold;
  MyUnfold1d myunfold1d;
    
  // calculate value for MC to be normalized to
  // (basically, the fit result of ttbar divided by selection efficiency)
  
  const double npresel = hpresel->Integral();
  
  TH2F* hnonselected = (TH2F*)preselfile->Get("nonselected_reweighted");
  const double nmigmatrixselected = myunfold.migmatrix->Integral();
  const double nmigmatrix = nmigmatrixselected + hnonselected->Integral();
  const double totalseleff = nmigmatrixselected / nmigmatrix;
  const double expected_ttbar_presel = (nttbarele /*+ nttbarmu*/) / totalseleff;
  
  
  // calculate "true" asymmetry of our sample
  const double truepos = hpresel->Integral(1, nbinsmafter, nbinsetaafter/2+1, nbinsetaafter);
  const double trueneg = hpresel->Integral(1, nbinsmafter, 1, nbinsetaafter/2);
  const double trueasy = (truepos-trueneg)/(truepos+trueneg);
  
  // ...and also true asymmetries for the individual mass bins
  double difftrueasys[nbinsmafter];
  for(int i=0; i<nbinsmafter; i++)
  {
    const double pos = hpresel->Integral(i+1, i+1, nbinsetaafter/2+1, nbinsetaafter);
    const double neg = hpresel->Integral(i+1, i+1, 1, nbinsetaafter/2);
    difftrueasys[i] = (pos-neg)/(pos+neg);
  }
 
  
  /// init all kinds of histograms for our results
  
  TH1F* asys = new TH1F("asys", "asys", 500, -1., 1.);
  TH1F* asys1d = new TH1F("asys1d", "asys1d", 500, -1., 1.);
  
  TH1F* asypull = new TH1F("asypull", "asypull", 500, -10.0, 10.0);
  TH1F* asypull1d = new TH1F("asypull1d", "asypull1d", 500, -10.0, 10.0);
  
  TH1F* asydiff = new TH1F("asydiff", "asydiff", 100, -0.1, 0.1);
  TH1F* asydiff1d = new TH1F("asydiff1d", "asydiff1d", 100, -0.1, 0.1);
    
  TH1F* ndiff = new TH1F("ndiff", "ndiff", 100, -1200000, 1200000);
  TH1F* ndiff1d = new TH1F("ndiff1d", "ndiff1d", 100, -1200000, 1200000);
  
  TH1F* nreldiff = new TH1F("nreldiff", "nreldiff", 100, -3, 3);
  TH1F* nreldiff1d = new TH1F("nreldiff1d", "nreldiff1d", 100, -3, 3);
  
  TH1F* log10taus = new TH1F("log10taus", "log10taus", 1000, -6, -1.5);
  TH1F* log10taus1d = new TH1F("log10taus1d", "log10taus1d", 1000, -6, -1.5);

  
  
//  TH1F* asyerrors = new TH1F("asyerrors", "asyerrors", 100, 0.008, 0.014);
//  TH1F* asyerrors1d = new TH1F("asyerrors1d", "asyerrors1d", 100, 0.008, 0.014);

  TH1F* asyerrors = new TH1F("asyerrors", "asyerrors", 1000, 0., 0.5);
  TH1F* asyerrors1d = new TH1F("asyerrors1d", "asyerrors1d", 1000, 0., 0.5);

  
  TH2F* asy1dvs2d = new TH2F("asy1dvs2d", "asy1dvs2d", 50, -0.1, 0.1, 50, -0.1, 0.1);
  TH1F* asydifference = new TH1F("asydifference", "asydifference", 100, -0.25, 0.25);  // differences for 1d vs 2d unfolding
  
  // init histos for relative differences bin-by-bin
  TH1F* reldiffs[nbinsafter];
  TH1F* pulls[nbinsafter];
  for(int i=0; i<nbinsafter; i++)
  {
    char name[20];
    sprintf(name, "reldiff%d", i);
    reldiffs[i] = new TH1F(name, name, 100, -0.5, 0.5);
    sprintf(name, "pull%d", i);
    pulls[i] = new TH1F(name, name, 100, -5, 5);
  }
  
  TH1F* reldiffs1d[nbinsafter1d];
  TH1F* pulls1d[nbinsafter1d];
  for(int i=0; i<nbinsafter1d; i++)
  {
    char name[20];
    sprintf(name, "1dreldiff%d", i);
    reldiffs1d[i] = new TH1F(name, name, 100, -0.5, 0.5);
    sprintf(name, "1dpull%d", i);
    pulls1d[i] = new TH1F(name, name, 100, -5, 5);
  }
  
  // histos for asymmetries in the different mass bins
  TH1F* diffasys[nbinsmafter];
  for(int i=0; i<nbinsmafter; i++)
  {
    char name[20];
    sprintf(name, "diffasybin%d", i);
    diffasys[i] = new TH1F(name, name, 100, -0.35, 0.35);
  }
  
  TH1F* diffasydiffs[nbinsmafter];
  TH1F* diffasypulls[nbinsmafter];
  TH1F* diffasyerrors[nbinsmafter];
  for(int i=0; i<nbinsmafter; i++)
  {
    char name[20];
    sprintf(name, "diffasydiffbin%d", i);
    diffasydiffs[i] = new TH1F(name, name, 100, -0.35, 0.35);
    sprintf(name, "diffasypullbin%d", i);
    diffasypulls[i] = new TH1F(name, name, 500, -10., 10.);
    
    sprintf(name, "diffasyerrorbin%d", i);
    double lowedge = 0.;
    double highedge = 0.5;
/*
    if(!TString(nameOfXAxisVarShort).CompareTo("m"))
    {
      if(i==0) {lowedge=0.025; highedge=0.035;}
      else if(i==1) {lowedge=0.01; highedge=0.02;}
      else {lowedge=0.015; highedge=0.025;}
    }
    else if(!TString(nameOfXAxisVarShort).CompareTo("pt"))
    {
      if(i==0) {lowedge=0.02; highedge=0.03;}
      else if(i==1) {lowedge=0.01; highedge=0.02;}
      else {lowedge=0.015; highedge=0.025;}
    }
    else if(!TString(nameOfXAxisVarShort).CompareTo("y"))
    {
      if(i==0) {lowedge=0.015; highedge=0.025;}
      else if(i==1) {lowedge=0.01; highedge=0.02;}
      else {lowedge=0.015; highedge=0.025;}
    }
*/
    diffasyerrors[i] = new TH1F(name, name, 1000, lowedge, highedge);
  }
  
  // init array with true presel bin contents for comparison
  double truecontents[nbinsafter];
  for(int i=0; i<nbinsafter; i++)
  {
    truecontents[i] = hpreselunwrapped->GetBinContent(i+1) * expected_ttbar_presel / npresel; 
  }
  
  double truecontents1d[nbinsafter1d];
  for(int i=0; i<nbinsafter1d; i++)
  {
    truecontents1d[i] = hpreselunwrapped1d->GetBinContent(i+1) * expected_ttbar_presel / npresel; 
  }
  

  
  TH1::AddDirectory(kFALSE);
  
  /// perform actual experiments
  for(int n=0; n<numexperiments; n++)
  {
    if(n%1000 == 0)
      cout << "Experiment number " << n << endl;
      
    // create pseudo data histo
    TH2F hrec("hrec","hrec", nbinsm, 0.5, nbinsm+0.5, nbinseta, 0.5, nbinseta+0.5);
    TH2F hrec1d("hrec1d","hrec1d", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d); 
    
    // loop over backgrounds and signal to compose pseudo-data-sample
    for(int i=0; i<4; i++)
    {
      TH2F* h = 0; 
      TH2F* h1d = 0; 
      double nfit = 0;
      double fiterror = 0;
      
      
      
      switch (i)
      {
        case 0:  h=myunfold.httbar_ele; h1d=myunfold1d.httbar_ele; nfit=nttbarele; fiterror=0; break; 
	  //        case 1:  h=myunfold.httbar_mu; h1d=myunfold1d.httbar_mu; nfit=nttbarmu; fiterror=0; break;
        
        default: const int tempi = i-1;
                 h=myunfold.bgwrappedlist[tempi];
                 h1d=myunfold1d.bgwrappedlist[tempi];
                 nfit=myunfold.bgintegrals[tempi];
                 fiterror=myunfold.bgabserrors[tempi];
                 break;
        
        //         case 2:  h=myunfold.hwjets_neg_ele; h1d=myunfold1d.hwjets_neg_ele; nfit=nwjetsnegele; fiterror=nwjetsnegeleerror; break;
        //         case 3:  h=myunfold.hzjets_ele; h1d=myunfold1d.hzjets_ele; nfit=nzjetsele; fiterror=nzjetseleerror; break;
        //         case 4:  h=myunfold.hstt_ele; h1d=myunfold1d.hstt_ele; nfit=nstele; fiterror=nsteleerror; break; //it's named stt, but it's whole single t
        //         case 5:  h=myunfold.hqcd_ele; h1d=myunfold1d.hqcd_ele; nfit=nqcdele; fiterror=nqcdeleerror; break;
        //         
        //         case 6:  h=myunfold.httbar_mu; h1d=myunfold1d.httbar_mu; nfit=nttbarmu; fiterror=0; break;
        //         case 7:  h=myunfold.hwjets_pos_mu; h1d=myunfold1d.hwjets_pos_mu; nfit=nwjetsposmu; fiterror=nwjetsposmuerror; break;
        //         case 8:  h=myunfold.hwjets_neg_mu; h1d=myunfold1d.hwjets_neg_mu; nfit=nwjetsnegmu; fiterror=nwjetsnegmuerror; break;
        //         case 9:  h=myunfold.hzjets_mu; h1d=myunfold1d.hzjets_mu; nfit=nzjetsmu; fiterror=nzjetsmuerror; break;
        //         case 10: h=myunfold.hstt_mu; h1d=myunfold1d.hstt_mu; nfit=nstmu; fiterror=nstmuerror; break; //it's named stt, but it's whole single t
        //         case 11: h=myunfold.hqcd_mu; h1d=myunfold1d.hqcd_mu; nfit=nqcdmu; fiterror=nqcdmuerror; break;
      }
      
      // WARNING assumes fiterror scales with lumi - pessimistic
      const int ndraw = random->Poisson(drawfactor * ( fiterror != 0 ? random->Gaus(nfit, fiterror) : nfit ));
            
      // vary all distributions using bin-by-bin errors.
      // this accounts for the limited MC statistics
      {
        TH2F* hclone = (TH2F*) h->Clone();
        TH2F* h1dclone = (TH2F*) h1d->Clone();
        
        const int nbinsx = hclone->GetNbinsX();
        const int nbinsx1d = h1dclone->GetNbinsX();
        const int nbinsy = hclone->GetNbinsY();
        const int nbinsy1d = h1dclone->GetNbinsY();
        
        for(int x=0; x<nbinsx+1; x++)
        {
          for(int y=0; y<nbinsy+1; y++)
          {
            const double val = hclone->GetBinContent(x, y);
            const double err = hclone->GetBinError(x, y);
            hclone->SetBinContent(x, y, random->Gaus(val, err));
          }
        }
        
        for(int x=0; x<nbinsx1d+1; x++)
        {
          for(int y=0; y<nbinsy1d+1; y++)
          {
            const double val = h1dclone->GetBinContent(x, y);
            const double err = h1dclone->GetBinError(x, y);
            h1dclone->SetBinContent(x, y, random->Gaus(val, err));
          }
        }
        
        h = hclone;
        h1d = h1dclone;
      }
      
      Double_t mass, sensvar;
      
      for(int j=0; j<ndraw; j++)
      {
        h->GetRandom2(mass, sensvar); // writes into the two variables!
        fill_nooverflow_2d(&hrec, mass, sensvar, 1);
        
        h1d->GetRandom2(mass, sensvar); // writes into the two variables!
        fill_nooverflow_2d(&hrec1d, mass, sensvar, 1);
      }
      
    }
      
    // unfold it



    TH2F* hunfold = myunfold.unfoldHisto(&hrec, true, drawfactor); // we'll treat the pseudo experiments as if they were real data.
    log10taus->Fill(TMath::Log10(myunfold.lastTau));
    
    //hrec1d.Sumw2(); // recalc errors, hopefully this works

    
    TH2F* hunfold1d = myunfold1d.unfoldHisto(&hrec1d, true, drawfactor);
    log10taus1d->Fill(TMath::Log10(myunfold1d.lastTau));


    
    // calc asys and other values needed for histos
    const double pos = hunfold->Integral(1, nbinsmafter, nbinsetaafter/2+1, nbinsetaafter);
    const double neg = hunfold->Integral(1, nbinsmafter, 1, nbinsetaafter/2);
    double asy = (pos-neg)/(pos+neg);
/*    
    if(correct_for_lincheck)
    {
      asy -= lincheck_offset; // correct for offset of linearity check
      asy /= lincheck_slope;  // correct for linearity check not yielding slope of 1
    }
*/    
    const double asy_error = asymmetryerror_afterunfolding_2d(myunfold.errormatrix, nbinsafter, pos, neg, nbinsetaafter);
    
//    cout << "asy_error: " << asy_error << endl;
//    cout << "asydiff: " <<  asy-trueasy << endl;


    asys->Fill(asy);
    asyerrors->Fill(asy_error);
    asypull->Fill((asy-trueasy)/asy_error);
    asydiff->Fill(asy-trueasy);
    
    for(int i=0; i<nbinsafter; i++)
    {
      const double reldiff = ( myunfold.hunfoldunwrapped->GetBinContent(i+1) - truecontents[i] * drawfactor ) / ( truecontents[i] *drawfactor );
      reldiffs[i]->Fill(reldiff);
      const double pull = ( myunfold.hunfoldunwrapped->GetBinContent(i+1) - truecontents[i] * drawfactor ) / myunfold.hunfoldunwrapped->GetBinError(i+1);
      pulls[i]->Fill(pull);
    }
    
    for(int i=0; i<nbinsafter1d; i++)
    {
      const double reldiff = ( myunfold1d.hunfoldunwrapped->GetBinContent(i+1) - truecontents1d[i] * drawfactor ) / ( truecontents1d[i] *drawfactor );
      reldiffs1d[i]->Fill(reldiff);
      const double pull = ( myunfold1d.hunfoldunwrapped->GetBinContent(i+1) - truecontents1d[i] * drawfactor ) / myunfold1d.hunfoldunwrapped->GetBinError(i+1);
      pulls1d[i]->Fill(pull);
    }
    
    
    for(int i=0; i<nbinsmafter; i++)
    {
      const double pos = hunfold->Integral(i+1, i+1, nbinsetaafter/2+1, nbinsetaafter);
      const double neg = hunfold->Integral(i+1, i+1, 1, nbinsetaafter/2);
      double asy = (pos-neg)/(pos+neg);
/*      
      if(correct_for_lincheck)
      {
        asy -= lincheck_diffoffset[i]; // correct for linearity check having offset
        asy /= lincheck_diffslope[i]; // correct for linearity check not yielding slope of 1
      }
*/      
      diffasys[i]->Fill(asy);
      
      diffasydiffs[i]->Fill(asy-difftrueasys[i]);
      
      const double differr = asymmetryerror_afterunfolding_2d_onexbin(myunfold.errormatrix, nbinsafter, pos, neg, nbinsetaafter, i); 
      diffasypulls[i]->Fill((asy-difftrueasys[i])/differr);
      
      diffasyerrors[i]->Fill(differr);
    }
    
    
    const double pos1d = hunfold1d->Integral(1, nbinsmafter1d, nbinsetaafter1d/2+1, nbinsetaafter1d);
    const double neg1d = hunfold1d->Integral(1, nbinsmafter1d, 1, nbinsetaafter1d/2);
    double asy1d = (pos1d-neg1d)/(pos1d+neg1d);
    const double asy_error1d = asymmetryerror_afterunfolding_1d(myunfold1d.errormatrix, myunfold1d.hunfoldunwrapped);
/*    
    if(correct_for_lincheck)
    {
      asy1d -= lincheck_offset1d; // correct for offset of linearity check
      asy1d /= lincheck_slope1d;  // correct for linearity check not yielding slope of 1
    }
*/    

//    cout << "asy_error1d: " << asy_error1d << endl;
//    cout << "asydiff1d: " <<  asy1d-trueasy << endl;



    asys1d->Fill(asy1d);
    asyerrors1d->Fill(asy_error1d);
    asydifference->Fill(asy - asy1d);
    asy1dvs2d->Fill(asy1d, asy);
    asypull1d->Fill((asy1d-trueasy)/asy_error1d);
    asydiff1d->Fill(asy1d-trueasy);
    
    ndiff->Fill(pos+neg-myunfold.expectedNUnfolded);
    ndiff1d->Fill(pos1d+neg1d-myunfold1d.expectedNUnfolded);
    
    nreldiff->Fill((pos+neg-myunfold.expectedNUnfolded)/myunfold.expectedNUnfolded);
    nreldiff1d->Fill((pos1d+neg1d-myunfold1d.expectedNUnfolded)/myunfold1d.expectedNUnfolded);
  }

  TH1::AddDirectory(kTRUE);

  /*
  canv1->Divide(2);
  canv1->cd(1);
  asys->Draw("HIST");
  canv1->cd(2);
  asyerrors->Draw("HIST");
  */
  
  #ifndef FINDTAU
    TString outfilename = TString("output/")+fileSuffix+"/"+"pseudo.root";
//    if(correct_for_lincheck)
//      outfilename = TString("output/")+fileSuffix+"/"+"pseudo_corrected.root";
  #else
    TString outfilename = TString("output/")+fileSuffix+"/"+"pseudo_tau.root";
  #endif
  
  TFile out(outfilename,"recreate");
  
  asys->Write();
  asyerrors->Write();
  asys1d->Write();
  asyerrors1d->Write();
  asypull->Write();
  asypull1d->Write();
  asydiff->Write();
  asydiff1d->Write();
  log10taus->Write();
  log10taus1d->Write();
  cout << " log10tau mean: " << log10taus->GetMean() << endl;
  cout << " log10tau1d mean: " << log10taus1d->GetMean() << endl;
  ndiff->Write();
  ndiff1d->Write();
  nreldiff->Write();
  nreldiff1d->Write();
  
  // text-based output of reldiff offsets
  TString reweightname = "noreweight";
  mkdir(TString("output/")+fileSuffix+"/reweightings/",0777);
  ofstream resfile(TString("output/")+fileSuffix+"/reweightings/"+reweightname+".txt");
  resfile << reweightname << " ";
  resfile << asydiff->GetMean() << " " << asydiff1d->GetMean() << " ";
  for(int i=0; i<nbinsmafter; i++)
  {
    resfile << diffasydiffs[i]->GetMean() << " ";
  }
  resfile << endl;
  resfile.close();
  
  
  // from here on output of pulls and reldiffs
  mkdir(TString("output/")+fileSuffix+"/pulls_and_reldiffs/",0777);
  
  for(int i=0; i<nbinsafter; i++)
  {
    reldiffs[i]->Write();
    reldiffs[i]->Draw("HIST");
    saveAs(canv1, TString("output/")+fileSuffix+"/pulls_and_reldiffs/"+reldiffs[i]->GetName());
    
    pulls[i]->Fit("gaus","");
    pulls[i]->Draw("HIST");
    saveAs(canv1, TString("output/")+fileSuffix+"/pulls_and_reldiffs/"+pulls[i]->GetName());
    pulls[i]->Write();
  }

  for(int i=0; i<nbinsafter1d; i++)
  {
    reldiffs1d[i]->Write();
    reldiffs1d[i]->Draw("HIST");
    saveAs(canv1, TString("output/")+fileSuffix+"/pulls_and_reldiffs/"+reldiffs1d[i]->GetName());
    
    pulls1d[i]->Fit("gaus","");
    pulls1d[i]->Draw("HIST");
    saveAs(canv1, TString("output/")+fileSuffix+"/pulls_and_reldiffs/"+pulls1d[i]->GetName());
    pulls1d[i]->Write();
  }
  
  double reldiff_deviations = 0;
  double pull_deviations = 0;
  cout << "\n\n---- pulls and reldiffs ----\n";
  for(int i=0; i<nbinsafter1d; i++)
  {
    cout << "reldiff1d[" << i << "]  = " << reldiffs1d[i]->GetMean() << endl;
    cout << "       [" << i << "] +- " << reldiffs1d[i]->GetMeanError() << endl;
    cout << "pull1d[" << i << "]  = " << pulls1d[i]->GetRMS() << endl;
    cout << "    [" << i << "] +- " << pulls1d[i]->GetRMSError() << endl;
    
    reldiff_deviations += fabs( reldiffs1d[i]->GetMean() / reldiffs1d[i]->GetMeanError() );
    pull_deviations += fabs( (pulls1d[i]->GetRMS() -1) / pulls1d[i]->GetRMSError() ) ;
  }
  cout << "Mean reldiff deviation (sigma): " << reldiff_deviations/nbinsafter << endl;
  cout << "Mean pull deviation (sigma): " << pull_deviations/nbinsafter << endl;
  cout << "\n----------------------------\n";
  
  for(int i=0; i<nbinsmafter; i++)
  {
    diffasys[i]->Write();
    diffasydiffs[i]->Write();
    diffasypulls[i]->Write();
    diffasyerrors[i]->Write();
  }
  
  asydifference->Write();
  asy1dvs2d->Write();

  TNtuple* asytuple = new TNtuple("trueasys", "trueasys", "trueasy");
  asytuple->Fill(trueasy);
  for(int i=0; i<nbinsafter; i++)
  {
    asytuple->Fill(difftrueasys[i]);
  }
  asytuple->Write();
  
  
  out.Close();
  
  cout << "Finished successfully" << endl;

  return 0;
}



int main()
{
  return pseudoexp();
}
