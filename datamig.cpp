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


#include "cmsstyle.hpp"
#include "samplelocations.hpp"
#include "fitresults.hpp"
#include "helpers.hpp"

#include "binning.hpp"
#include "specialhelpers.hpp"

// WARNING make sure variables in binning.hpp are set correctly


#include "myunfold_class.hpp"
#include "myunfold_class1d.hpp"

const int numexperiments = 5000; // should be about 50,000

int datamig()
{
  TH1::SetDefaultSumw2(true);

  
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
  TFile* preselfile = new TFile("output/seleff.root");
  TH2F* hpresel = (TH2F*)preselfile->Get("preN");
  TH1F* hpreselunwrapped = new TH1F("hpreselunwrapped","hpreselunwrapped", nbinsafter, 0.5, nbinsafter+0.5); ;
  unwrap2dhisto(hpresel, hpreselunwrapped);
  cmsstyle->setup_style_2D(hpresel, labelOfXAxisVar, labelOfSensVar);
  
  // same for "1d" case where there's only 1 m-bin
  TFile* preselfile1d = new TFile("output/1dseleff.root");
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
 
  
  /// init all kinds of histograms for our results
  
  TH1F* asys = new TH1F("asys", "asys", 1000, -1., 1.);
  TH1F* asys1d = new TH1F("asys1d", "asys1d", 1000, -1., 1.);
 
  TH1F* asydiff = new TH1F("asydiff", "asydiff", 100, -0.1, 0.1);
  TH1F* asydiff1d = new TH1F("asydiff1d", "asydiff1d", 100, -0.1, 0.1);
  
  TH1F* asyerrors = new TH1F("asyerrors", "asyerrors", 100, 0.0, 0.02);
  TH1F* asyerrors1d = new TH1F("asyerrors1d", "asyerrors1d", 100, 0.00, 0.02);
 
  
  // histos for asymmetries in the different mass bins
  TH1F* diffasys[nbinsmafter];
  for(int i=0; i<nbinsmafter; i++)
  {
    char name[20];
    sprintf(name, "diffasybin%d", i);
    diffasys[i] = new TH1F(name, name, 1000, -1., 1.);
  }
  
  TH1F* diffasydiffs[nbinsmafter];
  for(int i=0; i<nbinsmafter; i++)
  {
    char name[20];
    sprintf(name, "diffasydiffbin%d", i);
    diffasydiffs[i] = new TH1F(name, name, 100, -0.1, 0.1);
  }
  
  // init array with true presel bin contents for comparison
  double truecontents[nbinsafter];
  for(int i=0; i<nbinsafter; i++)
  {
    truecontents[i] = hpreselunwrapped->GetBinContent(i+1) * expected_ttbar_presel / npresel; 
  }
  
  TH1::AddDirectory(kFALSE);
  
  /// create data histo
  TH2F* hrec;
  TH2F* hrec1d;  
  H2rec myh2rec("hrec", xedgesrec, yedgesrec);
  H2rec1d myh2rec1d("hrec1d", xedgesrec1d, yedgesrec1d);  
  TFile* inputfile = new TFile(sl_data_comb);
  
  /// read out tuples
  //  TNtuple* tuple = (TNtuple*)inputfile->Get("TTbarSelHyppsiprime/BBB_UnfoldingTuple");
  TNtuple* tuple = (TNtuple*)inputfile->Get("Tuple");
  
  float mttbargen, mttbarrec, diffabsetagen, diffabsetarec, weight;
  
  tuple->SetBranchAddress(nameOfXAxisVarGen, &mttbargen);
  tuple->SetBranchAddress(nameOfXAxisVarRec, &mttbarrec);
  tuple->SetBranchAddress(nameOfSensVarGen, &diffabsetagen);
  tuple->SetBranchAddress(nameOfSensVarRec, &diffabsetarec);
  tuple->SetBranchAddress("weight", &weight);
  
  /// fill histos containing actual data
  const int nentries = tuple->GetEntries();
  
  for(int i=0; i<nentries; i++)
  {
    tuple->GetEntry(i);
    myh2rec.fill(fabs(mttbarrec), diffabsetarec, weight);
    myh2rec1d.fill(fabs(mttbarrec), diffabsetarec, weight);
  }
  
  hrec = myh2rec.toTH2F();
  hrec1d = myh2rec1d.toTH2F();
  
  vary_migmatrix_for_each_unfolding = true;
  
  /// perform actual experiments
  for(int n=0; n<numexperiments; n++)
  {
    if(n%100 == 0)
      cout << "Experiment number " << n << endl;
      
    // unfold it
    TH2F* hunfold = myunfold.unfoldHisto(hrec, true); // we'll treat the pseudo experiments as if they were real data.
        
    hrec1d->Sumw2(); // recalc errors, hopefully this works
    
    TH2F* hunfold1d = myunfold1d.unfoldHisto(hrec1d, true);
    
    // calc asys and other values needed for histos
    const double pos = hunfold->Integral(1, nbinsmafter, nbinsetaafter/2+1, nbinsetaafter);
    const double neg = hunfold->Integral(1, nbinsmafter, 1, nbinsetaafter/2);
    double asy = (pos-neg)/(pos+neg);
    

    const double asy_error = asymmetryerror_afterunfolding_2d(myunfold.errormatrix, nbinsafter, pos, neg, nbinsetaafter);
    
    asys->Fill(asy);
    asyerrors->Fill(asy_error);
    asydiff->Fill(asy-trueasy);
    
    
    
    for(int i=0; i<nbinsmafter; i++)
    {
      const double pos = hunfold->Integral(i+1, i+1, nbinsetaafter/2+1, nbinsetaafter);
      const double neg = hunfold->Integral(i+1, i+1, 1, nbinsetaafter/2);
      double asy = (pos-neg)/(pos+neg);
      

      diffasys[i]->Fill(asy);
      
      const double truepos = hpresel->Integral(i+1, i+1, nbinsetaafter/2+1, nbinsetaafter);
      const double trueneg = hpresel->Integral(i+1, i+1, 1, nbinsetaafter/2);
      const double trueasy = (truepos-trueneg)/(truepos+trueneg);
      diffasydiffs[i]->Fill(asy-trueasy);
      
      const double differr = asymmetryerror_afterunfolding_2d_onexbin(myunfold.errormatrix, nbinsafter, pos, neg, nbinsetaafter, i); 
    }
    
    
    const double pos1d = hunfold1d->Integral(1, nbinsmafter1d, nbinsetaafter1d/2+1, nbinsetaafter1d);
    const double neg1d = hunfold1d->Integral(1, nbinsmafter1d, 1, nbinsetaafter1d/2);
    double asy1d = (pos1d-neg1d)/(pos1d+neg1d);
    const double asy_error1d = asymmetryerror_afterunfolding_1d(myunfold1d.errormatrix, myunfold1d.hunfoldunwrapped);
    

    
    asys1d->Fill(asy1d);
    asyerrors1d->Fill(asy_error1d);
    asydiff1d->Fill(asy1d-trueasy);
    
  }

  TH1::AddDirectory(kTRUE);

  
  
  TString outfilename = TString("output/datamig.root");
  
  
  TFile out(outfilename,"recreate");
  
  asys->Write();
  asyerrors->Write();
  asys1d->Write();
  asyerrors1d->Write();
  asydiff->Write();
  asydiff1d->Write();

  
 
  for(int i=0; i<nbinsmafter; i++)
  {
    diffasys[i]->Write();
    diffasydiffs[i]->Write();
  }
  
  
  cout << "Finished successfully" << endl;

  return 0;
}



int main()
{
  return datamig();
}
