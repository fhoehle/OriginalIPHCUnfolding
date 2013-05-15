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
#include "fitresults.hpp"
#include "helpers.hpp"
#include "samplelocations.hpp"


#include "binning.hpp"

// WARNING make sure variables in binning.hpp are set correctly


#include "myunfold_class.hpp"
#include "myunfold_class1d.hpp"


int lincheck()
{
  TH1::SetDefaultSumw2(true);

  const int numexperiments = 600; // experiments per reweighting factor was 600
  const int onlygenbin = -1; // -1 to use all bins
  
  #ifdef FINDTAU 
    cout << "ERR: FINDTAU still enabled" << endl; // trigger compiler error, cause we probably didn't want to run like this
  #endif
  
  TRandom3* random = new TRandom3();
  random->SetSeed();
  
  //style stuff
  CMSStyle* cmsstyle = new CMSStyle();
  //gStyle->SetOptStat("");
  TCanvas* canv1 = new TCanvas("canv1","canv1",300,0,800,600);
  canv1->SetBatch();
  canv1->SetRightMargin(0.04); // this would be cms recommendation, but we need...
  canv1->SetRightMargin(0.14); //  this - to make z-axis stuff show completely in 2d plots 
  
  // for true generated values from before selection
  // (also reweighted distributions)
//  TFile* preselfile = new TFile(TString("output/")+fileSuffix+"/seleff.root");
//  TFile* preselfile1d = new TFile(TString("output/")+fileSuffix+"/1dseleff.root");
  TFile* preselfile = new TFile(TString("output/")+"seleff.root");
  TFile* preselfile1d = new TFile(TString("output/")+"1dseleff.root");

  
  MyUnfold myunfold;
  MyUnfold1d myunfold1d;
  
  ofstream rewfile(TString(TString("output/")+fileSuffix+"/lincheck_reweightedasys.txt"));
  rewfile << "sample\ttrue asy\tunf asy\n";
  
  
  /// init all kinds of histograms we need as input
  
  #ifdef MADGRAPH
//    TFile* file_ttbar_ele = new TFile("/portal/ekpams2/home/froscher/outputroot/2dunfolding/ttbar_ele_histos.root");
    TFile* file_ttbar_ele = new TFile("/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_TTbarAllMCatNLO_PileUpReweighted.root");
//    TFile* file_ttbar_ele = new TFile("/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_TTbarSignalTuneZ2_7TeV-madgraph-tauola_PileUpReweighted.root");
    //    TFile* file_ttbar_mu = new TFile("/portal/ekpams2/home/froscher/outputroot/2dunfolding/ttbar_mu_histos.root");
  #endif
  #ifdef POWHEG
    TFile* file_ttbar_ele;
    //    TFile* file_ttbar_mu;
    if(!random_split)
    {
      file_ttbar_ele = new TFile(sl_ttbarele);
      //      file_ttbar_mu = new TFile(sl_ttbarmu);
    }
    else
    {
//      file_ttbar_ele = new TFile("FIXME/portal/ekpams2/home/froscher/outputroot/powhegsample/ttbar_ele_combined_subsample2_histos.root");
      file_ttbar_ele = new TFile("/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_TTbarAllMCatNLO_PileUpReweighted.root");
//      file_ttbar_ele = new TFile("/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_TTbarSignalTuneZ2_7TeV-madgraph-tauola_PileUpReweighted.root");
      //      file_ttbar_mu = new TFile("FIXME/portal/ekpams2/home/froscher/outputroot/powhegsample/ttbar_mu_combined_subsample2_histos.root");      
    }
      
  #endif
  
//  TNtuple* ttbartuple_ele = (TNtuple*)file_ttbar_ele->Get("TTbarSelHyppsiprime/BBB_UnfoldingTuple");
//  TNtuple* ttbartuple_mu = (TNtuple*)file_ttbar_mu->Get("TTbarSelHyppsiprime/BBB_UnfoldingTuple");
  TNtuple* ttbartuple_ele = (TNtuple*)file_ttbar_ele->Get("Tuple");
  //  TNtuple* ttbartuple_mu = (TNtuple*)file_ttbar_mu->Get("Tuple");
  
  TH2F* hrecttbar_rew0_ele;
  TH2F* hrecttbar_rewp05_ele;
  TH2F* hrecttbar_rewp10_ele;
  TH2F* hrecttbar_rewp15_ele;
  TH2F* hrecttbar_rewp20_ele;
  TH2F* hrecttbar_rewp25_ele;
  TH2F* hrecttbar_rewm05_ele;
  TH2F* hrecttbar_rewm10_ele;
  TH2F* hrecttbar_rewm15_ele;
  TH2F* hrecttbar_rewm20_ele;
  TH2F* hrecttbar_rewm25_ele;
  /*
  TH2F* hrecttbar_rew0_mu;
  TH2F* hrecttbar_rewp05_mu;
  TH2F* hrecttbar_rewp10_mu;
  TH2F* hrecttbar_rewp15_mu;
  TH2F* hrecttbar_rewp20_mu;
  TH2F* hrecttbar_rewp25_mu;
  TH2F* hrecttbar_rewm05_mu;
  TH2F* hrecttbar_rewm10_mu;
  TH2F* hrecttbar_rewm15_mu;
  TH2F* hrecttbar_rewm20_mu;
  TH2F* hrecttbar_rewm25_mu;  
  */
  TH2F* hrecttbar_rew0_ele_1d;
  TH2F* hrecttbar_rewp05_ele_1d;
  TH2F* hrecttbar_rewp10_ele_1d;
  TH2F* hrecttbar_rewp15_ele_1d;
  TH2F* hrecttbar_rewp20_ele_1d;
  TH2F* hrecttbar_rewp25_ele_1d;
  TH2F* hrecttbar_rewm05_ele_1d;
  TH2F* hrecttbar_rewm10_ele_1d;
  TH2F* hrecttbar_rewm15_ele_1d;
  TH2F* hrecttbar_rewm20_ele_1d;
  TH2F* hrecttbar_rewm25_ele_1d;
  /*  
  TH2F* hrecttbar_rew0_mu_1d;
  TH2F* hrecttbar_rewp05_mu_1d;
  TH2F* hrecttbar_rewp10_mu_1d;
  TH2F* hrecttbar_rewp15_mu_1d;
  TH2F* hrecttbar_rewp20_mu_1d;
  TH2F* hrecttbar_rewp25_mu_1d;
  TH2F* hrecttbar_rewm05_mu_1d;
  TH2F* hrecttbar_rewm10_mu_1d;
  TH2F* hrecttbar_rewm15_mu_1d;
  TH2F* hrecttbar_rewm20_mu_1d;
  TH2F* hrecttbar_rewm25_mu_1d;  
  */
  
  H2gen myh2gen("bla", xedgesgen, yedgesgen); // only to check binning
  float mttbarrec, diffabsetarec, mttbargen, diffabsetagen, weight;
  
  for(int i=0; i<11; i++) // yes, this is a little slow - but not bad enough to motivate the effort of fixing it
  {
    cout << "Initializing histo " << i << endl;
    
    TNtuple* t = 0;
    TH2F** h = 0; 
    TH2F** h1d = 0; 
    double rew = 0; // reweighting factor "k"
    TString name = "";


    
    switch (i)
    {
      case 0: t=ttbartuple_ele; h=&hrecttbar_rew0_ele; h1d=&hrecttbar_rew0_ele_1d; name="hrecttbar_rew0_ele"; rew=0; break;
      case 1: t=ttbartuple_ele; h=&hrecttbar_rewp05_ele; h1d=&hrecttbar_rewp05_ele_1d; name="hrecttbar_rewp05_ele"; rew=0.05; break;
      case 2: t=ttbartuple_ele; h=&hrecttbar_rewp10_ele; h1d=&hrecttbar_rewp10_ele_1d; name="hrecttbar_rewp10_ele"; rew=0.10; break;
      case 3: t=ttbartuple_ele; h=&hrecttbar_rewp15_ele; h1d=&hrecttbar_rewp15_ele_1d; name="hrecttbar_rewp15_ele"; rew=0.15; break;
      case 4: t=ttbartuple_ele; h=&hrecttbar_rewp20_ele; h1d=&hrecttbar_rewp20_ele_1d; name="hrecttbar_rewp20_ele"; rew=0.20; break;
      case 5: t=ttbartuple_ele; h=&hrecttbar_rewp25_ele; h1d=&hrecttbar_rewp25_ele_1d; name="hrecttbar_rewp25_ele"; rew=0.25; break;
      case 6: t=ttbartuple_ele; h=&hrecttbar_rewm05_ele; h1d=&hrecttbar_rewm05_ele_1d; name="hrecttbar_rewm05_ele"; rew=-0.05; break;
      case 7: t=ttbartuple_ele; h=&hrecttbar_rewm10_ele; h1d=&hrecttbar_rewm10_ele_1d; name="hrecttbar_rewm10_ele"; rew=-0.10; break;
      case 8: t=ttbartuple_ele; h=&hrecttbar_rewm15_ele; h1d=&hrecttbar_rewm15_ele_1d; name="hrecttbar_rewm15_ele"; rew=-0.15; break;
      case 9: t=ttbartuple_ele; h=&hrecttbar_rewm20_ele; h1d=&hrecttbar_rewm20_ele_1d; name="hrecttbar_rewm20_ele"; rew=-0.20; break;
      case 10: t=ttbartuple_ele; h=&hrecttbar_rewm25_ele; h1d=&hrecttbar_rewm25_ele_1d; name="hrecttbar_rewm25_ele"; rew=-0.25; break;
	/*      
      case 11: t=ttbartuple_mu; h=&hrecttbar_rew0_mu; h1d=&hrecttbar_rew0_mu_1d; name="hrecttbar_rew0_mu"; rew=0; break;
      case 12: t=ttbartuple_mu; h=&hrecttbar_rewp05_mu; h1d=&hrecttbar_rewp05_mu_1d; name="hrecttbar_rewp05_mu"; rew=0.05; break;
      case 13: t=ttbartuple_mu; h=&hrecttbar_rewp10_mu; h1d=&hrecttbar_rewp10_mu_1d; name="hrecttbar_rewp10_mu"; rew=0.10; break;
      case 14: t=ttbartuple_mu; h=&hrecttbar_rewp15_mu; h1d=&hrecttbar_rewp15_mu_1d; name="hrecttbar_rewp15_mu"; rew=0.15; break;
      case 15: t=ttbartuple_mu; h=&hrecttbar_rewp20_mu; h1d=&hrecttbar_rewp20_mu_1d; name="hrecttbar_rewp20_mu"; rew=0.20; break;
      case 16: t=ttbartuple_mu; h=&hrecttbar_rewp25_mu; h1d=&hrecttbar_rewp25_mu_1d; name="hrecttbar_rewp25_mu"; rew=0.25; break;
      case 17: t=ttbartuple_mu; h=&hrecttbar_rewm05_mu; h1d=&hrecttbar_rewm05_mu_1d; name="hrecttbar_rewm05_mu"; rew=-0.05; break;
      case 18: t=ttbartuple_mu; h=&hrecttbar_rewm10_mu; h1d=&hrecttbar_rewm10_mu_1d; name="hrecttbar_rewm10_mu"; rew=-0.10; break;
      case 19: t=ttbartuple_mu; h=&hrecttbar_rewm15_mu; h1d=&hrecttbar_rewm15_mu_1d; name="hrecttbar_rewm15_mu"; rew=-0.15; break;
      case 20: t=ttbartuple_mu; h=&hrecttbar_rewm20_mu; h1d=&hrecttbar_rewm20_mu_1d; name="hrecttbar_rewm20_mu"; rew=-0.20; break;
      case 21: t=ttbartuple_mu; h=&hrecttbar_rewm25_mu; h1d=&hrecttbar_rewm25_mu_1d; name="hrecttbar_rewm25_mu"; rew=-0.25; break;
	*/
    }

    H2rec myh2(name, xedgesrec, yedgesrec);
    
    char name1d[50] = "";
    sprintf(name1d, "%s_1d", name.Data());
    *h1d = new TH2F(name1d, name1d, nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
    


    t->SetBranchAddress(nameOfXAxisVarRec, &mttbarrec);
    t->SetBranchAddress(nameOfSensVarRec, &diffabsetarec);
    t->SetBranchAddress(nameOfXAxisVarGen, &mttbargen);
    t->SetBranchAddress(nameOfSensVarGen, &diffabsetagen);
    t->SetBranchAddress("weight", &weight);
    weight *= lumi_factor;
    
    const int nentries = t->GetEntries();
    
    int nentries_toskip = split_training_test_pseudo && !force_no_linchecksamplesplit ? int(nentries*splitfactor_training_test) : 0; 
    int nentries_toconsider = nentries;
    
    if(invert_split && split_training_test_pseudo && !force_no_linchecksamplesplit)
    {
      nentries_toconsider = nentries_toskip;
      nentries_toskip = 0;
    }
    
    if(modulo_split)
    {
      nentries_toskip = 0;
      nentries_toconsider = nentries;
    }
    
    for(int j=nentries_toskip; j<nentries_toconsider; j++)
    {
      if(modulo_split && !force_no_linchecksamplesplit)
        if(j%2 == (invert_split xor look_at_samemod_split))
          continue;
        
      //if(gRandom->Integer(2) != invert_split)
      //  continue;
      
      t->GetEntry(j);
      
      if(onlygenbin != -1)
        if(myh2gen.findXBin(mttbargen) != onlygenbin)
          continue;

          myh2.fill(fabs(mttbarrec), diffabsetarec, weight*(1+rew*diffabsetagen));
         fill_nooverflow_2d(*h1d, fabs(mttbarrec), diffabsetarec, weight*(1+rew*diffabsetagen));
    }
    
    *h = myh2.toTH2F();
  }
  
  
  ///init histos for results
  
  TH2F* gen_vs_unf = new TH2F("gen_vs_unf","gen_vs_unf", 100, -2.4, 2.4, 100, -2.4, 2.4);
  TH2F* gen_vs_unf1d = new TH2F("gen_vs_unf1d","gen_vs_unf1d", 100, -2.4, 2.4, 100, -2.4, 2.4);
  
  TH1F* hresult_rew0 = new TH1F("hresult_rew0","hresult_rew0", 100, -2.4, 2.4);
  TH1F* hresult_rewp05 = new TH1F("hresult_rewp05","hresult_rewp05", 100, -2.4, 2.4);
  TH1F* hresult_rewp10 = new TH1F("hresult_rewp10","hresult_rewp10", 100, -2.4, 2.4);
  TH1F* hresult_rewp15 = new TH1F("hresult_rewp15","hresult_rewp15", 100, -2.4, 2.4);
  TH1F* hresult_rewp20 = new TH1F("hresult_rewp20","hresult_rewp20", 100, -2.4, 2.4);
  TH1F* hresult_rewp25 = new TH1F("hresult_rewp25","hresult_rewp25", 100, -2.4, 2.4);
  TH1F* hresult_rewm05 = new TH1F("hresult_rewm05","hresult_rewm05", 100, -2.4, 2.4);
  TH1F* hresult_rewm10 = new TH1F("hresult_rewm10","hresult_rewm10", 100, -2.4, 2.4);
  TH1F* hresult_rewm15 = new TH1F("hresult_rewm15","hresult_rewm15", 100, -2.4, 2.4);
  TH1F* hresult_rewm20 = new TH1F("hresult_rewm20","hresult_rewm20", 100, -2.4, 2.4);
  TH1F* hresult_rewm25 = new TH1F("hresult_rewm25","hresult_rewm25", 100, -2.4, 2.4);  
  
  TH1F* hresult_rew0_1d = new TH1F("hresult_rew0_1d","hresult_rew0_1d", 100, -2.4, 2.4);
  TH1F* hresult_rewp05_1d = new TH1F("hresult_rewp05_1d","hresult_rewp05_1d", 100, -2.4, 2.4);
  TH1F* hresult_rewp10_1d = new TH1F("hresult_rewp10_1d","hresult_rewp10_1d", 100, -2.4, 2.4);
  TH1F* hresult_rewp15_1d = new TH1F("hresult_rewp15_1d","hresult_rewp15_1d", 100, -2.4, 2.4);
  TH1F* hresult_rewp20_1d = new TH1F("hresult_rewp20_1d","hresult_rewp20_1d", 100, -2.4, 2.4);
  TH1F* hresult_rewp25_1d = new TH1F("hresult_rewp25_1d","hresult_rewp25_1d", 100, -2.4, 2.4);
  TH1F* hresult_rewm05_1d = new TH1F("hresult_rewm05_1d","hresult_rewm05_1d", 100, -2.4, 2.4);
  TH1F* hresult_rewm10_1d = new TH1F("hresult_rewm10_1d","hresult_rewm10_1d", 100, -2.4, 2.4);
  TH1F* hresult_rewm15_1d = new TH1F("hresult_rewm15_1d","hresult_rewm15_1d", 100, -2.4, 2.4);
  TH1F* hresult_rewm20_1d = new TH1F("hresult_rewm20_1d","hresult_rewm20_1d", 100, -2.4, 2.4);
  TH1F* hresult_rewm25_1d = new TH1F("hresult_rewm25_1d","hresult_rewm25_1d", 100, -2.4, 2.4);  
  
  TH1F* hlincheck = new TH1F("hlincheck", "hlincheck", 1000, -2.4, 2.4);
  TH1F* hlincheck1d = new TH1F("hlincheck1d", "hlincheck1d", 1000, -2.4, 2.4);
  
  TH1F* hlincheckdiffasys[nbinsmafter];
  for(int i=0; i<nbinsmafter; i++)
  {
    char name[20];
    sprintf(name, "lincheck_diffasybin%d", i);
    hlincheckdiffasys[i] = new TH1F(name, name, 1000, -2.4, 2.4);
  }
  
  TH1::AddDirectory(kFALSE);
  
  
  /// loop over different asymmetry reweighting factors
  for(int ireweight; ireweight <11; ireweight++)
  {
    cout << "Run with reweight factor: " << ireweight << endl;
    
    TString prename;
    TH2F* ttbarrec_histo_ele; 
    //    TH2F* ttbarrec_histo_mu;
    TH2F* ttbarrec_histo_ele_1d; 
    //    TH2F* ttbarrec_histo_mu_1d;
    TH1F* hresult;
    TH1F* hresult1d;
    
    TH1F* diffasys[nbinsmafter]; // contain asymmetry values for a specific m bin
    for(int i=0; i<nbinsmafter; i++)
    {
      char name[20];
      sprintf(name, "asyrecm%d", i);
      diffasys[i] = new TH1F(name, name, 100, -2.4, 2.4);
    }
    
    switch (ireweight)
    {
    case 0:  prename = "preNkm05"; ttbarrec_histo_ele=hrecttbar_rewm05_ele; /*ttbarrec_histo_mu=hrecttbar_rewm05_mu;*/
      ttbarrec_histo_ele_1d=hrecttbar_rewm05_ele_1d; /*ttbarrec_histo_mu_1d=hrecttbar_rewm05_mu_1d;*/
               hresult=hresult_rewm05; hresult1d=hresult_rewm05_1d; break;
    case 1:  prename = "preNkm10"; ttbarrec_histo_ele=hrecttbar_rewm10_ele; /*ttbarrec_histo_mu=hrecttbar_rewm10_mu;*/
      ttbarrec_histo_ele_1d=hrecttbar_rewm10_ele_1d; /*ttbarrec_histo_mu_1d=hrecttbar_rewm10_mu_1d;*/
               hresult=hresult_rewm10; hresult1d=hresult_rewm10_1d; break;
    case 2:  prename = "preNkm15"; ttbarrec_histo_ele=hrecttbar_rewm15_ele; /*ttbarrec_histo_mu=hrecttbar_rewm15_mu;*/
      ttbarrec_histo_ele_1d=hrecttbar_rewm15_ele_1d; /*ttbarrec_histo_mu_1d=hrecttbar_rewm15_mu_1d;*/
               hresult=hresult_rewm15; hresult1d=hresult_rewm15_1d; break;
    case 3:  prename = "preNkm20"; ttbarrec_histo_ele=hrecttbar_rewm20_ele; /*ttbarrec_histo_mu=hrecttbar_rewm20_mu;*/
      ttbarrec_histo_ele_1d=hrecttbar_rewm20_ele_1d; /*ttbarrec_histo_mu_1d=hrecttbar_rewm20_mu_1d;*/
               hresult=hresult_rewm20; hresult1d=hresult_rewm20_1d; break;
    case 4:  prename = "preNkm25"; ttbarrec_histo_ele=hrecttbar_rewm25_ele; /*ttbarrec_histo_mu=hrecttbar_rewm25_mu;*/
      ttbarrec_histo_ele_1d=hrecttbar_rewm25_ele_1d; /*ttbarrec_histo_mu_1d=hrecttbar_rewm25_mu_1d;*/
               hresult=hresult_rewm25; hresult1d=hresult_rewm25_1d; break;
    case 5:  prename = "preNkp05"; ttbarrec_histo_ele=hrecttbar_rewp05_ele; /*ttbarrec_histo_mu=hrecttbar_rewp05_mu;*/
      ttbarrec_histo_ele_1d=hrecttbar_rewp05_ele_1d; /*ttbarrec_histo_mu_1d=hrecttbar_rewp05_mu_1d;*/
               hresult=hresult_rewp05; hresult1d=hresult_rewp05_1d; break;
    case 6:  prename = "preNkp10"; ttbarrec_histo_ele=hrecttbar_rewp10_ele; /*ttbarrec_histo_mu=hrecttbar_rewp10_mu;*/
      ttbarrec_histo_ele_1d=hrecttbar_rewp10_ele_1d; /*ttbarrec_histo_mu_1d=hrecttbar_rewp10_mu_1d;*/
               hresult=hresult_rewp10; hresult1d=hresult_rewp10_1d; break;
    case 7:  prename = "preNkp15"; ttbarrec_histo_ele=hrecttbar_rewp15_ele; /*ttbarrec_histo_mu=hrecttbar_rewp15_mu;*/
      ttbarrec_histo_ele_1d=hrecttbar_rewp15_ele_1d; /*ttbarrec_histo_mu_1d=hrecttbar_rewp15_mu_1d;*/
               hresult=hresult_rewp15; hresult1d=hresult_rewp15_1d; break;
    case 8:  prename = "preNkp20"; ttbarrec_histo_ele=hrecttbar_rewp20_ele; /*ttbarrec_histo_mu=hrecttbar_rewp20_mu;*/
      ttbarrec_histo_ele_1d=hrecttbar_rewp20_ele_1d; /*ttbarrec_histo_mu_1d=hrecttbar_rewp20_mu_1d;*/
               hresult=hresult_rewp20; hresult1d=hresult_rewp20_1d; break;
    case 9:  prename = "preNkp25"; ttbarrec_histo_ele=hrecttbar_rewp25_ele; /*ttbarrec_histo_mu=hrecttbar_rewp25_mu;*/
      ttbarrec_histo_ele_1d=hrecttbar_rewp25_ele_1d; /*ttbarrec_histo_mu_1d=hrecttbar_rewp25_mu_1d;*/
               hresult=hresult_rewp25; hresult1d=hresult_rewp25_1d; break;     
    case 10: prename = "preN"; ttbarrec_histo_ele=hrecttbar_rew0_ele; /*ttbarrec_histo_mu=hrecttbar_rew0_mu;*/
      ttbarrec_histo_ele_1d=hrecttbar_rew0_ele_1d; /*ttbarrec_histo_mu_1d=hrecttbar_rew0_mu_1d;*/
               hresult=hresult_rew0; hresult1d=hresult_rew0_1d; break;
    }
    
    TH2F* hpre = (TH2F*)preselfile->Get(prename);
    TH2F* hpre1d = (TH2F*)preselfile1d->Get(prename);
    
    
    
    // calculate "true" asymmetry of our sample
    const double truepos = hpre->Integral(1, nbinsmafter, nbinsetaafter/2+1, nbinsetaafter);
    const double trueneg = hpre->Integral(1, nbinsmafter, 1, nbinsetaafter/2);
    const double trueasy = (truepos-trueneg)/(truepos+trueneg);
    
    
    // ...and also true asymmetries for the individual mass bins
    double difftrueasys[nbinsmafter];
    for(int i=0; i<nbinsmafter; i++)
    {
      const double pos = hpre->Integral(i+1, i+1, nbinsetaafter/2+1, nbinsetaafter);
      const double neg = hpre->Integral(i+1, i+1, 1, nbinsetaafter/2);
      difftrueasys[i] = (pos-neg)/(pos+neg);
    }
    
    /// perform actual experiments
    for(int n=0; n<numexperiments; n++)
    {
      if(n%1000 == 0)
        cout << "Experiment number " << n << endl;
        
      // create pseudo data histo
      //TH2F hrec("hrec","hrec", nbinsm, mbinedges, nbinseta, etabinedges);
      TH2F hrec("hrec","hrec", nbinsm, 0.5, nbinsm+0.5, nbinseta, 0.5, nbinseta+0.5);
      TH2F hrec1d("hrec1d","hrec1d", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d); 

      // loop over backgrounds and signal to compose pseudo-data-sample
      for(int i=0; i<5; i++)
      {
        TH2F* h = 0; 
        TH2F* h1d = 0; 
        double nfit = 0;
        double fiterror = 0;
        
        switch(i)
        {
	case 0: h=ttbarrec_histo_ele; h1d=ttbarrec_histo_ele_1d; nfit=nttbarele; fiterror=0; break;
	    //          case 1: h=ttbarrec_histo_mu; h1d=ttbarrec_histo_mu_1d; nfit=nttbarmu; fiterror=0; break;
          
          default: const int tempi = i-1;
            h=myunfold.bgwrappedlist[tempi];
            h1d=myunfold1d.bgwrappedlist[tempi];
            nfit=myunfold.bgintegrals[tempi];
            fiterror=myunfold.bgabserrors[tempi];
            break;
        
          //         if(i==1) {h=myunfold.hwjets_pos_ele; h1d=myunfold1d.hwjets_pos_ele; nfit=nwjetsposele; fiterror=nwjetsposeleerror;}
          //         if(i==2) {h=myunfold.hwjets_neg_ele; h1d=myunfold1d.hwjets_neg_ele; nfit=nwjetsnegele; fiterror=nwjetsnegeleerror;}
          //         if(i==3) {h=myunfold.hzjets_ele; h1d=myunfold1d.hzjets_ele; nfit=nzjetsele; fiterror=nzjetseleerror;}
          //         if(i==4) {h=myunfold.hstt_ele; h1d=myunfold1d.hstt_ele; nfit=nstele; fiterror=nsteleerror;} //it's named stt, but it's whole single t
          //         if(i==5) {h=myunfold.hqcd_ele; h1d=myunfold1d.hqcd_ele; nfit=nqcdele; fiterror=nqcdeleerror;}
          //         
          //         if(i==7) {h=myunfold.hwjets_pos_mu; h1d=myunfold1d.hwjets_pos_mu; nfit=nwjetsposmu; fiterror=nwjetsposmuerror;}
          //         if(i==8) {h=myunfold.hwjets_neg_mu; h1d=myunfold1d.hwjets_neg_mu; nfit=nwjetsnegmu; fiterror=nwjetsnegmuerror;}
          //         if(i==9) {h=myunfold.hzjets_mu; h1d=myunfold1d.hzjets_mu; nfit=nzjetsmu; fiterror=nzjetsmuerror;}
          //         if(i==10) {h=myunfold.hstt_mu; h1d=myunfold1d.hstt_mu; nfit=nstmu; fiterror=nstmuerror;} //it's named stt, but it's whole single t
          //         if(i==11) {h=myunfold.hqcd_mu; h1d=myunfold1d.hqcd_mu; nfit=nqcdmu; fiterror=nqcdmuerror;}
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
        
      } // end pseudo sample generation
    
      // unfold it
      TH2F* hunfold = myunfold.unfoldHisto(&hrec, true, drawfactor); // we'll treat the pseudo experiments as if they were real data.
      
      hrec1d.Sumw2(); // recalc errors, hopefully this works
      
      TH2F* hunfold1d = myunfold1d.unfoldHisto(&hrec1d, true, drawfactor);
      
      // calc asy
      
      const double pos = hunfold->Integral(1, nbinsmafter, nbinsetaafter/2+1, nbinsetaafter);
      const double neg = hunfold->Integral(1, nbinsmafter, 1, nbinsetaafter/2);
      const double asy = (pos-neg)/(pos+neg);
      
       //const double asy_error = asymmetryerror_afterunfolding_2d(myunfold.errormatrix, nbinsafter, pos, neg, nbinsetaafter);
      
      const double pos1d = hunfold1d->Integral(1, nbinsmafter1d, nbinsetaafter1d/2+1, nbinsetaafter1d);
      const double neg1d = hunfold1d->Integral(1, nbinsmafter1d, 1, nbinsetaafter1d/2);
      const double asy1d = (pos1d-neg1d)/(pos1d+neg1d);
       //const double asy_error1d = asymmetryerror_afterunfolding_1d(myunfold1d.errormatrix, myunfold1d.hunfoldunwrapped);
      
      // fill histograms
      gen_vs_unf->Fill(trueasy, asy);
      gen_vs_unf1d->Fill(trueasy, asy1d);
      hresult->Fill(asy);
      hresult1d->Fill(asy1d);
      
      for(int i=0; i<nbinsmafter; i++)
      {
        const double pos = hunfold->Integral(i+1, i+1, nbinsetaafter/2+1, nbinsetaafter);
        const double neg = hunfold->Integral(i+1, i+1, 1, nbinsetaafter/2);
        const double asy = (pos-neg)/(pos+neg);
        diffasys[i]->Fill(asy);
      }
      
    } // end experiment loop
    
    // generate money plots
    
    
    #define GETBINERROR GetRMS // i'm very sorry... use this to switch GetRMS and GetMeanError
   
    hlincheck->SetBinContent(hlincheck->FindFixBin(trueasy), hresult->GetMean());
    hlincheck->SetBinError(hlincheck->FindFixBin(trueasy), hresult->GETBINERROR());
    
    hlincheck1d->SetBinContent(hlincheck1d->FindFixBin(trueasy), hresult1d->GetMean());
    hlincheck1d->SetBinError(hlincheck1d->FindFixBin(trueasy), hresult1d->GETBINERROR());
    
    cout << hresult->GetName() << "\t" << trueasy << "\t" << hresult->GetMean() << "\t" << hresult->GETBINERROR() << endl;
    cout << hresult1d->GetName() << "\t" << trueasy << "\t" << hresult1d->GetMean() << "\t" << hresult1d->GETBINERROR() << endl;
    
    rewfile << hresult->GetName() << "\t" << trueasy << "\t" << hresult->GetMean() << "\t" << hresult->GETBINERROR() << endl;
    rewfile << hresult1d->GetName() << "\t" << trueasy << "\t" << hresult1d->GetMean() << "\t" << hresult1d->GETBINERROR() << endl;
    for(int i=0; i<nbinsmafter; i++)
    {
      rewfile << hresult->GetName() << " bin " << i << "\t" << difftrueasys[i] << "\t" << diffasys[i]->GetMean() << "\t" << diffasys[i]->GETBINERROR() << endl;
    }
    
    // differential money plots
    
    for(int i=0; i<nbinsmafter; i++)
    {
      hlincheckdiffasys[i]->SetBinContent(hlincheckdiffasys[i]->FindFixBin(difftrueasys[i]), diffasys[i]->GetMean());
      hlincheckdiffasys[i]->SetBinError(hlincheckdiffasys[i]->FindFixBin(difftrueasys[i]), diffasys[i]->GETBINERROR());
    }
    
    // cleanup
    for(int i=0; i<nbinsmafter; i++)
    {
      delete diffasys[i];
    }
  } // end reweight factor loop

  TH1::AddDirectory(kTRUE);

  /// do fits and save results
  
  TF1* fitfunc = new TF1("fitfunc", "[0]*x+[1]", -2.4, 2.4);
  
  TString genbinstring = "";
  if(onlygenbin!=-1)
  {
    genbinstring += "bin";
    genbinstring += onlygenbin;
  }
  
  ofstream file(TString(TString("output/")+fileSuffix+"/lincheck_fitresult"+genbinstring+".hpp"));
  
  file << "#ifndef LINCHECK_FITRESULT_HPP\n#define LINCHECK_FITRESULT_HPP\n\n";
  
  hlincheck->Fit(fitfunc);
  
  file << "const double lincheck_slope +/- lincheck_slope_error = " << fitfunc->GetParameter(0) << "+/-" << fitfunc->GetParError(0) << ";" << endl;
  file << "const double lincheck_offset +/- lincheck_offset = " << fitfunc->GetParameter(1) << "+/-" << fitfunc->GetParError(1) << ";" << endl;
  
  hlincheck1d->Fit(fitfunc);
  file << "const double lincheck_slope1d +/- lincheck_slope1d_error = " << fitfunc->GetParameter(0) << "+/-" << fitfunc->GetParError(0) << ";" << endl;
  file << "const double lincheck_offset1d +/- lincheck_offset1d_error = " << fitfunc->GetParameter(1) << "+/-" << fitfunc->GetParError(1) << ";" << endl;
  
  double slopelist[nbinsmafter];  
  double offsetlist[nbinsmafter];  
  double slopelist_error[nbinsmafter];  
  double offsetlist_error[nbinsmafter];  
    
  for(int i=0; i<nbinsmafter; i++)
  {
    hlincheckdiffasys[i]->Fit(fitfunc);
    slopelist[i] = fitfunc->GetParameter(0);
    offsetlist[i] = fitfunc->GetParameter(1);
    slopelist_error[i] = fitfunc->GetParError(0);
    offsetlist_error[i] = fitfunc->GetParError(1);

  }

  file << "const double lincheck_diffslope[] +/- lincheck_diffslope_error[] = {";
  for(int i=0; i<nbinsmafter; i++)
  {
    file << slopelist[i];    

    file << "+/-";    

    file << slopelist_error[i];    
    
    if(i!=nbinsmafter-1)
      file << ",";    
  }
  file << "};\n"; 
  
  file << "const double lincheck_diffoffset[] +/- lincheck_diffoffset[] = {";
  for(int i=0; i<nbinsmafter; i++)
  {
    file << offsetlist[i];

    file << "+/-";    

    file << offsetlist_error[i];    
    
    if(i!=nbinsmafter-1)
      file << ",";    
  }
  file << "};\n\n"; 
  
  file << "#endif // LINCHECK_FITRESULT_HPP\n";
  
  file.close();
  
  /// write our result histograms into a file
  TFile out(TString("output/")+fileSuffix+"/"+"lincheck"+genbinstring+".root","recreate");
  
  gen_vs_unf->Write();
  gen_vs_unf1d->Write();
  
  hresult_rew0->Write();
  hresult_rewp05->Write();
  hresult_rewp10->Write();
  hresult_rewp15->Write();
  hresult_rewp20->Write();
  hresult_rewp25->Write();
  hresult_rewm05->Write();
  hresult_rewm10->Write();
  hresult_rewm15->Write();
  hresult_rewm20->Write();
  hresult_rewm25->Write();
  
  hresult_rew0_1d->Write();
  hresult_rewp05_1d->Write();
  hresult_rewp10_1d->Write();
  hresult_rewp15_1d->Write();
  hresult_rewp20_1d->Write();
  hresult_rewp25_1d->Write();
  hresult_rewm05_1d->Write();
  hresult_rewm10_1d->Write();
  hresult_rewm15_1d->Write();
  hresult_rewm20_1d->Write();
  hresult_rewm25_1d->Write();
    
  hlincheck->Write();
  hlincheck1d->Write();
  
  for(int i=0; i<nbinsmafter; i++)
  {
    hlincheckdiffasys[i]->Write();
  }
  
  out.Close();
  
  rewfile.close();
  
  cout << "Finished successfully" << endl;

  return 0;
}



int main()
{
  return lincheck();
}
