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
#include <iomanip>


#include <sys/stat.h> // for mkdir
#include <sys/types.h>

#include "../cmsstyle.hpp"
#include "../fitresults.hpp"
#include "../helpers.hpp"
#include "../samplelocations.hpp"


#include "../binning.hpp"

// WARNING make sure variables in binning.hpp are set correctly


#include "../myunfold_class.hpp"
#include "../myunfold_class1d.hpp"


const int numexperiments = 5000; // should be about 50,000


int reweighted_pseudoexp()
{
  bool inverted = false;
  bool centered = false;
  bool negative = false;
  bool zprime = false; // doesnt combine with rest
  bool contsquare = false; // doesnt combine with rest
  
  string env_invert = getEnv("PSEUDOWEIGHT_INVERT");
  if(env_invert.find("TRUE") != string::npos) inverted = true;
  else if(env_invert.find("FALSE") != string::npos) inverted = false;
  else
  {
    cerr << 1 << "ENVVAR WRONG, EXITING" << endl;
    return 1; 
  }
  
  string env_negative = getEnv("PSEUDOWEIGHT_NEGATIVE");
  if(env_negative.find("TRUE") != string::npos) negative = true;
  else if(env_negative.find("FALSE") != string::npos) negative = false;
  else
  {
    cerr << 2 << "ENVVAR WRONG, EXITING" << endl;
    return 1; 
  }
  
  string env_center = getEnv("PSEUDOWEIGHT_CENTER");
  if(env_center.find("TRUE") != string::npos) centered = true;
  else if(env_center.find("FALSE") != string::npos) centered = false;
  else
  {
    cerr << 3 << "ENVVAR WRONG, EXITING" << endl;
    return 1; 
  }
  
  string env_zprime = getEnv("PSEUDOWEIGHT_ZPRIME");
  if(env_zprime.find("TRUE") != string::npos) zprime = true;
  else if(env_zprime.find("FALSE") != string::npos) zprime = false;
  else
  {
    cerr << 4 << "ENVVAR WRONG, EXITING" << endl;
    return 1; 
  }
  
  string env_square = getEnv("PSEUDOWEIGHT_SQUARE");
  if(env_square.find("TRUE") != string::npos) contsquare = true;
  else if(env_square.find("FALSE") != string::npos) contsquare = false;
  else
  {
    cerr << 5 << "ENVVAR WRONG, EXITING" << endl;
    return 1; 
  }
  
  
  TString rewvar = "";
  string env_rewvar = getEnv("REWEIGHT_VAR");
  if(env_rewvar.find("M") != string::npos) rewvar = "m";
  else if(env_rewvar.find("PT") != string::npos) rewvar = "pt";
  else if(env_rewvar.find("Y") != string::npos) rewvar = "y";
  else
  {
    cerr << 6 << "ENVVAR WRONG, EXITING" << endl;
    return 1; 
  }
  
  bool donotreweight;
  string env_notrew = getEnv("DONOTREWEIGHT");
  if(env_notrew.find("TRUE") != string::npos) donotreweight = true;
  else if(env_notrew.find("FALSE") != string::npos) donotreweight = false;
  else
  {
    cerr << 7 << "ENVVAR WRONG, EXITING" << endl;
    return 1; 
  }
  
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
//  TFile* preselfile1d = new TFile(TString("output/")+fileSuffix+"/1dseleff.root");
  TFile* preselfile = new TFile("seleff.root");
  TFile* preselfile1d = new TFile("1dseleff.root");
  
  TH2F* hpresel = 0;
  TH2F* hpresel1d = 0;
  
  TH2F* preselected_nonflat = 0;
  
  
  TH2F* zprime_nonflat = 0;
  TH1F* zprimeposplotz = 0;
  TH1F* zprimenegplotz = 0;
  TH1F* zprimeposplot = 0;
  TH1F* zprimenegplot = 0;
  double intfrac = 1;
  
  
  TString preselname = "";

  {  
    if(contsquare)
      preselname = "preNkcontsquare";
    else if(inverted && !centered && !negative)
      preselname = "preNkcontinv";
    else if(!inverted && !centered && !negative)
      preselname = "preNkcont";
    else if(inverted && !centered && negative)
      preselname = "preNkcontneginv";
    else if(!inverted && !centered && negative)
      preselname = "preNkcontneg";
    else if(!inverted && centered)
      preselname = "preNkcontcent";
    else if(inverted && centered)
      preselname = "preNkcontcentinv";
    
    preselname += "_" + rewvar;
    
    if(donotreweight)
      preselname = "preN";
    
    hpresel = (TH2F*)preselfile->Get(preselname);
    hpresel1d = (TH2F*)preselfile1d->Get(preselname);
  }

  
  double difftrueasys[nbinsmafter];
  for(int i=0; i<nbinsmafter; i++)
  {
    const double truepos = hpresel->Integral(i+1, i+1, nbinsetaafter/2+1, nbinsetaafter);
    const double trueneg = hpresel->Integral(i+1, i+1, 1, nbinsetaafter/2);
    difftrueasys[i] = (truepos-trueneg)/(truepos+trueneg);
  }
  
  TH1F* hpreselunwrapped = new TH1F("hpreselunwrapped","hpreselunwrapped", nbinsafter, 0.5, nbinsafter+0.5); ;
  unwrap2dhisto(hpresel, hpreselunwrapped);
  cmsstyle->setup_style_2D(hpresel, labelOfXAxisVar, labelOfSensVar);
  
  TH1F* hpreselunwrapped1d = new TH1F("hpreselunwrapped1d","hpreselunwrapped1d", nbinsafter1d, 0.5, nbinsafter1d+0.5); ;
  unwrap2dhisto(hpresel1d, hpreselunwrapped1d);
  cmsstyle->setup_style_2D(hpresel1d, labelOfXAxisVar, labelOfSensVar);
  
  
  MyUnfold myunfold;
  MyUnfold myunfold1d;
  
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
  
  TH1F* asys = new TH1F("asys", "asys", 100, -0.3, 0.3);
  TH1F* asys1d = new TH1F("asys1d", "asys1d", 100, -0.3, 0.3);
  
  TH1F* asypull = new TH1F("asypull", "asypull", 100, -3, 3);
  TH1F* asypull1d = new TH1F("asypull1d", "asypull1d", 100, -3, 3);
  
  TH1F* asydiff = new TH1F("asydiff", "asydiff", 100, -0.1, 0.1);
  TH1F* asydiff1d = new TH1F("asydiff1d", "asydiff1d", 100, -0.1, 0.1);
  
  TH1F* ndiff = new TH1F("ndiff", "ndiff", 100, -120000, 120000);
  TH1F* ndiff1d = new TH1F("ndiff1d", "ndiff1d", 100, -120000, 120000);

  
  TH1F* asyerrors = new TH1F("asyerrors", "asyerrors", 100, 0.02, 0.05);
  TH1F* asyerrors1d = new TH1F("asyerrors1d", "asyerrors1d", 100, 0.02, 0.05);
  
  
  // init histos for relative differences bin-by-bin
  TH1F* reldiffs[nbinsafter];
  TH1F* pulls[nbinsafter];
  for(int i=0; i<nbinsafter; i++)
  {
    char name[20];
    sprintf(name, "reldiff%d", i);
    reldiffs[i] = new TH1F(name, name, 100, -0.6, 0.6);
    sprintf(name, "pull%d", i);
    pulls[i] = new TH1F(name, name, 100, -5, 5);
  }
  
  TH1F* reldiffs1d[nbinsafter1d];
  TH1F* pulls1d[nbinsafter1d];
  for(int i=0; i<nbinsafter1d; i++)
  {
    char name[20];
    sprintf(name, "1dreldiff%d", i);
    reldiffs1d[i] = new TH1F(name, name, 100, -0.6, 0.6);
    sprintf(name, "1dpull%d", i);
    pulls1d[i] = new TH1F(name, name, 100, -5, 5);
  }
  
  // histos for asymmetries in the different mass bins
  TH1F* diffasys[nbinsmafter];
  for(int i=0; i<nbinsmafter; i++)
  {
    char name[20];
    sprintf(name, "diffasybin%d", i);
    diffasys[i] = new TH1F(name, name, 100, -0.5, 0.5);
  }
  
  TH1F* diffasydiffs[nbinsmafter];
  TH1F* diffasypulls[nbinsmafter];
  for(int i=0; i<nbinsmafter; i++)
  {
    char name[20];
    sprintf(name, "diffasydiffbin%d", i);
    diffasydiffs[i] = new TH1F(name, name, 100, -0.35, 0.35);
    sprintf(name, "diffasypullbin%d", i);
    diffasypulls[i] = new TH1F(name, name, 100, -3.5, 3.5);
  }
  
  // histos for KS tests between data and reweighted reco
  TH1F* dataks_m = new TH1F("dataks_m", "dataks_m", 100, 0, 1500);
  TH1F* dataks_pt = new TH1F("dataks_pt", "dataks_pt", 100, 0, 400);
  TH1F* dataks_y = new TH1F("dataks_y", "dataks_y", 100, -2.5, 2.5);
  TH1F* dataks_diffabsy = new TH1F("dataks_diffabsy", "dataks_diffabsy", 100, -4, 4);
  
  TH1F* mcks_m_ele = new TH1F("mcks_m_ele", "mcks_m_ele", 100, 0, 1500);
  TH1F* mcks_pt_ele = new TH1F("mcks_pt_ele", "mcks_pt_ele", 100, 0, 400);
  TH1F* mcks_y_ele = new TH1F("mcks_y_ele", "mcks_y_ele", 100, -2.5, 2.5);
  TH1F* mcks_diffabsy_ele = new TH1F("mcks_diffabsy_ele", "mcks_diffabsy_ele", 100, -4, 4);
  /*  TH1F* mcks_m_mu = new TH1F("mcks_m_mu", "mcks_m_mu", 100, 0, 1500);
  TH1F* mcks_pt_mu = new TH1F("mcks_pt_mu", "mcks_pt_mu", 100, 0, 400);
  TH1F* mcks_y_mu = new TH1F("mcks_y_mu", "mcks_y_mu", 100, -2.5, 2.5);
  TH1F* mcks_diffabsy_mu = new TH1F("mcks_diffabsy_mu", "mcks_diffabsy_mu", 100, -4, 4);
  */
  TH1F *bg_m, *bg_pt, *bg_y, *bg_diffabsy;
  
  {
    TFile* fwjetsposele = new TFile(sl_wjetsposele);
    //TFile* fwjetsnegele = new TFile(sl_wjetsnegele);
    TFile* fzjetsele = new TFile(sl_zjetsele);
    //    TFile* fsinglet_tchan_t_ele = new TFile(sl_singlet_tchan_t_ele);
    //    TFile* fsinglet_tchan_tbar_ele = new TFile(sl_singlet_tchan_tbar_ele);
    TFile* fsinglet_twchan_t_ele = new TFile(sl_singlet_twchan_t_ele);
    TFile* fsinglet_twchan_tbar_ele = new TFile(sl_singlet_twchan_tbar_ele);
    //TFile* fqcdele = new TFile(sl_qcdele);
        
    /*    TFile* fwjetsposmu = new TFile(sl_wjetsposmu);
    TFile* fwjetsnegmu = new TFile(sl_wjetsnegmu);
    TFile* fzjetsmu = new TFile(sl_zjetsmu);
    TFile* fsinglet_tchan_t_mu = new TFile(sl_singlet_tchan_t_mu);
    TFile* fsinglet_tchan_tbar_mu = new TFile(sl_singlet_tchan_tbar_mu);
    TFile* fsinglet_twchan_t_mu = new TFile(sl_singlet_twchan_t_mu);
    TFile* fsinglet_twchan_tbar_mu = new TFile(sl_singlet_twchan_tbar_mu);
    TFile* fqcdmu = new TFile(sl_qcdmu);
    */
    TH1F** h;
    TString histoname;
    
    for(int i=0; i<4; i++)
    {
      if(i==0) { histoname = "TTbarSelHyppsiprime/m_ttbar"; h=&bg_m;}
      if(i==1) { histoname = "TTbarSelHyppsiprime/pt_ttbar"; h=&bg_pt;}
      if(i==2) { histoname = "TTbarSelHyppsiprime/y_ttbar"; h=&bg_y;}
      if(i==3) { histoname = "TTbarSelHyppsiprime/AAE_diffabsy"; h=&bg_diffabsy;}
      
      TH1F* hwjetsele = (TH1F*)fwjetsposele->Get(histoname);
      hwjetsele->Scale(pred_pos_wjets_ele_met*betawjetsposele/hwjetsele->Integral());
      //TH1F* hwjetsnegele = (TH1F*)fwjetsnegele->Get(histoname);
      //hwjetsnegele->Scale(pred_neg_wjets_ele_met*betawjetsnegele/hwjetsnegele->Integral());
      //hwjetsele->Add(hwjetsnegele);
      
      
      TH1F* hzjetsele = (TH1F*)fzjetsele->Get(histoname);
      hzjetsele->Scale((pred_pos_zjets_ele_met/*+pred_neg_zjets_ele_met*/)*betazjetsele/hzjetsele->Integral());
      
      //      TH1F* hsinglet_tchanele = (TH1F*)fsinglet_tchan_t_ele->Get(histoname);
      //      hsinglet_tchanele->Scale((pred_pos_st_t_top_ele_met/*+pred_neg_st_t_top_ele_met*/)*betastele/hsinglet_tchanele->Integral());
      //      TH1F* hsinglet_tchanele_tbar = (TH1F*)fsinglet_tchan_tbar_ele->Get(histoname);
      
      //      hsinglet_tchanele_tbar->Scale((pred_pos_st_t_atop_ele_met/*+pred_neg_st_t_atop_ele_met*/)*betastele/hsinglet_tchanele_tbar->Integral());
      //      hsinglet_tchanele->Add(hsinglet_tchanele_tbar);
      
      
      TH1F* hsinglet_twchanele = (TH1F*)fsinglet_twchan_t_ele->Get(histoname);
      hsinglet_twchanele->Scale((pred_pos_st_tw_top_ele_met/*+pred_neg_st_tw_top_ele_met*/)*betastele/hsinglet_twchanele->Integral());
      TH1F* hsinglet_twchanele_tbar = (TH1F*)fsinglet_twchan_tbar_ele->Get(histoname);
      hsinglet_twchanele_tbar->Scale((pred_pos_st_tw_top_ele_met/*+pred_neg_st_tw_atop_ele_met*/)*betastele/hsinglet_twchanele_tbar->Integral());
      hsinglet_twchanele->Add(hsinglet_twchanele_tbar);
    
      //TH1F* hqcdele = (TH1F*)fqcdele->Get(histoname);
      //hqcdele->Scale((pred_pos_qcd_ele_met+pred_neg_qcd_ele_met)*betaqcdele/hqcdele->Integral());
      
      
      
      // mu 
      
      //TH1F* hwjets = (TH1F*)fwjetsposele->Get(histoname);
      //hwjets->Scale(pred_pos_wjets_mu_met*betawjetsposmu/hwjets->Integral());
      //TH1F* hwjetsneg = (TH1F*)fwjetsnegmu->Get(histoname);
      //hwjetsneg->Scale(pred_neg_wjets_mu_met*betawjetsnegmu/hwjetsneg->Integral());
      //hwjets->Add(hwjetsneg);
      
      
      //TH1F* hzjets = (TH1F*)fzjetsele->Get(histoname);
      //hzjets->Scale((pred_pos_zjets_mu_met+pred_neg_zjets_mu_met)*betazjetsmu/hzjets->Integral());
      
      
      //TH1F* hsinglet_tchan = (TH1F*)fsinglet_tchan_t_ele->Get(histoname);
      //hsinglet_tchan->Scale((pred_pos_st_t_top_mu_met+pred_neg_st_t_top_mu_met)*betastmu/hsinglet_tchan->Integral());
      //TH1F* hsinglet_tchan_tbar = (TH1F*)fsinglet_tchan_tbar_ele->Get(histoname);
      //hsinglet_tchan_tbar->Scale((pred_pos_st_t_atop_mu_met+pred_neg_st_t_atop_mu_met)*betastmu/hsinglet_tchan_tbar->Integral());
      //hsinglet_tchan->Add(hsinglet_tchan_tbar);
      
      
      //TH1F* hsinglet_twchan = (TH1F*)fsinglet_twchan_t_ele->Get(histoname);
      //hsinglet_twchan->Scale((pred_pos_st_tw_top_mu_met+pred_neg_st_tw_top_mu_met)*betastmu/hsinglet_twchan->Integral());
      //TH1F* hsinglet_twchan_tbar = (TH1F*)fsinglet_twchan_tbar_ele->Get(histoname);
      //hsinglet_twchan_tbar->Scale((pred_pos_st_tw_atop_mu_met+pred_neg_st_tw_atop_mu_met)*betastmu/hsinglet_twchan_tbar->Integral());
      //hsinglet_twchan->Add(hsinglet_twchan_tbar);        
      
      //TH1F* hqcd = (TH1F*)fqcdele->Get(histoname);
      //hqcd->Scale((pred_pos_qcd_mu_met+pred_neg_qcd_mu_met)*betaqcdmu/hqcd->Integral());
      
      // combine
      //      hwjets->Add(hwjetsele);
      //      hzjets->Add(hzjetsele);
      //      hsinglet_tchan->Add(hsinglet_tchanele);
      //      hsinglet_twchan->Add(hsinglet_twchanele);
      //      hqcd->Add(hqcdele);
      /*    
      hsinglet_tchan->Add(hsinglet_twchan);
      
      hzjets->Add(hqcd);
      hwjets->Add(hzjets);
      hsinglet_tchan->Add(hwjets);
      */

      //hsinglet_tchanele->Add(hsinglet_twchanele);
      
      //hzjets->Add(hqcd);
      hwjetsele->Add(hzjetsele);
      //      hsinglet_tchanele->Add(hwjetsele);
      hsinglet_twchanele->Add(hwjetsele);



      *h = hsinglet_twchanele;
    }
    
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
  
  
  // init reweighted reconstructed histos
  
  #ifdef MADGRAPH
//    TFile* file_ttbar_ele = new TFile("/portal/ekpams2/home/froscher/outputroot/2dunfolding/ttbar_ele_histos.root");
    TFile* file_ttbar_ele = new TFile("Tuples_and_Results_23_April_need_fix_MCatNLO_et_QCD/TupleFileSkimmed_TTbarAllMCatNLO.root");
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
      file_ttbar_ele = new TFile("Tuples_and_Results_23_April_need_fix_MCatNLO_et_QCD/TupleFileSkimmed_TTbarAllMCatNLO.root");
      //      file_ttbar_mu = new TFile("FIXME/portal/ekpams2/home/froscher/outputroot/powhegsample/ttbar_mu_combined_subsample2_histos.root");
    }
  #endif
  
//  TNtuple* ttbartuple_ele = (TNtuple*)file_ttbar_ele->Get("TTbarSelHyppsiprime/BBB_UnfoldingTuple");
//  TNtuple* ttbartuple_mu = (TNtuple*)file_ttbar_mu->Get("TTbarSelHyppsiprime/BBB_UnfoldingTuple");
  TNtuple* ttbartuple_ele = (TNtuple*)file_ttbar_ele->Get("Tuple");
  //  TNtuple* ttbartuple_mu = (TNtuple*)file_ttbar_mu->Get("Tuple");
  
  TH2F* hrecttbar_rew_ele;
  //  TH2F* hrecttbar_rew_mu;
  TH2F* hrecttbar_rew_ele1d;
  //  TH2F* hrecttbar_rew_mu1d;
  
  TH2F* rew_selected = new TH2F("rew_selected","rew_selected", 100, xedgesgen[0], 2* xedgesgen[nbinsmafter],
                                100, 2* yedgesgen[0][0], 2* yedgesgen[0][nbinsetaafter]); 
  TH2F* rew_reco = new TH2F("rew_reco","rew_reco", 100, xedgesgen[0], 2* xedgesgen[nbinsmafter],
                                100, 2* yedgesgen[0][0], 2* yedgesgen[0][nbinsetaafter]);

                                
  H2rec reco_rew("reco_rew", xedgesrec, yedgesrec);
  H2rec reco_orig("reco_orig", xedgesrec, yedgesrec);
  H2rec reco_data("reco_data", xedgesrec, yedgesrec);
                                

  H2gen myh2gen("myh2gen", xedgesgen, yedgesgen);
  
  float *mttbarrec, diffabsetarec, mttbargen, diffabsetagen, weight;
  float mrec, ptrec, yrec;
  float xgen;
  
  for(int i=0; i<1; i++)
  {
    cout << "Initializing histo " << i << endl;
    
    TNtuple* t = 0;
    TH2F** h = 0; 
    TH2F** h1d = 0; 
    TString name = "";
        
    switch (i)
    {
      case 0: t=ttbartuple_ele; h=&hrecttbar_rew_ele; h1d=&hrecttbar_rew_ele1d; name="hrecttbar_rew_ele"; break;
	//      case 1: t=ttbartuple_mu; h=&hrecttbar_rew_mu; h1d=&hrecttbar_rew_mu1d; name="hrecttbar_rew_mu"; break;
    }
    
    H2rec myh2(name, xedgesrec, yedgesrec);   
    H2rec1d myh21d(name+"1d", xedgesrec1d, yedgesrec1d);   
    
    t->SetBranchAddress(nameOfSensVarRec, &diffabsetarec);
    t->SetBranchAddress(nameOfXAxisVarGen, &mttbargen);
    t->SetBranchAddress(nameOfSensVarGen, &diffabsetagen);
    
    double reweighting_loweredge, reweighting_upperedge;
    
    if(!rewvar.CompareTo("m"))
    {
      t->SetBranchAddress("mttbar_gen", &xgen);
      reweighting_loweredge = 340;
      reweighting_upperedge = 650;
    }
    if(!rewvar.CompareTo("pt"))
    {
      t->SetBranchAddress("pt_ttbar_gen", &xgen);
      reweighting_loweredge = 0;
      reweighting_upperedge = 75;
    }  
    if(!rewvar.CompareTo("y"))
    {
      t->SetBranchAddress("y_ttbar_gen", &xgen);
      reweighting_loweredge = 0;
      reweighting_upperedge = 1.1;
    }
    
    
    
    if(!TString(nameOfXAxisVarRec).CompareTo("mttbar_rec")) mttbarrec = &mrec;
    if(!TString(nameOfXAxisVarRec).CompareTo("pt_ttbar")) mttbarrec = &ptrec;
    if(!TString(nameOfXAxisVarRec).CompareTo("y_ttbar_rec")) mttbarrec = &yrec;
    
    t->SetBranchAddress("mttbar_rec", &mrec);
    t->SetBranchAddress("pt_ttbar", &ptrec);
    t->SetBranchAddress("y_ttbar_rec", &yrec);
    t->SetBranchAddress("weight", &weight);
    
    int nentries = t->GetEntries();
    int nentries_toskip = split_training_test_pseudo ? int(nentries*splitfactor_training_test) : 0; 
    int nentries_toconsider = nentries;
    
    if(invert_split && split_training_test_pseudo)
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
      if(modulo_split)
        if(j%2 == invert_split)
          continue;
        
      t->GetEntry(j);        
      
      double valx = fabs(mttbargen);
      const double valy = diffabsetagen;
      double contk = 0 ;
      double contweight = 0;
      
      {
        if(discont)
          valx = myh2gen.getXBinCenter(myh2gen.findXBin(valx));
        
        const double wfactor = 0.125; // used to be 0.25
        if(contsquare)
        {
          contk = 2*wfactor * (fabs(xgen)-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge)
                            * (fabs(xgen)-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge);
          if(contk < 0) contk = 0;
          if(contk > wfactor) contk = wfactor;
        }          
        else
        {
          contk = wfactor * (fabs(xgen)-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge);
          if(contk < 0) contk = 0;
          if(contk > wfactor) contk = wfactor;
          
          if(inverted && !centered && !negative)
            contk = wfactor - contk; // invert reweighting to be bigger for smaller M
          if(inverted && !centered && negative)
            contk = -1*(wfactor - contk); 
          if(!inverted && !centered && negative)
            contk = -1*contk; 
            
          if(!inverted && centered)
            contk = contk-(wfactor/2); // center it
          if(inverted && centered)
            contk = -1*(contk-(wfactor/2)); // center it, invert
        }
        
        contweight = (1+contk*diffabsetagen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
      }
      
      if(donotreweight)
        contweight = 1.0;
      
      myh2.fill(fabs(*mttbarrec), diffabsetarec, weight*contweight);
      myh21d.fill(fabs(*mttbarrec), diffabsetarec, weight*contweight);      
      
      reco_rew.fill(fabs(*mttbarrec), diffabsetarec, weight*contweight);
      reco_orig.fill(fabs(*mttbarrec), diffabsetarec, weight);
      
      fill_nooverflow_2d(rew_selected, mttbargen, diffabsetagen, weight*contweight);
      fill_nooverflow_2d(rew_reco, *mttbarrec, diffabsetarec, weight*contweight);
      
      if(!i)
      {
        fill_nooverflow_1d(mcks_m_ele, mrec, weight*contweight);
        fill_nooverflow_1d(mcks_pt_ele, ptrec, weight*contweight);
        fill_nooverflow_1d(mcks_y_ele, yrec, weight*contweight);
        fill_nooverflow_1d(mcks_diffabsy_ele, diffabsetarec, weight*contweight);
      }
      /*      else
      {
        fill_nooverflow_1d(mcks_m_mu, mrec, weight*contweight);
        fill_nooverflow_1d(mcks_pt_mu, ptrec, weight*contweight);
        fill_nooverflow_1d(mcks_y_mu, yrec, weight*contweight);
        fill_nooverflow_1d(mcks_diffabsy_mu, diffabsetarec, weight*contweight);
	}*/
    }
    
    *h = myh2.toTH2F(); 
    *h1d = myh21d.toTH2F(); 
  }
  
  scale_to(mcks_m_ele, nttbarele);
  scale_to(mcks_pt_ele, nttbarele);
  scale_to(mcks_y_ele, nttbarele);
  scale_to(mcks_diffabsy_ele, nttbarele);
  /*  scale_to(mcks_m_mu, nttbarmu);
  scale_to(mcks_pt_mu, nttbarmu);
  scale_to(mcks_y_mu, nttbarmu);
  scale_to(mcks_diffabsy_mu, nttbarmu);
  */
  /// get data tuple. only needed to see how much reconstructed values differ from MC expectation
  {
    TFile* inputfile = new TFile(sl_data_comb);
//    TNtuple* tuple = (TNtuple*)inputfile->Get("TTbarSelHyppsiprime/BBB_UnfoldingTuple");
    TNtuple* tuple = (TNtuple*)inputfile->Get("Tuple");
    
    float mrec, ptrec, yrec, diffabsetarec, weight;
    
    tuple->SetBranchAddress("mttbar_rec", &mrec);
    tuple->SetBranchAddress("pt_ttbar", &ptrec);
    tuple->SetBranchAddress("y_ttbar_rec", &yrec);
    tuple->SetBranchAddress(nameOfSensVarRec, &diffabsetarec);
    tuple->SetBranchAddress("weight", &weight);
    
    float* xvar = 0;
    if(!TString(nameOfXAxisVarRec).CompareTo("mttbar_rec")) xvar = &mrec;
    if(!TString(nameOfXAxisVarRec).CompareTo("pt_ttbar")) xvar = &ptrec;
    if(!TString(nameOfXAxisVarRec).CompareTo("y_ttbar_rec")) xvar = &yrec;
    
    const int nentries = tuple->GetEntries();
    
    for(int i=0; i<nentries; i++)
    {
      tuple->GetEntry(i);

      weight = 1.;      
      reco_data.fill(fabs(*xvar), diffabsetarec, weight);
      
      fill_nooverflow_1d(dataks_m, mrec, weight);
      fill_nooverflow_1d(dataks_pt, ptrec, weight);
      fill_nooverflow_1d(dataks_y, yrec, weight);
      fill_nooverflow_1d(dataks_diffabsy, diffabsetarec, weight);
    }   
    
  }



  
  TH1::AddDirectory(kFALSE);
  
  /// perform actual experiments
  for(int n=0; n<numexperiments; n++)
  {
    if(n%1000 == 0)
      cout << "Experiment number " << n << endl;
      
    // create pseudo data histo
    //TH2F hrec("hrec","hrec", nbinsm, mbinedges, nbinseta, etabinedges);
    TH2F hrec("hrec","hrec", nbinsm, 0.5, nbinsm+0.5, nbinseta, 0.5, nbinseta+0.5);
    TH2F hrec1d("hrec1d","hrec1d", nbinsm1d, 0.5, nbinsm1d+0.5, nbinseta1d, 0.5, nbinseta1d+0.5);

    // loop over backgrounds and signal to compose pseudo-data-sample
    for(int i=0; i<12; i++)
    {
      TH2F* h = 0; 
      TH2F* h1d = 0; 
      double nfit = 0;
      double fiterror = 0;
      
      //if(i!=0 && i!=6) continue; // to turn off bg
      
      switch (i)
      {
        case 0:  h=hrecttbar_rew_ele; h1d=hrecttbar_rew_ele1d; nfit=nttbarele; fiterror=0; break; 
	  //        case 1:  h=hrecttbar_rew_mu; h1d=hrecttbar_rew_mu1d; nfit=nttbarmu; fiterror=0; break;      
        
        default: const int tempi = i-2;
          h=myunfold.bgwrappedlist[tempi];
          h1d=myunfold1d.bgwrappedlist[tempi];
          nfit=myunfold.bgintegrals[tempi];
          fiterror=myunfold.bgabserrors[tempi];
          break;
          
        //         case 1:  h=myunfold.hwjets_pos_ele; h1d=myunfold1d.hwjets_pos_ele; nfit=nwjetsposele; fiterror=nwjetsposeleerror; break;
        //         case 2:  h=myunfold.hwjets_neg_ele; h1d=myunfold1d.hwjets_neg_ele; nfit=nwjetsnegele; fiterror=nwjetsnegeleerror; break;
        //         case 3:  h=myunfold.hzjets_ele; h1d=myunfold1d.hzjets_ele; nfit=nzjetsele; fiterror=nzjetseleerror; break;
        //         case 4:  h=myunfold.hstt_ele; h1d=myunfold1d.hstt_ele; nfit=nstele; fiterror=nsteleerror; break; //it's named stt, but it's whole single t
        //         case 5:  h=myunfold.hqcd_ele; h1d=myunfold1d.hqcd_ele; nfit=nqcdele; fiterror=nqcdeleerror; break;
        //   
        //         case 7:  h=myunfold.hwjets_pos_mu; h1d=myunfold1d.hwjets_pos_mu; nfit=nwjetsposmu; fiterror=nwjetsposmuerror; break;
        //         case 8:  h=myunfold.hwjets_neg_mu; h1d=myunfold1d.hwjets_neg_mu; nfit=nwjetsnegmu; fiterror=nwjetsnegmuerror; break;
        //         case 9:  h=myunfold.hzjets_mu; h1d=myunfold1d.hzjets_mu; nfit=nzjetsmu; fiterror=nzjetsmuerror; break;
        //         case 10: h=myunfold.hstt_mu; h1d=myunfold1d.hstt_mu; nfit=nstmu; fiterror=nstmuerror; break; //it's named stt, but it's whole single t
        //         case 11: h=myunfold.hqcd_mu; h1d=myunfold1d.hqcd_mu; nfit=nqcdmu; fiterror=nqcdmuerror; break;
      }
      
      const int ndraw = random->Poisson(drawfactor * ( fiterror != 0 ? random->Gaus(nfit, fiterror) : nfit ));
      // WARNING this above assumes that fiterror scales with lumi - probably not the case, so we're overly pessimistic
      
      
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
    TH2F* hunfold1d = myunfold1d.unfoldHisto(&hrec1d, true, drawfactor); // we'll treat the pseudo experiments as if they were real data.
    
    
    // calc asy
    const double pos = hunfold->Integral(1, nbinsmafter, nbinsetaafter/2+1, nbinsetaafter);
    const double neg = hunfold->Integral(1, nbinsmafter, 1, nbinsetaafter/2);
    double asy = (pos-neg)/(pos+neg);
    /*    
    if(correct_for_lincheck)
    {
      asy -= lincheck_offset; // correct for linearity check having offset
      asy /= lincheck_slope;  // correct for linearity check not yielding slope of 1
    }
    */
    const double asy_error = asymmetryerror_afterunfolding_2d(myunfold.errormatrix, nbinsafter, pos, neg, nbinsetaafter);
    
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
    }
    
    // 1d stuff
    
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
    asys1d->Fill(asy1d);
    asyerrors1d->Fill(asy_error1d);
    asypull1d->Fill((asy1d-trueasy)/asy_error1d);
    asydiff1d->Fill(asy1d-trueasy);
  
    ndiff->Fill(pos+neg-expected_ttbar_presel);
    ndiff1d->Fill(pos1d+neg1d-expected_ttbar_presel);
    

  }

  TH1::AddDirectory(kTRUE);

  /*
  canv1->Divide(2);
  canv1->cd(1);
  asys->Draw("HIST");
  canv1->cd(2);
  asyerrors->Draw("HIST");
  */
  
  TString reweightname = "reweightedpseudo_noninv";
  
  if(zprime)
    reweightname = "reweightedpseudo_zprime";
  else if(contsquare)
    reweightname = "reweightedpseudo_contsquare";
  else
  {
    if(!inverted && !centered && negative)
      reweightname = "reweightedpseudo_neg";
    if(inverted && !centered && !negative)
      reweightname = "reweightedpseudo_inv";
    if(inverted && !centered && negative)
      reweightname = "reweightedpseudo_neginv";
    if(inverted && centered)
      reweightname = "reweightedpseudo_centered_inv";
    if(!inverted && centered)
      reweightname = "reweightedpseudo_centered_noninv";
  }
  if(donotreweight)
    reweightname = "noreweight";

  TString outname = TString("output/")+fileSuffix+"/reweightings/reweight"+rewvar+"/"+reweightname;
  /*  
  if(correct_for_lincheck)
    outname += "_corrected";
  */
  outname += ".root";
  
  TFile out(outname,"recreate");
  
  asys->Write();
  asyerrors->Write();
  asys1d->Write();
  asyerrors1d->Write();
  asypull->Write();
  asypull1d->Write();
  asydiff->Write();
  asydiff1d->Write();
  ndiff->Write();
  ndiff1d->Write();
  
  rew_selected->Write();
  rew_reco->Write();

  /// text-based output of asymmetry values
  mkdir(TString("output/")+fileSuffix+"/reweightings/",0777);
  ofstream resfile(TString("output/")+fileSuffix+"/reweightings/reweight"+rewvar+"/"+reweightname+".txt");
  resfile << reweightname << " ";
  resfile << asys->GetMean() << " " << asys1d->GetMean() << " ";
  for(int i=0; i<nbinsmafter; i++)
  {
    resfile << diffasys[i]->GetMean() << " ";
  }
  
  resfile << "\ntrue values\t";
  resfile << trueasy << " " << trueasy << " ";
  for(int i=0; i<nbinsmafter; i++)
  {
    resfile << difftrueasys[i] << " ";
  }
  
  {
    resfile << endl << endl;
    
    double asy_orig[nbinsm];
    double asy_rew[nbinsm];
    double asy_data[nbinsm];
    
    for(int i=0; i<nbinsm; i++)
    {
      TH2F* h_orig = reco_orig.toTH2F();
      TH2F* h_rew = reco_rew.toTH2F();
      TH2F* h_data = reco_data.toTH2F();
      double pos, neg;
      
      pos = h_orig->Integral(i+1, i+1, nbinseta/2+1, nbinseta);
      neg = h_orig->Integral(i+1, i+1, 1, nbinseta/2);
      asy_orig[i] = (pos-neg)/(pos+neg);
      
      pos = h_rew->Integral(i+1, i+1, nbinseta/2+1, nbinseta);
      neg = h_rew->Integral(i+1, i+1, 1, nbinseta/2);
      asy_rew[i] = (pos-neg)/(pos+neg);
      
      pos = h_data->Integral(i+1, i+1, nbinseta/2+1, nbinseta);
      neg = h_data->Integral(i+1, i+1, 1, nbinseta/2);
      asy_data[i] = (pos-neg)/(pos+neg);
          
    }
    
    resfile << "ori  ";
    for(int i=0; i<nbinsm; i++)
      resfile << asy_orig[i] << "  ";
    resfile << "\nrew  ";
    for(int i=0; i<nbinsm; i++)
      resfile << asy_rew[i] << "  ";
    resfile << "\ndat  ";
    for(int i=0; i<nbinsm; i++)
      resfile << asy_data[i] << "  ";
    
    resfile << endl;
    
    resfile << "\ndiffrew  ";
    for(int i=0; i<nbinsm; i++)
      resfile << setiosflags(ios::fixed) << setprecision(3) << asy_rew[i] - asy_orig[i] << "  ";
    
    resfile << "\ndiffdat  ";
    for(int i=0; i<nbinsm; i++)
      resfile << setiosflags(ios::fixed) << setprecision(3) << asy_data[i] - asy_orig[i] << "  ";
    
    
    // ks tests
    {
      /*
      mcks_m_ele->Add(mcks_m_mu);
      mcks_pt_ele->Add(mcks_pt_mu);
      mcks_y_ele->Add(mcks_y_mu);
      mcks_diffabsy_ele->Add(mcks_diffabsy_mu);
      */
      mcks_m_ele->Add(bg_m);
      mcks_pt_ele->Add(bg_pt);
      mcks_y_ele->Add(bg_y);
      mcks_diffabsy_ele->Add(bg_diffabsy);
      
      resfile << "\n\nKS values:  ";
      resfile << mcks_m_ele->KolmogorovTest(dataks_m) << "  ";
      resfile << mcks_pt_ele->KolmogorovTest(dataks_pt) << "  ";
      resfile << mcks_y_ele->KolmogorovTest(dataks_y) << "  ";
      resfile << mcks_diffabsy_ele->KolmogorovTest(dataks_diffabsy) << "\n";
      
      mcks_m_ele->Write("mcks_m");
      mcks_pt_ele->Write("mcks_pt");
      mcks_y_ele->Write("mcks_y");
      mcks_diffabsy_ele->Write("mcks_diffabsy");
      
      dataks_m->Write("dataks_m");
      dataks_pt->Write("dataks_pt");
      dataks_y->Write("dataks_y");
      dataks_diffabsy->Write("dataks_diffabsy");
    }
  }
  
  resfile << endl;
  resfile.close();
  
  
  
  /// draw/save plots
  for(int i=0; i<nbinsafter; i++)
  {
    reldiffs[i]->Write();
    reldiffs[i]->Draw("HIST");
    //canv1->SaveAs((string("output/")+reldiffs[i]->GetName()+".png").c_str());
    
    pulls[i]->Fit("gaus","");
    pulls[i]->Draw("HIST");
    //canv1->SaveAs((string("output/")+pulls[i]->GetName()+".png").c_str());
    pulls[i]->Write();

  }
  
  for(int i=0; i<nbinsafter1d; i++)
  {
    reldiffs1d[i]->Write();
    reldiffs1d[i]->Draw("HIST");
    //canv1->SaveAs((string("output/")+reldiffs[i]->GetName()+".png").c_str());
    
    pulls1d[i]->Fit("gaus","");
    pulls1d[i]->Draw("HIST");
    //canv1->SaveAs((string("output/")+pulls[i]->GetName()+".png").c_str());
    pulls1d[i]->Write();
  }
  
  for(int i=0; i<nbinsmafter; i++)
  {
    diffasys[i]->Write();
    diffasydiffs[i]->Write();
    diffasypulls[i]->Write();
  }
  
  
  out.Close();
  
  cout << "Finished successfully" << endl;

  return 0;
}



int main()
{
  return reweighted_pseudoexp();
}
