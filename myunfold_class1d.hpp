#ifndef MYUNFOLD_CLASS1D_HPP
#define MYUNFOLD_CLASS1D_HPP

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
#include <TNtuple.h>
#include <TUnfold.h>
#include <TUnfoldSys.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <sstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <TVector.h>
#include <TMatrix.h>
#include <TPave.h>
#include <TPaveText.h>
#include <TVirtualPad.h>
#include <TClass.h>
#include "TApplication.h"
#include "TRandom3.h"

#include "binning.hpp"

#include "cmsstyle.hpp"
#include "fitresults.hpp"
#include "helpers.hpp"
#include "samplelocations.hpp"

#include "usesystematics.hpp" // gives us whether to use systematics or not



// WARNING make sure variables in binning.hpp are set correctly


// This class is basically the same as the regular 2d MyUnfold, but it explicitly chooses
//  migmatrix1d.root and binning constants that contain only one bin for mttbar or equivalent
//  measurables.
// The intention is to be able to do both 2d unfolding and 1d unfolding at the same time
//  by using both classes.
class MyUnfold1d
{
  public:
    MyUnfold1d();
    
    TH2F* unfoldHisto(TH2F* data, bool isRealData, double bgFactor=1);
  
  
  // everything after this line should be protected, but I'm too lazy to do stuff the right way
  // - so it'll stay public and I can use the histograms in the actual programs
    TH2F* httbar_ele;
    TH2F* hwjets_pos_ele;
  //    TH2F* hwjets_neg_ele;
    TH2F* hzjets_ele;
    TH2F* hdiboson_ele;
    TH2F* hstt_ele;
  //    TH2F* hstt_top_ele;
  //    TH2F* hstt_tbar_ele;
    TH2F* hsttw_ele;
    TH2F* hsttw_top_ele;
    TH2F* hsttw_atop_ele;
  //    TH2F* hsttw_tbar_ele;
  //    TH2F* hqcd_ele;
  /*    TH2F* httbar_mu;
    TH2F* hwjets_pos_mu;
    TH2F* hwjets_neg_mu;
    TH2F* hzjets_mu;
    TH2F* hstt_mu;
    TH2F* hstt_top_mu;
    TH2F* hstt_tbar_mu;
    TH2F* hsttw_mu;
    TH2F* hsttw_top_mu;
    TH2F* hsttw_tbar_mu;
    TH2F* hqcd_mu;
  */
    TH1F* hwjets_pos_eleunwrapped;
  //    TH1F* hwjets_neg_eleunwrapped;
    TH1F* hzjets_eleunwrapped;
    TH1F* hdiboson_eleunwrapped;
    TH1F* hst_eleunwrapped;
  //    TH1F* hqcd_eleunwrapped;
  /*    TH1F* hwjets_pos_muunwrapped;
    TH1F* hwjets_neg_muunwrapped;
    TH1F* hzjets_muunwrapped;
    TH1F* hst_muunwrapped;
    TH1F* hqcd_muunwrapped;
  */
    TH2F* migmatrix;
    TH1F* preNunwrapped;
    
    TH2D* correlationmatrix;
    TH2D* errormatrix;
    TGraph* lcurve;
    TH2F* hunfold;
    TH1D* hunfoldunwrapped;
    
    double lastTau;
    double seleff;
    double scaleBias;
    double expectedNUnfolded;  
    

    void minimizeRhoAverage(TUnfoldSys* unfold, TH1F* hdata, int nsteps, double log10min, double log10max);
    
    std::vector<double> bgrelerrors;
    std::vector<double> bgabserrors;
    std::vector<double> bgintegrals;
    std::vector<TH1F*> bglist;
    std::vector<TH2F*> bgwrappedlist;
    std::vector<TH2F*> bgwrappedlist_simplebins;
    
};

MyUnfold1d::MyUnfold1d()
  : httbar_ele(NULL),
    hwjets_pos_ele(NULL),
    //    hwjets_neg_ele(NULL),
    hzjets_ele(NULL),
    hdiboson_ele(NULL),
    hstt_ele(NULL),
    //hstt_top_ele(NULL),
    //    hstt_tbar_ele(NULL),
    hsttw_ele(NULL),
    hsttw_top_ele(NULL),
    hsttw_atop_ele(NULL),
    //    hsttw_tbar_ele(NULL),
    //    hqcd_ele(NULL),
    /*    httbar_mu(NULL),
    hwjets_pos_mu(NULL),
    hwjets_neg_mu(NULL),
    hzjets_mu(NULL),
    hstt_mu(NULL),
    hstt_top_mu(NULL),
    hstt_tbar_mu(NULL),
    hsttw_mu(NULL),
    hsttw_top_mu(NULL),
    hsttw_tbar_mu(NULL),
    hqcd_mu(NULL),
    */
    preNunwrapped(NULL),
    correlationmatrix(NULL),
    errormatrix(NULL),
    lcurve(NULL),
    hunfoldunwrapped(NULL)
{
  // files containing mc information
  #ifdef POWHEG
    TFile* file_ttbar_ele;
    file_ttbar_ele = new TFile(sl_ttbarele);
  #endif
  TFile* file_wjets_pos_ele = new TFile(sl_wjetsposele);
  //  TFile* file_wjets_neg_ele = new TFile(sl_wjetsnegele);
  TFile* file_zjets_ele = new TFile(sl_zjetsele);
  TFile* file_diboson_ele = new TFile(sl_dibosonele);
  //  TFile* file_stt_top_ele = new TFile(sl_singlet_tchan_t_ele);
  //  TFile* file_stt_tbar_ele = new TFile(sl_singlet_tchan_tbar_ele);
  TFile* file_sttw_top_ele = new TFile(sl_singlet_twchan_t_ele);
  TFile* file_sttw_atop_ele = new TFile(sl_singlet_twchan_tbar_ele);
  //  TFile* file_sttw_tbar_ele = new TFile(sl_singlet_twchan_tbar_ele);
  //  TFile* file_qcd_ele = new TFile(sl_qcdele);
  /*  
  #ifdef POWHEG
    TFile* file_ttbar_mu;
    file_ttbar_mu = new TFile(sl_ttbarmu);
  #endif
  TFile* file_wjets_pos_mu = new TFile(sl_wjetsposmu);
  TFile* file_wjets_neg_mu = new TFile(sl_wjetsnegmu);
  TFile* file_zjets_mu = new TFile(sl_zjetsmu);
  TFile* file_stt_top_mu = new TFile(sl_singlet_tchan_t_mu);
  TFile* file_stt_tbar_mu = new TFile(sl_singlet_tchan_tbar_mu);
  TFile* file_sttw_top_mu = new TFile(sl_singlet_twchan_t_mu);
  TFile* file_sttw_tbar_mu = new TFile(sl_singlet_twchan_tbar_mu);
  TFile* file_qcd_mu = new TFile(sl_qcdmu);
  */  
  // tuples for filling histos
//  TString tuplename("TTbarSelHyppsiprime/BBB_UnfoldingTuple");
  TString tuplename("Tuple");
  
  TNtuple* ttbartuple_ele = (TNtuple*)file_ttbar_ele->Get(tuplename);
  TNtuple* wjetspostuple_ele = (TNtuple*)file_wjets_pos_ele->Get(tuplename);
  //  TNtuple* wjetsnegtuple_ele = (TNtuple*)file_wjets_neg_ele->Get(tuplename);
  TNtuple* zjetstuple_ele = (TNtuple*)file_zjets_ele->Get(tuplename);
  TNtuple* dibosontuple_ele = (TNtuple*)file_diboson_ele->Get(tuplename);
  //  TNtuple* stt_top_tuple_ele = (TNtuple*)file_stt_top_ele->Get(tuplename);
  //  TNtuple* stt_tbar_tuple_ele = (TNtuple*)file_stt_tbar_ele->Get(tuplename);
  TNtuple* sttw_top_tuple_ele = (TNtuple*)file_sttw_top_ele->Get(tuplename);
  TNtuple* sttw_atop_tuple_ele = (TNtuple*)file_sttw_atop_ele->Get(tuplename);
  //  TNtuple* sttw_tbar_tuple_ele = (TNtuple*)file_sttw_tbar_ele->Get(tuplename);
  //  TNtuple* qcdtuple_ele = (TNtuple*)file_qcd_ele->Get(tuplename);
  /*  
  TNtuple* ttbartuple_mu = (TNtuple*)file_ttbar_mu->Get(tuplename);
  TNtuple* wjetspostuple_mu = (TNtuple*)file_wjets_pos_mu->Get(tuplename);
  TNtuple* wjetsnegtuple_mu = (TNtuple*)file_wjets_neg_mu->Get(tuplename);
  TNtuple* zjetstuple_mu = (TNtuple*)file_zjets_mu->Get(tuplename);
  TNtuple* stt_top_tuple_mu = (TNtuple*)file_stt_top_mu->Get(tuplename);
  TNtuple* stt_tbar_tuple_mu = (TNtuple*)file_stt_tbar_mu->Get(tuplename);
  TNtuple* sttw_top_tuple_mu = (TNtuple*)file_sttw_top_mu->Get(tuplename);
  TNtuple* sttw_tbar_tuple_mu = (TNtuple*)file_sttw_tbar_mu->Get(tuplename);
  TNtuple* qcdtuple_mu = (TNtuple*)file_qcd_mu->Get(tuplename);
  */  
  // mc histograms, to be filled out of tuples
  httbar_ele = new TH2F("httbar_ele","httbar_ele", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  hwjets_pos_ele = new TH2F("hwjets_pos_ele","hwjets_pos_ele", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  //  hwjets_neg_ele = new TH2F("hwjets_neg_ele","hwjets_neg_ele", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  hzjets_ele = new TH2F("hzjets_ele","hzjets_ele", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  //  hstt_top_ele = new TH2F("hstt_top_ele","hstt_top_ele", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  //  hstt_tbar_ele = new TH2F("hstt_tbar_ele","hstt_tbar_ele", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  hsttw_top_ele = new TH2F("hsttw_top_ele","hsttw_top_ele", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  hsttw_atop_ele = new TH2F("hsttw_atop_ele","hsttw_atop_ele", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  hdiboson_ele = new TH2F("hdiboson_ele","hdiboson_ele", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  //  hsttw_tbar_ele = new TH2F("hsttw_tbar_ele","hsttw_tbar_ele", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  //  hqcd_ele = new TH2F("hqcd_ele","hqcd_ele", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  /*
  httbar_mu = new TH2F("httbar_mu","httbar_mu", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  hwjets_pos_mu = new TH2F("hwjets_pos_mu","hwjets_pos_mu", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  hwjets_neg_mu = new TH2F("hwjets_neg_mu","hwjets_neg_mu", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  hzjets_mu = new TH2F("hzjets_mu","hzjets_mu", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  hstt_top_mu = new TH2F("hstt_top_mu","hstt_top_mu", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  hstt_tbar_mu = new TH2F("hstt_tbar_mu","hstt_tbar_mu", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  hsttw_top_mu = new TH2F("hsttw_top_mu","hsttw_top_mu", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  hsttw_tbar_mu = new TH2F("hsttw_tbar_mu","hsttw_tbar_mu", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  hqcd_mu = new TH2F("hqcd_mu","hqcd_mu", nbinsm1d, mbinedges1d, nbinseta1d, etabinedges1d);
  */
  // result histo
  hunfold = new TH2F("hunfold","hunfold", nbinsmafter1d, mbinedgesafter1d, nbinsetaafter1d, etabinedges1dafter);
  
  // fill histos for mc samples
  float mttbarrec, diffabsetarec, weight;
  
  for(int i=0; i<6; i++)
  {
    TNtuple* t = 0;
    TH2F* h = 0; 
    double nfit = 0;
    
    if(i==0) {t=ttbartuple_ele; h=httbar_ele; }
    if(i==1) {t=wjetspostuple_ele; h=hwjets_pos_ele; }
    //    if(i==2) {t=wjetsnegtuple_ele; h=hwjets_neg_ele; }
    if(i==2) {t=zjetstuple_ele; h=hzjets_ele; }
    //    if(i==4) {t=stt_top_tuple_ele; h=hstt_top_ele; }
    //    if(i==5) {t=stt_tbar_tuple_ele; h=hstt_tbar_ele; }
    if(i==3) {t=sttw_top_tuple_ele; h=hsttw_top_ele; }
    if(i==4) {t=sttw_atop_tuple_ele; h=hsttw_atop_ele; }
    if(i==5) {t=dibosontuple_ele; h=hdiboson_ele; }
    //    if(i==4) {t=sttw_tbar_tuple_ele; h=hsttw_tbar_ele; }
    //    if(i==8) {t=qcdtuple_ele; h=hqcd_ele; }
    /*    
    if(i==9) {t=ttbartuple_mu; h=httbar_mu; }
    if(i==10) {t=wjetspostuple_mu; h=hwjets_pos_mu; }
    if(i==11) {t=wjetsnegtuple_mu; h=hwjets_neg_mu; }
    if(i==12) {t=zjetstuple_mu; h=hzjets_mu; }
    if(i==13) {t=stt_top_tuple_mu; h=hstt_top_mu; }
    if(i==14) {t=stt_tbar_tuple_mu; h=hstt_tbar_mu; }
    if(i==15) {t=sttw_top_tuple_mu; h=hsttw_top_mu; }
    if(i==16) {t=sttw_tbar_tuple_mu; h=hsttw_tbar_mu; }
    if(i==17) {t=qcdtuple_mu; h=hqcd_mu; }
    */
    
    t->SetBranchAddress(nameOfXAxisVarRec, &mttbarrec);
    t->SetBranchAddress(nameOfSensVarRec, &diffabsetarec);
    t->SetBranchAddress("weight", &weight);
    weight *= lumi_factor;
    
    const int nentries = t->GetEntries();
    
    const bool is_ttbar = (t == ttbartuple_ele) /*|| (t == ttbartuple_mu)*/;
    int nentries_toskip = 0; 
    int nentries_toconsider = nentries;
    

    
    for(int j=nentries_toskip; j<nentries_toconsider; j++)
    {
        
      t->GetEntry(j);
      
      fill_nooverflow_2d(h, fabs(mttbarrec), diffabsetarec, weight);
    }
    
  }

  // combine singletop channels
  /*
  scale_to(hstt_top_ele, pred_pos_st_t_top_ele_met + pred_neg_st_t_top_ele_met);
  scale_to(hstt_tbar_ele, pred_pos_st_t_atop_ele_met + pred_neg_st_t_atop_ele_met);
  hstt_ele = (TH2F*)hstt_top_ele->Clone("hstt_ele");
  hstt_ele->Add(hstt_tbar_ele); 
  */

  /*
  scale_to(hstt_top_mu, pred_pos_st_t_top_mu_met + pred_neg_st_t_top_mu_met);
  scale_to(hstt_tbar_mu, pred_pos_st_t_atop_mu_met + pred_neg_st_t_atop_mu_met);
  hstt_mu = (TH2F*)hstt_top_mu->Clone("hstt_mu");
  hstt_mu->Add(hstt_tbar_mu); 
  */  

  hsttw_top_ele->Add(hsttw_atop_ele); 
  scale_to(hsttw_top_ele, pred_pos_st_tw_top_ele_met /*+ pred_neg_st_tw_top_ele_met*/);
  //scale_to(hsttw_tbar_ele, pred_pos_st_tw_atop_ele_met /*+ pred_neg_st_tw_atop_ele_met*/);
  hstt_ele = (TH2F*)hsttw_top_ele->Clone("hstt_ele");
  //hstt_ele->Add(hsttw_tbar_ele); 

  /*
  scale_to(hsttw_top_mu, pred_pos_st_tw_top_mu_met + pred_neg_st_tw_top_mu_met);
  scale_to(hsttw_tbar_mu, pred_pos_st_tw_atop_mu_met + pred_neg_st_tw_atop_mu_met);
  hsttw_mu = (TH2F*)hsttw_top_mu->Clone("hsttw_mu");
  hsttw_mu->Add(hsttw_tbar_mu); 
  */
  //  hstt_ele->Add(hsttw_ele);
  //  hstt_mu->Add(hsttw_mu);

  
  
  // unwrap bg histograms to 1d
  hwjets_pos_eleunwrapped = new TH1F("hwjets_pos_eleunwrapped","hwjets_pos_eleunwrapped", nbins1d, 0.5, nbins1d+0.5); 
  unwrap2dhisto(hwjets_pos_ele, hwjets_pos_eleunwrapped);
  //  hwjets_neg_eleunwrapped = new TH1F("hwjets_neg_eleunwrapped","hwjets_neg_eleunwrapped", nbins1d, 0.5, nbins1d+0.5); 
  //  unwrap2dhisto(hwjets_neg_ele, hwjets_neg_eleunwrapped);
  hzjets_eleunwrapped = new TH1F("hzjets_eleunwrapped","hzjets_eleunwrapped", nbins1d, 0.5, nbins1d+0.5); 
  unwrap2dhisto(hzjets_ele, hzjets_eleunwrapped);
  hst_eleunwrapped = new TH1F("hst_eleunwrapped","hst_eleunwrapped", nbins1d, 0.5, nbins1d+0.5); 
  unwrap2dhisto(hstt_ele, hst_eleunwrapped);  
  hdiboson_eleunwrapped = new TH1F("hdiboson_eleunwrapped","hdiboson_eleunwrapped", nbins1d, 0.5, nbins1d+0.5); 
  unwrap2dhisto(hdiboson_ele, hdiboson_eleunwrapped);
  //  hqcd_eleunwrapped = new TH1F("hqcd_eleunwrapped","hqcd_eleunwrapped", nbins1d, 0.5, nbins1d+0.5); 
  //  unwrap2dhisto(hqcd_ele, hqcd_eleunwrapped);
  /*  hwjets_pos_muunwrapped = new TH1F("hwjets_pos_muunwrapped","hwjets_pos_muunwrapped", nbins1d, 0.5, nbins1d+0.5); 
  unwrap2dhisto(hwjets_pos_mu, hwjets_pos_muunwrapped);
  hwjets_neg_muunwrapped = new TH1F("hwjets_neg_muunwrapped","hwjets_neg_muunwrapped", nbins1d, 0.5, nbins1d+0.5); 
  unwrap2dhisto(hwjets_neg_mu, hwjets_neg_muunwrapped);
  hzjets_muunwrapped = new TH1F("hzjets_muunwrapped","hzjets_muunwrapped", nbins1d, 0.5, nbins1d+0.5); 
  unwrap2dhisto(hzjets_mu, hzjets_muunwrapped);
  hst_muunwrapped = new TH1F("hst_muunwrapped","hst_muunwrapped", nbins1d, 0.5, nbins1d+0.5); 
  unwrap2dhisto(hstt_mu, hst_muunwrapped);  
  hqcd_muunwrapped = new TH1F("hqcd_muunwrapped","hqcd_muunwrapped", nbins1d, 0.5, nbins1d+0.5); 
  unwrap2dhisto(hqcd_mu, hqcd_muunwrapped);
  */
  
  // get migration matrixfile
//  TFile* matrixfile = new TFile(TString("output/")+fileSuffix+"/1dmigmatrix.root");
  TFile* matrixfile = new TFile(TString("output/")+"1dmigmatrix.root");
  migmatrix = (TH2F*)matrixfile->Get("migmatrixcombined");
  
  seleff = migmatrix->Integral(1, nbinsafter1d, 1, nbins1d)
            / migmatrix->Integral(0, nbinsafter1d+1, 0, nbins1d+1);
  
  // get bias distribution
//  TFile* selefffile = new TFile(TString("output/")+fileSuffix+"/1dseleff.root");
  TFile* selefffile = new TFile(TString("output/")+"1dseleff.root");
  preNunwrapped = new TH1F("preNunwrapped", "preNunwrapped", nbinsafter1d, 0.5, nbinsafter1d+0.5);
  unwrap2dhisto((TH2F*)selefffile->Get("preN"), preNunwrapped);
  

  /// create rotated background using cov matrix 
  // WARNING needs changing when bg fit variables change in any way
  cerr << "FIXME bgfactor is ignored now in some parts of the code, not in others" << endl;
  const double bgFactor = 1;
  
  std::vector<TH1F*> bglist_nonrot;
  if(!swapwjets)
  {
    bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hwjets_pos_eleunwrapped->Clone(), bgFactor * pred_pos_wjets_ele_met * betawjetsposele));
    //    bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hwjets_neg_eleunwrapped->Clone(), bgFactor * pred_neg_wjets_ele_met * betawjetsnegele));
    //    bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hwjets_pos_muunwrapped->Clone(), bgFactor * pred_pos_wjets_mu_met * betawjetsposmu));
    //    bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hwjets_neg_muunwrapped->Clone(), bgFactor * pred_neg_wjets_mu_met * betawjetsnegmu));
  }
  else
  {
    //    bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hwjets_neg_eleunwrapped->Clone(), bgFactor * pred_pos_wjets_ele_met * betawjetsposele));
//    bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hwjets_pos_eleunwrapped->Clone(), bgFactor * pred_neg_wjets_ele_met * betawjetsnegele));
    //    bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hwjets_neg_muunwrapped->Clone(), bgFactor * pred_pos_wjets_mu_met * betawjetsposmu));
    //    bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hwjets_pos_muunwrapped->Clone(), bgFactor * pred_neg_wjets_mu_met * betawjetsnegmu));
  }
  bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hzjets_eleunwrapped->Clone(), bgFactor * nzjetsele));
  bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hst_eleunwrapped->Clone(), bgFactor * nstele));
  bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hdiboson_eleunwrapped->Clone(), bgFactor * ndibosonele));
  //  bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hzjets_muunwrapped->Clone(), bgFactor * nzjetsmu));
  //  bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hst_muunwrapped->Clone(), bgFactor * nstmu));
  //  bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hqcd_eleunwrapped->Clone(), bgFactor * nqcdele));
  //  bglist_nonrot.push_back((TH1F*)scale_to((TH1F*)hqcd_muunwrapped->Clone(), bgFactor * nqcdmu));
  
  std::vector<double> bgpredlist;
  bgpredlist.push_back(bgFactor * pred_pos_wjets_ele_met);
  //  bgpredlist.push_back(bgFactor * pred_neg_wjets_ele_met);
  //  bgpredlist.push_back(bgFactor * pred_pos_wjets_mu_met);
  //  bgpredlist.push_back(bgFactor * pred_neg_wjets_mu_met);
  bgpredlist.push_back(bgFactor * (pred_pos_zjets_ele_met /*+ pred_neg_zjets_ele_met*/));
  bgpredlist.push_back(bgFactor * pred_st_ele_met);
  bgpredlist.push_back(bgFactor * (pred_pos_diboson_ele_met /*+ pred_neg_zjets_ele_met*/));
  //  bgpredlist.push_back(bgFactor * (pred_pos_zjets_mu_met + pred_neg_zjets_mu_met));
  //  bgpredlist.push_back(bgFactor * pred_st_mu_met);
  //  bgpredlist.push_back(bgFactor * (pred_pos_qcd_ele_met + pred_neg_qcd_ele_met));
  //  bgpredlist.push_back(bgFactor * (pred_pos_qcd_mu_met + pred_neg_qcd_mu_met));
  

//  TFile* covfile = new TFile("covarianz_matrix.root");
//  TH1D* covmatrix_flat = (TH1D*) covfile->Get("covarianz_matrix0");
  
  
  const int bg_process_n = bglist_nonrot.size();
//  TMatrixD covmatrix(bg_process_n,bg_process_n);
  
  // starting both at 2 so that ttbar is taken out
/*
  for(int i=2; i<bg_process_n+2; i++)
  {    
    for(int j=2; j<bg_process_n+2; j++)
    {
      covmatrix[i-2][j-2] = covmatrix_flat->GetBinContent(i*(bg_process_n+2)+j+1);      
      covmatrix[i-2][j-2] *= bgpredlist[i-2] * bgpredlist[j-2]; 
    }    
  }
  
  TVectorD eigenvalues(bg_process_n);
  TMatrixD eigenvectors = covmatrix.EigenVectors(eigenvalues);
*/  
  for(int i=0; i<bg_process_n; i++)
  {    
    TH1F* current = (TH1F*)bglist_nonrot[i]->Clone();
    current->Reset();
/*    
    for(int j=0; j<bg_process_n; j++)
    {
//      current->Add(bglist_nonrot[j], pow(eigenvectors(j,i),2));
      current->Add(bglist_nonrot[j], 1.);
    }
*/
    current->Add(bglist_nonrot[i]);
   
    bglist.push_back(current);
/*
    bgrelerrors.push_back(sqrt(eigenvalues[i]) / current->Integral());
    bgabserrors.push_back(sqrt(eigenvalues[i]));
    bgintegrals.push_back(current->Integral());
*/

    bgrelerrors.push_back(1. / sqrt(current->Integral()));
    bgabserrors.push_back(sqrt(current->Integral()));
    bgintegrals.push_back(current->Integral());

    
    // save wrapped version as well, for bg subtraction etc
    TH2F* wrapped = (TH2F*) httbar_ele->Clone();
    TH2F* wrapped_simplebins = new TH2F("","",nbinsm1d, 0.5, nbinsm1d+0.5, nbinseta1d, 0.5, nbinseta1d+0.5);
    int n = 0;
    for(int x=1; x<=wrapped->GetXaxis()->GetNbins(); x++)
    {
      for(int y=1; y<=wrapped->GetYaxis()->GetNbins(); y++)
      {
        n++;
        wrapped->SetBinContent(x, y, current->GetBinContent(n));
        wrapped->SetBinError(x, y, current->GetBinError(n));
        wrapped_simplebins->SetBinContent(x, y, current->GetBinContent(n));
        wrapped_simplebins->SetBinError(x, y, current->GetBinError(n));
      }
    }
    bgwrappedlist.push_back(wrapped);
    bgwrappedlist_simplebins.push_back(wrapped_simplebins);
  }
  
}



TH2F* MyUnfold1d::unfoldHisto(TH2F* hdata, bool isRealData, double bgFactor)
{
  // concatenate rows to make 1d histo out of the 2d histo
  TH1F* hdataunwrapped = new TH1F("hdataunwrapped","hdataunwrapped", nbins1d, 0.5, nbins1d+0.5); 
  unwrap2dhisto(hdata, hdataunwrapped);
    
  /// unfold the sucker
  TH2F* migmat = migmatrix;
  
  // check with random migmatrix variations
  if(vary_migmatrix_for_each_unfolding)
  {
    migmat = (TH2F*) migmatrix->Clone();
    for(int x=1; x<nbinsafter1d; x++)
      for(int y=1; y<nbins1d; y++)
      {
        // if bin error is zero we use average weight of events 1.1 as error
//        migmat->SetBinContent(x,y, gRandom->Gaus(migmatrix->GetBinContent(x,y),migmatrix->GetBinError(x,y)?migmatrix->GetBinError(x,y):1.1));
        migmat->SetBinContent(x,y, gRandom->Gaus(migmatrix->GetBinContent(x,y),migmatrix->GetBinError(x,y)?migmatrix->GetBinError(x,y):0.01325));
      }
  }  
  
  TUnfoldSys unfold(migmat, TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeNone);
  
  unfold.SetInput(hdataunwrapped); 
  
  // set bias
  if(use_bias_distribution)
    unfold.SetBias(preNunwrapped); 
  
  if(isRealData)
//    expectedNUnfolded = (hdataunwrapped->Integral() /*- (n_bg_ele*/ /*+n_bg_mu*/ /*)*bgFactor */) / seleff;
    expectedNUnfolded = (hdataunwrapped->Integral() - (n_bg_ele /*+n_bg_mu*/ )*bgFactor ) / seleff;
  else
    expectedNUnfolded = hdataunwrapped->Integral() / seleff; // object attribute! will be used outside this method.
  
  scaleBias = expectedNUnfolded / preNunwrapped->Integral(); // object attribute! will be used outside this method.
  
  
  // bg subtraction  
  if(isRealData)
  {
    for(int i=0; i<bglist.size(); i++)
    {
      unfold.SubtractBackground(bglist[i], TString::Format("rotated bg%d", i).Data(), 1, bgrelerrors[i]);  
//      unfold.SubtractBackground(bglist[i], TString::Format("rotated bg%d", i).Data(), 1, bgrelerrors[i]);  
    }

  }
  
  unfold.RegularizeBins2D(1,1,nbinsetaafter1d,nbinsetaafter1d,nbinsmafter1d, TUnfold::kRegModeCurvature); 
  
  const double tau = pow(10, logtau_for_binning_1d);
  
  #ifndef FINDTAU
    unfold.DoUnfold(tau, hdataunwrapped, scaleBias);
  #else
    minimizeRhoAverage(&unfold, hdataunwrapped, 200, -6, -2);
  #endif
  
  lastTau = unfold.GetTau();
  
  delete hunfoldunwrapped;
  hunfoldunwrapped = unfold.GetOutput("hunfoldunwrapped", "unfolded", 0.5, nbinsafter1d+0.5);

  delete correlationmatrix;
  correlationmatrix = unfold.GetRhoIJ("correlation", "correlation matrix");
  
  delete errormatrix;
  errormatrix = unfold.GetEmatrix("errormatrix", "error matrix"); // only error types d and e... FIXME
  
  /// extract 2d histogram out of unfolded unwrapped histo

  int n3 = 0;
  for(int x=1; x<=hunfold->GetXaxis()->GetNbins(); x++)
  {
    for(int y=1; y<=hunfold->GetYaxis()->GetNbins(); y++)
    {
      n3++;
      hunfold->SetBinContent(x, y, hunfoldunwrapped->GetBinContent(n3));
      hunfold->SetBinError(x, y, hunfoldunwrapped->GetBinError(n3));
    }
  }
  
  delete hdataunwrapped;

  
  return hunfold;
}


void MyUnfold1d::minimizeRhoAverage(TUnfoldSys* unfold, TH1F* hdata, int nsteps, double log10min, double log10max)
{
  double bestlogtau = -1000;
  double bestavg = 1e5;
  
  const double step = (log10max-log10min)/nsteps;
  
  for(double logtau=log10min; logtau<=log10max; logtau+=step)
  {
    unfold->DoUnfold(pow(10, logtau), hdata, scaleBias);    
    
    const double avg = unfold->GetRhoAvg();
    
    if(avg < bestavg)
    {
      bestavg = avg;
      bestlogtau = logtau;
    }
    
  }
  
  unfold->DoUnfold(pow(10, bestlogtau), hdata, scaleBias); 
  
}


#endif // MYUNFOLD_CLASS1D_HPP
