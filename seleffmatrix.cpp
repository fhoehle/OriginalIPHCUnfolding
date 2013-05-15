#include <iostream>
#include <assert.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>
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


#include "myh2/myh2.hpp"

#include "binning.hpp"

// WARNING make sure variables in binning.hpp are set correctly


int seleffmatrix()
{
  TH1::SetDefaultSumw2(true);

  CMSStyle* cmsstyle = new CMSStyle();
  gStyle->SetOptStat("");
  TCanvas* canv1 = new TCanvas("canv1","canv1",300,0,800,600);
  canv1->SetRightMargin(0.04);
  
  // to make z-axis stuff show completely in 2d plots
  canv1->SetRightMargin(0.14);
  
  H2gen* seleff;
  H2gen* seleffele;
  //  H2gen* seleffmu;
  TH2F* selefffine;
  H2gen* nonselected;
  H2gen* nonselectedtotal;
  H2gen* nonselectedele;
  //  H2gen* nonselectedmu;
  H2gen* preN = new H2gen("preN","preN", xedgesgen, yedgesgen); 
  
  // these are going to be reweighted in the sensitive variable - for the lincheck
  H2gen* preNkm05 = new H2gen("preNkm05","preNkm05", xedgesgen, yedgesgen); 
  H2gen* preNkm10 = new H2gen("preNkm10","preNkm10", xedgesgen, yedgesgen); 
  H2gen* preNkm15 = new H2gen("preNkm15","preNkm15", xedgesgen, yedgesgen); 
  H2gen* preNkm20 = new H2gen("preNkm20","preNkm20", xedgesgen, yedgesgen); 
  H2gen* preNkm25 = new H2gen("preNkm25","preNkm25", xedgesgen, yedgesgen); 
  H2gen* preNkp05 = new H2gen("preNkp05","preNkp05", xedgesgen, yedgesgen); 
  H2gen* preNkp10 = new H2gen("preNkp10","preNkp10", xedgesgen, yedgesgen); 
  H2gen* preNkp15 = new H2gen("preNkp15","preNkp15", xedgesgen, yedgesgen); 
  H2gen* preNkp20 = new H2gen("preNkp20","preNkp20", xedgesgen, yedgesgen); 
  H2gen* preNkp25 = new H2gen("preNkp25","preNkp25", xedgesgen, yedgesgen); 
  
  // these are samples with a gradient reweighting
  H2gen* preNkcont_m = new H2gen("preNkcont_m","preNkcont_m", xedgesgen, yedgesgen); 
  H2gen* preNkcontinv_m = new H2gen("preNkcontinv_m","preNkcontinv_m", xedgesgen, yedgesgen); 
  H2gen* preNkcontneg_m = new H2gen("preNkcontneg_m","preNkcontneg_m", xedgesgen, yedgesgen); 
  H2gen* preNkcontneginv_m = new H2gen("preNkcontneginv_m","preNkcontneginv_m", xedgesgen, yedgesgen); 
  H2gen* preNkcontcent_m = new H2gen("preNkcontcent_m","preNkcontcent_m", xedgesgen, yedgesgen); 
  H2gen* preNkcontcentinv_m = new H2gen("preNkcontcentinv_m","preNkcontcentinv_m", xedgesgen, yedgesgen); 
  H2gen* preNkcontsquare_m = new H2gen("preNkcontsquare_m","preNkcontsquare_m", xedgesgen, yedgesgen); 
  
  H2gen* preNkcont_pt = new H2gen("preNkcont_pt","preNkcont_pt", xedgesgen, yedgesgen); 
  H2gen* preNkcontinv_pt = new H2gen("preNkcontinv_pt","preNkcontinv_pt", xedgesgen, yedgesgen); 
  H2gen* preNkcontneg_pt = new H2gen("preNkcontneg_pt","preNkcontneg_pt", xedgesgen, yedgesgen); 
  H2gen* preNkcontneginv_pt = new H2gen("preNkcontneginv_pt","preNkcontneginv_pt", xedgesgen, yedgesgen); 
  H2gen* preNkcontcent_pt = new H2gen("preNkcontcent_pt","preNkcontcent_pt", xedgesgen, yedgesgen); 
  H2gen* preNkcontcentinv_pt = new H2gen("preNkcontcentinv_pt","preNkcontcentinv_pt", xedgesgen, yedgesgen); 
  H2gen* preNkcontsquare_pt = new H2gen("preNkcontsquare_pt","preNkcontsquare_pt", xedgesgen, yedgesgen); 
  
  H2gen* preNkcont_y = new H2gen("preNkcont_y","preNkcont_y", xedgesgen, yedgesgen); 
  H2gen* preNkcontinv_y = new H2gen("preNkcontinv_y","preNkcontinv_y", xedgesgen, yedgesgen); 
  H2gen* preNkcontneg_y = new H2gen("preNkcontneg_y","preNkcontneg_y", xedgesgen, yedgesgen); 
  H2gen* preNkcontneginv_y = new H2gen("preNkcontneginv_y","preNkcontneginv_y", xedgesgen, yedgesgen); 
  H2gen* preNkcontcent_y = new H2gen("preNkcontcent_y","preNkcontcent_y", xedgesgen, yedgesgen); 
  H2gen* preNkcontcentinv_y = new H2gen("preNkcontcentinv_y","preNkcontcentinv_y", xedgesgen, yedgesgen); 
  H2gen* preNkcontsquare_y = new H2gen("preNkcontsquare_y","preNkcontsquare_y", xedgesgen, yedgesgen); 
  
  
  H2gen* postNges = new H2gen("postNges","postNges", xedgesgen, yedgesgen);
  H2gen* postNele = new H2gen("postNele","postNele", xedgesgen, yedgesgen);
  //  H2gen* postNmu = new H2gen("postNmu","postNmu", xedgesgen, yedgesgen);

  #ifdef BUILD_1D_BASIS
    TH1F* preN1d = new TH1F("preN1d","preN1d", nbinsetaafter1d, etabinedges1dafter); 
    TH1F* postN1d = new TH1F("postN1d","postN1d", nbinsetaafter1d, etabinedges1dafter); 
  #endif
  
  TH2F* preselected_nonflat = new TH2F("preselected_nonflat","preselectednonflat", 100, xedgesgen[0], xedgesgen[nbinsmafter], 100, yedgesgen[0][0], yedgesgen[0][nbinsetaafter]); 
  
  TH2F* preselected_nonflat_fewbins = new TH2F("preselected_nonflat_fewbins","preselectednonflat_fewbins", 8, xedgesgen[0], xedgesgen[nbinsmafter], 10, yedgesgen[0][0], yedgesgen[0][nbinsetaafter]); 
  TH1F* preselected_nonflat_asys = new TH1F("preselected_nonflat_asys","preselectednonflat_asys", 8, xedgesgen[0], xedgesgen[nbinsmafter]); 
  
  TH2F* preNfine = new TH2F("preNfine","preNfine", 25, xedgesgen[0], xedgesgen[nbinsmafter], 20, yedgesgen[0][0], yedgesgen[0][nbinsetaafter]); 
  TH2F* postNfine = new TH2F("postNfine","postNfine", 25, xedgesgen[0], xedgesgen[nbinsmafter], 20, yedgesgen[0][0], yedgesgen[0][nbinsetaafter]); 
  
  TH1F* sensgenselplots[nbinsmafter];
  for(int i=0; i<nbinsmafter; i++)
  {
    sensgenselplots[i] = new TH1F("","", 100, yedgesgen[0][0], yedgesgen[0][nbinsetaafter]);
  }
  TH1F* senspreselplotinc = new TH1F("","", 100, yedgesgen[0][0], yedgesgen[0][nbinsetaafter]);
  TH1F* senspreselplots[nbinsmafter];
  for(int i=0; i<nbinsmafter; i++)
  {
    senspreselplots[i] = new TH1F("","", 100, yedgesgen[0][0], yedgesgen[0][nbinsetaafter]);
  }
  
  TH2F* preselcorrelationplot = new TH2F("","", 100, xedgesgen[0], xedgesgen[nbinsmafter], 100, yedgesgen[0][0], yedgesgen[0][nbinsetaafter]);
  TH2F* genselcorrelationplot = new TH2F("","", 100, xedgesgen[0], xedgesgen[nbinsmafter], 100, yedgesgen[0][0], yedgesgen[0][nbinsetaafter]);
  
  TH1F* kinematicpreselplot = new TH1F("kinematicpreselplot","kinematicpreselplot", 100, pretty_kin_loweredge, pretty_kin_upperedge);
  
  
  float *mttbargen, diffabsy_gen, weight;
  float diffysquare_gen, pt_ttbar_gen, pz_ttbar_gen;
  float eta_ttbar_gen, y_ttbar_gen;
      

  
  int nentries, nentries_toconsider, nentries_toskip;


  #ifdef POWHEG
    int nPreSamples = 1;
  #endif

  /// fill pre histo
  for(int j=0; j<nPreSamples; j++)
  {
    TFile* preselfile;
    
    cout << "Starting with presel file #" << j << endl;
    
    double sampleweight = 1;
    

    #ifdef POWHEG
      if(j==0)  preselfile = new TFile(sl_preselfile0);
/*      if(j==1)  preselfile = new TFile(sl_preselfile1);
      if(j==2)  preselfile = new TFile(sl_preselfile2);
      if(j==3)  preselfile = new TFile(sl_preselfile3);
      if(j==4)  preselfile = new TFile(sl_preselfile4);
      if(j==5)  preselfile = new TFile(sl_preselfile5);
      if(j==6)  preselfile = new TFile(sl_preselfile6);
      if(j==7)  preselfile = new TFile(sl_preselfile7);
      if(j==8)  preselfile = new TFile(sl_preselfile8);
      if(j==9)  preselfile = new TFile(sl_preselfile9);
      if(j==10)  preselfile = new TFile(sl_preselfile10);
      if(j==11)  preselfile = new TFile(sl_preselfile11);
      if(j==12)  preselfile = new TFile(sl_preselfile12);
      if(j==13)  preselfile = new TFile(sl_preselfile13);
      if(j==14)  preselfile = new TFile(sl_preselfile14);
      if(j==15)  preselfile = new TFile(sl_preselfile15);
      if(j==16)  preselfile = new TFile(sl_preselfile16);    
*/
    #endif

    
//    TNtuple* pretuple = (TNtuple*)preselfile->Get("DataTuple"); 
//    TNtuple* pretuple = (TNtuple*)preselfile->Get("TTbarSelHyppsiprime/BBB_UnfoldingTuple"); 
    TNtuple* pretuple = (TNtuple*)preselfile->Get("Tuple"); 

    pretuple->SetBranchAddress(nameOfSensVarGen, &diffabsy_gen);
    

    float mgen, ptgen, ygen;
    pretuple->SetBranchAddress("mttbar_gen", &mgen);
    pretuple->SetBranchAddress("pt_ttbar_gen", &ptgen);
    pretuple->SetBranchAddress("y_ttbar_gen", &ygen);

    if(!TString(nameOfXAxisVarRec).CompareTo("mttbar_rec")) mttbargen = &mgen;
    if(!TString(nameOfXAxisVarRec).CompareTo("pt_ttbar")) mttbargen = &ptgen;
    if(!TString(nameOfXAxisVarRec).CompareTo("y_ttbar_rec")) mttbargen = &ygen;

    pretuple->SetBranchAddress("weight", &weight);
    weight *= lumi_factor;

    nentries = pretuple->GetEntries();

    nentries_toconsider = nentries; 
    nentries_toskip = 0;

    for(int i=nentries_toskip; i<nentries_toconsider; i++)
    {

      pretuple->GetEntry(i);
      weight *= sampleweight;

      
      const double valx = fabs(*mttbargen);
      const double valy = diffabsy_gen;
      
      const double reweighting_loweredge = 0;
      const double reweighting_upperedge = 75;
      
      //cout << valx << "\t" << valy << "\t" << weight << "\n";
      
      fill_nooverflow_2d(preNfine, valx, valy, weight);
      preselected_nonflat->Fill(valx, valy, weight);
      preselected_nonflat_fewbins->Fill(valx, valy, weight);
      
      kinematicpreselplot->Fill(valx, weight);
      
      preN->fill(valx, valy, weight);

      
      #ifdef BUILD_1D_BASIS
        fill_nooverflow_1d(preN1d, valy, weight); 
      #endif
      
      preNkm05->fill(valx, valy, weight*(1-0.05*diffabsy_gen));
      preNkm10->fill(valx, valy, weight*(1-0.10*diffabsy_gen));
      preNkm15->fill(valx, valy, weight*(1-0.15*diffabsy_gen));
      preNkm20->fill(valx, valy, weight*(1-0.20*diffabsy_gen));
      preNkm25->fill(valx, valy, weight*(1-0.25*diffabsy_gen));
      
      preNkp05->fill(valx, valy, weight*(1+0.05*diffabsy_gen));
      preNkp10->fill(valx, valy, weight*(1+0.10*diffabsy_gen));
      preNkp15->fill(valx, valy, weight*(1+0.15*diffabsy_gen));
      preNkp20->fill(valx, valy, weight*(1+0.20*diffabsy_gen));
      preNkp25->fill(valx, valy, weight*(1+0.25*diffabsy_gen));
      
      {
        double valx2 = fabs(mgen);
        const double reweighting_loweredge = 340;
        const double reweighting_upperedge = 650;
        
        double contk, contweight;
        
        // below: lots of reweighted histograms are filled for checks of the unfolding's reliability
        
//        const double wfactor = 0.125; // used to be 0.25
        const double wfactor = 0.375; // used to be 0.25
        
        // noninverted, squared
        contk = 2*wfactor * (valx2-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge)
                          * (valx2-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge);
        if(contk < 0) contk = 0;
        if(contk > wfactor) contk = wfactor;
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontsquare_m->fill(valx, valy, weight*contweight);
        
        // noninverted
        contk = wfactor * (valx2-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge);
        if(contk < 0) contk = 0;
        if(contk > wfactor) contk = wfactor;
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcont_m->fill(valx, valy, weight*contweight);
        
        //... and negative
        contk *= -1;
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontneg_m->fill(valx, valy, weight*contweight);
        contk *= -1;
        
        // inverted
        contk = wfactor - contk; // invert reweighting to be bigger for smaller M
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontinv_m->fill(valx, valy, weight*contweight);
        
        // ... and negative
        contk *= -1;
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontneginv_m->fill(valx, valy, weight*contweight);
        contk *= -1;
        
        // centered
        contk = wfactor * (valx2-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge);
        if(contk < 0) contk = 0;
        if(contk > wfactor) contk = wfactor;
        contk -= wfactor/2; // center
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontcent_m->fill(valx, valy, weight*contweight);
        
        // centered, inverted
        contk = -1 * contk; // invert
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontcentinv_m->fill(valx, valy, weight*contweight);      
      }
      
      // same as above, but reweighted in pt(ttbar)
      {
        double valx2 = fabs(ptgen);
        const double reweighting_loweredge = 0;
        const double reweighting_upperedge = 75;
        
        double contk, contweight;
        
        const double wfactor = 0.125; // used to be 0.25
        
        // noninverted, squared
        contk = 2*wfactor * (valx2-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge)
        * (valx2-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge);
        if(contk < 0) contk = 0;
        if(contk > wfactor) contk = wfactor;
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontsquare_pt->fill(valx, valy, weight*contweight);
        
        // noninverted
        contk = wfactor * (valx2-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge);
        if(contk < 0) contk = 0;
        if(contk > wfactor) contk = wfactor;
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcont_pt->fill(valx, valy, weight*contweight);
        
        //... and negative
        contk *= -1;
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontneg_pt->fill(valx, valy, weight*contweight);
        contk *= -1;
        
        // inverted
        contk = wfactor - contk; // invert reweighting to be bigger for smaller M
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontinv_pt->fill(valx, valy, weight*contweight);
        
        // ... and negative
        contk *= -1;
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontneginv_pt->fill(valx, valy, weight*contweight);
        contk *= -1;
        
        // centered
        contk = wfactor * (valx2-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge);
        if(contk < 0) contk = 0;
        if(contk > wfactor) contk = wfactor;
        contk -= wfactor/2; // center
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontcent_pt->fill(valx, valy, weight*contweight);
        
        // centered, inverted
        contk = -1 * contk; // invert
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontcentinv_pt->fill(valx, valy, weight*contweight);      
      }
      
      // same as above, but reweighted in y(ttbar)
      {
        double valx2 = fabs(ygen);
        const double reweighting_loweredge = 0;
        const double reweighting_upperedge = 1.1;
        
        double contk, contweight;
        
        const double wfactor = 0.125; // used to be 0.25
        
        // noninverted, squared
        contk = 2*wfactor * (valx2-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge)
        * (valx2-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge);
        if(contk < 0) contk = 0;
        if(contk > wfactor) contk = wfactor;
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontsquare_y->fill(valx, valy, weight*contweight);
        
        // noninverted
        contk = wfactor * (valx2-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge);
        if(contk < 0) contk = 0;
        if(contk > wfactor) contk = wfactor;
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcont_y->fill(valx, valy, weight*contweight);
        
        //... and negative
        contk *= -1;
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontneg_y->fill(valx, valy, weight*contweight);
        contk *= -1;
        
        // inverted
        contk = wfactor - contk; // invert reweighting to be bigger for smaller M
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontinv_y->fill(valx, valy, weight*contweight);
        
        // ... and negative
        contk *= -1;
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontneginv_y->fill(valx, valy, weight*contweight);
        contk *= -1;
        
        // centered
        contk = wfactor * (valx2-reweighting_loweredge)/(reweighting_upperedge-reweighting_loweredge);
        if(contk < 0) contk = 0;
        if(contk > wfactor) contk = wfactor;
        contk -= wfactor/2; // center
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontcent_y->fill(valx, valy, weight*contweight);
        
        // centered, inverted
        contk = -1 * contk; // invert
        contweight = (1+contk*diffabsy_gen);
        if(contweight < 0) contweight = 0;
        if(contweight > 2) contweight = 2;
        preNkcontcentinv_y->fill(valx, valy, weight*contweight);      
      }
      
      
      const int genbinm = preN->findXBin(valx);
      senspreselplots[genbinm]->Fill(valy, weight);
      senspreselplotinc->Fill(valy, weight);
      preselcorrelationplot->Fill(valx, valy, weight);
    }
  }
  

  
  #ifdef POWHEG
    int nPostSamples = 1;
  #endif
  
  /// fill post histo
  for(int j=0; j<nPostSamples; j++)
  {
    cout << "Starting with postsel sample #" << j << ".\n";
    TFile* postselfile = 0;
    H2gen* postN;
    //double sampleweight = 1;
    
    #ifdef POWHEG
      if(j==0)
      {
        postN = postNele;
        postselfile = new TFile(sl_ttbarele); 
      }
      /*      if(j==1)
      {
        postN = postNmu;
        postselfile = new TFile(sl_ttbarmu); 
      }
      */
    #endif
    
//    TNtuple* posttuple = (TNtuple*)postselfile->Get("TTbarSelHyppsiprime/BBB_UnfoldingTuple"); 
    TNtuple* posttuple = (TNtuple*)postselfile->Get("Tuple"); 
    float mttbargen;
    posttuple->SetBranchAddress(nameOfXAxisVarGen, &mttbargen);
    posttuple->SetBranchAddress(nameOfSensVarGen, &diffabsy_gen);
    posttuple->SetBranchAddress("weight", &weight);
    weight *= lumi_factor;

    nentries = posttuple->GetEntries();

    nentries_toconsider = nentries; 
    nentries_toskip = 0;
    

    
    //for(int i=nentries_toconsider; i<nentries; i++)
    for(int i=nentries_toskip; i<nentries_toconsider; i++)
    {

          
      posttuple->GetEntry(i);
      //weight *= sampleweight;
      
      postN->fill(fabs(mttbargen), diffabsy_gen, weight);
      postNges->fill(fabs(mttbargen), diffabsy_gen, weight);
      #ifdef BUILD_1D_BASIS
        fill_nooverflow_1d(postN1d, diffabsy_gen, weight); 
      #endif     
        
      
      // total histo for visualization, so non-symmetric
      fill_nooverflow_2d(postNfine, fabs(mttbargen), diffabsy_gen, weight);
      
      
      const int genbinm = preN->findXBin(fabs(mttbargen));
      sensgenselplots[genbinm]->Fill(diffabsy_gen, weight);
      genselcorrelationplot->Fill(fabs(mttbargen), diffabsy_gen, weight);
      
    }
    
  
  }

  /// finish
  
  nonselected = preN->clone();
  nonselected->setName("nonselected_reweighted");

  nonselected->add(postNele, -betattbarele);
  //  nonselected->add(postNmu, -betattbarmu); 
  
  nonselectedtotal = preN->clone();
  nonselectedtotal->setName("nonselected_total");
  nonselectedtotal->add(postNele,-1);
  //  nonselectedtotal->add(postNmu, -1); 
  
  nonselectedele = preN->clone();
  nonselectedele->setName("nonselectedele");
  nonselectedele->add(postNele,-1);
  /*
  nonselectedmu = preN->clone();
  nonselectedmu->setName("nonselectedmu");
  nonselectedmu->add(postNmu,-1);
  */
  seleff = postNges->clone();
  seleff->setName("seleff");
  seleff->divide(preN);
  
  seleffele = postNele->clone();
  seleffele->setName("seleffele");
  seleffele->divide(preN);
  /*
  seleffmu = postNmu->clone();
  seleffmu->setName("seleffmu");
  seleffmu->divide(preN);
  */
  selefffine = (TH2F*)postNfine->Clone();
  selefffine->SetName("seleff_fine");
  selefffine->Divide(preNfine);
  

  /// output
  
  TH2F* nonselected_t = nonselected->toVisTH2F();
  cmsstyle->setup_style_2D(nonselected_t, labelOfXAxisVar, labelOfSensVar);
  nonselected_t->SetMinimum(0);
  nonselected_t->Draw("COLZ");
  drawGenBinEdges();
  saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring+"nonselected"));
  
  TH2F* preselected_t = preN->toVisTH2F();
  cmsstyle->setup_style_2D(preselected_t, labelOfXAxisVar, labelOfSensVar);
  preselected_t->Draw("COLZ");
  drawGenBinEdges();
  saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring+"preselection_automin"));
  
  preselected_t->SetMinimum(0);
  preselected_t->Draw("COLZ");
  //preselected_t->GetZaxis()->SetRangeUser(2460.,2500.);
  drawGenBinEdges();
  saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring+"preselection"));
  
  cmsstyle->setup_style_2D(preselected_nonflat, labelOfXAxisVar, labelOfSensVar);
  preselected_nonflat->SetMinimum(0);
  preselected_nonflat->Draw("COLZ");
  saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring+"preselection_nonflat_nolines"));
  drawGenBinEdges();
  saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring+"preselection_nonflat"));
  
  {
    for(int i=0; i<preselected_nonflat_fewbins->GetNbinsX(); i++)
    {
      double poserr, negerr; // new from Frank
//      double pos = preselected_nonflat_fewbins->Integral(i+1, i+1, 6, 10);
//      double neg = preselected_nonflat_fewbins->Integral(i+1, i+1, 1, 5);
      double pos = preselected_nonflat_fewbins->IntegralAndError(i+1, i+1, 6, 10, poserr); // new from Frank
      double neg = preselected_nonflat_fewbins->IntegralAndError(i+1, i+1, 1, 5, negerr); // new from Frank
      if(pos+neg == 0)
      {
        pos = neg = 1;
      }
      const double asy = (pos-neg)/(pos+neg);
//      const double err = asymmetry_error_naive(pos, neg);
      const double err = asymmetry_error_custom_errors(pos, poserr, neg, negerr); // new from Frank
      preselected_nonflat_asys->SetBinContent(i+1, asy);
      preselected_nonflat_asys->SetBinError(i+1, err);
    }
    cmsstyle->setup_style(preselected_nonflat_asys, labelOfXAxisVar, "A_{C}", kGreen+2, 0, 0);  
    double marg = canv1->GetRightMargin();
    canv1->SetRightMargin(0.04);
    preselected_nonflat_asys->SetMaximum(preselected_nonflat_asys->GetMaximum()*1.4);
    preselected_nonflat_asys->Draw("EHIST");
    TLegend leg(0.80, 0.85, 0.92, 0.92);
    leg.AddEntry(preselected_nonflat_asys,"t#bar{t}","el")->SetTextAlign(13);
    leg.SetBorderSize(0);
    leg.SetLineStyle(0);
    leg.SetTextFont(42);
    leg.SetFillStyle(0);
    leg.Draw();
    cmsstyle->CMSSimulation();  
    saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring+"preselection_asys_"+nameOfXAxisVarShort));
    canv1->SetRightMargin(marg);
    
  }
    
  preselected_nonflat->Draw("COLZ");
  drawGenBinEdges(2);
  saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring+"preselection_nonflat_thicklines"));
  
  TH2F* postNges_t = postNges->toVisTH2F();
  cmsstyle->setup_style_2D(postNges_t, labelOfXAxisVar, labelOfSensVar);
  postNges_t->SetMinimum(0);
  postNges_t->Draw("COLZ");
  drawGenBinEdges();
  saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring+"selected"));
  
  TH2F* seleff_t = seleff->toVisTH2F();
  cmsstyle->setup_style_2D(seleff_t, labelOfXAxisVar, labelOfSensVar);
  seleff_t->SetZTitle("selection efficiency");
  Float_t prevmargin = canv1->GetRightMargin();
  canv1->SetRightMargin(0.19);
  seleff_t->SetMinimum(0);
  seleff_t->Draw("COLZ");
  drawGenBinEdges();
  saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring+"seleff_"+nameOfXAxisVarShort));
  canv1->SetRightMargin(prevmargin);
  
  cmsstyle->setup_style_2D(selefffine, labelOfXAxisVar, labelOfSensVar);
  selefffine->SetMinimum(0);
  selefffine->SetMaximum(0.25);
  selefffine->Draw("COLZ");
  drawGenBinEdges();
  saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring+"selefffine"));
  
  #ifdef BUILD_1D_BASIS
    TH1F* seleff1d = (TH1F*) postN1d->Clone();
    seleff1d->Divide(preN1d);
    
    cmsstyle->setup_style(seleff1d, labelOfSensVar, "selection efficiency", 1, 0, 0);
    seleff1d->SetMinimum(0);
    seleff1d->SetMaximum(seleff1d->GetMaximum()*1.5);
    seleff1d->Draw("HIST");
    saveAs(canv1, (TString("output")+fileSuffix+"/"+dimstring+"seleff_1dplot"));
  #endif
  
  cmsstyle->setup_style(senspreselplotinc, labelOfSensVar, "N", 1, 0, 0);
  senspreselplotinc->Draw("HIST");
  saveAs(gPad, (TString("output")+fileSuffix+"/"+ dimstring + "senspreselplotinc" ));
  senspreselplotinc->SetMaximum(senspreselplotinc->GetMaximum()*1.2);
  senspreselplotinc->SetMinimum(0);
  draw1dGenBinEdges(senspreselplotinc, 2);
  saveAs(gPad, (TString("output")+fileSuffix+"/"+ dimstring + "senspreselplotinc_withlines" ));
  
  {
    double marg = canv1->GetRightMargin();
    canv1->SetRightMargin(0.04);
    cmsstyle->setup_style(kinematicpreselplot, labelOfXAxisVar, "N", 1, 0, 0);
    kinematicpreselplot->SetMaximum(kinematicpreselplot->GetMaximum()*1.2);  
    kinematicpreselplot->Draw("HIST");
    TLegend leg(0.80, 0.85, 0.92, 0.92);
    leg.AddEntry(kinematicpreselplot,"t#bar{t}","l");
    leg.SetBorderSize(0);
    leg.SetLineStyle(0);
    leg.SetTextFont(42);
    leg.SetFillStyle(0);
    leg.Draw();
    cmsstyle->CMSSimulation();
    saveAs(gPad, (TString("output")+fileSuffix+"/"+ dimstring + "kinematicpreselplot_"+nameOfXAxisVarShort));
    canv1->SetRightMargin(marg);
  }
    
  for(int i=0; i<nbinsmafter; i++)
  {
    cmsstyle->setup_style(senspreselplots[i], labelOfSensVar, "N", 1, 0, 0);
    senspreselplots[i]->Draw("E");
    saveAs(gPad, (TString("output")+fileSuffix+"/"+ dimstring + TString::Format("senspreselplotbin%d",i) ));
  }
  for(int i=0; i<nbinsmafter; i++)
  {
    cmsstyle->setup_style(sensgenselplots[i], labelOfSensVar, "N", 1, 0, 0);
    sensgenselplots[i]->Draw("E");
    saveAs(gPad, (TString("output")+fileSuffix+"/"+ dimstring + TString::Format("sensgenselplotbin%d",i) ));
  }
  
  
  cmsstyle->setup_style_2D(preselcorrelationplot, labelOfXAxisVar, labelOfSensVar);
  preselcorrelationplot->Draw("COLZ");
  drawGenBinEdges();
  saveAs(gPad, (TString("output")+fileSuffix+"/"+dimstring  + "preselcorrelationplot" ));
  
  cmsstyle->setup_style_2D(genselcorrelationplot, labelOfXAxisVar, labelOfSensVar);
  genselcorrelationplot->Draw("COLZ");
  drawGenBinEdges();
  saveAs(gPad, (TString("output")+fileSuffix+"/"+dimstring + "genselcorrelationplot" ));
  
////////////  TFile* file = new TFile(TString("output/")+fileSuffix+"/"+dimstring+"seleff.root","recreate");
  TFile* file = new TFile(TString("output/")+dimstring+"seleff.root","recreate");
  file->cd();
  seleff->write();
  seleffele->write();
  //  seleffmu->write();
  selefffine->Write();
  nonselected->write();
  nonselectedtotal->write();
  nonselectedele->write();
  //  nonselectedmu->write();
  preN->write();
  postNges->write();
  postNele->write();
  //  postNmu->write();
  
  preNkm05->write();
  preNkm10->write();
  preNkm15->write();
  preNkm20->write();
  preNkm25->write();
  
  preNkp05->write();
  preNkp10->write();
  preNkp15->write();
  preNkp20->write();
  preNkp25->write();
  
  preNkcont_m->write();
  preNkcontinv_m->write();
  preNkcontneg_m->write();
  preNkcontneginv_m->write();
  preNkcontcent_m->write();
  preNkcontcentinv_m->write();
  preNkcontsquare_m->write();
  
  preNkcont_pt->write();
  preNkcontinv_pt->write();
  preNkcontneg_pt->write();
  preNkcontneginv_pt->write();
  preNkcontcent_pt->write();
  preNkcontcentinv_pt->write();
  preNkcontsquare_pt->write();
  
  preNkcont_y->write();
  preNkcontinv_y->write();
  preNkcontneg_y->write();
  preNkcontneginv_y->write();
  preNkcontcent_y->write();
  preNkcontcentinv_y->write();
  preNkcontsquare_y->write();
  
  preselected_nonflat->Write();

  kinematicpreselplot->Write();
  
  file->Close();
  
  return 0;
}


int main()
{
  return seleffmatrix();
}
