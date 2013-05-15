///////////////////////
///// IMPORTANT
/////////////////////
// To find binning for "gen" use unselected/generated tuple
// To find binning for "rec" use selected/reconstructed tuple
//////////////////

#include "helpers.hpp"


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
#include <TPave.h>
#include <TPaveText.h>
#include <TVirtualPad.h>
#include <TClass.h>
#include "TApplication.h"
#include "TRandom3.h"
#include "TChain.h"

#include "binning.hpp"


/// we try to determine a binning for a variable that is flat in the reconstructed ttbar spectrum or something similar


/// finds the right binning for the named variables in the named tuple, with the named number of bins.
/// the output will be array definitions, named (x|y)edges[nameExt]
void getBinning(TNtuple* tuple, const char* nameExt, const char* xVarName, const int nXBins, const char* yVarName, const int nYBins, bool forceNoOverlaps)
{
  TH1::AddDirectory(kFALSE);
  
  const int nentries = tuple->GetEntries();
  
  float varx, vary, weight;
  tuple->SetBranchAddress(xVarName, &varx);
  tuple->SetBranchAddress(yVarName, &vary);
  tuple->SetBranchAddress("weight", &weight);
  weight *= lumi_factor;
  
  /// get min and max values
  double xmax = 0;
  double xmin = 2000;
  double ymax = 0;
  double ymin = 2000;
  for(int j=0; j<nentries; j++)
  {
    tuple->GetEntry(j);
    
    const float absvarx = fabs(varx);
      
    if(absvarx < xmin)
      xmin = absvarx;
    if(absvarx > xmax)
      xmax = absvarx;    
    
    if(vary < ymin)
      ymin = vary;
    if(vary > ymax)
      ymax = vary;

  }

  /// fill 1d histo with x values, lots of bins
  TH1F* histo = new TH1F("histo","histo", 10000, xmin, xmax);
  for(int j=0; j<nentries; j++)
  {
    tuple->GetEntry(j);
    
    const float var = fabs(varx);
        
    fill_nooverflow_1d(histo, var, weight);
  }
  
  /// split histogram into properly-sized bits, output as array code
  const double totint = histo->Integral();
  const double desiredBinIntegral = totint/nXBins;
  
  cout << "const double xedges"<<nameExt<<"[] = {" << xmin << ", ";
  
  double xedges[nXBins+1];
  xedges[0] = xmin;
  int iedge = 1;
  
  int prevBinRangeEnd = 0;
  for(int i=1; i<=10000; i++)
  {
    const double integral = histo->Integral(prevBinRangeEnd+1, i);
    if(integral>desiredBinIntegral)
    {
      const double edge = histo->GetXaxis()->GetBinUpEdge(i);
      
      cout << edge << ", ";
      xedges[iedge++] = edge;
      
      prevBinRangeEnd = i;
    }
  }
  
  const double lastedge = histo->GetXaxis()->GetBinUpEdge(10000);
  cout << lastedge << "};" << endl;
  xedges[nXBins] = lastedge; 
  
  /// now for the y binnings
  // fill new histogram, with known x binning and lots of bins for y
  TH2F* histo2 = new TH2F("histo2","histo2", nXBins, xedges , 10000, ymin, ymax);
  for(int j=0; j<nentries; j++)
  {
    tuple->GetEntry(j);
    
    const float absvarx = fabs(varx);
    
    fill_nooverflow_2d(histo2, absvarx, vary, weight);
  }
  
  // find the right y binning for each bin of x
  //cout << "const double yedges"<<nameExt<<"[]["<< nYBins+1<<"] = { " ;

  double yedges[nXBins][nYBins+1];
  
  for(int x=0; x<nXBins; x++)
  {
    int xlowbin = x+1;
    int xhighbin = x+1;
    
    if(forceNoOverlaps)
    {
      xlowbin = 1;
      xhighbin = nXBins;
    }
    
    const int ntestbins = 10000;
    const double totint = histo2->Integral(xlowbin, xhighbin, 1, ntestbins);
    const double desiredBinIntegral = totint/nYBins;
    
    //if(x!=0)
    //  cout << ",\n\t\t";
    //cout << " {" << ymin << ", ";
    
    yedges[x][0] = ymin;
    int iedge = 1;
    
    int prevBinRangeEnd = 0;
    for(int i=1; i<=ntestbins; i++)
    {
      const double integral = histo2->Integral(xlowbin, xhighbin, prevBinRangeEnd+1, i);
      if(integral>desiredBinIntegral)
      {
        double newedge = histo2->GetYaxis()->GetBinUpEdge(i);
        
        if(fabs(newedge)<0.01) 
          newedge = 0; // so that we don't make errors calculating asys later
        
        //cout << newedge << ", ";
        yedges[x][iedge++] = newedge;
        
        prevBinRangeEnd = i;
      }
    }
    
    const double lastedge = histo2->GetYaxis()->GetBinUpEdge(ntestbins);
    //cout << lastedge << "}";
    yedges[x][nYBins] = lastedge; 
    
  }
  //cout << "};" << endl;
  
  // normal output is finished, now do symmetrized binning
  
  
  cout << "const double yedges"<<nameExt<<"[]["<< nYBins+1<<"] = { " ;
   
  for(int x=0; x<nXBins; x++)
  {
    if(x!=0)
      cout << ",\n\t\t";
    
    
    double yedges_sym[nYBins+1];
    yedges_sym[nYBins/2] = 0;
    for(int i=0; i<nYBins/2; i++)
    {
      yedges_sym[i] =        -( fabs(yedges[x][i]) + fabs(yedges[x][nYBins-i]) ) /2.;
      yedges_sym[nYBins-i] = +( fabs(yedges[x][i]) + fabs(yedges[x][nYBins-i]) ) /2.;
    }
    
    cout << " {" << yedges_sym[0] << ", ";
    for(int i=1; i<nYBins; i++)
    {
      cout << yedges_sym[i] << ", ";
    }
    cout << yedges_sym[nYBins] << "}";
    
  }
  cout << "};" << endl;

}

void findbinning()
{
  TFile* file_ttbar;
  TNtuple* tuple;
  
  /// rec
  #ifdef MADGRAPH
////////    file_ttbar = new TFile("/portal/ekpams2/home/froscher/outputroot/2dunfolding/ttbar_histos.root");
//////    file_ttbar = new TFile("ttbar_mu_histos.root");




//    file_ttbar = new TFile("/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileNonSkimmed_TTbarAllMCatNLO.root");
    file_ttbar = new TFile("/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_TTbarSignalTuneZ2_7TeV-madgraph-tauola_PileUpReweighted.root");
//      file_ttbar = new TFile("/afs/cern.ch/work/c/cardaci/Tuples_and_Results_24_May/TupleFileSkimmed_TTbarAllMCatNLO_PileUpReweighted.root");



//    file_ttbar = new TFile("/afs/cern.ch/work/c/cardaci/precut_merged_nlo.root");
//    tuple = (TNtuple*)file_ttbar->Get("TTbarSelHyppsiprime/BBB_UnfoldingTuple");
    tuple = (TNtuple*)file_ttbar->Get("Tuple");
  #endif
  #ifdef POWHEG
//    TChain ch("TTbarSelHyppsiprime/BBB_UnfoldingTuple");

    TChain ch("Tuple");
//    TChain ch("DataTuple");

    // WARNING we should actually reweight here
///////////    ch.Add("/portal/ekpams2/home/froscher/outputroot/powhegsample/ttbar_ele_combined_histos.root");
//////////    ch.Add("/portal/ekpams2/home/froscher/outputroot/powhegsample/ttbar_mu_combined_histos.root");
///////    ch.Add("ttbar_mu_histos.root");
////////    ch.Add("ttbar_mu_histos.root");




//    ch.Add("/afs/cern.ch/work/c/cardaci/Tuples_and_Results_29_June_1btag_selection/TupleFileNonSkimmed_TTbarAllMCatNLO.root");
    ch.Add("/afs/cern.ch/work/c/cardaci/Tuples_and_Results_29_June_1btag_selection/TupleFileSkimmed_TTbarSignalTuneZ2_7TeV-madgraph-tauola_PileUpReweighted.root");
//    ch.Add("/afs/cern.ch/work/c/cardaci/Tuples_and_Results_29_June_1btag_selection/TupleFileSkimmed_TTbarAllMCatNLO_PileUpReweighted.root");




//    ch.Add("/afs/cern.ch/work/c/cardaci/precut_merged_nlo.root");

  //  ch.Add("TupleFileSkimmed.root");
    tuple = (TNtuple*)&ch;
  #endif
  
  
  //getBinning(tuple, "rec", nameOfXAxisVarRec, 6, nameOfSensVarRec, 16, false); //2d
  getBinning(tuple, "rec", nameOfXAxisVarRec, 1, nameOfSensVarRec, 16, true); //1d
  //getBinning(tuple, "gen", nameOfXAxisVarGen, 3, nameOfSensVarGen, 8, false); //2d
  //getBinning(tuple, "gen", nameOfXAxisVarGen, 1, nameOfSensVarGen, 8, true); //1d
  
  /// presel-gen && unfolded  
/*  #ifdef MADGRAPH
    file_ttbar = new TFile("/storage/8/froscher/outputroot/2dunfolding/presel/ttbar_3_pre_ttbar.root");
    tuple = (TNtuple*)file_ttbar->Get("DataTuple"); 
  #endif
  #ifdef POWHEG
    TChain ch2("DataTuple");
    // WARNING we should actually reweight here
    ch2.Add("/storage/8/froscher/outputroot/powheg_genasy/dilep_0_pre_ttbar.root");
    ch2.Add("/storage/8/froscher/outputroot/powheg_genasy/semilep_1_pre_ttbar.root");
    ch2.Add("/storage/8/froscher/outputroot/powheg_genasy/semilep_2_pre_ttbar.root");
    ch2.Add("/storage/8/froscher/outputroot/powheg_genasy/semilep_3_pre_ttbar.root");
    ch2.Add("/storage/8/froscher/outputroot/powheg_genasy/semilep_4_pre_ttbar.root");
    ch2.Add("/storage/8/froscher/outputroot/powheg_genasy/semilep_5_pre_ttbar.root");
    tuple = (TNtuple*)&ch2;
  #endif
  
  getBinning(tuple, "gen", nameOfXAxisVarGen, 3, nameOfSensVarGen, 8, false); //2d
  //getBinning(tuple, "gen", nameOfXAxisVarGen, 1, nameOfSensVarGen, 8, true); //1d
*/


  return; 
}

int main()
{
  findbinning();
  return 0;
}

