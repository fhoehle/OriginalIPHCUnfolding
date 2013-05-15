#ifndef SPECIALHELPERS_H
#define SPECIALHELPERS_H

#include <iostream>
#include <list>
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
#include <TList.h>
#include <TPaveText.h>
#include <TVirtualPad.h>
#include <TClass.h>
#include "TApplication.h"
#include "TRandom3.h"

#include "binning.hpp"

using namespace std;


/// this file contains helper functions which depend on binning.hpp or other special settings




void unwrap2dhisto(H2gen* h2, TH1* h1)
{
  int n=0;
  for(int x=0; x<nbinsmafter; x++)
  {
    for(int y=0; y<nbinsetaafter; y++)
    {
      n++;
      h1->SetBinContent(n, h2->getBinContent(x,y));
      h1->SetBinError(n, h2->getBinError(x,y));
    }
  }
  //h1->Sumw2();
}

void unwrap2dhisto(H2rec* h2, TH1* h1)
{
  int n=0;
  for(int x=0; x<nbinsmafter; x++)
  {
    for(int y=0; y<nbinsetaafter; y++)
    {
      n++;
      h1->SetBinContent(n, h2->getBinContent(x,y));
      h1->SetBinError(n, h2->getBinError(x,y));
    }
  }
  //h1->Sumw2();
}


void drawGenBinEdges(int width = 0)
{
  for(int i=1; i<nbinsmafter; i++)
  {
    TLine line(xedgesgen[i], yedgesgen[i][0],xedgesgen[i],yedgesgen[i][nbinsetaafter]);
    if(width)
      line.SetLineWidth(width);
    line.DrawClone("SAME");
  }
  for(int x=0; x<nbinsmafter; x++)
  {
    for(int y=1; y<nbinsetaafter; y++)
    { 
      TLine line(xedgesgen[x], yedgesgen[x][y],xedgesgen[x+1],yedgesgen[x][y]);
      if(width)
        line.SetLineWidth(width);
      line.DrawClone("SAME");
    }
  } 
}

void draw1dGenBinEdges(TH1* h, int width = 0)
{
  for(int i=1; i<nbinsetaafter1d; i++)
  {
    TLine line(etabinedges1dafter[i], h->GetMinimum(), etabinedges1dafter[i], h->GetMaximum());
    if(width)
      line.SetLineWidth(width);
    line.DrawClone("SAME");
  }

}

void draw1dRecBinEdges(TH1* h, int width = 0)
{
  for(int i=1; i<nbinseta1d; i++)
  {
    TLine line(etabinedges1d[i], h->GetMinimum(), etabinedges1d[i], h->GetMaximum());
    if(width)
      line.SetLineWidth(width);
    line.DrawClone("SAME");
  }
  
}

void drawRecBinEdges()
{
  for(int i=1; i<nbinsm; i++)
  {
    TLine line(xedgesrec[i], yedgesrec[i][0],xedgesrec[i],yedgesrec[i][nbinseta]);
    line.DrawClone("SAME");
  }
  for(int x=0; x<nbinsm; x++)
  {
    for(int y=1; y<nbinseta; y++)
    { 
      TLine line(xedgesrec[x], yedgesrec[x][y],xedgesrec[x+1],yedgesrec[x][y]);
      line.DrawClone("SAME");
    }
  } 
}


void drawMigMatBinEdges()
{
  for(int i=1; i<nbinsmafter; i++)
  {
    TLine line(0.5+nbinsetaafter*i, 0.5, 0.5+nbinsetaafter*i, 0.5+nbins);
    line.DrawClone("SAME");
  }
  for(int i=1; i<nbinsm; i++)
  { 
    TLine line(0.5, 0.5+nbinseta*i, 0.5+nbinsafter, 0.5+nbinseta*i);
    line.DrawClone("SAME");
  }
   
}

#endif // SPECIALHELPERS_H
