#include "cmsstyle.hpp"
#include "fitresults.hpp"
#include "helpers.hpp"
#include "samplelocations.hpp"


#include "binning.hpp"

// WARNING make sure variables in binning.hpp are set correctly


//#include "myunfold_class.hpp"
#include "myunfold_class1d.hpp"


int classy_unfold1d()
{
  TH1::SetDefaultSumw2(true);

  TFile* inputfile;
  
  
  inputfile = new TFile(sl_data_comb);

  
 
//  TFile* matrixfile = new TFile(TString("output/")+fileSuffix+"1dmigmatrix.root");
//  TFile* preselfile = new TFile(TString("output/")+fileSuffix+"1dseleff.root");
  TFile* matrixfile = new TFile(TString("output/")+"1dmigmatrix.root");
  TFile* preselfile = new TFile(TString("output/")+"1dseleff.root");
  MyUnfold1d myunfold;

  CMSStyle* cmsstyle = new CMSStyle();
  gStyle->SetOptStat("");
  TCanvas* canv1 = new TCanvas("canv1","canv1",300,0,800,600);
  canv1->SetRightMargin(0.04);
  // to make z-axis stuff show completely in 2d plots
  canv1->SetRightMargin(0.14);

  
  /// get migration matrix and presel distribution
  TH2F* migmatrix = (TH2F*)matrixfile->Get("migmatrixcombined");
  cmsstyle->setup_style_2D(migmatrix, "gen", "rec");
  
  TH2F* hpresel = (TH2F*)preselfile->Get("preN");
  
  
  /// create 2d histos
  TH2F* hrec;
  TH2F* hgen;
  
  H2rec1d myh2rec("hrec", xedgesrec1d, yedgesrec1d);
  H2gen1d myh2gen("hgen", xedgesgen1d, yedgesgen1d);
  
  TH1F* sensrecoplots[nbinsmafter1d];
  for(int i=0; i<nbinsmafter1d; i++)
  {
    sensrecoplots[i] = new TH1F("","", 100, yedgesgen1d[0][0], yedgesgen1d[0][nbinsetaafter1d]);
  }
  
  
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
  int nentries = tuple->GetEntries();
  
  // so that we're stat. independent from migmatrix
  int nentries_toskip = 0; 
                               
  int nentries_toconsider = nentries;

  for(int i=nentries_toskip; i<nentries_toconsider; i++)
  {
    
    tuple->GetEntry(i);
    
    double weightasymmetry = 1;
    
    //if(reweight) weightasymmetry = reweightk * diffabsetagen + 1;
    
    //fill_nooverflow_2d(hgen, fabs(mttbargen), diffabsetagen, weight*weightasymmetry);
    //fill_nooverflow_2d(hrec, fabs(mttbarrec), diffabsetarec, weight*weightasymmetry);
    myh2gen.fill(fabs(mttbargen), diffabsetagen, weight*weightasymmetry);
    myh2rec.fill(fabs(mttbarrec), diffabsetarec, weight*weightasymmetry);
    
    const int genbinm = myh2gen.findXBin(fabs(mttbarrec));
    sensrecoplots[genbinm]->Fill(diffabsetarec, weight*weightasymmetry);
  }
  
  hgen = myh2gen.toTH2F();
  hrec = myh2rec.toTH2F();
  

  /// subtract background manually (not relevant to unfolding, just for some additional numbers)
  TH2F* hrecwithoutbg = (TH2F*)hrec->Clone();
  
  if(workingondata)
  {
    cerr << "Warning: The following warnings can be ignored. " << endl;
    const int size = myunfold.bgwrappedlist.size();
    for(int i=0; i<size; i++)
      hrecwithoutbg->Add(myunfold.bgwrappedlist[i], -1);    
    
    cout << "Integral before bg subtraction: " << hrec->Integral() << endl;
    cout << "Integral after bg subtraction: " << hrecwithoutbg->Integral() << endl;
  }
  

  
  TH2F* hunfold = myunfold.unfoldHisto(hrec, workingondata); // where the magic happens...
  
   
  TString style = "COLZ";
  
  hrec->SetMinimum(0);
  hrec->SetMaximum(hrec->GetMaximum()*1.1);
  hrecwithoutbg->SetMinimum(0);
  hrecwithoutbg->SetMaximum(hrec->GetMaximum()*1.1);
  hunfold->SetMinimum(0);
  hunfold->SetMaximum(hunfold->GetMaximum()*1.1);
  
  cmsstyle->setup_style_2D(hunfold, labelOfXAxisVar, labelOfSensVar);
  cmsstyle->setup_style_2D(hrec, labelOfXAxisVar, labelOfSensVar);
  cmsstyle->setup_style_2D(hpresel, labelOfXAxisVar, labelOfSensVar);
  hunfold->SetNdivisions(505);
  hrec->SetNdivisions(505);
  hpresel->SetNdivisions(505);
  

  

  hrecwithoutbg->Draw(style);
  saveAs(gPad, (TString("output/")+fileSuffix+"1dclassyunfold_hrecwithoutbg"));

  hunfold->Draw(style);
  saveAs(gPad, (TString("output/")+fileSuffix+"1dclassyunfold_hunfold"));
  hrec->Draw(style);
  saveAs(gPad, (TString("output/")+fileSuffix+"1dclassyunfold_hrec"));
  if(myunfold.lcurve)
  {
    myunfold.lcurve->Draw("AP");
    saveAs(gPad, (TString("output/")+fileSuffix+"1dclassyunfold_lcurve"));
  }
  myunfold.correlationmatrix->Draw("COLZ");
  saveAs(gPad, (TString("output/")+fileSuffix+"1dclassyunfold_correlationmatrix"));

  myunfold.errormatrix->Draw("COLZ");
  saveAs(gPad, (TString("output/")+fileSuffix+"1dclassyunfold_errormatrix"));
  
  for(int i=0; i<nbinsmafter1d; i++)
  {
    cmsstyle->setup_style(sensrecoplots[i], labelOfSensVar, "N", 1, 0, 0);
    sensrecoplots[i]->Draw("E");
    saveAs(gPad, (TString("output/")+fileSuffix+ TString::Format("1dclassyunfold_sensrecoplotbin%d",i) ));
  }
  
  /// calc asy
  double pos, neg, asy;
  
  cout << endl << endl;
   

  
  pos = hrec->Integral(1, nbinsm1d, nbinseta1d/2+1, nbinseta1d);
  neg = hrec->Integral(1, nbinsm1d, 1, nbinseta1d/2);
  asy = (pos-neg)/(pos+neg);
  cout << "reconstructed raw asy: " << asy << endl;
  cout << "err for that: " << asymmetry_error_naive(pos, neg) << endl;
  
  { 
    double poserr, negerr;
    pos = hrecwithoutbg->IntegralAndError(1, nbinsm1d, nbinseta1d/2+1, nbinseta1d, poserr);
    neg = hrecwithoutbg->IntegralAndError(1, nbinsm1d, 1, nbinseta1d/2, negerr);
    asy = (pos-neg)/(pos+neg);
    cout << "bg-subtracted asy: " << asy << endl;
    cout << "err for that: " << asymmetry_error_custom_errors(pos, poserr, neg, negerr) << endl << endl;
  }
  
  
  pos = hunfold->Integral(1, nbinsmafter1d, nbinsetaafter1d/2+1, nbinsetaafter1d);
  neg = hunfold->Integral(1, nbinsmafter1d, 1, nbinsetaafter1d/2);
  asy = (pos-neg)/(pos+neg);
  double asy_error = asymmetryerror_afterunfolding_1d(myunfold.errormatrix, myunfold.hunfoldunwrapped);
  {
    cout << "shift on the unfolded asy: " << asy - 0.0167293 << endl;
    cout << "unfolded asy: " << asy << endl;
    cout << "err for that: " << asy_error << endl;

  }
  
  
  
  {
    canv1->SetRightMargin(0.04);
    TH1F* hunfold1d = new TH1F("hunfold1d", "hunfold1d", nbinsetaafter1d, etabinedges1dafter);
    TH1F* hsim1d = new TH1F("hsim1d", "hsim1d", nbinsetaafter1d, etabinedges1dafter);
    
    TH1F* hpreselunwrapped = (TH1F*) myunfold.hunfoldunwrapped->Clone();
    unwrap2dhisto(hpresel, hpreselunwrapped);
    
    cmsstyle->setup_style(hunfold1d, labelOfSensVar, "1/#sigma d#sigma/d(|y_{t}|-|y_{#bar{t}}|)", 1, 0, 0);
    cmsstyle->setup_style(hsim1d, labelOfSensVar, "1/#sigma d#sigma/d(|y_{t}|-|y_{#bar{t}}|)", 2, 0, 0);
    
    
    for(int i=0; i< nbinsetaafter1d; i++)
    {
      hunfold1d->SetBinContent(i+1, myunfold.hunfoldunwrapped->GetBinContent(i+1));
      hunfold1d->SetBinError(i+1, myunfold.hunfoldunwrapped->GetBinError(i+1));
      
      hsim1d->SetBinContent(i+1, hpreselunwrapped->GetBinContent(i+1));
      hsim1d->SetBinError(i+1, hpreselunwrapped->GetBinError(i+1));
    }
    
    scale_to(hunfold1d, 1.0);
    scale_to(hsim1d, 1.0);
    
    for(int i=0; i< nbinsetaafter1d; i++)
    {
      hunfold1d->SetBinContent(i+1, hunfold1d->GetBinContent(i+1)/hunfold1d->GetBinWidth(i+1));
      hunfold1d->SetBinError(i+1, hunfold1d->GetBinError(i+1)/hunfold1d->GetBinWidth(i+1));
      
      hsim1d->SetBinContent(i+1, hsim1d->GetBinContent(i+1)/hsim1d->GetBinWidth(i+1));
      hsim1d->SetBinError(i+1, hsim1d->GetBinError(i+1)/hsim1d->GetBinWidth(i+1));
    }
    
    hunfold1d->SetMaximum();
    hunfold1d->SetMinimum(0);
    const double offset = hunfold1d->GetMaximum() - hunfold1d->GetMinimum();
    hunfold1d->SetMaximum(hunfold1d->GetMaximum()+offset*1.1);
    
    char resulttext[40] = "";
    sprintf(resulttext, "A_{C}=%.3f #pm %.3f", asy, asy_error);
    

    hunfold1d->SetMarkerStyle(23);
    hunfold1d->SetMarkerSize(1.5);
    hunfold1d->Draw("Ehist");

    hsim1d->SetFillStyle(3353);
    hsim1d->SetFillColor(kGreen-2);
    hsim1d->SetLineColor(kGreen-2);
    hsim1d->SetMarkerStyle(0);


    hsim1d->Draw("Ehist,same");
    
    TLegend leg(0.65, 0.75, 0.96, 0.92);
    leg.AddEntry(hunfold1d,"(Data - BG) Unfolded ","lp");
    leg.AddEntry(hsim1d,"MCatNLO Simulation","f");
    leg.SetBorderSize(0);
    leg.SetLineStyle(0);
    leg.SetTextFont(42);
    leg.SetFillStyle(0);
    leg.Draw();
    cmsstyle->CMSPreliminary();
    cmsstyle->ShowLumi(lumiString);
    cmsstyle->TextforCMS(resulttext, 2.5);
    cmsstyle->TextforCMS("dilepton", 3.5);
    saveAs(gPad, (TString("output/")+fileSuffix+"1dclassyunfold_unfolded"));
    
  }
  
  
  
  
  return 0;
}


int main()
{
  classy_unfold1d();
  return 0;
}
