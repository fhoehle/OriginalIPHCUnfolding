#include "cmsstyle.hpp"
#include "fitresults.hpp"
#include "helpers.hpp"
#include "samplelocations.hpp"


#include "binning.hpp"

// WARNING make sure variables in binning.hpp are set correctly


#include "myunfold_class.hpp"
//#include "myunfold_class1d.hpp"


int classy_unfold()
{
  TH1::SetDefaultSumw2(true);

  TFile* inputfile;
  


  inputfile = new TFile(sl_data_comb);


  
//  TFile* matrixfile = new TFile(TString("output/")+fileSuffix+"migmatrix.root");
//  TFile* preselfile = new TFile(TString("output/")+fileSuffix+"seleff.root");
  TFile* matrixfile = new TFile(TString("output/")+"migmatrix.root");
  TFile* preselfile = new TFile(TString("output/")+"seleff.root");

  MyUnfold myunfold;



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
  
  H2rec myh2rec("hrec", xedgesrec, yedgesrec);
  H2gen myh2gen("hgen", xedgesgen, yedgesgen);
  
  TH1F* sensrecoplots[nbinsmafter];
  for(int i=0; i<nbinsmafter; i++)
  {
    sensrecoplots[i] = new TH1F("","", 100, yedgesgen[0][0], yedgesgen[0][nbinsetaafter]);
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
  

  int nentries_toskip = 0; 
                               
  int nentries_toconsider = nentries;
  

  for(int i=nentries_toskip; i<nentries_toconsider; i++)
  {

    
    tuple->GetEntry(i);
    
    double weightasymmetry = 1;
  
    myh2gen.fill(fabs(mttbargen), diffabsetagen, weight*weightasymmetry);
    myh2rec.fill(fabs(mttbarrec), diffabsetarec, weight*weightasymmetry);
    
    const int genbinm = myh2gen.findXBin(fabs(mttbarrec));
    sensrecoplots[genbinm]->Fill(diffabsetarec, weight*weightasymmetry);
  }
  
  hgen = myh2gen.toTH2F();
  hrec = myh2rec.toTH2F();
  


  

  
  TH2F* hunfold = myunfold.unfoldHisto(hrec, workingondata); // where the magic happens...
  
   
  TString style = "COLZ";
  
  hrec->SetMinimum(0);
  hunfold->SetMinimum(0);
  hunfold->SetMaximum(hunfold->GetMaximum()*1.1);
  
  cmsstyle->setup_style_2D(hunfold, labelOfXAxisVar, labelOfSensVar);
  cmsstyle->setup_style_2D(hrec, labelOfXAxisVar, labelOfSensVar);
  cmsstyle->setup_style_2D(hpresel, labelOfXAxisVar, labelOfSensVar);
  hunfold->SetNdivisions(505);
  hrec->SetNdivisions(505);
  hpresel->SetNdivisions(505);
  

  


  //hrecwithoutbg->Draw(style);
  //saveAs(gPad, (TString("output/")+fileSuffix+"classyunfold_hrecwithoutbg"));


  hunfold->Draw(style);
  saveAs(gPad, (TString("output/")+fileSuffix+"classyunfold_hunfold"));
  hrec->Draw(style);
  saveAs(gPad, (TString("output/")+fileSuffix+"classyunfold_hrec"));
  if(myunfold.lcurve)
  {
    myunfold.lcurve->Draw("AP");
    saveAs(gPad, (TString("output/")+fileSuffix+"classyunfold_lcurve"));
  }

  myunfold.correlationmatrix->Draw("COLZ");
  saveAs(gPad, (TString("output/")+fileSuffix+"classyunfold_correlationmatrix"));

  myunfold.errormatrix->Draw("COLZ");
  saveAs(gPad, (TString("output/")+fileSuffix+"classyunfold_errormatrix"));
  
  for(int i=0; i<nbinsmafter; i++)
  {
    cmsstyle->setup_style(sensrecoplots[i], labelOfSensVar, "N", 1, 0, 0);
    sensrecoplots[i]->Draw("E");
    saveAs(gPad, (TString("output/")+fileSuffix+ TString::Format("classyunfold_sensrecoplotbin%d",i) ));
  }
  
  /// calc asy
  double pos, neg, asy;
  double posrec, negrec, asyrec;
  
  cout << endl << endl;
   
  
  pos = hrec->Integral(1, nbinsm, nbinseta/2+1, nbinseta);
  neg = hrec->Integral(1, nbinsm, 1, nbinseta/2);
  asy = (pos-neg)/(pos+neg);
  cout << "reconstructed raw asy: " << asy << endl;
  cout << "err for that: " << asymmetry_error_naive(pos, neg) << endl;
  
  
  
  //diffvis for unfolded data
  TH1F* diffvis = new TH1F("diffvis", "diffvis", nbinsmafter, xedgesgen);
  TH1F* diffvistrue = new TH1F("diffvistrue", "diffvistrue", nbinsmafter, xedgesgen);
  TH1F* mvis = new TH1F("mvis", "mvis", nbinsmafter, xedgesgen);
  cmsstyle->setup_style(diffvis, labelOfXAxisVar, "A_{C}", 1, 0, 0);
  cmsstyle->setup_style(diffvistrue, labelOfXAxisVar, "A_{C}", 2, 0, 0);
  cmsstyle->setup_style(mvis, labelOfXAxisVar, "N/GeV", 2, 0, 0);
  
  for(int i=0; i<nbinsmafter; i++)
  {
    pos = hunfold->Integral(i+1, i+1, nbinsetaafter/2+1, nbinsetaafter);
    neg = hunfold->Integral(i+1, i+1, 1, nbinsetaafter/2);
    asy = (pos-neg)/(pos+neg);

    posrec = hrec->Integral(i+1, i+1, nbinsetaafter/2+1, nbinsetaafter);
    negrec = hrec->Integral(i+1, i+1, 1, nbinsetaafter/2);
    asyrec = (posrec-negrec)/(posrec+negrec);



    //const double err = asymmetry_error_naive(pos, neg);
    const double err = asymmetryerror_afterunfolding_2d_onexbin(myunfold.errormatrix, nbinsafter, pos, neg, nbinsetaafter, i);
    
    const double truepos = hpresel->Integral(i+1, i+1, nbinsetaafter/2+1, nbinsetaafter);
    const double trueneg = hpresel->Integral(i+1, i+1, 1, nbinsetaafter/2);
    const double trueasy = (truepos-trueneg)/(truepos+trueneg);
    const double trueerr = asymmetry_error_naive(truepos, trueneg);
    

    cout << "reconstructed raw asy in mass bin " << i << ": " << asyrec << endl;
    cout << "err for that: " << asymmetry_error_naive(posrec, negrec) << endl;


    cout << "unfolded asy in mass bin " << i << ": " << asy << endl;

    //if (i == 0) cout << "shift on the unfolded asy in mass bin " << i << ": " << asy + 0.0934648 << endl;
    //if (i == 1) cout << "shift on the unfolded asy in mass bin " << i << ": " << asy - 0.105143 << endl;
    //if (i == 2) cout << "shift on the unfolded asy in mass bin " << i << ": " << asy + 0.166374 << endl;

    cout << "err for that: " << err << endl ;

    //if (i == 0) cout << "stat. error on syst. error in mass bin " << i << ": " << err + 0.115749 << endl;
    //if (i == 1) cout << "stat. error on syst. error in mass bin " << i << ": " << err + 0.132972 << endl;
    //if (i == 2) cout << "stat. error on syst. error in mass bin " << i << ": " << err + 0.0824862 << endl;
       
    const double nmass = hunfold->Integral(i+1, i+1, 1, nbinsetaafter);
    
    const int nbin = i+1; 
    diffvis->SetBinContent(nbin, asy);
    diffvis->SetBinError(nbin, err);
    diffvistrue->SetBinContent(nbin, trueasy);
    diffvistrue->SetBinError(nbin, trueerr);
    mvis->SetBinContent(nbin, nmass/(xedgesgen[i+1]-xedgesgen[i]));
    mvis->SetBinError(nbin, sqrt(nmass/(xedgesgen[i+1]-xedgesgen[i])));
    
    cout << endl;
  }
  
  {
    double marg = canv1->GetRightMargin();
    canv1->SetRightMargin(0.04);
    
    
   
    // theory histos
    TH1F* theory = 0;
    TH1F* exo = 0;
    
    if(!TString(nameOfXAxisVarShort).CompareTo("m"))
    {
      theory = (TH1F*)diffvis->Clone("theorydiffvis");
      cmsstyle->setup_style(theory, labelOfXAxisVar, "A_{C}", 4, 0, 0);
      theory->SetBinContent(1, 0.01032);
      theory->SetBinContent(2, 0.01222);
      theory->SetBinContent(3, 0.01519);
      theory->SetBinError(1, 0.00035);
      theory->SetBinError(2, 0.0006);
      theory->SetBinError(3, 0.00037);
/*
      theory = (TH1F*)diffvis->Clone("theorydiffvis");
      cmsstyle->setup_style(theory, labelOfXAxisVar, "A_{C}", 4, 0, 0);
      theory->SetBinContent(1, 0.0077);
      theory->SetBinContent(2, 0.0112);
      theory->SetBinContent(3, 0.0157);
      theory->SetBinError(1, 0.0003);
      theory->SetBinError(2, 0.0004);
      theory->SetBinError(3, 0.0006);
*/

      exo = (TH1F*)diffvis->Clone("exodiffvis");
      cmsstyle->setup_style(exo, labelOfXAxisVar, "A_{C}", 8, 0, 0);
      exo->SetBinContent(1, 0.0032);
      exo->SetBinContent(2, 0.0090);
      exo->SetBinContent(3, 0.0490);
      exo->SetBinError(1, 0.0032*0.07);
      exo->SetBinError(2, 0.0090*0.09);
      exo->SetBinError(3, 0.0490*0.15);
    }

    
    if(!TString(nameOfXAxisVarShort).CompareTo("y"))
    {
      exo = (TH1F*)diffvis->Clone("exodiffvis");
      cmsstyle->setup_style(exo, labelOfXAxisVar, "A_{C}", 8, 0, 0);
      exo->SetBinContent(1, 0.022);
      exo->SetBinContent(2, 0.014);
      exo->SetBinContent(3, 0.021);
      exo->SetBinError(1, 0.022*0.10);
      exo->SetBinError(2, 0.014*0.12);
      exo->SetBinError(3, 0.021*0.12);

      theory = (TH1F*)diffvis->Clone("theorydiffvis");
      cmsstyle->setup_style(theory, labelOfXAxisVar, "A_{C}", 4, 0, 0);
      theory->SetBinContent(1, 0.01089);
      theory->SetBinContent(2, 0.00922);
      theory->SetBinContent(3, 0.02276);
      theory->SetBinError(1, 0.00046);
      theory->SetBinError(2, 0.0003);
      theory->SetBinError(3, 0.00090);

/*
      theory = (TH1F*)diffvis->Clone("theorydiffvis");
      cmsstyle->setup_style(theory, labelOfXAxisVar, "A_{C}", 4, 0, 0);
      theory->SetBinContent(1, 0.0030);
      theory->SetBinContent(2, 0.0086);
      theory->SetBinContent(3, 0.0235);
      theory->SetBinError(1, 0.0002);
      theory->SetBinError(2, 0.0004);
      theory->SetBinError(3, 0.0010);
*/

    }
    

    if(!TString(nameOfXAxisVarShort).CompareTo("pt"))
    {
      theory = (TH1F*)diffvis->Clone("theorydiffvis");
      cmsstyle->setup_style(theory, labelOfXAxisVar, "A_{C}", 4, 0, 0);
      theory->SetBinContent(1, 0.01626);
      theory->SetBinContent(2, -0.00585);
      theory->SetBinContent(3, -0.00334);
      theory->SetBinError(1, 0.00072);
      theory->SetBinError(2, 0.00039);
      theory->SetBinError(3, 0.00022);
    }






    // asymmetric syst. error graphs
    TGraphAsymmErrors* syserrs = 0;   
    float sysup[nbinsmafter];
    float sysdown[nbinsmafter];
    
    if(!TString(nameOfXAxisVarShort).CompareTo("m"))
    {
      sysup[0] = 0.073;
      sysup[1] = 0.093;
      sysup[2] = 0.096;
      sysdown[0] = 0.036;
      sysdown[1] = 0.056;
      sysdown[2] = 0.052;
/*
      sysup[0] = 0.017;
      sysup[1] = 0.016;
      sysup[2] = 0.024;
      sysdown[0] = 0.017;
      sysdown[1] = 0.016;
      sysdown[2] = 0.024;
*/

    }
    if(!TString(nameOfXAxisVarShort).CompareTo("pt"))
    {
      sysup[0] = 0.051;
      sysup[1] = 0.133;
      sysup[2] = 0.052;
      sysdown[0] = 0.022;
      sysdown[1] = 0.075;
      sysdown[2] = 0.066;
/*
      sysup[0] = 0.020;
      sysup[1] = 0.010;
      sysup[2] = 0.019;
      sysdown[0] = 0.020;
      sysdown[1] = 0.010;
      sysdown[2] = 0.019;
*/
    }
    if(!TString(nameOfXAxisVarShort).CompareTo("y"))
    {
      sysup[0] = 0.043;
      sysup[1] = 0.137;
      sysup[2] = 0.085;
      sysdown[0] = 0.110;
      sysdown[1] = 0.138;
      sysdown[2] = 0.008;
/*
      sysup[0] = 0.010;
      sysup[1] = 0.010;
      sysup[2] = 0.022;
      sysdown[0] = 0.010;
      sysdown[1] = 0.010;
      sysdown[2] = 0.022;
*/
    }
      
      
    {
      syserrs = new TGraphAsymmErrors(diffvis);
      syserrs->SetLineWidth(2);
      syserrs->SetMarkerStyle(0);
      
      for(int i=0; i<nbinsmafter; i++)
      {
       const double staterr = diffvis->GetBinError(i+1); 
        syserrs->SetPointEYhigh(i, sqrt( pow(staterr,2) + pow(sysup[i],2) ));
        syserrs->SetPointEYlow(i, sqrt( pow(staterr,2) + pow(sysdown[i],2) ));                
      }      
      
    }
    
    
    {
      diffvis->SetMaximum();
      diffvis->SetMinimum();
      double max = diffvis->GetMaximum();
      double min = diffvis->GetMinimum();
      
      if(theory)
      {
        theory->SetMaximum();
        theory->SetMinimum();
        max = max > theory->GetMaximum() ? max : theory->GetMaximum();
        min = min < theory->GetMinimum() ? min : theory->GetMinimum();
        theory->SetMarkerStyle(0);
      }
      else
      {
        diffvistrue->SetMaximum();
        diffvistrue->SetMinimum();
        max = max > diffvistrue->GetMaximum() ? max : diffvistrue->GetMaximum();
        min = min < diffvistrue->GetMinimum() ? min : diffvistrue->GetMinimum();
        diffvistrue->SetMarkerStyle(0);
      }
      
      if(exo)
      {
        exo->SetMaximum();
        exo->SetMinimum();
        max = max > exo->GetMaximum() ? max : exo->GetMaximum();
        min = min < exo->GetMinimum() ? min : exo->GetMinimum();
        exo->SetMarkerStyle(0);
      }
      
      const double offset = max - min;
      diffvis->SetMaximum(diffvis->GetMaximum()+offset*1.15);
      
      
      if(!TString(nameOfXAxisVarShort).CompareTo("y"))
      {
        diffvis->SetMaximum(diffvis->GetMaximum()+offset*0.6);
        diffvis->SetMinimum(diffvis->GetMinimum()-offset*0.8);
      }
      
      if(!TString(nameOfXAxisVarShort).CompareTo("pt"))
      {
        diffvis->SetMinimum(diffvis->GetMinimum()-offset*0.5);
      }
      
      
    }
    
    
    diffvis->Draw("E1hist");

    if(!theory)
      diffvistrue->Draw("Ehist,same");
    if(exo)
      exo->Draw("Ehist,same");
    if(theory)
      theory->Draw("Ehist,same");

    if(syserrs)
      syserrs->Draw("p");

    diffvis->Draw("E1hist,same"); // second time to overwrite syst. error bars

    if(theory)
      theory->Draw("Ehist,same"); // to see the tiny error bars where they can be seen
    if(!theory)
      diffvistrue->Draw("Ehist,same");
    if(exo)
      exo->Draw("Ehist,same");

    
    double leg_right = 0.96;
    double leg_left = 0.65;
    double leg_bottom = 0.75;
    if(!TString(nameOfXAxisVarShort).CompareTo("m"))
    {
      leg_right = 0.90;
    }
    if(!TString(nameOfXAxisVarShort).CompareTo("pt"))
    {
      leg_left = 0.55;
      leg_bottom = 0.75;
    }
    if(!TString(nameOfXAxisVarShort).CompareTo("y"))
    {
      leg_left = 0.55;
      leg_bottom = 0.72;
    }
    TLegend leg(leg_left, leg_bottom, leg_right, 0.92);
    leg.AddEntry(diffvis,"Data","lp");
    if(!theory)
      leg.AddEntry(diffvistrue,"MCatNLO Simulation","l");
    if(exo)
      leg.AddEntry(exo,"EFT","l");
    if(theory)
      leg.AddEntry(theory,"NLO prediction","l");
    leg.SetBorderSize(0);
    leg.SetLineStyle(0);
    leg.SetTextFont(42);
    leg.SetFillStyle(0);
    leg.Draw();
    cmsstyle->CMSPreliminary();
    cmsstyle->ShowLumi(lumiString);
    cmsstyle->TextforCMS("dilepton", 2.5);
    saveAs(gPad, (TString("output/")+fileSuffix+"classyunfold_diffvis"));
    //mvis->SetMaximum(20000);
    //mvis->SetMinimum(-5000);
    mvis->Draw("Ehist");
    saveAs(gPad, (TString("output/")+fileSuffix+"classyunfold_mvis"));
    
    canv1->SetRightMargin(marg);
  }
  
  
  
  pos = hunfold->Integral(1, nbinsmafter, nbinsetaafter/2+1, nbinsetaafter);
  neg = hunfold->Integral(1, nbinsmafter, 1, nbinsetaafter/2);
  asy = (pos-neg)/(pos+neg);

  {
    double asy_error = asymmetryerror_afterunfolding_2d(myunfold.errormatrix, nbinsafter, pos, neg, nbinsetaafter);
    cout << "unfolded asy: " << asy << endl;
    cout << "err for that: " << asy_error << endl;

  }






  // some correlation plots to follow: change palette temporarily for those
  
  Int_t MyPalette[100];
  Double_t r[]    = {0.0, 1.0, 1.0};
  Double_t g[]    = {0.0, 1.0, 0.0};
  Double_t b[]    = {1.0, 1.0, 0.0};
  Double_t stop[] = {.0, 0.5, 1.0};
  Int_t FI = TColor::CreateGradientColorTable(3, stop, r, g, b, 100);
  for (int i=0;i<100;i++)
    MyPalette[i] = FI+i;
  
  gStyle->SetPalette(100, MyPalette); 
  myunfold.correlationmatrix->SetMaximum(1);
  myunfold.correlationmatrix->SetMinimum(-1);
  myunfold.correlationmatrix->Draw("COLZ");
  saveAs(gPad, (TString("output/")+fileSuffix+"/"+"classyunfold_cormatrix"));
  
  gStyle->SetPalette(1); 
  
  
  /// find correlaton between bins of unfolded asymmetry (in linear approximation)
  /// (see http://en.wikipedia.org/wiki/Propagation_of_uncertainty)
  // construct jacobian of the asymmetry functions
  // each row corresponds to an asymmetry bin; the columns contain the partial derivatives of this asymmetry function
  TMatrix jacobian(nbinsmafter, nbinsafter);
  
  for(int row=0; row<nbinsmafter; row++)
  {
    for(int col=0; col<nbinsafter; col++)
    {
      const int secondarybin = col/nbinsetaafter;
      
      if(secondarybin != row)
        continue; // partial derivative is zero anyway
      
      const int sensvarbin = col - secondarybin*nbinsetaafter;
      
      const bool is_pos_sensvar = sensvarbin >= nbinsetaafter/2;
      
      const double npos = hunfold->Integral(row+1, row+1, nbinsetaafter/2+1, nbinsetaafter);
      const double nneg = hunfold->Integral(row+1, row+1, 1, nbinsetaafter/2);
      
      if(is_pos_sensvar)
      {
        jacobian[row][col] = 2 * nneg / pow(npos+nneg, 2);
      }
      else
      {
        jacobian[row][col] = -2 * npos / pow(npos+nneg, 2);
      }
      
    }    
  }
  
  // get covmatrix as TMatrix
  TMatrix cov(nbinsafter, nbinsafter);
  
  for(int row=0; row<nbinsafter; row++)
  {
    for(int col=0; col<nbinsafter; col++)
    {
      cov[row][col] = myunfold.errormatrix->GetBinContent(row+1, col+1);
    }
  }
  
  // transform covmatrix using jacobian, yielding cov for asy bins
  TMatrix temp(nbinsafter, nbinsmafter);
  temp.MultT(cov, jacobian);
  
  TMatrix asycovmat(nbinsmafter, nbinsmafter);
  asycovmat.Mult(jacobian, temp);
    
  // transform back to th2, calculate correlation matrix, draw both
  TH2F* hasycov = new TH2F("hasycov", "hasycov", nbinsmafter, 0.5, nbinsmafter+0.5, nbinsmafter, 0.5, nbinsmafter+0.5);
  TH2F* hasycor = new TH2F("hasycor", "hasycor", nbinsmafter, 0.5, nbinsmafter+0.5, nbinsmafter, 0.5, nbinsmafter+0.5);
  
  for(int row=0; row<nbinsmafter; row++)
  {
    for(int col=0; col<nbinsmafter; col++)
    {
      hasycov->SetBinContent(row+1, col+1, asycovmat[row][col]);
      
      const double cor = asycovmat[row][col] / sqrt(asycovmat[row][row]*asycovmat[col][col]);
      hasycor->SetBinContent(row+1, col+1, cor);
    }
  }
  
  TString painttextformat = gStyle->GetPaintTextFormat();
  
  gStyle->SetPaintTextFormat("1.5f");
  hasycov->SetMarkerSize(1.3);
  hasycov->Draw("colz text");
  saveAs(gPad, (TString("output/")+fileSuffix+"/"+"classyunfold_asycovmatrix"));
  
  gStyle->SetPalette(100, MyPalette); 
  gStyle->SetPaintTextFormat("1.2f");
  hasycor->SetMarkerSize(2);
  hasycor->SetMaximum(1);
  hasycor->SetMinimum(-1);
  hasycor->Draw("colz text");
  saveAs(gPad, (TString("output/")+fileSuffix+"/"+"classyunfold_asycormatrix"));
  
  gStyle->SetPaintTextFormat(painttextformat.Data());
  gStyle->SetPalette(1);   









  
  return 0;
}


int main()
{
  classy_unfold();
  return 0;
}
