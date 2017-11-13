#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#ifndef __CINT__
#include "ROMETreeInfo.h"
#endif
#define Npm 4760
#define NSiPM 4092
#define NRow 93
#define NLine 44

Double_t rangemax4SiPM;
Double_t range4SiPM;
Double_t rangemax4PMT;
Double_t range4PMT;
Double_t GainMax;
Double_t rangemin;
Double_t rangemax;
Int_t rangemode;
std::vector<Int_t> PMdata;
Int_t Nvalidpm;

//Int_t runarray[8]={308583,308584,308585,308586,308587,308588,308589,308590};
//Int_t runarray[]={308540,308541,308542,308543,308544,308545,308546,308547};
Int_t runarray[]={308525,308526,308527,308528,308529,308530,308531,308532};
Int_t runnum=(sizeof(runarray)/sizeof(runarray[0]));

void Draw(Int_t *runarray,Bool_t PMTanalysis, Bool_t SiPManalysis);
void Detail(Int_t *runarray, int ch);
void getPMdata(Int_t run,Bool_t PMTanalysis, Bool_t SiPManalysis);
Double_t vargaus(Double_t *x,Double_t *par);
Int_t arraysearch(std::vector<Int_t> array, Int_t value);
void Analysing(Int_t *runarray);
void ModeDeclear(Bool_t PMTanalysis,Bool_t SiPManalysis);

void strongLED_online(){
  rangemode=0;
  rangemax4SiPM=2.0;
  range4SiPM=16.0;
  rangemax4PMT=0.5;
  range4PMT=5.0;
  GainMax=5.0*pow(10,8);
  Draw(runarray,true,false);
  //Detail(runarray,4141);
}

void Draw(Int_t *runarray,Bool_t PMTanalysis,Bool_t SiPManalysis) {
  Analysing(runarray);
  ModeDeclear(PMTanalysis,SiPManalysis);
  getPMdata(runarray[0],PMTanalysis,SiPManalysis);
  TCanvas *canvas1 = new TCanvas("canvas1","Gain",600,600);
  TGraphErrors *chgain = new TGraphErrors(Nvalidpm);
  TClonesArray *linearGarr =new TClonesArray("TGraphErrors",Nvalidpm);
  TClonesArray *linFarr = new TClonesArray("TF1",Nvalidpm);
  for(int i=0;i<Nvalidpm;i++){
    int ch = PMdata[i];
    new((*linearGarr)[i]) TGraphErrors(runnum);
    new((*linFarr)[i]) TF1(Form("linear%d",ch),"[0]*x+[1]");
  }
  for(int run = 0; run < runnum; run++) {
    //file read
    TFile *frec = new TFile(Form("$(MEG2SYS)/analyzer/recfiles/rec%06d.root", runarray[run]),"READ");
    TTree *rec = (TTree*)frec->Get("rec");
    TBranch* br2 = rec->GetBranch("xecwfcl");
    TClonesArray* XECWaveformAnalysisResult = new TClonesArray("MEGXECWaveformAnalysisResult");
    br2->SetAddress(&XECWaveformAnalysisResult);

    //create arrays
    TClonesArray *chargeHarr = new TClonesArray("TH1D",Nvalidpm);
    TClonesArray *fitFarr = new TClonesArray("TF1",Nvalidpm);
    for(int i=0;i<Nvalidpm;i++){
      int ch= PMdata[i];
      if (ch<NSiPM) {
        rangemin=rangemax4SiPM-range4SiPM;
        rangemax=rangemax4SiPM;
      }else{
        rangemin=rangemax4PMT-range4PMT;
        rangemax=rangemax4PMT;
      }
      new((*chargeHarr)[i]) TH1D(Form("charge histogram %d;charge[10^{9}e];events",ch),Form("charge histogram title %d",ch),500,rangemin,rangemax);
      new((*fitFarr)[i]) TF1(Form("gausvar%d",ch),vargaus,-20,0.5,3);
    }

    //Filling phase
    Int_t Nevent=br2->GetEntries();
    for (int eve = 0; eve < Nevent; ++eve) {
      br2->GetEntry(eve);
      for(int i=0;i<Nvalidpm;i++){
        double charge=0;
        int ch = PMdata[i];
        charge = ((MEGXECWaveformAnalysisResult*)(XECWaveformAnalysisResult->At(ch)))->GetchargeAt(rangemode);
        ((TH1D*)((*chargeHarr)[i]))->Fill(charge);
      }
    }
    //Fitting phase
    for(int i=0;i<Nvalidpm;i++){
      int ch=PMdata[i];
      Double_t defvar=((TH1D*)((*chargeHarr)[i]))->GetRMS()*((TH1D*)((*chargeHarr)[i]))->GetRMS();
      ((TF1*)(*fitFarr)[i])->SetParameters(100,defvar,((TH1D*)((*chargeHarr)[i]))->GetMean());
      ((TH1D*)((*chargeHarr)[i]))->Fit(Form("gausvar%d",ch),"NQ");
      Double_t Mean=-((TF1*)((*fitFarr)[i]))->GetParameter(2);
      Double_t Meanerr=-((TF1*)((*fitFarr)[i]))->GetParError(2);
      Double_t Variance=((TF1*)((*fitFarr)[i]))->GetParameter(1);
      Double_t Varerr=((TF1*)((*fitFarr)[i]))->GetParError(1);
      ((TGraphErrors*)((*linearGarr)[i]))->SetPoint(run,Mean,Variance);
      ((TGraphErrors*)((*linearGarr)[i]))->SetPointError(run,Meanerr,Varerr);
    }
  }

  for(int i=0;i<Nvalidpm;i++){
    chNum=PMdata[i];
    ((TGraphErrors*)((*linearGarr)[i]))->Fit(Form("linear%d",chNum),"Q");
    Double_t tmpgain=((TF1*)((*linFarr)[i]))->GetParameter(0)*pow(10,9);
    Double_t tmpgainerr=((TF1*)((*linFarr)[i]))->GetParError(0)*pow(10,9);
    Double_t tmpnoisevar=((TF1*)((*linFarr)[i]))->GetParameter(1);
    if(tmpgain<=0||tmpgainerr>tmpgain){
      std::cout<<"bad id:    "<<i<<"  ch:    "<<chNum<<"  gain error:    "<<tmpgainerr<<std::endl;
    }else{
      chgain->SetPoint(i,chNum,tmpgain);
      chgain->SetPointError(i,0,tmpgainerr);
    }
  }
  canvas1->cd();
  chgain->SetMarkerStyle(22);
  chgain->SetMarkerColor(2);
  chgain->SetMarkerSize(1);
  chgain->SetMaximum(GainMax);
  chgain->SetMinimum(0);
  chgain->SetTitle("Gain of each channel;channel;Gain");
  chgain->Draw("ap");
}

void Detail(Int_t *runarray, int ch) {
  TCanvas *canvas1 = new TCanvas("canvas1","Charge distribution",600,400);
  TCanvas *canvas2 = new TCanvas("canvas2","Linear Fitting",600,400);
  getPMdata(runarray[0],true,true);
  if(arraysearch(PMdata,ch)==-1){
    std::cerr<<"Channel not found"<<std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::cout<< "channel:    "<<ch<<std::endl;
  Analysing(runarray);
  TGraphErrors *gr =new TGraphErrors(runnum);
  TClonesArray *chargeHarr = new TClonesArray("TH1D",runnum);
  TClonesArray *fitFarr = new TClonesArray("TF1",runnum);
  Double_t MeanMax = 0;
  Double_t VarMax=0;
  canvas1->cd();
  for(int run = 0; run < runnum; run++) {
    new((*fitFarr)[run]) TF1(Form("gausvar%d",run),vargaus,-10,0.5,3);
    if (ch<NSiPM) {
      rangemin=rangemax4SiPM-range4SiPM;
      rangemax=rangemax4SiPM;
    }else{
      rangemin=rangemax4PMT-range4PMT;
      rangemax=rangemax4PMT;
    }
    new((*chargeHarr)[run]) TH1D(Form("charge histogram %d",ch),Form("charge histogram title %d;Charge[10^{9}e];Events",ch),200,rangemin,rangemax);
    TFile *frec = new TFile(Form("$(MEG2SYS)/analyzer/recfiles/rec%06d.root", runarray[run]),"READ");

    TTree *rec = (TTree*)frec->Get("rec");
    TBranch* br2 = rec->GetBranch("xecwfcl");
    TClonesArray* XECWaveformAnalysisResult = new TClonesArray("MEGXECWaveformAnalysisResult");
    br2->SetAddress(&XECWaveformAnalysisResult);

    Int_t Nevent=br2->GetEntries();
    for (int eve = 0; eve < Nevent; ++eve) {
      br2->GetEntry(eve);
      double charge=0;
      charge = ((MEGXECWaveformAnalysisResult*)(XECWaveformAnalysisResult->At(ch)))->GetchargeAt(rangemode);
      ((TH1D*)((*chargeHarr)[run]))->Fill(charge);
    }
    Double_t defvar=((TH1D*)((*chargeHarr)[run]))->GetRMS()*((TH1D*)((*chargeHarr)[run]))->GetRMS();
    ((TF1*)(*fitFarr)[run])->SetParameters(100,defvar,((TH1D*)((*chargeHarr)[run]))->GetMean());
    ((TH1D*)((*chargeHarr)[run]))->SetLineColor(run+2);
    ((TH1D*)((*chargeHarr)[run]))->Fit(Form("gausvar%d",run),"NQ");
    if(run==0){
      ((TH1D*)((*chargeHarr)[run]))->Draw();
      ((TF1*)((*fitFarr)[run]))->Draw("same");
    }else{
      ((TH1D*)((*chargeHarr)[run]))->Draw("same");
      ((TF1*)((*fitFarr)[run]))->Draw("same");
    }
    Double_t Mean=-((TF1*)((*fitFarr)[run]))->GetParameter(2);
    Double_t Meanerr=-((TF1*)((*fitFarr)[run]))->GetParError(2);
    Double_t Variance=((TF1*)((*fitFarr)[run]))->GetParameter(1);
    Double_t Varerr=((TF1*)((*fitFarr)[run]))->GetParError(1);
    if(run==runnum-1){
      MeanMax=Mean;
      VarMax=Variance;
    }
    gr->SetPoint(run,Mean,Variance);
    gr->SetPointError(run,Meanerr,Varerr);
  }

  canvas2->cd();
  gr->SetTitle("Mean vs Variance;Mean[10^{9}e];Variance");
  gr->GetXaxis()->SetLimits(0,MeanMax*1.2);
  gr->SetMaximum(VarMax*1.2);
  gr->SetMinimum(0);
  TF1 *linF= new TF1("linearfit","[0]*x+[1]",0,5);
  gr->Fit("linearfit","Q");
  Double_t Gain=linF->GetParameter(0);
  Double_t GainErr=linF->GetParError(0);
  Double_t NoiseVariance=linF->GetParameter(1);
  std::cout<<"Gain:    "<<Gain<<"+-"<<GainErr<<std::endl;
  std::cout<<"Noise Variance:    "<<NoiseVariance<<std::endl;
  gr->SetMarkerStyle(22);
  gr->SetMarkerColor(2);
  gr->SetMarkerSize(1);
  gr->Draw("ap");
}



void getPMdata(Int_t run,Bool_t PMTanalysis, Bool_t SiPManalysis) {
  TFile *frec = new TFile(Form("$(MEG2SYS)/analyzer/recfiles/rec%06d.root", run),"READ");
  TClonesArray* pmrhArray = (TClonesArray*)frec->Get("XECPMRunHeader");
  MEGXECPMRunHeader *pmrh = 0;
  for (int iPM = 1; iPM < Npm; iPM++) {
    pmrh = (MEGXECPMRunHeader*)(pmrhArray->At(iPM));
    Int_t DCid=pmrh->GetDRSChipID();
    if (DCid>=0) {
      if (PMTanalysis==true&&iPM>=NSiPM) {
        PMdata.push_back(iPM);
      }else if(SiPManalysis==true&&iPM<NSiPM){
        PMdata.push_back(iPM);
      }
    }
  }
  Nvalidpm=PMdata.size();
}

Double_t vargaus(Double_t *x,Double_t *par){
  Float_t xx=x[0];
  Double_t f = par[0]/sqrt(2*TMath::Pi()*par[1])*TMath::Exp(-pow(xx-par[2],2)/(2*par[1]));
  return f;
}

Int_t arraysearch(std::vector<Int_t> array, Int_t value){
  Int_t sizearr=array.size();
  Bool_t found=false;
  int i=0;
  for (i = 0; i < sizearr; i++) {
    if (array[i]==value) {
      found=true;
      break;
    }
  }
   if(found==false){
      i=-1;
   }
  return i;
}

void Analysing(Int_t *runarray){
  std::cout<<"Analysing: "<<std::endl;
  for (int run = 0; run < runnum; run++) {
std::cout<<"Run "<<runarray[run]<<std::endl;
  }
}

void ModeDeclear(Bool_t PMTanalysis,Bool_t SiPManalysis){
  if (PMTanalysis==true) {
    std::cout<<"PMT analysis enabled"<<std::endl;
  }else{
    std::cout<<"PMT analysis disabled"<<std::endl;
  }
  if (SiPManalysis==true) {
    std::cout<<"SiPM analysis enabled"<<std::endl;
  }else{
    std::cout<<"SiPM analysis disabled"<<std::endl;
  }
}
