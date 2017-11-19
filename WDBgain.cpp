#include <string>
#include <iostream>
#include <fstream>
#include <sstream>  // Required for stringstream
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TBranch.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <math.h>
#include <TSystem.h>
#include <TMath.h>

void getPMdata(Int_t run);
Double_t vargaus(Double_t *x,Double_t *par);
Int_t arraysearch(std::vector<Int_t> array, Int_t value);

#define Npm 4760
#define Nmppc 4092
#define NRow 93
#define NLine 44

std::vector<Int_t> PMdata;
std::vector<std::vector<Double_t>> ChargeMean;

Int_t Nvalidpm;
Int_t Nvalidmppc=0;

Int_t rangemode=0;
Int_t runarray[3]={308580,308581,308582};//PMT 9,Nov
Double_t WDBGainarr[3]={100,10,1};

Int_t runnum=(sizeof(runarray)/sizeof(runarray[0]));

//Double_t varrange[2]={0,0.0001};//HV 47V

//Double_t posoff=1.0;
//Double_t chargerange[2]={-6.0+posoff,posoff};
Double_t rangemin;
Double_t rangemax;
//Double_t GainMax=5.0*pow(10,8);

void WDBgain() {

  getPMdata(runarray[0]);
  ChargeMean.resize(Nvalidpm);
  TCanvas *canvas1 = new TCanvas("canvas1","Gain",600,600);
  TCanvas *canvas2 = new TCanvas("canvas2","Ch Gain",600,600);
  TCanvas *canvas3 = new TCanvas("canvas3","WDB gain",600,600);

  TGraph* grWDBgain = new TGraph();
  TGraph* grWDBgain1st = new TGraph();
  TGraph* grWDBgain2nd = new TGraph();
  TH1D* hGain1st= new TH1D("gain 1st","gain 1st;Gain;channels",100,0,20);
  TH1D* hGain2nd= new TH1D("gain 2nd","gain 2nd;Gain;channels",100,0,20);
  TClonesArray *linearGarr =new TClonesArray("TGraphErrors",Nvalidpm);
  TClonesArray *linFarr = new TClonesArray("TF1",Nvalidpm);
  for(int i=0;i<Nvalidpm;i++){
    ChargeMean[i].resize(runnum);
    int ch= PMdata[i];
    new((*linearGarr)[i]) TGraphErrors(runnum);
  }
  TH1D* testH=new TH1D("test","test title",100,-0.1,0);
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
      rangemin=-0.2;
      switch (run) {
        case 0:
        rangemax=20.0;
        break;
        case 1:
        rangemax=2.0;
        break;
        case 2:
        rangemax=0.2;
        break;
        default:
        rangemax=20.0;
        break;
      }
      new((*chargeHarr)[i]) TH1D(Form("charge histogram %d;charge[10^{9}e];events",ch),Form("charge histogram title %d",ch),400,rangemin,rangemax);
      new((*fitFarr)[i]) TF1(Form("gausvar%d",ch),vargaus,-20,0.5,3);
    }

    //Filling phase
    Int_t Nevent=br2->GetEntries();
    for (int eve = 0; eve < Nevent; ++eve) {
      br2->GetEntry(eve);
      // Double_t testcharge = ((MEGXECWaveformAnalysisResult*)(XECWaveformAnalysisResult->At(352)))->GetchargeAt(rangemode);
      // //std::cout<<"test charge: "<<testcharge<<std::endl;
      // if (run==0) {
      // testH->Fill(testcharge);
      // }

      for(int i=0;i<Nvalidpm;i++){
        double charge=0;
        int ch=PMdata[i];
        charge = -((MEGXECWaveformAnalysisResult*)(XECWaveformAnalysisResult->At(ch)))->GetchargeAt(rangemode);
        ((TH1D*)((*chargeHarr)[i]))->Fill(charge);
      }
    }

    //Fitting region
    for(int i=0;i<Nvalidpm;i++){
      int ch=PMdata[i];
      Double_t defvar=((TH1D*)((*chargeHarr)[i]))->GetRMS()*((TH1D*)((*chargeHarr)[i]))->GetRMS();
      ((TF1*)(*fitFarr)[i])->SetParameters(100,defvar,((TH1D*)((*chargeHarr)[i]))->GetMean());
      ((TH1D*)((*chargeHarr)[i]))->Fit(Form("gausvar%d",ch),"NQ","",0,20);
      // if(i==0&&run==0){
      //   ((TH1D*)((*chargeHarr)[i]))->Draw();
      // }
      Double_t Mean=((TF1*)((*fitFarr)[i]))->GetParameter(2);
      Double_t Meanerr=((TF1*)((*fitFarr)[i]))->GetParError(2);
      // Double_t Mean=((TH1D*)((*chargeHarr)[i]))->GetMean();
      // Double_t Meanerr=((TH1D*)((*chargeHarr)[i]))->GetMeanError();
      ((TGraphErrors*)((*linearGarr)[i]))->SetPoint(run,WDBGainarr[run],Mean);
      ((TGraphErrors*)((*linearGarr)[i]))->SetPointError(run,0,Meanerr);
      if(run==1){
        grWDBgain->SetPoint(grWDBgain->GetN(),ch,Mean);
      }
      ChargeMean[i][run]=Mean;
    }
  }
  for (int i = 0; i < Nvalidpm; i++) {
    int ch=PMdata[i];
    if(ChargeMean[i][0]!=0&&ChargeMean[i][1]!=0&&ch<Nmppc){
      Double_t Gainfirst= ChargeMean[i][1]/ChargeMean[i][2];
      Double_t Gainsecond= ChargeMean[i][0]/ChargeMean[i][1];
      hGain1st->Fill(Gainfirst);
      hGain2nd->Fill(Gainsecond);
      grWDBgain1st->SetPoint(grWDBgain1st->GetN(),ch,Gainfirst);
      grWDBgain2nd->SetPoint(grWDBgain2nd->GetN(),ch,Gainsecond);
      if(Gainfirst>12||Gainsecond>12||Gainfirst<0||Gainsecond<0){
        std::cout<<"outlier: "<<ch<<std::endl;
      }
    }

  }
  canvas1->cd();
  hGain2nd->SetLineColor(kBlue);
  hGain2nd->Fit("gaus","N");
  hGain2nd->Draw();
  hGain1st->SetLineColor(kRed);
  hGain1st->Fit("gaus","N");
  hGain1st->Draw("same");

  //testH->Draw();
  canvas2->cd();
  grWDBgain->SetMarkerStyle(20);
  grWDBgain->SetMarkerColor(kRed);
  grWDBgain->Draw("ap");
  canvas3->cd();
  grWDBgain1st->SetTitle("Gain;channels");
  grWDBgain1st->SetMarkerStyle(20);
  grWDBgain1st->SetMarkerColor(kRed);
  grWDBgain1st->SetMaximum(20);
  grWDBgain1st->SetMinimum(0);
  grWDBgain1st->Draw("ap");
  grWDBgain2nd->SetMarkerStyle(20);
  grWDBgain2nd->SetMarkerColor(kBlue);
  grWDBgain2nd->Draw("same p");
}

void WDBgain(int ch) {
  getPMdata(runarray[0]);
  if(arraysearch(PMdata,ch)==-1){
    std::cerr<<"Channel not found"<<std::endl;
    std::exit(EXIT_FAILURE);
  }

  TClonesArray *chargeHarr = new TClonesArray("TH1D",runnum);
  TClonesArray *fitFarr = new TClonesArray("TF1",Nvalidpm);
  TCanvas *canvas1 = new TCanvas("canvas1","Gain",600,600);
  TCanvas *canvas2 = new TCanvas("canvas2","Ch Gain",600,600);

  TH1D* testH=new TH1D("test","test title",100,-0.1,0);
  TGraphErrors* gr=new TGraphErrors();
  for(int run = 0; run < runnum; run++) {
    //file read
    TFile *frec = new TFile(Form("$(MEG2SYS)/analyzer/recfiles/rec%06d.root", runarray[run]),"READ");
    TTree *rec = (TTree*)frec->Get("rec");
    TBranch* br2 = rec->GetBranch("xecwfcl");
    TClonesArray* XECWaveformAnalysisResult = new TClonesArray("MEGXECWaveformAnalysisResult");
    br2->SetAddress(&XECWaveformAnalysisResult);

    rangemin=-0.2;
    switch (run) {
      case 0:
      rangemax=20.0;
      break;
      case 1:
      rangemax=2.0;
      break;
      case 2:
      rangemax=0.2;
      break;
      default:
      rangemax=20.0;
      break;
    }
    new((*chargeHarr)[run]) TH1D(Form("charge histogram %d;charge[10^{9}e];events",ch),Form("charge histogram title %d",ch),400,rangemin,rangemax);
    new((*fitFarr)[run]) TF1(Form("gausvar%d",ch),vargaus,-20,0.5,3);

    //Filling phase
    Int_t Nevent=br2->GetEntries();
    for (int eve = 0; eve < Nevent; ++eve) {
      br2->GetEntry(eve);


      Double_t charge = -((MEGXECWaveformAnalysisResult*)(XECWaveformAnalysisResult->At(ch)))->GetchargeAt(rangemode);
      ((TH1D*)((*chargeHarr)[run]))->Fill(charge);

    }
    Double_t defvar=((TH1D*)((*chargeHarr)[run]))->GetRMS()*((TH1D*)((*chargeHarr)[run]))->GetRMS();
    ((TF1*)(*fitFarr)[run])->SetParameters(100,defvar,((TH1D*)((*chargeHarr)[run]))->GetMean());
    ((TH1D*)((*chargeHarr)[run]))->Fit(Form("gausvar%d",ch),"NQ","",-10,0);
    Double_t Mean=((TF1*)((*fitFarr)[run]))->GetParameter(2);
    Double_t Meanerr=((TF1*)((*fitFarr)[run]))->GetParError(2);
    gr->SetPoint(run,WDBGainarr[run],Mean);
    gr->SetPointError(run,0,Meanerr);
  }
  canvas1->cd();
  canvas1->SetLogy();
  canvas1->SetLogx();
  gr->SetTitle("Charge;WDB Gain;|Q|[10^{9}e]");
  gr->SetMarkerColor(kRed);
  gr->SetMarkerStyle(20);
  gr->Draw("ap");


  //testH->Draw();
  canvas2->cd();

  for (int i = 0; i < runnum; i++) {
    ((TH1D*)((*chargeHarr)[i]))->SetLineColor(i+2);
    if(i==0){
      ((TH1D*)((*chargeHarr)[i]))->SetTitle("charge histogram;|Q|[10^{9}e];events");
      ((TH1D*)((*chargeHarr)[i]))->Draw();
    }else{
      ((TH1D*)((*chargeHarr)[i]))->Draw("same");
    }

  }

}


void getPMdata(Int_t run) {
  TFile *frec = new TFile(Form("$(MEG2SYS)/analyzer/recfiles/rec%06d.root", run),"READ");
  TClonesArray* pmrhArray = (TClonesArray*)frec->Get("XECPMRunHeader");
  MEGXECPMRunHeader *pmrh = 0;
  for (int iPM = 1; iPM < Npm; iPM++) {
    pmrh = (MEGXECPMRunHeader*)(pmrhArray->At(iPM));
    Int_t DCid=pmrh->GetDRSChipID();
    if (DCid>=0) {
      PMdata.push_back(iPM);
      if(iPM<Nmppc){
        Nvalidmppc+=1;
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
