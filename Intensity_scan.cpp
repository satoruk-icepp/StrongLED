#include <string>
#include <iostream>
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
#include <Intensity.h>

Int_t rangemode=0;
//Int_t runarray[6]={305670,305671,305672,305673,305674,305675};
//Int_t runarray[6]={306596,306597,306598,306599,306600,306601};//HV54, Config I 
//Int_t runarray[6]={306611,306612,306613,306614,306615,306616};//HV52,Config I
//Int_t runarray[6]={306617,306618,306619,306620,306621,306622};//HV51,Config I
//Int_t runarray[6]={306623,306624,306625,306626,306627,306628};//HV50,Config I
//Int_t runarray[6]={306629,306630,306631,306632,306633,306634};//HV49,Config I
//Int_t runarray[5]={306636,306637,306638,306639,306640};//HV48,Config I
//Int_t runarray[5]={306664,306665,306666,306667,306668};
//Int_t runarray[5]={306679,306680,306681,306682,306683};
//Int_t runarray[5]={306688,306689,306690,306691,306692};
Int_t runarray[6]={306902,306903,306905,306906,306907,306908};//OV7V,Config A
//Int_t runarray[6]={306909,306910,306911,306912,306913,306914};//OV5V,Config A
//Int_t runarray[6]={306922,306923,306924,306926,306927,306928};//OV1V,Config A
//Int_t runarray[5]={306674,306675,306676,306677,306678};//HV Ultrastrong
//Int_t runarray[5]={306701,306702,306703,306704,306705};
//Int_t runarray[3]={306703,306704,306705};
//Int_t runarray[6]={306902,306903,306905,306906,306907,306908};

Int_t runnum=(sizeof(runarray)/sizeof(runarray[0]));

//Double_t varrange[2]={0,0.0001};//HV 47V

//Double_t posoff=1.0;
//Double_t chargerange[2]={-6.0+posoff,posoff};
Double_t rangemax4SiPM=2.0;
Double_t range4SiPM=16.0;
Double_t rangemax4PMT=0.5;
Double_t range4PMT=5.0;
Double_t rangemin;
Double_t rangemax;
Double_t GainMax=5.0*pow(10,8);

void Intensity_scan() {
   TString foutname;
   foutname.Form("gain_%d-%d_%d.root",runarray[0],runarray[runnum-1],rangemode);
   TFile* fout= new TFile(foutname,"recreate");
   TTree* tout = new TTree("tree","tree");
   Int_t chNum;
   Double_t onegain;
   Double_t onenoise;
   Double_t onegainerr;
   tout->Branch("chNum",&chNum);
   tout->Branch("gain",&onegain);
   tout->Branch("noise",&onenoise);
   tout->Branch("gainerr",&onegainerr);
   getPMdata(runarray[0]);
   TCanvas *canvas1 = new TCanvas("canvas1","Gain Histogram",600,400);
   TCanvas *canvas2 = new TCanvas("canvas2","Ch Gain",600,400);
   TCanvas *canvas3 = new TCanvas("canvas3","Gain Map",400,600);
   TCanvas *canvas4 = new TCanvas("canvas4","Noise Map",400,600);
   std::vector<Double_t> Gain;
   std::vector<Double_t> GainErr;
   std::vector<Double_t> NoiseVariance;
   Gain.resize(Nvalidpm);
   GainErr.resize(Nvalidpm);
   NoiseVariance.resize(Nvalidpm);
   TH1D *gainhist=new TH1D("gain hist", "gain hist title;Gain;# of Channel",100,0,5.0*pow(10,8));
   TGraphErrors *chgain = new TGraphErrors(Nvalidpm);
   TClonesArray *linearGarr =new TClonesArray("TGraphErrors",Nvalidpm);
   TClonesArray *linFarr = new TClonesArray("TF1",Nvalidpm);
   for(int i=0;i<Nvalidpm;i++){
      int ch= PMdata[i];
      new((*linearGarr)[i]) TGraphErrors(runnum);
      new((*linFarr)[i]) TF1(Form("linear%d",ch),"[0]*x+[1]");
   }
   for(int run = 0; run < runnum; run++) {
      //file read
      TFile *frec = new TFile(Form("$(MEG2SYS)/analyzer/rec%06d.root", runarray[run]),"READ");
      TTree *rec = (TTree*)frec->Get("rec");
      TBranch* br2 = rec->GetBranch("xecwfcl");
      TClonesArray* XECWaveformAnalysisResult = new TClonesArray("MEGXECWaveformAnalysisResult");
      br2->SetAddress(&XECWaveformAnalysisResult);

      //create arrays
      TClonesArray *chargeHarr = new TClonesArray("TH1D",Nvalidpm);
      TClonesArray *fitFarr = new TClonesArray("TF1",Nvalidpm);

      for(int i=0;i<Nvalidpm;i++){
         int ch= PMdata[i];
         if (ch<Nmppc) {
            rangemin=rangemax4SiPM-range4SiPM;
            rangemax=rangemax4SiPM;
         }else{
            rangemin=rangemax4PMT-range4PMT;
            rangemax=rangemax4PMT;
         }
         new((*chargeHarr)[i]) TH1D(Form("charge histogram %d;charge[10^{9}e];events",ch),Form("charge histogram title %d",ch),100,rangemin,rangemax);
         new((*fitFarr)[i]) TF1(Form("gausvar%d",ch),vargaus,-10,0.5,3);
      }

      //Filling phase
      Int_t Nevent=br2->GetEntries();
      for (int eve = 0; eve < Nevent; ++eve) {
         br2->GetEntry(eve);
         for(int i=0;i<Nvalidpm;i++){
            double charge=0;
            int ch=PMdata[i];
            charge = ((MEGXECWaveformAnalysisResult*)(XECWaveformAnalysisResult->At(ch)))->GetchargeAt(rangemode);
            ((TH1D*)((*chargeHarr)[i]))->Fill(charge);
         }
      }
      //Fitting region
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
   TString outfilename;
   outfilename.Form("gain_%d-%d_%d.txt",runarray[0],runarray[runnum-1],rangemode);
   ofstream writing_file;
   writing_file.open(outfilename, std::ios::out);
   for(int i=0;i<Nvalidpm;i++){
      chNum=PMdata[i];
      ((TGraphErrors*)((*linearGarr)[i]))->Fit(Form("linear%d",chNum),"Q");
      Gain[i]=((TF1*)((*linFarr)[i]))->GetParameter(0)*pow(10,9);
      GainErr[i]=((TF1*)((*linFarr)[i]))->GetParError(0)*pow(10,9);
      NoiseVariance[i]=((TF1*)((*linFarr)[i]))->GetParameter(1);
      chgain->SetPoint(i,chNum,Gain[i]);
      chgain->SetPointError(i,0,GainErr[i]);
      gainhist->Fill(Gain[i]);
//      writing_file << ch << ","<< Gain[i]<<","<<GainErr[i]<<","<<NoiseVariance[i]<< endl;
      if(Gain[i]<0||GainErr[i]>Gain[i]){
         std::cout<<"bad id:    "<<i<<"  ch:    "<<chNum<<"  gain error:    "<<GainErr[i]<<std::endl;
      }/*else if(Gain[i]>5*pow(10,7)){
         std::cout<<"high gain:    "<<i<<"     ch    "<<ch<<std::endl;
         }*/
      onegain=Gain[i];
      onegainerr=GainErr[i];
      onenoise=NoiseVariance[i];
      tout->Fill();
   }
   fout->cd();
   tout->Write();
   fout->Close();
   canvas1->cd();
   gainhist->Draw();
   for(int i=0;i<Nmppc;i++){
      GainAllSiPM[i]=ExtendAllSiPM(Gain,PMdata,i);
      NoiseAllSiPM[i]=ExtendAllSiPM(NoiseVariance,PMdata,i);
   }
   canvas3->cd();
   InnerGeometry(GainAllSiPM,GainMax);
   canvas4->cd();
   InnerGeometry(NoiseAllSiPM,1);
   canvas2->cd();
   chgain->SetMarkerStyle(22);
   chgain->SetMarkerColor(2);
   chgain->SetMarkerSize(1);
   chgain->SetMaximum(GainMax);
   chgain->SetMinimum(0);
   chgain->SetTitle("Gain of each channel;channel;Gain");
   chgain->Draw("ap");
   canvas2->Print(Form("run%d-%d_chgain.pdf",runarray[0],runarray[runnum-1]));

   //gain histogram
   //gainhist->Write();

   /*char gainhistpng[64];
     sprintf(gainhistpng,"gainhist_%d-%d_%d_ultra.png",runarray[0],runarray[runnum-1],rangemode);
     canvas2->Print(gainhistpng);*/

}


void Intensity_scan(int ch) {
  getPMdata(runarray[0]);
  if(arraysearch(PMdata,ch)==-1){
     std::cerr<<"Channel not found"<<std::endl;
  }
  TGraphErrors *gr =new TGraphErrors(runnum);

  gr->SetMarkerStyle(22);
  gr->SetMarkerColor(2);
  gr->SetMarkerSize(1);
  std::cout<< "channel:    "<<ch<<std::endl;

  TCanvas *canvas1 = new TCanvas("canvas1","fit",600,400);
  TCanvas *canvas6 = new TCanvas("canvas6","chhist",600,400);
  TClonesArray *chargeHarr = new TClonesArray("TH1D",runnum);
  TClonesArray *fitFarr = new TClonesArray("TF1",runnum);
  Double_t MeanMax = 0;
  Double_t VarMax=0;
  for(int run = 0; run < runnum; run++) {
    new((*fitFarr)[run]) TF1(Form("gausvar%d",run),vargaus,-10,0.5,3);
    if (ch<Nmppc) {
      rangemin=rangemax4SiPM-range4SiPM;
      rangemax=rangemax4SiPM;
    }else{
      rangemin=rangemax4PMT-range4PMT;
      rangemax=rangemax4PMT;
    }
    new((*chargeHarr)[run]) TH1D(Form("charge histogram %d",ch),Form("charge histogram title %d;Charge[10^{9}e];Events",ch),100,rangemin,rangemax);
    TFile *frec = new TFile(Form("$(MEG2SYS)/analyzer/rec%06d.root", runarray[run]),"READ");

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
    if(run==0){
      ((TH1D*)((*chargeHarr)[run]))->Fit(Form("gausvar%d",run),"NQ");
      ((TH1D*)((*chargeHarr)[run]))->Draw();
      ((TF1*)((*fitFarr)[run]))->Draw("same");

    }else{
      ((TH1D*)((*chargeHarr)[run]))->Fit(Form("gausvar%d",run),"NQ");
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

  canvas1->cd();
  gr->SetTitle("Mean vs Variance;Mean[10^{9}e];Variance");
  gr->GetXaxis()->SetLimits(0,MeanMax*1.2);
  gr->SetMaximum(VarMax*1.2);
  gr->SetMinimum(0);
  TF1 *linF= new TF1("linearfit","[0]*x+[1]",0,5);
  gr->Fit("linearfit","Q");
  Double_t Gain=linF->GetParameter(0);
  Double_t GainErr=linF->GetParError(0);
  Double_t NoiseVariance=linF->GetParameter(1);
  std::cout<<"Gain:    "<<Gain<<std::endl;
  std::cout<<"Gain Error:    "<<GainErr<<std::endl;
  std::cout<<"Noise Variance:    "<<NoiseVariance<<std::endl;

  gr->Draw("ap");
}
