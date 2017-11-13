TString makefoutname(int runini, int runend,int rangemode){
  TString foutname;
  foutname.Form("$(MEG2SYS)/analyzer/macros/xec/sLED/gaindata/gain_%d-%d_%d.root",runini,runend,rangemode);
  return foutname;
}

void PMTcor(){
  TCanvas* canvas1 =new TCanvas("canvas1","Gain Difference",600,600);
  TCanvas* canvas2 =new TCanvas("canvas2","Gain Comparison",600,600);

  //Int_t runone[2]={306902,306908};
  //Int_t runone[2]={306909,306914};
  //Int_t runtwo[2]={306922,306928};
  //Int_t runtwo[2]={306701,306705};
  //Int_t runone[2]={306998,307004};
  //Int_t runtwo[2]={308540,308547};
  Int_t runone[2]={308540,308547};
  Int_t runtwo[2]={308583,308590};

  Int_t rangemode=0;
  TFile* fin[2];
  fin[0]= new TFile(makefoutname(runone[0],runone[1],rangemode),"read");
  fin[1]= new TFile(makefoutname(runtwo[0],runtwo[1],rangemode),"read");

  Int_t chNum[2];
  Double_t gain[2];
  Double_t noise[2];
  Double_t gainerr[2];
  TTree* tin[2];
  Int_t Nevent[2];
  for(int i=0;i<2;i++){
    tin[i] = (TTree*)fin[i]->Get("tree");
    tin[i]->SetBranchAddress("chNum", &chNum[i]);
    tin[i]->SetBranchAddress("gain", &gain[i]);
    tin[i]->SetBranchAddress("noise", &noise[i]);
    tin[i]->SetBranchAddress("gainerr", &gainerr[i]);
    Nevent[i]=tin[i]->GetEntries();
  }
  TString histtitle;
  histtitle.Form("gain difference: run%d-%d vs run %d-%d",runone[0],runone[1],runtwo[0],runtwo[1]);
  histtitle=histtitle+";Gain_{new}-Gain_{old};Channels";
  TH1D* gaindiffhist = new TH1D("Gain_{new}-Gain_{old}",histtitle,200,-2.5*pow(10,7),2.5*pow(10,7));
  TF1* fitgaus = new TF1("f1","[0]*TMath::Gaus(x,[1],[2])");
  TF1* corfit = new TF1("f2","[0]*x+[1]");
  fitgaus->SetParameters(100,0,8*pow(10,5));
  TGraphErrors* GainComp[2];
  TGraphErrors* GainCor=new TGraphErrors();

  for(int i=0;i<2;i++){
    GainComp[i]  = new TGraphErrors();
    GainComp[i]->SetMarkerStyle(22);
  }
  Int_t ich[2]={};
  Int_t idg=0;
  while(ich[0]<Nevent[0]&&ich[1]<Nevent[1]){
    tin[0]->GetEntry(ich[0]);
    tin[1]->GetEntry(ich[1]);

    if(chNum[0]==chNum[1]){
      if(chNum[0]>4092){

        Bool_t quality=true;
        for (int run = 0; run < 2; run++) {
          if(gainerr[run]>0.2*gain[run]||gain[run]<0||gain[run]>pow(10,7)){
            quality=false;
            std::cout<<"bad channel: "<<chNum[0]<<std::endl;
            break;
          }
        }
        if(quality==true){
          for(int i=0;i<2;i++){
            GainComp[i]->SetPoint(idg,chNum[i],gain[i]);
            GainComp[i]->SetPointError(idg,0,gainerr[i]);
          }
          Int_t id_cor=GainCor->GetN();
          GainCor->SetPoint(id_cor,gain[0],gain[1]);
          GainCor->SetPointError(id_cor,gainerr[0],gainerr[1]);
          idg++;
          Double_t GainDiff=gain[1]-gain[0];
          gaindiffhist->Fill(GainDiff);
        }
      }
      ich[0]++;
      ich[1]++;
    }else if(chNum[0]<chNum[1]){
      ich[0]++;
    }else{
      ich[1]++;
    }
  }
  canvas1->cd();
  fitgaus->SetLineColor(2);
  gaindiffhist->Fit("f1");
  gaindiffhist->Draw();

  canvas2->cd();
  GainCor->SetTitle("Gain Correlation;Gain_{latest};Gain_{MEG I}");
  GainCor->GetXaxis()->SetLimits(0,15*pow(10,6));
  GainCor->SetMaximum(15*pow(10,6));
  GainCor->SetMinimum(0);
  GainCor->Fit("f2");
  GainCor->Draw("ap");
  canvas2->Print(Form("Gain_correlation_PMT_run%d-%dvsrun%d-%d.png",runone[0],runone[1],runtwo[0],runtwo[1]));
  // GainComp[0]->SetMarkerColor(2);
  // GainComp[0]->Draw("ap");
  // GainComp[1]->Draw("p");

}
