#define Npm 4760
#define Nmppc 4092
#define NRow 93
#define NLine 44

std::vector<Int_t> PMdata;
Int_t Nvalidpm;
Int_t Nvalidmppc=0;

Double_t GainAllSiPM[Nmppc];
Double_t NoiseAllSiPM[Nmppc];

void getPMdata(Int_t run) {
  TFile *frec = new TFile(Form("$(MEG2SYS)/analyzer/rec%06d.root", run),"READ");
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

Int_t arraysearch(std::vector<Double_t> array, Double_t value){
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

int colordtm(Double_t value,Double_t vmin, Double_t vmax){
  if (value<vmin) {
    value=vmin;
  }
  if (value>vmax) {
    value=vmax;
  }
  double ratio=(value-vmin)/(vmax-vmin);
  double crange=49;
  double cmin=51;
  double dbcolor=cmin+crange*ratio;
  int color=floor(dbcolor);
  return color;
}

TString safeName(TString name){
  TObject* old = gROOT->FindObject(name);
  if(old) delete old;
  return name;
}

Double_t vargaus(Double_t *x,Double_t *par){
  Float_t xx=x[0];
  Double_t f = par[0]/sqrt(2*TMath::Pi()*par[1])*TMath::Exp(-pow(xx-par[2],2)/(2*par[1]));
  return f;
}

Double_t ExtendAllSiPM(std::vector<Double_t> Base,std::vector<Int_t> PMdata,int i){
   Double_t value=0;
   if(arraysearch(PMdata,i)!=-1){
      value=Base[arraysearch(PMdata,i)];
   }
   return value;
}

void InnerGeometry(Double_t PropertyAllSiPM[Nmppc],Double_t GainMax){
   Double_t CenterBoxWidth=0.8;
   Double_t CenterBoxHeight=0.8;
   Double_t CenterBoxOriginX=(1+CenterBoxWidth)/2;
   Double_t CenterBoxOriginY=(1+CenterBoxHeight)/2;
   Double_t SiPMBoxWidth=CenterBoxWidth/NLine;
   Double_t SiPMBoxHeight=CenterBoxHeight/NRow;
   for(int i = 0; i<NRow; i++){
      for(int j = 0; j<NLine; j++){
         Double_t OneBoxXmax=CenterBoxOriginX-j*SiPMBoxWidth;
         Double_t OneBoxYmax=CenterBoxOriginY-i*SiPMBoxHeight;
         Double_t OneBoxXmin=OneBoxXmax-SiPMBoxWidth;
         Double_t OneBoxYmin=OneBoxYmax-SiPMBoxHeight;
         TBox *b= new TBox(OneBoxXmin,OneBoxYmin,OneBoxXmax,OneBoxYmax);
         TString chnum;
         chnum.Form("%d",i*NRow+j);
         TText *t=new TText(OneBoxXmin,OneBoxYmin,chnum);
         b->SetFillColor(colordtm(PropertyAllSiPM[i*NLine+j],0,GainMax));
         b->Draw();
         //t->Draw();
      }
   }
}
