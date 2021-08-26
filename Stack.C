#include "RooTFnBinding.h"

using namespace RooFit ;

double GetLegacyResBarrel(const int& year, double p){

  // This is only for Muon Barrel

  double res = -1.;

  if (year == 2016) {
    res = 0.006152416185063
      + 0.000100116444003*p
      - 0.000000103718255*pow(p,2)
      + 0.000000000056979*pow(p,3)
      - 0.000000000000011*pow(p,4);
  } else if (year == 2017) {
    res = 0.005265784542684
      + 0.000114819686533*p
      - 0.000000130882947*pow(p,2)
      + 0.000000000069409*pow(p,3)
      - 0.000000000000013*pow(p,4);
  } else if (year == 2018) {
    res = 0.006204084326376
      + 0.000095866688875*p
      - 0.000000097138211*pow(p,2)
      + 0.000000000048821*pow(p,3)
      - 0.000000000000009*pow(p,4);
  }

  return res;

}

TGraph* GetLegacyPlot(const int& year, TGraph* gNew){
  const int np = gNew->GetN();
  TGraph* g = new TGraph(np);
  double p,res;
  for (int i = 0; i < np; ++i) {
    gNew->GetPoint(i,p,res);
    double_t resLegacy = GetLegacyResBarrel(year,p);
    std::cout << year << "\t" << p << "\t" << res << "\t" << resLegacy << "\n" ;
    if (resLegacy < 0.) throw;
    g->SetPoint(i,p,resLegacy);
  }
  return g;
}

double DSCB(double *x, double *par){

  double alpha_l = par[0];
  double alpha_h = par[1];
  double n_l     = par[2];
  double n_h     = par[3];
  double mean    = par[4];
  double sigma   = par[5];
  double N       = par[6];
  double t = (x[0]-mean)/sigma;
  double result;

  double fAlphaL = alpha_l/n_l;
  double f2AlphaL = (n_l/alpha_l) - alpha_l - t;
  double fAlphaH = alpha_h/n_h;
  double f2AlphaH = (n_h/alpha_h) - alpha_h + t;

  if (-alpha_l <= t && alpha_h >= t){
    result = exp(-0.5*t*t);
  } else if (t < -alpha_l) {
    result = exp(-0.5*alpha_l*alpha_l)*pow(fAlphaL*f2AlphaL, -n_l);
  } else if (t > alpha_h) {
    result = exp(-0.5*alpha_h*alpha_h)*pow(fAlphaH*f2AlphaH, -n_h);
  }

  return N*result;

}

TGraphAsymmErrors* GetResolutionGraph(const int& year, const int& etaBins_, bool mergeTracker = false) {

  std::cout << "=============\t" << year << "\t=============\t" << etaBins_ << "\t=============\n";

  TFile* f1 = TFile::Open("MuonResolution_SameBins.root","READ");

  std::vector<std::string> etaBins = {
    "HPResB_G", "HPResE_G",
    "HPResB_T", "HPResE_T",
  };

  std::unordered_map<int, float> luminosity = {
    {2016, 35.92},
    {2017, 41.43},
    {2018, 59.74}
  };

  // std::vector<std::string> etaBins = {
  //   "HPtResB_G", "HPtResE_G",
  //   "HPtResB_T", "HPtResE_T",
  // };

  std::pair<float,float> MuonB = { 0., 1.2 };
  std::pair<float,float> MuonE = { 1.2, 2.4 };

  std::vector<std::string> titleEtaBins = {
    Form("%.1f <=|#eta|<= %.1f [globalHighPtId]", MuonB.first, MuonB.second),  //B_G
    Form("%.1f <|#eta|<= %.1f [globalHighPtId]", MuonE.first, MuonE.second),   //E_G
    Form("%.1f <=|#eta|<= %.1f [trackerHighPtId]", MuonB.first, MuonB.second), //B_T
    Form("%.1f <|#eta|<= %.1f [trackerHighPtId]", MuonE.first, MuonE.second),  //O_T
    Form("%.1f <=|#eta|<= %.1f [global + trackerHighPtId]", MuonB.first, MuonB.second), //B_G
    Form("%.1f <|#eta|<= %.1f [global + trackerHighPtId]", MuonE.first, MuonE.second),  //E_G
  };

  std::unordered_map<int,std::vector<std::pair<std::string,Double_t>>> samples =
    {
      {
        2016,
        {
          {"ZToMuMu_NNPDF30_13TeV-powheg_M_50_120", 2.116e+03},
          {"ZToMuMu_NNPDF30_13TeV-powheg_M_120_200", 2.058e+01},
          {"ZToMuMu_NNPDF30_13TeV-powheg_M_200_400", 2.890e+00},
          {"ZToMuMu_NNPDF30_13TeV-powheg_M_400_800", 2.515e-01},
          {"ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400", 1.709e-02},
          {"ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300", 1.370e-03},
          {"ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500", 8.282e-05},
          {"ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500", 3.414e-06},
          {"ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000", 3.650e-07},
          {"ZToMuMu_NNPDF30_13TeV-powheg_M_6000_Inf", 2.526e-08},
        },
      },
      {
        2017,
        {
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_50_120", 2.116e+03},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_120_200", 2.058e+01},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_200_400", 2.890e+00},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_400_800", 2.515e-01},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_800_1400", 1.709e-02},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_1400_2300", 1.370e-03},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_2300_3500", 8.282e-05},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_3500_4500", 3.414e-06},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_4500_6000", 3.650e-07},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_6000_Inf", 2.526e-08}
        }
      },
      {
        2018,
        {
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_50_120", 2.116e+03},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_120_200", 2.058e+01},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_200_400", 2.890e+00},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_400_800", 2.515e-01},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_800_1400", 1.709e-02},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_1400_2300", 1.370e-03},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_2300_3500", 8.282e-05},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_3500_4500", 3.414e-06},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_4500_6000", 3.650e-07},
          { "ZToMuMu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_M_6000_Inf", 2.526e-08},
        }
      }
    };

  auto rhpath = [&](const int& i){
    std::string hpath = Form("%d/%s/",year,samples[year][i].first.c_str());
    std::cout << hpath << std::endl;
    return hpath;
  };


  TCanvas* cps = new TCanvas("cps","cps",5*500,3*500);
  cps->Divide(5,3);

  TCanvas* cpt = new TCanvas("cpt","cpt",5*500,3*500);
  cpt->Divide(5,3);

  Int_t nParams = 7;

  Double_t xmin = -0.3;
  Double_t xmax = 0.3;

  std::vector<Double_t> sigmas_;
  std::vector<Double_t> ptBins_;
  std::vector<Double_t> sigmaErrors_;

  Double_t prevSigma_ = 0.01; // Bottom limit 1%

  auto GetHisto = [&] (const int& nSample, const int& nPbin, const int& nEtaBin){
    auto h3 =
      static_cast<TH3D*>(f1->Get(Form("%d/%s/%s",year,samples[year][nSample].first.c_str(),etaBins[nEtaBin].c_str())));
    auto h2 =
      static_cast<TH2D*>(h3->Project3D("zx"));
    auto h1 =
      static_cast<TH1D*>(h2->ProjectionY(Form("%s_h_%d",samples[year][nSample].first.c_str(),nPbin),nPbin)->Clone());
    TH1F* hCutFlow = static_cast<TH1F*>(f1->Get(Form("%s/HCutFlow",rhpath(nSample).c_str())));
    Double_t sumGenWeight = hCutFlow->GetBinContent(hCutFlow->GetXaxis()->FindBin("genWeight"));
    h1->Scale(luminosity[year]*1e3*samples[year][nSample].second/sumGenWeight);

    float xmin = h2->GetXaxis()->GetBinLowEdge(nPbin);
    float xmax = h2->GetXaxis()->GetBinLowEdge(nPbin+1);

    return std::make_pair(h1,std::make_pair(xmin,xmax));
  };



  auto GetMergedHisto = [&] (const int& initialBin, const int& endBin){
    THStack* hsMerged = new THStack(Form("hsMerged_%d_%d_%d",etaBins_,initialBin,endBin),"hsMerged");
    TH1D* hMerged;
    Double_t pBinMin, pBinMax;
    for(int j = initialBin; j <= endBin /*P Bins*/; ++j) {
      for (int i = 0; i < samples[year].size(); ++i) {
        auto h3 =
          static_cast<TH3D*>(f1->Get(Form("%d/%s/%s",year,samples[year][i].first.c_str(),etaBins[etaBins_].c_str())));
        auto h2 =
          static_cast<TH2D*>(h3->Project3D("zx"));
        pBinMin = h2->GetXaxis()->GetBinLowEdge(initialBin);
        pBinMax = h2->GetXaxis()->GetBinLowEdge(endBin+1);
        auto h1 =
          static_cast<TH1D*>(h2->ProjectionY(Form("%s_h_%d",samples[year][i].first.c_str(),j),j));

        TH1F* hCutFlow = static_cast<TH1F*>(f1->Get(Form("%s/HCutFlow",rhpath(i).c_str())));
        Double_t sumGenWeight = hCutFlow->GetBinContent(hCutFlow->GetXaxis()->FindBin("genWeight"));
        if( j == initialBin ) {
          hMerged = static_cast<TH1D*>(h1->Clone());
        } else {
          if (h1->GetEntries()>1e2)
            hMerged->Add(h1,luminosity[year]*1e3*samples[year][i].second/sumGenWeight);
        }
      }
    }
    hsMerged->SetTitle(Form("[%.1f:%.1f] GeV %s [%d];(1/p-1/p^{GEN})/(1/p^{GEN});Event Count",pBinMin,pBinMax,titleEtaBins[etaBins_].c_str(),year));
    hsMerged->Add(hMerged,"F");
    return hsMerged;
  };

  auto GetModifiedXSec = [] (const double& xnorm, const double& xsec) {
    const double factor = 1/xnorm;
    return factor*xsec;
  };

  auto GetHistoLimits = [] (const TH1* h) {

    float lowerLimit = 0;
    float upperLimit = 0;
    for (int i = h->GetNbinsX(); i != 0; --i) {
      if ( h->GetBinContent(i) != 0. )
        lowerLimit = h->GetBinContent(i);
    }
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      if ( h->GetBinContent(i) != 0. )
        upperLimit = h->GetBinContent(i);
    }

    return std::make_pair(lowerLimit, upperLimit);

  };


  int nPBins = static_cast<TH3D*>
    (f1->Get(Form("%d/%s/%s",year,samples[year][0].first.c_str(),
                  etaBins[etaBins_].c_str())))->GetNbinsX();

  for(int j = 1; j <= nPBins; ++j) {

    Double_t yMax = 0.;
    THStack* hs = new THStack("hs",Form("%s;P Resolution;Event Count",etaBins[etaBins_].c_str()));
    THStack* hs_ = static_cast<THStack*>(hs->Clone());
    hs->SetName(Form("hs_%s_%d",etaBins[etaBins_].c_str(),j));
    hs_->SetName(Form("hs__%s_%d",etaBins[etaBins_].c_str(),j));
    Double_t ptBinMin = 0;
    Double_t ptBinMax = 0.;

    for(int i = 0; i < samples[year].size(); ++i){

      std::cout << Form("%d/%s/%s\n",year,samples[year][i].first.c_str(),etaBins[etaBins_].c_str()) ;

      std::pair<TH1D*,std::pair<float,float>> hpack = GetHisto(i,j,etaBins_);

      auto h1 = hpack.first ;

      ptBinMin = hpack.second.first;
      ptBinMax = hpack.second.second;

      std::string sTitle = Form("[%.1f:%.1f] GeV %s [%d];(1/p-1/p^{GEN})/(1/p^{GEN});Event Count",ptBinMin,ptBinMax,titleEtaBins[etaBins_].c_str(),year);

      if (mergeTracker and etaBins_ < 2) {
        std::pair<TH1D*,std::pair<float,float>> hpackTracker = GetHisto(i,j,etaBins_ + 2);
        auto hh = hpackTracker.first;
        bool success = h1->Add(hh);
        if(!success) throw;
        sTitle = Form("[%.1f:%.1f] GeV %s [%d];(1/p-1/p^{GEN})/(1/p^{GEN});Event Count",ptBinMin,ptBinMax,titleEtaBins[etaBins_ + 4].c_str(),year);
      }


      hs->SetTitle(sTitle.c_str());
      hs_->SetTitle(sTitle.c_str());
      h1->SetTitle(Form("%s;P Resolution;EventCount",samples[year][i].first.c_str()));
      //h1->GetYaxis()->SetRangeUser(0.,yMax);
      h1->SetFillColor(i+1);
      cps->cd(i+1);
      h1->Draw("HIST");
      hs->Add(h1,"F");
      if( h1->GetEntries() > 1e2){
        hs_->Add(h1,"F");
      }
    }
    ptBins_.emplace_back(ptBinMin);

    if( !mergeTracker and etaBins_ > 1 and ptBinMin > 700) {
      const int nPBinsTmp = nPBins;
      hs_ = GetMergedHisto(j,nPBins);
      nPBins = j;
      j = nPBinsTmp+1;
      ptBinMax = 3000;
    }


    cps->cd(11);
    hs->Draw("HIST");
    cps->cd(12);
    hs_->Draw("HIST");
    TH1* h = static_cast<TH1D*>(hs_->GetStack()->Last());
    //h->SetTitle(Form("[%.1f:%.1f] GeV %s [%d];(1/p-1/p^{GEN})/(1/p^{GEN});Event Count",ptBinMin,ptBinMax,titleEtaBins[etaBins_].c_str(),year));

    if (etaBins_ == 0 and ptBinMin == 800.) {
      xmin = -0.2;
      xmax = 0.2;
    } else if (etaBins_ == 0 and ptBinMin == 1.2e3) {
      xmin = -0.2;
      xmax = 0.2;
    } else if (etaBins_ == 0 and ptBinMin == 72.) {
      xmin = -0.1;
      xmax = 0.1;
    } else if (etaBins_ == 3 and ptBinMin < 73.) {
      xmin = -0.2;
      xmax = 0.2;
    } else if (year == 2018 and etaBins_ == 3 and ptBinMin == 800) {
      xmin = -0.35;
      xmax = 0.35;
    } else {
      xmin = -0.3;
      xmax = 0.3;
    }



    TF1 *fxDCB = new TF1(Form("fxDCB_%d_%s_%.0f",year,etaBins[etaBins_].c_str(),ptBinMin),
                         DSCB, xmin, xmax, nParams);
    fxDCB->SetParNames("#alpha_{low}","#alpha_{high}","n_{low}", "n_{high}", "#mu", "#sigma", "N");
    fxDCB->SetParameters(1., 1., 10, 10, h->GetMean(), h->GetRMS(), h->GetMaximum());
    fxDCB->SetParLimits(5,prevSigma_,0.2); // Require higher sigma for higher Pt


    h->SetName(Form("%d_%s_%d",year,etaBins[etaBins_].c_str(),j));

    std::string fitOption = "QMWL";

    h->Fit(fxDCB,fitOption.c_str(),"",xmin,xmax);

    Double_t errorTmp = fxDCB->GetParError(5);
    Double_t sigmaTmp = fxDCB->GetParameter(5);

    int counter = 1;
    const double goodError = 0.011;

    auto flipFitOption  = [&] (){
      if(fitOption.compare("QMWL") == 0) {
        fitOption = "QME";
      } else if(fitOption.compare("QME") == 0) {
        fitOption = "QMWL";
      }
      std::cout << "Fliping fitOption:\n\t" << fitOption << "\n";
      counter = 1;
    };

    auto Fit = [&] () {
      std::cout << counter << "\n";
      if( counter > 10 ) flipFitOption();
      h->Fit(fxDCB,fitOption.c_str(),"",xmin,xmax);
      errorTmp = fxDCB->GetParError(5);
      sigmaTmp = fxDCB->GetParameter(5);
      ++counter;
    };

    while( errorTmp > goodError) {
      Fit();
    }

    if( isnan(errorTmp) ) {
      Fit();
      counter = 1;
      while( errorTmp > goodError){
        Fit();
      }
    }

    prevSigma_ = fxDCB->GetParameter(5);
    sigmas_.emplace_back(fxDCB->GetParameter(5));
    sigmaErrors_.emplace_back(fxDCB->GetParError(5));
    fxDCB->Draw("SAME");

    RooRealVar pres("pres","P Residual",xmin,xmax);
    pres.setBins(10000);
    RooAbsPdf* dcb = RooFit::bindPdf(fxDCB,pres);
    RooDataHist dh1("dh1","dh1",pres,h);
    std::string sTitle = Form("[%.1f:%.1f] GeV %s [%d];(1/p-1/p^{GEN})/(1/p^{GEN}) [#sigma %.2f];Event Count",ptBinMin,ptBinMax,titleEtaBins[etaBins_].c_str(),year,fxDCB->GetParameter(5));
    if(mergeTracker)
      sTitle = Form("[%.1f:%.1f] GeV %s [%d];(1/p-1/p^{GEN})/(1/p^{GEN}) [#sigma %.2f];Event Count",ptBinMin,ptBinMax,titleEtaBins[etaBins_+4].c_str(),year,fxDCB->GetParameter(5));
    RooPlot* frame = pres.frame(Title(sTitle.c_str()));
    frame->SetTitle(sTitle.c_str());

    dcb->fitTo(dh1,SumW2Error(true));
    dh1.plotOn(frame);
    dcb->plotOn(frame);

    cpt->cd(j);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    //hs_->Draw("HIST");
    //h->Draw();
    //h->GetXaxis()->SetRangeUser(-0.4,0.4);
    //fxDCB->Draw("SAME");
    frame->Draw();
    std::string filename = Form("Bin_%d_%s_%.0f_.png",year,etaBins[etaBins_].c_str(),ptBinMin);
    if(mergeTracker)
      filename = Form("Bin_%d_%s_MergedIDs_%.0f_.png",year,etaBins[etaBins_].c_str(),ptBinMin);
    cps->Print(filename.c_str());
  }

  ptBins_.emplace_back(3000); // Last limit


  //if (etaBins_ > 1) nPoints_ = 11;
  TGraphAsymmErrors* g_ = new TGraphAsymmErrors(nPBins);

  for(int i = 0; i < nPBins; ++i){
    Double_t mid = (ptBins_[i]+ptBins_[i+1])/2.;
    Double_t dx = mid - ptBins_[i];
    g_->SetPoint(i,mid,sigmas_[i]);
    g_->SetPointError(i,dx,dx,sigmaErrors_[i],sigmaErrors_[i]);
    //    g_->SetPointError(i,dx,dx,0,0);
  }

  std::string filename = Form("%d_%s_.png",year,etaBins[etaBins_].c_str());
  if(mergeTracker)
    filename = Form("%d_%s_Merged_.png",year,etaBins[etaBins_].c_str());
  cpt->Print(filename.c_str());

  f1->Close();

  return g_;

}

int Stack() {

  std::vector<int> etaBins = { 0, 1, 2, 3};

  TFile* fOut = TFile::Open("SummaryPlots.root","UPDATE");

  TCanvas* c = new TCanvas("c","c",3*500,3*500);
  c->Divide(3,3);

  TCanvas* cRatio = new TCanvas("cRatio","cRatio",3*500,1*500);
  cRatio->Divide(3,1);

  std::vector<int>  Years = { 2016, 2017, 2018 };

  std::pair<float,float> MuonB = { 0., 1.2 };
  std::pair<float,float> MuonE = { 1.2, 2.4 };

  std::vector<std::string> titleEtaBins = {
    Form("%.1f <=|#eta|<= %.1f [globalHighPtId]; P [GeV]; Resolution [%%]", MuonB.first, MuonB.second),
    Form("%.1f <|#eta|<= %.1f [globalHighPtId]; P [GeV]; Resolution [%%]", MuonE.first, MuonE.second),
    Form("%.1f <=|#eta|<= %.1f [trackerHighPtId]; P [GeV]; Resolution [%%]", MuonB.first, MuonB.second),
    Form("%.1f <|#eta|<= %.1f [trackerHighPtId]; P [GeV]; Resolution [%%]", MuonE.first, MuonE.second),
    Form("%.1f <=|#eta|<= %.1f [global+trackerHighPtId]; P [GeV]; Resolution [%%]", MuonB.first, MuonB.second),
    Form("%.1f <|#eta|<= %.1f [global+trackerHighPtId]; P [GeV]; Resolution [%%]", MuonE.first, MuonE.second),
  };

  std::vector<Int_t> colorEtaBins = {
    kRed, kBlack, kRed, kBlack
  };

  std::vector<std::string> legendEtaBins = {
    Form("%.1f <=|#eta|<= %.1f", MuonB.first, MuonB.second),
    Form("%.1f <|#eta|<= %.1f", MuonE.first, MuonE.second),
    Form("%.1f <=|#eta|<= %.1f", MuonB.first, MuonB.second),
    Form("%.1f <|#eta|<= %.1f", MuonE.first, MuonE.second),
  };

  auto FitToP4 = [&] (TGraph* gr, const int yr, const int etaBins_){
    gr->SetName(Form("MResolution_%d_%d",yr,etaBins_));
    std::string fxName = Form("p4_%d_%d",yr,etaBins_);
    TF1 *p4 = new TF1(fxName.c_str(),"pol4",53,2500);
    p4->SetParameters(0.016995688776817,
                      0.000083569126602,
                      0.000000002611754,
                      -0.000000000023200,
                      0.000000000000008);
    p4->SetParLimits(0,-0.1,0.1);
    p4->SetParLimits(1,-0.1,0.1);
    p4->SetParLimits(2,-0.1,0.1);
    p4->SetParLimits(3,-0.1,0.1);
    p4->SetParLimits(4,-0.1,0.1);

    for(int i = 0; i<8; ++i){
      gr->Fit(fxName.c_str(),"MWRE");
    }

    fOut->cd();
    p4->Write();
    gr->Write();
    p4->Print();

    return p4;
  };

  for(auto yr: Years){

    TMultiGraph *mgG = new TMultiGraph();
    mgG->SetName(Form("mgG_%d",yr));
    mgG->SetTitle(Form("P Resolution [%d] [globalHighPt]; P Reco; P Resolution [#sigma]",yr));

    TMultiGraph *mgT = new TMultiGraph();
    mgT->SetName(Form("mgT_%d",yr));
    mgT->SetTitle(Form("P Resolution [%d] [trackerHighPt]; P Reco; P Resolution [#sigma]",yr));

    TMultiGraph *mgGT = new TMultiGraph();
    mgGT->SetName(Form("mgGT_%d",yr));
    mgGT->SetTitle(Form("P Resolution [%d] [global+trackerHighPt]; P Reco; P Resolution [#sigma]",yr));

    TLegend *lG = new TLegend(0.7,0.2,0.95,0.3);
    lG->SetName(Form("lG_%d",yr));

    TLegend *lT = new TLegend(0.7,0.2,0.95,0.3);
    lG->SetName(Form("lT_%d",yr));

    TLegend *lGT = new TLegend(0.7,0.2,0.95,0.3);
    lGT->SetName(Form("lGT_%d",yr));

    TGraph* gRatioLegacy;
    TGraph* gNew;
    TGraph* gMergedId;

    for(auto etaBins_: etaBins){
      TGraphAsymmErrors *g = GetResolutionGraph(yr,etaBins_);
      FitToP4(g,yr,etaBins_);
      g->SetLineColor(colorEtaBins[etaBins_]);
      g->SetMarkerColor(colorEtaBins[etaBins_]);
      g->SetMarkerStyle(23);
      if(etaBins_ < 2){
        TGraphAsymmErrors *gg = GetResolutionGraph(yr,etaBins_,true);
	FitToP4(gg,yr,etaBins_);
        gg->SetLineColor(colorEtaBins[etaBins_]);
        gg->SetMarkerColor(colorEtaBins[etaBins_]);
        gg->SetMarkerStyle(23);
        lGT->AddEntry(gg,legendEtaBins[etaBins_].c_str());
        mgGT->Add(gg,"P");
        if (etaBins_ == 0) gNew = g;
        lG->AddEntry(g,legendEtaBins[etaBins_].c_str());
        mgG->Add(g,"P");
      } else {
        mgT->Add(g,"P");
        lT->AddEntry(g,legendEtaBins[etaBins_].c_str());
      }
    }

    // Print Ratio
    TMultiGraph *mgR = new TMultiGraph(Form("mgR_%d",yr),
       Form("UL[Red] vs Legacy[Black] %d;P[GeV]; Resolution [%%]",yr));
    mgR->Add(gNew,"P");
    std::cout << "nRatio: " << gNew->GetN() << "\n\n\n";
    gRatioLegacy = GetLegacyPlot(yr,gNew);
    gRatioLegacy->SetMarkerStyle(23);
    gRatioLegacy->SetMarkerColor(kBlack);
    std::cout << "nLegacy: " << gRatioLegacy->GetN() << "\n\n\n";
    gRatioLegacy->Print();
    mgR->Add(gRatioLegacy,"P");
    cRatio->cd(yr%2015);
    mgR->Draw("AP");
    mgR->GetYaxis()->SetRangeUser(0.,0.18);
    cRatio->Print("ComparissonToLegacy.png");

    c->cd(yr%2015);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    mgG->Draw("AP");
    lG->Draw();
    mgG->GetXaxis()->SetRangeUser(0,3000);
    mgG->GetYaxis()->SetRangeUser(0.,0.18);

    c->cd((yr%2015)+3);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    mgT->Draw("AP");
    lT->Draw();
    mgT->GetXaxis()->SetRangeUser(0,3000);
    mgT->GetYaxis()->SetRangeUser(0.,0.18);

    c->cd((yr%2015)+6);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    mgGT->Draw("AP");
    lGT->Draw();
    mgGT->GetXaxis()->SetRangeUser(0,3000);
    mgGT->GetYaxis()->SetRangeUser(0.,0.18);
  }

  c->Print(Form("ResolutionMeasurement.png"));

  return 0;
}
