#include "RooTFnBinding.h"
#include "RooCrystalBall.cxx"

using namespace RooFit ;

TGraphAsymmErrors* plotFits(Int_t year, std::string hname, Bool_t isData = false){

  std::cout << Form("\n\n=============== %d %s ==============\n\n", year, hname.c_str());

  TFile *f1 = TFile::Open("MuonStudies_ZPeakResolution.root");

  int nBins = 7;

  TCanvas* c1 = new TCanvas("c1","c1",2000,1000);
  c1->Divide(4,2);

  TCanvas* cPull = new TCanvas("cPull","cPull",4*500,2*500);
  cPull->Divide(4,2);

  std::unordered_map<int, float> luminosity = {
    {2016, 35.92},
    {2017, 41.43},
    {2018, 59.74}
  };

  std::unordered_map<int,std::vector<std::pair<std::string,Double_t>>> samples =
    {
      {
        2016,
        {
          {"DYJetsToMuMu_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM_pdfwgt_F-pythia8", 6.178e+00},
          // {"DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", 1.702e+02},
          // {"DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", 1.475e+02},
          // {"DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", 4.104e+01},
          // {"DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", 5.669e+00},
          // {"DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",1.360e+00},
          // {"DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",6.227e-01},
          // {"DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", 1.512e-01},
          // {"DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", 3.659e-03},
          // {"DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8", 5.128e+03},
          // {"DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8", 9.556e+02},
          // {"DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8", 3.600e+02},
        },
      },
      {
        2017,
        {
          {"DYJetsToMuMu_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM_pdfwgt_F-pythia8",6.182e+00},
          // {"DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 1.395e+02},
          // {"DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 1.404e+02},
          // {"DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 3.836e+01},
          // {"DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 5.220e+00},
          // {"DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 1.264e+00},
          // {"DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 5.685e-01},
          // {"DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 1.330e-01},
          // {"DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 2.981e-03},
          // {"DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8", 5.128e+03},
          // {"DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8", 9.556e+02},
          // {"DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8", 3.600e+02},
        }
      },
      {
        2018,
        {
          {"DYJetsToMuMu_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM_pdfwgt_F-pythia8",6.178e+00},
          // {"DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 1.398e+02},
          // {"DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 1.406e+02},
          // {"DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 3.840e+01},
          // {"DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 5.213e+00},
          // {"DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 1.265e+00},
          // {"DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 5.683e-01},
          // {"DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 1.330e-01},
          // {"DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", 2.987e-03},
          // {"DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8", 5.128e+03},
          // {"DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8", 9.556e+02},
          // {"DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8", 3.600e+02},
        }
      }
    };

  auto GetLumiFactor = [&](const int& year_, const std::string& sample) {
    TH1F* hCutFlow = static_cast<TH1F*>(f1->Get(Form("%d/%s/HCutFlow",year_,sample.c_str())));
    Double_t sumGenWeight = hCutFlow->GetBinContent(hCutFlow->GetXaxis()->FindBin("genWeight"));
    return luminosity[year_]*1e3/sumGenWeight;
  };

  std::string histopath = "";

  if(isData)
    histopath = Form("%d/ULSingleMuon/%s",year,hname.c_str());

  std::cout << histopath << "\n";

  std::vector<Float_t> sigmas;
  std::vector<Float_t> sigmaErrors;
  std::vector<Float_t> ptBins;


  TH2F* h2;

  if(isData)
    h2 = static_cast<TH2F*>(f1->Get(histopath.c_str())->Clone());

  for(auto sample: samples[year]){
    if(isData) break; // No need to stack up for data
    if(histopath.size()==0){
      histopath = Form("%d/%s/%s",year,sample.first.c_str(),hname.c_str());
      h2 = static_cast<TH2F*>(f1->Get(histopath.c_str())->Clone());
      h2->SetName(histopath.c_str());
      h2->Scale(sample.second*GetLumiFactor(year,sample.first));
    } else {
      h2->Add(static_cast<TH2F*>(f1->Get(histopath.c_str())),sample.second*GetLumiFactor(year,sample.first));
    }
    std::clog << Form("Adding sample:\t%s\n",sample.first.c_str());
  }

  Double_t YLimit = 0.;

  Bool_t breakLoop = false;

  TH1D *h;

  const float ptLimitLowStats_global = 400.;
  const float ptLimitLowStats_tracker = 140.;

  Double_t prevMean = 0.;
  Double_t prevSigma = 1.5;
  Double_t prevAlphaL = 2.;
  Double_t prevAlphaR = 2.;
  Double_t prevNL = 1.5;
  Double_t prevNR = 1.5;

  for(int k = 1; k <= nBins; ++k){

    if(breakLoop)
      break;

    c1->cd(k);

    Float_t ptBinLow = h2->GetXaxis()->GetBinLowEdge(k);
    ptBins.emplace_back(ptBinLow);
    Float_t ptBinHigh = h2->GetXaxis()->GetBinLowEdge(k+1);

    if( k == nBins )
      ptBins.emplace_back(ptBinHigh);


    if ( ((hname.find("HMassZPt_A_G") != std::string::npos)
          or
          (hname.find("HMassZPt_B_G") != std::string::npos))
         and ptBinLow > ptLimitLowStats_global ) {
      nBins = k;
      ptBins.emplace_back(ptBinHigh);
      break;
    } else if ( ((hname.find("HMassZPt_A_T") != std::string::npos)
                 or
                 (hname.find("HMassZPt_B_T") != std::string::npos)
                 or
                 (hname.find("HMassZPt_A_GT") != std::string::npos)
                 or
                 (hname.find("HMassZPt_B_GT") != std::string::npos)
                 )
                and ptBinLow > ptLimitLowStats_tracker) {
      nBins = k;
      ptBins.emplace_back(ptBinHigh);
      break;
    }

    h = static_cast<TH1D*>(h2->ProjectionY(Form("%s_%.0f",hname.c_str(),ptBinHigh),k));
    h->SetTitle(Form("%s [%.0f:%.0f];Dimuon Mass [GeV];Event Count", hname.c_str(), ptBinLow, ptBinHigh));

    if (k == 1)
      YLimit = h->GetMaximum()*1.1;

    //std::pair<float,float> MassWindow = {-150., 150.};

    RooRealVar *mass = new RooRealVar("mass","m_{Z} (GeV)",-100.,250.);
    mass->setBins(2000,"cache");

    RooRealVar *massZdpg = new RooRealVar("massZdpg","DPG Mass Z", 91.1855, 91.1897);
    RooRealVar *widthZdpg =  new RooRealVar("widthZdpg","DPG Width Z", 2.4929, 2.4975);
    RooBreitWigner *breitW = new RooBreitWigner("breitW","Fit PDF",*mass,*massZdpg,*widthZdpg);

    double zero = 1e-3;

    RooRealVar* dcbMean  = new RooRealVar("dcbMean","dcbMean",prevMean, -1.5, 1.5);
    RooRealVar* dcbSigma = new RooRealVar("dcbSigma","dcbSigma",prevSigma, 1.1, 10.);
    RooRealVar* dcbAlphaL = new RooRealVar("dcbAlphaL","dcbAlphaL",prevAlphaL, zero, 25.);
    RooRealVar* dcbNL = new RooRealVar("dcbNL","dcbNL",prevNL, zero, 25.);
    RooRealVar* dcbAlphaR = new RooRealVar("dcbAlphaR","dcbAlphaR",prevAlphaR, zero, 25.);
    RooRealVar* dcbNR = new RooRealVar("dcbNR","dcbNR",prevNR, zero, 25.);
    RooCrystalBall* dcb = new RooCrystalBall("dcb","dcb",*mass,*dcbMean,*dcbSigma,
                                          *dcbAlphaL, *dcbNL,
                                          *dcbAlphaR, *dcbNR);

    RooDataHist dh1("dh1","dh1",*mass,h);

    RooFFTConvPdf* bwdcb = new RooFFTConvPdf("bwdcb","BreitWigner DCB",
                                             *mass, *dcb, *breitW);
    bwdcb->setBufferFraction(0);
    RooFitResult* fit = bwdcb->fitTo(dh1,Range(75.,105.),
				     Save(true),PrintLevel(-1),
				     SumW2Error(false),Minos(false));

    double sigmaErrorConv =
      static_cast<RooRealVar*>(fit->floatParsFinal().find("dcbSigma"))->getError();

    prevMean =
      static_cast<RooRealVar*>(fit->floatParsFinal().find("dcbMean"))->getVal();
    prevSigma =
      static_cast<RooRealVar*>(fit->floatParsFinal().find("dcbSigma"))->getVal();
    prevAlphaL =
      static_cast<RooRealVar*>(fit->floatParsFinal().find("dcbAlphaL"))->getVal();
    prevAlphaR =
      static_cast<RooRealVar*>(fit->floatParsFinal().find("dcbAlphaR"))->getVal();
    prevNL =
      static_cast<RooRealVar*>(fit->floatParsFinal().find("dcbNL"))->getVal();
    prevNR =
      static_cast<RooRealVar*>(fit->floatParsFinal().find("dcbNR"))->getVal();

    std::cout << Form("\n\n\t%s\n\tSigma:%.4f\tError:%.4f\n\n",hname.c_str(),prevSigma,sigmaErrorConv);
    std::cout << Form("\tMean:%.4f\tAlphaL:%.4f\tAlphaR:%.4f\tNL:%.4f\tNR:%.4f\n\n",
                      prevMean,prevAlphaL,prevAlphaR,prevNL,prevNR);

    std::string parLabel
      = Form("#sigma=%.3f|#mu=%.3f|#alpha_{L}=%.3f|#alpha_{R}=%.3f|N_{L}=%.3f|N_{R}=%.3f",
	     prevSigma,prevMean,prevAlphaL,prevAlphaR,prevNL,prevNR);

    std::string title = Form("%s [%.0f:%.0f] MC [%d];#splitline{Pt [GeV]}{%s};Event Count",
                             hname.c_str(),ptBinLow, ptBinHigh,year,parLabel.c_str());
    if(isData)
      title = Form("%s [%.0f:%.0f] Data [%d];#splitline{Pt [GeV] }{%s};Event Count",
                   hname.c_str(),ptBinLow, ptBinHigh,year,parLabel.c_str());

    RooPlot* frame = mass->frame(Title(title.c_str()));

    dh1.plotOn(frame);

    bwdcb->plotOn(frame,LineColor(kYellow));

    sigmas.emplace_back(prevSigma);
    sigmaErrors.emplace_back(sigmaErrorConv);

    RooPlot* framePull = mass->frame(Title(Form("Pull %s",title.c_str())));
    framePull->SetName(Form("fPull_%.0f_%.0f_%s_%d",ptBinLow,ptBinHigh,hname.c_str(),year));
    RooHist* hpull = frame->pullHist();
    hpull->SetName(Form("hPull_%.0f_%.0f_%s_%d",ptBinLow,ptBinHigh,hname.c_str(),year));
    framePull->addPlotable(hpull,"P");

    gStyle->SetOptFit(1111);
    frame->Draw();
    frame->GetYaxis()->SetRangeUser(0.,YLimit);
    if( ptBinLow >= 150. ){
      frame->GetYaxis()->SetRangeUser(0.1,YLimit*10.);
      gPad->SetLogy();
      gPad->Modified();
      gPad->Update();
    }

    cPull->cd(k);
    framePull->Draw();

    //fxDCB->Draw("SAME");
  }

  std::string fname = Form("%d_%s_Fits.png",year,hname.c_str());

  if(isData)
    fname = Form("%d_%s_Fits_Data.png",year,hname.c_str());

  c1->Print(fname.c_str());

  if(!isData){
    cPull->Print(Form("%d_%s_Pull_Fits.png",year,hname.c_str()));
  }else{
    cPull->Print(Form("%d_%s_Pull_Fits_Data.png",year,hname.c_str()));
  }

  const Int_t nPoints = nBins - 1; // - underflow bin

  TGraphAsymmErrors* g = new TGraphAsymmErrors(nPoints);
  g->SetName(Form("%d_%s_g",year,hname.c_str()));

  for(int i = 0; i < nPoints; ++i){
    Double_t mid = (ptBins[i]+ptBins[i+1])/2.;
    Double_t dx = mid - ptBins[i];
    g->SetPoint(i,mid,sigmas[i]);
    g->SetPointError(i,dx,dx,sigmaErrors[i],sigmaErrors[i]);
  }

  c1->Clear();
  g->Print();
  g->SetMarkerColor(4);
  g->SetMarkerStyle(21);
  g->SetTitle(Form("Mass Resolution [%d]; Pt; Mass Resolution",year));
  g->Draw("AP");
  g->Print();
  fname = Form("%d_%s.png",year,hname.c_str());
  if(isData)
    fname = Form("%d_%s_Data.png",year,hname.c_str());
  c1->Print(fname.c_str());
  delete c1;

  return g;

}

int StackMRes(){

  std::function<TH1*(TGraphAsymmErrors*)>
    convertToHisto = [] (TGraphAsymmErrors* g) {

      const int nPoints = g->GetN();

      Double_t* x = g->GetX();
      Double_t* y = g->GetY();

      std::vector<double> xBins;
      std::vector<double> values;
      std::vector<double> errors;
      for(int i = 0; i < nPoints; ++i){
        xBins.emplace_back(x[i] - g->GetErrorXlow(i));
        values.emplace_back(y[i]);
        errors.emplace_back(g->GetErrorYhigh(i));
        //std::cout << Form("[%d] Bin: %.1f\t%.1f\t%.1f\n",i,xBins[i],values[i],errors[i]);
      }
      xBins.emplace_back(x[nPoints-1] + g->GetErrorXhigh(nPoints-1));
      //std::cout << Form("[%d] Bin: %.1f\n",xBins.size()-1,xBins[xBins.size()-1]);
      values.emplace_back(y[nPoints]);
      errors.emplace_back(g->GetErrorYhigh(nPoints));

      TH1F* h = new TH1F(Form("hfr_%s",g->GetName()),"Ratio",nPoints,&xBins[0]);

      for (int i = 0; i < nPoints; ++i) {
        h->SetBinContent(i+1,values[i]);
        h->SetBinError(i+1,errors[i]);
      }

      return h;
  };

  std::function<TGraphAsymmErrors*(TGraphAsymmErrors*,TGraphAsymmErrors*)>
    getRatio  = [&] (TGraphAsymmErrors* gNum, TGraphAsymmErrors* gDen){
      TGraphAsymmErrors* gr = new TGraphAsymmErrors(convertToHisto(gNum),convertToHisto(gDen),"pois");
      return gr;
    };

  std::vector<std::string> Histos2D =
    { "HMassZPt_A_G", "HMassZPt_B_G",
      "HMassZPt_A_T","HMassZPt_B_T",
      "HMassZPt_A_GG_L1","HMassZPt_B_GG_L1",
      "HMassZPt_A_GT_L1","HMassZPt_B_GT_L1",
      "HDPMassZ_A", "HDPMassZ_B", "HDPMassZ_C"};

  std::unordered_map<std::string,std::string> hTitles = {
    { "HMassZPt_A_G", "Z Peak Resolution [globalHighPt] [ 0. < |#eta| < 1.2 ]" },
    { "HMassZPt_B_G", "Z Peak Resolution [globalHighPt] [ 1.2 < |#eta| < 2.4 ]" },
    { "HMassZPt_A_T", "Z Peak Resolution [trackerHighPt] [ 0. < |#eta| < 1.2 ]" },
    { "HMassZPt_B_T", "Z Peak Resolution [trackerHighPt] [ 1.2 < |#eta| < 2.4 ]" },
    { "HMassZPt_A_GG_L1", "Z Peak Resolution [global + global] [ 0. < |#eta| < 1.2 ]"},
    { "HMassZPt_B_GG_L1", "Z Peak Resolution [global + global] [ 1.2 < |#eta| < 2.4 ]"},
    { "HMassZPt_A_GT_L1", "Z Peak Resolution [global + tracker] [ 0. < |#eta| < 1.2 ]"},
    { "HMassZPt_B_GT_L1", "Z Peak Resolution [global + tracker] [ 1.2 < |#eta| < 2.4 ]"},
    { "HDPMassZ_A", "Both Muons |#eta|<1.2;Leading #mu P_{T} (GeV);#sigma_{#mu#mu} at Z peak (GeV)"},
    { "HDPMassZ_B", "Leading Muon 1.2<|#eta|<1.6;Leading #mu P_{T} (GeV);#sigma_{#mu#mu} at Z peak (GeV)"},
    { "HDPMassZ_C", "Leading Muon |#eta|>1.6;Leading #mu P_{T} (GeV);#sigma_{#mu#mu} at Z peak (GeV)"}
  };


  TCanvas* cp1 = new TCanvas("cp1","cp1", 4*500, 3*500);
  cp1->Divide(4,3);

  std::vector<Int_t> years = {2016,2017,2018};

  for(auto year: years){
    Int_t k = 1;
    for(auto hname : Histos2D){

      const Float_t leftMargin = 0.12;
      const Float_t rightMargin = 0.12;
      const Float_t topMargin = 0.12;
      const Float_t bottomMargin = 0.5;

      cp1->cd(k);
      auto mainPad = new TPad(Form("mainPad_%s",hname.c_str()),"mainPad",0.,0.25,1.,1.);
      mainPad->Draw();
      mainPad->SetBottomMargin(1e-3);
      mainPad->SetLeftMargin(leftMargin);
      mainPad->SetRightMargin(rightMargin);
      mainPad->SetTickx();
      mainPad->SetTicky();

      TMultiGraph *mg = new TMultiGraph();
      mg->SetName(Form("%d_mg_%s",year,hname.c_str()));
      mg->SetTitle(Form("%s [%d];Pt [GeV];Mass Resolution at Z Peak [GeV]",hTitles[hname].c_str(),year));

      TGraphAsymmErrors* gmc = plotFits(year, hname);
      gmc->SetLineColor(kRed);
      gmc->SetMarkerColor(kRed);
      std::cout << Form("GMC Printout %s:\n\n\n",gmc->GetName());
      gmc->Print();

      TGraphAsymmErrors* gdata = plotFits(year, hname, true);
      gdata->SetLineColor(kBlack);
      gdata->SetMarkerColor(kBlack);
      std::cout << Form("GDATA Printout %s:\n\n\n",gdata->GetName());
      gdata->Print();

      mainPad->cd();
      mg->Add(gmc,"P");
      mg->Add(gdata,"P");
      mg->Draw("AP");
      mg->GetYaxis()->SetRangeUser(0.,6.);

      TLegend *l = new TLegend(0.6,0.7,0.9,0.9);
      l->SetName(Form("%d",year));
      l->AddEntry(gmc,"MC");
      l->AddEntry(gdata,Form("%d Data",year));
      l->Draw();

      auto subPad = new TPad(Form("subPad_%s",hname.c_str()),"subPad",0.,0.,1.,0.25);
      subPad->SetLeftMargin(leftMargin);
      subPad->SetRightMargin(rightMargin);
      subPad->SetTopMargin(1e-3);
      subPad->SetBottomMargin(bottomMargin);
      subPad->SetGrid();
      subPad->SetFrameLineWidth(1);

      const int font = 43;
      const float fontSize = 15.;
      const float labelSize = 0.15;
      TGraphAsymmErrors* gRatio = getRatio(gdata,gmc);
      gRatio->GetXaxis()->SetTitleFont(font);
      gRatio->GetXaxis()->SetTitleSize(fontSize);
      gRatio->GetXaxis()->SetLabelSize(labelSize);
      gRatio->GetXaxis()->SetTitleOffset(12.0);
      gRatio->GetYaxis()->SetTitleOffset(4.0);
      gRatio->GetYaxis()->SetNdivisions(6,3,0);
      gRatio->GetYaxis()->SetTitle("Data/MC");
      gRatio->GetXaxis()->SetTitle("Leading Muon Pt [GeV]");
      gRatio->GetYaxis()->SetTitleFont(font);
      gRatio->GetYaxis()->SetTitleSize(fontSize);
      gRatio->GetYaxis()->SetLabelSize(labelSize);
      gRatio->SetMinimum(0.5);
      gRatio->SetMaximum(1.5);

      cp1->cd(k);
      subPad->Draw();
      subPad->cd();
      gRatio->Draw("AP");
      gRatio->GetXaxis()->SetRangeUser(mg->GetXaxis()->GetXmin(),
				       mg->GetXaxis()->GetXmax());


      std::cout << k << "\t" << hname << "\n";
      ++k;
    }
    cp1->Print(Form("%d_.png",year));
  }

  return 0;
}
