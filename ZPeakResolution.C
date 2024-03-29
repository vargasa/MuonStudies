#include "TFile.h"
#include "ZPeakResolution.h"
#include "TMath.h"
#include "TLegend.h"
#include "TError.h"
#include "TVector2.h"
#include <functional>
#include <bitset>

#define FillS(xx) Fill(xx,1.)

ZPeakResolution::ZPeakResolution(TTree *)
{

  HCutFlow = 0;

  HMassZPt_A_G = 0;
  HMassZPt_B_G = 0;
  HMassZPt_A_T = 0;
  HMassZPt_B_T = 0;

  HMassZPt_A_GG = 0;
  HMassZPt_A_GT = 0;
  HMassZPt_B_GG = 0;
  HMassZPt_B_GT = 0;

  HDPMassZ_A = 0;
  HDPMassZ_B = 0;
  HDPMassZ_C = 0;

  HMuonPtl1 = 0;
  HMuonPtl2 = 0;

  HMassZ = 0;

  HIp3dl1 = 0;
  HIp3dl2 = 0;
}

void ZPeakResolution::Begin(TTree *tree) {

  if (fInput->FindObject("SampleName")) {
    // Lesson: TString can't be in TCollection
    TNamed *p = dynamic_cast<TNamed *>(fInput->FindObject("SampleName"));
    SampleName = p->GetTitle();
  }
}

void ZPeakResolution::SlaveBegin(TTree *tree) {

  TH1::SetDefaultSumw2();

  HCutFlow = new TH1D("HCutFlow","",50,0.,50.);  /* Limits are meaningless here */
  fOutput->Add(HCutFlow);

  const Double_t PBins[16] = {
    52.,72.,100.,150.,200.,
    300.,400.,600.,800.,1000.,
    1300.,1600.,2000.,2500,3100,
    4500
  };

  const Double_t PResBins[198] = {
    -0.99,-0.98,-0.97,-0.96,-0.95,-0.94,
    -0.93,-0.92,-0.91,-0.9,-0.89,-0.88,
    -0.87,-0.86,-0.85,-0.84,-0.83,-0.82,
    -0.81,-0.8,-0.79,-0.78,-0.77,-0.76,
    -0.75,-0.74,-0.73,-0.72,-0.71,-0.7,
    -0.69,-0.68,-0.67,-0.66,-0.65,-0.64,
    -0.63,-0.62,-0.61,-0.6,-0.59,-0.58,
    -0.57,-0.56,-0.55,-0.54,-0.53,-0.52,
    -0.51,-0.5,-0.49,-0.48,-0.47,-0.46,
    -0.45,-0.44,-0.43,-0.42,-0.41,-0.4,
    -0.39,-0.38,-0.37,-0.36,-0.35,-0.34,
    -0.33,-0.32,-0.31,-0.3,-0.29,-0.28,
    -0.27,-0.26,-0.25,-0.24,-0.23,-0.22,
    -0.21,-0.2,-0.19,-0.18,-0.17,-0.16,
    -0.15,-0.14,-0.13,-0.12,-0.11,-0.1,
    -0.09,-0.08,-0.07,-0.06,-0.05,-0.04,
    -0.03,-0.02,-0.01,7.5287e-16,0.01,0.02,
    0.03,0.04,0.05,0.06,0.07,0.08,
    0.09,0.1,0.11,0.12,0.13,0.14,
    0.15,0.16,0.17,0.18,0.19,0.2,
    0.21,0.22,0.23,0.24,0.25,0.26,
    0.27,0.28,0.29,0.3,0.31,0.32,
    0.33,0.34,0.35,0.36,0.37,0.38,
    0.39,0.4,0.41,0.42,0.43,0.44,
    0.45,0.46,0.47,0.48,0.49,0.5,
    0.51,0.52,0.53,0.54,0.55,0.56,
    0.57,0.58,0.59,0.6,0.61,0.62,
    0.63,0.64,0.65,0.66,0.67,0.68,
    0.69,0.7,0.71,0.72,0.73,0.74,
    0.75,0.76,0.77,0.78,0.79,0.8,
    0.81,0.82,0.83,0.84,0.85,0.86,
    0.87,0.88,0.89,0.9,0.91,0.92,
    0.93,0.94,0.95,0.96,0.97,0.98,
  };

  HMuonPtl1 = new TH1F("HMuonPtl1","",50,0,5000);
  fOutput->Add(HMuonPtl1);

  HMuonPtl2 = new TH1F("HMuonPtl2","",50,0,5000);
  fOutput->Add(HMuonPtl2);


  const Double_t PtBins_MRes[8] = {
    53,75,100,150,200,300,450,800
  };

  const Double_t ZMassBins[19] = {
    70., 72., 74., 76., 78.,
    80., 82., 84., 86., 88.,
    90., 92., 94., 96., 98.,
    100., 102., 104., 106.
  };

  HMassZPt_A_G = new TH2F("HMassZPt_A_G","", 7, PtBins_MRes, 18, ZMassBins);
  fOutput->Add(HMassZPt_A_G);

  HMassZPt_B_G = static_cast<TH2F*>(HMassZPt_A_G->Clone());
  HMassZPt_B_G->SetName("HMassZPt_B_G");
  fOutput->Add(HMassZPt_B_G);

  HMassZPt_A_T = static_cast<TH2F*>(HMassZPt_A_G->Clone());
  HMassZPt_A_T->SetName("HMassZPt_A_T");
  fOutput->Add(HMassZPt_A_T);

  HMassZPt_B_T = static_cast<TH2F*>(HMassZPt_A_G->Clone());
  HMassZPt_B_T->SetName("HMassZPt_B_T");
  fOutput->Add(HMassZPt_B_T);

  HDPMassZ_A = static_cast<TH2F*>(HMassZPt_A_G->Clone());
  HDPMassZ_A->SetName("HDPMassZ_A");
  fOutput->Add(HDPMassZ_A);

  HDPMassZ_B = static_cast<TH2F*>(HMassZPt_A_G->Clone());
  HDPMassZ_B->SetName("HDPMassZ_B");
  fOutput->Add(HDPMassZ_B);

  HDPMassZ_C = static_cast<TH2F*>(HMassZPt_A_G->Clone());
  HDPMassZ_C->SetName("HDPMassZ_C");
  fOutput->Add(HDPMassZ_C);

  HMassZ = new TH1F("HMassZ","", 100,50,10000);
  fOutput->Add(HMassZ);

  HMassZPt_A_GG = new TH3F("HMassZPt_A_GG","", 7, PtBins_MRes, 7, PtBins_MRes, 18, ZMassBins);
  fOutput->Add(HMassZPt_A_GG);

  HMassZPt_B_GG = new TH3F("HMassZPt_B_GG","", 7, PtBins_MRes, 7, PtBins_MRes, 18, ZMassBins);
  fOutput->Add(HMassZPt_B_GG);

  HMassZPt_A_GT = new TH3F("HMassZPt_A_GT","", 7, PtBins_MRes, 7, PtBins_MRes, 18, ZMassBins);
  fOutput->Add(HMassZPt_A_GT);

  HMassZPt_B_GT = new TH3F("HMassZPt_B_GT","", 7, PtBins_MRes, 7, PtBins_MRes, 18, ZMassBins);
  fOutput->Add(HMassZPt_B_GT);

  HIp3dl1 = new TH1F("HIp3dl1","",400,-0.2,0.2);
  fOutput->Add(HIp3dl1);

  HIp3dl2 = static_cast<TH1F*>(HIp3dl1->Clone());
  HIp3dl2->SetName("HIp3dl2");
  fOutput->Add(HIp3dl2);

  HMuonPt = new TH2F("HMuonPt","HMuonPt", 90., 0., 450., 90., 0., 450.);
  fOutput->Add(HMuonPt);

  HMuonEta = new TH2F("HMuonEta","HMuonEta", 24, -2.4, 2.4, 24, -2.4, 2.4);
  fOutput->Add(HMuonEta);

  HMuonPhi = new TH2F("HMuonPhi","HMuonPhi", 32, -3.15, 3.15, 32, -3.15, 3.15);
  fOutput->Add(HMuonPhi);
}

void ZPeakResolution::SortByDescPt(std::vector<UInt_t>& GoodIdx, const Leptons& l){
   std::function<bool(const int&, const int&)> sortDescPt =
    [&](const Int_t& i, const Int_t& j){
      return !(l.pt[i] < l.pt[j]);
    };

   std::sort(GoodIdx.begin(),GoodIdx.end(),
            sortDescPt);
}

std::vector<UInt_t> ZPeakResolution::GetGoodMuon(const Muons& Mu){
  std::vector<UInt_t> GoodIndex = {};
  if(!MuonTest()) return GoodIndex;
  if(*Mu.n == 0) return GoodIndex;

  leadMuIdx = LeadingIdx(Mu);

  if(Mu.pt[leadMuIdx] < 53.)  /* HLT_Mu50_OR_HLT_TkMu50 lower pt limit from SFDB*/
    return GoodIndex;
  const Float_t MaxEta = 2.4;
  const Float_t MinPt = 53.;
  if(abs(Mu.eta[leadMuIdx]) >= MaxEta){
    HCutFlow->FillS("LeadingMuOut");
    return GoodIndex;
  }
  GoodIndex.reserve(10);
  const int MuonPdgId = 13;
  for (UInt_t i=0; i<*Mu.n;++i){
    if( Muon_highPurity[i] and Muon_highPtId[i] >=1 and
        abs(Mu.eta[i]) < MaxEta and
        Mu.pt[i] > MinPt and
        Muon_tkRelIso[i] < 5e-2
#ifndef CMSDATA
        and abs(Mu.pdgId[i]) == MuonPdgId
#endif
        )
      GoodIndex.emplace_back(i);
  }

  SortByDescPt(GoodIndex,Mu);

  return GoodIndex;
}

Int_t ZPeakResolution::LeadingIdx(const Leptons& l) {

  Int_t idx = -1;
  Double_t pt = -1.;

  for(uint i = 0; i < *l.n; ++i){
    if( l.pt[i] > pt ){
      idx = i;
      pt = l.pt[i];
    }
  }

  return idx;

}



std::vector<std::pair<UInt_t,UInt_t>> ZPeakResolution::GetMuonPairs(const Muons& m) const{

  // return leading pair sorted by Pt

  std::vector<std::pair<UInt_t,UInt_t>> pairs;

  for(uint i = 0; i < GoodMuon.size(); ++i){
    for(uint j = i+1; j < GoodMuon.size(); ++j){
      if (m.charge[i] != m.charge[j]) {
        if(m.pt[i] > m.pt[j]) {
          pairs.push_back(std::make_pair(GoodMuon[i],GoodMuon[j]));
        } else {
          pairs.push_back(std::make_pair(GoodMuon[j],GoodMuon[i]));
        }
      }
    }
  }

  return pairs;

};

float ZPeakResolution::MassRecoZ(const float& pt1, const float& eta1, const float& phi1, const float& m1,
                             const float& pt2, const float& eta2, const float& phi2, const float& m2) const{

  ROOT::Math::PtEtaPhiMVector l1(pt1, eta1, phi1, m1);
  ROOT::Math::PtEtaPhiMVector l2(pt2, eta2, phi2, m2);

  return (l1+l2).M();
};

ZPairInfo ZPeakResolution::FindZ(const Muons& mu) const{
  return FindZ(ZPeakResolution::GetMuonPairs(mu), mu, GoodMuon);
}

ZPairInfo ZPeakResolution::FindZ(const std::vector<std::pair<UInt_t,UInt_t>>& Pairs,
                             const Leptons& l,
                             const std::vector<UInt_t>& GoodLepton) const{

  ZPairInfo z1;

  if(Pairs.size() == 0)
    return z1;

  std::vector<ZPairInfo> PairsWithMass;


  for(const auto& pair: Pairs){
    UInt_t i = pair.first;
    UInt_t j = pair.second;
    Float_t m = ZPeakResolution::MassRecoZ(l.pt[i],l.eta[i],l.phi[i],l.mass,
                                       l.pt[j],l.eta[j],l.phi[j],l.mass);
    z1.Mass = m;
    z1.Pair = pair;
    PairsWithMass.emplace_back(z1);
  }


  return PairsWithMass[0]; /*Pair with smallest delta*/

}

#ifndef CMSDATA
std::pair<Int_t,Int_t> ZPeakResolution::GetMother(std::pair<Int_t,Int_t> Daughter) const{
  return ZPeakResolution::GetMother(Daughter.first,Daughter.second);
}
#endif
#ifndef CMSDATA
std::pair<Int_t,Int_t> ZPeakResolution::GetMother(Int_t GenPartIdx, Int_t PdgId /*\Mu or e*/) const{

  if(GenPartIdx == -1 ) return std::make_pair(GenPartIdx,9999);

  std::pair<Int_t,Int_t> Mother;
  Mother.first = GenPart_genPartIdxMother[GenPartIdx];
  Mother.second = GenPart_pdgId[Mother.first];

  if(abs(Mother.second) == abs(PdgId)){
    Mother = ZPeakResolution::GetMother(Mother);
  }

  return Mother;
}
#endif

Bool_t ZPeakResolution::CheckMuonPair(const std::pair<UInt_t,UInt_t>& p) const{

  const Float_t MinSubleadPt = 10.;
  if (Muon_pt[p.second] < MinSubleadPt) return kFALSE;
  const Float_t farFromPV = 1e-2;
  if (Muon_ip3d[p.first] > farFromPV or Muon_ip3d[p.second] > farFromPV) return kFALSE;

  Bool_t GlobalHighPtl1 = (Muon_highPtId[p.first] == 2);
  Bool_t GlobalHighPtl2 = (Muon_highPtId[p.second] == 2);

  if (!(GlobalHighPtl1 or GlobalHighPtl2)) {
    HCutFlow->FillS("FailZGlbHighPt");
    return kFALSE;
  } else {
    if (GlobalHighPtl1 and !Muon_isPFcand[p.first]) {
      HCutFlow->FillS("FailZGlbHighPt&PFCand");
      return kFALSE;
    } else if (GlobalHighPtl2 and !Muon_isPFcand[p.second]) {
      HCutFlow->FillS("FailZGlbHighPt&PFCand");
      return kFALSE;
    }
  }
  return kTRUE;
}


Bool_t ZPeakResolution::Process(Long64_t entry) {

  l1 = l2 = -1;
  PairMu = false;
  PairZMass = -1.;

  ReadEntry(entry);

  double genW = 1.;

#ifndef CMSDATA
  genW = *genWeight;
  HCutFlow->Fill("genWeight",genW);
#endif

  HCutFlow->FillS("NoCuts");

  if (!MuonTest()){
    HCutFlow->FillS("FailMuonHLTs");
    return kFALSE;
  }

  if (!FlagsTest()){
    HCutFlow->FillS("FailFlags");
    return kFALSE;
  }

#ifndef CMSDATA
  Muons Mus(nMuon,Muon_tunepRelPt,Muon_pt,Muon_eta,Muon_phi,
            Muon_charge,Muon_dxy,Muon_dz,Muon_pfRelIso03_all,
            Muon_ip3d,Muon_sip3d,
            Muon_tightId, Muon_genPartIdx, Muon_pdgId);
#endif
#ifdef CMSDATA
  Muons Mus(nMuon,Muon_tunepRelPt,Muon_pt,Muon_eta,Muon_phi,
            Muon_charge,Muon_dxy,Muon_dz,Muon_pfRelIso03_all,
            Muon_ip3d,Muon_sip3d,
            Muon_tightId);
#endif

  GoodMuon = ZPeakResolution::GetGoodMuon(Mus);

  if( (GoodMuon.size()) < 2 ){
    HCutFlow->FillS("goodMu<2");
    return kFALSE;
  }

  ////////////// Define Z //////////////
  const Float_t MinZMass = 70.;
  const Float_t MaxZMass = 111.;


  ZPairInfo ztmu;

  ztmu = FindZ(Mus);
  if(!ztmu.empty() && CheckMuonPair(ztmu.Pair)) PairMu = true;

  if(PairMu){
     HCutFlow->FillS("PairMu");
  } else {
    HCutFlow->FillS("NoPair");
    return kFALSE;
  }

   l1 = ztmu.Pair.first;
   l2 = ztmu.Pair.second;

   lep1 = PtEtaPhiMVector(Mus.pt[l1],Mus.eta[l1],
                          Mus.phi[l1],Muons::mass);
   lep2 = PtEtaPhiMVector(Mus.pt[l2],Mus.eta[l2],
                          Mus.phi[l2],Muons::mass);

   zb   = lep1 + lep2;

   Double_t w = 1.;

   HMuonPtl1->Fill(lep1.Pt(), genW);
   HMuonPtl2->Fill(lep2.Pt(), genW);
   HMassZ->Fill(zb.M(), genW);

   HCutFlow->FillS("ZCandidate");

   auto printEventInfo = [&](std::string histoName){
     std::cout << Form("%s\tRun: %u\tEvent: %llu\tLumiblock: %u\tmll: %.4f\n",
                       histoName.c_str(), *run, *event, *luminosityBlock, zb.M());

   };

   auto FillHistos53 = [&](TH2F* A_G, TH2F* B_G,
                           TH2F* A_T, TH2F* B_T) {
     TH2F *hl1, *hl2 = NULL;

     std::pair<Float_t,Float_t> etaLimit = {1.2, 2.4 };
     const Int_t globalId = 2;

     if( abs(lep1.Eta()) <= etaLimit.first ) { // A
       if(Muon_highPtId[l1] == globalId){
         hl1 = A_G;
       } else {
         hl1 = A_T;
       }
       if(Muon_highPtId[l2] == globalId){
         hl2 = A_G;
       } else {
         hl2 = A_T;
       }
     } else if ( abs(lep1.Eta()) > etaLimit.first and
                 abs(lep1.Eta()) <= etaLimit.second ) { //B
       if(Muon_highPtId[l1] == globalId){
         hl1 = B_G;
       } else {
         hl1 = B_T;
       }
       if(Muon_highPtId[l2] == globalId){
         hl2 = B_G;
       } else {
         hl2 = B_T;
       }
     }


     hl1->Fill(lep1.Pt(), zb.M(), genW);
     hl2->Fill(lep2.Pt(), zb.M(), genW);

#ifdef CMSDATA
     if(lep1.Pt() > 300. and lep1.Pt() < 450.)
       printEventInfo(hl1->GetName());
#endif

     HIp3dl1->Fill(Muon_ip3d[l1]);
     HIp3dl2->Fill(Muon_ip3d[l2]);

     if ( Muon_highPtId[l1] == globalId
          and Muon_highPtId[l2] == globalId ) {
       if ( abs(lep1.Eta()) > etaLimit.first and
            abs(lep1.Eta()) <= etaLimit.second ) {
         HMassZPt_B_GG->Fill(lep1.Pt(),lep2.Pt(),zb.M(), genW);
       } else {
         HMassZPt_A_GG->Fill(lep1.Pt(),lep2.Pt(),zb.M(), genW);
       }
     } else {
       if ( abs(lep1.Eta()) > etaLimit.first and
            abs(lep1.Eta()) <= etaLimit.second ) {
         HMassZPt_B_GT->Fill(lep1.Pt(),lep2.Pt(),zb.M(), genW);
       } else {
         HMassZPt_A_GT->Fill(lep1.Pt(),lep2.Pt(),zb.M(), genW);
       }
     }

   };

   auto FillHistosDP = [&](TH2F* A, TH2F* B, TH2F* C){

     std::pair<float,float> etaLimits = { 1.2, 1.6 };
     const int globalId = 2;

     TH2F* theOne;

     if ( Muon_highPtId[l1] == globalId
          and Muon_highPtId[l2] == globalId) {
       if ( abs(lep1.Eta()) < etaLimits.first
            and abs(lep2.Eta()) < etaLimits.first){
         theOne = A;
       } else if ( abs(lep1.Eta()) > etaLimits.first
                   and abs(lep2.Eta()) < etaLimits.second ){
         theOne = B;
       } else if ( abs(lep1.Eta()) > etaLimits.second ){
         theOne = C;
       }
       theOne->Fill(lep1.Pt(),zb.M(),genW);
       theOne->Fill(lep2.Pt(),zb.M(),genW);
     }
   };


   std::pair<float,float> MassWindow = { 75., 110. };

   if ( zb.M() > MassWindow.first and zb.M() < MassWindow.second ){
     FillHistos53(HMassZPt_A_G, HMassZPt_B_G,
                  HMassZPt_A_T, HMassZPt_B_T);

     FillHistosDP(HDPMassZ_A,HDPMassZ_B,HDPMassZ_C);

     HMuonPt->Fill(lep1.Pt(),lep2.Pt());
     HMuonEta->Fill(lep1.Eta(),lep2.Pt());
     HMuonPhi->Fill(lep1.Phi(),lep2.Phi());
   }

   return kTRUE;
}

void ZPeakResolution::Terminate() {

#ifdef Y2016
  const Int_t Year_ = 2016;
#elif defined(Y2017)
  const Int_t Year_ = 2017;
#elif defined(Y2018)
  const Int_t Year_ = 2018;
#endif

  std::unique_ptr<TFile> fOut(TFile::Open("MuonStudies_ZPeakResolution_MoreHistos.root","UPDATE"));
  fOut->mkdir(Form("%d",Year_));
  fOut->mkdir(Form("%d/%s",Year_,SampleName.Data()));
  fOut->cd(Form("%d/%s",Year_,SampleName.Data()));

  TObjLink *lnk = fOutput->FirstLink();
  while (lnk) {
    lnk->GetObject()->Write();
    lnk = lnk->Next();
  }

  fOut->Write();
  fOut->Close();

}
