#pragma once

#include "NuHepMC/HepMC3Features.hxx"

#include "NuHepMC/Constants.hxx"
#include "NuHepMC/EventUtils.hxx"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

#include "TH2.h"

#include <iostream>
#include <vector>

enum Classification {
  k1p_only = 0,
  k1n_only,
  k1pi0_1p,
  k1piplus_1p,
  k1p_1n,
  k1p_any_n,
  k1p_1gamma,
  k2p_only,
  k2p_any_n,
  k3p_any_n,
  kany_n_only,
  k1pi0_any_p,
  k1pi0_any_n,
  k1pi0_any_np,
  kother,
  kNumClass
};

std::vector<Classification> pclasses = {k1p_only, k1n_only, k1pi0_1p,
                                        k1piplus_1p};

std::ostream &operator<<(std::ostream &os, Classification c) {
  switch (c) {
  case k1p_only:
    return os << "1p";
  case k1p_1n:
    return os << "1p, 1n";
  case k1p_any_n:
    return os << "1p, 2+n";
  case k1p_1gamma:
    return os << "1p, 1gamma";
  case k2p_only:
    return os << "2p";
  case k2p_any_n:
    return os << "2p, 1+n";
  case k3p_any_n:
    return os << "3p, 1+n";
  case k1n_only:
    return os << "1n";
  case kany_n_only:
    return os << "1+n";
  case k1pi0_1p:
    return os << "1pi0, 1p";
  case k1piplus_1p:
    return os << "1pi+, 1p";
  case k1pi0_any_p:
    return os << "1pi0, 1+p";
  case k1pi0_any_n:
    return os << "1pi0, 1+n";
  case k1pi0_any_np:
    return os << "1pi0, 1+n, 1+p";
  case kother:
    return os << "other";
  }
  throw;
}

std::string to_string(Classification c) {
  switch (c) {
  case k1p_only:
    return "k1p_only";
  case k1p_1n:
    return "k1p_1n";
  case k1p_1gamma:
    return "k1p_1gamma";
  case k2p_only:
    return "k2p_only";
  case k1n_only:
    return "k1n_only";
  case k1pi0_1p:
    return "k1pi0_1p";
  case k1piplus_1p:
    return "k1piplus_1p";
  case k1pi0_any_p:
    return "k1pi0_any_p";
  case k1pi0_any_n:
    return "k1pi0_any_n";
  case k1pi0_any_np:
    return "k1pi0_any_np";
  case k1p_any_n:
    return "k1p_any_n";
  case k2p_any_n:
    return "k2p_any_n";
  case k3p_any_n:
    return "k3p_any_n";
  case kany_n_only:
    return "kany_n_only";
  case kother:
    return "kother";
  }
  throw;
}

std::pair<std::string, std::string> TransparencyName(Classification c,
                                                     std::string suffix = "") {

  std::string primpart;
  switch (c) {
  case k1p_only: {
    primpart = "proton";
    break;
  }
  case k1n_only: {
    primpart = "neutron";
    break;
  }
  case k1pi0_1p: {
    primpart = "pi0";
    break;
  }
  case k1piplus_1p: {
    primpart = "piplus";
    break;
  }
  default:
    throw;
  }

  return {to_string(c) + "_" + primpart + "_transp" + suffix,
          to_string(c) + "_" + primpart + "_all" + suffix};
}

Classification
GetClassification(std::vector<HepMC3::ConstGenParticlePtr> const &particles,
                  double ToGeV) {

  int nprotons = 0;
  int nneutrons = 0;
  int npi0 = 0;
  int npip = 0;
  int nlep = 0;
  int nother = 0;

  for (auto const &pt : particles) {
    switch (pt->pid()) {
    case 2212: {
      nprotons++;
      break;
    }
    case 2112: {
      nneutrons++;
      break;
    }
    case 111: {
      npi0++;
      break;
    }
    case 211: {
      npip++;
      break;
    }
    case 11:
    case -11:
    case 12:
    case -12:
    case 13:
    case -13:
    case 14:
    case -14: {
      nlep++;
      break;
    }
    case 22: { // ignore gammas < 15 MeV
      if ((pt->momentum().e() * ToGeV) > 0.015) {
        return kother;
      }
      break;
    }
    default: {
      if (pt->pid() < 1E6) {
        return kother;
      }
    }
    }
  }

  if ((nlep > 1) || (npi0 > 1) || (npip > 1)) {
    return kother;
  }

  if (npi0 != 0) {
    if ((nneutrons + nprotons) == 0) {
      return kother;
    }
    if (npip != 0) {
      return kother;
    }
    return nneutrons ? (nprotons ? k1pi0_any_np : k1pi0_any_n)
                     : (nprotons == 1 ? k1pi0_1p : k1pi0_any_p);
  }
  if (npip != 0) {
    if ((nprotons == 1) && (nneutrons == 0)) {
      return k1piplus_1p;
    }
    return kother;
  }

  switch (nprotons) {
  case 0: {
    return (nneutrons == 1) ? k1n_only : kany_n_only;
  }
  case 1: {

    return nneutrons ? ((nneutrons == 1) ? k1p_1n : k1p_any_n) : k1p_only;
  }
  case 2: {
    return nneutrons ? k2p_any_n : k2p_only;
  }
  case 3: {
    return k3p_any_n;
  }
  default: {
    return kother;
  }
  }
}

std::vector<HepMC3::ConstGenParticlePtr>
GetPreFSIParticles(HepMC3::GenEvent &evt, bool isGENIE = false) {
  if (isGENIE) {
    std::vector<HepMC3::ConstGenParticlePtr> prefsiparts;
    for (auto const &pt : evt.particles()) {
      if (pt->status() == 26) {
        prefsiparts.push_back(pt);
      }
    }
    for (auto const &pt :
         NuHepMC::Event::GetPrimaryVertex(evt)->particles_out()) {
      if ((std::abs(pt->pid()) >= 11) && (std::abs(pt->pid()) <= 16)) {
        prefsiparts.push_back(pt);
      }
    }
    return prefsiparts;
  } else {
    return NuHepMC::Event::GetPrimaryVertex(evt)->particles_out();
  }
}

Classification PrimaryClassification(HepMC3::GenEvent &evt, double ToGeV,
                                     bool isGENIE = false) {
  return GetClassification(GetPreFSIParticles(evt, isGENIE), ToGeV);
}

Classification FSClassification(HepMC3::GenEvent &evt, double ToGeV) {
  std::vector<HepMC3::ConstGenParticlePtr> fsparts;
  for (auto const &pt : evt.particles()) {
    if (pt->status() == NuHepMC::ParticleStatus::UndecayedPhysical) {
      fsparts.push_back(pt);
    }
  }
  return GetClassification(fsparts, ToGeV);
}

void RowNormTH2(TH2 *h2) {
  for (int j = 0; j < h2->GetYaxis()->GetNbins(); ++j) {
    double sum = 0;
    for (int i = 0; i < h2->GetXaxis()->GetNbins(); ++i) {
      h2->SetBinContent(i + 1, j + 1, h2->GetBinContent(i + 1, j + 1));
      h2->SetBinError(i + 1, j + 1, h2->GetBinError(i + 1, j + 1));
      sum += h2->GetBinContent(i + 1, j + 1);
    }

    if (sum) {
      for (int i = 0; i < h2->GetXaxis()->GetNbins(); ++i) {
        h2->SetBinContent(i + 1, j + 1, h2->GetBinContent(i + 1, j + 1) / sum);
        h2->SetBinError(i + 1, j + 1, h2->GetBinError(i + 1, j + 1) / sum);
      }
    }
  }
}

TH1D *CutOffZeroBin(TH1 *h, int rebinx = 1) {
  std::vector<double> xbins;
  xbins.push_back(h->GetXaxis()->GetBinLowEdge(2));
  for (int j = 1; j < h->GetXaxis()->GetNbins(); ++j) {
    xbins.push_back(h->GetXaxis()->GetBinUpEdge(j + 1));
  }

  TH1D *hc = new TH1D((std::string(h->GetName()) + "_nozero").c_str(),
                      (std::string(";") + h->GetXaxis()->GetTitle() + ";" +
                       h->GetYaxis()->GetTitle())
                          .c_str(),
                      xbins.size() - 1, xbins.data());

  for (int i = 0; i < hc->GetXaxis()->GetNbins(); ++i) {
    hc->SetBinContent(i + 1, h->GetBinContent(i + 2));
    hc->SetBinError(i + 1, h->GetBinError(i + 2));
  }
  if (rebinx > 1) {
    hc->RebinX(rebinx);
  }
  hc->SetDirectory(nullptr);
  return hc;
}