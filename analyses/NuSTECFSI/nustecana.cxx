// Leave this at the top to enable features detected at build time in headers in
// HepMC3
#include "NuHepMC/HepMC3Features.hxx"

#include "commonana.hxx"

#include "HepMC3/ReaderFactory.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"
#include "HepMC3/Setup.h"

#include "NuHepMC/Constants.hxx"
#include "NuHepMC/CrossSectionUtils.hxx"
#include "NuHepMC/EventUtils.hxx"
#include "NuHepMC/ReaderUtils.hxx"

#include <iostream>
#include <sstream>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

double ToGeV = 1;
double fatx = 1;

NuHepMC::StatusCodeDescriptors vtxstatus;
NuHepMC::StatusCodeDescriptors partstatus;
NuHepMC::StatusCodeDescriptors proc_ids;

// plots
std::unique_ptr<TH2D> TotalNeutronKE_1p_only;
std::unique_ptr<TH2D> TotalNeutralE_1p_only;
std::unique_ptr<TH2D> TotalPi0E_1piplus_1p;
std::unique_ptr<TH2D> TotalNeutralE_1piplus_1p;

// transparency
std::map<Classification,
         std::pair<std::unique_ptr<TH1D>, std::unique_ptr<TH1D>>>
    Transparency;
std::map<Classification,
         std::pair<std::unique_ptr<TH1D>, std::unique_ptr<TH1D>>>
    Transparency_5deg;

std::unique_ptr<TH2D> PrimaryToFinalStateSmearing;

std::pair<std::unique_ptr<TH1D>, std::unique_ptr<TH1D>>
TransparencyFact(Classification c, std::string suffix = "") {

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
  case k1pi0_1p:
  case k1pi0_any_n:
  case k1pi0_any_np: {
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

  auto tn = TransparencyName(c, suffix);

  return {std::make_unique<TH1D>(
              tn.first.c_str(),
              (";Primary " + primpart +
               " KE (GeV); Nuclear transparency (#theta_{deflect} < 5^{#circ})")
                  .c_str(),
              50, 0, 1),
          std::make_unique<TH1D>(
              tn.second.c_str(),
              (";Primary " + primpart +
               " KE (GeV); Nuclear transparency (#theta_{deflect} < 5^{#circ})")
                  .c_str(),
              50, 0, 1)};
}

std::pair<double, double> GetNeutronNeutralEnergy(HepMC3::GenEvent &evt) {
  std::pair<double, double> NeutronNeutralEnergy{0, 0};
  for (auto const &pt : evt.particles()) {
    if (pt->status() != NuHepMC::ParticleStatus::UndecayedPhysical) {
      continue;
    }

    switch (std::abs(pt->pid())) {
    case 111: {
      NeutronNeutralEnergy.second += pt->momentum().e();
    }
    case 2112: {
      double Tneut = (pt->momentum().e() - pt->momentum().m());
      NeutronNeutralEnergy.first += Tneut;
      NeutronNeutralEnergy.second += Tneut;
    }
    }
  }
  return NeutronNeutralEnergy;
}

HepMC3::ConstGenParticlePtr
GetPrimaryParticle(Classification c,
                   std::vector<HepMC3::ConstGenParticlePtr> const &parts) {
  for (auto const &part : parts) {
    switch (c) {
    case k1p_only: {
      if (part->pid() == 2212) {
        return part;
      }
    }
    case k1n_only: {
      if (part->pid() == 2112) {
        return part;
      }
    }
    case k1pi0_1p: {
      if (part->pid() == 111) {
        return part;
      }
    }
    case k1piplus_1p: {
      if (part->pid() == 211) {
        return part;
      }
    }
    }
  }
  return nullptr;
}

void ProcessEvent(HepMC3::GenEvent &evt) {

  auto pc_pos =
      std::find(pclasses.begin(), pclasses.end(), PrimaryClassification(evt));

  if (pc_pos == pclasses.end()) {
    return;
  }
  auto pclass = *pc_pos;

  auto fsclass = FSClassification(evt);
  PrimaryToFinalStateSmearing->Fill(fsclass, pclass);

  double w = evt.weights()[0];

  HepMC3::ConstGenParticlePtr primpart = GetPrimaryParticle(
      pclass, NuHepMC::Event::GetPrimaryVertex(evt)->particles_out());

  double pKE = (primpart->momentum().e() - primpart->momentum().m()) * ToGeV;

  if (fsclass == pclass) {
    HepMC3::ConstGenParticlePtr fspart = GetPrimaryParticle(
        pclass, NuHepMC::Event::GetParticles_All(
                    evt, NuHepMC::ParticleStatus::UndecayedPhysical));

    auto fs_mom = fspart->momentum();
    auto prim_mom = primpart->momentum();

    double costheta = (fs_mom.x() * prim_mom.x() + fs_mom.y() * prim_mom.y() +
                       fs_mom.z() * prim_mom.z()) /
                      (fs_mom.length() * prim_mom.length());

    double theta = std::acos((costheta > 1) ? 1 : costheta) * 180.0 / M_PI;

    if (theta < 5) {
      Transparency_5deg[pclass].first->Fill(pKE, w);
    }
    Transparency[pclass].first->Fill(pKE, w);
  } // end if topo stayed the same

  Transparency_5deg[pclass].second->Fill(pKE, w);
  Transparency[pclass].second->Fill(pKE, w);

  switch (pclass) {
  case k1p_only: {
    auto NeutronNeutralEnergy = GetNeutronNeutralEnergy(evt);
    TotalNeutronKE_1p_only->Fill(NeutronNeutralEnergy.first * ToGeV, pKE, w);
    TotalNeutralE_1p_only->Fill(NeutronNeutralEnergy.second * ToGeV, pKE, w);
    break;
  }
  case k1piplus_1p: {
    auto NeutronNeutralEnergy = GetNeutronNeutralEnergy(evt);
    TotalPi0E_1piplus_1p->Fill(
        (NeutronNeutralEnergy.second - NeutronNeutralEnergy.first) * ToGeV, pKE,
        w);
    TotalNeutralE_1piplus_1p->Fill(NeutronNeutralEnergy.second * ToGeV, pKE, w);
    break;
  }
  }
}

int main(int argc, char const *argv[]) {

  TH1::SetDefaultSumw2(true);

  if (argc < 2) {
    std::cout << "[RUNLIKE]: " << argv[0] << " <infile.hepmc3>" << std::endl;
    return 1;
  }

  std::string inf = argv[1];
  std::string out = argv[2];
  std::string dir = (argc > 3) ? argv[3] : "";

  auto rdr = HepMC3::deduce_reader(inf);
  if (!rdr) {
    std::cout << "Failed to instantiate HepMC3::Reader from " << inf
              << std::endl;
    return 1;
  }

  HepMC3::GenEvent evt;

  PrimaryToFinalStateSmearing = std::make_unique<TH2D>(
      "PrimaryToFinalStateSmearing",
      ";Final State topo.;Post-Hard Scatter topo.;Row-normalized", kNumClass, 0,
      kNumClass, pclasses.size(), 0, pclasses.size());

  TotalNeutronKE_1p_only = std::make_unique<TH2D>(
      "TotalNeutronKE_1p_only",
      ";#sum T_{neutron};T_{prot}^{preFSI};Row-normalized", 81, 0, 1, 50, 0, 1);
  TotalNeutralE_1p_only = std::make_unique<TH2D>(
      "TotalNeutralE_1p_only",
      ";#sum E_{neutral};T_{prot}^{preFSI};Row-normalized", 81, 0, 1, 50, 0, 1);
  TotalPi0E_1piplus_1p = std::make_unique<TH2D>(
      "TotalPi0E_1piplus_1p",
      ";#sum E_{#pi^{0}};T_{#pi+}^{preFSI};Row-normalized", 81, 0, 1, 50, 0, 1);
  TotalNeutralE_1piplus_1p = std::make_unique<TH2D>(
      "TotalNeutralE_1piplus_1p",
      ";#sum E_{neutral};T_{#pi+}^{preFSI};Row-normalized", 81, 0, 1, 50, 0, 1);

  Transparency[k1p_only] = TransparencyFact(k1p_only);
  Transparency[k1n_only] = TransparencyFact(k1n_only);
  Transparency[k1pi0_1p] = TransparencyFact(k1pi0_1p);
  Transparency[k1piplus_1p] = TransparencyFact(k1piplus_1p);

  Transparency_5deg[k1p_only] = TransparencyFact(k1p_only, "_lt5deg");
  Transparency_5deg[k1n_only] = TransparencyFact(k1n_only, "_lt5deg");
  Transparency_5deg[k1pi0_1p] = TransparencyFact(k1pi0_1p, "_lt5deg");
  Transparency_5deg[k1piplus_1p] = TransparencyFact(k1piplus_1p, "_lt5deg");

  size_t NEvents = 0;
  while (!rdr->failed()) {
    rdr->read_event(evt);

    if (!NEvents) { // can only reliably read run_info after reading an event,
                    // so do it on the first one
      ToGeV = NuHepMC::Event::ToMeVFactor(evt) * 1E-3;

      proc_ids = NuHepMC::GR4::ReadProcessIdDefinitions(evt.run_info());
      vtxstatus = NuHepMC::GR5::ReadVertexStatusIdDefinitions(evt.run_info());
      partstatus =
          NuHepMC::GR6::ReadParticleStatusIdDefinitions(evt.run_info());

      std::cout << "Input file reports that it contains "
                << NuHepMC::GC2::ReadExposureNEvents(evt.run_info())
                << " events" << std::endl;
      std::cout << "Processed " << NEvents << " events";
    }

    if (NEvents && !(NEvents % 10000)) {
      std::cout << "\r                                                ";
      std::cout << "\rProcessed " << NEvents << " events" << std::flush;
    }

    if (!rdr->failed()) {
      NEvents++;
    } else {
      break;
    }

    // if (NEvents > 1E5) {
    //   break;
    // }

    ProcessEvent(evt);
  }
  std::cout << "Processed " << NEvents << " events" << std::endl;

  for (int i = 0; i < PrimaryToFinalStateSmearing->GetXaxis()->GetNbins();
       ++i) {
    std::stringstream ss;
    ss << Classification(i);
    PrimaryToFinalStateSmearing->GetXaxis()->SetBinLabel(i + 1,
                                                         ss.str().c_str());
  }

  for (int i = 0; i < PrimaryToFinalStateSmearing->GetYaxis()->GetNbins();
       ++i) {
    std::stringstream ss;
    ss << pclasses[i];
    PrimaryToFinalStateSmearing->GetYaxis()->SetBinLabel(i + 1,
                                                         ss.str().c_str());
  }

  TFile fout(out.c_str(), "RECREATE");

  TDirectory *dout = &fout;
  if (dir.length()) {
    dout = fout.mkdir(dir.c_str());
  }

  RowNormTH2(PrimaryToFinalStateSmearing.get());
  dout->WriteObject(PrimaryToFinalStateSmearing.release(),
                    "PrimaryToFinalStateSmearing");

  RowNormTH2(TotalNeutronKE_1p_only.get());
  dout->WriteObject(CutOffZeroBin(TotalNeutronKE_1p_only.get()),
                    "TotalNeutronKE_1p_only_nozero");
  dout->WriteObject(TotalNeutronKE_1p_only.release(), "TotalNeutronKE_1p_only");
  RowNormTH2(TotalNeutralE_1p_only.get());
  dout->WriteObject(CutOffZeroBin(TotalNeutralE_1p_only.get()),
                    "TotalNeutralE_1p_only_nozero");
  dout->WriteObject(TotalNeutralE_1p_only.release(), "TotalNeutralE_1p_only");
  RowNormTH2(TotalPi0E_1piplus_1p.get());
  dout->WriteObject(CutOffZeroBin(TotalPi0E_1piplus_1p.get()),
                    "TotalPi0E_1piplus_1p_nozero");
  dout->WriteObject(TotalPi0E_1piplus_1p.release(), "TotalPi0E_1piplus_1p");
  RowNormTH2(TotalNeutralE_1piplus_1p.get());
  dout->WriteObject(CutOffZeroBin(TotalNeutralE_1piplus_1p.get()),
                    "TotalNeutralE_1piplus_1p_nozero");
  dout->WriteObject(TotalNeutralE_1piplus_1p.release(),
                    "TotalNeutralE_1piplus_1p");

  for (auto &a : Transparency) {
    if (a.second.first) {
      a.second.first->Divide(a.second.second.get());
      std::string name = a.second.first->GetName();
      dout->WriteObject(a.second.first.release(), name.c_str());
      name = a.second.second->GetName();
      dout->WriteObject(a.second.second.release(), name.c_str());
    }
  }
  for (auto &a : Transparency_5deg) {
    if (a.second.first) {
      a.second.first->Divide(a.second.second.get());
      std::string name = a.second.first->GetName();
      dout->WriteObject(a.second.first.release(), name.c_str());
      name = a.second.second->GetName();
      dout->WriteObject(a.second.second.release(), name.c_str());
    }
  }
}