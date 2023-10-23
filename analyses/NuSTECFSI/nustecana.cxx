// Leave this at the top to enable features detected at build time in headers in
// HepMC3
#include "NuHepMC/HepMC3Features.hxx"

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

enum Classification {
  k1p_only = 0,
  k1n_only,
  k1pi0_1p,
  k1p_1n,
  k1p_any_n,
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
std::vector<Classification> pclasses = {k1p_only, k1n_only, k1pi0_1p};

std::ostream &operator<<(std::ostream &os, Classification c) {
  switch (c) {
  case k1p_only:
    return os << "1p";
  case k1p_1n:
    return os << "1p, 1n";
  case k1p_any_n:
    return os << "1p, 2+n";
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
  case k2p_only:
    return "k2p_only";
  case k1n_only:
    return "k1n_only";
  case k1pi0_1p:
    return "k1pi0_1p";
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

Classification
GetClassification(std::vector<HepMC3::ConstGenParticlePtr> const &particles,
                  HepMC3::GenEvent &evt) {

  int nprotons = 0;
  int nneutrons = 0;
  int npi0 = 0;
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
    default: {
      if (pt->pid() < 1E6) {
        return {
            kother,
        };
      }
    }
    }
  }

  if ((nlep > 1) || (npi0 > 1)) {
    return kother;
  }
  if (npi0) {
    if (!(nneutrons + nprotons)) {
      return kother;
    }
    return nneutrons ? (nprotons ? k1pi0_any_np : k1pi0_any_n)
                     : (nprotons == 1 ? k1pi0_1p : k1pi0_any_p);
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

Classification PrimaryClassification(HepMC3::GenEvent &evt) {
  return GetClassification(
      NuHepMC::Event::GetPrimaryVertex(evt)->particles_out(), evt);
}

Classification FSClassification(HepMC3::GenEvent &evt) {
  std::vector<HepMC3::ConstGenParticlePtr> fsparts;
  for (auto const &pt : evt.particles()) {
    if (pt->status() == NuHepMC::ParticleStatus::UndecayedPhysical) {
      fsparts.push_back(pt);
    }
  }
  return GetClassification(fsparts, evt);
}

// binned in prim 1part kinematics, lead/sublead for 2212,2112,111
std::map<Classification, std::map<int, std::vector<std::unique_ptr<TH2D>>>> KE;
std::map<Classification, std::unique_ptr<TH2D>> EHadVis;
std::map<Classification, std::map<int, std::unique_ptr<TH2D>>> Multiplicity;

// transparency
std::map<Classification,
         std::map<int, std::pair<std::unique_ptr<TH1D>, std::unique_ptr<TH1D>>>>
    Transparency;

std::map<Classification, size_t> primtopocount;
std::unique_ptr<TH2D> PrimaryToFinalStateSmearing;

std::unique_ptr<TH2D> KEHistFact(Classification c, int pdg, int l) {

  std::string pdgs;
  switch (pdg) {
  case 2112: {
    pdgs = "neutron";
    break;
  }
  case 2212: {
    pdgs = "proton";
    break;
  }
  case 111: {
    pdgs = "pi0";
    break;
  }
  }

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
  default:
    throw;
  }

  std::string ld;
  switch (l) {
  case 0: {
    ld = "lead";
    break;
  }
  case 1: {
    ld = "subl";
    break;
  }
  case 2: {
    ld = "subsubl";
    break;
  }
  }

  return std::make_unique<TH2D>(
      (to_string(c) + "_" + ld + "_" + pdgs + "_KE").c_str(),
      (";KE " + pdgs + " (GeV);Primary " + primpart + " KE (GeV); P.D.F")
          .c_str(),
      50,  1E-6, 1, 9, 0.1, 1);
}

std::unique_ptr<TH2D> EHadVisHistFact(Classification c) {

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
  default:
    throw;
  }

  return std::make_unique<TH2D>(
      (to_string(c) + "_EHadVis").c_str(),
      (";EHadVis (GeV);Primary " + primpart + " KE (GeV); P.D.F").c_str(), 50, 1E-6,
      2, 9, 0.1, 1);
}

std::unique_ptr<TH2D> MultHistFact(Classification c, int pdg) {

  std::string pdgs;
  switch (pdg) {
  case 2112: {
    pdgs = "neutron";
    break;
  }
  case 2212: {
    pdgs = "proton";
    break;
  }
  case 111: {
    pdgs = "pi0";
    break;
  }
  }

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
  default:
    throw;
  }

  return std::make_unique<TH2D>(
      (to_string(c) + "_" + pdgs + "_mult").c_str(),
      (";" + pdgs + " Multiplicity;Primary " + primpart + " KE (GeV); P.D.F")
          .c_str(),
      10, 0, 10, 9, 0.1, 1);
}

std::pair<std::unique_ptr<TH1D>, std::unique_ptr<TH1D>>
TransparencyFact(Classification c) {

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
  default:
    throw;
  }

  return {std::make_unique<TH1D>(
              (to_string(c) + "_" + primpart + "_transp").c_str(),
              (";Primary " + primpart +
               " KE (GeV); Nuclear transparency (#theta_{deflect} < 5^{#circ})")
                  .c_str(),
              50, 0, 1),
          std::make_unique<TH1D>(
              (to_string(c) + "_" + primpart + "_all").c_str(),
              (";Primary " + primpart +
               " KE (GeV); Nuclear transparency (#theta_{deflect} < 5^{#circ})")
                  .c_str(),
              50, 0, 1)};
}

void RowNormTH2(TH2 *h2) {
  for (int j = 0; j < h2->GetYaxis()->GetNbins(); ++j) {
    double sum = 0;
    for (int i = 0; i < h2->GetXaxis()->GetNbins(); ++i) {
      h2->SetBinContent(i + 1, j + 1,
                        h2->GetBinContent(i + 1, j + 1) /
                            (h2->GetYaxis()->GetBinWidth(j + 1) *
                             h2->GetXaxis()->GetBinWidth(i + 1)));
      h2->SetBinError(i + 1, j + 1,
                      h2->GetBinError(i + 1, j + 1) /
                          (h2->GetYaxis()->GetBinWidth(j + 1) *
                           h2->GetXaxis()->GetBinWidth(i + 1)));
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

double GetEHadVis(HepMC3::GenEvent &evt) {
  double EHadVis = 0;
  for (auto const &pt : evt.particles()) {
    if (pt->status() != NuHepMC::ParticleStatus::UndecayedPhysical) {
      continue;
    }

    switch (std::abs(pt->pid())) {
    case 211:
    case 2212: {
      EHadVis += (pt->momentum().e() - pt->momentum().m());
      break;
    }
    case 111: {
      EHadVis += pt->momentum().e();
    }
    case NuHepMC::ParticleNumber::NuclearRemnant:
    case 13:
    case 14:
    case 2112: {
      break;
    }
    default: {
      std::cout << pt->pid() << std::endl;
    }
    }
  }
  return EHadVis;
}

void ProcessEvent(HepMC3::GenEvent &evt) {

  auto pc_pos =
      std::find(pclasses.begin(), pclasses.end(), PrimaryClassification(evt));

  if (pc_pos == pclasses.end()) {
    return;
  }
  auto pclass = *pc_pos;

  // HepMC3::Print::listing(evt);
  // std::cout << "prim class: " << pclass << std::endl;

  auto fsclass = FSClassification(evt);
  PrimaryToFinalStateSmearing->Fill(fsclass,
                                    std::distance(pclasses.begin(), pc_pos));

  double w = evt.weights()[0];

  primtopocount[pclass]++;

  constexpr std::array<int, 3> pids = {111, 2112, 2212};
  std::array<std::vector<double>, 3> KEs{};
  std::array<size_t, 3> Mult{};

  for (auto &part : evt.particles()) {
    if (part->status() != NuHepMC::ParticleStatus::UndecayedPhysical) {
      continue;
    }

    for (size_t pid_i = 0; pid_i < pids.size(); ++pid_i) {
      if (part->pid() != pids[pid_i]) {
        continue;
      }
      Mult[pid_i]++;
      KEs[pid_i].push_back((part->momentum().e() - part->momentum().m()) *
                           ToGeV);
    }
  }

  HepMC3::ConstGenParticlePtr primpart;
  for (auto const &part :
       NuHepMC::Event::GetPrimaryVertex(evt)->particles_out()) {
    switch (pclass) {
    case k1p_only: {
      if (part->pid() == 2212) {
        primpart = part;
      }
      break;
    }
    case k1n_only: {
      if (part->pid() == 2112) {
        primpart = part;
      }
      break;
    }
    case k1pi0_1p: {
      if (part->pid() == 111) {
        primpart = part;
      }
      break;
    }
    }
  }

  double pKE = (primpart->momentum().e() - primpart->momentum().m()) * ToGeV;

  if (!EHadVis.count(pclass)) {
    EHadVis[pclass] = EHadVisHistFact(pclass);
    EHadVis[pclass]->SetDirectory(nullptr);
  }
  EHadVis[pclass]->Fill(GetEHadVis(evt) * ToGeV, pKE, w);

  if (!(Transparency.count(pclass) &&
        Transparency[pclass].count(primpart->pid()))) {
    Transparency[pclass][primpart->pid()] = TransparencyFact(pclass);
  }

  if (fsclass == (pclass)) {
    HepMC3::ConstGenParticlePtr fspart;
    for (auto &part : evt.particles()) {
      if (part->status() != NuHepMC::ParticleStatus::UndecayedPhysical) {
        continue;
      }
      switch (pclass) {
      case k1p_only: {
        if (part->pid() == 2212) {
          fspart = part;
        }
        break;
      }
      case k1n_only: {
        if (part->pid() == 2112) {
          fspart = part;
        }
        break;
      }
      case k1pi0_1p: {
        if (part->pid() == 111) {
          fspart = part;
        }
        break;
      }
      }
    }

    auto fs_mom = fspart->momentum();
    auto prim_mom = primpart->momentum();

    double costheta = (fs_mom.x() * prim_mom.x() + fs_mom.y() * prim_mom.y() +
                       fs_mom.z() * prim_mom.z()) /
                      (fs_mom.length() * prim_mom.length());

    double theta = std::acos((costheta > 1) ? 1 : costheta) * 180.0 / M_PI;

    if (theta < 5) {
      Transparency[pclass][primpart->pid()].first->Fill(pKE, w);
    }
  } // end if topo stayed the same

  Transparency[pclass][primpart->pid()].second->Fill(pKE, w);

  for (size_t i = 0; i < pids.size(); ++i) {
    std::stable_sort(KEs[i].rbegin(), KEs[i].rend());

    if (!(Multiplicity.count(pclass) && Multiplicity[pclass].count(pids[i]))) {
      Multiplicity[pclass][pids[i]] = MultHistFact(pclass, pids[i]);
      Multiplicity[pclass][pids[i]]->SetDirectory(nullptr);
    }
    Multiplicity[pclass][pids[i]]->Fill(Mult[i], pKE, w);

    for (size_t j = 0; j < 3; ++j) {
      if (KEs[i].size() <= j) {
        break;
      }

      if (!(KE.count(pclass) && KE[pclass].count(pids[i]) &&
            (KE[pclass][pids[i]].size() >= (j + 1)))) {
        KE[pclass][pids[i]].resize(j + 1);
        KE[pclass][pids[i]][j] = KEHistFact(pclass, pids[i], j);
        KE[pclass][pids[i]][j]->SetDirectory(nullptr);
      }

      KE[pclass][pids[i]][j]->Fill(KEs[i][j], pKE, w);
    }
  }
}

void DumpH(TH1D const &h) {
  std::cout << "# name: " << h.GetName() << std::endl;
  std::cout << "# xtitle: " << h.GetXaxis()->GetTitle() << std::endl;
  std::cout << "# binlow, binup, bincontent, binerror" << std::endl;
  for (int i = 0; i < h.GetXaxis()->GetNbins(); ++i) {
    std::cout << h.GetXaxis()->GetBinLowEdge(i + 1) << ", "
              << h.GetXaxis()->GetBinUpEdge(i + 1) << ", "
              << h.GetBinContent(i + 1) << ", " << h.GetBinError(i + 1)
              << std::endl;
  }
  std::cout << std::endl;
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
    ss << Classification(pclasses[i]);
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

  for (auto &a : KE) {
    std::string name = EHadVis[a.first]->GetName();
    std::cout << "hist: " << name << " " << double(EHadVis[a.first]->GetEntries()) << "/"
              << double(primtopocount[a.first]) << " for this topology = "
              << double(EHadVis[a.first]->GetEntries() * 100) / double(primtopocount[a.first])
              << "%" << std::endl;
    if (double(EHadVis[a.first]->GetEntries()) / double(primtopocount[a.first]) <
        0.05) {
      continue;
    }
    RowNormTH2(EHadVis[a.first].get());
    dout->WriteObject(EHadVis[a.first].release(), name.c_str());

    for (auto &b : a.second) {

      if (Transparency[a.first][b.first].first) {

        Transparency[a.first][b.first].first->Divide(
            Transparency[a.first][b.first].second.get());
        std::string name = Transparency[a.first][b.first].first->GetName();
        dout->WriteObject(Transparency[a.first][b.first].first.release(),
                          name.c_str());
        name = Transparency[a.first][b.first].second->GetName();
        dout->WriteObject(Transparency[a.first][b.first].second.release(),
                          name.c_str());
      }

      if (double(Multiplicity[a.first][b.first]->GetEntries()) /
              double(primtopocount[a.first]) >
          0.05) {
        name = Multiplicity[a.first][b.first]->GetName();
        RowNormTH2(Multiplicity[a.first][b.first].get());
        dout->WriteObject(Multiplicity[a.first][b.first].release(),
                          name.c_str());
      }

      for (auto &c : b.second) {

        name = c->GetName();
        std::cout << "hist: " << name << " " << double(c->GetEntries()) << "/"
                  << double(primtopocount[a.first]) << " for this topology = "
                  << double(c->GetEntries() * 100) /
                         double(primtopocount[a.first])
                  << "%" << std::endl;
        if (double(c->GetEntries()) / double(primtopocount[a.first]) < 0.05) {
          continue;
        }
        RowNormTH2(c.get());
        dout->WriteObject(c.release(), name.c_str());
      }
    }
  }
}