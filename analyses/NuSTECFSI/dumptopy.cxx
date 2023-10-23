#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include "fmt/core.h"

#include <iostream>
#include <sstream>

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

std::string KEHistName(Classification c, int pdg, int l) {

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

  return to_string(c) + "_" + ld + "_" + pdgs + "_KE";
}

std::string MultHistName(Classification c, int pdg) {

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

  return to_string(c) + "_" + pdgs + "_mult";
}

std::pair<std::string, std::string> TransparencyName(Classification c) {

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
  default:
    throw;
  }

  return {to_string(c) + "_" + primpart + "_transp",
          to_string(c) + "_" + primpart + "_all"};
}

std::string EHadVisHistName(Classification c) {

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

  return to_string(c) + "_EHadVis";
}

void DumpH(std::ostream &os, TH1D const *h, std::string key,
           std::string hname) {
  if (!h) {
    os << fmt::format(R"("{}": "{}-NOTFOUND")", key, hname) << std::endl;
    return;
  }
  os << fmt::format(R"("{}": {{ )", key) << std::endl;
  os << fmt::format(R"(  "name": "{}",)", h->GetName()) << std::endl;
  os << fmt::format(R"(  "xtitle": "{}",)", h->GetXaxis()->GetTitle())
     << std::endl;
  os << fmt::format(R"(  "ytitle": "{}",)", h->GetYaxis()->GetTitle())
     << std::endl;
  os << fmt::format(R"(  "v": """)", key) << std::endl;
  os << "    # binlow, binup, bincontent, binerror\n    #----------"
     << std::endl;
  for (int i = 0; i < h->GetXaxis()->GetNbins(); ++i) {
    os << fmt::format(R"(      {}, {}, {}, {})",
                      h->GetXaxis()->GetBinLowEdge(i + 1),
                      h->GetXaxis()->GetBinUpEdge(i + 1),
                      h->GetBinContent(i + 1), h->GetBinError(i + 1))
       << std::endl;
  }
  os << fmt::format(R"(    #----------
    """,
}})") << std::flush;
}

void DumpH(std::ostream &os, TH2D const *h, std::string key,
           std::string hname) {
  if (!h) {
    os << fmt::format(R"("{}": "{}-NOTFOUND")", key, hname) << std::endl;
    return;
  }
  os << fmt::format(R"("{}": [)", key) << std::endl;

  for (int j = 0; j < h->GetYaxis()->GetNbins(); ++j) {
    os << fmt::format(R"(  {{ 
    "ybin": ({},{}), )",
                      h->GetYaxis()->GetBinLowEdge(j + 1),
                      h->GetYaxis()->GetBinUpEdge(j + 1))
       << std::endl;
    os << fmt::format(R"(    "name": "{}",)", h->GetName()) << std::endl;
    os << fmt::format(R"(    "xtitle": "{}",)", h->GetXaxis()->GetTitle())
       << std::endl;
    os << fmt::format(R"(    "ytitle": "{}",)", h->GetYaxis()->GetTitle())
       << std::endl;
    os << fmt::format(R"(    "v": """)", key) << std::endl;
    os << "      # binlow, binup, bincontent, binerror\n      #----------"
       << std::endl;
    for (int i = 0; i < h->GetXaxis()->GetNbins(); ++i) {
      os << fmt::format(R"(        {}, {}, {}, {})",
                        h->GetXaxis()->GetBinLowEdge(i + 1),
                        h->GetXaxis()->GetBinUpEdge(i + 1),
                        h->GetBinContent(i + 1, j + 1),
                        h->GetBinError(i + 1, j + 1))
         << std::endl;
    }
    os << fmt::format(R"(      #----------
    """
  }},)")
       << std::endl;
  }
  os << "]" << std::flush;
}

void DumpHAN(std::ostream &os, TH2D const *h, std::string key,
             std::string hname) {
  if (!h) {
    os << fmt::format(R"("{}": "{}-NOTFOUND")", key, hname) << std::endl;
    return;
  }
  os << fmt::format(R"("{}": [)", key) << std::endl;

  for (int j = 0; j < h->GetYaxis()->GetNbins(); ++j) {
    os << fmt::format(R"(  {{ 
    "ybin": "{}", )",
                      h->GetYaxis()->GetBinLabel(j + 1))
       << std::endl;
    os << fmt::format(R"(    "name": "{}",)", h->GetName()) << std::endl;
    os << fmt::format(R"(    "xtitle": "{}",)", h->GetXaxis()->GetTitle())
       << std::endl;
    os << fmt::format(R"(    "ytitle": "{}",)", h->GetYaxis()->GetTitle())
       << std::endl;
    os << fmt::format(R"(    "v": """)", key) << std::endl;
    os << "      # binlabel | bincontent | binerror\n      #----------"
       << std::endl;
    for (int i = 0; i < h->GetXaxis()->GetNbins(); ++i) {
      os << fmt::format(
                R"(        {} | {} | {})", h->GetXaxis()->GetBinLabel(i + 1),
                h->GetBinContent(i + 1, j + 1), h->GetBinError(i + 1, j + 1))
         << std::endl;
    }
    os << fmt::format(R"(      #----------
    """
  }},)")
       << std::endl;
  }
  os << "]" << std::flush;
}

int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cout << "[RUNLIKE]: " << argv[0] << " <infile.hepmc3>" << std::endl;
    return 1;
  }

  std::string inf = argv[1];
  std::string dir = (argc > 2) ? argv[2] : "";

  TFile fin(inf.c_str(), "READ");

  TDirectory *dout = &fin;

  std::ostream &os = std::cout;

  os << dir << " = { " << std::endl;
  os << fmt::format(R"("label": "{}",)", dir) << std::endl;

  DumpHAN(os, dout->Get<TH2D>("PrimaryToFinalStateSmearing"),
          "PrimaryToFinalStateSmearing", "PrimaryToFinalStateSmearing");
  os << "," << std::endl;

  DumpH(os, dout->Get<TH1D>(TransparencyName(k1p_only).first.c_str()),
        "proton_transp", TransparencyName(k1p_only).first);
  os << "," << std::endl;
  DumpH(os, dout->Get<TH1D>(TransparencyName(k1n_only).first.c_str()),
        "neutron_transp", TransparencyName(k1n_only).first);
  os << "," << std::endl;
  DumpH(os, dout->Get<TH1D>(TransparencyName(k1pi0_1p).first.c_str()),
        "pi0_transp", TransparencyName(k1pi0_1p).first);
  os << "," << std::endl;

  for (auto s : {k1p_only, k1n_only}) {
    for (auto p : {2212, 2112}) {
      DumpH(os, dout->Get<TH2D>(MultHistName(s, p).c_str()), MultHistName(s, p),
            MultHistName(s, p));
      os << "," << std::endl;
    }
  }
  DumpH(os, dout->Get<TH2D>(MultHistName(k1pi0_1p, 111).c_str()),
        MultHistName(k1pi0_1p, 111), MultHistName(k1pi0_1p, 111));
  os << "," << std::endl;

  for (auto s : {k1p_only, k1n_only}) {
    for (auto p : {2212, 2112}) {
      for (auto l : {0, 1, 2}) {
        DumpH(os, dout->Get<TH2D>(KEHistName(s, p, l).c_str()),
              KEHistName(s, p, l), KEHistName(s, p, l));
        os << "," << std::endl;
      }
    }
  }
  DumpH(os, dout->Get<TH2D>(KEHistName(k1pi0_1p, 111, 0).c_str()),
        KEHistName(k1pi0_1p, 111, 0), KEHistName(k1pi0_1p, 111, 0));
  os << "," << std::endl;
  for (auto s : {k1p_only, k1n_only, k1pi0_1p}) {
    DumpH(os, dout->Get<TH2D>(EHadVisHistName(s).c_str()), EHadVisHistName(s),
          EHadVisHistName(s));
    os << "," << std::endl;
  }

  os << "} " << std::endl;
}