#include "commonana.hxx"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"

#include "fmt/core.h"

#include <iostream>
#include <sstream>

void DumpH(std::ostream &os, TH1D const *h, std::string hname) {
  if (!h) {
    os << fmt::format(R"("{}": "{}-NOTFOUND")", hname, hname) << std::endl;
    return;
  }
  os << fmt::format(R"("{}": {{ )", hname) << std::endl;
  os << fmt::format(R"(  "name": "{}",)", h->GetName()) << std::endl;
  os << fmt::format(R"(  "xtitle": "{}",)", h->GetXaxis()->GetTitle())
     << std::endl;
  os << fmt::format(R"(  "ytitle": "{}",)", h->GetYaxis()->GetTitle())
     << std::endl;
  os << fmt::format(R"(  "v": """)", hname) << std::endl;
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

struct THBlob {
  std::unique_ptr<TH2D> TotalNeutronKE_1p_only;
  std::unique_ptr<TH2D> TotalNeutralE_1p_only;
  std::unique_ptr<TH2D> TotalPi0E_1piplus_1p;
  std::unique_ptr<TH2D> TotalNeutralE_1piplus_1p;

  std::unique_ptr<TH1D> Transparency_1p_only;
  std::unique_ptr<TH1D> Transparency_1piplus_1p;

  std::unique_ptr<TH1D> TotalNeutronKE_1p_only_below200;
  std::unique_ptr<TH1D> TotalNeutronKE_1p_only_above200;
  std::unique_ptr<TH1D> TotalNeutralE_1p_only_below200;
  std::unique_ptr<TH1D> TotalNeutralE_1p_only_above200;

  std::unique_ptr<TH1D> TotalPi0E_1piplus_1p_below300;
  std::unique_ptr<TH1D> TotalPi0E_1piplus_1p_above300;
  std::unique_ptr<TH1D> TotalNeutralE_1piplus_1p_below300;
  std::unique_ptr<TH1D> TotalNeutralE_1piplus_1p_above300;

  THBlob(std::string const &fin_name) {
    TFile fin(fin_name.c_str(), "READ");

    TotalNeutronKE_1p_only =
        std::unique_ptr<TH2D>(fin.Get<TH2D>("TotalNeutronKE_1p_only_nozero"));
    TotalNeutronKE_1p_only->SetDirectory(nullptr);
    TotalNeutralE_1p_only =
        std::unique_ptr<TH2D>(fin.Get<TH2D>("TotalNeutralE_1p_only_nozero"));
    TotalNeutralE_1p_only->SetDirectory(nullptr);
    TotalPi0E_1piplus_1p =
        std::unique_ptr<TH2D>(fin.Get<TH2D>("TotalPi0E_1piplus_1p_nozero"));
    TotalPi0E_1piplus_1p->SetDirectory(nullptr);
    TotalNeutralE_1piplus_1p =
        std::unique_ptr<TH2D>(fin.Get<TH2D>("TotalNeutralE_1piplus_1p_nozero"));
    TotalNeutralE_1piplus_1p->SetDirectory(nullptr);

    Transparency_1p_only =
        std::unique_ptr<TH1D>(fin.Get<TH1D>("k1p_only_proton_transp"));
    Transparency_1p_only->SetDirectory(nullptr);
    Transparency_1piplus_1p =
        std::unique_ptr<TH1D>(fin.Get<TH1D>("k1piplus_1p_piplus_transp"));
    Transparency_1piplus_1p->SetDirectory(nullptr);

    YBin();
  }

  void Write(TDirectory *dir) {
    dir->WriteObject(TotalNeutronKE_1p_only->Clone(), "TotalNeutronKE_1p_only");
    dir->WriteObject(TotalNeutralE_1p_only->Clone(), "TotalNeutralE_1p_only");
    dir->WriteObject(TotalPi0E_1piplus_1p->Clone(), "TotalPi0E_1piplus_1p");
    dir->WriteObject(TotalNeutralE_1piplus_1p->Clone(),
                     "TotalNeutralE_1piplus_1p");

    dir->WriteObject(Transparency_1p_only->Clone(), "Transparency_1p_only");
    dir->WriteObject(Transparency_1piplus_1p->Clone(),
                     "Transparency_1piplus_1p");

    dir->WriteObject(TotalNeutronKE_1p_only_below200->Clone(),
                     "TotalNeutronKE_1p_only_below200");
    dir->WriteObject(TotalNeutronKE_1p_only_above200->Clone(),
                     "TotalNeutronKE_1p_only_above200");
    dir->WriteObject(TotalNeutralE_1p_only_below200->Clone(),
                     "TotalNeutralE_1p_only_below200");
    dir->WriteObject(TotalNeutralE_1p_only_above200->Clone(),
                     "TotalNeutralE_1p_only_above200");

    dir->WriteObject(TotalPi0E_1piplus_1p_below300->Clone(),
                     "TotalPi0E_1piplus_1p_below300");
    dir->WriteObject(TotalPi0E_1piplus_1p_above300->Clone(),
                     "TotalPi0E_1piplus_1p_above300");
    dir->WriteObject(TotalNeutralE_1piplus_1p_below300->Clone(),
                     "TotalNeutralE_1piplus_1p_below300");
    dir->WriteObject(TotalNeutralE_1piplus_1p_above300->Clone(),
                     "TotalNeutralE_1piplus_1p_above300");
  }

  void YBin() {

    double rebin = 4;

    TotalNeutronKE_1p_only_below200 =
        std::unique_ptr<TH1D>(TotalNeutronKE_1p_only->ProjectionX(
            "TotalNeutronKE_1p_only_below200", 1, 10));
    TotalNeutronKE_1p_only_below200->Scale(100.0 / 10.0, "");
    TotalNeutronKE_1p_only_below200->SetDirectory(nullptr);

    TotalNeutronKE_1p_only_above200 =
        std::unique_ptr<TH1D>(TotalNeutronKE_1p_only->ProjectionX(
            "TotalNeutronKE_1p_only_above200", 11, 50));
    TotalNeutronKE_1p_only_above200->RebinX(rebin);
    TotalNeutronKE_1p_only_above200->Scale(100.0 / (40.0), "");
    TotalNeutronKE_1p_only_above200->SetDirectory(nullptr);

    TotalNeutralE_1p_only_below200 =
        std::unique_ptr<TH1D>(TotalNeutralE_1p_only->ProjectionX(
            "TotalNeutralE_1p_only_below200", 1, 10));
    TotalNeutralE_1p_only_below200->Scale(100.0 / (10.0)), "";
    TotalNeutralE_1p_only_below200->SetDirectory(nullptr);

    TotalNeutralE_1p_only_above200 =
        std::unique_ptr<TH1D>(TotalNeutralE_1p_only->ProjectionX(
            "TotalNeutralE_1p_only_above200", 11, 50));
    TotalNeutralE_1p_only_above200->RebinX(rebin);
    TotalNeutralE_1p_only_above200->Scale(100.0 / (40.0), "");
    TotalNeutralE_1p_only_above200->SetDirectory(nullptr);

    TotalPi0E_1piplus_1p_below300 =
        std::unique_ptr<TH1D>(TotalPi0E_1piplus_1p->ProjectionX(
            "TotalPi0E_1piplus_1p_below300", 1, 15));
    TotalPi0E_1piplus_1p_below300->Scale(100.0 / (15.0), "");
    TotalPi0E_1piplus_1p_below300->SetDirectory(nullptr);

    TotalPi0E_1piplus_1p_above300 =
        std::unique_ptr<TH1D>(TotalPi0E_1piplus_1p->ProjectionX(
            "TotalPi0E_1piplus_1p_above300", 16, 50));
    TotalPi0E_1piplus_1p_above300->RebinX(rebin);
    TotalPi0E_1piplus_1p_above300->Scale(100.0 / (35.0), "");
    TotalPi0E_1piplus_1p_above300->SetDirectory(nullptr);

    TotalNeutralE_1piplus_1p_below300 =
        std::unique_ptr<TH1D>(TotalNeutralE_1piplus_1p->ProjectionX(
            "TotalNeutralE_1piplus_1p_below300", 1, 15));
    TotalNeutralE_1piplus_1p_below300->Scale(100.0 / (15.0), "");
    TotalNeutralE_1piplus_1p_below300->SetDirectory(nullptr);

    TotalNeutralE_1piplus_1p_above300 =
        std::unique_ptr<TH1D>(TotalNeutralE_1piplus_1p->ProjectionX(
            "TotalNeutralE_1piplus_1p_above300", 16, 50));
    TotalNeutralE_1piplus_1p_above300->RebinX(rebin);
    TotalNeutralE_1piplus_1p_above300->Scale(100.0 / (35.0), "");
    TotalNeutralE_1piplus_1p_above300->SetDirectory(nullptr);
  }
};

int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cout << "[RUNLIKE]: " << argv[0] << " <infile.hepmc3>" << std::endl;
    return 1;
  }

  std::vector<std::pair<std::string, THBlob>> infs;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];

    auto colpos = arg.find_first_of(':');
    std::string name = arg.substr(0, colpos);
    std::string fname = arg.substr(colpos + 1);

    infs.emplace_back(name, THBlob(fname));
    // std::cout << " Read " << name << " histos from " << fname << std::endl;
  }

  // TFile fout("dumpedhists.root", "RECREATE");

  // for (auto &hb : infs) {
  //   auto gdir = fout.mkdir(hb.first.c_str());

  //   std::ostream &os = std::cout;

  //   os << hb.first << " = { " << std::endl;
  //   os << fmt::format(R"("label": "{}",)", hb.first) << std::endl;

  //   DumpH(os, hb.second.Transparency_1p_only.get(), "transp_1p");
  //   os << "," << std::endl;
  //   DumpH(os, hb.second.Transparency_1piplus_1p.get(), "transp_1piplus_1p");
  //   os << "," << std::endl;
  //   DumpH(os, hb.second.TotalNeutronKE_1p_only_below200.get(),
  //         "neutronKE_1p_below200MeV");
  //   os << "," << std::endl;
  //   DumpH(os, hb.second.TotalNeutronKE_1p_only_above200.get(),
  //         "neutronKE_1p_above200MeV");
  //   os << "," << std::endl;
  //   DumpH(os, hb.second.TotalNeutralE_1p_only_below200.get(),
  //         "neutralE_1p_below200MeV");
  //   os << "," << std::endl;
  //   DumpH(os, hb.second.TotalNeutralE_1p_only_above200.get(),
  //         "neutralE_1p_above200MeV");
  //   os << "," << std::endl;
  //   DumpH(os, hb.second.TotalPi0E_1piplus_1p_below300.get(),
  //         "pi0E_1piplus_1p_below300MeV");
  //   os << "," << std::endl;
  //   DumpH(os, hb.second.TotalPi0E_1piplus_1p_above300.get(),
  //         "pi0E_1piplus_1p_above300MeV");
  //   os << "," << std::endl;
  //   DumpH(os, hb.second.TotalNeutralE_1piplus_1p_below300.get(),
  //         "neutralE_1piplus_1p_below300MeV");
  //   os << "," << std::endl;
  //   DumpH(os, hb.second.TotalNeutralE_1piplus_1p_above300.get(),
  //         "neutralE_1piplus_1p_above300MeV");
  //   os << "," << std::endl;

  //   os << "} " << std::endl;

  //   hb.second.Write(gdir);
  // }
  double fontsize = 0.075;

  TLatex ltx;
  ltx.SetTextSize(fontsize);
  ltx.SetTextFont(132);

  gStyle->SetOptStat(false);
  { // 1p->neutron
    TCanvas c1("c1", "", 1400, 800);

    TPad pleft("pleft", "", 0, 0, 0.5, 1);
    pleft.AppendPad();
    pleft.SetLeftMargin(0.28);
    pleft.SetRightMargin(0.03);
    pleft.SetTopMargin(0.03);
    pleft.SetBottomMargin(0.25);
    TPad pright("pright", "", 0.5, 0, 1, 1);
    pright.AppendPad();
    pright.SetLeftMargin(0.05);
    pright.SetRightMargin(0.26);
    pright.SetTopMargin(0.03);
    pright.SetBottomMargin(0.25);

    pleft.cd();
    int first = 0;

    TLegend *legendl = new TLegend(0.25, 0.5, 0.425, 0.8);
    legendl->SetTextFont(132);
    legendl->SetTextSize(fontsize);
    legendl->SetBorderSize(0);
    legendl->SetFillStyle(0);
    legendl->SetNColumns(1);

    double maxy = 12;

    int cols[] = {TColor::GetColor("#0077bb"), TColor::GetColor("#009988"),
                  TColor::GetColor("#ee7733")};

    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;
      if (hb.TotalNeutronKE_1p_only_below200->Integral() < 1E-8) {
        continue;
      }

      hb.TotalNeutronKE_1p_only_below200->GetYaxis()->SetRangeUser(0, maxy);
      hb.TotalNeutronKE_1p_only_below200->GetYaxis()->SetLabelSize(fontsize);
      hb.TotalNeutronKE_1p_only_below200->GetYaxis()->SetTitleSize(fontsize);
      hb.TotalNeutronKE_1p_only_below200->GetYaxis()->SetTitleFont(132);
      hb.TotalNeutronKE_1p_only_below200->GetYaxis()->SetLabelFont(132);
      hb.TotalNeutronKE_1p_only_below200->GetYaxis()->SetTitle(
          "P.D.F. x10^{2}");
      hb.TotalNeutronKE_1p_only_below200->GetYaxis()->SetNdivisions(505);

      hb.TotalNeutronKE_1p_only_below200->GetXaxis()->SetRangeUser(0, 0.2);
      hb.TotalNeutronKE_1p_only_below200->GetXaxis()->SetLabelSize(fontsize);
      hb.TotalNeutronKE_1p_only_below200->GetXaxis()->SetTitleSize(fontsize);
      hb.TotalNeutronKE_1p_only_below200->GetXaxis()->SetLabelFont(132);
      hb.TotalNeutronKE_1p_only_below200->GetXaxis()->SetTitleOffset(1);
      hb.TotalNeutronKE_1p_only_below200->GetXaxis()->SetNdivisions(505);
      hb.TotalNeutronKE_1p_only_below200->GetXaxis()->SetTitle("");

      hb.TotalNeutronKE_1p_only_below200->SetLineColor(cols[g_it]);
      hb.TotalNeutronKE_1p_only_below200->SetLineWidth(2);
      hb.TotalNeutronKE_1p_only_below200->Draw((!first++) ? "HIST"
                                                          : "HIST SAME");

      legendl->AddEntry(hb.TotalNeutronKE_1p_only_below200.get(),
                        infs[g_it].first.c_str(), "l");
    }

    pright.cd();
    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;
      if (hb.TotalNeutronKE_1p_only_above200->Integral() < 1E-8) {
        continue;
      }
      hb.TotalNeutronKE_1p_only_above200->GetYaxis()->SetRangeUser(0, maxy);
      hb.TotalNeutronKE_1p_only_above200->GetYaxis()->SetLabelOffset(1);
      hb.TotalNeutronKE_1p_only_above200->GetYaxis()->SetLabelSize(fontsize);
      hb.TotalNeutronKE_1p_only_above200->GetYaxis()->SetTitleSize(fontsize);
      hb.TotalNeutronKE_1p_only_above200->GetYaxis()->SetTitleFont(132);
      hb.TotalNeutronKE_1p_only_above200->GetYaxis()->SetLabelFont(132);
      hb.TotalNeutronKE_1p_only_above200->GetYaxis()->SetNdivisions(505);

      hb.TotalNeutronKE_1p_only_above200->GetXaxis()->SetRangeUser(0, 1);
      hb.TotalNeutronKE_1p_only_above200->GetXaxis()->SetLabelSize(fontsize);
      hb.TotalNeutronKE_1p_only_above200->GetXaxis()->SetTitleSize(fontsize);
      hb.TotalNeutronKE_1p_only_above200->GetXaxis()->SetLabelFont(132);
      hb.TotalNeutronKE_1p_only_above200->GetXaxis()->SetTitleOffset(1);
      hb.TotalNeutronKE_1p_only_above200->GetXaxis()->SetNdivisions(505);
      hb.TotalNeutronKE_1p_only_above200->GetXaxis()->SetTitle("");

      hb.TotalNeutronKE_1p_only_above200->SetLineColor(cols[g_it]);
      hb.TotalNeutronKE_1p_only_above200->SetLineWidth(2);
      hb.TotalNeutronKE_1p_only_above200->Draw((!first++) ? "HIST"
                                                          : "HIST SAME");
    }
    c1.cd();
    ltx.SetTextAlign(32);
    ltx.DrawLatexNDC(0.39, 0.08, "#sum");
    ltx.SetTextAlign(12);
    ltx.DrawLatexNDC(0.41, 0.08, "#it{T}_{neutron} (GeV)");

    ltx.DrawLatexNDC(0.2, 0.9, "#it{T}_{p}^{prim.} < 0.2 GeV");
    ltx.DrawLatexNDC(0.55, 0.9, "0.2 < #it{T}_{p}^{prim.} < 1 GeV");

    legendl->Draw();
    c1.Print("TotalNeutronKE_1p.pdf");
    c1.Print("allplots.pdf[");
    c1.Print("allplots.pdf");
  }

  { // 1p->neutron
    TCanvas c1("c1", "", 1400, 800);

    TPad pleft("pleft", "", 0, 0, 0.5, 1);
    pleft.AppendPad();
    pleft.SetLeftMargin(0.28);
    pleft.SetRightMargin(0.03);
    pleft.SetTopMargin(0.03);
    pleft.SetBottomMargin(0.25);
    TPad pright("pright", "", 0.5, 0, 1, 1);
    pright.AppendPad();
    pright.SetLeftMargin(0.05);
    pright.SetRightMargin(0.26);
    pright.SetTopMargin(0.03);
    pright.SetBottomMargin(0.25);

    pleft.cd();
    int first = 0;
    TLegend *legendl = new TLegend(0.25, 0.5, 0.425, 0.8);
    legendl->SetTextFont(132);
    legendl->SetTextSize(fontsize);
    legendl->SetBorderSize(0);
    legendl->SetFillStyle(0);
    legendl->SetNColumns(1);

    double maxy = 12;

    int cols[] = {TColor::GetColor("#0077bb"), TColor::GetColor("#009988"),
                  TColor::GetColor("#ee7733")};

    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;
      if (hb.TotalNeutralE_1p_only_below200->Integral() < 1E-8) {
        continue;
      }

      hb.TotalNeutralE_1p_only_below200->GetYaxis()->SetRangeUser(0, maxy);
      hb.TotalNeutralE_1p_only_below200->GetYaxis()->SetLabelSize(fontsize);
      hb.TotalNeutralE_1p_only_below200->GetYaxis()->SetTitleSize(fontsize);
      hb.TotalNeutralE_1p_only_below200->GetYaxis()->SetTitleFont(132);
      hb.TotalNeutralE_1p_only_below200->GetYaxis()->SetLabelFont(132);
      hb.TotalNeutralE_1p_only_below200->GetYaxis()->SetTitle("P.D.F. x10^{2}");
      hb.TotalNeutralE_1p_only_below200->GetYaxis()->SetNdivisions(505);

      hb.TotalNeutralE_1p_only_below200->GetXaxis()->SetRangeUser(0, 0.2);
      hb.TotalNeutralE_1p_only_below200->GetXaxis()->SetLabelSize(fontsize);
      hb.TotalNeutralE_1p_only_below200->GetXaxis()->SetTitleSize(fontsize);
      hb.TotalNeutralE_1p_only_below200->GetXaxis()->SetLabelFont(132);
      hb.TotalNeutralE_1p_only_below200->GetXaxis()->SetTitleOffset(1);
      hb.TotalNeutralE_1p_only_below200->GetXaxis()->SetNdivisions(505);
      hb.TotalNeutralE_1p_only_below200->GetXaxis()->SetTitle("");

      hb.TotalNeutralE_1p_only_below200->SetLineColor(cols[g_it]);
      hb.TotalNeutralE_1p_only_below200->SetLineWidth(2);
      hb.TotalNeutralE_1p_only_below200->Draw((!first++) ? "HIST"
                                                         : "HIST SAME");
      legendl->AddEntry(hb.TotalNeutralE_1p_only_below200.get(),
                        infs[g_it].first.c_str(), "l");
    }

    pright.cd();
    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;
      if (hb.TotalNeutralE_1p_only_above200->Integral() < 1E-8) {
        continue;
      }
      hb.TotalNeutralE_1p_only_above200->GetYaxis()->SetRangeUser(0, maxy);
      hb.TotalNeutralE_1p_only_above200->GetYaxis()->SetLabelOffset(1);
      hb.TotalNeutralE_1p_only_above200->GetYaxis()->SetLabelSize(fontsize);
      hb.TotalNeutralE_1p_only_above200->GetYaxis()->SetTitleSize(fontsize);
      hb.TotalNeutralE_1p_only_above200->GetYaxis()->SetTitleFont(132);
      hb.TotalNeutralE_1p_only_above200->GetYaxis()->SetLabelFont(132);
      hb.TotalNeutralE_1p_only_above200->GetYaxis()->SetNdivisions(505);

      hb.TotalNeutralE_1p_only_above200->GetXaxis()->SetRangeUser(0, 1);
      hb.TotalNeutralE_1p_only_above200->GetXaxis()->SetLabelSize(fontsize);
      hb.TotalNeutralE_1p_only_above200->GetXaxis()->SetTitleSize(fontsize);
      hb.TotalNeutralE_1p_only_above200->GetXaxis()->SetLabelFont(132);
      hb.TotalNeutralE_1p_only_above200->GetXaxis()->SetTitleOffset(1);
      hb.TotalNeutralE_1p_only_above200->GetXaxis()->SetNdivisions(505);
      hb.TotalNeutralE_1p_only_above200->GetXaxis()->SetTitle("");

      hb.TotalNeutralE_1p_only_above200->SetLineColor(cols[g_it]);
      hb.TotalNeutralE_1p_only_above200->SetLineWidth(2);
      hb.TotalNeutralE_1p_only_above200->Draw((!first++) ? "HIST"
                                                         : "HIST SAME");
    }
    c1.cd();
    ltx.SetTextAlign(32);
    ltx.DrawLatexNDC(0.39, 0.08, "#sum");
    ltx.SetTextAlign(12);
    ltx.DrawLatexNDC(0.41, 0.08, "#it{T}_{neutron} + #it{E}_{#pi^{0}} (GeV)");

    ltx.DrawLatexNDC(0.2, 0.9, "#it{T}_{p}^{prim.} < 0.2 GeV");
    ltx.DrawLatexNDC(0.55, 0.9, "0.2 < #it{T}_{p}^{prim.} < 1 GeV");

    legendl->Draw();
    c1.Print("TotalNeutralE_1p.pdf");
    c1.Print("allplots.pdf");
  }

  { // 1p->neutron
    TCanvas c1("c1", "", 1400, 800);

    TPad pleft("pleft", "", 0, 0, 0.5, 1);
    pleft.AppendPad();
    pleft.SetLeftMargin(0.28);
    pleft.SetRightMargin(0.03);
    pleft.SetTopMargin(0.03);
    pleft.SetBottomMargin(0.25);
    TPad pright("pright", "", 0.5, 0, 1, 1);
    pright.AppendPad();
    pright.SetLeftMargin(0.05);
    pright.SetRightMargin(0.26);
    pright.SetTopMargin(0.03);
    pright.SetBottomMargin(0.25);

    pleft.cd();
    int first = 0;
    TLegend *legendl = new TLegend(0.25, 0.5, 0.425, 0.8);
    legendl->SetTextFont(132);
    legendl->SetTextSize(fontsize);
    legendl->SetBorderSize(0);
    legendl->SetFillStyle(0);
    legendl->SetNColumns(1);

    double maxy = 3;

    int cols[] = {TColor::GetColor("#0077bb"), TColor::GetColor("#009988"),
                  TColor::GetColor("#ee7733")};

    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;
      if (hb.TotalPi0E_1piplus_1p_below300->Integral() < 1E-8) {
        continue;
      }

      hb.TotalPi0E_1piplus_1p_below300->GetYaxis()->SetRangeUser(0, maxy);
      hb.TotalPi0E_1piplus_1p_below300->GetYaxis()->SetLabelSize(fontsize);
      hb.TotalPi0E_1piplus_1p_below300->GetYaxis()->SetTitleSize(fontsize);
      hb.TotalPi0E_1piplus_1p_below300->GetYaxis()->SetTitleFont(132);
      hb.TotalPi0E_1piplus_1p_below300->GetYaxis()->SetLabelFont(132);
      hb.TotalPi0E_1piplus_1p_below300->GetYaxis()->SetTitle("P.D.F. x10^{2}");
      hb.TotalPi0E_1piplus_1p_below300->GetYaxis()->SetNdivisions(505);

      hb.TotalPi0E_1piplus_1p_below300->GetXaxis()->SetRangeUser(0.1, 0.3);
      hb.TotalPi0E_1piplus_1p_below300->GetXaxis()->SetLabelSize(fontsize);
      hb.TotalPi0E_1piplus_1p_below300->GetXaxis()->SetTitleSize(fontsize);
      hb.TotalPi0E_1piplus_1p_below300->GetXaxis()->SetLabelFont(132);
      hb.TotalPi0E_1piplus_1p_below300->GetXaxis()->SetTitleOffset(1);
      hb.TotalPi0E_1piplus_1p_below300->GetXaxis()->SetNdivisions(505);
      hb.TotalPi0E_1piplus_1p_below300->GetXaxis()->SetTitle("");

      hb.TotalPi0E_1piplus_1p_below300->SetLineColor(cols[g_it]);
      hb.TotalPi0E_1piplus_1p_below300->SetLineWidth(2);
      hb.TotalPi0E_1piplus_1p_below300->Draw((!first++) ? "HIST" : "HIST SAME");
      legendl->AddEntry(hb.TotalPi0E_1piplus_1p_below300.get(),
                        infs[g_it].first.c_str(), "l");
    }

    pright.cd();
    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;
      if (hb.TotalPi0E_1piplus_1p_above300->Integral() < 1E-8) {
        continue;
      }
      hb.TotalPi0E_1piplus_1p_above300->GetYaxis()->SetRangeUser(0, maxy);
      hb.TotalPi0E_1piplus_1p_above300->GetYaxis()->SetLabelOffset(1);
      hb.TotalPi0E_1piplus_1p_above300->GetYaxis()->SetLabelSize(fontsize);
      hb.TotalPi0E_1piplus_1p_above300->GetYaxis()->SetTitleSize(fontsize);
      hb.TotalPi0E_1piplus_1p_above300->GetYaxis()->SetTitleFont(132);
      hb.TotalPi0E_1piplus_1p_above300->GetYaxis()->SetLabelFont(132);
      hb.TotalPi0E_1piplus_1p_above300->GetYaxis()->SetNdivisions(505);

      hb.TotalPi0E_1piplus_1p_above300->GetXaxis()->SetRangeUser(0, 1);
      hb.TotalPi0E_1piplus_1p_above300->GetXaxis()->SetLabelSize(fontsize);
      hb.TotalPi0E_1piplus_1p_above300->GetXaxis()->SetTitleSize(fontsize);
      hb.TotalPi0E_1piplus_1p_above300->GetXaxis()->SetLabelFont(132);
      hb.TotalPi0E_1piplus_1p_above300->GetXaxis()->SetTitleOffset(1);
      hb.TotalPi0E_1piplus_1p_above300->GetXaxis()->SetNdivisions(505);
      hb.TotalPi0E_1piplus_1p_above300->GetXaxis()->SetTitle("");

      hb.TotalPi0E_1piplus_1p_above300->SetLineColor(cols[g_it]);
      hb.TotalPi0E_1piplus_1p_above300->SetLineWidth(2);
      hb.TotalPi0E_1piplus_1p_above300->Draw((!first++) ? "HIST" : "HIST SAME");
    }
    c1.cd();
    ltx.SetTextAlign(32);
    ltx.DrawLatexNDC(0.39, 0.08, "#sum");
    ltx.SetTextAlign(12);
    ltx.DrawLatexNDC(0.41, 0.08, "#it{E}_{#pi^{0}} (GeV)");

    ltx.DrawLatexNDC(0.2, 0.9, "#it{T}_{#pi^{+}}^{prim.} < 0.3 GeV");
    ltx.DrawLatexNDC(0.55, 0.9, "0.3 < #it{T}_{#pi^{+}}^{prim.} < 1 GeV");

    legendl->Draw();
    c1.Print("TotalPi0E_1piplus_1p.pdf");
    c1.Print("allplots.pdf");
  }

  { // 1p->neutron
    TCanvas c1("c1", "", 1400, 800);

    TPad pleft("pleft", "", 0, 0, 0.5, 1);
    pleft.AppendPad();
    pleft.SetLeftMargin(0.28);
    pleft.SetRightMargin(0.03);
    pleft.SetTopMargin(0.03);
    pleft.SetBottomMargin(0.25);
    TPad pright("pright", "", 0.5, 0, 1, 1);
    pright.AppendPad();
    pright.SetLeftMargin(0.05);
    pright.SetRightMargin(0.26);
    pright.SetTopMargin(0.03);
    pright.SetBottomMargin(0.25);

    pleft.cd();
    int first = 0;
    TLegend *legendl = new TLegend(0.25, 0.5, 0.425, 0.8);
    legendl->SetTextFont(132);
    legendl->SetTextSize(fontsize);
    legendl->SetBorderSize(0);
    legendl->SetFillStyle(0);
    legendl->SetNColumns(1);

    double maxy = 10;

    int cols[] = {TColor::GetColor("#0077bb"), TColor::GetColor("#009988"),
                  TColor::GetColor("#ee7733")};

    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;
      if (hb.TotalNeutralE_1piplus_1p_below300->Integral() < 1E-8) {
        continue;
      }

      hb.TotalNeutralE_1piplus_1p_below300->GetYaxis()->SetRangeUser(0, maxy);
      hb.TotalNeutralE_1piplus_1p_below300->GetYaxis()->SetLabelSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_below300->GetYaxis()->SetTitleSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_below300->GetYaxis()->SetTitleFont(132);
      hb.TotalNeutralE_1piplus_1p_below300->GetYaxis()->SetLabelFont(132);
      hb.TotalNeutralE_1piplus_1p_below300->GetYaxis()->SetTitle(
          "P.D.F. x10^{2}");
      hb.TotalNeutralE_1piplus_1p_below300->GetYaxis()->SetNdivisions(505);

      hb.TotalNeutralE_1piplus_1p_below300->GetXaxis()->SetRangeUser(0, 0.5);
      hb.TotalNeutralE_1piplus_1p_below300->GetXaxis()->SetLabelSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_below300->GetXaxis()->SetTitleSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_below300->GetXaxis()->SetLabelFont(132);
      hb.TotalNeutralE_1piplus_1p_below300->GetXaxis()->SetTitleOffset(1);
      hb.TotalNeutralE_1piplus_1p_below300->GetXaxis()->SetNdivisions(505);
      hb.TotalNeutralE_1piplus_1p_below300->GetXaxis()->SetTitle("");

      hb.TotalNeutralE_1piplus_1p_below300->SetLineColor(cols[g_it]);
      hb.TotalNeutralE_1piplus_1p_below300->SetLineWidth(2);
      hb.TotalNeutralE_1piplus_1p_below300->Draw((!first++) ? "HIST"
                                                            : "HIST SAME");
      legendl->AddEntry(hb.TotalNeutralE_1piplus_1p_below300.get(),
                        infs[g_it].first.c_str(), "l");
    }

    pright.cd();
    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;
      if (hb.TotalNeutralE_1piplus_1p_above300->Integral() < 1E-8) {
        continue;
      }
      hb.TotalNeutralE_1piplus_1p_above300->GetYaxis()->SetRangeUser(0, maxy);
      hb.TotalNeutralE_1piplus_1p_above300->GetYaxis()->SetLabelOffset(1);
      hb.TotalNeutralE_1piplus_1p_above300->GetYaxis()->SetLabelSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_above300->GetYaxis()->SetTitleSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_above300->GetYaxis()->SetTitleFont(132);
      hb.TotalNeutralE_1piplus_1p_above300->GetYaxis()->SetLabelFont(132);
      hb.TotalNeutralE_1piplus_1p_above300->GetYaxis()->SetNdivisions(505);

      hb.TotalNeutralE_1piplus_1p_above300->GetXaxis()->SetRangeUser(0, 1);
      hb.TotalNeutralE_1piplus_1p_above300->GetXaxis()->SetLabelSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_above300->GetXaxis()->SetTitleSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_above300->GetXaxis()->SetLabelFont(132);
      hb.TotalNeutralE_1piplus_1p_above300->GetXaxis()->SetTitleOffset(1);
      hb.TotalNeutralE_1piplus_1p_above300->GetXaxis()->SetNdivisions(505);
      hb.TotalNeutralE_1piplus_1p_above300->GetXaxis()->SetTitle("");

      hb.TotalNeutralE_1piplus_1p_above300->SetLineColor(cols[g_it]);
      hb.TotalNeutralE_1piplus_1p_above300->SetLineWidth(2);
      hb.TotalNeutralE_1piplus_1p_above300->Draw((!first++) ? "HIST"
                                                            : "HIST SAME");
    }
    c1.cd();
    ltx.SetTextAlign(32);
    ltx.DrawLatexNDC(0.39, 0.08, "#sum");
    ltx.SetTextAlign(12);
    ltx.DrawLatexNDC(0.41, 0.08, "#it{T}_{neutron} + #it{E}_{#pi^{0}} (GeV)");

    ltx.DrawLatexNDC(0.2, 0.9, "#it{T}_{#pi^{+}}^{prim.} < 0.3 GeV");
    ltx.DrawLatexNDC(0.55, 0.9, "0.3 < #it{T}_{#pi^{+}}^{prim.} < 1 GeV");

    legendl->Draw();
    c1.Print("TotalNeutralE_1piplus_1p.pdf");
    c1.Print("allplots.pdf");
    c1.Print("allplots.pdf]");
  }
}