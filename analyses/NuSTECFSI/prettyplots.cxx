#include "commonana.hxx"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"

#include "fmt/core.h"

#include <iostream>
#include <sstream>

bool reshape = true;

struct THBlob {
  std::unique_ptr<TH2D> PrimaryToFinalStateSmearing;

  std::unique_ptr<TH1D> PreFSIKinematics_1p;

  std::unique_ptr<TH2D> TotalNeutronKE_1p_only;

  std::unique_ptr<TH2D> PreFSIKinematics_1piplus_1p;

  std::unique_ptr<TH3D> TotalPi0E_1piplus_1p;
  std::unique_ptr<TH3D> TotalNeutralE_1piplus_1p;

  std::unique_ptr<TH1D> Transparency_1p_only;
  std::unique_ptr<TH1D> Transparency_1piplus_1p;

  std::unique_ptr<TH1D> TotalNeutronKE_1p_only_belowcut;
  std::unique_ptr<TH1D> TotalNeutronKE_1p_only_abovecut;

  std::unique_ptr<TH1D> TotalPi0E_1piplus_1p_belowcut;
  std::unique_ptr<TH1D> TotalPi0E_1piplus_1p_abovecut;
  std::unique_ptr<TH1D> TotalNeutralE_1piplus_1p_belowcut;
  std::unique_ptr<TH1D> TotalNeutralE_1piplus_1p_abovecut;

  THBlob(std::string const &fin_name, THBlob const *other) {
    TFile fin(fin_name.c_str(), "READ");

    PrimaryToFinalStateSmearing =
        std::unique_ptr<TH2D>(fin.Get<TH2D>("PrimaryToFinalStateSmearing"));
    PrimaryToFinalStateSmearing->SetDirectory(nullptr);
    RowNormTH2(PrimaryToFinalStateSmearing.get());

    PreFSIKinematics_1p =
        std::unique_ptr<TH1D>(fin.Get<TH1D>("PreFSIKinematics_1p"));
    PreFSIKinematics_1p->SetDirectory(nullptr);

    PreFSIKinematics_1piplus_1p =
        std::unique_ptr<TH2D>(fin.Get<TH2D>("PreFSIKinematics_1piplus_1p"));
    PreFSIKinematics_1piplus_1p->SetDirectory(nullptr);

    TotalNeutronKE_1p_only =
        std::unique_ptr<TH2D>(fin.Get<TH2D>("TotalNeutronKE_1p_only"));
    TotalNeutronKE_1p_only->SetDirectory(nullptr);

    if (other && reshape) {
      // reweight each row by the relevant preFSI kinematics shape.

      for (int i = 0; i < PreFSIKinematics_1p->GetXaxis()->GetNbins(); ++i) {
        double w = other->PreFSIKinematics_1p->GetBinContent(i + 1) /
                   (PreFSIKinematics_1p->GetBinContent(i + 1));
        if (std::isnormal(w)) {
          for (int j = 0; j < TotalNeutronKE_1p_only->GetXaxis()->GetNbins();
               ++j) {
            double bc = TotalNeutronKE_1p_only->GetBinContent(j + 1, i + 1);
            TotalNeutronKE_1p_only->SetBinContent(j + 1, i + 1, bc * w);
            double be = TotalNeutronKE_1p_only->GetBinError(j + 1, i + 1);
            TotalNeutronKE_1p_only->SetBinError(j + 1, i + 1, be * w);
          }
        } else {
          std::cout << "PreFSIKinematics_1p bin " << i << ": w = " << w
                    << " = ("
                    << other->PreFSIKinematics_1p->GetBinContent(i + 1) << "/"
                    << PreFSIKinematics_1p->GetBinContent(i + 1) << ")"
                    << std::endl;
        }
      }
    }

    TotalPi0E_1piplus_1p =
        std::unique_ptr<TH3D>(fin.Get<TH3D>("TotalPi0E_1piplus_1p"));
    TotalPi0E_1piplus_1p->SetDirectory(nullptr);

    TotalNeutralE_1piplus_1p =
        std::unique_ptr<TH3D>(fin.Get<TH3D>("TotalNeutralE_1piplus_1p"));
    TotalNeutralE_1piplus_1p->SetDirectory(nullptr);

    if (other && reshape) {
      // reweight each row by the relevant preFSI kinematics shape.

      for (int i = 0; i < PreFSIKinematics_1piplus_1p->GetXaxis()->GetNbins();
           ++i) {
        for (int j = 0; j < PreFSIKinematics_1piplus_1p->GetYaxis()->GetNbins();
             ++j) {
          double w =
              other->PreFSIKinematics_1piplus_1p->GetBinContent(i + 1, j + 1) /
              (PreFSIKinematics_1piplus_1p->GetBinContent(i + 1, j + 1));
          if (std::isnormal(w)) {
            for (int k = 0; k < TotalPi0E_1piplus_1p->GetXaxis()->GetNbins();
                 ++k) {
              double bc =
                  TotalPi0E_1piplus_1p->GetBinContent(k + 1, j + 1, i + 1);
              TotalPi0E_1piplus_1p->SetBinContent(k + 1, j + 1, i + 1, bc * w);
              double be =
                  TotalPi0E_1piplus_1p->GetBinError(k + 1, j + 1, i + 1);
              TotalPi0E_1piplus_1p->SetBinError(k + 1, j + 1, i + 1, be * w);
            }
          } else {
            std::cout << "PreFSIKinematics_1piplus_1p bin (" << i << "," << j
                      << "): w = " << w << " = ("
                      << other->PreFSIKinematics_1piplus_1p->GetBinContent(
                             i + 1, j + 1)
                      << "/"
                      << PreFSIKinematics_1piplus_1p->GetBinContent(i + 1,
                                                                    j + 1)
                      << ")" << std::endl;
          }
        }
      }
    }

    Transparency_1p_only =
        std::unique_ptr<TH1D>(fin.Get<TH1D>("k1p_only_proton_transp"));
    Transparency_1p_only->SetDirectory(nullptr);
    Transparency_1piplus_1p =
        std::unique_ptr<TH1D>(fin.Get<TH1D>("k1piplus_1p_piplus_transp"));
    Transparency_1piplus_1p->SetDirectory(nullptr);

    double rebin = 4;

    TotalNeutronKE_1p_only_belowcut =
        std::unique_ptr<TH1D>(CutOffZeroBin(TotalNeutronKE_1p_only->ProjectionX(
            "TotalNeutronKE_1p_only_belowcut", 2, 11)));
    TotalNeutronKE_1p_only_belowcut->Scale(
        1.0 / (other ? other->TotalNeutronKE_1p_only->Integral()
                     : TotalNeutronKE_1p_only->Integral()));
    TotalNeutronKE_1p_only_belowcut->SetDirectory(nullptr);

    TotalNeutronKE_1p_only_abovecut = std::unique_ptr<TH1D>(
        CutOffZeroBin(TotalNeutronKE_1p_only->ProjectionX(
                          "TotalNeutronKE_1p_only_abovecut", 12, 51),
                      rebin));
    TotalNeutronKE_1p_only_abovecut->Scale(
        1.0 /
        (double(rebin) * (other ? other->TotalNeutronKE_1p_only->Integral()
                                : TotalNeutronKE_1p_only->Integral())));
    TotalNeutronKE_1p_only_abovecut->SetDirectory(nullptr);

    TotalPi0E_1piplus_1p_belowcut =
        std::unique_ptr<TH1D>(CutOffZeroBin(TotalPi0E_1piplus_1p->ProjectionX(
            "TotalPi0E_1piplus_1p_belowcut", 2, 51, 2, 16)));
    TotalPi0E_1piplus_1p_belowcut->Scale(
        1.0 / (other ? other->TotalPi0E_1piplus_1p->Integral()
                     : TotalPi0E_1piplus_1p->Integral()));
    TotalPi0E_1piplus_1p_belowcut->SetDirectory(nullptr);

    TotalPi0E_1piplus_1p_abovecut = std::unique_ptr<TH1D>(
        CutOffZeroBin(TotalPi0E_1piplus_1p->ProjectionX(
                          "TotalPi0E_1piplus_1p_abovecut", 2, 51, 17, 31),
                      rebin));
    TotalPi0E_1piplus_1p_abovecut->Scale(
        1.0 / (double(rebin) * (other ? other->TotalPi0E_1piplus_1p->Integral()
                                      : TotalPi0E_1piplus_1p->Integral())));
    TotalPi0E_1piplus_1p_abovecut->SetDirectory(nullptr);

    TotalNeutralE_1piplus_1p_belowcut = std::unique_ptr<TH1D>(
        CutOffZeroBin(TotalNeutralE_1piplus_1p->ProjectionX(
            "TotalNeutralE_1piplus_1p_belowcut", 2, 51, 2, 16)));
    TotalNeutralE_1piplus_1p_belowcut->Scale(
        1.0 / (other ? other->TotalNeutralE_1piplus_1p->Integral()
                     : TotalNeutralE_1piplus_1p->Integral()));
    TotalNeutralE_1piplus_1p_belowcut->SetDirectory(nullptr);

    TotalNeutralE_1piplus_1p_abovecut = std::unique_ptr<TH1D>(
        CutOffZeroBin(TotalNeutralE_1piplus_1p->ProjectionX(
                          "TotalNeutralE_1piplus_1p_abovecut", 2, 51, 17, 31),
                      rebin));
    TotalNeutralE_1piplus_1p_abovecut->Scale(
        1.0 /
        (double(rebin) * (other ? other->TotalNeutralE_1piplus_1p->Integral()
                                : TotalNeutralE_1piplus_1p->Integral())));
    TotalNeutralE_1piplus_1p_abovecut->SetDirectory(nullptr);
  }
};

int cols[] = {TColor::GetColor("#DDAA33"), TColor::GetColor("#BB5566"),
              TColor::GetColor("#004488"), TColor::GetColor("#000000")};

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

    infs.emplace_back(
        name, THBlob(fname, infs.size() ? &infs.front().second : nullptr));
  }

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
    pleft.SetRightMargin(0.05);
    pleft.SetTopMargin(0.1);
    pleft.SetBottomMargin(0.25);
    TPad pright("pright", "", 0.5, 0, 1, 1);
    pright.AppendPad();
    pright.SetLeftMargin(0.02);
    pright.SetRightMargin(0.26);
    pright.SetTopMargin(0.1);
    pright.SetBottomMargin(0.25);

    pleft.cd();
    int first = 0;

    TLegend *legendl = new TLegend(0.62, 0.575, 0.7975, 0.875);
    legendl->SetTextFont(132);
    legendl->SetTextSize(fontsize);
    legendl->SetBorderSize(0);
    legendl->SetFillStyle(0);
    legendl->SetNColumns(1);

    double maxy = 5E-2;

    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;
      if (hb.TotalNeutronKE_1p_only_belowcut->Integral() < 1E-8) {
        continue;
      }

      hb.TotalNeutronKE_1p_only_belowcut->GetYaxis()->SetRangeUser(0, maxy);
      hb.TotalNeutronKE_1p_only_belowcut->GetYaxis()->SetLabelSize(fontsize);
      hb.TotalNeutronKE_1p_only_belowcut->GetYaxis()->SetTitleSize(fontsize);
      hb.TotalNeutronKE_1p_only_belowcut->GetYaxis()->SetTitleFont(132);
      hb.TotalNeutronKE_1p_only_belowcut->GetYaxis()->SetLabelFont(132);
      hb.TotalNeutronKE_1p_only_belowcut->GetYaxis()->SetTitle("P.D.F.");
      hb.TotalNeutronKE_1p_only_belowcut->GetYaxis()->SetNdivisions(505);

      hb.TotalNeutronKE_1p_only_belowcut->GetXaxis()->SetRangeUser(0, 0.18);
      hb.TotalNeutronKE_1p_only_belowcut->GetXaxis()->SetLabelSize(fontsize);
      hb.TotalNeutronKE_1p_only_belowcut->GetXaxis()->SetTitleSize(fontsize);
      hb.TotalNeutronKE_1p_only_belowcut->GetXaxis()->SetLabelFont(132);
      hb.TotalNeutronKE_1p_only_belowcut->GetXaxis()->SetTitleOffset(1);
      hb.TotalNeutronKE_1p_only_belowcut->GetXaxis()->SetNdivisions(505);
      hb.TotalNeutronKE_1p_only_belowcut->GetXaxis()->SetTitle("");

      hb.TotalNeutronKE_1p_only_belowcut->SetLineColor(cols[g_it]);
      hb.TotalNeutronKE_1p_only_belowcut->SetLineWidth(3);
      hb.TotalNeutronKE_1p_only_belowcut->Draw((!first++) ? "EHIST"
                                                          : "EHIST SAME");

      legendl->AddEntry(hb.TotalNeutronKE_1p_only_belowcut.get(),
                        infs[g_it].first.c_str(), "l");
    }

    pright.cd();
    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;
      if (hb.TotalNeutronKE_1p_only_abovecut->Integral() < 1E-8) {
        continue;
      }
      hb.TotalNeutronKE_1p_only_abovecut->Scale(3);
      hb.TotalNeutronKE_1p_only_abovecut->GetYaxis()->SetRangeUser(0, maxy);
      hb.TotalNeutronKE_1p_only_abovecut->GetYaxis()->SetLabelOffset(1);
      hb.TotalNeutronKE_1p_only_abovecut->GetYaxis()->SetLabelSize(fontsize);
      hb.TotalNeutronKE_1p_only_abovecut->GetYaxis()->SetTitleSize(fontsize);
      hb.TotalNeutronKE_1p_only_abovecut->GetYaxis()->SetTitleFont(132);
      hb.TotalNeutronKE_1p_only_abovecut->GetYaxis()->SetLabelFont(132);
      hb.TotalNeutronKE_1p_only_abovecut->GetYaxis()->SetNdivisions(505);

      hb.TotalNeutronKE_1p_only_abovecut->GetXaxis()->SetRangeUser(0, 1);
      hb.TotalNeutronKE_1p_only_abovecut->GetXaxis()->SetLabelSize(fontsize);
      hb.TotalNeutronKE_1p_only_abovecut->GetXaxis()->SetTitleSize(fontsize);
      hb.TotalNeutronKE_1p_only_abovecut->GetXaxis()->SetLabelFont(132);
      hb.TotalNeutronKE_1p_only_abovecut->GetXaxis()->SetTitleOffset(1);
      hb.TotalNeutronKE_1p_only_abovecut->GetXaxis()->SetNdivisions(505);
      hb.TotalNeutronKE_1p_only_abovecut->GetXaxis()->SetTitle("");

      hb.TotalNeutronKE_1p_only_abovecut->SetLineColor(cols[g_it]);
      hb.TotalNeutronKE_1p_only_abovecut->SetLineWidth(3);
      hb.TotalNeutronKE_1p_only_abovecut->Draw((!first++) ? "EHIST"
                                                          : "EHIST SAME");
    }
    c1.cd();
    ltx.DrawLatexNDC(0.75, 0.5, "#times 3");
    ltx.SetTextAlign(32);
    ltx.DrawLatexNDC(0.39, 0.08, "#sum");
    ltx.SetTextAlign(12);
    ltx.DrawLatexNDC(0.41, 0.08, "#it{T}_{neutron} (GeV)");

    ltx.DrawLatexNDC(0.2, 0.95, "#it{T}_{p}^{prim.} < 0.2 GeV");
    ltx.DrawLatexNDC(0.55, 0.95, "0.2 < #it{T}_{p}^{prim.} < 1 GeV");

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
    pleft.SetRightMargin(0.05);
    pleft.SetTopMargin(0.1);
    pleft.SetBottomMargin(0.25);
    TPad pright("pright", "", 0.5, 0, 1, 1);
    pright.AppendPad();
    pright.SetLeftMargin(0.02);
    pright.SetRightMargin(0.26);
    pright.SetTopMargin(0.1);
    pright.SetBottomMargin(0.25);

    pleft.cd();
    int first = 0;
    TLegend *legendl = new TLegend(0.7, 0.6, 0.8775, 0.9);
    legendl->SetTextFont(132);
    legendl->SetTextSize(fontsize);
    legendl->SetBorderSize(0);
    legendl->SetFillStyle(0);
    legendl->SetNColumns(1);

    double maxy = 1.05E-2;

    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;
      if (hb.TotalPi0E_1piplus_1p_belowcut->Integral() < 1E-8) {
        continue;
      }

      hb.TotalPi0E_1piplus_1p_belowcut->GetYaxis()->SetRangeUser(0, maxy);
      hb.TotalPi0E_1piplus_1p_belowcut->GetYaxis()->SetLabelSize(fontsize);
      hb.TotalPi0E_1piplus_1p_belowcut->GetYaxis()->SetTitleSize(fontsize);
      hb.TotalPi0E_1piplus_1p_belowcut->GetYaxis()->SetTitleFont(132);
      hb.TotalPi0E_1piplus_1p_belowcut->GetYaxis()->SetLabelFont(132);
      hb.TotalPi0E_1piplus_1p_belowcut->GetYaxis()->SetTitle("P.D.F.");
      hb.TotalPi0E_1piplus_1p_belowcut->GetYaxis()->SetNdivisions(505);

      hb.TotalPi0E_1piplus_1p_belowcut->GetXaxis()->SetRangeUser(0.13, 0.5);
      hb.TotalPi0E_1piplus_1p_belowcut->GetXaxis()->SetLabelSize(fontsize);
      hb.TotalPi0E_1piplus_1p_belowcut->GetXaxis()->SetTitleSize(fontsize);
      hb.TotalPi0E_1piplus_1p_belowcut->GetXaxis()->SetLabelFont(132);
      hb.TotalPi0E_1piplus_1p_belowcut->GetXaxis()->SetTitleOffset(1);
      hb.TotalPi0E_1piplus_1p_belowcut->GetXaxis()->SetNdivisions(505);
      hb.TotalPi0E_1piplus_1p_belowcut->GetXaxis()->SetTitle("");

      hb.TotalPi0E_1piplus_1p_belowcut->SetLineColor(cols[g_it]);
      hb.TotalPi0E_1piplus_1p_belowcut->SetLineWidth(3);
      hb.TotalPi0E_1piplus_1p_belowcut->Draw((!first++) ? "EHIST"
                                                        : "EHIST SAME");
      legendl->AddEntry(hb.TotalPi0E_1piplus_1p_belowcut.get(),
                        infs[g_it].first.c_str(), "l");
    }

    pright.cd();
    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;
      if (hb.TotalPi0E_1piplus_1p_abovecut->Integral() < 1E-8) {
        continue;
      }
      hb.TotalPi0E_1piplus_1p_abovecut->Scale(3);
      hb.TotalPi0E_1piplus_1p_abovecut->GetYaxis()->SetRangeUser(0, maxy);
      hb.TotalPi0E_1piplus_1p_abovecut->GetYaxis()->SetLabelOffset(1);
      hb.TotalPi0E_1piplus_1p_abovecut->GetYaxis()->SetLabelSize(fontsize);
      hb.TotalPi0E_1piplus_1p_abovecut->GetYaxis()->SetTitleSize(fontsize);
      hb.TotalPi0E_1piplus_1p_abovecut->GetYaxis()->SetTitleFont(132);
      hb.TotalPi0E_1piplus_1p_abovecut->GetYaxis()->SetLabelFont(132);
      hb.TotalPi0E_1piplus_1p_abovecut->GetYaxis()->SetNdivisions(505);

      hb.TotalPi0E_1piplus_1p_abovecut->GetXaxis()->SetRangeUser(0.13, 1);
      hb.TotalPi0E_1piplus_1p_abovecut->GetXaxis()->SetLabelSize(fontsize);
      hb.TotalPi0E_1piplus_1p_abovecut->GetXaxis()->SetTitleSize(fontsize);
      hb.TotalPi0E_1piplus_1p_abovecut->GetXaxis()->SetLabelFont(132);
      hb.TotalPi0E_1piplus_1p_abovecut->GetXaxis()->SetTitleOffset(1);
      hb.TotalPi0E_1piplus_1p_abovecut->GetXaxis()->SetNdivisions(505);
      hb.TotalPi0E_1piplus_1p_abovecut->GetXaxis()->SetTitle("");

      hb.TotalPi0E_1piplus_1p_abovecut->SetLineColor(cols[g_it]);
      hb.TotalPi0E_1piplus_1p_abovecut->SetLineWidth(3);
      hb.TotalPi0E_1piplus_1p_abovecut->Draw((!first++) ? "EHIST"
                                                        : "EHIST SAME");
    }
    c1.cd();
    ltx.DrawLatexNDC(0.75, 0.5, "#times 3");
    ltx.SetTextAlign(32);
    ltx.DrawLatexNDC(0.39, 0.08, "#sum");
    ltx.SetTextAlign(12);
    ltx.DrawLatexNDC(0.41, 0.08, "#it{E}_{#pi^{0}} (GeV)");

    ltx.DrawLatexNDC(0.2, 0.95, "#it{T}_{#pi^{+}}^{prim.} < 0.3 GeV");
    ltx.DrawLatexNDC(0.55, 0.95, "0.3 < #it{T}_{#pi^{+}}^{prim.} < 1 GeV");

    legendl->Draw();
    c1.Print("TotalPi0E_1piplus_1p.pdf");
    c1.Print("allplots.pdf");
  }

  { // 1p->neutron
    TCanvas c1("c1", "", 1400, 800);

    TPad pleft("pleft", "", 0, 0, 0.5, 1);
    pleft.AppendPad();
    pleft.SetLeftMargin(0.28);
    pleft.SetRightMargin(0.05);
    pleft.SetTopMargin(0.1);
    pleft.SetBottomMargin(0.25);
    TPad pright("pright", "", 0.5, 0, 1, 1);
    pright.AppendPad();
    pright.SetLeftMargin(0.02);
    pright.SetRightMargin(0.26);
    pright.SetTopMargin(0.1);
    pright.SetBottomMargin(0.25);

    pleft.cd();
    int first = 0;
    TLegend *legendl = new TLegend(0.7, 0.6, 0.8775, 0.9);
    legendl->SetTextFont(132);
    legendl->SetTextSize(fontsize);
    legendl->SetBorderSize(0);
    legendl->SetFillStyle(0);
    legendl->SetNColumns(1);

    double maxy = 12E-2;

    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;
      if (hb.TotalNeutralE_1piplus_1p_belowcut->Integral() < 1E-8) {
        continue;
      }

      hb.TotalNeutralE_1piplus_1p_belowcut->GetYaxis()->SetRangeUser(0, maxy);
      hb.TotalNeutralE_1piplus_1p_belowcut->GetYaxis()->SetLabelSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_belowcut->GetYaxis()->SetTitleSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_belowcut->GetYaxis()->SetTitleFont(132);
      hb.TotalNeutralE_1piplus_1p_belowcut->GetYaxis()->SetLabelFont(132);
      hb.TotalNeutralE_1piplus_1p_belowcut->GetYaxis()->SetTitle("P.D.F.");
      hb.TotalNeutralE_1piplus_1p_belowcut->GetYaxis()->SetNdivisions(505);

      hb.TotalNeutralE_1piplus_1p_belowcut->GetXaxis()->SetRangeUser(0, 0.6);
      hb.TotalNeutralE_1piplus_1p_belowcut->GetXaxis()->SetLabelSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_belowcut->GetXaxis()->SetTitleSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_belowcut->GetXaxis()->SetLabelFont(132);
      hb.TotalNeutralE_1piplus_1p_belowcut->GetXaxis()->SetTitleOffset(1);
      hb.TotalNeutralE_1piplus_1p_belowcut->GetXaxis()->SetNdivisions(505);
      hb.TotalNeutralE_1piplus_1p_belowcut->GetXaxis()->SetTitle("");

      hb.TotalNeutralE_1piplus_1p_belowcut->SetLineColor(cols[g_it]);
      hb.TotalNeutralE_1piplus_1p_belowcut->SetLineWidth(3);
      hb.TotalNeutralE_1piplus_1p_belowcut->Draw((!first++) ? "EHIST"
                                                            : "EHIST SAME");
      legendl->AddEntry(hb.TotalNeutralE_1piplus_1p_belowcut.get(),
                        infs[g_it].first.c_str(), "l");
    }

    pright.cd();
    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;
      if (hb.TotalNeutralE_1piplus_1p_abovecut->Integral() < 1E-8) {
        continue;
      }
      hb.TotalNeutralE_1piplus_1p_abovecut->Scale(7.5);
      hb.TotalNeutralE_1piplus_1p_abovecut->GetYaxis()->SetRangeUser(0, maxy);
      hb.TotalNeutralE_1piplus_1p_abovecut->GetYaxis()->SetLabelOffset(1);
      hb.TotalNeutralE_1piplus_1p_abovecut->GetYaxis()->SetLabelSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_abovecut->GetYaxis()->SetTitleSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_abovecut->GetYaxis()->SetTitleFont(132);
      hb.TotalNeutralE_1piplus_1p_abovecut->GetYaxis()->SetLabelFont(132);
      hb.TotalNeutralE_1piplus_1p_abovecut->GetYaxis()->SetNdivisions(505);

      hb.TotalNeutralE_1piplus_1p_abovecut->GetXaxis()->SetRangeUser(0, 1);
      hb.TotalNeutralE_1piplus_1p_abovecut->GetXaxis()->SetLabelSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_abovecut->GetXaxis()->SetTitleSize(fontsize);
      hb.TotalNeutralE_1piplus_1p_abovecut->GetXaxis()->SetLabelFont(132);
      hb.TotalNeutralE_1piplus_1p_abovecut->GetXaxis()->SetTitleOffset(1);
      hb.TotalNeutralE_1piplus_1p_abovecut->GetXaxis()->SetNdivisions(505);
      hb.TotalNeutralE_1piplus_1p_abovecut->GetXaxis()->SetTitle("");

      hb.TotalNeutralE_1piplus_1p_abovecut->SetLineColor(cols[g_it]);
      hb.TotalNeutralE_1piplus_1p_abovecut->SetLineWidth(3);
      hb.TotalNeutralE_1piplus_1p_abovecut->Draw((!first++) ? "EHIST"
                                                            : "EHIST SAME");
    }
    c1.cd();
    ltx.DrawLatexNDC(0.75, 0.5, "#times 7.5");
    ltx.SetTextAlign(32);
    ltx.DrawLatexNDC(0.39, 0.08, "#sum");
    ltx.SetTextAlign(12);
    ltx.DrawLatexNDC(0.41, 0.08, "#it{T}_{neutron} + #it{E}_{#pi^{0}} (GeV)");

    ltx.DrawLatexNDC(0.2, 0.95, "#it{T}_{#pi^{+}}^{prim.} < 0.3 GeV");
    ltx.DrawLatexNDC(0.55, 0.95, "0.3 < #it{T}_{#pi^{+}}^{prim.} < 1 GeV");

    legendl->Draw();
    c1.Print("TotalNeutralE_1piplus_1p.pdf");
    c1.Print("allplots.pdf");
  }

  { // transparency
    TCanvas c1("c1", "", 1200, 1200);

    c1.SetLeftMargin(0.25);
    c1.SetRightMargin(0.05);
    c1.SetTopMargin(0.05);
    c1.SetBottomMargin(0.25);

    int first = 0;
    TLegend *legendl = new TLegend(0.5, 0.65, 0.94, 0.93);
    legendl->SetTextFont(132);
    legendl->SetTextSize(fontsize);
    legendl->SetBorderSize(0);
    legendl->SetFillStyle(0);
    legendl->SetNColumns(1);

    double maxy = 1;

    first = 0;
    for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
      auto &hb = infs[g_it].second;

      hb.Transparency_1p_only->GetYaxis()->SetRangeUser(0, maxy);
      hb.Transparency_1p_only->GetYaxis()->SetLabelSize(fontsize);
      hb.Transparency_1p_only->GetYaxis()->SetTitleSize(fontsize);
      hb.Transparency_1p_only->GetYaxis()->SetTitleFont(132);
      hb.Transparency_1p_only->GetYaxis()->SetLabelFont(132);
      hb.Transparency_1p_only->GetYaxis()->SetTitle(
          "Proton Transparency, #theta_{scatter} < 5^{#circ}");
      hb.Transparency_1p_only->GetYaxis()->SetNdivisions(505);

      hb.Transparency_1p_only->GetXaxis()->SetRangeUser(0.02, 1);
      hb.Transparency_1p_only->GetXaxis()->SetLabelSize(fontsize);
      hb.Transparency_1p_only->GetXaxis()->SetTitleSize(fontsize);
      hb.Transparency_1p_only->GetXaxis()->SetLabelFont(132);
      hb.Transparency_1p_only->GetXaxis()->SetTitleFont(132);
      hb.Transparency_1p_only->GetXaxis()->SetTitleOffset(1);
      hb.Transparency_1p_only->GetXaxis()->SetNdivisions(505);
      hb.Transparency_1p_only->GetXaxis()->SetTitle("T_{Proton} GeV");

      hb.Transparency_1p_only->SetLineColor(cols[g_it]);
      hb.Transparency_1p_only->SetLineWidth(4);
      hb.Transparency_1p_only->Draw((!first++) ? "EHIST" : "EHIST SAME");
      legendl->AddEntry(hb.Transparency_1p_only.get(), infs[g_it].first.c_str(),
                        "l");
    }

    legendl->Draw();
    c1.Print("TotalNeutralE_1piplus_1p.pdf");
    c1.Print("allplots.pdf");
  }
  {

    TCanvas c1("c1", "", 1200, 1200);

    int first = 0;
    TLegend *legendl = new TLegend(0.05, 0.9, 0.95, 1);
    legendl->SetTextFont(132);
    legendl->SetTextSize(fontsize * 0.75);
    legendl->SetBorderSize(0);
    legendl->SetFillStyle(0);
    legendl->SetNColumns(4);

    std::vector<TPad *> pads;

    double ywidth = 0.25;

    for (int i = 0; i < 3; ++i) {
      pads.emplace_back(new TPad(("p_" + std::to_string(i)).c_str(), "", 0.,
                                 0.175 + ywidth * i, 1.,
                                 0.175 + ywidth * (i + 1)));
      pads.back()->AppendPad();
      pads.back()->SetLeftMargin(0.1);
      pads.back()->SetRightMargin(0.1);
      pads.back()->SetTopMargin(0.05);
      pads.back()->SetBottomMargin(0.05);
    }

    double maxy = 0.65;

    int indx[] = {1, 3, 4};

    for (int i = 0; i < 3; ++i) {
      pads[i]->cd();

      first = 0;
      for (size_t g_it = 0; g_it < infs.size(); ++g_it) {
        auto &hb = infs[g_it].second;

        auto proj = hb.PrimaryToFinalStateSmearing->ProjectionX(
            ("PrimaryToFinalStateSmearing_" + std::to_string(g_it) + "_" +
             std::to_string(i))
                .c_str(),
            indx[i], indx[i]);

        proj->GetYaxis()->SetRangeUser(0, maxy);
        proj->GetYaxis()->SetLabelSize(fontsize * 2);
        proj->GetYaxis()->SetTitleSize(fontsize * 2);
        proj->GetYaxis()->SetTitleFont(132);
        proj->GetYaxis()->SetLabelFont(132);
        proj->GetXaxis()->SetTitleOffset(0.01);
        proj->GetYaxis()->SetTitle("P.D.F");
        proj->GetYaxis()->SetNdivisions(505);

        proj->GetXaxis()->SetRangeUser(0, 1);
        proj->GetXaxis()->SetLabelSize(fontsize);
        proj->GetXaxis()->SetTitleSize(fontsize);
        proj->GetXaxis()->SetLabelFont(132);
        proj->GetXaxis()->SetTitleOffset(1);
        proj->GetXaxis()->SetLabelOffset(1);
        proj->GetXaxis()->SetNdivisions(505);

        proj->SetLineColor(cols[g_it]);
        proj->SetLineWidth(3);
        proj->Draw((!first++) ? "EHIST" : "EHIST SAME");
        if (!i) {
          legendl->AddEntry(proj, infs[g_it].first.c_str(), "l");
        }

        ltx.SetTextAlign(12);
        ltx.SetTextSize(0.15);
        ltx.DrawLatexNDC(
            0.3, 0.8,
            (std::string("Pre-FSI Topology: ") +
             hb.PrimaryToFinalStateSmearing->GetYaxis()->GetBinLabel(indx[i]))
                .c_str());
      }
    }

    c1.cd();

    for (int i = 0; i < infs.front()
                            .second.PrimaryToFinalStateSmearing->GetXaxis()
                            ->GetNbins();
         ++i) {
      double l = 0.1, r = 0.1;
      double w = 1 - (l + r);
      double bw = w / double(infs.front()
                                 .second.PrimaryToFinalStateSmearing->GetXaxis()
                                 ->GetNbins());
      ltx.SetTextSize(0.03);
      ltx.SetTextAngle(-45);
      ltx.DrawLatexNDC(l + bw * i + 0.5 * bw, 0.175,
                       infs.front()
                           .second.PrimaryToFinalStateSmearing->GetXaxis()
                           ->GetBinLabel(i + 1));
    }

    ltx.SetTextSize(fontsize);
    ltx.SetTextAngle(0);
    ltx.DrawLatexNDC(0.25, 0.05, "Final State Topology");

    legendl->Draw();
    c1.Print("TotalNeutralE_1piplus_1p.pdf");
    c1.Print("allplots.pdf");
    c1.Print("allplots.pdf]");
  }
}