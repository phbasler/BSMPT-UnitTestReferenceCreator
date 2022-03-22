// SPDX-FileCopyrightText: 2021 Philipp Basler
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <exception>
#include <iostream>
#include <stdlib.h>
#include <vector>

#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>

#include <fstream>
#include <map>

using std::exception;

int main(int argc, char *argv[])
try
{
  const std::vector<double> example_point_C2HDM{/* lambda_1 = */ 3.29771,
                                                /* lambda_2 = */ 0.274365,
                                                /* lambda_3 = */ 4.71019,
                                                /* lambda_4 = */ -2.23056,
                                                /* Re(lambda_5) = */ -2.43487,
                                                /* Im(lambda_5) = */ 0.124948,
                                                /* Re(m_{12}^2) = */ 2706.86,
                                                /* tan(beta) = */ 4.64487,
                                                /* Yukawa Type = */ 1};

  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);

  const double testVW = 0.1;

  const std::string ClassName{"Compare_C2HDM"};
  const std::string headerFileName{"C2HDM.h"};
  std::ofstream header(headerFileName);
  header
      << "// SPDX-FileCopyrightText: 2021 Philipp Basler \n"
      << "//\n"
      << "// SPDX-License-Identifier: GPL-3.0-or-later\n"
      << "#include <BSMPT/minimizer/Minimizer.h>\n"
      << "#include <map>\n"
      << "#include <vector>\n"
      << "class Compare_C2HDM\n "
      << "{\n"
      << "public:\n"
      << "  using Matrix3D = std::vector<std::vector<std::vector<double>>>;\n"
      << "  using Matrix2D = std::vector<std::vector<double>>;\n"
      << "  " << ClassName << "();\n"
      << "  Matrix3D CheckTripleCT;\n"
      << "  Matrix3D CheckTripleCW;\n"
      << "  Matrix3D CheckTripleTree;\n"
      << "  std::map<int, BSMPT::Minimizer::EWPTReturnType> EWPTPerSetting;\n"
      << "  std::map<int,double> LWPerSetting;\n"
      << "  std::map<int,std::vector<double>> vevSymmetricPerSetting;\n"
      << "  std::map<int,std::vector<double>> etaPerSetting;\n"
      << "  const double testVW = " << testVW << ";\n"
      << "};\n";
  header.close();

  int WhichMin;
  Minimizer::EWPTReturnType EWPT;

  std::ofstream source("C2HDM.cpp");
  source << "// SPDX-FileCopyrightText: 2021 Philipp Basler \n"
         << "//\n"
         << "// SPDX-License-Identifier: GPL-3.0-or-later\n"
         << "#include \"" << headerFileName << "\" \n"
         << ClassName << "::" << ClassName << "()\n"
         << "{\n"
         << "  std::size_t NHiggs = " << modelPointer->get_NHiggs() << ";\n"
         << "  CheckTripleTree = Matrix3D{NHiggs, Matrix2D{NHiggs, "
            "  std::vector<double>(NHiggs, 0)}};\n"
         << "  CheckTripleCW =   Matrix3D{NHiggs, Matrix2D{NHiggs, "
            "  std::vector<double>(NHiggs, 0)}};\n"
         << "  CheckTripleCT =   Matrix3D{NHiggs, Matrix2D{NHiggs, "
            "  std::vector<double>(NHiggs, 0)}};\n";

  std::map<int, Minimizer::EWPTReturnType> mdata;
  for (bool UseGSL : {false, true})
  {
    for (bool UseCMAES : {false, true})
    {
      for (bool UseNLopt : {false, true})
      {
        if (not UseGSL and not UseCMAES and not UseNLopt) continue;

        WhichMin = Minimizer::CalcWhichMinimizer(UseGSL, UseCMAES, UseNLopt);
        EWPT     = Minimizer::PTFinder_gen_all(modelPointer, 0, 300, WhichMin);
        mdata[WhichMin] = EWPT;
        source << "  EWPTPerSetting[" << WhichMin
               << "].Tc = " << mdata[WhichMin].Tc << ";" << std::endl
               << "  EWPTPerSetting[" << WhichMin
               << "].vc = " << mdata[WhichMin].vc << ";" << std::endl;
        for (const auto &el : EWPT.EWMinimum)
        {
          if (std::abs(el) > 1e-5)
            source << "  EWPTPerSetting[" << WhichMin
                   << "].EWMinimum.push_back(" << el << ");" << std::endl;
          else
            source << "  EWPTPerSetting[" << WhichMin
                   << "].EWMinimum.push_back(" << 0 << ");" << std::endl;
        }

        std::vector<double> vevsymmetricSolution, checksym, startpoint;
        for (const auto &el : EWPT.EWMinimum)
          startpoint.push_back(0.5 * el);
        vevsymmetricSolution = Minimizer::Minimize_gen_all(
            modelPointer, EWPT.Tc + 1, checksym, startpoint, WhichMin, true);

        for (const auto &el : vevsymmetricSolution)
        {
          const auto value = (std::abs(el) > 1e-5) ? el : 0;
          source << "  vevSymmetricPerSetting[" << WhichMin << "].push_back("
                 << value << ");" << std::endl;
        }

        auto config =
            std::pair<std::vector<bool>, int>{std::vector<bool>(5, true), 1};
        Baryo::CalculateEtaInterface EtaInterface(config);

        if (EWPT.vc / EWPT.Tc > 1)
        {
          auto eta = EtaInterface.CalcEta(testVW,
                                          EWPT.EWMinimum,
                                          vevsymmetricSolution,
                                          EWPT.Tc,
                                          modelPointer,
                                          Minimizer::WhichMinimizerDefault);

          source << "  LWPerSetting[" << WhichMin
                 << "] = " << EtaInterface.getLW() << ";" << std::endl;

          for (const auto &el : eta)
          {
            source << "  etaPerSetting[" << WhichMin << "].push_back(" << el
                   << ");" << std::endl;
          }
        }
      }
    }
  }

  modelPointer->Prepare_Triple();
  modelPointer->TripleHiggsCouplings();
  auto NHiggs = modelPointer->get_NHiggs();
  for (std::size_t i{0}; i < NHiggs; ++i)
  {
    for (std::size_t j{0}; j < NHiggs; ++j)
    {
      for (std::size_t k{0}; k < NHiggs; ++k)
      {
        auto value =
            modelPointer->get_TripleHiggsCorrectionsTreePhysical(i, j, k);
        if (value != 0)
        {
          source << "  CheckTripleTree.at(" << i << ").at(" << j << ").at(" << k
                 << ") = " << value << ";\n";
        }
        value = modelPointer->get_TripleHiggsCorrectionsCTPhysical(i, j, k);
        if (value != 0)
        {
          source << "  CheckTripleCT.at(" << i << ").at(" << j << ").at(" << k
                 << ") = " << value << ";\n";
        }
        value = modelPointer->get_TripleHiggsCorrectionsCWPhysical(i, j, k);
        if (value != 0)
        {
          source << "  CheckTripleCW.at(" << i << ").at(" << j << ").at(" << k
                 << ") = " << value << ";\n";
        }
      }
    }
  }

  source << "}\n";
  source.close();

  return EXIT_SUCCESS;
}
catch (int)
{
  return EXIT_SUCCESS;
}
catch (exception &e)
{
  std::cerr << e.what() << std::endl;
  return EXIT_FAILURE;
}
