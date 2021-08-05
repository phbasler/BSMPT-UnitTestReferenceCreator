// SPDX-FileCopyrightText: 2021 Philipp Basler
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <exception>
#include <iostream>
#include <stdlib.h>
#include <vector>

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>

#include <fstream>
#include <map>

using std::exception;

int main(int argc, char *argv[])
try
{
  const std::vector<double> example_point_CPINTHEDARK{
      /* m11s = */ -7823.7540500000005,
      /* m22s = */ 242571.64899822656,
      /* mSs = */ 109399.20176343,
      /* ReA = */ 93.784159581909734,
      /* ImA = */ 126.30387933116994,
      /* L1 = */ 0.25810698810286969,
      /* L2 = */ 4.6911643599657609,
      /* L3 = */ -0.21517372505705856,
      /* L4 = */ -0.42508424793839744,
      /* L5 = */ -0.13790431680607695,
      /* L6 = */ 15.075540949860104,
      /* L7 = */ 6.7788372529237835,
      /* L8 = */ -1.8651245632976341};

  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CPINTHEDARK);
  modelPointer->initModel(example_point_CPINTHEDARK);

  const std::string ClassName{"Compare_CPINTHEDARK"};
  const std::string headerFileName{"CPINTHEDARK.h"};
  std::ofstream header(headerFileName);
  header
      << "// SPDX-FileCopyrightText: 2021 Philipp Basler \n"
      << "//\n"
      << "// SPDX-License-Identifier: GPL-3.0-or-later\n"
      << "#include <BSMPT/minimizer/Minimizer.h>\n"
      << "#include <map>\n"
      << "#include <vector>\n"
      << "class Compare_CPINTHEDARK\n "
      << "{\n"
      << "public:\n"
      << "\tusing Matrix3D = std::vector<std::vector<std::vector<double>>>;\n"
      << "\tusing Matrix2D = std::vector<std::vector<double>>;\n"
      << "\t" << ClassName << "();\n"
      << "\tMatrix3D CheckTripleCT;\n"
      << "\tMatrix3D CheckTripleCW;\n"
      << "\tMatrix3D CheckTripleTree;\n"
      << "\tstd::map<int, BSMPT::Minimizer::EWPTReturnType> EWPTPerSetting;\n"
      << "};\n";
  header.close();

  int WhichMin;
  Minimizer::EWPTReturnType EWPT;

  std::ofstream source("CPINTHEDARK.cpp");
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
