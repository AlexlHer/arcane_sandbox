#include "arcane/launcher/ArcaneLauncher.h"
 
#include "arcane/utils/ITraceMng.h"
#include "arcane/utils/FatalErrorException.h"
#include "arcane/utils/Real3.h"
 
#include "arcane/core/MeshReaderMng.h"
#include "arcane/core/IMesh.h"
#include "arcane/core/ISubDomain.h"
#include "arcane/core/IParallelMng.h"
#include "arcane/core/ItemGroup.h"
#include "arcane/core/VariableTypes.h"
 
#include "arcane/utils/Exception.h"
 
#include <iostream>


using namespace Arcane;
 
void executeSample(const String& case_file)
{
  StandaloneSubDomain launcher{ ArcaneLauncher::createStandaloneSubDomain(case_file) };
  ISubDomain* sd = launcher.subDomain();

  ITraceMng* tm = launcher.traceMng();

  MeshReaderMng mrm(sd);

  IMesh* mesh = mrm.readMesh("AdditionalMesh", "plancher.msh", sd->parallelMng());
 
  Int32 nb_cell = mesh->nbCell();
  tm->info() << "NB_CELL=" << nb_cell;

  VariableNodeReal3& nodes_coordinates = mesh->nodesCoordinates();
  mesh->allCells().isAllItems();
  ENUMERATE_ (Cell, icell, mesh->allCells()) {
    icell.localId();
    Cell cell = *icell;
    Real3 cell_center;
    for (Node node : cell.nodes()) {
      cell_center += nodes_coordinates[node];
    }
    cell_center /= cell.nbNode();
    tm->info() << "Cell=" << cell.uniqueId() << " center=" << cell_center;
  }
}

int main(int argc, char* argv[])
{
  String case_file;
 
  auto func = [&] {
    std::cout << "Sample: StandaloneSubDomain\n";

    CommandLineArguments cmd_line_args(&argc, &argv);
    ArcaneLauncher::init(cmd_line_args);
    if (argc > 1)
      case_file = argv[argc - 1];
    executeSample(case_file);
  };
 
  return arcaneCallFunctionAndCatchException(func);
}

// int
// main(int argc,char* argv[])
// {
//   ArcaneLauncher::init(CommandLineArguments(&argc,&argv));
//   auto& app_build_info = ArcaneLauncher::applicationBuildInfo();
//   app_build_info.setCodeName("HelloWorld");
//   app_build_info.setCodeVersion(VersionInfo(1,0,0));
//   if(ArcaneLauncher::printHelp()){
//     return 0;
//   }
//   return ArcaneLauncher::run();
// }