#include "inc_gar_discover.h"
#include "../types.h"
#include <gflags/gflags.h>
#include <gflags/gflags_declare.h>
#include <glog/logging.h>

int main(int argc, char* argv[]) {
  FLAGS_stderrthreshold = 0;

  gflags::SetUsageMessage(
      "Usage: mpiexec [mpi_opts] ./inc_gar_discover [grape_opts]");
  if (argc == 1) {
    gflags::ShowUsageWithFlagsRestrict(argv[0], "inc_gar_discover");
    exit(1);
  }
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  gflags::ShutDownCommandLineFlags();

  google::InitGoogleLogging("inc_gar_discover");
  google::InstallFailureSignalHandler();

  grape::Init();

  grape::Run<int64_t, uint32_t, grape::Data, grape::EdgeData>();

  grape::Finalize();

  google::ShutdownGoogleLogging();
}
