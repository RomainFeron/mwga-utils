#include "maf_stats.h"

void run_mafstats(int argc, char* argv[]) {

    Parameters parameters = parse_args(argc, argv);
    (*commands[parameters.command])(parameters);  // Call the MafStats method corresponding to the specified CLI subcommand
}

