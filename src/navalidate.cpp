#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>
#include "validate.h"


int main(int argc, char **argv) {
    CCP4Program prog("navalidate", "0.0.1", "$Date: 2023/12/17");
    prog.set_termination_message("Failed");

    std::cout << std::endl << "Copyright 2023 Jordan Dialpuri and University of York." << std::endl << std::endl;
    prog.summary_beg();
    prog.summary_end();

    std::string pdb_in;

    CCP4CommandInput args(argc, argv, true);
    int arg = 0;
    while (++arg < args.size()) {
        if (args[arg] == "-pdbin") {
            if (arg++ < args.size()) {
                pdb_in = args[arg];
            }
        } else {
            std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
            args.clear();
        }
    }

    clipper::MiniMol input_model;
    clipper::MMDBfile mfile;
    mfile.read_file(pdb_in);
    mfile.import_minimol(input_model);

    Validate validator = Validate(input_model);
    validator.validate();


    prog.set_termination_message( "Normal termination" );

    return 0;
}
