# Parallelization of Sound Waves Propagation Using MPI

## Overview

This project demonstrates the parallelization of sound wave propagation using the Message Passing Interface (MPI). It includes code and resources to compile and run simulations of sound wave propagation in different scenarios.

## Files and Directories

- **main.c**: The main program file for the simulation.
- **main_Using_MpiPrint.c**: Alternative main file with additional MPI printing functions.
- **utils.c / utils.h**: Utility functions and headers.
- **visu.c**: Visualization functions.
- **parameters**: Directory containing parameter files.
- **Makefile**: Script to compile and manage the project.
- **mpSlurm.sh**: SLURM script for job submission.
- **commands.txt**: Command list for running the simulations.
- **project2022.pdf**: Detailed project documentation.
- **HSPC presentation.pptx**: Presentation of the project.

## Compilation and Execution

### Compilation

To compile the project, use the following command:
```sh
make
```

### Running the Simulation

To run the simulation, execute:
```sh
./a.out param.txt
```
You can modify the simulation parameters by editing the `param.txt` file.

### Cleaning Up

To clean up the generated files, use:
```sh
make delete
```

To clean the executables and object files, use:
```sh
make clean
```

### Generating a Movie

To generate a movie from the simulation output, use:
```sh
gen_movie
```

## Specifications

Simulation outputs are saved with step numbers in their filenames (e.g., step 4000 -> `out_p.dat003999`).

## Additional Information

For detailed information about the project, refer to the `project2022.pdf` file and the power point presentation.
