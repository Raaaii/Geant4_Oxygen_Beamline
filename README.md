# Geant4 Oxygen Beamline

## Overview
This repository contains a Geant4-based simulation tool for microdosimetric analysis in the context of a Zero Degree Beamline for Oxygen-16 in hadrontherapy. The simulation includes specific parameters configured in a macro file, making it easier to handle and customize.

## Features
- **Hadrontherapy Simulation:** Geant4-based tool for simulating the behavior of Oxygen-16 particles in a zero-degree beamline.
- **Microdosimetric Analysis:** Comprehensive analysis of microdosimetric parameters, including Linear Energy Transfer (LET) and dose distributions.
- **Macro File Configuration:** Easy customization and handling with a dedicated macro file containing specific parameters for the simulation.

## Structure
- **`/src`:** Source code for the Geant4 simulation tool.
- **`/macros`:** Contains the macro file with configurable parameters for the simulation.
- **`/analysis`:** Results and analysis files, including LET and dose distributions.
- **`/experimental_results`:** Experimental data and analysis results for comparison.

## Getting Started
To run the simulation:
1. Clone the repository: `git clone https://github.com/your-username/Geant4_Oxygen_Beamline.git`
2. Navigate to the source directory: `cd Geant4_Oxygen_Beamline/src`
3. Compile the code: `make`
4. Run the simulation with the desired macro file: `./your_executable your_macro_file.mac`

## Analysis
The repository includes pre-generated analysis results in the `/Expermeintal_analysis` directory. Additionally, experimental results and analysis are provided for comparison.



