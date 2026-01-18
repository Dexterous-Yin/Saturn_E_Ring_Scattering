# Code for "How Saturn's tenuous E ring sculpts its electron radiation belts"

This repository contains the simulation codes used in the study:
> **How Saturn's tenuous E ring sculpts its electron radiation belts**  
> Authors: Ze-Fan Yin et al.

---
## 📂 Overview
This project investigates the interaction between Saturn's E-ring ice grains and energetic electrons, using a quantitative model chain:
1. **Geant4 Simulations**: Modeling single-grain scattering effects.
2. **Monte Carlo Simulations**: Calculating the collective scattering effects (Green's function).
3. **Fokker-Planck Solver**: Modeling the equilibrium state between radial diffusion and E-ring scattering.

## ⚙️ Prerequisites
To run these codes, you will need:
* [**Geant4**](http://geant4.web.cern.ch/) (Toolkit for the simulation of the passage of particles through matter)  - *Required only if you wish to re-run single-grain interaction simulations.*
* **MATLAB** (Tested on version R2024b) - *Required for Monte Carlo simulations, Equilibrium solver, and Plotting.*
> **Note:** To reproduce the paper figures (Figs. 2-5), **only MATLAB is required**. We provide the necessary pre-calculated datasets (like Green's functions) to skip the computationally intensive Geant4 and Monte Carlo steps.

## 💼 Repository structure
```text
Saturn_E_Ring_Scattering/
│
├── Geant4_Simulation/          # Geant4 simulation core files, adapted from official B2a example
│   ├── src/                    # Source files (.cc)
│   ├── include/                # Header files (.hh)
|   ├── run1.mac & run2.mac     # Example macro file to execute the simulation
|   ├── simulation_full.ipynb   # Sample scripts for full simulation
│   └── sample_output/          # Example raw data used for Fig. 1 insets and Python notebook for preview
│
├── Monte_Carlo_Simulation/     # Monte Carlo test particle simulation for collective effect of E-ring ice grains
│   ├── calc_Green_Function_Full.m   # Main Monte Carlo logic (Requires Geant4 raw data)
│   ├── derive_equilibrium.m         # Fokker-Planck equilibrium solver
│   └── Data/                        # Contains pre-calculated Green's Function matrices
│
└── Plotting_Scripts/           # MATLAB scripts for visualization
    ├── plot_fig2_single_effect.m  # Script to reproduce Fig. 2
    ├── plot_fig3_green_function.m # Script to reproduce Fig. 3
    ├── plot_fig4_PA_evolution.m   # Script to reproduce Fig. 4
    └── plot_fig5_equilibrium.m    # Script to reproduce Fig. 5
```

## 📊 Reproducing paper figures

You can reproduce the figures directly using the provided scripts and data. To reproduce the results, please download the pre-calculated .mat files from the Releases Page and place it in the `Monte_Carlo_Simulation/Data/` folder.
- **Fig. 2** – Run `Plotting_Scripts/plot_fig2_single_effect.m`
-  **Fig. 3** – Run `Plotting_Scripts/plot_fig3_green_function.m`
-  **Fig. 4** – Run `Plotting_Scripts/plot_fig4_PA_evolution.m`. This script includes the simulation for temporal evolution of electron distributions.
-  **Fig. 5** – Run `Plotting_Scripts/plot_fig5_equilibrium.m`
## 📈 Detailed Usage

### 1. Geant4 Simulation

Navigate to the directory and compile:
```shell
cd Geant4_Simulation
mkdir build
cd build
cmake ../
make -j 4
```
Then follow `simulation_full.ipynb` to conduct the full simulation for different electron energies and grain radii. One can also uncomment the lines
```cpp
G4double randomX = 0; //r * cos(theta);
G4double randomY = 0; //r * sin(theta);
G4double randomZ = -worldZHalfLength/2;
```
in the souce file `src/B2PrimaryGeneratorAction.cc` to introduce a uniformly distributed electron source over the grain's circular cross-section.

### 2. Monte Carlo Simulation

These steps are for re-calculating the Green's function from scratch. They require the full Geant4 raw dataset (~10,000 files), which is not hosted here due to size limits.

(1) Run `Monte_Carlo_Simulation/get_K_map.m` to obtain the map for the second adiabatic invariant K.  
(2) Run `Monte_Carlo_Simulation/calc_Green_Function_Half.m` to set up simulation grids and calculate the corresponding physical quantities.  
(3) Run `Monte_Carlo_Simulation/calc_Green_Function_Full.m` to calculate the Green's function for input at each grid. *Note that this script reads raw Geant4 outputs. You must point the file path to your local Geant4 data directory.*

### 3. Equilibrium Solver
Run `Monte_Carlo_Simulation/derive_equilibrium.m` to solve the equilibrium between E-ring scattering and radial diffusion. 

## License
This code is released under the MIT License. See the `LICENSE` file for details. 


