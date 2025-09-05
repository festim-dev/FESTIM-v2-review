# FESTIM-v2.0-review

## Overview

This repository contains Python scripts for the FESTIM v2.0 review paper.

## Contents

1. [How to Use](#how-to-use) 
2. [Advection demo](#advection-demo)
3. [Performance testing](#performance-testing)
4. [Multiphysics coupling](#multi-physics-coupling)
    - [OpenFOAM](#OpenFOAM)
    - [OpenMC](#OpenMC)

## How to Use

1. Clone this repository to your local machine.
2. Create a Conda environment using the provided `environment.yml` file:

    ```bash
    conda env create -f environment_festim_v2.yml
    ```

   This will set up a Conda environment named `festim-v2-review-env` with all the required dependencies for running the FESTIM scripts and visualisation scripts.


3. Activate the Conda environment:

    ```bash
    conda activate festim-v2-review-env
    ```

4. Execute the Python scripts using the activated Conda environment and ensure compatibility with FESTIM requirements.

5. Navigate to the desired folder based on the simulation or comparison you are interested in.


## Advection demo

The `advection_demo` folder contains scripts demonstrating FESTIM v2.0's new advection capabilities for modeling hydrogen transport in flowing media. 
The demo showcases how to couple diffusion and advection processes, comparing results with pure diffusion cases to highlight the effects of velocity fields on hydrogen concentration distributions.

Running the test file:

```bash 
python advection_demo.py
```


## Performance testing

The `performance` folder includes scripts testing the computational performance of FESTIM v1.4 and v2.0 on a 2D multi-material diffusion case adapted from the [V&V book discontinuity verification case](https://festim-vv-report.readthedocs.io/en/latest/verification/mms/discontinuity.html). 
The test case runs a transient simulation on a 500×500 mesh to benchmark different numerical methods and compare runtimes between FESTIM versions. 

>[!NOTE]
>To run all the scripts in the [Performance testing](#performance-testing) section, a conda environment for festim v1.4 is also required, which can be created by:
>
>```bash
>conda env create -f environment_festim_v2.yml
>```
>
>Then activate it the same way as done previously.

Running the performance tests for festim v1 int he conda environemnt (`festim-v1-review-env`):

```bash 
python test_performance_v1.py
```

Running the performance tests for festim v2 int he conda environemnt (`festim-v2-review-env`):

```bash 
python test_performance_v2.py
```

Before finally plotting the results:
```bash 
python plot_runtimes.py
```

## Multi-physics coupling

The `coupling` folder includes scripts which demonstrate how external solvers such as OpenFOAM and OpenMC can be coupled with FESTIM.

### OpenFOAM

Scripts in the `coupling_cfd` subfolder are used for coupling FESTIM simulations with OpenFOAM.

>[!NOTE]
>OpenFOAM isn't included in this conda environment. In order to re-run the OpenFOAM scripts, refer to [OpenFOAM's installation instructions](https://www.openfoam.com/download/openfoam-installation-on-linux).


### OpenMC

Scripts in the `coupling_neutronics` subfolder are used for coupling FESTIM simulations with OpenMC.


## Contact

For any questions or issues, please contact darkj385@mit.edu