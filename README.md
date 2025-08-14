# FESTIM-v2.0-review

## Overview

This repository contains Python scripts for the FESTIM v2.0 review paper.

## Contents

1. [Coupling to external solvers](#coupling-to-external-solvers)
    - [OpenFOAM](#OpenFOAM)
    - [OpenMC](#OpenMC)
2. [V&V](#V&V)
3. [How to Use](#how-to-use)


## Coupling to external solvers

The `coupling` folder includes scripts for comparing FESTIM results with other models.

### OpenFOAM

Scripts in this subfolder are used for coupling FESTIM simulations with OpenFOAM.

### OpenMC

Scripts in this subfolder are used for coupling FESTIM simulations with OpenMC.

## V&V

The `v&v` folder contains scripts for verifying the accuracy of FESTIM using the method of manufactured solutions and validating FESTIM against experimental data.

## How to Use

1. Clone this repository to your local machine.
2. Create a Conda environment using the provided `environment.yml` file:

    ```bash
    conda env create -f environment.yml
    ```

   This will set up a Conda environment named `festim-v2-review-env` with all the required dependencies for running the FESTIM scripts and visualisation scripts.

3. Activate the Conda environment:

    ```bash
    conda activate festim-v2-review-env
    ```

4. Execute the Python scripts using the activated Conda environment and ensure compatibility with FESTIM requirements.

5. Navigate to the desired folder based on the simulation or comparison you are interested in.

**Note**: OpenFOAM isn't included in this conda environment. In order to run the OpenFOAM scripts, refer to [OpenFOAM's installation instructions](https://www.openfoam.com/download/openfoam-installation-on-linux).


## Contact

For any questions or issues, please contact remidm@mit.edu.