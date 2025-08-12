# Split-Edge Tracking in Pipe Flow

## Overview
This repository contains a **Python workflow** for **edge tracking in pipe flow** using **precompiled [OpenPipeFlow](http://openpipeflow.org/) simulations**.  
The code implements a **bisection procedure** between two supplied initial conditions to locate the edge-states on the **edge** separating trajectories that develope into a single puff from those that evolve into two puffs.

It is designed to:
- Automate the setup and execution of bisection steps.
- Classify initial conditions into single- or two-puff outcomes.
- Perform forward integration.
- Generate diagnostic plots for manual inspection.

---

## Features
- **Bisection procedure** to locate edge states between given initial conditions.
- Automated **classification** into single- or two-puff outcomes.
- **Forward integration** for dynamical evolution of the bounding states.
- Precompiled OpenPipeFlow simulation templates for fast execution.
- Plot generation for manual classification checks and L² divergence analysis.

---

## How It Works

### 1. Initial Conditions
- Place the two starting states in the `Main_Dir/IC_0000/` directory:
  - `s_1p.cdf.dat` — evolves into a single puff.
  - `s_2p.cdf.dat` — evolves into two puffs.
- Edit `Main_Dir/bisec.config` in the main directory to set parameters for the bisection (e.g., tolerance, iteration limits).

### 2. Bisection Procedure
- The code bisects between the bounding states until the desired state separation is reached.
- Each iteration:
  - Produces a new initial condition file: `s_{X}p_{n}.cdf.dat` where:
    - `X = 1` → single puff outcome.
    - `X = 2` → two puff outcome.
    - `n` → bisection iteration number.
  - Saves a classification figure: `bisec_sim_{n}_class.png` for manual verification.
- Progress is logged in a `status` file.

### 3. Forward Integration
- After bisection, selected states are integrated forward for a fixed duration `T`:
  - Creates two simulations: `fwdInt_1p/` and `fwdInt_2p/`.
  - Produces:
    - Classification figures (as in bisection).
    - `fwdInt_dist.png` — L² distance between states vs time.
    - `fwdInt_QU_states.png` — QU profiles before and after integration.
- Results are recorded in a `status` file, indicating the states chosen for the next bisection step.

---

## Directory Structure

### Directory Layout
```
├── bisec_algo_main.py # Main biesction script
├── bisec_algo_aux.py  # Functions used during the bisection algorithm 
├── opf_analysis.py    # Analysis scripts for pipe flow data
├── bisection.f90      # Fortran code for the avg_state.out utility
├── Bisection_run/
  │
  ├── utils/  
  │   ├── bisec_sim/           # Template simulation for bisection
  │   ├── fwdInt_sim/          # Template simulation for forward integration
  │   ├── p2m.out              # Extracts real-space velocity fields (need to compile seperetly)
  │   ├── avg_state.out        # Averages two states for bisection (need to compile seperetly)
  │   ├── sim_run.sh           # Run a simulation
  │   ├── avg_state.sh         # Average two states
  │   ├── mv_ic_ext.sh         # Extract a single state
  │   └── mv_series_ext.sh     # Extract a series of states
  │   
  │
  ├── IC_0000/                 # Initial condition directory
  ├── bisec.config             # Configuration file for bisection
  └── main.state               # Traker of the current bisection step
```

---

### At Runtime (Generated)
```
├── Bisection_run/
    IC_{n}/
       ├── s_1p.cdf.dat          # 1-puff initial condition
       └── s_2p.cdf.dat          # 2-puff initial condition
    
    Bisec_{num}/
       ├── bisec_sim_{n}         # n-th bisection integration
       ├── status                # Log of bisection iterations
       ├── s_{X}p_{n}.cdf.dat    # Intermediate states (X = 1 or 2 puff outcome)
       ├── bisec_sim_{n}_class.png  # Classification plot
       └── ...                      # Other outputs
    
    FWdInt_{num}/
       ├── fwdInt_1p/             # Forward integration for 1 puff
       ├── fwdInt_2p/             # Forward integration for 2 puffs
       ├── fwdInt_dist.png        # L² distance evolution
       ├── fwdInt_QU_states.png   # QU profiles before/after integration
       ├── classification figures
       └── status                 # Forward integration log

```


## Installation

### Prerequisites
- **Python**: 3.x  
- **OpenPipeFlow**: Precompiled binaries included in `utils/`
- **Dependencies**: `numpy`, `matplotlib`, `netCDF4`

---

## Usage

### Configuration File Reference — `bisec.config`
This file controls the **bisection** and **forward integration** procedures, as well as state classification and data conversion.  
It is divided into several sections:

#### [State distance]
Defines thresholds for state separation during bisection and forward integration.
- **`L_bs`** — Minimum $L^2$ distance between bounding states required to stop the bisection process.  
- **`L_fwd`** — Maximum $L^2$ distance allowed for the forward integration step before it is considered divergent.  
- **`force_change`** — Boolean flag to force a change in classification between iterations. If `True` ensures that bisection is performed untill both `1p` and `2p` states have been replaced.

#### [State definition]
Parameters used to classify states into **puff types** based on flow features.
- **`F_th`** — Threshold $q_{th}$ for $q$ variable consided as turbulent detection in the flow field.  
- **`min_gap`** —  Spatial separation between detected structures $w$ to be considered distinct (in '2p' state). If $w>$`min_gap` the state will be considered as '2p' (if tuebulent lenfgth condition is also satesfied)  
- **`small_gap`** — Spatial seperation between structures $w$ to not be considered a single puff. Thus only if gap $w<$`small_gap` the structure can be classified as '1p'. (This parameter is not used in practics, by seting 'small_gap=0'). 
- **`p0_Ft_lim`, `p1_Ft_lim`, `p2_Ft_lim`, `p3_Ft_lim`** — Ranges of turbulent length $l_{turb} defining classification categories for 0 - unclassided, 1 - one puff, 2 - two puffs, or 3 - three (or more) puffs.  
- **`min_class_time`** — Minimum simulation time required before a classification is accepted.

#### [Conversion settings]
Parameters for converting simulation data into physical fields for analysis or visualization.
- **`Lz`** — Domain length in the axial ($z$) direction.  
- **`nx`, `ny`, `nz`** — Number of grid points in the $x$, $y$, and axial direction $z$ respectively to which the data is interpolated for analysis.

#### [Misc]
Miscellaneous runtime controls.
- **`c_timer_min`** — Minimum allowed wall-clock time (in minutes) before considering a simulation step complete. The script will wait repeatedly checks every `c_timer_min` if the simulations where complete.  

### Bisection state file `main.state`:
Keeps track of the current state of the bisection:
### [State]
- `s` - The name of the bisection step. `IC` - analyzing inital conditions, `Bis` - bisection step, `Fwd` - forward  integration, `Adv` - advancing to the next step. When initating should be `IC`.
- `n` - The number of the current iteration. 


### Preperaing the bisection simulation
# Setup Instructions

1. **Clone the repository**  
2. **Compile two [OpenPipeFlow](http://openpipeflow.org/) template simulations:**  
   - **`bisec_sim/`** — Long simulation time limit $T_{\text{bisec}}$ with a low saving rate.  
     - If $T_{\text{bisec}}$ is too short to achieve a classification, another simulation will be initiated with the initial conditions taken from the last snapshot.  
     - This process repeats until a classification is achieved.  
     - The saving rate should be kept relatively low to reduce memory usage.  
   - **`fwdInt_sim/`** — Short simulation time limit $T_{\text{fwd}}$ with a high saving rate.  
     - $T_{\text{fwd}}$ sets the maximum integration time for the bounding trajectories and should be slightly longer than the expected separation time.  
     - The high saving rate ensures that the states can be used for spatio-temporal analysis of the edge states.
3. **Compile two utility programs:**  
   - **`p2m.out`** — Built-in openpiflow utility `prim2matlab.f90` for data extraction.  
   - **`avg_state.out`** — Custom utility for generating a bisection state between two input states, using the provided `bisection.f90` OpenPipeFlow utility code.
4. **Set the appropriate parameters** in `bisec.config`.
5. **Create the initial directory** `IC_0000/` containing:  
   - `s_1p.cdf.dat` — Flow state known to be (or evolve into) a one-puff state.  
   - `s_2p.cdf.dat` — Flow state known to be (or evolve into) a two-puff state.
6. **Initialize the state file** with:  
   ```
   s = IC
   n = 0
   ```

### Run the workflow
```bash
python3 bisec_algo_main.py '/Bisection_run/'
```

---

## Outputs
- **Flow fields**: ``Bisection_run/FwdInt_{num}/fwsInt_{1,2}p/state_{X}.cdf.dat`` samples of the edge-state (when converged)
- **Analysis plots**: For manual verification.
- **Forward integration plots**:  
- **Logs**: `status` files documenting all iterations.

---

## References
If you use this code in your research, please cite:
- *[[	Self-Replication of Turbulent Puffs: On the edge between chaotic saddles](https://doi.org/10.48550/arXiv.2505.05075)]*
- *[[Openpipeflow](10.1016/j.softx.2017.05.003)]*

---

## Contact
- **Anton Svirsky** — *[anton.sv@campus.technion.ac.il]*  
