# Edge Tracking in Pipe Flow

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

## ⚙️ How It Works

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

## 📂 Directory Structure

### At Runtime (Generated)
```
IC_{num}/
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

### Repository Layout
```
Edge_Tracking_Pipeflow/
│
├── utils/                   # Precompiled OpenPipeFlow tools
│   ├── bisec_sim            # Template simulation for bisection
│   ├── fwdInt_sim           # Template simulation for forward integration
│   ├── p2m.out              # Extracts real-space velocity fields
│   ├── avg_state.out        # Averages two states for bisection│
│   ├── sim_run.sh           # Run a simulation
│   ├── avg_state.sh         # Average two states
│   ├── mv_ic_ext.sh         # Extract a single state
│   ├── mv_series_ext.sh     # Extract a series of states
│   └── ...
│
├── IC_0000/                 # Initial condition directory
├── bisec.config             # Configuration file for bisection
```

---

## 🛠 Installation

### Prerequisites
- **Python**: 3.x  
- **OpenPipeFlow**: Precompiled binaries included in `utils/`
- **Dependencies**: `numpy`, `matplotlib`, 'netCDF4'

---

## 🚀 Usage

### Run the workflow
```bash
python3 bisec_algo_main.py 'Main_Dir'
```

---

## 📊 Outputs
-  **Flow fields** in FwdInt_{num}/fwsInt_{1,2}p/state_{X}.cdf.dat samples of the edge-state (if converged)
- **Bisection plots**: `bisec_sim_{n}_class.png` — classification verification.
- **Forward integration plots**:  
  - `fwdInt_dist.png` — L² distance evolution.  
  - `fwdInt_QU_states.png` — QU profiles.
- **Logs**: `status` files documenting all iterations.

---

## 📚 References
If you use this code in your research, please cite:
- *[[Your Paper / DOI here](https://doi.org/10.48550/arXiv.2505.05075)]*
- *[[Openpipeflow](10.1016/j.softx.2017.05.003)]*

---

## 📄 License
This project is licensed under the **[MIT / GPL / other]** License — see the [LICENSE](LICENSE) file for details.

---

## 📬 Contact
For questions or collaborations:
- **Anton Svirsky** — *[your-email@example.com]*  
