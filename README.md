# NAMD workflow for nanoparticle–membrane permeation free energies (SMD + ABF)

Scripts and input templates for computing the free-energy profile (potential of mean force, PMF) of a gold nanoparticle crossing a lipid bilayer, using **NAMD**, **VMD**, the **Colvars** module and **Python**. This is the workflow behind:

> Rajeev, A. V., Singh, G. K., Joshi, H., & Unni, H. N. (2026). Comparing the translocation of the bare and surface-functionalised gold nanoparticle across normal and cancer-mimicking cell membranes. *Physical Chemistry Chemical Physics*, 28, 19489–19503. https://doi.org/10.1039/D6CP01177B

The repository holds only scripts and templates. Trajectories, structures (PSF/PDB) and force-field files are **not** included.

## Workflow

```
01_equilibration  →  02_smd  →  03_abf_setup  →  04_hpc_run  →  05_analysis
 staged NVT/NPT       pull NP     cut SMD path      run each        convergence, PMF,
 with restraints      along z     into windows      ABF window      membrane & water metrics
```

| Folder | What it contains |
|---|---|
| `01_equilibration/` | NAMD inputs for staged NVT → NPT equilibration (`nvt1`, `npt1–3`, `npt_restart`). `gen_rest_*.tcl` (VMD) write reference PDBs that flag the nanoparticle and/or lipid atoms for harmonic positional restraints. |
| `02_smd/` | `gen_smd.tcl` marks nanoparticle atoms as SMD atoms; `smd.namd` pulls the nanoparticle through the membrane along −z at constant velocity; `smd_merge.sh` and `smd_analysis.tcl` extract the force–position profile from the NAMD log. |
| `03_abf_setup/` | `get_windows.tcl` (VMD + pbctools) reads the SMD trajectory, measures the nanoparticle–membrane z-distance per frame and writes starting coordinates and box files for each ABF window. `config_replicate.sh` / `colvar_replicate.sh` generate one NAMD config and one Colvars file per window; `*_move.sh` and `restart_files_move.sh` place them into `win<i>/` folders. `load.tcl` checks the nanoparticle z-position of every window's starting structure and `load_all_restart.tcl` loads them into VMD for a visual check. `merge/` holds the NAMD/Colvars inputs that merge the gradients and counts of all windows into a single PMF (Colvars `inputPrefix`). `templates/` has one example window config and Colvars file (`distanceZ` collective variable, ABF bias with harmonic walls). |
| `04_hpc_run/` | `run_monitor.sh` runs windows on a GPU workstation and starts `monitor_copy.py`, which saves time-stamped copies of the ABF output (`.grad`, `.count`, `.pmf`, Colvars state) every 1 ns so convergence can be checked afterwards. `slurm_gpu_loop.sh` / `slurm_cpu_mpi_loop.sh` run the windows on a SLURM cluster. |
| `05_analysis/abf_convergence/` | Collect time-resolved ABF files and plot gradient, sample-count and PMF convergence against simulation time for each system. |
| `05_analysis/pmf/` | Shift PMFs to a bulk-water reference and plot PMFs for different membranes/nanoparticles together. |
| `05_analysis/membrane/` | Lipid density profiles (MDAnalysis), area per lipid with block averaging, and deuterium order parameters (S<sub>CD</sub>) for the sn-1/sn-2 chains. |
| `05_analysis/water_pore/` | VMD scripts for water co-transport, pore/headgroup metrics and local membrane thinning around the nanoparticle, plus a Python plotting script. |
| `05_analysis/smd/` | Plot SMD force against nanoparticle position. |

## Requirements

- NAMD 2.14 or 3.0 (with Colvars)
- VMD 1.9.x (with the `pbctools` plugin)
- Python 3.9+ with `numpy`, `matplotlib`, `pandas`, `MDAnalysis` (`pip install -r requirements.txt`)
- CHARMM36 topology/parameter files (`toppar`), plus INTERFACE force-field parameters for gold. These are not redistributed here.
- Optional: SLURM for cluster runs, gnuplot for `plot_all_windows.sh`

## Using the scripts

1. **Force-field paths.** The NAMD inputs point to `/path/to/toppar/`. Replace this with the location of your CHARMM36 files, for example:
   `sed -i 's#/path/to/toppar/#/your/toppar/#' 01_equilibration/*.namd 02_smd/*.namd 03_abf_setup/*.sh 03_abf_setup/templates/*.conf`
2. **System-specific selections.** Atom serial numbers in `get_windows.tcl` and the Colvars files, the residue name of the nanoparticle (`IAUM`) and lipid names in the analysis scripts match the systems in the paper. Change them for a new system.
3. **Window ranges.** Loops such as `for ((i=2; i<=61; i++))` define which windows are generated or run. Edit them to match your window count.
4. **Run order.** Equilibration → SMD → `get_windows.tcl` → `config_replicate.sh` + `colvar_replicate.sh` → move scripts → run windows → `collect_grad_count_merge.sh` → merge all windows with `merge/merge.conf` → analysis.

Scripts that read files from fixed folders (for example the `SYSTEMS` dictionary in `convergence_count_grad_pmf.py`) use paths of the form `data/ABF/<system>/`. Point them at your own output folders.

## Third-party software

This repository contains only scripts and input files written for this project. NAMD, VMD, the Colvars module and the CHARMM36/INTERFACE force-field files are **not** included. Obtain them from their developers and use them under their own licences. The MIT licence below applies only to the files in this repository.

## Acknowledgements

The NAMD configuration files for equilibration and steered MD (`01_equilibration/*.namd`, `02_smd/smd.namd`) were adapted from templates provided by Dr. Himanshu Joshi (IIT Hyderabad). The ABF window and merge configurations (`03_abf_setup/templates/win0.conf`, `03_abf_setup/merge/merge.conf`, and the configs generated by `config_replicate.sh`) follow the layout of the NAMD tutorial configuration files from the Theoretical and Computational Biophysics Group, University of Illinois Urbana-Champaign. This work was carried out at IIT Hyderabad under the supervision of Dr. Harikrishnan N. Unni.

## Citation

If you use these scripts, please cite the paper above (see `CITATION.cff`).

## Author

Athul V. Rajeev, PhD scholar, Department of Biomedical Engineering, Indian Institute of Technology Hyderabad
ORCID: [0009-0008-4748-3931](https://orcid.org/0009-0008-4748-3931)

## Licence

MIT, see `LICENSE`.
