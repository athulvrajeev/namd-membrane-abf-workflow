import MDAnalysis as mda
from MDAnalysis.analysis import density
import matplotlib.pyplot as plt
import numpy as np

# Load the universe
u = mda.Universe('../gold_cancer_lipid_40A.psf', '../output_npt1/npt1.dcd')

lipid_types = ['POPS', 'POPE', 'POPC', 'PSM', 'CHL1']
density_profiles = {}

for lipid in lipid_types:
    lipid_atoms = u.select_atoms(f'resname {lipid}')
    print(f"{lipid}: {lipid_atoms.n_atoms} atoms selected")  # Debug selection

    if lipid_atoms.n_atoms > 0:
        D = density.DensityAnalysis(lipid_atoms, delta=1.0)
        D.run()
        print(f"Density calculated for {lipid}")  # Debug density calculation
        
        density_grid = D.results.density.grid
        if density_grid.ndim == 3:
            z_density = np.mean(density_grid, axis=(0, 1))
            z_coords = D.results.density.edges[2][:-1]
            density_profiles[lipid] = (z_coords, z_density)
            print(f"Density profile stored for {lipid}")
        else:
            print(f"Unexpected density grid dimensions for {lipid}: {density_grid.ndim}")
    else:
        print(f"No atoms found for lipid {lipid}. Skipping.")
        
#Set Times New Roman font
plt.rcParams['font.family'] = 'Times New Roman'

# Plot all density profiles
plt.figure(figsize=(10, 6))
for lipid, (z_coords, z_density) in density_profiles.items():
    plt.plot(z_coords, z_density, label=f'{lipid}')
plt.xlabel('Z-coordinate (Å)')
plt.ylabel('Density (kg/m³)')
plt.legend()
plt.title('Density Profiles of Lipid Bilayer along z-axis')
plt.savefig('debug_density_profiles_by_lipid_type_dcd.png', dpi=600)
plt.show()

# Save individual density profiles to CSV
for lipid, (z_coords, z_density) in density_profiles.items():
    csv_filename = f'{lipid}_density_profile.csv'
    try:
        np.savetxt(csv_filename, 
                   np.column_stack((z_coords, z_density)), 
                   delimiter=',', 
                   header='Z-coordinate,Density', 
                   comments='')
        print(f"CSV file for {lipid} saved as: {csv_filename}")
    except Exception as e:
        print(f"Error saving CSV file for {lipid}: {e}")


