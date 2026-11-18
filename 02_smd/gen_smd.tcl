mol load psf gold_normal_lipid_40A.psf pdb gold_normal_lipid_40A.pdb
set sel [atomselect top all]
set sel1 [atomselect top "resname IAUM"]
$sel set occupancy 0
$sel1 set occupancy 1
$sel writepdb gold_normal_lipid_40A_smd.pdb
quit
