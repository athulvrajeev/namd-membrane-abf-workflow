mol load psf gold_normal_lipid_40A.psf pdb gold_normal_lipid_40A.pdb
set sel [atomselect top all]
set sel1 [atomselect top "resname IAUM"]
$sel set beta 0
$sel1 set beta 1
$sel writepdb gold_normal_lipid_40A_gold_restrained.pdb
quit
