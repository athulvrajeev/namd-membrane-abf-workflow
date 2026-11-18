## ============================================================
## SCRIPT 3: Water Density Profile Along Z Near NP
## ============================================================
## Computes the number of water molecules in thin z-slabs
## centered on the NP position at each frame.
## This produces a 2D map (frame vs z-slab) showing water
## co-transport and pore formation as the NP moves through.
##
## Usage (in VMD):
##   vmd -psf your_system.psf -dcd your_trajectory.dcd
##   source 3_water_density_profile.tcl
##
## OUTPUT:
##   water_density_profile_<label>.dat  — full 2D data (frame x z-slab)
##   water_zprofile_summary_<label>.dat — mean water count per z-slab
## ============================================================

## ---- USER SETTINGS ----

set psf_file    "system.psf"         ;# PSF topology file
set dcd_file    "smd_trajectory.dcd"     ;# DCD trajectory file
set label       "bare_cancer"             ;# Output label (e.g. bare_normal, meth_cancer)

set np_sel_text  "segname NM"
set wat_sel_text "water and name OH2"
set mem_sel_text "name P"

## Z-slab parameters (relative to NP COM)
set z_min   -30.0     ;# Angstrom below NP COM
set z_max    30.0     ;# Angstrom above NP COM
set dz        2.0     ;# Slab thickness

## Lateral cylinder radius around NP COM for water counting
set cyl_radius  8.0   ;# Angstrom

## ---- END USER SETTINGS ----

mol new $psf_file type psf waitfor all
mol addfile $dcd_file type dcd waitfor all

set mol_id [molinfo top get id]
set nframes [molinfo $mol_id get numframes]

## Build z-slab bin centers
set bins {}
set z $z_min
while {$z <= $z_max} {
    lappend bins $z
    set z [expr {$z + $dz}]
}
set nbins [llength $bins]

## Open output files
set outfile  [open "water_density_profile_${label}.dat" w]
set sumfile  [open "water_zprofile_summary_${label}.dat" w]

## Write header for 2D file
set header "# Frame  Z_NP_rel"
foreach b $bins { append header "  z=[format %.1f $b]" }
puts $outfile $header

## Accumulators for summary
set sum_counts [lrepeat $nbins 0.0]

set np_sel  [atomselect $mol_id $np_sel_text]
set mem_sel [atomselect $mol_id $mem_sel_text]

puts "Processing $nframes frames for water density profile..."

for {set f 0} {$f < $nframes} {incr f} {

    $np_sel  frame $f
    $mem_sel frame $f

    set np_com  [measure center $np_sel weight mass]
    set np_x    [lindex $np_com 0]
    set np_y    [lindex $np_com 1]
    set np_z    [lindex $np_com 2]

    set mem_com [measure center $mem_sel weight mass]
    set mem_z   [lindex $mem_com 2]
    set np_z_rel [expr {$np_z - $mem_z}]

    set row "$f  [format %.3f $np_z_rel]"
    set bin_idx 0

    foreach bz $bins {
        ## Absolute z slab center
        set slab_z_abs [expr {$np_z + $bz}]
        set slab_lo    [expr {$slab_z_abs - $dz/2.0}]
        set slab_hi    [expr {$slab_z_abs + $dz/2.0}]

        ## Count waters in cylinder around NP COM, within this z-slab
        set wat_sel [atomselect $mol_id \
            "($wat_sel_text) and \
             z > $slab_lo and z < $slab_hi and \
             sqrt((x - $np_x)^2 + (y - $np_y)^2) < $cyl_radius" frame $f]
        set n_wat [$wat_sel num]
        $wat_sel delete

        append row "  $n_wat"

        ## Accumulate for summary
        lset sum_counts $bin_idx [expr {[lindex $sum_counts $bin_idx] + $n_wat}]
        incr bin_idx
    }

    puts $outfile $row

    if {[expr {$f % 100}] == 0} {
        puts "  Frame $f / $nframes  |  NP Z_rel = [format %.2f $np_z_rel] A"
    }
}

$np_sel delete
$mem_sel delete
close $outfile

## Write summary: mean water count per z-slab across all frames
puts $sumfile "# Z_slab_relative_to_NP(A)   Mean_Water_Count"
set bin_idx 0
foreach bz $bins {
    set mean_w [expr {[lindex $sum_counts $bin_idx] / double($nframes)}]
    puts $sumfile "[format %.1f $bz]   [format %.3f $mean_w]"
    incr bin_idx
}
close $sumfile

puts "Done."
puts "  2D profile : water_density_profile_${label}.dat"
puts "  Z-summary  : water_zprofile_summary_${label}.dat"
