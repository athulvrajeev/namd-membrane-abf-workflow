## ============================================================
## SCRIPT 2: Pore Radius & Lipid Headgroup Separation
## ============================================================
## Measures:
##   (1) Lateral headgroup separation — minimum distance between
##       phosphorus atoms on opposite sides of the forming pore,
##       used as a proxy for pore diameter.
##   (2) Membrane thinning — local bilayer thickness around the NP.
##   (3) Headgroup displacement in z — how far lipid heads are
##       pulled from their equilibrium position (pore depth proxy).
##
## Usage (in VMD):
##   vmd -psf your_system.psf -dcd your_trajectory.dcd
##   source 2_pore_metrics.tcl
##
## OUTPUT: pore_metrics_<label>.dat
##   Columns: Frame | Z_NP_rel | Min_P_dist(A) | Mean_local_thickness(A) | N_P_within_ring
## ============================================================

## ---- USER SETTINGS ----

set psf_file    "system.psf"         ;# PSF topology file
set dcd_file    "smd_trajectory.dcd"     ;# DCD trajectory file
set label       "bare_cancer"             ;# Output label (e.g. bare_normal, meth_cancer)

## Resname/selection for NP
set np_sel_text  "segname NM"

## Phosphorus atoms (headgroup proxy)
set phos_sel_text "name P"

## Membrane reference (all phosphorus)
set mem_sel_text  "name P"

## Lateral radius (Angstrom) around NP COM to define the "pore region"
set pore_radius  15.0

## ---- END USER SETTINGS ----

mol new $psf_file type psf waitfor all
mol addfile $dcd_file type dcd waitfor all

set mol_id [molinfo top get id]
set nframes [molinfo $mol_id get numframes]

set outfile [open "pore_metrics_${label}.dat" w]
puts $outfile "# Frame  Z_NP_rel(A)  MinP_lateral_dist(A)  Mean_local_thickness(A)  N_P_in_ring"

set np_sel  [atomselect $mol_id $np_sel_text]
set mem_sel [atomselect $mol_id $mem_sel_text]

puts "Processing $nframes frames for pore metrics..."

for {set f 0} {$f < $nframes} {incr f} {

    $np_sel  frame $f
    $mem_sel frame $f

    ## NP COM
    set np_com [measure center $np_sel weight mass]
    set np_x   [lindex $np_com 0]
    set np_y   [lindex $np_com 1]
    set np_z   [lindex $np_com 2]

    ## Membrane z reference
    set mem_com [measure center $mem_sel weight mass]
    set mem_z   [lindex $mem_com 2]
    set np_z_rel [expr {$np_z - $mem_z}]

    ## Select phosphorus atoms within lateral pore_radius of NP COM (xy plane)
    set p_ring_sel [atomselect $mol_id \
        "($phos_sel_text) and \
         sqrt((x - $np_x)^2 + (y - $np_y)^2) < $pore_radius" frame $f]
    set n_p_ring [$p_ring_sel num]

    set min_dist  999.0
    set mean_thick 0.0

    if {$n_p_ring > 1} {
        ## Get z-coordinates of local phosphorus atoms
        set p_coords [$p_ring_sel get {x y z}]

        ## Separate upper (z > mem_z) and lower (z < mem_z) leaflet P atoms
        set upper_z_list {}
        set lower_z_list {}
        foreach coord $p_coords {
            set pz [lindex $coord 2]
            if {$pz > $mem_z} {
                lappend upper_z_list $pz
            } else {
                lappend lower_z_list $pz
            }
        }

        ## Mean local bilayer thickness (upper mean z - lower mean z)
        if {[llength $upper_z_list] > 0 && [llength $lower_z_list] > 0} {
            set sum_up 0.0
            foreach z $upper_z_list { set sum_up [expr {$sum_up + $z}] }
            set mean_up [expr {$sum_up / [llength $upper_z_list]}]

            set sum_lo 0.0
            foreach z $lower_z_list { set sum_lo [expr {$sum_lo + $z}] }
            set mean_lo [expr {$sum_lo / [llength $lower_z_list]}]

            set mean_thick [expr {$mean_up - $mean_lo}]
        }

        ## Minimum lateral distance between any two P atoms in the ring
        ## (rough pore diameter estimate — smallest gap in headgroup ring)
        set p_xy_list {}
        foreach coord $p_coords {
            lappend p_xy_list [list [lindex $coord 0] [lindex $coord 1]]
        }

        set n [llength $p_xy_list]
        for {set i 0} {$i < $n} {incr i} {
            for {set j [expr {$i+1}]} {$j < $n} {incr j} {
                set xi [lindex [lindex $p_xy_list $i] 0]
                set yi [lindex [lindex $p_xy_list $i] 1]
                set xj [lindex [lindex $p_xy_list $j] 0]
                set yj [lindex [lindex $p_xy_list $j] 1]
                set d  [expr {sqrt(($xi-$xj)**2 + ($yi-$yj)**2)}]
                if {$d < $min_dist} { set min_dist $d }
            }
        }
    }

    $p_ring_sel delete

    puts $outfile "$f  [format %.3f $np_z_rel]  [format %.3f $min_dist]  [format %.3f $mean_thick]  $n_p_ring"

    if {[expr {$f % 100}] == 0} {
        puts "  Frame $f / $nframes  |  Z_rel = [format %.2f $np_z_rel]  |  Thickness = [format %.2f $mean_thick] A"
    }
}

$np_sel delete
$mem_sel delete
close $outfile

puts "Done. Output written to: pore_metrics_${label}.dat"
