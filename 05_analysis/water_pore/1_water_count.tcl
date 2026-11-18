# ============================================================
# Water co-transport in membrane core around NP translocation axis
# ============================================================

set psf_file    "system.psf"         ;# PSF topology file
set dcd_file    "smd_trajectory.dcd"     ;# DCD trajectory file
set label       "bare_cancer"             ;# Output label (e.g. bare_normal, meth_cancer)

set np_sel_text   "segname NM"          ;# change if needed
set wat_sel_text  "water and name OH2"
set mem_sel_text  "name P"

set cyl_radius    6.0                   ;# radial distance from NP axis, Angstrom
set core_half_z   15.0                  ;# membrane core: +/- 15 A from membrane COM

mol new $psf_file type psf waitfor all
mol addfile $dcd_file type dcd waitfor all

set mol_id [molinfo top get id]
set nframes [molinfo $mol_id get numframes]

set outfile [open "water_core_count_${label}.dat" w]
puts $outfile "# Frame  NP_Z_rel_A  Water_count_core_cylinder"

set np_sel  [atomselect $mol_id $np_sel_text]
set mem_sel [atomselect $mol_id $mem_sel_text]

for {set f 0} {$f < $nframes} {incr f} {

    $np_sel frame $f
    $mem_sel frame $f

    set np_com  [measure center $np_sel weight mass]
    set mem_com [measure center $mem_sel weight mass]

    set np_x [lindex $np_com 0]
    set np_y [lindex $np_com 1]
    set np_z [lindex $np_com 2]
    set mem_z [lindex $mem_com 2]

    set np_z_rel [expr {$np_z - $mem_z}]
    set z_low  [expr {$mem_z - $core_half_z}]
    set z_high [expr {$mem_z + $core_half_z}]

    set water_sel [atomselect $mol_id "
        ($wat_sel_text) and
        z > $z_low and z < $z_high and
        sqrt((x - $np_x)*(x - $np_x) + (y - $np_y)*(y - $np_y)) < $cyl_radius
    " frame $f]

    set nwater [$water_sel num]
    $water_sel delete

    puts $outfile "$f  [format %.3f $np_z_rel]  $nwater"

    if {[expr {$f % 100}] == 0} {
        puts "Frame $f / $nframes : Zrel = [format %.2f $np_z_rel] A, waters = $nwater"
    }
}

$np_sel delete
$mem_sel delete
close $outfile

puts "Done: water_core_count_${label}.dat"
