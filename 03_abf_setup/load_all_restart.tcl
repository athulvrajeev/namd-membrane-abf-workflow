# Load all ABF restart coordinate files into VMD as frames

set psf "gold_can_membrane.psf"
set prefix "Can_gold"
set inputdir "input"
set first 0
set last 100

mol new $psf

set out [open "loaded_frames.dat" w]
puts $out "# VMD_frame  ABF_window  file"

for {set i $first} {$i <= $last} {incr i} {

    set coor "${inputdir}/${prefix}-${i}.restart.coor"

    if {[file exists $coor]} {
        mol addfile $coor waitfor all
        set frame [expr [molinfo top get numframes] - 1]
        puts $out "$frame $i $coor"
        puts "Loaded window $i as VMD frame $frame"
    } else {
        puts "Missing: $coor"
    }
}

close $out

# Basic visualization
mol delrep 0 top

mol representation Lines
mol color Name
#mol selection "not water"
mol addrep top

mol representation VDW 1.0 12
mol color Name
mol selection "serial 65867 to 67551"
mol addrep top

animate goto 0
puts "Done. Use the frame slider to inspect all ABF windows."
puts "Frame-to-window mapping written to loaded_frames.dat"
