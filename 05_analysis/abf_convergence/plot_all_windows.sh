#!/bin/bash

mkdir -p plots

for i in {2..57}
do
    count_file="win${i}/output/win${i}.count"
    pmf_file="win${i}/output/win${i}.pmf"

    if [ -f "$count_file" ]; then
        gnuplot <<EOF
set terminal png size 1000,700
set output "plots/win${i}_count.png"
set title "win${i} count"
set xlabel "Reaction coordinate"
set ylabel "Count"
plot "${count_file}" using 1:2 with lines title "win${i} count"
EOF
        echo "Created plots/win${i}_count.png"
    else
        echo "Missing count file for win${i}"
    fi

    if [ -f "$pmf_file" ]; then
        gnuplot <<EOF
set terminal png size 1000,700
set output "plots/win${i}_pmf.png"
set title "win${i} PMF"
set xlabel "Reaction coordinate"
set ylabel "PMF"
plot "${pmf_file}" using 1:2 with lines title "win${i} PMF"
EOF
        echo "Created plots/win${i}_pmf.png"
    else
        echo "Missing PMF file for win${i}"
    fi

done
