#!/bin/bash

# Assign the filename
filename="serialcuda_pulse.job"

for ((endtime = 2000; endtime < 4000; endtime+=100))
do
        let a=$endtime+100
        if [[ $endtime != "" && $a != "" ]]; then
        sed -i "s/$endtime/$a/" $filename
        fi
	sbatch serialcuda_pulse.job
done

