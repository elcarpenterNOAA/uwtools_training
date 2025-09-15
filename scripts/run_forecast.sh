#!/bin/bash


config=$1
cycle=$2

module load nco

set -eux

# Prepare run directory
uw fv3 files_copied --cycle $cycle -c $config --key-path cycled_forecast --batch

rundir=$(echo 'cycle: !datetime '$cycle | uw config realize -i $config --output-format yaml \
--key-path cycled_forecast.fv3.rundir --update-format yaml)

for fn in $rundir/INPUT/*.nc ; do
  ncatted -a checksum,,d,, $fn
done
ncatted -O -a source,global,c,c,'FV3GFS GAUSSIAN NETCDF FILE' $rundir/INPUT/fv_core.res.tile1.nc

uw fv3 run --cycle $cycle -c $config --key-path cycled_forecast --batch
