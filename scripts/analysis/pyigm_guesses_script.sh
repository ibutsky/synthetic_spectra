#### change these parameters to specify which spectra to analyze #####
# Simulation time at which spectra was generated
time=11.2
# Redshift at which spectra was generated (or estimate)
redshift=0.2
# Name of model: 'stream' or 'anisd'
model="stream"
# Orientation of the disk (relative to observer/line-of-sight)
# Options: 'face', 'edge_theta0', 'edge_theta1.0', 'edge_theta1.5'
orient="edge_theta1.5"
# The impact paramter; 10 - 100, in intervals of 10
radius=10
# The desired number of impact parameters to analyze in this script
nsteps=10



#### Don't change below this line, unless you're sure you want to ####
### defining the end radius based on the desired number of steps
end_radius=$[$radius+$[10*$nsteps]]
# go into the working directory
cd ../../data/analyzed_spectra
while [ $radius -lt $end_radius ]; do
#    python stitch_fits.py $orient $model $time $radius
    # defining the base name that many files share
    basename="COS-FUV_"$orient$"_"$model"_"$time"Gyr_r"$radius"kpc"
    # make the directory where we'll store all relevant data for this spectrum
    mkdir $basename
    # navigate to spectrum-specific directory
    cd $basename
    # move unanalyzed fits files to this directory
    mv "../../unanalyzed_spectra/"$basename.fits .
    # make info file
    python ../../../scripts/analysis/create_info_file.py $orient $model $time $radius
    # set continuum fit 
    lt_continuumfit --redshift $redshift $basename".fits" $basename"_ibnorm.fits"
    # run pyigm_igmguesses 
    pyigm_igmguesses $basename"_ibnorm.fits" -o $basename"_lineids.json"
    pyigm_igmguesses $basename"_ibnorm.fits" -p $basename"_lineids.json" -o $basename"_lineids.json"
    # list all .joebvp files in a file named "flist"
    ls *.joebvp >>flist
    # done here.  move up one directory
    cd ../
    # Next radius.
    radius=$[$radius+10]
done
# go back to the directory where we started
cd ../../scripts/analysis

