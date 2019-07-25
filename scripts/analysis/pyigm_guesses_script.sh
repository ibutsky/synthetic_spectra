#### change these parameters to specify which spectra to analyze #####
# Name of model: 'P0' or 'tempest'
unanalyzed_spectra_dir="../../data/unanalyzed_spectra"
analyzed_spectra_dir="../../data/analyzed_spectra"
model="tempest"
#model="P0"
# Redshift at which spectra was generated (or estimate) 
redshift=0.20
# ray_id is the identifying number of the ray (the last number before the .fits)
ray_id=2
# The desired number of impact parameters to analyze in this script
nsteps=1


#### Don't change below this line, unless you're sure you want to ####
### defining the end radius based on the desired number of steps
end_ray_id=$[$ray_id+$nsteps]
# go into the working directory
#cd ../../data/analyzed_spectra
while [ $ray_id -lt $end_ray_id ]; do
    # defining the base name that many files share
    basename="COS-FUV_"$model$"_z"$redshift"_"$ray_id
    # make the directory where we'll store all relevant data for this spectrum
    mkdir $basename
    # navigate to spectrum-specific directory
    cd $basename
    # move unanalyzed fits files to this directory
    mv "../"$unanalyzed_spectra_dir"/"$basename.fits  $basename.fits
    # set continuum fit 
    lt_continuumfit --redshift $redshift $basename".fits" $basename"_ibnorm.fits"
    # run pyigm_igmguesses 
    pyigm_igmguesses $basename"_ibnorm.fits" -o $basename"_lineids.json"
    pyigm_igmguesses $basename"_ibnorm.fits" -p $basename"_lineids.json" -o $basename"_lineids.json"
    # list all .joebvp files in a file named "flist"
    ls *.joebvp >>flist
    # Next radius.
    cd ../
    mv $basename $analyzed_spectra_dir"/."
    ray_id=$[$ray_id+1]
done
# go back to the directory where we started


