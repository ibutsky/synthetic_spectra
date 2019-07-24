#### change these parameters to specify which spectra to analyze #####
unanalyzed_spectra_dir="../../data/unanalyzed_spectra"
analyzed_spectra_dir="../../data/analyzed_spectra"
model="tempest"
model="P0"
# Redshift at which spectra was generated (or estimate) 
redshift=0.25
# ray_id is the identifying number of the ray (the last number before the .fits)
ray_id=10


#### Don't change below this line, unless you're sure you want to ####
### defining the end radius based on the desired number of steps
end_ray_id=$[$ray_id+$nsteps]
# go into the working directory
cd ../../data/analyzed_spectra
# defining the base name that many files share
basename="COS-FUV_"$model$"_z"$redshift"_"$ray_id
# navigate to spectrum-specific directory
cd $basename
open -a Preview FitInspection.pdf &    
# run pyigm_igmguesses 
pyigm_igmguesses $basename"_ibnorm.fits" -p $basename"_lineids.json" -o $basename"_lineids.json"
# fix the zsys in the generated .joebvp files
python ../../../scripts/analysis/fix_joebvp_files.py $redshift
# list all .joebvp files in a file named "flist"
rm flist
ls *.joebvp >>flist

cd ../../../scripts/analysis



