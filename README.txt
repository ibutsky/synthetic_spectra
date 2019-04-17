File Name Explanation:

1. Instrument name: COS-G130M or COS-G160M. Once stitched together, 
   the instrument name becomes COS-FUV
2. Orientation of ray: "edge" or "face"
   If the orientation is "edge" then, there's a second, "theta"
   parameter. 
3. Theta (optional) possible values: 0, pi/4(0.7), pi/3(1.0), pi/2(1.5)
   (note: 0 is along the z-axis)
2. CR model: anisotropic diffusion or streaming
3. Age of galaxy in Gyr (Right now, they're all 11 Gyr)
4. Impact parameter in kpc (range from 10 to 200 kpc)



Description of Workflow: 

1. Generate synthetic spectra from simulations. This is done on Blue Waters using 
   the script "generate_spectra.py" in scripts/analysis/. The files are then copied 
   over into the folders TridentG130M and TridentG160M. Each "sightline" has a 
   corresponding spectrum generated with the COS G130M and G160M instrument line-spread 
   functions. Move all synthetic spectra files to data/unanalyzed_spectra/.

2. Navigate to the scripts/analysis/ folder. Stitch the G130M and G160M spectra together 
   into one fits file folder by running "stitch_fits.py". By default, this will scan the
   contents of the data/unanalyzed_spectra/ folder and stitch together spectra that have 
   both "COS-G130M" and "COS-G160M" counterparts (following the naming convention above). 
   It will stitch the spectra into a single file starting with "COS-UV" (and ending with
   the orientation, theta, crmodel, galaxy age, and impact parameter). The format of this
   new file will have the structure expected by the subsequent analysis tools. The 
   individual G130M and G160M will be moved to the data/unanalyzed_spectra/trash/ 
   directory, which will be periodically cleared.
 
3. The next analysis steps (generating initial guesses of the voigt-profile parameters) 
   are wrapped into a single bash script. Edit the top section of"pyigm_guesses_script.sh"
   to indicate the model name, orientation, time, redshift, and the projected radius. If
   you want, you can keep the automation going longer by changing nsteps to be greater 
   than 1. I recommend doing this after you've used the bash script a few times. After 
   editing the relevant parameters, you can run the script in your terminal by running:
   "bash pyigm_guesses_script.sh". After the script is done running, it will create a new
   folder (with a name following the conventions described above) in data/analyzed_spectra/. 
   All relevant data to this synthetic spectrum will be moved into this folder, including 
   the original synthetic spectrum file from data/unanalyzed_spectra/. Note: for instructions 
   on how to use line tools and pyigm_igmguesses, see 'veepernotes.txt'.

4. Combine all data into single h5 file using "combine_all_spectra.py". Note that this will check
   the analyzed data folder and combine the information for all of the spectra in that folder.
   It will save the data to "combined_spectra.h5". Remember to move / rename this file to 
   the expected working directory.



#### Old. will delete soon
3. Go into the SpectralAnalysis folder. We'll be here a little while. 

4. Fix the header files of all of the .fits files by running $python fix_fits_files.py

5. Look for absorption features in spectra. The file "pyigm_guesses_script.sh" bundles together
   several different steps of this process (#TODO: elaborate). You will need to edit this file to change
   the orientation, model, time, and redshift. The script then loops through a bunch of radii starting with
   the variable "radius" and ending with "end_radius". I recommend starting this with just one iteration 
   before moving on to many. 
