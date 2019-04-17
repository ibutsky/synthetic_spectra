import os
import sys
import numpy as np

def create_info_file(view, model, time, r, working_dir = '.'):

    if time == '11.2':
        redshift = '0.2'

    folder = '%s/COS-FUV_%s_%s_%sGyr_r%skpc'\
        %(working_dir, view, model, time, r)
    print(folder)
    if os.path.isdir(folder):
        outfile = open('%s/info.txt'%(folder), 'w')
        outfile.write("# view model time impact\n")
        outfile.write('%s %s %s %s %s\n'%(view, model, time, redshift, r))
        outfile.close()
               

#view = sys.argv[1]
#model = sys.argv[2]
#time = sys.argv[3]
#impact = sys.argv[4]
#working_dir = '../'

#create_info_file(view, model, time, impact, working_dir = working_dir)       

working_dir = '../../data/analyzed_spectra'
views = ['face', 'edge_theta0', 'edge_theta1.0', 'edge_theta1.5']
models = ['anisd', 'stream']
impacts = np.arange(10, 210, 10)
time = '11.2'

for view in views:
    for model in models:
        for impact in impacts:
            impact = str(impact)
            create_info_file(view, model, time, impact, working_dir = working_dir)

