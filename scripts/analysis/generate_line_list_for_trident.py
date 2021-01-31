from linetools.lists.linelist import LineList
import sys

def generate_line_list_file(ions = []):
    out = open('pyigm_line_list.txt', 'w')
    out.write('#Ion\tWavelength [A]\tgamma\t\tf_value\t\talt. name\n')

    llist = LineList('ISM')
    for ion in ions:
        ion_data = llist.all_transitions(ion.replace(" ", ""))
        print(ion)
        wls      = ion_data['wrest']
        gammas   = ion_data['gamma']
        fvals    = ion_data['f']
        altnames = ion_data['name']
        if wls.size == 1:
            sys.stdout.flush()
            wls      = [wls.value]
            gammas    = [gammas.value]
            fvals    = [fvals]
            altnames = [altnames]
        for wl, gamma, fval, altname in zip(wls, gammas, fvals, altnames):
            out.write("%s\t%0.6f\t%e\t%e\t%s\n"\
                          %(ion, wl, gamma, fval, altname))
            
    out.close()

ions = ['H I', 'O VI', 'C II', 'C III', 'C IV', \
                'Si II', 'Si III', 'Si IV', 'N V', 'Mg II']
generate_line_list_file(ions = ions)

