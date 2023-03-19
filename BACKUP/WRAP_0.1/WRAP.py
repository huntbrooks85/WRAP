#-----------------------------------------------------------------------#
# WRAP v0.1
# By Hunter Brooks, at NAU, Flagstaff, 2023/2/~
#
# Purpose: Gathers photometry and astrometry from various catalogs given a RA, DEC, and Radius
#-----------------------------------------------------------------------#
from catalog_module.importmodule import *

#Sets the different outputs for the 3 different tabs
ML_KEY_SINGLE = '-ML-' + sg.WRITE_ONLY_KEY
ML_KEY_MULTI = '-ML2-' + sg.WRITE_ONLY_KEY
ML_KEY_MODEL = '-ML3-' + sg.WRITE_ONLY_KEY

# Prints when a catalog search has started
def start(catalog, tab):
    sg.cprint('------------------------------------------------                                                                                                              ', c='wheat4', end='', key = tab)
    sg.cprint('               Started ' + str(catalog) + ' Search                                                                                                            ', c='wheat4', end='', key = tab)
    sg.cprint('------------------------------------------------                                                                                                              ', c='wheat4', end='', key = tab)

#Prints when a catalog search has finished and found an object
def found(catalog, tab):
    sg.cprint('------------------------------------------------                                                                                                             ', c='wheat4', end='', key = tab)
    sg.cprint('             Finished ' + str(catalog) + ' Search                                                                                                            ', c='wheat4', end='', key = tab)
    sg.cprint('------------------------------------------------                                                                                                             ', c='wheat4', end='', key = tab)

#Prints when a catalog search has finished and did not find an object
def not_found(catalog, tab):
    sg.cprint('------------------------------------------------                                                                                                              ', c='wheat4', end='', key = tab)
    sg.cprint('                 Object Not Found                                                                                                                             ', c='wheat4', end='', key = tab)
    sg.cprint('                      Finished ' + str(catalog) + ' Search                                                                                                    ', c='wheat4', end='', key = tab)
    sg.cprint('------------------------------------------------                                                                                                              ', c='wheat4', end='', key = tab)

def wrap_end(tab):
    sg.cprint('------------------------------------------------                                                                                                              ', c='wheat4', end='', key = tab)
    sg.cprint('         All Catalogs Have Been Searched                                                                                                              ', c='wheat4', end='', key = tab)
    sg.cprint('             Finished Running WRAP                                                                                                                   ', c='wheat4', end='', key = tab)
    sg.cprint('------------------------------------------------                                                                                                              ', c='wheat4', end='', key = tab)

#Sets the theme for WRAP
sg.theme('LightGrey4')
#Makes the drop down window for types of file in the multi-object search
file_type_list = ['CSV', 'FITS', 'ASCII', 'IPAC']

#Makes the layout of the WRAP program, providing a place to put the ra, dec, and radius the output file name and which catalogs needs to be used
layout1 = [ [sg.Image(filename = ('image_module/temp_logo_1_50.png'), size = (135, 95)),          sg.Text('WRAP', justification='center', size=(5, 1), font = ('Chalkduster', 52)),     sg.Image(filename = 'image_module/BYW_Logo.png', size = (205, 95))],
        
            [sg.Text('RA', font = ('Times New Roman', 22), size=(13, 1), justification='center'), sg.Text('DEC', font = ('Times New Roman', 22), size=(13, 1), justification='center'), sg.Text('RADIUS', font = ('Times New Roman', 22), size=(13, 1), justification='center')],
            [sg.Text('(Degrees)', font = ('Times New Roman', 20), size=(18, 1), justification='center'), sg.Text('(Degrees)', font = ('Times New Roman', 20), size=(11, 1), justification='center'), sg.Text('(Arcsecs)', font = ('Times New Roman', 20), size=(20, 1), justification='center')],
            [sg.InputText(size=(18, 1), key = 'RA', font = ('Times New Roman', 15)),              sg.InputText(size=(18, 2), key = 'DEC', font = ('Times New Roman', 15)),              sg.InputText(size=(18, 2), key = 'RADIUS', font = ('Times New Roman', 15))],
            [sg.Text('Output File Name', size=(50, 1), justification='center', font = ('Times New Roman', 22))],
            [sg.InputText(key = 'output', font = ('Times New Roman', 15), size = (70, 3), justification='center')],

            [sg.Text('Catalogs:', justification='center', size=(50, 1), font = ('Times New Roman', 25))],   
            [sg.Checkbox('AllWISE', key = 'AW', font = ('Times New Roman', 22), size = (10, 2)),        sg.Checkbox('CatWISE 2020', key = 'catwise', font = ('Times New Roman', 22), size = (14, 2)),     sg.Checkbox('Gaia', key = 'gaia', font = ('Times New Roman', 22), size = (9, 2))],
            [sg.Checkbox('VISTA', key = 'VSA', font = ('Times New Roman', 22), size = (9, 2)),          sg.Checkbox('WFACM', key = 'UKIDSS', font = ('Times New Roman', 22), size = (10, 2)),             sg.Checkbox('2MASS', key = '2MASS', font = ('Times New Roman', 22), size = (10, 2))],
            [sg.Checkbox('PanSTARRS', key = 'ps', font = ('Times New Roman', 22), size = (13, 2)),      sg.Checkbox('NSC', key = 'nsc', font = ('Times New Roman', 22), size = (8, 2)),                   sg.Checkbox('GALEX', key = 'galex', font = ('Times New Roman', 22), size = (10, 2))],

            [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'), sg.Button('Help', size = (17), button_color = '#F7CC7C'), sg.Button('Close WRAP', size = (17), button_color = '#E48671')], 

            [sg.Text("\u0332".join('Output'), size=(50, 1), justification='center', font = ('Times New Roman', 15))],
            [sg.Multiline(size=(85, 6), write_only=(True), key=ML_KEY_SINGLE, reroute_stdout=True, echo_stdout_stderr=True, reroute_cprint=True)]]

#Makes the layout of the WRAP program, providing a place to put a csv list of ra and decs and radius the output file name and which catalogs needs to be used
layout2 = [ [sg.Image(filename = ('image_module/temp_logo_1_50.png'), size = (135, 95)),          sg.Text('WRAP', justification='center', size=(5, 1), font = ('Chalkduster', 52)), sg.Image(filename = 'image_module/BYW_Logo.png', size = (205, 95))],
        
            [sg.Text('FILE DIRECTORY', font = ('Times New Roman', 22), size=(50, 1), justification='center')],
            [sg.Text('(CSV, FITS, ASCII, IPAC)', font = ('Times New Roman', 20), size=(50, 1), justification='center')],
            [sg.FileBrowse('File Browser', size = (80, 1), key = 'file')],
            [sg.Text('RADIUS', font = ('Times New Roman', 22), size=(17, 1), justification='center'), sg.Text('FILETYPE', font = ('Times New Roman', 22), size=(9, 1), justification='center'),   sg.Text('Output File Name', size=(25, 1), justification='center', font = ('Times New Roman', 22))],
            [sg.InputText(size=(22, 2), key = 'RADIUS_multi', font = ('Times New Roman', 15)),              sg.Combo(file_type_list, size = (13), font = ('Times New Roman', 15), key = 'type'),   sg.InputText(key = 'output', font = ('Times New Roman', 15), size = (22, 2), justification='center')],

            [sg.Text('Catalogs:', justification='center', size=(50, 1), font = ('Times New Roman', 25))],
            [sg.Checkbox('AllWISE', key = 'AW_multi', font = ('Times New Roman', 22), size = (10, 2)),        sg.Checkbox('CatWISE 2020', key = 'catwise_multi', font = ('Times New Roman', 22), size = (14, 2)),     sg.Checkbox('Gaia', key = 'gaia_multi', font = ('Times New Roman', 22), size = (9, 2))],
            [sg.Checkbox('VISTA', key = 'VSA_multi', font = ('Times New Roman', 22), size = (9, 2)),          sg.Checkbox('WFACM', key = 'UKIDSS_multi', font = ('Times New Roman', 22), size = (10, 2)),             sg.Checkbox('2MASS', key = '2MASS_multi', font = ('Times New Roman', 22), size = (10, 2))],
            [sg.Checkbox('PanSTARRS', key = 'ps_multi', font = ('Times New Roman', 22), size = (13, 2)),      sg.Checkbox('NSC', key = 'nsc_multi', font = ('Times New Roman', 22), size = (8, 2)),                   sg.Checkbox('GALEX', key = 'galex_multi', font = ('Times New Roman', 22), size = (10, 2))],

            [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'), sg.Button('Help', size = (17), button_color = '#F7CC7C'), sg.Button('Close WRAP', size = (17), button_color = '#E48671')], 

            [sg.Text("\u0332".join('Output'), size=(50, 1), justification='center', font = ('Times New Roman', 15))],
            [sg.Multiline(size=(85, 6), write_only=(True), key=ML_KEY_MULTI, reroute_stdout=True, echo_stdout_stderr=True, reroute_cprint=True)]]

#Makes the layour of the WRAP program, that allows the user to fit their data to multiple different models
layout3 = [[sg.Image(filename = ('image_module/temp_logo_1_50.png'), size = (135, 95)),           sg.Text('WRAP', justification='center', size=(5, 1), font = ('Chalkduster', 52)), sg.Image(filename = 'image_module/BYW_Logo.png', size = (205, 95))], 

            [sg.Text('FILE DIRECTORY', font = ('Times New Roman', 22), size=(50, 1), justification='center')],
            [sg.Text('(CSV, FITS, ASCII, IPAC)', font = ('Times New Roman', 20), size=(50, 1), justification='center')],
            [sg.FileBrowse(size = (80, 1), key = 'file')],

            [sg.Text('Output File Name', size=(50, 1), justification='center', font = ('Times New Roman', 22))],
            [sg.InputText(key = 'output', font = ('Times New Roman', 15), size = (70, 3), justification='center')],
            
            [sg.Text('Models:', justification='center', size=(50, 1), font = ('Times New Roman', 25))], 
            [sg.Checkbox('\u03A7\u00b2-Color', key = 'chi_color', font = ('Times New Roman', 22), size = (16, 2)),               sg.Checkbox('\u03A7\u00b2-Spectra', key = 'chi_spectra', font = ('Times New Roman', 22), size = (13, 2)),           sg.Checkbox('C-C Plots', key = 'cc_plots', font = ('Times New Roman', 22))],

            [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'), sg.Button('Help', size = (17), button_color = '#F7CC7C'), sg.Button('Close WRAP', size = (17), button_color = '#E48671')],

            [sg.Text("\u0332".join('Output'), size=(50, 1), justification='center', font = ('Times New Roman', 15))],
            [sg.Multiline(size=(85, 6), write_only=(True), key=ML_KEY_MODEL, reroute_stdout=True, echo_stdout_stderr=True, reroute_cprint=True)]]

tab_layout = [[sg.TabGroup([[sg.Tab('Single Object', layout1, title_color='#F9F8F3', background_color='#FBF5EE', element_justification= 'center', key = 'Single Obect Search'),
                        sg.Tab('Multi-Object', layout2, title_color='#F9F8F3',background_color='#FBF5EE', element_justification= 'center', key = 'Multi-Object Search'),
                        sg.Tab('Model Fitting', layout3, title_color='#F9F8F3',background_color='#FBF5EE', element_justification= 'center', key = 'Model Fitting')]], 
                        tab_location='centertop', title_color='Black', tab_background_color='#F9F8F3',selected_title_color='Black', selected_background_color='#979793', border_width = 6, font = ('Times New Roman', 18), enable_events = True, key = 'tab_group'), sg.Button('Close')]] 

#Makes the window in the layout from above
window = sg.Window('WRAP', tab_layout, size = (550, 610), grab_anywhere=True, finalize=True)

while True:
    #Opens the window
    event, values = window.read()
    group = values['tab_group']

    if group == 'Single Obect Search':    
    #Then the 'Run WRAP' button is clicked it checks to see if the ra, dec, and radius are filled
        if event in (None, 'Run WRAP'):
            if values['RA'] == '':
                sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_SINGLE)
                sg.cprint('        Please enter a RA!                                                                                        ', c='wheat4', end='', key = ML_KEY_SINGLE)
                sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_SINGLE)
            if values['DEC'] == '':
                sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_SINGLE)
                sg.cprint('        Please enter a DEC!                                                                                       ', c='wheat4', end='', key = ML_KEY_SINGLE)
                sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_SINGLE)
            if values['RADIUS'] == '':
                sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_SINGLE)
                sg.cprint('        Please enter a radius!                                                                                    ', c='wheat4', end='', key = ML_KEY_SINGLE)
                sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_SINGLE)

            #If the ra, dec and radius are filled it then finds the catalogs wanted to be used and then runs the catalog searches with the ra, dec, and radius provided
            if values['RA'] != '' and values['DEC'] != '' and values['RADIUS'] != '':
                sg.cprint('------------------------------------------------                                                                                                                             ', c='wheat4', end='', key = ML_KEY_SINGLE)
                sg.cprint('Running Search on RA (deg): ' + str(values['RA']) + '                                                                                                                        ', c='wheat4', end='', key = ML_KEY_SINGLE)
                sg.cprint('DEC (deg): ' + str(values['DEC']) + '                                                                                                                                        ', c='wheat4', end='', key = ML_KEY_SINGLE)
                sg.cprint('RADIUS (arcsec): ' + str(values['RADIUS']) + '                                                                                                                               ', c='wheat4', end='', key = ML_KEY_SINGLE)
                sg.cprint('------------------------------------------------                                                                                                                             ', c='wheat4', end='', key = ML_KEY_SINGLE)

                #Makes fake information for all of the bands, in case the catalog relating to it was not called
                ra_use = float(values['RA'])
                dec_use = float(values['DEC'])
                radius_use = int(values['RADIUS'])

                #Imports all the catalog functions
                from catalog_module.aw_search import allwise_image
                from catalog_module.vsa_search import vsa_image
                from catalog_module.ukidss_search import ukidss_image
                from catalog_module.ps_search import ps_image
                from catalog_module.twomass_search import twomass_image
                from catalog_module.gaia_search import gaia_image
                from catalog_module.catwise_search import catwise_image
                from catalog_module.nsc_search import nsc_image
                from catalog_module.galex_search import galex_image

                #Makes lists for each of the functions names and variables
                catalogs = ['AW', 'catwise', 'gaia', 'VSA', 'UKIDSS', '2MASS', 'ps', 'nsc', 'galex']
                catalog_names = ['AllWISE', 'CatWISE 2020', 'GAIA', 'VSA', 'UKIDSS', '2MASS', 'PanSTARRS', 'NSC', 'GALEX']
                catalog_functions = [allwise_image, catwise_image, gaia_image, vsa_image, ukidss_image, twomass_image, ps_image, nsc_image, galex_image]
                value_names = [[ 'w1', 'w2', 'w3', 'w4'],
                               ['cw_w1', 'cw_w2'], 
                               ['parallax', 'g_mag', 'rad_v'],
                               ['vsa_j', 'vsa_h', 'vsa_ks'],
                               ['uks_j', 'uks_h', 'uks_k'],
                               ['j', 'h', 'ks'],
                               ['g', 'r', 'i', 'z', 'y'],
                               ['g', 'r', 'i', 'z'], 
                               ['nuv', 'fuv']]

                #Creates fake list of the data and data name
                photometry = []
                photometry_name = []
                #Starts searching through each catalog
                for q in range(len(catalogs)): 

                    #Calls the catalog and records the data down
                    if values[catalogs[q]] == True: 
                        start(catalog_names[q], ML_KEY_SINGLE)
                        search_catalog = catalog_functions[q](ra_use, dec_use, radius_use)
                        photometry.append(search_catalog)
                        photometry_name.append(value_names[q])

                        #Tells the user if the object was found or not
                        if isinstance(sum(filter(None, search_catalog)), float) == True:
                            found(catalog_names[q], ML_KEY_SINGLE)
                        else: 
                            not_found(catalog_names[q], ML_KEY_SINGLE)

                    #Puts in null data if catalog is not called
                    elif values[catalogs[q]] == False:
                        empty = np.empty((1, len(value_names[q]), 1))
                        empty[:] = np.nan
                        empty = empty.tolist()
                        photometry.append(empty)
                        photometry_name.append(value_names[q])

                    #Writes all of the data gathered into a csv file
                    if q == (len(catalogs) - 1):
                        wrap_end(ML_KEY_SINGLE)
                        output = values['output']
                        myFile = open('Output/' + str(output) + '.csv', 'w')
                        writer = csv.writer(myFile)
                        writer.writerow(photometry_name)
                        writer.writerow(photometry)

        #Provides the user with the authors information if the 'Help' button is pressed
        if event in (None, 'Help'):
            sg.cprint('------------------------------------------------                                                     ', end='', c='wheat4', key = ML_KEY_SINGLE)
            sg.cprint('     Thank you for using WRAP!                                                                       ', end='', c='wheat4', key = ML_KEY_SINGLE)
            sg.cprint('   Authors Contact: hcb98@nau.edu                                                                    ', end='', c='wheat4', key = ML_KEY_SINGLE)
            sg.cprint('------------------------------------------------                                                     ', end='', c='wheat4', key = ML_KEY_SINGLE)

        #Closes WRAP if the 'Close WRAP' button is pressed
        if event in (None, 'Close WRAP'):
            break

    if group == 'Multi-Object Search':
        if event in (None, 'Run WRAP1'):
            if values['file'] == '': 
                sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_MULTI)
                sg.cprint('        Please enter a file!                                                                                        ', c='wheat4', end='', key = ML_KEY_MULTI)
                sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_MULTI)
            if values['RADIUS_multi'] == '':
                sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_MULTI)
                sg.cprint('        Please enter a radius!                                                                                        ', c='wheat4', end='', key = ML_KEY_MULTI)
                sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_MULTI)
            if values['type'] == '': 
                sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_MULTI)
                sg.cprint('        Please enter a file type!                                                                                        ', c='wheat4', end='', key = ML_KEY_MULTI)
                sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_MULTI)

            if values['file'] != '' and values['RADIUS_multi'] != '' and values['type'] != '':
                sg.cprint('------------------------------------------------                                                                                                                                   ', c='wheat4', end='', key = ML_KEY_MULTI)
                sg.cprint('Directory: ' + str(values['file']) + '                                                                                                                                             ', c='wheat4', end='', key = ML_KEY_MULTI)
                sg.cprint('File Type: ' + str(values['type']) + '                                                                                                                                             ', c='wheat4', end='', key = ML_KEY_MULTI)
                sg.cprint('RADIUS (arcsec): ' + str(values['RADIUS_multi']) + '                                                                                                                               ', c='wheat4', end='', key = ML_KEY_MULTI)
                sg.cprint('------------------------------------------------                                                                                                                                   ', c='wheat4', end='', key = ML_KEY_MULTI)

                if values['type'] == 'CSV':
                    csv_table = pd.read_csv(values['file'])
                    ra_list = csv_table['ra'].tolist()
                    dec_list = csv_table['dec'].tolist()
                if values['type'] == 'FITS':
                    pass
                if values['type'] == 'ASCII':
                    pass
                if values['type'] == 'IPAC':
                    pass

                #Writes all of the data gathered into a csv file 
                output = values['output']
                header = ['ra', 'dec', 'W1', 'W2', 'W3', 'W4', 'J', 'H', 'Ks', 'g', 'r', 'i', 'z', 'y', 'plx', 'g_gaia', 'rad_v', 'cw_w1', 'cw_w2', 'uks_j', 'uks_h', 'uks_k', 'vsa_j', 'vsa_h', 'vsa_ks']
                myFile = open('Output/' + str(output) + '.csv', 'w')
                writer = csv.writer(myFile)
                writer.writerow(header)

                #Imports all the catalog functions
                from catalog_module.aw_search import allwise_image
                from catalog_module.vsa_search import vsa_image
                from catalog_module.ukidss_search import ukidss_image
                from catalog_module.ps_search import ps_image
                from catalog_module.twomass_search import twomass_image
                from catalog_module.gaia_search import gaia_image
                from catalog_module.catwise_search import catwise_image
                from catalog_module.nsc_search import nsc_image
                from catalog_module.galex_search import galex_image

                #Makes lists for each of the functions names and variables
                catalogs = ['AW_multi', 'catwise_multi', 'gaia_multi', 'VSA_multi', 'UKIDSS_multi', '2MASS_multi', 'ps_multi', 'nsc_multi', 'galex_multi']
                catalog_names = ['AllWISE', 'CatWISE 2020', 'GAIA', 'VSA', 'UKIDSS', '2MASS', 'PanSTARRS', 'NSC', 'GALEX']
                catalog_functions = [allwise_image, catwise_image, gaia_image, vsa_image, ukidss_image, twomass_image, ps_image, nsc_image, galex_image]
                value_names = [[ 'w1', 'w2', 'w3', 'w4'],
                               ['cw_w1', 'cw_w2'], 
                               ['parallax', 'g_mag', 'rad_v'],
                               ['vsa_j', 'vsa_h', 'vsa_ks'],
                               ['uks_j', 'uks_h', 'uks_k'],
                               ['j', 'h', 'ks'],
                               ['g', 'r', 'i', 'z', 'y'],
                               ['g', 'r', 'i', 'z'], 
                               ['nuv', 'fuv']]

                for index in range(len(ra_list)): 
                    #Makes fake information for all of the bands, in case the catalog relating to it was not called
                    ra_use = ra_list[index]
                    dec_use = dec_list[index]
                    radius_use = int(values['RADIUS_multi'])

                    #Creates fake list of the data and data name
                    photometry = []
                    photometry_name = []
                    #Starts searching through each catalog
                    for q in range(len(catalogs)): 
                        if values[catalogs[q]] == True: 
                            #Calls the catalog and records the data down if the function is called
                            start(catalog_names[q], ML_KEY_MULTI)
                            search_catalog = catalog_functions[q](ra_use, dec_use, radius_use)
                            photometry.append(search_catalog)
                            photometry_name.append(value_names[q])

                            #Tells the user if the object was found or not
                            if isinstance(sum(filter(None, search_catalog)), float) == True:
                                found(catalog_names[q], ML_KEY_SINGLE)
                            else: 
                                not_found(catalog_names[q], ML_KEY_SINGLE)

                        elif values[catalogs[q]] == False:
                            #Puts in fake data if catalog is not called
                            empty = np.empty((1, len(value_names[q]), 1))
                            empty[:] = np.nan
                            empty = empty.tolist()
                            photometry.append(empty)
                            photometry_name.append(value_names[q])

                        if q == (len(catalogs) - 1):
                            wrap_end(ML_KEY_MULTI)
                            # Writes all of the data gathered into a csv file 
                            writer.writerow(photometry)

        #Provides the user with the authors information if the 'Help' button is pressed
        if event in (None, 'Help2'):
            sg.cprint('------------------------------------------------                                                     ', end='', c='wheat4', key = ML_KEY_MULTI)
            sg.cprint('     Thank you for using WRAP!                                                                       ', end='', c='wheat4', key = ML_KEY_MULTI)
            sg.cprint('   Authors Contact: hcb98@nau.edu                                                                    ', end='', c='wheat4', key = ML_KEY_MULTI)
            sg.cprint('------------------------------------------------                                                     ', end='', c='wheat4', key = ML_KEY_MULTI)

        #Closes WRAP if the 'Close WRAP' button is pressed
        if event in (None, 'Close WRAP3'):
            break

window.close()
