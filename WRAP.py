#-----------------------------------------------------------------------#
# WRAP v0.4
# By Hunter Brooks, at NAU, Flagstaff: May 5, 2023
#
# Purpose: Gathers photometry and astrometry from various 
#          ultra-violet, optical, and near-infrared catalogs 
#          given a RA, DEC, and Radius by the user
#-----------------------------------------------------------------------#

#Imports all needed packages
from catalog_module.importmodule import *

#Imports all the catalog scripts
from catalog_module.aw_search import allwise_image
from catalog_module.vsa_search import vsa_image
from catalog_module.ukidss_search import ukidss_image
from catalog_module.ps_search import ps_image
from catalog_module.twomass_search import twomass_image
from catalog_module.gaia_search import gaia_image
from catalog_module.catwise_search import catwise_image
from catalog_module.nsc_search import nsc_image
from catalog_module.galex_search import galex_image

def single_object_search(): 
  '''Performs the single object search by going to each catalog script and recording the data down.
  Once all the data is recorded it it put into a CSV file for the user.'''
  for q in range(len(catalogs)): 

    #Calls the catalog and records the data down
    if values[catalogs[q]] == True: 
      wrap_start(catalog_names[q], ML_KEY_SINGLE)
      search_catalog = catalog_functions[q](ra_use, dec_use, radius_use)
      photometry.append(search_catalog)
      photometry_name.append(value_names[q])

      #Tells the user if the object was found or not
      if search_catalog is None: 
        wrap_end(ML_KEY_SINGLE)
      else: 
        if isinstance(sum(filter(None, search_catalog[0:len(search_catalog) - 2])), float) == True:
          wrap_found(catalog_names[q], ML_KEY_SINGLE)
        else: 
          wrap_not_found(catalog_names[q], ML_KEY_SINGLE)

    #Puts in null data if catalog is not called
    elif values[catalogs[q]] == False:
      empty = [np.nan] * len(value_names[q])
      photometry.append(empty)
      photometry_name.append(value_names[q])

    #Writes all of the data gathered into a csv file 
    if q == (len(catalogs) - 1):
      wrap_end(ML_KEY_SINGLE)

      #Makes the output file name
      if values['output'] == '':
        output = 'WRAP_output'
      else: 
        output = values['output']

      #Writes the CSV file with the photometry and astrometry gathered
      myFile = open('Output/' + str(output) + '.csv', 'w')
      writer = csv.writer(myFile)
      flat_photometry_list = [item for sublist in photometry for item in sublist]
      flat_photometry_name_list = [item for sublist in photometry_name for item in sublist]
      writer.writerow(flat_photometry_name_list)
      writer.writerow(flat_photometry_list)
      myFile.close()

def single_tab_check(): 
  '''Checks if the User put in the correct formats for the RA, DEC, and Radius options.'''

  #Checks if the RA tab is entered with a number
  fake_list = []
  try: 
    float(values['RA'])
  except ValueError:
    sg.cprint('------------------------------------------------                                                                  ', c='wheat4', key = ML_KEY_SINGLE)
    sg.cprint('Please enter a Correct RA!                                                                                        ', c='wheat4', end='', key = ML_KEY_SINGLE)
    sg.cprint('------------------------------------------------                                                                  ', c='wheat4', key = ML_KEY_SINGLE)
    fake_list.append(1)

  #Checks if the DEC tab is entered with a number
  try: 
    float(values['DEC'])
  except ValueError:
    sg.cprint('------------------------------------------------                                                                  ', c='wheat4', key = ML_KEY_SINGLE)
    sg.cprint('Please enter a Correct DEC!                                                                                       ', c='wheat4', end='', key = ML_KEY_SINGLE)
    sg.cprint('------------------------------------------------                                                                  ', c='wheat4', key = ML_KEY_SINGLE)
    fake_list.append(2)

  #Checks if the RADIUS tab is entered with a number
  if values['RADIUS'].isnumeric() == False:
    sg.cprint('------------------------------------------------                                                                  ', c='wheat4', key = ML_KEY_SINGLE)
    sg.cprint('Please enter a Correct radius!                                                                                    ', c='wheat4', end='', key = ML_KEY_SINGLE)
    sg.cprint('------------------------------------------------                                                                  ', c='wheat4', key = ML_KEY_SINGLE)
    fake_list.append(3)
  return fake_list

def multi_object_search():
  '''Performs the Multi-object search by going to each catalog script and recording the data down.
  Once all the data is recorded it it put into a CSV file for the user.'''
  for index in range(len(ra_list)): 

    #Make variables for the RA, DEC, and RADIUS used
    ra_use = ra_list[index]
    dec_use = dec_list[index]
    radius_use = int(values['RADIUS_multi'])

    #Creates fake list for the data and data names
    photometry = []
    photometry_name = []

    #Starts searching through each catalog
    for q in range(len(catalogs)): 

      #Calls the catalog and records the data down
      if values[catalogs_multi[q]] == True: 
        wrap_start(catalog_names[q], ML_KEY_MULTI)
        search_catalog = catalog_functions[q](ra_use, dec_use, radius_use)
        photometry.append(search_catalog)  

        #Tells the user if the object was found or not
        if isinstance(sum(filter(None, search_catalog[0:len(search_catalog) - 2])), float) == True:
          wrap_found(catalog_names[q], ML_KEY_SINGLE)
        else: 
          wrap_not_found(catalog_names[q], ML_KEY_SINGLE)  

      #Puts in null data if catalog is not called
      elif values[catalogs_multi[q]] == False:
        empty = [np.nan] * len(value_names[q])
        photometry.append(empty)
        photometry_name.append(value_names[q])   

      #Writes all of the data gathered into a csv file
      if q == (len(catalogs) - 1):
        wrap_end(ML_KEY_MULTI)   

        #Writes a new list for the objects photometry and astrometry
        flat_photometry_list = [item for sublist in photometry for item in sublist]
        writer.writerow(flat_photometry_list)

    #Finishes the file once every object is searched
    if index == len(ra_list) - 1:
      myFile.close()

def multi_tab_check():
  '''Checks if the User put in anything for the File, FileType, and Radius options.'''
  if values['file'] == '': 
    sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_MULTI)
    sg.cprint('Please enter a correct file directory!                                                               ', c='wheat4', end='', key = ML_KEY_MULTI)
    sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_MULTI)

  #Checks if the RADIUS tab is entered
  if values['RADIUS_multi'].isnumeric() == False:
    sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_MULTI)
    sg.cprint('Please enter a correct radius value!                                                                 ', c='wheat4', end='', key = ML_KEY_MULTI)
    sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_MULTI)

  #Checks if the FILETYPE tab is entered 
  if values['type'] == '': 
    sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_MULTI)
    sg.cprint('Please enter a correct file type!                                                                    ', c='wheat4', end='', key = ML_KEY_MULTI)
    sg.cprint('------------------------------------------------                                                     ', c='wheat4', key = ML_KEY_MULTI)

def multi_tab_table(): 
  '''Reads in the file the user chose for the multi-object search'''

  #Reads in the file depending on the filetype
  if values['type'] == 'CSV':
    csv_table = pd.read_csv(values['file'])
    ra_list = csv_table['ra'].tolist()
    dec_list = csv_table['dec'].tolist()
  if values['type'] == 'FITS':
    fits_table = fits.open(values['file'])
    fits_data = fits_table[1].data
    ra_list = fits_data['ra'].tolist()
    dec_list = fits_data['dec'].tolist()
  if values['type'] == 'ASCII':
    ascii_table = ascii.read(values['file'])
    ra_list = ascii_table['ra'].tolist()
    dec_list = ascii_table['dec'].tolist()
  if values['type'] == 'IPAC':
    ipac_table = ascii.read(values['file'], format = 'ipac')
    ra_list = ipac_table['ra'].tolist()
    dec_list = ipac_table['dec'].tolist()
  return ra_list, dec_list

#Sets the different outputs for the 2 different tabs
ML_KEY_SINGLE = '-ML-'  + sg.WRITE_ONLY_KEY
ML_KEY_MULTI  = '-ML2-' + sg.WRITE_ONLY_KEY

#Prints when the catalog search is started
def wrap_start(catalog, tab):
  '''The print function for when a catalog search has started'''

  sg.cprint('------------------------------------------------                                                                                                               ', c='wheat4', end='', key = tab)
  sg.cprint('Started ' + str(catalog) + ' Search                                                                                                                            ', c='wheat4', end='', key = tab)
  sg.cprint('------------------------------------------------                                                                                                               ', c='wheat4', end='', key = tab)

#Prints when a catalog search has finished and found an object
def wrap_found(catalog, tab):
  '''Print function for when the catalog search was successful.'''

  sg.cprint('------------------------------------------------                                                                                                               ', c='wheat4', end='', key = tab)
  sg.cprint('Finished ' + str(catalog) + ' Search                                                                                                                           ', c='wheat4', end='', key = tab)
  sg.cprint('------------------------------------------------                                                                                                               ', c='wheat4', end='', key = tab)

#Prints when a catalog search has finished and did not find an object
def wrap_not_found(catalog, tab):
  '''Print function for when the catalog search was not successful.'''

  sg.cprint('------------------------------------------------                                                                                                               ', c='wheat4', end='', key = tab)
  sg.cprint('Object Not Found                                                                                                                                               ', c='wheat4', end='', key = tab)
  sg.cprint('Finished ' + str(catalog) + ' Search                                                                                                                           ', c='wheat4', end='', key = tab)
  sg.cprint('------------------------------------------------                                                                                                               ', c='wheat4', end='', key = tab)

#Prints when the catalog search has finished
def wrap_end(tab):
  '''The print function for when all catalogs have been searched'''

  sg.cprint('------------------------------------------------                                                                                                               ', c='wheat4', end='', key = tab)
  sg.cprint('All Catalogs Have Been Searched                                                                                                                                ', c='wheat4', end='', key = tab)
  sg.cprint('Finished Running WRAP                                                                                                                                          ', c='wheat4', end='', key = tab)
  sg.cprint('------------------------------------------------                                                                                                               ', c='wheat4', end='', key = tab)

#Sets the theme for WRAP
sg.theme('LightGrey4')

#Sets the layout if the user is not on a windows machine
if platform != 'win32':
  #Makes the layout of WRAP for the single object search, by providing a location for: ra, dec, radius, output file name, catalogs, and output
  layout_single = [ [sg.Image(filename = ('image_module/WRAP_Logo.png'), size = (135, 95)),                    sg.Text('WRAP', justification='center', size=(5, 1), font = ('Chalkduster', 52)),                 sg.Image(filename = 'image_module/BYW_Logo.png', size = (205, 95))],
          
                    [sg.Text('RA', font = ('Times New Roman', 22), size=(13, 1), justification='center'),           sg.Text('DEC', font = ('Times New Roman', 22), size=(13, 1), justification='center'),             sg.Text('RADIUS', font = ('Times New Roman', 22), size=(13, 1), justification='center')],
                    [sg.Text('(Degrees)', font = ('Times New Roman', 20), size=(18, 1), justification='center'),    sg.Text('(Degrees)', font = ('Times New Roman', 20), size=(11, 1), justification='center'),       sg.Text('(Arcsecs)', font = ('Times New Roman', 20), size=(20, 1), justification='center')],
                    [sg.InputText(size=(18, 1), key = 'RA', font = ('Times New Roman', 15)),                        sg.InputText(size=(18, 2), key = 'DEC', font = ('Times New Roman', 15)),                          sg.InputText(size=(18, 2), key = 'RADIUS', font = ('Times New Roman', 15))],
                    
                    [sg.Text('Output File Name', size=(50, 1), justification='center', font = ('Times New Roman', 22))],
                    [sg.InputText(key = 'output', font = ('Times New Roman', 15), size = (70, 3), justification='center')],

                    [sg.Text('Catalogs:', justification='center', size=(50, 1), font = ('Times New Roman', 25))],   
                    [sg.Checkbox('CatWISE 2020', key = 'catwise', font = ('Times New Roman', 22), size = (14, 2)),   sg.Checkbox('AllWISE', key = 'AW', font = ('Times New Roman', 22), size = (10, 2)),               sg.Checkbox('Gaia', key = 'gaia', font = ('Times New Roman', 22), size = (9, 2))],
                    [sg.Checkbox('VISTA', key = 'VSA', font = ('Times New Roman', 22), size = (9, 2)),               sg.Checkbox('WFCAM', key = 'UKIDSS', font = ('Times New Roman', 22), size = (10, 2)),             sg.Checkbox('2MASS', key = '2MASS', font = ('Times New Roman', 22), size = (10, 2))],
                    [sg.Checkbox('PanSTARRS', key = 'ps', font = ('Times New Roman', 22), size = (13, 2)),           sg.Checkbox('NSC', key = 'nsc', font = ('Times New Roman', 22), size = (8, 2)),                   sg.Checkbox('GALEX', key = 'galex', font = ('Times New Roman', 22), size = (10, 2))],

                    [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'),                                   sg.Button('Help', size = (17), button_color = '#F7CC7C'),                                         sg.Button('Close WRAP', size = (17), button_color = '#E48671')], 

                    [sg.Text("\u0332".join('Output'), size=(50, 1), justification='center', font = ('Times New Roman', 15))],
                    [sg.Multiline(size=(85, 6), write_only=(True), key=ML_KEY_SINGLE, reroute_stdout=True, echo_stdout_stderr=True, reroute_cprint=True)]]

  #Makes the drop down window for types of file in the multi-object search
  filetype_list = ['CSV', 'FITS', 'ASCII', 'IPAC']
  #Makes the layout of WRAP for the multi-object search, by providing a location for: file directory, radius, filetype, output file name, catalogs, and output
  layout_multi = [ [sg.Image(filename = ('image_module/WRAP_Logo.png'), size = (135, 95)),                      sg.Text('WRAP', justification='center', size=(5, 1), font = ('Chalkduster', 52)),                       sg.Image(filename = 'image_module/BYW_Logo.png', size = (205, 95))],
          
                  [sg.Text('FILE DIRECTORY', font = ('Times New Roman', 22), size=(50, 1), justification='center')],
                  [sg.Text('(CSV, FITS, ASCII, IPAC)', font = ('Times New Roman', 20), size=(50, 1), justification='center')],

                  [sg.FileBrowse('File Browser', size = (80, 1), key = 'file', file_types = [('CSV Files', '*.csv'), ('FITS Files', '*.fits'), ('ASCII Files', '*.txt'), ('IPAC Files', '*.txt')])],
                  [sg.Text('RADIUS', font = ('Times New Roman', 22), size=(17, 1), justification='center'),              sg.Text('FILETYPE', font = ('Times New Roman', 22), size=(9, 1), justification='center'),               sg.Text('Output File Name', size=(25, 1), justification='center', font = ('Times New Roman', 22))],
                  [sg.InputText(size=(22, 2), key = 'RADIUS_multi', font = ('Times New Roman', 15)),                     sg.Combo(filetype_list, size = (13), font = ('Times New Roman', 15), key = 'type'),                     sg.InputText(key = 'output2', font = ('Times New Roman', 15), size = (22, 2), justification='center')],

                  [sg.Text('Catalogs:', justification='center', size=(50, 1), font = ('Times New Roman', 25))],
                  [sg.Checkbox('CatWISE 2020', key = 'catwise_multi', font = ('Times New Roman', 22), size = (14, 2)),   sg.Checkbox('AllWISE', key = 'AW_multi', font = ('Times New Roman', 22), size = (10, 2)),               sg.Checkbox('Gaia', key = 'gaia_multi', font = ('Times New Roman', 22), size = (9, 2))],
                  [sg.Checkbox('VISTA', key = 'VSA_multi', font = ('Times New Roman', 22), size = (9, 2)),               sg.Checkbox('WFCAM', key = 'UKIDSS_multi', font = ('Times New Roman', 22), size = (10, 2)),             sg.Checkbox('2MASS', key = '2MASS_multi', font = ('Times New Roman', 22), size = (10, 2))],
                  [sg.Checkbox('PanSTARRS', key = 'ps_multi', font = ('Times New Roman', 22), size = (13, 2)),           sg.Checkbox('NSC', key = 'nsc_multi', font = ('Times New Roman', 22), size = (8, 2)),                   sg.Checkbox('GALEX', key = 'galex_multi', font = ('Times New Roman', 22), size = (10, 2))],

                  [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'),                                         sg.Button('Help', size = (17), button_color = '#F7CC7C'),                                               sg.Button('Close WRAP', size = (17), button_color = '#E48671')], 

                  [sg.Text("\u0332".join('Output'), size=(50, 1), justification='center', font = ('Times New Roman', 15))],
                  [sg.Multiline(size=(85, 6), write_only=(True), key=ML_KEY_MULTI, reroute_stdout=True, echo_stdout_stderr=True, reroute_cprint=True)]]

  #Makes the general layout for WRAP
  tab_layout = [[sg.TabGroup([[sg.Tab('Single Object',   layout_single,       title_color='#F9F8F3',          background_color='#FBF5EE',   element_justification= 'center',     key = 'Single Obect Search'),
                              sg.Tab('Multi-Object',    layout_multi,        title_color='#F9F8F3',          background_color='#FBF5EE',   element_justification= 'center',     key = 'Multi-Object Search')]], 
                              tab_location='centertop', title_color='Black', tab_background_color='#F9F8F3', selected_title_color='Black', selected_background_color='#979793', border_width = 6, font = ('Times New Roman', 18), enable_events = True, key = 'tab_group'), sg.Button('Close')]] 

  #Generates the window based off the layouts above
  window = sg.Window('WRAP', tab_layout, size = (550, 610), grab_anywhere=False, finalize=True, enable_close_attempted_event = True)
  
  #Makes lists for: catalog names, catalog functions, and fake variables
  catalogs = ['catwise', 'AW', 'gaia', 
              'VSA', 'UKIDSS', '2MASS', 
              'ps', 'nsc', 'galex']
  catalogs_multi = ['catwise_multi', 'AW_multi', 'gaia_multi', 
                    'VSA_multi', 'UKIDSS_multi', '2MASS_multi', 
                    'ps_multi', 'nsc_multi', 'galex_multi']
  catalog_names = ['CatWISE 2020', 'AllWISE', 'GAIA', 
                  'VSA', 'UKIDSS', '2MASS', 
                  'PanSTARRS', 'NSC', 'GALEX']
  catalog_functions = [catwise_image, allwise_image, gaia_image, 
                      vsa_image, ukidss_image, twomass_image, 
                      ps_image, nsc_image, galex_image]

  #Defines all of the variable names
  value_names = [['cw_ra', 'cw_ra_e', 'cw_dec', 'cw_dec_e', 'cw_w1', 'cw_w1_e', 'cw_w2', 'cw_w2_e', 'cw_pmra', 'cw_pmra_e', 'cw_pmdec', 'cw_pmdec_e', 'cw_mjd', 'cw_catalog', 'cw_notes'], 
                ['aw_ra', 'aw_ra_e', 'aw_dec', 'aw_dec_e', 'aw_w1', 'aw_w1_e', 'aw_w2', 'aw_w2_e', 'aw_w3', 'aw_w3_e', 'aw_w4', 'aw_w4_e', 'aw_pmra', 'aw_pmra_e', 'aw_pmdec', 'aw_pmdec_e', 'aw_catalog', 'aw_notes'], 
                ['gaia_ra', 'gaia_ra_e', 'gaia_dec', 'gaia_dec_e', 'gaia_parallax', 'gaia_parallax_e', 'gaia_radv', 'gaia_radv_e', 'gaia_pmra', 'gaia_pmra_e', 'gaia_pmdec', 'gaia_pmdec_e', 'gaia_g', 'gaia_g_e', 'gaia_bp', 'gaia_bp_e', 'gaia_rp', 'gaia_rp_e', 'gaia_year', 'gaia_catalog', 'gaia_notes'], 
                ['vsa_ra', 'vsa_dec', 'vsa_y', 'vsa_y_e', 'vsa_j', 'vsa_j_e', 'vsa_h', 'vsa_h_e', 'vsa_ks', 'vsa_ks_e', 'vsa_mjd_y', 'vsa_mjd_j', 'vsa_mjd_h', 'vsa_mjd_ks', 'vsa_catalog', 'vsa_notes'], 
                ['wfcam_ra', 'wfcam_ra_e', 'wfcam_dec', 'wfcam_dec_e', 'wfcam_y', 'wfcam_y_e', 'wfcam_j', 'wfcam_j_e', 'wfcam_h', 'wfcam_h_e', 'wfcam_k', 'wfcam_k_e', 'wfcam_pmra', 'wfcam_pmra_e', 'wfcam_pmdec', 'wfcam_pmdec_e', 'wfcam_epoch', 'wfcam_catalog', 'wfcam_notes'], 
                ['2mass_ra', '2mass_dec', '2mass_j', '2mass_j_e', '2mass_h', '2mass_h_e', '2mass_ks', '2mass_ks_e', '2mass_catalog', '2mass_notes'], 
                ['ps_ra', 'ps_ra_e', 'ps_dec', 'ps_dec_e', 'ps_g', 'ps_g_e', 'ps_r', 'ps_r_e', 'ps_i', 'ps_i_e', 'ps_z', 'ps_z_e', 'ps_y', 'ps_y_e', 'ps_mjd', 'ps_catalog', 'ps_notes'], 
                ['nsc_ra', 'nsc_ra_e', 'nsc_dec', 'nsc_dec_e', 'nsc_g', 'nsc_g_e', 'nsc_r', 'nsc_r_e', 'nsc_i', 'nsc_i_e', 'nsc_z', 'nsc_z_e', 'nsc_u', 'nsc_u_e', 'nsc_y', 'nsc_y_e', 'nsc_pmra', 'nsc_pmra_e', 'nsc_pmdec', 'nsc_pmdec_e', 'nsc_mjd', 'nsc_catalog', 'nsc_notes'],  
                ['galex_ra', 'galex_dec', 'galex_fuv', 'galex_fuv_e', 'galex_nuv', 'galex_nuv_e', 'galex_catalog', 'galex_notes']]

  #Defines the header for the CSV file
  header = ['cw_ra', 'cw_ra_e', 'cw_dec', 'cw_dec_e', 'cw_w1', 'cw_w1_e', 'cw_w2', 'cw_w2_e', 'cw_pmra', 'cw_pmra_e', 'cw_pmdec', 'cw_pmdec_e', 'cw_mjd', 'cw_catalog', 'cw_notes', 
            'aw_ra', 'aw_ra_e', 'aw_dec', 'aw_dec_e', 'aw_w1', 'aw_w1_e', 'aw_w2', 'aw_w2_e', 'aw_w3', 'aw_w3_e', 'aw_w4', 'aw_w4_e', 'aw_pmra', 'aw_pmra_e', 'aw_pmdec', 'aw_pmdec_e', 'aw_catalog', 'aw_notes', 
            'gaia_ra', 'gaia_ra_e', 'gaia_dec', 'gaia_dec_e', 'gaia_parallax', 'gaia_parallax_e', 'gaia_radv', 'gaia_radv_e', 'gaia_pmra', 'gaia_pmra_e', 'gaia_pmdec', 'gaia_pmdec_e', 'gaia_g', 'gaia_g_e', 'gaia_bp', 'gaia_bp_e', 'gaia_rp', 'gaia_rp_e', 'gaia_year', 'gaia_catalog', 'gaia_notes', 
            'vsa_ra', 'vsa_dec', 'vsa_y', 'vsa_y_e', 'vsa_j', 'vsa_j_e', 'vsa_h', 'vsa_h_e', 'vsa_ks', 'vsa_ks_e', 'vsa_mjd_y', 'vsa_mjd_j', 'vsa_mjd_h', 'vsa_mjd_ks', 'vsa_catalog', 'vsa_notes', 
            'wfcam_ra', 'wfcam_ra_e', 'wfcam_dec', 'wfcam_dec_e', 'wfcam_y', 'wfcam_y_e', 'wfcam_j', 'wfcam_j_e', 'wfcam_h', 'wfcam_h_e', 'wfcam_k', 'wfcam_k_e', 'wfcam_pmra', 'wfcam_pmra_e', 'wfcam_pmdec', 'wfcam_pmdec_e', 'wfcam_epoch', 'wfcam_catalog', 'wfcam_notes', 
            '2mass_ra', '2mass_dec', '2mass_j', '2mass_j_e', '2mass_h', '2mass_h_e', '2mass_ks', '2mass_ks_e', '2mass_catalog', '2mass_notes', 
            'ps_ra', 'ps_ra_e', 'ps_dec', 'ps_dec_e', 'ps_g', 'ps_g_e', 'ps_r', 'ps_r_e', 'ps_i', 'ps_i_e', 'ps_z', 'ps_z_e', 'ps_y', 'ps_y_e', 'ps_mjd', 'ps_catalog', 'ps_notes', 
            'nsc_ra', 'nsc_ra_e', 'nsc_dec', 'nsc_dec_e', 'nsc_g', 'nsc_g_e', 'nsc_r', 'nsc_r_e', 'nsc_i', 'nsc_i_e', 'nsc_z', 'nsc_z_e', 'nsc_u', 'nsc_u_e', 'nsc_y', 'nsc_y_e', 'nsc_pmra', 'nsc_pmra_e', 'nsc_pmdec', 'nsc_pmdec_e', 'nsc_mjd', 'nsc_catalog', 'nsc_notes',  
            'galex_ra', 'galex_dec', 'galex_fuv', 'galex_fuv_e', 'galex_nuv', 'galex_nuv_e', 'galex_catalog', 'galex_notes']

#Sets the layout if the user is on a windows machine
elif platform == 'win32':
  #Makes the layout of WRAP for the single object search, by providing a location for: ra, dec, radius, output file name, catalogs, and output
  layout_single = [ [sg.Image(filename = ('image_module/WRAP_Logo.png'), size = (135, 95)),                    sg.Text('WRAP', justification='center', size=(6, 1), font = ('Ink Free', 45)),                 sg.Image(filename = 'image_module/BYW_Logo.png', size = (205, 95))],
          
                    [sg.Text('RA', font = ('Times New Roman', 18), size=(13, 1), justification='center'),           sg.Text('DEC', font = ('Times New Roman', 18), size=(10, 1), justification='center'),             sg.Text('RADIUS', font = ('Times New Roman', 18), size=(13, 1), justification='center')],
                    [sg.Text('(Degrees)', font = ('Times New Roman', 16), size=(14, 1), justification='center'),    sg.Text('(Degrees)', font = ('Times New Roman', 16), size=(11, 1), justification='center'),       sg.Text('(Arcsecs)', font = ('Times New Roman', 16), size=(16, 1), justification='center')],
                    [sg.InputText(size=(15), key = 'RA', font = ('Times New Roman', 15)),                        sg.InputText(size=(15, 2), key = 'DEC', font = ('Times New Roman', 15)),                          sg.InputText(size=(15, 2), key = 'RADIUS', font = ('Times New Roman', 15))],
                    
                    [sg.Text('Output File Name', size=(50, 1), justification='center', font = ('Times New Roman', 18))],
                    [sg.InputText(key = 'output', font = ('Times New Roman', 15), size = (70, 3), justification='center')],

                    [sg.Text('Catalogs:', justification='center', size=(50, 1), font = ('Times New Roman', 22))],   
                    [sg.Checkbox('CatWISE 2020', key = 'catwise', font = ('Times New Roman', 15), size = (14, 2)),   sg.Checkbox('AllWISE', key = 'AW', font = ('Times New Roman', 15), size = (10, 2)),               sg.Checkbox('Gaia', key = 'gaia', font = ('Times New Roman', 15), size = (9, 2))],
                    [sg.Checkbox('VISTA', key = 'VSA', font = ('Times New Roman', 15), size = (9, 2)),               sg.Checkbox('WFCAM', key = 'UKIDSS', font = ('Times New Roman', 15), size = (10, 2)),             sg.Checkbox('2MASS', key = '2MASS', font = ('Times New Roman', 15), size = (10, 2))],
                    [sg.Checkbox('PanSTARRS', key = 'ps', font = ('Times New Roman', 15), size = (13, 2)),           sg.Checkbox('GALEX', key = 'galex', font = ('Times New Roman', 15), size = (10, 2))],

                    [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'),                                   sg.Button('Help', size = (17), button_color = '#F7CC7C'),                                         sg.Button('Close WRAP', size = (17), button_color = '#E48671')], 

                    [sg.Text("\u0332".join('Output'), size=(50, 1), justification='center', font = ('Times New Roman', 15))],
                    [sg.Multiline(size=(85, 6), write_only=(True), key=ML_KEY_SINGLE, reroute_stdout=True, echo_stdout_stderr=True, reroute_cprint=True)]]

  #Makes the drop down window for types of file in the multi-object search
  filetype_list = ['CSV', 'FITS', 'ASCII', 'IPAC']
  #Makes the layout of WRAP for the multi-object search, by providing a location for: file directory, radius, filetype, output file name, catalogs, and output
  layout_multi = [ [sg.Image(filename = ('image_module/WRAP_Logo.png'), size = (135, 95)),                      sg.Text('WRAP', justification='center', size=(6, 1), font = ('Ink Free', 45)),                       sg.Image(filename = 'image_module/BYW_Logo.png', size = (205, 95))],
          
                  [sg.Text('FILE DIRECTORY', font = ('Times New Roman', 18), size=(50, 1), justification='center')],
                  [sg.Text('(CSV, FITS, ASCII, IPAC)', font = ('Times New Roman', 16), size=(50, 1), justification='center')],
                  [sg.FileBrowse('File Browser', size = (80, 1), key = 'file', file_types = [('CSV Files', '*.csv'), ('FITS Files', '*.fits'), ('ASCII Files', '*.txt'), ('IPAC Files', '*.txt')])],

                  [sg.Text('RADIUS', font = ('Times New Roman', 16), size=(14, 1), justification='center'),              sg.Text('FILETYPE', font = ('Times New Roman', 16), size=(10, 1), justification='center'),               sg.Text('Output File Name', size=(26, 1), justification='center', font = ('Times New Roman', 16))],
                  [sg.InputText(size=(16, 2), key = 'RADIUS_multi', font = ('Times New Roman', 15)),                     sg.Combo(filetype_list, size = (13), font = ('Times New Roman', 15), key = 'type'),                     sg.InputText(key = 'output2', font = ('Times New Roman', 15), size = (22, 2), justification='center')],

                  [sg.Text('Catalogs:', justification='center', size=(50, 1), font = ('Times New Roman', 25))],
                  [sg.Checkbox('CatWISE 2020', key = 'catwise_multi', font = ('Times New Roman', 15), size = (14, 2)),   sg.Checkbox('AllWISE', key = 'AW_multi', font = ('Times New Roman', 15), size = (10, 2)),               sg.Checkbox('Gaia', key = 'gaia_multi', font = ('Times New Roman', 15), size = (9, 2))],
                  [sg.Checkbox('VISTA', key = 'VSA_multi', font = ('Times New Roman', 15), size = (9, 2)),               sg.Checkbox('WFCAM', key = 'UKIDSS_multi', font = ('Times New Roman', 15), size = (10, 2)),             sg.Checkbox('2MASS', key = '2MASS_multi', font = ('Times New Roman', 15), size = (10, 2))],
                  [sg.Checkbox('PanSTARRS', key = 'ps_multi', font = ('Times New Roman', 15), size = (13, 2)),           sg.Checkbox('GALEX', key = 'galex_multi', font = ('Times New Roman', 15), size = (10, 2))],

                  [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'),                                         sg.Button('Help', size = (17), button_color = '#F7CC7C'),                                               sg.Button('Close WRAP', size = (17), button_color = '#E48671')], 

                  [sg.Text("\u0332".join('Output'), size=(50, 1), justification='center', font = ('Times New Roman', 15))],
                  [sg.Multiline(size=(85, 6), write_only=(True), key=ML_KEY_MULTI, reroute_stdout=True, echo_stdout_stderr=True, reroute_cprint=True)]]

  #Makes the general layout for WRAP
  tab_layout = [[sg.TabGroup([[sg.Tab('Single Object',   layout_single,       title_color='#F9F8F3',          background_color='#FBF5EE',   element_justification= 'center',     key = 'Single Obect Search'),
                              sg.Tab('Multi-Object',    layout_multi,        title_color='#F9F8F3',          background_color='#FBF5EE',   element_justification= 'center',     key = 'Multi-Object Search')]], 
                              tab_location='centertop', title_color='Black', tab_background_color='#F9F8F3', selected_title_color='Black', selected_background_color='#979793', border_width = 6, font = ('Times New Roman', 18), enable_events = True, key = 'tab_group'), sg.Button('Close')]] 


  #Generates the window based off the layouts above
  window = sg.Window('WRAP', tab_layout, size = (550, 650), grab_anywhere=False, finalize=True, enable_close_attempted_event = True)

    #Makes lists for: catalog names, catalog functions, and fake variables
  catalogs = ['catwise', 'AW', 'gaia', 
              'VSA', 'UKIDSS', '2MASS', 
              'ps', 'galex']
  catalogs_multi = ['catwise_multi', 'AW_multi', 'gaia_multi', 
                    'VSA_multi', 'UKIDSS_multi', '2MASS_multi', 
                    'ps_multi', 'galex_multi']
  catalog_names = ['CatWISE 2020', 'AllWISE', 'GAIA', 
                  'VSA', 'UKIDSS', '2MASS', 
                  'PanSTARRS', 'GALEX']
  catalog_functions = [catwise_image, allwise_image, gaia_image, 
                      vsa_image, ukidss_image, twomass_image, 
                      ps_image, galex_image]

  #Defines all of the variable names
  value_names = [['cw_ra', 'cw_ra_e', 'cw_dec', 'cw_dec_e', 'cw_w1', 'cw_w1_e', 'cw_w2', 'cw_w2_e', 'cw_pmra', 'cw_pmra_e', 'cw_pmdec', 'cw_pmdec_e', 'cw_mjd', 'cw_catalog', 'cw_notes'], 
                ['aw_ra', 'aw_ra_e', 'aw_dec', 'aw_dec_e', 'aw_w1', 'aw_w1_e', 'aw_w2', 'aw_w2_e', 'aw_w3', 'aw_w3_e', 'aw_w4', 'aw_w4_e', 'aw_pmra', 'aw_pmra_e', 'aw_pmdec', 'aw_pmdec_e', 'aw_catalog', 'aw_notes'], 
                ['gaia_ra', 'gaia_ra_e', 'gaia_dec', 'gaia_dec_e', 'gaia_parallax', 'gaia_parallax_e', 'gaia_radv', 'gaia_radv_e', 'gaia_pmra', 'gaia_pmra_e', 'gaia_pmdec', 'gaia_pmdec_e', 'gaia_g', 'gaia_g_e', 'gaia_bp', 'gaia_bp_e', 'gaia_rp', 'gaia_rp_e', 'gaia_year', 'gaia_catalog', 'gaia_notes'], 
                ['vsa_ra', 'vsa_dec', 'vsa_y', 'vsa_y_e', 'vsa_j', 'vsa_j_e', 'vsa_h', 'vsa_h_e', 'vsa_ks', 'vsa_ks_e', 'vsa_mjd_y', 'vsa_mjd_j', 'vsa_mjd_h', 'vsa_mjd_ks', 'vsa_catalog', 'vsa_notes'], 
                ['wfcam_ra', 'wfcam_ra_e', 'wfcam_dec', 'wfcam_dec_e', 'wfcam_y', 'wfcam_y_e', 'wfcam_j', 'wfcam_j_e', 'wfcam_h', 'wfcam_h_e', 'wfcam_k', 'wfcam_k_e', 'wfcam_pmra', 'wfcam_pmra_e', 'wfcam_pmdec', 'wfcam_pmdec_e', 'wfcam_epoch', 'wfcam_catalog', 'wfcam_notes'], 
                ['2mass_ra', '2mass_dec', '2mass_j', '2mass_j_e', '2mass_h', '2mass_h_e', '2mass_ks', '2mass_ks_e', '2mass_catalog', '2mass_notes'], 
                ['ps_ra', 'ps_ra_e', 'ps_dec', 'ps_dec_e', 'ps_g', 'ps_g_e', 'ps_r', 'ps_r_e', 'ps_i', 'ps_i_e', 'ps_z', 'ps_z_e', 'ps_y', 'ps_y_e', 'ps_mjd', 'ps_catalog', 'ps_notes'], 
                ['galex_ra', 'galex_dec', 'galex_fuv', 'galex_fuv_e', 'galex_nuv', 'galex_nuv_e', 'galex_catalog', 'galex_notes']]

  #Defines the header for the CSV file
  header = ['cw_ra', 'cw_ra_e', 'cw_dec', 'cw_dec_e', 'cw_w1', 'cw_w1_e', 'cw_w2', 'cw_w2_e', 'cw_pmra', 'cw_pmra_e', 'cw_pmdec', 'cw_pmdec_e', 'cw_mjd', 'cw_catalog', 'cw_notes', 
            'aw_ra', 'aw_ra_e', 'aw_dec', 'aw_dec_e', 'aw_w1', 'aw_w1_e', 'aw_w2', 'aw_w2_e', 'aw_w3', 'aw_w3_e', 'aw_w4', 'aw_w4_e', 'aw_pmra', 'aw_pmra_e', 'aw_pmdec', 'aw_pmdec_e', 'aw_catalog', 'aw_notes', 
            'gaia_ra', 'gaia_ra_e', 'gaia_dec', 'gaia_dec_e', 'gaia_parallax', 'gaia_parallax_e', 'gaia_radv', 'gaia_radv_e', 'gaia_pmra', 'gaia_pmra_e', 'gaia_pmdec', 'gaia_pmdec_e', 'gaia_g', 'gaia_g_e', 'gaia_bp', 'gaia_bp_e', 'gaia_rp', 'gaia_rp_e', 'gaia_year', 'gaia_catalog', 'gaia_notes', 
            'vsa_ra', 'vsa_dec', 'vsa_y', 'vsa_y_e', 'vsa_j', 'vsa_j_e', 'vsa_h', 'vsa_h_e', 'vsa_ks', 'vsa_ks_e', 'vsa_mjd_y', 'vsa_mjd_j', 'vsa_mjd_h', 'vsa_mjd_ks', 'vsa_catalog', 'vsa_notes', 
            'wfcam_ra', 'wfcam_ra_e', 'wfcam_dec', 'wfcam_dec_e', 'wfcam_y', 'wfcam_y_e', 'wfcam_j', 'wfcam_j_e', 'wfcam_h', 'wfcam_h_e', 'wfcam_k', 'wfcam_k_e', 'wfcam_pmra', 'wfcam_pmra_e', 'wfcam_pmdec', 'wfcam_pmdec_e', 'wfcam_epoch', 'wfcam_catalog', 'wfcam_notes', 
            '2mass_ra', '2mass_dec', '2mass_j', '2mass_j_e', '2mass_h', '2mass_h_e', '2mass_ks', '2mass_ks_e', '2mass_catalog', '2mass_notes', 
            'ps_ra', 'ps_ra_e', 'ps_dec', 'ps_dec_e', 'ps_g', 'ps_g_e', 'ps_r', 'ps_r_e', 'ps_i', 'ps_i_e', 'ps_z', 'ps_z_e', 'ps_y', 'ps_y_e', 'ps_mjd', 'ps_catalog', 'ps_notes', 
            'galex_ra', 'galex_dec', 'galex_fuv', 'galex_fuv_e', 'galex_nuv', 'galex_nuv_e', 'galex_catalog', 'galex_notes']

#Keeps the window open
while True:

  #Reads all of the events and values, then reads which tab is currently in
  event, values = window.read()
  group = values['tab_group']

  #Runs the program if the user is in 'Single Object' tab
  if group == 'Single Obect Search':   

    #Checks if the 'Run WRAP' button is clicked
    if event in (None, 'Run WRAP'):

      #Calls the "single_tab_check" function
      fake_list = single_tab_check()
        
      #If the RA, DEC, and RADIUS tabs are filled then runs the search on the catalogs used
      if len(fake_list) == 0: 
        sg.cprint('------------------------------------------------                                                                                                                             ', c='wheat4', end='', key = ML_KEY_SINGLE)
        sg.cprint('Running Search on RA (deg): ' + str(values['RA']) + '                                                                                                                        ', c='wheat4', end='', key = ML_KEY_SINGLE)
        sg.cprint('DEC (deg): ' + str(values['DEC']) + '                                                                                                                                        ', c='wheat4', end='', key = ML_KEY_SINGLE)
        sg.cprint('RADIUS (arcsec): ' + str(values['RADIUS']) + '                                                                                                                               ', c='wheat4', end='', key = ML_KEY_SINGLE)
        sg.cprint('------------------------------------------------                                                                                                                             ', c='wheat4', end='', key = ML_KEY_SINGLE)

        #Make variables for the RA, DEC, and RADIUS used
        ra_use, dec_use, radius_use = float(values['RA']), float(values['DEC']), int(values['RADIUS'])
        
        #Creates fake list for the data and data names
        photometry, photometry_name = [], []

        #Calls the "single_object_search" function
        single_object_search()

    #Provides the user with the authors information if the 'Help' button is pressed
    if event in (None, 'Help'):
      sg.cprint('------------------------------------------------                                                                                    ', c='wheat4',         key = ML_KEY_SINGLE)
      sg.cprint('     Thank you for using WRAP!                                                                                                      ', c='wheat4', end='', key = ML_KEY_SINGLE)
      sg.cprint('   Authors Contact: hcb98@nau.edu                                                                                                   ', c='wheat4', end='', key = ML_KEY_SINGLE)
      sg.cprint('------------------------------------------------                                                                                    ', c='wheat4',         key = ML_KEY_SINGLE)

    #Closes WRAP if the 'Close WRAP' button is pressed
    if event in (None, 'Close WRAP'):
      break

  #Runs the program if the user is in 'Multi-Object' tab
  if group == 'Multi-Object Search':

    #Checks if the 'Run WRAP' button is clicked
    if event in (None, 'Run WRAP0'):

      #Calls the "multi_tab_check" function
      multi_tab_check()

      #If the RA, DEC, and RADIUS tabs are filled then runs the search on the catalogs used
      if values['file'] != '' and values['RADIUS_multi'].isnumeric() == True and (values['type'] == 'CSV' or values['type'] == 'FITS' or values['type'] == 'ASCII' or values['type'] == 'IPAC'):
        sg.cprint('------------------------------------------------                                                                                                                                   ', c='wheat4', end='', key = ML_KEY_MULTI)
        sg.cprint('Directory: ' + str(values['file']) + '                                                                                                                                             ', c='wheat4', end='', key = ML_KEY_MULTI)
        sg.cprint('File Type: ' + str(values['type']) + '                                                                                                                                             ', c='wheat4', end='', key = ML_KEY_MULTI)
        sg.cprint('RADIUS (arcsec): ' + str(values['RADIUS_multi']) + '                                                                                                                               ', c='wheat4', end='', key = ML_KEY_MULTI)
        sg.cprint('------------------------------------------------                                                                                                                                   ', c='wheat4', end='', key = ML_KEY_MULTI)

        #Calls the "multi_tab_table" function
        ra_list, dec_list = multi_tab_table()

        #Makes a csv file and writes the header
        output = values['output2']
        
        #Makes the output file name
        if values['output2'] == '':
          output = 'WRAP_output'
        else: 
          output = values['output2']

        #Makes the CSV file and writes the header
        myFile = open('Output/' + str(output) + '.csv', 'w')
        writer = csv.writer(myFile)
        writer.writerow(header)
        
        #Calls the "multi_object_search" function
        multi_object_search()

    #Provides the user with the authors information if the 'Help' button is pressed
    if event in (None, 'Help1'):
      sg.cprint('------------------------------------------------                                                     ', end='', c='wheat4', key = ML_KEY_MULTI)
      sg.cprint('     Thank you for using WRAP!                                                                       ', end='', c='wheat4', key = ML_KEY_MULTI)
      sg.cprint('   Authors Contact: hcb98@nau.edu                                                                    ', end='', c='wheat4', key = ML_KEY_MULTI)
      sg.cprint('------------------------------------------------                                                     ', end='', c='wheat4', key = ML_KEY_MULTI)

    #Closes WRAP if the 'Close WRAP' button is pressed
    if event in (None, 'Close WRAP2'):
      break

#Closes the window
window.close()