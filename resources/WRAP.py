#-----------------------------------------------------------------------#
# WRAP.py v2.0.0
# By Hunter Brooks, at NAU/UToledo, Flagstaff: April 23, 2025
#-----------------------------------------------------------------------#



# Import all needed packages.
# ------------------------------------------------------------- #
# Import WRAP Scripts
from catalog_api import *
from image_api import *
from plotting import *

# Relavent java and dataframe packages
import argparse
import pandas as pd
import json
import sys
# ------------------------------------------------------------- #



# WRAP
# ------------------------------------------------------------- #   
# Wipes previous WRAP memory usage
sys.stdout.flush()
def wrap(ra, dec, radius, catalog_list):
    '''
    Runs API catalog and image search using astroquery and ViZieR

        Parameters:
            ra (float): The Right Accension in Degrees
            dec (float): The Declination in Degrees
            radius (float): The search radius in arcseconds 
            catalog_list (list): The list of user selected catalogs

        Returns:
            df (pandas dataframe): The user selected object data from the selected catalogs
    '''
    
    # Master Dictionary of all relavent catalog information
    catalog_info = [
        {"name": "CatWISE"  , "table_id": "II/365/catwise", 'table_data': ['RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_DE_ICRS', 'W1mproPM', 'W2mproPM', 'e_W1mproPM', 'e_W2mproPM', 'pmRA', 'pmDE', 'e_pmRA', 'e_pmDE']                                          , 'image_id': ['WISE 3.4', 'WISE 4.6']                        , 'image_selection': [True, True]                    , 'image_names': ['W1', 'W2']             , 'table_header': ['CW_RA', 'CW_DEC', 'CW_RA_E', 'CW_DEC_E', 'CW_W1', 'CW_W2', 'CW_W1_E', 'CW_W2_E', 'CW_PMRA', 'CW_PMDEC', 'CW_PMRA_E', 'CW_PMDEC_E', 'CW_NOTES']},
        {"name": "AllWISE"  , "table_id": "II/328/allwise", 'table_data': ['RAJ2000', 'DEJ2000', 'W1mag', 'W2mag', 'W3mag', 'W4mag', 'e_W1mag', 'e_W2mag', 'e_W3mag', 'e_W4mag', 'pmRA', 'pmDE', 'e_pmRA', 'e_pmDE']                                        , 'image_id': ['WISE 3.4', 'WISE 4.6', 'WISE 12', 'WISE 22']  , 'image_selection': [True, True, False, False]      , 'image_names': ['W1', 'W2', 'W3', 'W4'] , 'table_header': ['AW_RA', 'AW_DEC', 'AW_W1', 'AW_W2', 'AW_W3', 'AW_W4', 'AW_W1_E', 'AW_W2_E', 'AW_W3_E', 'AW_W4_E', 'AW_PMRA', 'AW_PMDEC', 'AW_PMRA_E', 'AW_PMDEC_E', 'AW_NOTES']},
        {"name": "Gaia"     , "table_id": "I/350/gaiaedr3", 'table_data': ['RAJ2000', 'DEJ2000', 'e_RA_ICRS', 'e_DE_ICRS', 'Gmag', 'BPmag', 'RPmag', 'e_Gmag', 'e_BPmag', 'e_RPmag', 'pmRA', 'pmDE', 'e_pmRA', 'e_pmDE', 'Plx', 'RVDR2', 'e_Plx', 'e_RVDR2'], 'image_id': ['WISE 3.4', 'WISE 4.6', 'WISE 12', 'WISE 22']  , 'image_selection': [True, True, False, False]      , 'image_names': ['W1', 'W2', 'W3', 'W4'] , 'table_header': ['GAIA_RA', 'GAIA_DEC', 'GAIA_RA_E', 'GAIA_DEC_E', 'GAIA_G', 'GAIA_BP', 'GAIA_RP', 'GAIA_G_E', 'GAIA_BP_E', 'GAIA_RP_E', 'GAIA_PMRA', 'GAIA_PMDEC', 'GAIA_PMRA_E', 'GAIA_PMDEC_E', 'GAIA_PLX', 'GAIA_RV', 'GAIA_PLX_E', 'GAIA_RV_E', 'GAIA_NOTES']},
        {"name": "VHS"    , "table_id": "II/367/vhs_dr5", 'table_data': ['RAJ2000', 'DEJ2000', 'Jap3', 'Hap3', 'Ksap3', 'e_Jap3', 'e_Hap3', 'e_Ksap3']                                                                                                    , 'image_id': 'VHS'                                         , 'image_selection': [True, False, True]             , 'image_names': ['J', 'H', 'K']          , 'table_header': ['VSA_RA', 'VSA_DEC', 'VSA_J', 'VSA_H', 'VSA_K', 'VSA_J_E', 'VSA_H_E', 'VSA_K_E', 'VSA_NOTES']}, 
        {"name": "UKIDSS"    , "table_id": "II/319"        , 'table_data': ['RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', 'Ymag', 'Jmag1', 'Jmag2', 'Hmag', 'Kmag', 'e_Ymag', 'e_Jmag1', 'e_Jmag2', 'e_Hmag', 'e_Kmag', 'pmRA', 'pmDE', 'e_pmRA', 'e_pmDE'], 'image_id': ['UKIDSS-Y', 'UKIDSS-J', 'UKIDSS-H', 'UKIDSS-K'], 'image_selection': [False, True, False, True]      , 'image_names': ['Y', 'J', 'H', 'K']     , 'table_header': ['WFCAM_RA', 'WFCAM_DEC', 'WFCAM_RA_E', 'WFCAM_DEC_E', 'WFCAM_Y', 'WFCAM_J1', 'WFCAM_J2', 'WFCAM_H', 'WFCAM_K', 'WFCAM_Y_E', 'WFCAM_J1_E', 'WFCAM_J2_E', 'WFCAM_H_E', 'WFCAM_K_E', 'WFCAM_PMRA', 'WFCAM_PMDE', 'WFCAM_PMRA_E', 'WFCAM_PMDEC_E', 'WFCAM_NOTES']},
        {"name": "2MASS"    , "table_id": "II/246/out"    , 'table_data': ['RAJ2000', 'DEJ2000', 'Jmag', 'Hmag', 'Kmag', 'Jcmsig', 'Hcmsig', 'Kcmsig']                                                                                                      , 'image_id': ['2MASS-J', '2MASS-H', '2MASS-K']               , 'image_selection': [True, False, True]             , 'image_names': ['J', 'H', 'K']          , 'table_header': ['2MASS_RA', '2MASS_DEC', '2MASS_J', '2MASS_H', '2MASS_K', '2MASS_J_E', '2MASS_H_E', '2MASS_K_E', '2MASS_NOTES']},
        {"name": "PanSTARRS", "table_id": "II/349/ps1"    , 'table_data': ['RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', 'gmag', 'rmag', 'imag', 'zmag', 'ymag', 'e_gmag', 'e_rmag', 'e_imag', 'e_zmag', 'e_ymag']                                        , 'image_id': 'PanSTARRS'                                     , 'image_selection': [True, True, False, False]      , 'image_names': ['r', 'i', 'z', 'y']     , 'table_header': ['PS_RA', 'PS_DEC', 'PS_RA_E', 'PS_DEC_E', 'PS_G', 'PS_R', 'PS_I', 'PS_Z', 'PS_Y', 'PS_G_E', 'PS_R_E', 'PS_I_E', 'PS_Z_E', 'PS_Y_E', 'PS_NOTES']}, 
        {"name": "NSC"      , "table_id": "NSC"           , 'table_data': ['ra', 'dec', 'raerr', 'decerr', 'gmag', 'rmag', 'imag', 'zmag', 'umag', 'ymag', 'gerr', 'rerr', 'ierr', 'zerr', 'uerr', 'yerr', 'pmra', 'pmdec', 'pmraerr', 'pmdecerr', 'mjd']   , 'image_id': 'NSC'                                           , 'image_selection': [False, True, True, True, False], 'image_names': ['g', 'r', 'i', 'z', 'Y'], 'table_header': ['NSC_RA', 'NSC_DEC', 'NSC_RA_E', 'NSC_DEC_E', 'NSC_G', 'NSC_R', 'NSC_I', 'NSC_Z', 'NSC_U', 'NSC_Y', 'NSC_G_E', 'NSC_R_E', 'NSC_I_E', 'NSC_Z_E', 'NSC_U_E', 'NSC_Y_E', 'NSC_PMRA', 'NSC_PMDEC', 'NSC_PMRA_E', 'NSC_PMDEC_E', 'NSC_MJD', 'NSC_NOTES']}, 
        {"name": "GALEX"    , "table_id": "II/312/ais"    , 'table_data': ['RAJ2000', 'DEJ2000', 'FUV', 'NUV', 'e_FUV', 'e_NUV']                                                                                                                            , 'image_id': ['GALEX Near UV', 'GALEX Far UV']               , 'image_selection': [True, True]                    , 'image_names': ['NUV', 'FUV']           , 'table_header': ['GALEX_RA', 'GALEX_DEC', 'GALEX_FUV', 'GALEX_NUV', 'GALEX_FUV_E', 'GALEX_NUV_E', 'GALEX_NOTES']}
    ]
    
    # Empty header list
    photometry = []
    photometry_name = []
    photometry.append([ra, dec, radius])
    photometry_name.append(['input_ra', 'input_dec', 'input_radius'])

    # Runs the catalog and image query
    for q in range(len(catalog_info)): 
        if catalog_info[q]['name'] in catalog_list: 
            table = table_query(ra, dec, radius, catalog_info[q]["table_id"]) # Table Query 
            images, w = image_query(ra, dec, radius, catalog_info[q]["image_id"]) # Image Query
            click_data = image_plot(ra, dec, radius, catalog_info[q], table, images, w) # Plotting
            
            photometry.append(click_data)
            photometry_name.append(catalog_info[q]['table_header'])
            
    # Flattens list
    photometry = [item for sublist in photometry for item in sublist]
    photometry_name = [item for sublist in photometry_name for item in sublist]
    
    # Returns data frame
    total = photometry_name + photometry
    df = pd.DataFrame([], columns=total)
    df.to_csv(sys.stdout, index=False)
    return df
# ------------------------------------------------------------- #   

# ['VHS', 'UKIDSS']
# wrap(1, 1, 150, ['NSC'])
# Makes usable with Java frontend
# ------------------------------------------------------------- #   
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ra", type=float)
    parser.add_argument("dec", type=float)
    parser.add_argument("radius", type=float)
    parser.add_argument("catalogs", nargs="+")
    
    args = parser.parse_args()
    df = wrap(args.ra, args.dec, args.radius, args.catalogs)
# ------------------------------------------------------------- #   
