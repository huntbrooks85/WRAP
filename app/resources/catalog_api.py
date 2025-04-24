#-----------------------------------------------------------------------#
# WRAP.catalog_api v2.0.0
# By Hunter Brooks, at NAU/UToledo, Flagstaff: April 23, 2025
#-----------------------------------------------------------------------#



# Import all needed packages.
# ------------------------------------------------------------- #
# System Packages
from sys import platform

# WCS Orientation Packages
from astropy import units as u
from astropy.coordinates import SkyCoord

# Image and Table Quering Packages
from astroquery.vizier import Vizier

# NOIRLab Source Catatlog Package
if platform != 'win32':
  from dl.helpers.utils import convert
  from dl import authClient as ac, queryClient as qc
# ------------------------------------------------------------- #



# WRAP TABLE QUERY FUNCTIONS
# ------------------------------------------------------------- #
def table_query(ra, dec, radius, catalog):
  if catalog != 'NSC':
    try:
      # Set row limit for Vizier query
      coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='fk5')
      search_radius = radius * u.arcsec  # Convert radius to astropy units
      result = Vizier(columns=['**'], row_limit = 1e10).query_region(coord, width=[search_radius, search_radius], catalog=catalog, frame = 'galactic') # Query Vizier catalog with all columns
      if len(result) == 0:
        return 0  # Return 0 if no results are found
      return result
    except:
      return 0  # Return 0 in case of any exceptions

  elif catalog == 'NSC':
    try:
      # Construct SQL query to search NSC database within specified radius
      query = " \
      SELECT ra, dec, gmag, gerr, rmag, rerr, imag, ierr, zmag, zerr, umag, uerr, ymag, yerr, pmra, pmraerr, pmdec, pmdecerr, raerr, decerr, mjd  \
      FROM nsc_dr2.object \
      WHERE ra > " + str((ra - ((radius / 7200) - 0.0001))) + " and ra < " + str((ra + ((radius / 7200) - 0.0001))) + " " \
      "AND dec > " + str((dec - ((radius / 7200) - 0.0001))) + " and dec < " + str((dec + ((radius / 7200) - 0.0001))) + " "
      
      # Execute the query and convert the response to a pandas DataFrame
      response = qc.query(sql=query, format='csv')
      df = convert(response, 'pandas')
      return df.head(100000000000)  # Return the first 100000000000 rows of the DataFrame (adjust as necessary)
    except:
      return 0  # Return 0 in case of any exceptions
# ------------------------------------------------------------- #
