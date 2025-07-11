#-----------------------------------------------------------------------#
# WRAP.image_api v2.0.0
# By Hunter Brooks, at NAU/UToledo, Flagstaff: April 23, 2025
#-----------------------------------------------------------------------#



# Import all needed packages.
# ------------------------------------------------------------- #
# Math Packages
import numpy as np

# System Packages
import sys,os
import requests
from io import BytesIO
from sys import platform
from pyvo.dal import sia
from bs4 import BeautifulSoup

# WCS Orientation Packages
from astropy.wcs import WCS
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.utils.data import download_file

# Image and Table Quering Packages
from astroquery.vsa import Vsa
from astroquery.skyview import SkyView

# NOIRLab Source Catatlog Package
if platform != 'win32':
  from dl.helpers.utils import convert
  from dl import authClient as ac, queryClient as qc
# ------------------------------------------------------------- #



# WRAP QUERY FUNCTIONS
# ------------------------------------------------------------- #
def image_query(ra, dec, radius, catalog):
  if catalog not in ['VHS', 'PanSTARRS', 'NSC']:  
    try: 
      # Clear SkyView cache and query images for non-VHS, PanSTARRS, and NSC catalogs
      SkyView.clear_cache()
      radius_deg = radius * u.arcsec # Convert radius to astropy units
      SkyView.URL = 'https://skyview.gsfc.nasa.gov/current/cgi/basicform.pl'
      images = SkyView.get_images(position=f'{ra}d {dec}d', coordinates='J2000', survey=catalog, radius=radius_deg)
      fits_header = images[0][0].header # Extract FITS header from the first image
      w = WCS(fits_header) # Create WCS object from the FITS header
      
      image_cutout = []
      position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5')
      size = u.Quantity([radius, radius], u.arcsec)
      for image in images: 
        image_cutout.append((Cutout2D(image[0].data, position, size, fill_value = np.nan, wcs = w.celestial)).data)
        wcs_cropped = (Cutout2D(image[0].data, position, size, fill_value = np.nan, wcs = w.celestial)).wcs  
      return image_cutout, wcs_cropped
    except: 
      return 0, 0

  elif catalog == 'VHS':
    # try:
      position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5')
      size = u.Quantity([radius, radius], u.arcsec)
      
      # Vsa.clear_cache()
      Vsa.BASE_URL = 'http://vsa.roe.ac.uk:8080/vdfs/'
      url_J = Vsa.get_image_list(SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs'), image_width=radius * u.arcsec, waveband='J', database='VHSDR4')
      url_H = Vsa.get_image_list(SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs'), image_width=radius * u.arcsec, waveband='H', database='VHSDR4')
      url_K = Vsa.get_image_list(SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs'), image_width=radius * u.arcsec, waveband='Ks', database='VHSDR4')

      temp_J = url_J[0]
      file_vsa_J = download_file(temp_J, cache = False)
      data_vsa_J = fits.getdata(file_vsa_J)
      
      hdu_j = fits.open(file_vsa_J)[1]
      w = WCS(hdu_j.header)
      
      cutout_j = (Cutout2D(data_vsa_J, position, size, wcs = w.celestial)).data
      
      if len(url_H) == 0: 
        temp_H = ''
        cutout_h = np.zeros_like(cutout_j)
      else: 
        temp_H = url_H[0]
        file_vsa_H = download_file(temp_H, cache = False)
        data_vsa_H = fits.getdata(file_vsa_H)
        cutout_h = Cutout2D(data_vsa_H, position, size, wcs = w.celestial)
        
      if len(url_K) == 0: 
        temp_K = ''
        cutout_k = np.zeros_like(cutout_j)
      else: 
        temp_K = url_K[0]
        file_vsa_K = download_file(temp_K, cache = False)
        data_vsa_K = fits.getdata(file_vsa_K)
        cutout_ks = Cutout2D(data_vsa_K, position, size, wcs = w.celestial)

      return [cutout_j, cutout_h, cutout_ks], w
    # except: 
      # return 0, 0

  elif catalog == 'PanSTARRS':
    try:
      # Construct URLs for PanSTARRS image retrieval
      new_dec = f'+{dec}' if dec > 0 else str(dec)
      ps_image_url = f'http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos={ra}{new_dec}&filter=color&filter=g&filter=r&filter=i&filter=z&filter=y&filetypes=stack&auxiliary=data&size={int(radius)*4}&output_size=0&verbose=0&autoscale=99.500000&catlist='
      
      # Download metadata from the constructed URL
      allwise_metadata = requests.get(ps_image_url)
      script_dir = os.path.dirname(os.path.abspath(__file__))
      script_dir = script_dir.split('/WRAP.app/Contents/app/resources')[0]
      metadata_path = f'{script_dir}/ps_metadata.txt'
      open(metadata_path, 'wb').write(allwise_metadata.content)

      # Parse metadata to find image links
      r_finder = 'amp;'
      r_finder_list = []
      with open(metadata_path, 'r') as fp:
        lines = fp.readlines()
        r_finder_list = [line for line in lines if r_finder in line]
      
      # Extract image URLs from the metadata
      r_link = ('http:' + (r_finder_list[2]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)
      i_link = ('http:' + (r_finder_list[3]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)
      z_link = ('http:' + (r_finder_list[4]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)
      y_link = ('http:' + (r_finder_list[5]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)
      
      # Download FITS files and get data for r, i, z, and y bands
      file_ps_r, file_ps_i, file_ps_z, file_ps_y = download_file(r_link, cache=True), download_file(i_link, cache=True), download_file(z_link, cache=True), download_file(y_link, cache=True)
      data_ps_r, data_ps_i, data_ps_z, data_ps_y = fits.getdata(file_ps_r), fits.getdata(file_ps_i), fits.getdata(file_ps_z), fits.getdata(file_ps_y)
      
      hdu_r = fits.open(file_ps_r)[0]
      w = WCS(hdu_r.header)  # Create WCS object from the r band FITS header
      
      # Remove temporary metadata file
      os.remove(metadata_path)
      return [data_ps_r, data_ps_i, data_ps_z, data_ps_y], w
    except:
      os.remove(metadata_path)
      return 0, 0

  elif catalog == 'NSC':
    try: 
      #Defines the catalog that is searched
      DEF_ACCESS_URL = "https://datalab.noirlab.edu/sia/des_dr2"
      svc = sia.SIAService(DEF_ACCESS_URL)

      #Finds all of the image urls for the ra, dec, and radius given
      imgTable = svc.search((ra,dec), float(radius/3600)).to_table()

      #Tests if any images were found
      if len(imgTable) > 0:
        #Obtains all of the image url, bands, and types from the table
        image_urls = imgTable['access_url']
        band_passes = imgTable['obs_bandpass']
        type = imgTable['prodtype']

        #Obtains the r band image
        for p in range(len(band_passes)):
          if band_passes[p] == 'g' and type[p] == 'image':
            image_url_g = image_urls[p]
            break

        #Obtains the r band image
        for p in range(len(band_passes)):
          if band_passes[p] == 'r' and type[p] == 'image':
            image_url_r = image_urls[p]
            break

        #Obtains the i band image
        for p in range(len(band_passes)):
          if band_passes[p] == 'i' and type[p] == 'image':
            image_url_i = image_urls[p]
            break

        #Obtains the z band image
        for p in range(len(band_passes)):
          if band_passes[p] == 'z' and type[p] == 'image':
            image_url_z = image_urls[p]
            break

        #Obtains the z band image
        for p in range(len(band_passes)):
          if band_passes[p] == 'Y' and type[p] == 'image':
            image_url_Y = image_urls[p]
            break

        #Downloads the images
        file_allwise_g, file_allwise_r, file_allwise_i, file_allwise_z, file_allwise_Y = download_file(image_url_g, cache=True), download_file(image_url_r, cache=True), download_file(image_url_i, cache=True), download_file(image_url_z, cache=True), download_file(image_url_Y, cache=True)
        data_allwise_g, data_allwise_r, data_allwise_i, data_allwise_z, data_allwise_Y = fits.getdata(file_allwise_g), fits.getdata(file_allwise_r), fits.getdata(file_allwise_i), fits.getdata(file_allwise_z), fits.getdata(file_allwise_Y)

        #Loads the WCS from the i band image
        hdu_r = fits.open(file_allwise_r)[0]
        wcs = WCS(hdu_r.header)
        return [data_allwise_g, data_allwise_r, data_allwise_i, data_allwise_z, data_allwise_Y], wcs
      else: 
        return 0, 0
    except: 
      return 0, 0
# ------------------------------------------------------------- #
