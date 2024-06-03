#-----------------------------------------------------------------------#
# WRAP v1.5.0
# By Hunter Brooks, at NAU, Flagstaff: May 27, 2024
#
# Purpose: Gathers photometry and astrometry from various 
#          ultra-violet, optical, and near-infrared catalogs 
#          given a RA, DEC, and Radius by the user
#-----------------------------------------------------------------------#



# Import all needed packages.
# ------------------------------------------------------------- #
# GUI Packages
import PySimpleGUI as sg

# Math Packages
import math
import numpy as np
import pandas as pd

# System Packages
import csv
import sys,os
import requests
import webbrowser
from sys import platform
from pyvo.dal import sia
from bs4 import BeautifulSoup

# Image Plotting Packagaes
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
from matplotlib.widgets import Slider, Button, CheckButtons, TextBox

# WCS Orientation Packages
from astropy.wcs import WCS
from astropy.io import fits
from astropy.io import ascii
from astropy.time import Time
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.utils.data import download_file

# Image and Table Quering Packages
from astroquery.vsa import Vsa
from astroquery.vizier import Vizier
from astroquery.skyview import SkyView


# NOIRLab Source Catatlog Package
if platform != 'win32':
  from dl.helpers.utils import convert
  from dl import authClient as ac, queryClient as qc
# ------------------------------------------------------------- #



# PRINTING FUNCTIONS
# ------------------------------------------------------------- #
def blockPrint():
  '''Makes a function that blocks the printing function'''
  sys.stdout = open(os.devnull, 'w')

def enablePrint():
  '''Makes a function that allows the printing function'''
  sys.stdout = sys.__stdout__
# ------------------------------------------------------------- #



# WRAP PRINTING OUTPUT FUNCTIONS
# ------------------------------------------------------------- #
def wrap_start(catalog, tab):
  '''The print function for when a catalog search has started'''

  print('#------------------------------------------------#')
  print('Started ' + str(catalog) + ' Search')
  print('#------------------------------------------------#')

def wrap_found(catalog, tab):
  '''Print function for when the catalog search was successful.'''

  print('#------------------------------------------------#')
  print('Finished ' + str(catalog) + ' Search')
  print('#------------------------------------------------#')

def wrap_not_found(catalog, tab):
  '''Print function for when the catalog search was not successful.'''

  print('#------------------------------------------------#')
  print('Object Not Found')
  print('Finished ' + str(catalog) + ' Search')
  print('#------------------------------------------------#')

def wrap_end(tab):
  '''The print function for when all catalogs have been searched'''

  print('')
  print('#------------------------------------------------#')
  print('All Catalogs Have Been Searched')
  print('Finished Running WRAP')
  print('#------------------------------------------------#')
  print('')
# ------------------------------------------------------------- #



# WRAP QUERY FUNCTIONS
# ------------------------------------------------------------- #
def image_query(ra, dec, radius, catalog):
  blockPrint()  # Suppress print statements during execution

  if catalog not in ['VSA', 'PS2', 'NSC']:  
    try: 
      # Clear SkyView cache and query images for non-VHS, PS2, and NSC catalogs
      SkyView.clear_cache()
      radius_deg = radius * u.arcsec # Convert radius to astropy units
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

  elif catalog == 'VSA':
    try:
      position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5')
      size = u.Quantity([radius, radius], u.arcsec)
      
      Vsa.clear_cache()
      url_J = Vsa.get_image_list(SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs'), image_width=radius * u.arcsec, waveband='J', database='VHSDR5')
      url_H = Vsa.get_image_list(SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs'), image_width=radius * u.arcsec, waveband='H', database='VHSDR5')
      url_K = Vsa.get_image_list(SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs'), image_width=radius * u.arcsec, waveband='Ks', database='VHSDR5')
      
      temp_J = url_J[0].replace("http://horus.roe.ac.uk/wsa/cgi-bin/getFImage.cgi?file=", "http://vsa.roe.ac.uk/cgi-bin/getImage.cgi?file=")
      response_J = requests.get(temp_J)
      soup_J =  BeautifulSoup(response_J.content, 'html.parser')
      link_tag_J = soup_J.find('a', href=True, text="download FITS file")
      fits_link_J = link_tag_J['href']
      file_vsa_J = download_file(fits_link_J, cache = False)
      data_vsa_J = fits.getdata(file_vsa_J)
      
      hdu_j = fits.open(file_vsa_J)[1]
      w = WCS(hdu_j.header)
      
      cutout_j = (Cutout2D(data_vsa_J, position, size, wcs = w.celestial)).data
      
      if len(url_H) == 0: 
        temp_H = ''
        cutout_h = np.zeros_like(cutout_j)
      else: 
        temp_H = url_H[0].replace("http://horus.roe.ac.uk/wsa/cgi-bin/getFImage.cgi?file=", "http://vsa.roe.ac.uk/cgi-bin/getImage.cgi?file=")
        response_H = requests.get(temp_H)
        soup_H =  BeautifulSoup(response_H.content, 'html.parser')
        link_tag_H = soup_H.find('a', href=True, text="download FITS file")
        fits_link_H = link_tag_H['href']
        file_vsa_H = download_file(fits_link_H, cache = False)
        data_vsa_H = fits.getdata(file_vsa_H)
        cutout_h = Cutout2D(data_vsa_H, position, size, wcs = w.celestial)
        
      if len(url_K) == 0: 
        temp_K = ''
        cutout_k = np.zeros_like(cutout_j)
      else: 
        temp_K = url_K[0].replace("http://horus.roe.ac.uk/wsa/cgi-bin/getFImage.cgi?file=", "http://vsa.roe.ac.uk/cgi-bin/getImage.cgi?file=")
        response_K = requests.get(temp_K)
        soup_K =  BeautifulSoup(response_K.content, 'html.parser')
        link_tag_K = soup_K.find('a', href=True, text="download FITS file")
        fits_link_K = link_tag_K['href']
        file_vsa_K = download_file(fits_link_K, cache = False)
        data_vsa_K = fits.getdata(file_vsa_K)
        cutout_ks = Cutout2D(data_vsa_K, position, size, wcs = w.celestial)

      return [cutout_j, cutout_h, cutout_ks], w
    except: 
      return 0, 0

  elif catalog == 'PS2':
    try:
      # Construct URLs for PS2 image retrieval
      new_dec = f'+{dec}' if dec > 0 else str(dec)
      ps_image_url = f'http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos={ra}{new_dec}&filter=color&filter=g&filter=r&filter=i&filter=z&filter=y&filetypes=stack&auxiliary=data&size={radius * 4}&output_size=0&verbose=0&autoscale=99.500000&catlist='
      
      # Download metadata from the constructed URL
      allwise_metadata = requests.get(ps_image_url)
      metadata_path = f'{directory}Output/metadata/ps_metadata.txt'
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

      # Define function to convert MJD to ISO date
      def mjd_to_date(mjd):
        time = Time(mjd, format='mjd')
        iso_date = time.iso
        return iso_date
      
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
      imgTable = svc.search((ra,dec), (radius/3600)).to_table()
      
      print(imgTable)

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

def table_query(ra, dec, radius, catalog):
  blockPrint()  # Suppress print statements during execution

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



# WRAP PLOTTING FUNCTIONS
# ------------------------------------------------------------- #   
def image_plot(ra, dec, radius, catalog_info): 
  # Perform image and table queries
  images, w = image_query(ra, dec, radius, catalog_info['image_id'])
  table = table_query(ra, dec, radius, catalog_info['table_id'])
  enablePrint()  # Enable printing

  # Check if queries returned valid results
  if type(images) != int and type(w) != int and type(table) != int: 
    try: 
      # Configure plot appearance
      plt.rcParams['toolbar'] = 'None'
      plt.style.use('Solarize_Light2')
      plt.rcParams["figure.figsize"] = [8, 8]

      # Create a new figure with WCS projection
      fig_1, ax = plt.subplots(subplot_kw={'projection': w})
      fig_1.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1)

      # Hide axis ticks and labels
      plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
      plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
      plt.suptitle(f'{catalog_info["name"]} Search', fontsize=35, y=0.98, fontfamily='Times New Roman')
      
      plt.grid(linewidth = 0)
      figure = plt.gcf()
      figure.set_size_inches(4.75, 6)

      # Extract object coordinates from the table
      try: 
        object_ra = table[0][catalog_info['table_data'][0]].tolist()
        object_dec = table[0][catalog_info['table_data'][1]].tolist()
      except:
        object_ra = table[catalog_info['table_data'][0]].tolist()
        object_dec = table[catalog_info['table_data'][1]].tolist()      
      
      # Combine image data
      default = catalog_info['image_selection']
      if isinstance(images[0], np.ndarray): 
        total_data = np.zeros_like(images[0].data)
        for index in range(len(images)):
          if default[index] == True: 
            total_data += images[index].data
        max_shape = np.nanmax(total_data.shape)
      else:
        total_data = np.zeros_like(images[0][0].data)
        for image in images:
          index = images.index(image)
          if default[index] == True: 
            total_data += image[0].data
        max_shape = np.nanmax(total_data.shape)
              
      # Set initial circle size for scatter plot
      circle_size = (radius * 3)
      scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('icrs'), s=circle_size, edgecolor='#40E842', facecolor='none')

      # Set initial normalization for image display
      init_bot, init_top = 45, 95
      norm1_total = matplotlib.colors.Normalize(vmin=np.nanpercentile(total_data.data, init_bot), vmax=np.nanpercentile(total_data.data, init_top))
      ax.imshow(total_data.data, cmap='Greys', norm=norm1_total)
      
      if catalog_info['name'] == 'VSA': 
        plt.xlim(0, max_shape), plt.ylim(max_shape, 0) 
      elif catalog_info['name'] == 'NSC': 
        plt.xlim(max_shape, 0), plt.ylim(0, max_shape) 
      else: 
        plt.xlim(0, max_shape), plt.ylim(0, max_shape) 

      # Create a cursor for the plot
      cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
      
      # Create an annotation for mouse events
      annotation = ax.annotate('', xy=(0, 0), xytext=(20, -20), arrowprops=dict(arrowstyle='wedge'), fontsize=12, color='red')
      annotation.set_visible(False)

      # Create sliders for adjusting image stretch
      freq_top = plt.axes([0.25, 0.12, 0.65, 0.03])
      slider_top = Slider(ax=freq_top, label='Top Stretch:', valmin=50, valmax=100, valinit=init_top, color='#E48671')
      freq_bottom = plt.axes([0.25, 0.087, 0.65, 0.03])
      slider_bottom = Slider(ax=freq_bottom, label='Bottom Stretch:', valmin=0, valmax=50, valinit=init_bot, color='#E48671')

      # Create slider for adjusting circle size
      circle_slid_location = plt.axes([0.25, 0.055, 0.65, 0.03])
      circle_slider = Slider(ax=circle_slid_location, label='Circle Size:', valmin=(circle_size - 2.5*radius), valmax=(circle_size + 1*radius), valinit=circle_size, color='#E48671')

      # Create a text box for notes
      axbox = plt.axes([0.15, 0.02, 0.8, 0.03])
      text = ''
      text_box = TextBox(axbox, 'Notes:', initial=text, textalignment="center")

      # Create a button for "Object Not Found"
      axes_button = plt.axes([0.04, 0.855, 0.92, 0.04])
      close = Button(axes_button, 'Object Not Found', color='#E48671')
      
      #Make checkbuttons with all of the different image bands
      rax = plt.axes([0.045, 0.4, 0.105, 0.12])
      labels = catalog_info['image_names']
      real_data = []
      if isinstance(images[0], np.ndarray): 
        for image in images:
          real_data.append(image.data)
      else:
        for image in images:
          real_data.append(image[0].data)
      check = CheckButtons(rax, labels, default)

      # Track mouse click locations
      location = []
      def mouse_event(event):
        location.append(event.xdata)
        location.append(event.ydata)
        location.append(event.inaxes)
      cid = fig_1.canvas.mpl_connect('button_press_event', mouse_event)
      
      # Update image stretch based on slider values
      def update_slider_stretch(val):
        if isinstance(images[0], np.ndarray):
          total_data = np.zeros_like(images[0].data)
        else:
          total_data = np.zeros_like(images[0][0].data)       
          
        for tf_ind in default: 
          if tf_ind == False: 
            pass
          if tf_ind == True: 
            index = default.index(tf_ind)
            if isinstance(images[0], np.ndarray):
              total_data += images[index].data
            else:
              total_data += images[index][0].data
                    
        norm1_w1 = matplotlib.colors.Normalize(vmin=np.nanpercentile(total_data.data, slider_bottom.val), vmax=np.nanpercentile(total_data.data, slider_top.val))
        ax.imshow(total_data.data, cmap='Greys', norm=norm1_w1)

      # Update text from text box
      text_list = [text]
      def submit(expression):
        text = expression
        text_list.append(text)
        
      #Update the image depending on what the user chooses
      def update_button(label):
        '''Updates the list of activated images and updates the image the user can see'''

        total_data = 0
        for lab in labels:
          if lab == label:
            index = labels.index(lab)
            if default[index] == False:
              default[index] = True
            elif default[index] == True: 
              default[index] = False     
              
        if isinstance(images[0], np.ndarray):
          total_data = np.zeros_like(images[0].data)
        else:
          total_data = np.zeros_like(images[0][0].data)       
          
        for tf_ind in default: 
          if tf_ind == False: 
            pass
          if tf_ind == True: 
            index = default.index(tf_ind)
            if isinstance(images[0], np.ndarray):
              total_data += images[index].data
            else:
              total_data += images[index][0].data
        norm1_total = matplotlib.colors.Normalize(vmin=np.nanpercentile(total_data.data, init_bot), vmax=np.nanpercentile(total_data.data, init_top))
        ax.imshow(total_data.data, cmap='Greys', norm=norm1_total)

      # Connect sliders and text box to their update functions
      slider_top.on_changed(update_slider_stretch)
      slider_bottom.on_changed(update_slider_stretch)
      text_box.on_text_change(submit)
      check.on_clicked(update_button)

      # Main loop to handle user interactions
      n = -1
      while True: 
        press = plt.waitforbuttonpress()
        text_max = len(text_list) - 1
        
        if press == False:
          n += 3

          # Find which axes was clicked
          click_axes = str(location[n])
          click_axes = click_axes.split('WCSAxesSubplot', 2)[0]
          
          # Handle "Object Not Found" button click
          if click_axes == 'Axes(0.04,0.855;0.92x0.04)':
            next_window()
            return len(catalog_info['table_data'])*[np.nan] + [text_list[text_max]]
          
          # Handle circle size slider adjustment
          if click_axes == 'Axes(0.25,0.055;0.65x0.03)': 
            scatter.remove()
            scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('icrs'), s=circle_slider.val, edgecolor='#40E842', facecolor='none')
          
          # Handle main plot click
          if click_axes == '': 
            # Find the closest point to the location clicked
            coord = w.pixel_to_world_values(location[len(location) - 3], location[len(location) - 2])
            distance = []
            for i in range(len(object_ra)):
              distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))
            list_location = distance.index(np.nanmin(distance))
            table_data = []
            for col_name in catalog_info['table_data']: 
              try: 
                col_data = table[0][f'{col_name}'].tolist()
                table_data.append(col_data[list_location])
              except: 
                col_data = table[f'{col_name}'].tolist()
                table_data.append(col_data[list_location])              
            next_window()
            return table_data + [text_list[text_max]]
        
        elif press is None:
          next_window()
          return len(catalog_info['table_data'])*[np.nan] + [text_list[text_max]]
    except Exception as e: 
      print('#------------------------------------------------#')
      print("Please Report This Error to GitHub:", e)
      print('#------------------------------------------------#')
      return len(catalog_info['table_data'])*[np.nan] + ['Catalog Data Not Retrieved']
  else: 
    return len(catalog_info['table_data'])*[np.nan] + ['Catalog Data Not Retrieved']

def next_window(): 
  # Clear the current figure and close all plots
  plt.clf(), plt.close('all')
  
  # Create a new figure
  plt.figure(1)
  
  # Display a message in the plot
  plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
  
  # Set plot limits and disable the grid
  plt.xlim(0, 1), plt.ylim(0, 1), plt.grid(linewidth=0)
  
  # Get current axis and hide the tick labels and ticks
  ax = plt.gca()
  ax.xaxis.set_tick_params(labelbottom=False)
  ax.yaxis.set_tick_params(labelleft=False)
  ax.set_xticks([])
  ax.set_yticks([])
  
  # Get the current figure and set its size
  figure2 = plt.gcf()
  figure2.set_size_inches(5.2, 1)
  
  # Pause briefly to display the figure, then clear and close it
  plt.pause(0.1), plt.clf(), plt.close('all')
  return
# ------------------------------------------------------------- #



# WRAP RUNNING FUNCTIONS
# ------------------------------------------------------------- #
def wiseview_link(ra, dec, radius): 
  '''Opens WISEView link to aid the user in finding their object'''

  #Gets the link for SMDET
  link = 'http://byw.tools/wiseview#ra=' + str(ra) + '&dec=' + str(dec) + '&size=' + str(radius) + '&band=3&speed=164.06&minbright=-10.0000&maxbright=80&window=0.75&diff_window=1&linear=1&color=&zoom=8.5&border=1&gaia=1&invert=1&maxdyr=0&scandir=0&neowise=0&diff=0&outer_epochs=0&unique_window=1&smooth_scan=0&shift=0&pmra=0&pmdec=0&synth_a=0&synth_a_sub=0&synth_a_ra=&synth_a_dec=&synth_a_w1=&synth_a_w2=&synth_a_pmra=0&synth_a_pmdec=0&synth_a_mjd=&synth_b=0&synth_b_sub=0&synth_b_ra=&synth_b_dec=&synth_b_w1=&synth_b_w2=&synth_b_pmra=0&synth_b_pmdec=0&synth_b_mjd=&smdet_coadd_id=1863p620&smdet_mask_idx=3&smdet_obj_px=&smdet_obj_py='
  webbrowser.open(link)
  
def single_object_search(): 
  '''Performs the single object search by going to each catalog script and recording the data down.
  Once all the data is recorded it it put into a CSV file for the user.'''

  #Runs the wiseview_link function
  if values['wiseview'] == True: 
    wiseview_link(ra_use, dec_use, radius_use)

  #Creates fake list for the data and data names
  photometry = []
  photometry_name = []
  photometry.append([ra_use, dec_use, radius_use])
  photometry_name.append(['input_ra', 'input_dec', 'input_radius'])

  for q in range(len(catalog_info)): 

    #Calls the catalog and records the data down
    if values[f"SINGLE_{catalog_info[q]['name']}"] == True: 
      wrap_start(catalog_info[q]['name'], ML_KEY_SINGLE)
      catalog_result = image_plot(ra_use, dec_use, radius_use, catalog_info[q])
      photometry.append(catalog_result)
      photometry_name.append(catalog_info[q]['table_header'])

      #Tells the user if the object was found or not
      if catalog_result == 0: 
        wrap_end(ML_KEY_SINGLE)
      else: 
        wrap_found(catalog_info[q]['name'], ML_KEY_SINGLE)

    #Writes all of the data gathered into a csv file 
    if q == 8:
      #Makes the output file name
      if values['output'] == '':
        output = 'WRAP_output'
      else: 
        output = values['output']

      #Writes the CSV file with the photometry and astrometry gathered
      new_directory = (directory.split('Output/')[0]) + 'Output/'
      myFile = open(str(new_directory) + str(output) + '.csv', 'w')
      writer = csv.writer(myFile)
      flat_photometry_list = [item for sublist in photometry for item in sublist]
      flat_photometry_name_list = [item for sublist in photometry_name for item in sublist]
      writer.writerow(flat_photometry_name_list)
      writer.writerow(flat_photometry_list)
      myFile.close()
  wrap_end(ML_KEY_SINGLE)

def single_tab_check(): 
  '''Checks if the User put in the correct formats for the RA, DEC, and Radius options.'''

  #Checks if the RA tab is entered with a number
  fake_list = []
  try: 
    float(values['RA'])
  except ValueError:
    print('#------------------------------------------------#')
    print('Please enter a Correct RA!')
    print('#------------------------------------------------#')
    fake_list.append(1)

  #Checks if the DEC tab is entered with a number
  try: 
    float(values['DEC'])
  except ValueError:
    print('#------------------------------------------------#')
    print('Please enter a Correct DEC!')
    print('#------------------------------------------------#')
    fake_list.append(2)

  #Checks if the RADIUS tab is entered with a number
  if values['RADIUS'].isnumeric() == False or (int(values['RADIUS']) > 500) == True or (int(values['RADIUS']) < 100) == True:
    print('#------------------------------------------------#')
    print('Please enter a Correct radius!')
    print('#------------------------------------------------#')
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
    photometry.append([ra_use, dec_use, radius_use])
    photometry_name.append(['input_ra', 'input_dec', 'input_radius'])

    #Runs the wiseview_link function
    if values['wiseview_multi'] == True: 
      wiseview_link(ra_use, dec_use, radius_use)

    for q in range(len(catalog_info)): 

      #Calls the catalog and records the data down
      if values[f"MULTI_{catalog_info[q]['name']}"] == True: 
        wrap_start(catalog_info[q]['name'], ML_KEY_SINGLE)
        catalog_result = image_plot(ra_use, dec_use, radius_use, catalog_info[q])
        photometry.append(catalog_result)
        photometry_name.append(catalog_info[q]['table_header'])

        #Tells the user if the object was found or not
        if catalog_result == 0: 
          wrap_end(ML_KEY_SINGLE)
        else: 
          wrap_found(catalog_info[q]['name'], ML_KEY_SINGLE) 

      #Writes all of the data gathered into a csv file
      if q == 8:
        if index == 0: 
          #Makes the CSV file and writes the header
          flat_photometry_name_list = [item for sublist in photometry_name for item in sublist]
          myFile = open('Output/' + str(output) + '.csv', 'w')
          writer = csv.writer(myFile)
          writer.writerow(flat_photometry_name_list)
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
    print('#------------------------------------------------#')
    print('Please enter a correct file directory!')
    print('#------------------------------------------------#')

  #Checks if the RADIUS tab is entered
  if values['RADIUS_multi'].isnumeric() == False or (int(values['RADIUS_multi']) > 500) == True or (int(values['RADIUS_multi']) < 100) == True:
    print('#------------------------------------------------#')
    print('Please enter a correct radius value!')
    print('#------------------------------------------------#')

  #Checks if the FILETYPE tab is entered 
  if values['type'] == '': 
    print('#------------------------------------------------#')
    print('Please enter a correct file type!')
    print('#------------------------------------------------#')

def multi_tab_table(): 
  '''Reads in the file the user chose for the multi-object search'''
  try:
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
  except: 
    print('#------------------------------------------------#')
    print('Please enter a correct file type!')
    print('#------------------------------------------------#')
# ------------------------------------------------------------- #



# GENERAL WRAP INFORMATION 
# ------------------------------------------------------------- #
# Sets General WRAP Theme
sg.theme('LightBrown3')
ML_KEY_SINGLE = '-ML-'  + sg.WRITE_ONLY_KEY
ML_KEY_MULTI  = '-ML2-' + sg.WRITE_ONLY_KEY

# Relevent Catalog Information for Query
catalog_info = [
  {"name": "CatWISE"  , "table_id": "II/365/catwise", 'table_data': ['RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_DE_ICRS', 'W1mproPM', 'W2mproPM', 'e_W1mproPM', 'e_W2mproPM', 'pmRA', 'pmDE', 'e_pmRA', 'e_pmDE']                                          , 'image_id': ['WISE 3.4', 'WISE 4.6']                        , 'image_selection': [True, True]                    , 'image_names': ['W1', 'W2']             , 'table_header': ['CW_RA', 'CW_DEC', 'CW_RA_E', 'CW_DEC_E', 'CW_W1', 'CW_W2', 'CW_W1_E', 'CW_W2_E', 'CW_PMRA', 'CW_PMDEC', 'CW_PMRA_E', 'CW_PMDEC_E', 'CW_NOTES']},
  {"name": "AllWISE"  , "table_id": "II/328/allwise", 'table_data': ['RAJ2000', 'DEJ2000', 'W1mag', 'W2mag', 'W3mag', 'W4mag', 'e_W1mag', 'e_W2mag', 'e_W3mag', 'e_W4mag', 'pmRA', 'pmDE', 'e_pmRA', 'e_pmDE']                                        , 'image_id': ['WISE 3.4', 'WISE 4.6', 'WISE 12', 'WISE 22']  , 'image_selection': [True, True, False, False]      , 'image_names': ['W1', 'W2', 'W3', 'W4'] , 'table_header': ['AW_RA', 'AW_DEC', 'AW_W1', 'AW_W2', 'AW_W3', 'AW_W4', 'AW_W1_E', 'AW_W2_E', 'AW_W3_E', 'AW_W4_E', 'AW_PMRA', 'AW_PMDEC', 'AW_PMRA_E', 'AW_PMDEC_E', 'AW_NOTES']},
  {"name": "Gaia"     , "table_id": "I/350/gaiaedr3", 'table_data': ['RAJ2000', 'DEJ2000', 'e_RA_ICRS', 'e_DE_ICRS', 'Gmag', 'BPmag', 'RPmag', 'e_Gmag', 'e_BPmag', 'e_RPmag', 'pmRA', 'pmDE', 'e_pmRA', 'e_pmDE', 'Plx', 'RVDR2', 'e_Plx', 'e_RVDR2'], 'image_id': ['WISE 3.4', 'WISE 4.6', 'WISE 12', 'WISE 22']  , 'image_selection': [True, True, False, False]      , 'image_names': ['W1', 'W2', 'W3', 'W4'] , 'table_header': ['GAIA_RA', 'GAIA_DEC', 'GAIA_RA_E', 'GAIA_DEC_E', 'GAIA_G', 'GAIA_BP', 'GAIA_RP', 'GAIA_G_E', 'GAIA_BP_E', 'GAIA_RP_E', 'GAIA_PMRA', 'GAIA_PMDEC', 'GAIA_PMRA_E', 'GAIA_PMDEC_E', 'GAIA_PLX', 'GAIA_RV', 'GAIA_PLX_E', 'GAIA_RV_E', 'GAIA_NOTES']},
  {"name": "VSA"      , "table_id": "II/367/vhs_dr5", 'table_data': ['RAJ2000', 'DEJ2000', 'Jap3', 'Hap3', 'Ksap3', 'e_Jap3', 'e_Hap3', 'e_Ksap3']                                                                                                    , 'image_id': 'VSA'                                           , 'image_selection': [True, False, True]             , 'image_names': ['J', 'H', 'K']          , 'table_header': ['VSA_RA', 'VSA_DEC', 'VSA_J', 'VSA_H', 'VSA_K', 'VSA_J_E', 'VSA_H_E', 'VSA_K_E', 'VSA_NOTES']}, 
  {"name": "WFCAM"    , "table_id": "II/319"        , 'table_data': ['RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', 'Ymag', 'Jmag1', 'Jmag2', 'Hmag', 'Kmag', 'e_Ymag', 'e_Jmag1', 'e_Jmag2', 'e_Hmag', 'e_Kmag', 'pmRA', 'pmDE', 'e_pmRA', 'e_pmDE'], 'image_id': ['UKIDSS-Y', 'UKIDSS-J', 'UKIDSS-H', 'UKIDSS-K'], 'image_selection': [False, True, False, True]      , 'image_names': ['Y', 'J', 'H', 'K']     , 'table_header': ['WFCAM_RA', 'WFCAM_DEC', 'WFCAM_RA_E', 'WFCAM_DEC_E', 'WFCAM_Y', 'WFCAM_J1', 'WFCAM_J2', 'WFCAM_H', 'WFCAM_K', 'WFCAM_Y_E', 'WFCAM_J1_E', 'WFCAM_J2_E', 'WFCAM_H_E', 'WFCAM_K_E', 'WFCAM_PMRA', 'WFCAM_PMDE', 'WFCAM_PMRA_E', 'WFCAM_PMDEC_E', 'WFCAM_NOTES']},
  {"name": "2MASS"    , "table_id": "II/246/out"    , 'table_data': ['RAJ2000', 'DEJ2000', 'Jmag', 'Hmag', 'Kmag', 'Jcmsig', 'Hcmsig', 'Kcmsig']                                                                                                      , 'image_id': ['2MASS-J', '2MASS-H', '2MASS-K']               , 'image_selection': [True, False, True]             , 'image_names': ['J', 'H', 'K']          , 'table_header': ['2MASS_RA', '2MASS_DEC', '2MASS_J', '2MASS_H', '2MASS_K', '2MASS_J_E', '2MASS_H_E', '2MASS_K_E', '2MASS_NOTES']},
  {"name": "PanSTARRS", "table_id": "II/349/ps1"    , 'table_data': ['RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', 'gmag', 'rmag', 'imag', 'zmag', 'ymag', 'e_gmag', 'e_rmag', 'e_imag', 'e_zmag', 'e_ymag']                                        , 'image_id': 'PS2'                                           , 'image_selection': [True, True, False, False]      , 'image_names': ['r', 'i', 'z', 'y']     , 'table_header': ['PS_RA', 'PS_DEC', 'PS_RA_E', 'PS_DEC_E', 'PS_G', 'PS_R', 'PS_I', 'PS_Z', 'PS_Y', 'PS_G_E', 'PS_R_E', 'PS_I_E', 'PS_Z_E', 'PS_Y_E', 'PS_NOTES']}, 
  {"name": "NSC"      , "table_id": "NSC"           , 'table_data': ['ra', 'dec', 'raerr', 'decerr', 'gmag', 'rmag', 'imag', 'zmag', 'umag', 'ymag', 'gerr', 'rerr', 'ierr', 'zerr', 'uerr', 'yerr', 'pmra', 'pmdec', 'pmraerr', 'pmdecerr', 'mjd']   , 'image_id': 'NSC'                                           , 'image_selection': [False, True, True, True, False], 'image_names': ['g', 'r', 'i', 'z', 'Y'], 'table_header': ['NSC_RA', 'NSC_DEC', 'NSC_RA_E', 'NSC_DEC_E', 'NSC_G', 'NSC_R', 'NSC_I', 'NSC_Z', 'NSC_U', 'NSC_Y', 'NSC_G_E', 'NSC_R_E', 'NSC_I_E', 'NSC_Z_E', 'NSC_U_E', 'NSC_Y_E', 'NSC_PMRA', 'NSC_PMDEC', 'NSC_PMRA_E', 'NSC_PMDEC_E', 'NSC_MJD', 'NSC_NOTES']}, 
  {"name": "GALEX"    , "table_id": "II/312/ais"    , 'table_data': ['RAJ2000', 'DEJ2000', 'FUV', 'NUV', 'e_FUV', 'e_NUV']                                                                                                                            , 'image_id': ['GALEX Near UV', 'GALEX Far UV']               , 'image_selection': [True, True]                    , 'image_names': ['NUV', 'FUV']           , 'table_header': ['GALEX_RA', 'GALEX_DEC', 'GALEX_FUV', 'GALEX_NUV', 'GALEX_FUV_E', 'GALEX_NUV_E', 'GALEX_NOTES']}
]

#Makes the drop down window for types of file in the multi-object search
filetype_list = ['CSV', 'FITS', 'ASCII', 'IPAC']

# Gets the file location of the WRAP app
script_path = os.path.abspath(__file__)
directory = (script_path.split('WRAP.py', 2))[0]
# ------------------------------------------------------------- #
  
  
  
# WRAP GUI LAYOUT
# ------------------------------------------------------------- #
if platform != 'win32':
  #Makes the layout of WRAP for the single object search, by providing a location for: ra, dec, radius, output file name, catalogs, and output
  layout_single = [
    [sg.Image(filename = (directory + '/Output/metadata/WRAP_Logo.png'), size = (135, 95)),                         sg.Text('WRAP', justification='center', size=(5, 1), font = ('Chalkduster', 52)),                 sg.Image(filename = directory + '/Output/metadata/BYW_Logo.png', size = (205, 95))],
                   
    [sg.Text('RA', font = ('Times New Roman', 22), size=(13, 1), justification='center'),           sg.Text('DEC', font = ('Times New Roman', 22), size=(13, 1), justification='center'),             sg.Text('RADIUS', font = ('Times New Roman', 22), size=(13, 1), justification='center')],
    [sg.Text('(Degrees)', font = ('Times New Roman', 20), size=(18, 1), justification='center'),    sg.Text('(Degrees)', font = ('Times New Roman', 20), size=(11, 1), justification='center'),       sg.Text('(Arcsecs)', font = ('Times New Roman', 20), size=(20, 1), justification='center')],
    [sg.InputText(size=(18, 1), key = 'RA', font = ('Times New Roman', 15)),                        sg.InputText(size=(18, 2), key = 'DEC', font = ('Times New Roman', 15)),                          sg.InputText(size=(18, 2), key = 'RADIUS', font = ('Times New Roman', 15))],
                    
    [sg.Text('Output File Name', size=(50, 1), justification='center', font = ('Times New Roman', 22))],
    [sg.InputText(key = 'output', font = ('Times New Roman', 15), size = (70, 3), justification='center')],

    [sg.Text('Catalogs:', justification='center', size=(50, 1), font = ('Times New Roman', 25))],   
    [sg.Checkbox('CatWISE 2020', key = 'SINGLE_CatWISE', font = ('Times New Roman', 22), size = (14, 2)),   sg.Checkbox('AllWISE', key = 'SINGLE_AllWISE', font = ('Times New Roman', 22), size = (10, 2)),               sg.Checkbox('Gaia', key = 'SINGLE_Gaia', font = ('Times New Roman', 22), size = (9, 2))],
    [sg.Checkbox('VISTA', key = 'SINGLE_VSA', font = ('Times New Roman', 22), size = (9, 2)),               sg.Checkbox('WFCAM', key = 'SINGLE_WFCAM', font = ('Times New Roman', 22), size = (10, 2)),             sg.Checkbox('2MASS', key = 'SINGLE_2MASS', font = ('Times New Roman', 22), size = (10, 2))],
    [sg.Checkbox('PanSTARRS', key = 'SINGLE_PanSTARRS', font = ('Times New Roman', 22), size = (13, 2)),           sg.Checkbox('NSC', key = 'SINGLE_NSC', font = ('Times New Roman', 22), size = (8, 2)),                   sg.Checkbox('GALEX', key = 'SINGLE_GALEX', font = ('Times New Roman', 22), size = (10, 2))],
    
    [sg.Checkbox('Select All',   enable_events=True, key='Check_All'),                               sg.Checkbox('Deselect All', enable_events=True, key='Uncheck_All'),                               sg.Checkbox('WiseView', enable_events=True, key='wiseview')],

    [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'),                                   sg.Button('Help', size = (17), button_color = '#F7CC7C'),                                         sg.Button('Close WRAP', size = (17), button_color = '#E48671')]
                  ]

  #Makes the layout of WRAP for the multi-object search, by providing a location for: file directory, radius, filetype, output file name, catalogs, and output
  layout_multi = [
    [sg.Image(filename = (directory + '/Output/metadata/WRAP_Logo.png'), size = (135, 95)),                         sg.Text('WRAP', justification='center', size=(5, 1), font = ('Chalkduster', 52)),                 sg.Image(filename = directory + '/Output/metadata/BYW_Logo.png', size = (205, 95))],
                  
    [sg.Text('FILE DIRECTORY', font = ('Times New Roman', 22), size=(50, 1), justification='center')],
    [sg.Text('(CSV, FITS, ASCII, IPAC)', font = ('Times New Roman', 20), size=(50, 1), justification='center')],

    [sg.FileBrowse('File Browser', size = (80, 1), key = 'file', file_types = [('CSV Files', '*.csv'), ('FITS Files', '*.fits'), ('ASCII Files', '*.txt'), ('IPAC Files', '*.txt')])],
    [sg.Text('RADIUS', font = ('Times New Roman', 22), size=(17, 1), justification='center'),              sg.Text('FILETYPE', font = ('Times New Roman', 22), size=(9, 1), justification='center'),               sg.Text('Output File Name', size=(25, 1), justification='center', font = ('Times New Roman', 22))],
    [sg.InputText(size=(22, 2), key = 'RADIUS_multi', font = ('Times New Roman', 15)),                     sg.Combo(filetype_list, size = (13), font = ('Times New Roman', 15), key = 'type'),                     sg.InputText(key = 'output2', font = ('Times New Roman', 15), size = (22, 2), justification='center')],

    [sg.Text('Catalogs:', justification='center', size=(50, 1), font = ('Times New Roman', 25))],
    [sg.Checkbox('CatWISE 2020', key = 'MULTI_CatWISE', font = ('Times New Roman', 22), size = (14, 2)),          sg.Checkbox('AllWISE', key = 'MULTI_AllWISE', font = ('Times New Roman', 22), size = (10, 2)),               sg.Checkbox('Gaia', key = 'MULTI_Gaia', font = ('Times New Roman', 22), size = (9, 2))],
    [sg.Checkbox('VISTA', key = 'MULTI_VSA', font = ('Times New Roman', 22), size = (9, 2)),                      sg.Checkbox('WFCAM', key = 'MULTI_WFCAM', font = ('Times New Roman', 22), size = (10, 2)),                   sg.Checkbox('2MASS', key = 'MULTI_2MASS', font = ('Times New Roman', 22), size = (10, 2))],
    [sg.Checkbox('PanSTARRS', key = 'MULTI_PanSTARRS', font = ('Times New Roman', 22), size = (13, 2)),           sg.Checkbox('NSC', key = 'MULTI_NSC', font = ('Times New Roman', 22), size = (8, 2)),                        sg.Checkbox('GALEX', key = 'MULTI_GALEX', font = ('Times New Roman', 22), size = (10, 2))],
    
    [sg.Checkbox('Select All',   enable_events=True, key='Check_All_Multi'),                               sg.Checkbox('Deselect All', enable_events=True, key='Uncheck_All_Multi'),                               sg.Checkbox('WiseView', enable_events=True, key='wiseview_multi')],

    [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'),                                         sg.Button('Help', size = (17), button_color = '#F7CC7C'),                                               sg.Button('Close WRAP', size = (17), button_color = '#E48671')]
                  ]

  #Makes the general layout for WRAP
  tab_layout = [[sg.TabGroup([[
    sg.Tab('Single Object',   layout_single,       title_color='#F9F8F3',          background_color='#eeeccb',   element_justification= 'center',     key = 'Single Obect Search'),
    sg.Tab('Multi-Object',    layout_multi,        title_color='#F9F8F3',          background_color='#eeeccb',   element_justification= 'center',     key = 'Multi-Object Search')]], 
    tab_location='centertop', title_color='Black', tab_background_color='#F9F8F3', selected_title_color='Black', selected_background_color='#9C873E', border_width = 6, font = ('Times New Roman', 18), enable_events = True, key = 'tab_group'), sg.Button('Close')
              ]] 

  #Generates the window based off the layouts above
  window = sg.Window('WRAP', tab_layout, size = (550, 535), grab_anywhere=False, finalize=True, enable_close_attempted_event = True)
  
if platform == 'win32':  
  #Makes the layout of WRAP for the single object search, by providing a location for: ra, dec, radius, output file name, catalogs, and output
  layout_single = [
    [sg.Image(filename = ( directory + '/Output/metadata/WRAP_Logo.png'), size = (135, 95)),                         sg.Text('WRAP', justification='center', size=(5, 1), font = ('Chalkduster', 52)),                 sg.Image(filename = directory + '/Output/metadata/BYW_Logo.png', size = (205, 95))],
                   
    [sg.Text('RA', font = ('Times New Roman', 22), size=(13, 1), justification='center'),           sg.Text('DEC', font = ('Times New Roman', 22), size=(13, 1), justification='center'),             sg.Text('RADIUS', font = ('Times New Roman', 22), size=(13, 1), justification='center')],
    [sg.Text('(Degrees)', font = ('Times New Roman', 20), size=(18, 1), justification='center'),    sg.Text('(Degrees)', font = ('Times New Roman', 20), size=(11, 1), justification='center'),       sg.Text('(Arcsecs)', font = ('Times New Roman', 20), size=(20, 1), justification='center')],
    [sg.InputText(size=(18, 1), key = 'RA', font = ('Times New Roman', 15)),                        sg.InputText(size=(18, 2), key = 'DEC', font = ('Times New Roman', 15)),                          sg.InputText(size=(18, 2), key = 'RADIUS', font = ('Times New Roman', 15))],
                    
    [sg.Text('Output File Name', size=(50, 1), justification='center', font = ('Times New Roman', 22))],
    [sg.InputText(key = 'output', font = ('Times New Roman', 15), size = (70, 3), justification='center')],

    [sg.Text('Catalogs:', justification='center', size=(50, 1), font = ('Times New Roman', 25))],   
    [sg.Checkbox('CatWISE 2020', key = 'SINGLE_CatWISE', font = ('Times New Roman', 22), size = (14, 2)),   sg.Checkbox('AllWISE', key = 'SINGLE_AllWISE', font = ('Times New Roman', 22), size = (10, 2)),               sg.Checkbox('Gaia', key = 'SINGLE_Gaia', font = ('Times New Roman', 22), size = (9, 2))],
    [sg.Checkbox('VISTA', key = 'SINGLE_VSA', font = ('Times New Roman', 22), size = (9, 2)),               sg.Checkbox('WFCAM', key = 'SINGLE_WFCAM', font = ('Times New Roman', 22), size = (10, 2)),             sg.Checkbox('2MASS', key = 'SINGLE_2MASS', font = ('Times New Roman', 22), size = (10, 2))],
    [sg.Checkbox('PanSTARRS', key = 'SINGLE_PanSTARRS', font = ('Times New Roman', 22), size = (13, 2)),    sg.Checkbox('GALEX', key = 'SINGLE_GALEX', font = ('Times New Roman', 22), size = (10, 2))],
    
    [sg.Checkbox('Select All',   enable_events=True, key='Check_All'),                               sg.Checkbox('Deselect All', enable_events=True, key='Uncheck_All'),                               sg.Checkbox('WiseView', enable_events=True, key='wiseview')],

    [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'),                                   sg.Button('Help', size = (17), button_color = '#F7CC7C'),                                         sg.Button('Close WRAP', size = (17), button_color = '#E48671')]
                  ]

  #Makes the layout of WRAP for the multi-object search, by providing a location for: file directory, radius, filetype, output file name, catalogs, and output
  layout_multi = [
    [sg.Image(filename = (directory + '/Output/metadata/WRAP_Logo.png'), size = (135, 95)),                         sg.Text('WRAP', justification='center', size=(5, 1), font = ('Chalkduster', 52)),                 sg.Image(filename = directory + '/Output/metadata/BYW_Logo.png', size = (205, 95))],
                  
    [sg.Text('FILE DIRECTORY', font = ('Times New Roman', 22), size=(50, 1), justification='center')],
    [sg.Text('(CSV, FITS, ASCII, IPAC)', font = ('Times New Roman', 20), size=(50, 1), justification='center')],

    [sg.FileBrowse('File Browser', size = (80, 1), key = 'file', file_types = [('CSV Files', '*.csv'), ('FITS Files', '*.fits'), ('ASCII Files', '*.txt'), ('IPAC Files', '*.txt')])],
    [sg.Text('RADIUS', font = ('Times New Roman', 22), size=(17, 1), justification='center'),              sg.Text('FILETYPE', font = ('Times New Roman', 22), size=(9, 1), justification='center'),               sg.Text('Output File Name', size=(25, 1), justification='center', font = ('Times New Roman', 22))],
    [sg.InputText(size=(22, 2), key = 'RADIUS_multi', font = ('Times New Roman', 15)),                     sg.Combo(filetype_list, size = (13), font = ('Times New Roman', 15), key = 'type'),                     sg.InputText(key = 'output2', font = ('Times New Roman', 15), size = (22, 2), justification='center')],

    [sg.Text('Catalogs:', justification='center', size=(50, 1), font = ('Times New Roman', 25))],
    [sg.Checkbox('CatWISE 2020', key = 'MULTI_CatWISE', font = ('Times New Roman', 22), size = (14, 2)),          sg.Checkbox('AllWISE', key = 'MULTI_AllWISE', font = ('Times New Roman', 22), size = (10, 2)),               sg.Checkbox('Gaia', key = 'MULTI_Gaia', font = ('Times New Roman', 22), size = (9, 2))],
    [sg.Checkbox('VISTA', key = 'MULTI_VSA', font = ('Times New Roman', 22), size = (9, 2)),                      sg.Checkbox('WFCAM', key = 'MULTI_WFCAM', font = ('Times New Roman', 22), size = (10, 2)),                   sg.Checkbox('2MASS', key = 'MULTI_2MASS', font = ('Times New Roman', 22), size = (10, 2))],
    [sg.Checkbox('PanSTARRS', key = 'MULTI_PanSTARRS', font = ('Times New Roman', 22), size = (13, 2)),           sg.Checkbox('GALEX', key = 'MULTI_GALEX', font = ('Times New Roman', 22), size = (10, 2))],
    
    [sg.Checkbox('Select All',   enable_events=True, key='Check_All_Multi'),                               sg.Checkbox('Deselect All', enable_events=True, key='Uncheck_All_Multi'),                               sg.Checkbox('WiseView', enable_events=True, key='wiseview_multi')],

    [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'),                                         sg.Button('Help', size = (17), button_color = '#F7CC7C'),                                               sg.Button('Close WRAP', size = (17), button_color = '#E48671')]
                  ]

  #Makes the general layout for WRAP
  tab_layout = [[sg.TabGroup([[
    sg.Tab('Single Object',   layout_single,       title_color='#F9F8F3',          background_color='#eeeccb',   element_justification= 'center',     key = 'Single Obect Search'),
    sg.Tab('Multi-Object',    layout_multi,        title_color='#F9F8F3',          background_color='#eeeccb',   element_justification= 'center',     key = 'Multi-Object Search')]], 
    tab_location='centertop', title_color='Black', tab_background_color='#F9F8F3', selected_title_color='Black', selected_background_color='#9C873E', border_width = 6, font = ('Times New Roman', 18), enable_events = True, key = 'tab_group'), sg.Button('Close')
              ]] 

  #Generates the window based off the layouts above
  window = sg.Window('WRAP', tab_layout, size = (550, 535), grab_anywhere=False, finalize=True, enable_close_attempted_event = True)
# ------------------------------------------------------------- #

print('#------------------------------------------------#')
print('               Opening WRAP                       ')
print('#------------------------------------------------#')

# RUNNING WRAP GUI
# ------------------------------------------------------------- #
#Keeps the window open
while True:

  #Reads all of the events and values, then reads which tab is currently in
  event, values = window.read()
  group = values['tab_group']

  #Runs the program if the user is in 'Single Object' tab
  if group == 'Single Obect Search':  

    #Adds the functionality of the 'Select All' and 'Deselect All' Buttons
    if event == 'Check_All':
      window['SINGLE_CatWISE'].update(True),   window['SINGLE_AllWISE'].update(True), window['SINGLE_Gaia'].update(True)
      window['SINGLE_VSA'].update(True),       window['SINGLE_WFCAM'].update(True),   window['SINGLE_2MASS'].update(True)
      window['SINGLE_PanSTARRS'].update(True), window['SINGLE_GALEX'].update(True)
      window['Uncheck_All'].update(False)
      if platform != 'win32':
        window['SINGLE_NSC'].update(True)
    if event == 'Uncheck_All':
      window['SINGLE_CatWISE'].update(False),   window['SINGLE_AllWISE'].update(False), window['SINGLE_Gaia'].update(False)
      window['SINGLE_VSA'].update(False),       window['SINGLE_WFCAM'].update(False),   window['SINGLE_2MASS'].update(False)
      window['SINGLE_PanSTARRS'].update(False), window['SINGLE_GALEX'].update(False)
      window['Check_All'].update(False)
      if platform != 'win32':
        window['SINGLE_NSC'].update(False)

    #Checks if the 'Run WRAP' button is clicked
    if event in (None, 'Run WRAP'):

      #Calls the "single_tab_check" function
      fake_list = single_tab_check()
        
      #If the RA, DEC, and RADIUS tabs are filled then runs the search on the catalogs used
      if len(fake_list) == 0: 
        print('')
        print('#------------------------------------------------#')
        print('RA (deg): ' + str(values['RA']))
        print('DEC (deg): ' + str(values['DEC']))
        print('RADIUS (arcsec): ' + str(values['RADIUS']))
        print('#------------------------------------------------#')
        print('')

        #Make variables for the RA, DEC, and RADIUS used
        ra_use, dec_use, radius_use = float(values['RA']), float(values['DEC']), int(values['RADIUS'])
        
        #Creates fake list for the data and data names
        photometry, photometry_name = [], []

        #Calls the "single_object_search" function
        single_object_search()

    #Provides the user with the authors information if the 'Help' button is pressed
    if event in (None, 'Help'):
      print('#------------------------------------------------#')
      print('#     Thank you for using WRAP!                  #')
      print('#   Authors Contact: hcb98@nau.edu               #')
      print('#------------------------------------------------#')

    #Closes WRAP if the 'Close WRAP' button is pressed
    if event in (None, 'Close WRAP'):
      print('#------------------------------------------------#')
      print('               Closing WRAP                       ')
      print('#------------------------------------------------#')
      break
    
  #Runs the program if the user is in 'Multi-Object' tab
  if group == 'Multi-Object Search':

    #Adds the functionality of the 'Select All' and 'Deselect All' Buttons
    if event == 'Check_All_Multi':
      window['MULTI_CatWISE'].update(True),   window['MULTI_AllWISE'].update(True), window['MULTI_Gaia'].update(True)
      window['MULTI_VSA'].update(True),       window['MULTI_WFCAM'].update(True),   window['MULTI_2MASS'].update(True)
      window['MULTI_PanSTARRS'].update(True), window['MULTI_GALEX'].update(True)
      window['Uncheck_All_Multi'].update(False)
      if platform != 'win32':
        window['MULTI_NSC'].update(True)
    if event == 'Uncheck_All_Multi':
      window['MULTI_CatWISE'].update(False),   window['MULTI_AllWISE'].update(False), window['MULTI_Gaia'].update(False)
      window['MULTI_VSA'].update(False),       window['MULTI_WFCAM'].update(False),   window['MULTI_2MASS'].update(False)
      window['MULTI_PanSTARRS'].update(False), window['MULTI_GALEX'].update(False)
      window['Check_All_Multi'].update(False)
      if platform != 'win32':
        window['MULTI_NSC'].update(False)        

    #Checks if the 'Run WRAP' button is clicked
    if event in (None, 'Run WRAP0'):

      #Calls the "multi_tab_check" function
      multi_tab_check()

      #If the RA, DEC, and RADIUS tabs are filled then runs the search on the catalogs used
      if values['file'] != '' and values['RADIUS_multi'].isnumeric() == True and (values['type'] == 'CSV' or values['type'] == 'FITS' or values['type'] == 'ASCII' or values['type'] == 'IPAC'):
        print('')
        print('#------------------------------------------------#')
        print('Directory: ' + str(values['file']))
        print('File Type: ' + str(values['type']))
        print('RADIUS (arcsec): ' + str(values['RADIUS_multi']))
        print('#------------------------------------------------#')
        print('')

        #Calls the "multi_tab_table" function
        ra_list, dec_list = multi_tab_table()

        #Makes a csv file and writes the header
        output = values['output2']
        
        #Makes the output file name
        if values['output2'] == '':
          output = 'WRAP_output'
        else: 
          output = values['output2']

        #Calls the "multi_object_search" function
        multi_object_search()

    #Provides the user with the authors information if the 'Help' button is pressed
    if event in (None, 'Help1'):
      print('#------------------------------------------------#')
      print('#     Thank you for using WRAP!                  #')
      print('#   Authors Contact: hcb98@nau.edu               #')
      print('#------------------------------------------------#')

    #Closes WRAP if the 'Close WRAP' button is pressed
    if event in (None, 'Close WRAP2'):
      print('#------------------------------------------------#')
      print('               Closing WRAP                       ')
      print('#------------------------------------------------#')
      break

#Closes the window
window.close()
# ------------------------------------------------------------- #