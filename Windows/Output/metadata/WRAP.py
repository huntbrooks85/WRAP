#-----------------------------------------------------------------------#
# WRAP v1.0.0
# By Hunter Brooks, at NAU, Flagstaff: December 1, 2023
#
# Purpose: Gathers photometry and astrometry from various 
#          ultra-violet, optical, and near-infrared catalogs 
#          given a RA, DEC, and Radius by the user
#-----------------------------------------------------------------------#

# Import all needed packages.
# ------------------------------------------------------------- #
import cv2
import csv
import math
import time
import sys,os
import astropy
import requests
import truncate
import webbrowser
import matplotlib
import numpy as np
import pandas as pd
from astropy import wcs
from sys import platform
from pyvo.dal import sia
import PySimpleGUI as sg
import matplotlib.cm as cm
from astropy.wcs import WCS
from astropy.io import fits
from astropy.io import ascii
from astropy.time import Time
from matplotlib import pyplot
from bs4 import BeautifulSoup
from astroquery.vsa import Vsa
from astropy import units as u
from astropy.table import Table
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia
from matplotlib import transforms
from astropy.utils.data import conf
from astropy.nddata import Cutout2D
import astropy.coordinates as coord
from astroquery.ukidss import Ukidss
from matplotlib.widgets import Cursor
from astroquery.ipac.irsa import Irsa
from astropy.coordinates import SkyCoord
from matplotlib.transforms import Affine2D
from matplotlib.transforms import Affine2D
from astropy.utils.data import download_file
from astropy.visualization.wcsaxes import WCSAxes
from astropy.visualization import astropy_mpl_style
from astroquery.mast import Observations, MastMissions, Catalogs
from matplotlib.widgets import Slider, Button, CheckButtons, TextBox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, FigureManager
from astropy.visualization import (PercentileInterval, SinhStretch, ImageNormalize)
# ------------------------------------------------------------- #

# Gets the file location of the WRAP app
# ------------------------------------------------------------- #
executable_path = sys.executable
directory = os.path.dirname(executable_path)
directory = directory.split('WRAP.app', 2)[0]
# ------------------------------------------------------------- #

#Does not load these packages if the user is on a Windows machine
# ------------------------------------------------------------- #
if platform != 'win32':
    from dl.helpers.utils import convert
    from dl import authClient as ac, queryClient as qc

#Ignores all of the warnings from the packages above
# ------------------------------------------------------------- #
import warnings
warnings.filterwarnings("ignore")

#                      AllWISE SEARCH                           #
# ------------------------------------------------------------- #
def blockPrint():
  '''Makes a function that blocks the printing function'''
  sys.stdout = open(os.devnull, 'w')

def enablePrint():
  '''Makes a function that allows the printing function'''
  sys.stdout = sys.__stdout__

def allwise_image(ra, dec, radius): 
  ''' First, it gets the images from the AllWISE API from IRSA and downloads the images. 
  Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
  Third, calls the table function to get all the objects from the AllWISE source catalog. 
  Fourth, makes the window for the user to click the object with all settings. 
  Finally, finds the closest object to the click and records the data. '''

  #Makes outline for the window of the plot
  plt.rcParams['toolbar'] = 'None'
  matplotlib.use("TkAgg")
  plt.style.use('Solarize_Light2')
  blockPrint()

  #Finds all the metadata that relates to the ra and dec searched, mostly to find the APIs for the W1, W2, W3, and W4 images
  metadata_allwise_link = 'http://irsa.ipac.caltech.edu/ibe/sia/wise/allwise/p3am_cdd?POS=' + str(ra) + ',' + str(dec) + '&SIZE=' + str(radius/3600)
  try:
    allwise_metadata = requests.get(metadata_allwise_link)
  except: 
    ra_aw_e, dec_aw_e, w1, w1_sigma, w2, w2_sigma, w3, w3_sigma, w4, w4_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    ra_allwise = ra
    dec_allwise = dec
    plt.close('all')
    plt.figure().clear()
    return ra_allwise, ra_aw_e, dec_allwise, dec_aw_e, w1, w1_sigma, w2, w2_sigma, w3, w3_sigma, w4, w4_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma, 'AllWISE Source Catalog', 'AllWISE is Having Maintenance'
  open(str(directory) + '/Output/metadata/AllWISE_metadata.txt', 'wb').write(allwise_metadata.content)

  #With this metadata it finds the API link for the W1 and W2 images
  w1_finder, w2_finder, w3_finder, w4_finder = 'W1 Coadd', 'W2 Coadd', 'W3 Coadd', 'W4 Coadd'
  with open(str(directory) + '/Output/metadata/AllWISE_metadata.txt', 'r') as fp:
   lines = fp.readlines()
   for line in lines:
      if line.find(w1_finder) != -1:
        w1_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]
      elif line.find(w2_finder) != -1:
        w2_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]
      elif line.find(w3_finder) != -1:
        w3_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]
      elif line.find(w4_finder) != -1:
        w4_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]

  #Download the W1 and W2 images
  file_allwise_w1, file_allwise_w2 = download_file(w1_allwise_image_url, cache=True), download_file(w2_allwise_image_url, cache=True)
  file_allwise_w3, file_allwise_w4 = download_file(w3_allwise_image_url, cache=True), download_file(w4_allwise_image_url, cache=True) 

  data_allwise_w1, data_allwise_w2 =  fits.getdata(file_allwise_w1), fits.getdata(file_allwise_w2)
  data_allwise_w3, data_allwise_w4 = fits.getdata(file_allwise_w3), fits.getdata(file_allwise_w4)

  #Find the location of all the object found in AllWISE in the radius choosen by the user 
  location_data = allwise_table(ra, dec, radius)
  object_ra, object_dec = location_data['ra'].tolist(), location_data['dec'].tolist()
  object_ra_sigma, object_dec_sigma = location_data['sigra'].tolist(), location_data['sigdec'].tolist()
  w1_list, w2_list = location_data['w1mpro'].tolist(), location_data['w2mpro'].tolist()
  w3_list, w4_list = location_data['w3mpro'].tolist(), location_data['w4mpro'].tolist()
  w1_list_sigma, w2_list_sigma = location_data['w1sigmpro'].tolist(), location_data['w2sigmpro'].tolist()
  w3_list_sigma, w4_list_sigma = location_data['w3sigmpro'].tolist(), location_data['w4sigmpro'].tolist()
  pmra_list, pmdec_list = location_data['pmra'].tolist(), location_data['pmdec'].tolist()
  pmra_sigma_list, pmdec_sigma_list = location_data['sigpmra'].tolist(), location_data['sigpmdec'].tolist()

  #Obtains the headers for each image
  hdu_w1, hdu_w2 = fits.open(file_allwise_w1)[0], fits.open(file_allwise_w2)[0]
  hdu_w3, hdu_w4 = fits.open(file_allwise_w3)[0], fits.open(file_allwise_w4)[0]
  wcs = WCS(hdu_w1.header)

  #Make a cutout from the coadd image for the RA and DEC put in
  position, size = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0'), u.Quantity([radius, radius], u.arcsec)
  cutout_w1 = Cutout2D(data_allwise_w1, position, size, fill_value = np.nan, wcs = wcs.celestial)
  cutout_w2 = Cutout2D(data_allwise_w2, position, size, fill_value = np.nan, wcs = wcs.celestial)
  cutout_w3 = Cutout2D(data_allwise_w3, position, size, fill_value = np.nan, wcs = wcs.celestial)
  cutout_w4 = Cutout2D(data_allwise_w4, position, size, fill_value = np.nan, wcs = wcs.celestial)
  wcs_cropped = cutout_w1.wcs
  enablePrint()

  #Obtains the dates for each image
  date_w1, date_w2 = hdu_w1.header['MIDOBS'].split('T', 2)[0], hdu_w2.header['MIDOBS'].split('T', 2)[0]
  date_w3, date_w4 = hdu_w3.header['MIDOBS'].split('T', 2)[0], hdu_w4.header['MIDOBS'].split('T', 2)[0]

  #Defining a mouse click as an event on the plot
  location = []
  plt.rcParams["figure.figsize"] = [8, 8]
  plt.rcParams["figure.autolayout"] = True
  def mouse_event(event):
    '''Makes a list of the x, y, and axes the mouse click is.'''
    
    location.append(event.ydata)
    location.append(event.xdata)
    location.append(event.inaxes)
  plt.connect('button_press_event', mouse_event)

  #Sets the WCS coordinates for the plots
  total_data = cutout_w1.data + cutout_w2.data
  ax = plt.subplot(projection = wcs_cropped)

  #Plots the objects found in the radius
  circle_size = (radius*3)
  scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

  #Normalize the image and plots it
  init_top, init_bot = 95, 45
  norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, init_bot), vmax = np.nanpercentile(total_data.data, init_top))
  ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

  #Makes the figure look pretty
  plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
  fontdict = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
  plt.suptitle('AllWISE Search', fontsize = 35, y = 0.96, fontfamily = 'Times New Roman')
  ax.set_title('Dates: \n'
             + 'W1 Date: ' + str(date_w1) + ' (Y-M-D)  ' + '  W2 Date: ' + str(date_w2) + ' (Y-M-D)\n'
             + 'W3 Date: ' + str(date_w3) + ' (Y-M-D)  ' + '  W4 Date: ' + str(date_w4) + ' (Y-M-D)\n', fontdict = fontdict, y = 1.05)
  plt.grid(linewidth = 0)
  figure = plt.gcf()
  plt.xlim(0, max(total_data.shape))
  plt.ylim(0, max(total_data.shape))
  if platform != 'win32':
    figure.set_size_inches(4.75, 7.25)
  elif platform == 'win32':
    figure.set_size_inches(4.75, 7.25)
  # figure.canvas.set_window_title('AllWISE Search')
  mng = pyplot.get_current_fig_manager()
  mng.window.resizable(False, False)
  
  #Makes a cursor to aid the User
  cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
  annotation = ax.annotate('', xy=(0, 0), xytext=(20, -20), arrowprops=dict(arrowstyle='wedge'), fontsize=12, color='red')
  annotation.set_visible(False)

  #Make checkbuttons with all of the different image bands
  rax = plt.axes([0.045, 0.4, 0.115, 0.1])
  labels = ['W1', 'W2', 'W3', 'W4']
  real_data = [cutout_w1.data, cutout_w2.data, cutout_w3.data, cutout_w4.data]
  default = [True, True, False, False]
  check = CheckButtons(rax, labels, default)

  #Adds a slider for the scaling of the image
  freq_top = plt.axes([0.25, 0.12, 0.65, 0.03])
  slider_top = Slider(ax = freq_top, label = 'Top Stetch:', valmin = 50, valmax = 100, valinit = init_top, color = '#E48671')
  freq_bottom = plt.axes([0.25, 0.087, 0.65, 0.03])
  slider_bottom = Slider(ax = freq_bottom, label = 'Bottom Stetch:', valmin = 0, valmax = 50, valinit = init_bot, color = '#E48671')

  #Adds a slider for the circle size
  circle_slid_location = plt.axes([0.25, 0.055, 0.65, 0.03])
  circle_slider = Slider(ax = circle_slid_location, label = 'Circle Size:', valmin = (circle_size - 2.5*radius), valmax = (circle_size + 1*radius), valinit = circle_size, color = '#E48671')

  #Adds a notes section that the user can add notes about their data
  axbox = plt.axes([0.15, 0.02, 0.8, 0.03])
  text = ''
  text_box = TextBox(axbox, 'Notes:', initial = text, textalignment="center")

  #Make a button that can be clicked if no object is found
  axes_button = plt.axes([0.04, 0.775, 0.92, 0.04])
  close = Button(axes_button, 'Object Not Found', color = '#E48671')

  #Updates the image depending on what the user chooses
  def update_button(label):
    '''Updates the list of activated images and updates the image the user can see.'''

    total_data = 0
    for lab in labels:
      if lab == label:
        index = labels.index(lab)
        if default[index] == False:
          default[index] = True
        elif default[index] == True: 
          default[index] = False
    for d in range(len(default)):
      if default == [False, False, False, False]: 
        total_data = real_data[0]*0
      if default[d] == True: 
        total_data = total_data + real_data[d]
      else: 
        pass
    norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
    ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

  #Updates the scaling when the slider is changed
  def update_slider_stretch(val):
    '''Updates the stretch the user can see, based in percentiles'''

    total_data = 0
    for d in range(len(default)):
      if default[d] == True: 
        total_data = total_data + real_data[d]
      else: 
        pass
    norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
    ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

  #Updates the notes added by the user when there is an input
  text_list = [text]
  def submit(expression):
    '''Updates the list of types in the 'Notes' setting'''

    text = expression
    text_list.append(text)

  #Allows the sliders and buttons to be pressed
  check.on_clicked(update_button)
  slider_top.on_changed(update_slider_stretch)
  slider_bottom.on_changed(update_slider_stretch)
  text_box.on_text_change(submit)

  #Display image until it is clicked to find the object
  n = -1
  while True:
    press = plt.waitforbuttonpress()
    text_max = len(text_list) - 1

    #Checks that it was a mouse click
    if press == False:
      n += 3

      #Finds which axes was clicked
      click_axes = str(location[n])
      click_axes = click_axes.split('WCSAxesSubplot', 2)[0]

      #Checks if the image was clicked
      if click_axes == '':
          
        #Makes a pop-up window with success text
        plt.clf()
        plt.close('all')
        plt.figure(1)
        plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for AllWISE! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.grid(linewidth = 0)
        ax = plt.gca()
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.yaxis.set_tick_params(labelleft=False)
        ax.set_xticks([])
        ax.set_yticks([])
        figure2 = plt.gcf()
        figure2.set_size_inches(4.75, 1)
        # figure2.canvas.set_window_title('Successful AllWISE Search')
        mng2 = pyplot.get_current_fig_manager()
        mng2.window.resizable(False, False)
        plt.pause(0.1)
        plt.clf()
        plt.close('all')

        #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
        coord = wcs_cropped.pixel_to_world_values(location[n-4],location[n-5])
        distance = []
        for i in range(len(object_ra)):
          distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))
        list_location = distance.index(np.min(distance))

        #Obtains all of the clicked objects data
        ra_allwise, dec_allwise = object_ra[list_location], object_dec[list_location]
        ra_aw_e, dec_aw_e = object_ra_sigma[list_location], object_dec_sigma[list_location]
        w1, w2 = w1_list[list_location], w2_list[list_location]
        w3, w4 = w3_list[list_location], w4_list[list_location]
        w1_sigma, w2_sigma = w1_list_sigma[list_location], w2_list_sigma[list_location]
        w3_sigma, w4_sigma = w3_list_sigma[list_location], w4_list_sigma[list_location]
        pmra, pmdec = pmra_list[list_location], pmdec_list[list_location]
        pmra_sigma, pmdec_sigma = pmra_sigma_list[list_location], pmdec_sigma_list[list_location]
        return ra_allwise, ra_aw_e, dec_allwise, dec_aw_e, w1, w1_sigma, w2, w2_sigma, w3, w3_sigma, w4, w4_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma, 'AllWISE Source Catalog', text_list[text_max]
      
      #Checks if the "Object Not Found" button was clicked
      elif click_axes == 'Axes(0.04,0.775;0.92x0.04)':
        ra_aw_e, dec_aw_e, w1, w1_sigma, w2, w2_sigma, w3, w3_sigma, w4, w4_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        ra_allwise = ra
        dec_allwise = dec

        #Makes a pop-up window with success text
        plt.clf()
        plt.close('all')
        plt.figure(1)
        plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for AllWISE! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.grid(linewidth = 0)
        ax = plt.gca()
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.yaxis.set_tick_params(labelleft=False)
        ax.set_xticks([])
        ax.set_yticks([])
        figure2 = plt.gcf()
        figure2.set_size_inches(4.75, 1)
        # figure2.canvas.set_window_title('Successful AllWISE Search')
        mng2 = pyplot.get_current_fig_manager()
        mng2.window.resizable(False, False)
        plt.pause(0.1)
        plt.clf()
        plt.close('all')

        return ra_allwise, ra_aw_e, dec_allwise, dec_aw_e, w1, w1_sigma, w2, w2_sigma, w3, w3_sigma, w4, w4_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma, 'AllWISE Source Catalog', 'Object Not Found was Pressed'
      
      #Updates the circle size when slider is moved
      elif click_axes == 'Axes(0.25,0.055;0.65x0.03)':
        scatter.remove()
        scatter = ax.scatter(object_ra, object_dec, transform = ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')

    #Checks if the window was closed
    elif press is None:
      ra_aw_e, dec_aw_e, w1, w1_sigma, w2, w2_sigma, w3, w3_sigma, w4, w4_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
      ra_allwise = ra
      dec_allwise = dec
      plt.close('all')
      plt.figure().clear()
      return ra_allwise, ra_aw_e, dec_allwise, dec_aw_e, w1, w1_sigma, w2, w2_sigma, w3, w3_sigma, w4, w4_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma, 'AllWISE Source Catalog', text_list[text_max]
  
def allwise_table(ra, dec, radius): 
  '''Find all the objects in the radius defined by the user'''

  blockPrint()
 
  #Uses astroquery to find all objects in the radius
  location_data = Irsa.query_region(
    coord.SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='fk5'), catalog='allwise_p3as_psd', spatial='Cone', radius=(radius - 1) * u.arcsec)
  return location_data

#                      CATWISE SEARCH                           #
# ------------------------------------------------------------- #
def catwise_image(ra, dec, radius): 
  ''' First, it gets the images from the CatWISE 2020 API from IRSA and downloads the images. 
  Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
  Third, calls the table function to get all the objects from the CatWISE 2020 source catalog. 
  Fourth, makes the window for the user to click the object with all settings. 
  Finally, finds the closest object to the click and records the data. '''

  #Makes outline for the window of the plot
  plt.rcParams['toolbar'] = 'None'
  matplotlib.use("TkAgg")
  plt.style.use('Solarize_Light2')
  blockPrint()

  #Finds all the metadata that relates to the ra and dec searched, mostly to find the APIs for the W1, W2, W3, and W4 images
  metadata_allwise_link = 'http://irsa.ipac.caltech.edu/ibe/sia/wise/allwise/p3am_cdd?POS=' + str(ra) + ',' + str(dec) + '&SIZE=' + str(radius/3600)
  try:
    allwise_metadata = requests.get(metadata_allwise_link)
  except: 
    mjd, ra_cw_e, dec_cw_e, w1, w1_sigma, w2, w2_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    ra_catwise = ra
    dec_catwise = dec
    shape_x, shape_y = total_data.shape[0], total_data.shape[1]
    plt.close('all')
    plt.figure().clear()
    return ra_catwise, ra_cw_e, dec_catwise, dec_cw_e, w1, w1_sigma, w2, w2_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma, mjd, 'CatWISE 2020 Catalog', 'CatWISE is Having Maintenance'
  open(str(directory) + '/Output/metadata/catwise_metadata.txt', 'wb').write(allwise_metadata.content)

  #With this metadata it finds the API link for the W1 and W2 images
  w1_finder, w2_finder = 'W1 Coadd', 'W2 Coadd'
  with open(str(directory) + '/Output/metadata/catwise_metadata.txt', 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      if line.find(w1_finder) != -1:
        w1_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]
      elif line.find(w2_finder) != -1:
        w2_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]

  #Download the W1 and W2 images
  file_allwise_w1, file_allwise_w2 = download_file(w1_allwise_image_url, cache=True), download_file(w2_allwise_image_url, cache=True)
  data_allwise_w1, data_allwise_w2 = fits.getdata(file_allwise_w1), fits.getdata(file_allwise_w2)

  #Find the location of all the object found in CatWISE in the radius choosen by the user 
  location_data = catwise_table(ra, dec, radius)
  object_mjd = location_data['meanobsmjd'].tolist()
  object_ra, object_dec = location_data['ra'].tolist(), location_data['dec'].tolist()
  object_ra_sigma, object_dec_sigma = location_data['sigra'].tolist(), location_data['sigdec'].tolist()
  w1_list, w2_list = location_data['w1mpro'].tolist(), location_data['w2mpro'].tolist()
  w1_list_sigma, w2_list_sigma = location_data['w1sigmpro'].tolist(), location_data['w2sigmpro'].tolist()
  pmra_list, pmdec_list = location_data['pmra'].tolist(), location_data['pmdec'].tolist()
  pmra_sigma_list, pmdec_sigma_list = location_data['sigpmra'].tolist(), location_data['sigpmdec'].tolist()

  #Make a cutout from the coadd image for the RA and DEC put in
  hdu_w1, hdu_w2 = fits.open(file_allwise_w1)[0], fits.open(file_allwise_w2)[0]
  wcs = WCS(hdu_w1.header)
  position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
  size = u.Quantity([radius, radius], u.arcsec)
  cutout_w1 = Cutout2D(data_allwise_w1, position, size, fill_value = np.nan, wcs = wcs.celestial)
  cutout_w2 = Cutout2D(data_allwise_w2, position, size, fill_value = np.nan, wcs = wcs.celestial)
  wcs_cropped = cutout_w1.wcs
  enablePrint()

  #Obtains the dates for each image
  date_w1, date_w2 = hdu_w1.header['MIDOBS'].split('T', 2)[0], hdu_w2.header['MIDOBS'].split('T', 2)[0]

  #Defining a mouse click as an event on the plot
  location = []
  plt.rcParams["figure.figsize"] = [8, 8]
  plt.rcParams["figure.autolayout"] = True
  def mouse_event(event):
    '''Makes a list of the x, y, and axes the mouse click is.'''

    location.append(event.ydata)
    location.append(event.xdata)
    location.append(event.inaxes)
  plt.connect('button_press_event', mouse_event)

  #Sets the WCS coordinates for the plots
  total_data = cutout_w1.data + cutout_w2.data
  ax = plt.subplot(projection = wcs_cropped)

  #Plots the objects found in the radius
  circle_size = (radius*3)
  scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

  #Normalize the image and plots it
  init_top = 95
  init_bot = 45
  norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, init_bot), vmax = np.nanpercentile(total_data.data, init_top))
  ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

  #Makes the figure look pretty
  plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
  fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
  plt.suptitle('CatWISE Search', fontsize = 35, y = 0.96, fontfamily = 'Times New Roman')
  ax.set_title('Dates: \n'
             + 'W1 Date: ' + str(date_w1) + ' (Y-M-D)  \n' 
             + 'W2 Date: ' + str(date_w2) + ' (Y-M-D)  \n', fontdict = fontdict_1, y = 1.05)
  plt.grid(linewidth = 0)
  figure = plt.gcf()
  plt.xlim(0, max(total_data.shape))
  plt.ylim(0, max(total_data.shape))
  if platform != 'win32':
    figure.set_size_inches(4.75, 7.25)
  elif platform == 'win32':
    figure.set_size_inches(4.75, 7.25)
  # figure.canvas.set_window_title('CatWISE Search')
  mng = pyplot.get_current_fig_manager()
  mng.window.resizable(False, False)
  
  #Makes a cursor to aid the User
  cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
  annotation = ax.annotate('', xy=(0, 0), xytext=(20, -20), arrowprops=dict(arrowstyle='wedge'), fontsize=12, color='red')
  annotation.set_visible(False)

  #Make checkbuttons with all of the different image bands
  rax = plt.axes([0.045, 0.4, 0.115, 0.08])
  labels = ['W1', 'W2']
  real_data = [cutout_w1.data, cutout_w2.data]
  default = [True, True]
  check = CheckButtons(rax, labels, default)

  #Adds a slider for the scaling of the image
  freq_top = plt.axes([0.25, 0.12, 0.65, 0.03])
  slider_top = Slider(ax = freq_top, label = 'Top Stetch:', valmin = 50, valmax = 100, valinit = init_top, color = '#E48671')
  freq_bottom = plt.axes([0.25, 0.087, 0.65, 0.03])
  slider_bottom = Slider(ax = freq_bottom, label = 'Bottom Stetch:', valmin = 0, valmax = 50, valinit = init_bot, color = '#E48671')

  #Adds a slider for the circle size
  circle_slid_location = plt.axes([0.25, 0.055, 0.65, 0.03])
  circle_slider = Slider(ax = circle_slid_location, label = 'Circle Size:', valmin = (circle_size - 2.5*radius), valmax = (circle_size + 1*radius), valinit = circle_size, color = '#E48671')

  #Adds a notes section that the user can add notes about their data
  axbox = plt.axes([0.15, 0.02, 0.8, 0.03])
  text = ''
  text_box = TextBox(axbox, 'Notes:', initial = text, textalignment="center")

  #Make a button that can be clicked if no object is found
  axes_button = plt.axes([0.04, 0.775, 0.92, 0.04])
  close = Button(axes_button, 'Object Not Found', color = '#E48671')

  #Update the image depending on what the user chooses
  def update_button(label):
    '''Updates the list of activated images and updates the image the user can see.'''

    total_data = 0
    for lab in labels:
      if lab == label:
        index = labels.index(lab)
        if default[index] == False:
          default[index] = True
        elif default[index] == True: 
          default[index] = False
    for d in range(len(default)):
      if default == [False, False]: 
        total_data = real_data[0]*0
      if default[d] == True: 
        total_data = total_data + real_data[d]
      else: 
        pass
    norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
    ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

  #Updates the scaling when the slider is changed
  def update_slider_stretch(val):
    '''Updates the stretch the user can see, based in percentiles'''
    
    total_data = 0
    for d in range(len(default)):
      if default[d] == True: 
        total_data = total_data + real_data[d]
      else: 
        pass
    norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
    ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

  #Updates the notes added by the user when there is an input
  text_list = [text]
  def submit(expression):
    '''Updates the list of types in the 'Notes' setting'''
    
    text = expression
    text_list.append(text)

  #Allows the sliders and buttons to be pressed
  check.on_clicked(update_button)
  slider_top.on_changed(update_slider_stretch)
  slider_bottom.on_changed(update_slider_stretch)
  text_box.on_text_change(submit)

  #Display image until it is clicked to find the object
  n = -1
  while True:
    press = plt.waitforbuttonpress()
    text_max = len(text_list) - 1

    #Checks that it was a mouse click
    if press == False:
      n += 3

      #Finds which axes was clicked
      click_axes = str(location[n])
      click_axes = click_axes.split('WCSAxesSubplot', 2)[0]

      #Checks if the image was clicked
      if click_axes == '':
        #Makes a pop-up window with success text
        plt.clf()
        plt.close('all')
        plt.figure(1)
        plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for CatWISE! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.grid(linewidth = 0)
        ax = plt.gca()
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.yaxis.set_tick_params(labelleft=False)
        ax.set_xticks([])
        ax.set_yticks([])
        figure2 = plt.gcf()
        figure2.set_size_inches(4.85, 1)
        # figure2.canvas.set_window_title('Successful CatWISE Search')
        mng2 = pyplot.get_current_fig_manager()
        mng2.window.resizable(False, False)
        plt.pause(0.1)
        plt.clf()
        plt.close('all')

        #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
        coord = wcs_cropped.pixel_to_world_values(location[n-4],location[n-5])
        distance = []
        for i in range(len(object_ra)):
          distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))

        list_location = distance.index(np.min(distance))
        mjd = object_mjd[list_location]
        ra_catwise, dec_catwise = object_ra[list_location], object_dec[list_location]
        ra_cw_e, dec_cw_e = object_ra_sigma[list_location], object_dec_sigma[list_location]
        w1, w2 = w1_list[list_location], w2_list[list_location]
        w1_sigma, w2_sigma = w1_list_sigma[list_location], w2_list_sigma[list_location]
        pmra, pmdec = pmra_list[list_location], pmdec_list[list_location]
        pmra_sigma, pmdec_sigma = pmra_sigma_list[list_location], pmdec_sigma_list[list_location]
        return ra_catwise, ra_cw_e, dec_catwise, dec_cw_e, w1, w1_sigma, w2, w2_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma, mjd, 'CatWISE 2020 Catalog', text_list[text_max] 
      
      #Checks if the "Object Not Found" button was clicked
      elif click_axes == 'Axes(0.04,0.775;0.92x0.04)':
        mjd, ra_cw_e, dec_cw_e, w1, w1_sigma, w2, w2_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        ra_catwise = ra
        dec_catwise = dec

        #Makes a pop-up window with success text
        plt.clf()
        plt.close('all')
        plt.figure(1)
        plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for CatWISE! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.grid(linewidth = 0)
        ax = plt.gca()
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.yaxis.set_tick_params(labelleft=False)
        ax.set_xticks([])
        ax.set_yticks([])
        figure2 = plt.gcf()
        figure2.set_size_inches(4.85, 1)
        # figure2.canvas.set_window_title('Successful CatWISE Search')
        mng2 = pyplot.get_current_fig_manager()
        mng2.window.resizable(False, False)
        plt.pause(0.1)
        plt.clf()
        plt.close('all')

        return ra_catwise, ra_cw_e, dec_catwise, dec_cw_e, w1, w1_sigma, w2, w2_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma, mjd, 'CatWISE 2020 Catalog', 'Object Not Found was Pressed'
      
      #Updates the circle size when slider is moved
      elif click_axes == 'Axes(0.25,0.055;0.65x0.03)':
        scatter.remove()
        scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')
        
    #Checks if the window was closed
    elif press is None:
      mjd, ra_cw_e, dec_cw_e, w1, w1_sigma, w2, w2_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
      ra_catwise = ra
      dec_catwise = dec
      shape_x, shape_y = total_data.shape[0], total_data.shape[1]
      plt.close('all')
      plt.figure().clear()
      return ra_catwise, ra_cw_e, dec_catwise, dec_cw_e, w1, w1_sigma, w2, w2_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma, mjd, 'CatWISE 2020 Catalog', text_list[text_max]
  
def catwise_table(ra, dec, radius): 
  '''Find all the objects in the radius defined by the user'''

  enablePrint()

  #Uses astroquery to find all objects in the radius
  location_data = Irsa.query_region(
    coord.SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='fk5'), catalog='catwise_2020', spatial='Cone', radius=(radius - 1) * u.arcsec)
  return location_data

#                         GAIA SEARCH                           #
# ------------------------------------------------------------- #
def gaia_image(ra, dec, radius): 
  ''' First, it gets the images from the GAIA API from the GAIA Archive and downloads the images. 
  Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
  Third, calls the table function to get all the objects from the GAIA source catalog. 
  Fourth, makes the window for the user to click the object with all settings. 
  Finally, finds the closest object to the click and records the data. '''

  #Makes outline for the window of the plot
  plt.rcParams['toolbar'] = 'None'
  matplotlib.use("TkAgg")
  plt.style.use('Solarize_Light2')
  blockPrint()

  #Finds all the metadata that relates to the ra and dec searched, mostly to find the APIs for the W1 and W2 images
  metadata_allwise_link = 'http://irsa.ipac.caltech.edu/ibe/sia/wise/allwise/p3am_cdd?POS=' + str(ra) + ',' + str(dec) + '&SIZE=' + str(radius/3600)
  try:
    allwise_metadata = requests.get(metadata_allwise_link)
  except: 
    ra_gaia_e, dec_gaia_e, par, par_e, rad, rad_e, pmra, pmra_e, pmdec, pmdec_e, g, g_e, bp, bp_e, rp, rp_e, year = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    ra_gaia = ra
    dec_gaia = dec
    plt.close('all')
    plt.figure().clear()
    return ra_gaia, ra_gaia_e, dec_gaia, dec_gaia_e, par, par_e, rad, rad_e, pmra, pmra_e, pmdec, pmdec_e, g, g_e, bp, bp_e, rp, rp_e, year, 'GAIA DR3 Archive', 'GAIA is Having Maintenance'
  open(str(directory) + '/Output/metadata/gaia_metadata.txt', 'wb').write(allwise_metadata.content)

  #With this metadata it finds the API link for the W1 and W2 images
  w1_finder, w2_finder, w3_finder, w4_finder = 'W1 Coadd', 'W2 Coadd', 'W3 Coadd', 'W4 Coadd'
  with open(str(directory) + '/Output/metadata/gaia_metadata.txt', 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      if line.find(w1_finder) != -1:
        w1_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]
      elif line.find(w2_finder) != -1:
        w2_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]

  #Download the W1 and W2 images
  file_allwise_w1, file_allwise_w2 = download_file(w1_allwise_image_url, cache=True), download_file(w2_allwise_image_url, cache=True)
  data_allwise_w2, data_allwise_w1 = fits.getdata(file_allwise_w2), fits.getdata(file_allwise_w1)
  
  #Find the location of all the object found in GAIA in the radius choosen by the user 
  location_data = gaia_table(ra, dec, radius)
  object_year = location_data['ref_epoch'].tolist()
  object_ra, object_ra_sigma = location_data['ra'].tolist(), location_data['ra_error']
  object_dec, object_dec_sigma = location_data['dec'].tolist(), location_data['dec_error']
  parallax_list, parallax_list_sigma = location_data['parallax'].tolist(), location_data['parallax_error'].tolist()
  rad_v_list, rad_v_list_sigma = location_data['radial_velocity'].tolist(), location_data['radial_velocity_error'].tolist()
  pmra_list, pmra_list_sigma = location_data['pmra'].tolist(), location_data['pmra_error'].tolist()
  pmdec_list, pmdec_list_sigma = location_data['pmdec'].tolist(), location_data['pmdec_error'].tolist()
  g_list, g_list_e = location_data['phot_g_mean_mag'].tolist(), location_data['phot_g_mean_flux_error'].tolist()
  bp_list, bp_list_e = location_data['phot_bp_mean_mag'].tolist(), location_data['phot_bp_mean_flux_error'].tolist()
  rp_list, rp_list_e = location_data['phot_rp_mean_mag'].tolist(), location_data['phot_rp_mean_flux_error'].tolist()
  g_list_flux, bp_list_flux, rp_list_flux = location_data['phot_g_mean_flux'].tolist(), location_data['phot_bp_mean_flux'].tolist(), location_data['phot_rp_mean_flux'].tolist()

  #Obtains the headers for each image
  hdu_w1, hdu_w2 = fits.open(file_allwise_w1)[0], fits.open(file_allwise_w2)[0]
  wcs1_w1 = WCS(hdu_w1.header)

  #Make a cutout from the coadd image for the RA and DEC put in
  position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
  size = u.Quantity([radius, radius], u.arcsec)
  cutout_w1 = Cutout2D(data_allwise_w1, position, size, fill_value = np.nan, wcs = wcs1_w1.celestial)
  cutout_w2 = Cutout2D(data_allwise_w2, position, size, fill_value = np.nan, wcs = wcs1_w1.celestial)
  wcs_cropped_w1 = cutout_w1.wcs
  enablePrint()

  #Obtains the dates for each image
  date_w1, date_w2 = hdu_w1.header['MIDOBS'].split('T', 2)[0], hdu_w2.header['MIDOBS'].split('T', 2)[0]

  #Defining a mouse click as an event on the plot
  location = []
  plt.rcParams["figure.figsize"] = [8, 8]
  plt.rcParams["figure.autolayout"] = True
  def mouse_event(event):
    '''Makes a list of the x, y, and axes the mouse click is.'''

    location.append(event.ydata)
    location.append(event.xdata)
    location.append(event.inaxes)
  plt.connect('button_press_event', mouse_event)

  #Sets the WCS coordinates for the plots
  total_data = cutout_w1.data + cutout_w2.data
  ax = plt.subplot(projection = wcs_cropped_w1)

  #Plots the objects found in the radius
  circle_size = (radius*2)
  scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

  #Normalize the image and plots it
  init_top = 95
  init_bot = 45
  norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, init_bot), vmax = np.nanpercentile(total_data.data, init_top))
  ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

  #Makes the figure look pretty
  plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
  fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
  plt.suptitle('Gaia Search', fontsize = 35, y = 0.96, fontfamily = 'Times New Roman')
  ax.set_title('Dates: \n'
             + 'W1 Date: ' + str(date_w1) + ' (Y-M-D)  \n' 
             + 'W2 Date: ' + str(date_w2) + ' (Y-M-D)  \n', fontdict = fontdict_1, y = 1.05)
  plt.grid(linewidth = 0)
  figure = plt.gcf()
  plt.xlim(0, max(total_data.shape))
  plt.ylim(0, max(total_data.shape))
  if platform != 'win32': 
    figure.set_size_inches(4.75, 7.25)
  elif platform == 'win32': 
    figure.set_size_inches(4.75, 7.25)
  # figure.canvas.set_window_title('Gaia Search')
  mng = pyplot.get_current_fig_manager()
  mng.window.resizable(False, False)
  
  #Makes a cursor to aid the User
  cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
  annotation = ax.annotate('', xy=(0, 0), xytext=(20, -20), arrowprops=dict(arrowstyle='wedge'), fontsize=12, color='red')
  annotation.set_visible(False)

  #Make checkbuttons with all of the different image bands
  rax = plt.axes([0.045, 0.4, 0.115, 0.08])
  labels = ['W1', 'W2']
  real_data = [cutout_w1.data, cutout_w2.data]
  default = [True, True]
  check = CheckButtons(rax, labels, default)

  #Adds a slider for the scaling of the image
  freq_top = plt.axes([0.25, 0.12, 0.65, 0.03])
  slider_top = Slider(ax = freq_top, label = 'Top Stetch:', valmin = 50, valmax = 100, valinit = init_top, color = '#E48671')
  freq_bottom = plt.axes([0.25, 0.087, 0.65, 0.03])
  slider_bottom = Slider(ax = freq_bottom, label = 'Bottom Stetch:', valmin = 0, valmax = 50, valinit = init_bot, color = '#E48671')

  #Adds a slider for the circle size
  circle_slid_location = plt.axes([0.25, 0.055, 0.65, 0.03])
  circle_slider = Slider(ax = circle_slid_location, label = 'Circle Size:', valmin = 5, valmax = (circle_size + 1*radius), valinit = circle_size, color = '#E48671')

  #Adds a notes section that the user can add notes about their data
  axbox = plt.axes([0.15, 0.02, 0.8, 0.03])
  text = ''
  text_box = TextBox(axbox, 'Notes:', initial = text, textalignment="center")

  #Make a button that can be clicked if no object is found
  axes_button = plt.axes([0.04, 0.775, 0.92, 0.04])
  close = Button(axes_button, 'Object Not Found', color = '#E48671')

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
    for d in range(len(default)):
      if default == [False, False]: 
        total_data = real_data[0]*0
      if default[d] == True: 
        total_data = total_data + real_data[d]
      else: 
        pass
    norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
    ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

  #Updates the scaling when the slider is changed
  def update_slider_stretch(val):
    '''Updates the stretch the user can see, based in percentiles'''

    total_data = 0
    for d in range(len(default)):
      if default[d] == True: 
        total_data = total_data + real_data[d]
      else: 
        pass
    norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
    ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

  #Updates the notes added by the user when there is an input
  text_list = [text]
  def submit(expression):
    '''Updates the list of types in the 'Notes' setting'''

    text = expression
    text_list.append(text)

  #Allows the sliders and buttons to be pressed
  check.on_clicked(update_button)
  slider_top.on_changed(update_slider_stretch)
  slider_bottom.on_changed(update_slider_stretch)
  text_box.on_text_change(submit)

  #Display image until it is clicked to find the object
  n = -1
  while True:
    press = plt.waitforbuttonpress()
    text_max = len(text_list) - 1

    #Checks that it was a mouse click
    if press == False:
      n += 3

      #Finds which axes was clicked
      click_axes = str(location[n])
      click_axes = click_axes.split('WCSAxesSubplot', 2)[0]

      #Checks if the image was clicked
      if click_axes == '':
        #Makes a pop-up window with success text
        plt.clf()
        plt.close('all')
        plt.figure(1)
        plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for Gaia! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.grid(linewidth = 0)
        ax = plt.gca()
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.yaxis.set_tick_params(labelleft=False)
        ax.set_xticks([])
        ax.set_yticks([])
        figure2 = plt.gcf()
        figure2.set_size_inches(4.5, 1)
        # figure2.canvas.set_window_title('Successful Gaia Search')
        mng2 = pyplot.get_current_fig_manager()
        mng2.window.resizable(False, False)
        plt.pause(0.1)
        plt.clf()
        plt.close('all')

        #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
        coord = wcs_cropped_w1.pixel_to_world_values(location[n-4],location[n-5])
        distance = []
        for i in range(len(object_ra)):
          distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))

        list_location = distance.index(np.min(distance))
        year = object_year[list_location]
        ra_gaia, dec_gaia = object_ra[list_location], object_dec[list_location]
        ra_gaia_e, dec_gaia_e = object_ra_sigma[list_location], object_dec_sigma[list_location]
        par, par_e = parallax_list[list_location], parallax_list_sigma[list_location]
        rad, rad_e = rad_v_list[list_location], rad_v_list_sigma[list_location]
        pmra, pmra_e = pmra_list[list_location], pmra_list_sigma[list_location]
        pmdec, pmdec_e = pmdec_list[list_location], pmdec_list_sigma[list_location]
        g = g_list[list_location]
        bp = bp_list[list_location]
        rp = rp_list[list_location]
        if g_list_e[list_location] != None:
          g_e = np.sqrt(((2.5/np.log(10))*(g_list_e[list_location]/g_list_flux[list_location]))**2 + (0.0027553202**2))
        else: 
          g_e = np.nan
        if bp_list_e[list_location] != None:
          bp_e = np.sqrt(((2.5/np.log(10))*(bp_list_e[list_location]/bp_list_flux[list_location]))**2 + (0.0027901700**2))
        else: 
          bp_e = np.nan
        if rp_list_e[list_location] != None:
          rp_e = np.sqrt(((2.5/np.log(10))*(rp_list_e[list_location]/rp_list_flux[list_location]))**2 + (0.0037793818**2))
        else: 
          rp_e = np.nan
        return ra_gaia, ra_gaia_e, dec_gaia, dec_gaia_e, par, par_e, rad, rad_e, pmra, pmra_e, pmdec, pmdec_e, g, g_e, bp, bp_e, rp, rp_e, year, 'GAIA DR3 Archive', text_list[text_max]
      
      #Checks if the Object not Found button was clicked
      elif click_axes == 'Axes(0.04,0.775;0.92x0.04)':
        ra_gaia_e, dec_gaia_e, par, par_e, rad, rad_e, pmra, pmra_e, pmdec, pmdec_e, g, g_e, bp, bp_e, rp, rp_e, year = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        ra_gaia = ra
        dec_gaia = dec

        #Makes a pop-up window with success text
        plt.clf()
        plt.close('all')
        plt.figure(1)
        plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for Gaia! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.grid(linewidth = 0)
        ax = plt.gca()
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.yaxis.set_tick_params(labelleft=False)
        ax.set_xticks([])
        ax.set_yticks([])
        figure2 = plt.gcf()
        figure2.set_size_inches(4.5, 1)
        # figure2.canvas.set_window_title('Successful Gaia Search')
        mng2 = pyplot.get_current_fig_manager()
        mng2.window.resizable(False, False)
        plt.pause(0.1)
        plt.clf()
        plt.close('all')

        return ra_gaia, ra_gaia_e, dec_gaia, dec_gaia_e, par, par_e, rad, rad_e, pmra, pmra_e, pmdec, pmdec_e, g, g_e, bp, bp_e, rp, rp_e, year, 'GAIA DR3 Archive', 'Object Not Found was Pressed'
      
      #Updates the circle size when slider is moved
      elif click_axes == 'Axes(0.25,0.055;0.65x0.03)':
        scatter.remove()
        scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')
        
    #Checks if the window was closed
    elif press is None:
      ra_gaia_e, dec_gaia_e, par, par_e, rad, rad_e, pmra, pmra_e, pmdec, pmdec_e, g, g_e, bp, bp_e, rp, rp_e, year = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
      ra_gaia = ra
      dec_gaia = dec
      plt.close('all')
      plt.figure().clear()
      return ra_gaia, ra_gaia_e, dec_gaia, dec_gaia_e, par, par_e, rad, rad_e, pmra, pmra_e, pmdec, pmdec_e, g, g_e, bp, bp_e, rp, rp_e, year, 'GAIA DR3 Archive', text_list[text_max]

def gaia_table(ra, dec, radius): 
  '''Find all the objects in the radius defined by the user'''

  blockPrint()

  #Makes the SQL code to run it into the GAIA search
  query = "SELECT TOP 2000 \
  gaia_source.ra,gaia_source.dec,gaia_source.ra_error,gaia_source.phot_g_mean_flux_error,gaia_source.phot_bp_mean_flux_error,gaia_source.phot_rp_mean_flux_error,gaia_source.ref_epoch,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,gaia_source.phot_g_mean_mag,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.phot_g_mean_flux,gaia_source.phot_bp_mean_flux,gaia_source.phot_rp_mean_flux \
  FROM gaiadr3.gaia_source \
  WHERE \
  CONTAINS( \
  POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec), \
  CIRCLE('ICRS', " + str(ra) + "," + str(dec) + "," + str(((radius/2) * 0.000277778) - 0.000277778)+ ")" \
  ")=1"

  #Run this SQL quiery into the online GAIA database
  job = Gaia.launch_job(query)
  results = job.get_results()
  return results

#                        GALEX SEARCH                           #
# ------------------------------------------------------------- #
def galex_image(ra, dec, radius_use): 
    ''' First, it gets the images from the GALEX API from MAST and downloads the images. 
    Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
    Third, calls the table function to get all the objects from the GALEX source catalog. 
    Fourth, makes the window for the user to click the object with all settings. 
    Finally, finds the closest object to the click and records the data. '''

    #Makes outline for the window of the plot
    plt.rcParams['toolbar'] = 'None'
    plt.style.use('Solarize_Light2')
    matplotlib.use("TkAgg")

    blockPrint()
    #Obtains all of the observations from MAST
    try:
      obs_table = Observations.query_region(str(ra) + ' ' + str(dec), radius = (radius_use/1800))
    except: 
      fuv, fuv_e, nuv, nuv_e = np.nan, np.nan, np.nan, np.nan
      ra_galex, dec_galex = ra, dec
      plt.close('all')
      plt.figure().clear()
      return ra_galex, dec_galex, fuv, fuv_e, nuv, nuv_e, 'MAST GALEX Data Release 7', 'GALEX is Having Maintenance'
    obs_table = pd.DataFrame(data = np.array(obs_table))
    
    #Finds only the GALEX images
    telescope = obs_table['obs_collection'].tolist()
    for p in range(len(telescope)):
        if telescope[p] == 'GALEX':
            pass
        else: 
            obs_table = obs_table.drop(p)

    #Finds the real images from GALEX
    urls = obs_table['dataURL'].tolist()
    obs_table = pd.DataFrame(data = np.array(obs_table))
    for n in range(len(urls)):
        temp_url = urls[n].split('-', 7)
        if temp_url[len(temp_url) - 1] == 'int.fits.gz':
            pass
        else: 
            obs_table = obs_table.drop(n)
    
    #Finds all of the NUV band urls
    nuv_list_url = []
    urls = obs_table[27].tolist()
    bands = obs_table[5].tolist()
    for n in range(len(bands)):
        if bands[n] == 'NUV': 
            nuv_list_url.append(urls[n])
        
    #Runs the code if an image is found
    if len(nuv_list_url) > 0: 
        #Find the location of all the object found in GALEX in the radius choosen by the user 
        observation_table = galex_table(ra, dec, radius_use)
        object_ra, object_dec = observation_table['ra'].tolist(), observation_table['dec'].tolist()
        fuv_list, fuv_list_e = observation_table['fuv_mag'].tolist(), observation_table['fuv_magerr'].tolist()
        nuv_list, nuv_list_e = observation_table['nuv_mag'].tolist(), observation_table['nuv_magerr'].tolist()

        #Combines all of the real images
        total_data = 0
        for i in range(len(nuv_list_url)):
            #Downloads the images
            file_allwise_nuv = download_file(nuv_list_url[i], cache=True)
            data_allwise_nuv = fits.getdata(file_allwise_nuv)

            #Make a cutout from the coadd image for the RA and DEC put in
            hdu_w1 = fits.open(file_allwise_nuv)[0]
            wcs1_w1 = WCS(hdu_w1.header)

            #Makes a cutout for the images
            position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
            size = u.Quantity([radius_use, radius_use], u.arcsec)
            cutout_w1 = Cutout2D(data_allwise_nuv, position, size, fill_value = np.nan, wcs = wcs1_w1.celestial)

            #Combines the data
            total_data = cutout_w1.data + total_data

            #Gets the new WCS coordinates for the image cutout
            wcs_cropped_w1 = cutout_w1.wcs

        #Finds the date for the last image
        date_w1 = hdu_w1.header['OBS-DATE']
        enablePrint()

        #Defining a mouse click as an event on the plot
        location = []
        plt.rcParams["figure.figsize"] = [8, 8]
        plt.rcParams["figure.autolayout"] = True
        def mouse_event(event):
            '''Makes a list of the x, y, and axes the mouse click is.'''

            location.append(event.ydata)
            location.append(event.xdata)
            location.append(event.inaxes)
        plt.connect('button_press_event', mouse_event)

        #Sets the WCS coordinates for the plots
        ax = plt.subplot(projection = wcs_cropped_w1)

        #Plots the objects found in the radius
        circle_size = (radius_use*4)
        scatter = plt.scatter(object_ra, object_dec, transform = ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

        #Normalize the image and plots it
        init_top, init_bot = 95, 45
        norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data, init_bot), vmax = np.nanpercentile(total_data, init_top))
        plt.imshow(total_data, cmap = 'Greys', norm = norm1_w1)

        #Makes the figure look pretty
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        fontdict = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
        plt.suptitle('GALEX Search', fontsize = 35, y = 0.95, fontfamily = 'Times New Roman')
        ax.set_title('Dates: \n'
                   + 'FUV Date: ' + str(date_w1) + ' (Y-M-D)  ', fontdict = fontdict, y = 1.105)
        plt.grid(linewidth = 0)
        figure = plt.gcf()
        plt.xlim(0, max(total_data.shape))
        plt.ylim(0, max(total_data.shape))
        if platform != 'win32':
            figure.set_size_inches(4.75, 7.25)
        elif platform == 'win32':
            figure.set_size_inches(4.75, 7.25)
        # figure.canvas.set_window_title('GALEX Search')
        mng = pyplot.get_current_fig_manager()
        mng.window.resizable(False, False)
        
        #Makes a cursor to aid the User
        cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
        annotation = ax.annotate('', xy=(0, 0), xytext=(20, -20), arrowprops=dict(arrowstyle='wedge'), fontsize=12, color='red')
        annotation.set_visible(False)

        #Adds a slider for the scaling of the image
        freq_top = plt.axes([0.25, 0.12, 0.65, 0.03])
        slider_top = Slider(ax = freq_top, label = 'Top Stetch:', valmin = 50, valmax = 100, valinit = init_top, color = '#E48671')
        freq_bottom = plt.axes([0.25, 0.087, 0.65, 0.03])
        slider_bottom = Slider(ax = freq_bottom, label = 'Bottom Stetch:', valmin = 0, valmax = 50, valinit = init_bot, color = '#E48671')

        #Adds a slider for the circle size
        circle_slid_location = plt.axes([0.25, 0.055, 0.65, 0.03])
        circle_slider = Slider(ax = circle_slid_location, label = 'Circle Size:', valmin = (circle_size - 2.5*radius_use), valmax = (circle_size + 1*radius_use), valinit = circle_size, color = '#E48671')

        #Adds a notes section that the user can add notes about their data
        axbox = plt.axes([0.15, 0.02, 0.8, 0.03])
        text = ''
        text_box = TextBox(axbox, 'Notes:', initial = text, textalignment="center")

        #Make a button that can be clicked if no object is found
        axes_button = plt.axes([0.04, 0.775, 0.92, 0.04])
        close = Button(axes_button, 'Object Not Found', color = '#E48671')

        #Updates the scaling when the slider is changed
        def update_slider_stretch(val):
            '''Updates the stretch the user can see, based in percentiles'''

            norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
            ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

        #Updates the notes added by the user when there is an input
        text_list = [text]
        def submit(expression):
            '''Updates the list of types in the 'Notes' setting'''

            text = expression
            text_list.append(text)

        #Allows the sliders and buttons to be pressed
        slider_top.on_changed(update_slider_stretch)
        slider_bottom.on_changed(update_slider_stretch)
        text_box.on_text_change(submit)

        n = -1
        while True:
            press = plt.waitforbuttonpress()
            text_max = len(text_list) - 1

            #Checks that it was a mouse click
            if press == False:
                n += 3

                #Finds which axes was clicked
                click_axes = str(location[n])
                click_axes = click_axes.split('WCSAxesSubplot', 2)[0]

                #Checks if the image was clicked
                if click_axes == '':
                    #Makes a pop-up window with success text
                    plt.clf()
                    plt.close('all')
                    plt.figure(1)
                    plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for GALEX! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
                    plt.xlim(0, 1)
                    plt.ylim(0, 1)
                    plt.grid(linewidth = 0)
                    ax = plt.gca()
                    ax.xaxis.set_tick_params(labelbottom=False)
                    ax.yaxis.set_tick_params(labelleft=False)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    figure2 = plt.gcf()
                    figure2.set_size_inches(4.75, 1)
                    # figure2.canvas.set_window_title('Successful GALEX Search')
                    mng2 = pyplot.get_current_fig_manager()
                    mng2.window.resizable(False, False)
                    plt.pause(0.1)
                    plt.clf()
                    plt.close('all')

                    #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
                    coord = wcs_cropped_w1.pixel_to_world_values(location[n-4], location[n-5])
                    distance = []
                    for i in range(len(object_ra)):
                        distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))
                    list_location = distance.index(np.min(distance))
                    ra_galex, dec_galex = object_ra[list_location], object_dec[list_location]
                    fuv, fuv_e = fuv_list[list_location], fuv_list_e[list_location]
                    nuv, nuv_e = nuv_list[list_location], nuv_list_e[list_location]
                    return ra_galex, dec_galex, fuv, fuv_e, nuv, nuv_e, 'MAST GALEX Data Release 7', text_list[text_max]
                
                #Checks if the "Object Not Found" button was clicked
                elif click_axes == 'Axes(0.04,0.775;0.92x0.04)':
                    fuv, fuv_e, nuv, nuv_e = np.nan, np.nan, np.nan, np.nan
                    ra_galex, dec_galex = ra, dec

                    #Makes a pop-up window with success text
                    plt.clf()
                    plt.close('all')
                    plt.figure(1)
                    plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for GALEX! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
                    plt.xlim(0, 1)
                    plt.ylim(0, 1)
                    plt.grid(linewidth = 0)
                    ax = plt.gca()
                    ax.xaxis.set_tick_params(labelbottom=False)
                    ax.yaxis.set_tick_params(labelleft=False)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    figure2 = plt.gcf()
                    figure2.set_size_inches(4.75, 1)
                    # figure2.canvas.set_window_title('Successful GALEX Search')
                    mng2 = pyplot.get_current_fig_manager()
                    mng2.window.resizable(False, False)
                    plt.pause(0.1)
                    plt.clf()
                    plt.close('all')
        
                    return ra_galex, dec_galex, fuv, fuv_e, nuv, nuv_e, 'MAST GALEX Release Data 7', 'Object Not Found was Pressed'
                
                #Updates the circle size when slider is moved
                elif click_axes == 'Axes(0.25,0.055;0.65x0.03)':
                    scatter.remove()
                    scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')

            #Checks if the window was closed
            elif press is None:
                fuv, fuv_e, nuv, nuv_e = np.nan, np.nan, np.nan, np.nan
                ra_galex, dec_galex = ra, dec
                plt.close('all')
                plt.figure().clear()
                return ra_galex, dec_galex, fuv, fuv_e, nuv, nuv_e, 'MAST GALEX Data Release 7', text_list[text_max]
                
    #Returns null values if the images are not found    
    else: 
        fuv, fuv_e, nuv, nuv_e = np.nan, np.nan, np.nan, np.nan
        ra_galex, dec_galex = ra, dec
        return ra_galex, dec_galex, fuv, fuv_e, nuv, nuv_e, 'MAST GALEX Data Release 7', 'Image Not Found'

def galex_table(ra, dec, radius): 
    '''Find all the objects in the radius defined by the user'''
    
    blockPrint()

    #Gets the table of all of the data in the radius from the user
    catalog_data = Catalogs.query_region(str(ra) + ' ' + str(dec), radius = ((radius/7200) - 0.000277778), catalog = "Galex", table = "mean")
    return catalog_data

#                          NSC SEARCH                           #
# ------------------------------------------------------------- #
def nsc_image(ra, dec, radius): 
    ''' First, it gets the images from the NSC API from AstroLab and downloads the images. 
    Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
    Third, calls the table function to get all the objects from the NSC source catalog. 
    Fourth, makes the window for the user to click the object with all settings. 
    Finally, finds the closest object to the click and records the data. '''
    
    #Makes outline for the window of the plot
    plt.rcParams['toolbar'] = 'None'
    matplotlib.use("TkAgg")
    plt.style.use('Solarize_Light2')
    blockPrint()

    #Defines the catalog that is searched
    DEF_ACCESS_URL = "https://datalab.noirlab.edu/sia/des_dr1"
    try:
      svc = sia.SIAService(DEF_ACCESS_URL)
    except: 
      ra_nsc_e, dec_nsc_e, g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e, mjd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
      ra_nsc, dec_nsc = ra, dec
      return ra_nsc, ra_nsc_e, dec_nsc, dec_nsc_e, g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e, mjd, 'NoirLab Source Catalog Data Release 2', 'NSC is Having Maintenance'

    #Finds all of the image urls for the ra, dec, and radius given
    imgTable = svc.search((ra,dec), (radius/3600)).to_table()

    #Tests if any images were found
    if len(imgTable) > 0:
        #Runs the table function to obtain the g, r, i, and z bands
        location_data = nsc_table(ra, dec, radius)
        mjd_list = location_data['mjd'].tolist()
        object_ra, object_dec = location_data['ra'].tolist(), location_data['dec'].tolist()
        object_ra_e, object_dec_e = location_data['raerr'].tolist(), location_data['decerr'].tolist()
        g_list, g_list_e = location_data['gmag'].tolist(), location_data['gerr'].tolist()
        r_list, r_list_e = location_data['rmag'].tolist(), location_data['rerr'].tolist()
        i_list, i_list_e = location_data['imag'].tolist(), location_data['ierr'].tolist()
        z_list, z_list_e = location_data['zmag'].tolist(), location_data['zerr'].tolist()
        u_list, u_list_e = location_data['umag'].tolist(), location_data['uerr'].tolist()
        y_list, y_list_e = location_data['ymag'].tolist(), location_data['yerr'].tolist()
        pmra_list, pmra_list_e = location_data['pmra'].tolist(), location_data['pmraerr'].tolist()
        pmdec_list, pmdec_list_e = location_data['pmdec'].tolist(), location_data['pmdecerr'].tolist()

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
        file_allwise_g, file_allwise_r, file_allwise_i = download_file(image_url_g, cache=True), download_file(image_url_r, cache=True), download_file(image_url_i, cache=True)
        file_allwise_z, file_allwise_Y = download_file(image_url_z, cache=True), download_file(image_url_Y, cache=True)
        data_allwise_g, data_allwise_r, data_allwise_i = fits.getdata(file_allwise_g), fits.getdata(file_allwise_r), fits.getdata(file_allwise_i)
        data_allwise_z, data_allwise_Y = fits.getdata(file_allwise_z), fits.getdata(file_allwise_Y)

        #Loads the WCS from the i band image
        hdu_g, hdu_r, hdu_i = fits.open(file_allwise_r)[0], fits.open(file_allwise_i)[0], fits.open(file_allwise_g)[0]
        hdu_z, hdu_Y = fits.open(file_allwise_z)[0], fits.open(file_allwise_Y)[0]
        wcs = WCS(hdu_g.header)

        #Makes the cutouts
        position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
        size = u.Quantity([radius, radius], u.arcsec)
        cutout_g = Cutout2D(data_allwise_g, position, size, fill_value = np.nan, wcs = wcs.celestial)
        cutout_r = Cutout2D(data_allwise_r, position, size, fill_value = np.nan, wcs = wcs.celestial)
        cutout_i = Cutout2D(data_allwise_i, position, size, fill_value = np.nan, wcs = wcs.celestial)
        cutout_z = Cutout2D(data_allwise_z, position, size, fill_value = np.nan, wcs = wcs.celestial)
        cutout_Y = Cutout2D(data_allwise_Y, position, size, fill_value = np.nan, wcs = wcs.celestial)
        wcs_cropped_i = cutout_i.wcs
        enablePrint()

        #Obtains the dates for each image
        date_g, date_r, date_i = hdu_g.header[36].split('T', 2)[0], hdu_r.header[36].split('T', 2)[0], hdu_i.header[36].split('T', 2)[0]
        date_z, date_Y = hdu_z.header[36].split('T', 2)[0], hdu_Y.header[36].split('T', 2)[0]

        #Allows the user to click the image to find the object
        location = []
        plt.rcParams["figure.figsize"] = [8, 8]
        plt.rcParams["figure.autolayout"] = True
        def mouse_event(event):
            '''Makes a list of the x, y, and axes the mouse click is.'''

            location.append(event.ydata)
            location.append(event.xdata)
            location.append(event.inaxes)
        fig_1 = plt.figure()
        cid = fig_1.canvas.mpl_connect('button_press_event', mouse_event)

        #Sets the WCS coordinates for the plots
        total_data = cutout_i.data + cutout_z.data + cutout_r.data
        ax = plt.subplot(projection = wcs_cropped_i)
        
        #Plots the objects found in the radius
        circle_size = (radius*3)
        scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

        #Normalize the image and plots it
        init_bot, init_top = 45, 95
        norm1_total = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, init_bot), vmax = np.nanpercentile(total_data.data, init_top))
        ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_total)

        #Makes the figure look pretty
        pixel_radius = 3.785*radius
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
        plt.suptitle('NSC Search', fontsize = 35, y = 0.96, fontfamily = 'Times New Roman')
        ax.set_title('Dates: \n'
                   + 'g Date: ' + str(date_g) + ' (Y-M-D)  ' + '  r Date: ' + str(date_r) + ' (Y-M-D)\n'
                   + 'i Date: ' + str(date_i) + ' (Y-M-D)  ' + '  z Date: ' + str(date_z) + ' (Y-M-D)\n'
                   + 'Y Date: ' + str(date_Y) + ' (Y-M-D)\n', fontdict = fontdict_1, y = 1.04)
        plt.grid(linewidth = 0)
        figure = plt.gcf()
        plt.xlim(pixel_radius, 0)
        plt.ylim(0, pixel_radius)
        figure.set_size_inches(4.75, 7.35)
        # figure.canvas.set_window_title('NSC Search')
        mng = pyplot.get_current_fig_manager()
        mng.window.resizable(False, False)
        
        #Makes a cursor to aid the User
        cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
        annotation = ax.annotate('', xy=(0, 0), xytext=(20, -20), arrowprops=dict(arrowstyle='wedge'), fontsize=12, color='red')
        annotation.set_visible(False)

        #Make checkbuttons with all of the different image bands
        rax = plt.axes([0.045, 0.4, 0.12, 0.14])
        labels = ['g', 'r', 'i', 'z', 'Y']
        real_data = [cutout_g.data, cutout_r.data, cutout_i.data, cutout_z.data, cutout_Y.data]
        default = [False, True, True, True, False]
        check = CheckButtons(rax, labels, default)

        #Adds a slider for the scaling of the image
        freq_top = plt.axes([0.25, 0.12, 0.65, 0.03])
        slider_top = Slider(ax = freq_top, label = 'Top Stetch:', valmin = 50, valmax = 100, valinit = init_top, color = '#E48671')
        freq_bottom = plt.axes([0.25, 0.087, 0.65, 0.03])
        slider_bottom = Slider(ax = freq_bottom, label = 'Bottom Stetch:', valmin = 0, valmax = 50, valinit = init_bot, color = '#E48671')

        #Adds a slider for the circle size
        circle_slid_location = plt.axes([0.25, 0.055, 0.65, 0.03])
        circle_slider = Slider(ax = circle_slid_location, label = 'Circle Size:', valmin = (circle_size - 2.5*radius), valmax = (circle_size + 1*radius), valinit = circle_size, color = '#E48671')

        #Adds a notes section that the user can add notes about their data
        axbox = plt.axes([0.15, 0.02, 0.8, 0.03])
        text = ''
        text_box = TextBox(axbox, 'Notes:', initial = text, textalignment="center")

        #Make a button that can be clicked if no object is found
        axes_button = plt.axes([0.04, 0.767, 0.92, 0.04])
        close = Button(axes_button, 'Object Not Found', color = '#E48671')

        #Update the image depending on what the user chooses
        def update_button(label):
            '''Updates the list of activated images and updates the image the user can see.'''

            total_data = 0
            for lab in labels:
                if lab == label:
                    index = labels.index(lab)
                    if default[index] == False:
                        default[index] = True
                    elif default[index] == True: 
                        default[index] = False
            for d in range(len(default)):
                if default == [False, False, False, False, False]: 
                    total_data = real_data[0]*0
                if default[d] == True: 
                    total_data = total_data + real_data[d]
                else: 
                    pass
            norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
            ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

        #Updates the scaling when the slider is changed
        def update_slider_stretch(val):
            '''Updates the stretch the user can see, based in percentiles'''

            total_data = 0
            for d in range(len(default)):
                if default[d] == True: 
                    total_data = total_data + real_data[d]
                else: 
                    pass
            norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
            ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

        #Updates the notes added by the user when there is an input
        text_list = [text]
        def submit(expression):
            '''Updates the list of types in the 'Notes' setting'''

            text = expression
            text_list.append(text)
                
        #Allows the sliders and buttons to be pressed
        check.on_clicked(update_button)
        slider_top.on_changed(update_slider_stretch)
        slider_bottom.on_changed(update_slider_stretch)
        text_box.on_text_change(submit)

        #Display image until it is clicked to find the object
        n = -1
        while True:
            press = plt.waitforbuttonpress()
            text_max = len(text_list) - 1

            #Checks that it was a mouse click
            if press == False:
                n += 3
                
                #Finds which axes was clicked
                click_axes = str(location[n])
                click_axes = click_axes.split('WCSAxesSubplot', 2)[0]

                #Checks if the image was clicked
                if click_axes == '':
                    #Makes a pop-up window with success text
                    plt.clf()
                    plt.close('all')
                    plt.figure(1)
                    plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for NSC! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
                    plt.xlim(0, 1)
                    plt.ylim(0, 1)
                    plt.grid(linewidth = 0)
                    ax = plt.gca()
                    ax.xaxis.set_tick_params(labelbottom=False)
                    ax.yaxis.set_tick_params(labelleft=False)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    figure2 = plt.gcf()
                    figure2.set_size_inches(4.5, 1)
                    # figure2.canvas.set_window_title('Successful NSC Search')
                    mng2 = pyplot.get_current_fig_manager()
                    mng2.window.resizable(False, False)
                    plt.pause(0.1)
                    plt.clf()
                    plt.close('all')

                    #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
                    coord = wcs_cropped_i.pixel_to_world_values(location[n-4],location[n-5])
                    distance = []
                    for i in range(len(object_ra)):
                        distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))

                    list_location = distance.index(np.min(distance))
                    mjd = mjd_list[list_location]
                    ra_nsc, dec_nsc = object_ra[list_location], object_dec[list_location]
                    ra_nsc_e, dec_nsc_e = object_ra_e[list_location], object_dec_e[list_location]
                    g, g_e = g_list[list_location], g_list_e[list_location]
                    r, r_e = r_list[list_location], r_list_e[list_location]
                    i, i_e = i_list[list_location], i_list_e[list_location]
                    z, z_e = z_list[list_location], z_list_e[list_location]
                    u_mag, u_mag_e = u_list[list_location], u_list_e[list_location]
                    y, y_e = y_list[list_location], y_list_e[list_location]
                    pmra, pmra_e = pmra_list[list_location], pmra_list_e[list_location]
                    pmdec, pmdec_e = pmdec_list[list_location], pmdec_list_e[list_location]
                    return ra_nsc, ra_nsc_e, dec_nsc, dec_nsc_e, g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e, mjd, 'NoirLab Source Catalog Data Release 2', text_list[text_max]
                
                #Checks if the Object not Found button was clicked
                elif click_axes == 'Axes(0.04,0.767;0.92x0.04)':
                    ra_nsc_e, dec_nsc_e, g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e, mjd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                    ra_nsc, dec_nsc = ra, dec

                    #Makes a pop-up window with success text
                    plt.clf()
                    plt.close('all')
                    plt.figure(1)
                    plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for NSC! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
                    plt.xlim(0, 1)
                    plt.ylim(0, 1)
                    plt.grid(linewidth = 0)
                    ax = plt.gca()
                    ax.xaxis.set_tick_params(labelbottom=False)
                    ax.yaxis.set_tick_params(labelleft=False)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    figure2 = plt.gcf()
                    figure2.set_size_inches(4.5, 1)
                    # figure2.canvas.set_window_title('Successful NSC Search')
                    mng2 = pyplot.get_current_fig_manager()
                    mng2.window.resizable(False, False)
                    plt.pause(0.1)
                    plt.clf()
                    plt.close('all')
        
                    return ra_nsc, ra_nsc_e, dec_nsc, dec_nsc_e, g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e, mjd, 'NoirLab Source Catalog Data Release 2', 'Object Not Found was Pressed'
                
                #Changes the circle size if the slider is changed
                elif click_axes == 'Axes(0.25,0.055;0.65x0.03)':
                    scatter.remove()
                    scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')
                
            #Checks if the window was closed
            elif press is None:
                ra_nsc_e, dec_nsc_e, g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e, mjd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                ra_nsc, dec_nsc = ra, dec
                plt.close('all')
                plt.figure().clear()
                return ra_nsc, ra_nsc_e, dec_nsc, dec_nsc_e, g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e, mjd, 'NoirLab Source Catalog Data Release 2', text_list[text_max]
            
    #Returns null values if image is not found
    else: 
        ra_nsc_e, dec_nsc_e, g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e, mjd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        ra_nsc, dec_nsc = ra, dec
        return ra_nsc, ra_nsc_e, dec_nsc, dec_nsc_e, g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e, mjd, 'NoirLab Source Catalog Data Release 2', 'Image Not Found'

def nsc_table(ra, dec, radius): 
    '''Find all the objects in the radius defined by the user'''

    blockPrint()

    #Makes a SQL query using the ra, dec, and radius
    query = " \
    SELECT ra, dec, gmag, gerr, rmag, rerr, imag, ierr, zmag, zerr, umag, uerr, ymag, yerr, pmra, pmraerr, pmdec, pmdecerr, raerr, decerr, mjd  \
    FROM nsc_dr2.object \
    WHERE ra > " + str((ra - ((radius/7200) - 0.0001))) + " and ra < " + str((ra + ((radius/7200) - 0.0001))) + " " \
    "AND dec > " + str((dec - ((radius/7200) - 0.0001))) + " and dec < " + str((dec + ((radius/7200) - 0.0001))) + " " \
    ""

    #Run this SQL quiery into the online NSC database
    response = qc.query(sql=query,format='csv')
    df = convert(response,'pandas')
    return df.head(100000000000)

#                    PANSTARRS SEARCH                           #
# ------------------------------------------------------------- #
def ps_image(ra, dec, radius):
  ''' First, it gets the images from the PanSTARRS API from MAST and downloads the images. 
  Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
  Third, calls the table function to get all the objects from the PanSTARRS source catalog. 
  Fourth, makes the window for the user to click the object with all settings. 
  Finally, finds the closest object to the click and records the data. '''
  
  #Makes outline for the window of the plot
  plt.rcParams['toolbar'] = 'None'
  matplotlib.use("TkAgg")
  plt.style.use('Solarize_Light2')
  blockPrint()

  #Makes the dec usable for the image url API
  if dec > 0: 
    new_dec = '+' + str(dec)
  elif dec <= 0: 
    new_dec = str(dec)

  #Find the panstarr metadata for the image API
  ps_image_url = 'http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=' + str(ra) + str(new_dec)  + '&filter=color&filter=g&filter=r&filter=i&filter=z&filter=y&filetypes=stack&auxiliary=data&size=' + str(radius * 4) + '&output_size=0&verbose=0&autoscale=99.500000&catlist='
  try:
    allwise_metadata = requests.get(ps_image_url)
  except: 
    ps_ra_e, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    ps_ra, ps_dec = ra, dec
    plt.close('all')
    plt.figure().clear()
    return ps_ra, ps_ra_e, ps_dec, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd, 'PanSTARRS Data Release 2', 'PanSTARRS is having Maintenance'
  open(str(directory) + '/Output/metadata/ps_metadata.txt', 'wb').write(allwise_metadata.content)

  #With the metadata from panstarrs find the image API in r band
  r_finder= 'amp;'
  r_finder_list = []
  with open(str(directory) + '/Output/metadata/ps_metadata.txt', 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      if line.find(r_finder) != -1:
        r_finder_list.append((lines[lines.index(line)]))

  #Checks if the images are found
  if len(r_finder_list) > 0: 

    #Obtains the data table for panstarrs and then gets all of the relatvent data from the table
    panstarr_data = ps_table(ra, dec, radius)
    object_mjd = panstarr_data['epochMean'].tolist()
    ra_list, dec_list = panstarr_data['raStack'].tolist(), panstarr_data['decStack'].tolist()
    ra_list_e, dec_list_e = panstarr_data['raStackErr'].tolist(), panstarr_data['decStackErr'].tolist()
    g_list, g_list_e = panstarr_data['gApMag'].tolist(), panstarr_data['gApMagErr'].tolist()
    r_list, r_list_e = panstarr_data['rApMag'].tolist(), panstarr_data['rApMagErr'].tolist()
    i_list, i_list_e = panstarr_data['iApMag'].tolist(), panstarr_data['iApMagErr'].tolist()
    z_list, z_list_e = panstarr_data['zApMag'].tolist(), panstarr_data['zApMagErr'].tolist()
    y_list, y_list_e = panstarr_data['yApMag'].tolist(), panstarr_data['yApMagErr'].tolist()

    #Gets all of the image urls
    if len(r_finder_list) < 6: 
      ps_ra_e, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
      ps_ra, ps_dec = ra, dec
      return ps_ra, ps_ra_e, ps_dec, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd, 'PanSTARRS Data Release 2', 'Image Not Found'
    r_link = ('http:' + (r_finder_list[2]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)
    i_link = ('http:' + (r_finder_list[3]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)
    z_link = ('http:' + (r_finder_list[4]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)
    y_link = ('http:' + (r_finder_list[5]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)

    if r_link[0:7] == 'http://':

      #Downloads the panstarrs r, i, and z image
      file_ps_r, file_ps_i, file_ps_z, file_ps_y = download_file(r_link, cache=True), download_file(i_link, cache=True), download_file(z_link, cache=True), download_file(y_link, cache=True)
      data_ps_r, data_ps_i, data_ps_z, data_ps_y = fits.getdata(file_ps_r), fits.getdata(file_ps_i), fits.getdata(file_ps_z), fits.getdata(file_ps_y)

      #Removes the null RA and DEC values that PanSTARRS return
      for elem in list(ra_list):
        if elem == -999.0 or elem is None:
          ra_location = ra_list.index(elem)
          ra_list.remove(elem)
          dec_list.pop(ra_location)
          g_list.pop(ra_location)
          g_list_e.pop(ra_location)
          r_list.pop(ra_location)
          r_list_e.pop(ra_location)
          i_list.pop(ra_location)
          i_list_e.pop(ra_location)
          z_list.pop(ra_location)
          z_list_e.pop(ra_location)
          y_list.pop(ra_location)
          y_list_e.pop(ra_location)

      if len(ra_list) == 0: 
        ps_ra_e, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        ps_ra, ps_dec = ra, dec
        plt.close('all')
        plt.figure().clear()
        return ps_ra, ps_ra_e, ps_dec, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd, 'PanSTARRS Data Release 2', 'Image Found but Array is Null'

      #Gets the headers from the images
      hdu_r, hdu_i, hdu_z, hdu_y = fits.open(file_ps_r)[0], fits.open(file_ps_i)[0], fits.open(file_ps_z)[0], fits.open(file_ps_y)[0]
      wcs1 = WCS(hdu_r.header)

      #Make a cutout from the coadd image for the ra and dec put in
      if type(dec) == str: 
        dec = float(dec.replace('+', '', 1))
      position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
      size = u.Quantity([radius, radius], u.arcsec)
      cutout_r = Cutout2D(data_ps_r, position, size, fill_value = np.nan, wcs = wcs1.celestial)
      cutout_i = Cutout2D(data_ps_i, position, size, fill_value = np.nan, wcs = wcs1.celestial)
      cutout_z = Cutout2D(data_ps_z, position, size, fill_value = np.nan, wcs = wcs1.celestial)
      cutout_y = Cutout2D(data_ps_y, position, size, fill_value = np.nan, wcs = wcs1.celestial)
      wcs_cropped = cutout_r.wcs
      enablePrint()

      #Gets the dates that each image was taken
      def mjd_to_date(mjd):
        time = Time(mjd, format='mjd')
        iso_date = time.iso
        return iso_date
      
      date_r = (mjd_to_date(hdu_r.header['MJD-OBS'])).split(' ', 2)[0]
      date_i = (mjd_to_date(hdu_i.header['MJD-OBS'])).split(' ', 2)[0]
      date_z = (mjd_to_date(hdu_z.header['MJD-OBS'])).split(' ', 2)[0]
      date_y = (mjd_to_date(hdu_y.header['MJD-OBS'])).split(' ', 2)[0]
  
      #Defining a mouse click as an event on the plot
      location = []
      plt.rcParams["figure.figsize"] = [8, 8]
      plt.rcParams["figure.autolayout"] = True
      def mouse_event(event):
        '''Makes a list of the x, y, and axes the mouse click is.'''

        location.append(event.ydata)
        location.append(event.xdata)
        location.append(event.inaxes)
      fig_1 = plt.figure()
      cid = fig_1.canvas.mpl_connect('button_press_event', mouse_event)

      #Sets the WCS coordinates for the plots
      total_data = cutout_r.data + cutout_i.data
      ax = plt.subplot(projection = wcs_cropped)
      
      #Plots the objects found in the radius
      circle_size = (radius*3)
      scatter = ax.scatter(ra_list, dec_list, transform=ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

      #Normalize the image and plots it
      init_top = 95
      init_bot = 45
      norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, init_bot), vmax = np.nanpercentile(total_data.data, init_top))
      ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

      #Makes the figure look pretty
      plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
      plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
      fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
      plt.suptitle('PanSTARRS Search', fontsize = 35, y = 0.96, fontfamily = 'Times New Roman')
      ax.set_title('Dates: \n'
                + 'r Date: ' + str(date_r) + ' (Y-M-D)  ' + '  i Date: ' + str(date_i) + ' (Y-M-D)   \n'
                + 'z Date: ' + str(date_z) + ' (Y-M-D)  ' + '  y Date: ' + str(date_y) + ' (Y-M-D)   \n', fontdict = fontdict_1, y = 1.05)
      plt.grid(linewidth = 0)
      figure = plt.gcf()
      plt.xlim(0, max(total_data.shape))
      plt.ylim(0, max(total_data.shape))
      if platform != 'win32':
        figure.set_size_inches(4.75, 7.25)
      elif platform == 'win32':
        figure.set_size_inches(4.75, 7.25)
      # figure.canvas.set_window_title('PanSTARRS Search')
      mng = pyplot.get_current_fig_manager()
      mng.window.resizable(False, False)
      
      #Makes a cursor to aid the User
      cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
      annotation = ax.annotate('', xy=(0, 0), xytext=(20, -20), arrowprops=dict(arrowstyle='wedge'), fontsize=12, color='red')
      annotation.set_visible(False)

      #Make checkbuttons with all of the different image bands
      rax = plt.axes([0.045, 0.4, 0.115, 0.1])
      labels = ['r', 'i', 'z', 'y']
      real_data = [cutout_r.data, cutout_i.data, cutout_z.data, cutout_y.data]
      default = [True, True, False, False]
      check = CheckButtons(rax, labels, default)

      #Adds a slider for the scaling of the image
      freq_top = plt.axes([0.25, 0.12, 0.65, 0.03])
      slider_top = Slider(ax = freq_top, label = 'Top Stetch:', valmin = 50, valmax = 100, valinit = init_top, color = '#E48671')
      freq_bottom = plt.axes([0.25, 0.087, 0.65, 0.03])
      slider_bottom = Slider(ax = freq_bottom, label = 'Bottom Stetch:', valmin = 0, valmax = 50, valinit = init_bot, color = '#E48671')

      #Adds a slider for the circle size
      circle_slid_location = plt.axes([0.25, 0.055, 0.65, 0.03])
      circle_slider = Slider(ax = circle_slid_location, label = 'Circle Size:', valmin = 5, valmax = (circle_size + 1*radius), valinit = circle_size, color = '#E48671')

      #Adds a notes section that the user can add notes about their data
      axbox = plt.axes([0.15, 0.02, 0.8, 0.03])
      text = ''
      text_box = TextBox(axbox, 'Notes:', initial = text, textalignment = "center")

      #Make a button that can be clicked if no object is found
      axes_button = plt.axes([0.04, 0.775, 0.92, 0.04])
      close = Button(axes_button, 'Object Not Found', color = '#E48671')

      #Update the image depending on what the user chooses
      def update_button(label):
        '''Updates the list of activated images and updates the image the user can see.'''

        total_data = 0
        for lab in labels:
          if lab == label:
            index = labels.index(lab)
            if default[index] == False:
              default[index] = True
            elif default[index] == True: 
              default[index] = False
        for d in range(len(default)):
          if default == [False, False, False, False]: 
            total_data = real_data[0]*0
          if default[d] == True: 
            total_data = total_data + real_data[d]
          else: 
            pass
        norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
        ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

      #Updates the scaling when the slider is changed
      def update_slider_stretch(val):
        '''Updates the stretch the user can see, based in percentiles'''

        total_data = 0
        for d in range(len(default)):
          if default[d] == True: 
            total_data = total_data + real_data[d]
          else: 
            pass
        norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
        ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)
        
      #Updates the notes added by the user when there is an input
      text_list = [text]
      def submit(expression):
        '''Updates the list of types in the 'Notes' setting'''

        text = expression
        text_list.append(text)

      #Allows the sliders and buttons to be pressed
      check.on_clicked(update_button)
      slider_top.on_changed(update_slider_stretch)
      slider_bottom.on_changed(update_slider_stretch)
      text_box.on_text_change(submit)

      #Display image until it is clicked to find the object
      n = -1
      while True:
        press = plt.waitforbuttonpress()
        text_max = len(text_list) - 1

        #Checks that it was a mouse click
        if press == False:
          n += 3

          #Finds which axes was clicked
          click_axes = str(location[n])
          click_axes = click_axes.split('WCSAxesSubplot', 2)[0]

          #Checks if the image was clicked
          if click_axes == '':
            #Makes a pop-up window with success text
            plt.clf()
            plt.close('all')
            plt.figure(1)
            plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for PanSTARRS! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.grid(linewidth = 0)
            ax = plt.gca()
            ax.xaxis.set_tick_params(labelbottom=False)
            ax.yaxis.set_tick_params(labelleft=False)
            ax.set_xticks([])
            ax.set_yticks([])
            figure2 = plt.gcf()
            figure2.set_size_inches(5.05, 1)
            # figure2.canvas.set_window_title('Successful PanSTARRS Search')
            mng2 = pyplot.get_current_fig_manager()
            mng2.window.resizable(False, False)
            plt.pause(0.1)
            plt.clf()
            plt.close('all')

            #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
            coord = wcs_cropped.pixel_to_world_values(location[n-4], location[n-5])
            distance = []
            for i in range(len(ra_list)):
              distance.append(math.dist(coord, [float(ra_list[i]), float(dec_list[i])]))

            list_location = distance.index(np.min(distance))
            mjd = object_mjd[list_location]
            ps_ra, ps_dec = ra_list[list_location], dec_list[list_location]
            ps_ra_e, ps_dec_e = ra_list_e[list_location], dec_list_e[list_location]
            g, g_e = g_list[list_location], g_list_e[list_location]
            r, r_e = r_list[list_location], r_list_e[list_location]
            i, i_e = i_list[list_location], i_list_e[list_location]
            z, z_e = z_list[list_location], z_list_e[list_location]
            y, y_e = y_list[list_location], y_list_e[list_location]
            return ps_ra, ps_ra_e, ps_dec, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd, 'PanSTARRS Data Release 2', text_list[text_max]
          
          #Checks if the Object not Found button was clicked
          elif click_axes == 'Axes(0.04,0.775;0.92x0.04)':
            ps_ra_e, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
            ps_ra, ps_dec = ra, dec

            #Makes a pop-up window with success text
            plt.clf()
            plt.close('all')
            plt.figure(1)
            plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for PanSTARRS! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.grid(linewidth = 0)
            ax = plt.gca()
            ax.xaxis.set_tick_params(labelbottom=False)
            ax.yaxis.set_tick_params(labelleft=False)
            ax.set_xticks([])
            ax.set_yticks([])
            figure2 = plt.gcf()
            figure2.set_size_inches(5.05, 1)
            # figure2.canvas.set_window_title('Successful PanSTARRS Search')
            mng2 = pyplot.get_current_fig_manager()
            mng2.window.resizable(False, False)
            plt.pause(0.1)
            plt.clf()
            plt.close('all')

            return ps_ra, ps_ra_e, ps_dec, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd, 'PanSTARRS Data Release 2', 'Object Not Found was Pressed'
          
          #Updates the circle size when slider is moved
          elif click_axes == 'Axes(0.25,0.055;0.65x0.03)':
            scatter.remove()
            scatter = ax.scatter(ra_list, dec_list, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')

        #Checks if the window was closed
        elif press is None:
          ps_ra_e, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
          ps_ra, ps_dec = ra, dec
          plt.close('all')
          plt.figure().clear()
          return ps_ra, ps_ra_e, ps_dec, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd, 'PanSTARRS Data Release 2', text_list[text_max]
    
    #If the images were not found returns null values
    elif r_link[0:7] != 'http://':
      ps_ra_e, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
      ps_ra, ps_dec = ra, dec
      return ps_ra, ps_ra_e, ps_dec, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd, 'PanSTARRS Data Release 2', 'Image Not Found'

  #If the images were not found returns null values
  elif len(r_finder_list) == 0: 
    ps_ra_e, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    ps_ra, ps_dec = ra, dec
    return ps_ra, ps_ra_e, ps_dec, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd, 'PanSTARRS Data Release 2', 'Image Not Found'

#Find all the objects in the radius defined by the user
def ps_table(ra, dec, radius): 
  '''Find all the objects in the radius defined by the user'''
  
  blockPrint()

  # #Finds the table for the Panstarrs data around the ra and dec given by the user
  catalog_data = Catalogs.query_region(str(ra) + ' ' + str(dec), radius = ((radius/7200) - 0.000277778), catalog = "Panstarrs", table = "stack")
  return catalog_data

#                        2MASS SEARCH                           #
# ------------------------------------------------------------- #
def twomass_image(ra, dec, radius): 
  ''' First, it gets the images from the 2MASS API from the IRSA Archive and downloads the images. 
  Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
  Third, calls the table function to get all the objects from the 2MASS source catalog. 
  Fourth, makes the window for the user to click the object with all settings. 
  Finally, finds the closest object to the click and records the data. '''
    
  #Makes outline for the window of the plot
  plt.rcParams['toolbar'] = 'None'
  matplotlib.use("TkAgg")
  plt.style.use('Solarize_Light2')
  blockPrint()

  #Finds all the metadata that relates to the ra and dec searched, mostly to find the APIs for the W1, W2, W3, and W4 images
  metadata_2mass_link = 'https://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?POS=' + str(ra) + ',' + str(dec) + '&SIZE=' + str(radius/3600)
  try:
    twomass_metadata = requests.get(metadata_2mass_link)
  except: 
    j, j_e, h, h_e, ks, ks_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    ra_2mass, dec_2mass = ra, dec
    plt.close('all')
    plt.figure().clear()
    return ra_2mass, dec_2mass, j, j_e, h, h_e, ks, ks_e, '2MASS All-Sky Point Source Catalog', '2MASS is Having Maintenance'
  open(str(directory) + '/Output/metadata/TWOMASS_metadata.txt', 'wb').write(twomass_metadata.content)

  #With this metadata it finds the API link for the W1 and W2 images
  J_finder, h_finder, k_finder = 'All-Sky Release Survey J-Band Atlas Image', 'All-Sky Release Survey H-Band Atlas Image', 'All-Sky Release Survey K-Band Atlas Image'
  j_twomass_image_url_list, h_twomass_image_url_list, k_twomass_image_url_list = [], [], []
  with open(str(directory) + '/Output/metadata/TWOMASS_metadata.txt', 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      if line.find(J_finder) != -1:
        j_twomass_image_url_list.append(((lines[lines.index(line) + 1]).split('[')[2]).split(']')[0])
      elif line.find(h_finder) != -1:
        h_twomass_image_url_list.append(((lines[lines.index(line) + 1]).split('[')[2]).split(']')[0])
      elif line.find(k_finder) != -1:
        k_twomass_image_url_list.append(((lines[lines.index(line) + 1]).split('[')[2]).split(']')[0])

  #Finds the properly cropped image from 2MASS
  shape_difference = []
  for i in range(len(j_twomass_image_url_list)):
    j_twomass_image_url= j_twomass_image_url_list[i]

    #Download the W1 and W2 images
    file_allwise_j = download_file(j_twomass_image_url, cache = True)
    data_allwise_j = fits.getdata(file_allwise_j)

    #Gets the headers from the images
    hdu_j = fits.open(file_allwise_j)[0]
    wcs = WCS(hdu_j.header)

    #Make a cutout from the coadd image for the RA and DEC put in
    position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
    size = u.Quantity([radius, radius], u.arcsec)

    #Checks if the ra and dec is in the image wcs
    ra_dec_pixel = wcs.world_to_pixel_values(ra, dec)
    if ((0 <= int(ra_dec_pixel[1]) <= data_allwise_j.shape[0]) == False) or ((0 <= int(ra_dec_pixel[0]) <= data_allwise_j.shape[1]) == False): 
      shape_difference.append(10000000)
    else:
      cutout_j = Cutout2D(data_allwise_j, position, size, fill_value = np.nan, wcs = wcs.celestial)

      #Adds the difference between the height and width's absolute value to a list
      shape_difference.append(abs(cutout_j.shape[0] - cutout_j.shape[1]))

  #Checks if the  image is found
  if (all(x == shape_difference[0] for x in shape_difference) == False) and (len(shape_difference) > 0):

    #Finds the best url link
    shape_best = shape_difference.index(np.min(shape_difference))
    j_twomass_image_url, h_twomass_image_url, k_twomass_image_url = j_twomass_image_url_list[shape_best], h_twomass_image_url_list[shape_best], k_twomass_image_url_list[shape_best]

    #Download the W1 and W2 images for the best url link
    file_allwise_j, file_allwise_h, file_allwise_k = download_file(j_twomass_image_url, cache = True), download_file(h_twomass_image_url, cache = True), download_file(k_twomass_image_url, cache = True)
    data_allwise_j, data_allwise_h, data_allwise_k = fits.getdata(file_allwise_j), fits.getdata(file_allwise_h), fits.getdata(file_allwise_k)

    #Gets the headers from the images for the best url link
    hdu_j, hdu_h, hdu_k = fits.open(file_allwise_j)[0], fits.open(file_allwise_h)[0], fits.open(file_allwise_k)[0]
    wcs = WCS(hdu_j.header)

    #Make a cutout from the coadd image for the RA and DEC put in for the best url link
    cutout_j = Cutout2D(data_allwise_j, position, size, fill_value = np.nan, wcs = wcs.celestial)
    cutout_h = Cutout2D(data_allwise_h, position, size, fill_value = np.nan, wcs = wcs.celestial)
    cutout_k = Cutout2D(data_allwise_k, position, size, fill_value = np.nan, wcs = wcs.celestial)
    wcs_cropped = cutout_j.wcs

    #Find the location of all the object found in AllWISE in the radius choosen by the user 
    location_data = twomass_table(ra, dec, radius)
    object_ra, object_dec = location_data['ra'].tolist(), location_data['dec'].tolist()
    j_list, j_list_e = location_data['j_m'].tolist(), location_data['j_msigcom'].tolist()
    h_list, h_list_e = location_data['h_m'].tolist(), location_data['h_msigcom'].tolist()
    ks_list, ks_list_e = location_data['k_m'].tolist(), location_data['k_msigcom'].tolist()
    enablePrint()

    #Gets the dates of when the images were taken
    date_j, date_h, date_k = hdu_j.header[64], hdu_h.header[64], hdu_k.header[64]

    #Defining a mouse click as an event on the plot
    location = []
    plt.rcParams["figure.figsize"] = [8, 8]
    plt.rcParams["figure.autolayout"] = True
    def mouse_event(event):
      '''Makes a list of the x, y, and axes the mouse click is.'''

      location.append(event.ydata)
      location.append(event.xdata)
      location.append(event.inaxes)
    fig_1 = plt.figure()
    cid = fig_1.canvas.mpl_connect('button_press_event', mouse_event)

    #Sets the WCS coordinates for the plots
    total_data = cutout_j.data + cutout_k.data
    ax = plt.subplot(projection = wcs_cropped)

    #Plots the objects found in the radius
    circle_size = (radius*3)
    scatter = ax.scatter(object_ra, object_dec, transform = ax.get_transform('fk5'), s = circle_size, edgecolor = '#40E842', facecolor = 'none')

    #Normalize the image and plots it
    stretch = 99
    norm1_w1 = ImageNormalize(total_data.data, PercentileInterval(stretch), stretch = SinhStretch())
    ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

    #Makes the figure look pretty
    plt.tick_params(axis = 'x', which = 'both', bottom = False, top = False, labelbottom = False)
    plt.tick_params(axis = 'y', which = 'both', bottom = False, top = False, labelbottom = False)
    fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
    plt.suptitle('2MASS Search', fontsize = 35, y = 0.96, fontfamily = 'Times New Roman')
    ax.set_title('Dates: \n'
              + 'J Date: ' + str(date_j) + ' (YYMMDD)  ' + '  H Date: ' + str(date_h) + ' (YYMMDD)\n'
              + 'K Date: ' + str(date_k) + ' (YYMMDD)  \n', fontdict = fontdict_1, y = 1.04)
    plt.grid(linewidth = 0)
    shape_x = cutout_j.shape[1]
    shape_y = cutout_j.shape[0]
    plt.xlim(shape_x - radius, radius + (shape_x - radius))
    plt.ylim(shape_y - radius, radius + (shape_y - radius))
    
    figure = plt.gcf()
    if platform != 'win32':
      figure.set_size_inches(4.75, 6.95)
    elif platform == 'win32':
      figure.set_size_inches(4.75, 7.05)
    # figure.canvas.set_window_title('2MASS Search')
    mng = pyplot.get_current_fig_manager()
    mng.window.resizable(False, False)
    
    #Makes a cursor to aid the User
    cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
    annotation = ax.annotate('', xy=(0, 0), xytext=(20, -20), arrowprops=dict(arrowstyle='wedge'), fontsize=12, color='red')
    annotation.set_visible(False)

    #Make checkbuttons with all of the different image bands
    rax = plt.axes([0.045, 0.4, 0.105, 0.12])
    labels = ['J', 'H', 'K']
    real_data = [cutout_j.data, cutout_h.data, cutout_k.data]
    default = [True, False, True]
    check = CheckButtons(rax, labels, default)

    #Adds a slider for the scaling of the image
    freq_bottom = plt.axes([0.25, 0.1, 0.65, 0.03])
    slider_bottom = Slider(ax = freq_bottom, label = 'Stetch:', valmin = 80, valmax = 100, valinit = stretch, color = '#E48671')

    #Adds a slider for the circle size
    circle_slid_location = plt.axes([0.25, 0.068, 0.65, 0.03])
    circle_slider = Slider(ax = circle_slid_location, label = 'Circle Size:', valmin = (circle_size - 2.5*radius), valmax = (circle_size + 1*radius), valinit = circle_size, color = '#E48671')

    #Adds a notes section that the user can add notes about their data
    axbox = plt.axes([0.15, 0.025, 0.8, 0.03])
    text = ''
    text_box = TextBox(axbox, 'Notes:', initial = text, textalignment = "center")

    #Make a button that can be clicked if no object is found
    axes_button = plt.axes([0.04, 0.78, 0.92, 0.04])
    close = Button(axes_button, 'Object Not Found', color = '#E48671')

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
      for d in range(len(default)):
        if default != [True, True, True]: 
          total_data = real_data[1]*0
        if default[d] == True: 
          total_data = total_data + real_data[d]
        else: 
          pass
      norm1_w1 = ImageNormalize(total_data.data, PercentileInterval(slider_bottom.val), stretch = SinhStretch())
      ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

    #Updates the scaling when the slider is changed
    def update_slider_stretch(val):
      '''Updates the stretch the user can see, based in percentiles'''

      total_data = 0
      for d in range(len(default)):
        if default[d] == True: 
          total_data = total_data + real_data[d]
        else: 
          pass
      norm1_w1 = ImageNormalize(total_data.data, PercentileInterval(slider_bottom.val), stretch = SinhStretch())
      ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

    #Updates the notes added by the user when there is an input
    text_list = [text]
    def submit(expression):
      '''Updates the list of types in the 'Notes' setting'''

      text = expression
      text_list.append(text)   

    #Allows the sliders and buttons to be pressed
    check.on_clicked(update_button)
    slider_bottom.on_changed(update_slider_stretch)
    text_box.on_text_change(submit)

    #Display image until it is clicked to find the object
    n = -1
    while True:
      press = plt.waitforbuttonpress()
      text_max = len(text_list) - 1

      #Checks that it was a mouse click
      if press == False:
        n += 3

        #Finds which axes was clicked
        click_axes = str(location[n])
        click_axes = click_axes.split('WCSAxesSubplot', 2)[0]

        #Checks if the image was clicked
        if click_axes == '':
          #Makes a pop-up window with success text
          plt.clf()
          plt.close('all')
          plt.figure(1)
          plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for 2MASS! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
          plt.xlim(0, 1)
          plt.ylim(0, 1)
          plt.grid(linewidth = 0)
          ax = plt.gca()
          ax.xaxis.set_tick_params(labelbottom=False)
          ax.yaxis.set_tick_params(labelleft=False)
          ax.set_xticks([])
          ax.set_yticks([])
          figure2 = plt.gcf()
          figure2.set_size_inches(4.75, 1)
          # figure2.canvas.set_window_title('Successful 2MASS Search')
          mng2 = pyplot.get_current_fig_manager()
          mng2.window.resizable(False, False)
          plt.pause(0.1)
          plt.clf()
          plt.close('all')

          #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
          coord = wcs_cropped.pixel_to_world_values(location[n-4],location[n-5])
          distance = []
          for i in range(len(object_ra)):
            distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))

          list_location = distance.index(np.min(distance))
          ra_2mass, dec_2mass = object_ra[list_location], object_dec[list_location]
          j, j_e = j_list[list_location], j_list_e[list_location]
          h, h_e = h_list[list_location], h_list_e[list_location]
          ks, ks_e = ks_list[list_location], ks_list_e[list_location]
          return ra_2mass, dec_2mass, j, j_e, h, h_e, ks, ks_e, '2MASS All-Sky Point Source Catalog', text_list[text_max]
        
        #Checks if the Object not Found button was clicked
        elif click_axes == 'Axes(0.04,0.78;0.92x0.04)':
          j, j_e, h, h_e, ks, ks_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
          ra_2mass, dec_2mass = ra, dec

          #Makes a pop-up window with success text
          plt.clf()
          plt.close('all')
          plt.figure(1)
          plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for 2MASS! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
          plt.xlim(0, 1)
          plt.ylim(0, 1)
          plt.grid(linewidth = 0)
          ax = plt.gca()
          ax.xaxis.set_tick_params(labelbottom=False)
          ax.yaxis.set_tick_params(labelleft=False)
          ax.set_xticks([])
          ax.set_yticks([])
          figure2 = plt.gcf()
          figure2.set_size_inches(4.75, 1)
          # figure2.canvas.set_window_title('Successful 2MASS Search')
          mng2 = pyplot.get_current_fig_manager()
          mng2.window.resizable(False, False)
          plt.pause(0.1)
          plt.clf()
          plt.close('all')
          
          return ra_2mass, dec_2mass, j, j_e, h, h_e, ks, ks_e, '2MASS All-Sky Point Source Catalog', 'Object Not Found was Pressed'
        
        #Adds the functionality of the circle slider bar
        elif click_axes == 'Axes(0.25,0.068;0.65x0.03)':
          scatter.remove()
          scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')

      #Checks if the window was closed
      elif press is None:
        j, j_e, h, h_e, ks, ks_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        ra_2mass, dec_2mass = ra, dec
        plt.close('all')
        plt.figure().clear()
        return ra_2mass, dec_2mass, j, j_e, h, h_e, ks, ks_e, '2MASS All-Sky Point Source Catalog', text_list[text_max]
  else: 
      j, j_e, h, h_e, ks, ks_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
      ra_2mass, dec_2mass = ra, dec
      text_list = 'Image Not Found'
      plt.close('all')
      plt.figure().clear()
      return ra_2mass, dec_2mass, j, j_e, h, h_e, ks, ks_e, '2MASS All-Sky Point Source Catalog', text_list

def twomass_table(ra, dec, radius): 
  '''Find all the objects in the radius defined by the user'''
  
  blockPrint()

  #Uses astroquery to find all objects in the radius
  location_data = Irsa.query_region(coord.SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='fk5'), catalog='fp_psc', spatial='Cone', radius=(radius - 1) * u.arcsec)
  return location_data

#                        WFCAM SEARCH                           #
# ------------------------------------------------------------- #
def ukidss_image(ra, dec, radius): 
    ''' First, it gets the images from the WFCAM API from the WFCAM Archive and downloads the images. 
    Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
    Third, calls the table function to get all the objects from the WFCAM source catalog. 
    Fourth, makes the window for the user to click the object with all settings. 
    Finally, finds the closest object to the click and records the data. '''
      
    #Makes outline for the window of the plot
    plt.rcParams['toolbar'] = 'None'
    matplotlib.use("TkAgg")
    plt.style.use('Solarize_Light2')
    blockPrint()
    
    try:
      url_J = Ukidss.get_image_list(SkyCoord(1, 1, unit = (u.deg, u.deg), frame = 'fk5'), image_width = (150) * u.arcsec, image_height = (150) * u.arcsec, waveband = 'J', database = 'UHSDR1', programme_id = 'UHSDR1')
    except: 
      ra_wfcam_e, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
      ra_wfcam = ra
      dec_wfcam = dec
      return ra_wfcam, ra_wfcam_e, dec_wfcam, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch, 'WFCAM', 'WFCAM is Having Maintenance'
    
    #Obtains all of the urls in J, H, and K from UKIDSS
    database_list = ['UKIDSSDR11PLUS', 'UHSDR1']
    catalo = ['LAS', 'UHSDR1']
    for data in range(len(database_list)):
        url_J = Ukidss.get_image_list(SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'), image_width = (radius) * u.arcsec, image_height = (radius) * u.arcsec, waveband = 'J', database = database_list[data], programme_id = catalo[data])
        url_H = Ukidss.get_image_list(SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'), image_width = (radius) * u.arcsec, image_height = (radius) * u.arcsec, waveband = 'H', database = database_list[data], programme_id = catalo[data])
        url_K = Ukidss.get_image_list(SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'), image_width = (radius) * u.arcsec, image_height = (radius) * u.arcsec, waveband = 'K', database = database_list[data], programme_id = catalo[data])
        url_Y = Ukidss.get_image_list(SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'), image_width = (radius) * u.arcsec, image_height = (radius) * u.arcsec, waveband = 'Y', database = database_list[data], programme_id = catalo[data])
        if len(url_J) > 0: 
            data = database_list[data]
            break

    #Checking to see if the images exist
    if len(url_J) > 0: 

        #Calls the UKIDSS Table
        table = ukidss_table(ra, dec, radius)
        
        #Defines each variables depending on if the image was found in UKIDSS DR11 or UHS DR1
        if data == 'UHSDR1': 

            #Downloading the fits images
            file_ukidss_J = download_file(url_J[0], cache = True) 
            data_ukidss_J = fits.getdata(file_ukidss_J)

            #Obtains the headers from the images
            hdu_j = fits.open(file_ukidss_J)[1]
            wcs_j = WCS(hdu_j.header)

            #Obtains the dates for each image
            date_j = str(hdu_j.header['HISTORY']).split(' ', 2)[0]

            #Gets the columns from the table
            object_ra, object_dec = table['ra'].tolist(), table['dec'].tolist()
            J_list, J_list_e = table['jAperMag3'].tolist(), table['jAperMag3Err'].tolist()
        
        else:
          if len(url_J) > 0 and len(url_H) > 0 and len(url_K) > 0 and len(url_Y) > 0: 
            #Downloading the fits images
            file_ukidss_J, file_ukidss_H, file_ukidss_K, file_ukidss_Y = download_file(url_J[0], cache=True), download_file(url_H[0], cache=True), download_file(url_K[0], cache=True), download_file(url_Y[0], cache=True)
            data_ukidss_J, data_ukidss_H, data_ukidss_K, data_ukidss_Y = fits.getdata(file_ukidss_J), fits.getdata(file_ukidss_H), fits.getdata(file_ukidss_K), fits.getdata(file_ukidss_Y)

            #Obtains the headers from the images
            hdu_j, hdu_h, hdu_k, hdu_y = fits.open(file_ukidss_J)[1], fits.open(file_ukidss_H)[1], fits.open(file_ukidss_K)[1], fits.open(file_ukidss_Y)[1]
            wcs_j, wcs_h, wcs_k, wcs_y = WCS(hdu_j.header),  WCS(hdu_h.header), WCS(hdu_k.header), WCS(hdu_y.header)

            #Make a cutout from the coadd image for the RA and DEC put in
            position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5')
            size = u.Quantity([radius, radius], u.arcsec)
            cutout_j = Cutout2D(data_ukidss_J, position, size, wcs = wcs_j.celestial)
            cutout_h = Cutout2D(data_ukidss_H, position, size, wcs = wcs_h.celestial)
            cutout_k = Cutout2D(data_ukidss_K, position, size, wcs = wcs_k.celestial)
            cutout_y = Cutout2D(data_ukidss_Y, position, size, wcs = wcs_y.celestial)
            wcs_cropped = cutout_j.wcs

            #Resizes the H, K, and Y bands to that of J
            j_shape = (cutout_j.data).shape
            if j_shape[0] != j_shape[1]:
                cutout_j = cutout_j.data[0:min(j_shape) - 1, 0:min(j_shape) - 1]
                j_shape = (cutout_j.data).shape
            h_reshape, k_reshape, y_reshape = cv2.resize(cutout_h.data, dsize = j_shape, interpolation = cv2.INTER_NEAREST), cv2.resize(cutout_k.data, dsize = j_shape, interpolation = cv2.INTER_NEAREST), cv2.resize(cutout_y.data, dsize = j_shape, interpolation = cv2.INTER_NEAREST)
            total_data = cutout_j.data + y_reshape

            #Obtains the dates for each image
            date_y, date_j, date_h, date_k = str(hdu_y.header['HISTORY']).split(' ', 2)[0], str(hdu_j.header['HISTORY']).split(' ', 2)[0], str(hdu_h.header['HISTORY']).split(' ', 6)[2], str(hdu_k.header['HISTORY']).split(' ', 6)[2]

            #Gets the columns from the table
            object_epoch = table['epoch'].tolist()
            object_ra, object_dec = table['ra'].tolist(), table['dec'].tolist()
            object_ra_e, object_dec_e = table['sigRa'].tolist(), table['sigDec'].tolist()
            Y_list, Y_list_e = table['yAperMag3'].tolist(), table['yAperMag3Err'].tolist()
            J_list, J_list_e = table['jAperMag3'].tolist(), table['jAperMag3Err'].tolist()
            H_list, H_list_e = table['hAperMag3'].tolist(), table['hAperMag3Err'].tolist()
            K_list, K_list_e = table['kAperMag3'].tolist(), table['kAperMag3Err'].tolist()
            pmra_list, pmra_list_e = table['muRa'].tolist(), table['sigMuRa'].tolist()
            pmdec_list, pmdec_list_e = table['muDec'].tolist(), table['sigMuDec'].tolist()
            enablePrint()
          else: 
            J_mag, ra_wfcam, J_e, dec_wfcam = np.nan, np.nan, np.nan, np.nan
            text = 'Image Not Found'
            ra_wfcam_e, dec_wfcam_e, Y_mag, Y_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
            return ra_wfcam, ra_wfcam_e, dec_wfcam, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch, data, text

        #Defining a mouse click as an event on the plot
        location = []
        plt.rcParams["figure.figsize"] = [8, 8]
        plt.rcParams["figure.autolayout"] = True
        def mouse_event(event):
            '''Makes a list of the x, y, and axes the mouse click is'''

            location.append(event.ydata)
            location.append(event.xdata)
            location.append(event.inaxes)
        fig_1 = plt.figure()
        cid = fig_1.canvas.mpl_connect('button_press_event', mouse_event)

        #Sets the plot depending on if it was found in UHS or UKIDSS
        if data == 'UHSDR1':

            #Makes the subplot for the plot
            ax = plt.subplot()

            #Gets the cutout for the fits iamge
            wcs = WCS(hdu_j.header)
            position = wcs.world_to_pixel_values(ra, dec)
            size = u.Quantity([(radius * 2.5), (radius * 2.5)], u.pixel)
            cutout_UHS = Cutout2D(data_ukidss_J, position, size)

            #Finds the camera orientation
            cam_type = hdu_j.header['CAMNUM']
            print(cam_type)

            #Obtains the shape of the cutout and sets the circle size for the scatter plot
            shape = min(cutout_UHS.shape)
            circle_size = (radius*3)

            #Converts the ra
            ra_dec_pixel = wcs.world_to_pixel_values(object_ra, object_dec)
            
            # Finds the radius size in pixels
            pixel_radius = 2.5*radius
            if cam_type == 1: 
                #Makes the dec negative
                total_data = np.rot90(cutout_UHS.data, 3)
                if total_data.shape[1] < total_data.shape[0]:
                  minus_dec = [(-x + total_data.shape[1]) for x in ra_dec_pixel[1]]
                else:
                  minus_dec = [(-x + pixel_radius) for x in ra_dec_pixel[1]]
                
                #Plots the correctly orientated image
                scatter = ax.scatter(minus_dec, ra_dec_pixel[0], s = circle_size, edgecolor = '#40E842', facecolor = 'none')
                x_lower, x_upper = (total_data.shape[1] - pixel_radius) + pixel_radius, -(pixel_radius - total_data.shape[1])
                y_lower, y_upper = total_data.shape[0] - pixel_radius, total_data.shape[0]

            elif cam_type == 2:
                #Plots the correctly orientated image
                scatter = ax.scatter(ra_dec_pixel[0], ra_dec_pixel[1], s = circle_size, edgecolor = '#40E842', facecolor = 'none')
                total_data = cutout_UHS.data
                x_lower, x_upper = 0, total_data.shape[0] + (pixel_radius - total_data.shape[0])
                y_lower, y_upper = total_data.shape[1] + (total_data.shape[1] - pixel_radius), (total_data.shape[1] - pixel_radius)

            elif cam_type == 3:
                #Makes the ra negative
                total_data = np.rot90(cutout_UHS.data)
                if total_data.shape[0] < total_data.shape[1]:
                  minus_ra = [(-x + total_data.shape[0]) for x in ra_dec_pixel[0]]
                else:
                  minus_ra = [(-x + pixel_radius) for x in ra_dec_pixel[0]]

                #Plots the correctly orientated image
                scatter = ax.scatter(ra_dec_pixel[1], minus_ra, s = circle_size, edgecolor = '#40E842', facecolor = 'none')
                x_lower, x_upper = total_data.shape[1], -(pixel_radius - total_data.shape[1])
                y_lower, y_upper = total_data.shape[0] - pixel_radius, total_data.shape[0]
                
            elif cam_type == 4: 
                #Plots the correctly orientated image
                scatter = ax.scatter(ra_dec_pixel[0], ra_dec_pixel[1], s = circle_size, edgecolor = '#40E842', facecolor = 'none')
                total_data = cutout_UHS.data
                x_lower, x_upper = pixel_radius, 0
                y_lower, y_upper = 0, pixel_radius

            # Normalize the image and plots it
            init_top, init_bot = 95, 45
            norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(cutout_UHS.data, init_bot), vmax = np.nanpercentile(cutout_UHS.data, init_top))
            print(total_data.shape)
            ax.imshow(total_data, cmap = 'Greys', norm = norm1_w1, origin = 'upper')
            plt.xlim(x_lower, x_upper)
            plt.ylim(y_lower, y_upper)

            # Formats the window correctly
            fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
            figure = plt.gcf()
            ax.set_title('Dates: \n' 'J Date: ' + str(date_j) + ' (YYYYMMDD) \n', fontdict=fontdict_1, y=1.05)
            if platform != 'win32':
                figure.set_size_inches(4.75, 7.15)
            elif platform == 'win32': 
                figure.set_size_inches(4.75, 7.15)
        else:
            #Sets the WCS coordinates for the plots
            ax = plt.subplot(projection = wcs_cropped)

            #Plots the objects found in the radius
            circle_size = (radius*3)
            scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

            cam_type = ''

            #Normalize the image and plots it
            init_top, init_bot = 95, 45
            norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, init_bot), vmax = np.nanpercentile(total_data.data, init_top))
            ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

            #Formats the window correctly
            fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
            figure = plt.gcf()
            ax.set_title('Dates: \n'
                    + 'Y Date: ' + str(date_y) + ' (YYYYMMDD)   ' + 'J Date: ' + str(date_j) + ' (YYYYMMDD) \n'
                    + 'H Date: ' + str(date_h) + ' (Y/M/D)   ' + 'K Date: ' + str(date_k) + ' (Y/M/D) \n', fontdict = fontdict_1, y = 1.05)
            if platform != 'win32':
                figure.set_size_inches(4.75, 7.25)
            elif platform == 'win32':
                figure.set_size_inches(4.75, 7.25)
            plt.xlim(max(total_data.shape), 0)
            plt.ylim(0, max(total_data.shape))
            plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)

        #Finishes formatting the problem correctly
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.yaxis.set_tick_params(labelleft=False)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.suptitle('WFCAM Search', fontsize = 35, y = 0.96, fontfamily = 'Times New Roman')
        plt.grid(linewidth = 0)
        # figure.canvas.set_window_title('WFCAM Search')
        mng = pyplot.get_current_fig_manager()
        mng.window.resizable(False, False)
        
        #Makes a cursor to aid the User
        cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
        annotation = ax.annotate('', xy=(0, 0), xytext=(20, -20), arrowprops=dict(arrowstyle='wedge'), fontsize=12, color='red')
        annotation.set_visible(False)

        #Make checkbuttons with all of the different image bands
        default = [True]
        real_data = [total_data]
        if data != 'UHSDR1':
            shape = max(cutout_j.shape)
            rax = plt.axes([0.045, 0.4, 0.115, 0.1])
            labels = ['Y', 'J', 'H', 'K']
            real_data = [y_reshape, cutout_j.data, h_reshape, k_reshape]
            default = [True, True, False, False]
            check = CheckButtons(rax, labels, default)

        #Adds a slider for the scaling of the image
        freq_top = plt.axes([0.25, 0.12, 0.65, 0.03])
        slider_top = Slider(ax = freq_top, label = 'Top Stetch:', valmin = 50, valmax = 100, valinit = init_top, color = '#E48671')
        freq_bottom = plt.axes([0.25, 0.087, 0.65, 0.03])
        slider_bottom = Slider(ax = freq_bottom, label = 'Bottom Stetch:', valmin = 0, valmax = 50, valinit = init_bot, color = '#E48671')

        #Adds a slider for the circle size
        circle_slid_location = plt.axes([0.25, 0.055, 0.65, 0.03])
        circle_slider = Slider(ax = circle_slid_location, label = 'Circle Size:', valmin = (circle_size - 2.5*radius), valmax = (circle_size + 1*radius), valinit = circle_size, color = '#E48671')

        #Adds a notes section that the user can add notes about their data
        axbox = plt.axes([0.15, 0.02, 0.8, 0.03])
        text = ''
        text_box = TextBox(axbox, 'Notes:', initial = text, textalignment="center")

        #Make a button that can be clicked if no object is found
        axes_button = plt.axes([0.04, 0.775, 0.92, 0.04])
        close = Button(axes_button, 'Object Not Found', color = '#E48671')

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
            for d in range(len(default)):
                if default == [False, False, False, False]: 
                    total_data = real_data[0]*0
                if default[d] == True: 
                    total_data = total_data + real_data[d]
                else: 
                    pass
            norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
            ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

        #Updates the scaling when the slider is changed
        def update_slider_stretch(val):
            '''Updates the stretch the user can see, based in percentiles'''

            total_data = 0
            for d in range(len(default)):
                if default[d] == True: 
                    total_data = total_data + real_data[d]
                else: 
                    pass
            norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
            if data == 'UHSDR1':
                ax.imshow(total_data, cmap = 'Greys', norm = norm1_w1)
            else:
                ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

        #Updates the notes added by the user when there is an input
        text_list = [text]
        def submit(expression):
            '''Updates the list of types in the 'Notes' setting'''

            text = expression
            text_list.append(text)

        #Allows the sliders and buttons to be pressed
        if data != 'UHSDR1':
            check.on_clicked(update_button)
        slider_top.on_changed(update_slider_stretch)
        slider_bottom.on_changed(update_slider_stretch)
        text_box.on_text_change(submit)

        #Display image until it is clicked to find the object
        n = -1
        while True:
            press = plt.waitforbuttonpress()
            text_max = len(text_list) - 1
            enablePrint()

            #Checks that it was a mouse click
            if press == False:
                n += 3

                #Finds which axes was clicked
                click_axes = str(location[n])
                if data != 'UHSDR1':
                    click_axes = click_axes.split('WCSAxesSubplot', 2)[0]
                else: 
                    click_axes = click_axes.split('AxesSubplot', 2)[0]

                #Checks if the image was clicked
                if click_axes == '':
                    #Makes a pop-up window with success text
                    plt.clf()
                    plt.close('all')
                    plt.figure(1)
                    plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for WFCAM! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
                    plt.xlim(0, 1)
                    plt.ylim(0, 1)
                    plt.grid(linewidth = 0)
                    ax = plt.gca()
                    ax.xaxis.set_tick_params(labelbottom=False)
                    ax.yaxis.set_tick_params(labelleft=False)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    figure2 = plt.gcf()
                    figure2.set_size_inches(4.75, 1)
                    # figure2.canvas.set_window_title('Successful WFCAM Search')
                    mng2 = pyplot.get_current_fig_manager()
                    mng2.window.resizable(False, False)
                    plt.pause(0.1)
                    plt.clf()
                    plt.close('all')

                    #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
                    if data != 'UHSDR1':
                        coord = wcs_cropped.pixel_to_world_values(location[n-4],location[n-5])
                        distance = []
                        for i in range(len(object_ra)):
                            distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))
                        list_location = distance.index(np.min(distance))
                    else: 
                        if (cam_type == 2) or (cam_type == 4):
                            coord = wcs.pixel_to_world_values(location[n-4],location[n-5])
                            distance = []
                            for i in range(len(object_ra)):
                                distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))
                            list_location = distance.index(np.min(distance))
                        elif cam_type == 3:
                            coord = [location[n - 4], location[n - 5]]
                            distance = []
                            for i in range(len(minus_ra)):
                                distance.append(math.dist(coord, [float(ra_dec_pixel[1][i]), float(minus_ra[i])]))
                            list_location = distance.index(np.min(distance))
                        elif cam_type == 1: 
                            coord = [location[n - 4], location[n - 5]]
                            distance = []
                            for i in range(len(minus_dec)):
                                distance.append(math.dist(coord, [float(minus_dec[i]), float(ra_dec_pixel[0][i])]))
                            list_location = distance.index(np.min(distance))

                    if data != 'UHSDR1':
                        epoch = object_epoch[list_location]
                        ra_wfcam, dec_wfcam = object_ra[list_location], object_dec[list_location]
                        ra_wfcam_e, dec_wfcam_e = object_ra_e[list_location], object_dec_e[list_location]
                        Y_mag, Y_e = Y_list[list_location], Y_list_e[list_location]
                        J_mag, J_e = J_list[list_location], J_list_e[list_location]
                        H_mag, H_e = H_list[list_location], H_list_e[list_location]
                        K_mag, K_e = K_list[list_location], K_list_e[list_location]
                        pmra, pmra_e = pmra_list[list_location], pmra_list_e[list_location]
                        pmdec, pmdec_e = pmdec_list[list_location], pmdec_list_e[list_location]
                        return ra_wfcam, ra_wfcam_e, dec_wfcam, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch, data, text_list[text_max] 
                    else: 
                        ra_wfcam, dec_wfcam = object_ra[list_location], object_dec[list_location]
                        J_mag, J_e = J_list[list_location], J_list_e[list_location]
                        ra_wfcam_e, dec_wfcam_e, Y_mag, Y_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                        return ra_wfcam, ra_wfcam_e, dec_wfcam, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch, data, text_list[text_max] 
                
                #Checks if the "Object Not Found" button was clicked
                elif click_axes == 'Axes(0.04,0.775;0.92x0.04)':
                    ra_wfcam_e, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                    ra_wfcam = ra
                    dec_wfcam = dec
                    
                    #Makes a pop-up window with success text
                    plt.clf()
                    plt.close('all')
                    plt.figure(1)
                    plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for WFCAM! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
                    plt.xlim(0, 1)
                    plt.ylim(0, 1)
                    plt.grid(linewidth = 0)
                    ax = plt.gca()
                    ax.xaxis.set_tick_params(labelbottom=False)
                    ax.yaxis.set_tick_params(labelleft=False)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    figure2 = plt.gcf()
                    figure2.set_size_inches(4.75, 1)
                    # figure2.canvas.set_window_title('Successful WFCAM Search')
                    mng2 = pyplot.get_current_fig_manager()
                    mng2.window.resizable(False, False)
                    plt.pause(0.1)
                    plt.clf()
                    plt.close('all')

                    return ra_wfcam, ra_wfcam_e, dec_wfcam, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch, data, 'Object Not Found was Pressed'
                
                #Updates the circle size when slider is moved
                elif click_axes == 'Axes(0.25,0.055;0.65x0.03)':
                    scatter.remove()
                    if data != 'UHSDR1':
                        scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')
                    else:
                        if (cam_type == 2) or (cam_type == 4):
                            scatter = ax.scatter(ra_dec_pixel[0], ra_dec_pixel[1], s = circle_slider.val, edgecolor = '#40E842', facecolor = 'none')
                        elif cam_type == 3:
                            scatter = ax.scatter(ra_dec_pixel[1], minus_ra, s = circle_slider.val, edgecolor = '#40E842', facecolor = 'none')
                        elif cam_type == 1:
                            scatter = ax.scatter(minus_dec, ra_dec_pixel[0], s = circle_slider.val, edgecolor = '#40E842', facecolor = 'none')

            #Checks if the window was closed
            elif press is None:
                ra_wfcam_e, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                ra_wfcam = ra
                dec_wfcam = dec
                plt.close('all')
                plt.figure().clear()
                return ra_wfcam, ra_wfcam_e, dec_wfcam, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch, data, text_list[text_max] 

    #If the image is not found return null values
    else: 
        ra_wfcam_e, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        ra_wfcam = ra
        dec_wfcam = dec
        text_list = 'Image Not Found'
        return ra_wfcam, ra_wfcam_e, dec_wfcam, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch, 'WFCAM', text_list
    
def ukidss_table(ra, dec, radius): 
    '''Find all the objects in the radius defined by the user'''

    blockPrint()

    #Find the table of all the objects found in UKIDSS in the radius choosen by the user
    program_list = Ukidss.list_catalogs()
    for prom in program_list:
        if prom != 'UHS':
            table = Ukidss.query_region(
                    SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'),
                    radius = (radius/2) * u.arcsec, 
                    programme_id = prom)
            if len(table) > 0: 
                return table
        elif prom == 'UHS':
            table = Ukidss.query_region(
                    SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'),
                    radius = ((radius/2) - 1) * u.arcsec, 
                    programme_id = prom, database = 'UHSDR1')
            if len(table) > 0: 
                return table
        
#                          VSA SEARCH                           #
# ------------------------------------------------------------- #
def vsa_image(ra, dec, radius): 
  ''' First, it gets the images from the VSA API from the VSA Archive and downloads the images. 
  Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
  Third, calls the table function to get all the objects from the VSA source catalog. 
  Fourth, makes the window for the user to click the object with all settings. 
  Finally, finds the closest object to the click and records the data. '''
  
  #Makes outline for the window of the plot
  plt.rcParams['toolbar'] = 'None'
  matplotlib.use("TkAgg")
  plt.style.use('Solarize_Light2')
  blockPrint()
  
  try:
    url_J = Vsa.get_image_list(SkyCoord(1, 1, unit = (u.deg, u.deg), frame = 'fk5'), image_width = 150 * u.arcsec, waveband = 'J', database = 'VHSDR6')
  except: 
    ymjd, jmjd, hmjd, ksmjd, y, y_e, j, j_e, h, h_e, ks, ks_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    ra_vsa, dec_vsa = ra, dec
    return ra_vsa, dec_vsa, y, y_e, j, j_e, h, h_e, ks, ks_e, ymjd, jmjd, hmjd, ksmjd, 'VSA', 'VSA is Having Maintenance'

  #Obtains all of the urls in J, H, and K from VSA
  database_list = ['VHSDR6', 'VVVDR5', 'VMCDR6', 'VIDEODR6', 'ULTRAVISTADR4']
  for img in database_list:
    url_J, url_H, url_Ks = [Vsa.get_image_list(
      SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'),
                image_width = radius * u.arcsec,
                waveband = band,
                database = img)
      for band in ['J', 'H', 'Ks']]
    if len(url_J) == 0:
      pass
    else: 
      break
  
  #Checking to see if the images exist
  if len(url_J) > 0 and len(url_H) > 0 and len(url_Ks) > 0: 

    #Downloading the fits images
    file_vsa_J, file_vsa_H, file_vsa_Ks = download_file(url_J[0], cache=True), download_file(url_H[0], cache=True), download_file(url_Ks[0], cache=True)
    data_vsa_J, data_vsa_H, data_vsa_Ks = fits.getdata(file_vsa_J), fits.getdata(file_vsa_H), fits.getdata(file_vsa_Ks)

    table = vsa_table(ra, dec, radius)
    
    #Obtains the photometry and astrometry from the catalog
    object_mjd_y, object_mjd_j = table['yMjd'].tolist(), table['jMjd'].tolist()
    object_mjd_h, object_mjd_ks = table['hMjd'].tolist(), table['ksMjd'].tolist()
    object_ra, object_dec = table['ra'].tolist(), table['dec'].tolist()
    Y_list, Y_list_e = table['yAperMag3'].tolist(), table['yAperMag3Err'].tolist()
    J_list, J_list_e = table['jAperMag3'].tolist(), table['jAperMag3Err'].tolist()
    H_list, H_list_e = table['hAperMag3'].tolist(), table['hAperMag3Err'].tolist()
    Ks_list, Ks_list_e = table['ksAperMag3'].tolist(), table['ksAperMag3Err'].tolist()

    #Reads in the header from the image
    hdu_j, hdu_h, hdu_k = fits.open(file_vsa_J)[1], fits.open(file_vsa_H)[1], fits.open(file_vsa_Ks)[1]
    wcs1_j = WCS(hdu_j.header)

    #Make a cutout from the coadd image for the RA and DEC put in
    position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5')
    size = u.Quantity([radius, radius], u.arcsec)
    cutout_j = Cutout2D(data_vsa_J, position, size, wcs = wcs1_j.celestial)
    cutout_h = Cutout2D(data_vsa_H, position, size, wcs = wcs1_j.celestial)
    cutout_ks = Cutout2D(data_vsa_Ks, position, size, wcs = wcs1_j.celestial)
    wcs_cropped = cutout_j.wcs
    enablePrint()

    #Finds the dates the images were taken
    date_j, date_h, date_k = hdu_j.header['DATE-OBS'].split('T', 2)[0], hdu_h.header['DATE-OBS'].split('T', 2)[0], hdu_k.header['DATE-OBS'].split('T', 2)[0]

    #Defining a mouse click as an event on the plot
    location = []
    plt.rcParams["figure.figsize"] = [8, 8]
    plt.rcParams["figure.autolayout"] = True
    def mouse_event(event):
      '''Makes a list of the x, y, and axes the mouse click is.'''

      location.append(event.ydata)
      location.append(event.xdata)
      location.append(event.inaxes)
    fig_1 = plt.figure()
    cid = fig_1.canvas.mpl_connect('button_press_event', mouse_event)
    
    #Sets the WCS coordinates for the plots
    ks_shape = cutout_ks.shape
    j_reshape, h_reshape = cv2.resize(cutout_j.data, (ks_shape[1], ks_shape[0]), interpolation = cv2.INTER_NEAREST), cv2.resize(cutout_h.data, (ks_shape[1], ks_shape[0]), interpolation = cv2.INTER_NEAREST)
    total_data = j_reshape + cutout_ks.data
    ax = plt.subplot(projection = wcs_cropped)

    #Plots the objects found in the radius
    circle_size = (radius*3)
    scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

    #Normalize the image and plots it
    init_top = 95
    init_bot = 45
    norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, init_bot), vmax = np.nanpercentile(total_data.data, init_top))
    ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

    #Makes the figure look pretty
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
    fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
    plt.suptitle('VISTA Search', fontsize = 35, y = 0.96, fontfamily = 'Times New Roman')
    ax.set_title('Dates: \n'
               + 'J Date: ' + str(date_j) + ' (Y-M-D)  ' + '  H Date: ' + str(date_h) + ' (Y-M-D)\n'
               + 'K Date: ' + str(date_k) + ' (Y-M-D) \n', fontdict = fontdict_1, y = 1.05)
    plt.grid(linewidth = 0)
    shape = min(cutout_ks.shape)
    plt.xlim(0, shape)
    plt.ylim(shape, 0)
    figure = plt.gcf()
    if platform != 'win32': 
      figure.set_size_inches(4.75, 7.25)
    elif platform == 'win32': 
      figure.set_size_inches(4.75, 7.25)
    # figure.canvas.set_window_title('VISTA Search')
    mng = pyplot.get_current_fig_manager()
    mng.window.resizable(False, False)
    
    #Makes a cursor to aid the User
    cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
    annotation = ax.annotate('', xy=(0, 0), xytext=(20, -20), arrowprops=dict(arrowstyle='wedge'), fontsize=12, color='red')
    annotation.set_visible(False)

    #Make checkbuttons with all of the different image bands
    rax = plt.axes([0.045, 0.4, 0.105, 0.12])
    labels = ['J', 'H', 'Ks']
    real_data = [j_reshape, h_reshape, cutout_ks.data]
    default = [True, False, True]
    check = CheckButtons(rax, labels, default)

    #Adds a slider for the scaling of the image
    freq_top = plt.axes([0.25, 0.12, 0.65, 0.03])
    slider_top = Slider(ax = freq_top, label = 'Top Stetch:', valmin = 50, valmax = 100, valinit = init_top, color = '#E48671')
    freq_bottom = plt.axes([0.25, 0.087, 0.65, 0.03])
    slider_bottom = Slider(ax = freq_bottom, label = 'Bottom Stetch:', valmin = 0, valmax = 50, valinit = init_bot, color = '#E48671')

    #Adds a slider for the circle size
    circle_slid_location = plt.axes([0.25, 0.055, 0.65, 0.03])
    circle_slider = Slider(ax = circle_slid_location, label = 'Circle Size:', valmin = 5, valmax = (circle_size + 1*radius), valinit = circle_size, color = '#E48671')

    #Adds a notes section that the user can add notes about their data
    axbox = plt.axes([0.15, 0.02, 0.8, 0.03])
    text = ''
    text_box = TextBox(axbox, 'Notes:', initial = text, textalignment="center")

    #Make a button that can be clicked if no object is found
    axes_button = plt.axes([0.04, 0.775, 0.92, 0.04])
    close = Button(axes_button, 'Object Not Found', color = '#E48671')

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
      for d in range(len(default)):
        if default == [False, False, False]:
          total_data = real_data[0]*0
        if default[d] == True: 
          total_data = total_data + real_data[d]
        else: 
          pass
      norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
      ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

    #Updates the scaling when the slider is changed
    def update_slider_stretch(val):
      '''Updates the stretch the user can see, based in percentiles'''

      total_data = 0
      for d in range(len(default)):
        if default[d] == True: 
          total_data = total_data + real_data[d]
        else: 
          pass
      norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
      ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

    #Updates the notes added by the user when there is an input
    text_list = [text]
    def submit(expression):
      '''Updates the list of types in the 'Notes' setting'''

      text = expression
      text_list.append(text)

    #Allows the sliders and buttons to be pressed
    check.on_clicked(update_button)
    slider_top.on_changed(update_slider_stretch)
    slider_bottom.on_changed(update_slider_stretch)
    text_box.on_text_change(submit)

    #Display image until it is clicked to find the object
    n = -1
    while True:
      press = plt.waitforbuttonpress()
      text_max = len(text_list) - 1

      #Checks that it was a mouse click
      if press == False:
        n += 3

        #Finds which axes was clicked
        click_axes = str(location[n])
        click_axes = click_axes.split('WCSAxesSubplot', 2)[0]

        #Checks if the image was clicked
        if click_axes == '':

          #Makes a pop-up window with success text
          plt.clf()
          plt.close('all')
          plt.figure(1)
          plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for VISTA! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
          plt.xlim(0, 1)
          plt.ylim(0, 1)
          plt.grid(linewidth = 0)
          ax = plt.gca()
          ax.xaxis.set_tick_params(labelbottom=False)
          ax.yaxis.set_tick_params(labelleft=False)
          ax.set_xticks([])
          ax.set_yticks([])
          figure2 = plt.gcf()
          figure2.set_size_inches(4.75, 1)
          # figure2.canvas.set_window_title('Successful 2MASS Search')
          mng2 = pyplot.get_current_fig_manager()
          mng2.window.resizable(False, False)
          plt.pause(0.1)
          plt.clf()
          plt.close('all')

          #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
          coord = wcs_cropped.pixel_to_world_values(location[n-4],location[n-5])
          distance = []
          for i in range(len(object_ra)):
            distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))

          list_location = distance.index(np.min(distance))
          ymjd, jmjd = object_mjd_y[list_location], object_mjd_j[list_location]
          hmjd, ksmjd = object_mjd_h[list_location], object_mjd_ks[list_location]
          ra_vsa, dec_vsa = object_ra[list_location], object_dec[list_location]
          y, y_e = Y_list[list_location], Y_list_e[list_location]
          j, j_e = J_list[list_location], J_list_e[list_location]
          h, h_e = H_list[list_location], H_list_e[list_location]
          ks, ks_e = Ks_list[list_location], Ks_list_e[list_location]
          return ra_vsa, dec_vsa, y, y_e, j, j_e, h, h_e, ks, ks_e, ymjd, jmjd, hmjd, ksmjd, img, text_list[text_max]
        
        #Checks if the Object not Found button was clicked
        elif click_axes == 'Axes(0.04,0.775;0.92x0.04)':
          ymjd, jmjd, hmjd, ksmjd, y, y_e, j, j_e, h, h_e, ks, ks_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
          ra_vsa, dec_vsa = ra, dec
          shape_x, shape_y = total_data.shape[0], total_data.shape[1]

          #Makes a pop-up window with success text
          plt.clf()
          plt.close('all')
          plt.figure(1)
          plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for VISTA! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
          plt.xlim(0, 1)
          plt.ylim(0, 1)
          plt.grid(linewidth = 0)
          ax = plt.gca()
          ax.xaxis.set_tick_params(labelbottom=False)
          ax.yaxis.set_tick_params(labelleft=False)
          ax.set_xticks([])
          ax.set_yticks([])
          figure2 = plt.gcf()
          figure2.set_size_inches(4.75, 1)
          # figure2.canvas.set_window_title('Successful 2MASS Search')
          mng2 = pyplot.get_current_fig_manager()
          mng2.window.resizable(False, False)
          plt.pause(0.1)
          plt.clf()
          plt.close('all')
          
          return ra_vsa, dec_vsa, y, y_e, j, j_e, h, h_e, ks, ks_e, ymjd, jmjd, hmjd, ksmjd, img, 'Object Not Found was Pressed'
        
        #Allows updates for the circle size slider bar
        elif click_axes == 'Axes(0.25,0.055;0.65x0.03)':
          scatter.remove()
          scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')

      #Checks if the window was closed
      elif press is None:
        ymjd, jmjd, hmjd, ksmjd, y, y_e, j, j_e, h, h_e, ks, ks_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        ra_vsa, dec_vsa = ra, dec
        plt.clf()
        plt.close('all')
        return ra_vsa, dec_vsa, y, y_e, j, j_e, h, h_e, ks, ks_e, ymjd, jmjd, hmjd, ksmjd, img, text_list[text_max]
      
  #Returns null values if the images are not found
  else: 
    ymjd, jmjd, hmjd, ksmjd, y, y_e, j, j_e, h, h_e, ks, ks_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    ra_vsa, dec_vsa = ra, dec
    text_list = 'Image Not Found'
    return ra_vsa, dec_vsa, y, y_e, j, j_e, h, h_e, ks, ks_e, ymjd, jmjd, hmjd, ksmjd, 'VSA', text_list
  
def vsa_table(ra, dec, radius): 
  '''Find all the objects in the radius defined by the user'''

  blockPrint()
    
  #Find the location of all the object found in VSA in the radius choosen by the user
  catalog_list = ['VHS', 'VVV', 'VMC', 'VIKING', 'VIDEO', 'UltraVISTA']
  database_list = ['VHSDR6', 'VVVDR5', 'VMCDR6', 'VIKINGDR5', 'VIDEODR6', 'ULTRAVISTADR4']
  for i in range(len(catalog_list)): 
    table = Vsa.query_region(
        SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'),
                radius = ((radius/2) - 1) * u.arcsec,
                programme_id = catalog_list[i],
                database = database_list[i])
    if len(table) > 0:
      return table
  
#                       WISEVIEW LINK                           #
# ------------------------------------------------------------- #
def wiseview_link(ra, dec, radius): 
  '''Opens WISEView link to aid the user in finding their object'''

  #Gets the link for SMDET
  link = 'http://byw.tools/wiseview#ra=' + str(ra) + '&dec=' + str(dec) + '&size=' + str(radius) + '&band=3&speed=164.06&minbright=-10.0000&maxbright=80&window=0.75&diff_window=1&linear=1&color=&zoom=8.5&border=1&gaia=1&invert=1&maxdyr=0&scandir=0&neowise=0&diff=0&outer_epochs=0&unique_window=1&smooth_scan=0&shift=0&pmra=0&pmdec=0&synth_a=0&synth_a_sub=0&synth_a_ra=&synth_a_dec=&synth_a_w1=&synth_a_w2=&synth_a_pmra=0&synth_a_pmdec=0&synth_a_mjd=&synth_b=0&synth_b_sub=0&synth_b_ra=&synth_b_dec=&synth_b_w1=&synth_b_w2=&synth_b_pmra=0&synth_b_pmdec=0&synth_b_mjd=&smdet_coadd_id=1863p620&smdet_mask_idx=3&smdet_obj_px=&smdet_obj_py='
  webbrowser.open(link)

#                       GUI FUNCTIONS                           #
# ------------------------------------------------------------- #
def single_object_search(): 
  '''Performs the single object search by going to each catalog script and recording the data down.
  Once all the data is recorded it it put into a CSV file for the user.'''

  #Runs the wiseview_link function
  wiseview_link(ra_use, dec_use, radius_use)

  #Creates fake list for the data and data names
  photometry = []
  photometry_name = []
  photometry.append([ra_use, dec_use, radius_use])
  photometry_name.append(['input_ra', 'input_dec', 'input_radius'])

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
      new_directory = (directory.split('Output/')[0]) + '/Output/'
      myFile = open(str(new_directory) + str(output) + '.csv', 'w')
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
  if values['RADIUS'].isnumeric() == False or (int(values['RADIUS']) > 500) == True or (int(values['RADIUS']) < 100) == True:
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
    photometry.append([ra_use, dec_use, radius_use])
    photometry_name.append(['input_ra', 'input_dec', 'input_radius'])

    #Runs the wiseview_link function
    wiseview_link(ra_use, dec_use, radius_use)

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
  if values['RADIUS_multi'].isnumeric() == False or (int(values['RADIUS_multi']) > 500) == True or (int(values['RADIUS_multi']) < 100) == True:
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
sg.theme('LightBrown3')

#                        GUI LAYOUT                             #
# ------------------------------------------------------------- #
#Sets the layout if the user is not on a windows machine
if platform != 'win32':
  #Makes the layout of WRAP for the single object search, by providing a location for: ra, dec, radius, output file name, catalogs, and output
  layout_single = [[sg.Image(filename = (str(directory) + '/Output/metadata/WRAP_Logo.png'), size = (135, 95)),                    sg.Text('WRAP', justification='center', size=(6, 1), font = ('Chalkduster', 45)),                 sg.Image(filename = str(directory) + '/Output/metadata/BYW_Logo.png', size = (205, 95))],
                   
                    [sg.Text('RA', font = ('Times New Roman', 22), size=(13, 1), justification='center'),           sg.Text('DEC', font = ('Times New Roman', 22), size=(13, 1), justification='center'),             sg.Text('RADIUS', font = ('Times New Roman', 22), size=(13, 1), justification='center')],
                    [sg.Text('(Degrees)', font = ('Times New Roman', 20), size=(18, 1), justification='center'),    sg.Text('(Degrees)', font = ('Times New Roman', 20), size=(11, 1), justification='center'),       sg.Text('(Arcsecs)', font = ('Times New Roman', 20), size=(20, 1), justification='center')],
                    [sg.InputText(size=(18, 1), key = 'RA', font = ('Times New Roman', 15)),                        sg.InputText(size=(18, 2), key = 'DEC', font = ('Times New Roman', 15)),                          sg.InputText(size=(18, 2), key = 'RADIUS', font = ('Times New Roman', 15))],
                    
                    [sg.Text('Output File Name', size=(50, 1), justification='center', font = ('Times New Roman', 22))],
                    [sg.InputText(key = 'output', font = ('Times New Roman', 15), size = (70, 3), justification='center')],

                    [sg.Text('Catalogs:', justification='center', size=(50, 1), font = ('Times New Roman', 25))],   
                    [sg.Checkbox('CatWISE 2020', key = 'catwise', font = ('Times New Roman', 22), size = (14, 2)),   sg.Checkbox('AllWISE', key = 'AW', font = ('Times New Roman', 22), size = (10, 2)),               sg.Checkbox('Gaia', key = 'gaia', font = ('Times New Roman', 22), size = (9, 2))],
                    [sg.Checkbox('VISTA', key = 'VSA', font = ('Times New Roman', 22), size = (9, 2)),               sg.Checkbox('WFCAM', key = 'UKIDSS', font = ('Times New Roman', 22), size = (10, 2)),             sg.Checkbox('2MASS', key = '2MASS', font = ('Times New Roman', 22), size = (10, 2))],
                    [sg.Checkbox('PanSTARRS', key = 'ps', font = ('Times New Roman', 22), size = (13, 2)),           sg.Checkbox('NSC', key = 'nsc', font = ('Times New Roman', 22), size = (8, 2)),                   sg.Checkbox('GALEX', key = 'galex', font = ('Times New Roman', 22), size = (10, 2))],
                    [sg.Checkbox('Select All',   enable_events=True, key='Check_All'),                               sg.Checkbox('Deselect All', enable_events=True, key='Uncheck_All')],

                    [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'),                                   sg.Button('Help', size = (17), button_color = '#F7CC7C'),                                         sg.Button('Close WRAP', size = (17), button_color = '#E48671')], 

                    [sg.Text("\u0332".join('Output'), size=(50, 1), justification='center', font = ('Times New Roman', 15))],
                    [sg.Multiline(size=(85, 6), write_only=(True), key=ML_KEY_SINGLE, reroute_stdout=True, echo_stdout_stderr=True, reroute_cprint=True)]]

  #Makes the drop down window for types of file in the multi-object search
  filetype_list = ['CSV', 'FITS', 'ASCII', 'IPAC']
  #Makes the layout of WRAP for the multi-object search, by providing a location for: file directory, radius, filetype, output file name, catalogs, and output
  layout_multi = [[sg.Image(filename = (str(directory) + '/Output/metadata/WRAP_Logo.png'), size = (135, 95)),                    sg.Text('WRAP', justification='center', size=(6, 1), font = ('Chalkduster', 45)),                 sg.Image(filename = str(directory) + '/Output/metadata/BYW_Logo.png', size = (205, 95))],
                  
                  [sg.Text('FILE DIRECTORY', font = ('Times New Roman', 22), size=(50, 1), justification='center')],
                  [sg.Text('(CSV, FITS, ASCII, IPAC)', font = ('Times New Roman', 20), size=(50, 1), justification='center')],

                  [sg.FileBrowse('File Browser', size = (80, 1), key = 'file', file_types = [('CSV Files', '*.csv'), ('FITS Files', '*.fits'), ('ASCII Files', '*.txt'), ('IPAC Files', '*.txt')])],
                  [sg.Text('RADIUS', font = ('Times New Roman', 22), size=(17, 1), justification='center'),              sg.Text('FILETYPE', font = ('Times New Roman', 22), size=(9, 1), justification='center'),               sg.Text('Output File Name', size=(25, 1), justification='center', font = ('Times New Roman', 22))],
                  [sg.InputText(size=(22, 2), key = 'RADIUS_multi', font = ('Times New Roman', 15)),                     sg.Combo(filetype_list, size = (13), font = ('Times New Roman', 15), key = 'type'),                     sg.InputText(key = 'output2', font = ('Times New Roman', 15), size = (22, 2), justification='center')],

                  [sg.Text('Catalogs:', justification='center', size=(50, 1), font = ('Times New Roman', 25))],
                  [sg.Checkbox('CatWISE 2020', key = 'catwise_multi', font = ('Times New Roman', 22), size = (14, 2)),   sg.Checkbox('AllWISE', key = 'AW_multi', font = ('Times New Roman', 22), size = (10, 2)),               sg.Checkbox('Gaia', key = 'gaia_multi', font = ('Times New Roman', 22), size = (9, 2))],
                  [sg.Checkbox('VISTA', key = 'VSA_multi', font = ('Times New Roman', 22), size = (9, 2)),               sg.Checkbox('WFCAM', key = 'UKIDSS_multi', font = ('Times New Roman', 22), size = (10, 2)),             sg.Checkbox('2MASS', key = '2MASS_multi', font = ('Times New Roman', 22), size = (10, 2))],
                  [sg.Checkbox('PanSTARRS', key = 'ps_multi', font = ('Times New Roman', 22), size = (13, 2)),           sg.Checkbox('NSC', key = 'nsc_multi', font = ('Times New Roman', 22), size = (8, 2)),                   sg.Checkbox('GALEX', key = 'galex_multi', font = ('Times New Roman', 22), size = (10, 2))],
                  [sg.Checkbox('Select All',   enable_events=True, key='Check_All_Multi'),                               sg.Checkbox('Deselect All', enable_events=True, key='Uncheck_All_Multi')],

                  [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'),                                         sg.Button('Help', size = (17), button_color = '#F7CC7C'),                                               sg.Button('Close WRAP', size = (17), button_color = '#E48671')], 

                  [sg.Text("\u0332".join('Output'), size=(50, 1), justification='center', font = ('Times New Roman', 15))],
                  [sg.Multiline(size=(85, 6), write_only=(True), key=ML_KEY_MULTI, reroute_stdout=True, echo_stdout_stderr=True, reroute_cprint=True)]]

  #Makes the general layout for WRAP
  tab_layout = [[sg.TabGroup([[sg.Tab('Single Object',   layout_single,       title_color='#F9F8F3',          background_color='#eeeccb',   element_justification= 'center',     key = 'Single Obect Search'),
                              sg.Tab('Multi-Object',    layout_multi,        title_color='#F9F8F3',          background_color='#eeeccb',   element_justification= 'center',     key = 'Multi-Object Search')]], 
                              tab_location='centertop', title_color='Black', tab_background_color='#F9F8F3', selected_title_color='Black', selected_background_color='#9C873E', border_width = 6, font = ('Times New Roman', 18), enable_events = True, key = 'tab_group'), sg.Button('Close')]] 

  #Generates the window based off the layouts above
  window = sg.Window('WRAP', tab_layout, size = (550, 600), grab_anywhere=False, finalize=True, enable_close_attempted_event = True)
  
  #Makes lists for: catalog names, catalog functions, and fake variables
  catalogs = ['catwise', 'AW',     'gaia', 
              'VSA',     'UKIDSS', '2MASS', 
              'ps',      'nsc',    'galex']
  catalogs_multi = ['catwise_multi',  'AW_multi',     'gaia_multi', 
                    'VSA_multi',      'UKIDSS_multi', '2MASS_multi', 
                    'ps_multi',       'nsc_multi',    'galex_multi']
  catalog_names = ['CatWISE 2020', 'AllWISE', 'GAIA', 
                  'VSA',           'UKIDSS',  '2MASS', 
                  'PanSTARRS',     'NSC',     'GALEX']
  catalog_functions = [catwise_image,   allwise_image,  gaia_image, 
                      vsa_image,        ukidss_image,   twomass_image, 
                      ps_image,         nsc_image,      galex_image]

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
  header = ['input_ra', 'input_dec', 'input_radius', 
            'cw_ra', 'cw_ra_e', 'cw_dec', 'cw_dec_e', 'cw_w1', 'cw_w1_e', 'cw_w2', 'cw_w2_e', 'cw_pmra', 'cw_pmra_e', 'cw_pmdec', 'cw_pmdec_e', 'cw_mjd', 'cw_catalog', 'cw_notes', 
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
  layout_single = [[sg.Image(filename = (str(directory) + '/Output/metadata/WRAP_Logo.png'), size = (135, 95)),                    sg.Text('WRAP', justification='center', size=(6, 1), font = ('Ink Free', 45)),                 sg.Image(filename = str(directory) + '/Output/metadata/BYW_Logo.png', size = (205, 95))],
                    [sg.Text('RA', font = ('Times New Roman', 18), size=(13, 1), justification='center'),           sg.Text('DEC', font = ('Times New Roman', 18), size=(10, 1), justification='center'),             sg.Text('RADIUS', font = ('Times New Roman', 18), size=(13, 1), justification='center')],
                    [sg.Text('(Degrees)', font = ('Times New Roman', 16), size=(14, 1), justification='center'),    sg.Text('(Degrees)', font = ('Times New Roman', 16), size=(11, 1), justification='center'),       sg.Text('(Arcsecs)', font = ('Times New Roman', 16), size=(16, 1), justification='center')],
                    [sg.InputText(size=(15), key = 'RA', font = ('Times New Roman', 15)),                        sg.InputText(size=(15, 2), key = 'DEC', font = ('Times New Roman', 15)),                          sg.InputText(size=(15, 2), key = 'RADIUS', font = ('Times New Roman', 15))],
                    
                    [sg.Text('Output File Name', size=(50, 1), justification='center', font = ('Times New Roman', 18))],
                    [sg.InputText(key = 'output', font = ('Times New Roman', 15), size = (70, 3), justification='center')],

                    [sg.Text('Catalogs:', justification='center', size=(50, 1), font = ('Times New Roman', 22))],   
                    [sg.Checkbox('CatWISE 2020', key = 'catwise', font = ('Times New Roman', 15), size = (14, 2)),   sg.Checkbox('AllWISE', key = 'AW', font = ('Times New Roman', 15), size = (10, 2)),               sg.Checkbox('Gaia', key = 'gaia', font = ('Times New Roman', 15), size = (9, 2))],
                    [sg.Checkbox('VISTA', key = 'VSA', font = ('Times New Roman', 15), size = (9, 2)),               sg.Checkbox('WFCAM', key = 'UKIDSS', font = ('Times New Roman', 15), size = (10, 2)),             sg.Checkbox('2MASS', key = '2MASS', font = ('Times New Roman', 15), size = (10, 2))],
                    [sg.Checkbox('PanSTARRS', key = 'ps', font = ('Times New Roman', 15), size = (13, 2)),           sg.Checkbox('GALEX', key = 'galex', font = ('Times New Roman', 15), size = (10, 2))],
                    [sg.Checkbox('Select All',   enable_events=True, key='Check_All'),                               sg.Checkbox('Deselect All', enable_events=True, key='Uncheck_All')],

                    [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'),                                   sg.Button('Help', size = (17), button_color = '#F7CC7C'),                                         sg.Button('Close WRAP', size = (17), button_color = '#E48671')], 

                    [sg.Text("\u0332".join('Output'), size=(50, 1), justification='center', font = ('Times New Roman', 15))],
                    [sg.Multiline(size=(85, 6), write_only=(True), key=ML_KEY_SINGLE, reroute_stdout=True, echo_stdout_stderr=True, reroute_cprint=True)]]

  #Makes the drop down window for types of file in the multi-object search
  filetype_list = ['CSV', 'FITS', 'ASCII', 'IPAC']
  #Makes the layout of WRAP for the multi-object search, by providing a location for: file directory, radius, filetype, output file name, catalogs, and output
  layout_multi = [[sg.Image(filename = (str(directory) + '/Output/metadata/WRAP_Logo.png'), size = (135, 95)),                    sg.Text('WRAP', justification='center', size=(6, 1), font = ('Ink Free', 45)),                 sg.Image(filename = str(directory) + '/Output/metadata/BYW_Logo.png', size = (205, 95))],
                  
                  [sg.Text('FILE DIRECTORY', font = ('Times New Roman', 18), size=(50, 1), justification='center')],
                  [sg.Text('(CSV, FITS, ASCII, IPAC)', font = ('Times New Roman', 16), size=(50, 1), justification='center')],
                  [sg.FileBrowse('File Browser', size = (80, 1), key = 'file', file_types = [('CSV Files', '*.csv'), ('FITS Files', '*.fits'), ('ASCII Files', '*.txt'), ('IPAC Files', '*.txt')])],

                  [sg.Text('RADIUS', font = ('Times New Roman', 16), size=(14, 1), justification='center'),              sg.Text('FILETYPE', font = ('Times New Roman', 16), size=(10, 1), justification='center'),               sg.Text('Output File Name', size=(26, 1), justification='center', font = ('Times New Roman', 16))],
                  [sg.InputText(size=(16, 2), key = 'RADIUS_multi', font = ('Times New Roman', 15)),                     sg.Combo(filetype_list, size = (13), font = ('Times New Roman', 15), key = 'type'),                     sg.InputText(key = 'output2', font = ('Times New Roman', 15), size = (22, 2), justification='center')],

                  [sg.Text('Catalogs:', justification='center', size=(50, 1), font = ('Times New Roman', 25))],
                  [sg.Checkbox('CatWISE 2020', key = 'catwise_multi', font = ('Times New Roman', 15), size = (14, 2)),   sg.Checkbox('AllWISE', key = 'AW_multi', font = ('Times New Roman', 15), size = (10, 2)),               sg.Checkbox('Gaia', key = 'gaia_multi', font = ('Times New Roman', 15), size = (9, 2))],
                  [sg.Checkbox('VISTA', key = 'VSA_multi', font = ('Times New Roman', 15), size = (9, 2)),               sg.Checkbox('WFCAM', key = 'UKIDSS_multi', font = ('Times New Roman', 15), size = (10, 2)),             sg.Checkbox('2MASS', key = '2MASS_multi', font = ('Times New Roman', 15), size = (10, 2))],
                  [sg.Checkbox('PanSTARRS', key = 'ps_multi', font = ('Times New Roman', 15), size = (13, 2)),           sg.Checkbox('GALEX', key = 'galex_multi', font = ('Times New Roman', 15), size = (10, 2))],
                  [sg.Checkbox('Select All',   enable_events=True, key='Check_All_Multi'),                               sg.Checkbox('Deselect All', enable_events=True, key='Uncheck_All_Multi')],

                  [sg.Button('Run WRAP', size = (17), button_color = '#95D49B'),                                         sg.Button('Help', size = (17), button_color = '#F7CC7C'),                                               sg.Button('Close WRAP', size = (17), button_color = '#E48671')], 

                  [sg.Text("\u0332".join('Output'), size=(50, 1), justification='center', font = ('Times New Roman', 15))],
                  [sg.Multiline(size=(85, 6), write_only=(True), key=ML_KEY_MULTI, reroute_stdout=True, echo_stdout_stderr=True, reroute_cprint=True)]]

  #Makes the general layout for WRAP
  tab_layout = [[sg.TabGroup([[sg.Tab('Single Object',   layout_single,       title_color='#F9F8F3',          background_color='#eeeccb',   element_justification= 'center',     key = 'Single Obect Search'),
                              sg.Tab('Multi-Object',    layout_multi,        title_color='#F9F8F3',          background_color='#eeeccb',   element_justification= 'center',     key = 'Multi-Object Search')]], 
                              tab_location='centertop', title_color='Black', tab_background_color='#F9F8F3', selected_title_color='Black', selected_background_color='#9C873E', border_width = 6, font = ('Times New Roman', 18), enable_events = True, key = 'tab_group'), sg.Button('Close')]] 


  #Generates the window based off the layouts above
  window = sg.Window('WRAP', tab_layout, size = (550, 680), grab_anywhere=False, finalize=True, enable_close_attempted_event = True)

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
  header = ['input_ra', 'input_dec', 'input_radius', 
            'cw_ra', 'cw_ra_e', 'cw_dec', 'cw_dec_e', 'cw_w1', 'cw_w1_e', 'cw_w2', 'cw_w2_e', 'cw_pmra', 'cw_pmra_e', 'cw_pmdec', 'cw_pmdec_e', 'cw_mjd', 'cw_catalog', 'cw_notes', 
            'aw_ra', 'aw_ra_e', 'aw_dec', 'aw_dec_e', 'aw_w1', 'aw_w1_e', 'aw_w2', 'aw_w2_e', 'aw_w3', 'aw_w3_e', 'aw_w4', 'aw_w4_e', 'aw_pmra', 'aw_pmra_e', 'aw_pmdec', 'aw_pmdec_e', 'aw_catalog', 'aw_notes', 
            'gaia_ra', 'gaia_ra_e', 'gaia_dec', 'gaia_dec_e', 'gaia_parallax', 'gaia_parallax_e', 'gaia_radv', 'gaia_radv_e', 'gaia_pmra', 'gaia_pmra_e', 'gaia_pmdec', 'gaia_pmdec_e', 'gaia_g', 'gaia_g_e', 'gaia_bp', 'gaia_bp_e', 'gaia_rp', 'gaia_rp_e', 'gaia_year', 'gaia_catalog', 'gaia_notes', 
            'vsa_ra', 'vsa_dec', 'vsa_y', 'vsa_y_e', 'vsa_j', 'vsa_j_e', 'vsa_h', 'vsa_h_e', 'vsa_ks', 'vsa_ks_e', 'vsa_mjd_y', 'vsa_mjd_j', 'vsa_mjd_h', 'vsa_mjd_ks', 'vsa_catalog', 'vsa_notes', 
            'wfcam_ra', 'wfcam_ra_e', 'wfcam_dec', 'wfcam_dec_e', 'wfcam_y', 'wfcam_y_e', 'wfcam_j', 'wfcam_j_e', 'wfcam_h', 'wfcam_h_e', 'wfcam_k', 'wfcam_k_e', 'wfcam_pmra', 'wfcam_pmra_e', 'wfcam_pmdec', 'wfcam_pmdec_e', 'wfcam_epoch', 'wfcam_catalog', 'wfcam_notes', 
            '2mass_ra', '2mass_dec', '2mass_j', '2mass_j_e', '2mass_h', '2mass_h_e', '2mass_ks', '2mass_ks_e', '2mass_catalog', '2mass_notes', 
            'ps_ra', 'ps_ra_e', 'ps_dec', 'ps_dec_e', 'ps_g', 'ps_g_e', 'ps_r', 'ps_r_e', 'ps_i', 'ps_i_e', 'ps_z', 'ps_z_e', 'ps_y', 'ps_y_e', 'ps_mjd', 'ps_catalog', 'ps_notes', 
            'nsc_ra', 'nsc_ra_e', 'nsc_dec', 'nsc_dec_e', 'nsc_g', 'nsc_g_e', 'nsc_r', 'nsc_r_e', 'nsc_i', 'nsc_i_e', 'nsc_z', 'nsc_z_e', 'nsc_u', 'nsc_u_e', 'nsc_y', 'nsc_y_e', 'nsc_pmra', 'nsc_pmra_e', 'nsc_pmdec', 'nsc_pmdec_e', 'nsc_mjd', 'nsc_catalog', 'nsc_notes',  
            'galex_ra', 'galex_dec', 'galex_fuv', 'galex_fuv_e', 'galex_nuv', 'galex_nuv_e', 'galex_catalog', 'galex_notes']
  
#                             GUI                               #
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
        window['catwise'].update(True), window['AW'].update(True), window['gaia'].update(True)
        window['VSA'].update(True), window['UKIDSS'].update(True), window['2MASS'].update(True)
        window['ps'].update(True), window['galex'].update(True)
        window['Uncheck_All'].update(False)
        if platform != 'win32':
          window['nsc'].update(True)
    elif event == 'Uncheck_All':
        window['catwise'].update(False), window['AW'].update(False), window['gaia'].update(False)
        window['VSA'].update(False), window['UKIDSS'].update(False), window['2MASS'].update(False)
        window['ps'].update(False), window['galex'].update(False)
        window['Check_All'].update(False) 
        if platform != 'win32':
           window['nsc'].update(False)

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
      sg.cprint('------------------------------------------------                                                                                                      ', c='wheat4',         key = ML_KEY_SINGLE)
      sg.cprint('     Thank you for using WRAP!                                                                                                                        ', c='wheat4', end='', key = ML_KEY_SINGLE)
      sg.cprint('   Authors Contact: hcb98@nau.edu                                                                                                                     ', c='wheat4', end='', key = ML_KEY_SINGLE)
      sg.cprint('------------------------------------------------                                                                                                      ', c='wheat4',         key = ML_KEY_SINGLE)

    #Closes WRAP if the 'Close WRAP' button is pressed
    if event in (None, 'Close WRAP'):
      break

  #Runs the program if the user is in 'Multi-Object' tab
  if group == 'Multi-Object Search':

    #Adds the functionality of the 'Select All' and 'Deselect All' Buttons
    if event == 'Check_All_Multi':
        window['catwise_multi'].update(True), window['AW_multi'].update(True), window['gaia_multi'].update(True)
        window['VSA_multi'].update(True), window['UKIDSS_multi'].update(True), window['2MASS_multi'].update(True)
        window['ps_multi'].update(True), window['galex_multi'].update(True)
        window['Uncheck_All_Multi'].update(False)
        if platform != 'win32':
          window['nsc_multi'].update(True)
    elif event == 'Uncheck_All_Multi':
        window['catwise_multi'].update(False), window['AW_multi'].update(False), window['gaia_multi'].update(False)
        window['VSA_multi'].update(False), window['UKIDSS_multi'].update(False), window['2MASS_multi'].update(False)
        window['ps_multi'].update(False), window['galex_multi'].update(False)
        window['Check_All_Multi'].update(False) 
        if platform != 'win32':
          window['nsc_multi'].update(False)

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
      sg.cprint('------------------------------------------------                                                                                              ', end='', c='wheat4', key = ML_KEY_MULTI)
      sg.cprint('     Thank you for using WRAP!                                                                                                                ', end='', c='wheat4', key = ML_KEY_MULTI)
      sg.cprint('   Authors Contact: hcb98@nau.edu                                                                                                             ', end='', c='wheat4', key = ML_KEY_MULTI)
      sg.cprint('------------------------------------------------                                                                                              ', end='', c='wheat4', key = ML_KEY_MULTI)

    #Closes WRAP if the 'Close WRAP' button is pressed
    if event in (None, 'Close WRAP2'):
      break

#Closes the window
window.close()