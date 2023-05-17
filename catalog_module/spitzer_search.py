#Import all of the packages
from catalog_module.importmodule import *

def blockPrint():
  '''Makes a function that blocks the printing function'''
  sys.stdout = open(os.devnull, 'w')

def enablePrint():
  '''Makes a function that allows the printing function'''
  sys.stdout = sys.__stdout__

def spitzer_image(ra, dec, radius): 
  ''' First, it gets the images from the Spitzer API from IRSA and downloads the images. 
  Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
  Third, calls the table function to get all the objects from the Spitzer source catalog. 
  Fourth, makes the window for the user to click the object with all settings. 
  Finally, finds the closest object to the click and records the data. '''
  pos_t2 = sha.query(ra = ra, dec = dec, size = (radius/3600))
  return pos_t2


def allwise_table(ra, dec, radius): 
  '''Find all the objects in the radius defined by the user'''

  blockPrint()
 
  #Uses astroquery to find all objects in the radius
  location_data = Irsa.query_region(coord.SkyCoord(ra, dec, unit = (u.deg,u.deg), frame = 'fk5'), catalog = 'allwise_p3as_psd', spatial = 'Box', width = (radius - 1) * u.arcsec)
  return location_data