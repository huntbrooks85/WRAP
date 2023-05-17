#Import all of the packages
from catalog_module.importmodule import *

def blockPrint():
  '''Makes a function that blocks the printing function'''
  sys.stdout = open(os.devnull, 'w')

def wiseview_link(ra, dec, radius): 
  '''Opens WISEView link to aid the user in finding their object'''

  #Gets the link for SMDET
  link = 'http://byw.tools/wiseview#ra=' + str(ra) + '&dec=' + str(dec) + '&size=' + str(radius) + '&band=3&speed=164.06&minbright=-10.0000&maxbright=80&window=0.75&diff_window=1&linear=1&color=&zoom=8.5&border=1&gaia=1&invert=1&maxdyr=0&scandir=0&neowise=0&diff=0&outer_epochs=0&unique_window=1&smooth_scan=0&shift=0&pmra=0&pmdec=0&synth_a=0&synth_a_sub=0&synth_a_ra=&synth_a_dec=&synth_a_w1=&synth_a_w2=&synth_a_pmra=0&synth_a_pmdec=0&synth_a_mjd=&synth_b=0&synth_b_sub=0&synth_b_ra=&synth_b_dec=&synth_b_w1=&synth_b_w2=&synth_b_pmra=0&synth_b_pmdec=0&synth_b_mjd=&smdet_coadd_id=1863p620&smdet_mask_idx=3&smdet_obj_px=&smdet_obj_py='
  webbrowser.open(link)
