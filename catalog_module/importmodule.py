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
import matplotlib
import numpy as np
import pandas as pd
from pyvo.dal import sia
import PySimpleGUI as sg
from astropy.wcs import WCS
from astropy.io import fits
from astropy.io import ascii
from bs4 import BeautifulSoup
from astroquery.vsa import Vsa
from astropy import units as u
from astropy.table import Table
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia
from astropy.nddata import Cutout2D
import astropy.coordinates as coord
from dl.helpers.utils import convert
from astroquery.ukidss import Ukidss
from astroquery.ipac.irsa import Irsa
from astropy.coordinates import SkyCoord
from astropy.utils.data import download_file
from dl import authClient as ac, queryClient as qc
from astroquery.mast import Observations, MastMissions, Catalogs
from matplotlib.widgets import Slider, Button, CheckButtons, TextBox
from astropy.visualization import (PercentileInterval, SinhStretch, ImageNormalize)

#Checks if the users system is greater than Python 3.8
# ------------------------------------------------------------- #
ver = sys.version_info
if ver.major >= 3.8:
    sys.exit('Python 3.8, or newer, is required to run WRAP!')

#Ignores all of the warnings from the packages above
# ------------------------------------------------------------- #
import warnings
warnings.filterwarnings("ignore")
