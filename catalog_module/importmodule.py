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
from bs4 import BeautifulSoup
from astroquery.vsa import Vsa
from astropy import units as u
from astropy.table import Table
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia
from matplotlib import transforms
from astropy.nddata import Cutout2D
import astropy.coordinates as coord
# from astroquery.ipac.irsa import sha
from astroquery.ukidss import Ukidss
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

#Does not load these packages if the user is on a Windows machine
# ------------------------------------------------------------- #
if platform != 'win32':
    from dl.helpers.utils import convert
    from dl import authClient as ac, queryClient as qc

#Checks if the users system is greater than Python 3.8 and is not a Linux machine
# ------------------------------------------------------------- #
ver = sys.version_info
if ver.major < 3:
    sys.exit('Python 3.X is required to run WRAP!')

if platform == "linux" or platform == "linux2":
    sys.exit('WRAP is not supported on Linux, please use another operating system!')

#Ignores all of the warnings from the packages above
# ------------------------------------------------------------- #
import warnings
warnings.filterwarnings("ignore")