Purpose: 
Gathers photometry and astrometry from various ultra-violet, optical, and near-infrared catalogs given a RA, DEC, and Radius by the user.

Catalog Full Names and References: 
 CatWISE2020 & AllWISE -> Wide-Field Infrared Survey Explorer, https://irsa.ipac.caltech.edu/Missions/wise.html
 Gaia -> GAIA ESA, https://gea.esac.esa.int/archive/
 VISTA -> Visible and Infrared Survey Telescope for Astronomy, http://horus.roe.ac.uk/vsa/index.html
 WFCAM -> Wide Field Camera, http://wsa.roe.ac.uk/index.html
 2MASS -> Two Micron All Sky Survey, https://irsa.ipac.caltech.edu/Missions/2mass.html
 PanSTARRS -> Panoramic Survey Telescope and Rapid Response System, https://outerspace.stsci.edu/display/PANSTARRS/
 NSC -> NOIRLab Source Catalog, https://datalab.noirlab.edu/nscdr2/index.php
 GALEX -> Galaxy Evolution Explorer, http://www.galex.caltech.edu

Python: 
-> Python 3.8.8 is needed for WRAP
-> Python 3.8.8 link: https://www.python.org/downloads/release/python-388/

Packages Needed: 
-> csv
-> math
-> sys
-> os
-> astropy
-> requests
-> truncate
-> matplotlib
-> numpy = 1.22.0
-> pandas
-> pyvo
-> PySimpleGUI
-> bs4
-> astroquery
-> astro-datalab

How to Start Program: 
-> Only tested on MacOS Ventura 13.2.1, problems may occur for older versions of MacOS and versions for Windows and Linux will be made in the future. 
-> Note that the window close button has been disabled, to close WRAP please click the red "Close WRAP" button. 

1) Open terminal 
2) Go to the file directory that contains "WRAP.py" (example: "cd /Users/hunter_brooks8/Desktop/WRAP_Project/WRAP")
3) Type "python WRAP.py" into the terminal

How to use: 
1) Select if you are find photometry and astrometry for a single object or multiple objects.

	For Single Object Search
2) Input the RA and DEC (in degrees) of the object being searched for
3) Input the radius (in arcsecs) around the RA and DEC that is being searched for
4) Put in the output file name that you want, do not put a file type at the end (it always saves to a csv file)
5) Select the catalogs that will be searched through
6) Click the green "Run Wrap" button
7) Click on the your object when the image appears (if your object is not there click the red "Object Not Found" button)
8) Repeat until every catalog image is shown
9) All data is saved to a csv file with your file output name in the "Output" folder (if no output name was put it will default to "WRAP_output.csv")

	For Multi-Object Search
2) Click the "Browse" button and click on you file with a list of RA ("ra") and DEC's ("dec") (in degrees)
3) Input the radius (in arcsecs) around the RA and DEC that is being searched for
4) Click the "Filetype" browser and select the file type being used
5) Put in the output file name that you want, do not put a file type at the end
6) Select the catalogs that will be searched through
7) Click the green "Run Wrap" button
8) Click on the your object when the image appears (if your object is not there click the red "Object Not Found" button)
9) Repeat until every catalog image is shown and every object has been run though
10) All data is saved to a csv file with your file output name in the "Output" folder (if no output name was put it will default to "WRAP_output.csv")

Author Contact: 
Hunter Brooks -> hcb98@nau.edu
