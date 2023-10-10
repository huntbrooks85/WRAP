<h1 align="center" id="title">✨ WRAP ✨</h1>

<p id="description"> The Wide-field Retrieval of Astrodata Program (WRAP) is a tool designed to aid the retrieval of photometric and astrometric data from a diverse range of ultraviolet, optical, and near-infrared catalogs. Users provide WRAP with Right Ascension (RA), Declination (DEC) coordinates and a search radius for their object of interest. The program then allows users to interactively select their object by clicking on it within catalog images. WRAP subsequently compiles and presents all relevant photometric and astrometric information gathered from the chosen catalogs in a CSV file, streamlining the data retrieval process.</p>

<h2>🛠️ Installation Steps:</h2>

<p>1. Python </p>

```
-> Python 3.8 or newer is needed for WRAP 
```
* Python 3.8 link: https://www.python.org/downloads/release/python-388/
* Python Installation Guide: https://wiki.python.org/moin/BeginnersGuide/Download

<p>2. Packages </p>

```
-> numpy
-> pandas
-> pyvo==1.4 
-> astropy==5.2.2
-> truncate==0.11 
-> requests==2.28.1 
-> astroquery==0.4.6 
-> matplotlib==3.4.2
-> PySimpleGUI==4.60.5
-> beautifulsoup4==4.11.1
-> OpenCV-python==4.7.0.72
```
* None Windows Users
```
-> astro-datalab== 2.20.1 
```

* PIP Installation Guide: https://pip.pypa.io/en/stable/installation/
* How to Use PIP: https://packaging.python.org/en/latest/tutorials/installing-packages/

<p>3. Fixing Astroquery Module </p>

1) Open the "replacement_module" directory
2) Copy the "core.py" file onto your Desktop
3) Go to the file directory containing the UKIDSS astroquery files (example: "/User/anaconda3/lib/python3.8/site-packages/astroquery/ukidss")
4) If you are struggling to find your file directory use these steps:
   1) Type: ```python3```, in your terminal
   2) Type: ```import astroquery```, in your terminal
   3) Type: ```astroquery.__file__```, in your terminal
   4) This will provide you with your astroquery file directory location
   5) Leave the ```python3``` window by pressing "control + z"
5) Open the "ukidss" folder in this directory, so that you can see the existing "core.py" file. 
6) Replace the "core.py" file in the astroquery UKIDSS directory with the "core.py" file copied on your Desktop. 

<h2> 🏆 Opening Application: </h2>

*   Note 1: Only tested on MacOS Ventura 13.X and Windows 11, problems may occur for older versions of MacOS and Windows (needs to be tested on Linux)
*   Note 2: Windows does not support the astro-datalab package, therefore Windows does not have the Noirlab Source Catalog option. 
*   Note 3: The window close button has been disabled, to close WRAP please click the red "Close WRAP" button at the bottom. 
*   Note 4: The orientation for all of the catalogs is North pointed up and East pointing left. 
*   Note 5: 2MASS can have strange imaging cropping, this is a warning that it may happen and do not be alarmed by it. 

<p> How to Start Program </p>

1) Open terminal 
2) Go to the file directory that contains "WRAP.py" (example: "cd /Users/hunter/Desktop/WRAP")
3) Type ```python3 WRAP.py``` (MacOS) or ```python3 .\WRAP.py``` (Windows) into the terminal

<p> How to Use WRAP </p>

<pp> For Single Object Search </pp>
1) Input the RA and DEC (in degrees) of the object being searched for
2) Input the radius (in arcsecs) around the RA and DEC that is being searched for
3) Put in the output file name that you want, do not put a file type at the end (it always saves to a csv file)
4) Select the catalogs that will be searched through
5) Click the green "Run Wrap" button
6) Your standard web browser will open with WISEView to aid you in finding your object
7) Click on the your object when the image appears (if your object is not there click the red "Object Not Found" button)
8) Repeat until every catalog image is shown
9) All data is saved to a csv file with your file output name in the "Output" folder (if no output name was put it will default to "WRAP_output.csv")

<pp> For Multi-Object Search </pp>
1) Click the "Browse" button and click on you file with a list of RA ("ra") and DEC's ("dec") (in degrees)
2) Input the radius (in arcsecs) around the RA and DEC that is being searched for
3) Click the "Filetype" browser and select the file type being used
4) Put in the output file name that you want, do not put a file type at the end
5) Select the catalogs that will be searched through
6) Click the green "Run Wrap" button
7) Your standard web browser will open with WISEView to aid you in finding your object (it will do this for every new object)
8) Click on the your object when the image appears (if your object is not there click the red "Object Not Found" button)
9) Repeat until every catalog image is shown and every object has been run though
10) All data is saved to a csv file with your file output name in the "Output" folder (if no output name was put it will default to "WRAP_output.csv")

<h2> 📞 Support </h2>

Mr. Hunter Brooks -> hcb98@nau.edu
Mr. Dan Caselden -> dancaselden@gmail.com

Slack Community: Email Mr. Hunter Brooks for a link to the slack community. 

<h2> 📖 Acknowledgments </h2>

1) If you were to publish any of the data gathered by WRAP please make sure that you are correctly acknowledging where the data comes from. Below is a thank you for all of the hard scientists that made all of the catalogs in WRAP. 
2)  Alonside using the correct acknowledgments for each catalog please cite Brook et al. (2023), in prep. when using WRAP for any publication. 
