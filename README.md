<h1 align="center" id="title">✨ WRAP ✨</h1>

<p id="description"> The Wide-field Retrieval of Astrodata Program (WRAP) is a tool designed to aid the retrieval of photometric and astrometric data from a diverse range of ultraviolet, optical, and near-infrared catalogs. Users provide WRAP with Right Ascension (RA), Declination (DEC) coordinates and a search radius for their object of interest. The program then allows users to interactively select their object by clicking on it within catalog images. WRAP subsequently compiles and presents all relevant photometric and astrometric information gathered from the chosen catalogs in a CSV file, streamlining the data retrieval process.</p>

<h2> 🛠️ Installation Steps 🛠️ </h2>

```
-> Python 3.8.8 is Required
```
* Python 3.8.8 link: https://www.python.org/downloads/release/python-388/
* Python Installation Guide: https://wiki.python.org/moin/BeginnersGuide/Download

<h2> 🏆 Opening Application 🏆 </h2>
<p><b>Program Setup</b></p>

1) Open terminal 
2) Go to the file directory that contains "WRAP.py" (example: "cd /Users/hunter/Desktop/WRAP")
3) Activate relavent virtual environment by typing ```source wrap_mac/bin/activate``` (MacOS) or ```wrap_windows\Scripts\activate``` (Windows)
4) Windows may return the error: ```\wrap_windows\Scripts\activate cannot be
loaded because running scripts is disabled on this system. For more information, see about_Execution_Policies at```, run the command ```Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass``` to fix the error
1) Type ```python3 WRAP.py``` (MacOS) or ```python3 .\WRAP.py``` (Windows) into the terminal

<p><b>How to Use WRAP</b></p>

<pp><i>Single Object Search</i></pp>
1) Input the RA and DEC (in degrees) of the object being searched for
2) Input the radius (in arcsecs) around the RA and DEC that is being searched for
3) Put in the output file name that you want, do not put a file type at the end (it always saves to a csv file)
4) Select the catalogs that will be searched through
5) Click the green "Run Wrap" button
6) Your standard web browser will open with WISEView to aid you in finding your object
7) Click on the your object when the image appears (if your object is not there click the red "Object Not Found" button)
8) Repeat until every catalog image is shown
9) All data is saved to a csv file with your file output name in the "Output" folder (if no output name was put it will default to "WRAP_output.csv")

<pp><i>Multi-Object Search</i></pp>
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

<p><b>Significant Details</b></p>

*   Note 1: Only tested on MacOS >13 and Windows >10, problems may occur for older versions of MacOS and Windows (IS NOT SUPPORTED ON LINUX)
*   Note 2: Windows does not support the astro-datalab package, therefore Windows does not have the Noirlab Source Catalog option. 
*   Note 3: The window close button has been disabled, to close WRAP please click the red "Close WRAP" button at the bottom. 
*   Note 4: The orientation for all of the catalogs is North pointed up and East pointing left. 
*   Note 5: 2MASS may have strange imaging cropping, this is a warning that it may happen and do not be alarmed by it. 


<h2> 🌓 Sky Coverage 🌗 </h2>
<p align="center">
    <a href="https://ibb.co/1LGq00T"><img src="/Output/metadata/Coverage_Map-1.png" width="75%"></a> <br>
    Sky-Coverage for the Catalogs Included in WRAP
</p>

<h2> 📞 Support 📞 </h2>

Mr. Hunter Brooks -> hcb98@nau.edu
Mr. Dan Caselden -> dancaselden@gmail.com

Slack Community: Email Mr. Hunter Brooks for a link to the slack community. 

<h2> 📖 Acknowledgments 📖 </h2>

1) If you were to publish any of the data gathered by WRAP please make sure that you are correctly acknowledging where the data comes from. Below is a thank you for all of the hard scientists that made all of the catalogs in WRAP. 
2)  Alonside using the correct acknowledgments for each catalog please cite Brook et al. (2023), in prep. when using WRAP for any publication. 
