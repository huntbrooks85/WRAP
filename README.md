<p align="center">
    <a href="https://ibb.co/9WGmwZV"><img src="/Output/metadata/LOGO.png" width="50%"></a> <br>
</p>

<h1 align="center" id="title"> 💫✨ WRAP ✨💫 </h1>

<p id="description"> The Wide-field Retrieval of Astrodata Program (WRAP) is a tool designed to aid the retrieval of photometric and astrometric data from a diverse range of ultraviolet, optical, and near-infrared catalogs. Users provide WRAP with Right Ascension (RA), Declination (DEC) coordinates and a search radius for their object of interest. The program then allows users to interactively select their object by clicking on it within catalog images. WRAP subsequently compiles and presents all relevant photometric and astrometric information gathered from the chosen catalogs in a CSV file, streamlining the data retrieval process.</p>

<div align="center">
  <h2>🛠️ Installation 🛠️</h2>
</div>

<div align="center">
<p><b>Program Setup</b></p>
</div>

<pp><b>🐍 Anaconda Environment 🐍</b><pp>
**Note**: Making use of an Anaconda (Conda) environment is highly recommended due to potential compatibility issues with older Python and pip packages when using WRAP.

1. **Download Anaconda Navigator**: To get started, please visit [Anaconda's official website](https://www.anaconda.com) and download the version of Anaconda Navigator that matches your operating system.

2. **Access Anaconda Environments**: After downloading and installing Anaconda Navigator, launch the application. In the Navigator's main window, you will find an "Environments" option on the left. Click on it to access your Conda environments.

3. **Create a New Conda Environment**: In the "Environments" section, you will see a plus sign labeled "Create." Click on it to create a new Conda environment.

4. **Customize Environment Settings**: Give your new environment a name of your choice (e.g., "WRAP"). Next, select a Python version that is closest to Python 3.8 (e.g., 3.8.18). Then, click the "Create" button to create the environment.

5. **Environment Creation**: The environment creation process may take a while, as Conda will install the necessary packages and dependencies for WRAP.

6. **Activate Your Conda Environment**: Once the environment is successfully created, you can activate it in your terminal. Open your terminal and enter the following command, replacing `*conda name*` with the name you chose for your environment:
   ```bash
   conda activate *conda name*

<pp><b>PIP Installations</b><pp>
How to use PIP

<div align="center">
<h2> 🏆 Opening Application 🏆 </h2>
</div>

Type ```python3 WRAP.py``` (MacOS) or ```python3 .\WRAP.py``` (Windows) into the terminal

<div align="center">
<p><b>How to Use WRAP</b></p>
</div>

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

<div align="center">
<p><b>Significant Details</b></p>
</div>

*   Note 1: Only tested on MacOS >13 and Windows >10, problems may occur for older versions of MacOS and Windows (IS NOT SUPPORTED ON LINUX)
*   Note 2: Windows does not support the astro-datalab package, therefore Windows does not have the Noirlab Source Catalog option. 
*   Note 3: The window close button has been disabled, to close WRAP please click the red "Close WRAP" button at the bottom. 
*   Note 4: The orientation for all of the catalogs is North pointed up and East pointing left. 
*   Note 5: 2MASS may have strange imaging cropping, this is a warning that it may happen and do not be alarmed by it. 


<h2> 🌓 Sky Coverage 🌗 </h2>
</div>
<p align="center">
    <a href="https://ibb.co/GQkGgNN"><img src="/Output/metadata/Coverage_Map-1.png" width="100%"></a> <br>
    Sky-Coverage for the Catalogs Included in WRAP
</p>

<div align="center">
<h2> 📞 Support 📞 </h2>
</div>

Mr. Hunter Brooks -> hcb98@nau.edu
Mr. Dan Caselden -> dancaselden@gmail.com

Slack Community: Email Mr. Hunter Brooks for a link to the slack community. 

<div align="center">
<h2> 📖 Acknowledgments 📖 </h2>
</div>

1) If you were to publish any of the data gathered by WRAP please make sure that you are correctly acknowledging where the data comes from. Below is a thank you for all of the hard scientists that made all of the catalogs in WRAP. 
2)  Alonside using the correct acknowledgments for each catalog please cite Brook et al. (2023), in prep. when using WRAP for any publication. 
