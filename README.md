<p align="center">
    <a href="https://ibb.co/9WGmwZV"><img src="/Source_Code/Output/metadata/LOGO.png" width="50%"></a> <br>
</p>

<h1 align="center" id="title"> 💫✨ WRAP ✨💫 </h1>
<div align="center">

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10359982.svg)](https://doi.org/10.5281/zenodo.10359982)

<div align="center">
  <b><a href="https://hcb985.wixsite.com/wrap-byw">WRAP Website</a></b>
  <br>
</div>

<p id="description"> The <b>Wide-field Retrieval of Astrodata Program (WRAP)</b> is a versatile tool designed to simplify the retrieval of photometric and astrometric data from a wide range of ultraviolet, optical, and near-infrared catalogs. Users input Right Ascension (RA) and Declination (DEC) coordinates, along with a search radius for their object of interest. The program then enables users to interactively select their object by simply clicking on it within catalog images. WRAP efficiently compiles and presents all pertinent photometric and astrometric information obtained from the selected catalogs, making it accessible in a convenient CSV file. This streamlines the data retrieval process, enhancing the user's experience</p>
</div>

<div align="center">
  <h2>🛠️ Installation For Source Code 🛠️</h2>
</div>

<div align="center">
  <b>THIS IS ONLY FOR THE SOURCE CODE, NOT THE MACOS/WINDOWS APPLICATION</b>
</div>

<div align="center">
<pp><b>🦾 Installing Git 🦾</b><pp>
</div>
<div align="center">
<pp><b>-----------------------------------------</b><pp>
</div>

1. **Download Git**: To get started, you'll need to install Git, which is a version control system. You can download the official Git installer for your operating system from [Git's official website](https://git-scm.com/downloads).

2. **Verification**: After installation, you can verify that Git is properly installed by opening a terminal or command prompt and running the following command:  
```bash
   git --version
   ```

<div align="center">
<pp><b>🐍 Anaconda Environment 🐍</b><pp>
</div>
<div align="center">
<pp><b>-----------------------------------------</b><pp>
</div>

**Note**: Making use of an Anaconda (Conda) environment is highly recommended due to potential compatibility issues with older Python and pip packages when using WRAP.

1. **Download Anaconda Navigator**: To get started, please visit [Anaconda's official website](https://www.anaconda.com) and download the version of Anaconda Navigator that matches your operating system.

2. **Access Anaconda Environments**: After downloading and installing Anaconda Navigator, launch the application. In the Navigator's main window, you will find an "Environments" option on the left. Click on it to access your Conda environments.

3. **Create a New Conda Environment**: In the "Environments" section, you will see a plus sign labeled "Create." Click on it to create a new Conda environment.

4. **Customize Environment Settings**: Give your new environment a name of your choice (e.g., "WRAP"). Next, select a Python version that is closest to Python 3.8 (e.g., 3.8.18). Then, click the "Create" button to create the environment.

5. **Environment Creation**: The environment creation process may take a while, as Conda will install the necessary packages and dependencies for WRAP.

6. **Activate Your Conda Environment**: Once the environment is successfully created, you can activate it in your terminal. Open your terminal and enter the following command, replacing `*conda name*` with the name you chose for your environment:
   ```bash
   conda activate *conda name*

<div align="center">
  <p><b>⬇️ PIP Installation ⬇️</b></p>
</div>
<div align="center">
<pp><b>-----------------------------------------</b><pp>
</div>

<p><i>MacOS</i></p>

1. Once the conda environment is activated, go to your directory containing the WRAP contents (e.g., "cd Documents/GitHub/WRAP").
2. Enter the "pip_module" directory (e.g., "cd pip_module")
3. Run the command:
   ```bash
   pip install -r requirements-mac.txt

<ppp><i>Windows</i><ppp>
1. Once the conda environment is activated, go to your directory containing the WRAP contents (e.g., "cd Documents/GitHub/WRAP").
2. Enter the "pip_module" directory (e.g., "cd pip_module")
3. Run the command:
   ```bash
   pip install -r requirements-windows.txt


<div align="center">
<pp><b> Opening WRAP </b><pp>
</div>
<div align="center">
<pp><b>-----------------------------------------</b><pp>
</div>

To open the application, type the following command into your terminal:

- For MacOS: `python3 WRAP.py`
- For Windows: `python3 .\WRAP.py`

<div align="center">
  <h2>🏆 Using WRAP 🏆</h2>
</div>

<div align="center">
  <b><a href="https://hcb985.wixsite.com/wrap-byw">WRAP Website</a></b>
</div>

<div align="center">
  <p><b>How to Use WRAP</b></p>
</div>
<div align="center">
<pp><b>-----------------------------------------</b><pp>
</div>

**Single Object Search**

1. Input the Right Ascension (RA) and Declination (DEC) of the object being searched for in degrees.
2. Input the search radius in arcseconds around the RA and DEC.
3. Specify the desired output file name (without a file type). WRAP will save results to a CSV file by default.
4. Select the catalogs to be searched through.
5. Click the green "Run Wrap" button.
6. Your standard web browser will open with WISEView to assist you in finding your object.
7. Click on your object when the image appears (if your object is not there, click the red "Object Not Found" button).
8. Repeat these steps until every catalog image is shown.
9. All data is saved to a CSV file with your specified output name in the "Output" folder. If no output name was provided, it will default to "WRAP_output.csv."

**Multi-Object Search**

1. Click the "Browse" button and select your file with a list of RA ("ra") and DEC ("dec") coordinates in degrees.
2. Input the search radius in arcseconds around the RA and DEC.
3. Choose the file type used.
4. Specify the output file name (without a file type).
5. Select the catalogs to be searched through.
6. Click the green "Run Wrap" button.
7. Your standard web browser will open with WISEView to assist you in finding your object (it will do this for every new object).
8. Click on your object when the image appears (if your object is not there, click the red "Object Not Found" button).
9. Repeat these steps until every catalog image is shown, and every object has been processed.
10. All data is saved to a CSV file with your specified output name in the "Output" folder. If no output name was provided, it will default to "WRAP_output.csv."

<div align="center">
  <p><b>⚠️ Significant Details ⚠️</b></p>
  <p>-----------------------------------------</p>
</div>

- **Note 1**: Only tested on MacOS >11 and Windows >8; problems may occur for older versions of MacOS and Windows ***(NOT SUPPORTED ON LINUX)***.
- **Note 2**: Windows does not support the astro-datalab package; therefore, Windows does not have the Noirlab Source Catalog option.
- **Note 3**: The window close button has been disabled; to close WRAP, please click the red "Close WRAP" button at the bottom.
- **Note 4**: The orientation for all of the catalogs is North pointed up and East pointing left.
- **Note 5**: 2MASS may have strange imaging cropping; this is a warning that it may happen, so do not be alarmed by it.



<div align="center">
  <h2>🌓 Sky Coverage 🌗</h2>
</div>

<p align="center">
  <a href="https://ibb.co/GQkGgNN"><img src="/Source_Code/Output/metadata/Coverage_Map-1.png" width="100%"></a> <br>
  Sky-Coverage for the Catalogs Included in WRAP
</p>

<div align="center">
  <h2>📞 Support 📞</h2>
</div>

- **Mr. Hunter Brooks**
  - Email: hcb98@nau.edu

<div align="center">
  <h2>📖 Acknowledgments 📖</h2>
</div>

1. If you intend to publish any of the data gathered by WRAP, please ensure that you correctly acknowledge the sources of the data. 

2. In addition to using the correct acknowledgments for each catalog, please cite "Brook et al. (2023), in prep." when using WRAP for any publication.

