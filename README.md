<p align="center">
    <a href="https://ibb.co/9WGmwZV"><img src="/Icons/LOGO.png" width="50%"></a> <br>
</p>

<h1 align="center" id="title"> 💫✨ WRAP ✨💫 </h1>
<div align="center">

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10359982.svg)](https://doi.org/10.5281/zenodo.10359982)

<p id="description"> The <b>Wide-field Retrieval of Astrodata Program (WRAP)</b> is a versatile tool designed to simplify the retrieval of photometric and astrometric data from a wide range of ultraviolet, optical, and near-infrared catalogs. Users input Right Ascension (RA) and Declination (DEC) coordinates, along with a search radius for their object of interest. The program then enables users to interactively select their object by simply clicking on it within catalog images. WRAP efficiently compiles and presents all pertinent photometric and astrometric information obtained from the selected catalogs, making it accessible in a convenient CSV file. This streamlines the data retrieval process, enhancing the user's experience</p>
</div>

<p align="center">
<div align="center">
  <h2>⬇️ Installing WRAP ⬇️</h2>
</div>

1.  **Download Python3.8.8**  
   Visit the [Python3.8.8](https://www.python.org/downloads/release/python-388/) webpage and download Python3.8.8. 

2.  **Download Java**  
   Visit the [Java SE Development Kit](https://www.oracle.com/java/technologies/downloads/) webpage and download Java SE Development Kit.

3. **Download WRAP**  
   Visit the [Releases](https://github.com/huntbrooks85/WRAP/releases) section on the right side of the GitHub page and download the latest version of WRAP.

4. **Open Terminal and Navigate to the WRAP Directory**  
   Use the terminal to move into the directory where WRAP was downloaded. For example:

   ```bash
   cd ~/Documents/GitHub/WRAP
   ```

5. **Run the Build Script**  
   Execute the build script with the following command:

   ```bash
   ./build_wrap.sh
   ```

   > **Note:** During this process, you may be prompted to install dependencies such as Python 3.8 or Homebrew. Please follow the instructions provided in the terminal if prompted.

6. **Install the WRAP Application**  
   After the script completes, a file named `WRAP-v2.1.0.pkg` will be created in the WRAP directory. Double-click this package to begin the installation process.

7. **Launch WRAP**  
   Once installed, WRAP will appear in your **Applications** folder. You can now launch it like any standard macOS application.
</p>


<div align="center">
  <h2>🏆 Using WRAP 🏆</h2>
</div>

<div align="center">
  <p><b>How to Use WRAP</b></p>
</div>
<div align="center">
</div>

<div align="left">
  <p><b>Single-Query:</b></p>
</div>

1. Enter the **Right Ascension (R.A.)** and **Declination (Decl.)** of the object in **decimal degrees**.
2. Specify the **search radius** in arcseconds around the R.A./Decl.
3. Select the **catalogs** you wish to query.
4. Click the **"Run"** button.
5. If the **WiseView** option is selected, your default browser will open WiseView to assist in identifying the object.
6. When the image appears, **click on your object**.  
   If it’s not visible, click the **red "Object Not Found"** button.
7. Repeat this process until all catalog images are reviewed.
8. Your query results will appear in the **"Table List"** panel on the left.
9. Click **"Save"** to export the results as a CSV file.

<div align="left">
  <p><b>Multi-Query:</b></p>
</div>

1. Check the **"Multi-Query"** checkbox.
2. Upload a **CSV file**:
   - R.A. column must be named `input_ra`
   - Decl. column must be named `input_dec`
3. Enter the **search radius** in arcseconds.
4. Select the **catalogs** to search.
5. Click the **"Run"** button.
6. If the **WiseView** option is selected, your default browser will open it.
7. When the image appears, **click on your object**, or use the **"Object Not Found"** button if it's missing.
8. Repeat the process for each catalog image.
9. Query results will be listed in the **"Table List"** panel. Click **"Save"** to export them as a CSV file.

<div align="center">
  <p><b>⚠️ Significant Details ⚠️</b></p>
</div>

- **Note 1**: Only tested on MacOS >11 and Windows >8; problems may occur for older versions of MacOS and Windows ***(NOT SUPPORTED ON LINUX)***.
- **Note 2**: Windows does not support the astro-datalab package; therefore, Windows does not have the Noirlab Source Catalog option.
- **Note 3**: The orientation for all of the catalogs is North pointed up and East pointing left.
- **Note 4**: Catalog imaging may have strange croppings; this is a warning that it may happen, do not be alarmed by it this is simplfy cropping at the edge of the detector.

<div align="center">
  <h2>📞 Support 📞</h2>
</div>

- **Mr. Hunter Brooks**
  - Email: hcb98@nau.edu

<div align="center">
  <h2>📖 Acknowledgments 📖</h2>
</div>

1. If you intend to publish any of the data gathered by WRAP, please ensure that you correctly acknowledge the sources of the data. 

2. In addition to using the correct acknowledgments for each catalog, please cite [Brook et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023RNAAS...7..272B/abstract) when using WRAP for any publication.

