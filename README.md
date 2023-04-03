<h1 align="center" id="title">✨ WRAP ✨</h1>

<p id="description">Gathers photometry and astrometry from various ultra-violet, optical, and near-infrared catalogs given a RA, DEC, and Radius by the user. Completed by the user clicking their object overplotted on the catalog images and gathering all of the photometry and astrometry from those various catalogs.</p>
  
<h2>🔭 Catalogs </h2>

*   CatWISE2020 & AllWISE -> Wide-Field Infrared Survey Explorer https://irsa.ipac.caltech.edu/Missions/wise.html
*   Gaia ->GAIA ESA https://gea.esac.esa.int/archive/
*   VISTA -> Visible and Infrared Survey Telescope for Astronomy http://horus.roe.ac.uk/vsa/index.html
*   WFCAM -> Wide Field Camera http://wsa.roe.ac.uk/index.html
*   2MASS -> Two Micron All Sky Survey https://irsa.ipac.caltech.edu/Missions/2mass.html
*   PanSTARRS -> Panoramic Survey Telescope and Rapid Response System https://outerspace.stsci.edu/display/PANSTARRS/
*   NSC -> NOIRLab Source Catalog https://datalab.noirlab.edu/nscdr2/index.php
*   GALEX -> Galaxy Evolution Explorer http://www.galex.caltech.edu

<h2>🛠️ Installation Steps:</h2>

<p>1. Python </p>

```
-> Python 3.8.8 is needed for WRAP 
-> Python 3.8.8 link: https://www.python.org/downloads/release/python-388/
```

<p>2. Packages </p>

```
-> pyvo == 1.4 
-> numpy = 1.22.0
-> pandas == 1.5.3 
-> astropy == 5.2.1
-> truncate == 0.11 
-> requests == 2.28.1 
-> astroquery == 0.4.6 
-> matplotlib == 3.5.0 
-> PySimpleGUI == 4.60.4 
-> OpenCV-python == 1.22.0 
-> astro-datalab == 2.20.1 
-> beautifulsoup4 == 4.11.1
```

<h2> 🏆 Opening Application: </h2>
*   Only tested on MacOS Ventura 13.2.1, problems may occur for older versions of MacOS and versions for Windows and Linux will be made in the future. 
*   Note that the window close button has been disabled, to close WRAP please click the red "Close WRAP" button at the buttom. 

<p> How to Start Program </p>

1) Open terminal 
2) Go to the file directory that contains "WRAP.py" (example: "cd /Users/hunter_brooks8/Desktop/WRAP_Project/WRAP")
3) Type "python WRAP.py" into the terminal

<p> How to Use WRAP </p>

<pp> For Single Object Search </pp>
1) Input the RA and DEC (in degrees) of the object being searched for
2) Input the radius (in arcsecs) around the RA and DEC that is being searched for
3) Put in the output file name that you want, do not put a file type at the end (it always saves to a csv file)
4) Select the catalogs that will be searched through
5) Click the green "Run Wrap" button
6) Click on the your object when the image appears (if your object is not there click the red "Object Not Found" button)
7) Repeat until every catalog image is shown
8) All data is saved to a csv file with your file output name in the "Output" folder (if no output name was put it will default to "WRAP_output.csv")

<pp> For Multi-Object Search </pp>
1) Click the "Browse" button and click on you file with a list of RA ("ra") and DEC's ("dec") (in degrees)
2) Input the radius (in arcsecs) around the RA and DEC that is being searched for
3) Click the "Filetype" browser and select the file type being used
4) Put in the output file name that you want, do not put a file type at the end
5) Select the catalogs that will be searched through
6) Click the green "Run Wrap" button
7) Click on the your object when the image appears (if your object is not there click the red "Object Not Found" button)
8) Repeat until every catalog image is shown and every object has been run though
9) All data is saved to a csv file with your file output name in the "Output" folder (if no output name was put it will default to "WRAP_output.csv")


<h2> Variable Names </h2> 
*   Photometry Filter Profiles: http://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse

 <p> CatWISE2020: </p>

```
	-> cw_ra: The Right Accenceion (RA) from CatWISE2020
	-> cw_dec: The Declination (DEC) from CatWISE2020
	-> cw_w1 & cw_w1_e: The W1 band from CatWISE2020 and its uncertainties
	-> cw_w2 & cw_w2_e: The W2 band from CatWISE2020 and its uncertainties
	-> cw_pmra & cw_pmra_e: The proper motion in RA from CatWISE2020 and its uncertainties
	-> cw_pmdec & cw_pmdec_e: The proper motion in DEC from CatWISE2020 and its uncertainties
	-> cw_notes: The notes written in the CatWISE2020 by the user

```
 <p> AllWISE: </p>
 ```
	-> aw_ra: The Right Accenceion (RA) from AllWISE
	-> aw_dec: The Declination (DEC) from AllWISE
	-> aw_w1 & aw_w1_e: The W1 band from AllWISE and its uncertainties
	-> aw_w2 & aw_w2_e: The W2 band from AllWISE and its uncertainties
	-> aw_w3 & aw_w3_e: The W3 band from AllWISE and its uncertainties
	-> aw_w4 & aw_w4_e: The W4 band from AllWISE and its uncertainties
	-> aw_pmra & aw_pmra_e: The proper motion in RA from AllWISE and its uncertainties
	-> aw_pmdec & aw_pmdec_e: The proper motion in DEC from AllWISE and its uncertainties
	-> aw_notes: The notes written in the AllWISE by the user
```
 <p> Gaia: </p>
 ```
	-> gaia_ra: The Right Accenceion (RA) from Gaia
	-> gaia_dec: The Declination (DEC) from Gaia
	-> gaia_parallax & gaia_parallax_e: The parallax from Gaia and its uncertainties
	-> gaia_radv & gaia_radv_e: The radial velocity from Gaia and its uncertainties
	-> gaia_pmra & gaia_pmra_e: The proper motion in RA from Gaia and its uncertainties
	-> gaia_pmdec & gaia_pmdec_e: The proper motion in DEC from Gaia and its uncertainties
	-> gaia_g: The g band from Gaia and its uncertainties
	-> gaia_bp: The bp band from Gaia and its uncertainties
	-> gaia_rp: The rp band from Gaia and its uncertainties
	-> gaia_notes: The notes written in the Gaia by the user
```
 <p> VISTA: </p>
 ```
	-> vsa_ra: The Right Accenceion (RA) from VISTA
	-> vsa_dec: The Declination (DEC) from VISTA
	-> vsa_y & vsa_y_e: The Y band from VISTA and its uncertainties
	-> vsa_j & vsa_j_e: The J band from VISTA and its uncertainties
	-> vsa_h & vsa_h_e: The H band from VISTA and its uncertainties
	-> vsa_ks & vsa_ks_e: The Ks band from VISTA and its uncertainties
	-> vsa_notes: The notes written in the VISTA by the user
```
 <p> WFCAM: </p>
 ```
	-> wfcam_ra: The Right Accenceion (RA) from WFCAM
	-> wfcam_dec: The Declination (DEC) from WFCAM
	-> wfcam_y & wfcam_y_e: The Y band from WFCAM and its uncertainties
	-> wfcam_j & wfcam_j_e: The J band from WFCAM and its uncertainties
	-> wfcam_h & wfcam_h_e: The H band from WFCAM and its uncertainties
	-> wfcam_k & wfcam_k_e: The K band from WFCAM and its uncertainties
	-> wfcam_pmra & wfcam_pmra_e: The proper motion in RA from WFCAM and its uncertainties
	-> wfcam_pmdec & wfcam_pmdec_e: The proper motion in DEC from WFCAM and its uncertainties
	-> wfcam_notes: The notes written in the WFCAM by the user
```
 <p> 2MASS: </p>
 ```
	-> 2mass_ra: The Right Accenceion (RA) from 2MASS
	-> 2mass_dec: The Declination (DEC) from 2MASS
	-> 2mass_j & 2mass_j_e: The J band from 2MASS and its uncertainties
	-> 2mass_h & 2mass_h_e: The J band from 2MASS and its uncertainties
	-> 2mass_ks & 2mass_ks_e: The Ks band from 2MASS and its uncertainties
	-> 2mass_notes: The notes written in the 2MASS by the user
```
 <p> PanSTARRS: </p>
 ```
	-> ps_ra: The Right Accenceion (RA) from PanSTARRS
	-> ps_dec: The Declination (DEC) from PanSTARRS
	-> ps_g & ps_g_e: The g band from PanSTARRS and its uncertainties
	-> ps_r & ps_r_e: The r band from PanSTARRS and its uncertainties
	-> ps_i & ps_i_e: The i band from PanSTARRS and its uncertainties
	-> ps_z & ps_z_e: The z band from PanSTARRS and its uncertainties
	-> ps_y & ps_y_e: The y band from PanSTARRS and its uncertainties
	-> ps_notes: The notes written in the PanSTARRS by the user
```
 <p> NSC: </p>
 ```
	-> nsc_ra: The Right Accenceion (RA) from NSC
	-> nsc_dec: The Declination (DEC) from NSC
	-> nsc_g & nsc_g_e: The g band from NSC and its uncertainties
	-> nsc_r & nsc_r_e: The r band from NSC and its uncertainties
	-> nsc_i & nsc_i_e: The i band from NSC and its uncertainties
	-> nsc_z & nsc_z_e: The z band from NSC and its uncertainties
	-> nsc_u & nsc_u_e: The u band from NSC and its uncertainties
	-> nsc_pmra & nsc_pmra_e: The proper motion in RA from NSC and its uncertainties
	-> nsc_pmedc & nsc_pmdec_e: The proper motion in DEC from NSC and its uncertainties
	-> nsc_notes: The notes written in the NSC by the user
```
 <p> GALEX: </p>
 ```
	-> galex_ra: The Right Accenceion (RA) from GALEX
	-> galex_dec: The Declination (DEC) from GALEX
	-> galex_fuv & galex_fuv_e: The FUV band from GALEX and its uncertainties
	-> galex_nuv & galex_nuv_e: The NUV band from GALEX and its uncertainties
	-> galex_notes: The notes written in the GALEX by the user
```
<h2> 📞 Support </h2>

Hunter Brooks -> hcb98@nau.edu