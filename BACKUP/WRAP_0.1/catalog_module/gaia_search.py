from catalog_module.importmodule import *

#Makes a function that blocks the printing function
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

#Makes a function that allows the printing function
def enablePrint():
    sys.stdout = sys.__stdout__

#Does the GAIA search
def gaia_image(ra, dec, radius): 
    plt.rcParams['toolbar'] = 'None'
    plt.style.use('Solarize_Light2')
    blockPrint()
    #Finds all the metadata that relates to the ra and dec searched, mostly to find the APIs for the W1 and W2 images
    metadata_allwise_link = 'http://irsa.ipac.caltech.edu/ibe/sia/wise/allwise/p3am_cdd?POS=' + str(ra) + ',' + str(dec) + '&SIZE=' + str(radius/3600)
    allwise_metadata = requests.get(metadata_allwise_link)
    open('Output/metadata/gaia_metadata.txt', 'wb').write(allwise_metadata.content)

    #With this metadata it finds the API link for the W1 and W2 images
    w1_finder, w2_finder = 'W1 Coadd', 'W2 Coadd'
    with open('Output/metadata/gaia_metadata.txt', 'r') as fp:
        lines = fp.readlines()
        for line in lines:
            if line.find(w1_finder) != -1:
                w1_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]
            elif line.find(w2_finder) != -1:
                w2_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]

    #Download the W1, W2, W3, and W4 image
    file_allwise_w1 = download_file(w1_allwise_image_url, cache=True)
    data_allwise_w1 = fits.getdata(file_allwise_w1)
    file_allwise_w2 = download_file(w2_allwise_image_url, cache=True)
    data_allwise_w2 = fits.getdata(file_allwise_w2)

    #Find the location of all the object found in GAIA in the radius choosen by the user 
    location_data = gaia_table(ra, dec, radius)
    object_ra = location_data['ra'].tolist()
    object_dec = location_data['dec'].tolist()
    parallax_list = location_data['parallax'].tolist()
    g_list = location_data['phot_g_mean_mag'].tolist()
    rad_v_list = location_data['radial_velocity'].tolist()

    #Make a cutout from the coadd image for the RA and DEC put in
    hdu = fits.open(file_allwise_w1)[0]
    wcs1_w1 = WCS(hdu.header)

    position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
    size = u.Quantity([radius, radius], u.arcsec)
    cutout_w1 = Cutout2D(data_allwise_w1, position, size, fill_value = np.nan, wcs = wcs1_w1.celestial)
    wcs_cropped_w1 = cutout_w1.wcs
    cutout_w2 = Cutout2D(data_allwise_w2, position, size, fill_value = np.nan, wcs = wcs1_w1.celestial)

    #Defining a mouse click as an event on the plot
    location = []
    plt.rcParams["figure.figsize"] = [8, 8]
    plt.rcParams["figure.autolayout"] = True
    def mouse_event(event):
        location.append(event.ydata)
        location.append(event.xdata)
        location.append(event.inaxes)
    fig_1 = plt.figure()
    cid = fig_1.canvas.mpl_connect('button_press_event', mouse_event)

    #Sets the WCS coordinates for the plots
    w1_w2 = cutout_w1.data + cutout_w2.data
    ax = plt.subplot(projection = wcs_cropped_w1)

    #Plots the objects found in the radius
    plt.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = (radius*3), edgecolor='#40E842', facecolor='none')

    #Add a mouse hovering ability
    mplcursors.cursor(hover = 2, highlight = True)

    #Normalize the image and plots it
    norm1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(w1_w2.data, 45), vmax = np.nanpercentile(w1_w2.data, 98))
    plt.imshow(w1_w2.data, cmap = 'Greys', norm = norm1)

    #Makes the figure look pretty
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
    fontdict_1 = {'family':'Times New Roman','color':'k','size':20}
    plt.title('GAIA Search', fontdict = fontdict_1)
    plt.grid(linewidth = 0)
    figure = plt.gcf()
    figure.set_size_inches(5, 5.1) 

    #Make a button that can be clicked if no object is found
    axes = plt.axes([0.04, 0.012, 0.92, 0.04])
    close = Button(axes, 'Object Not Found', color = '#E48671')

    #Display image until it is clicked to find the object
    plt.waitforbuttonpress(0)
    plt.close('all')
    plt.figure().clear()

    #Find the closest point to the location clicked to obtain GAIA data
    if len(location) > 0: 
        if str(location[2]) == 'Axes(0.04,0.012;0.92x0.04)': 
            parallax, g_mag, rad_v = np.nan, np.nan, np.nan
            return parallax, g_mag, rad_v
        else: 
            coord = wcs_cropped_w1.pixel_to_world_values(location[0],location[1])
            distance = []
            for i in range(len(object_ra)):
                distance.append(math.dist(coord, [object_ra[i], object_dec[i]]))
            list_location = distance.index(np.min(distance))
            parallax = parallax_list[list_location]
            g_mag = g_list[list_location]
            rad_v = rad_v_list[list_location]
            return parallax, g_mag, rad_v

    elif len(location) == 0: 
        parallax, g_mag, rad_v = np.nan, np.nan, np.nan
        return parallax, g_mag, rad_v

def gaia_table(ra, dec, radius): 
    blockPrint()
    #Makes the SQL code to run it into the GAIA search
    query = "SELECT TOP 2000 \
    gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.parallax,gaia_source.pmra,gaia_source.pmdec,gaia_source.ruwe,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.phot_variable_flag,gaia_source.non_single_star,gaia_source.has_xp_continuous,gaia_source.has_xp_sampled,gaia_source.has_rvs,gaia_source.has_epoch_photometry,gaia_source.has_epoch_rv,gaia_source.has_mcmc_gspphot,gaia_source.has_mcmc_msc,gaia_source.teff_gspphot,gaia_source.logg_gspphot,gaia_source.mh_gspphot,gaia_source.distance_gspphot,gaia_source.azero_gspphot,gaia_source.ag_gspphot,gaia_source.ebpminrp_gspphot \
    FROM gaiadr3.gaia_source \
    WHERE \
    CONTAINS( \
	POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec), \
	CIRCLE('ICRS', " + str(ra) + "," + str(dec) + "," + str((radius/2) * 0.000277778)+ ")" \
    ")=1"

    #Run this SQL quiery into the online GAIA database
    job = Gaia.launch_job_async(query)
    results = job.get_results()
    return results