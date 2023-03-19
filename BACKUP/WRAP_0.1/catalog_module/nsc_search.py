from catalog_module.importmodule import *

#Makes a function that blocks the printing function
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

#Makes a function that allows the printing function
def enablePrint():
    sys.stdout = sys.__stdout__

#Does the NSC search
def nsc_image(ra, dec, radius): 
    plt.rcParams['toolbar'] = 'None'
    plt.style.use('Solarize_Light2')
    blockPrint()

    #Defines the catalog that is searched
    DEF_ACCESS_URL = "https://datalab.noirlab.edu/sia/des_dr1"
    svc = sia.SIAService(DEF_ACCESS_URL)

    #Finds all of the image urls for the ra, dec, and radius given
    imgTable = svc.search((ra,dec), (radius/3600)).to_table()

    #Tests if any images were found
    if len(imgTable) > 0:
        #Runs the table function to obtain the g, r, i, and z bands
        location_data = nsc_table(ra, dec, radius)
        object_ra = location_data['ra'].tolist()
        object_dec = location_data['dec'].tolist()
        g_list = location_data['gmag'].tolist()
        r_list = location_data['rmag'].tolist()
        i_list = location_data['imag'].tolist()
        z_list = location_data['zmag'].tolist()

        #Obtains all of the image url, bands, and types from the table
        image_urls = imgTable['access_url']
        band_passes = imgTable['obs_bandpass']
        type = imgTable['prodtype']

        #Obtains the i band image
        for p in range(len(band_passes)):
            if band_passes[p] == 'i' and type[p] == 'image':
                image_url_i = image_urls[p]
                break

        #Obtains the r band image
        for p in range(len(band_passes)):
            if band_passes[p] == 'r' and type[p] == 'image':
                image_url_r = image_urls[p]
                break

        #Obtains the z band image
        for p in range(len(band_passes)):
            if band_passes[p] == 'z' and type[p] == 'image':
                image_url_z = image_urls[p]
                break

        #Downloads the i, z, and r bands images
        file_allwise_i = download_file(image_url_i, cache=True)
        data_allwise_i = fits.getdata(file_allwise_i)
        file_allwise_z = download_file(image_url_z, cache=True)
        data_allwise_z = fits.getdata(file_allwise_z)
        file_allwise_r = download_file(image_url_r, cache=True)
        data_allwise_r = fits.getdata(file_allwise_r)

        #Loads the WCS from the i band image
        hdu = fits.open(file_allwise_i)[0]
        wcs1_i = WCS(hdu.header)

        #Makes the cutouts
        position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
        size = u.Quantity([radius, radius], u.arcsec)
        cutout_i = Cutout2D(data_allwise_i, position, size, fill_value = np.nan, wcs = wcs1_i.celestial)
        cutout_z = Cutout2D(data_allwise_z, position, size, fill_value = np.nan, wcs = wcs1_i.celestial)
        cutout_r = Cutout2D(data_allwise_r, position, size, fill_value = np.nan, wcs = wcs1_i.celestial)
        wcs_cropped_i = cutout_i.wcs

        enablePrint()
        #Allows the user to click the image to find the object
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
        r_i_z = cutout_i.data + cutout_z.data + cutout_r.data
        ax = plt.subplot(projection = wcs_cropped_i)
        
        #Plots the objects found in the radius
        plt.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = (radius*3), edgecolor='#40E842', facecolor='none')

        #Add a mouse hovering ability
        mplcursors.cursor(hover = 2, highlight = True)

        #Normalize the image and plots it
        norm1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(r_i_z.data, 45), vmax = np.nanpercentile(r_i_z.data, 98))
        plt.imshow(r_i_z.data, cmap = 'Greys', norm = norm1)

        #Makes the figure look pretty
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        fontdict_1 = {'family':'Times New Roman','color':'k','size':20}
        plt.title('NOIRLab Source Catalog Search', fontdict = fontdict_1)
        plt.grid(linewidth = 0)
        ax.invert_xaxis()
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
                g, r, i, z = np.nan, np.nan, np.nan, np.nan
                return g, r, i, z
            else: 
                coord = wcs_cropped_i.pixel_to_world_values(location[0],location[1])
                distance = []
                for i in range(len(object_ra)):
                    distance.append(math.dist(coord, [object_ra[i], object_dec[i]]))
                list_location = distance.index(np.min(distance))
                g = g_list[list_location]
                r = r_list[list_location]
                i = i_list[list_location]
                z = z_list[list_location]
                return g, r, i, z

        elif len(location) == 0: 
            g, r, i, z = np.nan, np.nan, np.nan, np.nan
            return g, r, i, z
    else: 
        g, r, i, z = np.nan, np.nan, np.nan, np.nan
        return g, r, i, z

def nsc_table(ra, dec, radius): 
    blockPrint()
    #Makes a SQL query using the ra, dec, and radius
    query = " \
    SELECT ra, dec, gmag, rmag, imag, zmag \
    FROM nsc_dr1.object \
    WHERE ra > " + str((ra - (radius/3600))) + " and ra < " + str((ra + (radius/3600))) + " " \
    "AND dec > " + str((dec - (radius/3600))) + " and dec < " + str((dec + (radius/3600))) + " " \
    ""

    #Run this SQL quiery into the online NSC database
    response = qc.query(sql=query,format='csv')
    df = convert(response,'pandas')
    return df.head(100000000000)