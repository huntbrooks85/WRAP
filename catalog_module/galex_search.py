from catalog_module.importmodule import *

#Makes a function that blocks the printing function
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

#Makes a function that allows the printing function
def enablePrint():
    sys.stdout = sys.__stdout__

#Does the GALEX search
def galex_image(ra, dec, radius_use): 
    plt.rcParams['toolbar'] = 'None'
    plt.style.use('Solarize_Light2')

    blockPrint()
    #Obtains all of the observations from MAST
    obs_table = Observations.query_region(str(ra) + ' ' + str(dec), radius=(radius_use/3600))
    obs_table = pd.DataFrame(data = np.array(obs_table))

    #Finds only the GALEX images
    project = obs_table['project'].tolist()
    for p in range(len(project)):
        if project[p] == 'AIS':
            pass
        else: 
            obs_table = obs_table.drop(p)

    #Finds the real images from GALEX
    nuv_urls = obs_table['dataURL'].tolist()
    for n in range(len(nuv_urls)):
        temp_url = nuv_urls[n].split('-', 7)
        if temp_url[len(temp_url) - 1] == 'int.fits.gz':
            image_url = nuv_urls[n]
            break
        else: 
            image_url = ''
            pass

    if image_url != '': 
        #Find the location of all the object found in GALEX in the radius choosen by the user 
        observation_table = galex_table(ra, dec, radius_use)
        object_ra = observation_table['ra'].tolist()
        object_dec = observation_table['dec'].tolist()
        fuv_list = observation_table['fuv_mag'].tolist()
        nuv_list = observation_table['nuv_mag'].tolist()

        #Downloads the images
        file_allwise_nuv = download_file(image_url, cache=True)
        data_allwise_nuv = fits.getdata(file_allwise_nuv)

        #Make a cutout from the coadd image for the RA and DEC put in
        hdu_w1 = fits.open(file_allwise_nuv)[0]
        wcs1_w1 = WCS(hdu_w1.header)

        position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
        size = u.Quantity([radius_use, radius_use], u.arcsec)
        cutout_w1 = Cutout2D(data_allwise_nuv, position, size, fill_value = np.nan, wcs = wcs1_w1.celestial)
        wcs_cropped_w1 = cutout_w1.wcs

        enablePrint()
        #Defining a mouse click as an event on the plot
        location = []
        plt.rcParams["figure.figsize"] = [8, 8]
        plt.rcParams["figure.autolayout"] = True
        def mouse_event(event):
            location.append(event.ydata)
            location.append(event.xdata)
            location.append(event.inaxes)
        plt.connect('button_press_event', mouse_event)

        #Sets the WCS coordinates for the plots
        ax = plt.subplot(projection = wcs_cropped_w1)

        #Plots the objects found in the radius
        plt.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = (radius_use*3), edgecolor='#40E842', facecolor='none')

        #Add a mouse hovering ability
        mplcursors.cursor(hover = 2, highlight = True)

        #Normalize the image and plots it
        norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(cutout_w1.data, 45), vmax = np.nanpercentile(cutout_w1.data, 95))
        plt.imshow(cutout_w1.data, cmap = 'Greys', norm = norm1_w1)

        #Makes the figure look pretty
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        fontdict_1 = {'family':'Times New Roman','color':'k','size':20}
        plt.title('GALEX Search', fontdict = fontdict_1)
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

        #Outputs null values if not clicked
        if len(location) == 0: 
            nuv, fuv = np.nan, np.nan
            plt.close('all')
            plt.figure().clear()
            return nuv, fuv

        #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
        if str(location[2]) == 'Axes(0.04,0.012;0.92x0.04)': 
            nuv, fuv = np.nan, np.nan
            return nuv, fuv
        else: 
            coord = wcs_cropped_w1.pixel_to_world_values(location[0],location[1])
            distance = []
            for i in range(len(object_ra)):
                distance.append(math.dist(coord, [object_ra[i], object_dec[i]]))
            list_location = distance.index(np.min(distance))
            fuv = fuv_list[list_location]
            nuv = nuv_list[list_location]
            return nuv, fuv
    else: 
        nuv, fuv = np.nan, np.nan
        return nuv, fuv

def galex_table(ra, dec, radius): 
    blockPrint()
    #Gets the table of all of the data in the radius from the user
    catalog_data = Catalogs.query_region(str(ra) + ' ' + str(dec), radius=(radius/3600), catalog="Galex", table="mean")
    return catalog_data