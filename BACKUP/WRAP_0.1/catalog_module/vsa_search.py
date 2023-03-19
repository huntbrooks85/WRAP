from catalog_module.importmodule import *

#Makes a function that blocks the printing function
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

#Makes a function that allows the printing function
def enablePrint():
    sys.stdout = sys.__stdout__

def vsa_image(ra, dec, radius_use): 
    plt.rcParams['toolbar'] = 'None'
    plt.style.use('Solarize_Light2')
    blockPrint()
    #Obtains all of the urls in J, H, and K from VSA
    url_J, url_H, url_Ks = [
        Vsa.get_image_list(
            SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'),
                     image_width = radius_use * u.arcsec,
                     waveband = band,
                     database = 'VHSDR6')
        for band in ['J', 'H', 'Ks']]

    #Checking to see if the images exist
    if len(url_J) > 0 and len(url_H) > 0 and len(url_Ks) > 0: 
        #Downloading the fits images
        file_vsa_J = download_file(url_J[0], cache=True)
        data_vsa_J = fits.getdata(file_vsa_J)
        file_vsa_H = download_file(url_H[0], cache=True)
        data_vsa_H = fits.getdata(file_vsa_H)
        file_vsa_Ks = download_file(url_Ks[0], cache=True)
        data_vsa_Ks = fits.getdata(file_vsa_Ks)

        #Find the location of all the object found in VSA in the radius choosen by the user
        catalog_list = ['VHS', 'VVV', 'VMC', 'VIKING', 'VIDEO', 'UltraVISTA']
        for cat in catalog_list: 
            table = Vsa.query_region(
                SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'),
                     radius = (radius_use/2) * u.arcsec,
                     programme_id = cat,
                     database = 'VHSDR6')
            if len(table) > 0:
                break
        if len(table) == 0: 
            j, h, k = np.nan, np.nan, np.nan
            return j, h, k

        object_ra = table['ra'].tolist()
        object_dec = table['dec'].tolist()
        J_list = table['jAperMag3'].tolist()
        H_list = table['hAperMag3'].tolist()
        Ks_list = table['ksAperMag3'].tolist()

        #Make a cutout from the coadd image for the RA and DEC put in
        hdu_j = fits.open(file_vsa_J)[1]
        wcs1_j = WCS(hdu_j.header)

        position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5')
        size = u.Quantity([radius_use, radius_use], u.arcsec)
        cutout_j = Cutout2D(data_vsa_J, position, size, wcs = wcs1_j.celestial)
        cutout_h = Cutout2D(data_vsa_H, position, size, wcs = wcs1_j.celestial)
        cutout_ks = Cutout2D(data_vsa_Ks, position, size, wcs = wcs1_j.celestial)
        wcs_cropped_j = cutout_j.wcs

        enablePrint()
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
        j_h_k = cutout_j.data + cutout_h.data + cutout_ks.data
        ax = plt.subplot(projection = wcs_cropped_j)

        #Plots the objects found in the radius
        plt.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = (radius_use*3), edgecolor='#40E842', facecolor='none')

        #Add a mouse hovering ability
        mplcursors.cursor(hover = 2, highlight = True)

        #Normalize the image and plots it
        norm1_j = matplotlib.colors.Normalize(vmin = np.nanpercentile(j_h_k.data, 45), vmax = np.nanpercentile(j_h_k.data, 98))
        plt.imshow(j_h_k.data, cmap = 'Greys', norm = norm1_j)

        #Makes the figure look pretty
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        fontdict_1 = {'family':'Times New Roman','color':'k','size':20}
        plt.title('VSA Search', fontdict = fontdict_1)
        plt.grid(linewidth = 0)
        ax.invert_yaxis()
        figure = plt.gcf()
        figure.set_size_inches(5, 5.1) 

        #Make a button that can be clicked if no object is found
        axes = plt.axes([0.04, 0.012, 0.92, 0.04])
        close = Button(axes, 'Object Not Found', color = '#E48671')

        #Display image until it is clicked to find the object
        plt.waitforbuttonpress(0)
        plt.close('all')
        plt.figure().clear()

        #Find the closest point to the location clicked to obtain J, H, K photometry
        if len(location) > 0:
            if str(location[2]) == 'Axes(0.04,0.012;0.92x0.04)': 
                j, h, k = np.nan, np.nan, np.nan
                return j, h, k
            else: 
                coord = wcs1_j.pixel_to_world_values(location[0],location[1])
                distance = []
                for i in range(len(object_ra)):
                    distance.append(math.dist(coord, [object_ra[i], object_dec[i]]))
                list_location = distance.index(np.min(distance))
                j = J_list[list_location] 
                h = H_list[list_location]
                k = Ks_list[list_location]
                return j, h, k

        elif len(location) == 0: 
            j, h, k = np.nan, np.nan, np.nan
            return j, h, k

    else: 
            j, h, k = np.nan, np.nan, np.nan
            return j, h, k
