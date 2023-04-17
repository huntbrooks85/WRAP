from catalog_module.importmodule import *

def blockPrint():
    '''Makes a function that blocks the printing function'''
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    '''Makes a function that allows the printing function'''
    sys.stdout = sys.__stdout__

def nsc_image(ra, dec, radius): 
    ''' First, it gets the images from the NSC API from AstroLab and downloads the images. 
    Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
    Third, calls the table function to get all the objects from the NSC source catalog. 
    Fourth, makes the window for the user to click the object with all settings. 
    Finally, finds the closest object to the click and records the data. '''
    
    #Makes outline for the window of the plot
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
        object_ra, object_dec = location_data['ra'].tolist(), location_data['dec'].tolist()
        g_list, g_list_e = location_data['gmag'].tolist(), location_data['gerr'].tolist()
        r_list, r_list_e = location_data['rmag'].tolist(), location_data['rerr'].tolist()
        i_list, i_list_e = location_data['imag'].tolist(), location_data['ierr'].tolist()
        z_list, z_list_e = location_data['zmag'].tolist(), location_data['zerr'].tolist()
        u_list, u_list_e = location_data['umag'].tolist(), location_data['uerr'].tolist()
        y_list, y_list_e = location_data['ymag'].tolist(), location_data['yerr'].tolist()
        pmra_list, pmra_list_e = location_data['pmra'].tolist(), location_data['pmraerr'].tolist()
        pmdec_list, pmdec_list_e = location_data['pmdec'].tolist(), location_data['pmdecerr'].tolist()

        #Obtains all of the image url, bands, and types from the table
        image_urls = imgTable['access_url']
        band_passes = imgTable['obs_bandpass']
        type = imgTable['prodtype']

        #Obtains the r band image
        for p in range(len(band_passes)):
            if band_passes[p] == 'g' and type[p] == 'image':
                image_url_g = image_urls[p]
                break

        #Obtains the r band image
        for p in range(len(band_passes)):
            if band_passes[p] == 'r' and type[p] == 'image':
                image_url_r = image_urls[p]
                break

        #Obtains the i band image
        for p in range(len(band_passes)):
            if band_passes[p] == 'i' and type[p] == 'image':
                image_url_i = image_urls[p]
                break

        #Obtains the z band image
        for p in range(len(band_passes)):
            if band_passes[p] == 'z' and type[p] == 'image':
                image_url_z = image_urls[p]
                break

        #Obtains the z band image
        for p in range(len(band_passes)):
            if band_passes[p] == 'Y' and type[p] == 'image':
                image_url_Y = image_urls[p]
                break

        #Downloads the images
        file_allwise_g, file_allwise_r, file_allwise_i = download_file(image_url_g, cache=True), download_file(image_url_r, cache=True), download_file(image_url_i, cache=True)
        file_allwise_z, file_allwise_Y = download_file(image_url_z, cache=True), download_file(image_url_Y, cache=True)
        data_allwise_g, data_allwise_r, data_allwise_i = fits.getdata(file_allwise_g), fits.getdata(file_allwise_r), fits.getdata(file_allwise_i)
        data_allwise_z, data_allwise_Y = fits.getdata(file_allwise_z), fits.getdata(file_allwise_Y)

        #Loads the WCS from the i band image
        hdu_g, hdu_r, hdu_i = fits.open(file_allwise_r)[0], fits.open(file_allwise_i)[0], fits.open(file_allwise_g)[0]
        hdu_z, hdu_Y = fits.open(file_allwise_z)[0], fits.open(file_allwise_Y)[0]
        wcs = WCS(hdu_g.header)

        #Makes the cutouts
        position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
        size = u.Quantity([radius, radius], u.arcsec)
        cutout_g = Cutout2D(data_allwise_g, position, size, fill_value = np.nan, wcs = wcs.celestial)
        cutout_r = Cutout2D(data_allwise_r, position, size, fill_value = np.nan, wcs = wcs.celestial)
        cutout_i = Cutout2D(data_allwise_i, position, size, fill_value = np.nan, wcs = wcs.celestial)
        cutout_z = Cutout2D(data_allwise_z, position, size, fill_value = np.nan, wcs = wcs.celestial)
        cutout_Y = Cutout2D(data_allwise_Y, position, size, fill_value = np.nan, wcs = wcs.celestial)
        wcs_cropped_i = cutout_i.wcs
        enablePrint()

        #Obtains the dates for each image
        date_g, date_r, date_i = hdu_g.header[36].split('T', 2)[0], hdu_r.header[36].split('T', 2)[0], hdu_i.header[36].split('T', 2)[0]
        date_z, date_Y = hdu_z.header[36].split('T', 2)[0], hdu_Y.header[36].split('T', 2)[0]

        #Allows the user to click the image to find the object
        location = []
        plt.rcParams["figure.figsize"] = [8, 8]
        plt.rcParams["figure.autolayout"] = True
        def mouse_event(event):
            '''Makes a list of the x, y, and axes the mouse click is.'''

            location.append(event.ydata)
            location.append(event.xdata)
            location.append(event.inaxes)
        fig_1 = plt.figure()
        cid = fig_1.canvas.mpl_connect('button_press_event', mouse_event)

        #Sets the WCS coordinates for the plots
        total_data = cutout_i.data + cutout_z.data + cutout_r.data
        ax = plt.subplot(projection = wcs_cropped_i)
        
        #Plots the objects found in the radius
        circle_size = (radius*3)
        scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

        #Normalize the image and plots it
        init_bot, init_top = 45, 95
        norm1_total = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, init_bot), vmax = np.nanpercentile(total_data.data, init_top))
        ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_total)

        #Makes the figure look pretty
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
        plt.suptitle('NSC Search', fontsize = 35, y = 0.96, fontfamily = 'Times New Roman')
        ax.set_title('Dates: \n'
                   + 'g Date: ' + str(date_g) + ' (Y-M-D)  ' + '  r Date: ' + str(date_r) + ' (Y-M-D)\n'
                   + 'i Date: ' + str(date_i) + ' (Y-M-D)  ' + '  z Date: ' + str(date_z) + ' (Y-M-D)\n'
                   + 'Y Date: ' + str(date_Y) + ' (Y-M-D)\n', fontdict = fontdict_1, y = 1.04)
        plt.grid(linewidth = 0)
        figure = plt.gcf()
        plt.xlim(len(total_data[0]), 0)
        figure.set_size_inches(4.75, 7.05)
        figure.canvas.set_window_title('NSC Search')

        #Make checkbuttons with all of the different image bands
        rax = plt.axes([0.045, 0.4, 0.12, 0.14])
        labels = ['g', 'r', 'i', 'z', 'Y']
        real_data = [cutout_g.data, cutout_r.data, cutout_i.data, cutout_z.data, cutout_Y.data]
        default = [False, True, True, True, False]
        check = CheckButtons(rax, labels, default)

        #Adds a slider for the scaling of the image
        freq_top = plt.axes([0.25, 0.12, 0.65, 0.03])
        slider_top = Slider(ax = freq_top, label = 'Top Stetch:', valmin = 50, valmax = 100, valinit = init_top, color = '#E48671')
        freq_bottom = plt.axes([0.25, 0.087, 0.65, 0.03])
        slider_bottom = Slider(ax = freq_bottom, label = 'Bottom Stetch:', valmin = 0, valmax = 50, valinit = init_bot, color = '#E48671')

        #Adds a slider for the circle size
        circle_slid_location = plt.axes([0.25, 0.055, 0.65, 0.03])
        circle_slider = Slider(ax = circle_slid_location, label = 'Circle Size:', valmin = (circle_size - 2.5*radius), valmax = (circle_size + 1*radius), valinit = circle_size, color = '#E48671')

        #Adds a notes section that the user can add notes about their data
        axbox = plt.axes([0.15, 0.02, 0.8, 0.03])
        text = ''
        text_box = TextBox(axbox, 'Notes:', initial = text, textalignment="center")

        #Make a button that can be clicked if no object is found
        axes_button = plt.axes([0.04, 0.767, 0.92, 0.04])
        close = Button(axes_button, 'Object Not Found', color = '#E48671')

        #Update the image depending on what the user chooses
        def update_button(label):
            '''Updates the list of activated images and updates the image the user can see.'''

            total_data = 0
            for lab in labels:
                if lab == label:
                    index = labels.index(lab)
                    if default[index] == False:
                        default[index] = True
                    elif default[index] == True: 
                        default[index] = False
            for d in range(len(default)):
                if default == [False, False, False, False, False]: 
                    total_data = real_data[0]*0
                if default[d] == True: 
                    total_data = total_data + real_data[d]
                else: 
                    pass
            norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
            ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

        #Updates the scaling when the slider is changed
        def update_slider_stretch(val):
            '''Updates the stretch the user can see, based in percentiles'''

            total_data = 0
            for d in range(len(default)):
                if default[d] == True: 
                    total_data = total_data + real_data[d]
                else: 
                    pass
            norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
            ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

        #Updates the notes added by the user when there is an input
        text_list = [text]
        def submit(expression):
            '''Updates the list of types in the 'Notes' setting'''

            text = expression
            text_list.append(text)
                
        #Allows the sliders and buttons to be pressed
        check.on_clicked(update_button)
        slider_top.on_changed(update_slider_stretch)
        slider_bottom.on_changed(update_slider_stretch)
        text_box.on_text_change(submit)

        #Display image until it is clicked to find the object
        n = -1
        while True:
            press = plt.waitforbuttonpress()
            text_max = len(text_list) - 1

            #Checks that it was a mouse click
            if press == False:
                n += 3
                
                #Finds which axes was clicked
                click_axes = str(location[n])
                click_axes = click_axes.split('WCSAxesSubplot', 2)[0]

                #Checks if the image was clicked
                if click_axes == '':
                    plt.close('all')
                    plt.figure().clear()

                    #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
                    coord = wcs_cropped_i.pixel_to_world_values(location[n-5],location[n-4])
                    distance = []
                    for i in range(len(object_ra)):
                        distance.append(math.dist(coord, [object_ra[i], object_dec[i]]))
                    list_location = distance.index(np.min(distance))
                    ra_nsc, dec_nsc = object_ra[list_location], object_dec[list_location]
                    g, g_e = g_list[list_location], g_list_e[list_location]
                    r, r_e = r_list[list_location], r_list_e[list_location]
                    i, i_e = i_list[list_location], i_list_e[list_location]
                    z, z_e = z_list[list_location], z_list_e[list_location]
                    u_mag, u_mag_e = u_list[list_location], u_list_e[list_location]
                    y, y_e = y_list[list_location], y_list_e[list_location]
                    pmra, pmra_e = pmra_list[list_location], pmra_list_e[list_location]
                    pmdec, pmdec_e = pmdec_list[list_location], pmdec_list_e[list_location]
                    return ra_nsc, dec_nsc, g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e, text_list[text_max]
                
                #Checks if the Object not Found button was clicked
                elif click_axes == 'Axes(0.04,0.012;0.92x0.04)':
                    g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                    ra_nsc, dec_nsc = ra, dec
                    plt.close('all')
                    plt.figure().clear()
                    return ra_nsc, dec_nsc, g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e, 'Object Not Found was Pressed'
                
                #Changes the circle size if the slider is changed
                elif click_axes == 'Axes(0.25,0.095;0.65x0.03)':
                    scatter.remove()
                    scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')
                
            #Checks if the window was closed
            elif press is None:
                g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                ra_nsc, dec_nsc = ra, dec
                plt.close('all')
                plt.figure().clear()
                return ra_nsc, dec_nsc, g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e, text_list[text_max]
            
    #Returns null values if image is not found
    else: 
        g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        ra_nsc, dec_nsc = ra, dec
        return ra_nsc, dec_nsc, g, g_e, r, r_e, i, i_e, z, z_e, u_mag, u_mag_e, y, y_e, pmra, pmra_e, pmdec, pmdec_e, 'Image Not Found'

def nsc_table(ra, dec, radius): 
    '''Find all the objects in the radius defined by the user'''

    blockPrint()

    #Makes a SQL query using the ra, dec, and radius
    query = " \
    SELECT ra, dec, gmag, gerr, rmag, rerr, imag, ierr, zmag, zerr, umag, uerr, ymag, yerr, pmra, pmraerr, pmdec, pmdecerr  \
    FROM nsc_dr2.object \
    WHERE ra > " + str((ra - (radius/7200))) + " and ra < " + str((ra + (radius/7200))) + " " \
    "AND dec > " + str((dec - (radius/7200))) + " and dec < " + str((dec + (radius/7200))) + " " \
    ""

    #Run this SQL quiery into the online NSC database
    response = qc.query(sql=query,format='csv')
    df = convert(response,'pandas')
    return df.head(100000000000)