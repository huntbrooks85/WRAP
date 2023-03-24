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

        #Downloads the i, z, and r bands images
        file_allwise_g = download_file(image_url_g, cache=True)
        data_allwise_g = fits.getdata(file_allwise_g)
        file_allwise_r = download_file(image_url_r, cache=True)
        data_allwise_r = fits.getdata(file_allwise_r)
        file_allwise_i = download_file(image_url_i, cache=True)
        data_allwise_i = fits.getdata(file_allwise_i)
        file_allwise_z = download_file(image_url_z, cache=True)
        data_allwise_z = fits.getdata(file_allwise_z)
        file_allwise_Y = download_file(image_url_Y, cache=True)
        data_allwise_Y = fits.getdata(file_allwise_Y)

        #Loads the WCS from the i band image
        hdu_g = fits.open(file_allwise_g)[0]
        wcs1_g = WCS(hdu_g.header)
        hdu_r = fits.open(file_allwise_r)[0]
        hdu_i = fits.open(file_allwise_i)[0]
        hdu_z = fits.open(file_allwise_z)[0]
        hdu_Y = fits.open(file_allwise_Y)[0]

        #Makes the cutouts
        position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
        size = u.Quantity([radius, radius], u.arcsec)
        cutout_g = Cutout2D(data_allwise_g, position, size, fill_value = np.nan, wcs = wcs1_g.celestial)
        cutout_r = Cutout2D(data_allwise_r, position, size, fill_value = np.nan, wcs = wcs1_g.celestial)
        cutout_i = Cutout2D(data_allwise_i, position, size, fill_value = np.nan, wcs = wcs1_g.celestial)
        cutout_z = Cutout2D(data_allwise_z, position, size, fill_value = np.nan, wcs = wcs1_g.celestial)
        cutout_Y = Cutout2D(data_allwise_Y, position, size, fill_value = np.nan, wcs = wcs1_g.celestial)
        wcs_cropped_i = cutout_i.wcs

        enablePrint()
        date_g = hdu_g.header[36].split('T', 2)[0]
        time_g = hdu_g.header[36].split('T', 2)[1]
        date_r = hdu_r.header[36].split('T', 2)[0]
        time_r = hdu_r.header[36].split('T', 2)[1]
        date_i = hdu_i.header[36].split('T', 2)[0]
        time_i = hdu_i.header[36].split('T', 2)[1]
        date_z = hdu_z.header[36].split('T', 2)[0]
        time_z = hdu_z.header[36].split('T', 2)[1]
        date_Y = hdu_Y.header[36].split('T', 2)[0]
        time_Y = hdu_Y.header[36].split('T', 2)[1]
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
        total_data = cutout_i.data + cutout_z.data + cutout_r.data
        ax = plt.subplot(projection = wcs_cropped_i)
        
        #Plots the objects found in the radius
        circle_size = (radius*3)
        scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

        #Normalize the image and plots it
        init_top = 95
        init_bot = 45
        norm1_total = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, init_bot), vmax = np.nanpercentile(total_data.data, init_top))
        ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_total)

        #Makes the figure look pretty
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
        ax.set_title(     'Dates: \n'
                        + 'g Date: ' + str(date_g) + ' (Y/M/D)  ' + '  r Date: ' + str(date_r) + ' (Y/M/D)\n'
                        + 'i Date: ' + str(date_i) + ' (Y/M/D)  '+ '  z Date: ' + str(date_z) + ' (Y/M/D)\n'
                        + 'Y Date: ' + str(date_Y) + ' (Y/M/D)\n'
                        + 'Times: \n'
                        + 'g Time: ' + str(time_g) + ' (UT)  ' + '  r Time: ' + str(time_r) + ' (UT)\n'
                        + 'i Time: ' + str(time_i) + ' (UT)  '+ '  z Time: ' + str(time_z) + ' (UT)\n'
                        + 'Y Time: ' + str(time_Y), fontdict = fontdict_1, y = 1)
        plt.grid(linewidth = 0)
        figure = plt.gcf()
        figure.set_size_inches(4.75, 6.95)
        figure.canvas.set_window_title('NSC Search')

    #Make checkbuttons with all of the different image bands
    rax = plt.axes([0.045, 0.4, 0.105, 0.12])
    labels = ['g', 'r', 'i', 'z', 'Y']
    real_data = [cutout_g.data, cutout_r.data, cutout_i.data, cutout_z.data, cutout_Y.data]
    default = [False, True, True, True, False]
    check = CheckButtons(rax, labels, default)

    #Adds a slider for the scaling of the image
    freq_top = plt.axes([0.25, 0.155, 0.65, 0.03])
    slider_top = Slider(ax = freq_top, label = 'Top Stetch:', valmin = 50, valmax = 100, valinit = init_top, color = '#E48671')
    freq_bottom = plt.axes([0.25, 0.125, 0.65, 0.03])
    slider_bottom = Slider(ax = freq_bottom, label = 'Bottom Stetch:', valmin = 0, valmax = 50, valinit = init_bot, color = '#E48671')

    #Adds a slider for the circle size
    circle_slid_location = plt.axes([0.25, 0.095, 0.65, 0.03])
    circle_slider = Slider(ax = circle_slid_location, label = 'Circle Size:', valmin = (circle_size - 2.5*radius), valmax = (circle_size + 1*radius), valinit = circle_size, color = '#E48671')

    #Adds a notes section that the user can add notes about their data
    axbox = plt.axes([0.25, 0.06, 0.65, 0.03])
    text = ''
    text_box = TextBox(axbox, 'Notes:', initial = text, textalignment="center")

    #Make a button that can be clicked if no object is found
    axes_button = plt.axes([0.04, 0.012, 0.92, 0.04])
    close = Button(axes_button, 'Object Not Found', color = '#E48671')

    #Update the image depending on what the user chooses
    def update_button(label):
        total_data = 0
        for lab in labels:
            if lab == label:
                index = labels.index(lab)
                if default[index] == False:
                    default[index] = True
                elif default[index] == True: 
                    default[index] = False
        for d in range(len(default)):
            if default[d] == True: 
                total_data = total_data + real_data[d]
            else: 
                pass
        norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
        ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

    #Updates the scaling when the slider is changed
    def update_slider_stretch(val):
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
                g = g_list[list_location]
                r = r_list[list_location]
                i = i_list[list_location]
                z = z_list[list_location]
                return g, r, i, z, text_list[text_max]
            
            #Checks if the Object not Found button was clicked
            elif click_axes == 'Axes(0.04,0.012;0.92x0.04)':
                g, r, i, z = np.nan, np.nan, np.nan, np.nan
                plt.close('all')
                plt.figure().clear()
                return g, r, i, z, text_list[text_max]
            
            elif click_axes == 'Axes(0.25,0.095;0.65x0.03)':
                scatter.remove()
                scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')
            
        #Checks if the window was closed
        elif press is None:
            g, r, i, z = np.nan, np.nan, np.nan, np.nan
            plt.close('all')
            plt.figure().clear()
            return g, r, i, z, text_list[text_max]

def nsc_table(ra, dec, radius): 
    blockPrint()
    #Makes a SQL query using the ra, dec, and radius
    query = " \
    SELECT ra, dec, gmag, rmag, imag, zmag \
    FROM nsc_dr1.object \
    WHERE ra > " + str((ra - (radius/7200))) + " and ra < " + str((ra + (radius/7200))) + " " \
    "AND dec > " + str((dec - (radius/7200))) + " and dec < " + str((dec + (radius/7200))) + " " \
    ""

    #Run this SQL quiery into the online NSC database
    response = qc.query(sql=query,format='csv')
    df = convert(response,'pandas')
    return df.head(100000000000)