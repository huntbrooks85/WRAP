#Import all of the packages
from catalog_module.importmodule import *

def blockPrint():
    '''Makes a function that blocks the printing function'''
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    '''Makes a function that allows the printing function'''
    sys.stdout = sys.__stdout__

def galex_image(ra, dec, radius_use): 
    ''' First, it gets the images from the GALEX API from MAST and downloads the images. 
    Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
    Third, calls the table function to get all the objects from the GALEX source catalog. 
    Fourth, makes the window for the user to click the object with all settings. 
    Finally, finds the closest object to the click and records the data. '''

    #Makes outline for the window of the plot
    plt.rcParams['toolbar'] = 'None'
    plt.style.use('Solarize_Light2')
    matplotlib.use("TkAgg")

    blockPrint()
    #Obtains all of the observations from MAST
    obs_table = Observations.query_region(str(ra) + ' ' + str(dec), radius = (radius_use/1800))
    obs_table = pd.DataFrame(data = np.array(obs_table))

    #Finds only the GALEX images
    telescope = obs_table['obs_collection'].tolist()
    for p in range(len(telescope)):
        if telescope[p] == 'GALEX':
            pass
        else: 
            obs_table = obs_table.drop(p)

    #Finds the real images from GALEX
    urls = obs_table['dataURL'].tolist()
    obs_table = pd.DataFrame(data = np.array(obs_table))
    for n in range(len(urls)):
        temp_url = urls[n].split('-', 7)
        if temp_url[len(temp_url) - 1] == 'int.fits.gz':
            pass
        else: 
            obs_table = obs_table.drop(n)
    
    #Finds all of the NUV band urls
    nuv_list_url = []
    urls = obs_table[27].tolist()
    bands = obs_table[5].tolist()
    for n in range(len(bands)):
        if bands[n] == 'NUV': 
            nuv_list_url.append(urls[n])
        
    #Runs the code if an image is found
    if len(nuv_list_url) > 0: 
        #Find the location of all the object found in GALEX in the radius choosen by the user 
        observation_table = galex_table(ra, dec, radius_use)
        object_ra, object_dec = observation_table['ra'].tolist(), observation_table['dec'].tolist()
        fuv_list, fuv_list_e = observation_table['fuv_mag'].tolist(), observation_table['fuv_magerr'].tolist()
        nuv_list, nuv_list_e = observation_table['nuv_mag'].tolist(), observation_table['nuv_magerr'].tolist()

        #Combines all of the real images
        total_data = 0
        for i in range(len(nuv_list_url)):
            #Downloads the images
            file_allwise_nuv = download_file(nuv_list_url[i], cache=True)
            data_allwise_nuv = fits.getdata(file_allwise_nuv)

            #Make a cutout from the coadd image for the RA and DEC put in
            hdu_w1 = fits.open(file_allwise_nuv)[0]
            wcs1_w1 = WCS(hdu_w1.header)

            #Makes a cutout for the images
            position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
            size = u.Quantity([radius_use, radius_use], u.arcsec)
            cutout_w1 = Cutout2D(data_allwise_nuv, position, size, fill_value = np.nan, wcs = wcs1_w1.celestial)

            #Combines the data
            total_data = cutout_w1.data + total_data

            #Gets the new WCS coordinates for the image cutout
            wcs_cropped_w1 = cutout_w1.wcs

        #Finds the date for the last image
        date_w1 = hdu_w1.header['OBS-DATE']
        enablePrint()

        #Defining a mouse click as an event on the plot
        location = []
        plt.rcParams["figure.figsize"] = [8, 8]
        plt.rcParams["figure.autolayout"] = True
        def mouse_event(event):
            '''Makes a list of the x, y, and axes the mouse click is.'''

            location.append(event.ydata)
            location.append(event.xdata)
            location.append(event.inaxes)
        plt.connect('button_press_event', mouse_event)

        #Sets the WCS coordinates for the plots
        ax = plt.subplot(projection = wcs_cropped_w1)

        #Plots the objects found in the radius
        circle_size = (radius_use*4)
        scatter = plt.scatter(object_ra, object_dec, transform = ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

        #Normalize the image and plots it
        init_top, init_bot = 95, 45
        norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data, init_bot), vmax = np.nanpercentile(total_data, init_top))
        plt.imshow(total_data, cmap = 'Greys', norm = norm1_w1)

        #Makes the figure look pretty
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        fontdict = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
        plt.suptitle('GALEX Search', fontsize = 35, y = 0.95, fontfamily = 'Times New Roman')
        ax.set_title('Dates: \n'
                   + 'FUV Date: ' + str(date_w1) + ' (Y-M-D)  ', fontdict = fontdict, y = 1.105)
        plt.grid(linewidth = 0)
        figure = plt.gcf()
        plt.xlim(0, max(total_data.shape))
        plt.ylim(0, max(total_data.shape))
        if platform != 'win32':
            figure.set_size_inches(4.75, 7.25)
        elif platform == 'win32':
            figure.set_size_inches(4.75, 7.25)
        figure.canvas.set_window_title('GALEX Search')
        mng = pyplot.get_current_fig_manager()
        mng.window.resizable(False, False)

        #Adds a slider for the scaling of the image
        freq_top = plt.axes([0.25, 0.12, 0.65, 0.03])
        slider_top = Slider(ax = freq_top, label = 'Top Stetch:', valmin = 50, valmax = 100, valinit = init_top, color = '#E48671')
        freq_bottom = plt.axes([0.25, 0.087, 0.65, 0.03])
        slider_bottom = Slider(ax = freq_bottom, label = 'Bottom Stetch:', valmin = 0, valmax = 50, valinit = init_bot, color = '#E48671')

        #Adds a slider for the circle size
        circle_slid_location = plt.axes([0.25, 0.055, 0.65, 0.03])
        circle_slider = Slider(ax = circle_slid_location, label = 'Circle Size:', valmin = (circle_size - 2.5*radius_use), valmax = (circle_size + 1*radius_use), valinit = circle_size, color = '#E48671')

        #Adds a notes section that the user can add notes about their data
        axbox = plt.axes([0.15, 0.02, 0.8, 0.03])
        text = ''
        text_box = TextBox(axbox, 'Notes:', initial = text, textalignment="center")

        #Make a button that can be clicked if no object is found
        axes_button = plt.axes([0.04, 0.775, 0.92, 0.04])
        close = Button(axes_button, 'Object Not Found', color = '#E48671')

        #Updates the scaling when the slider is changed
        def update_slider_stretch(val):
            '''Updates the stretch the user can see, based in percentiles'''

            norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, slider_bottom.val), vmax = np.nanpercentile(total_data.data, slider_top.val))
            ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

        #Updates the notes added by the user when there is an input
        text_list = [text]
        def submit(expression):
            '''Updates the list of types in the 'Notes' setting'''

            text = expression
            text_list.append(text)

        #Allows the sliders and buttons to be pressed
        slider_top.on_changed(update_slider_stretch)
        slider_bottom.on_changed(update_slider_stretch)
        text_box.on_text_change(submit)

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
                    #Makes a pop-up window with success text
                    plt.clf()
                    plt.close('all')
                    plt.figure(1)
                    plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for GALEX! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
                    plt.xlim(0, 1)
                    plt.ylim(0, 1)
                    plt.grid(linewidth = 0)
                    ax = plt.gca()
                    ax.xaxis.set_tick_params(labelbottom=False)
                    ax.yaxis.set_tick_params(labelleft=False)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    figure2 = plt.gcf()
                    figure2.set_size_inches(4.75, 1)
                    figure2.canvas.set_window_title('Successful GALEX Search')
                    mng2 = pyplot.get_current_fig_manager()
                    mng2.window.resizable(False, False)
                    plt.pause(0.1)
                    plt.clf()
                    plt.close('all')

                    #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
                    coord = wcs_cropped_w1.pixel_to_world_values(location[n-4], location[n-5])
                    distance = []
                    for i in range(len(object_ra)):
                        distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))
                    list_location = distance.index(np.min(distance))
                    ra_galex, dec_galex = object_ra[list_location], object_dec[list_location]
                    fuv, fuv_e = fuv_list[list_location], fuv_list_e[list_location]
                    nuv, nuv_e = nuv_list[list_location], nuv_list_e[list_location]
                    return ra_galex, dec_galex, fuv, fuv_e, nuv, nuv_e, 'MAST GALEX Data Release 7', text_list[text_max]
                
                #Checks if the "Object Not Found" button was clicked
                elif click_axes == 'Axes(0.04,0.775;0.92x0.04)':
                    fuv, fuv_e, nuv, nuv_e = np.nan, np.nan, np.nan, np.nan
                    ra_galex, dec_galex = ra, dec

                    #Makes a pop-up window with success text
                    plt.clf()
                    plt.close('all')
                    plt.figure(1)
                    plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for GALEX! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
                    plt.xlim(0, 1)
                    plt.ylim(0, 1)
                    plt.grid(linewidth = 0)
                    ax = plt.gca()
                    ax.xaxis.set_tick_params(labelbottom=False)
                    ax.yaxis.set_tick_params(labelleft=False)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    figure2 = plt.gcf()
                    figure2.set_size_inches(4.75, 1)
                    figure2.canvas.set_window_title('Successful GALEX Search')
                    mng2 = pyplot.get_current_fig_manager()
                    mng2.window.resizable(False, False)
                    plt.pause(0.1)
                    plt.clf()
                    plt.close('all')
        
                    return ra_galex, dec_galex, fuv, fuv_e, nuv, nuv_e, 'MAST GALEX Release Data 7', 'Object Not Found was Pressed'
                
                #Updates the circle size when slider is moved
                elif click_axes == 'Axes(0.25,0.055;0.65x0.03)':
                    scatter.remove()
                    scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')

            #Checks if the window was closed
            elif press is None:
                fuv, fuv_e, nuv, nuv_e = np.nan, np.nan, np.nan, np.nan
                ra_galex, dec_galex = ra, dec
                plt.close('all')
                plt.figure().clear()
                return ra_galex, dec_galex, fuv, fuv_e, nuv, nuv_e, 'MAST GALEX Data Release 7', text_list[text_max]
                
    #Returns null values if the images are not found    
    else: 
        fuv, fuv_e, nuv, nuv_e = np.nan, np.nan, np.nan, np.nan
        ra_galex, dec_galex = ra, dec
        return ra_galex, dec_galex, fuv, fuv_e, nuv, nuv_e, 'MAST GALEX Data Release 7', 'Image Not Found'

def galex_table(ra, dec, radius): 
    '''Find all the objects in the radius defined by the user'''
    
    blockPrint()

    #Gets the table of all of the data in the radius from the user
    catalog_data = Catalogs.query_region(str(ra) + ' ' + str(dec), radius = ((radius/7200) - 0.000277778), catalog = "Galex", table = "mean")
    return catalog_data