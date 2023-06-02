#Import all of the packages
from catalog_module.importmodule import *

def blockPrint():
    '''Makes a function that blocks the printing function'''
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    '''Makes a function that allows the printing function'''
    sys.stdout = sys.__stdout__

def ukidss_image(ra, dec, radius): 
    ''' First, it gets the images from the WFCAM API from the WFCAM Archive and downloads the images. 
    Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
    Third, calls the table function to get all the objects from the WFCAM source catalog. 
    Fourth, makes the window for the user to click the object with all settings. 
    Finally, finds the closest object to the click and records the data. '''
      
    #Makes outline for the window of the plot
    plt.rcParams['toolbar'] = 'None'
    matplotlib.use("TkAgg")
    plt.style.use('Solarize_Light2')
    blockPrint()
    
    #Obtains all of the urls in J, H, and K from UKIDSS
    database_list = ['UKIDSSDR11PLUS', 'UHSDR1']
    catalo = ['LAS', 'UHSDR1']
    for data in range(len(database_list)):
        url_J = Ukidss.get_image_list(SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'), image_width = (radius) * u.arcsec, image_height = (radius) * u.arcsec, waveband = 'J', database = database_list[data], programme_id = catalo[data])
        url_H = Ukidss.get_image_list(SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'), image_width = (radius) * u.arcsec, image_height = (radius) * u.arcsec, waveband = 'H', database = database_list[data], programme_id = catalo[data])
        url_K = Ukidss.get_image_list(SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'), image_width = (radius) * u.arcsec, image_height = (radius) * u.arcsec, waveband = 'K', database = database_list[data], programme_id = catalo[data])
        url_Y = Ukidss.get_image_list(SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'), image_width = (radius) * u.arcsec, image_height = (radius) * u.arcsec, waveband = 'Y', database = database_list[data], programme_id = catalo[data])
        if len(url_J) > 0: 
            data = database_list[data]
            break

    #Checking to see if the images exist
    if len(url_J) > 0: 

        #Calls the UKIDSS Table
        table = ukidss_table(ra, dec, radius)
        
        #Defines each variables depending on if the image was found in UKIDSS DR11 or UHS DR1
        if data == 'UHSDR1': 

            #Downloading the fits images
            file_ukidss_J = download_file(url_J[0], cache = True) 
            data_ukidss_J = fits.getdata(file_ukidss_J)

            #Obtains the headers from the images
            hdu_j = fits.open(file_ukidss_J)[1]
            wcs_j = WCS(hdu_j.header)

            #Obtains the dates for each image
            date_j = str(hdu_j.header['HISTORY']).split(' ', 2)[0]

            #Gets the columns from the table
            object_ra, object_dec = table['ra'].tolist(), table['dec'].tolist()
            J_list, J_list_e = table['jAperMag3'].tolist(), table['jAperMag3Err'].tolist()
        
        else:
            #Downloading the fits images
            file_ukidss_J, file_ukidss_H, file_ukidss_K, file_ukidss_Y = download_file(url_J[0], cache=True), download_file(url_H[0], cache=True), download_file(url_K[0], cache=True), download_file(url_Y[0], cache=True)
            data_ukidss_J, data_ukidss_H, data_ukidss_K, data_ukidss_Y = fits.getdata(file_ukidss_J), fits.getdata(file_ukidss_H), fits.getdata(file_ukidss_K), fits.getdata(file_ukidss_Y)

            #Obtains the headers from the images
            hdu_j, hdu_h, hdu_k, hdu_y = fits.open(file_ukidss_J)[1], fits.open(file_ukidss_H)[1], fits.open(file_ukidss_K)[1], fits.open(file_ukidss_Y)[1]
            wcs_j, wcs_h, wcs_k, wcs_y = WCS(hdu_j.header),  WCS(hdu_h.header), WCS(hdu_k.header), WCS(hdu_y.header)

            #Make a cutout from the coadd image for the RA and DEC put in
            position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5')
            size = u.Quantity([radius, radius], u.arcsec)
            cutout_j = Cutout2D(data_ukidss_J, position, size, wcs = wcs_j.celestial)
            cutout_h = Cutout2D(data_ukidss_H, position, size, wcs = wcs_h.celestial)
            cutout_k = Cutout2D(data_ukidss_K, position, size, wcs = wcs_k.celestial)
            cutout_y = Cutout2D(data_ukidss_Y, position, size, wcs = wcs_y.celestial)
            wcs_cropped = cutout_j.wcs

            #Resizes the H, K, and Y bands to that of J
            j_shape = (cutout_j.data).shape
            if j_shape[0] != j_shape[1]:
                cutout_j = cutout_j.data[0:min(j_shape) - 1, 0:min(j_shape) - 1]
                j_shape = (cutout_j.data).shape
            h_reshape, k_reshape, y_reshape = cv2.resize(cutout_h.data, dsize = j_shape, interpolation = cv2.INTER_NEAREST), cv2.resize(cutout_k.data, dsize = j_shape, interpolation = cv2.INTER_NEAREST), cv2.resize(cutout_y.data, dsize = j_shape, interpolation = cv2.INTER_NEAREST)
            total_data = cutout_j.data + y_reshape

            #Obtains the dates for each image
            date_y, date_j, date_h, date_k = str(hdu_y.header['HISTORY']).split(' ', 2)[0], str(hdu_j.header['HISTORY']).split(' ', 2)[0], str(hdu_h.header['HISTORY']).split(' ', 6)[2], str(hdu_k.header['HISTORY']).split(' ', 6)[2]

            #Gets the columns from the table
            object_epoch = table['epoch'].tolist()
            object_ra, object_dec = table['ra'].tolist(), table['dec'].tolist()
            object_ra_e, object_dec_e = table['sigRa'].tolist(), table['sigDec'].tolist()
            Y_list, Y_list_e = table['yAperMag3'].tolist(), table['yAperMag3Err'].tolist()
            J_list, J_list_e = table['jAperMag3'].tolist(), table['jAperMag3Err'].tolist()
            H_list, H_list_e = table['hAperMag3'].tolist(), table['hAperMag3Err'].tolist()
            K_list, K_list_e = table['kAperMag3'].tolist(), table['kAperMag3Err'].tolist()
            pmra_list, pmra_list_e = table['muRa'].tolist(), table['sigMuRa'].tolist()
            pmdec_list, pmdec_list_e = table['muDec'].tolist(), table['sigMuDec'].tolist()
            enablePrint()

        #Defining a mouse click as an event on the plot
        location = []
        plt.rcParams["figure.figsize"] = [8, 8]
        plt.rcParams["figure.autolayout"] = True
        def mouse_event(event):
            '''Makes a list of the x, y, and axes the mouse click is'''

            location.append(event.ydata)
            location.append(event.xdata)
            location.append(event.inaxes)
        fig_1 = plt.figure()
        cid = fig_1.canvas.mpl_connect('button_press_event', mouse_event)

        #Sets the plot depending on if it was found in UHS or UKIDSS
        if data == 'UHSDR1':

            #Makes the subplot for the plot
            ax = plt.subplot()

            #Gets the cutout for the fits iamge
            wcs = WCS(hdu_j.header)
            position = wcs.world_to_pixel_values(ra, dec)
            size = u.Quantity([(radius * 2.5), (radius * 2.5)], u.pixel)
            cutout_UHS = Cutout2D(data_ukidss_J, position, size)

            #Finds the camera orientation
            cam_type = hdu_j.header['CAMNUM']

            #Obtains the shape of the cutout and sets the circle size for the scatter plot
            shape = max(cutout_UHS.shape)
            circle_size = (radius*3)

            #Converts the ra
            ra_dec_pixel = wcs.world_to_pixel_values(object_ra, object_dec)

            if cam_type == 1: 
                #Makes the dec negative
                minus_dec = [(-x + shape) for x in ra_dec_pixel[1]]
                
                #Plots the correctly orientated image
                scatter = ax.scatter(minus_dec, ra_dec_pixel[0], s = circle_size, edgecolor = '#40E842', facecolor = 'none')
                total_data = np.rot90(cutout_UHS.data, 3)
                plt.xlim(min(cutout_UHS.shape), (min(cutout_UHS.shape) - max(cutout_UHS.shape)))
                plt.ylim(0, shape)

            elif cam_type == 2:
                #Plots the correctly orientated image
                scatter = ax.scatter(ra_dec_pixel[0], ra_dec_pixel[1], s = circle_size, edgecolor = '#40E842', facecolor = 'none')
                total_data = cutout_UHS.data
                plt.xlim(0, shape)
                plt.ylim(min(cutout_UHS.shape), (min(cutout_UHS.shape) - max(cutout_UHS.shape)))

            elif cam_type == 3:
                #Makes the ra negative
                minus_ra = [(-x + shape) for x in ra_dec_pixel[0]]

                #Plots the correctly orientated image
                scatter = ax.scatter(ra_dec_pixel[1], minus_ra, s = circle_size, edgecolor = '#40E842', facecolor = 'none')
                total_data = np.rot90(cutout_UHS.data)
                plt.xlim(min(cutout_UHS.shape), (min(cutout_UHS.shape) - max(cutout_UHS.shape)))
                plt.ylim(0, shape)
                
            elif cam_type == 4: 
                #Plots the correctly orientated image
                scatter = ax.scatter(ra_dec_pixel[0], ra_dec_pixel[1], s = circle_size, edgecolor = '#40E842', facecolor = 'none')
                total_data = cutout_UHS.data
                plt.xlim(min(cutout_UHS.shape), (min(cutout_UHS.shape) - max(cutout_UHS.shape)))
                plt.ylim(0, shape)

            # Normalize the image and plots it
            init_top, init_bot = 95, 45
            norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(cutout_UHS.data, init_bot), vmax = np.nanpercentile(cutout_UHS.data, init_top))
            ax.imshow(total_data, cmap = 'Greys', norm = norm1_w1, origin='lower')

            # Formats the window correctly
            fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
            figure = plt.gcf()
            ax.set_title('Dates: \n' 'J Date: ' + str(date_j) + ' (YYYYMMDD) \n', fontdict=fontdict_1, y=1.05)
            if platform != 'win32':
                figure.set_size_inches(4.75, 7.15)
            elif platform == 'win32': 
                figure.set_size_inches(4.75, 7.15)
        else:
            #Sets the WCS coordinates for the plots
            ax = plt.subplot(projection = wcs_cropped)

            #Plots the objects found in the radius
            circle_size = (radius*3)
            scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

            cam_type = ''

            #Normalize the image and plots it
            init_top, init_bot = 95, 45
            norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, init_bot), vmax = np.nanpercentile(total_data.data, init_top))
            ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

            #Formats the window correctly
            fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
            figure = plt.gcf()
            ax.set_title('Dates: \n'
                    + 'Y Date: ' + str(date_y) + ' (YYYYMMDD)   ' + 'J Date: ' + str(date_j) + ' (YYYYMMDD) \n'
                    + 'H Date: ' + str(date_h) + ' (Y/M/D)   ' + 'K Date: ' + str(date_k) + ' (Y/M/D) \n', fontdict = fontdict_1, y = 1.05)
            if platform != 'win32':
                figure.set_size_inches(4.75, 7.25)
            elif platform == 'win32':
                figure.set_size_inches(4.75, 7.25)
            plt.xlim(max(total_data.shape), 0)
            plt.ylim(0, max(total_data.shape))
            plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)

        #Finishes formatting the problem correctly
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.yaxis.set_tick_params(labelleft=False)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.suptitle('WFCAM Search', fontsize = 35, y = 0.96, fontfamily = 'Times New Roman')
        plt.grid(linewidth = 0)
        figure.canvas.set_window_title('WFCAM Search')
        mng = pyplot.get_current_fig_manager()
        mng.window.resizable(False, False)

        #Make checkbuttons with all of the different image bands
        default = [True]
        real_data = [total_data]
        if data != 'UHSDR1':
            shape = max(cutout_j.shape)
            rax = plt.axes([0.045, 0.4, 0.115, 0.1])
            labels = ['Y', 'J', 'H', 'K']
            real_data = [y_reshape, cutout_j.data, h_reshape, k_reshape]
            default = [True, True, False, False]
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
        axes_button = plt.axes([0.04, 0.775, 0.92, 0.04])
        close = Button(axes_button, 'Object Not Found', color = '#E48671')

         #Update the image depending on what the user chooses
        def update_button(label):
            '''Updates the list of activated images and updates the image the user can see'''

            total_data = 0
            for lab in labels:
                if lab == label:
                    index = labels.index(lab)
                    if default[index] == False:
                        default[index] = True
                    elif default[index] == True: 
                        default[index] = False
            for d in range(len(default)):
                if default == [False, False, False, False]: 
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
            if data == 'UHSDR1':
                ax.imshow(total_data, cmap = 'Greys', norm = norm1_w1)
            else:
                ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

        #Updates the notes added by the user when there is an input
        text_list = [text]
        def submit(expression):
            '''Updates the list of types in the 'Notes' setting'''

            text = expression
            text_list.append(text)

        #Allows the sliders and buttons to be pressed
        if data != 'UHSDR1':
            check.on_clicked(update_button)
        slider_top.on_changed(update_slider_stretch)
        slider_bottom.on_changed(update_slider_stretch)
        text_box.on_text_change(submit)

        #Display image until it is clicked to find the object
        n = -1
        while True:
            press = plt.waitforbuttonpress()
            text_max = len(text_list) - 1
            enablePrint()

            #Checks that it was a mouse click
            if press == False:
                n += 3

                #Finds which axes was clicked
                click_axes = str(location[n])
                if data != 'UHSDR1':
                    click_axes = click_axes.split('WCSAxesSubplot', 2)[0]
                else: 
                    click_axes = click_axes.split('AxesSubplot', 2)[0]

                #Checks if the image was clicked
                if click_axes == '':
                    #Makes a pop-up window with success text
                    plt.clf()
                    plt.close('all')
                    plt.figure(1)
                    plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for WFCAM! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
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
                    figure2.canvas.set_window_title('Successful WFCAM Search')
                    mng2 = pyplot.get_current_fig_manager()
                    mng2.window.resizable(False, False)
                    plt.pause(0.1)
                    plt.clf()
                    plt.close('all')

                    #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
                    if data != 'UHSDR1':
                        coord = wcs_cropped.pixel_to_world_values(location[n-4],location[n-5])
                        distance = []
                        for i in range(len(object_ra)):
                            distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))
                        list_location = distance.index(np.min(distance))
                    else: 
                        if (cam_type == 2) or (cam_type == 4):
                            coord = wcs.pixel_to_world_values(location[n-4],location[n-5])
                            distance = []
                            for i in range(len(object_ra)):
                                distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))
                            list_location = distance.index(np.min(distance))
                        elif cam_type == 3:
                            coord = [location[n - 4], location[n - 5]]
                            distance = []
                            for i in range(len(minus_ra)):
                                distance.append(math.dist(coord, [float(ra_dec_pixel[1][i]), float(minus_ra[i])]))
                            list_location = distance.index(np.min(distance))
                        elif cam_type == 1: 
                            coord = [location[n - 4], location[n - 5]]
                            distance = []
                            for i in range(len(minus_dec)):
                                distance.append(math.dist(coord, [float(minus_dec[i]), float(ra_dec_pixel[0][i])]))
                            list_location = distance.index(np.min(distance))

                    if data != 'UHSDR1':
                        epoch = object_epoch[list_location]
                        ra_wfcam, dec_wfcam = object_ra[list_location], object_dec[list_location]
                        ra_wfcam_e, dec_wfcam_e = object_ra_e[list_location], object_dec_e[list_location]
                        Y_mag, Y_e = Y_list[list_location], Y_list_e[list_location]
                        J_mag, J_e = J_list[list_location], J_list_e[list_location]
                        H_mag, H_e = H_list[list_location], H_list_e[list_location]
                        K_mag, K_e = K_list[list_location], K_list_e[list_location]
                        pmra, pmra_e = pmra_list[list_location], pmra_list_e[list_location]
                        pmdec, pmdec_e = pmdec_list[list_location], pmdec_list_e[list_location]
                        return ra_wfcam, ra_wfcam_e, dec_wfcam, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch, data, text_list[text_max] 
                    else: 
                        ra_wfcam, dec_wfcam = object_ra[list_location], object_dec[list_location]
                        J_mag, J_e = J_list[list_location], J_list_e[list_location]
                        ra_wfcam_e, dec_wfcam_e, Y_mag, Y_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                        return ra_wfcam, ra_wfcam_e, dec_wfcam, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch, data, text_list[text_max] 
                
                #Checks if the "Object Not Found" button was clicked
                elif click_axes == 'Axes(0.04,0.775;0.92x0.04)':
                    ra_wfcam_e, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                    ra_wfcam = ra
                    dec_wfcam = dec
                    
                    #Makes a pop-up window with success text
                    plt.clf()
                    plt.close('all')
                    plt.figure(1)
                    plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for WFCAM! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
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
                    figure2.canvas.set_window_title('Successful WFCAM Search')
                    mng2 = pyplot.get_current_fig_manager()
                    mng2.window.resizable(False, False)
                    plt.pause(0.1)
                    plt.clf()
                    plt.close('all')

                    return ra_wfcam, ra_wfcam_e, dec_wfcam, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch, data, 'Object Not Found was Pressed'
                
                #Updates the circle size when slider is moved
                elif click_axes == 'Axes(0.25,0.055;0.65x0.03)':
                    scatter.remove()
                    if data != 'UHSDR1':
                        scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')
                    else:
                        if (cam_type == 2) or (cam_type == 4):
                            scatter = ax.scatter(ra_dec_pixel[0], ra_dec_pixel[1], s = circle_slider.val, edgecolor = '#40E842', facecolor = 'none')
                        elif cam_type == 3:
                            scatter = ax.scatter(ra_dec_pixel[1], minus_ra, s = circle_slider.val, edgecolor = '#40E842', facecolor = 'none')
                        elif cam_type == 1:
                            scatter = ax.scatter(minus_dec, ra_dec_pixel[0], s = circle_slider.val, edgecolor = '#40E842', facecolor = 'none')

            #Checks if the window was closed
            elif press is None:
                ra_wfcam_e, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                ra_wfcam = ra
                dec_wfcam = dec
                plt.close('all')
                plt.figure().clear()
                return ra_wfcam, ra_wfcam_e, dec_wfcam, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch, data, text_list[text_max] 

    #If the image is not found return null values
    else: 
        ra_wfcam_e, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        ra_wfcam = ra
        dec_wfcam = dec
        text_list = 'Image Not Found'
        return ra_wfcam, ra_wfcam_e, dec_wfcam, dec_wfcam_e, Y_mag, Y_e, J_mag, J_e, H_mag, H_e, K_mag, K_e, pmra, pmra_e, pmdec, pmdec_e, epoch, 'WFCAM', text_list
    
def ukidss_table(ra, dec, radius): 
    '''Find all the objects in the radius defined by the user'''

    blockPrint()

    #Find the table of all the objects found in UKIDSS in the radius choosen by the user
    program_list = Ukidss.list_catalogs()
    for prom in program_list:
        if prom != 'UHS':
            table = Ukidss.query_region(
                    SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'),
                    radius = (radius/2) * u.arcsec, 
                    programme_id = prom)
            if len(table) > 0: 
                return table
        elif prom == 'UHS':
            table = Ukidss.query_region(
                    SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'),
                    radius = ((radius/2) - 1) * u.arcsec, 
                    programme_id = prom, database = 'UHSDR1')
            if len(table) > 0: 
                return table
        