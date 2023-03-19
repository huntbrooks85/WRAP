from catalog_module.importmodule import *

#Makes a function that blocks the printing function
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

#Makes a function that allows the printing function
def enablePrint():
    sys.stdout = sys.__stdout__

#Does an Allwise search
def twomass_image(ra, dec, radius): 
    plt.rcParams['toolbar'] = 'None'
    plt.style.use('Solarize_Light2')
    blockPrint()
    #Finds all the metadata that relates to the ra and dec searched, mostly to find the APIs for the W1, W2, W3, and W4 images
    metadata_2mass_link = 'https://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?POS=' + str(ra) + ',' + str(dec) + '&SIZE=' + str(radius/3600)
    twomass_metadata = requests.get(metadata_2mass_link)
    open('Output/metadata/TWOMASS_metadata.txt', 'wb').write(twomass_metadata.content)

    #With this metadata it finds the API link for the W1 and W2 images
    J_finder, h_finder, k_finder = 'All-Sky Release Survey J-Band Atlas Image', 'All-Sky Release Survey H-Band Atlas Image', 'All-Sky Release Survey K-Band Atlas Image'
    with open('Output/metadata/TWOMASS_metadata.txt', 'r') as fp:
        lines = fp.readlines()
        for line in lines:
            if line.find(J_finder) != -1:
                j_twomass_image_url = ((lines[lines.index(line) + 1]).split('[')[2]).split(']')[0]
            elif line.find(h_finder) != -1:
                h_twomass_image_url = ((lines[lines.index(line) + 1]).split('[')[2]).split(']')[0]
            elif line.find(k_finder) != -1:
                k_twomass_image_url = ((lines[lines.index(line) + 1]).split('[')[2]).split(']')[0]

    #Download the W1 and W2 images
    file_allwise_j = download_file(j_twomass_image_url, cache=True)
    data_allwise_j = fits.getdata(file_allwise_j)
    file_allwise_h = download_file(h_twomass_image_url, cache=True)
    data_allwise_h = fits.getdata(file_allwise_h)
    file_allwise_k = download_file(k_twomass_image_url, cache=True)
    data_allwise_k = fits.getdata(file_allwise_k)

    #Find the location of all the object found in AllWISE in the radius choosen by the user 
    location_data = twomass_table(ra, dec, radius)
    object_ra = location_data['ra'].tolist()
    object_dec = location_data['dec'].tolist()
    j_list = location_data['j_m'].tolist()
    h_list = location_data['h_m'].tolist()
    ks_list = location_data['k_m'].tolist()

    #Make a cutout from the coadd image for the ra and dec put in
    hdu_j = fits.open(file_allwise_j)[0]
    hdu_h = fits.open(file_allwise_h)[0]
    hdu_k = fits.open(file_allwise_k)[0]
    wcs1_w1 = WCS(hdu_j.header)

    position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
    size = u.Quantity([radius, radius], u.arcsec)
    cutout_j = Cutout2D(data_allwise_j, position, size, fill_value = np.nan, wcs = wcs1_w1.celestial)
    cutout_h = Cutout2D(data_allwise_h, position, size, fill_value = np.nan, wcs = wcs1_w1.celestial)
    cutout_k = Cutout2D(data_allwise_k, position, size, fill_value = np.nan, wcs = wcs1_w1.celestial)
    wcs_cropped = cutout_j.wcs

    enablePrint()
    date_j = hdu_j.header[64]
    time_j = hdu_j.header[65]
    date_h = hdu_h.header[64]
    time_h = hdu_h.header[65]
    date_k = hdu_k.header[64]
    time_k = hdu_k.header[65]
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
    total_data = cutout_j.data + cutout_k.data
    ax = plt.subplot(projection = wcs_cropped)

    #Plots the objects found in the radius
    circle_size = (radius*3)
    scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

    #Normalize the image and plots it
    init_top = 95
    init_bot = 45
    norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, init_bot), vmax = np.nanpercentile(total_data.data, init_top))
    ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

    #Makes the figure look pretty
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
    fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
    ax.set_title(     'Dates: \n'
                    + 'J Date: ' + str(date_j) + ' (Y/M/D)  ' + '  H Date: ' + str(date_h) + ' (Y/M/D)\n'
                    + 'K Date: ' + str(date_k) + ' (Y/M/D)  \n'
                    + 'Times: \n'
                    + 'J Time: ' + str(time_j) + ' (UT)  ' + '  H Time: ' + str(time_h) + ' (UT)\n'
                    + 'K Time: ' + str(time_k) + ' (UT) \n', fontdict = fontdict_1, y = 1)
    plt.grid(linewidth = 0)
    figure = plt.gcf()
    figure.set_size_inches(4.75, 6.95)
    figure.canvas.set_window_title('2MASS Search')

    #Make checkbuttons with all of the different image bands
    rax = plt.axes([0.045, 0.4, 0.105, 0.12])
    labels = ['J', 'H', 'K']
    real_data = [cutout_j.data, cutout_h.data, cutout_k.data]
    default = [True, False, True]
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
                coord = wcs_cropped.pixel_to_world_values(location[n-5],location[n-4])
                distance = []
                for i in range(len(object_ra)):
                    distance.append(math.dist(coord, [object_ra[i], object_dec[i]]))
                list_location = distance.index(np.min(distance))
                j = j_list[list_location]
                h = h_list[list_location]
                ks = ks_list[list_location]
                return j, h, ks, text_list[text_max]
            
            #Checks if the Object not Found button was clicked
            elif click_axes == 'Axes(0.04,0.012;0.92x0.04)':
                j, h, ks= np.nan, np.nan, np.nan
                plt.close('all')
                plt.figure().clear()
                return j, h, ks, text_list[text_max]
            
            elif click_axes == 'Axes(0.25,0.095;0.65x0.03)':
                scatter.remove()
                scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')
            
        #Checks if the window was closed
        elif press is None:
            j, h, ks= np.nan, np.nan, np.nan
            plt.close('all')
            plt.figure().clear()
            return j, h, ks, text_list[text_max]

#Find all the objects in the radius defined by the user
def twomass_table(ra, dec, radius): 
    blockPrint()
    if dec < 0:
        ra_hours, dec_hours = int(ra*(24/360)), int(dec)
        dec_hours = '+' + str(dec_hours)
        ra_min, dec_min = int((ra*(24/360))%1 * 60), int((int(dec) - dec) * 60)
        ra_sec, dec_sec = (((ra*(24/360))%1 * 60)%1 * 60), (((int(dec) - dec) * 60)%1 * 60)
        ra_tot, dec_tot = (str(ra_hours) + 'h+' + str(ra_min) + 'm+' + str(ra_sec) + 's'), (str(dec_hours) + 'd+' + str(dec_min) + 'm+' + str(dec_sec) + 's')
    elif dec > 0: 
        ra_hours, dec_hours = int(ra*(24/360)), int(dec)
        dec_hours = '+' + str(dec_hours)
        ra_min, dec_min = int((ra*(24/360))%1 * 60), int((dec)%1 * 60)
        ra_sec, dec_sec = (((ra*(24/360))%1 * 60)%1 * 60), (((dec)%1 * 60)%1 * 60)
        ra_tot, dec_tot = (str(ra_hours) + 'h+' + str(ra_min) + 'm+' + str(ra_sec) + 's'), (str(dec_hours) + 'd+' + str(dec_min) + 'm+' + str(dec_sec) + 's')

    #Obtains the url for the data table relating to the location chosen 
    twomass_table_url = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?spatial=box&catalog=fp_psc&objstr=' + str(ra_tot) + str(dec_tot) + '&size=' + str(radius) + '&outfmt=1'
    twomass_table = ascii.read(twomass_table_url, format = 'ipac')

    return twomass_table