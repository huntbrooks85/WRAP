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
    hdu = fits.open(file_allwise_j)[0]
    wcs1_w1 = WCS(hdu.header)

    position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
    size = u.Quantity([radius, radius], u.arcsec)
    cutout_j = Cutout2D(data_allwise_j, position, size, fill_value = np.nan, wcs = wcs1_w1.celestial)
    cutout_h = Cutout2D(data_allwise_h, position, size, fill_value = np.nan, wcs = wcs1_w1.celestial)
    cutout_k = Cutout2D(data_allwise_k, position, size, fill_value = np.nan, wcs = wcs1_w1.celestial)
    wcs_cropped = cutout_j.wcs

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
    j_h_k = cutout_j.data + cutout_h.data  + cutout_k.data
    ax = plt.subplot(projection = wcs_cropped)

    #Plots the objects found in the radius
    plt.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = (radius*3), edgecolor='#40E842', facecolor='none')

    #Add a mouse hovering ability
    mplcursors.cursor(hover = 2, highlight = True)

    #Normalize the image and plots it
    norm1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(j_h_k.data, 45), vmax = np.nanpercentile(j_h_k.data, 98))
    plt.imshow(j_h_k.data, cmap = 'Greys', norm = norm1)

    #Makes the figure look pretty
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
    fontdict_1 = {'family':'Times New Roman','color':'k','size':20}
    plt.title('2MASS Search', fontdict = fontdict_1)
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

    #Find the closest point to the location clicked to obtain W1, W2, W3, and W4
    if len(location) > 0:
        if str(location[2]) == 'Axes(0.04,0.012;0.92x0.04)': 
            j, h, ks = np.nan, np.nan, np.nan
            return j, h, ks
        else: 
            coord = wcs_cropped.pixel_to_world_values(location[0],location[1])
            distance = []
            for i in range(len(object_ra)):
                distance.append(math.dist(coord, [object_ra[i], object_dec[i]]))
            list_location = distance.index(np.min(distance))
            j = j_list[list_location] 
            h = h_list[list_location]
            ks = ks_list[list_location]
            return j, h, ks
    elif len(location) == 0:
        j, h, ks = np.nan, np.nan, np.nan
        return j, h, ks

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