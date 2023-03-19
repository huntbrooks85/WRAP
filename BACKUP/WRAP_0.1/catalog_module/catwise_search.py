from catalog_module.importmodule import *

#Makes a function that blocks the printing function
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

#Makes a function that allows the printing function
def enablePrint():
    sys.stdout = sys.__stdout__

def catwise_image(ra, dec, radius): 
    plt.rcParams['toolbar'] = 'None'
    plt.style.use('Solarize_Light2')
    blockPrint()
    #Finds all the metadata that relates to the ra and dec searched, mostly to find the APIs for the W1, W2, W3, and W4 images
    metadata_allwise_link = 'http://irsa.ipac.caltech.edu/ibe/sia/wise/allwise/p3am_cdd?POS=' + str(ra) + ',' + str(dec) + '&SIZE=' + str(radius/3600)
    allwise_metadata = requests.get(metadata_allwise_link)
    open('Output/metadata/catwise_metadata.txt', 'wb').write(allwise_metadata.content)

    #With this metadata it finds the API link for the W1 and W2 images
    w1_finder, w2_finder = 'W1 Coadd', 'W2 Coadd'
    with open('Output/metadata/catwise_metadata.txt', 'r') as fp:
        lines = fp.readlines()
        for line in lines:
            if line.find(w1_finder) != -1:
                w1_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]
            elif line.find(w2_finder) != -1:
                w2_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]

    #Download the W1 and W2 images
    file_allwise_w1 = download_file(w1_allwise_image_url, cache=True)
    data_allwise_w1 = fits.getdata(file_allwise_w1)
    file_allwise_w2 = download_file(w2_allwise_image_url, cache=True)
    data_allwise_w2 = fits.getdata(file_allwise_w2)

    #Find the location of all the object found in CatWISE in the radius choosen by the user 
    location_data = catwise_table(ra, dec, radius)
    object_ra = location_data['ra'].tolist()
    object_dec = location_data['dec'].tolist()
    w1_list = location_data['w1mpro'].tolist()
    w2_list = location_data['w2mpro'].tolist()
    
    #Make a cutout from the coadd image for the RA and DEC put in
    hdu_w1 = fits.open(file_allwise_w1)[0]
    wcs1_w1 = WCS(hdu_w1.header)

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
    norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(w1_w2.data, 45), vmax = np.nanpercentile(w1_w2.data, 95))
    plt.imshow(w1_w2.data, cmap = 'Greys', norm = norm1_w1)
    
    #Makes the figure look pretty
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
    fontdict_1 = {'family':'Times New Roman','color':'k','size':20}
    plt.title('CatWISE 2020 Search', fontdict = fontdict_1)
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

    #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
    if len(location) > 0:
        if str(location[2]) == 'Axes(0.04,0.012;0.92x0.04)': 
            w1, w2 = np.nan, np.nan
            return w1, w2
        else: 
            coord = wcs_cropped_w1.pixel_to_world_values(location[0],location[1])
            distance = []
            for i in range(len(object_ra)):
                distance.append(math.dist(coord, [object_ra[i], object_dec[i]]))
            list_location = distance.index(np.min(distance))
            w1 = w1_list[list_location] 
            w2 = w2_list[list_location]
            return w1, w2

    elif len(location) == 0: 
        w1, w2= np.nan, np.nan
        return w1, w2

def catwise_table(ra, dec, radius): 
    blockPrint()
    #Converts degree to HH:MM:SS and DD:MM:SS
    if dec < 0:
        ra_hours, dec_hours = int(ra*(24/360)), int(dec)
        ra_min, dec_min = int((ra*(24/360))%1 * 60), int((int(dec) - dec) * 60)
        ra_sec, dec_sec = (((ra*(24/360))%1 * 60)%1 * 60), (((int(dec) - dec) * 60)%1 * 60)
        ra_tot, dec_tot = (str(ra_hours) + 'h+' + str(ra_min) + 'm+' + str(ra_sec) + 's'), (str(dec_hours) + 'd+' + str(dec_min) + 'm+' + str(dec_sec) + 's')
    elif dec > 0: 
        ra_hours, dec_hours = int(ra*(24/360)), int(dec)
        ra_min, dec_min = int((ra*(24/360))%1 * 60), int((dec)%1 * 60)
        ra_sec, dec_sec = (((ra*(24/360))%1 * 60)%1 * 60), (((dec)%1 * 60)%1 * 60)
        ra_tot, dec_tot = (str(ra_hours) + 'h+' + str(ra_min) + 'm+' + str(ra_sec) + 's'), (str(dec_hours) + 'd+' + str(dec_min) + 'm+' + str(dec_sec) + 's')

    #Obtains the url for the data table relating to the location chosen 
    location_url = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?catalog=catwise_2020&spatial=cone&radius=' + str(radius/2) + '&radunits=arcsec&objstr=' + str(ra_tot) + '+' + str(dec_tot) + '&size=' + str(radius/2) + '&outfmt=1&selcols=ra,dec,w1mpro,w2mpro'
    location_data = ascii.read(location_url, format = 'ipac')
    return location_data