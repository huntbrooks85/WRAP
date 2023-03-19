from catalog_module.importmodule import *

#Makes a function that blocks the printing function
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

#Makes a function that allows the printing function
def enablePrint():
    sys.stdout = sys.__stdout__

#Does an PanSTARRS search
def ps_image(ra, dec, radius):
    plt.rcParams['toolbar'] = 'None'
    plt.style.use('Solarize_Light2')
    blockPrint()
    #Makes the dec usable for the image url API
    if dec > 0: 
        dec = '+' + str(dec)

    #Find the panstarr metadata for the image API
    ps_image_url = 'http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=' + str(ra) + str(dec)  + '&filter=color&filter=g&filter=r&filter=i&filter=z&filter=y&filetypes=stack&auxiliary=data&size=' + str(radius * 4) + '&output_size=0&verbose=0&autoscale=99.500000&catlist='
    allwise_metadata = requests.get(ps_image_url)
    open('Output/metadata/ps_metadata.txt', 'wb').write(allwise_metadata.content)

    #With the metadata from panstarrs find the image API in r band
    r_finder= 'amp;'
    r_finder_list = []
    with open('Output/metadata/ps_metadata.txt', 'r') as fp:
        lines = fp.readlines()
        for line in lines:
            if line.find(r_finder) != -1:
                r_finder_list.append((lines[lines.index(line)]))
    if len(r_finder_list) > 0: 
        r_link = ('http:' + (r_finder_list[2]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)
        i_link = ('http:' + (r_finder_list[3]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)
        z_link = ('http:' + (r_finder_list[4]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)

        #Downloads the panstarrs r, i, and z image
        file_ps_r = download_file(r_link, cache=True)
        data_ps_r = fits.getdata(file_ps_r)
        file_ps_i = download_file(i_link, cache=True)
        data_ps_i = fits.getdata(file_ps_i)
        file_ps_z = download_file(z_link, cache=True)
        data_ps_z = fits.getdata(file_ps_z)

        #Obtains the data table for panstarrs and then gets all of the relatvent data from the table
        panstarr_data = ps_table(ra, dec, radius)
        ra_list = panstarr_data['raStack'].tolist()
        dec_list = panstarr_data['decStack'].tolist()
        g_list = panstarr_data['gMeanApMag'].tolist()
        r_list = panstarr_data['rMeanApMag'].tolist()
        i_list = panstarr_data['iMeanApMag'].tolist()
        z_list = panstarr_data['zMeanApMag'].tolist()
        y_list = panstarr_data['yMeanApMag'].tolist()

        #Removes the null RA and DEC values that PanSTARRS return
        for elem in list(ra_list):
            if elem == -999.0:
                ra_location = ra_list.index(elem)
                dec_list.pop(ra_location)
                g_list.pop(ra_location)
                r_list.pop(ra_location)
                i_list.pop(ra_location)
                z_list.pop(ra_location)
                y_list.pop(ra_location)
                ra_list.remove(elem)

        #Make a cutout from the coadd image for the ra and dec put in
        hdu = fits.open(file_ps_r)[0]
        wcs1 = WCS(hdu.header)

        if type(dec) == str: 
            dec = float(dec.replace('+', '', 1))

        position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
        size = u.Quantity([radius, radius], u.arcsec)
        cutout_r = Cutout2D(data_ps_r, position, size, fill_value = np.nan, wcs = wcs1.celestial)
        cutout_i = Cutout2D(data_ps_i, position, size, fill_value = np.nan, wcs = wcs1.celestial)
        cutout_z = Cutout2D(data_ps_z, position, size, fill_value = np.nan, wcs = wcs1.celestial)
        wcs_cropped = cutout_r.wcs

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
        r_i_z = cutout_r.data + cutout_i.data + cutout_z.data
        ax = plt.subplot(projection = wcs_cropped)

        #Plots the objects found in the radius
        plt.scatter(ra_list, dec_list, transform=ax.get_transform('fk5'), s = (radius*3), edgecolor='#40E842', facecolor='none') 

        #Add a mouse hovering ability
        mplcursors.cursor(hover = 2, highlight = True)

        #Normalize the image and plots it
        norm1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(r_i_z.data, 45), vmax = np.nanpercentile(r_i_z.data, 99))
        plt.imshow(r_i_z.data, cmap = 'Greys', norm = norm1)
        
        #Makes the figure look pretty
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        fontdict_1 = {'family':'Times New Roman','color':'k','size':20}
        plt.title('PanSTARRS Search', fontdict = fontdict_1)
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

        #Find the closest point to the location clicked to obtain g, r, i, z, and y
        if len(location) > 0:
            if str(location[2]) == 'Axes(0.04,0.012;0.92x0.04)': 
                g, r, i, z, y = np.nan, np.nan, np.nan, np.nan, np.nan
                return g, r, i, z, y
            else: 
                coord = wcs_cropped.pixel_to_world_values(location[0],location[1])
                distance = []
                for p in range(len(ra_list)):
                    distance.append(math.dist(coord, [ra_list[p], dec_list[p]]))
                list_location = distance.index(np.min(distance))
                g = g_list[list_location] 
                r = r_list[list_location] 
                i = i_list[list_location] 
                z = z_list[list_location] 
                y = y_list[list_location] 
                return g, r, i, z, y

        elif len(location) == 0: 
            g, r, i, z, y = np.nan, np.nan, np.nan, np.nan, np.nan
            return g, r, i, z, y
    
    elif len(r_finder_list) == 0: 
        g, r, i, z, y = np.nan, np.nan, np.nan, np.nan, np.nan
        return g, r, i, z, y

#Find all the objects in the radius defined by the user
def ps_table(ra, dec, radius): 
    blockPrint()
    #Finds the table for the Panstarrs data around the ra and dec given by the user
    url_panstarr = 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr1/mean?ra=' + str(ra) + '&dec=' + str(dec) + '&radius=' + str(radius/7200) + '&nDetections.gte=1&pagesize=5001&format=csv'
    test_url_panstarr = requests.get(url_panstarr)

    #Checks to see if the table exists and then obtains all of the relavent data to the search
    if test_url_panstarr.status_code == 200:
        panstarr_data = pd.read_csv(url_panstarr)
    return panstarr_data