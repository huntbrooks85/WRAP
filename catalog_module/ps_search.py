#Import all of the packages
from catalog_module.importmodule import *

def blockPrint():
  '''Makes a function that blocks the printing function'''
  sys.stdout = open(os.devnull, 'w')

def enablePrint():
  '''Makes a function that allows the printing function'''
  sys.stdout = sys.__stdout__

def ps_image(ra, dec, radius):
  ''' First, it gets the images from the PanSTARRS API from MAST and downloads the images. 
  Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
  Third, calls the table function to get all the objects from the PanSTARRS source catalog. 
  Fourth, makes the window for the user to click the object with all settings. 
  Finally, finds the closest object to the click and records the data. '''
  
  #Makes outline for the window of the plot
  plt.rcParams['toolbar'] = 'None'
  plt.style.use('Solarize_Light2')
  blockPrint()

  #Makes the dec usable for the image url API
  if dec > 0: 
    new_dec = '+' + str(dec)
  elif dec <= 0: 
    new_dec = str(dec)

  #Find the panstarr metadata for the image API
  ps_image_url = 'http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=' + str(ra) + str(new_dec)  + '&filter=color&filter=g&filter=r&filter=i&filter=z&filter=y&filetypes=stack&auxiliary=data&size=' + str(radius * 4) + '&output_size=0&verbose=0&autoscale=99.500000&catlist='
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

  #Checks if the images are found
  if len(r_finder_list) > 0: 

    #Obtains the data table for panstarrs and then gets all of the relatvent data from the table
    panstarr_data = ps_table(ra, dec, radius)
    object_mjd = panstarr_data['epochMean'].tolist()
    ra_list, dec_list = panstarr_data['raStack'].tolist(), panstarr_data['decStack'].tolist()
    ra_list_e, dec_list_e = panstarr_data['raStackErr'].tolist(), panstarr_data['decStackErr'].tolist()
    g_list, g_list_e = panstarr_data['gApMag'].tolist(), panstarr_data['gApMagErr'].tolist()
    r_list, r_list_e = panstarr_data['rApMag'].tolist(), panstarr_data['rApMagErr'].tolist()
    i_list, i_list_e = panstarr_data['iApMag'].tolist(), panstarr_data['iApMagErr'].tolist()
    z_list, z_list_e = panstarr_data['zApMag'].tolist(), panstarr_data['zApMagErr'].tolist()
    y_list, y_list_e = panstarr_data['yApMag'].tolist(), panstarr_data['yApMagErr'].tolist()

    #Gets all of the image urls
    r_link = ('http:' + (r_finder_list[2]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)
    i_link = ('http:' + (r_finder_list[3]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)
    z_link = ('http:' + (r_finder_list[4]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)
    y_link = ('http:' + (r_finder_list[5]).split('href="', 3)[3].split('"', 1)[0]).replace('amp;', '', 6)

    #Downloads the panstarrs r, i, and z image
    file_ps_r, file_ps_i, file_ps_z, file_ps_y = download_file(r_link, cache=True), download_file(i_link, cache=True), download_file(z_link, cache=True), download_file(y_link, cache=True)
    data_ps_r, data_ps_i, data_ps_z, data_ps_y = fits.getdata(file_ps_r), fits.getdata(file_ps_i), fits.getdata(file_ps_z), fits.getdata(file_ps_y)

    #Removes the null RA and DEC values that PanSTARRS return
    for elem in list(ra_list):
      if elem == -999.0 or elem is None:
        ra_location = ra_list.index(elem)
        ra_list.remove(elem)
        dec_list.pop(ra_location)
        g_list.pop(ra_location)
        g_list_e.pop(ra_location)
        r_list.pop(ra_location)
        r_list_e.pop(ra_location)
        i_list.pop(ra_location)
        i_list_e.pop(ra_location)
        z_list.pop(ra_location)
        z_list_e.pop(ra_location)
        y_list.pop(ra_location)
        y_list_e.pop(ra_location)

    #Gets the headers from the images
    hdu_r, hdu_i, hdu_z, hdu_y = fits.open(file_ps_r)[0], fits.open(file_ps_i)[0], fits.open(file_ps_z)[0], fits.open(file_ps_y)[0]
    wcs1 = WCS(hdu_r.header)

    #Make a cutout from the coadd image for the ra and dec put in
    if type(dec) == str: 
      dec = float(dec.replace('+', '', 1))
    position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
    size = u.Quantity([radius, radius], u.arcsec)
    cutout_r = Cutout2D(data_ps_r, position, size, fill_value = np.nan, wcs = wcs1.celestial)
    cutout_i = Cutout2D(data_ps_i, position, size, fill_value = np.nan, wcs = wcs1.celestial)
    cutout_z = Cutout2D(data_ps_z, position, size, fill_value = np.nan, wcs = wcs1.celestial)
    cutout_y = Cutout2D(data_ps_y, position, size, fill_value = np.nan, wcs = wcs1.celestial)
    wcs_cropped = cutout_r.wcs
    enablePrint()

    #Gets the dates that each image was taken
    date_r, date_i, date_z, date_y = hdu_r.header[8].split('T', 2)[0], hdu_i.header[8].split('T', 2)[0], hdu_z.header[8].split('T', 2)[0], hdu_y.header[8].split('T', 2)[0]

    #Defining a mouse click as an event on the plot
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
    total_data = cutout_r.data + cutout_i.data
    ax = plt.subplot(projection = wcs_cropped)
    
    #Plots the objects found in the radius
    circle_size = (radius*3)
    scatter = ax.scatter(ra_list, dec_list, transform=ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

    #Normalize the image and plots it
    init_top = 95
    init_bot = 45
    norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, init_bot), vmax = np.nanpercentile(total_data.data, init_top))
    ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

    #Makes the figure look pretty
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
    fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
    plt.suptitle('PanSTARRS Search', fontsize = 35, y = 0.96, fontfamily = 'Times New Roman')
    ax.set_title('Dates: \n'
              + 'r Date: ' + str(date_r) + ' (Y-M-D)  ' + '  i Date: ' + str(date_i) + ' (Y-M-D)   \n'
              + 'z Date: ' + str(date_z) + ' (Y-M-D)  ' + '  y Date: ' + str(date_y) + ' (Y-M-D)   \n', fontdict = fontdict_1, y = 1.05)
    plt.grid(linewidth = 0)
    figure = plt.gcf()
    if platform != 'win32':
      figure.set_size_inches(4.75, 6.95)
    elif platform == 'win32':
      figure.set_size_inches(4.75, 7.25)
    figure.canvas.set_window_title('PanSTARRS Search')

    #Make checkbuttons with all of the different image bands
    rax = plt.axes([0.045, 0.4, 0.115, 0.1])
    labels = ['r', 'i', 'z', 'y']
    real_data = [cutout_r.data, cutout_i.data, cutout_z.data, cutout_y.data]
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
    text_box = TextBox(axbox, 'Notes:', initial = text, textalignment = "center")

    #Make a button that can be clicked if no object is found
    axes_button = plt.axes([0.04, 0.775, 0.92, 0.04])
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
          coord = wcs_cropped.pixel_to_world_values(location[n-4],location[n-5])
          distance = []
          for i in range(len(ra_list)):
            distance.append(math.dist(coord, [ra_list[i], dec_list[i]]))

          list_location = distance.index(np.min(distance))
          mjd = object_mjd[list_location]
          ps_ra, ps_dec = ra_list[list_location], dec_list[list_location]
          ps_ra_e, ps_dec_e = ra_list_e[list_location], dec_list_e[list_location]
          g, g_e = g_list[list_location], g_list_e[list_location]
          r, r_e = r_list[list_location], r_list_e[list_location]
          i, i_e = i_list[list_location], i_list_e[list_location]
          z, z_e = z_list[list_location], z_list_e[list_location]
          y, y_e = y_list[list_location], y_list_e[list_location]
          return ps_ra, ps_ra_e, ps_dec, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd, text_list[text_max]
        
        #Checks if the Object not Found button was clicked
        elif click_axes == 'Axes(0.04,0.775;0.92x0.04)':
          ps_ra_e, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
          ps_ra, ps_dec = ra, dec
          plt.close('all')
          plt.figure().clear()
          return ps_ra, ps_ra_e, ps_dec, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd, 'Object Not Found was Pressed'
        
        #Updates the circle size when slider is moved
        elif click_axes == 'Axes(0.25,0.055;0.65x0.03)':
          scatter.remove()
          scatter = ax.scatter(ra_list, dec_list, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')

      #Checks if the window was closed
      elif press is None:
        ps_ra_e, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        ps_ra, ps_dec = ra, dec
        plt.close('all')
        plt.figure().clear()
        return ps_ra, ps_ra_e, ps_dec, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd, text_list[text_max]
    
  #If the images were not found returns null values
  elif len(r_finder_list) == 0: 
    ps_ra_e, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    ps_ra, ps_dec = ra, dec
    return ps_ra, ps_ra_e, ps_dec, ps_dec_e, g, g_e, r, r_e, i, i_e, z, z_e, y, y_e, mjd, 'Image Not Found'

#Find all the objects in the radius defined by the user
def ps_table(ra, dec, radius): 
  '''Find all the objects in the radius defined by the user'''
  
  blockPrint()

  # #Finds the table for the Panstarrs data around the ra and dec given by the user
  catalog_data = Catalogs.query_region(str(ra) + ' ' + str(dec), radius = ((radius/7200) - 0.000277778), catalog = "Panstarrs", table = "stack")
  return catalog_data