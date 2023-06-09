#Import all of the packages
from catalog_module.importmodule import *

def blockPrint():
  '''Makes a function that blocks the printing function'''
  sys.stdout = open(os.devnull, 'w')

def enablePrint():
  '''Makes a function that allows the printing function'''
  sys.stdout = sys.__stdout__

def vsa_image(ra, dec, radius): 
  ''' First, it gets the images from the VSA API from the VSA Archive and downloads the images. 
  Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
  Third, calls the table function to get all the objects from the VSA source catalog. 
  Fourth, makes the window for the user to click the object with all settings. 
  Finally, finds the closest object to the click and records the data. '''
  
  #Makes outline for the window of the plot
  plt.rcParams['toolbar'] = 'None'
  matplotlib.use("TkAgg")
  plt.style.use('Solarize_Light2')
  blockPrint()

  #Obtains all of the urls in J, H, and K from VSA
  database_list = ['VHSDR6', 'VVVDR5', 'VMCDR6', 'VIKINGDR5', 'VIDEODR6', 'ULTRAVISTADR4']
  for img in database_list:
    url_J, url_H, url_Ks = [Vsa.get_image_list(
      SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'),
                image_width = radius * u.arcsec,
                waveband = band,
                database = img)
      for band in ['J', 'H', 'Ks']]
    if len(url_J) == 0:
      pass
    else: 
      break
  
  #Checking to see if the images exist
  if len(url_J) > 0 and len(url_H) > 0 and len(url_Ks) > 0: 

    #Downloading the fits images
    file_vsa_J, file_vsa_H, file_vsa_Ks = download_file(url_J[0], cache=True), download_file(url_H[0], cache=True), download_file(url_Ks[0], cache=True)
    data_vsa_J, data_vsa_H, data_vsa_Ks = fits.getdata(file_vsa_J), fits.getdata(file_vsa_H), fits.getdata(file_vsa_Ks)

    table = vsa_table(ra, dec, radius)
    
    #Obtains the photometry and astrometry from the catalog
    object_mjd_y, object_mjd_j = table['yMjd'].tolist(), table['jMjd'].tolist()
    object_mjd_h, object_mjd_ks = table['hMjd'].tolist(), table['ksMjd'].tolist()
    object_ra, object_dec = table['ra'].tolist(), table['dec'].tolist()
    Y_list, Y_list_e = table['yAperMag3'].tolist(), table['yAperMag3Err'].tolist()
    J_list, J_list_e = table['jAperMag3'].tolist(), table['jAperMag3Err'].tolist()
    H_list, H_list_e = table['hAperMag3'].tolist(), table['hAperMag3Err'].tolist()
    Ks_list, Ks_list_e = table['ksAperMag3'].tolist(), table['ksAperMag3Err'].tolist()

    #Reads in the header from the image
    hdu_j, hdu_h, hdu_k = fits.open(file_vsa_J)[1], fits.open(file_vsa_H)[1], fits.open(file_vsa_Ks)[1]
    wcs1_j = WCS(hdu_j.header)

    #Make a cutout from the coadd image for the RA and DEC put in
    position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5')
    size = u.Quantity([radius, radius], u.arcsec)
    cutout_j = Cutout2D(data_vsa_J, position, size, wcs = wcs1_j.celestial)
    cutout_h = Cutout2D(data_vsa_H, position, size, wcs = wcs1_j.celestial)
    cutout_ks = Cutout2D(data_vsa_Ks, position, size, wcs = wcs1_j.celestial)
    wcs_cropped = cutout_j.wcs
    enablePrint()

    #Finds the dates the images were taken
    date_j, date_h, date_k = hdu_j.header['DATE-OBS'].split('T', 2)[0], hdu_h.header['DATE-OBS'].split('T', 2)[0], hdu_k.header['DATE-OBS'].split('T', 2)[0]

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
    ks_shape = cutout_ks.shape
    j_reshape, h_reshape = cv2.resize(cutout_j.data, (ks_shape[1], ks_shape[0]), interpolation = cv2.INTER_NEAREST), cv2.resize(cutout_h.data, (ks_shape[1], ks_shape[0]), interpolation = cv2.INTER_NEAREST)
    total_data = j_reshape + cutout_ks.data
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
    plt.suptitle('VISTA Search', fontsize = 35, y = 0.96, fontfamily = 'Times New Roman')
    ax.set_title('Dates: \n'
               + 'J Date: ' + str(date_j) + ' (Y-M-D)  ' + '  H Date: ' + str(date_h) + ' (Y-M-D)\n'
               + 'K Date: ' + str(date_k) + ' (Y-M-D) \n', fontdict = fontdict_1, y = 1.05)
    plt.grid(linewidth = 0)
    shape = min(cutout_ks.shape)
    plt.xlim(0, shape)
    plt.ylim(shape, 0)
    figure = plt.gcf()
    if platform != 'win32': 
      figure.set_size_inches(4.75, 7.25)
    elif platform == 'win32': 
      figure.set_size_inches(4.75, 7.25)
    figure.canvas.set_window_title('VISTA Search')
    mng = pyplot.get_current_fig_manager()
    mng.window.resizable(False, False)

    #Make checkbuttons with all of the different image bands
    rax = plt.axes([0.045, 0.4, 0.105, 0.12])
    labels = ['J', 'H', 'Ks']
    real_data = [j_reshape, h_reshape, cutout_ks.data]
    default = [True, False, True]
    check = CheckButtons(rax, labels, default)

    #Adds a slider for the scaling of the image
    freq_top = plt.axes([0.25, 0.12, 0.65, 0.03])
    slider_top = Slider(ax = freq_top, label = 'Top Stetch:', valmin = 50, valmax = 100, valinit = init_top, color = '#E48671')
    freq_bottom = plt.axes([0.25, 0.087, 0.65, 0.03])
    slider_bottom = Slider(ax = freq_bottom, label = 'Bottom Stetch:', valmin = 0, valmax = 50, valinit = init_bot, color = '#E48671')

    #Adds a slider for the circle size
    circle_slid_location = plt.axes([0.25, 0.055, 0.65, 0.03])
    circle_slider = Slider(ax = circle_slid_location, label = 'Circle Size:', valmin = 5, valmax = (circle_size + 1*radius), valinit = circle_size, color = '#E48671')

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
        if default == [False, False, False]:
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

          #Makes a pop-up window with success text
          plt.clf()
          plt.close('all')
          shape_x, shape_y = total_data.shape[0], total_data.shape[1]
          plt.figure(1)
          plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for VISTA! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
          plt.xlim(0, 1)
          plt.ylim(0, 1)
          plt.grid(linewidth = 0)
          ax = plt.gca()
          ax.xaxis.set_tick_params(labelbottom=False)
          ax.yaxis.set_tick_params(labelleft=False)
          ax.set_xticks([])
          ax.set_yticks([])
          figure2 = plt.gcf()
          figure2.set_size_inches(4.6, 1)
          figure2.canvas.set_window_title('Successful VISTA Search')
          mng2 = pyplot.get_current_fig_manager()
          mng2.window.resizable(False, False)
          plt.pause(0.1)
          plt.clf()
          plt.close('all')

          #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
          coord = wcs_cropped.pixel_to_world_values(location[n-4],location[n-5])
          distance = []
          for i in range(len(object_ra)):
            distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))

          list_location = distance.index(np.min(distance))
          ymjd, jmjd = object_mjd_y[list_location], object_mjd_j[list_location]
          hmjd, ksmjd = object_mjd_h[list_location], object_mjd_ks[list_location]
          ra_vsa, dec_vsa = object_ra[list_location], object_dec[list_location]
          y, y_e = Y_list[list_location], Y_list_e[list_location]
          j, j_e = J_list[list_location], J_list_e[list_location]
          h, h_e = H_list[list_location], H_list_e[list_location]
          ks, ks_e = Ks_list[list_location], Ks_list_e[list_location]
          return ra_vsa, dec_vsa, y, y_e, j, j_e, h, h_e, ks, ks_e, ymjd, jmjd, hmjd, ksmjd, img, text_list[text_max]
        
        #Checks if the Object not Found button was clicked
        elif click_axes == 'Axes(0.04,0.775;0.92x0.04)':
          ymjd, jmjd, hmjd, ksmjd, y, y_e, j, j_e, h, h_e, ks, ks_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
          ra_vsa, dec_vsa = ra, dec
          shape_x, shape_y = total_data.shape[0], total_data.shape[1]

          #Makes a pop-up window with success text
          plt.clf()
          plt.close('all')
          shape_x, shape_y = total_data.shape[0], total_data.shape[1]
          plt.figure(1)
          plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for VISTA! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
          plt.xlim(0, 1)
          plt.ylim(0, 1)
          plt.grid(linewidth = 0)
          ax = plt.gca()
          ax.xaxis.set_tick_params(labelbottom=False)
          ax.yaxis.set_tick_params(labelleft=False)
          ax.set_xticks([])
          ax.set_yticks([])
          figure2 = plt.gcf()
          figure2.set_size_inches(4.6, 1)
          figure2.canvas.set_window_title('Successful VISTA Search')
          mng2 = pyplot.get_current_fig_manager()
          mng2.window.resizable(False, False)
          plt.pause(0.1)
          plt.clf()
          plt.close('all')
          
          return ra_vsa, dec_vsa, y, y_e, j, j_e, h, h_e, ks, ks_e, ymjd, jmjd, hmjd, ksmjd, img, 'Object Not Found was Pressed'
        
        #Allows updates for the circle size slider bar
        elif click_axes == 'Axes(0.25,0.055;0.65x0.03)':
          scatter.remove()
          scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')

      #Checks if the window was closed
      elif press is None:
        ymjd, jmjd, hmjd, ksmjd, y, y_e, j, j_e, h, h_e, ks, ks_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        ra_vsa, dec_vsa = ra, dec
        plt.clf()
        plt.close('all')
        return ra_vsa, dec_vsa, y, y_e, j, j_e, h, h_e, ks, ks_e, ymjd, jmjd, hmjd, ksmjd, img, text_list[text_max]
      
  #Returns null values if the images are not found
  else: 
    ymjd, jmjd, hmjd, ksmjd, y, y_e, j, j_e, h, h_e, ks, ks_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    ra_vsa, dec_vsa = ra, dec
    text_list = 'Image Not Found'
    return ra_vsa, dec_vsa, y, y_e, j, j_e, h, h_e, ks, ks_e, ymjd, jmjd, hmjd, ksmjd, 'VSA', text_list
  
def vsa_table(ra, dec, radius): 
  '''Find all the objects in the radius defined by the user'''

  blockPrint()
    
  #Find the location of all the object found in VSA in the radius choosen by the user
  catalog_list = ['VHS', 'VVV', 'VMC', 'VIKING', 'VIDEO', 'UltraVISTA']
  database_list = ['VHSDR6', 'VVVDR5', 'VMCDR6', 'VIKINGDR5', 'VIDEODR6', 'ULTRAVISTADR4']
  for i in range(len(catalog_list)): 
    table = Vsa.query_region(
        SkyCoord(ra, dec, unit = (u.deg, u.deg), frame = 'fk5'),
                radius = ((radius/2) - 1) * u.arcsec,
                programme_id = catalog_list[i],
                database = database_list[i])
    if len(table) > 0:
      return table