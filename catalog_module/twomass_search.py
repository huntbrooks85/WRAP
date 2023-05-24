#Import all of the packages
from catalog_module.importmodule import *

def blockPrint():
  '''Makes a function that blocks the printing function'''
  sys.stdout = open(os.devnull, 'w')

def enablePrint():
  '''Makes a function that allows the printing function'''
  sys.stdout = sys.__stdout__

def twomass_image(ra, dec, radius): 
  ''' First, it gets the images from the 2MASS API from the IRSA Archive and downloads the images. 
  Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
  Third, calls the table function to get all the objects from the 2MASS source catalog. 
  Fourth, makes the window for the user to click the object with all settings. 
  Finally, finds the closest object to the click and records the data. '''
    
  #Makes outline for the window of the plot
  plt.rcParams['toolbar'] = 'None'
  plt.style.use('Solarize_Light2')
  blockPrint()

  #Finds all the metadata that relates to the ra and dec searched, mostly to find the APIs for the W1, W2, W3, and W4 images
  metadata_2mass_link = 'https://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?POS=' + str(ra) + ',' + str(dec) + '&SIZE=' + str(radius/3600)
  twomass_metadata = requests.get(metadata_2mass_link)
  open('Output/metadata/TWOMASS_metadata.txt', 'wb').write(twomass_metadata.content)

  #With this metadata it finds the API link for the W1 and W2 images
  J_finder, h_finder, k_finder = 'All-Sky Release Survey J-Band Atlas Image', 'All-Sky Release Survey H-Band Atlas Image', 'All-Sky Release Survey K-Band Atlas Image'
  j_twomass_image_url_list, h_twomass_image_url_list, k_twomass_image_url_list = [], [], []
  with open('Output/metadata/TWOMASS_metadata.txt', 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      if line.find(J_finder) != -1:
        j_twomass_image_url_list.append(((lines[lines.index(line) + 1]).split('[')[2]).split(']')[0])
      elif line.find(h_finder) != -1:
        h_twomass_image_url_list.append(((lines[lines.index(line) + 1]).split('[')[2]).split(']')[0])
      elif line.find(k_finder) != -1:
        k_twomass_image_url_list.append(((lines[lines.index(line) + 1]).split('[')[2]).split(']')[0])

  #Finds the properly cropped image from 2MASS
  shape_difference = []
  for i in range(len(j_twomass_image_url_list)):
    j_twomass_image_url= j_twomass_image_url_list[i]

    #Download the W1 and W2 images
    file_allwise_j = download_file(j_twomass_image_url, cache = True)
    data_allwise_j = fits.getdata(file_allwise_j)

    #Gets the headers from the images
    hdu_j = fits.open(file_allwise_j)[0]
    wcs = WCS(hdu_j.header)

    #Make a cutout from the coadd image for the RA and DEC put in
    position = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0')
    size = u.Quantity([radius, radius], u.arcsec)

    #Checks if the ra and dec is in the image wcs
    ra_dec_pixel = wcs.world_to_pixel_values(ra, dec)
    if ((0 <= int(ra_dec_pixel[1]) <= data_allwise_j.shape[0]) == False) or ((0 <= int(ra_dec_pixel[0]) <= data_allwise_j.shape[1]) == False): 
      shape_difference.append(10000000)
    else:
      cutout_j = Cutout2D(data_allwise_j, position, size, fill_value = np.nan, wcs = wcs.celestial)

      #Adds the difference between the height and width's absolute value to a list
      shape_difference.append(abs(cutout_j.shape[0] - cutout_j.shape[1]))

  #Checks if the  image is found
  if (all(x == shape_difference[0] for x in shape_difference) == False) and (len(shape_difference) > 0):

    #Finds the best url link
    shape_best = shape_difference.index(np.min(shape_difference))
    j_twomass_image_url, h_twomass_image_url, k_twomass_image_url = j_twomass_image_url_list[shape_best], h_twomass_image_url_list[shape_best], k_twomass_image_url_list[shape_best]

    #Download the W1 and W2 images for the best url link
    file_allwise_j, file_allwise_h, file_allwise_k = download_file(j_twomass_image_url, cache = True), download_file(h_twomass_image_url, cache = True), download_file(k_twomass_image_url, cache = True)
    data_allwise_j, data_allwise_h, data_allwise_k = fits.getdata(file_allwise_j), fits.getdata(file_allwise_h), fits.getdata(file_allwise_k)

    #Gets the headers from the images for the best url link
    hdu_j, hdu_h, hdu_k = fits.open(file_allwise_j)[0], fits.open(file_allwise_h)[0], fits.open(file_allwise_k)[0]
    wcs = WCS(hdu_j.header)

    #Make a cutout from the coadd image for the RA and DEC put in for the best url link
    cutout_j = Cutout2D(data_allwise_j, position, size, fill_value = np.nan, wcs = wcs.celestial)
    cutout_h = Cutout2D(data_allwise_h, position, size, fill_value = np.nan, wcs = wcs.celestial)
    cutout_k = Cutout2D(data_allwise_k, position, size, fill_value = np.nan, wcs = wcs.celestial)
    wcs_cropped = cutout_j.wcs

    #Find the location of all the object found in AllWISE in the radius choosen by the user 
    location_data = twomass_table(ra, dec, radius)
    object_ra, object_dec = location_data['ra'].tolist(), location_data['dec'].tolist()
    j_list, j_list_e = location_data['j_m'].tolist(), location_data['j_msigcom'].tolist()
    h_list, h_list_e = location_data['h_m'].tolist(), location_data['h_msigcom'].tolist()
    ks_list, ks_list_e = location_data['k_m'].tolist(), location_data['k_msigcom'].tolist()
    enablePrint()

    #Gets the dates of when the images were taken
    date_j, date_h, date_k = hdu_j.header[64], hdu_h.header[64], hdu_k.header[64]

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
    total_data = cutout_j.data + cutout_k.data
    ax = plt.subplot(projection = wcs_cropped)

    #Plots the objects found in the radius
    circle_size = (radius*3)
    scatter = ax.scatter(object_ra, object_dec, transform = ax.get_transform('fk5'), s = circle_size, edgecolor = '#40E842', facecolor = 'none')

    #Normalize the image and plots it
    stretch = 99
    norm1_w1 = ImageNormalize(total_data.data, PercentileInterval(stretch), stretch = SinhStretch())
    ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

    #Makes the figure look pretty
    plt.tick_params(axis = 'x', which = 'both', bottom = False, top = False, labelbottom = False)
    plt.tick_params(axis = 'y', which = 'both', bottom = False, top = False, labelbottom = False)
    fontdict_1 = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
    plt.suptitle('2MASS Search', fontsize = 35, y = 0.96, fontfamily = 'Times New Roman')
    ax.set_title('Dates: \n'
              + 'J Date: ' + str(date_j) + ' (YYMMDD)  ' + '  H Date: ' + str(date_h) + ' (YYMMDD)\n'
              + 'K Date: ' + str(date_k) + ' (YYMMDD)  \n', fontdict = fontdict_1, y = 1.04)
    plt.grid(linewidth = 0)
    shape = max(cutout_j.shape)
    plt.xlim(0, shape)
    plt.ylim(0, shape)
    figure = plt.gcf()
    if platform != 'win32':
      figure.set_size_inches(4.75, 6.75)
    elif platform == 'win32':
      figure.set_size_inches(4.75, 7.05)
    figure.canvas.set_window_title('2MASS Search')

    #Make checkbuttons with all of the different image bands
    rax = plt.axes([0.045, 0.4, 0.105, 0.12])
    labels = ['J', 'H', 'K']
    real_data = [cutout_j.data, cutout_h.data, cutout_k.data]
    default = [True, False, True]
    check = CheckButtons(rax, labels, default)

    #Adds a slider for the scaling of the image
    freq_bottom = plt.axes([0.25, 0.1, 0.65, 0.03])
    slider_bottom = Slider(ax = freq_bottom, label = 'Stetch:', valmin = 80, valmax = 100, valinit = stretch, color = '#E48671')

    #Adds a slider for the circle size
    circle_slid_location = plt.axes([0.25, 0.068, 0.65, 0.03])
    circle_slider = Slider(ax = circle_slid_location, label = 'Circle Size:', valmin = (circle_size - 2.5*radius), valmax = (circle_size + 1*radius), valinit = circle_size, color = '#E48671')

    #Adds a notes section that the user can add notes about their data
    axbox = plt.axes([0.15, 0.025, 0.8, 0.03])
    text = ''
    text_box = TextBox(axbox, 'Notes:', initial = text, textalignment = "center")

    #Make a button that can be clicked if no object is found
    axes_button = plt.axes([0.04, 0.78, 0.92, 0.04])
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
        if default != [True, True, True]: 
          total_data = real_data[1]*0
        if default[d] == True: 
          total_data = total_data + real_data[d]
        else: 
          pass
      norm1_w1 = ImageNormalize(total_data.data, PercentileInterval(slider_bottom.val), stretch = SinhStretch())
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
      norm1_w1 = ImageNormalize(total_data.data, PercentileInterval(slider_bottom.val), stretch = SinhStretch())
      ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

    #Updates the notes added by the user when there is an input
    text_list = [text]
    def submit(expression):
      '''Updates the list of types in the 'Notes' setting'''

      text = expression
      text_list.append(text)   

    #Allows the sliders and buttons to be pressed
    check.on_clicked(update_button)
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
          shape_x, shape_y = total_data.shape[0], total_data.shape[1]
          ax.text(shape_x/20, shape_y/5, 'Your Click Has Been Successfully Recorded for 2MASS! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
          plt.pause(0.1)
          plt.clf()
          plt.close('all')

          #Find the closest point to the location clicked to obtain W1, W2, W3, and W4 photometry
          coord = wcs_cropped.pixel_to_world_values(location[n-4],location[n-5])
          distance = []
          for i in range(len(object_ra)):
            distance.append(math.dist(coord, [object_ra[i], object_dec[i]]))

          list_location = distance.index(np.min(distance))
          ra_2mass, dec_2mass = object_ra[list_location], object_dec[list_location]
          j, j_e = j_list[list_location], j_list_e[list_location]
          h, h_e = h_list[list_location], h_list_e[list_location]
          ks, ks_e = ks_list[list_location], ks_list_e[list_location]
          return ra_2mass, dec_2mass, j, j_e, h, h_e, ks, ks_e, '2MASS All-Sky Point Source Catalog', text_list[text_max]
        
        #Checks if the Object not Found button was clicked
        elif click_axes == 'Axes(0.04,0.78;0.92x0.04)':
          j, j_e, h, h_e, ks, ks_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
          ra_2mass, dec_2mass = ra, dec
          shape_x, shape_y = total_data.shape[0], total_data.shape[1]
          ax.text(shape_x/20, shape_y/5, 'Your Click Has Been Successfully Recorded for 2MASS! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
          plt.pause(0.1)
          plt.clf()
          plt.close('all')
          return ra_2mass, dec_2mass, j, j_e, h, h_e, ks, ks_e, '2MASS All-Sky Point Source Catalog', 'Object Not Found was Pressed'
        
        #Adds the functionality of the circle slider bar
        elif click_axes == 'Axes(0.25,0.068;0.65x0.03)':
          scatter.remove()
          scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')

      #Checks if the window was closed
      elif press is None:
        j, j_e, h, h_e, ks, ks_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        ra_2mass, dec_2mass = ra, dec
        plt.close('all')
        plt.figure().clear()
        return ra_2mass, dec_2mass, j, j_e, h, h_e, ks, ks_e, '2MASS All-Sky Point Source Catalog', text_list[text_max]
  else: 
      j, j_e, h, h_e, ks, ks_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
      ra_2mass, dec_2mass = ra, dec
      text_list = 'Image Not Found'
      plt.close('all')
      plt.figure().clear()
      return ra_2mass, dec_2mass, j, j_e, h, h_e, ks, ks_e, '2MASS All-Sky Point Source Catalog', text_list

def twomass_table(ra, dec, radius): 
  '''Find all the objects in the radius defined by the user'''
  
  blockPrint()

  #Uses astroquery to find all objects in the radius
  location_data = Irsa.query_region(coord.SkyCoord(ra, dec, unit = (u.deg,u.deg), frame = 'fk5'), catalog = 'fp_psc', spatial = 'Box', width = (radius - 1) * u.arcsec)
  return location_data