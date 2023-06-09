#Import all of the packages
from catalog_module.importmodule import *

def blockPrint():
  '''Makes a function that blocks the printing function'''
  sys.stdout = open(os.devnull, 'w')

def enablePrint():
  '''Makes a function that allows the printing function'''
  sys.stdout = sys.__stdout__

def allwise_image(ra, dec, radius): 
  ''' First, it gets the images from the AllWISE API from IRSA and downloads the images. 
  Second, crops the images, around the RA and DEC from the user and grabs relavent data to the image. 
  Third, calls the table function to get all the objects from the AllWISE source catalog. 
  Fourth, makes the window for the user to click the object with all settings. 
  Finally, finds the closest object to the click and records the data. '''

  #Makes outline for the window of the plot
  plt.rcParams['toolbar'] = 'None'
  matplotlib.use("TkAgg")
  plt.style.use('Solarize_Light2')
  blockPrint()

  #Finds all the metadata that relates to the ra and dec searched, mostly to find the APIs for the W1, W2, W3, and W4 images
  metadata_allwise_link = 'http://irsa.ipac.caltech.edu/ibe/sia/wise/allwise/p3am_cdd?POS=' + str(ra) + ',' + str(dec) + '&SIZE=' + str(radius/3600)
  allwise_metadata = requests.get(metadata_allwise_link)
  open('Output/metadata/AllWISE_metadata.txt', 'wb').write(allwise_metadata.content)

  #With this metadata it finds the API link for the W1 and W2 images
  w1_finder, w2_finder, w3_finder, w4_finder = 'W1 Coadd', 'W2 Coadd', 'W3 Coadd', 'W4 Coadd'
  with open('Output/metadata/AllWISE_metadata.txt', 'r') as fp:
   lines = fp.readlines()
   for line in lines:
      if line.find(w1_finder) != -1:
        w1_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]
      elif line.find(w2_finder) != -1:
        w2_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]
      elif line.find(w3_finder) != -1:
        w3_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]
      elif line.find(w4_finder) != -1:
        w4_allwise_image_url = ((lines[lines.index(line) + 1]).split('>', 1)[1]).split('</', 1)[0]

  #Download the W1 and W2 images
  file_allwise_w1, file_allwise_w2 = download_file(w1_allwise_image_url, cache=True), download_file(w2_allwise_image_url, cache=True)
  file_allwise_w3, file_allwise_w4 = download_file(w3_allwise_image_url, cache=True), download_file(w4_allwise_image_url, cache=True) 

  data_allwise_w1, data_allwise_w2 =  fits.getdata(file_allwise_w1), fits.getdata(file_allwise_w2)
  data_allwise_w3, data_allwise_w4 = fits.getdata(file_allwise_w3), fits.getdata(file_allwise_w4)

  #Find the location of all the object found in AllWISE in the radius choosen by the user 
  location_data = allwise_table(ra, dec, radius)
  object_ra, object_dec = location_data['ra'].tolist(), location_data['dec'].tolist()
  object_ra_sigma, object_dec_sigma = location_data['sigra'].tolist(), location_data['sigdec'].tolist()
  w1_list, w2_list = location_data['w1mpro'].tolist(), location_data['w2mpro'].tolist()
  w3_list, w4_list = location_data['w3mpro'].tolist(), location_data['w4mpro'].tolist()
  w1_list_sigma, w2_list_sigma = location_data['w1sigmpro'].tolist(), location_data['w2sigmpro'].tolist()
  w3_list_sigma, w4_list_sigma = location_data['w3sigmpro'].tolist(), location_data['w4sigmpro'].tolist()
  pmra_list, pmdec_list = location_data['pmra'].tolist(), location_data['pmdec'].tolist()
  pmra_sigma_list, pmdec_sigma_list = location_data['sigpmra'].tolist(), location_data['sigpmdec'].tolist()

  #Obtains the headers for each image
  hdu_w1, hdu_w2 = fits.open(file_allwise_w1)[0], fits.open(file_allwise_w2)[0]
  hdu_w3, hdu_w4 = fits.open(file_allwise_w3)[0], fits.open(file_allwise_w4)[0]
  wcs = WCS(hdu_w1.header)

  #Make a cutout from the coadd image for the RA and DEC put in
  position, size = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5', equinox = 'J2000.0'), u.Quantity([radius, radius], u.arcsec)
  cutout_w1 = Cutout2D(data_allwise_w1, position, size, fill_value = np.nan, wcs = wcs.celestial)
  cutout_w2 = Cutout2D(data_allwise_w2, position, size, fill_value = np.nan, wcs = wcs.celestial)
  cutout_w3 = Cutout2D(data_allwise_w3, position, size, fill_value = np.nan, wcs = wcs.celestial)
  cutout_w4 = Cutout2D(data_allwise_w4, position, size, fill_value = np.nan, wcs = wcs.celestial)
  wcs_cropped = cutout_w1.wcs
  enablePrint()

  #Obtains the dates for each image
  date_w1, date_w2 = hdu_w1.header['MIDOBS'].split('T', 2)[0], hdu_w2.header['MIDOBS'].split('T', 2)[0]
  date_w3, date_w4 = hdu_w3.header['MIDOBS'].split('T', 2)[0], hdu_w4.header['MIDOBS'].split('T', 2)[0]

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
  total_data = cutout_w1.data + cutout_w2.data
  ax = plt.subplot(projection = wcs_cropped)

  #Plots the objects found in the radius
  circle_size = (radius*3)
  scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('fk5'), s = circle_size, edgecolor='#40E842', facecolor='none')

  #Normalize the image and plots it
  init_top, init_bot = 95, 45
  norm1_w1 = matplotlib.colors.Normalize(vmin = np.nanpercentile(total_data.data, init_bot), vmax = np.nanpercentile(total_data.data, init_top))
  ax.imshow(total_data.data, cmap = 'Greys', norm = norm1_w1)

  #Makes the figure look pretty
  plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
  plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
  fontdict = {'family':'Times New Roman','color':'k','size':11, 'style':'italic'}
  plt.suptitle('AllWISE Search', fontsize = 35, y = 0.96, fontfamily = 'Times New Roman')
  ax.set_title('Dates: \n'
             + 'W1 Date: ' + str(date_w1) + ' (Y-M-D)  ' + '  W2 Date: ' + str(date_w2) + ' (Y-M-D)\n'
             + 'W3 Date: ' + str(date_w3) + ' (Y-M-D)  ' + '  W4 Date: ' + str(date_w4) + ' (Y-M-D)\n', fontdict = fontdict, y = 1.05)
  plt.grid(linewidth = 0)
  figure = plt.gcf()
  plt.xlim(0, max(total_data.shape))
  plt.ylim(0, max(total_data.shape))
  if platform != 'win32':
    figure.set_size_inches(4.75, 7.25)
  elif platform == 'win32':
    figure.set_size_inches(4.75, 7.25)
  figure.canvas.set_window_title('AllWISE Search')
  mng = pyplot.get_current_fig_manager()
  mng.window.resizable(False, False)

  #Make checkbuttons with all of the different image bands
  rax = plt.axes([0.045, 0.4, 0.115, 0.1])
  labels = ['W1', 'W2', 'W3', 'W4']
  real_data = [cutout_w1.data, cutout_w2.data, cutout_w3.data, cutout_w4.data]
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

  #Updates the image depending on what the user chooses
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
          
        #Makes a pop-up window with success text
        plt.clf()
        plt.close('all')
        plt.figure(1)
        plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for AllWISE! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
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
        figure2.canvas.set_window_title('Successful AllWISE Search')
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

        #Obtains all of the clicked objects data
        ra_allwise, dec_allwise = object_ra[list_location], object_dec[list_location]
        ra_aw_e, dec_aw_e = object_ra_sigma[list_location], object_dec_sigma[list_location]
        w1, w2 = w1_list[list_location], w2_list[list_location]
        w3, w4 = w3_list[list_location], w4_list[list_location]
        w1_sigma, w2_sigma = w1_list_sigma[list_location], w2_list_sigma[list_location]
        w3_sigma, w4_sigma = w3_list_sigma[list_location], w4_list_sigma[list_location]
        pmra, pmdec = pmra_list[list_location], pmdec_list[list_location]
        pmra_sigma, pmdec_sigma = pmra_sigma_list[list_location], pmdec_sigma_list[list_location]
        return ra_allwise, ra_aw_e, dec_allwise, dec_aw_e, w1, w1_sigma, w2, w2_sigma, w3, w3_sigma, w4, w4_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma, 'AllWISE Source Catalog', text_list[text_max]
      
      #Checks if the "Object Not Found" button was clicked
      elif click_axes == 'Axes(0.04,0.775;0.92x0.04)':
        ra_aw_e, dec_aw_e, w1, w1_sigma, w2, w2_sigma, w3, w3_sigma, w4, w4_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        ra_allwise = ra
        dec_allwise = dec

        #Makes a pop-up window with success text
        plt.clf()
        plt.close('all')
        plt.figure(1)
        plt.text(0.06, 0.25, 'Your Click Has Been Successfully Recorded for AllWISE! \n              Please Wait for the Next Catalog to Load!', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
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
        figure2.canvas.set_window_title('Successful AllWISE Search')
        mng2 = pyplot.get_current_fig_manager()
        mng2.window.resizable(False, False)
        plt.pause(0.1)
        plt.clf()
        plt.close('all')

        return ra_allwise, ra_aw_e, dec_allwise, dec_aw_e, w1, w1_sigma, w2, w2_sigma, w3, w3_sigma, w4, w4_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma, 'AllWISE Source Catalog', 'Object Not Found was Pressed'
      
      #Updates the circle size when slider is moved
      elif click_axes == 'Axes(0.25,0.055;0.65x0.03)':
        scatter.remove()
        scatter = ax.scatter(object_ra, object_dec, transform = ax.get_transform('fk5'), s = circle_slider.val, edgecolor='#40E842', facecolor='none')

    #Checks if the window was closed
    elif press is None:
      ra_aw_e, dec_aw_e, w1, w1_sigma, w2, w2_sigma, w3, w3_sigma, w4, w4_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
      ra_allwise = ra
      dec_allwise = dec
      plt.close('all')
      plt.figure().clear()
      return ra_allwise, ra_aw_e, dec_allwise, dec_aw_e, w1, w1_sigma, w2, w2_sigma, w3, w3_sigma, w4, w4_sigma, pmra, pmra_sigma, pmdec, pmdec_sigma, 'AllWISE Source Catalog', text_list[text_max]
  
def allwise_table(ra, dec, radius): 
  '''Find all the objects in the radius defined by the user'''

  blockPrint()
 
  #Uses astroquery to find all objects in the radius
  location_data = Irsa.query_region(coord.SkyCoord(ra, dec, unit = (u.deg,u.deg), frame = 'fk5'), catalog = 'allwise_p3as_psd', spatial = 'Box', width = (radius - 1) * u.arcsec)
  return location_data
