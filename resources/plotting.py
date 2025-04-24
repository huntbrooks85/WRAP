#-----------------------------------------------------------------------#
# WRAP.plotting v2.0.0
# By Hunter Brooks, at NAU/UToledo, Flagstaff: April 23, 2025
#-----------------------------------------------------------------------#



# Import all needed packages.
# ------------------------------------------------------------- #
# Math Packages
import math
import numpy as np

# System Packages
from sys import platform

# Image Plotting Packagaes
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
from matplotlib.widgets import Slider, Button, CheckButtons, TextBox

# NOIRLab Source Catatlog Package
if platform != 'win32':
  from dl.helpers.utils import convert
  from dl import authClient as ac, queryClient as qc
# ------------------------------------------------------------- #



# WRAP PLOTTING FUNCTIONS
# ------------------------------------------------------------- #   
def image_plot(ra, dec, radius, catalog_info, table, images, w): 
  # Check if queries returned valid results
  if type(images) != int and type(w) != int and type(table) != int: 
    try: 
      # Configure plot appearance
      plt.rcParams['toolbar'] = 'None'
      plt.style.use('bmh')
      plt.rcParams["figure.figsize"] = [8, 8]

      # Create a new figure with WCS projection
      fig_1, ax = plt.subplots(subplot_kw={'projection': w})
      fig_1.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1)

      # Hide axis ticks and labels
      plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
      plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
      plt.suptitle(f'{catalog_info["name"]} Search', fontsize=35, y=0.98, fontfamily='Times New Roman')
      
      plt.grid(linewidth = 0)
      figure = plt.gcf()
      figure.set_size_inches(4.75, 6)

      # Extract object coordinates from the table
      try: 
        object_ra = table[0][catalog_info['table_data'][0]].tolist()
        object_dec = table[0][catalog_info['table_data'][1]].tolist()
      except:
        object_ra = table[catalog_info['table_data'][0]].tolist()
        object_dec = table[catalog_info['table_data'][1]].tolist()      
      
      # Combine image data
      default = catalog_info['image_selection']
      if isinstance(images[0], np.ndarray): 
        total_data = np.zeros_like(images[0].data)
        for index in range(len(images)):
          if default[index] == True: 
            total_data += images[index].data
        max_shape = np.nanmax(total_data.shape)
      else:
        total_data = np.zeros_like(images[0][0].data)
        for image in images:
          index = images.index(image)
          if default[index] == True: 
            total_data += image[0].data
        max_shape = np.nanmax(total_data.shape)
              
      # Set initial circle size for scatter plot
      circle_size = (radius * 3)
      scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('icrs'), s=circle_size, linewidths=2, edgecolor='#40E842', facecolor='none')

      # Set initial normalization for image display
      init_bot, init_top = 45, 95
      norm1_total = matplotlib.colors.Normalize(vmin=np.nanpercentile(total_data.data, init_bot), vmax=np.nanpercentile(total_data.data, init_top))
      ax.imshow(total_data.data, cmap='Greys', norm=norm1_total)
      
      if catalog_info['name'] == 'VISTA': 
        plt.xlim(0, max_shape), plt.ylim(max_shape, 0) 
      elif catalog_info['name'] == 'NSC': 
        plt.xlim(max_shape, 0), plt.ylim(0, max_shape) 
      else: 
        plt.xlim(0, max_shape), plt.ylim(0, max_shape) 

      # Create a cursor for the plot
      cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
      
      # Create an annotation for mouse events
      annotation = ax.annotate('', xy=(0, 0), xytext=(20, -20), arrowprops=dict(arrowstyle='wedge'), fontsize=12, color='red')
      annotation.set_visible(False)

      # Create sliders for adjusting image stretch
      freq_top = plt.axes([0.25, 0.12, 0.65, 0.03])
      slider_top = Slider(ax=freq_top, label='Top Stretch:', valmin=50, valmax=100, valinit=init_top, color='#E48671')
      freq_bottom = plt.axes([0.25, 0.087, 0.65, 0.03])
      slider_bottom = Slider(ax=freq_bottom, label='Bottom Stretch:', valmin=0, valmax=50, valinit=init_bot, color='#E48671')

      # Create slider for adjusting circle size
      circle_slid_location = plt.axes([0.25, 0.055, 0.65, 0.03])
      circle_slider = Slider(ax=circle_slid_location, label='Circle Size:', valmin=(circle_size - 2.5*radius), valmax=(circle_size + 1*radius), valinit=circle_size, color='#E48671')

      # Create a text box for notes
      axbox = plt.axes([0.15, 0.02, 0.8, 0.03])
      text = 'No Inputted Text'
      text_box = TextBox(axbox, 'Notes:', initial=text, textalignment="center")

      # Create a button for "Object Not Found"
      axes_button = plt.axes([0.04, 0.855, 0.92, 0.04])
      close = Button(axes_button, 'Object Not Found', color='#E48671')
      
      #Make checkbuttons with all of the different image bands
      rax = plt.axes([0.045, 0.4, 0.105, 0.12])
      labels = catalog_info['image_names']
      real_data = []
      if isinstance(images[0], np.ndarray): 
        for image in images:
          real_data.append(image.data)
      else:
        for image in images:
          real_data.append(image[0].data)
      check = CheckButtons(rax, labels, default)

      # Track mouse click locations
      location = []
      def mouse_event(event):
        location.append(event.xdata)
        location.append(event.ydata)
        location.append(event.inaxes)
      cid = fig_1.canvas.mpl_connect('button_press_event', mouse_event)
      
      # Update image stretch based on slider values
      def update_slider_stretch(val):
        if isinstance(images[0], np.ndarray):
          total_data = np.zeros_like(images[0].data)
        else:
          total_data = np.zeros_like(images[0][0].data)       
          
        for tf_ind in default: 
          if tf_ind == False: 
            pass
          if tf_ind == True: 
            index = default.index(tf_ind)
            if isinstance(images[0], np.ndarray):
              total_data += images[index].data
            else:
              total_data += images[index][0].data
                    
        norm1_w1 = matplotlib.colors.Normalize(vmin=np.nanpercentile(total_data.data, slider_bottom.val), vmax=np.nanpercentile(total_data.data, slider_top.val))
        ax.imshow(total_data.data, cmap='Greys', norm=norm1_w1)

      # Update text from text box
      text_list = [text]
      def submit(expression):
        text = expression
        text_list.append(text)
        
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
              
        if isinstance(images[0], np.ndarray):
          total_data = np.zeros_like(images[0].data)
        else:
          total_data = np.zeros_like(images[0][0].data)       
          
        for tf_ind in default: 
          if tf_ind == False: 
            pass
          if tf_ind == True: 
            index = default.index(tf_ind)
            if isinstance(images[0], np.ndarray):
              total_data += images[index].data
            else:
              total_data += images[index][0].data
        norm1_total = matplotlib.colors.Normalize(vmin=np.nanpercentile(total_data.data, init_bot), vmax=np.nanpercentile(total_data.data, init_top))
        ax.imshow(total_data.data, cmap='Greys', norm=norm1_total)

      # Connect sliders and text box to their update functions
      slider_top.on_changed(update_slider_stretch)
      slider_bottom.on_changed(update_slider_stretch)
      text_box.on_text_change(submit)
      check.on_clicked(update_button)

      # Main loop to handle user interactions
      n = -1
      while True: 
        press = plt.waitforbuttonpress()
        text_max = len(text_list) - 1
        
        if press == False:
          n += 3

          # Find which axes was clicked
          click_axes = str(location[n])
          click_axes = click_axes.split('WCSAxesSubplot', 2)[0]
          
          # Handle "Object Not Found" button click
          if click_axes == 'Axes(0.04,0.855;0.92x0.04)':
            next_window()
            return len(catalog_info['table_data'])*[np.nan] + [text_list[text_max]]
          
          # Handle circle size slider adjustment
          if click_axes == 'Axes(0.25,0.055;0.65x0.03)': 
            scatter.remove()
            scatter = ax.scatter(object_ra, object_dec, transform=ax.get_transform('icrs'), s=circle_slider.val, linewidth=2, edgecolor='#40E842', facecolor='none')
          
          # Handle main plot click
          if click_axes == '': 
            # Find the closest point to the location clicked
            coord = w.pixel_to_world_values(location[len(location) - 3], location[len(location) - 2])
            distance = []
            for i in range(len(object_ra)):
              distance.append(math.dist(coord, [float(object_ra[i]), float(object_dec[i])]))
            list_location = distance.index(np.nanmin(distance))
            table_data = []
            for col_name in catalog_info['table_data']: 
              try: 
                col_data = table[0][f'{col_name}'].tolist()
                table_data.append(col_data[list_location])
              except: 
                col_data = table[f'{col_name}'].tolist()
                table_data.append(col_data[list_location]) 
            next_window()             
            return table_data + [text_list[text_max]]
        
        elif press is None:
          next_window()
          return len(catalog_info['table_data'])*[np.nan] + [text_list[text_max]]
    except Exception as e: 
      return len(catalog_info['table_data'])*[np.nan] + ['Catalog Data Not Retrieved']
  else: 
    return len(catalog_info['table_data'])*[np.nan] + ['Catalog Data Not Retrieved']
# ------------------------------------------------------------- #



# Pop-up window when plot is clicked
# ------------------------------------------------------------- #
def next_window(): 
  # Clear the current figure and close all plots
  plt.clf(), plt.close('all')
  
  # Create a new figure
  plt.figure(1)
  
  # Display a message in the plot
  plt.text(0.04, 0.4, 'Your Click Has Been Recorded \n       Loading Next Catalog', style='oblique', bbox={'facecolor': '#40E842', 'alpha': 1, 'pad': 10})
  
  # Set plot limits and disable the grid
  plt.xlim(0, 1), plt.ylim(0, 1), plt.grid(linewidth=0)
  
  # Get current axis and hide the tick labels and ticks
  ax = plt.gca()
  ax.xaxis.set_tick_params(labelbottom=False)
  ax.yaxis.set_tick_params(labelleft=False)
  ax.set_xticks([])
  ax.set_yticks([])
  
  # Get the current figure and set its size
  figure2 = plt.gcf()
  figure2.set_size_inches(3, 0.5)
  
  # Pause briefly to display the figure, then clear and close it
  plt.pause(0.1), plt.clf(), plt.close('all')
  return
# ------------------------------------------------------------- #
