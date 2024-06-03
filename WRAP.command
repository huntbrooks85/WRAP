#!/bin/bash
# Get the directory of the current script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"

# Function to close terminal window
close_terminal_window() {
    osascript -e "tell application \"Terminal\" to close (every window whose id is $1)"
}

# Check if Python 3.8 is installed
if ! command -v python3.8 &> /dev/null; then
    echo "Python 3.8 is not installed. Installing..."
    # You might need to provide instructions or a link to a standalone installer for macOS
    echo "Please install Python 3.8 using a standalone installer or a package manager like Homebrew."
    exit 1
fi

# Install the required Python packages
pip3.8 install astro_datalab==2.20.1
pip3.8 install astropy==5.2.2
pip3.8 install astroquery==0.4.7
pip3.8 install beautifulsoup4==4.11.1
pip3.8 install matplotlib==3.5.0
pip3.8 install numpy==1.22.0
pip3.8 install opencv_python==4.7.0.72
pip3.8 install pandas==1.5.3
pip3.8 install PySimpleGUI==4.60.5
pip3.8 install pyvo==1.4
pip3.8 install Requests==2.32.3
pip3.8 install numpy==1.22.0


# Get the ID of the terminal window
WINDOW_ID=$(osascript -e 'tell application "Terminal" to id of window 1')

# Clear the terminal
clear

# Run the Python script located in the same directory
python3.8 "$SCRIPT_DIR/WRAP.py"

# Close the terminal window after the script finishes
close_terminal_window $WINDOW_ID
