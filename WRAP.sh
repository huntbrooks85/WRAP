#!/bin/bash

# Prints out opening
echo "Opening WRAP"

CONDA_ENV_NAME="WRAP"
PIP_MODULE="Output/metadata/pip_module/requirements-mac.txt"

# Set the path to Miniconda installation directory
MINICONDA_DIR="$HOME/miniconda"

# Function to install Miniconda
install_miniconda() {
    echo "Installing Miniconda..."
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS
        curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
        bash Miniconda3-latest-MacOSX-x86_64.sh -b -p $MINICONDA_DIR
        rm Miniconda3-latest-MacOSX-x86_64.sh
        source "$MINICONDA_DIR/etc/profile.d/conda.sh"
    elif [[ "$OSTYPE" == "msys" || "$OSTYPE" == "win32" ]]; then
        # Windows
        powershell.exe -Command "& {Invoke-WebRequest https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe -OutFile Miniconda3-latest-Windows-x86_64.exe}"
        Miniconda3-latest-Windows-x86_64.exe /S /D=$MINICONDA_DIR
        rm Miniconda3-latest-Windows-x86_64.exe
        source "$MINICONDA_DIR/Scripts/activate"
    else
        echo "Unsupported operating system: $OSTYPE"
        exit 1
    fi
}

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    install_miniconda
else
    echo "Conda is already installed."
fi

# Create and activate Conda environment
if ! conda env list | grep -q "$CONDA_ENV_NAME"; then
    echo "Creating Conda environment: $CONDA_ENV_NAME"
    conda create --name $CONDA_ENV_NAME python=3.8 git -y
fi

echo "Activating Conda environment: $CONDA_ENV_NAME"
source activate $CONDA_ENV_NAME

# Install pip_module if the environment doesn't already have it
if ! python -c "import $PIP_MODULE" &> /dev/null; then
    echo "Installing $PIP_MODULE"
    pip install -r $PIP_MODULE
else
    echo "$PIP_MODULE is already installed in the Conda environment."
fi

# Run your script or command here
clear
echo "Opened WRAP"
python Output/metadata/WRAP.py
echo "Closing WRAP"
