#!/bin/bash

CONDA_ENV_NAME="WRAP"

if [[ "$OSTYPE" == "darwin"* ]]; then
    PIP_MODULE="Output/metadata/pip_module/requirements-mac.txt"
    echo "Running on macOS"
elif [[ "$OSTYPE" == "msys" || "$OSTYPE" == "win32" ]]; then
    PIP_MODULE="Output/metadata/pip_module/requirements-windows.txt"
    echo "Running on Windows"
else
    echo "Unsupported operating system: $OSTYPE"
    exit 1
fi

echo "Activating Conda environment: $CONDA_ENV_NAME"
source activate $CONDA_ENV_NAME

pip install -r $PIP_MODULE

# Run your script or command here
clear
echo "Opened WRAP"
python Output/metadata/WRAP.py
echo "Closing WRAP"
