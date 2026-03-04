#!/bin/bash

# Check if python3.8 exists
if command -v python3.8 &> /dev/null; then
    echo "✅ Python 3.8 is already installed: $(python3.8 --version)"
else
    echo "🔍 Python 3.8 not found. Please install Python 3.8.8 at https://www.python.org/downloads/release/python-388/"
fi

set -e  # Exit on error

# Create virtual environment
echo "🔧 Creating Python virtual environment at ./resources/myenv ..."
python3.8 -m venv ./resources/myenv

# Activate virtual environment and install requirements
echo "📦 Installing Python packages from requirements.txt ..."
source ./resources/myenv/bin/activate
pip install --upgrade pip
pip install -r ./resources/requirements.txt
deactivate

# Clean previous Java builds
echo "🧹 Removing old .class files and WRAP.jar if exist..."
rm -f *.class WRAP.jar

# Compile Java
echo "🧵 Compiling WRAP.java ..."
javac WRAP.java

# Build jar
echo "📦 Creating WRAP.jar ..."
jar cmfv MainClass.txt WRAP.jar *.class

# Package app using jpackage
echo "📦 Running jpackage ..."
jpackage \
  --name WRAP \
  --input . \
  --main-jar WRAP.jar \
  --resource-dir resources \
  --type pkg \
  --icon Icons/WRAPLogo.icns \
  --app-version 2.1.1

echo "✅ Build and package complete!"
