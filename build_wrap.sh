#!/bin/bash

# Check if python3.8 exists
if command -v python3.8 &> /dev/null; then
    echo "✅ Python 3.8 is already installed: $(python3.8 --version)"
else
    echo "🔍 Python 3.8 not found. Installing with Homebrew..."
    
    # Check if Homebrew is installed
    if ! command -v brew &> /dev/null; then
        echo "❌ Homebrew is not installed. Please install Homebrew first: https://brew.sh/"
        exit 1
    fi

    # Install Python 3.8 using brew
    brew install python@3.8

    # Add it to PATH (if not already linked)
    if ! command -v python3.8 &> /dev/null; then
        echo "🔧 Adding Python 3.8 to PATH..."
        echo 'export PATH="/usr/local/opt/python@3.8/bin:$PATH"' >> ~/.bash_profile
        export PATH="/usr/local/opt/python@3.8/bin:$PATH"
    fi

    echo "✅ Python 3.8 installed successfully: $(python3.8 --version)"
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
  --app-version 2.1.0

echo "✅ Build and package complete!"
