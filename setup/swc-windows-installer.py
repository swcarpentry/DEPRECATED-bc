#!/usr/bin/env python

"""Software Carpentry Windows Installer

Helps mimic a *nix environment on Windows with as little work as possible.

The script:
* Provides standard ipython operation for msysgit
* Provides standard nosetests behavior for msysgit
* Installs nano and makes it accessible from msysgit

To use:

1. Install Anaconda CE Python distribution
   http://continuum.io/anacondace.html
2. Install msysgit
   http://code.google.com/p/msysgit/downloads/list?q=full+installer+official+git
3. Run swc_windows_installer.py
   You should be able to simply double click the file in Windows
"""

import shutil
import zipfile
import os

import requests


def install_nano(python_scripts_directory):
    """Download and install the nano text editor"""
    url = "http://www.nano-editor.org/dist/v2.2/NT/nano-2.2.6.zip"
    r = requests.get(url)
    output_file = open('nano.zip', 'wb')
    output_file.write(r.content)
    output_file.close()
    nano_zip = zipfile.ZipFile('nano.zip')
    nano_files = ['nano.exe', 'cygwin1.dll', 'cygintl-8.dll',
                  'cygiconv-2.dll', 'cyggcc_s-1.dll']
    for file_name in nano_files:
        nano_zip.extract(file_name, '.')
        shutil.move(file_name, python_scripts_directory)
    os.remove('nano.zip')

def create_ipython_entry_point(python_scripts_directory):
    """Creates a terminal-based IPython entry point for msysgit"""
    output_file = open(python_scripts_directory + 'ipython', 'w')
    file_contents = """#!/usr/bin/env python
import sys
from IPython.frontend.html.notebook.notebookapp import launch_new_instance as launch_notebook
from IPython.frontend.terminal.ipapp import launch_new_instance as launch_ipython

def main():
    if len(sys.argv) > 1 and sys.argv[1] == 'notebook':
        sys.exit(launch_notebook())
    else:
        sys.exit(launch_ipython())

if __name__ == '__main__':
    main()
"""

    output_file.write(file_contents)

def create_nosetests_entry_point(python_scripts_directory):
    """Creates a terminal-based nosetests entry point for msysgit"""
    output_file = open(python_scripts_directory + 'nosetests', 'w')
    file_contents = """#!/usr/bin/env/ python
import sys
import nose

if __name__ == '__main__':
    sys.exit(nose.core.main())
"""
    output_file.write(file_contents)

def main():
    python_scripts_directory = "C:\\Anaconda\\Scripts\\"
    #python_scripts_directory = "./scripts/"
    create_ipython_entry_point(python_scripts_directory)
    create_nosetests_entry_point(python_scripts_directory)
    install_nano(python_scripts_directory)


if __name__ == '__main__':
    main()
