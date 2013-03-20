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

try:  # Python 3
    from io import BytesIO as _BytesIO
except ImportError:  # Python 2
    from StringIO import StringIO as _BytesIO
import os.path
try:  # Python 3
    from urllib.request import urlopen as _urlopen
except ImportError:  # Python 2
    from urllib2 import urlopen as _urlopen
import zipfile


def install_nano(install_directory):
    """Download and install the nano text editor"""
    url = "http://www.nano-editor.org/dist/v2.2/NT/nano-2.2.6.zip"
    r = _urlopen(url)
    nano_zip_content = _BytesIO(r.read())
    nano_zip = zipfile.ZipFile(nano_zip_content)
    nano_files = ['nano.exe', 'cygwin1.dll', 'cygintl-8.dll',
                  'cygiconv-2.dll', 'cyggcc_s-1.dll']
    for file_name in nano_files:
        nano_zip.extract(file_name, install_directory)

def create_ipython_entry_point(python_scripts_directory):
    """Creates a terminal-based IPython entry point for msysgit"""
    contents = '\n'.join([
            '#!/usr/bin/env python',
            'from IPython.frontend.terminal.ipapp import launch_new_instance',
            'launch_new_instance()',
            '',
            ])
    with open(os.path.join(python_scripts_directory, 'ipython'), 'w') as f:
        f.write(contents)

def create_nosetests_entry_point(python_scripts_directory):
    """Creates a terminal-based nosetests entry point for msysgit"""
    contents = '\n'.join([
            '#!/usr/bin/env/ python',
            'import sys',
            'import nose',
            "if __name__ == '__main__':",
            '    sys.exit(nose.core.main())',
            '',
            ])
    with open(os.path.join(python_scripts_directory, 'nosetests'), 'w') as f:
        f.write(contents)


def main():
    python_scripts_directory = "C:\\Anaconda\\Scripts\\"
    #python_scripts_directory = "./scripts/"
    create_ipython_entry_point(python_scripts_directory)
    create_nosetests_entry_point(python_scripts_directory)
    install_nano(python_scripts_directory)


if __name__ == '__main__':
    main()
