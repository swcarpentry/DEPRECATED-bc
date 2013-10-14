#!/usr/bin/env python

"""Software Carpentry Nano Installer for Windows

Installs nano and makes it the default editor in msysgit

To use:

1. Install Python
2. Install msysgit
   http://code.google.com/p/msysgit/downloads/list?q=full+installer+official+git
3. Run swc_nano_installer.py
   You should be able to simply double click the file in Windows

This is a stripped down version of swc_windows_installer.py
originally written by Ethan White and W. Trevor Price.

"""

try:  # Python 3
    from io import BytesIO as _BytesIO
except ImportError:  # Python 2
    from StringIO import StringIO as _BytesIO
import os
try:  # Python 3
    from urllib.request import urlopen as _urlopen
except ImportError:  # Python 2
    from urllib2 import urlopen as _urlopen
import zipfile


def install_nano(install_dir):
    """Download and install the nano text editor"""
    url = "http://www.nano-editor.org/dist/v2.2/NT/nano-2.2.6.zip"
    r = _urlopen(url)
    nano_zip_content = _BytesIO(r.read())
    nano_zip = zipfile.ZipFile(nano_zip_content)
    nano_files = ['nano.exe', 'cygwin1.dll', 'cygintl-8.dll',
                  'cygiconv-2.dll', 'cyggcc_s-1.dll']
    for file_name in nano_files:
        nano_zip.extract(file_name, install_dir)

def make_bash_profile(home_dir, nano_dir):
    """Creates a .bash_profile file for nano setup

    Adds nano to the path and sets the default editor to nano

    """

    nano_path = make_posix_path(nano_dir)
    contents = '\n'.join(['',
                          '# Add nano to path and set as default editor',
                          '# Added by the Software Carpentry nano installer',
                          'export PATH=$PATH:%s' % nano_path,
                          'export EDITOR=nano',
                          ''])
    with open(os.path.join(home_dir, '.bash_profile'), 'a') as f:
        f.write(contents)

def make_posix_path(windows_path):
    """Convert a Windows path to a posix path"""
    return windows_path.replace('\\', '/').replace('C:', '/c')

def main():
    home_dir = os.path.expanduser("~")
    nano_dir = os.path.join(home_dir, '.nano')
    if not os.path.exists(nano_dir):
        os.makedirs(nano_dir)
    install_nano(nano_dir)
    make_bash_profile(home_dir, nano_dir)

if __name__ == '__main__':
    main()
