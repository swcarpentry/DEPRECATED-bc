#!/usr/bin/env python

"""Software Carpentry Windows Installer

Helps mimic a *nix environment on Windows with as little work as possible.

The script:
* Installs nano and makes it accessible from msysgit
* Provides standard nosetests behavior for msysgit

To use:

1. Install Python, IPython, and Nose.  An easy way to do this is with
   the Anaconda CE Python distribution
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
import os
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
    os.makedirs(install_directory)
    for file_name in nano_files:
        nano_zip.extract(file_name, install_directory)


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


def update_bash_profile(extra_paths=()):
    """Create or append to a .bash_profile for Software Carpentry

    Adds nano to the path, sets the default editor to nano, and adds
    additional paths for other executables.
    """
    lines = [
        '',
        '# Add paths for Software-Carpentry-installed scripts and executables',
        'export PATH=$PATH:{}'.format(':'.join(
            make_posix_path(path) for path in extra_paths),),
        '',
        '# Make nano the default editor',
        'export EDITOR=nano',
        '',
        ]
    config_path = os.path.join(os.path.expanduser('~'), '.bash_profile')
    with open(config_path, 'a') as f:
        f.write('\n'.join(lines))

def make_posix_path(windows_path):
    """Convert a Windows path to a posix path"""
    return windows_path.replace('\\', '/').replace('C:', '/c')


def main():
    swc_dir = os.path.join(os.path.expanduser('~'), '.swc')
    bin_dir = os.path.join(swc_dir, 'bin')
    create_nosetests_entry_point(python_scripts_directory=bin_dir)
    nano_dir = os.path.join(swc_dir, 'lib', 'nano')
    install_nano(installation_directory=nano_dir)
    update_bash_profile(extra_paths=(nano_dir, bin_dir))


if __name__ == '__main__':
    main()
