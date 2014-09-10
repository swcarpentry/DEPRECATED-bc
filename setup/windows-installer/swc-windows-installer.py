#!/usr/bin/env python

"""Software Carpentry Windows Installer

Helps mimic a *nix environment on Windows with as little work as possible.

The script:
* Installs nano and makes it accessible from msysgit
* Installs sqlite3 and makes it accessible from msysGit
* Creates ~/nano.rc with links to syntax highlighting configs
* Provides standard nosetests behavior for msysgit

To use:

1. Install Python, IPython, and Nose.  An easy way to do this is with
   the Anaconda CE Python distribution
   http://continuum.io/anacondace.html
2. Install msysgit
   http://code.google.com/p/msysgit/downloads/list
3. Run swc_windows_installer.py
   You should be able to simply double click the file in Windows

"""

import hashlib
try:  # Python 3
    from io import BytesIO as _BytesIO
except ImportError:  # Python 2
    from StringIO import StringIO as _BytesIO
import os
import re
import sys
import tarfile
try:  # Python 3
    from urllib.request import urlopen as _urlopen
except ImportError:  # Python 2
    from urllib2 import urlopen as _urlopen
import zipfile


if sys.version_info >= (3, 0):  # Python 3
    open3 = open
else:
    def open3(file, mode='r', newline=None):
        if newline:
            if newline != '\n':
                raise NotImplementedError(newline)
            f = open(file, mode + 'b')
        else:
            f = open(file, mode)
        return f


def download(url, sha1):
    """Download a file and verify it's hash"""
    r = _urlopen(url)
    byte_content = r.read()
    download_sha1 = hashlib.sha1(byte_content).hexdigest()
    if download_sha1 != sha1:
        raise ValueError(
            'downloaded {!r} has the wrong SHA1 hash: {} != {}'.format(
                url, download_sha1, sha1))
    return byte_content


def splitall(path):
    """Split a path into a list of components

    >>> splitall('nano-2.2.6/doc/Makefile.am')
    ['nano-2.2.6', 'doc', 'Makefile.am']
    """
    parts = []
    while True:
        head, tail = os.path.split(path)
        if tail:
            parts.insert(0, tail)
        elif head:
            parts.insert(0, head)
            break
        else:
            break
        path = head
    return parts


def transform(tarinfo, strip_components=0):
    """Transform TarInfo objects for extraction"""
    path_components = splitall(tarinfo.name)
    try:
        tarinfo.name = os.path.join(*path_components[strip_components:])
    except TypeError:
        if len(path_components) <= strip_components:
            return None
        raise
    return tarinfo


def tar_install(url, sha1, install_directory, compression='*',
                strip_components=0):
    """Download and install a tar bundle"""
    if not os.path.isdir(install_directory):
        tar_bytes = download(url=url, sha1=sha1)
        tar_io = _BytesIO(tar_bytes)
        filename = os.path.basename(url)
        mode = 'r:{}'.format(compression)
        tar_file = tarfile.open(filename, mode, tar_io)
        os.makedirs(install_directory)
        members = [
            transform(tarinfo=tarinfo, strip_components=strip_components)
            for tarinfo in tar_file]
        tar_file.extractall(
            path=install_directory,
            members=[m for m in members if m is not None])


def zip_install(url, sha1, install_directory):
    """Download and install a zipped bundle"""
    if not os.path.isdir(install_directory):
        zip_bytes = download(url=url, sha1=sha1)
        zip_io = _BytesIO(zip_bytes)
        zip_file = zipfile.ZipFile(zip_io)
        os.makedirs(install_directory)
        zip_file.extractall(install_directory)


def install_nano(install_directory):
    """Download and install the nano text editor"""
    zip_install(
        url='http://www.nano-editor.org/dist/v2.2/NT/nano-2.2.6.zip',
        sha1='f5348208158157060de0a4df339401f36250fe5b',
        install_directory=install_directory)


def install_nanorc(install_directory):
    """Download and install nano syntax highlighting"""
    tar_install(
        url='http://www.nano-editor.org/dist/v2.2/nano-2.2.6.tar.gz',
        sha1='f2a628394f8dda1b9f28c7e7b89ccb9a6dbd302a',
        install_directory=install_directory,
        strip_components=1)
    home = os.path.expanduser('~')
    nanorc = os.path.join(home, 'nano.rc')
    if not os.path.isfile(nanorc):
        syntax_dir = os.path.join(install_directory, 'doc', 'syntax')
        with open3(nanorc, 'w', newline='\n') as f:
            for filename in os.listdir(syntax_dir):
                if filename.endswith('.nanorc'):
                    path = os.path.join(syntax_dir, filename)
                    rel_path = os.path.relpath(path, home)
                    include_path = make_posix_path(os.path.join('~', rel_path))
                    f.write('include {}\n'.format(include_path))


def install_sqlite(install_directory):
    """Download and install the sqlite3 shell"""
    zip_install(
        url='https://sqlite.org/2014/sqlite-shell-win32-x86-3080403.zip',
        sha1='1a8ab0ca9f4c51afeffeb49bd301e1d7f64741bb',
        install_directory=install_directory)


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
    if not os.path.isdir(python_scripts_directory):
        os.makedirs(python_scripts_directory)
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
        'export PATH=\"$PATH:{}\"'.format(':'.join(
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
    for regex, sub in [
            (re.compile(r'\\'), '/'),
            (re.compile('^[Cc]:'), '/c'),
            ]:
        windows_path = regex.sub(sub, windows_path)
    return windows_path


def main():
    swc_dir = os.path.join(os.path.expanduser('~'), '.swc')
    bin_dir = os.path.join(swc_dir, 'bin')
    nano_dir = os.path.join(swc_dir, 'lib', 'nano')
    nanorc_dir = os.path.join(swc_dir, 'share', 'nanorc')
    sqlite_dir = os.path.join(swc_dir, 'lib', 'sqlite')
    create_nosetests_entry_point(python_scripts_directory=bin_dir)
    install_nano(install_directory=nano_dir)
    install_nanorc(install_directory=nanorc_dir)
    install_sqlite(install_directory=sqlite_dir)
    update_bash_profile(extra_paths=(nano_dir, sqlite_dir, bin_dir))


if __name__ == '__main__':
    print("Preparing your Software Carpentry awesomeness!")
    main()
