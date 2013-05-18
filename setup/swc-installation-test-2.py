#!/usr/bin/env python

"""Test script to check for required functionality.

Execute this code at the command line by typing:

  python swc-installation-test-2.py

Run the script and follow the instructions it prints at the end.

This script requires at least Python 2.6.  You can check the version
of Python that you have installed with 'swc-installation-test-1.py'.

By default, this script will test for all the dependencies your
instructor thinks you need.  If you want to test for a different set
of packages, you can list them on the command line.  For example:

  python swc-installation-test-2.py git virtual-editor

This is useful if the original test told you to install a more recent
version of a particular dependency, and you just want to re-test that
dependency.
"""

from __future__ import print_function  # for Python 2.6 compatibility

import distutils.ccompiler as _distutils_ccompiler
import fnmatch as _fnmatch
try:  # Python 2.7 and 3.x
    import importlib as _importlib
except ImportError:  # Python 2.6 and earlier
    class _Importlib (object):
        """Minimal workarounds for functions we need
        """
        @staticmethod
        def import_module(name):
            module = __import__(name)
            for n in name.split('.')[1:]:
                module = getattr(module, n)
            return module
    _importlib = _Importlib()
import logging as _logging
import os as _os
import platform as _platform
import re as _re
import shlex as _shlex
import subprocess as _subprocess
import sys as _sys
try:  # Python 3.x
    import urllib.parse as _urllib_parse
except ImportError:  # Python 2.x
    import urllib as _urllib_parse  # for quote()


if not hasattr(_shlex, 'quote'):  # Python versions older than 3.3
    # Use the undocumented pipes.quote()
    import pipes as _pipes
    _shlex.quote = _pipes.quote


__version__ = '0.1'

# Comment out any entries you don't need
CHECKS = [
# Shell
    'virtual-shell',
# Editors
    'virtual-editor',
# Browsers
    'virtual-browser',
# Version control
    'git',
    'hg',              # Command line tool
    #'mercurial',       # Python package
    'EasyMercurial',
# Build tools and packaging
    'make',
    'virtual-pypi-installer',
    'setuptools',
    #'xcode',
# Testing
    'nosetests',       # Command line tool
    'nose',            # Python package
# SQL
    'sqlite3',         # Command line tool
    'sqlite3-python',  # Python package
# Python
    'python',
    'ipython',         # Command line tool
    'IPython',         # Python package
    'argparse',        # Useful for utility scripts
    'numpy',
    'scipy',
    'matplotlib',
    'pandas',
    'sympy',
    'Cython',
    'networkx',
    'mayavi.mlab',
    ]

CHECKER = {}

_ROOT_PATH = _os.sep
if _platform.system() == 'win32':
    _ROOT_PATH = 'c:\\'


class InvalidCheck (KeyError):
    def __init__(self, check):
        super(InvalidCheck, self).__init__(check)
        self.check = check

    def __str__(self):
        return self.check


class DependencyError (Exception):
    _default_url = 'http://software-carpentry.org/setup/'
    _setup_urls = {  # (system, version, package) glob pairs
        ('*', '*', 'Cython'): 'http://docs.cython.org/src/quickstart/install.html',
        ('Linux', '*', 'EasyMercurial'): 'http://easyhg.org/download.html#download-linux',
        ('Darwin', '*', 'EasyMercurial'): 'http://easyhg.org/download.html#download-mac',
        ('Windows', '*', 'EasyMercurial'): 'http://easyhg.org/download.html#download-windows',
        ('*', '*', 'EasyMercurial'): 'http://easyhg.org/download.html',
        ('*', '*', 'argparse'): 'https://pypi.python.org/pypi/argparse#installation',
        ('*', '*', 'ash'): 'http://www.in-ulm.de/~mascheck/various/ash/',
        ('*', '*', 'bash'): 'http://www.gnu.org/software/bash/manual/html_node/Basic-Installation.html#Basic-Installation',
        ('Linux', '*', 'chromium'): 'http://code.google.com/p/chromium/wiki/LinuxBuildInstructions',
        ('Darwin', '*', 'chromium'): 'http://code.google.com/p/chromium/wiki/MacBuildInstructions',
        ('Windows', '*', 'chromium'): 'http://www.chromium.org/developers/how-tos/build-instructions-windows',
        ('*', '*', 'chromium'): 'http://www.chromium.org/developers/how-tos',
        ('Windows', '*', 'emacs'): 'http://www.gnu.org/software/emacs/windows/Installing-Emacs.html',
        ('*', '*', 'emacs'): 'http://www.gnu.org/software/emacs/#Obtaining',
        ('*', '*', 'firefox'): 'http://www.mozilla.org/en-US/firefox/new/',
        ('Linux', '*', 'gedit'): 'http://www.linuxfromscratch.org/blfs/view/svn/gnome/gedit.html',
        ('*', '*', 'git'): 'http://git-scm.com/downloads',
        ('*', '*', 'google-chrome'): 'https://www.google.com/intl/en/chrome/browser/',
        ('*', '*', 'hg'): 'http://mercurial.selenic.com/',
        ('*', '*', 'mercurial'): 'http://mercurial.selenic.com/',
        ('*', '*', 'IPython'): 'http://ipython.org/install.html',
        ('*', '*', 'ipython'): 'http://ipython.org/install.html',
        ('*', '*', 'jinja'): 'http://jinja.pocoo.org/docs/intro/#installation',
        ('*', '*', 'kate'): 'http://kate-editor.org/get-it/',
        ('*', '*', 'make'): 'http://www.gnu.org/software/make/',
        ('Darwin', '*', 'matplotlib'): 'http://matplotlib.org/users/installing.html#building-on-osx',
        ('Windows', '*', 'matplotlib'): 'http://matplotlib.org/users/installing.html#installing-on-windows',
        ('*', '*', 'matplotlib'): 'http://matplotlib.org/users/installing.html#installing',
        ('*', '*', 'mayavi.mlab'): 'http://docs.enthought.com/mayavi/mayavi/installation.html',
        ('*', '*', 'nano'): 'http://www.nano-editor.org/dist/latest/faq.html#3',
        ('*', '*', 'networkx'): 'http://networkx.github.com/documentation/latest/install.html#installing',
        ('*', '*', 'nose'): 'https://nose.readthedocs.org/en/latest/#installation-and-quick-start',
        ('*', '*', 'nosetests'): 'https://nose.readthedocs.org/en/latest/#installation-and-quick-start',
        ('*', '*', 'notepad++'): 'http://notepad-plus-plus.org/download/v6.3.html',
        ('*', '*', 'numpy'): 'http://docs.scipy.org/doc/numpy/user/install.html',
        ('*', '*', 'pandas'): 'http://pandas.pydata.org/pandas-docs/stable/install.html',
        ('*', '*', 'pip'): 'http://www.pip-installer.org/en/latest/installing.html',
        ('*', '*', 'python'): 'http://www.python.org/download/releases/2.7.3/#download',
        ('*', '*', 'pyzmq'): 'https://github.com/zeromq/pyzmq/wiki/Building-and-Installing-PyZMQ',
        ('Linux', '*', 'scipy'): 'http://www.scipy.org/Installing_SciPy/Linux',
        ('Darwin', '*', 'scipy'): 'http://www.scipy.org/Installing_SciPy/Mac_OS_X',
        ('Windows', '*', 'scipy'): 'http://www.scipy.org/Installing_SciPy/Windows',
        ('*', '*', 'scipy'): 'http://www.scipy.org/Installing_SciPy',
        ('*', '*', 'setuptools'): 'https://pypi.python.org/pypi/setuptools#installation-instructions',
        ('*', '*', 'sqlite3'): 'http://www.sqlite.org/download.html',
        ('*', '*', 'sublime-text'): 'http://www.sublimetext.com/2',
        ('*', '*', 'sympy'): 'http://docs.sympy.org/dev/install.html',
        ('Darwin', '*', 'textmate'): 'http://macromates.com/',
        ('Darwin', '*', 'textwrangler'): 'http://www.barebones.com/products/textwrangler/download.html',
        ('*', '*', 'tornado'): 'http://www.tornadoweb.org/',
        ('*', '*', 'vim'): 'http://www.vim.org/download.php',
        ('Darwin', '*', 'xcode'): 'https://developer.apple.com/xcode/',
        ('*', '*', 'xemacs'): 'http://www.us.xemacs.org/Install/',
        ('*', '*', 'zsh'): 'http://www.zsh.org/',
        }

    def _get_message(self):
        return self._message
    def _set_message(self, message):
        self._message = message
    message = property(_get_message, _set_message)

    def __init__(self, checker, message, causes=None):
        super(DependencyError, self).__init__(message)
        self.checker = checker
        self.message = message
        if causes is None:
            causes = []
        self.causes = causes

    def get_url(self):
        system = _platform.system()
        version = None
        for pversion in (
            'linux_distribution',
            'mac_ver',
            'win32_ver',
            ):
            value = getattr(_platform, pversion)()
            if value[0]:
                version = value[0]
                break
        package = self.checker.name
        for (s,v,p),url in self._setup_urls.items():
            if (_fnmatch.fnmatch(system, s) and
                    _fnmatch.fnmatch(version, v) and
                    _fnmatch.fnmatch(package, p)):
                return url
        return self._default_url

    def __str__(self):
        url = self.get_url()
        lines = [
            'check for {0} failed:'.format(self.checker.full_name()),
            '  ' + self.message,
            '  For instructions on installing an up-to-date version, see',
            '  ' + url,
            ]
        if self.causes:
            lines.append('  causes:')
            for cause in self.causes:
                lines.extend('  ' + line for line in str(cause).splitlines())
        return '\n'.join(lines)


def check(checks=None):
    successes = []
    failures = []
    if not checks:
        checks = CHECKS
    for check in checks:
        try:
            checker = CHECKER[check]
        except KeyError as e:
            raise InvalidCheck(check)# from e
        _sys.stdout.write('check {0}...\t'.format(checker.full_name()))
        try:
            version = checker.check()
        except DependencyError as e:
            failures.append(e)
            _sys.stdout.write('fail\n')
        else:
            _sys.stdout.write('pass\n')
            successes.append((checker, version))
    if successes:
        print('\nSuccesses:\n')
        for checker,version in successes:
            print('{0} {1}'.format(
                    checker.full_name(),
                    version or 'unknown'))
    if failures:
        print('\nFailures:')
        printed = []
        for failure in failures:
            if failure not in printed:
                print()
                print(failure)
                printed.append(failure)
        return False
    return True


class Dependency (object):
    def __init__(self, name, long_name=None, minimum_version=None,
                 version_delimiter='.', and_dependencies=None,
                 or_dependencies=None):
        self.name = name
        self.long_name = long_name or name
        self.minimum_version = minimum_version
        self.version_delimiter = version_delimiter
        if not and_dependencies:
            and_dependencies = []
        self.and_dependencies = and_dependencies
        if not or_dependencies:
            or_dependencies = []
        self.or_dependencies = or_dependencies
        self._check_error = None

    def __str__(self):
        return '<{0} {1}>'.format(type(self).__name__, self.name)

    def full_name(self):
        if self.name == self.long_name:
            return self.name
        else:
            return '{0} ({1})'.format(self.long_name, self.name)

    def check(self):
        if self._check_error:
            raise self._check_error
        try:
            self._check_dependencies()
            return self._check()
        except DependencyError as e:
            self._check_error = e  # cache for future calls
            raise

    def _check_dependencies(self):
        for dependency in self.and_dependencies:
            if not hasattr(dependency, 'check'):
                dependency = CHECKER[dependency]
            try:
                dependency.check()
            except DependencyError as e:
                raise DependencyError(
                    checker=self,
                    message=(
                        'some dependencies for {0} were not satisfied'
                        ).format(self.full_name()),
                    causes=[e])
        self.or_pass = None
        or_errors = []
        for dependency in self.or_dependencies:
            if not hasattr(dependency, 'check'):
                dependency = CHECKER[dependency]
            try:
                version = dependency.check()
            except DependencyError as e:
                or_errors.append(e)
            else:
                self.or_pass = {
                    'dependency': dependency,
                    'version': version,
                    }
                break  # no need to test other dependencies
        if self.or_dependencies and not self.or_pass:
            raise DependencyError(
                checker=self,
                message=(
                    '{0} requires at least one of the following dependencies'
                    ).format(self.full_name()),
                    causes=or_errors)

    def _check(self):
        version = self._get_version()
        parsed_version = None
        if hasattr(self, '_get_parsed_version'):
            parsed_version = self._get_parsed_version()
        if self.minimum_version:
            self._check_version(version=version, parsed_version=parsed_version)
        return version

    def _get_version(self):
        raise NotImplementedError(self)

    def _minimum_version_string(self):
        return self.version_delimiter.join(
            str(part) for part in self.minimum_version)

    def _check_version(self, version, parsed_version=None):
        if not parsed_version:
            parsed_version = self._parse_version(version=version)
        if not parsed_version or parsed_version < self.minimum_version:
            raise DependencyError(
                checker=self,
                message='outdated version of {0}: {1} (need >= {2})'.format(
                    self.full_name(), version, self._minimum_version_string()))

    def _parse_version(self, version):
        if not version:
            return None
        parsed_version = []
        for part in version.split(self.version_delimiter):
            try:
                parsed_version.append(int(part))
            except ValueError as e:
                raise DependencyError(
                    checker=self,
                    message=(
                        'unparsable {0!r} in version {1} of {2}, (need >= {3})'
                        ).format(
                        part, version, self.full_name(),
                        self._minimum_version_string()))# from e
        return tuple(parsed_version)


class PythonDependency (Dependency):
    def __init__(self, name='python', long_name='Python version',
                 minimum_version=(2, 6), **kwargs):
        super(PythonDependency, self).__init__(
            name=name, long_name=long_name, minimum_version=minimum_version,
            **kwargs)

    def _get_version(self):
        return _sys.version

    def _get_parsed_version(self):
        return _sys.version_info


CHECKER['python'] = PythonDependency()


class CommandDependency (Dependency):
    exe_extension = _distutils_ccompiler.new_compiler().exe_extension

    def __init__(self, command, paths=None, version_options=('--version',),
                 stdin=None, version_regexp=None, version_stream='stdout',
                 **kwargs):
        if 'name' not in kwargs:
            kwargs['name'] = command
        super(CommandDependency, self).__init__(**kwargs)
        self.command = command
        self.paths = paths
        self.version_options = version_options
        self.stdin = None
        if not version_regexp:
            regexp = r'([\d][\d{0}]*[\d])'.format(self.version_delimiter)
            version_regexp = _re.compile(regexp)
        self.version_regexp = version_regexp
        self.version_stream = version_stream

    def _get_command_version_stream(self, command=None, stdin=None,
                                    expect=(0,)):
        if command is None:
            command = self.command + (self.exe_extension or '')
        if not stdin:
            stdin = self.stdin
        if stdin:
            popen_stdin = _subprocess.PIPE
        else:
            popen_stdin = None
        try:
            p = _subprocess.Popen(
                [command] + list(self.version_options), stdin=popen_stdin,
                stdout=_subprocess.PIPE, stderr=_subprocess.PIPE,
                universal_newlines=True)
        except OSError as e:
            raise DependencyError(
                checker=self,
                message="could not find '{0}' executable".format(command),
                )# from e
        stdout,stderr = p.communicate(stdin)
        status = p.wait()
        if status not in expect:
            lines = [
                "failed to execute: {0} {1}".format(
                    command,
                    ' '.join(_shlex.quote(arg)
                             for arg in self.version_options)),
                'status: {0}'.format(status),
                ]
            for name,string in [('stdout', stdout), ('stderr', stderr)]:
                if string:
                    lines.extend([name + ':', string])
            raise DependencyError(checker=self, message='\n'.join(lines))
        for name,string in [('stdout', stdout), ('stderr', stderr)]:
            if name == self.version_stream:
                if not string:
                    raise DependencyError(
                        checker=self,
                        message='empty version stream on {0} for {1}'.format(
                            self.version_stream, command))
                return string
        raise NotImplementedError(self.version_stream)

    def _get_version_stream(self, **kwargs):
        paths = [self.command + (self.exe_extension or '')]
        if self.exe_extension:
            paths.append(self.command)  # also look at the extension-less path
        if self.paths:
            paths.extend(self.paths)
        or_errors = []
        for path in paths:
            try:
                return self._get_command_version_stream(command=path, **kwargs)
            except DependencyError as e:
                or_errors.append(e)
        raise DependencyError(
            checker=self,
            message='errors finding {0} version'.format(
                self.full_name()),
            causes=or_errors)

    def _get_version(self):
        version_stream = self._get_version_stream()
        match = self.version_regexp.search(version_stream)
        if not match:
            raise DependencyError(
                checker=self,
                message='no version string in output:\n{0}'.format(
                    version_stream))
        return match.group(1)


def _program_files_paths(*args):
    "Utility for generating MS Windows search paths"
    pf = _os.environ.get('ProgramFiles', '/usr/bin')
    pfx86 = _os.environ.get('ProgramFiles(x86)', pf)
    paths = [_os.path.join(pf, *args)]
    if pfx86 != pf:
        paths.append(_os.path.join(pfx86, *args))
    return paths


for command,long_name,minimum_version,paths in [
        ('sh', 'Bourne Shell', None, None),
        ('ash', 'Almquist Shell', None, None),
        ('bash', 'Bourne Again Shell', None, None),
        ('csh', 'C Shell', None, None),
        ('ksh', 'KornShell', None, None),
        ('dash', 'Debian Almquist Shell', None, None),
        ('tcsh', 'TENEX C Shell', None, None),
        ('zsh', 'Z Shell', None, None),
        ('git', 'Git', (1, 7, 0), None),
        ('hg', 'Mercurial', (2, 0, 0), None),
        ('EasyMercurial', None, (1, 3), None),
        ('pip', None, None, None),
        ('sqlite3', 'SQLite 3', None, None),
        ('nosetests', 'Nose', (1, 0, 0), None),
        ('ipython', 'IPython script', (0, 13), None),
        ('emacs', 'Emacs', None, None),
        ('xemacs', 'XEmacs', None, None),
        ('vim', 'Vim', None, None),
        ('vi', None, None, None),
        ('nano', 'Nano', None, None),
        ('gedit', None, None, None),
        ('kate', 'Kate', None, None),
        ('notepad++', 'Notepad++', None,
         _program_files_paths('Notepad++', 'notepad++.exe')),
        ('firefox', 'Firefox', None,
         _program_files_paths('Mozilla Firefox', 'firefox.exe')),
        ('google-chrome', 'Google Chrome', None,
         _program_files_paths('Google', 'Chrome', 'Application', 'chrome.exe')
         ),
        ('chromium', 'Chromium', None, None),
        ]:
    if not long_name:
        long_name = command
    CHECKER[command] = CommandDependency(
        command=command, paths=paths, long_name=long_name,
        minimum_version=minimum_version)
del command, long_name, minimum_version, paths  # cleanup namespace


class MakeDependency (CommandDependency):
    makefile = '\n'.join([
            'all:',
            '\t@echo "MAKE_VERSION=$(MAKE_VERSION)"',
            '\t@echo "MAKE=$(MAKE)"',
            '',
            ])

    def _get_version(self):
        try:
            return super(MakeDependency, self)._get_version()
        except DependencyError as e:
            version_options = self.version_options
            self.version_options = ['-f', '-']
            try:
                stream = self._get_version_stream(stdin=self.makefile)
                info = {}
                for line in stream.splitlines():
                    try:
                        key,value = line.split('=', 1)
                    except ValueError as ve:
                        raise e# from NotImplementedError(stream)
                    info[key] = value
                if info.get('MAKE_VERSION', None):
                    return info['MAKE_VERSION']
                elif info.get('MAKE', None):
                    return None
                raise e
            finally:
                self.version_options = version_options


CHECKER['make'] = MakeDependency(command='make', minimum_version=None)


class EasyInstallDependency (CommandDependency):
    def _get_version(self):
        try:
            return super(EasyInstallDependency, self)._get_version()
        except DependencyError as e:
            version_stream = self.version_stream
            try:
                self.version_stream = 'stderr'
                stream = self._get_version_stream(expect=(1,))
                if 'option --version not recognized' in stream:
                    return 'unknown (possibly Setuptools?)'
            finally:
                self.version_stream = version_stream


CHECKER['easy_install'] = EasyInstallDependency(
    command='easy_install', long_name='Setuptools easy_install',
    minimum_version=None)


class PathCommandDependency (CommandDependency):
    """A command that doesn't support --version or equivalent options

    On some operating systems (e.g. OS X), a command's executable may
    be hard to find, or not exist in the PATH.  Work around that by
    just checking for the existence of a characteristic file or
    directory.  Since the characteristic path may depend on OS,
    installed version, etc., take a list of paths, and succeed if any
    of them exists.
    """
    def _get_command_version_stream(self, *args, **kwargs):
        raise NotImplementedError()

    def _get_version_stream(self, *args, **kwargs):
        raise NotImplementedError()

    def _get_version(self):
        for path in self.paths:
            if _os.path.exists(path):
                return None
        raise DependencyError(
            checker=self,
            message=(
                'nothing exists at any of the expected paths for {0}:\n    {1}'
                ).format(
                self.full_name(),
                '\n    '.join(p for p in self.paths)))


for paths,name,long_name in [
        ([_os.path.join(_ROOT_PATH, 'Applications', 'Sublime Text 2.app')],
         'sublime-text', 'Sublime Text'),
        ([_os.path.join(_ROOT_PATH, 'Applications', 'TextMate.app')],
         'textmate', 'TextMate'),
        ([_os.path.join(_ROOT_PATH, 'Applications', 'TextWrangler.app')],
         'textwrangler', 'TextWrangler'),
        ([_os.path.join(_ROOT_PATH, 'Applications', 'Safari.app')],
         'safari', 'Safari'),
        ([_os.path.join(_ROOT_PATH, 'Applications', 'Xcode.app'),  # OS X >=1.7
          _os.path.join(_ROOT_PATH, 'Developer', 'Applications', 'Xcode.app'
                        )  # OS X 1.6,
          ],
         'xcode', 'Xcode'),
        ]:
    if not long_name:
        long_name = name
    CHECKER[name] = PathCommandDependency(
        command=None, paths=paths, name=name, long_name=long_name)
del paths, name, long_name  # cleanup namespace


class PythonPackageDependency (Dependency):
    def __init__(self, package, **kwargs):
        if 'name' not in kwargs:
            kwargs['name'] = package
        if 'and_dependencies' not in kwargs:
            kwargs['and_dependencies'] = []
        if 'python' not in kwargs['and_dependencies']:
            kwargs['and_dependencies'].append('python')
        super(PythonPackageDependency, self).__init__(**kwargs)
        self.package = package

    def _get_version(self):
        package = self._get_package(self.package)
        return self._get_version_from_package(package)

    def _get_package(self, package):
        try:
            return _importlib.import_module(package)
        except ImportError as e:
            raise DependencyError(
                checker=self,
                message="could not import the '{0}' package for {1}".format(
                    package, self.full_name()),
                )# from e

    def _get_version_from_package(self, package):
        try:
            version = package.__version__
        except AttributeError:
            version = None
        return version


for package,name,long_name,minimum_version,and_dependencies in [
        ('nose', None, 'Nose Python package',
         CHECKER['nosetests'].minimum_version, None),
        ('jinja2', 'jinja', 'Jinja', (2, 6), None),
        ('zmq', 'pyzmq', 'PyZMQ', (2, 1, 4), None),
        ('IPython', None, 'IPython Python package',
         CHECKER['ipython'].minimum_version, ['jinja', 'tornado', 'pyzmq']),
        ('argparse', None, 'Argparse', None, None),
        ('numpy', None, 'NumPy', None, None),
        ('scipy', None, 'SciPy', None, None),
        ('matplotlib', None, 'Matplotlib', None, None),
        ('pandas', None, 'Pandas', (0, 8), None),
        ('sympy', None, 'SymPy', None, None),
        ('Cython', None, None, None, None),
        ('networkx', None, 'NetworkX', None, None),
        ('mayavi.mlab', None, 'MayaVi', None, None),
        ('setuptools', None, 'Setuptools', None, None),
        ]:
    if not name:
        name = package
    if not long_name:
        long_name = name
    kwargs = {}
    if and_dependencies:
        kwargs['and_dependencies'] = and_dependencies
    CHECKER[name] = PythonPackageDependency(
        package=package, name=name, long_name=long_name,
        minimum_version=minimum_version, **kwargs)
# cleanup namespace
del package, name, long_name, minimum_version, and_dependencies, kwargs


class MercurialPythonPackage (PythonPackageDependency):
    def _get_version(self):
        try:  # mercurial >= 1.2
            package = _importlib.import_module('mercurial.util')
        except ImportError as e:  # mercurial <= 1.1.2
            package = self._get_package('mercurial.version')
            return package.get_version()
        else:
            return package.version()


CHECKER['mercurial'] = MercurialPythonPackage(
    package='mercurial.util', name='mercurial',
    long_name='Mercurial Python package',
    minimum_version=CHECKER['hg'].minimum_version)


class TornadoPythonPackage (PythonPackageDependency):
    def _get_version_from_package(self, package):
        return package.version

    def _get_parsed_version(self):
        package = self._get_package(self.package)
        return package.version_info


CHECKER['tornado'] = TornadoPythonPackage(
    package='tornado', name='tornado', long_name='Tornado', minimum_version=(2, 0))


class SQLitePythonPackage (PythonPackageDependency):
    def _get_version_from_package(self, package):
        return _sys.version

    def _get_parsed_version(self):
        return _sys.version_info


CHECKER['sqlite3-python'] = SQLitePythonPackage(
    package='sqlite3', name='sqlite3-python',
    long_name='SQLite Python package',
    minimum_version=CHECKER['sqlite3'].minimum_version)


class UserTaskDependency (Dependency):
    "Prompt the user to complete a task and check for success"
    def __init__(self, prompt, **kwargs):
        super(UserTaskDependency, self).__init__(**kwargs)
        self.prompt = prompt

    def _check(self):
        if _sys.version_info >= (3, ):
            result = input(self.prompt)
        else:  # Python 2.x
            result = raw_input(self.prompt)
        return self._check_result(result)

    def _check_result(self, result):
        raise NotImplementedError()


class EditorTaskDependency (UserTaskDependency):
    def __init__(self, **kwargs):
        self.path = _os.path.expanduser(_os.path.join(
                '~', 'swc-installation-test.txt'))
        self.contents = 'Hello, world!'
        super(EditorTaskDependency, self).__init__(
            prompt=(
                'Open your favorite text editor and create the file\n'
                '  {0}\n'
                'containing the line:\n'
                '  {1}\n'
                'Press enter here after you have done this.\n'
                'You may remove the file after you have finished testing.'
                ).format(self.path, self.contents),
            **kwargs)

    def _check_result(self, result):
        message = None
        try:
            with open(self.path, 'r') as f:
                contents = f.read()
        except IOError as e:
            raise DependencyError(
                checker=self,
                message='could not open {0!r}: {1}'.format(self.path, e)
                )# from e
        if contents.strip() != self.contents:
            raise DependencyError(
                checker=self,
                message=(
                    'file contents ({0!r}) did not match the expected {1!r}'
                    ).format(contents, self.contents))


CHECKER['other-editor'] = EditorTaskDependency(
    name='other-editor', long_name='')


class VirtualDependency (Dependency):
    def _check(self):
        return '{0} {1}'.format(
            self.or_pass['dependency'].full_name(),
            self.or_pass['version'])


for name,long_name,dependencies in [
        ('virtual-shell', 'command line shell', (
            'bash',
            'dash',
            'ash',
            'zsh',
            'ksh',
            'csh',
            'tcsh',
            'sh',
            )),
        ('virtual-editor', 'text/code editor', (
            'emacs',
            'xemacs',
            'vim',
            'vi',
            'nano',
            'gedit',
            'kate',
            'notepad++',
            'sublime-text',
            'textmate',
            'textwrangler',
            'other-editor',  # last because it requires user interaction
            )),
        ('virtual-browser', 'web browser', (
            'firefox',
            'google-chrome',
            'chromium',
            'safari',
            )),
        ('virtual-pypi-installer', 'PyPI installer', (
            'easy_install',
            'pip',
            )),
        ]:
    CHECKER[name] = VirtualDependency(
        name=name, long_name=long_name, or_dependencies=dependencies)
del name, long_name, dependencies  # cleanup namespace


def _print_info(key, value, indent=19):
    print('{0}{1}: {2}'.format(key, ' '*(indent-len(key)), value))

def print_system_info():
    print("If you do not understand why the above failures occurred,")
    print("copy and send the *entire* output (all info above and summary")
    print("below) to the instructor for help.")
    print()
    print('==================')
    print('System information')
    print('==================')
    _print_info('os.name', _os.name)
    _print_info('os.uname', _platform.uname())
    _print_info('platform', _sys.platform)
    _print_info('platform+', _platform.platform())
    for pversion in (
            'linux_distribution',
            'mac_ver',
            'win32_ver',
            ):
        value = getattr(_platform, pversion)()
        if value[0]:
            _print_info(pversion, value)
    _print_info('prefix', _sys.prefix)
    _print_info('exec_prefix', _sys.exec_prefix)
    _print_info('executable', _sys.executable)
    _print_info('version_info', _sys.version_info)
    _print_info('version', _sys.version)
    _print_info('environment', '')
    for key,value in sorted(_os.environ.items()):
        print('  {0}={1}'.format(key, value))
    print('==================')

def print_suggestions(instructor_fallback=True):
    print()
    print('For suggestions on installing missing packages, see')
    print('http://software-carpentry.org/setup/')
    print('')
    print('For instructings on installing a particular package,')
    print('see the failure message for that package printed above.')
    if instructor_fallback:
        print('')
        print('For help, email the *entire* output of this script to')
        print('your instructor.')


if __name__ == '__main__':
    import optparse as _optparse

    parser = _optparse.OptionParser(usage='%prog [options] [check...]')
    epilog = __doc__
    parser.format_epilog = lambda formatter: '\n' + epilog
    parser.add_option(
        '-v', '--verbose', action='store_true',
        help=('print additional information to help troubleshoot '
              'installation issues'))
    options,args = parser.parse_args()
    try:
        passed = check(args)
    except InvalidCheck as e:
        print("I don't know how to check for {0!r}".format(e.check))
        print('I do know how to check for:')
        for key,checker in sorted(CHECKER.items()):
            if checker.long_name != checker.name:
                print('  {0} {1}({2})'.format(
                        key, ' '*(20-len(key)), checker.long_name))
            else:
                print('  {0}'.format(key))
        _sys.exit(1)
    if not passed:
        if options.verbose:
            print()
            print_system_info()
            print_suggestions(instructor_fallback=True)
        _sys.exit(1)
