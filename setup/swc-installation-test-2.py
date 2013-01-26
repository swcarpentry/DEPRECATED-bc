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
# Build tools and packaging
    'make',
    'virtual-pypi-installer',
    'setuptools',
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
    'numpy',
    'scipy',
    'matplotlib',
    'sympy',
    'Cython',
    'networkx',
    'mayavi.mlab',
    ]

CHECKER = {}


class InvalidCheck (KeyError):
    def __init__(self, check):
        super(InvalidCheck, self).__init__(check)
        self.check = check

    def __str__(self):
        return self.check


class DependencyError (Exception):
    def _get_message(self):
        return self._message
    def _set_message(self, message):
        self._message = message
    message = property(_get_message, _set_message)

    def __init__(self, checker, message):
        super(DependencyError, self).__init__(message)
        self.checker = checker
        self.message = message

    def __str__(self):
        url = 'http://software-carpentry.org/setup/'  # TODO: per-package URL
        return 'check for {0} failed:\n{1}\n{2}\n{3}'.format(
            self.checker.full_name(), self.message,
            'For instructions on installing an up-to-date version, see',
            url)


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
            dependency.check()
        self.or_pass = or_error = None
        for dependency in self.or_dependencies:
            if not hasattr(dependency, 'check'):
                dependency = CHECKER[dependency]
            try:
                version = dependency.check()
            except DependencyError as e:
                or_error = e
            else:
                self.or_pass = {
                    'dependency': dependency,
                    'version': version,
                    }
                break  # no need to test other dependencies
        if self.or_dependencies and not self.or_pass:
            raise or_error

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

    def __init__(self, command, version_options=('--version',), stdin=None,
                 version_regexp=None, version_stream='stdout', **kwargs):
        if 'name' not in kwargs:
            kwargs['name'] = command
        super(CommandDependency, self).__init__(**kwargs)
        self.command = command
        self.version_options = version_options
        self.stdin = None
        if not version_regexp:
            regexp = r'([\d][\d{0}]*[\d])'.format(self.version_delimiter)
            version_regexp = _re.compile(regexp)
        self.version_regexp = version_regexp
        self.version_stream = version_stream

    def _get_version_stream(self, stdin=None, expect=(0,)):
        if not stdin:
            stdin = self.stdin
        if stdin:
            popen_stdin = _subprocess.PIPE
        else:
            popen_stdin = None
        command = self.command + (self.exe_extension or '')
        try:
            p = _subprocess.Popen(
                [command] + list(self.version_options), stdin=popen_stdin,
                stdout=_subprocess.PIPE, stderr=_subprocess.PIPE,
                close_fds=True, shell=False, universal_newlines=True)
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
                return string
        raise NotImplementedError(self.version_stream)

    def _get_version(self):
        version_stream = self._get_version_stream()
        match = self.version_regexp.search(version_stream)
        if not match:
            raise DependencyError(
                checker=self,
                message='no version string in output:\n{0}'.format(
                    version_stream))
        return match.group(1)


for command,long_name,minimum_version in [
        ('sh', 'Bourne Shell', None),
        ('ash', 'Almquist Shell', None),
        ('bash', 'Bourne Again Shell', None),
        ('csh', 'C Shell', None),
        ('ksh', 'KornShell', None),
        ('dash', 'Debian Almquist Shell', None),
        ('tcsh', 'TENEX C Shell', None),
        ('zsh', 'Z Shell', None),
        ('git', 'Git', (1, 7, 0)),
        ('hg', 'Mercurial', (2, 0, 0)),
        ('pip', None, None),
        ('sqlite3', 'SQLite 3', None),
        ('nosetests', 'Nose', (1, 0, 0)),
        ('ipython', 'IPython script', (0, 13)),
        ('emacs', 'Emacs', None),
        ('xemacs', 'XEmacs', None),
        ('vim', 'Vim', None),
        ('vi', None, None),
        ('nano', 'Nano', None),
        ('kate', 'Kate', None),
        ('notepad++', 'Notepad++', None),
        ('firefox', 'Firefox', None),
        ('google-chrome', 'Google Chrome', None),
        ('chromium', 'Chromium', None),
        ]:
    if not long_name:
        long_name = command
    CHECKER[command] = CommandDependency(
        command=command, long_name=long_name, minimum_version=minimum_version)
del command, long_name, minimum_version  # cleanup namespace


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


for package,name,long_name,minimum_version in [
        ('nose', None, 'Nose Python package',
         CHECKER['nosetests'].minimum_version),
        ('IPython', None, 'IPython Python package',
         CHECKER['ipython'].minimum_version),
        ('numpy', None, 'NumPy', None),
        ('scipy', None, 'SciPy', None),
        ('matplotlib', None, 'Matplotlib', None),
        ('sympy', None, 'SymPy', None),
        ('Cython', None, None, None),
        ('networkx', None, 'NetworkX', None),
        ('mayavi.mlab', None, 'MayaVi', None),
        ('setuptools', None, 'Setuptools', None),
        ]:
    if not name:
        name = package
    if not long_name:
        long_name = name
    CHECKER[name] = PythonPackageDependency(
        package=package, name=name, long_name=long_name,
        minimum_version=minimum_version)
del package, name, long_name, minimum_version  # cleanup namespace


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


class SQLitePythonPackage (PythonPackageDependency):
    def _get_version_from_package(self, package):
        return _sys.version

    def _get_parsed_version(self):
        return _sys.version_info


CHECKER['sqlite3-python'] = SQLitePythonPackage(
    package='sqlite3', name='sqlite3-python',
    long_name='SQLite Python package',
    minimum_version=CHECKER['sqlite3'].minimum_version)


class VirtualDependency (Dependency):
    def _check(self):
        return '{0} {1}'.format(
            self.or_pass['dependency'].full_name(),
            self.or_pass['version'])


for name,dependencies in [
        ('virtual-shell', (
            'bash',
            'dash',
            'ash',
            'zsh',
            'ksh',
            'csh',
            'tcsh',
            'sh',
            )),
        ('virtual-editor', (
            'emacs',
            'xemacs',
            'vim',
            'vi',
            'nano',
            'kate',
            'notepad++',
            )),
        ('virtual-browser', (
            'firefox',
            'google-chrome',
            'chromium',
            )),
        ('virtual-pypi-installer', (
            'easy_install',
            'pip',
            )),
        ]:
    CHECKER[name] = VirtualDependency(
        name=name, long_name=name, or_dependencies=dependencies)
del name, dependencies  # cleanup namespace


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
    try:
        passed = check(_sys.argv[1:])
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
        print()
        print_system_info()
        print_suggestions(instructor_fallback=True)
        _sys.exit(1)
