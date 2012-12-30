#!/usr/bin/env python

"""Test script to check for required functionality.

Execute this code at the command line by typing:

  python swc-installation-test-2.py

How to get a command line:

- On OSX run this with the Terminal application.

- On Windows, go to the Start menu, select 'Run' and type 'cmd'
(without the quotes) to run the 'cmd.exe' Windows Command Prompt.

- On Linux, either use your login shell directly, or run one of a
  number of graphical terminals (e.g. 'xterm', 'gnome-terminal', ...).

Run the script and follow the instructions it prints at the end.

This script requires at least Python 2.6.  You can check the version
of Python that you have installed with 'swc-installation-test-1.py'.
"""

import distutils.ccompiler as _distutils_ccompiler
import importlib as _importlib
import logging as _logging
import os as _os
import platform as _platform
import re as _re
import subprocess as _subprocess
import sys as _sys


__version__ = '0.1'

# Comment out any entries you don't need
CHECKS = [
# Shell
    'bash',
# Editors
    'virtual-editor',
# Browsers
    'virtual-browser',
# Version control
    'git',
    'hg',              # Command line tool
    'mercurial',       # Python package
# Build tools and packaging
    'make',
    'easy_install',
    'setuptools',
# Testing
    'nosetests',       # Command line tool
    'nose',            # Python package
# SQL
    'sqlite3',         # Command line tool
    'sqlite3-python',  # Python package
# Python
    'python',
    'IPython',
    'numpy',
    'scipy',
    'matplotlib',
    'sympy',
    'Cython',
    'networkx',
    'mayavi.mlab',
    ]

CHECKER = {}


class DependencyError (Exception):
    def __init__(self, checker, message):
        self.checker = checker
        self.message = message

    def __str__(self):
        return 'check for {0} failed:\n{1}'.format(
            self.checker.full_name(), self.message)


def check(checks=None):
    failures = []
    if not checks:
        checks = CHECKS
    for check in checks:
        checker = CHECKER[check]
        _sys.stdout.write('check {0}...\t'.format(checker.full_name()))
        try:
            checker.check()
        except DependencyError as e:
            failures.append(e)
            _sys.stdout.write('fail\n')
        else:
            _sys.stdout.write('pass\n')
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
            self._check()
        except DependencyError as e:
            self._check_error = e  # cache for future calls
            raise

    def _check_dependencies(self):
        for dependency in self.and_dependencies:
            if not hasattr(dependency, 'check'):
                dependency = CHECKER[dependency]
            dependency.check()
        or_pass = not bool(self.or_dependencies)
        or_error = None
        for dependency in self.or_dependencies:
            if not hasattr(dependency, 'check'):
                dependency = CHECKER[dependency]
            try:
                dependency.check()
            except DependencyError as e:
                or_error = e
            else:
                or_pass = True
        if not or_pass:
            print(or_error)
            raise or_error

    def _check(self):
        version = self._get_version()
        parsed_version = None
        if hasattr(self, '_get_parsed_version'):
            parsed_version = self._get_parsed_version()
        if self.minimum_version:
            self._check_version(version=version, parsed_version=parsed_version)

    def _get_version(self):
        raise NotImplementedError(self)

    def _check_version(self, version, parsed_version=None):
        if not parsed_version:
            parsed_version = self._parse_version(version=version)
        if parsed_version < self.minimum_version:
            raise DependencyError(
                checker=self,
                message='outdated version of {0}: {1} (need >= {2})'.format(
                    self.full_name(), version,
                    self.version_delimiter.join(
                        str(part) for part in self.minimum_version)))

    def _parse_version(self, version):
        parsed_version = []
        for part in version.split(self.version_delimiter):
            try:
                parsed_version.append(int(part))
            except ValueError as e:
                raise NotImplementedError((version, part)) from e
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

    def __init__(self, command, version_option='--version',
                 version_regexp=None, version_stream='stdout', **kwargs):
        if 'name' not in kwargs:
            kwargs['name'] = command
        super(CommandDependency, self).__init__(**kwargs)
        self.command = command
        self.version_option = version_option
        if not version_regexp:
            regexp = r'([\d][\d{0}]*[\d])'.format(self.version_delimiter)
            version_regexp = _re.compile(regexp)
        self.version_regexp = version_regexp
        self.version_stream = version_stream

    def _get_version_stream(self):
        command = self.command + (self.exe_extension or '')
        try:
            p = _subprocess.Popen(
                [command, self.version_option],
                stdout=_subprocess.PIPE, stderr=_subprocess.PIPE,
                close_fds=True, shell=False, universal_newlines=True)
        except OSError as e:
            raise DependencyError(
                checker=self,
                message="could not find '{0}' executable".format(command),
                ) from e
        stdout,stderr = p.communicate()
        status = p.wait()
        if status:
            lines = [
                "failed to execute '{0} {1}':".format(
                    command, self.version_option),
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
        ('bash', 'Bourne Again Shell', (4, 0)),
        ('easy_install', 'Setuptools easy_install', None),
        ('git', 'Git', (1, 8, 0)),
        ('hg', 'Mercurial', (2, 0, 0)),
        ('make', None, None),
        ('sqlite3', 'SQLite 3', None),
        ('nosetests', 'Nose', (1, 0, 0)),
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
        try:
            package = _importlib.import_module(self.package)
        except ImportError as e:
            raise DependencyError(
                checker=self,
                message="could not import the '{0}' package".format(
                    self.package),
                ) from e
        try:
            version = package.__version__
        except AttributeError:
            version = 'unknown'
        return version


for package,name,long_name,minimum_version in [
        ('mercurial', None, 'Mercurial Python package',
         CHECKER['hg'].minimum_version),
        ('nose', None, 'Nose Python package',
         CHECKER['nosetests'].minimum_version),
        ('sqlite3', 'sqlite3-python', 'SQLite Python package',
         CHECKER['sqlite3'].minimum_version),
        ('IPython', None, None, None),
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


class VirtualDependency (Dependency):
    def _check(self):
        pass


for name,dependencies in [
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
        ]:
    CHECKER[name] = VirtualDependency(
        name=name, long_name=name, or_dependencies=dependencies)
del name, dependencies  # cleanup namespace


def print_system_info():
    print("If you do not understand why the above failures occurred,")
    print("copy and send the *entire* output (all info above and summary")
    print("below) to the instructor for help.")
    print()
    print('==================')
    print('System information')
    print('==================')
    print('os.name      : {0}'.format(_os.name))
    try:
        print('os.uname     : {0}'.format(_os.uname()))
    except:
        pass
    print('platform     : {0}'.format(_sys.platform))
    print('platform+    : {0}'.format(_platform.platform()))
    print('prefix       : {0}'.format(_sys.prefix))
    print('exec_prefix  : {0}'.format(_sys.exec_prefix))
    print('executable   : {0}'.format(_sys.executable))
    print('version_info : {0}'.format(_sys.version_info))
    print('version      : {0}'.format(_sys.version))
    print('environment  :')
    for key,value in sorted(_os.environ.items()):
        print('  {0}={1}'.format(key, value))
    print('==================')


if __name__ == '__main__':
    if not check(_sys.argv[1:]):
        print()
        print_system_info()
