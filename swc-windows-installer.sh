#!/bin/sh
#
# Software Carpentry Windows Installer
#
# Helps mimic a *nix environment on Windows with as little work as possible.
#
# The script:
# * Installs nano and makes it accessible from msysGit
# * Provides standard nosetests behavior (if Python and the nose
#   module are installed) for msysGit
#
# To use:
#
# 1. Install msysGit
#    http://code.google.com/p/msysgit/downloads/list?q=full+installer+official+git
# 2. Run swc_windows_installer.sh
#    You should be able to simply double click the file in Windows

msg()
{
	MESSAGE="${1}"
	echo "${1}" >&2
}

die()
{
	msg "${@}"
	exit 1
}

zip_install()
{
	URL="${1}"
	EXPECTED_MD5="${2}"
	INSTALL_DIRECTORY="${3}"
	ZIP=$(basename "${URL}")
	if test ! -d "${INSTALL_DIRECTORY}"
	then
		msg "download ${URL} and install into ${INSTALL_DIRECTORY}"
		curl -o "${ZIP}" "${URL}" || die "couldn't download ${URL}"
		MD5=$(md5sum "${ZIP}" | cut -d ' ' -f 1) || die "couldn't hash ${ZIP}"
		if test "${MD5}" != "${EXPECTED_MD5}"
		then
			die "downloaded ${URL} has the wrong MD5 hash: ${MD5} != ${EXPECTED_MD5}"
		fi
		mkdir -p "${INSTALL_DIRECTORY}" ||
			die "couldn't create ${INSTALL_DIRECTORY}"
		unzip "${ZIP}" -d "${INSTALL_DIRECTORY}" ||
			die "couldn't unzip ${ZIP} into ${INSTALL_DIRECTORY}"
	fi
}

install_nano()
{
	INSTALL_DIRECTORY="${1}"
	zip_install 'http://www.nano-editor.org/dist/v2.2/NT/nano-2.2.6.zip' \
		'4d8987b64a6be0f8de29d51ab5dc5a7a' \
		"${INSTALL_DIRECTORY}"
}

# Creates a terminal-based nosetests entry point for msysGit
create_nosetests_entry_point()
{
	INSTALL_DIRECTORY="${1}"
	if python -c 'import nose' 2> /dev/null
	then
		msg "create nosetests entry point in ${INSTALL_DIRECTORY}"
		if test ! -f "${INSTALL_DIRECTORY}"
		then
			mkdir -p "${INSTALL_DIRECTORY}"
		fi
		cat > "${INSTALL_DIRECTORY}/nosetests" <<-EOF
			#!/usr/bin/env/ python
			import sys
			import nose

			if __name__ == '__main__':
			    sys.exit(nose.core.main())
			EOF
	fi
}

add_swc_paths()
{
	FIRST_PATH="${1}"
	if ! grep "PATH.*:${FIRST_PATH}" ~/.bash_profile > /dev/null 2>&1
	then
		msg "add SWC paths to PATH in ~/.bash_profile"
		_PATH="\${PATH}"
		for EXTRA_PATH in "${@}"
		do
			_PATH="${_PATH}:${EXTRA_PATH}"
		done
		cat >> ~/.bash_profile <<-EOF

			# Add paths for Software-Carpentry-installed scripts and executables
			export PATH="${_PATH}"
			EOF
	fi
}

set_default_editor()
{
	EDITOR="${1}"
	if ! grep EDITOR ~/.bash_profile > /dev/null
	then
		msg "set the default EDITOR to '${EDITOR}'"
		cat >> ~/.bash_profile <<-EOF

			# Set the default editor
			export EDITOR='${EDITOR}'
			EOF
	fi
}

main()
{
	SWC_DIR=~/.swc
	BIN_DIR="${SWC_DIR}/bin"
	NANO_DIR="${SWC_DIR}/lib/nano"
	create_nosetests_entry_point "${BIN_DIR}" || die "couldn't create nosetests"
	install_nano "${NANO_DIR}" || die "couldn't install nano"
	add_swc_paths "${NANO_DIR}" "${BIN_DIR}" || die "couldn't add SWC paths"
	set_default_editor nano || die "couldn't set EDITOR"
}

main
