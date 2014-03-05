#!/bin/bash

# Setup new bootcamp
#
# This script follow the steps in https://vimeo.com/87241285.
#
# Usage:
#
#     $ ./setup-bootcamp.sh bootcamp_id owner [instructor [instructor ...]]
#
# Note:
#
# This script will not store any password for security reasons.

# Functions

function create_repo {
    # TODO: Add description
    url=https://api.github.com/user/repos
    ghpages=https://${OWNER}.github.io/${BOOTCAMPID}
    curl -u "${OWNER}:${PASSWORD}" -i \
        -d "{\"name\":\"${BOOTCAMPID}\",\"homepage\":\"${ghpages}\"${REPO_DEFAULTS}}" \
        ${url}
}

function set_team {
    for instructor in ${INSTRUCTORS};
    do
        # This don't work.
        #
        # More information at http://stackoverflow.com/q/21698009/1802726.
        url=https://api.github.com/repos/${OWNER}/${BOOTCAMPID}/collaborators/${instructor}
        curl -u "${OWNER}:${PASSWORD}" -i -X PUT \
            ${url}
    done
}

function clone_and_push {
    git clone -b gh-pages -o bc \
        https://github.com/swcarpentry/bc.git ../${BOOTCAMPID}
    cd ../${BOOTCAMPID}
    git remote add origin \
        ${BOOTCAMPURL/https:\/\//https:\/\/${OWNER}@}.git
    # This is only supported by git >= 1.7.9
    #
    # More information: http://stackoverflow.com/a/5343146/1802726
    git config credential.username ${OWNER}
    # Try to avoid user be asked for password. Don't work properly.
    #
    # git push \
    #      --force-with-lease=gh-pages:gh-pages \
    #     --repo=${BOOTCAMPURL/https:\/\//https:\/\/${OWNER}:${PASSWORD}@}.git
    #
    # Try to use "[Here
    # string](http://www.gnu.org/software/bash/manual/bashref.html#Here-Strings)"
    # but it don't work either.
    #
    # git push origin gh-pages <<< ${PASSWORD}
    git push origin gh-pages
}

function end_message {
    echo "1. Repository create at ${BOOTCAMPURL}"
    echo "2. Repository cloned at ../${BOOTCAMPID}"
    echo ""
    echo "TODO:"
    echo ""
    echo "1. Update ../${BOOTCAMPID}/index.html"
    echo "2. Check the information at ../${BOOTCAMPID}/index.html"
    echo "3. Update the repository:"
    echo ""
    echo "    $ git push origin gh-pages"
}

function main {
    create_repo
    set_team
    clone_and_push
    end_message
}

# Check the number of input arguments
if test $# -lt 2
then
    echo "You must provide at least the bootcamp id and"
    echo "your GitHub username."
    exit 1
fi

# Setup environments variables
BOOTCAMPID=$1
OWNER=$2
BOOTCAMPURL=https://github.com/${OWNER}/${BOOTCAMPID}
INSTRUCTORS=${@:3}
REPO_DEFAULTS=,\"has_issues\":false,\"has_wiki\":false,\"has_downloads\":false

# Ask for password
echo "Enter the password for $OWNER:"
# Read the password from user.
read -s PASSWORD

# Call main function
main
