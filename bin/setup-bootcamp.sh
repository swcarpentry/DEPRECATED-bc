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
# This script will ask the users for their GitHub password twice. The first time
# it will store the password in memory and use it when using GitHub API. The
# second time will be git requesting it to push the repository.

# Functions

function create_repo {
    # TODO: Add description
    #
    # Description can be a little complex to handle in the command line, because
    # of that the script don't support it right now.
    url=https://api.github.com/user/repos
    ghpages=https://${OWNER}.github.io/${BOOTCAMPID}
    curl -u "${OWNER}:${PASSWORD}" -i \
        -d "{\"name\":\"${BOOTCAMPID}\",\"homepage\":\"${ghpages}\"${REPO_DEFAULTS}}" \
        ${url}
}

function set_team {
    for instructor in ${INSTRUCTORS};
    do
        url=https://api.github.com/repos/${OWNER}/${BOOTCAMPID}/collaborators/${instructor}
        curl -u "${OWNER}:${PASSWORD}" -i -X PUT -d "{}" \
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
    git push -u origin gh-pages
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
    echo "    $ git push gh-pages"
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
