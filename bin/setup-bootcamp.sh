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
        url=https://api.github.com/repos/${OWNER}/${BOOTCAMPID}/collaborators/${instructor}
        curl -u "${OWNER}:${PASSWORD}" -i -X PUT \
            ${url}
    done
}

function clone_and_push {
    git clone -b gh-pages -o bc \
        https://github.com/swcarpentry/bc ../${BOOTCAMPID}
    cd ../${BOOTCAMPID}
    git remote add origin ${BOOTCAMPURL}
    git push \
        --repo ${BOOTCAMPURL/https:\/\//https:\/\/${OWNER}:${PASSWORD}@} \
        gh-pages
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
