#!/bin/bash

## Setup new workshop
##
## This script follow the steps in https://vimeo.com/87241285.
##
## Usage:
##
##     $ ./setup-workshop.sh workshop_id description owner [instructor [instructor ...]]
##
## Parameters:
##
## - `workshop_id`: is the identifier provide by Software Carpentry for the
##   workshop that you will be running. It will be used for the name of the Git
##   repository and is a name like YYYY-MM-DD-site, e.g., `2014-03-31-esu`.
## - `description`: is the description of the workshop that you will be running.
##   It should be multi-word and must be enclose with single or double quotes.
## - `owner`: this is your GitHub username.
## - `instructors`: this is the GitHub username of folks that will be teaching
##    with you.
##
## Note:
##
## This script will ask the users for their GitHub password twice. The first time
## it will store the password in memory and use it when using GitHub API. The
## second time will be git requesting it to push the repository.

# Functions

function check_pwd {
    if test $(git remote --verbose 2>&1 | grep -c swcarpentry/bc) -ge 1
    then
        echo "Look like that you are inside bc repository."
        echo "Since cloning a repository inside another can create some problems,"
        echo "please move to outside bc repository and call this script again."
        exit 2
    fi
}

function create_repo {
    url=https://api.github.com/user/repos
    ghpages=https://${OWNER}.github.io/${WORKSHOPID}
    curl -f -u "${OWNER}:${PASSWORD}" -i \
        -d "{\"name\":\"${WORKSHOPID}\",\"description\":\"${DESCRIPTION}\",\"homepage\":\"${ghpages}\"${REPO_DEFAULTS}}" \
        ${url}
    if test $? -ne 0
    then
        echo "Can't create the remote repository. Aborting."
        exit $?
    fi
}

function set_team {
    for instructor in ${INSTRUCTORS};
    do
        url=https://api.github.com/repos/${OWNER}/${WORKSHOPID}/collaborators/${instructor}
        curl -f -u "${OWNER}:${PASSWORD}" -i -X PUT -d "{}" \
            ${url}
        if test $? -ne 0
        then
            echo "WARNING: can't add ${instructor} as collaborators."
        fi
    done
}

function clone_and_push {
    git clone -b gh-pages -o bc \
        https://github.com/swcarpentry/bc.git ${REPOSITORYPATH}
    cd ${REPOSITORYPATH}
    git remote add origin \
        ${WORKSHOPURL/https:\/\//https:\/\/${OWNER}@}.git
    # This is only supported by git >= 1.7.9
    #
    # More information: http://stackoverflow.com/a/5343146/1802726
    git config credential.username ${OWNER}
    # Try to avoid user be asked for password. Don't work properly.
    #
    # git push \
    #      --force-with-lease=gh-pages:gh-pages \
    #     --repo=${WORKSHOPURL/https:\/\//https:\/\/${OWNER}:${PASSWORD}@}.git
    #
    # Try to use "[Here
    # string](http://www.gnu.org/software/bash/manual/bashref.html#Here-Strings)"
    # but it don't work either.
    #
    # git push origin gh-pages <<< ${PASSWORD}
    git push -u origin gh-pages
}

function end_message {
    echo "LOG:"
    echo ""
    echo "1. Repository create at ${WORKSHOPURL}"
    echo "2. Repository cloned at ${REPOSITORYPATH}"
    echo ""
    echo "TODO:"
    echo ""
    echo "1. Update ${REPOSITORYPATH}/index.html"
    echo "2. Check the information at ${REPOSITORYPATH}/index.html"
    echo "3. Update the repository:"
    echo ""
    echo "    $ git push gh-pages"
}

function main {
    check_pwd
    create_repo
    set_team
    clone_and_push
    end_message
}

# Check the number of input arguments
if test $# -lt 3
then
    grep '^##' $0 | sed -r 's/## ?//'
    exit 1
fi

# Setup environments variables
WORKSHOPID=$1
REPOSITORYPATH=./$1
DESCRIPTION="$2"
OWNER=$3
WORKSHOPURL=https://github.com/${OWNER}/${WORKSHOPID}
INSTRUCTORS=${@:4}
REPO_DEFAULTS=,\"has_issues\":false,\"has_wiki\":false,\"has_downloads\":false

# Ask for password
echo "Enter the password for $OWNER:"
# Read the password from user.
read -s PASSWORD

# Call main function
main
