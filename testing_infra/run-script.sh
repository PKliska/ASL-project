#!/usr/bin/env bash


#### !nix-shell --packages python3Full -i "bash"

set -e  # abort script at first error, when a command exits with non-zero status

printf "Activating python virtual env... (venv needs to be in root dir of project)\n"
source ../.venv/bin/activate

if [ -z "$1" ]; then # first arg is not set
    printf "Usage:\n    cd testing_infra\n    run-script.sh <PATH_TO_XONSH_SCRIPT>\n"
    exit 1
fi

xonsh "$1" "${@:2}" # 1st arg is the program to be run and pass the other args into the program
