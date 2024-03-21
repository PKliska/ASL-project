#!/usr/bin/env bash


#### !nix-shell --packages python3Full -i "bash"

set -e  # abort script at first error, when a command exits with non-zero status

source ../.venv/bin/activate

if [ -z "$1" ]; then # first arg is not set
    printf "Usage:\n    run-script.sh <PATH_TO_XONSH_SCRIPT>\n"
    exit 1
fi

printf "\n\n"

xonsh "$1"