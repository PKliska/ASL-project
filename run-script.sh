#!/usr/bin/env nix-shell
#!nix-shell --packages python3Full -i "bash"

source ./.venv/bin/activate

if [ -z "$1" ]; then # first arg is not set
    printf "Usage:\n    run-script.sh <PATH_TO_XONSH_SCRIPT>\n"
    exit 1
fi

printf "\n\n"

xonsh "$1"