#!/usr/bin/env bash


#### !nix-shell --packages python3Full -i "bash"

set -e  # abort script at first error, when a command exits with non-zero status

printf "Activating python virtual env... (venv needs to be in root dir of project)\n"
source ../.venv/bin/activate

if [ -z "$1" ]; then # first arg is not set
    printf "Usage:\n    cd testing_infra\n    run-script.sh <PATH_TO_XONSH_SCRIPT>\n"
    exit 1
fi

pueue clean

printf "❗️ Already running tasks:\n"
pueue status status=running

printf "\n✅ Adding task to run queue\n"
task_id=$(pueue add --print-task-id "xonsh $1 ${@:2}") # 1st arg is the program to be run and pass the other args into the program

printf "❗️ Queue now (first 5 tasks):\n"
pueue status first 5
printf "\n"

kill_task() {
    pueue kill "$task_id" || pueue remove "$task_id"
    sleep 1
}

trap kill_task SIGINT

printf "❗️ Waiting for your task to start executing...\n"
pueue follow "$task_id"

