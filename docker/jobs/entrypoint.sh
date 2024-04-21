#!/bin/sh

until mkdir ${WORK_HOME}/.ssh; do
    echo "creating ${WORK_HOME}/.ssh directory..."
done

until
    ssh-keyscan ${IFS_SERVER_HOST:-domain.example} >>${WORK_HOME}/.ssh/known_hosts
do
    echo "adding IFS server to known hosts"
done

python cgan_ui/jobs.py
