#!/bin/bash


PIDs=`ps aux | grep "from multiprocessing" | grep -v grep | grep $USER | awk '{ print $2 }'`
echo $PIDs

for PID in $PIDs; do
   ps aux | grep $PID | grep -v grep
   kill -9 $PID
done
