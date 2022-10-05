#!/bin/bash

PORT=`cat PORT`

echo "lsof -i :$PORT"
echo
lsof -i :$PORT
echo
echo

PIDs=`lsof -i :$PORT | awk '(NR>1) {print $2}'`

echo "now killing following processes:"
echo
for PID in $PIDs; do
   ps aux | grep $PID | grep -v grep
   kill -9 $PID
done
