#!/bin/bash

echo "lsof -i :8000"
echo
lsof -i :8000
echo
echo

PIDs=`lsof -i :8000 | awk '(NR>1) {print $2}'`

echo "now killing following processes:"
echo
for PID in $PIDs; do
   ps aux | grep $PID | grep -v grep
   kill -9 $PID
done
