

KEYWORDS="process_worker_pool HTEX parsl multiprocessing"

for KEY in $KEYWORDS; do
    PIDs=`ps aux | grep $KEY | grep -v grep | awk '(NR>1) {print $2}'`
    echo "Grepping and killing $KEY processes"
    echo $PIDs

    for PID in $PIDs; do
        kill $PID
    done
    sleep 0.1
    echo
done
