
for i in `seq 10`
do
    echo $i
    cmsRun TunePythia84jForUE_cfg.py nFile=$i > /dev/null 2>&1 &
done
