#! /bin/bash

read_mbytes=()
rx_mbytes=()
wtime=()
eff=()
exe=()
read_rate=()
rx_rate=()

metrics() {
    local prmon=$1/prmon.json
    local output=$1/rdf_ttbar.out
    if [[ ! -f ${prmon} || ! -f ${output} ]] ; then
	echo "Job output not found. Exiting..." 1>&2
	exit 1
    fi
    grep '^processing took' ${output} &> /dev/null
    if [[ ! $? ]] ; then
	echo "Job output incomplete. Exiting..."1>&2
	exit 1
    fi
    local workers=$2
    rmb=$(grep -A 18 Max ${prmon} | awk '/read_bytes/ {sub(/,/, "", $2); print ($2 / 1024 ** 2)}')
    read_mbytes+=($rmb)
    rxmb=$(grep -A 18 Max ${prmon} | awk '/rx_bytes/ {sub(/,/, "", $2); print ($2 / 1024 ** 2)}')
    rx_mbytes+=($rxmb)
    u=$(grep -A 18 Max ${prmon} | awk '/utime/ {sub(/,/, "", $2); print $2}')
    s=$(grep -A 18 Max ${prmon} | awk '/stime/ {sub(/,/, "", $2); print $2}')
    w=$(grep -A 18 Max ${prmon} | awk '/wtime/ {sub(/,/, "", $2); print $2}')
    wtime+=($w)
    eff+=($(awk "BEGIN {print(100*($u+$s)/$w/$workers)}"))
    e=$(grep '^processing took' ${output} | awk '{print $3}')
    exe+=($e)
    read_rate+=($(awk "BEGIN {print($rmb / $e)}"))
    rx_rate+=($(awk "BEGIN {print($rxmb / $e)}"))

}

base=$1

jobs="${base} ${base}.*"
workers=$(echo $base | awk -F_ '{print($NF)}')
locstor=$(echo ${base} | grep 'local')
njobs=0
for job in $jobs; do
    if [ -d "$job" ] ; then
	njobs=$((njobs+1))
	metrics $job $workers
    fi
done
if [ "${njobs}" -lt 2 ] ; then
    echo "Only ${njobs} job found, cannot estimate errors. Exiting..." 1>&2
    exit 1
fi

echo "${njobs} jobs for ${workers} workers" 1>&2
#echo "wtime,wtime_err,CPU_eff,CPU_eff_err,rate,rate_err,tput,tput_err"
echo -n "${workers},"
echo -n ${wtime[@]} | awk '{for (i=1;i<=NF;i++)sum+=$i; avg=(sum/NF); for (i=1;i<=NF;i++)sum2+=($i-avg)^2;err=sqrt(sum2/(NF-1)); printf "%.1f,%.1f,", avg, err}'
echo -n ${eff[@]} | awk '{for (i=1;i<=NF;i++)sum+=$i; avg=(sum/NF); for (i=1;i<=NF;i++)sum2+=($i-avg)^2;err=sqrt(sum2/(NF-1)); printf "%.2f,%.2f,", avg, err}'
if [ -n "$locstor" ] ; then
    echo -n ${read_rate[@]} | awk '{for (i=1;i<=NF;i++)sum+=$i; avg=(sum/NF); for (i=1;i<=NF;i++)sum2+=($i-avg)^2;err=sqrt(sum2/(NF-1)); printf "%.1f,%.1f\n", avg, err}'
else
    echo -n ${rx_rate[@]} | awk '{for (i=1;i<=NF;i++)sum+=$i; avg=(sum/NF); for (i=1;i<=NF;i++)sum2+=($i-avg)^2;err=sqrt(sum2/(NF-1)); printf "%.1f,%.1f\n", avg, err}'
fi
