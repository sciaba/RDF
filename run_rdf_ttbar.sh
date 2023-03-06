#! /bin/bash

# Default values
afname="unl"
workers=4
fraction=50

usage() {
    echo "Usage : run_rdf_ttbar.sh [-n NFILES] -a AFNAME -w NWORKERS [-m] -h]"
    echo
    echo "        -n NFILES: specify the maximum number of files per dataset"
    echo "        -a AFNAME: specify the data source. Possible choices:"
    echo "                   cern-xrootd: from EOSCMS via xrootd"
    echo "                   cern-http: from EOSCMS via HTTPS"
    echo "                   cernbox-xrootd: from CERNBOX via xrootd"
    echo "                   cern-local: from the local filesystem of iopef01"
    echo "                   unl: from Mebraska via HTTP"
    echo "        -w NWORKERS: specify how many workers should be used"
    echo "        -m: use of the merged dataset should be used instead of the original"
    echo "        -h: this help message"
}

while getopts "n:a:w:mh" arg; do
    case $arg in
	n)
	    nfiles=$OPTARG
	    ;;
	a)
	    afname=$OPTARG
	    ;;
	w)
	    workers=$OPTARG
	    ;;
	m)
	    merged=1
	    ;;
	?|h)
	    usage
	    exit 1
	    ;;
    esac
done

sudo sysctl vm.drop_caches=3

if [ -z "${nfiles}" ] ; then
    nfiles='all'
    nf='-1'
else
    nf=${nfiles}
fi

if [ -z "${merged}" ] ; then
    WDIR="rdf_ttbar_orig_run_"
    datasets="ntuples.json"
else
    WDIR="rdf_ttbar_merged_run_"
    datasets="ntuples_merged.json"
fi
WDIR="${WDIR}${nfiles}_${afname}_${workers}"
if [ -d ${WDIR} ] ; then
    echo "Directory already exit. Exiting..."
    exit 1
fi
mkdir ${WDIR}
cat rdf_ttbar_template.py | sed "s/_NFILES_/${nf}/" | sed "s/_AFNAME_/${afname}/" | sed "s/_WORKERS_/${workers}/" | sed "s/_DATASETS_/${datasets}/" > ${WDIR}/rdf_ttbar_tmp.py
cd ${WDIR}
ln -s ../utils utils
ln -s ../helper.cpp .
ln -s ../ntuples.json .
ln -s ../ntuples_merged.json .
export EXTRA_CLING_ARGS="-O2"
prmon -i 5 -- python3 rdf_ttbar_tmp.py &> rdf_ttbar.out
