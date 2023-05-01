#! /bin/bash

# Default values
maxjobs=3
afname="unl"
workers=4

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
    echo "        -x hdd/ssd: read through the CERN X-Cache instance"
    echo "        -c: disable the xrootd connection multiplexing"
    echo "        -p: add export XRD_PARALLELEVTLOOP=10 to optimise xrootd I/O"
    echo "        -h: this help message"
}

while getopts "n:a:w:mx:cph" arg; do
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
	x)
	    xcache=$OPTARG
	    ;;
	m)
	    merged=1
	    ;;
	c)
	    mpoff=1
	    ;;
	p)
	    xrdopt=1
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

if [ -n "${xrdopt}" ] ; then
    export XRD_PARALLELEVTLOOP=10
    WDIR="${WDIR}xrd_"
fi

if [ -z "${mpoff}" ] ; then
    mpoff=0
else
    WDIR="${WDIR}nomp_"
fi

if [ -z "${xcache}" ] ; then
    xcache='False'
elif [ "${xcache}" == 'hdd' ] ; then
    WDIR="${WDIR}xcache_"
elif [ "${xcache}" == 'ssd' ] ; then
    WDIR="${WDIR}xcachessd_"
fi

WDIR="${WDIR}${nfiles}_${afname}_${workers}"
if [[ -d "${WDIR}.${maxjobs}" ]] ; then
    echo "Enough jobs already. Exiting..."
    exit 1
fi
pre=$(ls -1d ${WDIR}.* 2> /dev/null | tail -1)
if [ -n "$pre" ] ; then
    pre=${pre: -1}
else
    pre=0
fi
new=$((pre + 1))
if [ -d $WDIR ] ; then
    echo "Moving $WDIR to $WDIR.$new"
    mv $WDIR $WDIR.$new
fi
mkdir ${WDIR}
if [ $? != 0 ] ; then
    echo "Could not create test dir. Exiting..."
    exit 1
fi
cat rdf_ttbar_template.py | sed "s/_NFILES_/${nf}/" | sed "s/_AFNAME_/${afname}/" | sed "s/_WORKERS_/${workers}/" | sed "s/_DATASETS_/${datasets}/" | sed "s/_XCACHE_/${xcache}/" | sed "s/_MPOFF_/${mpoff}/" > ${WDIR}/rdf_ttbar_tmp.py
cd ${WDIR}
ln -s ../utils utils
ln -s ../helper.cpp .
ln -s ../ntuples.json .
ln -s ../ntuples_merged.json .
export EXTRA_CLING_ARGS="-O2"
export XRD_APPNAME="AGCRDF"
env > env.out
prmon -i 5 -- python3 rdf_ttbar_tmp.py &> rdf_ttbar.out
