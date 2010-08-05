#!/bin/bash

# A script to ensure that useq ChIPSeq application is called with the correct arguments.



function usage
{
    echo ""
    echo "usage: useq [[[-t treatment_dir -c control_dir -o output_dir -g genome_ver -f output_format ] | [-h]]"
    echo ""
    echo -e "       treatment_dir: \tTreatment alignment file directories"
    echo -e "       control_dir: \tControl alignment file directories"
    echo -e "       genome_ver: \tGenome version (e.g. H_sapiens_Feb_2009, M_musculus_Jul_2007)"
    echo -e "       output_format: \t'useq' or 'bar'; 'useq' is default"
    echo ""
}

# main 

if [ $# -lt 1 ]; then
    usage
    exit
fi

genome_ver='M_musculus_Jul_2007'
output_format='useq'



while [ "$1" != "" ]; do
    case $1 in
        -t | --treatment_dir )      shift
                                    treatment_dir=$1
                                    ;;
        -c | --control_dir )        shift
                                    control_dir=$1
                                    ;;
        -o | --output_dir )         shift
                                    output_dir=$1
                                    ;;
        -g | --genome_ver )         shift
                                    genome_ver=$1
                                    ;;
        -f | --output_format )      shift
                                    output_format=$1
                                    ;;
        -h | --help )               usage
                                    exit
                                    ;;
        * )                         usage
                                    exit 1
    esac
    shift
done

if [ -z "$treatment_dir" ]; then
    echo "No treatment directory specified"
    exit 1
fi

if [ -z "$control_dir" ]; then
    echo "No control directory specified"
    exit 1
fi

if [ -z "$output_dir" ]; then
    echo "No output directory specified"
    exit 1
fi

if [ "$output_format" = "bar" ]; then
    java -Xmx2G -jar /usr/local/bin/useq/Apps/ChIPSeq -y sam -v "$genome_ver" -t "$treatment_dir" -c "$control_dir" -s "$output_dir" -r /usr/bin/R64
else
    java -Xmx2G -jar /usr/local/bin/useq/Apps/ChIPSeq -y sam -v "$genome_ver" -t "$treatment_dir" -c "$control_dir" -s "$output_dir" -r /usr/bin/R64 -u
fi

