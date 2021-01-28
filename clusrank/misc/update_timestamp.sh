#!/bin/bash

# Note: this script is should be sourced from the project root directory

set -e

if [ "$(uname)" == "Darwin" ] || [ "$(uname)" == "Linux" ];
then
    printf "Updating date, version, and copyright year.\n"

    # define some variables
    yr=$(date +%Y)
    dt=$(date +%Y-%m-%d)
    cprt_R=misc/copyright.R
    cprt_cpp=misc/copyright.cpp
    citation=inst/CITATION
    version=$(grep "Version" DESCRIPTION | awk '{print $NF}')

    # update copyright year in the template headers
    regexp1="s/Copyright \(C\) 2015-[0-9]+/Copyright \(C\) 2015-$yr/"
    sed -i.bak -E "$regexp1" $cprt_R && rm $cprt_R.bak
    sed "s_#_/_g" $cprt_R > $cprt_cpp

    # update copyright year in all R scripts
    for Rfile in R/*.R
    do
        if ! grep -q 'Copyright (C)' $Rfile;
        then
            if [ $Rfile != "R/RcppExports.R" ];
            then
                cat $cprt_R $Rfile > .tmp
                mv .tmp $Rfile
            fi
        fi
        sed -i.bak -E "$regexp1" $Rfile && rm $Rfile.bak
    done

    # update copyright year in all C++ scripts
    for cppfile in src/*.cpp
    do
        if ! grep -q 'Copyright (C)' $cppfile;
        then
            if [ $cppfile != "src/RcppExports.cpp" ];
            then
                cat $cprt_cpp $cppfile > .tmp
                mv .tmp $cppfile
            fi
        fi
        sed -i.bak -E "$regexp1" $cppfile && rm $cppfile.bak
    done
    rm $cprt_cpp

    # update date in DESCRIPTION
    regexp2="s/Date: [0-9]{4}-[0-9]{1,2}-[0-9]{1,2}/Date: $dt/"
    sed -i.bak -E "$regexp2" DESCRIPTION && rm DESCRIPTION.bak

    # done
    printf "All updated.\n"
fi
