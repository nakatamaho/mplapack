#!/bin/bash

ROUTINEs="Rdot"

MP1s="_Float128 _Float64x dd"
MP2s="mpfr gmp qd"

for _routine in $ROUTINEs; do
    pos=-1
    if [ $_routine = "Rgemm" ]; then
        pos=4
    fi
    if [ $_routine = "Rgemv" ]; then
        pos=3
    fi
    if [ $_routine = "Rgetrf" ]; then
        pos=3
    fi
    if [ $_routine = "Rsyrk" ]; then
        pos=3
    fi
    if [ $_routine = "Rpotrf" ]; then
        pos=2
    fi
    if [ $_routine = "Raxpy" ]; then
        pos=2
    fi

    if [ $_routine = "Rdot" ]; then
        pos=2
    fi

    if [ $pos = -1 ]; then
        echo "error"
        exit
    fi
    echo "\subsection{{\tt $_routine} benchmarks}"
    while read line
    do
        IFS=, list=(${line})
        abbriv=${list[0]}
        corename=${list[1]}
        _dirname=${list[2]}
        mkdir -p ~/Desktop/eps
        (pushd 2022/$_dirname/logs ; grep -v "'double'" ${_routine}1.plt  | grep -v "set title" > ll ; sed -i 's/th 1/th 6/g' ll ; sed -i 's/set terminal pdf/set terminal eps/g' ll ; sed -i 's/plot/set key above\nplot/' ll ; gnuplot ll > ~/Desktop/eps/${_routine}1.${abbriv}.eps ; popd) >& /dev/null
        (pushd 2022/$_dirname/logs ; grep -v "'double'" ${_routine}2.plt  | grep -v "set title" > ll ; sed -i 's/th 1/th 6/g' ll ; sed -i 's/set terminal pdf/set terminal eps/g' ll ; sed -i 's/plot/set key above\nplot/' ll ; gnuplot ll > ~/Desktop/eps/${_routine}2.${abbriv}.eps ; popd) >& /dev/null
        echo "%%%%% $corename %%%%"
        echo "\subsubsection{{\tt $_routine} on ${corename}}"
        echo "In Figure~\ref{${_routine}1.${abbriv}}, we show the result of {\tt ${_routine}} performance for {\tt \_Float128}, {\tt \_Float64x} and {\tt double-double}, and in Figure~\ref{${_routine}2.${abbriv}} we show the result of {\tt ${_routine}} performance for {\tt MPFR}, {\tt GMP} and {\tt quad-double} on ${corename}."
        IFS=" "
        PEAKS=""
        for _mp1 in $MP1s; do
            _PEAK=`sort -k${pos}n 2022/$_dirname/logs/log.${_routine}.${_mp1} | tail -1 | awk "{print \\$${pos}}"`
            if [ `echo "$_PEAK < 100" | bc` = 1 ]; then
                PEAK=`printf "%3.1f\n" ${_PEAK}`
            else
                PEAK=`printf "%3.0f\n" ${_PEAK}`
            fi
            PEAKS="$PEAKS ${PEAK} MFlops,"
        done
        echo "The peak performances of the reference {\tt ${_routine}}s of {\tt \_Float128}, {\tt \_Float64x} and {\tt double-double} are$PEAKS respectively."

        PEAKS=""
        for _mp1 in $MP1s; do
            _PEAK=`sort -k${pos}n 2022/$_dirname/logs/log.${_routine}.${_mp1}_opt | tail -1 | awk "{print \\$${pos}}"`
            if [ `echo "$_PEAK < 100" | bc` = 1 ]; then
                PEAK=`printf "%3.1f\n" ${_PEAK}`
            else
                PEAK=`printf "%3.0f\n" ${_PEAK}`
            fi
            PEAKS="$PEAKS ${PEAK} MFlops,"
        done
        echo "The peak performances of simple OpenMP parallelized {\tt ${_routine}}s of {\tt \_Float128}, {\tt \_Float64x} and {\tt double-double} are$PEAKS respectively."

        PEAKS=""
        for _mp2 in $MP2s; do
            _PEAK=`sort -k${pos}n 2022/$_dirname/logs/log.${_routine}.${_mp2} | tail -1 | awk "{print \\$${pos}}"`
            if [ `echo "$_PEAK < 100" | bc` = 1 ]; then
                PEAK=`printf "%3.1f\n" ${_PEAK}`
            else
                PEAK=`printf "%3.0f\n" ${_PEAK}`
            fi
            PEAKS="$PEAKS ${PEAK} MFlops,"
        done
        echo "The peak performances of the reference {\tt ${_routine}}s of {\tt MPFR 512bit}, {\tt GMP 512bit} and {\tt quad-double} are$PEAKS respectively."

        PEAKS=""
        for _mp2 in $MP2s; do
            _PEAK=`sort -k${pos}n 2022/$_dirname/logs/log.${_routine}.${_mp2}_opt | tail -1 | awk "{print \\$${pos}}"`
            if [ `echo "$_PEAK < 100" | bc` = 1 ]; then
                PEAK=`printf "%3.1f\n" ${_PEAK}`
            else
                PEAK=`printf "%3.0f\n" ${_PEAK}`
            fi
            PEAKS="$PEAKS ${PEAK} MFlops,"
        done
        echo "The peak performances of simple OpenMP parallelized {\tt MPFR 512bit}, {\tt GMP 512bit} and {\tt quad-double} are$PEAKS respectively.\\\\"


        echo "\begin{figure}"
        echo "\caption{ {\tt ${_routine}} performance on ${corename} for {\tt \_Float128}, {\tt \_Float64x} and {\tt double-double} with/without simple OpenMP acceleration. }" 
        echo "\label{${_routine}1.${abbriv}}"
        echo "\begin{center}"
        echo "\includegraphics{${_routine}1.${abbriv}.eps}"
        echo "\end{center}"
        echo "\end{figure}"

        echo "\begin{figure}"
        echo "\caption{{\tt ${_routine}} performance on ${corename} for {\tt MPFR 512bit}, {\tt GMP 512bit} and {\tt quad-double} with/without simple OpenMP acceleration. }"
        echo "\label{${_routine}2.${abbriv}}"
        echo "\begin{center}"
        echo "\includegraphics{${_routine}2.${abbriv}.eps}"
        echo "\end{center}"
        echo "\end{figure}"
        echo ""
        echo ""
    done < table_abbriv.txt
done
