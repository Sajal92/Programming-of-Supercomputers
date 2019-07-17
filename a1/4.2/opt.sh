#!/bin/bash
g_or_i="z"
base=1
comb_file_name=""
h_or_p_or_c="z"

while [[ "$g_or_i" != "g" && "$g_or_i" != "i" ]]
do
    echo "GNU or Intel? (g/i): "
    read g_or_i
done
if [ "$g_or_i" = "g" ]; then
    comb_file_name="opts-gnu.txt"
    base=380.30533
else
    comb_file_name="opts-intel.txt"
    base=498.27622
fi

while [[ "$h_or_p_or_c" != "h" && "$h_or_p_or_c" != "p" && "$h_or_p_or_c" != "c" ]]
do
    echo "Harvest or plant or clean? (h/p/c): "
    read h_or_p_or_c
done

if [ "$h_or_p_or_c" = "h" ]; then
    max_fom=$base
    max_id="00"
    max_opts="00"
    max_speedup=1
    while IFS='' read -r line || [[ -n "$line" ]]; do
        IFS=' ' read -r id opts <<< "$line"

        fom=$(cat "$id"/ser_"$id".out | grep "FOM" | cut -f2 -d'=' | sed -e 's/^[[:space:]]*//' | cut -f1 -d' ')
        speedup=$(bc<<<"scale=2;$fom/$base")
        echo "Speed-up: $speedup" >> "$id"/ser_"$id".out

        if (( $(echo "$fom > $max_fom" | bc -l) )); then
            max_fom=$fom
            max_id=$id
            max_opts=$opts
            max_speedup=$speedup
        fi
    done < "$comb_file_name"
    echo "$max_id $max_opts $max_fom $max_speedup"
elif [ "$h_or_p_or_c" = "p" ]; then
    while IFS='' read -r line || [[ -n "$line" ]]; do
        IFS=' ' read -r id opts <<< "$line"
    
        mkdir $id
        cat Makefile | sed 's/^\(CXX\|LD\)\(FLAGS =\) \(.*\)/\1\2 '"$opts"' \3/' > "$id"/Makefile
        cat ser.sh | sed 's/\(output\|error\) = .*$/\1 = ser_'"$id"'.out/' > "$id"/ser.sh
        cp lul* $id
        ( cd $id && make && llsubmit ser.sh ) 2>&1 | tee -a output.log 
    done < "$comb_file_name"
else
    rm -r 0*
    rm -r 1*
    rm -r 2*
    rm -r 3*
    if [ "$g_or_i" = "g" ]; then
        rm -r 4*
        rm -r 5*
        rm -r 6*
    fi
fi
