#!/bin/bash

block_file=$(find . -type f -name "blockMeshDict" | head -n 1)

if [[ ! -f "$block_file" ]]; then
    echo "Cannot find blockMeshDict."
    exit 1
fi

read a b c < <(grep -oP 'domain\s*\(\K[^\)]+' "$block_file")
echo "Extracted domain from $block_file: $a $b $c"

decomp_file=$(find . -type f -name "decomposeParDict" | head -n 1)

if [[ ! -f "$decomp_file" ]]; then
    echo "Cannot find decomposeParDict."
    exit 1
fi

read d e f < <(awk '/simpleCoeffs/,/}/ {if ($1=="n") {gsub(/[()]/,""); print $2, $3, $4}}' "$decomp_file")
echo "Extracted n from $decomp_file: $d $e $f"

cp_file=$(find . -type f -name "couplingProperties" | head -n 1)

if [[ ! -f "$cp_file" ]]; then
    echo "Cannot find couplingProperties."
    exit 1
fi

cp "$cp_file" "$cp_file.bak"
echo "Backup created: $cp_file.bak"

sed -i -E "s/(meshx\s*)[0-9.eE+-]+;/\1$a;/" "$cp_file"
sed -i -E "s/(meshy\s*)[0-9.eE+-]+;/\1$b;/" "$cp_file"
sed -i -E "s/(meshz\s*)[0-9.eE+-]+;/\1$c;/" "$cp_file"

sed -i -E "s/(decx\s*)[0-9.eE+-]+;/\1$d;/" "$cp_file"
sed -i -E "s/(decy\s*)[0-9.eE+-]+;/\1$e;/" "$cp_file"
sed -i -E "s/(decz\s*)[0-9.eE+-]+;/\1$f;/" "$cp_file"

echo "Updated meshx/meshy/meshz to $a / $b / $c"
echo "Updated decx/decy/decz to $d / $e / $f"

control_file=$(find . -type f -name "controlDict" | head -n 1)

if [[ ! -f "$control_file" ]]; then
    echo "Cannot find controlDict."
    exit 1
fi

max_x=$(awk '/vertices/,/^\)/ {
    gsub(/[()]/,"");
    if ($0 ~ /^[[:space:]]*[0-9.+-]/) {
        count++;
        if (count==3) print $1;
    }
}' "$block_file")

convert_to_meters=$(awk '/convertToMeters/ {print $2}' "$block_file" | tr -d ';')

if [[ -z "$max_x" || -z "$a" || -z "$convert_to_meters" ]]; then
    echo "Failed to extract max_x, domain x, or convertToMeters."
    exit 1
fi

length_ratio=$(echo "$max_x / $a * $convert_to_meters" | bc -l)

echo "Extracted max_x=$max_x, domain_x=$a, convertToMeters=$convert_to_meters"
echo "Computed x_plane=$length_ratio"

cp "$control_file" "$control_file.bak"
echo "Backup created: $control_file.bak"

sed -i -E "s/^\s*(x_plane\s+)[0-9.eE+-]+;/\1$length_ratio;/" "$control_file"

echo "Updated x_plane in controlDict to $length_ratio"
