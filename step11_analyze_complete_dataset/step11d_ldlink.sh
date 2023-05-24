#!/bin/bash

expected_line="QueryGWASTraitRSNumberPosition(GRCh37)AllelesR2D'RiskAlleleEffectSize(95%CI)BetaorORP-value"
expected_error_message1="NoentriesintheGWASCatalogareidentifiedusingtheLDtraitsearchcriteria."
expected_error_message2="InputSNPlistdoesnotcontainanyvalidRSnumbersorcoordinates."


snps_dir="snps"
for file in $snps_dir/*.txt
do

  snps=$(cat "$file")
  output_filename="${file##*/}"
  output_filename="${output_filename%.txt}"
  output_filename="snps_output/${output_filename}_output.json"

  if [ -e "$output_filename" ]; then
    output_filename_exists=1
  else
    output_filename_exists=0
  fi

  retry=1

  while [ $retry -eq 1 ]
  do

    if [[ $output_filename_exists == 0 ]]; then
      response=$(curl -k -H "Content-Type: application/json" -X POST -d '{"snps": "'"$snps"'", "pop": "GBR", "r2_d": "r2", "r2_d_threshold": "0.8", "window": "500000", "genome_build": "grch37"}' 'https://ldlink.nci.nih.gov/LDlinkRest/ldtrait?token=69091934f7d8')

      error_message=$(echo "$response" | grep -oP '(?<="error": ")[^"]+')
      error_message=$(echo "$error_message" | tr -d '[:space:]')

      first_line=$(echo "$response" | head -n 1)
      first_line=$(echo "$first_line" | tr -d '[:space:]')
    fi

    if [[ $output_filename_exists == 1 ]]; then
      retry=0
    elif [[ $first_line == $expected_line ]]; then
      echo "$response" > "$output_filename"
      retry=0
    elif [[ $error_message == $expected_error_message1 ]]; then
      echo "$error_message"
      echo "nothing found" > "$output_filename"
      echo "$snps"
      retry=0
    elif [[ $error_message == $expected_error_message2 ]]; then
      echo "$error_message"
      echo "invalid name" > "$output_filename"
      echo "$snps"
      retry=0
    else
      echo "$error_message"
      retry=1
    fi
  done
done

#check boolean operations:
#zz="Input SNP list does not contain any valid RS numbers or coordinates. "
#zz2="Input SNP list does not contain any valid RS numbers or coordinates. " 
#if [ "$zz" == "$zz2" ]; then
#  zz3=1
#else
#  zz3=0
#fi
#echo "$zz3"

if [ "$first_line" == "$expected_line" ]; then
  zz3=1
else
  zz3=0
fi
echo "$zz3"

