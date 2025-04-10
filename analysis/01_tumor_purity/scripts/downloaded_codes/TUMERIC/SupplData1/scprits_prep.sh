#!/bin/bash
for file in ./*.R
do
    chmod 755 $file
    sed -i 's///g' $file
done


