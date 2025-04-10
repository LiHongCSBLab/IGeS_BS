#!/bin/bash
for file in ./*.py
do
    chmod 755 $file
    sed -i 's///g' $file
done


