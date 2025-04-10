#!/bin/bash


for file in ./*
do
    chmod 755 $file
    sed -i 's/\r$//g' $file
done

