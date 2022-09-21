#!/bin/bash

infile=$1
outfile=$2

mv $infile $outfile
mv "${infile}.tbi" "${outfile}.tbi"
