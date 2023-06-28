#!/bin/bash
run_dolphot() {
# args: epoch_num flt1, flt2, flt3, flt4

# Mask data
#wfc3mask -ncombine=1 ${1}_sci.fits;  wfc3mask -ncombine=1 ${2}_sci.fits;
wfc3mask -ncombine=1 ${1}_flc.fits;  wfc3mask -ncombine=1 ${2}_flc.fits; wfc3mask -ncombine=1 ${3}_flc.fits; wfc3mask -ncombine=1 ${4}_flc.fits;

# Get the sky
#calcsky ${1}_sci 10 25 -64 2.25 2.00; calcsky ${2}_sci 10 25 -64 2.25 2.00;
calcsky ${1}_flc 10 25 -64 2.25 2.00; calcsky ${2}_flc 10 25 -64 2.25 2.00; calcsky ${3}_flc 10 25 -64 2.25 2.00; calcsky ${4}_flc 10 25 -64 2.25 2.00; 

# Run dolphot
dolphot F555W.phot -pphot.param | tee -a dolphot.log

}

run_dolphot if2901p8q if2901p8q if2901p9q if2901paq
