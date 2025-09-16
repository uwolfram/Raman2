# General

Code evaluates individual Raman spectra. Note, that it is a finger exercise where I tried to code some modules.

Code was originally developed 2018 for [2] and based on [1] so cite these papers when you use it: 
    [1] Mirzaali, M., Schwiedrzik, J., Thaiwichai, S., Best, J., Michler, J., Zysset, P. & Wolfram, U. 
        Mechanical properties of cortical bone and their relationships with age, gender, composition and 
        microindentation properties in the elderly. Bone 93, 196–211 (2016).
    [2] Hennige, S., Wolfram, U., Wickes, L., Murray, F., Roberts, J., Kamenos, N., Schofield, S., Groetsch, A., 
        Spiesz, E., Aubin-Tam, M. & Etnoyer, P. Crumbling Reefs and Cold-Water Coral Habitat Loss in a Future Ocean: 
        Evidence of “Coralporosis” as an Indicator of Habitat Integrity. Frontiers in Marine Science 7, 1–16 (2020).

To evaluate spectra bands need to be specified. To call the analyses use for example:
    python3 ramanAnal.py -f CINMS_001_live_1_Copy_mod_Copy.txt -p T -lb 1050 -ub 1120 -b CaCO3
    python3 ramanAnal.py -f AC1_P01.txt -p T -lb 1000 -ub 1150 -b v1PO4

    python3 ramanAnal.py -f AC1_P01.txt -p T -lb 1620 -ub 1700 -b amid1

## Evaluated spectra
Spectra to be evaluated could be in bone tisste:
    v1PO4 
    v2PO4
    Amide1
    Amide3
    PYD
    Lipids
    
Check [1] for baclground info! Bands can be chosen as lower bound - upper bound:
    v1PO4   930 - 980, final version: 900 - 1000
    v2PO4   410 - 460, final version: 370 - 500 # note: upper bound for v1PO4 was changed from 460 to 459 to ease detection
    amid1   1620 - 1700, final version: 1620 - 1720
    amid3   1215 - 1300, final version: 1150 - 1350
    CH3     1365 - 1390
    pyd     1660
    lipid   1298

# Input

file : 
    Spectral file as ascii txt
    
plot :
    Whether control plots should be fiven
    
lb:
    lower bound of peak
    
ub:
    upper bound of peak
    

# Returns
* prints FWHM, Integral, Peak, and Pos for the specified band
* the background corrected spectra
* writes CINMS_001_live_1_*.pdf to illustrate wavelet (Lorentzian function) used for analyses


        
# Example run

see above


# Misc

note ramanAnal.py can call two different background corrections that do essentially the same
backgroundCorrection() # own code
backgroundCorrection2() # imports baselineCorrection


if you wish to execute for many spectra, you could make an executable text file with the following content:

#! /bin/bash

for item in *.txt
do
    ./ramanAnal.py $item plot
done
