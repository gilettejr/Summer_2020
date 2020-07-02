import os

def crossmatch(galaxy):
    
    os.system('./topcat -stilts tskymatch2 ifmt1=votable ifmt2=csv out=crossmatching/gaia/' + galaxy +' ofmt=csv ra1=ra dec1=dec ra2=RA dec2=DEC error=1.0 in1=crossmatching/gaia_pre/' + galaxy +' in2=crossmatching/ukirt_pre/' + galaxy)
    

    
