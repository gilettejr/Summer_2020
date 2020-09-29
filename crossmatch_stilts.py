import os

def crossmatch(galaxy,path='undefined'):
    
    if path=='cls_cut':
        
        os.system('./topcat -stilts tskymatch2 ifmt1=votable ifmt2=csv out=crossmatching/gaia/cls_cut' + galaxy +' ofmt=csv ra1=ra dec1=dec ra2=RA dec2=DEC error=1.0 in1=initial_data/gaia_pre/' + galaxy +' in2=crossmatching/ukirt_pre/cls_cut' + galaxy)
        
    else:
     
        os.system('./topcat -stilts tskymatch2 ifmt1=votable ifmt2=csv out=crossmatching/gaia/' + galaxy +' ofmt=csv ra1=ra dec1=dec ra2=RA dec2=DEC error=1.0 in1=initial_data/gaia_pre/' + galaxy +' in2=crossmatching/ukirt_pre/' + galaxy)
    
