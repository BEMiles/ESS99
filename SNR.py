#!/usr/bin/env python

#Brittany E. Miles bemiles@ucla.edu

usage= '''
        Usage: SNR.py rtt diam [-p pspin] [-d dec] [-a alb]

        Mandatory Arguments
        ------------------------------
        rtt- Round trip time (seconds)
        diam- Diameter of your target (meters)
    
        Optional Arguments
        ------------------------------
        pspin- Spin peroid of you target (hours)
        dec- Declination of your target (degrees)
        alb- Albedo of your target      (adimenisional)

	      The argument alb is to be used if the diameter of the 
        asteroid is not known. If used, no diam is input is needed.

        In absence of pspin and dec inputs the following assumptions
        will be made: pspin=2.4 hours  dec=18.0 degrees

        This program will output info for both Arecibo and Goldstone.
 
        Telescope
        ------------------------------       
        -Signal to Noise Ratio per Run
        -Number of Runs you can make
        -Signal to Noise Ratio Daily
        -Target Bandwidth

        or the program will say that you cannot observe
        that target based on input declination
        
     '''
import sys,math


try:
    sys.argv[2]
except IndexError:
    print usage
    sys.exit()

rtt=float(sys.argv[1]) 		            #round trip time (seconds)

diam=sys.argv[2]                            #diameter (meters)


                     
#optional block
if len(sys.argv)>3:                     
    count=0
    while count<len(sys.argv)-1:

        if '-p' in sys.argv[count]:
            pspin=(float(sys.argv[count+1])*3600.0)                  #spin period (seconds)
        
        if '-d' in sys.argv[count]:
            dec=(float(sys.argv[count+1])*(math.pi/180.0))           #declination (in radians) 
        
        if '-a' in sys.argv[count]:
            alb=float(sys.argv[count+1])                             #albedo (no dimension)
        count=count+1

try:
    pspin                                  
except NameError:
    pspin=2.4*3600                        

try:
    dec                                     
except NameError:
    dec=(18.0*(math.pi/180.0))               

#Physical Constants and Derived Variables

C=299792458.0                               #meters per second

AU=149597871000.0                           #meters

k=1.3806488E-23                             #Boltzmann constant

ObjectDistance=((rtt*0.5)*C)                #meters 

DTOR=math.pi/180                            #degree to radians

V_SUN=-26.762                               #magnitude of sun

try:
    alb
except NameError:
    diam=float(diam)
else:
    #ALBEDO TO DIAMETER
    AbsoluteMag=5.51  #WRONG this is Iris 7's H
    diam=2*AU*(10**((AbsoluteMag-V_SUN)/-5.0))*alb**(-0.5)



area=(math.pi)*((0.5*diam)**2)              #meters squared
   
sigma=0.1*(area)                            #assumption of asteroid cross section



#Rounds off Quantities
def Round(Quantity,place):
            Quantity=str(Quantity)

            Split=Quantity.split('.')

            Integer=Split[0] 
            Fractional=Split[1] 

            if len(Split[1])>place:
    
                if int(Fractional[place:place+1])>=5:
                    end=str(int(Fractional[:place])+1)
                    Quantity=Split[0]+'.'+end
                    return Quantity
                elif int(Fractional[place:place+1])<5:
                    end=str(int(Fractional[:place]))
                    Quantity=Split[0]+'.'+end
                    return Quantity
            else:
                print float(Quantity)

    
#Telescope Stats Typical settings for Arecibo and Goldstone
#Things like frequency can be changed in real life.

class Arecibo:
    Freq=2380.0*(10**6)   #hertz
    Lambd=C/Freq          #meters
    Ae=27600              #meters^2
    Tsys=25               #Kelvin           
    Pt=900000             #Watts
    Down=5                #seconds
    Lat=18.23             #degrees
    ViewRange=19.69       #degrees
    
class Goldstone:
    Freq=8560.0*(10**6)   #hertz
    Lambd=C/Freq          #meters
    Ae=2770               #meters^2
    Tsys=25               #Kelvin
    Pt=450000             #Watts
    Down=5                #seconds
    Lat=35.24             #degrees
    ViewRange=90          #degrees

Telescope=[Arecibo,Goldstone]


count=0
while count<=len(Telescope)-1:
    
    #Track Time Calculations.(Rise and set of an asteroid)
    cos_hn=math.cos(Telescope[count].ViewRange*DTOR)-math.sin(dec)*math.sin(Telescope[count].Lat*DTOR)
    cos_hd=math.cos(dec)*math.cos(Telescope[count].Lat*DTOR)
    
    cos_h=(cos_hn/cos_hd)

    if cos_h>1 or cos_h<-1:
        print ' '
        print Telescope[count].__name__
        print '------------------------------------'
        print 'can\'t see your target at '+Telescope[count].__name__
    else:
        hourangle=math.acos(cos_h)
        TrackTime=((hourangle/DTOR)/15.0)*2.0*3600.0               #seconds
        
        if 'Goldstone' in Telescope[count].__name__:    
            TrackTime=6*3600                                       #seconds


        def SNRCalculations():
            
            IT=rtt-Telescope[count].Down #Integration Time
            
            runs=0.5*(TrackTime/rtt)                                      
            runs=int(runs)
    

            Bandwidth=(4*math.pi*diam)/(Telescope[count].Lambd*pspin)               #Hertz
       
            PRN=(Telescope[count].Pt*(Telescope[count].Ae**2)*sigma)               
            PRD=((Telescope[count].Lambd**2)*4*math.pi*(ObjectDistance**4))         
        
            PowerReceived=(PRN/PRD)                                                 #Watts

            PowerNoise=k*Telescope[count].Tsys*Bandwidth                            #Watts   
    
            SNRperRun=(PowerReceived/PowerNoise)*((Bandwidth*IT)**0.5)              #Adimensional

            SNRDaily=(SNRperRun*(runs**0.5))                                        #Adimensional
        

            print ' '
            print Telescope[count].__name__
            print '------------------------------------'
            print 'SNR Per Run: '+str(Round(SNRperRun,1))
            print 'No. of Runs: '+str(runs)
            print 'SNR Daily: '+str(Round(SNRDaily,1))
            print 'Bandwidth: '+str(Round(Bandwidth,1))+' Hertz'

        SNRCalculations()
    count=count+1
