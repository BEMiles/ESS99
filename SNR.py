#!/usr/bin/env python

#Brittany E. Miles bemiles@ucla.edu 2013


usage= '''
        Usage: SNR.py rtt diam [-p pspin] [-d dec] [-a alb] [-h abmg]

        Mandatory Arguments
        ------------------------------
        rtt- Round trip time (seconds)
        diam- Diameter of your target (meters)
    
        Optional Arguments
        ------------------------------
        pspin- Spin peroid of your target (hours)
        dec- Declination of your target (degrees)
        alb- Albedo of your target (adimensional)
        abmg- Absolute magnitude of your target (adimenisional)

	If the diameter of an asteroid is not known you must
        provide an albedo and absolute magnitude. If used, no diam is 
        input is needed. In the absence of pspin and dec inputs the 
        following assumptions will be made: 
        pspin=2.4 hours  dec=18.0 degrees

        This program will output the following info for 
        both Arecibo and Goldstone telescopes.
 
        Telescope Name
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

rtt=float(sys.argv[1]) #round trip time (seconds)

diam=sys.argv[2] #diameter (meters)


                     
#optional arguments are found by searching through the command line

if len(sys.argv)>3:                     
    count=0
    while count<len(sys.argv)-1:

        if '-p' in sys.argv[count]:
            pspin=(float(sys.argv[count+1])*3600.0) #spin period (seconds)
        
        if '-d' in sys.argv[count]:
            dec=(float(sys.argv[count+1])*(math.pi/180.0)) #declination (in radians) 
        
        if '-h' in sys.argv[count]:
            abmg=float(sys.argv[count+1]) #absolute magnitude (no dimension)

        if '-a' in sys.argv[count]:                        
            alb=float(sys.argv[count+1]) #albedo (no dimension)
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

C=299792458.0 #meters per second

AU=149597871000.0 #meters

k=1.3806488E-23 #Boltzmann constant

ObjectDistance=((rtt*0.5)*C) #meters 

DTOR=math.pi/180 #degree to radians

try: 
    float(diam)
except ValueError:
    try:
        alb
    except NameError:
        alb=input('What is the albedo of your target?:')
    
    try:
        abmg
    except NameError:
        abmg=input('What is the absolute magnitude of your target?:')

    #Absolute Magnitude and albedo calculation for diameter
    diam=1329*(alb**-0.5)*(10**(-0.2*abmg))*10**3

else:
    diam=float(diam)

area=(math.pi)*((0.5*diam)**2) #meters squared
   
sigma=0.1*(area) #assumption of asteroid cross section




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

    
#Typical settings for Arecibo and Goldstone
#Things like frequency can be changed in real life.
#
#Freq-Frequency (hertz)
#Lambd-Wavelength (meters)
#Ae-Effective Area (meters^2)
#Tsys-Temperature of System (Kelvin)
#Pt-Transmitted Power (Watts)
#Down-Telescope Down Time between runs (seconds)
#Lat-Latitude (degrees)
#ViewRange-half of the part of a great circle that
#          is mapped by diurnal motion that the
#          telescope can follow. Or the half of
#          of the swath of a stationary telescope. (degrees)

class Arecibo:
    Freq=2380.0*(10**6)   
    Lambd=C/Freq          
    Ae=27600             
    Tsys=25                        
    Pt=900000             
    Down=5                
    Lat=18.23             
    ViewRange=19.69       
    
class Goldstone:
    Freq=8560.0*(10**6)  
    Lambd=C/Freq          
    Ae=2770               
    Tsys=25               
    Pt=450000             
    Down=5               
    Lat=35.24             
    ViewRange=90          

Telescope=[Arecibo,Goldstone]


count=0
while count<=len(Telescope)-1:
    
    #This block calculates the cosine of the the hourangle(cos(h))
    #If the object is visible, -1<cos(h)<1 a Track Time is calculated
    #or assigned and Signal to Noise Calculations are made. If an object
    #is not visible, the program will tell you that and won't do any
    #calculations.

 
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
        TrackTime=((hourangle/DTOR)/15.0)*2.0*3600.0 #seconds
        
        if 'Goldstone' in Telescope[count].__name__:    
            TrackTime=6*3600 #seconds


        

        def SNRCalculations():
            
            IT=rtt-Telescope[count].Down #Integration Time (seconds)
            
            runs=0.5*(TrackTime/rtt)                                      
            runs=int(runs)
    

            Bandwidth=(4*math.pi*diam)/(Telescope[count].Lambd*pspin) #Hertz
       
            PRN=(Telescope[count].Pt*(Telescope[count].Ae**2)*sigma)               
            PRD=((Telescope[count].Lambd**2)*4*math.pi*(ObjectDistance**4))         
        
            PowerReceived=(PRN/PRD) #Watts

            PowerNoise=k*Telescope[count].Tsys*Bandwidth #Watts   
    
            SNRperRun=(PowerReceived/PowerNoise)*((Bandwidth*IT)**0.5) #Adimensional

            SNRDaily=(SNRperRun*(runs**0.5)) #Adimensional
        

            print ' '
            print Telescope[count].__name__
            print '------------------------------------'
            print 'SNR Per Run: '+str(Round(SNRperRun,1))
            print 'No. of Runs: '+str(runs)
            print 'SNR Daily: '+str(Round(SNRDaily,1))
            print 'Bandwidth: '+str(Round(Bandwidth,1))+' Hertz'
            

        SNRCalculations()
    count=count+1
