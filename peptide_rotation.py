

#This script is used to vary the phi and psi angles systematically at one degree interval, from -180 to +180 degrees, within an L-ALA two-linked peptide unit.
#It also classifies each of the conformers as disallowed, partially allowed or fully allowed, based on a set of short contact criteria
import math
import numpy as np
import sys

#Function to calculate dihedral angle formed by 4 points in space
def FindDihedralAngle(A,B,C,D):
    BCi=C[0]-B[0]
    BCj=C[1]-B[1]
    BCk=C[2]-B[2]

    BAi=A[0]-B[0]
    BAj=A[1]-B[1]
    BAk=A[2]-B[2]

    CDi=D[0]-C[0]
    CDj=D[1]-C[1]
    CDk=D[2]-C[2]

    Q1i=(BCj*BAk)-(BCk*BAj)
    Q1j=(BCk*BAi)-(BCi*BAk)
    Q1k=(BCi*BAj)-(BCj*BAi)

    Q2i=(BCj*CDk)-(BCk*CDj)
    Q2j=(BCk*CDi)-(BCi*CDk)
    Q2k=(BCi*CDj)-(BCj*CDi)
    magQ1=math.sqrt((Q1i*Q1i)+(Q1j*Q1j)+(Q1k*Q1k))
    Q1i=Q1i/magQ1
    Q1j=Q1j/magQ1
    Q1k=Q1k/magQ1

    magQ2=math.sqrt((Q2i*Q2i)+(Q2j*Q2j)+(Q2k*Q2k))
    Q2i=Q2i/magQ2
    Q2j=Q2j/magQ2
    Q2k=Q2k/magQ2

    Q1dotQ2=(Q1i*Q2i)+(Q1j*Q2j)+(Q1k*Q2k)
    chi=math.acos(Q1dotQ2)
    chinew=math.degrees(chi)
    
    Q1=np.array([Q1i,Q1j,Q1k])
    Q2=np.array([Q2i,Q2j,Q2k])
    Q1crossQ2=np.cross(Q1,Q2)
    magBC=math.sqrt((BCi*BCi)+(BCj*BCj)+(BCk*BCk))
    unitBCi=BCi/magBC
    unitBCj=BCj/magBC
    unitBCk=BCk/magBC
    unitBC=np.array([unitBCi,unitBCj,unitBCk])
    anglesign=np.dot(Q1crossQ2,unitBC)
    if anglesign<0:
        chinew=chinew*-1
        
    return int(round(chinew))

#Function to find the coordinates of a point (x,y,z) after rotation by <angle> degrees about an axis with direction cosines (l,m,n), passing through point (a,b,c)
def resultantpoint(l,m,n,a,b,c,x,y,z,angle):
    angle=math.radians(angle)
    term1=(a*(m*m+n*n)-l*(b*m+c*n-l*x-m*y-n*z))*(1-math.cos(angle))
    term2=x*math.cos(angle)
    term3=(b*n-c*m-n*y+m*z)*math.sin(angle)
    x_new=term1+term2+term3
    term1=(b*(l*l+n*n)-m*(a*l+c*n-l*x-m*y-n*z))*(1-math.cos(angle))
    term2=y*math.cos(angle)
    term3=(c*l-a*n+n*x-l*z)*math.sin(angle)
    y_new=term1+term2+term3
    term1=(c*(l*l+m*m)-n*(a*l+b*m-l*x-m*y-n*z))*(1-math.cos(angle))
    term2=z*math.cos(angle)
    term3=(a*m-b*l-m*x+l*y)*math.sin(angle)
    z_new=term1+term2+term3
    newcoords=[0 for z in range(3)]
    newcoords[0]=x_new
    newcoords[1]=y_new
    newcoords[2]=z_new
    return newcoords

#Function to calculate Euclidian distance between points a1 and a2
def FindDistance(a1,a2):
    return round(math.sqrt(((a2[0]-a1[0])*(a2[0]-a1[0]))+((a2[1]-a1[1])*(a2[1]-a1[1]))+((a2[2]-a1[2])*(a2[2]-a1[2]))),3)

#Function to compare the interatomic distance between atoms a1 and a2 with their normal and outer contact limit d1 and d2 respectively
def CheckShortContacts(a1,a2,d1,d2):
    dist=FindDistance(a1,a2)
    if dist<=d1 and dist>=d2:
        return "partially allowed: "+str(dist)
    elif dist>d1:
        return "fully allowed: "+str(dist)
    elif dist<d2:
        return "disallowed: "+str(dist)

#Function to calculation the direction cosines of a line passing through points A and B    
def FindDirCosines(A,B):
    delx=B[0]-A[0]
    dely=B[1]-A[1]
    delz=B[2]-A[2]
    denom=math.sqrt(delx*delx+dely*dely+delz*delz)
    l=delx/denom
    m=dely/denom
    n=delz/denom
    a=[l,m,n]
    return a
    
#Function to bring the phi and psi values to zero in a two-linked peptide unit
def BringToZero(tempatoms):
    phi=FindDihedralAngle(tempatoms['C1'],tempatoms['N2'],tempatoms['CA2'],tempatoms['C2'])
    psi=FindDihedralAngle(tempatoms['N2'],tempatoms['CA2'],tempatoms['C2'],tempatoms['N3'])
    #bring phi to 0
    dircos=FindDirCosines(tempatoms['N2'], tempatoms['CA2'])
    atoms['CB']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['CB'][0],tempatoms['CB'][1],tempatoms['CB'][2],-phi)
    atoms['CA2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],-phi)
    atoms['HA']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['HA'][0],tempatoms['HA'][1],tempatoms['HA'][2],-phi)
    atoms['C2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['C2'][0],tempatoms['C2'][1],tempatoms['C2'][2],-phi)
    tempatoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['O2'][0],tempatoms['O2'][1],tempatoms['O2'][2],-phi)
    tempatoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['N3'][0],tempatoms['N3'][1],tempatoms['N3'][2],-phi)
    tempatoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['H3'][0],tempatoms['H3'][1],tempatoms['H3'][2],-phi)
    tempatoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['CA3'][0],tempatoms['CA3'][1],tempatoms['CA3'][2],-phi)
    #bring psi to 0
    dircos=FindDirCosines(tempatoms['CA2'],tempatoms['C2'])
    atoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['O2'][0],tempatoms['O2'][1],tempatoms['O2'][2],-psi)
    atoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['N3'][0],tempatoms['N3'][1],tempatoms['N3'][2],-psi)
    atoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['H3'][0],tempatoms['H3'][1],tempatoms['H3'][2],-psi)
    atoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['CA3'][0],tempatoms['CA3'][1],tempatoms['CA3'][2],-psi)


#Function to generate the different conformers of a two-linked peptide unit with phi and psi ranging from -180 to +180    
def MapsGenerator(structname):
    shrt_contact_file=open(structname+"_short_conts.txt","w") #This file stores all the short-contact information for partially allowed and disallowed conformers
    mapfile=open(structname+"_rmap.csv","w") #This file stores information about each integral phi-psi combination - whether it is fully allowed (recorded as 1), partially allowed (recorded as 2) or disallowed (recorded as 0)
    mapfile.write("phi,psi,a/d/p\n")
    BringToZero(atoms)
    
    for i in range(-180,181):
        for j in range(-180,181):
            finalres=''
            #psi rotation
            dircos=FindDirCosines(atoms['CA2'],atoms['C2'])
            atoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['O2'][0],atoms['O2'][1],atoms['O2'][2],1)
            atoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['N3'][0],atoms['N3'][1],atoms['N3'][2],1)
            atoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['H3'][0],atoms['H3'][1],atoms['H3'][2],1)
            atoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['CA3'][0],atoms['CA3'][1],atoms['CA3'][2],1)
            phi=FindDihedralAngle(atoms['C1'],atoms['N2'],atoms['CA2'],atoms['C2'])
            psi=FindDihedralAngle(atoms['N2'],atoms['CA2'],atoms['C2'],atoms['N3'])
            shrt_contact_file.writelines("\nphi="+str(phi)+"  psi="+str(psi)+"\n")
            
            #start checking for short contacts among non-bonded atoms within two-linked peptide unit
            res=CheckShortContacts(atoms['CA1'],atoms['H3'],2.40,2.20)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('CA1...H3 : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['HA'],2.40,2.20)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('C1...HA : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['C2'],3.00,2.90)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('C1...C2 : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['O2'],2.80,2.70)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('C1...O2 : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['N3'],2.90,2.80)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('C1...N3 : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['H3'],2.40,2.20)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('C1...H3 : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['CB'],3.20,3.00)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('C1...CB : '+res+'\n')
            
            res=CheckShortContacts(atoms['O1'],atoms['CB'],2.80,2.70)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('O1...CB : '+res+'\n')
            res=CheckShortContacts(atoms['O1'],atoms['HA'],2.40,2.20)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('O1...HA : '+res+'\n')
            res=CheckShortContacts(atoms['O1'],atoms['C2'],2.80,2.70)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('O1...C2 : '+res+'\n')    
            res=CheckShortContacts(atoms['O1'],atoms['O2'],2.70,2.60)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('O1...O2 : '+res+'\n')  
            res=CheckShortContacts(atoms['O1'],atoms['N3'],2.70,2.60)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('O1...N3 : '+res+'\n')  
            res=CheckShortContacts(atoms['O1'],atoms['H3'],2.40,2.20)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('O1...H3 : '+res+'\n') 
            res=CheckShortContacts(atoms['O1'],atoms['CA3'],2.80,2.70)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('O1...CA3 : '+res+'\n')  

            res=CheckShortContacts(atoms['N2'],atoms['O2'],2.70,2.60)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('N2...O2 : '+res+'\n')
            res=CheckShortContacts(atoms['N2'],atoms['N3'],2.70,2.60)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('N2...N3 : '+res+'\n')  
            res=CheckShortContacts(atoms['N2'],atoms['H3'],2.40,2.20)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('N2...H3 : '+res+'\n')
            
            res=CheckShortContacts(atoms['H2'],atoms['HA'],2.00,1.90)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('H2...HA : '+res+'\n')   
            res=CheckShortContacts(atoms['H2'],atoms['C2'],2.40,2.20)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('H2...C2 : '+res+'\n')
            res=CheckShortContacts(atoms['H2'],atoms['O2'],2.40,2.20)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('H2...O2 : '+res+'\n')
            res=CheckShortContacts(atoms['H2'],atoms['N3'],2.40,2.20)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('H2...N3 : '+res+'\n')
            res=CheckShortContacts(atoms['H2'],atoms['H3'],2.00,1.90)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('H2...H3 : '+res+'\n')
            res=CheckShortContacts(atoms['H2'],atoms['CB'],2.40,2.20)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('H2...CB : '+res+'\n')
            
            res=CheckShortContacts(atoms['CB'],atoms['O2'],2.80,2.70)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('CB...O2 : '+res+'\n')
            res=CheckShortContacts(atoms['CB'],atoms['N3'],2.90,2.80)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('CB...N3 : '+res+'\n')
            res=CheckShortContacts(atoms['CB'],atoms['H3'],2.40,2.20)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('CB...H3 : '+res+'\n')
            
            res=CheckShortContacts(atoms['HA'],atoms['O2'],2.40,2.20)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('HA...O2 : '+res+'\n')
            
            res=CheckShortContacts(atoms['HA'],atoms['N3'],2.40,2.20)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('HA...N3 : '+res+'\n')
            res=CheckShortContacts(atoms['HA'],atoms['H3'],2.00,1.90)
            finalres=finalres+res[0]
            if res[0]=='p' or res[0]=='d':
                shrt_contact_file.write('HA...H3 : '+res+'\n')
            
            if 'd' not in finalres and 'p' not in finalres:
                shrt_contact_file.write("This conformation is fully allowed\n--------------------------------------------------------\n")
                mapfile.write(str(phi)+" , "+str(psi)+" , 1\n")
            elif 'd' in finalres:
                shrt_contact_file.write("This conformation is disallowed\n--------------------------------------------------------\n")
                mapfile.write(str(phi)+" , "+str(psi)+" , 0\n")
            else:
                shrt_contact_file.write("This conformation is partially allowed\n--------------------------------------------------------\n")
                mapfile.write(str(phi)+" , "+str(psi)+" , 2\n") 
        
        #phi rotation    
        dircos=FindDirCosines(atoms['N2'], atoms['CA2'])
        atoms['CB']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['CB'][0],atoms['CB'][1],atoms['CB'][2],1)
        atoms['CA2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],1)
        atoms['HA']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['HA'][0],atoms['HA'][1],atoms['HA'][2],1)
        atoms['C2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['C2'][0],atoms['C2'][1],atoms['C2'][2],1)
        atoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['O2'][0],atoms['O2'][1],atoms['O2'][2],1)
        atoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['N3'][0],atoms['N3'][1],atoms['N3'][2],1)
        atoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['H3'][0],atoms['H3'][1],atoms['H3'][2],1)
        atoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['CA3'][0],atoms['CA3'][1],atoms['CA3'][2],1)   
            
    shrt_contact_file.close()
    mapfile.close()


#Execution starts from here

atoms=dict()

#L-ALA two-linked peptide unit is given as input through command line. It should be in PDB format

#Please Note:
''' The naming of atoms within the two-linked peptide unit should strictly follow the convention given below:

     O1      CB       H3
     ||       |        |
CA1--C1--N2--CA2--C2--N3--CA3
          |   |   || 
         H2  HA   O2 
'''
#All of the above atoms should be present
input_structure=sys.argv[1]
ipfile=open(input_structure,"r")
for line in ipfile:
    if "ATOM" in line:
        lineparts=line.split()
        xcoord=float(lineparts[6])
        ycoord=float(lineparts[7])
        zcoord=float(lineparts[8])
        atoms[lineparts[2]]=[xcoord,ycoord,zcoord]  
ipfile.close()
MapsGenerator()
    
        
            
            
        
        