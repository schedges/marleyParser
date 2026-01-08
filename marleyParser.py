#Usage: python parseMarley.py <inputFileName.ascii> <outputFileName.root>
#
#
#Hists:
# - Neutrino energies
# - Excitation energy
# - Lepton energy
# - Total EM energy
# - Total neutron energy
# - Single neutron event energy
# - Double neutron event energy (energies of individual neutrons in these events)
# - Triple neutron event energy (energies of individual neutrons in these events)
# - Gamma energies
# - Number of neutrons
#
#TTree:
# - Neutrino energy
# - Excitation energy
# - Lepton energy
# - EM energy
# - Neutron energy
# - nGammas
# - nNeutrons
#
# The general output format from MARLEY is:
#   Num initial particles, num final particles, excitation energy transferred to nucleus
#   pdf code, Four momentum 0-4, mass, charge (initial particles, each on own line)
#   pdf code, Four momentum 0-4, mass, charge (final particles, each on own line)

import sys
import ROOT
import math
import array
import numpy

oldFormat=0

if not len(sys.argv)==3:
  print("Usage error: python parseMarley.py <inputFile.ascii> <outputName.root>")
  sys.exit(0)

#Constants, in MeV
restMass_e=0.511
restMass_n=939.56
restMass_p=938.27
restMass_nu=0.
restMass_alpha=3727.38

#Histogram settings
maxEnergy=250 #MeV
nBins=1000

#Default in MARLEY
neutrinoDirection=(0,0,1)

#Helper function to calculate the angle between two vectors from
#https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/
def unit_vector(vector):
    return vector / numpy.linalg.norm(vector)

def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return numpy.arccos(numpy.clip(numpy.dot(v1_u, v2_u), -1.0, 1.0))

#https://arxiv.org/pdf/1412.4417.pdf for Xe
def getQuenchedEnergy(E0):
  Z=54
  kappa = 0.1394
  E0_keV=E0*1000
  epsilon = 11.5*(E0_keV)*math.pow(Z,-7./3.)
  g = 3*math.pow(epsilon,0.15) + 0.7*math.pow(epsilon,0.6) + epsilon
  L = kappa * g / (1 + kappa*g) #quenched energy
  return L*E0
  
#Histograms
#Input to scattering
neutrinoEnergyHist=ROOT.TH1D("neutrinoEnergyHist",";Incident neutrino energy (MeV);Counts",nBins,0,maxEnergy)
#Nuclear excitations
excitationEnergyHist=ROOT.TH1D("excitationEnergyHist",";Excitation energy (MeV);Counts",nBins,0,maxEnergy)
#Output particles
leptonEnergyHist=ROOT.TH1D("leptonEnergyHist",";Outgoing electron energy (MeV);Counts",nBins,0,maxEnergy)
leptonEnergyOneNeutronHist=ROOT.TH1D("leptonEnergyOneNeutronHist",";Outgoing electron energy (MeV);Counts",nBins,0,maxEnergy)
leptonEnergyTwoNeutronHist=ROOT.TH1D("leptonEnergyTwoNeutronHist",";Outgoing electron energy (MeV);Counts",nBins,0,maxEnergy)
leptonEnergyThreeNeutronHist=ROOT.TH1D("leptonEnergyThreeNeutronHist",";Outgoing electron energy (MeV);Counts",nBins,0,maxEnergy)
finalNeutrinoEnergyHist=ROOT.TH1D("finalNeutrinoEnergyHist",";Outgoing neutrino energy (MeV);Counts",nBins,0,maxEnergy)
emEnergyHist=ROOT.TH1D("emEnergyHist",";Total electromagnetic energy (MeV);Counts",nBins,0,maxEnergy)
visibleEnergyHist=ROOT.TH1D("visibleEnergyHist",";Total visible energy (MeV);Counts",nBins,0,maxEnergy)
neutronEnergyHist=ROOT.TH1D("neutronEnergyHist",";Neutron energy (MeV);Counts",nBins,0,maxEnergy)
singleNeutronEnergyHist=ROOT.TH1D("singleNeutronEnergyHist",";Neutron energy (MeV);Counts",nBins,0,maxEnergy)
doubleNeutronEnergyHist=ROOT.TH1D("doubleNeutronEnergyHist",";Neutron energy (MeV);Counts",nBins,0,maxEnergy)
tripleNeutronEnergyHist=ROOT.TH1D("tripleNeutronEnergyHist",";Neutron energy (MeV);Counts",nBins,0,maxEnergy)
quadrupleNeutronEnergyHist=ROOT.TH1D("quadrupleNeutronEnergyHist",";Neutron energy (MeV);Counts",nBins,0,maxEnergy)
gammaEnergyHist=ROOT.TH1D("gammaEnergyHist",";Gamma energy (MeV);Counts",nBins*10,0,maxEnergy)
boundGammaEnergyHist=ROOT.TH1D("boundGammaEnergyHist",";Gamma energy (MeV);Counts",nBins*10,0,maxEnergy)
nNeutronsHist=ROOT.TH1D("nNeutronsHist",";Number of neutrons;Counts",10,0,10)
nGammasHist=ROOT.TH1D("nGammasHist",";Number of gammas;Counts",20,0,20)
alphaEnergyHist=ROOT.TH1D("alphaEnergyHist",";Alpha energy (MeV);Counts",nBins,0,maxEnergy)
protonEnergyHist=ROOT.TH1D("protonEnergyHist",";Proton energy (MeV);Counts",nBins,0,maxEnergy)
doubleNeutronEnergies2D = ROOT.TH2D("doubleNeutronEnergies2D",";Energy of first neutron (MeV);Energy of second neutron (MeV);Counts",150,0,15,150,0,15)
#Final particles
finalNucleiHist=ROOT.TH1D("finalNucleiHist",";PDG Code;Counts",1000000,0,1000000)
finalNuclearRecoilHist=ROOT.TH1D("finalNuclearRecoilHist",";Recoil Energy (keVnr);Counts",1000,0,1)
#Direction hists
leptonAngleHist=ROOT.TH1D("leptonAngleHist",";Outgoing electron angle (degrees);Counts",90,0,180)
leptonAngleNeutrinoEnergyHist=ROOT.TH2D("leptonAngleNeutrinoEnergyHist",";Outgoing electron angle (degrees);Neutrino energy (MeV)",90,0,180,53,0,53)
leptonAngleLeptonEnergyHist=ROOT.TH2D("leptonAngleLeptonEnergyHist",";Outgoing electron angle (degrees);Outgoing electron energy (MeV)",90,0,180,53,0,53)
leptonAngleLeptonEnergyHistRadians=ROOT.TH2D("leptonAngleLeptonEnergyHistRadians",";cos#theta;E_{#varepsilon} (MeV);Counts",20,-1,1,53,0,53)
leptonAngleExcitationHistRadians=ROOT.TH2D("leptonAngleExcitationHistRadians",";cos#theta;#omega (MeV);Counts",50,1,-1,53,0,53)
neutrinoAngleExcitationHist=ROOT.TH2D("neutrinoAngleExcitationHist","#theta (deg);#omega (MeV);Counts",45,0,180,53,0,53)
neutrinoAngleNeutrinoEnergyHist=ROOT.TH2D("neutrinoAngleNeutrinoEnergyHist",";#theta (deg);#Neutrino Energy (MeV);Counts",45,0,180,53,0,53)
neutrinoAngleHist=ROOT.TH1D("neutrinoAngleHist",";Outgoing neutrino angle (degrees);Counts",90,0,180)
singleNeutronAngleHist=ROOT.TH1D("singleNeutronAngleHist",";Neutron angle (degrees);Counts",90,0,180)
singleNeutronAngleNeutrinoEnergyHist=ROOT.TH2D("singleNeutronAngleNeutrinoEnergyHist",";Neutron angle (degrees);Neutrino energy (MeV)",90,0,180,53,0,53)
singleNeutronAngleNeutronEnergyHist=ROOT.TH2D("singleNeutronAngleNeutronEnergyHist",";Neutron a (degrees);Neutron energy (MeV)",90,0,180,10,0,10)
twoNeutronAngleHist=ROOT.TH1D("twoNeutronAngleHist",";Angle N_{1} v. N_{2}(degrees);Counts",90,0,180)
twoNeutronAngleNeutrinoEnergyHist=ROOT.TH2D("twoNeutronAngleNeutrinoEnergyHist",";Angle N_{1} v. N_{2} (degrees);Neutrino energy (MeV)",90,0,180,53,0,53)
twoNeutronAngleNeutronEnergyHist=ROOT.TH2D("twoNeutronAngleNeutronEnergyHist",";Angle N_{1} v. N_{2} (degrees);Neutron paireEnergy",90,0,180,20,0,20)
threeDirectionNeutronDalitzHist=ROOT.TH2D("threeDirectionNeutronDalitzHist",";Angle N_{1} v. N_{2}(degrees);Angle N_{2} v. N_{3}",90,0,180,90,0,180)

#TTree
marleyTree=ROOT.TTree("marleyTree","Basic event information")

neutrinoEnergy=array.array('d',[0])
finalNeutrinoEnergy=array.array('d',[0])
neutrinoOpeningAngle=array.array('d',[0])
excitationEnergy=array.array('d',[0])
leptonEnergy=array.array('d',[0])
gammaEnergy=array.array('d',[0])
leptonOpeningAngle=array.array('d',[0])
emEnergy=array.array('d',[0])
visibleEnergy=array.array('d',[0])
neutronEnergy=array.array('d',[0])
protonEnergy=array.array('d',[0])
alphaEnergy=array.array('d',[0])
nGammas=array.array('i',[0])
nNeutrons=array.array('i',[0])
nProtons=array.array('i',[0])
nAlphas=array.array('i',[0])
finalNuclearState=array.array('i',[0])
initialNuclearState=array.array('i',[0])
finalNuclearRecoilEnergy=array.array('d',[0])

marleyTree.Branch('neutrinoEnergy',neutrinoEnergy,'neutrinoEnergy/D')
marleyTree.Branch('finalNeutrinoEnergy',finalNeutrinoEnergy,'finalNeutrinoEnergy/D')
marleyTree.Branch('neutrinoOpeningAngle',neutrinoOpeningAngle,'neutrinoOpeningAngle/D')
marleyTree.Branch('excitationEnergy',excitationEnergy,'excitationEnergy/D')
marleyTree.Branch('leptonEnergy',leptonEnergy,'leptonEnergy/D')
marleyTree.Branch('gammaEnergy',gammaEnergy,'gammaEnergy/D')
marleyTree.Branch('leptonOpeningAngle',leptonOpeningAngle,'leptonOpeningAngle/D')
marleyTree.Branch('emEnergy',emEnergy,'emEnergy/D')
marleyTree.Branch('visibleEnergy',visibleEnergy,'visibleEnergy/D')
marleyTree.Branch('neutronEnergy',neutronEnergy,'neutronEnergy/D')
marleyTree.Branch('protonEnergy',protonEnergy,'protonEnergy/D')
marleyTree.Branch('alphaEnergy',alphaEnergy,'alphaEnergy/D')
marleyTree.Branch('nGammas',nGammas,'nGammas/I')
marleyTree.Branch('nNeutrons',nNeutrons,'nNeutrons/I')
marleyTree.Branch('nProtons',nProtons,'nProtons/I')
marleyTree.Branch('nAlphas',nAlphas,'nAlphas/I')
marleyTree.Branch('finalNuclearState',finalNuclearState,'finalNuclearState/I')
marleyTree.Branch('initialNuclearState',initialNuclearState,'initialNuclearState/I')
marleyTree.Branch('finalNuclearRecoilEnergy',finalNuclearRecoilEnergy,'finalNuclearRecoilEnergy/D')

#Get input file 
inpFile=open(sys.argv[1],"r")

#Keeps track of how many particles to look for in event 
nInitialParticles=0
nFinalParticles=0
nThrows=0
nNeutronEvents=0
nGammaEvents=0

#Loop through input file 
for i,line in enumerate(inpFile):
  if i%100000==0:
    print("On entry "+str(i))
  #print(line)
  #Strip newline char
  line=line.strip("\n")
  #Split by space
  lineParts=line.split(" ")

  #Cross section first
  if i==0 and oldFormat==0:
    crossSection=float(lineParts[0])*3.89379290*math.pow(10,-22)*math.pow(10,40)
  #New event
  elif nInitialParticles==0 and nFinalParticles==0:
  
    #Num initial particles, num final particles, excitation energy transferred to nucleus
    nInitialParticles=int(lineParts[0])
    nFinalParticles=int(lineParts[1])
    nThrows+=1
    
    excitationEnergy[0]=float(lineParts[2])
    excitationEnergyHist.Fill(excitationEnergy[0])
    
    #Reset quantities
    emEnergy[0]=0
    gammaEnergy[0]=0
    neutronEnergy[0]=0
    protonEnergy[0]=0
    alphaEnergy[0]=0
    finalNeutrinoEnergy[0]=0
    finalNuclearRecoilEnergy[0]=0
    finalNuclearState[0]=0
    initialNuclearState[0]=0
    nNeutrons[0]=0
    nGammas[0]=0
    nProtons[0]=0
    nAlphas[0]=0
    visibleEnergy[0]=0
    neutrinoOpeningAngle[0]=0
    
    #Reset list
    neutronEnergyList=[]
    neutronDirectionList=[]
    
  #See if we're still parsing the initial files 
  elif nInitialParticles>0:
  
    #Reduce # of initial particles by 1
    nInitialParticles-=1
  
    #Load up data on initial particles
    #pdf code, Four momentum 0-4, mass, charge
    pdgCode=int(lineParts[0]) #12=neutrino, 2112=neutron, 22=gamma, 11=electron
    totalEnergy=float(lineParts[1])
    px=float(lineParts[2])
    py=float(lineParts[3])
    pz=float(lineParts[4])
    mass=float(lineParts[5])
    charge=float(lineParts[6])
    
    if pdgCode>10000:
      initialNuclearState[0]=pdgCode
    
    #Get initial neutrino energy if this particle corresponds to the neutrino
    if abs(pdgCode)==12 or abs(pdgCode)==14 or abs(pdgCode)==16:
      neutrinoEnergy[0]=totalEnergy-restMass_nu
      neutrinoEnergyHist.Fill(neutrinoEnergy[0])
      
  #See if we're parsing the product particles 
  elif nFinalParticles>0:
  
    #Reduce # of final particles by 1
    nFinalParticles-=1

    #Read quantities
    totalEnergy=float(lineParts[1])
    px=float(lineParts[2])
    py=float(lineParts[3])
    pz=float(lineParts[4])
    mass=float(lineParts[5])
    charge=float(lineParts[6])
  
    #Get info on the particle
    #pdf code, Four momentum 0-4, mass, charge
    pdgCode=int(lineParts[0]) #12=neutrino, 2112=neutron, 22=gamma, 11=electron
    if pdgCode>1000000000 and pdgCode!=1000020040:
      finalNucleiHist.Fill(pdgCode-1000000000)
      finalNuclearState[0]=(pdgCode-1000000000)
      finalNuclearRecoilEnergy[0]=totalEnergy-mass
      finalNuclearRecoilHist.Fill(finalNuclearRecoilEnergy[0])
      #quenchedEnergy=getQuenchedEnergy(finalNuclearRecoilEnergy[0])
      #visibleEnergy[0]+=quenchedEnergy
      
    #Leptons
    if abs(pdgCode)==11:
      energy=totalEnergy-restMass_e
      leptonDirection=(px,py,pz)
      leptonOpeningAngle[0]=angle_between(leptonDirection,neutrinoDirection)*180./math.pi
      
      leptonEnergyHist.Fill(energy)
      leptonAngleHist.Fill(leptonOpeningAngle[0])
      leptonAngleNeutrinoEnergyHist.Fill(leptonOpeningAngle[0],neutrinoEnergy[0])
      leptonAngleLeptonEnergyHist.Fill(leptonOpeningAngle[0],energy)
      leptonAngleLeptonEnergyHistRadians.Fill(numpy.cos(leptonOpeningAngle[0]*math.pi/180.),energy)
      leptonAngleExcitationHistRadians.Fill(numpy.cos(leptonOpeningAngle[0]*math.pi/180.),excitationEnergy[0])
      
      leptonEnergy[0]=energy
      emEnergy[0]+=energy
      visibleEnergy[0]+=energy
    #Neutrons
    elif pdgCode==2112:
      energy=totalEnergy-restMass_n
      neutronDirection=(px,py,pz)
      
      neutronEnergyList.append(energy)
      
      neutronEnergy[0]+=energy
      neutronDirectionList.append(neutronDirection)
      nNeutrons[0]+=1
      #quenchedEnergy=getQuenchedEnergy(energy)
      #visibleEnergy[0]+=quenchedEnergy
    #Gammas
    elif pdgCode==22:
      energy=totalEnergy
      
      gammaEnergyHist.Fill(totalEnergy)
      emEnergy[0]+=energy
      gammaEnergy[0]+=energy
      visibleEnergy[0]+=energy
      nGammas[0]+=1
    #Alpha
    elif pdgCode==1000020040:
      energy=totalEnergy-restMass_alpha
      alphaEnergy[0]+=energy
      nAlphas[0]+=1
      alphaEnergyHist.Fill(energy)
      #quenchedEnergy=getQuenchedEnergy(energy)
      #visibleEnergy[0]+=quenchedEnergy
    #Proton
    elif pdgCode==2212:
      energy=totalEnergy-restMass_p
      protonEnergyHist.Fill(energy)
      protonEnergy[0]+=energy
      nProtons[0]+=1
      #quenchedEnergy=getQuenchedEnergy(energy)
      #visibleEnergy[0]+=quenchedEnergy
    elif pdgCode==12 or pdgCode==14 or pdgCode==16:
      finalNeutrinoEnergyHist.Fill(totalEnergy-restMass_nu)
      finalNeutrinoEnergy[0]=totalEnergy-restMass_nu
      neutrinoOutgoingDirection=(px,py,pz)
      neutrinoOpeningAngle[0]=angle_between(neutrinoOutgoingDirection,neutrinoDirection)*180./math.pi
    
      
  #See if this is the last final particle. If so, we can update hists
  if nFinalParticles==0:
    if nProtons[0]==0 and nNeutrons[0]==0 and nAlphas[0]==0:
      boundGammaEnergyHist.Fill(gammaEnergy[0])
      
      neutrinoAngleHist.Fill(neutrinoOpeningAngle[0])
      neutrinoAngleExcitationHist.Fill(neutrinoOpeningAngle[0],excitationEnergy[0])
      neutrinoAngleNeutrinoEnergyHist.Fill(neutrinoOpeningAngle[0],neutrinoEnergy[0])
      
      
    marleyTree.Fill()
    finalNuclearRecoilEnergy[0]=0
    finalNuclearState[0]=0
    
    emEnergyHist.Fill(emEnergy[0])
    if nNeutrons[0]>0:
      for nrg in neutronEnergyList:
        neutronEnergyHist.Fill(nrg)
      nNeutronEvents+=1
    nNeutronsHist.Fill(nNeutrons[0])
    if nGammas[0]>0:
      nGammasHist.Fill(nGammas[0])
      nGammaEvents+=1
    
    #1n events
    if nNeutrons[0]==1:
      singleNeutronEnergyHist.Fill(neutronEnergyList[0])
      
      singleNeutronAngle=angle_between(neutronDirectionList[0],neutrinoDirection)*180./math.pi
      
      singleNeutronAngleHist.Fill(singleNeutronAngle)
      singleNeutronAngleNeutrinoEnergyHist.Fill(singleNeutronAngle,neutrinoEnergy[0])
      singleNeutronAngleNeutronEnergyHist.Fill(singleNeutronAngle,neutronEnergyList[0])
      
      leptonEnergyOneNeutronHist.Fill(leptonEnergy[0])
    #2n events
    elif nNeutrons[0]==2:
      for nrg in neutronEnergyList:
        doubleNeutronEnergyHist.Fill(nrg)
      doubleNeutronEnergies2D.Fill(neutronEnergyList[0],neutronEnergyList[1])
      
      twoNeutronAngle=angle_between(neutronDirectionList[0],neutronDirectionList[1])*180./math.pi
      
      twoNeutronAngleHist.Fill(twoNeutronAngle)
      twoNeutronAngleNeutrinoEnergyHist.Fill(twoNeutronAngle,neutrinoEnergy[0])
      twoNeutronAngleNeutronEnergyHist.Fill(twoNeutronAngle,neutronEnergyList[0]+neutronEnergyList[1])
      
      leptonEnergyTwoNeutronHist.Fill(leptonEnergy[0])
    #3n events
    elif nNeutrons[0]==3:
      for nrg in neutronEnergyList:
        tripleNeutronEnergyHist.Fill(nrg)
      angle_12 = angle_between(neutronDirectionList[0],neutronDirectionList[1])*180./math.pi
      angle_23 = angle_between(neutronDirectionList[1],neutronDirectionList[2])*180./math.pi
      threeDirectionNeutronDalitzHist.Fill(angle_12,angle_23)
      
      leptonEnergyThreeNeutronHist.Fill(leptonEnergy[0])
    elif nNeutrons[0]==4:
      for nrg in neutronEnergyList:
        quadrupleNeutronEnergyHist.Fill(nrg)
    
#Close the input file
inpFile.close()

#Create the output file
outFile=ROOT.TFile(sys.argv[2],"RECREATE")

#Write the TTRee
marleyTree.Write()

#Write the histograms
#Incident neutrino hists
neutrinoEnergyHist.Write()
#Nuclear excitation hist
excitationEnergyHist.Write()
#Final nucleus hists
finalNucleiHist.Write()
#Outgoing particle energy hists
leptonEnergyHist.Write()
finalNeutrinoEnergyHist.Write()
emEnergyHist.Write()
neutronEnergyHist.Write()
singleNeutronEnergyHist.Write()
doubleNeutronEnergyHist.Write()
tripleNeutronEnergyHist.Write()
quadrupleNeutronEnergyHist.Write()
gammaEnergyHist.Write()
boundGammaEnergyHist.Write()
alphaEnergyHist.Write()
protonEnergyHist.Write()
nNeutronsHist.Write()
nGammasHist.Write()
doubleNeutronEnergies2D.Write()
#Outgoing angle hists
leptonAngleHist.Write()
leptonAngleNeutrinoEnergyHist.Write()
leptonAngleLeptonEnergyHist.Write()
singleNeutronAngleHist.Write()
singleNeutronAngleNeutrinoEnergyHist.Write()
singleNeutronAngleNeutronEnergyHist.Write()
twoNeutronAngleHist.Write()
twoNeutronAngleNeutrinoEnergyHist.Write()
twoNeutronAngleNeutronEnergyHist.Write()
threeDirectionNeutronDalitzHist.Write()



#NORMALIZED CROSS SECTION HISTS
if oldFormat==0:
  leptonCrossSection=1*crossSection
  normalizedLeptonEnergyHist=leptonEnergyHist.Clone("normalizedLeptonEnergyHist")
  if normalizedLeptonEnergyHist.GetEntries()>0:
    normalizedLeptonEnergyHist.Scale(leptonCrossSection/normalizedLeptonEnergyHist.GetEntries())
    normalizedLeptonEnergyHist.GetYaxis().SetTitle("d#sigma/dE_{e} (x 10^{-40}cm^{2})")
    normalizedLeptonEnergyHist.Write()

  emCrossSection=1*crossSection
  normalizedEmEnergyHist=emEnergyHist.Clone("normalizedEmEnergyHist")
  if normalizedEmEnergyHist.GetEntries()>0:
    normalizedEmEnergyHist.Scale(emCrossSection/normalizedEmEnergyHist.GetEntries())
    normalizedEmEnergyHist.GetYaxis().SetTitle("d#sigma/dE_{em} (x 10^{-40}cm^{2})")
    normalizedEmEnergyHist.Write()

  gammaCrossSection=nGammaEvents*1./nThrows*crossSection
  normalizedGammaEnergyHist=gammaEnergyHist.Clone("normalizedGammaEnergyHist")
  if normalizedGammaEnergyHist.GetEntries()>0:
    normalizedGammaEnergyHist.Scale(gammaCrossSection/normalizedGammaEnergyHist.GetEntries())
    normalizedGammaEnergyHist.GetYaxis().SetTitle("d#sigma/dE_{#gamma} (x 10^{-40}cm^{2})")
    normalizedGammaEnergyHist.Write()

  singleNeutronCrossSection=singleNeutronEnergyHist.GetEntries()*1./nThrows*crossSection
  normalizedSingleNeutronEnergyHist=singleNeutronEnergyHist.Clone("normalizedSingleNeutronEnergyHist")
  if normalizedSingleNeutronEnergyHist.GetEntries()>0:
    normalizedSingleNeutronEnergyHist.Scale(singleNeutronCrossSection/normalizedSingleNeutronEnergyHist.GetEntries())
    normalizedSingleNeutronEnergyHist.GetYaxis().SetTitle("d#sigma/dE_{n} (x 10^{-40}cm^{2})")
    normalizedSingleNeutronEnergyHist.Write()

  doubleNeutronCrossSection=doubleNeutronEnergyHist.GetEntries()*1./nThrows*crossSection*1./2.
  normalizedDoubleNeutronEnergyHist=doubleNeutronEnergyHist.Clone("normalizedDoubleNeutronEnergyHist")
  if doubleNeutronEnergyHist.GetEntries()>0:
    normalizedDoubleNeutronEnergyHist.Scale(doubleNeutronCrossSection/doubleNeutronEnergyHist.GetEntries())
    normalizedDoubleNeutronEnergyHist.GetYaxis().SetTitle("d#sigma/dE_{n} (x 10^{-40}cm^{2})")
    normalizedDoubleNeutronEnergyHist.Write()

  tripleNeutronCrossSection=tripleNeutronEnergyHist.GetEntries()*1./nThrows*crossSection*1./3.
  normalizedTripleNeutronEnergyHist=tripleNeutronEnergyHist.Clone("normalizedTripleNeutronEnergyHist")
  if tripleNeutronEnergyHist.GetEntries()>0:
    normalizedTripleNeutronEnergyHist.Scale(tripleNeutronCrossSection/tripleNeutronEnergyHist.GetEntries())
    normalizedTripleNeutronEnergyHist.GetYaxis().SetTitle("d#sigma/dE_{n} (x 10^{-40}cm^{2})")
    normalizedTripleNeutronEnergyHist.Write()

  quadrupleNeutronCrossSection=quadrupleNeutronEnergyHist.GetEntries()*1./nThrows*crossSection*1./4.
  normalizedQuadrupleNeutronEnergyHist=quadrupleNeutronEnergyHist.Clone("normalizedQuadrupleNeutronEnergyHist")
  if quadrupleNeutronEnergyHist.GetEntries()>0:
    normalizedQuadrupleNeutronEnergyHist.Scale(quadrupleNeutronCrossSection/quadrupleNeutronEnergyHist.GetEntries())
    normalizedQuadrupleNeutronEnergyHist.GetYaxis().SetTitle("d#sigma/dE_{n} (x 10^{-40}cm^{2})")
    normalizedQuadrupleNeutronEnergyHist.Write()
    
  neutronCrossSection=singleNeutronCrossSection+doubleNeutronCrossSection+tripleNeutronCrossSection+quadrupleNeutronCrossSection
  normalizedNeutronEnergyHist=neutronEnergyHist.Clone("normalizedNeutronEnergyHist")
  if normalizedNeutronEnergyHist.GetEntries()>0:
    normalizedNeutronEnergyHist.Scale(neutronCrossSection/normalizedNeutronEnergyHist.GetEntries())
    normalizedNeutronEnergyHist.GetYaxis().SetTitle("d#sigma/dE_{n} (x 10^{-40}cm^{2})")
    normalizedNeutronEnergyHist.Write()


  alphaCrossSection=alphaEnergyHist.GetEntries()*1./nThrows*crossSection
  normalizedAlphaEnergyHist=alphaEnergyHist.Clone("normalizedAlphaEnergyHist")
  if normalizedAlphaEnergyHist.GetEntries()>0:
    normalizedAlphaEnergyHist.Scale(alphaCrossSection/alphaEnergyHist.GetEntries())
    normalizedAlphaEnergyHist.GetYaxis().SetTitle("d#sigma/dE_{#alpha} (x 10^{-40}cm^{2})")
    normalizedAlphaEnergyHist.Write()

  protonCrossSection=protonEnergyHist.GetEntries()*1./nThrows*crossSection
  normalizedProtonEnergyHist=protonEnergyHist.Clone("normalizedProtonEnergyHist")
  if protonEnergyHist.GetEntries()>0:
    normalizedProtonEnergyHist.Scale(protonCrossSection/protonEnergyHist.GetEntries())
    normalizedProtonEnergyHist.GetYaxis().SetTitle("d#sigma/dE_{p} (x 10^{-40}cm^{2})")
    normalizedProtonEnergyHist.Write()

  normalizedLeptonAngleLeptonEnergyHistRadians=leptonAngleLeptonEnergyHistRadians.Clone("normalizedLeptonAngleLeptonEnergyHistRadians")
  if normalizedLeptonAngleLeptonEnergyHistRadians.GetEntries()>0:
    normalizedLeptonAngleLeptonEnergyHistRadians.Scale(leptonCrossSection/normalizedLeptonAngleLeptonEnergyHistRadians.GetEntries())
    normalizedLeptonAngleLeptonEnergyHistRadians.Scale(100.) #So 10^-42 not 10^-40
    normalizedLeptonAngleLeptonEnergyHistRadians.GetZaxis().SetTitle("d#sigma/d(cos#theta)dE (x 10^{-42]cm^{2}/MeV #cdot 0.157 rad)")
    normalizedLeptonAngleLeptonEnergyHistRadians.Write()
  normalizedLeptonAngleExcitationEnergyRadiansHist=leptonAngleExcitationHistRadians.Clone("normalizedLeptonAngleExcitationEnergyRadiansHist")
  if normalizedLeptonAngleExcitationEnergyRadiansHist.GetEntries()>0:
    normalizedLeptonAngleExcitationEnergyRadiansHist.Scale(leptonCrossSection/normalizedLeptonAngleExcitationEnergyRadiansHist.GetEntries())
    normalizedLeptonAngleExcitationEnergyRadiansHist.Scale(100.) #So 10^-42 not 10^-40
    normalizedLeptonAngleExcitationEnergyRadiansHist.GetZaxis().SetTitle("d#sigma/d(cos#theta)dE (x 10^{-42]cm^{2}/MeV)")
    normalizedLeptonAngleExcitationEnergyRadiansHist.Write()

  normalizedFinalNucleiHist=finalNucleiHist.Clone("normalizedFinalNucleiHist")
  if normalizedFinalNucleiHist.GetEntries()>0:
    normalizedFinalNucleiHist.Scale(leptonCrossSection/normalizedFinalNucleiHist.GetEntries())
    normalizedFinalNucleiHist.GetYaxis().SetTitle("Cross Section (x 10^{-40}cm^{2})")
    normalizedFinalNucleiHist.Write()
    
  normalizedNNeutronsHist=nNeutronsHist.Clone("normalizedNNeutronsHist")
  if normalizedNNeutronsHist.GetEntries()>0:
    normalizedNNeutronsHist.Scale(crossSection/normalizedNNeutronsHist.GetEntries())
    normalizedNNeutronsHist.GetYaxis().SetTitle("Cross Section (x 10^{-40}cm^{2})")
    normalizedNNeutronsHist.Write()
      
  normalizedFinalNeutrinoEnergyHist=finalNeutrinoEnergyHist.Clone("normalizedFinalNeutrinoEnergyHist")
  if normalizedFinalNeutrinoEnergyHist.GetEntries()>0:
    normalizedFinalNeutrinoEnergyHist.Scale(crossSection/finalNeutrinoEnergyHist.GetEntries())
    normalizedFinalNeutrinoEnergyHist.GetYaxis().SetTitle("Cross Section (x 10^{-40}cm^{2})")
    normalizedFinalNeutrinoEnergyHist.Write()
    
  normalizedExcitationEnergyHist=excitationEnergyHist.Clone("normalizedExcitationEnergyHist")
  if normalizedExcitationEnergyHist.GetEntries()>0:
    normalizedExcitationEnergyHist.Scale(crossSection/excitationEnergyHist.GetEntries())
    normalizedExcitationEnergyHist.GetYaxis().SetTitle("Cross Section (x 10^{-40}cm^{2})")
    normalizedExcitationEnergyHist.Write()
  
  normalizedFinalNuclearRecoilHist=finalNuclearRecoilHist.Clone("normalizedFinalNuclearRecoilHist")
  if normalizedFinalNuclearRecoilHist.GetEntries()>0:
    normalizedFinalNuclearRecoilHist.Scale(crossSection/normalizedFinalNuclearRecoilHist.GetEntries())
    normalizedFinalNuclearRecoilHist.GetYaxis().SetTitle("Cross Section (x 10^{-40}cm^{2})")
    normalizedFinalNuclearRecoilHist.Write()
    
  normalizedLeptonEnergyOneNeutronHist=leptonEnergyOneNeutronHist.Clone("normalizedLeptonEnergyOneNeutronHist")
  if normalizedLeptonEnergyOneNeutronHist.GetEntries()>0:
    normalizedLeptonEnergyOneNeutronHist.Scale(singleNeutronCrossSection/normalizedLeptonEnergyOneNeutronHist.GetEntries())
    normalizedLeptonEnergyOneNeutronHist.GetYaxis().SetTitle("Cross Section(x 10^{-40}cm^{2})")
    normalizedLeptonEnergyOneNeutronHist.Write()
    
  normalizedLeptonEnergyTwoNeutronHist=leptonEnergyTwoNeutronHist.Clone("normalizedLeptonEnergyTwoNeutronHist")
  if normalizedLeptonEnergyTwoNeutronHist.GetEntries()>0:
    normalizedLeptonEnergyTwoNeutronHist.Scale(doubleNeutronCrossSection/normalizedLeptonEnergyTwoNeutronHist.GetEntries())
    normalizedLeptonEnergyTwoNeutronHist.GetYaxis().SetTitle("Cross Section(x 10^{-40}cm^{2})")
    normalizedLeptonEnergyTwoNeutronHist.Write()
    
  normalizedLeptonEnergyTwoNeutronHist=leptonEnergyTwoNeutronHist.Clone("normalizedLeptonEnergyTwoNeutronHist")
  if normalizedLeptonEnergyTwoNeutronHist.GetEntries()>0:
    normalizedLeptonEnergyTwoNeutronHist.Scale(doubleNeutronCrossSection/normalizedLeptonEnergyTwoNeutronHist.GetEntries())
    normalizedLeptonEnergyTwoNeutronHist.GetYaxis().SetTitle("Cross Section(x 10^{-40}cm^{2})")
    normalizedLeptonEnergyTwoNeutronHist.Write()
    
  normalizedLeptonEnergyThreeNeutronHist=leptonEnergyThreeNeutronHist.Clone("normalizedLeptonEnergyThreeNeutronHist")
  if normalizedLeptonEnergyThreeNeutronHist.GetEntries()>0:
    normalizedLeptonEnergyThreeNeutronHist.Scale(tripleNeutronCrossSection/normalizedLeptonEnergyThreeNeutronHist.GetEntries())
    normalizedLeptonEnergyThreeNeutronHist.GetYaxis().SetTitle("Cross Section(x 10^{-40}cm^{2})")
    normalizedLeptonEnergyThreeNeutronHist.Write()
    
  boundGammaCrossSection=boundGammaEnergyHist.GetEntries()*1./nThrows*crossSection
  normalizedBoundGammaEnergyHist=boundGammaEnergyHist.Clone("normalizedBoundGammaEnergyHist")
  if normalizedBoundGammaEnergyHist.GetEntries()>0:
    normalizedBoundGammaEnergyHist.Scale(boundGammaCrossSection/normalizedBoundGammaEnergyHist.GetEntries())
    normalizedBoundGammaEnergyHist.GetYaxis().SetTitle("Cross Section(x 10^{-40}cm^{2})")
    normalizedBoundGammaEnergyHist.Write()
    
  normalizedNeutrinoAngleHist=neutrinoAngleHist.Clone("normalizedNeutrinoAngleHist")
  if normalizedNeutrinoAngleHist.GetEntries()>0:
    normalizedNeutrinoAngleHist.Scale(crossSection/normalizedNeutrinoAngleHist.GetEntries())
    normalizedNeutrinoAngleHist.GetYaxis().SetTitle("Cross Section(x 10^{-40}cm^{2})")
    normalizedNeutrinoAngleHist.Write()
   
  normalizedNeutrinoAngleExcitationHist=neutrinoAngleExcitationHist.Clone("normalizedNeutrinoAngleExcitationHist")
  if normalizedNeutrinoAngleExcitationHist.GetEntries()>0:
    normalizedNeutrinoAngleExcitationHist.Scale(crossSection/normalizedNeutrinoAngleExcitationHist.GetEntries())
    normalizedNeutrinoAngleExcitationHist.GetZaxis().SetTitle("Cross Section(x 10^{-42}cm^{2})")
    normalizedNeutrinoAngleExcitationHist.Scale(100) #so 10^-42
    normalizedNeutrinoAngleExcitationHist.Write()
    
  normalizedNeutrinoAngleNeutrinoEnergyHist=neutrinoAngleNeutrinoEnergyHist.Clone("normalizedNeutrinoAngleNeutrinoEnergyHist")
  if normalizedNeutrinoAngleNeutrinoEnergyHist.GetEntries()>0:
    normalizedNeutrinoAngleNeutrinoEnergyHist.Scale(crossSection/normalizedNeutrinoAngleNeutrinoEnergyHist.GetEntries())
    normalizedNeutrinoAngleNeutrinoEnergyHist.GetZaxis().SetTitle("Cross Section(x 10^{-42}cm^{2})")
    normalizedNeutrinoAngleNeutrinoEnergyHist.Scale(100) #so 10^-42
    normalizedNeutrinoAngleNeutrinoEnergyHist.Write() 
       
    
    
    
  nDigits=3
  print('\nTotal bound gamma cross section: {0:.{1}e} x 10^-42 cm^2'.format(boundGammaCrossSection*100, nDigits))
  print('Total single neutron cross section: {0:.{1}e} x 10^-42 cm^2'.format(singleNeutronCrossSection*100, nDigits))
  print('Total proton cross section: {0:.{1}e} x 10^-42 cm^2'.format(protonCrossSection*100, nDigits))
  print('Total alpha cross section: {0:.{1}e} x 10^-42 cm^2'.format(alphaCrossSection*100, nDigits))
  print('Total cross section: {0:.{1}e} x 10^-42 cm^2\n'.format(crossSection*100, nDigits))
  
  
  print('\nTotal single neutron cross section: {0:.{1}e} x 10^-40 cm^2'.format(singleNeutronCrossSection*100, nDigits))
  print('Total double neutron cross section: {0:.{1}e} x 10^-40 cm^2\n'.format(doubleNeutronCrossSection*100, nDigits))
  print('Total triple neutron cross section: {0:.{1}e} x 10^-40 cm^2\n'.format(tripleNeutronCrossSection*100, nDigits))
  print('Total quadruple neutron cross section: {0:.{1}e} x 10^-40 cm^2\n'.format(quadrupleNeutronCrossSection*100, nDigits))
  '''
  print("Total lepton cross section: "+str(leptonCrossSection)+" x 10^-40 cm^2")
  print("Total gamma cross section: "+str(gammaCrossSection)+" x 10^-40 cm^2")
  print("Total neutron cross section: "+str(neutronCrossSection)+" x 10^-40 cm^2")
  print("Total double neutron cross section: "+str(doubleNeutronCrossSection)+" x 10^-40 cm^2")
  print("Total triple neutron cross section: "+str(tripleNeutronCrossSection)+" x 10^-40 cm^2")
  print("Total quadruple neutron cross section: "+str(quadrupleNeutronCrossSection)+" x 10^-40 cm^2")
  '''
#Close the output file
outFile.Close()