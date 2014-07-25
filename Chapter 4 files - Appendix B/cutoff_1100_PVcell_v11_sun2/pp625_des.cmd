
File {
	Grid = "n610_msh.tdr"
	Plot    = "n625_des.tdr"	*Commented out to limit memory usage on system.  Add back in for individual experiments, not long optimization simulations.
 	Current = "n625_des.plt"
	Output  = "n625_des.log"
	Parameter = "pp625_des.par"
	IlluminationSpectrum = "n621_am15d.txt"
      }
	
Electrode {
	  {Name = "front"  Voltage = 0 }
	  {Name = "back" Voltage = 0}
                }
                     
Physics {
  	AreaFactor = 2000000000.0
  	Fermi
  	HeteroInterface

  	EffectiveIntrinsicDensity(NoBandGapNarrowing)
  	
  	Recombination(
   	SRH
    	Auger
    	Radiative
  	)
  	Mobility(
   	 DopingDep
  	)

 	 ThermionicEmission
 	  
  	OpticalGeneration(
   	 QuantumYield= 1 * generated carriers/photon, default: 1, for lambda < 450, > 1
    	ComputeFromSpectrum
  	) * end OpticalGeneration
  	
  	Optics (
   	     TMM (
    	      NodesPerWavelength = 100
     	     Excitation (
      	        Theta = 0
       	       Polarization = 0
        	  )

         	 Stripe (
          	   Left = -24.0
          	   Right = 25.0
          	           )  *end stripe
       	          )  *end TMM
           	     ) * end Optics  
      	  	
  	eBarrierTunneling "TD_NLM1" (
    	  Band2Band
   	  TwoBand
 	 ) 

  	hBarrierTunneling "TD_NLM1"(
      	Band2Band
      	TwoBand
  	)
  	
  	eBarrierTunneling "TD_NLM2" (
    	  Band2Band
   	  TwoBand
 	 ) 

  	hBarrierTunneling "TD_NLM2"(
      	Band2Band
      	TwoBand
  	)
  	
	} 

		
Physics	(Region="GaInP_top_window"){
	MoleFraction(xFraction = .8)}
		
Physics	(Region = "GaInP_region1") {
	MoleFraction(xFraction =.75)}
	
Physics	(Region = "GaInP_region2") {
	MoleFraction(xFraction =.75)}
	
Physics	(Region="GaInP_top_BSF"){
	MoleFraction(xFraction =.8)}
	
Physics	(Region = "InAlAs_region3") {
	MoleFraction(xFraction =.4)}
	
Physics	(Region = "InAlAs_region4") {
	MoleFraction(xFraction =.3)}
	
Physics	(Region = "AlGaAs_mid_window") {
	MoleFraction(xFraction =.1)}
		
Physics	(Region = "AlGaAs_region5") {
	MoleFraction(xFraction =.06)}
	
Physics	(Region = "AlGaAs_region6") {
	MoleFraction(xFraction =.06)}
	
Physics	(Region = "AlGaAs_mid_BSF") {
	MoleFraction(xFraction = .15)}
	
Physics	(Region = "AlInP_region7") {
	MoleFraction(xFraction =.06)}
	
Physics	(Region = "InGaAs_region8") {
	MoleFraction(xFraction = .3)}
	

Physics	(Region = "InGaAs_bot_window") {
	MoleFraction(xFraction =.8)}
		
Physics	(Region = "InGaAs_region9") {
	MoleFraction(xFraction =.80)}*.775
	
Physics	(Region = "InGaAs_region10") {
	MoleFraction(xFraction =.80)}
	
Physics	(Region = "InGaAs_bot_BSF") {
	MoleFraction(xFraction =.8)}

Plot {
	eCurrent 
	hCurrent
	Potential 
	ElectricField
	eMobility 
	hMobility 
	eVelocity 
	hVelocity
	OpticalIntensity
	OpticalGeneration
	OpticalAbsorption
	OpticalField
	NDopantConcentration
	PDopantConcentration
	xMoleFraction Doping DonorConcentration AcceptorConcentration 
 	eEffectiveStateDensity hEffectiveStateDensity EffectiveIntrinsicDensity IntrinsicDensity
  	eDensity hDensity SpaceCharge
  	eQuasiFermiPotential hQuasiFermiPotential BandGap ConductionBandEnergy ValenceBandEnergy ElectronAffinity
  	ElectricField ElectricField/vector ElectrostaticPotential
  	eLifetime hLifetime SRH Auger TotalRecombination SurfaceRecombination RadiativeRecombination
  	eCurrent/Vector hCurrent/Vector current/vector
  	eMobility hMobility eVelocity hVelocity
  	SRH Auger TotalRecombination SurfaceRecombination RadiativeRecombination
  	RefractiveIndexRealVertexComponent0
  	RefractiveIndexImagVertexComponent0
	NonLocal   
        }
                
Math  {
	 Extrapolate 
 	 Derivatives   
 	 RelErrControl 
 	 Digits=11  
 	 ExtendedPrecision            
 	Iterations=20  
	Notdamped=100 
  	
	Method = Super
	ExitOnFailure
	Transient = BE
	RhsMax = 1e20
	
  	Number_of_Threads = maximum
	
	ErrRef(electron) = 1E7
  	ErrRef(hole) = 1E7
	
	BreakCriteria{
	Current (Contact = "front" maxval = .00001)
	}
	  
	  NonLocal "TD_NLM1" (
   	 RegionInterface = "InAlAs_region3/InAlAs_region4"
  	 Length=15e-7         
 	 Permeation = 15e-7
	  )
	  
	  NonLocal "TD_NLM2" (
   	 RegionInterface = "AlInP_region7/InGaAs_region8"
  	 Length=15e-7         
 	 Permeation = 15e-7
	  )
          }
          
            


Solve   {
	 Coupled (Iterations=20) {Poisson}
  
	 * Use Transient to help achieve initial convergence
	 Transient (
	   Initialtime=0 FinalTime=0.1 
	   InitialStep=1e-10 MinStep=1e-14
	 ){ Coupled {Poisson Electron Hole } }
	 * Plot (FilePrefix = "n625_sc")
    	
    	*CurrentPlot()
     	*System("rm -f n625_banddgm_des.tdr")
     	
     	Quasistationary ( 
   		 InitialStep=5e-2 MaxStep =1e-1 MinStep = 1e-3 
    		Goal{ voltage = 2.5 Name = "front" }
 	   ){ Coupled {Poisson Electron Hole  } } 
    
    
    	Plot (FilePrefix = "n625_debug")
    	
  	 Quasistationary ( 
   		 InitialStep=5e-10 MaxStep =1e-2 MinStep = 1e-10 
    		Goal{ voltage = 4.5 Name = "front" }
 	   ){ Coupled {Poisson Electron Hole  } }  
            }








