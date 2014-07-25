
File {
	Grid = "n91_msh.tdr"
	Plot    = "n391_des.tdr"	*Commented out to limit memory usage on system.  Add back in for individual experiments, not long optimization simulations.
 	Current = "n391_des.plt"
	Output  = "n391_des.log"
	Parameter = "pp391_des.par"
	IlluminationSpectrum = "n135_am15d.txt"
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
  	
	} 

		
Physics	(Region="Al.2Ga.8Sb_top_window"){
	MoleFraction(xFraction = .83)}
		
Physics	(Region = "Al.2Ga.8Sb_region1") {
	MoleFraction(xFraction =.83)}
	
Physics	(Region = "Al.2Ga.8Sb_region2") {
	MoleFraction(xFraction =.83)}
	
Physics	(Region="Al.2Ga.8Sb_top_BSF"){
	MoleFraction(xFraction =.8)}
	
Physics	(Region = "InAlAs_region3") {
	MoleFraction(xFraction =.1)}
	
Physics	(Region = "InAlAs_region4") {
	MoleFraction(xFraction =.3)}
	

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
 	 Digits=5              
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

          }
          
            


Solve   {
	 Coupled (Iterations=20) {Poisson}
  
	 * Use Transient to help achieve initial convergence
	 Transient (
	   Initialtime=0 FinalTime=0.1 
	   InitialStep=1e-6 MinStep=1e-10
	 ){ Coupled {Poisson Electron Hole } }
	 * Plot (FilePrefix = "n391_sc")
    	
    	*CurrentPlot()
     	*System("rm -f n391_banddgm_des.tdr")
     	
     	Quasistationary ( 
   		 InitialStep=5e-2 MaxStep =1e-2 MinStep = 1e-3 
    		Goal{ voltage = .5 Name = "front" }
 	   ){ Coupled {Poisson Electron Hole  } } 
    
  	 Quasistationary ( 
   		 InitialStep=5e-2 MaxStep =1e-3 MinStep = 1e-3 
    		Goal{ voltage = 2.3 Name = "front" }
 	   ){ Coupled {Poisson Electron Hole  } }  
            }





