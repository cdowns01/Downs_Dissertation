
File {
	Grid = "n71_msh.tdr"
	Plot    = "n76_des.tdr"	*Commented out to limit memory usage on system.  Add back in for individual experiments, not long optimization simulations.
 	Current = "n76_des.plt"
	Output  = "n76_des.log"
	Parameter = "pp76_des.par"
	IlluminationSpectrum = "n73_am15d.txt"
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
      	    	
  	eBarrierTunneling "TD_NLM2" (
    	  Band2Band
   	  TwoBand
 	 ) 

  	hBarrierTunneling "TD_NLM2"(
      	Band2Band
      	TwoBand
  	)
  	
	} 

Physics	(Region = "Al.2Ga.8Sb_mid_window") {
	MoleFraction(xFraction =.8)}
		
Physics	(Region = "Al.2Ga.8Sb_region5") {
	MoleFraction(xFraction =.8)}
	
Physics	(Region = "Al.2Ga.8Sb_region6") {
	MoleFraction(xFraction =.8)}
	
Physics	(Region = "Al.2Ga.8Sb_mid_BSF") {
	MoleFraction(xFraction =.8)}
	
Physics	(Region = "InAlAs_region7") {
	MoleFraction(xFraction =.2)}
	
*Physics	(Region = "InAlAs_region8") {
	*MoleFraction(xFraction = .375)}

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
 	 Digits=9              
 	Iterations=20  
	Notdamped=100 
	ExtendedPrecision
  	
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

	  
	   NonLocal "TD_NLM2" (
   	 RegionInterface = "InAlAs_region7/InAlAs_region8"
  	 Length=15e-7         
 	 Permeation = 15e-7
	  )
          }
          
            


Solve   {
	 Coupled (Iterations=20) {Poisson}
  
	 * Use Transient to help achieve initial convergence
	Transient (
	   Initialtime=0 FinalTime=0.1 
	   InitialStep=1e-2 MinStep=1e-7
	 ){ Coupled {Poisson Electron Hole } }


		* fast voltage controlled ramping until mpp
		Quasistationary ( 
		  InitialStep=1e-2 MaxStep =1e-3 MinStep = 1e-2 Increment=1.4 DoZero
		  Goal{ voltage = 1.5 Name="front" }
		){ Coupled {Poisson Electron Hole  } 
		 }  
 
            }






