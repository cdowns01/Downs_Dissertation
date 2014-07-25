#!MC 1100
#setdep @previous|-1@

## Delete the top frame
$!FRAMECONTROL DELETETOP

## Load data in Tecplot using the SWB-Loader
$!READDATASET "n@node|-2@_des.tdr" DATASETREADER = "SWB-Loader"  ##add node number for sdevice to plot appropriate response

## Fit the entire plot to the grid area
$!VIEW FIT

##---------------------------------------------------------------------------
## Tecplot macro for plotting the band diagram
## Function: Makes a plot of the banddiagram along a 1D cut perpendicular to X-axis at position |CutPos| 
## Created by Gergoe Letay, Synopsys LLC 2005
## Modified by Sameer Shah on 26/02/08
##---------------------------------------------------------------------------

## Turning off all graphics during preliminary portions of a macro file can greatly increase the efficiency of the macro.
$!DRAWGRAPHICS FALSE

## Specify axis for normal cut
## i.e. X => X-normal cut
$!VARSET |CutAxis| = 'X'

## Specify coordinate for cut
$!VARSET |CutPos| = @cutpos@

##------- Frame 2 ------------
## The orthogonal slicer merges zones
$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
	COMMAND = 'SLICEMODE MergedZones' 
	
## Create an orthogonal slice 	
$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
 	COMMAND = 'ORTHOSLICE |CutAxis| |CutPos| Active'   
 	
## Set the name of the new frame 	 
$!FRAMENAME  = 'Slice at |CutAxis| = |CutPos|' 

## Set the position, border and background attributes of the frame
$!FRAMELAYOUT
  HEADERCOLOR = CUSTOM8
  XYPOS
    {
    X = 0.19685
    Y = 0.19685
    }
  WIDTH = 11.299
  HEIGHT = 3.837
  
## Assigne attributes for axes in an XY line plot  
$!XYLINEAXIS GRIDAREA
{
DRAWBORDER=YES
}
  
$!XYLINEAXIS YDETAIL 1
{
COORDSCALE = LINEAR
TITLE
	{
	SHOWONAXISLINE = YES
	TITLEMODE = USETEXT 
	TEXT = 'ENERGY [eV]'
	}
}       

$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
	COMMAND = 'SET_PREFERENCE ListCoordinateVariables 1'
	## Sets a preference option
$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
  	COMMAND = 'GET_VARNUM_BY_TDR_QUANTITY "Y" varnrY' 
  	## Obtains the variable number associated with the specified name of the TDR dataset and assigns it to a Tecplot macro variable.
$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
  	COMMAND = 'GET_VARNUM_BY_TDR_QUANTITY "ConductionBandEnergy" varnrCB' 
$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
  	COMMAND = 'GET_VARNUM_BY_TDR_QUANTITY "ValenceBandEnergy" varnrVB' 
####$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
#####	COMMAND = 'GET_VARNUM_BY_TDR_QUANTITY "eQuasiFermiPotential" varnrElecEf' 
####$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
 ####	COMMAND = 'GET_VARNUM_BY_TDR_QUANTITY "hQuasiFermiPotential" varnrHoleEf' 

## Convert Quasi Fermi Potentials [V] to Quasi Fermi Energies [eV]
####$!ALTERDATA
####	EQUATION = "V|varnrElecEf|=-V|varnrElecEf|"
####$!ALTERDATA
####	EQUATION = "V|varnrHoleEf|=-V|varnrHoleEf|"

## Delete all line maps
$!DELETELINEMAPS   

## Create new line mappings
$!CREATELINEMAP
$!CREATELINEMAP
#$!CREATELINEMAP
#$!CREATELINEMAP

## Assign attributes for line mappings
$!LINEMAP [1]
  NAME     = '&DV&'
  ASSIGN
    {
    ZONE     = |NUMZONES|
    XAXISVAR = |varnrY|
    YAXISVAR = |varnrCB|
    YAXIS    = 1
    }
  LINES
    {
    COLOR         = BLACK
    LINETHICKNESS = 0.4
    }
$!LINEMAP [2]
  NAME     = '&DV&'
  ASSIGN
    {
    ZONE     = |NUMZONES|
    XAXISVAR = |varnrY|
    YAXISVAR = |varnrVB|
    YAXIS    = 1
    }
  LINES
    {
    COLOR         = 	BLACK
    LINETHICKNESS = 0.4
    }
#$!LINEMAP [3]
#  NAME     = '&DV&'
#  ASSIGN
#    {
#    ZONE     = |NUMZONES|
#    XAXISVAR = |varnrY|
#    YAXISVAR = |varnrElecEf|
#    YAXIS    = 1
#    }
#  LINES
#    {
#    COLOR         = RED
#    LINETHICKNESS = 0.4
#    LINEPATTERN   = DASHDOT
#    }
#$!LINEMAP [4]
#  NAME     = '&DV&'
#  ASSIGN
#    {
#    ZONE     = |NUMZONES|
#    XAXISVAR = |varnrY|
#    YAXISVAR = |varnrHoleEf|
#    YAXIS    = 1
#    }
#  LINES
#    {
#    COLOR         = BLUE
#    LINETHICKNESS = 0.4
#    LINEPATTERN   = DASHDOT
#    }
    
## Specify the set of line-mappings considered for plotting.    
$!ACTIVELINEMAPS = [1-2]

## Fit the entire plot to the grid area
$!VIEW FIT

## Reset the range on Y axis so that it equals the minimum and maximum of the data being plotted
$!VIEW AXISNICEFIT
  AXIS    = 'Y'
  AXISNUM = 1
  
## Compute diffusion lengths in [um]
$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
  	COMMAND = 'GET_VARNUM_BY_TDR_QUANTITY "eLifetime" varnreLifetime' 
$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
  	COMMAND = 'GET_VARNUM_BY_TDR_QUANTITY "eMobility" varnreMobility' 
$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
  	COMMAND = 'GET_VARNUM_BY_TDR_QUANTITY "hLifetime" varnrhLifetime' 
$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
  	COMMAND = 'GET_VARNUM_BY_TDR_QUANTITY "hMobility" varnrhMobility' 


$!ALTERDATA
	EQUATION = "{eDiffusionLength_um}=SQRT(V|varnreLifetime|*V|varnreMobility|*0.0258520)*1e4"
$!ALTERDATA
	EQUATION = "{hDiffusionLength_um}=SQRT(V|varnrhLifetime|*V|varnrhMobility|*0.0258520)*1e4"

## Compute Diffusion coefficients in [cm^2/s]
$!ALTERDATA
	EQUATION = "{eDiffusionCoeff_cm2s-1}=V|varnreMobility|*0.0258520"
$!ALTERDATA
	EQUATION = "{hDiffusionCoeff_cm2s-1}=V|varnrhMobility|*0.0258520"

  





