#!MC 1100

$!FRAMECONTROL DELETETOP

$!READDATASET "n76_des.tdr" DATASETREADER = "SWB-Loader"  ##add node number for sdevice to plot appropriate response

$!VIEW FIT


$!DRAWGRAPHICS FALSE

$!VARSET |CutAxis| = 'X'

$!VARSET |CutPos| = -25

$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
	COMMAND = 'SLICEMODE MergedZones' 
	
$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
 	COMMAND = 'ORTHOSLICE |CutAxis| |CutPos| Active'   
 	
$!FRAMENAME  = 'Slice at |CutAxis| = |CutPos|' 

$!FRAMELAYOUT
  HEADERCOLOR = CUSTOM8
  XYPOS
    {
    X = 0.19685
    Y = 0.19685
    }
  WIDTH = 11.299
  HEIGHT = 3.837
  
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
$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
  	COMMAND = 'GET_VARNUM_BY_TDR_QUANTITY "Y" varnrY' 
$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
  	COMMAND = 'GET_VARNUM_BY_TDR_QUANTITY "ConductionBandEnergy" varnrCB' 
$!ADDONCOMMAND ADDONID = 'Sentaurus Workbench Add-on' 
  	COMMAND = 'GET_VARNUM_BY_TDR_QUANTITY "ValenceBandEnergy" varnrVB' 


$!DELETELINEMAPS   

$!CREATELINEMAP
$!CREATELINEMAP

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
    
$!ACTIVELINEMAPS = [1-2]

$!VIEW FIT

$!VIEW AXISNICEFIT
  AXIS    = 'Y'
  AXISNUM = 1
  
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

$!ALTERDATA
	EQUATION = "{eDiffusionCoeff_cm2s-1}=V|varnreMobility|*0.0258520"
$!ALTERDATA
	EQUATION = "{hDiffusionCoeff_cm2s-1}=V|varnrhMobility|*0.0258520"

  






