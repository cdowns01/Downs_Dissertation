

set N     521
set i     3

proj_load  n393_des.plt PLT_JV($N)

set COLORS  [list orange blue red violet brown orange magenta]
set NCOLORS [llength $COLORS]
set color   [lindex  $COLORS [expr $i%$NCOLORS]]

cv_createDS J($N)  "PLT_JV($N) front OuterVoltage" "PLT_JV($N) front TotalCurrent"
cv_inv J($N) y 

cv_create V($N)  "PLT_JV($N) front OuterVoltage" "PLT_JV($N) front OuterVoltage"

cv_createWithFormula P($N)  "<V($N)>*<J($N)>" A A 


cv_setCurveAttr J($N)  "JV" $color solid  2 none 3 defcolor 1 defcolor


set Jsc($N) [cv_compute "vecvaly(<J($N)>,0)" A A A A]

gr_setAxisAttr X {Voltage (V)}  16 0 1 black 1 14 0 5 0
gr_setAxisAttr Y {Current Density (mA/cm^2)} 16 0 20 black 1 14 0 5 0

ft_scalar Jsc [format %.2f $Jsc($N)]

set Jmin [cv_compute "vecmin(<J($N)>)" A A A A]
if {$Jmin <= 0} {
	set Voc($N) [expr [cv_compute "veczero(<J($N)>)" A A A A]]
} elseif {$Jmin <= 1e-6} {
	set Voc($N) [expr [cv_compute "vecvalx(<J($N>,$Jmin)" A A A A]]
}

ft_scalar Voc [format %.4f $Voc($N)]

set Ps 10000 ;#Incident light power density for AM1.5d radiation in mW/cm^2
if {$Voc($N) > 0} {
	set Pm($N) [cv_compute "vecmax(<P($N)>)" A A A A]
	set FF($N) [expr $Pm($N)/($Voc($N)*$Jsc($N))*100]
	set Eff($N) [expr $Pm($N)/$Ps*100]
}

ft_scalar Pm [format %.2f $Pm($N)]
ft_scalar FF  [format %.2f $FF($N)]
ft_scalar Eff  [format %.2f $Eff($N)]





