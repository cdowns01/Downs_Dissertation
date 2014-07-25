set fileinput am15v2.txt
set fileoutput n137_am15d.txt
set fileinID  [open $fileinput r]
set fileoutID [open $fileoutput w]
set scale 100.0
set done_header 0
while { (![eof $fileinID]) } {
   gets $fileinID line
   if { [llength $line] > 0 } {
      puts $line
      if { $done_header == 0 } {
         puts $fileoutID $line
      } else {
         set wavelength [lindex $line 0]
         set intensity  [expr $scale*[lindex $line 1]]
         puts $fileoutID "$wavelength [format %1.6e $intensity]"
      }
      if { [string first "intensity in W/m2" $line] > -1 } {
         set done_header 1
      }  
   }
}

close $fileinID
close $fileoutID

