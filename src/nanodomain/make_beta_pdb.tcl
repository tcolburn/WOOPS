set psf [lindex $argv 0]
set pdb [lindex $argv 1]
set out [lindex $argv 2]

mol new $psf
mol addfile $pdb

set all  [atomselect top "all"]
set prot [atomselect top "protein"]
set wat  [atomselect top "water"]

$all set beta 0
$prot set beta 1
$wat set beta 2

$all writepdb $out

quit
