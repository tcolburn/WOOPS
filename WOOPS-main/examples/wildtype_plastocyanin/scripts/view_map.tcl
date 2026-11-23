# Quick VMD demonstration. Must first generate the visual PDBs:
#   python3 scripts/make_vmd_maps.py
# Then, from any directory:
#   vmd -e /path/to/examples/wildtype_plastocyanin/scripts/view_map.tcl
# Optional: WOOPS_MAP=p2 vmd -e scripts/view_map.tcl
set example [file normalize [file join [file dirname [info script]] ..]]
set map p1
if {[info exists env(WOOPS_MAP)]} { set map $env(WOOPS_MAP) }
if {[lsearch -exact {p1 p2 p21} $map] == -1} {
    puts stderr "WOOPS_MAP must be p1, p2, or p21"
    quit
}
set psf [file join $example input wildtype.psf]
set pdb [file join $example derived ${map}_map.pdb]
if {![file exists $pdb]} {
    puts stderr "Missing $pdb (run python3 scripts/make_vmd_maps.py first)"
    quit
}
mol new $psf type psf waitfor all
mol addfile $pdb type pdb waitfor all
mol delrep 0 top
mol representation Licorice 0.2 12 12
mol color Beta
mol selection "not water"
mol material Opaque
mol addrep top
color scale method BWR
display projection Orthographic
display resetview
puts "Showing WOOPS $map: values stored in visualization PDB Beta column."
puts "Do not use this derived PDB as WOOPS input: its beta selection flags are replaced."
