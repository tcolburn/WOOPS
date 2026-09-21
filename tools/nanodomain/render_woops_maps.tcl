# render_woops_maps.tcl
#
# Run from the simulation-specific WOOPS directory after:
#
#     python process_woops_legs.py
#
# Example:
#
#     cd /path/to/WOOPS/012
#     vmd -dispdev text -e render_woops_maps.tcl
#
# The Python script writes:
#
#     woops_processed/woops_render_config.tcl
#
# which supplies all input paths and map specifications.

proc woops_fail {message} {
    puts stderr ""
    puts stderr "WOOPS RENDER ERROR: $message"
    puts stderr ""
    quit
}

proc require_file {path description} {
    if {![file isfile $path]} {
        woops_fail "$description not found: $path"
    }
}

proc rotation_from_to {from_vector to_vector} {
    set from_length [veclength $from_vector]
    set to_length [veclength $to_vector]

    if {$from_length < 1.0e-12 || $to_length < 1.0e-12} {
        woops_fail "Cannot construct a rotation from a zero-length vector."
    }

    set a [vecnorm $from_vector]
    set b [vecnorm $to_vector]

    set cosine [vecdot $a $b]

    if {$cosine > 1.0} {
        set cosine 1.0
    }

    if {$cosine < -1.0} {
        set cosine -1.0
    }

    # Vectors are already parallel.
    if {[expr {abs($cosine - 1.0)}] < 1.0e-10} {
        return [transidentity]
    }

    # Vectors are antiparallel.
    if {[expr {abs($cosine + 1.0)}] < 1.0e-10} {
        set trial {1.0 0.0 0.0}

        if {[expr {abs([vecdot $a $trial])}] > 0.9} {
            set trial {0.0 0.0 1.0}
        }

        set axis [vecnorm [veccross $a $trial]]

        return [transabout $axis 180.0 deg]
    }

    set axis [vecnorm [veccross $a $b]]

    set pi [expr {acos(-1.0)}]
    set angle_degrees [expr {acos($cosine) * 180.0 / $pi}]

    return [transabout $axis $angle_degrees deg]
}

proc compute_selection_dipole {selection} {
    set coordinates [$selection get {x y z}]
    set charges [$selection get charge]

    if {[llength $coordinates] == 0} {
        woops_fail "The dipole atom selection is empty."
    }

    set dipole {0.0 0.0 0.0}

    foreach coordinate $coordinates charge $charges {
        set contribution [vecscale $charge $coordinate]
        set dipole [vecadd $dipole $contribution]
    }

    if {[veclength $dipole] < 1.0e-10} {
        woops_fail "The selected solute dipole is numerically zero."
    }

    return $dipole
}

proc apply_orientation {
    molid
    translate_matrix
    align_matrix
    face_matrix
} {
    set all [atomselect $molid all]

    $all move $translate_matrix
    $all move $align_matrix
    $all move $face_matrix

    $all delete
}

# Draw a single-headed arrow from start_point to end_point.
#
# The cylindrical shaft begins at start_point.
# The only conical arrowhead terminates at end_point.
proc draw_single_headed_arrow {
    molid
    start_point
    end_point
    color_name
    shaft_radius
    head_radius
    head_length
} {
    set arrow_vector [vecsub $end_point $start_point]
    set arrow_length [veclength $arrow_vector]

    if {$arrow_length <= $head_length} {
        woops_fail "Arrow length must exceed arrow-head length."
    }

    set arrow_unit [vecnorm $arrow_vector]

    set head_base [
        vecsub \
            $end_point \
            [vecscale $head_length $arrow_unit]
    ]

    graphics $molid color $color_name

    graphics $molid cylinder \
        $start_point \
        $head_base \
        radius $shaft_radius \
        resolution 40 \
        filled yes

    graphics $molid cone \
        $head_base \
        $end_point \
        radius $head_radius \
        resolution 40
}

proc convert_tga_to_png {tga_path png_path} {
    set magick_command [auto_execok magick]
    set convert_command [auto_execok convert]

    if {$magick_command ne ""} {
        if {[catch {exec magick $tga_path $png_path} message]} {
            puts stderr "WARNING: ImageMagick conversion failed: $message"
            puts stderr "         Keeping Tachyon output: $tga_path"
            return
        }

        file delete -force $tga_path
        puts "Wrote $png_path"

        return
    }

    if {$convert_command ne ""} {
        if {[catch {exec convert $tga_path $png_path} message]} {
            puts stderr "WARNING: ImageMagick conversion failed: $message"
            puts stderr "         Keeping Tachyon output: $tga_path"
            return
        }

        file delete -force $tga_path
        puts "Wrote $png_path"

        return
    }

    puts stderr "WARNING: Neither 'magick' nor 'convert' was found."
    puts stderr "         Keeping Tachyon output as $tga_path"
}


# ============================================================================
# Read the processing configuration
# ============================================================================

set launch_root [file normalize [pwd]]

set config_path [
    file join \
        $launch_root \
        woops_processed \
        woops_render_config.tcl
]

if {![file isfile $config_path]} {
    woops_fail "Run process_woops_legs.py first; missing $config_path"
}

source $config_path

require_file $woops_original_psf "Original PSF"
require_file $woops_coordinate_pdb "Coordinate PDB"

foreach map_spec $woops_map_specs {
    lassign \
        $map_spec \
        map_name \
        map_psf \
        map_min \
        map_max \
        map_midpoint

    require_file $map_psf "Mapped $map_name PSF"
}

file mkdir $woops_image_dir

puts ""
puts "WOOPS rendering configuration"
puts "  Root:       $woops_root"
puts "  Simulation: $woops_simulation_id"
puts "  PSF:        $woops_original_psf"
puts "  PDB:        $woops_coordinate_pdb"
puts "  Images:     $woops_image_dir"
puts ""


# ============================================================================
# Display and Tachyon settings
#
# VMD screen coordinates:
#
#     x = horizontal screen direction
#     y = vertical screen direction
#     z = direction normal to the screen
#
# The protein dipole is aligned with +y.
#
# The protein is then rotated around y until PRO16 is as close to +z as
# possible while preserving the vertical dipole alignment.
# ============================================================================

color Display Background white

axes location Off

display projection Orthographic
display resize $woops_image_width $woops_image_height
display antialias on
display depthcue off

display shadows on
display ambientocclusion on
display aoambient 0.80
display aodirect 0.40

light 0 on
light 1 off
light 2 off
light 3 off

set material_name WOOPS_AOChalky
if {[lsearch -exact [material list] $material_name] < 0} {
    material add $material_name copy AOChalky
}

material change opacity   $material_name 1.0
material change ambient   $material_name 0.28
material change diffuse   $material_name 0.52
material change specular  $material_name 0.08
material change shininess $material_name 0.12

catch {
    render aasamples TachyonInternal 16
}

catch {
    render aosamples TachyonInternal 16
}

color scale method RWB

# Load VMD's bundled color-scale-bar plugin. No other display, orientation,
# label, arrow, surface, or rendering settings are changed below.
if {[catch {package require colorscalebar} colorbar_error]} {
    woops_fail "Could not load VMD colorscalebar plugin: $colorbar_error"
}

# The QuickSurf manifold is intentionally fully opaque.
set material_name AOChalky
material change opacity $material_name 1.00


# ============================================================================
# Load and orient the original physical structure
# ============================================================================

set reference_molid [
    mol new \
        $woops_original_psf \
        type psf \
        waitfor all
]

mol addfile \
    $woops_coordinate_pdb \
    type pdb \
    waitfor all \
    molid $reference_molid


# ----------------------------------------------------------------------
# Center the solute at its mass-weighted center
# ----------------------------------------------------------------------

set dipole_selection [
    atomselect \
        $reference_molid \
        $woops_dipole_selection
]

if {[$dipole_selection num] == 0} {
    woops_fail "Dipole selection is empty: $woops_dipole_selection"
}

set solute_center [
    measure center \
        $dipole_selection \
        weight mass
]

set translate_matrix [
    transoffset [
        vecscale -1.0 $solute_center
    ]
]

set reference_all [
    atomselect \
        $reference_molid \
        all
]

$reference_all move $translate_matrix


# ----------------------------------------------------------------------
# Align the physical protein dipole with the vertical +y axis
# ----------------------------------------------------------------------

set solute_dipole [
    compute_selection_dipole \
        $dipole_selection
]

set align_matrix [
    rotation_from_to \
        $solute_dipole \
        {0.0 1.0 0.0}
]

$reference_all move $align_matrix


# ----------------------------------------------------------------------
# Rotate around y until PRO16 faces the camera
# ----------------------------------------------------------------------

set label_selection [
    atomselect \
        $reference_molid \
        $woops_label_selection
]

if {[$label_selection num] != 1} {
    woops_fail \
        "Expected exactly one PRO16 label atom for '$woops_label_selection'; found [$label_selection num]."
}

set label_coordinate [
    lindex \
        [$label_selection get {x y z}] \
        0
]

set label_x [lindex $label_coordinate 0]
set label_z [lindex $label_coordinate 2]

set face_matrix [transidentity]

set radial_length [
    expr {
        sqrt(
            $label_x*$label_x +
            $label_z*$label_z
        )
    }
]

if {$radial_length > 1.0e-10} {
    set pi [expr {acos(-1.0)}]

    set face_angle [
        expr {
            atan2(
                -$label_x,
                $label_z
            ) * 180.0 / $pi
        }
    ]

    set face_matrix [
        transaxis y $face_angle deg
    ]

    $reference_all move $face_matrix
}


# ----------------------------------------------------------------------
# Confirm the final orientation numerically
# ----------------------------------------------------------------------

set oriented_dipole [
    compute_selection_dipole \
        $dipole_selection
]

set oriented_dipole_unit [
    vecnorm \
        $oriented_dipole
]

set label_coordinate [
    lindex \
        [$label_selection get {x y z}] \
        0
]

puts "COM-referenced solute dipole after orientation:"
puts "  $oriented_dipole"

puts "Normalized dipole direction:"
puts "  $oriented_dipole_unit"

puts "Target normalized direction:"
puts "  0.0 1.0 0.0"

puts "PRO16 carbonyl-O coordinate after orientation:"
puts "  $label_coordinate"


# ============================================================================
# Determine the oriented protein bounds
# ============================================================================

set protein_bounds_selection [
    atomselect \
        $reference_molid \
        $woops_surface_selection
]

if {[$protein_bounds_selection num] == 0} {
    woops_fail "Protein surface selection is empty: $woops_surface_selection"
}

set protein_bounds [
    measure minmax \
        $protein_bounds_selection
]

set protein_min [
    lindex \
        $protein_bounds \
        0
]

set protein_max [
    lindex \
        $protein_bounds \
        1
]

set protein_xmin [lindex $protein_min 0]
set protein_ymin [lindex $protein_min 1]
set protein_zmin [lindex $protein_min 2]

set protein_xmax [lindex $protein_max 0]
set protein_ymax [lindex $protein_max 1]
set protein_zmax [lindex $protein_max 2]

puts "Oriented protein bounds:"
puts "  minimum = $protein_min"
puts "  maximum = $protein_max"


# ============================================================================
# Underlying physical structure
# ============================================================================

mol modselect \
    0 \
    $reference_molid \
    $woops_base_selection

mol modstyle \
    0 \
    $reference_molid \
    VDW \
    0.600000 \
    12.000000

mol modcolor \
    0 \
    $reference_molid \
    Name

mol modmaterial \
    0 \
    $reference_molid \
    Opaque


# ============================================================================
# PRO16 marker, leader, and label
#
# The text and horizontal part of the leader are placed outside the projected
# protein boundary and in front of the protein along +z.
# ============================================================================

graphics $reference_molid materials off

set label_x [
    lindex \
        $label_coordinate \
        0
]

set label_y [
    lindex \
        $label_coordinate \
        1
]

set graphics_front_z [
    expr {
        $protein_zmax + 6.0
    }
]

set label_front_point [
    list \
        $label_x \
        $label_y \
        $graphics_front_z
]

set label_elbow [
    list \
        [expr {$protein_xmax + 2.0}] \
        $label_y \
        $graphics_front_z
]

set label_line_end [
    list \
        [expr {$protein_xmax + 4.5}] \
        $label_y \
        $graphics_front_z
]

set label_text_position [
    list \
        [expr {$protein_xmax + 4.9}] \
        [expr {$label_y - 0.25}] \
        $graphics_front_z
]

# Orange marker and leader.
graphics $reference_molid color orange

graphics $reference_molid sphere \
    $label_coordinate \
    radius 0.48 \
    resolution 30

graphics $reference_molid line \
    $label_coordinate \
    $label_front_point \
    width 4 \
    style solid

graphics $reference_molid line \
    $label_front_point \
    $label_elbow \
    width 4 \
    style solid

graphics $reference_molid line \
    $label_elbow \
    $label_line_end \
    width 4 \
    style solid

# Bold black text on the white background.
graphics $reference_molid color black

graphics $reference_molid text \
    $label_text_position \
    "PRO16" \
    size 1.65 \
    thickness 3.0


# ============================================================================
# Fixed-direction protein-dipole arrow
#
# The arrow axis passes through the protein COM:
#
#     x = 0
#     z = 0
#
# because the protein was translated so that its mass-weighted center is at
# the origin.
#
# The total axial extent is centered at y = 0. The arrow length is chosen so:
#
#   1. it is never shorter than 48 Angstrom,
#   2. the lower shaft protrudes beyond the bottom of the protein,
#   3. the complete arrowhead protrudes beyond the top of the protein.
#
# The QuickSurf manifold is opaque, so the central part of the arrow should
# be depth-occluded where it passes through the protein.
# ============================================================================

# Three times the previous dimensions.
set requested_arrow_length 36.0
set dipole_shaft_radius 0.48
set dipole_head_radius 1.44
set dipole_head_length 5.10

# Clearance beyond the molecular bounds.
set bottom_clearance 4.0
set top_clearance [
    expr {
        $dipole_head_length + 2.0
    }
]

# Because the arrow is centered at the COM, the half-length must be large
# enough to clear both the negative-y and positive-y protein extents.
set half_length_for_bottom [
    expr {
        abs($protein_ymin) + $bottom_clearance
    }
]

set half_length_for_top [
    expr {
        abs($protein_ymax) + $top_clearance
    }
]

set requested_half_length [
    expr {
        0.5 * $requested_arrow_length
    }
]

set dipole_arrow_half_length [
    expr {
        max(
            $requested_half_length,
            $half_length_for_bottom,
            $half_length_for_top
        )
    }
]

set dipole_arrow_length [
    expr {
        2 * $dipole_arrow_half_length
    }
]

# The arrow lies on the dipole/COM axis, not in front of the protein.
set dipole_arrow_x 0.0
set dipole_arrow_z 0.0

set dipole_arrow_start [
    list \
        $dipole_arrow_x \
        [expr {-$dipole_arrow_half_length + 2}] \
        $dipole_arrow_z
]

set dipole_arrow_end [
    list \
        $dipole_arrow_x \
        [expr {$dipole_arrow_half_length + 2}] \
        $dipole_arrow_z
]

puts "Dipole-arrow geometry:"
puts "  center          = 0.0 0.0 0.0"
puts "  direction       = +y"
puts "  total length    = $dipole_arrow_length"
puts "  shaft radius    = $dipole_shaft_radius"
puts "  head radius     = $dipole_head_radius"
puts "  head length     = $dipole_head_length"
puts "  lower endpoint  = $dipole_arrow_start"
puts "  upper tip       = $dipole_arrow_end"

draw_single_headed_arrow \
    $reference_molid \
    $dipole_arrow_start \
    $dipole_arrow_end \
    black \
    $dipole_shaft_radius \
    $dipole_head_radius \
    $dipole_head_length


# Clean up reference selections.
$reference_all delete
$dipole_selection delete
$label_selection delete
$protein_bounds_selection delete


# ============================================================================
# Render each WOOPS map
# ============================================================================

foreach map_spec $woops_map_specs {
    lassign \
        $map_spec \
        map_name \
        map_psf \
        map_min \
        map_max \
        map_midpoint

    puts ""
    puts "Rendering $map_name from $map_psf"

    set map_molid [
        mol new \
            $map_psf \
            type psf \
            waitfor all
    ]

    mol addfile \
        $woops_coordinate_pdb \
        type pdb \
        waitfor all \
        molid $map_molid

    apply_orientation \
        $map_molid \
        $translate_matrix \
        $align_matrix \
        $face_matrix

    # The modified PSF stores the WOOPS value in the charge field.
    mol modselect \
        0 \
        $map_molid \
        $woops_surface_selection

    mol modstyle \
        0 \
        $map_molid \
        Surf \
        1.800000 \
        0.000000 \

    mol modcolor \
        0 \
        $map_molid \
        Charge

    mol modmaterial \
        0 \
        $map_molid \
        $material_name

    mol scaleminmax \
        $map_molid \
        0 \
        $map_min \
        $map_max

    color scale midpoint $map_midpoint

    display resetview
    scale by 1.15
    display update

    # Add a fixed color bar using this map's explicit scale limits.
    # Six intervals produce seven labels, including zero for the WOOPS ranges
    # used here: P1 [-1,1], P2 [-0.5,1] or symmetric difference limits,
    # and P21 symmetric limits. The active RWB scale and map_midpoint are used.
    if {[string first "_minus_" $woops_simulation_id] >= 0} {
        switch -- $map_name {
            p1      { set colorbar_title "Delta P1" }
            p2      { set colorbar_title "Delta P2" }
            p21     { set colorbar_title "Delta P21" }
            default { set colorbar_title $map_name }
        }
    } else {
        switch -- $map_name {
            p1      { set colorbar_title "P1" }
            p2      { set colorbar_title "P2" }
            p21     { set colorbar_title "P21" }
            default { set colorbar_title $map_name }
        }
    }

    ::ColorScaleBar::color_scale_bar \
        0.80 0.055 \
        0 1 \
        $map_min $map_max \
        6 black 0 \
        -1.40 0.80 \
        1 $map_molid 0 \
        1 $colorbar_title

    display update

    set basename "${woops_simulation_id}_${map_name}"

    set tga_path [
        file join \
            $woops_image_dir \
            "${basename}.tga"
    ]

    set png_path [
        file join \
            $woops_image_dir \
            "${basename}.png"
    ]

    render TachyonInternal $tga_path

    convert_tga_to_png \
        $tga_path \
        $png_path

    mol delete $map_molid
}

puts ""
puts "WOOPS rendering complete."
puts "Image directory: $woops_image_dir"
puts ""

quit
