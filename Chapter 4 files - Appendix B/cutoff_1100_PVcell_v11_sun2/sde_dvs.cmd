(sde:clear)
(define node "@node@")

(define MATERIALS (list "GaInP" "GaInP" "GaInP" "GaInP" "InAlAs" "InAlAs" "AlGaAs" "AlGaAs" "AlGaAs" "AlGaAs" "AlInP" "InGaAs" "InGaAs" "InGaAs" "InGaAs" "InGaAs" "GaAs" "Gold"))

(define WIDTHS (list @width@ @width@  @width@ @width@ @width@ @width@ @width@ @width@ @width@ @width@ @width@ @width@ @width@ @width@ @width@ @width@ @width@ @width@))
(define HEIGHTS (list .04 .3 @top_depth@ .05 @tunnel_depth@ @tunnel_depth@ .1 .5 1.25 .1 @tunnel_depth@ @tunnel_depth@ .1 .3 .7 .1 .5 .5))
(define REGION_LABELS (list "_top_window" "_region1" "_region2" "_top_BSF"  "_region3"  "_region4" "_mid_window" "_region5" "_region6" "_mid_BSF" "_region7" "_region8" "_bot_window" "_region9" "_region10" "_bot_BSF" "bot_contact" "region_gold"))

(define DOPING_CONCENTRATIONS (list 1e17 5e16 5e16 1e18 1e19 1e19 5e17 5e16 5e16 1e18 1e19 1e19 5e17 1e17 5e16 @doping@ @contact_doping@ 1e10))
(define DOPANT_SPECIES (list "PhosphorusActiveConcentration" "PhosphorusActiveConcentration" "BoronActiveConcentration" "BoronActiveConcentration" "BoronActiveConcentration" "PhosphorusActiveConcentration" "PhosphorusActiveConcentration" "PhosphorusActiveConcentration" "BoronActiveConcentration" "BoronActiveConcentration" "BoronActiveConcentration" "PhosphorusActiveConcentration" "PhosphorusActiveConcentration" "PhosphorusActiveConcentration" "BoronActiveConcentration" "BoronActiveConcentration" "BoronActiveConcentration" "PhosphorusActiveConcentration"))

(define PROFILE_NAMES (list "top_window" "n_GaInP" "p_GaInP" "top_BSF" "p_tunnel1" "n_tunnel1" "mid_window" "n_InGaAs" "p_InGaAs" "mid_BSF" "p_tunnel2" "n_tunnel2" "bot_window" "n_Ge" "p_Ge" "bot_BSF" "bot_contact" "gold"))

(define Y 0)
(define X 0)

(sdegeo:create-rectangle (position @<-width/2>@ 0 0) (position @<-width/2+contact_width>@ @<-ARC1_depth-ARC2_depth>@ 0) "Gold" "Gold_region_1")
(sdegeo:create-rectangle (position @<-width/2+contact_width>@ 0 0) (position @<width/2>@ @<-ARC1_depth>@ 0) "TiOx" "TiOx_region_1")
(sdegeo:create-rectangle (position @<-width/2+contact_width>@ @<-ARC1_depth>@ 0) (position @<width/2>@ @<-ARC1_depth-ARC2_depth>@ 0) "MgF" "MgF_region_1")

(for-each
	(lambda (MATERIAL WIDTH HEIGHT DOPING_CONCENTRATION DOPANT_SPECIE PROFILE_NAME REGION_LABEL) ;names of local variables
		(begin
			(define REGION (string-append MATERIAL REGION_LABEL))
			
			;device structure created
			(sdegeo:create-rectangle (position (- 0 (/ WIDTH 2)) Y 0.0) (position (/ WIDTH 2) (+ Y HEIGHT) 0.0) MATERIAL REGION)
			(set! Y (+ Y HEIGHT))
			
			;doping profiles applied
			(sdedr:define-constant-profile (string-append REGION "doping") DOPANT_SPECIE DOPING_CONCENTRATION)
			(sdedr:define-constant-profile-region PROFILE_NAME (string-append REGION "doping") REGION)			
			
			(define DEFINITION (string-append PROFILE_NAME "definition"))
			
			;refinement regions applied
			(sdedr:define-refinement-size  DEFINITION (/ WIDTH 20) (/ HEIGHT 20)  (/ WIDTH 40) (/ HEIGHT 40))
			(sdedr:define-refinement-region  (string-append DEFINITION "placement") DEFINITION REGION)

		)
	) MATERIALS WIDTHS HEIGHTS DOPING_CONCENTRATIONS DOPANT_SPECIES PROFILE_NAMES REGION_LABELS;lists of global variables
)

(sdegeo:define-contact-set "front" 4.0 (color:rgb 1 1 0) "##")
(sdegeo:define-contact-set "back" 4.0 (color:rgb 1 0 0) "##")

(sdegeo:define-2d-contact (find-edge-id (position @<-width/2+.1>@ @<-ARC1_depth-ARC2_depth>@ 0 )) "back")
(sdegeo:define-2d-contact (find-edge-id (position 0.1 Y 0 )) "front")

(sde:build-mesh "snmesh" "" (string-append "n" node "_msh") )

