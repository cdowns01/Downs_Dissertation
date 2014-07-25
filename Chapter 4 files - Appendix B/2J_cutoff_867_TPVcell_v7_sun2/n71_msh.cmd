Title "Untitled"

Controls {
}

Definitions {
	Constant "Al.2Ga.8Sb_mid_windowdoping" {
		Species = "PhosphorusActiveConcentration"
		Value = 1e+17
	}
	Constant "Al.2Ga.8Sb_region5doping" {
		Species = "PhosphorusActiveConcentration"
		Value = 5e+16
	}
	Constant "Al.2Ga.8Sb_region6doping" {
		Species = "BoronActiveConcentration"
		Value = 5e+16
	}
	Constant "Al.2Ga.8Sb_mid_BSFdoping" {
		Species = "BoronActiveConcentration"
		Value = 1e+18
	}
	Constant "InAlAs_region7doping" {
		Species = "BoronActiveConcentration"
		Value = 1e+19
	}
	Constant "InAlAs_region8doping" {
		Species = "PhosphorusActiveConcentration"
		Value = 1e+19
	}
	Constant "GaSb_bot_windowdoping" {
		Species = "PhosphorusActiveConcentration"
		Value = 5e+17
	}
	Constant "GaSb_region9doping" {
		Species = "PhosphorusActiveConcentration"
		Value = 5e+16
	}
	Constant "GaSb_region10doping" {
		Species = "BoronActiveConcentration"
		Value = 5e+16
	}
	Constant "GaSb_bot_BSFdoping" {
		Species = "BoronActiveConcentration"
		Value = 1e+18
	}
	Constant "Goldregion_golddoping" {
		Species = "PhosphorusActiveConcentration"
		Value = 1e+10
	}
	Refinement "mid_windowdefinition" {
		MaxElementSize = ( 2.5 0.002 )
		MinElementSize = ( 1.25 0.001 )
	}
	Refinement "n_InGaAsdefinition" {
		MaxElementSize = ( 2.5 0.015 )
		MinElementSize = ( 1.25 0.0075 )
	}
	Refinement "p_InGaAsdefinition" {
		MaxElementSize = ( 2.5 0.035 )
		MinElementSize = ( 1.25 0.0175 )
	}
	Refinement "mid_BSFdefinition" {
		MaxElementSize = ( 2.5 0.0025 )
		MinElementSize = ( 1.25 0.00125 )
	}
	Refinement "p_tunnel2definition" {
		MaxElementSize = ( 2.5 0.0015 )
		MinElementSize = ( 1.25 0.00075 )
	}
	Refinement "n_tunnel2definition" {
		MaxElementSize = ( 2.5 0.0015 )
		MinElementSize = ( 1.25 0.00075 )
	}
	Refinement "bot_windowdefinition" {
		MaxElementSize = ( 2.5 0.005 )
		MinElementSize = ( 1.25 0.0025 )
	}
	Refinement "n_Gedefinition" {
		MaxElementSize = ( 2.5 0.015 )
		MinElementSize = ( 1.25 0.0075 )
	}
	Refinement "p_Gedefinition" {
		MaxElementSize = ( 2.5 0.1 )
		MinElementSize = ( 1.25 0.05 )
	}
	Refinement "bot_BSFdefinition" {
		MaxElementSize = ( 2.5 0.025 )
		MinElementSize = ( 1.25 0.0125 )
	}
	Refinement "golddefinition" {
		MaxElementSize = ( 2.5 0.025 )
		MinElementSize = ( 1.25 0.0125 )
	}
}

Placements {
	Constant "mid_window" {
		Reference = "Al.2Ga.8Sb_mid_windowdoping"
		EvaluateWindow {
			Element = region ["Al.2Ga.8Sb_mid_window"]
		}
	}
	Constant "n_InGaAs" {
		Reference = "Al.2Ga.8Sb_region5doping"
		EvaluateWindow {
			Element = region ["Al.2Ga.8Sb_region5"]
		}
	}
	Constant "p_InGaAs" {
		Reference = "Al.2Ga.8Sb_region6doping"
		EvaluateWindow {
			Element = region ["Al.2Ga.8Sb_region6"]
		}
	}
	Constant "mid_BSF" {
		Reference = "Al.2Ga.8Sb_mid_BSFdoping"
		EvaluateWindow {
			Element = region ["Al.2Ga.8Sb_mid_BSF"]
		}
	}
	Constant "p_tunnel2" {
		Reference = "InAlAs_region7doping"
		EvaluateWindow {
			Element = region ["InAlAs_region7"]
		}
	}
	Constant "n_tunnel2" {
		Reference = "InAlAs_region8doping"
		EvaluateWindow {
			Element = region ["InAlAs_region8"]
		}
	}
	Constant "bot_window" {
		Reference = "GaSb_bot_windowdoping"
		EvaluateWindow {
			Element = region ["GaSb_bot_window"]
		}
	}
	Constant "n_Ge" {
		Reference = "GaSb_region9doping"
		EvaluateWindow {
			Element = region ["GaSb_region9"]
		}
	}
	Constant "p_Ge" {
		Reference = "GaSb_region10doping"
		EvaluateWindow {
			Element = region ["GaSb_region10"]
		}
	}
	Constant "bot_BSF" {
		Reference = "GaSb_bot_BSFdoping"
		EvaluateWindow {
			Element = region ["GaSb_bot_BSF"]
		}
	}
	Constant "gold" {
		Reference = "Goldregion_golddoping"
		EvaluateWindow {
			Element = region ["Goldregion_gold"]
		}
	}
	Refinement "mid_windowdefinitionplacement" {
		Reference = "mid_windowdefinition"
		RefineWindow = region ["Al.2Ga.8Sb_mid_window"]
	}
	Refinement "n_InGaAsdefinitionplacement" {
		Reference = "n_InGaAsdefinition"
		RefineWindow = region ["Al.2Ga.8Sb_region5"]
	}
	Refinement "p_InGaAsdefinitionplacement" {
		Reference = "p_InGaAsdefinition"
		RefineWindow = region ["Al.2Ga.8Sb_region6"]
	}
	Refinement "mid_BSFdefinitionplacement" {
		Reference = "mid_BSFdefinition"
		RefineWindow = region ["Al.2Ga.8Sb_mid_BSF"]
	}
	Refinement "p_tunnel2definitionplacement" {
		Reference = "p_tunnel2definition"
		RefineWindow = region ["InAlAs_region7"]
	}
	Refinement "n_tunnel2definitionplacement" {
		Reference = "n_tunnel2definition"
		RefineWindow = region ["InAlAs_region8"]
	}
	Refinement "bot_windowdefinitionplacement" {
		Reference = "bot_windowdefinition"
		RefineWindow = region ["GaSb_bot_window"]
	}
	Refinement "n_Gedefinitionplacement" {
		Reference = "n_Gedefinition"
		RefineWindow = region ["GaSb_region9"]
	}
	Refinement "p_Gedefinitionplacement" {
		Reference = "p_Gedefinition"
		RefineWindow = region ["GaSb_region10"]
	}
	Refinement "bot_BSFdefinitionplacement" {
		Reference = "bot_BSFdefinition"
		RefineWindow = region ["GaSb_bot_BSF"]
	}
	Refinement "golddefinitionplacement" {
		Reference = "golddefinition"
		RefineWindow = region ["Goldregion_gold"]
	}
}

