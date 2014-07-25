Title "Untitled"

Controls {
}

Definitions {
	Constant "AlInP_top_windowdoping" {
		Species = "PhosphorusActiveConcentration"
		Value = 5e+17
	}
	Constant "AlInP_region1doping" {
		Species = "PhosphorusActiveConcentration"
		Value = 5e+16
	}
	Constant "AlInP_region2doping" {
		Species = "BoronActiveConcentration"
		Value = 5e+16
	}
	Constant "AlInP_top_BSFdoping" {
		Species = "BoronActiveConcentration"
		Value = 1e+18
	}
	Constant "InAlAs_region3doping" {
		Species = "BoronActiveConcentration"
		Value = 1e+19
	}
	Constant "InAlAs_region4doping" {
		Species = "PhosphorusActiveConcentration"
		Value = 1e+19
	}
	Constant "GaInP_mid_windowdoping" {
		Species = "PhosphorusActiveConcentration"
		Value = 5e+17
	}
	Constant "GaInP_region5doping" {
		Species = "PhosphorusActiveConcentration"
		Value = 1e+17
	}
	Constant "GaInP_region6doping" {
		Species = "BoronActiveConcentration"
		Value = 1e+17
	}
	Constant "GaInP_mid_BSFdoping" {
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
	Constant "GaAs_bot_windowdoping" {
		Species = "PhosphorusActiveConcentration"
		Value = 1e+17
	}
	Constant "GaAs_region9doping" {
		Species = "PhosphorusActiveConcentration"
		Value = 1e+16
	}
	Constant "GaAs_region10doping" {
		Species = "BoronActiveConcentration"
		Value = 1e+16
	}
	Constant "GaAs_bot_BSFdoping" {
		Species = "BoronActiveConcentration"
		Value = 1e+18
	}
	Constant "Goldregion_golddoping" {
		Species = "PhosphorusActiveConcentration"
		Value = 1e+10
	}
	Refinement "top_windowdefinition" {
		MaxElementSize = ( 2.5 0.002 )
		MinElementSize = ( 1.25 0.001 )
	}
	Refinement "n_GaInPdefinition" {
		MaxElementSize = ( 2.5 0.015 )
		MinElementSize = ( 1.25 0.0075 )
	}
	Refinement "p_GaInPdefinition" {
		MaxElementSize = ( 2.5 0.025 )
		MinElementSize = ( 1.25 0.0125 )
	}
	Refinement "top_BSFdefinition" {
		MaxElementSize = ( 2.5 0.0025 )
		MinElementSize = ( 1.25 0.00125 )
	}
	Refinement "p_tunnel1definition" {
		MaxElementSize = ( 2.5 0.0015 )
		MinElementSize = ( 1.25 0.00075 )
	}
	Refinement "n_tunnel1definition" {
		MaxElementSize = ( 2.5 0.0015 )
		MinElementSize = ( 1.25 0.00075 )
	}
	Refinement "mid_windowdefinition" {
		MaxElementSize = ( 2.5 0.005 )
		MinElementSize = ( 1.25 0.0025 )
	}
	Refinement "n_InGaAsdefinition" {
		MaxElementSize = ( 2.5 0.025 )
		MinElementSize = ( 1.25 0.0125 )
	}
	Refinement "p_InGaAsdefinition" {
		MaxElementSize = ( 2.5 0.15 )
		MinElementSize = ( 1.25 0.075 )
	}
	Refinement "mid_BSFdefinition" {
		MaxElementSize = ( 2.5 0.005 )
		MinElementSize = ( 1.25 0.0025 )
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
		MaxElementSize = ( 2.5 0.15 )
		MinElementSize = ( 1.25 0.075 )
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
	Constant "top_window" {
		Reference = "AlInP_top_windowdoping"
		EvaluateWindow {
			Element = region ["AlInP_top_window"]
		}
	}
	Constant "n_GaInP" {
		Reference = "AlInP_region1doping"
		EvaluateWindow {
			Element = region ["AlInP_region1"]
		}
	}
	Constant "p_GaInP" {
		Reference = "AlInP_region2doping"
		EvaluateWindow {
			Element = region ["AlInP_region2"]
		}
	}
	Constant "top_BSF" {
		Reference = "AlInP_top_BSFdoping"
		EvaluateWindow {
			Element = region ["AlInP_top_BSF"]
		}
	}
	Constant "p_tunnel1" {
		Reference = "InAlAs_region3doping"
		EvaluateWindow {
			Element = region ["InAlAs_region3"]
		}
	}
	Constant "n_tunnel1" {
		Reference = "InAlAs_region4doping"
		EvaluateWindow {
			Element = region ["InAlAs_region4"]
		}
	}
	Constant "mid_window" {
		Reference = "GaInP_mid_windowdoping"
		EvaluateWindow {
			Element = region ["GaInP_mid_window"]
		}
	}
	Constant "n_InGaAs" {
		Reference = "GaInP_region5doping"
		EvaluateWindow {
			Element = region ["GaInP_region5"]
		}
	}
	Constant "p_InGaAs" {
		Reference = "GaInP_region6doping"
		EvaluateWindow {
			Element = region ["GaInP_region6"]
		}
	}
	Constant "mid_BSF" {
		Reference = "GaInP_mid_BSFdoping"
		EvaluateWindow {
			Element = region ["GaInP_mid_BSF"]
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
		Reference = "GaAs_bot_windowdoping"
		EvaluateWindow {
			Element = region ["GaAs_bot_window"]
		}
	}
	Constant "n_Ge" {
		Reference = "GaAs_region9doping"
		EvaluateWindow {
			Element = region ["GaAs_region9"]
		}
	}
	Constant "p_Ge" {
		Reference = "GaAs_region10doping"
		EvaluateWindow {
			Element = region ["GaAs_region10"]
		}
	}
	Constant "bot_BSF" {
		Reference = "GaAs_bot_BSFdoping"
		EvaluateWindow {
			Element = region ["GaAs_bot_BSF"]
		}
	}
	Constant "gold" {
		Reference = "Goldregion_golddoping"
		EvaluateWindow {
			Element = region ["Goldregion_gold"]
		}
	}
	Refinement "top_windowdefinitionplacement" {
		Reference = "top_windowdefinition"
		RefineWindow = region ["AlInP_top_window"]
	}
	Refinement "n_GaInPdefinitionplacement" {
		Reference = "n_GaInPdefinition"
		RefineWindow = region ["AlInP_region1"]
	}
	Refinement "p_GaInPdefinitionplacement" {
		Reference = "p_GaInPdefinition"
		RefineWindow = region ["AlInP_region2"]
	}
	Refinement "top_BSFdefinitionplacement" {
		Reference = "top_BSFdefinition"
		RefineWindow = region ["AlInP_top_BSF"]
	}
	Refinement "p_tunnel1definitionplacement" {
		Reference = "p_tunnel1definition"
		RefineWindow = region ["InAlAs_region3"]
	}
	Refinement "n_tunnel1definitionplacement" {
		Reference = "n_tunnel1definition"
		RefineWindow = region ["InAlAs_region4"]
	}
	Refinement "mid_windowdefinitionplacement" {
		Reference = "mid_windowdefinition"
		RefineWindow = region ["GaInP_mid_window"]
	}
	Refinement "n_InGaAsdefinitionplacement" {
		Reference = "n_InGaAsdefinition"
		RefineWindow = region ["GaInP_region5"]
	}
	Refinement "p_InGaAsdefinitionplacement" {
		Reference = "p_InGaAsdefinition"
		RefineWindow = region ["GaInP_region6"]
	}
	Refinement "mid_BSFdefinitionplacement" {
		Reference = "mid_BSFdefinition"
		RefineWindow = region ["GaInP_mid_BSF"]
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
		RefineWindow = region ["GaAs_bot_window"]
	}
	Refinement "n_Gedefinitionplacement" {
		Reference = "n_Gedefinition"
		RefineWindow = region ["GaAs_region9"]
	}
	Refinement "p_Gedefinitionplacement" {
		Reference = "p_Gedefinition"
		RefineWindow = region ["GaAs_region10"]
	}
	Refinement "bot_BSFdefinitionplacement" {
		Reference = "bot_BSFdefinition"
		RefineWindow = region ["GaAs_bot_BSF"]
	}
	Refinement "golddefinitionplacement" {
		Reference = "golddefinition"
		RefineWindow = region ["Goldregion_gold"]
	}
}

