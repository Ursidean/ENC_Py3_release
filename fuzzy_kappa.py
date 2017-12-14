"""
Module for the calculation of Fuzzy Kappa and Fuzzy Kappa Simulation. Note that
this implementation includes no masking maps.
"""

# Import the Map Comparison Library for calculation of Fuzzy Kappa and FKS.
import mcl


def fuzzy_kappa(amap_path, smap_path, fuzzy_coefficients):
    # Create the analysis id for the Fuzzy Kappa.
    analysis_id_fk = mcl.createAnalysis()
    # Load the actual map for the analysis of FK.
    mcl.loadMapActual(analysis_id_fk, amap_path)
    # Load the fuzzy weights for Fuzzy Kappa.
    mcl.loadFuzzyWeights(analysis_id_fk, fuzzy_coefficients)
    # Read in the map for analysis of fuzzy kappa.
    mcl.loadMapSimulated(analysis_id_fk, smap_path)
    # Calculate Fuzzy Kappa.
    fuzzy_kappa = mcl.getFuzzyKappa(analysis_id_fk)
    # Clear the analysis for Fuzzy Kappa (memory dump).
    mcl.clear(analysis_id_fk)
    return fuzzy_kappa


def fks(omap_path, amap_path, smap_path, fuzzy_trans_coefficients):
    # Create the analysis id for the FKS.
    analysis_id_fks = mcl.createAnalysis()
    # Load the original map for the analysis of FK.
    mcl.loadOriginalMap(analysis_id_fks, omap_path)
    # Load the actual map for the analysis.
    mcl.loadMapActual(analysis_id_fks, amap_path)
    # Load the fuzzy transition weights for FKS.
    mcl.loadTransitionFuzzyWeights(analysis_id_fks, fuzzy_trans_coefficients)
    # Read in the map for analysis of fuzzy Kappa Simulation.
    mcl.loadMapSimulated(analysis_id_fks, smap_path)
    # Calculate Fuzzy Kappa Simulation.
    fuzzy_kappa_sim = mcl.getFuzzyKappaSim(analysis_id_fks)
    # Clear the analysis for Fuzzy Kappa Simulation (memory dump).
    mcl.clear(analysis_id_fks)
    return fuzzy_kappa_sim
