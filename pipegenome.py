#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 10:23:31 2026

@author: igor
"""

import os
import shutil
import subprocess
import pandas as pd
import cobra
from escher import Builder

def choose_file(prompt_title="Select a file", file_types=[("All files", "*.*")]):
    """
    Opens a native file dialog to select a file, mimicking R's choose.file()
    """
    import tkinter as tk
    from tkinter import filedialog
    
    root = tk.Tk()
    root.withdraw()       
    root.lift()           
    root.attributes('-topmost', True) 
    
    file_path = filedialog.askopenfilename(title=prompt_title, filetypes=file_types)
    root.destroy()        
    return file_path

def run_bioinformatics_pipeline(faa_input, emapper_annotations=None, output_prefix="microalgae", carve_executable_path="carve"):
    """
    Automated pipeline for microalgae metabolic reconstruction, 
    targeted pathway extraction (C, N, and Calcification), and Escher visualization.
    """
    sbml_model_path = f"{output_prefix}_model.xml"
    csv_output_path = f"{output_prefix}_metabolic_pathways.csv"
    escher_output_path = f"{output_prefix}_metabolic_map.html"
    
    # -------------------------------------------------------------------------
    # PRE-FLIGHT CHECKS: Validate environment dependencies and input files
    # -------------------------------------------------------------------------
    print("[INFO] Performing pre-flight checks...")
    
    if not shutil.which(carve_executable_path):
        print(f"[CRITICAL ERROR] The CarveMe executable '{carve_executable_path}' was not found.")
        print("-> Please ensure CarveMe is installed and the correct Conda environment is activated.")
        return

    if not os.path.exists(faa_input):
        print(f"[CRITICAL ERROR] Input FASTA file not found at: '{faa_input}'")
        return

    if emapper_annotations and not os.path.exists(emapper_annotations):
        print(f"[WARNING] eggNOG-mapper file not found at: '{emapper_annotations}'")
        print("-> Pipeline will proceed with standard CarveMe de novo reconstruction.")
        emapper_annotations = None

    # -------------------------------------------------------------------------
    # STEP 1: Genome-scale Metabolic Reconstruction via CarveMe
    # -------------------------------------------------------------------------
    print("[INFO] Launching CarveMe for metabolic reconstruction...")
    
    # FIX: Removed '--protein' and '--eukaryote' as they are not recognized by CarveMe CLI.
    # CarveMe treats protein FASTA as the default input type.
    carve_cmd = [carve_executable_path, faa_input, "-o", sbml_model_path]
    
    # Eukaryotic annotation guidance is supplied here via the '-i' (init) flag
    if emapper_annotations:
        print(f"[INFO] Integrating functional annotations from eggNOG-mapper...")
        carve_cmd.extend(["-i", emapper_annotations])
    
    try:
        # Run the corrected command
        subprocess.run(carve_cmd, check=True)
        print(f"[SUCCESS] SBML model successfully generated: {sbml_model_path}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] CarveMe execution failed during runtime: {e}")
        return

    # -------------------------------------------------------------------------
    # STEP 2: Smart Parsing of eggNOG-mapper Annotations
    # -------------------------------------------------------------------------
    egg_data = {}
    if emapper_annotations:
        print("[INFO] Parsing eggNOG-mapper file dynamically to extract KO and EC numbers...")
        try:
            header_columns = None
            with open(emapper_annotations, 'r') as f:
                for line in f:
                    if line.startswith('#query'):
                        header_columns = line.lstrip('#').strip().split('\t')
                        break
            
            if header_columns:
                df_egg = pd.read_csv(emapper_annotations, sep='\t', comment='#', names=header_columns)
                for _, row in df_egg.iterrows():
                    gene_id = str(row['query'])
                    egg_data[gene_id] = {
                        'kos': str(row['KEGG_ko']) if 'KEGG_ko' in row and pd.notna(row['KEGG_ko']) else "",
                        'ecs': str(row['EC']) if 'EC' in row and pd.notna(row['EC']) else ""
                    }
        except Exception as e:
            print(f"[WARNING] Failed to parse eggNOG file dynamically: {e}.")

    # -------------------------------------------------------------------------
    # STEP 3: COBRApy Model Loading & Feature Extraction
    # -------------------------------------------------------------------------
    print("[INFO] Loading SBML model into COBRApy...")
    try:
        model = cobra.io.read_sbml_model(sbml_model_path)
    except Exception as e:
        print(f"[ERROR] Failed to read the generated SBML file: {e}")
        return
    
    pathway_data = []
    print("[INFO] Extracting metabolic pathways and screening for target enzymes...")
    
    for reaction in model.reactions:
        bigg_id = reaction.id
        pathway_name = reaction.subsystem if reaction.subsystem else "Unclassified"
        
        associated_genes = [gene.id for gene in reaction.genes]
        genes_str = "; ".join(associated_genes)
        
        ec_numbers = reaction.annotation.get("ec-code", "")
        if isinstance(ec_numbers, list): ec_numbers = "; ".join(ec_numbers)
        
        ko_codes = reaction.annotation.get("kegg.orthology", "")
        if isinstance(ko_codes, list): ko_codes = "; ".join(ko_codes)
        
        if associated_genes:
            for g in associated_genes:
                if g in egg_data:
                    if not ko_codes or ko_codes == "-": ko_codes = egg_data[g]['kos']
                    if not ec_numbers or ec_numbers == "-": ec_numbers = egg_data[g]['ecs']

        has_carbonic_anhydrase = "NÃO"
        has_urease = "NÃO"
        
        rxn_name_lower = reaction.name.lower()
        bigg_id_upper = bigg_id.upper()
        
        if "4.2.1.1" in str(ec_numbers) or "carbonic anhydrase" in rxn_name_lower or "CA" in bigg_id_upper.split("_"):
            has_carbonic_anhydrase = "SIM"
            
        if "3.5.1.5" in str(ec_numbers) or "urease" in rxn_name_lower or "UREA" in bigg_id_upper:
            has_urease = "SIM"
            
        pathway_data.append({
            "Nome da Via": pathway_name,
            "ID da Reacao (BiGG)": bigg_id,
            "Codigo KO (KOfam)": ko_codes,
            "Numero EC": ec_numbers,
            "Genes Associados": genes_str,
            "Presenca de Anidrase Carbonica": has_carbonic_anhydrase,
            "Presenca de Urease": has_urease
        })
        
    df_output = pd.DataFrame(pathway_data)
    df_output.to_csv(csv_output_path, index=False, sep=",")
    print(f"[SUCCESS] Analytical CSV generated: {csv_output_path}")

    # -------------------------------------------------------------------------
    # STEP 4: Escher Map Generation (Filtered by Genome Presence)
    # -------------------------------------------------------------------------
    print("[INFO] Generating eukaryotic metabolic map in Escher...")
    reaction_presence = {rxn.id: 1 for rxn in model.reactions}
    
    try:
        builder = Builder(
            map_name="RECON1.Metabolic Core", 
            reaction_data=reaction_presence
        )
        builder.reaction_styles = ['color', 'size', 'abs']
        builder.reaction_scale = [
            {'type': 'min', 'color': '#e0e0e0', 'size': 3},
            {'type': 'max', 'color': '#2ca02c', 'size': 10}
        ]
        builder.save_html(escher_output_path)
        print(f"[SUCCESS] Interactive Escher map generated: {escher_output_path}")
    except Exception as e:
        print(f"[ERROR] Escher visualization failed: {e}")
        
    print("="*70)
    print("[FINISHED] Pipeline successfully processed Carbon, Nitrogen, and Calcification screening.")

# --- Execution Block ---
if __name__ == "__main__":
    print("[INFO] Please select your input files in the pop-up windows...")
    
    input_faa_file = choose_file(
        prompt_title="Select the Microalgae FASTA file (.faa)",
        file_types=[("FASTA protein files", "*.faa"), ("All files", "*.*")]
    )
    
    eggnog_file = choose_file(
        prompt_title="Select eggNOG-mapper Annotations (Cancel if none)",
        file_types=[("Annotation files", "*.annotations"), ("TSV files", "*.tsv"), ("All files", "*.*")]
    )
    
    if not input_faa_file:
        print("[CRITICAL ERROR] No input FASTA file selected. Pipeline aborted.")
    else:
        if not eggnog_file:
            eggnog_file = None
            
        # If 'carve' is globally accessible in this conda env, keep it as "carve"
        # Otherwise, replace with the absolute path to your carve binary
        target_carve = "carve" 
        
        run_bioinformatics_pipeline(
            faa_input=input_faa_file, 
            emapper_annotations=eggnog_file,
            output_prefix="algae_metabolism_results",
            carve_executable_path=target_carve
        )