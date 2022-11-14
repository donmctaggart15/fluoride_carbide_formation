# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 16:15:04 2022

@author: donmc
"""

from simmate.database import connect
from simmate.database.third_parties import MatprojStructure
import pandas as pd
import copy
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.ext.matproj import MPRester  # This version needed for interface rxns
mpr = MPRester("eXE8IJUp6lZSxLAo")

#%% Code sections

# 1. Get fluoride structures
# 2. Get carbide structures
# 3. Calculate potentials
# 4. Get formation enthalpies of carbides (all carbides independent of corresponding fluoride structures)
# 5. Create table with all data

#%% 1: GET F STRUCUTRES 

# Make empty property lists
reduced_form_f_list = []
reduced_form_c_list = []
reduced_form_m_list = []
mpid_f_list = []
mpid_c_list = []
e_per_atom_f_list = []
e_per_atom_c_list = []
e_per_atom_m_list = []
e_hull_f_list = []
e_hull_c_list = []
bandgap_f_list = []
bandgap_c_list = []
scale_factor_list = []
spacegroup_f_list = []
spacegroup_c_list = []
atoms_per_fu_f_list = []
atoms_per_fu_c_list = []
num_f_list = []
num_c_list = []
num_m_list = []
num_m_list_c = []

# Get all 2-element structures that contain F and order by increasing hull energy
all_f_strucs = MatprojStructure.objects.filter(
    elements__icontains='"F"',
    nelements=2,
).order_by('energy_above_hull').all()

# Convert to dataframe to delete repeats
df_all_f_strucs = all_f_strucs.to_dataframe()

# Keep only lowest energy structure of repeat structures, first occurance will be lowest energy, delete all others
df_filter_f_strucs = df_all_f_strucs.drop_duplicates("formula_reduced")
df_filter_f_strucs = df_filter_f_strucs.reset_index(drop=True)

# Convert to "Structure" objects to replace F easily
# Note: this makes list of lists so entries must be called with x[0] (not just x) 
ids = df_filter_f_strucs['id'].tolist()
f_strucs_list = []
f_strucs_list_orig = []
for j,x in enumerate(ids):
    struc = MatprojStructure.objects.filter(
        id=x
    ).to_toolkit()
    f_strucs_list.append(struc)
    
#%% 2: GET CARBIDE STRUCTURES (~10 min to run)

# Remove F from each structure, do another filtering where "C-element" is chemical_system
for i,struc in enumerate(f_strucs_list):
    # Get properties for F structure
    mpid_f = df_filter_f_strucs.id[i]
    energy_per_atom_f = df_filter_f_strucs.energy_per_atom[i]
    spacegroup_f = df_filter_f_strucs.spacegroup[i]
    e_hull_f = df_filter_f_strucs.energy_above_hull[i]
    bandgap_f = df_filter_f_strucs.band_gap[i]
    red_form_f = df_filter_f_strucs.formula_reduced[i]
    
    # Use to retrieve non-F elements
    data_dict_f = struc[0].composition.to_data_dict
    # Use to get number of each non-F atom in reduced form
    reduced_dict_f = struc[0].composition.to_reduced_dict
    # Get first element in F structures
    first_element_f = data_dict_f['elements'][0]
    # Get second element in F structures 
    f_element = data_dict_f['elements'][1]
    # Some strucs have F as first element, if so, use the other element
    if first_element_f == 'F':
        f_element = first_element_f
        first_element_f = data_dict_f['elements'][1]
    # Get the number of the first element in reduced form
    num_first_element_f = reduced_dict_f[first_element_f]
    # Get number of F atoms in reduced form
    num_f = reduced_dict_f[f_element]
    
    # Get atoms per formula unit for voltage calculations
    formula_units_f = struc[0].composition.get_reduced_composition_and_factor()[1]
    atoms_per_fu_f = struc[0].num_sites / formula_units_f
    
    # Remove F in each entry in f_strucs_list to get M energy
    struc_m = copy.deepcopy(struc[0])
    struc_m.remove_species("F")
    
    # Replace C for F in each entry in f_strucs_list
    struc2 = copy.deepcopy(struc[0])
    struc2.replace_species({"F": "C"})

    # For each entry check if there is another MP entry with same chemical composition
    # identical to the C-replaced structure
    possible_all = MatprojStructure.objects.filter(
        chemical_system = struc2.composition.chemical_system
    ).order_by('energy_above_hull').all()
    
    # If carbide structures found perform more operations
    if len(possible_all) > 0:
        
        # Get energy of metal
        metal = MatprojStructure.objects.filter(
            chemical_system = struc_m.composition.chemical_system
        ).order_by('energy_above_hull').to_dataframe()
        metal_e_per_atom = metal.energy_per_atom[0]
        red_form_m = metal.formula_reduced[0]
        
        # Convert to dataframe to delete repeats
        df_possible_all = possible_all.to_dataframe()
        
        # Keep only most stable structure for repeat compositions
        df_possible_all_filter = df_possible_all.drop_duplicates("formula_reduced")
        df_possible_all_filter = df_possible_all_filter.reset_index(drop=True)
        
        # Print index of F struc, number of carbide strucs before and after filtering, and chemical system 
        print("F struc index:", i, "--> # carbide strucs:", len(df_possible_all_filter), " --> chem sys:", struc[0].composition.chemical_system)
        
        # Convert to "Structure" objects to replace F easily
        ids2 = df_possible_all_filter['id'].tolist()
        carbide_list = []
        for j,x in enumerate(ids2):
            comp = MatprojStructure.objects.filter(
                id=x
            ).to_toolkit()
            carbide_list.append(comp)
        
        # For each carbide, get properties and replace C with F to match with original F strucs 
        for k,carb in enumerate(carbide_list):
            carb2 = copy.deepcopy(carb[0])
            carb2.replace_species({"C": "F"})
            chem_sys_c = carb2.composition.chemical_system
            chem_sys_f = struc[0].composition.chemical_system
            red_form_c = df_possible_all_filter.formula_reduced[k]
            mpid_c = df_possible_all_filter.id[k]
            energy_per_atom_c = df_possible_all_filter.energy_per_atom[k]
            spacegroup_c = df_possible_all_filter.spacegroup[k]
            e_hull_c = df_possible_all_filter.energy_above_hull[k]
            bandgap_c = df_possible_all_filter.band_gap[k]
        
            if chem_sys_c == chem_sys_f:
                
                # Create scale factor for carbide so that chemical equation will be balanced  (Cr3C and CrF5 --> 0.33 Cr3C and CrF5)
                # Use to retrieve non-C elements
                data_dict_c = carb[0].composition.to_data_dict
                # Use to get number of each non-C atom in reduced form
                reduced_dict_c = carb[0].composition.to_reduced_dict
                # Get first element in C structures
                first_element_c = data_dict_c['elements'][0]
                # Get second element in strucs 
                c_element = data_dict_c['elements'][1]
                # Some strucs have C as first element, if so, use the other element
                if first_element_c == 'C':
                    c_element = first_element_c
                    first_element_c = data_dict_c['elements'][1]
                # Get the number of the first element in reduced form
                num_first_element_c = reduced_dict_c[first_element_c]
                # Get number of C atoms in reduced form
                num_c = reduced_dict_c[c_element]
                # Divide number of first element for F by number of first element for C
                scale_factor= num_first_element_f / num_first_element_c 
                
                # Get atoms per formula unit for voltage calculations
                formula_units_c = carb[0].composition.get_reduced_composition_and_factor()[1]
                atoms_per_fu_c = carb[0].num_sites / formula_units_c
                
                # Add properties to lists
                reduced_form_f_list.append(red_form_f)
                reduced_form_c_list.append(red_form_c)
                reduced_form_m_list.append(red_form_m)
                mpid_f_list.append(mpid_f)
                mpid_c_list.append(mpid_c)
                e_per_atom_f_list.append(energy_per_atom_f)
                e_per_atom_c_list.append(energy_per_atom_c)
                e_per_atom_m_list.append(metal_e_per_atom)
                e_hull_f_list.append(e_hull_f)
                e_hull_c_list.append(e_hull_c)
                bandgap_f_list.append(bandgap_f)
                bandgap_c_list.append(bandgap_c)
                scale_factor_list.append(scale_factor)
                spacegroup_f_list.append(spacegroup_f)
                spacegroup_c_list.append(spacegroup_c)
                atoms_per_fu_f_list.append(atoms_per_fu_f)
                atoms_per_fu_c_list.append(atoms_per_fu_c)
                num_f_list.append(num_f)
                num_c_list.append(num_c)
                num_m_list.append(num_first_element_f)
                num_m_list_c.append(num_first_element_c)
                
                print("     Carbide found for index", i, "(", red_form_f, ") -->", red_form_c, "--> scale factor:", round(scale_factor,2))
       
                
#%% 3: Potentials

# Rxn 1: MFx + xLi + yC --> M + xLiF + yC (fluoride/metal potential)

# Rxn 2: MFx + xLi + yC --> MCy + xLiF (fluoride/carbide potential)
        
voltage_metal_list = []
voltage_carbide_list = [] 
voltage_diff_list = []

# Get energies for Li, LiF, and graphite
li = MatprojStructure.objects.filter(id = "mp-1018134").to_dataframe()
lif = MatprojStructure.objects.filter(id = "mp-1138").to_dataframe()
carbon = MatprojStructure.objects.filter(id = "mp-569304").to_dataframe()
li_e_per_atom = li.energy_per_atom[0]
lif_e_per_atom = lif.energy_per_atom[0]
carbon_e_per_atom = carbon.energy_per_atom[0]

# Rxn 1 (metal formation): MFx + xLi + yC --> M + xLiF + yC
for i,x in enumerate(reduced_form_f_list):
    products = num_m_list[i] * e_per_atom_m_list[i] + num_f_list[i] * 2 * lif_e_per_atom + scale_factor_list[i] * num_c_list[i] * carbon_e_per_atom
    reactants = atoms_per_fu_f_list[i] * e_per_atom_f_list[i] + num_f_list[i] * li_e_per_atom + scale_factor_list[i] * num_c_list[i] * carbon_e_per_atom
    voltage = (products - reactants) / -num_f_list[i]
    voltage_metal_list.append(voltage)

# Rxn 2 (carbide formation): MFx + xLi + yC --> MCy + xLiF 
for i,x in enumerate(reduced_form_f_list):
    products = scale_factor_list[i] * atoms_per_fu_c_list[i] * e_per_atom_c_list[i] + num_f_list[i] * 2 * lif_e_per_atom
    reactants = atoms_per_fu_f_list[i] * e_per_atom_f_list[i] + num_f_list[i] * li_e_per_atom + scale_factor_list[i] * num_c_list[i] * carbon_e_per_atom
    voltage = (products - reactants) / -num_f_list[i]
    voltage_carbide_list.append(voltage)  

# Difference in potentials (Rxn 2 - Rxn 1)
for i,x in enumerate(voltage_carbide_list):
    voltage_diff = voltage_carbide_list[i] - voltage_metal_list[i]
    voltage_diff_list.append(voltage_diff)

#%% 4: CARBIDE FORMATION ENTHAPLY: M + yC --> MCy
# Use given MP carbide formation energies

# Get all C structures with 2 elements and order by hull energy
all_c_strucs = MatprojStructure.objects.filter(
    elements__icontains='"C"',
    nelements=2,
).order_by('energy_above_hull').all()

# Convert to dataframe
df_all_c_strucs = all_c_strucs.to_dataframe()
# Keep only most stable structure for repeat compositions
df_filter_c_strucs = df_all_c_strucs.drop_duplicates("formula_reduced")
df_filter_c_strucs = df_filter_c_strucs.reset_index(drop=True)

# Get reduced fromula, formation energy, and formation energy per atom for carbides
carb_ref_form = df_filter_c_strucs['formula_reduced'].tolist()
carb_form_e = df_filter_c_strucs['formation_energy'].tolist()
carb_form_e_per_atom = df_filter_c_strucs['formation_energy_per_atom'].tolist()

# Make table for carbides 
carbide_enthaply_table = pd.DataFrame({
    'carbide': carb_ref_form,
    'mp_formation_enthalpy': carb_form_e,
    'mp_formation_enthalpy_per_atom': carb_form_e_per_atom,
    })

carbide_enthaply_table.to_csv('carbide_mp_enthaply_table_final.csv')

#%% 5: Make table of fluoride/carbide pairs

fluoride_carbide_table = pd.DataFrame({
    'f_form': reduced_form_f_list,
    'c_form': reduced_form_c_list,
    'voltage_metal': voltage_metal_list,
    'voltage_carbide': voltage_carbide_list,
    'voltage_diff': voltage_diff_list,
    'e_per_atom_f': e_per_atom_f_list,
    'f_e_hull': e_hull_f_list,
    'c_e_hull': e_hull_c_list,
    'f_spacegroup': spacegroup_f_list,
    'c_spacegroup': spacegroup_c_list,
    'f_bandgap': bandgap_f_list,
    'c_bandgap': bandgap_c_list,
    'f_id': mpid_f_list,
    'c_id': mpid_c_list,
    'scale_factor': scale_factor_list,
    })

fluoride_carbide_table.to_csv('fluoride_carbide_table_final.csv')
