import json
from chembl_webresource_client.utils import utils as chembl_utils
from chembl_webresource_client.new_client import new_client
from fastmcp import FastMCP  # Changed import
import os
import sys
from env_loader import load_env_vars

load_env_vars()


# Create an MCP server
chembl_host = os.getenv("CHEMBL_HOST", "0.0.0.0")
chembl_port = int(os.getenv("CHEMBL_PORT", "8051"))

mcp = FastMCP(
    "ChEMBL MCP Server",
    description="MCP server for querying the ChEMBL database."
    # Removed dependencies from here, ensure it's installed in the environment if needed by tools
    # Removed host/port from settings, will be passed to run()
)

# --- Molecule Tools ---

@mcp.tool()
def get_molecule_by_chembl_id(chembl_id: str) -> str:
    """Retrieve a specific molecule by its unique ChEMBL ID (e.g., 'CHEMBL192').
    This is the most direct way to get a molecule if its ChEMBL ID is known.
    Args:
        chembl_id: The ChEMBL ID of the molecule to retrieve.
    Returns:
        A JSON string containing the molecule data, or None if not found.

    Example Return Schema (example for chembl_id CHEMBL192):
    [{'molecule_chembl_id': 'CHEMBL192',
        'molecule_structures': {'canonical_smiles': 'CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12',
        'molfile': '\n     RDKit          2D\n\n 33 36  0  0  0  0  0  0  0  0999 V2000\n    2.1000   -0.0042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.1000    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.5375   -0.0042    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4917   -0.3667    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8792   -0.0042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8042    0.9083    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4917    1.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8792    0.6833    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.2042    0.3458    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8042   -0.2417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.2875   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1583   -0.3750    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9333   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3208   -0.0333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1875    0.6083    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8958    0.6083    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.3958   -1.0917    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.7833   -0.0042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1583   -1.0917    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.2875   -1.1125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4917    1.7708    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9333   -1.1125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3208   -1.4542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.3958   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.7833   -1.4417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0750    1.5750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8042   -0.9500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8792   -1.4542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.9958   -1.4292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4958   -1.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4167   -1.3125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.1125   -1.4500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.0375   -0.9542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  2  0\n  3 13  1  0\n  4  1  1  0\n  5  4  2  0\n  6  2  1  0\n  7  2  1  0\n  8  5  1  0\n  9 10  2  0\n 10  1  1  0\n 11  5  1  0\n 12  3  1  0\n 13 14  2  0\n 14 11  1  0\n 15  3  2  0\n 16  3  2  0\n 17 25  1  0\n 18 12  1  0\n 19 12  1  0\n 20 11  2  0\n 21  7  2  0\n 22 23  2  0\n 23 20  1  0\n 24 18  1  0\n 25 19  1  0\n 26  6  1  0\n 27 10  1  0\n 28 20  1  0\n 29 17  1  0\n 30 28  1  0\n 31 27  1  0\n 32 30  1  0\n 33 31  1  0\n  9  6  1  0\n  8  7  1  0\n 22 13  1  0\n 17 24  1  0\nM  END\n\n> <chembl_id>\nCHEMBL192\n\n> <chembl_pref_name>\nSILDENAFIL\n\n',
        'standard_inchi': 'InChI=1S/C22H30N6O4S/c1-5-7-17-19-20(27(4)25-17)22(29)24-21(23-19)16-14-15(8-9-18(16)32-6-2)33(30,31)28-12-10-26(3)11-13-28/h8-9,14H,5-7,10-13H2,1-4H3,(H,23,24,29)',
        'standard_inchi_key': 'BNRNXUUZRGQAQC-UHFFFAOYSA-N'},
        'pref_name': 'SILDENAFIL'}]
    """
    try:
        molecule_client = new_client.molecule
        results = molecule_client.filter(chembl_id=chembl_id)
        return json.dumps({"data": results[0] if results else None})
    except Exception as e:
        return json.dumps({"error": f"Failed to get molecule by ChEMBL ID '{chembl_id}'.", "details": str(e)})


@mcp.tool()
def find_molecule_by_pref_name(pref_name: str) -> str:
    """Search for molecules by their preferred name (e.g., 'Aspirin').
    Useful when the exact ChEMBL ID is not known but the common name is.
    Args:
        pref_name: The preferred name to search for.
    Returns:
        A JSON string containing a list of dictionaries, each representing a matching molecule.

    Example Return Schema (example for pref_name aspirin):
    [{'atc_classifications': ['B01AC06',
        'N02BA01',
        'N02BA51',
        'A01AD05',
        'N02BA71'],
        'availability_type': 2,
        'biotherapeutic': None,
        'black_box_warning': 0,
        'chebi_par_id': 15365,
        'chirality': 2,
        'cross_references': [{'xref_id': 'aspirin',
            'xref_name': 'aspirin',
            'xref_src': 'DailyMed'},
        {'xref_id': '144203627',
            'xref_name': 'SID: 144203627',
            'xref_src': 'PubChem'},
        {'xref_id': '144209315',
            'xref_name': 'SID: 144209315',
            'xref_src': 'PubChem'},
        {'xref_id': '144210466',
            'xref_name': 'SID: 144210466',
            'xref_src': 'PubChem'},
        {'xref_id': '170465039',
            'xref_name': 'SID: 170465039',
            'xref_src': 'PubChem'},
        {'xref_id': '17389202',
            'xref_name': 'SID: 17389202',
            'xref_src': 'PubChem'},
        {'xref_id': '17390036',
            'xref_name': 'SID: 17390036',
            'xref_src': 'PubChem'},
        {'xref_id': '174007205',
            'xref_name': 'SID: 174007205',
            'xref_src': 'PubChem'},
        {'xref_id': '26747283',
            'xref_name': 'SID: 26747283',
            'xref_src': 'PubChem'},
        {'xref_id': '26752858',
            'xref_name': 'SID: 26752858',
            'xref_src': 'PubChem'},
        {'xref_id': '47193676',
            'xref_name': 'SID: 47193676',
            'xref_src': 'PubChem'},
        {'xref_id': '50105490',
            'xref_name': 'SID: 50105490',
            'xref_src': 'PubChem'},
        {'xref_id': '85230910',
            'xref_name': 'SID: 85230910',
            'xref_src': 'PubChem'},
        {'xref_id': '87798', 'xref_name': 'SID: 87798', 'xref_src': 'PubChem'},
        {'xref_id': '90340586',
            'xref_name': 'SID: 90340586',
            'xref_src': 'PubChem'},
        {'xref_id': '14', 'xref_name': 'aspirin', 'xref_src': 'TG-GATEs'},
        {'xref_id': 'Aspirin', 'xref_name': None, 'xref_src': 'Wikipedia'}],
        'dosed_ingredient': True,
        'first_approval': 1950,
        'first_in_class': 0,
        'helm_notation': None,
        'indication_class': 'Analgesic; Antirheumatic; Antipyretic',
        'inorganic_flag': 0,
        'max_phase': 4,
        'molecule_chembl_id': 'CHEMBL25',
        'molecule_hierarchy': {'molecule_chembl_id': 'CHEMBL25',
        'parent_chembl_id': 'CHEMBL25'},
        'molecule_properties': {'alogp': '1.31',
        'aromatic_rings': 1,
        'cx_logd': '-2.16',
        'cx_logp': '1.24',
        'cx_most_apka': '3.41',
        'cx_most_bpka': None,
        'full_molformula': 'C9H8O4',
        'full_mwt': '180.16',
        'hba': 3,
        'hba_lipinski': 4,
        'hbd': 1,
        'hbd_lipinski': 1,
        'heavy_atoms': 13,
        'molecular_species': 'ACID',
        'mw_freebase': '180.16',
        'mw_monoisotopic': '180.0423',
        'num_lipinski_ro5_violations': 0,
        'num_ro5_violations': 0,
        'psa': '63.60',
        'qed_weighted': '0.55',
        'ro3_pass': 'N',
        'rtb': 2},
        'molecule_structures': {'canonical_smiles': 'CC(=O)Oc1ccccc1C(=O)O',
        'molfile': '\n     RDKit          2D\n\n 13 13  0  0  0  0  0  0  0  0999 V2000\n    8.8810   -2.1206    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    8.8798   -2.9479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    9.5946   -3.3607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   10.3110   -2.9474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   10.3081   -2.1170    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    9.5928   -1.7078    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   11.0210   -1.7018    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   11.7369   -2.1116    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   11.0260   -3.3588    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   11.0273   -4.1837    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   11.7423   -4.5949    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   10.3136   -4.5972    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   11.0178   -0.8769    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  5  7  1  0\n  3  4  2  0\n  7  8  2  0\n  4  9  1  0\n  4  5  1  0\n  9 10  1  0\n  2  3  1  0\n 10 11  1  0\n  5  6  2  0\n 10 12  2  0\n  6  1  1  0\n  7 13  1  0\nM  END\n\n> <chembl_id>\nCHEMBL25\n\n> <chembl_pref_name>\nASPIRIN\n\n',
        'standard_inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
        'standard_inchi_key': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'},
        'molecule_synonyms': [{'molecule_synonym': '8-hour bayer',
            'syn_type': 'TRADE_NAME',
            'synonyms': '8-HOUR BAYER'},
        {'molecule_synonym': 'Acetosalic Acid',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'Acetosalic Acid'},
        {'molecule_synonym': 'Acetylsalic acid',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'ACETYLSALIC ACID'},
        {'molecule_synonym': 'Acetylsalicylic Acid',
            'syn_type': 'INN',
            'synonyms': 'Acetylsalicylic Acid'},
        {'molecule_synonym': 'Acetylsalicylic Acid',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'Acetylsalicylic Acid'},
        {'molecule_synonym': 'Acetylsalicylic acid',
            'syn_type': 'ATC',
            'synonyms': 'ACETYLSALICYLIC ACID'},
        {'molecule_synonym': 'Acetylsalicylic acid',
            'syn_type': 'OTHER',
            'synonyms': 'ACETYLSALICYLIC ACID'},
        {'molecule_synonym': 'Alka rapid',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'ALKA RAPID'},
        {'molecule_synonym': 'Anadin all night',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'ANADIN ALL NIGHT'},
        {'molecule_synonym': 'Angettes 75',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'ANGETTES 75'},
        {'molecule_synonym': 'Aspirin', 'syn_type': 'USAN', 'synonyms': 'Aspirin'},
        {'molecule_synonym': 'Aspirin', 'syn_type': 'BAN', 'synonyms': 'ASPIRIN'},
        {'molecule_synonym': 'Aspirin', 'syn_type': 'BNF', 'synonyms': 'ASPIRIN'},
        {'molecule_synonym': 'Aspirin', 'syn_type': 'FDA', 'synonyms': 'ASPIRIN'},
        {'molecule_synonym': 'Aspirin', 'syn_type': 'JAN', 'synonyms': 'ASPIRIN'},
        {'molecule_synonym': 'Aspirin',
            'syn_type': 'MERCK_INDEX',
            'synonyms': 'ASPIRIN'},
        {'molecule_synonym': 'Aspirin', 'syn_type': 'OTHER', 'synonyms': 'ASPIRIN'},
        {'molecule_synonym': 'Aspirin',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'ASPIRIN'},
        {'molecule_synonym': 'Aspirin', 'syn_type': 'USP', 'synonyms': 'ASPIRIN'},
        {'molecule_synonym': 'Aspro clr',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'ASPRO CLR'},
        {'molecule_synonym': 'BAY1019036',
            'syn_type': 'RESEARCH_CODE',
            'synonyms': 'BAY1019036'},
        {'molecule_synonym': 'Bayer extra strength aspirin for migraine pain',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'BAYER EXTRA STRENGTH ASPIRIN FOR MIGRAINE PAIN'},
        {'molecule_synonym': 'Danamep',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'DANAMEP'},
        {'molecule_synonym': 'Disprin cv',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'DISPRIN CV'},
        {'molecule_synonym': 'Disprin direct',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'DISPRIN DIRECT'},
        {'molecule_synonym': 'Durlaza',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'DURLAZA'},
        {'molecule_synonym': 'Ecotrin',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'Ecotrin'},
        {'molecule_synonym': 'Enprin',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'ENPRIN'},
        {'molecule_synonym': 'Equi-Prin',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'Equi-Prin'},
        {'molecule_synonym': 'Gencardia',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'GENCARDIA'},
        {'molecule_synonym': 'Levius',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'LEVIUS'},
        {'molecule_synonym': 'Max strgh aspro clr',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'MAX STRGH ASPRO CLR'},
        {'molecule_synonym': 'Measurin',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'MEASURIN'},
        {'molecule_synonym': 'Micropirin ec',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'MICROPIRIN EC'},
        {'molecule_synonym': 'NSC-27223',
            'syn_type': 'RESEARCH_CODE',
            'synonyms': 'NSC-27223'},
        {'molecule_synonym': 'NSC-406186',
            'syn_type': 'RESEARCH_CODE',
            'synonyms': 'NSC-406186'},
        {'molecule_synonym': 'Nu-seals 300',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'NU-SEALS 300'},
        {'molecule_synonym': 'Nu-seals 600',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'NU-SEALS 600'},
        {'molecule_synonym': 'Nu-seals 75',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'NU-SEALS 75'},
        {'molecule_synonym': 'Nu-seals cardio 75',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'NU-SEALS CARDIO 75'},
        {'molecule_synonym': 'Paynocil',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'PAYNOCIL'},
        {'molecule_synonym': 'Platet',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'PLATET'},
        {'molecule_synonym': 'Platet 300',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'PLATET 300'},
        {'molecule_synonym': 'Postmi 300',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'POSTMI 300'},
        {'molecule_synonym': 'Postmi 75',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'POSTMI 75'},
        {'molecule_synonym': 'Salicylic Acid Acetate',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'Salicylic Acid Acetate'},
        {'molecule_synonym': 'Vazalore',
            'syn_type': 'TRADE_NAME',
            'synonyms': 'VAZALORE'}],
        'molecule_type': 'Small molecule',
        'natural_product': 0,
        'oral': True,
        'parenteral': False,
        'polymer_flag': False,
        'pref_name': 'ASPIRIN',
        'prodrug': 0,
        'structure_type': 'MOL',
        'therapeutic_flag': True,
        'topical': False,
        'usan_stem': None,
        'usan_stem_definition': None,
        'usan_substem': None,
        'usan_year': None,
        'withdrawn_class': None,
        'withdrawn_country': None,
        'withdrawn_flag': False,
        'withdrawn_reason': None,
        'withdrawn_year': None}]
        The above values are not from query, but an example of what the output might look like.
    """
    try:
        molecule_client = new_client.molecule
        query = molecule_client.filter(pref_name__icontains=pref_name)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to find molecule by preferred name '{pref_name}'.", "details": str(e)})


@mcp.tool()
def find_molecule_by_synonym(synonym: str) -> str:
    """Find molecules by their synonyms (e.g., 'Viagra' for Sildenafil).
    Helpful if a molecule is known by an alternative name or brand name.
    Args:
        synonym: The synonym to search for (case-insensitive exact match).
    Returns:
        A JSON string containing a list of dictionaries, each representing a matching molecule.

    Example Return Schema (example for synonym viagra):
    [{'molecule_chembl_id': 'CHEMBL192'}, {'molecule_chembl_id': 'CHEMBL1737'}] 
    """
    try:
        molecule_client = new_client.molecule
        query = molecule_client.filter(
            molecule_synonyms__molecule_synonym__iexact=synonym)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to find molecule by synonym '{synonym}'.", "details": str(e)})


@mcp.tool()
def get_molecules_by_chembl_ids(chembl_ids: list[str]) -> str:
    """Retrieve multiple molecules by providing a list of their ChEMBL IDs.
    Efficient for batch lookups.
    Args:
        chembl_ids: A list of ChEMBL IDs (e.g., ['CHEMBL25', 'CHEMBL192']).
    Returns:
        A JSON string containing a list of dictionaries, each representing a molecule.
    
    Example Return Schema:
    [{'molecule_chembl_id': 'CHEMBL25', 'pref_name': 'ASPIRIN'}, {'molecule_chembl_id': 'CHEMBL27', 'pref_name': 'PROPRANOLOL'}, {'molecule_chembl_id': 'CHEMBL192', 'pref_name': 'SILDENAFIL'}] 
    """
    try:
        molecule_client = new_client.molecule
        query = molecule_client.filter(molecule_chembl_id__in=chembl_ids)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get molecules by ChEMBL IDs: {chembl_ids}.", "details": str(e)})


# @mcp.tool()
# def get_molecule_image_svg(chembl_id: str) -> str:
#     """Get the 2D chemical structure image of a molecule as an SVG string.
#     Useful for visualizing a molecule when its ChEMBL ID is known.
#     Args:
#         chembl_id: The ChEMBL ID of the molecule.
#     Returns:
#         A JSON string containing an SVG string representing the molecule's image, or None if not found/available.
    
#     """
#     try:
#         image_client = new_client.image
#         image_client.set_format('svg')
#         svg_data = image_client.get(chembl_id)
#         return json.dumps({"svg_image": svg_data})
#     except Exception as e:
#         return json.dumps({"error": f"Failed to get molecule image for ChEMBL ID '{chembl_id}'.", "details": str(e)})


@mcp.tool()
def find_similar_molecules_by_smiles(smiles: str, similarity_threshold: int = 70) -> str:
    """Find molecules structurally similar to a given molecule represented by its SMILES string.
    Uses Tanimoto similarity. A threshold of 70 means 70% similar or more.
    Args:
        smiles: The SMILES string of the query molecule (e.g., 'CCO').
        similarity_threshold: The minimum similarity percentage (0-100, default 70).
    Returns:
        A JSON string containing a list of dictionaries, each representing a similar molecule, including its similarity score.
    
    Example Return Schema:
    [
        {'molecule_chembl_id': 'CHEMBL477888', 'similarity': '85.4166686534881591796875'}
        {'molecule_chembl_id': 'CHEMBL477889', 'similarity': '85.4166686534881591796875'}
        {'molecule_chembl_id': 'CHEMBL478779', 'similarity': '85.4166686534881591796875'}
        {'molecule_chembl_id': 'CHEMBL2304268', 'similarity': '70.1754391193389892578125'}
    ]
    """
    try:
        similarity_client = new_client.similarity
        query = similarity_client.filter(
            smiles=smiles, similarity=similarity_threshold)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to find similar molecules by SMILES '{smiles}'.", "details": str(e)})


@mcp.tool()
def find_similar_molecules_by_chembl_id(chembl_id: str, similarity_threshold: int = 70) -> str:
    """Find molecules structurally similar to a given molecule identified by its ChEMBL ID.
    Uses Tanimoto similarity. A threshold of 70 means 70% similar or more.
    Args:
        chembl_id: The ChEMBL ID of the query molecule.
        similarity_threshold: The minimum similarity percentage (0-100, default 70).
    Returns:
        A JSON string containing a list of dictionaries, each representing a similar molecule, including its similarity score.
    
    Example Return Schema:
    [{'molecule_chembl_id': 'CHEMBL2296002', 'pref_name': None, 'similarity': '100'}, {'molecule_chembl_id': 'CHEMBL1697753', 'pref_name': 'ASPIRIN DL-LYSINE', 'similarity': '100'}, {'molecule_chembl_id': 'CHEMBL3833325', 'pref_name': 'CARBASPIRIN CALCIUM', 'similarity': '88.8888895511627197265625'}, {'molecule_chembl_id': 'CHEMBL3833404', 'pref_name': 'CARBASPIRIN', 'similarity': '88.8888895511627197265625'}, '...(remaining elements truncated)...']
    """
    try:
        similarity_client = new_client.similarity
        query = similarity_client.filter(
            chembl_id=chembl_id, similarity=similarity_threshold)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to find similar molecules by ChEMBL ID '{chembl_id}'.", "details": str(e)})


@mcp.tool()
def get_approved_drugs(order_by_mw: bool = False) -> str:
    """Retrieve all drugs that have reached the maximum clinical trial phase (approved drugs).
    Args:
        order_by_mw: If True, sorts the results by molecular weight (ascending).
    Returns:
        A JSON string containing a list of dictionaries, each representing an approved drug.
    
    Example Return Schema:
    [
        {'atc_classifications': ['V03AN03'], 'availability_type': 1, 'biotherapeutic': None, 'black_box_warning': 0, 'chebi_par_id': 30217, 'chirality': 2, 'cross_references': [], 'dosed_ingredient': True, 'first_approval': 2015, 'first_in_class': 0, 'helm_notation': None, 'indication_class': 'Gases, Diluent for', 'inorganic_flag': 1, 'max_phase': 4, 'molecule_chembl_id': 'CHEMBL1796997', 'molecule_hierarchy': {'molecule_chembl_id': 'CHEMBL1796997', 'parent_chembl_id': 'CHEMBL1796997'}, 'molecule_properties': {'alogp': None, 'aromatic_rings': None, 'cx_logd': None, 'cx_logp': None, 'cx_most_apka': None, 'cx_most_bpka': None, 'full_molformula': 'He', 'full_mwt': '4.00', 'hba': None, 'hba_lipinski': None, 'hbd': None, 'hbd_lipinski': None, 'heavy_atoms': None, 'molecular_species': None, 'mw_freebase': '4.00', 'mw_monoisotopic': '4.0026', 'num_lipinski_ro5_violations': None, 'num_ro5_violations': None, 'psa': None, 'qed_weighted': None, 'ro3_pass': None, 'rtb': None}, 'molecule_structures': {'canonical_smiles': '[He]', 'molfile': '\n     RDKit          2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   -3.1917    1.0958    0.0000 He  0  0  0  0  0 15  0  0  0  0  0  0\nM  END\n\n> <chembl_id>\nCHEMBL1796997\n\n> <chembl_pref_name>\nHELIUM\n\n', 'standard_inchi': 'InChI=1S/He', 'standard_inchi_key': 'SWQJXJOGLNCZEY-UHFFFAOYSA-N'}, 'molecule_synonyms': [{'molecule_synonym': 'E939', 'syn_type': 'E_NUMBER', 'synonyms': 'E939'}, {'molecule_synonym': 'E-939', 'syn_type': 'RESEARCH_CODE', 'synonyms': 'E-939'}, {'molecule_synonym': 'Helium', 'syn_type': 'ATC', 'synonyms': 'HELIUM'}, {'molecule_synonym': 'Helium', 'syn_type': 'FDA', 'synonyms': 'HELIUM'}, {'molecule_synonym': 'Helium', 'syn_type': 'MERCK_INDEX', 'synonyms': 'HELIUM'}, {'molecule_synonym': 'Helium', 'syn_type': 'OTHER', 'synonyms': 'HELIUM'}, {'molecule_synonym': 'Helium', 'syn_type': 'USP', 'synonyms': 'HELIUM'}, {'molecule_synonym': 'INS-939', 'syn_type': 'RESEARCH_CODE', 'synonyms': 'INS-939'}, {'molecule_synonym': 'INS NO.939', 'syn_type': 'RESEARCH_CODE', 'synonyms': 'INS NO.939'}], 'molecule_type': 'Small molecule', 'natural_product': 0, 'oral': False, 'parenteral': False, 'polymer_flag': False, 'pref_name': 'HELIUM', 'prodrug': 0, 'structure_type': 'MOL', 'therapeutic_flag': False, 'topical': True, 'usan_stem': '-ium', 'usan_stem_definition': 'quaternary ammonium derivatives', 'usan_substem': '-ium', 'usan_year': None, 'withdrawn_class': None, 'withdrawn_country': None, 'withdrawn_flag': False, 'withdrawn_reason': None, 'withdrawn_year': None},
        ...
    ]
    """
    try:
        molecule_client = new_client.molecule
        query = molecule_client.filter(max_phase=4)
        if order_by_mw:
            query = query.order_by('molecule_properties__mw_freebase')
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": "Failed to get approved drugs.", "details": str(e)})

# --- Activity Tools ---


@mcp.tool()
def get_activities_for_target(target_chembl_id: str, standard_type: str = "IC50") -> str:
    """Fetch bioactivity data (e.g., IC50, Ki) for a specific biological target.
    Useful for finding out which compounds are active against a particular target.
    Args:
        target_chembl_id: The ChEMBL ID of the target (e.g., 'CHEMBL240' for EGFR).
        standard_type: The type of bioactivity measurement (e.g., 'IC50', 'Ki', 'EC50'). Default is 'IC50'. 
                       Pass None to get all activity types.
    Returns:
        A JSON string containing a list of dictionaries, each representing an activity record.
    
    Example Return Schema (for target CHEMBL255):
    13200 <int>
    """
    try:
        activity_client = new_client.activity
        query = activity_client.filter(target_chembl_id=target_chembl_id)
        if standard_type:
            query = query.filter(standard_type__iexact=standard_type)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get activities for target ChEMBL ID '{target_chembl_id}'.", "details": str(e)})


@mcp.tool()
def get_activities_for_molecule(molecule_chembl_id: str, pchembl_value_exists: bool = True) -> str:
    """Retrieve all recorded bioactivities for a specific molecule.
    Useful for understanding the biological profile of a compound.
    Args:
        molecule_chembl_id: The ChEMBL ID of the molecule.
        pchembl_value_exists: If True (default), only returns activities that have a pChEMBL value (a standardized measure of potency).
    Returns:
        A JSON string containing a list of dictionaries, each representing an activity record.

    Example Return Schema (for molecule CHEMBL25):
    138 <int>
    """
    try:
        activity_client = new_client.activity
        query = activity_client.filter(molecule_chembl_id=molecule_chembl_id)
        if pchembl_value_exists:
            query = query.filter(pchembl_value__isnull=False)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get activities for molecule ChEMBL ID '{molecule_chembl_id}'.", "details": str(e)})

# --- Target Tools ---


@mcp.tool()
def find_target_by_gene_name(gene_name: str) -> str:
    """Search for biological targets (e.g., proteins, protein families) by a gene name or symbol.
    This tool searches within target synonyms, so it can find targets even if the gene name is not the preferred name.
    Args:
        gene_name: The gene name or symbol to search for (e.g., 'BRCA1', 'EGFR'). Case-insensitive contains match.
    Returns:
        A JSON string containing a list of dictionaries, each representing a matching target.
    Example Return Schema (for gene name 'BRD4'):
    [
        {'organism': 'Homo sapiens', 'pref_name': 'Bromodomain-containing protein 4', 'target_type': 'SINGLE PROTEIN'}
        {'organism': 'Mus musculus', 'pref_name': 'Bromodomain-containing protein 4', 'target_type': 'SINGLE PROTEIN'}
        {'organism': 'Homo sapiens', 'pref_name': 'BRD4/HDAC1', 'target_type': 'PROTEIN COMPLEX'}
        {'organism': 'Homo sapiens', 'pref_name': 'Cereblon/Cullin-4A/Bromodomain-containing protein 4', 'target_type': 'PROTEIN-PROTEIN INTERACTION'}
        {'organism': 'Homo sapiens', 'pref_name': 'Cereblon/Bromodomain-containing protein 4', 'target_type': 'PROTEIN-PROTEIN INTERACTION'}
        {'organism': 'Homo sapiens', 'pref_name': 'von Hippel-Lindau disease tumor suppressor/Bromodomain-containing protein 4', 'target_type': 'PROTEIN-PROTEIN INTERACTION'}
        {'organism': 'Homo sapiens', 'pref_name': 'Cereblon/DNA damage-binding protein 1/Bromodomain-containing protein 4', 'target_type': 'PROTEIN-PROTEIN INTERACTION'}
        {'organism': 'Homo sapiens', 'pref_name': 'von Hippel-Lindau disease tumor suppressor/Elongin-B/Elongin-C/Bromodomain-containing protein 4', 'target_type': 'PROTEIN-PROTEIN INTERACTION'}
        {'organism': 'Homo sapiens', 'pref_name': 'Bromodomain and extra-terminal motif (BET)', 'target_type': 'PROTEIN FAMILY'}
        {'organism': 'Homo sapiens', 'pref_name': 'BRD4/E3 ubiquitin-protein ligase Mdm2', 'target_type': 'PROTEIN-PROTEIN INTERACTION'}
        {'organism': 'Homo sapiens', 'pref_name': 'BRD4/E3 ubiquitin-protein ligase XIAP', 'target_type': 'PROTEIN-PROTEIN INTERACTION'}
    ]
    """
    try:
        target_client = new_client.target
        query = target_client.filter(target_synonym__icontains=gene_name)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to find target by gene name '{gene_name}'.", "details": str(e)})

# --- Generic Filter Tools ---


@mcp.tool()
def get_molecules_by_filter(filters: dict[str, str], order_by: list[str] = ['molecule_properties__mw_freebase']) -> str:
    """A flexible tool to retrieve molecules based on a custom set of filter conditions.
    Use this for complex queries not covered by more specific tools.
    Args:
        filters: A dictionary where keys are ChEMBL molecule fields (e.g., 'molecule_properties__mw_freebase') 
                 and values are the conditions (e.g., {'molecule_properties__mw_freebase__gte': 200, 'max_phase': 4}).
                 Supports Django-style lookups (e.g., '__gte', '__icontains').
        order_by: Optional list of fields to sort the results by (e.g., ['molecule_properties__mw_freebase']). 
                  Prefix with '-' for descending order (e.g., ['-num_ro5_violations']).
    Returns:
        A JSON string containing a list of dictionaries, each representing a molecule matching the criteria.
    Example Usage:
        To get molecules with molecular weight >= 200, ALOGP <= 5, and max_phase = 4:
        get_molecules_by_filter(filters={"molecule_properties__mw_freebase__gte": 200, 
                                       "molecule_properties__alogp__lte": 5, 
                                       "max_phase": 4})
    """
    try:
        molecule_client = new_client.molecule
        query = molecule_client.filter(**filters)
        if order_by:
            query = query.order_by(*order_by)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get molecules by filter: {filters}.", "details": str(e)})


@mcp.tool()
def get_activities_by_filter(filters: dict[str, str], order_by: list[str] = ['molecule_properties__mw_freebase']) -> str:
    """A flexible tool to retrieve bioactivity data based on a custom set of filter conditions.
    Use this for complex queries on activity data.
    Args:
        filters: A dictionary of filter conditions for activity fields (e.g., {'pchembl_value__gte': 6.0, 'standard_type__iexact': 'IC50'}).
        order_by: Optional list of fields to sort the results by (e.g., ['-pchembl_value']).
    Returns:
        A JSON string containing a list of dictionaries, each representing an activity record.
    Example Usage:
        To get IC50 activities with pChEMBL value >= 6.0 for target CHEMBL240:
        get_activities_by_filter(filters={"pchembl_value__gte": 6.0, 
                                          "standard_type__iexact": "IC50", 
                                          "target_chembl_id": "CHEMBL240"})
    """
    try:
        activity_client = new_client.activity
        query = activity_client.filter(**filters)
        if order_by:
            query = query.order_by(*order_by)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get activities by filter: {filters}.", "details": str(e)})


@mcp.tool()
def get_targets_by_filter(filters: dict[str, str], order_by: list[str] = ['molecule_properties__mw_freebase']) -> str:
    """A flexible tool to retrieve biological targets based on a custom set of filter conditions.
    Use this for complex queries on target data.
    Args:
        filters: A dictionary of filter conditions for target fields (e.g., {'pref_name__startswith': 'Cytochrome P450', 'target_type': 'SINGLE PROTEIN'}).
        order_by: Optional list of fields to sort the results by (e.g., ['pref_name']).
    Returns:
        A JSON string containing a list of dictionaries, each representing a target matching the criteria.
    Example Usage:
        To get human single protein targets whose preferred name starts with 'Cytochrome P450':
        get_targets_by_filter(filters={"pref_name__startswith": "Cytochrome P450", 
                                       "target_type": "SINGLE PROTEIN", 
                                       "organism__istartswith": "Homo sapiens"})
    """
    try:
        target_client = new_client.target
        query = target_client.filter(**filters)
        if order_by:
            query = query.order_by(*order_by)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get targets by filter: {filters}.", "details": str(e)})

# --- Utility Tools ---


@mcp.tool()
def smiles_to_ctab(smiles: str) -> str:
    """Convert a SMILES (Simplified Molecular Input Line Entry System) string to a CTAB (Chemical Table) block.
    CTAB is a text-based format for representing chemical structures.
    Args:
        smiles: The SMILES string to convert.
    Returns:
        A JSON string containing the molecule in CTAB format.
    """
    try:
        ctab_data = chembl_utils.smiles2ctab(smiles)
        return json.dumps({"ctab": ctab_data})
    except Exception as e:
        return json.dumps({"error": f"Failed to convert SMILES '{smiles}' to CTAB.", "details": str(e)})


@mcp.tool()
def compute_molecular_descriptors(smiles: str) -> str:
    """Calculate a set of physicochemical properties and descriptors for a molecule from its SMILES string.
    Examples of descriptors include Molecular Weight, AlogP, Number of Hydrogen Bond Donors/Acceptors, etc.
    Args:
        smiles: The SMILES string of the molecule.
    Returns:
        A JSON string containing various calculated molecular descriptors.

    Example Return Schema:
    {'qed': 0.5501217966938848,
    'MolWt': 180.15899999999996,
    'TPSA': 63.60000000000001,
    'HeavyAtomCount': 13,
    'NumAromaticRings': 1,
    'NumHAcceptors': 3,
    'NumHDonors': 1,
    'NumRotatableBonds': 2,
    'MolLogP': 1.3100999999999998,
    'MolecularFormula': 'C9H8O4',
    'Ro3Pass': 0,
    'NumRo5': 0,
    'MonoisotopicMolWt': 180.042258736}
    """
    try:
        ctab = chembl_utils.smiles2ctab(smiles)
        descriptors_json_str = chembl_utils.chemblDescriptors(ctab)
        descriptors_list = json.loads(descriptors_json_str)
        return json.dumps({"data": descriptors_list[0] if descriptors_list else {}})
    except Exception as e:
        return json.dumps({"error": f"Failed to compute molecular descriptors for SMILES '{smiles}'.", "details": str(e)})


@mcp.tool()
def compute_structural_alerts(smiles: str) -> str:
    """Identify known structural alerts (e.g., toxicophores, groups associated with reactivity) within a molecule from its SMILES string.
    Useful for early-stage assessment of potential liabilities.
    Args:
        smiles: The SMILES string of the molecule.
    Returns:
        A JSON string containing a list of dictionaries, where each dictionary represents a structural alert found in the molecule.
    Example Return Schema:
    [
        {'alert_id': 1030, 'alert_name': 'Ester', 'set_name': 'MLSMR', 'smarts': '[#6]-C(=O)O-[#6]'}
        {'alert_id': 1069, 'alert_name': 'vinyl michael acceptor1', 'set_name': 'MLSMR', 'smarts': '[#6]-[CH1]=C-C(=O)[#6,#7,#8]'}
    ]
    """
    try:
        ctab = chembl_utils.smiles2ctab(smiles)
        alerts_json_str = chembl_utils.structuralAlerts(ctab)
        alerts_list = json.loads(alerts_json_str)
        return json.dumps({"data": alerts_list[0] if alerts_list else []})
    except Exception as e:
        return json.dumps({"error": f"Failed to compute structural alerts for SMILES '{smiles}'.", "details": str(e)})


@mcp.tool()
def standardize_molecule_from_smiles(smiles: str) -> str:
    """Standardize a molecular structure provided as a SMILES string.
    This process typically involves neutralizing charges, removing salts, and applying a consistent representation.
    Args:
        smiles: The SMILES string of the molecule to standardize.
    Returns:
        A JSON string containing the standardized molecular information, usually including the standardized SMILES.

    Example Return Schema:
    [{'standard_molblock': '\n     RDKit          2D\n\n 19 17  0  0  0  0  0  0  0  0999 V2000\n    0.0000   -3.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.5170   -1.9538    0.0000 Na  0  0  0  0  0 15  0  0  0  0  0  0\n   -2.9244   -0.4442    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0602    0.0590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0638    1.0590    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1924   -0.4380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3282    0.0652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5396   -0.4318    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4038    0.0714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4002    1.0714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2644    1.5744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2608    2.5744    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5324    1.5682    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3318    1.0652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.0649    2.8712    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.8623    1.8920    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8683    1.7822    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4567    2.6936    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    4.1963    3.3666    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n  3  4  1  0\n  4  5  2  0\n  4  6  1  0\n  6  7  1  0\n  7  8  2  0\n  8  9  1  0\n  9 10  2  0\n 10 11  1  0\n 11 12  1  0\n 10 13  1  0\n 13 14  2  0\n 14  7  1  0\n 15 16  2  0\n 16 17  1  0\n 17 18  2  0\n 18 19  1  0\n 19 15  1  0\nM  CHG  2   2   1   3  -1\nM  END\n'}]
    """
    try:
        ctab = chembl_utils.smiles2ctab(smiles)
        standardized_json_str = chembl_utils.standardize(ctab)
        return json.dumps({"data": json.loads(standardized_json_str)})
    except Exception as e:
        return json.dumps({"error": f"Failed to standardize molecule from SMILES '{smiles}'.", "details": str(e)})


@mcp.tool()
def get_parent_molecule_from_smiles(smiles: str) -> str:
    """For a molecule that might be a salt or a mixture (represented by its SMILES string), this tool attempts to identify and return its parent structure.
    This often means removing counter-ions or selecting the largest covalent component.
    Args:
        smiles: The SMILES string of the molecule.
    Returns:
        A JSON string containing information about the parent molecule, usually including its SMILES.
    Example Return Schema:
    [{'parent_molblock': '\n     RDKit          2D\n\n 18 18  0  0  0  0  0  0  0  0999 V2000\n   -5.5170   -1.9538    0.0000 Na  0  0  0  0  0  1  0  0  0  0  0  0\n   -2.9244   -0.4442    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0602    0.0590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0638    1.0590    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1924   -0.4380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3282    0.0652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5396   -0.4318    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4038    0.0714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4002    1.0714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2644    1.5744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2608    2.5744    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5324    1.5682    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3318    1.0652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.0649    2.8712    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.8623    1.8920    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8683    1.7822    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4567    2.6936    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    4.1963    3.3666    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  2  0\n  3  5  1  0\n  5  6  1  0\n  6  7  2  0\n  7  8  1  0\n  8  9  2  0\n  9 10  1  0\n 10 11  1  0\n  9 12  1  0\n 12 13  2  0\n 13  6  1  0\n 14 15  2  0\n 15 16  1  0\n 16 17  2  0\n 17 18  1  0\n 18 14  1  0\nM  END\n',
      'exclude': False}]
    """
    try:
        ctab = chembl_utils.smiles2ctab(smiles)
        parent_json_str = chembl_utils.getParent(ctab)
        return json.dumps({"data": json.loads(parent_json_str)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get parent molecule from SMILES '{smiles}'.", "details": str(e)})


# --- Main execution for direct run or mcp dev ---
if __name__ == "__main__":
    print(f"Starting ChEMBL MCP Server...")
    print(f"Server Name: {mcp.name}")
    sys.stdout.flush()

    # For stdio transport (if needed for local testing, though fastmcp focuses on web transports)
    # if os.getenv("MCP_TRANSPORT") == "stdio":
    #     print(f"Running ChEMBL MCP Server with stdio transport")
    #     sys.stdout.flush()
    #     mcp.run(transport="stdio")
    # else:
    # Default to SSE transport with host and port passed directly
    print(
        f"Attempting to run ChEMBL MCP Server with FastMCP SSE transport on host {chembl_host}, port {chembl_port}")
    sys.stdout.flush()
    mcp.run(transport="sse", host=chembl_host, port=chembl_port)
