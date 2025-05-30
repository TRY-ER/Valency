from fastmcp import FastMCP  # Changed import
import os
import json
import requests
import sys
from env_loader import load_env_vars

load_env_vars()


# --- Environment Variables for ADMET MCP Server ---
# Host and port for this MCP server
admet_mcp_host = os.getenv("ADMET_MCP_HOST", "0.0.0.0")
admet_mcp_port = int(os.getenv("ADMET_MCP_PORT", "8055"))  # Port for uvicorn

# Host and port for the ADMET FastAPI application
admet_api_server_host = os.getenv("ADMET_SERVER_HOST", "0.0.0.0")
admet_api_server_port = int(os.getenv("ADMET_SERVER_PORT", "7373"))
admet_api_base_url = f"http://{admet_api_server_host}:{admet_api_server_port}"

mcp = FastMCP(
    "ADMET MCP Server",
    description="MCP server for retrieving ADMET predictions for a given SMILES string.",
    # Configure host and port
    settings={"host": admet_mcp_host, "port": admet_mcp_port}
    # The 'dependencies' argument might not be directly supported in jlowin/fastmcp's constructor in the same way.
    # If needed, refer to jlowin/fastmcp documentation for managing tool dependencies.
)

# --- ADMET API Tool ---

@mcp.tool()
def get_admet_prediction(smiles: str) -> str:
    """Get ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity)
    predictions for a given SMILES string.
    Args:
        smiles: The SMILES string representing the molecule.
    Returns:
        A JSON string containing the ADMET predictions, or an error message.
    Returns Example for sample input SMILES "CCO:
    {'data': {'molecular_weight': 46.069,
        'logP': -0.0014000000000000123,
        'hydrogen_bond_acceptors': 1.0,
        'hydrogen_bond_donors': 1.0,
        'Lipinski': 4.0,
        'QED': 0.40680796565539457,
        'stereo_centers': 0.0,
        'tpsa': 20.23,
        'AMES': 0.04434225065633655,
        'BBB_Martins': 0.9718312621116638,
        'Bioavailability_Ma': 0.8777502179145813,
        'CYP1A2_Veith': 0.0006845777155831456,
        'CYP2C19_Veith': 0.0010726468644861598,
        'CYP2C9_Substrate_CarbonMangels': 0.09108298048377036,
        'CYP2C9_Veith': 0.0004959922305715736,
        'CYP2D6_Substrate_CarbonMangels': 0.09421363705769181,
        'CYP2D6_Veith': 4.968681432728772e-05,
        'CYP3A4_Substrate_CarbonMangels': 0.20037514865398406,
        'CYP3A4_Veith': 7.176195385483197e-07,
        'Carcinogens_Lagunin': 0.4459096074104309,
        'ClinTox': 0.000478909959201701,
        'DILI': 0.09767487635836006,
        'HIA_Hou': 0.9967907905578614,
        'NR-AR-LBD': 0.00019546153271221556,
        'NR-AR': 0.004891549109015614,
        'NR-AhR': 0.00011306895430607256,
        'NR-Aromatase': 0.0002000639091420453,
        'NR-ER-LBD': 0.00016521390789421276,
        'NR-ER': 0.013651646859943866,
        'NR-PPAR-gamma': 0.00010114995820913464,
        'PAMPA_NCATS': 0.6344833254814148,
        'Pgp_Broccatelli': 0.004326826200122013,
        'SR-ARE': 0.005861052265390754,
        'SR-ATAD5': 4.132302001380594e-06,
        'SR-HSE': 0.0007515676756156608,
        'SR-MMP': 0.0003461426138528623,
        'SR-p53': 4.428979354997864e-05,
        'Skin_Reaction': 0.21598085761070251,
        'hERG': 0.010593873774632811,
        'Caco2_Wang': -4.069674962390725,
        'Clearance_Hepatocyte_AZ': 24.540657720877952,
        'Clearance_Microsome_AZ': -13.152814685501244,
        'Half_Life_Obach': 19.515910397116773,
        'HydrationFreeEnergy_FreeSolv': -4.625918485103303,
        'LD50_Zhu': 1.1119359795279147,
        'Lipophilicity_AstraZeneca': -0.024667792846933788,
        'PPBR_AZ': 3.99051964480106,
        'Solubility_AqSolDB': 1.2631092878304067,
        'VDss_Lombardo': 4.972023524356644,
        'molecular_weight_drugbank_approved_percentile': 1.1244668476153548,
        'logP_drugbank_approved_percentile': 23.148507173322994,
        'hydrogen_bond_acceptors_drugbank_approved_percentile': 6.649864288483908,
        'hydrogen_bond_donors_drugbank_approved_percentile': 36.68088406359054,
        'Lipinski_drugbank_approved_percentile': 63.80379992245056,
        'QED_drugbank_approved_percentile': 34.625823962776266,
        'stereo_centers_drugbank_approved_percentile': 22.489336952307095,
        'tpsa_drugbank_approved_percentile': 9.732454439705311,
        'AMES_drugbank_approved_percentile': 17.75882124854595,
        'BBB_Martins_drugbank_approved_percentile': 82.62892594028693,
        'Bioavailability_Ma_drugbank_approved_percentile': 66.96393951143854,
        'CYP1A2_Veith_drugbank_approved_percentile': 9.112058937572701,
        'CYP2C19_Veith_drugbank_approved_percentile': 2.365257851880574,
        'CYP2C9_Substrate_CarbonMangels_drugbank_approved_percentile': 33.15238464521132,
        'CYP2C9_Veith_drugbank_approved_percentile': 5.661108956960062,
        'CYP2D6_Substrate_CarbonMangels_drugbank_approved_percentile': 47.26638231872819,
        'CYP2D6_Veith_drugbank_approved_percentile': 0.6979449398991857,
        'CYP3A4_Substrate_CarbonMangels_drugbank_approved_percentile': 23.730127956572314,
        'CYP3A4_Veith_drugbank_approved_percentile': 2.2101589763474214,
        'Carcinogens_Lagunin_drugbank_approved_percentile': 85.53702985653354,
        'ClinTox_drugbank_approved_percentile': 2.055060100814269,
        'DILI_drugbank_approved_percentile': 20.08530438154323,
        'HIA_Hou_drugbank_approved_percentile': 53.23768902675455,
        'NR-AR-LBD_drugbank_approved_percentile': 2.9856533540131833,
        'NR-AR_drugbank_approved_percentile': 13.105854982551376,
        'NR-AhR_drugbank_approved_percentile': 1.6285381930981,
        'NR-Aromatase_drugbank_approved_percentile': 6.785575804575417,
        'NR-ER-LBD_drugbank_approved_percentile': 1.7060876308646762,
        'NR-ER_drugbank_approved_percentile': 3.9550213260953857,
        'NR-PPAR-gamma_drugbank_approved_percentile': 9.422256688639006,
        'PAMPA_NCATS_drugbank_approved_percentile': 41.644048080651416,
        'Pgp_Broccatelli_drugbank_approved_percentile': 14.191547111283443,
        'SR-ARE_drugbank_approved_percentile': 5.234587049243893,
        'SR-ATAD5_drugbank_approved_percentile': 1.7060876308646762,
        'SR-HSE_drugbank_approved_percentile': 8.026366808840635,
        'SR-MMP_drugbank_approved_percentile': 6.785575804575417,
        'SR-p53_drugbank_approved_percentile': 2.2101589763474214,
        'Skin_Reaction_drugbank_approved_percentile': 22.7995347033734,
        'hERG_drugbank_approved_percentile': 8.76308646762311,
        'Caco2_Wang_drugbank_approved_percentile': 97.75106630476928,
        'Clearance_Hepatocyte_AZ_drugbank_approved_percentile': 37.689026754556025,
        'Clearance_Microsome_AZ_drugbank_approved_percentile': 13.920124079100425,
        'Half_Life_Obach_drugbank_approved_percentile': 66.61496704148895,
        'HydrationFreeEnergy_FreeSolv_drugbank_approved_percentile': 90.5389685924777,
        'LD50_Zhu_drugbank_approved_percentile': 2.093834819697557,
        'Lipophilicity_AstraZeneca_drugbank_approved_percentile': 25.979061651803022,
        'PPBR_AZ_drugbank_approved_percentile': 2.365257851880574,
        'Solubility_AqSolDB_drugbank_approved_percentile': 99.49592865451726,
        'VDss_Lombardo_drugbank_approved_percentile': 74.40868553702985},
        'error': None,
        'warning': None}
    """
    if not smiles or not isinstance(smiles, str):
        return json.dumps({"error": "Invalid input. SMILES string must be a non-empty string."})

    predict_url = f"{admet_api_base_url}/predict/"
    try:
        payload = {"smiles": smiles}
        # Ensure correct content type
        headers = {"Content-Type": "application/json"}

        # Using a timeout for the request
        response = requests.post(
            predict_url, json=payload, headers=headers, timeout=60)  # 60 seconds timeout
        response.raise_for_status()  # Raises an HTTPError for bad responses (4XX or 5XX)

        return json.dumps(response.json())

    except requests.exceptions.HTTPError as http_err:
        error_detail = "No response content"
        if http_err.response is not None:
            try:
                error_detail = http_err.response.json()  # Try to get JSON error from FastAPI
            except json.JSONDecodeError:
                error_detail = http_err.response.text  # Fallback to text
        return json.dumps({"error": f"HTTP error occurred: {http_err}", "details": error_detail})
    except requests.exceptions.ConnectionError as conn_err:
        return json.dumps({"error": f"Connection error: Could not connect to ADMET service at {predict_url}.", "details": str(conn_err)})
    except requests.exceptions.Timeout as timeout_err:
        return json.dumps({"error": f"Request timed out while contacting ADMET service at {predict_url}.", "details": str(timeout_err)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get ADMET prediction for SMILES '{smiles}'.", "details": str(e)})


# --- Main execution for direct run or mcp dev ---
if __name__ == "__main__":
    print(f"Starting ADMET MCP Server...")
    print(f"Server Name: {mcp.name}")
    print(f"ADMET FastAPI app expected at: {admet_api_base_url}")
    # The host and port for SSE are now passed to mcp.run()
    print(
        f"Attempting to run ADMET MCP Server with FastMCP SSE transport on host {admet_mcp_host}, port {admet_mcp_port}")
    sys.stdout.flush()

    # mcp.run() from jlowin/fastmcp will start a FastAPI server.
    # For SSE, specify transport, host, and port directly.
    mcp.run(transport="sse", host=admet_mcp_host, port=admet_mcp_port)
