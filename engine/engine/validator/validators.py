from engine.validator.base import Validator
from rdkit import Chem
from typing import Tuple, List
from engine.discriminator.similarity_w_rcsb import RCSBSearcher, RCSBQuery
import re


class MolValidator(Validator):

    def __init__(self):
        self.type = "MOL"

    def validate(self, source_str: str):
        try:
            mol = Chem.MolFromSmiles(source_str)
            if mol is None:
                return False
            return True
        except:
            return False


class ProtValidator(Validator):

    def __init__(self):
        self.type = "PROT"

    def validate(self, source_str: str):
        try:
            query = RCSBQuery(entry_id=source_str, rows=1)
            searcher = RCSBSearcher()
            results = searcher.search(query)
            print("results >>", results)
            print("comparable >>", results["result_set"][0]["identifier"])
            print("source str>>", source_str) 
            if source_str == results["result_set"][0]["identifier"]:
                return True
            else:
                return False
        except Exception as e:
            return False


class PolymerValidator(Validator):

    def __init__(self):
        self.type = "POLY"

    def validate(self, source_str: str):
        # check if the source str contains at least 2 "*" characters
        if source_str.count("*") >= 2:
            mol = Chem.MolFromSmiles(source_str)
            if mol is None:
                return False
            return True
        else:
            return False


class PSMILESValidator(Validator):

    def __init__(self):
        self.type = "PSMILES"

    def validate(self, smiles):
        """
    Checks if a SMILES string represents a valid molecule.

    Args:
        smiles: The SMILES string to validate.

    Returns:
        True if the SMILES is valid, False otherwise.
    """
        smiles = smiles.strip().replace(" ", "").replace("[*]", "*")
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                return True
            else:
                return False
        except:
            return False


class PolyWDGStringValidator():
    """
    Class for validating the string format of weighted directed graph
    data.This class provides methods to validate the format of a
    datapoint string. This is a specific string format that is used
    for storing weighted directed polymer data in a parsable format.

    The format is as follows:
        [monomer1].[monomer2]|[fraction_of_monomer1]|[fraction_of_monomer2]|<[polymer_rule1]<[polymer_rule2]
    The polymer rule has an own format in it. Which is as follows:
        [[atom_index1]-[atom_index2]]:[fraction_of_bond_between_atom1_to_atom2]:[fraction_of_bond_between_atom2_to_atom1]

    This format is explicitly used for formatting the input for
    Weighted Directed Message Passing Neural Networks (wD-MPNN).
    The input format holds a SMART notation and regular expression
    formatting to keep molecular data with corresponding bonds and
    weights. Irrespective of this explicit usecase, the formatting
    can allow featurization of same data for other graph based neural
    networks.

    The validate method validates the proper formatting for monomer
    molecules, proper value of the fractions and valid atom indicies
    and corresponding weights in the polymer rules.

    Example
    -------
    >>> from deepchem.utils import PolyWDGStringValidator
    >>> validator = PolyWDGStringValidator()
    >>> validator.validate("[1*]C.C[2*]|0.5|0.5|<1-2:0.5:0.5")
    True

    References
    ----------
    .. [1] Aldeghi, Matteo, and Connor W. Coley. "A graph representation of molecular
        ensembles for polymer property prediction." Chemical Science 13.35 (2022): 10486-10498.
    """

    def __init__(self):
        self.type = "WDG"

    @staticmethod
    def get_parsed_vals(datapoint: str) -> Tuple[str, list, str]:
        """
        This static method parses the datapoint string into 3 parts:
        1. Monomer molecules
        2. Fragments
        3. Polymer rules

        Parameters
        ----------
        datapoint : str
            The datapoint string to parse

        Returns
        -------
        Tuple[str, list, str]
            A tuple containing the 3 parts of the datapoint string
        """
        base_parsed = datapoint.split("|")
        if len(base_parsed) < 3:
            raise ValueError(
                f"Invalid datapoint format: At least 3 splits should be there but found {len(base_parsed)} no. of splits"
            )
        monomer_mols = base_parsed[0]
        polymer_rules = base_parsed[-1]
        fragments = base_parsed[1:-1]
        return monomer_mols, fragments, polymer_rules

    @staticmethod
    def get_polymer_rules(rules_str: str) -> List[str]:
        """
        This static method parses the polymer rules string into a list of rules.

        Parameters
        ----------
        rules_str : str
            The polymer rules string to parse

        Returns
        -------
        List[str]
            A list containing the parsed rule strings
        """
        if len(rules_str.split("<")) == 1:
            raise ValueError(
                "Invalid rules string: The rule string must contain '<' as a separator for rules !"
            )
        return rules_str.split("<")[1:]

    def _validate_fragments(self, datapoint: str):
        """
        This method validate the number of fragments match
        the number of monomers.

        Parameters
        ----------
        datapoint : str
            The datapoint string to validate

        Raises
        ------
        ValueError
            If the number of fragments does not match the number of monomers
        """
        monomer_mols, fragments, _ = self.get_parsed_vals(datapoint)
        if len(fragments) != len(monomer_mols.split(".")):
            raise ValueError(
                f"Number of fragments and number of monomers should match. Mismatch -> No. of Fragments {len(fragments)} , No. of Monomers{len(monomer_mols.split('.'))}"
            )

    def _get_all_wildcards(self, text: str) -> List[str]:
        """
        This method returns all the wildcards present in the given string
        representation by using regular expression to detect digits after
        '*'.

        Parameters
        ----------
        text : str

        Returns
        -------
        List[str]
            A list of all wildcards present in the text
        """
        matches = re.findall(r"\d+(?=\*)", text)
        return matches

    def _validate_wildcards(self, datapoint: str):
        """
        This method validates the presence of wildcards in the polymer
        molecules string and ensures that the sequence of the wildcard
        notation is proper.

        Parameters
        ----------
        datapoint : str
            The datapoint string to validate

        Raises
        ------
        ValueError
            If the wildcards are not present in the sequce the maximum
            wildcard value, ValueError is raised.
        """
        monomer_mols, _, _ = self.get_parsed_vals(datapoint)
        max_wildcard = max(
            [int(x) for x in self._get_all_wildcards(monomer_mols)])
        for wildcard in range(1, max_wildcard + 1):
            if str(wildcard) + "*" not in monomer_mols:
                raise ValueError(
                    f"Invalid wildcard format: The wildcard {wildcard} is not present in the monomer molecules string  as per the sequence of the maximum {max_wildcard}!"
                )

    def _validate_polymer_rules(self, datapoint: str):
        """
        This method validates the format of the polymer rules string
        by checking for the presence of the '-' separator between the
        atom indexes, the correct number of splits in the rule string,
        and the validity of the atom indexes present in the monomer
        SMILES. It also checks if the atom indexes are in the correct
        correct count for a valid bond formation.

        Parameters
        ----------
        datapoint : str
            The datapoint string to validate

        Raises
        ------
        ValueError
            If the polymer rules string is invalid, ValueError is raised
            with appropriate error messages
        """
        monomer_mols, _, polymer_rules = self.get_parsed_vals(datapoint)
        polymer_rule_list = self.get_polymer_rules(polymer_rules)
        for rules in polymer_rule_list:
            splits = rules.split(":")
            if len(splits) != 3:
                raise ValueError(
                    f"Invalid polymer rule format: The rule must contain exactly 3 splits ! but found {len(splits)} splits"
                )
            if "-" not in splits[0]:
                raise ValueError(
                    f"Invalid polymer rule format: The bond string between two wildcard index must be seprated by '-', got invalid data {splits[0]}"
                )
            elif len(splits[0].split("-")) != 2 and any(
                    elem == "" for elem in splits[0].split("-")):
                raise ValueError(
                    f"Invalid polymer rule format: The first split must contain exactly 2 splits to depict connection between atom indexes! but found {len(splits[0].split('-'))} splits"
                )
            else:
                for wild_card_index in splits[0].split("-"):
                    if not wild_card_index.isdigit():
                        raise ValueError(
                            f"Invalid polymer rule format: The first split must contain only digits! but found {wild_card_index}"
                        )
                    if wild_card_index not in monomer_mols:
                        raise ValueError(
                            f"Invalid polymer rule format: The first split must contain only valid wild card indexes! but found {wild_card_index} which is not in {monomer_mols}"
                        )

    def _validate(self, datapoint: str):
        """
        This method validates the string format of weighted
        directed graph data. To validate the string format
        it checks for following conditions:

        1. The number of fragments and the number of monomer
           molecules should match.
        2. The wild card indexes should be present in the monomer
           molecules string and should be in the correct sequence.
        3. The polymer rules should be in the correct format.
        4. The atom indexes in the polymer rules should be valid
           and present in the monomer molecules string.

        It raises ValueError if the string format is invalid.

        Parameters
        ----------
        datapoint : str
            The datapoint string to validate

        Returns
        -------
        bool
            True if the string format is valid, None otherwise
            (Error will be raised otherwise)
        """
        self._validate_fragments(datapoint)
        self._validate_wildcards(datapoint)
        self._validate_polymer_rules(datapoint)
        return True

    def validate(self, datapoint: str):
        try:
            return self._validate(datapoint)
        except:
          return False


if __name__ == "__main__":
    validator = MolValidator()
    print("polymer valid >>", validator.validate("CC"))
    print("polymer invalid >>", validator.validate("CNCNFe"))
    validator = ProtValidator()
    print("protein valid >>", validator.validate("1MO8"))  # Valid pdb id
    print("protein invalid >>", validator.validate("XB123"))  # Invalid pdb id
    validator = PolymerValidator()
    print("polymer valid >>", validator.validate("[*]CC[*]"))  # Valid pdb id
    print("protein invalid >>", validator.validate("[*]CC"))  # Invalid pdb id
