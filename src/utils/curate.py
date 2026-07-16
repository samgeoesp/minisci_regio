import json

from rdkit import Chem
from rdkit.Chem import Mol


class HeteroarenePull:

    def __init__(self, track: list, alldescriptors: str, mol: Mol):
        self.track = track
        self.mol = mol
        self.descriptors = {}
        # Load in extracted descriptors.
        with open(alldescriptors) as f:
            self.alldescriptors = json.load(f)["actions"][-1]["result"][
                "descriptors"
            ]

    def pull_features(self) -> None:
        """
        Pull the correct features from the data based upon previous analysis.

        In order to write the code in a more readable manner, it has three
        loops through the reactions. This could be resolved but the overall
        computational cost and efficiency at this point of the extraction
        is not key.

        First it extracts the atom based features, followed by the bond
        based features. It will then extract the target before saving.
        """
        # Atom Descriptors.
        # Loop through every reaction centre.
        atom_descriptors = {}
        for rxn_site in self.track:
            tag = f"0_{rxn_site}"
            for atom_feat in self.alldescriptors["atomic"].keys():
                if atom_feat.startswith(f"{rxn_site}_"):
                    # Check that the atom is a Carbon.
                    assert atom_feat.endswith("C")
                    atom_descriptors[tag] = self.alldescriptors["atomic"][
                        atom_feat
                    ]
        # Bond Descriptors.
        # Loop through only reactions pulled in first step.
        bond_descriptors = {}
        for tag in atom_descriptors.keys():
            bf = {}
            reaction_tag = tag.split("_")
            rc = reaction_tag[-1]
            hm = Chem.AddHs(self.mol)
            atom = hm.GetAtomWithIdx(int(rc))
            bonded_atoms = [neigh.GetIdx() for neigh in atom.GetNeighbors()]
            atom_symb = [
                hm.GetAtomWithIdx(atom).GetSymbol() for atom in bonded_atoms
            ]
            combined = [f"{i}_{j}" for i, j in zip(bonded_atoms, atom_symb)]
            # Loop through connected atoms.
            rxn_features = self.alldescriptors
            pre_sort, distances = {}, []
            for second in combined:
                search1 = f"{rc}_C:{second}"
                search2 = f"{second}:{rc}_C"
                if search1 in rxn_features["bonding"].keys():
                    bond_f = rxn_features["bonding"][search1]
                    bond_f["search"] = search1
                    pre_sort[rxn_features["bonding"][search1]["dist"]] = bond_f
                    distances.append(rxn_features["bonding"][search1]["dist"])
                elif search2 in rxn_features["bonding"].keys():
                    bond_f = rxn_features["bonding"][search2]
                    bond_f["search"] = search2
                    pre_sort[rxn_features["bonding"][search2]["dist"]] = bond_f
                    distances.append(rxn_features["bonding"][search2]["dist"])

            # First pull the hydrogen bond.
            for dist in distances:
                if "_H" in pre_sort[dist]["search"]:
                    for feat in pre_sort[dist].keys():
                        bf[f"C_H_bond_{feat}"] = pre_sort[dist][feat]
                    distances.remove(dist)

            # Now add other bond descriptors based on distances.
            for ind, dist in enumerate(sorted(distances)):
                for feat in pre_sort[dist].keys():
                    bf[f"C_{ind}_bond_{feat}"] = pre_sort[dist][feat]

            # Remove unwanted columns (bond_search column).
            for col in [
                "C_H_bond_search",
                "C_0_bond_search",
                "C_1_bond_search",
            ]:
                if col in bf.keys():
                    del bf[col]

            # Populate dictionary.
            bond_descriptors[tag] = bf

        # Combine descriptors and return.
        for tag in bond_descriptors.keys():
            self.descriptors[tag] = {
                "descriptors": {
                    **atom_descriptors[tag],
                    **bond_descriptors[tag],
                }
            }


class RadicalPull:

    def __init__(self, track: list, alldescriptors: str, mol: Mol):
        self.track = track
        self.mol = mol
        self.descriptors = {}
        self.features = {}
        # Load in extracted descriptors.
        with open(alldescriptors) as f:
            self.alldescriptors = json.load(f)["actions"][-1]["result"][
                "descriptors"
            ]

    def pull_features(self) -> None:
        """
        Pull the correct features from the data based upon previous analysis.

        Functionality is less complex than for heteroarene species due to only
        pulling the atom based descriptors.
        """
        # Loop through the reaction sites.
        atom_descriptors = {}
        rxn_features = self.alldescriptors
        for rxn_site in self.track:
            tag = f"1_{rxn_site}"
            # Pull the atomic and molecular descriptors for the rxn_site.
            for atom_feat in self.alldescriptors["atomic"].keys():
                if atom_feat.startswith(f"{rxn_site}_"):
                    atom_descriptors[tag] = {
                        **rxn_features["atomic"][atom_feat],
                        **rxn_features["molecular"],
                    }
            # Combine descriptors and return.
            for tag in atom_descriptors.keys():
                self.descriptors[tag] = {
                    "descriptors": {
                        **atom_descriptors[tag],
                    }
                }


def collate_descriptors(hetdescriptors: dict, raddescriptors: dict) -> dict:
    """
    Combine the calculated heteroarene descriptors with the descriptors for
    the radical species.

    Args:
        hetdescriptors (dict): The heteroarene descriptors.
        raddescriptors (dict): The radical descriptors.
    Returns:
        features (dict): This is all the collated descriptors for the entire
            system.
    """
    # Unwanted features.
    features = {}
    remove = [
        "het_esp_var",
        "rad_alie_min",
        "rad_alie_max",
        "rad_alie_avg",
        "rad_alie_var",
        "rad_lea_min",
        "rad_lea_max",
        "rad_lea_avg",
        "rad_lea_var",
        "rad_esp_avg",
        "rad_esp_var",
        "rad_area_tot",
        "rad_area_pos",
        "rad_area_neg",
        "rad_esp_min",
        "rad_esp_max",
    ]

    # Create the combined descriptors.
    for rxncentre, d in raddescriptors.items():
        rad = {f"rad_{desc}": val for desc, val in d["descriptors"].items()}
    for rxncentre, d in hetdescriptors.items():
        het = {f"het_{desc}": val for desc, val in d["descriptors"].items()}
        features[rxncentre] = {"descriptors": {**het, **rad}}
        # Remove unwanted features.
        for frem in remove:
            if frem in features[rxncentre]["descriptors"].keys():
                del features[rxncentre]["descriptors"][frem]
    return features
