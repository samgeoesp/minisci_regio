from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import HybridizationType as hbt


class FindRXNCs:

    def __init__(self, molblock: str):
        self.molblock = molblock
        self.nohmol = Chem.MolFromMolBlock(molblock)
        self.mol = Chem.AddHs(self.nohmol)
        AllChem.EmbedMolecule(self.mol)
        self.track = {"heteroarene_rc": []}

    def find_CH(self) -> None:
        """
        Determine every SP2 Hybridized ring C-H bond. To find the CH bonds it
        takes these SP2 reaction centres that are in a heteroarene ring and
        checks that there is 1 H atom attached.
        """
        # Find every C ring atom.
        self.track["ra"] = [
            atom.GetIdx()
            for atom in self.mol.GetAtoms()
            if atom.IsInRing() and atom.GetSymbol() == "C"
        ]
        # Get atoms that have SP2 hybridization.
        self.track["sp2_ra"] = [
            at
            for at in self.track["ra"]
            if self.mol.GetAtomWithIdx(at).GetHybridization() is hbt.SP2
        ]
        # Determine if atom has a C-H bond.
        self.track["CH_centres"] = []
        for atom in self.track["sp2_ra"]:
            nbs = [
                nb.GetIdx()
                for nb in self.mol.GetAtomWithIdx(atom).GetNeighbors()
            ]
            if (
                len(
                    [
                        nb
                        for nb in nbs
                        if self.mol.GetAtomWithIdx(nb).GetSymbol() == "H"
                    ]
                )
                == 1
            ):
                self.track["CH_centres"].append(atom)

    def find_het(self) -> None:
        """
        Most reactants fall within this category. They are an aromatic ring and
        have at least one heteroatom.
        """
        # Pull all the aromatic rings.
        self.track["all_rings"] = [
            r for r in self.mol.GetRingInfo().AtomRings()
        ]
        for ring in self.track["all_rings"]:
            if all(self.mol.GetAtomWithIdx(at).GetIsAromatic() for at in ring):
                if any(
                    self.mol.GetAtomWithIdx(at).GetSymbol() != "C"
                    for at in ring
                ):
                    self.track["heteroarene_rc"].extend(
                        set(self.track["CH_centres"]) & set(ring)
                    )

    def find_fused_het(self) -> None:
        """
        find_het only looks for single rings that have a heteroatom. This
        instead finds every instance where a heteroaromatic ring is fused to
        another. It will also find rings that could have tautomeric rings that
        are aromatic and therefore, could react.
        """
        self.track["fused_rings"] = []
        self.track["fused_atoms"] = []
        # Loop through all rings in the system.
        for id1, ring1 in enumerate(self.track["all_rings"]):
            # Check aromaticity and heteroarene.
            arom1 = all(
                self.mol.GetAtomWithIdx(at).GetIsAromatic() for at in ring1
            ) and any(
                self.mol.GetAtomWithIdx(at).GetSymbol() != "C" for at in ring1
            )
            # Cross check all rings.
            for id2 in range(id1 + 1, len(self.track["all_rings"])):
                ring2 = self.track["all_rings"][id2]
                arom2 = all(
                    self.mol.GetAtomWithIdx(at).GetIsAromatic() for at in ring2
                ) and any(
                    self.mol.GetAtomWithIdx(at).GetSymbol() != "C"
                    for at in ring2
                )
                if len(set(ring1) & set(ring2)) != 1:
                    if arom1 or arom2:
                        self.track["fused_rings"].append(
                            list(set(ring1) ^ set(ring2))
                        )
                        self.track["fused_atoms"].extend(
                            list(set(ring1) & set(ring2))
                        )

        # Loop through all fused rings and determine if SP2
        for fused in self.track["fused_rings"]:
            if all(
                self.mol.GetAtomWithIdx(at).GetHybridization() is hbt.SP2
                for at in fused
            ):
                self.track["heteroarene_rc"].extend(
                    set(self.track["CH_centres"]) & set(fused)
                )
                self.track["heteroarene_rc"] = list(
                    set(self.track["heteroarene_rc"])
                )

    def removesym(self) -> None:
        """
        Evaluate the possible reaction centres to determine if they have
        symmetrical equivalents. If so, remove the lowest integer value.
        """
        # Create rank for molecule.
        rank = Chem.CanonicalRankAtoms(self.nohmol, breakTies=False)
        syms = {}
        for idx, r in enumerate(rank):
            syms.setdefault(r, []).append(idx)

        rem = []
        for sym in syms.values():
            rem += sym[: len(sym) - 1]

        # Remove the symmetrical sites.
        for symat in rem:
            if symat in self.track["heteroarene_rc"]:
                self.track["heteroarene_rc"].remove(symat)

    def determine(self):
        self.find_CH()
        self.find_het()
        self.find_fused_het()
        self.removesym()
