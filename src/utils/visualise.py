from __future__ import annotations

import tempfile

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


def _weight_to_pastel_colour(weight: float) -> tuple[float, float, float]:
    """Map a weight in [0, 1] from pastel blue (0) to pastel purple (1)."""
    w = max(0.0, min(1.0, float(weight)))
    blue = (0.67, 0.84, 0.95)
    purple = (0.545, 0.271, 0.820)
    return tuple(blue[i] + w * (purple[i] - blue[i]) for i in range(3))


WEIGHT_CMAP = LinearSegmentedColormap.from_list(
    "weight_pastel", [(0.67, 0.84, 0.95), (0.545, 0.271, 0.820)]
)


def plot_results(
    smiles: str, results: list[dict[str, str | float]], savepath: str
) -> None:
    """
    Take a SMILES string and plot the result onto the depiction of the given
    molecule.

    Arguments:
        smiles (str): The molecule in question.
        results (list(dict)): The list of reaction centres. Each entry in the list is
            a dictionary containing the keys "atom_idx" and "weight".
        savepath (str): Where to save the file.
    """
    # Generate the Mol object.
    mol = Chem.MolFromSmiles(smiles)

    # Loop through results and prepare per-atom colour highlights.
    highlight_atoms = []
    highlight_atom_colours = {}
    highlight_atom_radii = {}
    for res in results:
        atomidx = int(res["atom_idx"])
        weight = float(res["weight"])
        atom = mol.GetAtomWithIdx(atomidx)
        atom.SetProp("atomNote", f"{weight:.2f}")

        highlight_atoms.append(atomidx)
        highlight_atom_colours[atomidx] = _weight_to_pastel_colour(weight)
        highlight_atom_radii[atomidx] = 0.25

    drawer = rdMolDraw2D.MolDraw2DCairo(720, 720)
    opts = drawer.drawOptions()
    opts.fillHighlights = True
    opts.atomHighlightsAreCircles = True
    opts.annotationFontScale = 0.7
    opts.bondLineWidth = 3.0
    opts.baseFontSize = 3.0
    # opts.centreMoleculesBeforeDrawing = True

    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        mol,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=highlight_atom_colours,
        highlightAtomRadii=highlight_atom_radii,
    )
    drawer.FinishDrawing()

    with tempfile.NamedTemporaryFile(
        mode="wb", suffix=".png", delete=True
    ) as f:
        f.write(drawer.GetDrawingText())
        f.flush()

        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        img = plt.imread(f.name)
        ax.imshow(img)
        ax.axis("off")
        # ax.set_xlim(100, 550)

        norm = Normalize(vmin=0.0, vmax=1.0)
        sm = ScalarMappable(norm=norm, cmap=WEIGHT_CMAP)
        sm.set_array([])
        cbar = fig.colorbar(
            sm, ax=ax, fraction=0.1, aspect=15, location="bottom", pad=0.03
        )
        cbar.set_label("\nPredicted Values")

        plt.tight_layout()
        plt.savefig(f"{savepath}/prediction.png", dpi=300)
