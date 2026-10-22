# Data

This repository contains the reactions used to train and test the model.


`literature_reactions.json` and `patent_reactions.json` are the reactions used to train the model. `hte_reactions.json` is the collection of reactions in the HTE set.

Each reaction has a corresponding dictionary with `heteroarene`, `radical`, and `product` keys. The `patent_reactions.json` has the additional `patent_id` key whereas the `literature_reactions.json` has the `paper` key containing the DOI in which the reaction information was pulled from.

For every reaction, there is a `heteroarene_rc` which corresponds to the ground truth atom position where the reaction occurs on the heteroarene. Additionally, there is `other_heteroarene_rc` which denotes additional possible reactivity sites in that given molecule.

### HTE

The reactions forming the HTE set curated for this work can be found in  `hte_reactions.json`. Within this file there is two additional columns: `observed_regioisomer` and `observed_ratio`. The former is all regioisomers observed in running the experiments and the latter is the ratio between the observed regioisomers.

### Example

Below is an annotated example from the HTE set:

```
 "600019": {
    "heteroarene": "ClC1=C(C=CN2C)C2=NC=N1", 
    "radical": "[H]C1C([H])([H])N(C(=O)OC(C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H])C1([H])[H]",
    "product": "ClC1=C(C=C(C(C2)CN2C(OC(C)(C)C)=O)N3C)C3=NC=N1",
    "heteroarene_rc": 4,
    "other_heteroarene_rc": [
      3,
      9
    ],
    "observed_regioisomers": [
      "ClC1=C(C=C(C(C2)CN2C(OC(C)(C)C)=O)N3C)C3=NC=N1",
      "ClC1=C(C=CN2C)C2=NC(C3CN(C3)C(OC(C)(C)C)=O)=N1"
    ],
    "observed_ratio": [
      0.75,
      0.25
    ]
```

| Key | Explanation |
| --- | --- |
| `heteroarene` | The SMILES for the heteroarene. This SMILES string corresponds to the atom numberings found in `heteroarene_rc` and `other_heteroarene_rc`. |
| `radical` | The SMILES of the radical structure. |
| `product` | The SMILES of the major regioisomer. This SMILES is also in `observed_regioisomers` as the first entry in the list. |
| `heteroarene_rc` | The atom number of the reaction centre on the SMILES in `heteroarene`. |
| `other_heteroarene_rc` | The atom number of other possible reaction centres on the SMILES in `heteroarene`. |
| `observed_regioisomers` | The different experimentally observed regioisomers. Entry `[0]` corresponds to the same SMILES ad `product`. |
| `observed_ratio` | The experimentally observed ratios for the regioisomers in `observed_regioisomers`. Indexes match between these two columns. This ratio is shown as a sum normalisation. |