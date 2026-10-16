# Data

This repository contains the reactions used to train and test the model.


`literature_reactions.json` and `patent_reactions.json` are the reactions used to train the model. `experimental_reactions.json` is the collection of reactions in the HTE set.

Each reaction has a corresponding dictionary with `heteroarene`, `radical`, and `product` keys. The `patent_reactions.json` has the additional `patent_id` key whereas the `literature_reactions.json` has the `paper` key containing the DOI in which the reaction information was pulled from.

For every reaction, there is a `heteroarene_rc` which corresponds to the ground truth atom position where the reaction occurs on the heteroarene. Additionally, there is `other_heteroarene_rc` which denotes additional possible reactivity sites in that given molecule.