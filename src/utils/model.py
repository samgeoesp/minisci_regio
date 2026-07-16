import os
import pickle

import numpy as np


class Model:

    def __init__(self):
        # Load in the model using filename path.
        pth = os.path.dirname(__file__)
        with open(f"{pth}/model.pkl", "rb") as m:
            self.model = pickle.load(m)
        # Extract names and build results list.
        self.names = self.model[:-1].get_feature_names_out()
        self.indices = [int(name[1:]) for name in self.names]
        self.resdata = []

    def predict(self, descriptors: dict) -> None:
        """
        Make the prediction on the molecule.

        Args:
            descriptors (dict): The descriptors for each reaction site on a
                molecule.
        """
        preds = {}
        # Loop through each reaction site in the descriptors.
        for rxnsite, d in descriptors.items():

            # Check length of extracted descriptors.
            if len(list(d["descriptors"].values())) != len(self.names):
                raise ValueError(
                    f"Incorrect descriptor length for reaction site {rxnsite}."
                )

            Xarr = np.array([i for i in d["descriptors"].values()]).reshape(
                1, -1
            )
            # Make prediciton.
            ypred = self.model.predict(Xarr)[0]
            preds[rxnsite.split("_")[1]] = ypred

        # Save the predictions.
        for rxnsite, ypred in preds.items():
            self.resdata.append({"atom_idx": rxnsite, "weight": ypred})
