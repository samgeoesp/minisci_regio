import logging
import os
from contextlib import redirect_stdout

import qm_descriptors as qmd
from rdkit import Chem

from src.utils.curate import HeteroarenePull, RadicalPull, collate_descriptors
from src.utils.load import FindRXNCs
from src.utils.model import Model
from src.utils.visualise import plot_results


class Minisci:
    """
    Main class for running Minisci ML predictions using DFT derived
    descriptors from qm_descriptors package.
    """

    def __init__(
        self, workdir=None, hetarene=None, radical=None, plot=True
    ) -> None:
        self.descriptors = {"heteroarene": None, "radical": None}
        self.features = None
        self.radical = radical
        self.hetarene = hetarene
        self.root = os.getcwd()
        self.workdir = workdir
        self.plot = plot

    def _preamble(self) -> None:
        """
        Display information about MinisML.
        """
        sep = "=" * 78 + "\n"
        blank = "===" + " " * 72 + "===\n"
        with open("minisML.log", "a") as f:
            f.write(sep)
            f.write(blank)
            f.write(
                "==="
                + "MinisML - Minisci Regioselectivity Model".center(72)
                + "===\n"
            )
            f.write(blank)
            f.write("===" + "Written by Sam Espley".center(72) + "===\n")
            f.write(blank)
            f.write(sep)

    def _verifyinputs(self) -> None:
        """
        Load in the specified inputs and output any associated warnings/errors.
        """
        # Determine if the heteroarene is valid
        if self.hetarene is None:
            self.logger.error("No heteroarene provided.")
            raise MolInputError("No heteroarene provided.")
        elif Chem.MolFromSmiles(self.hetarene) is None:
            self.logger.error(
                f"Heteroarene ({self.hetarene}) is not a valid molecule."
            )
            raise MolInputError(
                f"Heteroarene ({self.hetarene}) is not a valid molecule."
            )

        # Deal with the radical and provide information to the user.
        if self.radical is None:
            self.logger.error("No radical provided.")
            raise MolInputError("No radical provided.")
        elif Chem.MolFromSmiles(self.radical) is None:
            self.logger.error(
                f"Radical ({self.radical}) is not a valid molecule."
            )
            raise MolInputError(
                f"Radical ({self.radical}) is not a valid molecule."
            )

        # Final check if valid structures loaded.
        self.logger.log(23, "Checking inputs ...")
        self.logger.info(f"Successfully loaded heteroarene - {self.hetarene}")
        self.logger.info(f"Successfully loaded radical - {self.radical}")

    def _generateworkdir(self) -> None:
        """
        Multiple calculations are needed - this is done by creating a
        working directory for each calculation as qm_descriptors needs an
        individual directory for each species.

        This is written to check for a 'calcs_minisci' directory before
        creating a job directory which is limited at 10,000 total.

        Within this working directory, 'heteroarene' and 'radical' directories
        are created.
        """
        os.chdir(self.workdir)
        try:
            os.makedirs("minisci", exist_ok=True)
            self.workdir = "{}/minisci".format(os.getcwd())
            os.chdir("minisci")
            # Make species specific directories.
            os.mkdir("het")
            os.mkdir("rad")
        except Exception as err:
            os.chdir(self.root)
            raise err

    def _generateradical(self) -> None:
        """
        Take the radical input and convert it to a SMILES containing the
        explicit hydrogens and the correct valency around the reaction centre.
        """
        mol = Chem.MolFromSmiles(self.radical)
        mol = Chem.AddHs(mol)
        rwmol = Chem.RWMol(mol)

        # Find the user defined radical centre.
        for atom in rwmol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomicNum(6)

        # Determine how many hydrogens are required for this radical species.
        add = 0
        for atom in rwmol.GetAtoms():
            if atom.GetSymbol() == "C" and atom.GetExplicitValence() != 4:
                if atom.GetExplicitValence() == 2:
                    add = 1
                    atomidx = atom.GetIdx()
                elif atom.GetExplicitValence() == 1:
                    add = 2
                    atomidx = atom.GetIdx()

        # Add the hydrogens to the radical species.
        if add == 0:
            pass
        elif add == 1:
            newH = rwmol.AddAtom(Chem.Atom(1))
            rwmol.AddBond(atomidx, newH, order=Chem.rdchem.BondType.SINGLE)
        elif add == 2:
            newH = rwmol.AddAtom(Chem.Atom(1))
            rwmol.AddBond(atomidx, newH, order=Chem.rdchem.BondType.SINGLE)
            newH = rwmol.AddAtom(Chem.Atom(1))
            rwmol.AddBond(atomidx, newH, order=Chem.rdchem.BondType.SINGLE)

        # Get updated molecule
        updated_mol = rwmol.GetMol()
        self.radical = Chem.MolToSmiles(updated_mol)
        self.radicalrc = [
            atom.GetIdx()
            for atom in updated_mol.GetAtoms()
            if atom.GetSymbol() == "C" and atom.GetExplicitValence() != 4
        ]

    def extraction(self) -> None:
        """
        Preprocess the inputs and run the calculations using qm_descriptors.

        The output from qm_descriptors is directed to a qm_descriptors.log file
        to prevent clutter to the terminal.

        All other logging will be stored in minisML.log.
        """
        # Create overall working directory.
        self._generateworkdir()

        # Set up the logger and create custom level for processes.
        self.logger = loggersetup()

        self._preamble()
        self.logger.info(f"Working Directory - {os.getcwd().split('/')[-1]}")

        # Verify inputs are valid.
        self._verifyinputs()

        # Clean the radical structure.
        self._generateradical()

        # Build mol objects.
        self.molrad = Chem.MolFromSmiles(self.radical)
        self.molhet = Chem.MolFromSmiles(self.hetarene)

        self.logger.log(23, "Finding reaction centres ...")

        # Determine reaction centres.
        frc = FindRXNCs(Chem.MolToMolBlock(self.molhet))
        frc.determine()
        frc.track["radical_rc"] = self.radicalrc

        # Run the qm_descriptors calculation and extraction.
        self.logger.log(23, "Starting calculations ...")
        for species, smi, mol in zip(
            ["het", "rad"],
            [self.hetarene, self.radical],
            [self.molhet, self.molrad],
        ):
            os.chdir(f"{self.workdir}/{species}")
            try:
                if species == "rad":
                    self.logger.info("      Radical calculation starting...")
                    # Run radical workflow.
                    redirect(
                        qmd.cdft.run_radical_cdft,
                        "qm_descriptors.log",
                        smiles=smi,
                        n_threads=16,
                        xc="b3lyp",
                        basis="6-311G(d,p)",
                    )
                    self.logger.info("      Radical calculation done!")
                    # Extract the descriptors.
                    radpull = RadicalPull(
                        frc.track["radical_rc"],
                        "workflow.json",
                        mol,
                    )
                    radpull.pull_features()
                    self.descriptors[species] = radpull.descriptors

                    os.chdir(self.workdir)
                else:
                    # Run standard workflow.
                    self.logger.info(
                        "      Heteroarene calculation starting..."
                    )
                    redirect(
                        qmd.cdft.run,
                        "qm_descriptors.log",
                        smiles=smi,
                        n_threads=16,
                        xc="b3lyp",
                        basis="6-311G(d,p)",
                    )
                    self.logger.info("      Heteroarene calculation done!")

                    # # Extract the descriptors.
                    hetpull = HeteroarenePull(
                        frc.track["heteroarene_rc"],
                        "workflow.json",
                        mol,
                    )
                    hetpull.pull_features()
                    self.descriptors[species] = hetpull.descriptors
                    os.chdir(self.workdir)
            except Exception as err:
                os.chdir(self.workdir)
                raise err
            os.chdir(self.workdir)

        self.features = collate_descriptors(
            self.descriptors["het"], self.descriptors["rad"]
        )
        os.chdir(self.root)

    def predict(self) -> None:
        """
        Running prediction on given molecules. Results are stored in
        model.resdata.

        A plot of the results will be generated in this block if self.plot=True
        was passed by the user.

        """
        # Make prediction.
        self.logger.log(23, "Making prediction ...")
        try:
            self.logger.info(f"Heteroarene --> {self.hetarene}")
            model = Model()
            model.predict(self.features)
            # Check whether the user wants to plot the results
            if self.plot:
                plot_results(self.hetarene, model.resdata, self.workdir)
            return model.resdata
        except Exception as err:
            self.logger.critical("Failed to predict for given reaction.")
            raise err


class MolInputError(Exception):
    """
    Custom error if either of the inputs are not suitable.
    """

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return f"Error with input structure - {self.message}"


def redirect(function, file, **kwargs):
    """
    Take the output of a given function and redirect it to a file instead.

    Args:
        function (callable): The function to redirect the output of.
        file (str): The name of the file to redirect output to.
        **kwargs: Additional keyword arguments for the function that you want
            to redirect the output for.
    """
    # Check the mode.
    mode = "a" if os.path.isfile(file) else "w"
    # Open file and run funciton with keyword arguments.
    with open(file, mode) as f:
        with redirect_stdout(f):
            function(**kwargs)


def loggersetup():
    """
    Define the logger and clear existing handles if the script is ran in a
    loop.
    """
    # Get logger.
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Check for existing handlers, if so remove.
    if logger.hasHandlers():
        logger.handlers.clear()

    # Create file handler.
    fh = logging.FileHandler("minisML.log")
    fh.setLevel(logging.INFO)

    # Adjust format.
    form = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%d-%m-%Y %H:%M:%S",
    )
    fh.setFormatter(form)
    logger.addHandler(fh)

    # Add custom logging level.
    logging.addLevelName(23, "PROCESS")

    return logger
