from src.minisci import Minisci
from src.utils.visualise import plot_results

# Extract features.
mml = Minisci(
    workdir=".",
    # hetarene="FC1=CC2=CC=C(Br)C=C2N=C1",
    # radical="*N(C(OC(C)(C)C)=O)[H]",
    hetarene="c1cnccc1",
    radical="C*C",
    plot=False,
)
mml.extraction()

# Make prediction.
data = mml.predict()

# Visualise result.
plot_results(
    # smiles="FC1=CC2=CC=C(Br)C=C2N=C1",
    smiles="c1cnccc1",
    results=data,
    savepath=".",
)
