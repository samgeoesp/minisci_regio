import pytest
from rdkit import Chem

from src.utils.load import FindRXNCs

toypyr = Chem.MolToMolBlock(Chem.MolFromSmiles("c1cnccc1"))
toyquin = Chem.MolToMolBlock(Chem.MolFromSmiles("c12ncccc1cccc2"))


@pytest.mark.parametrize(
    "molblock, expected_result",
    [
        (toypyr, True),
    ],
)
def test___init___FindRXNCs(molblock, expected_result):
    """
    Test the __init__ of FindRXNCs.
    """
    frc = FindRXNCs(molblock)
    assert len(frc.track) == 1
    assert frc.mol.GetConformer()


@pytest.mark.parametrize(
    "molblock, expected_result",
    [
        (toypyr, True),
    ],
)
def test_find_CH_FindRXNCs(molblock, expected_result):
    """
    Test the find_CH function of FindRXNCs.
    """
    frc = FindRXNCs(molblock)
    frc.find_CH()
    assert frc.track["CH_centres"]


@pytest.mark.parametrize(
    "molblock, expected_result",
    [
        (toypyr, True),
    ],
)
def test_find_het_FindRXNCs(molblock, expected_result):
    """
    Test the find_het function of FindRXNCs.
    """
    frc = FindRXNCs(molblock)
    frc.find_CH()
    frc.find_het()
    assert len(frc.track["heteroarene_rc"]) == 5  # Undecorated pyridine


@pytest.mark.parametrize(
    "molblock, expected_result",
    [
        (toyquin, True),
    ],
)
def test_find_fused_het_FindRXNCs(molblock, expected_result):
    """
    Test the find_het function of FindRXNCs.
    """
    frc = FindRXNCs(molblock)
    frc.find_CH()
    frc.find_het()
    frc.find_fused_het()
    assert len(frc.track["heteroarene_rc"]) == 7  # Undecorated isoquinoline


@pytest.mark.parametrize(
    "molblock, expected_result",
    [
        (toypyr, ["greater", 5, 3]),
        (toyquin, ["equal", 7, 7]),
    ],
)
def test_remove_sym_FindRXNCs(molblock, expected_result):
    """
    Test the find_het function of FindRXNCs.
    """
    frc = FindRXNCs(molblock)
    frc.find_CH()
    frc.find_het()
    frc.find_fused_het()
    before = len(frc.track["heteroarene_rc"])
    frc.removesym()
    after = len(frc.track["heteroarene_rc"])
    if expected_result[0] == "greater":
        assert before > after
    else:
        assert before == after
    assert before == expected_result[1]
    assert after == expected_result[2]
