from unittest.mock import patch

import pytest
from rdkit import Chem

from src.utils.curate import HeteroarenePull, RadicalPull, collate_descriptors

toyfeatures = {
    "actions": [
        {
            "result": {
                "descriptors": {
                    "atomic": {
                        "0_C": {"qn": 0},
                    },
                    "bonding": {
                        "0_C:1_H": {"dist": 1.4},
                        "1_C:0_C": {"dist": 1.4},
                    },
                    "molecular": {"ie_vertical": 8.4},
                }
            }
        }
    ]
}


@pytest.mark.parametrize(
    "track, all_descriptors, mol, expected_result",
    [
        (None, "None", None, FileNotFoundError),
        (None, None, None, 1),
    ],
)
def test___init___HeteroarenePull(
    track, all_descriptors, mol, expected_result
):
    """
    Test the initialisation of the HeteroarenePull class.
    """
    if expected_result == 1:
        with (
            patch("src.utils.curate.open") as patchopen,
            patch(
                "src.utils.curate.json.load",
                return_value={"actions": [{"result": {"descriptors": []}}]},
            ) as patchjson,
        ):
            HeteroarenePull(track, all_descriptors, mol)
            assert patchopen.call_count == 1
            assert patchjson.call_count == 1
    else:
        with pytest.raises(expected_result):
            HeteroarenePull(track, all_descriptors, mol)


@pytest.mark.parametrize(
    "track, all_descriptors, mol, _features, expected_result",
    [
        ([0], None, Chem.MolFromSmiles("C"), toyfeatures, True),
        ([0], None, Chem.MolFromSmiles("CC"), toyfeatures, True),
    ],
)
def test_pull_features_HeteroarenePull(
    track, all_descriptors, mol, _features, expected_result
):
    """
    Test the initialisation of the HeteroarenePull class.
    """
    if expected_result:
        with (
            patch("src.utils.curate.open") as patchopen,
            patch(
                "src.utils.curate.json.load", return_value=_features
            ) as patchjson,
        ):
            hp = HeteroarenePull(track, all_descriptors, mol)
            hp.pull_features()
            assert patchopen.call_count == 1
            assert patchjson.call_count == 1


@pytest.mark.parametrize(
    "track, all_descriptors, mol, expected_result",
    [
        (None, "None", None, FileNotFoundError),
        (None, None, None, 1),
    ],
)
def test___init___RadicalPull(track, all_descriptors, mol, expected_result):
    """
    Test the initialisation of the RadicalPull class.
    """
    if expected_result == 1:
        with (
            patch("src.utils.curate.open") as patchopen,
            patch(
                "src.utils.curate.json.load",
                return_value={"actions": [{"result": {"descriptors": []}}]},
            ) as patchjson,
        ):
            RadicalPull(track, all_descriptors, mol)
            assert patchopen.call_count == 1
            assert patchjson.call_count == 1
    else:
        with pytest.raises(expected_result):
            RadicalPull(track, all_descriptors, mol)


@pytest.mark.parametrize(
    "track, all_descriptors, mol, _features, expected_result",
    [
        ([0], None, Chem.MolFromSmiles("C*C"), toyfeatures, True),
    ],
)
def test_pull_features_RadicalPull(
    track, all_descriptors, mol, _features, expected_result
):
    """
    Test the initialisation of the RadicalPull class.
    """
    with (
        patch("src.utils.curate.open") as patchopen,
        patch(
            "src.utils.curate.json.load", return_value=_features
        ) as patchjson,
    ):
        rp = RadicalPull(track, all_descriptors, mol)
        rp.pull_features()
        assert patchopen.call_count == 1
        assert patchjson.call_count == 1


@pytest.mark.parametrize(
    "hetdescriptors, raddescriptors, expected_result",
    [
        (
            {0: {"descriptors": {"esp_var": 1, "dummy_desc": 1}}},
            {0: {"descriptors": {"esp_var": 1, "dummy_desc": 1}}},
            True,
        )
    ],
)
def test_collate_descriptors(hetdescriptors, raddescriptors, expected_result):
    """
    Test the collation of descriptors code.
    """
    descriptors = collate_descriptors(hetdescriptors, raddescriptors)
    assert len(descriptors[0]["descriptors"]) == 2, descriptors
