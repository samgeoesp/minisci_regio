from unittest.mock import patch

import pytest

from src.utils.model import Model

toydecr = {"0_0": {"descriptors": {"het_dummy": 22, "rad_dummy": 23}}}


@pytest.mark.parametrize("expected_result", [(True)])
def test___init___Model(expected_result):
    """
    Test the __init__ of Model.
    """

    with (
        patch("src.utils.model.open"),
        patch("src.utils.model.pickle.load"),
    ):
        model = Model()
        assert model.model
        assert model.resdata == []


@pytest.mark.parametrize("descriptors, expected_result", [(toydecr, True)])
def test_predict_Model(descriptors, expected_result):
    """
    Test the __init__ of Model.
    """

    with (
        patch("src.utils.model.open"),
        patch("src.utils.model.pickle.load"),
    ):
        model = Model()
        before = len(model.resdata)
        model.predict(descriptors)
        after = len(model.resdata)
        assert before < after
