from unittest.mock import MagicMock, mock_open, patch

import pytest

from src.minisci import Minisci, MolInputError, loggersetup, redirect
from src.utils.load import FindRXNCs


@pytest.mark.parametrize(
    "inputs, expected_workdir, expected_hetarene, expected_radical, expected_plot",
    [
        (
            {
                "workdir": "example",
                "hetarene": "c1cnccc1",
                "radical": "C*C",
                "plot": False,
            },
            "example",
            "c1cnccc1",
            "C*C",
            False,
        ),
        (
            {
                "workdir": None,
                "hetarene": None,
                "radical": None,
                "plot": True,
            },
            None,
            None,
            None,
            True,
        ),
    ],
)
def test__init__Minisci(
    inputs,
    expected_workdir,
    expected_hetarene,
    expected_radical,
    expected_plot,
):
    """
    Code to test the initialization of the Minisci class.
    """
    mock_logger = MagicMock()
    with patch("src.minisci.loggersetup", return_value=mock_logger):
        mml = Minisci(**inputs)
        assert mml.workdir == expected_workdir
        assert mml.hetarene == expected_hetarene
        assert mml.radical == expected_radical
        assert mml.plot == expected_plot
        assert mml.descriptors == {"heteroarene": None, "radical": None}
        assert mml.features is None


@pytest.mark.parametrize(
    "inputs, expected_result",
    [
        (
            {
                "workdir": "example",
                "hetarene": "c1cnccc1",
                "radical": "C*C",
                "plot": False,
            },
            True,
        ),
        (
            {
                "workdir": "example",
                "hetarene": "c1cnccc1",
                "radical": "C*C",
                "plot": False,
            },
            False,
        ),
    ],
)
def test_extraction_Minisci(inputs, expected_result):
    """
    Code to test the extraction of descriptors for the Minisci predictions.
    """
    mock_logger = MagicMock()

    with (
        patch.object(Minisci, "_generateworkdir"),
        patch.object(Minisci, "_preamble"),
        patch.object(Minisci, "_verifyinputs"),
        patch("src.minisci.loggersetup", return_value=mock_logger),
        patch.object(
            FindRXNCs,
            "determine",
            return_value={"radical_rc": 0, "heteroarene_rc": 0},
        ),
        patch("src.minisci.RadicalPull"),
        patch("src.minisci.os.chdir"),
        patch("src.minisci.redirect"),
    ):
        mml = Minisci(**inputs)
        # Use HeteroarenePull as the side effect creator.
        if expected_result:
            with patch("src.minisci.HeteroarenePull"):
                mml.extraction()
                assert mml
        else:
            with patch(
                "src.minisci.HeteroarenePull", side_effect=Exception("Error")
            ):
                with pytest.raises(Exception):
                    mml.extraction()


@pytest.mark.parametrize(
    "inputs, expected_result",
    [
        (
            {
                "workdir": "example",
                "hetarene": "c1cnccc1",
                "radical": "C*C",
                "plot": False,
            },
            True,
        ),
        (
            {
                "workdir": "example",
                "hetarene": "c1cnccc1",
                "radical": "C*C",
                "plot": False,
            },
            False,
        ),
    ],
)
def test_predict_Minisci(inputs, expected_result):
    """
    Code to test the prediction function for the Minisci predictions.
    """
    mock_logger = MagicMock()
    with (
        patch("src.minisci.loggersetup", return_value=mock_logger),
        patch.object(Minisci, "extraction"),
        patch("src.minisci.Model"),
    ):
        mml = Minisci(**inputs)
        mml.logger = mock_logger
        if expected_result:
            with patch("src.minisci.Model"):
                mml.predict()
                assert mml.hetarene
        else:
            with patch("src.minisci.Model", side_effect=Exception("Error")):
                with pytest.raises(Exception):
                    mml.predict()


@pytest.mark.parametrize(
    "inputs, expected_result",
    [
        (
            {
                "workdir": "example",
                "hetarene": "c1cnccc1",
                "radical": "C*C",
                "plot": False,
            },
            True,
        )
    ],
)
def test__preamble_Minisci(inputs, expected_result):
    """
    Code to test the preamble print for the Minisci predictions.
    """
    mockopen = mock_open()
    with (patch("src.minisci.open", mockopen),):
        mml = Minisci(**inputs)
        mml._preamble()
        assert mockopen().write.call_count == 7


@pytest.mark.parametrize(
    "inputs, expected_result",
    [
        (
            {
                "workdir": "example",
                "hetarene": None,
                "radical": "C*C",
                "plot": False,
            },
            MolInputError,
        ),
        (
            {
                "workdir": "example",
                "hetarene": "None",
                "radical": "C*C",
                "plot": False,
            },
            MolInputError,
        ),
        (
            {
                "workdir": "example",
                "hetarene": "c1cnccc1",
                "radical": None,
                "plot": False,
            },
            MolInputError,
        ),
        (
            {
                "workdir": "example",
                "hetarene": "c1cnccc1",
                "radical": "None",
                "plot": False,
            },
            MolInputError,
        ),
        (
            {
                "workdir": "example",
                "hetarene": "c1cnccc1",
                "radical": "C*C",
                "plot": False,
            },
            True,
        ),
    ],
)
def test__verifyinputs_Minisci(inputs, expected_result):
    """
    Code to test the verification of inputs for the Minisci predictions.
    """
    mock_logger = MagicMock()
    with (patch("src.minisci.loggersetup", return_value=mock_logger),):
        mml = Minisci(**inputs)
        mml.logger = mock_logger
        if expected_result is MolInputError:
            with pytest.raises(expected_result):
                mml._verifyinputs()
            assert mock_logger.error.call_count == 1
        else:
            mml._verifyinputs()
            assert mock_logger.log.call_count == 1
            assert mock_logger.info.call_count == 2


@pytest.mark.parametrize(
    "inputs, expected_result",
    [
        (
            {
                "workdir": "example",
                "hetarene": "c1cnccc1",
                "radical": "C*C",
                "plot": False,
            },
            True,
        ),
        (
            {
                "workdir": "example",
                "hetarene": "c1cnccc1",
                "radical": "C*C",
                "plot": False,
            },
            False,
        ),
    ],
)
def test__generateworkdir_Minisci(inputs, expected_result):
    """
    Code to test the generation of working directories for the Minisci
    predictions.
    """
    mockopen = mock_open()
    with (
        patch("src.minisci.open", mockopen),
        patch("src.minisci.os.getcwd") as patchgetcwd,
        patch("src.minisci.os.chdir") as patchchdir,
        patch("src.minisci.os.mkdir") as patchmkdir,
    ):
        mml = Minisci(**inputs)
        if expected_result:
            with patch("src.minisci.os.makedirs") as patchmakedirs:
                mml._generateworkdir()
                assert patchmakedirs.call_count == 1
                assert patchgetcwd.call_count == 2  # Once in __init__...
                assert patchchdir.call_count == 2
                assert patchmkdir.call_count == 2
        else:
            with patch("src.minisci.os.makedirs", side_effect=Exception):
                with pytest.raises(Exception):
                    mml._generateworkdir()
                assert patchchdir.call_count == 2


@pytest.mark.parametrize(
    "inputs, expected_result",
    [
        (
            {
                "workdir": "example",
                "hetarene": "c1cnccc1",
                "radical": "c1c*ccc1",
                "plot": False,
            },
            2,
        ),
        (
            {
                "workdir": "example",
                "hetarene": "c1cnccc1",
                "radical": "C*C",
                "plot": False,
            },
            1,
        ),
        (
            {
                "workdir": "example",
                "hetarene": "c1cnccc1",
                "radical": "*C",
                "plot": False,
            },
            0,
        ),
    ],
)
def test__generateradical_Minisci(inputs, expected_result):
    """
    Code to test the generation of the radical for the Minisci predictions.

    The different tests go through the different types of radicals based on
    hybridisation.
    """
    mml = Minisci(**inputs)
    mml._generateradical()
    assert mml.radicalrc == [expected_result]


@pytest.mark.parametrize(
    "message, expected_result",
    [
        ("This is an error", True),
    ],
)
def test___str___MolInputError(message, expected_result):
    """
    Code to test the MolInputError class.
    """
    molinperr = MolInputError(message)
    assert molinperr.__str__() == f"Error with input structure - {message}"


@pytest.mark.parametrize(
    "file_exists, expected_mode",
    [
        (False, "w"),
        (True, "a"),
    ],
)
def test_redirect(file_exists, expected_mode):
    """
    Code to test the redirect function for redirecting stdout to a file.
    """
    mock_function = MagicMock()
    mockopen = mock_open()

    with (
        patch("src.minisci.os.path.isfile", return_value=file_exists),
        patch("src.minisci.open", mockopen),
    ):
        redirect(mock_function, "test.log", arg1="value1", arg2="value2")

        # Check file was opened with correct mode
        mockopen.assert_called_once_with("test.log", expected_mode)
        # Check function was called with kwargs
        mock_function.assert_called_once_with(arg1="value1", arg2="value2")


@pytest.mark.parametrize(
    "has_handlers, expected_clear_count",
    [
        (False, 0),
        (True, 1),
    ],
)
def test_loggersetup(has_handlers, expected_clear_count):
    """
    Code to test the loggersetup function for creating the logger.
    """
    mock_logger = MagicMock()
    mock_logger.hasHandlers.return_value = has_handlers
    mock_handler = MagicMock()

    with (
        patch("src.minisci.logging.getLogger", return_value=mock_logger),
        patch("src.minisci.logging.FileHandler", return_value=mock_handler),
        patch("src.minisci.logging.Formatter"),
        patch("src.minisci.logging.addLevelName") as mock_addlevel,
    ):
        result = loggersetup()

        # Check logger was configured correctly
        mock_logger.setLevel.assert_called_once()
        mock_logger.hasHandlers.assert_called_once()
        assert mock_logger.handlers.clear.call_count == expected_clear_count

        # Check handler was added
        mock_logger.addHandler.assert_called_once_with(mock_handler)

        # Check custom level was added
        mock_addlevel.assert_called_once_with(23, "PROCESS")

        # Check logger is returned
        assert result == mock_logger
