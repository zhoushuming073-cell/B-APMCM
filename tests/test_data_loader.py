import pandas as pd

from apmcm.common.paths import ARCHIVE, DATA, PROCESS_INPUT, config
from apmcm.data.verify import verify


def test_corrected_archive_and_workbook_hashes():
    assert verify() == 10


def test_process_input_agrees_with_corrected_audit():
    settings = config("data.yaml")
    raw = pd.read_csv(PROCESS_INPUT)
    audit = pd.read_csv(DATA / "metadata/audit.csv")
    assert len(raw) == settings["expected_records"]
    assert raw.groupby("heat").size().to_dict() == audit.set_index("heat").rows.to_dict()
    assert raw[raw.attachment == 1][["CO2", "CO"]].max().tolist() == [21.74, 59.04]
