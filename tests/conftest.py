from pathlib import Path
import pytest
from mxspots.synthetic import get_cached_synthetic_frame, SyntheticFrame


@pytest.fixture(scope="session")
def test_data_dir() -> Path:
    return Path(__file__).resolve().parent / "data"


@pytest.fixture(scope="session")
def clean_frame(test_data_dir: Path) -> SyntheticFrame:
    return get_cached_synthetic_frame(str(test_data_dir / "clean.yaml"), max_spots=100)


@pytest.fixture(scope="session")
def insulin_frame(test_data_dir: Path) -> SyntheticFrame:
    return get_cached_synthetic_frame(str(test_data_dir / "insulin.yaml"), max_spots=100)


@pytest.fixture(scope="session")
def ice_frame(test_data_dir: Path) -> SyntheticFrame:
    return get_cached_synthetic_frame(str(test_data_dir / "ice.yaml"), max_spots=100)


@pytest.fixture(scope="session")
def lyso_split_frame(test_data_dir: Path) -> SyntheticFrame:
    return get_cached_synthetic_frame(str(test_data_dir / "lyso-split.yaml"), max_spots=100)
