import pytest

pytestmark = pytest.mark.skip(
    reason=(
        "Server/client cache removed; cache-specific tests are obsolete and " "skipped."
    )
)


def test_cache_removed_placeholder():
    assert True
