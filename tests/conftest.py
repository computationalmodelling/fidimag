import os
import pytest
import tempfile

@pytest.fixture(scope="session", autouse=True)
def config():
    tmpdir = tempfile.mkdtemp(suffix='fidimag-tests')
    os.chdir(tmpdir)
