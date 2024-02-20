'''
from dash.testing.application_runners import import_app
import pytest

@pytest.fixture(scope='module')
def test_client():
    app = import_app("tests.app")

    return app
'''