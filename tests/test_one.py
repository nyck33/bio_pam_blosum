from dash.testing.application_runners import import_app

def test_one(dash_duo):
    app = import_app("tests.test_app")
    dash_duo.start_server(app)
    dash_duo.wait_for_text_to_equal("h2", "Needleman-Wunsch", timeout=4)

