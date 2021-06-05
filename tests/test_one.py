from dash.testing.application_runners import import_app

def test_one(dash_duo):
    app = import_app("tests.app")
    dash_duo.start_server(app)
    dash_duo.wait_for_text_to_equal("h2", "Needleman-Wunsch, Smith-Waterman and Entrez", timeout=4)
    #assert dash_duo.find_element("h3").text == "Objective"
    assert dash_duo.find_element("h3").text == "PAM BLOSUM + NCBI"
    assert dash_duo.get_logs() == [], "Browser console should contain no errors"

    return None