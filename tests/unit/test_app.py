from main import app

def test_home_page():
    """
    Given a Flask app configured for testing
    when the '/' page is requested (GET)
    Then check that response is valid
    :return:
    """
    with app.test_client() as test_client:
        response = test_client.get('/')
        assert response.status_code == 200
        assert b"Needleman-Wunsch" in response.data

