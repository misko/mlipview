from fairchem_local_server import server


def test_global_calculator_singleton():
    # The calculator and predict unit should be singletons
    assert server._CALCULATOR is server._CALCULATOR
    assert server._PREDICT_UNIT is server._PREDICT_UNIT


def test_health_endpoint_structure():
    h = server.health()
    assert h.get('status') == 'ok'
    assert isinstance(h.get('device'), str)
    assert isinstance(h.get('model'), str)