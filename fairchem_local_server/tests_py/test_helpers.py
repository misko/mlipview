from fairchem_local_server import server


def test_health_endpoint_structure():
    h = server.health()
    assert h.get("status") == "ok"
    assert isinstance(h.get("device"), str)
    assert isinstance(h.get("model"), str)
