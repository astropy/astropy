def assert_celestial_component(component, index):
    assert component[:2] == ("celestial", index)
    assert callable(component[2])
