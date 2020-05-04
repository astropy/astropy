from .. import tests


def pytest_configure(config):
    tests.RUNNING = True
    print('SET RUNNING TO TRUE')


def pytest_unconfigure(config):
    tests.RUNNING = False
    print('SET RUNNING TO FALSE')
