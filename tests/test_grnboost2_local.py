import numpy as np

from edge_causality.grnboost2_local import EarlyStopMonitor


def test_early_stop_monitor_waits_for_full_window() -> None:
    class Model:
        oob_improvement_ = np.full(30, -1.0)

    monitor = EarlyStopMonitor(window_length=25)
    assert monitor(23, Model(), None) is False
    assert monitor(24, Model(), None) is True
