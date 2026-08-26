from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import sys
import tempfile
import types
import unittest

COMPARISON = Path(__file__).resolve().parents[1]


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, COMPARISON / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {filename}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class ComparisonFitResultsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        model = types.ModuleType("model")
        model.MASS_EDGES = tuple(round(1.0 + 0.1 * index, 1) for index in range(11))
        model.POLARIZATION = 0.3519
        model.POLARIZATION_ANGLE = 0.0
        model.bin_model = None
        model.bin_wave_models = None
        model.channel = None
        model.truth_model = None
        model.truth_wave_models = None
        sys.modules["model"] = model
        sys.modules.setdefault("laddu", types.ModuleType("laddu"))
        sys.modules.setdefault("numpy", types.ModuleType("numpy"))
        matplotlib = sys.modules.setdefault("matplotlib", types.ModuleType("matplotlib"))
        pyplot = types.ModuleType("matplotlib.pyplot")
        matplotlib.pyplot = pyplot
        sys.modules.setdefault("matplotlib.pyplot", pyplot)
        cls.prepare = load_module("comparison_prepare", "prepare.py")
        cls.plot_compare = load_module("comparison_plot_compare", "plot_compare.py")

    def test_amptools_floats_the_s_wave_magnitude_as_the_real_reference(self) -> None:
        config = self.prepare.amptools_config(0, Path("generated/amptools"))

        self.assertIn(
            "initialize ksks::PosRe::S0+ cartesian 1.0 0.0 real", config
        )
        self.assertNotIn(
            "initialize ksks::PosRe::S0+ cartesian 1.0 0.0 fixed", config
        )

    def test_amptools_reader_returns_the_complex_d_over_s_ratio(self) -> None:
        fit = """+++ Parameter Values and Errors +++
4
ksks::PosRe::S0+_re 3.0 0.0
ksks::PosRe::S0+_im 0.0 0.0
ksks::PosRe::D2+_re 0.0 0.0
ksks::PosRe::D2+_im 6.0 0.0
"""
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "zlm-bin-05.fit").write_text(fit, encoding="utf-8")
            points = self.plot_compare.read_amptools(directory)

        self.assertEqual(len(points), 1)
        self.assertAlmostEqual(points[0].magnitude, 2.0)
        self.assertAlmostEqual(points[0].phase, math.pi / 2.0)

    def test_pawian_floats_the_channel_normalization(self) -> None:
        defaults = """Reflectivity_f0(1500)_J0_M0_R+_K0Mag 0 0
Reflectivity_f2(1270)_J2_M2_R+_K0Mag 0 0
Reflectivity_f2(1270)_J2_M2_R+_K0Phi 0 0
channelType4protonKshort1Kshort2Scaling 1 0.01
"""
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            defaults_path = directory / "defaults.dat"
            parameters_path = directory / "start.dat"
            defaults_path.write_text(defaults, encoding="utf-8")
            fixed = self.prepare.prepare_parameters(defaults_path, parameters_path)
            config = self.prepare.pawian_config(0, directory, parameters_path, fixed)
            parameters = parameters_path.read_text(encoding="utf-8")

        self.assertIn(
            "Reflectivity_f2(1270)_J2_M2_R+_K0Mag\t0.5\t0.10000000000000001\t0\t5",
            parameters,
        )
        self.assertIn(
            "channelType4protonKshort1Kshort2Scaling\t1\t0.01\t0\t20",
            parameters,
        )
        self.assertNotIn("channelType4protonKshort1Kshort2Scaling", fixed)
        self.assertIn("fixPrimaryChannelScaling = false", config)
        self.assertIn("ratioMcToData = 1000", config)


if __name__ == "__main__":
    unittest.main()
