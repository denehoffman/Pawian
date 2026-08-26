from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

COMPARISON = Path(__file__).resolve().parents[1]
REPOSITORY = COMPARISON.parents[2]
REFERENCE = COMPARISON / "reference"
SCRIPT = COMPARISON / "amptools_reference.py"
AMPTOOLS_REVISION = json.loads(
    (REFERENCE / "compatibility-artifacts.json").read_text(encoding="utf-8")
)["reference"]["revision"]


class AmpToolsReferenceWorkflowTest(unittest.TestCase):
    def valid_result(self) -> dict[str, object]:
        sample_path = REFERENCE / "events.json"
        sample = json.loads(sample_path.read_text(encoding="utf-8"))
        events = []
        for source in sample["data"]:
            particles = [
                source["particles"][name] for name in sample["particle_order"]
            ]
            final_state_sum = [
                sum(momentum[component] for momentum in particles[1:])
                for component in range(4)
            ]
            events.append(
                {
                    "event_id": source["event_id"],
                    "weight": source["weight"],
                    "particles": particles,
                    "final_state_sum": final_state_sum,
                    "intensity": 1.0 / (8.0 * math.pi),
                }
            )
        return {
            "format_version": 1,
            "framework": {
                "name": "AmpTools",
                "revision": AMPTOOLS_REVISION,
            },
            "sample_sha256": hashlib.sha256(
                sample_path.read_bytes()
            ).hexdigest(),
            "particle_order": sample["particle_order"],
            "polarization_fraction": 0.0,
            "normalization": {
                "generated_before_evaluation": True,
                "generated_events": len(sample["normalization"]),
                "accepted_weight_sum": sum(
                    event["weight"] for event in sample["normalization"]
                ),
            },
            "normalization_source": {
                stream: {
                    "event_count": len(sample["normalization"]),
                    "events": [
                        {
                            "event_id": event["event_id"],
                            "weight": event["weight"],
                        }
                        for event in sample["normalization"]
                    ],
                }
                for stream in ("generated", "accepted")
            },
            "event_count": len(events),
            "weight_sum": sum(event["weight"] for event in sample["data"]),
            "events": events,
        }

    def test_compatibility_artifacts_match_recorded_digests(self) -> None:
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "verify-artifacts",
                "--manifest",
                str(REFERENCE / "compatibility-artifacts.json"),
                "--root",
                str(REPOSITORY),
            ],
            check=True,
            cwd=COMPARISON,
            text=True,
            capture_output=True,
        )

        self.assertIn("verified 4 compatibility artifacts", completed.stdout)

    def test_top_level_recipe_does_not_require_nix(self) -> None:
        completed = subprocess.run(
            ["just", "--dry-run", "amptools-reference"],
            check=True,
            cwd=REPOSITORY,
            text=True,
            capture_output=True,
        )

        rendered_recipe = completed.stdout + completed.stderr
        self.assertNotIn("nix", rendered_recipe)

    def test_build_comparison_recipe_does_not_require_nix(self) -> None:
        completed = subprocess.run(
            ["just", "--dry-run", "build-comparison"],
            check=True,
            cwd=REPOSITORY,
            text=True,
            capture_output=True,
        )

        rendered_recipe = completed.stdout + completed.stderr
        self.assertNotIn("nix", rendered_recipe)

    def test_modified_compatibility_artifact_fails_visibly(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in ("Zlm.cc", "Zlm.h", "wignerD.cc", "wignerD.h"):
                destination = root / "Examples/gammap/comparison" / name
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(COMPARISON / name, destination)
            with (root / "Examples/gammap/comparison/Zlm.cc").open(
                "a", encoding="utf-8"
            ) as artifact:
                artifact.write("\n// accidental edit\n")

            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "verify-artifacts",
                    "--manifest",
                    str(REFERENCE / "compatibility-artifacts.json"),
                    "--root",
                    str(root),
                ],
                check=False,
                cwd=COMPARISON,
                text=True,
                capture_output=True,
            )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Zlm.cc", completed.stderr)
        self.assertIn("digest mismatch", completed.stderr)

    def test_prepare_writes_auditable_fixed_s_wave_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "prepare",
                    "--sample",
                    str(REFERENCE / "events.json"),
                    "--output",
                    str(output),
                ],
                check=True,
                cwd=COMPARISON,
                text=True,
                capture_output=True,
            )

            metadata = json.loads(
                (output / "input-metadata.json").read_text(encoding="utf-8")
            )
            prepare_config = (output / "prepare-normalization.cfg").read_text(
                encoding="utf-8"
            )
            evaluation_config = (output / "reference.cfg").read_text(
                encoding="utf-8"
            )

            self.assertEqual(
                metadata["particle_order"], ["beam", "recoil", "ks1", "ks2"]
            )
            self.assertEqual(
                metadata["data"]["event_ids"], [1001, 1002, 1003, 1004]
            )
            self.assertEqual(metadata["data"]["event_count"], 4)
            self.assertEqual(metadata["normalization"]["event_count"], 6)
            self.assertAlmostEqual(
                metadata["normalization"]["weight_sum"], 6.0
            )
            self.assertEqual(metadata["beam_energy_gev"], 8.5)
            self.assertEqual(metadata["polarization_fraction"], 0.0)
            self.assertTrue((output / "data.root").is_file())
            self.assertTrue((output / "normalization.root").is_file())
            self.assertIn("normintfile reference", prepare_config)
            self.assertNotIn(" input", prepare_config)
            self.assertIn("normalization.ni input", evaluation_config)
            self.assertIn(
                "amplitude reference::PositiveReal::S0+ Zlm 0 0 1 1",
                evaluation_config,
            )
            self.assertIn(
                "initialize reference::PositiveReal::S0+ cartesian 1 0 fixed",
                evaluation_config,
            )

    def test_prepare_is_reproducible(self) -> None:
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            first_path = Path(first)
            second_path = Path(second)
            for output in (first_path, second_path):
                subprocess.run(
                    [
                        sys.executable,
                        str(SCRIPT),
                        "prepare",
                        "--sample",
                        str(REFERENCE / "events.json"),
                        "--output",
                        str(output),
                    ],
                    check=True,
                    cwd=COMPARISON,
                    text=True,
                    capture_output=True,
                )

            first_metadata = (first_path / "input-metadata.json").read_bytes()
            second_metadata = (
                second_path / "input-metadata.json"
            ).read_bytes()
            first_digest = hashlib.sha256(first_metadata).hexdigest()
            second_digest = hashlib.sha256(second_metadata).hexdigest()

        self.assertEqual(first_digest, second_digest)

    def test_validate_accepts_the_known_unpolarized_s_wave(self) -> None:
        result = self.valid_result()
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "result.json"
            path.write_text(json.dumps(result), encoding="utf-8")
            completed = subprocess.run(
                [sys.executable, str(SCRIPT), "validate", "--result", str(path)],
                check=True,
                cwd=COMPARISON,
                text=True,
                capture_output=True,
            )

        self.assertIn(
            "validated AmpTools unpolarized S-wave reference", completed.stdout
        )

    def test_validate_rejects_a_nonconstant_s_wave_intensity(self) -> None:
        result = self.valid_result()
        result["events"][-1]["intensity"] += 1e-4
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "result.json"
            path.write_text(json.dumps(result), encoding="utf-8")
            completed = subprocess.run(
                [sys.executable, str(SCRIPT), "validate", "--result", str(path)],
                check=False,
                cwd=COMPARISON,
                text=True,
                capture_output=True,
            )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("intensity", completed.stderr)

    def test_validate_rejects_a_stale_sample_digest(self) -> None:
        result = self.valid_result()
        result["sample_sha256"] = "0" * 64
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "result.json"
            path.write_text(json.dumps(result), encoding="utf-8")
            completed = subprocess.run(
                [sys.executable, str(SCRIPT), "validate", "--result", str(path)],
                check=False,
                cwd=COMPARISON,
                text=True,
                capture_output=True,
            )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("sample digest", completed.stderr)

    def test_validate_rejects_a_changed_event_weight(self) -> None:
        result = self.valid_result()
        result["events"][0]["weight"] = 2.0
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "result.json"
            path.write_text(json.dumps(result), encoding="utf-8")
            completed = subprocess.run(
                [sys.executable, str(SCRIPT), "validate", "--result", str(path)],
                check=False,
                cwd=COMPARISON,
                text=True,
                capture_output=True,
            )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("event 1001 weight", completed.stderr)

    def test_validate_rejects_changed_particle_ordering(self) -> None:
        result = self.valid_result()
        result["events"][0]["particles"][1:3] = reversed(
            result["events"][0]["particles"][1:3]
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "result.json"
            path.write_text(json.dumps(result), encoding="utf-8")
            completed = subprocess.run(
                [sys.executable, str(SCRIPT), "validate", "--result", str(path)],
                check=False,
                cwd=COMPARISON,
                text=True,
                capture_output=True,
            )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("momentum or ordering changed", completed.stderr)

    def test_validate_rejects_compensating_normalization_weight_changes(
        self,
    ) -> None:
        result = self.valid_result()
        events = result["normalization_source"]["generated"]["events"]
        events[0]["weight"] += 0.25
        events[1]["weight"] -= 0.25
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "result.json"
            path.write_text(json.dumps(result), encoding="utf-8")
            completed = subprocess.run(
                [sys.executable, str(SCRIPT), "validate", "--result", str(path)],
                check=False,
                cwd=COMPARISON,
                text=True,
                capture_output=True,
            )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn(
            "normalization generated event 2001 weight", completed.stderr
        )


if __name__ == "__main__":
    unittest.main()
