import tempfile
import unittest
import zipfile
from pathlib import Path
from test_design import ROOT
from manage_skill import assets, forbidden, install, package, validate


class ReleaseTests(unittest.TestCase):
    def test_manifest_and_cli_validation(self):
        self.assertTrue(validate(ROOT)["passed"])

    def test_package_members_clean(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "skill.zip"
            result = package(ROOT, path)
            self.assertGreater(result["files"], 20)
            with zipfile.ZipFile(path) as z:
                self.assertIsNone(z.testzip())
                self.assertFalse(any(forbidden(Path(n)) for n in z.namelist()))
                self.assertIn("ecd-snipr-harness/SKILL.md", z.namelist())

    def test_private_paths_forbidden(self):
        for name in ("local_data/a.json", "private_output/results.tsv", "examples/a.xlsx", "dist/release.zip", "cache/uniprot.json", "__pycache__/a.pyc", "real_data_output_installed/table.tsv"):
            self.assertTrue(forbidden(Path(name)), name)

    def test_install_backup_outside_skill_discovery(self):
        with tempfile.TemporaryDirectory() as temp:
            target = Path(temp) / "skills" / "ecd-snipr-harness"
            first = install(ROOT, target)
            self.assertIsNone(first["backup_path"])
            second = install(ROOT, target)
            self.assertEqual(Path(second["backup_path"]).parent, Path(temp) / "skill-backups")
            self.assertTrue((Path(second["backup_path"]) / "SKILL.md").exists())
            self.assertFalse((target / "outputs").exists())


if __name__ == "__main__":
    unittest.main()
