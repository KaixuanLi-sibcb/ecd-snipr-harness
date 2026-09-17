#!/usr/bin/env python3
"""Validate, package and install only allowlisted workflow assets."""

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from ecd_snipr.common import file_hash, read_json
from ecd_snipr.contracts import validate_project

FORBIDDEN = {"local_data", "real_data_output", "smoke_test_output", "live_validation_output", "coverage_output", "private_output", "outputs", "cache", ".git", "__pycache__", ".pytest_cache", "dist", "backup"}


def forbidden(path):
    return any(p in FORBIDDEN or p.endswith(".egg-info") or p.startswith(("real_data_output", "smoke_test_output")) for p in path.parts) or path.suffix.lower() in {".xlsx", ".xls", ".zip", ".pyc", ".pem", ".key"} or path.name in {".DS_Store", ".env"}


def assets(root):
    manifest = read_json(root / "manifest.json")
    result = []
    for relative in manifest["package_paths"]:
        path = root / relative
        if not path.exists() or not path.resolve().is_relative_to(root.resolve()):
            raise ValueError("Missing/outside manifest path: " + relative)
        items = sorted(path.rglob("*")) if path.is_dir() else [path]
        for item in items:
            rel = item.relative_to(root)
            if item.is_symlink():
                raise ValueError("Symlinks are not packaged: " + str(rel))
            if item.is_file() and not forbidden(rel):
                result.append(rel)
    return sorted(set(result))


def privacy(root):
    files = assets(root)
    tracked = subprocess.run(["git", "-C", str(root), "ls-files", "-z"], capture_output=True)
    violations = []
    scope = "package_allowlist_only_no_git_repository"
    if tracked.returncode == 0:
        scope = "tracked_paths_and_package_allowlist"
        violations += [p for p in tracked.stdout.decode().split("\0") if p and forbidden(Path(p))]
    secret_pattern = re.compile(r"(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{50,}|-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----)")
    for relative in files:
        text = (root / relative).read_text(encoding="utf-8")
        if secret_pattern.search(text):
            violations.append("credential_pattern:" + str(relative))
    return {"passed": not violations, "scope": scope, "packaged_files": len(files), "violations": violations,
            "limitation": "Path/credential-pattern audit, not a guarantee against arbitrary biological information in text; public release still requires human diff review."}


def validate(root):
    checks = {}
    paths = assets(root)
    skill = (root / "SKILL.md").read_text()
    checks["frontmatter"] = skill.startswith("---\nname: ecd-snipr-harness\ndescription:")
    checks["manifest_paths"] = bool(paths)
    for path in root.joinpath("scripts").glob("*.py"):
        command = subprocess.run([sys.executable, str(path), "--help"], capture_output=True)
        checks[path.name + "_help"] = command.returncode == 0
    validate_project(read_json(root / "examples/synthetic_project.json"))
    checks["fixture_contract"] = True
    checks["privacy"] = privacy(root)["passed"]
    # Check markdown's relative file links, ignoring anchors, URLs and code paths.
    broken = []
    for relative in paths:
        path = root / relative
        if path.suffix == ".md":
            for target in re.findall(r"\]\(([^)]+)\)", path.read_text()):
                if not target.startswith(("https://", "http://", "#")):
                    candidate = path.parent / target.split("#")[0]
                    if not candidate.exists():
                        broken.append(str(relative) + ":" + target)
    checks["relative_document_links"] = not broken
    return {"passed": all(checks.values()), "checks": checks, "broken_links": broken}


def package(root, output):
    if not validate(root)["passed"]:
        raise ValueError("Validation failed")
    output = Path(output)
    if output.exists():
        raise FileExistsError("Package destination exists; choose a new output name")
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED) as archive:
        for relative in assets(root):
            archive.write(root / relative, "ecd-snipr-harness/" + relative.as_posix())
    with zipfile.ZipFile(output) as archive:
        if archive.testzip() or any(forbidden(Path(n)) for n in archive.namelist()):
            raise ValueError("Package privacy/integrity failure")
    return {"path": str(output.resolve()), "sha256": file_hash(output), "files": len(assets(root))}


def install(root, destination, compatibility=None):
    destination = Path(destination).expanduser().absolute()
    if destination.resolve() == root.resolve() or destination.is_symlink():
        raise ValueError("Destination must not be the source or an existing symlink")
    if compatibility:
        compatibility = Path(compatibility).expanduser().absolute()
        if compatibility.exists() and not compatibility.is_symlink():
            raise ValueError("Refusing to replace real compatibility directory")
    if not validate(root)["passed"]:
        raise ValueError("Source validation failed")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".ecd-snipr-install-", dir=destination.parent))
    backup = None
    try:
        for relative in assets(root):
            target = staging / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(root / relative, target)
        if not validate(staging)["passed"]:
            raise ValueError("Staged installation validation failed")
        if destination.exists():
            timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
            backup_root = destination.parent.parent / "skill-backups"
            backup_root.mkdir(parents=True, exist_ok=True)
            backup = backup_root / (destination.name + ".backup-" + timestamp)
            destination.rename(backup)
        staging.rename(destination)
    except Exception:
        if backup and not destination.exists():
            backup.rename(destination)
        if staging.exists():
            shutil.rmtree(staging)
        raise
    if compatibility:
        compatibility.parent.mkdir(parents=True, exist_ok=True)
        if compatibility.is_symlink():
            compatibility.unlink()
        compatibility.symlink_to(destination, target_is_directory=True)
    return {"installed_path": str(destination), "backup_path": str(backup) if backup else None, "compatibility_symlink": str(compatibility) if compatibility else None, "validation": validate(destination)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--skill-dir", default=str(ROOT))
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("validate")
    sub.add_parser("privacy-check")
    cmd = sub.add_parser("package")
    cmd.add_argument("--output", required=True)
    cmd = sub.add_parser("install")
    cmd.add_argument("--destination", required=True)
    cmd.add_argument("--compatibility")
    args = p.parse_args()
    root = Path(args.skill_dir).resolve()
    if args.command == "validate":
        result = validate(root)
    elif args.command == "privacy-check":
        result = privacy(root)
    elif args.command == "package":
        result = package(root, args.output)
    else:
        result = install(root, args.destination, args.compatibility)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result.get("passed", True) else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError) as exc:
        print(json.dumps({"passed": False, "error": str(exc)}), file=sys.stderr)
        raise SystemExit(1)
