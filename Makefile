PYTHON ?= python3
OUTPUT ?= outputs
PACKAGE ?= dist/ecd-snipr-harness-$(shell date -u +%Y%m%dT%H%M%S).zip

.PHONY: validate smoke test privacy-check package all-checks-offline install-user
validate:
	$(PYTHON) scripts/manage_skill.py validate
smoke:
	$(PYTHON) scripts/ecd_snipr_cli.py smoke --outdir $(OUTPUT)/fixture --resume
test:
	$(PYTHON) -m unittest discover -s tests -p 'test_*.py' -v
privacy-check:
	$(PYTHON) scripts/manage_skill.py privacy-check
package:
	$(PYTHON) scripts/manage_skill.py package --output "$(PACKAGE)"
all-checks-offline: validate smoke test privacy-check package
install-user:
	$(PYTHON) scripts/manage_skill.py install --destination "$(HOME)/.agents/skills/ecd-snipr-harness" --compatibility "$(HOME)/.codex/skills/ecd-snipr-harness"
