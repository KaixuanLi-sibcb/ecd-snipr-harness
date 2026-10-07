PYTHON ?= python3
OUTPUT ?= outputs
PACKAGE ?= dist/ecd-snipr-harness-$(shell date -u +%Y%m%dT%H%M%S).zip

.PHONY: validate smoke test privacy-check package all-checks-offline install-user lab-evidence-fixture receiver-risk-fixture core-evidence-fixture local-predictors-fixture hpc-fixture
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
lab-evidence-fixture:
	$(PYTHON) scripts/ecd_snipr_cli.py lab-evidence --input examples/synthetic_lab.tsv --outdir $(OUTPUT)/lab-fixture --resume
receiver-risk-fixture:
	$(PYTHON) -m unittest discover -s tests -p test_receiver_risk.py -v
core-evidence-fixture:
	$(PYTHON) -m unittest discover -s tests -p test_core_evidence.py -v

local-predictors-fixture:
	$(PYTHON) -m unittest discover -s tests -p test_local_predictors.py -v

hpc-fixture:
	$(PYTHON) -m unittest discover -s tests -p test_hpc_workflow.py -v

all-checks-offline: validate smoke test lab-evidence-fixture receiver-risk-fixture core-evidence-fixture local-predictors-fixture hpc-fixture privacy-check package
install-user:
	$(PYTHON) scripts/manage_skill.py install --destination "$(HOME)/.agents/skills/ecd-snipr-harness" --compatibility "$(HOME)/.codex/skills/ecd-snipr-harness"
