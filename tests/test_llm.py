import json
from unittest.mock import patch
from urllib.error import HTTPError

import pytest

from polymer_scatt import llm


def example():
    return {"n_curves": 4639, "curves": [
        {"label": "T0", "status": "rejected", "q_unit": "A^-1",
         "candidate": {"rg_nm": 7.89, "qmax_rg": 3.09},
         "reasons": ["qmax*Rg exceeds the screening limit"],
         "q": [1.0, 2.0], "intensity": [9912.34, 8192.12], "sigma": [1, 2]},
    ], "api_key": "NEVER_SEND_THIS", "raw_data": "RAW_PRIVATE_DATA"}


def test_offline_never_contacts_endpoint():
    with patch.object(llm, "_post_json", side_effect=AssertionError("network")):
        result = llm.interpret(example(), endpoint="https://example.org/v1")
    assert result["executed"] is False
    assert result["model"] is None
    assert "no external or local language model was run" in result["narrative"]
    assert "T0: rejected" in result["narrative"]


def test_summary_redacts_raw_vectors_and_credentials():
    prompt = llm.build_prompt(example())
    serialized = json.dumps(prompt)
    assert "NEVER_SEND_THIS" not in serialized
    assert "RAW_PRIVATE_DATA" not in serialized
    assert "9912.34" not in serialized
    assert '"q":' not in serialized
    assert "7.89" in serialized
    assert "3.09" in serialized
    assert "MUST remain unresolved/rejected" in serialized


def test_large_series_is_explicitly_truncated():
    compact = llm.compact_summary({"n_curves": 1000, "curves": [example()["curves"][0]] * 1000})
    assert len(compact["curves"]) == 24
    assert compact["curves_omitted_records"] == 976
    assert compact["n_curves"] == 1000


def test_cli_summary_counts_and_group_outcomes_are_preserved():
    summary = {
        "title": "PolymerScattLab sample analysis", "n_curves": 4639,
        "n_accepted": 0, "n_rejected": 4639, "n_invalid": 0,
        "groups": [{"label": "MAXS1", "n_curves": 4, "n_accepted": 0,
                    "n_rejected": 4, "n_invalid": 0, "q_unit": "A^-1",
                    "rejection_reasons": ["outside_guinier_qrg_range"]}],
        "sample_state": "unknown", "chain_rg_status": "not established",
        "warnings": ["Archive q units unknown"],
    }
    result = llm.interpret(summary)
    assert result["prompt"]["compact_summary"] == summary
    narrative = result["narrative"]
    assert "4639 curves; 0 accepted apparent Guinier fits; 4639 rejected fits; 0 invalid curves" in narrative
    assert "MAXS1: 4 curves; 0 accepted apparent Guinier fits; 4 rejected fits; 0 invalid curves" in narrative
    assert "Sample state: unknown" in narrative
    assert "polymer-chain Rg status: not established" in narrative
    assert "MAXS1: unresolved" not in narrative


def test_ollama_mock_contract_and_no_result_mutation():
    data = example()
    with patch.object(llm, "_post_json", return_value={"message": {"content": "Rejected fit; no valid Rg."}}) as post:
        result = llm.interpret(data, provider="ollama", model="local-model")
    endpoint, payload, key = post.call_args.args
    assert endpoint == "http://localhost:11434/api/chat"
    assert payload["stream"] is False
    assert payload["model"] == "local-model"
    assert result["executed"] is True
    assert data["curves"][0]["status"] == "rejected"
    assert result["numerical_results_modified"] is False
    assert "9912.34" not in json.dumps(payload)


def test_compatible_mock_contract(monkeypatch):
    monkeypatch.setenv("POLYMERSCATT_API_KEY", "test-only-key")
    with patch.object(llm, "_post_json", return_value={"choices": [{"message": {"content": "Unknown units."}}]}) as post:
        result = llm.interpret(example(), provider="openai-compatible", model="configured-model", endpoint="https://example.org/v1")
    endpoint, payload, key = post.call_args.args
    assert endpoint == "https://example.org/v1/chat/completions"
    assert key == "test-only-key"
    assert "test-only-key" not in json.dumps(payload)
    assert "test-only-key" not in json.dumps(result)
    assert result["model"] == "configured-model"


@pytest.mark.parametrize("endpoint", ["http://example.org/v1", "https://user:secret@example.org/v1", "https://example.org/v1?key=secret", "file:///tmp/data"])
def test_rejects_insecure_or_embedded_secret_endpoints(endpoint):
    with pytest.raises(llm.LLMError) as exc:
        llm.interpret(example(), provider="openai", model="m", endpoint=endpoint)
    assert "secret" not in str(exc.value)


def test_configuration_is_explicit(monkeypatch):
    monkeypatch.delenv("POLYMERSCATT_LLM_MODEL", raising=False)
    monkeypatch.delenv("POLYMERSCATT_LLM_ENDPOINT", raising=False)
    with pytest.raises(llm.LLMError, match="explicit model"):
        llm.interpret(example(), provider="ollama")
    with pytest.raises(llm.LLMError, match="explicit OpenAI-compatible endpoint"):
        llm.interpret(example(), provider="openai", model="m")


def test_http_errors_do_not_leak_keys_or_provider_body():
    class FailingOpener:
        def open(self, req, timeout):
            raise HTTPError("https://secret-key@example.org", 401, "SECRET_PROVIDER_BODY", {}, None)
    with patch.object(llm.request, "build_opener", return_value=FailingOpener()):
        with pytest.raises(llm.LLMError) as exc:
            llm._post_json("https://example.org/v1/chat/completions", {}, "secret-key")
    assert "401" in str(exc.value)
    assert "secret-key" not in str(exc.value)
    assert "SECRET_PROVIDER_BODY" not in str(exc.value)


def test_malformed_model_response_has_controlled_error():
    with patch.object(llm, "_post_json", return_value={"choices": []}):
        with pytest.raises(llm.LLMError, match="assistant text"):
            llm.interpret(example(), provider="openai", model="m", endpoint="https://example.org/v1")
