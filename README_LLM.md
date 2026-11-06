# Optional LLM interpretation

The numerical fitting engine determines all parameters and rejection decisions.
The LLM explains a compact summary of those fixed results; its output cannot
change fitted values or validation gates. This is an inference integration,
not a model trained on the uploaded scattering data.

The default mode is **offline**. It creates a reviewable structured prompt and
a deterministic interpretation template. Its output explicitly records
`executed: false`; no language model is called. The sample analysis is usable
without credentials, an account, or an LLM server.

From the project folder, the command-line entry point accepts a generated
`summary.json`:

```bash
python -m polymer_scatt explain analysis_output/summary.json --provider offline
python -m polymer_scatt explain analysis_output/summary.json --provider ollama --model YOUR_INSTALLED_MODEL
```

The local browser interface (`python -m polymer_scatt serve`) also offers an
explicit interpretation action after a curve has been fitted. Local Ollama
must already be running and the chosen model must already be installed.

```python
from polymer_scatt.llm import interpret

# Fit summaries only, no raw curve arrays needed.
offline = interpret(summary)
print(offline["narrative"])
print(offline["prompt"]["messages"])  # inspect the exact proposed prompt

# Explicitly call an already-installed local Ollama model.
local = interpret(summary, provider="ollama", model="YOUR_INSTALLED_MODEL")

# Explicitly call a configured OpenAI-compatible endpoint.
remote = interpret(
    summary,
    provider="openai-compatible",
    model="YOUR_MODEL_NAME",
    endpoint="https://YOUR_PROVIDER/v1/chat/completions",
)
```

Set credentials in `POLYMERSCATT_API_KEY` (or `OPENAI_API_KEY` for compatible
mode), never in a data file or endpoint URL. Optional environment settings are
`POLYMERSCATT_LLM_MODEL` and `POLYMERSCATT_LLM_ENDPOINT`. Selecting a provider
is always explicit. No model is downloaded or installed by this extension.
Compatible providers must support the Chat Completions messages format.

Only allowlisted summary fields, compact sample metadata and the question are
sent. Raw q, intensity and error vectors, arbitrary file contents and API keys
are excluded from the prompt. Lists are capped at 24 representative records,
with an explicit omitted-record count; pass full-set counts separately.
Review the returned prompt before choosing a remote endpoint if sample labels
or metadata are confidential. Online use is subject to the selected provider's
data handling and costs. Remote requests require HTTPS, while a local server
may use HTTP on localhost. Redirects are refused.

Online output records `executed: true`, provider and requested model. Prose is
labelled for review and is not parsed back into numerical results. A model can
still make an incorrect statement; the deterministic report remains the
authoritative source for Rg and validity. The curated references are method
reading, not a live literature search.

API contracts were checked against official documentation:

- https://docs.ollama.com/api/chat
- https://developers.openai.com/api/reference/resources/chat/subresources/completions/methods/create

Mock transport tests exercise both request/response formats, offline operation,
data minimization and sanitized errors. No live model response is claimed.
