from lsdtt_qgis.naming import extract_json_payload

BANNER_AND_JSON = """\
=========================================================
|| Welcome to lsdtt-channel-extraction                  ||
|| Part of lsdtt3 by Simon M. Mudd                      ||
=========================================================
|| Version : 0.2.0
=========================================================

{
  "tool": "lsdtt-channel-extraction",
  "version": "0.2.0",
  "citation": "https://example.org",
  "features": ["fftw", "eigen"],
  "parameters": [
    {"name": "n_threads", "type": "int", "default": "0",
     "description": "Number of OpenMP threads.", "guidance": "",
     "group": "Run control", "role": "io", "parent": ""}
  ]
}
"""


def test_extract_json_skips_banner():
    payload = extract_json_payload(BANNER_AND_JSON)
    assert payload is not None
    assert payload["tool"] == "lsdtt-channel-extraction"
    assert payload["features"] == ["fftw", "eigen"]
    assert payload["parameters"][0]["name"] == "n_threads"


def test_extract_json_returns_none_without_json():
    assert extract_json_payload("Usage: tool <parameter_file>\n") is None
    assert extract_json_payload("") is None


def test_extract_json_returns_none_for_broken_json():
    assert extract_json_payload("banner\n{ not json at all") is None
