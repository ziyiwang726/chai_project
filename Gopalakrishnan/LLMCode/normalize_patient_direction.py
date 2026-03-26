#!/usr/bin/env python3
import os
import re

import pandas as pd

from llm_provider_client import normalize_provider, with_provider_tag


LLM_PROVIDER = normalize_provider()
FILLED_FILE = os.getenv("LLM_FILLED_FILE", with_provider_tag("taxon_article_matrix_filled.csv"))
OUTPUT_FILE = os.getenv("LLM_OUTPUT_FILE", FILLED_FILE)

CELL_RE = re.compile(r"^\s*([+-])(\s*[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?\s*([ABC])(?:.*)?)$")
BRACKET_RE = re.compile(r"\[(.*)\]\s*$")
HR_RE = re.compile(
    r"\b(?:subdistribution\s+)?(?:cause[- ]specific\s+)?(?:hazard\s+ratio|hr|shr|cshr|cs[- ]?hr)\s*[=:]?\s*([<>]?\s*-?\d*\.?\d+(?:e[-+]?\d+)?)",
    flags=re.IGNORECASE,
)
COEF_RE = re.compile(
    r"\b(?:coef(?:ficient)?|beta|β)\s*[=:]?\s*([-+]?\d*\.?\d+(?:e[-+]?\d+)?)",
    flags=re.IGNORECASE,
)

RESPONSE_POSITIVE = (
    " responder ",
    " responders ",
    " response ",
    " responded ",
    " responsive ",
    " clinical benefit ",
    " disease control ",
    " objective response ",
    " partial response ",
    " complete response ",
    " tumor shrinkage ",
)
RESPONSE_NEGATIVE = (
    " non-responder ",
    " nonresponders ",
    " non-responder ",
    " non-responders ",
    " nonresponder ",
    " nonresponders ",
    " no response ",
    " lack of response ",
    " resistance ",
    " resistant ",
    " progressive disease ",
    " progression ",
    " worse response ",
    " poor response ",
)
SURVIVAL_POSITIVE = (
    " longer survival ",
    " longer os ",
    " longer pfs ",
    " longer dfs ",
    " longer rfs ",
    " longer ttp ",
    " prolonged survival ",
    " prolonged os ",
    " prolonged pfs ",
    " improved survival ",
    " improved os ",
    " improved pfs ",
    " better survival ",
    " survival benefit ",
    " lower hazard ",
    " reduced hazard ",
    " lower risk of progression ",
    " reduced risk of progression ",
    " lower risk of recurrence ",
    " reduced risk of recurrence ",
    " lower risk of death ",
    " reduced risk of death ",
    " without recurrence ",
    " recurrence-free ",
    " disease-free ",
)
SURVIVAL_NEGATIVE = (
    " shorter survival ",
    " shorter os ",
    " shorter pfs ",
    " shorter dfs ",
    " shorter rfs ",
    " reduced pfs ",
    " reduced os ",
    " worse survival ",
    " poorer survival ",
    " poor survival ",
    " higher hazard ",
    " increased hazard ",
    " higher risk of progression ",
    " increased risk of progression ",
    " higher risk of recurrence ",
    " increased risk of recurrence ",
    " higher risk of death ",
    " increased risk of death ",
    " with recurrence ",
    " recurrence ",
)
GENERAL_POSITIVE = (
    " beneficial ",
    " benefit ",
    " favorable ",
    " improved outcome ",
    " better outcome ",
)
GENERAL_NEGATIVE = (
    " harmful ",
    " harm ",
    " unfavorable ",
    " adverse outcome ",
    " worse outcome ",
)


def extract_bracket_text(cell_value):
    match = BRACKET_RE.search(str(cell_value or ""))
    return (match.group(1) or "").strip() if match else ""


def norm_text(text):
    return f" {str(text or '').lower()} "


def parse_numeric(token):
    token = str(token or "").strip()
    token = token.lstrip("<>").strip()
    if not token:
        return None
    try:
        return float(token)
    except ValueError:
        return None


def collect_signal(values, positive_sign="+", negative_sign="-"):
    cleaned = [v for v in values if v is not None]
    if not cleaned:
        return None
    if all(v < 1 for v in cleaned):
        return positive_sign
    if all(v > 1 for v in cleaned):
        return negative_sign
    return None


def collect_coef_signal(values):
    cleaned = [v for v in values if v is not None]
    if not cleaned:
        return None
    if all(v < 0 for v in cleaned):
        return "+"
    if all(v > 0 for v in cleaned):
        return "-"
    return None


def has_any(text, patterns):
    return any(pattern in text for pattern in patterns)


def infer_sign(marker, effect_text):
    text = norm_text(effect_text)
    if not text.strip():
        return None

    signals = []

    if has_any(text, RESPONSE_NEGATIVE):
        signals.append("-")
    if has_any(text, RESPONSE_POSITIVE):
        signals.append("+")

    hr_values = [parse_numeric(token) for token in HR_RE.findall(effect_text)]
    hr_signal = collect_signal(hr_values)
    if hr_signal is not None:
        signals.append(hr_signal)

    hazard_context = any(term in text for term in (" cox ", " fine-gray ", " fine gray ", " cause-specific ", " hazard "))
    coef_values = [parse_numeric(token) for token in COEF_RE.findall(effect_text)]
    coef_signal = collect_coef_signal(coef_values) if hazard_context else None
    if coef_signal is not None:
        signals.append(coef_signal)

    if has_any(text, SURVIVAL_NEGATIVE):
        signals.append("-")
    if has_any(text, SURVIVAL_POSITIVE):
        signals.append("+")

    if marker == "C":
        if has_any(text, GENERAL_NEGATIVE):
            signals.append("-")
        if has_any(text, GENERAL_POSITIVE):
            signals.append("+")

    signals = [signal for signal in signals if signal in {"+", "-"}]
    if not signals:
        return None
    if len(set(signals)) != 1:
        return None
    return signals[0]


def normalize_cell(cell_value):
    match = CELL_RE.match(str(cell_value or ""))
    if not match:
        return cell_value, False

    current_sign = match.group(1)
    remainder = match.group(2)
    marker = match.group(3)
    effect_text = extract_bracket_text(cell_value)
    new_sign = infer_sign(marker, effect_text)

    if new_sign is None or new_sign == current_sign:
        return cell_value, False
    return f"{new_sign}{remainder}", True


def main():
    print(f"[CHECKPOINT] Loading {FILLED_FILE}")
    df = pd.read_csv(FILLED_FILE, dtype=str).fillna("")

    changed = 0
    reviewed = 0
    for col in df.columns[1:]:
        for row_idx in df.index:
            value = df.at[row_idx, col]
            if not str(value).strip():
                continue
            reviewed += 1
            new_value, did_change = normalize_cell(value)
            if did_change:
                df.at[row_idx, col] = new_value
                changed += 1

    df.to_csv(OUTPUT_FILE, index=False)
    print(
        f"[CHECKPOINT] normalize_patient_direction complete for provider={LLM_PROVIDER}: "
        f"reviewed={reviewed}, changed={changed}, saved={OUTPUT_FILE}"
    )


if __name__ == "__main__":
    main()
