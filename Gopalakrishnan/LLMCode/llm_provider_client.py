#!/usr/bin/env python3
import json
import os
import re
import requests

try:
    from openai import OpenAI
except Exception:
    OpenAI = None

try:
    from google import genai
except Exception:
    genai = None


PROVIDER_SIDECOV_DIRS = {
    "openai": "openAIGenerated",
    "claude": "claudeGenerated",
    "gemini": "geminiGenerated",
}

DEFAULT_MODELS = {
    "openai": "gpt-5.4",
    "claude": "claude-opus-4-6",
    "gemini": "gemini-3.1-pro-preview",
}

PROXY_DEFAULT_MODELS = {
    "openai": "openai.gpt-5.4",
    "claude": "anthropic.claude-4.6-opus",
    "gemini": "google.gemini-3.1-pro-preview",
}

MODEL_PRICING = {
    "openai": {
        "gpt-5.4-pro": {"input": 30.00, "cached_input": None, "output": 180.00, "long_context_threshold": 272000},
        "gpt-5.4": {"input": 2.50, "cached_input": 0.25, "output": 15.00, "long_context_threshold": 272000},
        "gpt-5.2-chat-latest": {"input": 1.75, "cached_input": 0.175, "output": 14.00, "long_context_threshold": None},
        "gpt-5.2": {"input": 1.75, "cached_input": 0.175, "output": 14.00, "long_context_threshold": None},
        "gpt-5.1-chat-latest": {"input": 1.25, "cached_input": 0.125, "output": 10.00, "long_context_threshold": None},
        "gpt-5.1": {"input": 1.25, "cached_input": 0.125, "output": 10.00, "long_context_threshold": None},
        "gpt-5-chat-latest": {"input": 1.25, "cached_input": 0.125, "output": 10.00, "long_context_threshold": None},
        "gpt-5": {"input": 1.25, "cached_input": 0.125, "output": 10.00, "long_context_threshold": None},
    },
    "claude": {
        "claude-opus-4-1": {"input": 15.00, "cached_input": 1.50, "output": 75.00, "long_context_threshold": None},
        "claude-opus-4": {"input": 15.00, "cached_input": 1.50, "output": 75.00, "long_context_threshold": None},
        "claude-sonnet-4": {"input": 3.00, "cached_input": 0.30, "output": 15.00, "long_context_threshold": None},
    },
    "gemini": {
        "gemini-2.5-pro": {"input": 1.25, "cached_input": 0.125, "output": 10.00, "long_context_threshold": 200000},
        "gemini-2.5-flash": {"input": 0.30, "cached_input": 0.03, "output": 2.50, "long_context_threshold": 200000},
    },
}

_openai_client = None
_gemini_client = None


def _request_timeout_seconds(default_seconds=300):
    raw = os.getenv("LLM_REQUEST_TIMEOUT_SECONDS", os.getenv("OPENAI_TIMEOUT_SECONDS", str(default_seconds)))
    try:
        return float(raw)
    except (TypeError, ValueError):
        return float(default_seconds)


def normalize_provider(provider=None):
    provider = (provider or os.getenv("LLM_PROVIDER", "openai")).strip().lower()
    aliases = {
        "anthropic": "claude",
        "claude": "claude",
        "gemini": "gemini",
        "google": "gemini",
        "openai": "openai",
    }
    if provider not in aliases:
        raise ValueError(f"Unsupported LLM provider: {provider}")
    return aliases[provider]


def get_proxy_api_key():
    return (
        os.getenv("OPENAI_API_KEY")
        or os.getenv("API_KEY")
        or os.getenv("CORNELL_API_KEY")
        or os.getenv("Cornell_API_Key")
    )


def get_provider_base_url(provider=None):
    raw = (
        os.getenv("OPENAI_BASE_URL")
        or os.getenv("API_BASE_URI")
        or os.getenv("API_BASE_URL")
        or os.getenv("API_Base_URI")
        or os.getenv("API_Base_URL")
        or os.getenv("CORNELL_API_BASE_URL")
        or os.getenv("Cornell_API_Base_URL")
    )
    if raw:
        return raw.rstrip("/")
    return None


def using_openai_compatible_proxy(provider=None):
    normalize_provider(provider)
    return bool(get_provider_base_url(provider) and get_proxy_api_key())


def get_provider_model(provider=None):
    provider = normalize_provider(provider)
    explicit = os.getenv("LLM_MODEL", "").strip()
    if explicit:
        return explicit

    if provider == "openai":
        default_model = PROXY_DEFAULT_MODELS[provider] if using_openai_compatible_proxy(provider) else DEFAULT_MODELS[provider]
        return os.getenv("OPENAI_MODEL", default_model)
    if provider == "claude":
        default_model = PROXY_DEFAULT_MODELS[provider] if using_openai_compatible_proxy(provider) else DEFAULT_MODELS[provider]
        return os.getenv("ANTHROPIC_MODEL", os.getenv("CLAUDE_MODEL", default_model))
    if provider == "gemini":
        default_model = PROXY_DEFAULT_MODELS[provider] if using_openai_compatible_proxy(provider) else DEFAULT_MODELS[provider]
        return os.getenv("GEMINI_MODEL", default_model)
    raise ValueError(f"Unsupported LLM provider: {provider}")


def get_provider_api_key(provider=None):
    provider = normalize_provider(provider)
    if using_openai_compatible_proxy(provider):
        return get_proxy_api_key()
    if provider == "openai":
        return get_proxy_api_key()
    if provider == "claude":
        return os.getenv("ANTHROPIC_API_KEY")
    if provider == "gemini":
        return os.getenv("GEMINI_API_KEY")
    return None


def get_provider_sidecov_dir(provider=None):
    provider = normalize_provider(provider)
    return os.getenv("SIDECOV_PROVIDER_DIR", PROVIDER_SIDECOV_DIRS[provider])


def get_provider_file_tag(provider=None):
    provider = normalize_provider(provider)
    explicit = os.getenv("LLM_FILE_TAG", "").strip()
    if explicit:
        return explicit
    return provider


def with_provider_tag(path, provider=None):
    tag = get_provider_file_tag(provider)
    if not tag:
        return path
    root, ext = os.path.splitext(path)
    return f"{root}_{tag}{ext}"


def get_script_dir():
    return os.path.dirname(os.path.abspath(__file__))


def get_runtime_dir():
    runtime_dir = os.getenv("LLM_RUNTIME_DIR", os.path.join(get_script_dir(), ".runtime"))
    os.makedirs(runtime_dir, exist_ok=True)
    return runtime_dir


def runtime_path(*parts, is_dir=False):
    path = os.path.join(get_runtime_dir(), *parts) if parts else get_runtime_dir()
    target = path if is_dir else os.path.dirname(path)
    if target:
        os.makedirs(target, exist_ok=True)
    return path


def ensure_parent_dir(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    return path


def assert_provider_ready(provider=None):
    provider = normalize_provider(provider)
    if using_openai_compatible_proxy(provider):
        if not get_proxy_api_key():
            raise RuntimeError("Missing Cornell/OpenAI-compatible proxy API key.")
        return

    api_key = get_provider_api_key(provider)
    if not api_key:
        raise RuntimeError(f"Missing API key for provider '{provider}'")
    if provider == "openai" and OpenAI is None:
        raise RuntimeError("Python package `openai` is not installed in this environment.")
    if provider == "gemini" and genai is None:
        raise RuntimeError("Python package `google-genai` is not installed in this environment.")


def usage_get(obj, key, default=0):
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def pricing_for_model(provider=None, model_name=None):
    provider = normalize_provider(provider)
    model_name = (model_name or get_provider_model(provider)).lower()
    pricing_table = MODEL_PRICING.get(provider, {})

    candidates = [model_name]
    if "." in model_name:
        candidates.append(model_name.split(".", 1)[1])

    for candidate in candidates:
        for key in sorted(pricing_table, key=len, reverse=True):
            if (
                candidate == key
                or candidate.startswith(f"{key}-")
                or candidate.startswith(f"{key}.")
            ):
                return dict(pricing_table[key], known=True)

    return {
        "input": 0.0,
        "cached_input": 0.0,
        "output": 0.0,
        "long_context_threshold": None,
        "known": False,
    }


def summarize_usage(usage, provider=None, model_name=None):
    provider = normalize_provider(provider)
    model_name = model_name or get_provider_model(provider)

    if provider == "openai":
        prompt_tokens = int(
            usage_get(usage, "prompt_tokens", usage_get(usage, "input_tokens", 0)) or 0
        )
        completion_tokens = int(
            usage_get(usage, "completion_tokens", usage_get(usage, "output_tokens", 0)) or 0
        )
        details = usage_get(
            usage,
            "input_tokens_details",
            usage_get(usage, "prompt_tokens_details", None),
        )
        cached_prompt_tokens = int(usage_get(details, "cached_tokens", 0) or 0)
        billable_prompt_tokens = max(prompt_tokens - cached_prompt_tokens, 0)
    elif provider == "claude":
        input_tokens = int(usage_get(usage, "input_tokens", 0) or 0)
        cache_creation_tokens = int(usage_get(usage, "cache_creation_input_tokens", 0) or 0)
        cache_read_tokens = int(usage_get(usage, "cache_read_input_tokens", 0) or 0)
        prompt_tokens = input_tokens + cache_creation_tokens + cache_read_tokens
        cached_prompt_tokens = cache_read_tokens
        completion_tokens = int(usage_get(usage, "output_tokens", 0) or 0)
        billable_prompt_tokens = input_tokens + cache_creation_tokens
    elif provider == "gemini":
        prompt_tokens = int(usage_get(usage, "prompt_token_count", 0) or 0)
        cached_prompt_tokens = int(usage_get(usage, "cached_content_token_count", 0) or 0)
        completion_tokens = int(
            (usage_get(usage, "candidates_token_count", 0) or 0) +
            (usage_get(usage, "thoughts_token_count", 0) or 0)
        )
        billable_prompt_tokens = max(prompt_tokens - cached_prompt_tokens, 0)
    else:
        prompt_tokens = 0
        cached_prompt_tokens = 0
        completion_tokens = 0
        billable_prompt_tokens = 0

    pricing = pricing_for_model(provider, model_name)
    input_rate = pricing["input"]
    cached_input_rate = pricing["cached_input"]
    output_rate = pricing["output"]
    long_context_applied = False
    threshold = pricing.get("long_context_threshold")

    if threshold and prompt_tokens > threshold:
        input_rate *= 2.0
        output_rate *= 1.5
        if cached_input_rate is not None:
            cached_input_rate *= 2.0
        long_context_applied = True

    estimated_cost_usd = (billable_prompt_tokens / 1_000_000.0) * input_rate
    if cached_input_rate is not None and cached_prompt_tokens:
        estimated_cost_usd += (cached_prompt_tokens / 1_000_000.0) * cached_input_rate
    estimated_cost_usd += (completion_tokens / 1_000_000.0) * output_rate

    return {
        "prompt_tokens": prompt_tokens,
        "cached_prompt_tokens": cached_prompt_tokens,
        "completion_tokens": completion_tokens,
        "estimated_cost_usd": estimated_cost_usd,
        "long_context_applied": long_context_applied,
        "pricing_known": pricing.get("known", False),
    }


def _parse_json_text(text):
    if not text:
        return {}
    try:
        return json.loads(text)
    except Exception:
        pass

    text2 = re.sub(r"^```(?:json)?", "", text.strip(), flags=re.IGNORECASE)
    text2 = re.sub(r"```$", "", text2.strip())
    try:
        return json.loads(text2)
    except Exception:
        return {}


def _get_openai_client():
    global _openai_client
    if _openai_client is None:
        assert_provider_ready("openai")
        kwargs = {"api_key": get_provider_api_key("openai")}
        base_url = get_provider_base_url("openai")
        if base_url:
            kwargs["base_url"] = base_url
        _openai_client = OpenAI(**kwargs)
    return _openai_client


def _get_gemini_client():
    global _gemini_client
    if _gemini_client is None:
        assert_provider_ready("gemini")
        _gemini_client = genai.Client(api_key=get_provider_api_key("gemini"))
    return _gemini_client


def _extract_proxy_choice_text(data):
    choices = data.get("choices") or []
    if not choices:
        return ""
    message = (choices[0] or {}).get("message") or {}
    content = message.get("content")
    if isinstance(content, str):
        return content.strip()
    if isinstance(content, list):
        texts = []
        for item in content:
            if isinstance(item, dict) and item.get("type") in {"text", "output_text"}:
                texts.append(item.get("text", ""))
        return "".join(texts).strip()
    return ""


def _call_proxy_json(system_prompt, user_prompt, model_name, reasoning_effort):
    base_url = get_provider_base_url("openai")
    headers = {
        "x-litellm-api-key": get_proxy_api_key(),
        "content-type": "application/json",
    }
    payload = {
        "model": model_name,
        "temperature": 0,
        "max_tokens": int(
            os.getenv(
                "LLM_PROXY_MAX_TOKENS",
                os.getenv(
                    "ANTHROPIC_MAX_TOKENS",
                    os.getenv("GEMINI_MAX_OUTPUT_TOKENS", "4096"),
                ),
            )
        ),
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
    }
    if not model_name.lower().startswith("anthropic."):
        payload["response_format"] = {"type": "json_object"}
    if model_name.lower().startswith("openai."):
        payload["reasoning_effort"] = reasoning_effort

    resp = requests.post(
        f"{base_url}/v1/chat/completions",
        headers=headers,
        json=payload,
        timeout=_request_timeout_seconds(),
    )
    resp.raise_for_status()
    data = resp.json()
    text = _extract_proxy_choice_text(data)
    return _parse_json_text(text), data.get("usage", {})


def _call_openai_json(system_prompt, user_prompt, model_name, reasoning_effort):
    if using_openai_compatible_proxy("openai"):
        return _call_proxy_json(system_prompt, user_prompt, model_name, reasoning_effort)

    client = _get_openai_client()
    timeout_seconds = _request_timeout_seconds()
    try:
        resp = client.responses.create(
            model=model_name,
            reasoning={"effort": reasoning_effort},
            input=[
                {"role": "system", "content": [{"type": "input_text", "text": system_prompt}]},
                {"role": "user", "content": [{"type": "input_text", "text": user_prompt}]},
            ],
            text={"format": {"type": "json_object"}},
            timeout=timeout_seconds,
        )
        text = getattr(resp, "output_text", "") or ""
        return _parse_json_text(text), getattr(resp, "usage", None)
    except Exception:
        resp = client.chat.completions.create(
            model=model_name,
            temperature=0,
            reasoning_effort=reasoning_effort,
            response_format={"type": "json_object"},
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt},
            ],
            timeout=timeout_seconds,
        )
        text = (resp.choices[0].message.content or "").strip()
        return _parse_json_text(text), getattr(resp, "usage", None)


def _call_claude_json(system_prompt, user_prompt, model_name):
    if using_openai_compatible_proxy("claude"):
        return _call_proxy_json(system_prompt, user_prompt, model_name, reasoning_effort="low")

    assert_provider_ready("claude")
    headers = {
        "x-api-key": get_provider_api_key("claude"),
        "anthropic-version": "2023-06-01",
        "content-type": "application/json",
    }
    payload = {
        "model": model_name,
        "system": system_prompt,
        "messages": [{"role": "user", "content": user_prompt}],
        "temperature": 0,
        "max_tokens": int(os.getenv("ANTHROPIC_MAX_TOKENS", "4096")),
    }
    resp = requests.post(
        "https://api.anthropic.com/v1/messages",
        headers=headers,
        json=payload,
        timeout=300,
    )
    resp.raise_for_status()
    data = resp.json()
    text = "".join(
        block.get("text", "")
        for block in (data.get("content") or [])
        if isinstance(block, dict) and block.get("type") == "text"
    )
    return _parse_json_text(text), data.get("usage", {})


def _call_gemini_json(system_prompt, user_prompt, model_name):
    if using_openai_compatible_proxy("gemini"):
        return _call_proxy_json(system_prompt, user_prompt, model_name, reasoning_effort="low")

    client = _get_gemini_client()
    config = genai.types.GenerateContentConfig(
        system_instruction=system_prompt,
        temperature=0,
        max_output_tokens=int(os.getenv("GEMINI_MAX_OUTPUT_TOKENS", "4096")),
        response_mime_type="application/json",
    )
    resp = client.models.generate_content(
        model=model_name,
        contents=user_prompt,
        config=config,
    )
    text = getattr(resp, "text", "") or ""
    return _parse_json_text(text), getattr(resp, "usage_metadata", None)


def call_llm_json(system_prompt, user_prompt, provider=None, reasoning_effort="high"):
    provider = normalize_provider(provider)
    model_name = get_provider_model(provider)

    if provider == "openai":
        parsed, usage = _call_openai_json(system_prompt, user_prompt, model_name, reasoning_effort)
    elif provider == "claude":
        parsed, usage = _call_claude_json(system_prompt, user_prompt, model_name)
    elif provider == "gemini":
        parsed, usage = _call_gemini_json(system_prompt, user_prompt, model_name)
    else:
        raise ValueError(f"Unsupported LLM provider: {provider}")

    return parsed, usage, model_name
